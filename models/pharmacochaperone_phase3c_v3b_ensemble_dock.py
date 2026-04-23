#!/usr/bin/env python3
"""
Phase 3c v3b — Ensemble docking of the Phase 3c v3 (fenamic expansion +
Phase 6b reversible covalent warheads) library against K1141 on Phase 5a
MD-sampled receptor conformers.

Identical pipeline to Phase 3c v2b but targets the 12k+ v3 library:
  Stage 1 — all ligands × snap_008 (global pocket rep), exh 8 × 3 modes.
  Stage 2 — top 50 × 5 k-means-selected conformers, exh 16 × 5 modes.
  Ensemble mean ΔG → Kd → f_PC.

Expected runtime: ~3h (2.5h Stage 1 on 12k ligands + 30 min Stage 2).

Survives terminal close via nohup; output JSON updated incrementally
(Stage 1 results checkpointed to JSON before Stage 2 starts).
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/egorlyfar/Brain/research/strc/models")
from pharmacochaperone_phase5b_ensemble_redock import (  # type: ignore
    strip_water_ions_to_pdb,
    pdb_to_pdbqt,
    BOX_CENTRE,
    BOX_SIZE,
    RT_KCAL,
    L_THERAPEUTIC_uM_HIGH,
    RESCUE_ETA_CONSERVATIVE,
)
from pharmacochaperone_phase3c_v2b_ensemble_dock import (  # type: ignore
    pocket_ca_coords,
    kmeans_1d_pocket,
    K1141_POCKET_RESIDUES,
)

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
SNAPSHOT_DIR = WORK / "artifacts" / "phase5a_snapshots"
LIB_JSON = WORK / "pharmacochaperone_phase3c_v3_fenamic_covalent_library.json"
DOCK_DIR = WORK / "docking_runs" / "3c_v3"
STAGE1_DIR = DOCK_DIR / "stage1"
STAGE2_DIR = DOCK_DIR / "stage2"
RECEPTOR_DIR = DOCK_DIR / "receptors"
OUT_JSON = WORK / "pharmacochaperone_phase3c_v3b_ensemble_dock.json"
CHECKPOINT_JSON = WORK / "pharmacochaperone_phase3c_v3b_ensemble_dock_checkpoint.json"

SHORTLIST_N = 50
N_CLUSTERS = 5
STAGE1_EXHAUSTIVENESS = 8
STAGE2_EXHAUSTIVENESS = 16


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


def select_receptors() -> tuple[list[int], list[Path]]:
    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)
    snapshots = sorted(SNAPSHOT_DIR.glob("snap_???.pdb"))
    log(f"Loading pocket Cα from {len(snapshots)} snapshots...")
    frames_pocket = []
    stripped_pdbs: list[Path] = []
    for snap in snapshots:
        prot = RECEPTOR_DIR / f"{snap.stem}_prot.pdb"
        if not prot.exists():
            strip_water_ions_to_pdb(snap, prot)
        stripped_pdbs.append(prot)
        frames_pocket.append(pocket_ca_coords(prot))
    frames_pocket = np.stack(frames_pocket, axis=0)
    mean = frames_pocket.mean(axis=0)
    dists_mean = np.linalg.norm(
        frames_pocket.reshape(len(frames_pocket), -1) - mean.flatten(), axis=1)
    rep_global = int(dists_mean.argmin())
    log(f"Global rep: snap_{rep_global:03d}")
    reps = kmeans_1d_pocket(frames_pocket, k=N_CLUSTERS, seed=42)
    log(f"K-means reps: {[f'snap_{r:03d}' for r in reps]}")
    rep_pdbqts = []
    for r in reps:
        pdbqt = RECEPTOR_DIR / f"snap_{r:03d}_receptor.pdbqt"
        if not pdbqt.exists():
            pdb_to_pdbqt(stripped_pdbs[r], pdbqt)
        rep_pdbqts.append(pdbqt)
    global_pdbqt = RECEPTOR_DIR / f"snap_{rep_global:03d}_receptor.pdbqt"
    if not global_pdbqt.exists():
        pdb_to_pdbqt(stripped_pdbs[rep_global], global_pdbqt)
    if rep_global not in reps:
        reps = [rep_global] + reps[:-1]
        rep_pdbqts = [global_pdbqt] + rep_pdbqts[:-1]
    return reps, rep_pdbqts


def vina_dock_one(receptor: Path, ligand: Path, out_dir: Path,
                  exhaustiveness: int, num_modes: int,
                  timeout: int = 180) -> float | None:
    cx, cy, cz = BOX_CENTRE
    sx, sy, sz = BOX_SIZE
    out_pdbqt = out_dir / f"{ligand.stem}__{receptor.stem}.pdbqt"
    log_txt = out_dir / f"{ligand.stem}__{receptor.stem}.log"
    cmd = [
        "vina", "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--cpu", "8",
        "--out", str(out_pdbqt),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    log_txt.write_text(r.stdout + "\n----STDERR----\n" + r.stderr)
    for line in r.stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def stage1(library: list[dict], rep_pdbqt: Path) -> list[tuple[str, str, float | None, bool]]:
    STAGE1_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    t0 = time.time()
    for i, entry in enumerate(library):
        lig_path = Path(entry["pdbqt_path"])
        if not lig_path.exists():
            continue
        dG = vina_dock_one(rep_pdbqt, lig_path, STAGE1_DIR,
                            STAGE1_EXHAUSTIVENESS, num_modes=3, timeout=120)
        results.append((entry["name"], entry["smiles"], dG,
                        entry.get("is_covalent", False)))
        if (i + 1) % 100 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta_min = (len(library) - i - 1) / rate / 60
            log(f"  stage1 {i+1}/{len(library)} ({rate:.2f} lig/s, ETA {eta_min:.0f} min)")
            # Checkpoint periodically
            CHECKPOINT_JSON.write_text(json.dumps({
                "stage": "1",
                "progress": i + 1,
                "total": len(library),
                "results_so_far": [
                    {"name": n, "smiles": s, "dG": dg, "is_covalent": cov}
                    for n, s, dg, cov in results
                ],
            }, indent=2, default=str))
    results.sort(key=lambda r: (r[2] if r[2] is not None else 0.0))
    return results


def stage2(shortlist: list, library_map: dict, receptors: list[Path]) -> list[dict]:
    STAGE2_DIR.mkdir(parents=True, exist_ok=True)
    final = []
    t0 = time.time()
    for i, (name, smi, dG_s1, is_cov) in enumerate(shortlist):
        lig_path = Path(library_map[name]["pdbqt_path"])
        per_rec = []
        for rec in receptors:
            dG = vina_dock_one(rec, lig_path, STAGE2_DIR,
                                STAGE2_EXHAUSTIVENESS, num_modes=5, timeout=240)
            if dG is not None:
                per_rec.append(dG)
        if not per_rec:
            continue
        arr = np.array(per_rec)
        mean_dG = float(arr.mean())
        std_dG = float(arr.std())
        Kd_uM = math.exp(mean_dG / RT_KCAL) * 1e6
        theta = L_THERAPEUTIC_uM_HIGH / (L_THERAPEUTIC_uM_HIGH + Kd_uM)
        f_PC = theta * RESCUE_ETA_CONSERVATIVE
        final.append({
            "name": name,
            "smiles": smi,
            "is_covalent": is_cov,
            "mean_dG": round(mean_dG, 3),
            "std_dG": round(std_dG, 3),
            "per_receptor_dG": [round(x, 3) for x in per_rec],
            "Kd_uM_from_mean": round(Kd_uM, 2),
            "bound_fraction_at_10uM": round(theta, 4),
            "f_PC_estimate": round(f_PC, 4),
            "MW": library_map[name].get("MW"),
            "logP": library_map[name].get("logP"),
        })
        elapsed = time.time() - t0
        rate = (i + 1) / elapsed
        log(f"  stage2 {i+1}/{len(shortlist)} {name[:45]:45s} "
            f"ΔG={mean_dG:.2f}±{std_dG:.2f} Kd={Kd_uM:.1f} μM f_PC={f_PC:.3f} "
            f"({rate*60:.1f} lig/min)")
    final.sort(key=lambda x: x["mean_dG"])
    return final


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shortlist-n", type=int, default=SHORTLIST_N)
    args = parser.parse_args()

    lib = json.loads(LIB_JSON.read_text())
    library = lib["ligands"]
    library_map = {e["name"]: e for e in library}
    log(f"=== Phase 3c v3b ensemble dock ===")
    log(f"Library size: {len(library)} ligands")

    rep_snaps, rep_pdbqts = select_receptors()
    log(f"Receptors: {[p.stem for p in rep_pdbqts]}")

    log("=== Stage 1 ===")
    s1 = stage1(library, rep_pdbqts[0])
    log(f"Stage 1 complete; top 10 by ΔG:")
    for r in s1[:10]:
        cov = "COV" if r[3] else "   "
        log(f"  {r[2]} kcal/mol {cov}  {r[0]}")

    log(f"=== Stage 2: top {args.shortlist_n} ===")
    shortlist = s1[:args.shortlist_n]
    s2 = stage2(shortlist, library_map, rep_pdbqts)
    log("Stage 2 final ranking:")
    for i, r in enumerate(s2[:20]):
        cov = "COV" if r["is_covalent"] else "   "
        log(f"  #{i+1} {cov} {r['name'][:50]:50s} "
            f"ΔG={r['mean_dG']:.2f}±{r['std_dG']:.2f} "
            f"Kd={r['Kd_uM_from_mean']:.1f} μM "
            f"f_PC={r['f_PC_estimate']:.3f}")

    green = [r for r in s2 if r["f_PC_estimate"] >= 0.50]
    yellow = [r for r in s2 if 0.25 <= r["f_PC_estimate"] < 0.50]
    red = [r for r in s2 if r["f_PC_estimate"] < 0.25]
    cov_green = [r for r in green if r["is_covalent"]]
    nc_green = [r for r in green if not r["is_covalent"]]
    log("=== VERDICT ===")
    log(f"GREEN (f_PC ≥ 0.50): {len(green)} — {len(cov_green)} covalent, {len(nc_green)} non-covalent")
    for r in green[:5]:
        c = "COV" if r["is_covalent"] else "NC "
        log(f"  → {c} {r['name'][:50]:50s} f_PC={r['f_PC_estimate']:.3f} Kd={r['Kd_uM_from_mean']:.1f} μM")
    log(f"YELLOW (0.25 ≤ f_PC < 0.50): {len(yellow)}")
    for r in yellow[:5]:
        c = "COV" if r["is_covalent"] else "NC "
        log(f"  → {c} {r['name'][:50]:50s} f_PC={r['f_PC_estimate']:.3f} Kd={r['Kd_uM_from_mean']:.1f} μM")
    log(f"RED: {len(red)}")

    OUT_JSON.write_text(json.dumps({
        "phase": "3c-v3b + 6b",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "library_size": len(library),
        "rep_snapshots": rep_snaps,
        "stage1_results_top100": [
            {"name": n, "smiles": s, "dG": dG, "is_covalent": cov}
            for (n, s, dG, cov) in s1[:100]
        ],
        "stage2_results": s2,
        "verdict": {
            "n_green": len(green),
            "n_yellow": len(yellow),
            "n_red": len(red),
            "n_cov_green": len(cov_green),
            "n_nc_green": len(nc_green),
            "top_green": green[:10],
            "top_yellow": yellow[:10],
        },
        "parameters": {
            "stage1_exhaustiveness": STAGE1_EXHAUSTIVENESS,
            "stage2_exhaustiveness": STAGE2_EXHAUSTIVENESS,
            "shortlist_n": args.shortlist_n,
            "n_clusters": N_CLUSTERS,
            "pocket_residues": K1141_POCKET_RESIDUES,
            "box_centre": BOX_CENTRE,
            "box_size": BOX_SIZE,
            "L_therapeutic_uM": L_THERAPEUTIC_uM_HIGH,
            "rescue_eta": RESCUE_ETA_CONSERVATIVE,
        },
    }, indent=2, default=str))
    log(f"Written {OUT_JSON}")


if __name__ == "__main__":
    main()
