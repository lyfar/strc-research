#!/usr/bin/env python3
"""
Phase 3c v2b — Two-stage ensemble docking of Phase 3c v2a library against
K1141 pocket with MD-sampled receptor conformers.

Stage 1 (breadth): dock all 667 ligands against ONE representative receptor
(snap closest to global k-means centroid) at Vina exhaustiveness 8.
Rank by ΔG, shortlist top 30.

Stage 2 (depth): dock shortlist × 5 k-means-selected receptor conformers at
exhaustiveness 16. Aggregate mean ± std ΔG → Kd (μM) → bound fraction at
therapeutic [L]=10 μM → f_PC estimate. Verdict:
  • f_PC ≥ 0.50 → GREEN (NORMAL rescue viable)
  • 0.25 ≤ f_PC < 0.50 → YELLOW (MILD-MODERATE rescue, stack lever)
  • f_PC < 0.25 → RED (insufficient, library expansion further needed)

Phase 5c established K1141 pocket is stable in MD (Cα RMSF 0.62 Å, void vol
720-850 Å³ across 20 frames), so k-means ensemble is tight but still reduces
noise from pose selection.

Deps: numpy, scipy, rdkit (implicit via 3a2a manifest), obabel, AutoDock Vina.
Runtime ~90 min on local Mac (55 min Stage 1 + 30 min Stage 2).
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Reuse hybrid36-safe stripper from Phase 5b
sys.path.insert(0, "/Users/egorlyfar/Brain/research/strc/models")
from pharmacochaperone_phase5b_ensemble_redock import (  # type: ignore
    strip_water_ions_to_pdb,
    pdb_to_pdbqt,
    vina_dock,
    BOX_CENTRE,
    BOX_SIZE,
    RT_KCAL,
    L_THERAPEUTIC_uM_HIGH,
    RESCUE_ETA_CONSERVATIVE,
)

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
SNAPSHOT_DIR = WORK / "artifacts" / "phase5a_snapshots"
LIB_JSON = WORK / "pharmacochaperone_phase3c_v2a_library_build.json"
DOCK_DIR = WORK / "docking_runs" / "3c_v2"
STAGE1_DIR = DOCK_DIR / "stage1"
STAGE2_DIR = DOCK_DIR / "stage2"
RECEPTOR_DIR = DOCK_DIR / "receptors"
OUT_JSON = WORK / "pharmacochaperone_phase3c_v2b_ensemble_dock.json"

# Phase 5c pocket-lining residues (9 residues; stable RMSF 0.62 Å)
K1141_POCKET_RESIDUES = [66, 67, 68, 538, 541, 542, 570, 571, 572]

SHORTLIST_N = 30
N_CLUSTERS = 5
STAGE1_EXHAUSTIVENESS = 8
STAGE2_EXHAUSTIVENESS = 16


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


def pocket_ca_coords(pdb_path: Path) -> np.ndarray:
    """Extract Cα coordinates of the 9 K1141 pocket residues.
    Uses OpenMM (hybrid36-safe). Returns (9, 3) array indexed by residue
    position in K1141_POCKET_RESIDUES."""
    from openmm.app import PDBFile
    pdb = PDBFile(str(pdb_path))
    # Map sequential residue index → Cα coords (after strip_water_ions
    # gives sequential residue IDs 1..N)
    # NB: pdb_path here is ALREADY the stripped protein PDB (1..701 residues)
    ca: dict[int, np.ndarray] = {}
    res_counter = 0
    for chain in pdb.topology.chains():
        if chain.id != "A":
            continue
        for res in chain.residues():
            res_counter += 1
            for atom in res.atoms():
                if atom.name == "CA":
                    p = pdb.positions[atom.index]
                    ca[res_counter] = np.array([
                        p.x.value_in_unit(p.unit) * 10,
                        p.y.value_in_unit(p.unit) * 10,
                        p.z.value_in_unit(p.unit) * 10,
                    ]) if hasattr(p, "x") else np.array([
                        p[0].value_in_unit(p[0].unit) * 10,
                        p[1].value_in_unit(p[1].unit) * 10,
                        p[2].value_in_unit(p[2].unit) * 10,
                    ])
                    break
    coords = np.stack([ca[r] for r in K1141_POCKET_RESIDUES if r in ca], axis=0)
    return coords


def kmeans_1d_pocket(frames_pocket: np.ndarray, k: int, seed: int = 42) -> list[int]:
    """Lightweight k-means on flattened pocket coords (F frames × 9*3).
    Returns: one index per cluster (the frame closest to each centroid)."""
    rng = np.random.default_rng(seed)
    X = frames_pocket.reshape(len(frames_pocket), -1)
    # Initialize by k-means++
    centroids = [X[rng.integers(0, len(X))]]
    for _ in range(k - 1):
        dists = np.array([min(np.sum((x - c) ** 2) for c in centroids) for x in X])
        probs = dists / dists.sum()
        idx = rng.choice(len(X), p=probs)
        centroids.append(X[idx])
    C = np.stack(centroids, axis=0)
    for _ in range(50):
        d = np.linalg.norm(X[:, None, :] - C[None, :, :], axis=-1)  # (F, k)
        assign = d.argmin(axis=1)
        new_C = np.stack([
            X[assign == i].mean(axis=0) if (assign == i).any() else C[i]
            for i in range(k)
        ], axis=0)
        if np.allclose(new_C, C, atol=1e-4):
            break
        C = new_C
    # Closest frame to each centroid
    reps = []
    for i in range(k):
        d_i = np.linalg.norm(X - C[i], axis=1)
        reps.append(int(d_i.argmin()))
    return reps


def select_receptors() -> tuple[list[int], list[Path]]:
    """Load 20 snapshots, extract K1141 pocket Cα, k-means → 5 reps.
    Also produce PDBQT for each rep + a 'mid' rep (snap closest to global
    centroid of all 20)."""
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
    frames_pocket = np.stack(frames_pocket, axis=0)  # (20, 9, 3)
    log(f"Pocket Cα tensor: {frames_pocket.shape}")

    # Global representative: frame closest to mean
    mean = frames_pocket.mean(axis=0)
    dists_mean = np.linalg.norm(
        frames_pocket.reshape(len(frames_pocket), -1) - mean.flatten(),
        axis=1,
    )
    rep_global = int(dists_mean.argmin())
    log(f"Global rep (closest to mean pocket): snap_{rep_global:03d}")

    # K-means for 5 representatives
    reps = kmeans_1d_pocket(frames_pocket, k=N_CLUSTERS, seed=42)
    log(f"K-means reps: {[f'snap_{r:03d}' for r in reps]}")

    # Prepare receptor PDBQTs
    rep_pdbqts = []
    for r in reps:
        pdbqt = RECEPTOR_DIR / f"snap_{r:03d}_receptor.pdbqt"
        if not pdbqt.exists():
            pdb_to_pdbqt(stripped_pdbs[r], pdbqt)
        rep_pdbqts.append(pdbqt)

    global_pdbqt = RECEPTOR_DIR / f"snap_{rep_global:03d}_receptor.pdbqt"
    if not global_pdbqt.exists():
        pdb_to_pdbqt(stripped_pdbs[rep_global], global_pdbqt)

    # Put global rep first in stage2 list if not already
    if rep_global not in reps:
        reps = [rep_global] + reps[:-1]
        rep_pdbqts = [global_pdbqt] + rep_pdbqts[:-1]
    return reps, rep_pdbqts


def vina_dock_one(receptor: Path, ligand: Path, out_dir: Path,
                  exhaustiveness: int, num_modes: int,
                  timeout: int = 300) -> float | None:
    cx, cy, cz = BOX_CENTRE
    sx, sy, sz = BOX_SIZE
    out_pdbqt = out_dir / f"{ligand.stem}__{receptor.stem}.pdbqt"
    log_txt = out_dir / f"{ligand.stem}__{receptor.stem}.log"
    cmd = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
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
        log(f"    TIMEOUT {ligand.stem}")
        return None
    log_txt.write_text(r.stdout + "\n----STDERR----\n" + r.stderr)
    for line in r.stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def stage1_breadth_screen(
    library: list[dict], rep_pdbqt: Path
) -> list[tuple[str, str, float | None]]:
    """Dock all ligands against one receptor. Return list of
    (name, smiles, dG) sorted by dG ascending (None at end)."""
    STAGE1_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    t0 = time.time()
    for i, entry in enumerate(library):
        lig_path = Path(entry["pdbqt_path"])
        if not lig_path.exists():
            continue
        dG = vina_dock_one(rep_pdbqt, lig_path, STAGE1_DIR,
                           STAGE1_EXHAUSTIVENESS, num_modes=3, timeout=120)
        results.append((entry["name"], entry["smiles"], dG))
        if (i + 1) % 25 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (len(library) - i - 1) / rate
            log(f"  stage1 {i+1}/{len(library)} done "
                f"({rate:.1f} lig/s, ETA {eta/60:.1f} min)")
    results.sort(key=lambda r: (r[2] if r[2] is not None else 0.0))
    return results


def stage2_ensemble_shortlist(
    shortlist: list[tuple[str, str, float | None]],
    library_map: dict[str, dict],
    receptors: list[Path],
) -> list[dict]:
    """Dock each shortlist ligand against all 5 receptor conformers at
    exhaustiveness 16. Aggregate mean ± std, convert to Kd + f_PC."""
    STAGE2_DIR.mkdir(parents=True, exist_ok=True)
    final = []
    t0 = time.time()
    for i, (name, smi, _dG_stage1) in enumerate(shortlist):
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
        # Kd (μM) from ensemble mean
        Kd_uM = math.exp(mean_dG / RT_KCAL) * 1e6
        # Bound fraction at [L] = 10 μM
        theta = L_THERAPEUTIC_uM_HIGH / (L_THERAPEUTIC_uM_HIGH + Kd_uM)
        f_PC = theta * RESCUE_ETA_CONSERVATIVE
        final.append({
            "name": name,
            "smiles": smi,
            "mean_dG": round(mean_dG, 3),
            "std_dG": round(std_dG, 3),
            "per_receptor_dG": [round(x, 3) for x in per_rec],
            "Kd_uM_from_mean": round(Kd_uM, 2),
            "bound_fraction_at_10uM": round(theta, 4),
            "f_PC_estimate": round(f_PC, 4),
            "MW": library_map[name].get("MW"),
            "logP": library_map[name].get("logP"),
            "HBD": library_map[name].get("HBD"),
            "HBA": library_map[name].get("HBA"),
            "TPSA": library_map[name].get("TPSA"),
        })
        elapsed = time.time() - t0
        rate = (i + 1) / elapsed
        log(f"  stage2 {i+1}/{len(shortlist)} {name} "
            f"ΔG={mean_dG:.2f}±{std_dG:.2f} "
            f"Kd={Kd_uM:.1f} μM f_PC={f_PC:.3f} "
            f"({rate*60:.1f} lig/min)")
    final.sort(key=lambda x: x["mean_dG"])
    return final


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["1", "2", "both"], default="both")
    parser.add_argument("--shortlist-n", type=int, default=SHORTLIST_N)
    args = parser.parse_args()

    lib = json.loads(LIB_JSON.read_text())
    library = lib["ligands"]
    library_map = {e["name"]: e for e in library}
    log(f"=== Phase 3c v2b ensemble dock ===")
    log(f"Library size: {len(library)} ligands")

    # Select receptors
    rep_snaps, rep_pdbqts = select_receptors()
    log(f"Receptors: {[p.stem for p in rep_pdbqts]}")

    stage1_results = []
    if args.stage in ("1", "both"):
        log("=== Stage 1: breadth screen ===")
        stage1_results = stage1_breadth_screen(library, rep_pdbqts[0])
        log(f"Stage 1 complete; top 10 by ΔG:")
        for r in stage1_results[:10]:
            log(f"  {r[2]} kcal/mol  {r[0]}")

    stage2_results = []
    if args.stage in ("2", "both"):
        log(f"=== Stage 2: ensemble dock on top {args.shortlist_n} ===")
        shortlist = stage1_results[:args.shortlist_n] if stage1_results else []
        if not shortlist and args.stage == "2":
            # Load prior stage1 results if present
            pass
        stage2_results = stage2_ensemble_shortlist(
            shortlist, library_map, rep_pdbqts,
        )
        log("Stage 2 final ranking:")
        for i, r in enumerate(stage2_results[:15]):
            log(f"  #{i+1} {r['name']:40s} "
                f"ΔG={r['mean_dG']:.2f}±{r['std_dG']:.2f} "
                f"Kd={r['Kd_uM_from_mean']:.1f} μM "
                f"f_PC={r['f_PC_estimate']:.3f}")

    # Verdict
    green = [r for r in stage2_results if r["f_PC_estimate"] >= 0.50]
    yellow = [r for r in stage2_results
              if 0.25 <= r["f_PC_estimate"] < 0.50]
    red = [r for r in stage2_results if r["f_PC_estimate"] < 0.25]
    log("=== VERDICT ===")
    log(f"GREEN (f_PC ≥ 0.50, NORMAL-rescue viable): {len(green)}")
    for r in green[:5]:
        log(f"  → {r['name']:40s} f_PC={r['f_PC_estimate']:.3f} "
            f"Kd={r['Kd_uM_from_mean']:.1f} μM")
    log(f"YELLOW (0.25 ≤ f_PC < 0.50, MILD/MODERATE rescue): {len(yellow)}")
    for r in yellow[:5]:
        log(f"  → {r['name']:40s} f_PC={r['f_PC_estimate']:.3f} "
            f"Kd={r['Kd_uM_from_mean']:.1f} μM")
    log(f"RED: {len(red)}")

    payload = {
        "phase": "3c-v2b",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "library_size": len(library),
        "rep_snapshots": rep_snaps,
        "stage1_results_top50": [
            {"name": n, "smiles": s, "dG": dG}
            for (n, s, dG) in stage1_results[:50]
        ],
        "stage2_results": stage2_results,
        "verdict": {
            "n_green": len(green),
            "n_yellow": len(yellow),
            "n_red": len(red),
            "top_green": green[:5],
            "top_yellow": yellow[:5],
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
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str))
    log(f"Written {OUT_JSON}")


if __name__ == "__main__":
    main()
