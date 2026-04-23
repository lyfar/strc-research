#!/usr/bin/env python3
"""
Phase 5e — Ensemble re-docking of Phase 3c v2 hits + Phase 4b shortlist on
Phase 5d E1659A MUTANT MD snapshots.

WHY
════
Phase 5b ran ensemble re-dock on WT Ultra-Mini (Phase 5a snapshots, no
mutation present). Phase 5d runs MD on the AF3-predicted full-length
E1659A mutant (job3-mutant.cif). This phase re-scores the same ligand
set against the mutant ensemble and compares f_PC(mut) vs f_PC(WT).

QUESTIONS ANSWERED
══════════════════
 1. Does E1659A shift K1141 pocket geometry enough to change predicted
    Kd for our Phase 3c v2 best hits?
 2. Is the pocket chemically more / less accessible in the mutant?
 3. Can the Phase 3c v2 fenamic hits (niflumic, flufenamic, sulfasalazine)
    and tafamidis-analog probe reach MILD-MODERATE f_PC thresholds on
    the actual disease target?

IDENTICAL PIPELINE TO Phase 5b for direct comparability — same box
coordinates, same Vina options, same ligand prep, same aggregation.

NOTE on chain A numbering
═════════════════════════
Phase 5a snapshots used Ultra-Mini sequential residue numbering 1-701.
Phase 5d snapshots use FULL-LENGTH AF3 numbering 1-1775. The K1141 pocket
centre (−22.027, −18.547, 2.215) is defined in Phase 5a's coordinate frame.
After MD on the mutant full-length structure, the absolute pocket centre
will drift from equilibration — we handle this by (a) aligning the mutant
prepped PDB's K1141 Cα to the Phase 5a pocket centre via rigid body fit
before docking, or (b) re-deriving pocket centre per snapshot from K1141
Cα + neighboring pocket residues (538-572 homologue in full-length).

The pragmatic choice for Phase 5e is (b): derive pocket box centre per
snapshot from K1141 Cα + adjacent Lys-pocket ring residues, keeping the
18 Å box size fixed.

Runtime: ~20 snapshots × 10 ligands × 8 CPU Vina @ 5 s = ~15 min wall.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
PHASE5D_JSON = WORK / "pharmacochaperone_phase5d_e1659a_md.json"
PHASE5B_JSON = WORK / "pharmacochaperone_phase5b_ensemble_redock.json"
DOCKING_DIR = WORK / "docking_runs" / "5e"
LIG_DIR_4B = WORK / "docking_runs" / "4b" / "ligands"
LIG_DIR_3CV2 = WORK / "docking_runs" / "3c_v2" / "ligands"
OUT_JSON = WORK / "pharmacochaperone_phase5e_mutant_ensemble_redock.json"

# Pocket box size identical to Phase 5b for direct comparability
BOX_SIZE = (18.0, 18.0, 18.0)

# K1141 pocket ring residues in FULL-LENGTH AF3 numbering.
# K1141 itself + 6 pocket-lining residues inferred from Phase 4a pocket
# reproducibility analysis (residue window 1135-1175 contains the lysine
# pocket structural motif in full-length). Box centre = centroid of these
# residues' Cα in each mutant snapshot.
POCKET_RING_RESNUMS_FULL = [1141, 1135, 1137, 1165, 1167, 1175]

# Ligands to re-dock
# Phase 4b top-5 (legacy shortlist; test whether old leads reach mutant)
LEGACY_LEADS = [
    "indole-3-acetic-acid",
    "naphthalene-2-carboxylic",
    "cyclopropane-phenyl-COOH",
    "salicylic-acid",
    "nicotinic-acid",
]
# Phase 3c v2 top fenamic hits + positive controls
V2_HITS = [
    "diflunisal",
    "niflumic-acid",
    "flufenamic-acid",
    "meclofenamic-acid",
    "sulfasalazine",
    "tafamidis-analog",
]

# Thermodynamics (identical to Phase 5b)
RT_KCAL = 0.5921  # kcal/mol at 298 K
L_THERAPEUTIC_uM_LOW = 1.0
L_THERAPEUTIC_uM_HIGH = 10.0
RESCUE_ETA_CONSERVATIVE = 0.5


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


STANDARD_AA_3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
}


def strip_water_ions_to_pdb(pdb_in: Path, pdb_out: Path) -> tuple[int, np.ndarray]:
    """Keep only standard protein residues on chain A. Returns
    (n_atoms, CA coord table indexed by full-length resnum)."""
    from openmm.app import PDBFile, Topology
    from openmm import unit

    pdb = PDBFile(str(pdb_in))
    top = pdb.topology
    positions = pdb.positions

    # Build CA lookup: resnum → (x, y, z) in Angstroms
    pos_ang_iter = positions.value_in_unit(unit.angstrom)
    ca_by_resnum: dict[int, np.ndarray] = {}

    new_top = Topology()
    new_chain = new_top.addChain("A")
    seen_atom_to_new = {}
    res_counter = 0
    kept_atom_indices = []
    for chain in top.chains():
        if chain.id != "A":
            continue
        for residue in chain.residues():
            if residue.name.strip().upper() not in STANDARD_AA_3:
                continue
            res_counter += 1
            # OpenMM residue.id preserves hybrid36 full-length number
            try:
                full_res_num = int(residue.id)
            except (ValueError, TypeError):
                full_res_num = res_counter
            new_res = new_top.addResidue(
                residue.name, new_chain, id=str(res_counter))
            for atom in residue.atoms():
                new_atom = new_top.addAtom(atom.name, atom.element, new_res)
                seen_atom_to_new[atom.index] = new_atom
                kept_atom_indices.append(atom.index)
                if atom.name.strip() == "CA":
                    p = pos_ang_iter[atom.index]
                    if hasattr(p, "x"):
                        ca_by_resnum[full_res_num] = np.array([p.x, p.y, p.z])
                    else:
                        ca_by_resnum[full_res_num] = np.array([p[0], p[1], p[2]])

    for b in top.bonds():
        a1, a2 = b[0], b[1]
        if a1.index in seen_atom_to_new and a2.index in seen_atom_to_new:
            new_top.addBond(seen_atom_to_new[a1.index],
                            seen_atom_to_new[a2.index])

    new_positions = [positions[i] for i in kept_atom_indices]
    new_positions_q = unit.Quantity(
        [p.value_in_unit(unit.nanometer) for p in new_positions],
        unit.nanometer)
    with open(pdb_out, "w") as fh:
        PDBFile.writeFile(new_top, new_positions_q, fh, keepIds=True)

    return len(kept_atom_indices), ca_by_resnum


def derive_box_centre(ca_by_resnum: dict[int, np.ndarray]) -> tuple[float, float, float]:
    """Compute K1141 pocket box centre = mean of pocket ring Cα."""
    found = [ca_by_resnum[r] for r in POCKET_RING_RESNUMS_FULL if r in ca_by_resnum]
    if len(found) < 3:
        raise RuntimeError(
            f"Only {len(found)} pocket ring residues found in snapshot "
            f"(need at least 3). Have: {[r for r in POCKET_RING_RESNUMS_FULL if r in ca_by_resnum]}")
    centre = np.mean(np.stack(found, axis=0), axis=0)
    return (float(centre[0]), float(centre[1]), float(centre[2]))


def pdb_to_pdbqt(pdb_in: Path, pdbqt_out: Path) -> None:
    cmd = ["obabel", str(pdb_in), "-O", str(pdbqt_out),
           "-xr", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if r.returncode != 0 or not pdbqt_out.exists():
        raise RuntimeError(f"obabel failed: {r.stderr}")


def vina_dock(receptor_pdbqt: Path, ligand_pdbqt: Path,
              box_centre: tuple[float, float, float],
              out_pdbqt: Path, log_txt: Path,
              cpu: int = 8) -> float | None:
    cx, cy, cz = box_centre
    sx, sy, sz = BOX_SIZE
    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", "16",
        "--num_modes", "5",
        "--cpu", str(cpu),
        "--out", str(out_pdbqt),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    log_txt.write_text(r.stdout + "\n----STDERR----\n" + r.stderr)
    for line in r.stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def resolve_ligand_pdbqt(name: str) -> Path | None:
    for d in (LIG_DIR_4B, LIG_DIR_3CV2):
        p = d / f"{name}.pdbqt"
        if p.exists():
            return p
    return None


def ensemble_redock(snapshots: list[Path], cpu: int = 8) -> dict:
    DOCKING_DIR.mkdir(parents=True, exist_ok=True)
    all_leads = LEGACY_LEADS + V2_HITS
    per_ligand: dict[str, list[float]] = {lig: [] for lig in all_leads}
    per_snapshot: dict[str, dict] = {}

    for i, snap in enumerate(snapshots):
        snap_id = f"snap_{i:03d}"
        log(f"=== {snap_id} / {len(snapshots)} ===")
        with tempfile.TemporaryDirectory() as tmpd:
            tmp = Path(tmpd)
            protein_pdb = tmp / f"{snap_id}_prot.pdb"
            n_atoms, ca_map = strip_water_ions_to_pdb(snap, protein_pdb)
            log(f"  chain A: {n_atoms} atoms, {len(ca_map)} Cα resnums")
            box_centre = derive_box_centre(ca_map)
            log(f"  pocket centre: ({box_centre[0]:.2f}, {box_centre[1]:.2f}, {box_centre[2]:.2f})")
            protein_pdbqt = DOCKING_DIR / f"{snap_id}_receptor.pdbqt"
            pdb_to_pdbqt(protein_pdb, protein_pdbqt)
            per_snapshot[snap_id] = {
                "box_centre": list(box_centre),
                "n_atoms": n_atoms,
                "results": {},
            }
            for lig in all_leads:
                lig_pdbqt = resolve_ligand_pdbqt(lig)
                if lig_pdbqt is None:
                    log(f"  MISSING ligand {lig}.pdbqt in 4b/ or 3c_v2/")
                    per_snapshot[snap_id]["results"][lig] = None
                    continue
                pose_out = DOCKING_DIR / f"{lig}__{snap_id}.pdbqt"
                log_out = DOCKING_DIR / f"{lig}__{snap_id}.log"
                aff = vina_dock(protein_pdbqt, lig_pdbqt,
                                box_centre, pose_out, log_out, cpu=cpu)
                per_snapshot[snap_id]["results"][lig] = aff
                if aff is not None:
                    per_ligand[lig].append(aff)
                log(f"  {lig:28s} ΔG={aff}")

    return {"per_ligand": per_ligand, "per_snapshot": per_snapshot}


def aggregate(per_ligand: dict[str, list[float]]) -> dict:
    out = {}
    for lig, vals in per_ligand.items():
        if not vals:
            out[lig] = {"n": 0, "valid": False}
            continue
        vals_sorted = sorted(vals)
        mean = sum(vals) / len(vals)
        if len(vals) > 1:
            variance = sum((x - mean) ** 2 for x in vals) / (len(vals) - 1)
            std = math.sqrt(variance)
        else:
            std = 0.0
        best = min(vals)
        median = vals_sorted[len(vals_sorted) // 2]
        kd_M = math.exp(mean / RT_KCAL)
        kd_uM = kd_M * 1e6
        theta_low = L_THERAPEUTIC_uM_LOW / (kd_uM + L_THERAPEUTIC_uM_LOW)
        theta_high = L_THERAPEUTIC_uM_HIGH / (kd_uM + L_THERAPEUTIC_uM_HIGH)
        f_PC_low = theta_low * RESCUE_ETA_CONSERVATIVE
        f_PC_high = theta_high * RESCUE_ETA_CONSERVATIVE
        out[lig] = {
            "n": len(vals),
            "valid": True,
            "mean_dG_kcal_mol": round(mean, 3),
            "std_dG_kcal_mol": round(std, 3),
            "median_dG_kcal_mol": round(median, 3),
            "best_dG_kcal_mol": round(best, 3),
            "Kd_uM": round(kd_uM, 3),
            "bound_fraction_at_1uM": round(theta_low, 3),
            "bound_fraction_at_10uM": round(theta_high, 3),
            "f_PC_at_1uM_conservative": round(f_PC_low, 3),
            "f_PC_at_10uM_conservative": round(f_PC_high, 3),
        }
    return out


def compare_with_wt(agg_mut: dict) -> dict:
    """Compare mutant f_PC to WT Phase 5b f_PC for same ligand."""
    if not PHASE5B_JSON.exists():
        return {"note": "Phase 5b JSON missing — skipping WT vs mut comparison"}
    p5b = json.loads(PHASE5B_JSON.read_text())
    wt = p5b.get("aggregated_per_ligand", {})
    comp = {}
    for lig, mut_stats in agg_mut.items():
        if not mut_stats.get("valid"):
            continue
        wt_stats = wt.get(lig)
        if not wt_stats or not wt_stats.get("valid"):
            comp[lig] = {"wt_available": False}
            continue
        dmean = mut_stats["mean_dG_kcal_mol"] - wt_stats["mean_dG_kcal_mol"]
        kd_ratio = mut_stats["Kd_uM"] / wt_stats["Kd_uM"] if wt_stats["Kd_uM"] > 0 else None
        comp[lig] = {
            "wt_available": True,
            "wt_mean_dG": wt_stats["mean_dG_kcal_mol"],
            "mut_mean_dG": mut_stats["mean_dG_kcal_mol"],
            "ddG_mut_minus_wt": round(dmean, 3),
            "kd_ratio_mut_over_wt": round(kd_ratio, 3) if kd_ratio else None,
            "wt_f_PC_10uM": wt_stats.get("f_PC_estimate_at_10uM_conservative"),
            "mut_f_PC_10uM": mut_stats["f_PC_at_10uM_conservative"],
            "f_PC_delta_mut_minus_wt": round(
                mut_stats["f_PC_at_10uM_conservative"]
                - (wt_stats.get("f_PC_estimate_at_10uM_conservative") or 0), 3),
            "verdict": (
                "mutant tighter" if dmean < -0.3
                else "mutant weaker" if dmean > 0.3
                else "no meaningful shift"),
        }
    return comp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase5d", default=str(PHASE5D_JSON))
    ap.add_argument("--cpu", type=int, default=8)
    args = ap.parse_args()

    log("Phase 5e mutant ensemble re-docking")
    p5d = json.loads(Path(args.phase5d).read_text())
    snapshots = [Path(p) for p in p5d["result"]["snapshot_pdb_paths"]]
    log(f"Loaded {len(snapshots)} Phase 5d mutant snapshots")

    dock_results = ensemble_redock(snapshots, cpu=args.cpu)
    aggregated = aggregate(dock_results["per_ligand"])
    wt_mut_compare = compare_with_wt(aggregated)

    ranked = sorted(
        [(lig, stats) for lig, stats in aggregated.items() if stats.get("valid")],
        key=lambda kv: kv[1]["mean_dG_kcal_mol"])

    summary = {
        "batch": "pharmacochaperone_phase5e_mutant_ensemble_redock",
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_snapshots": len(snapshots),
        "phase5d_source": str(args.phase5d),
        "leads_docked": LEGACY_LEADS + V2_HITS,
        "box_size": list(BOX_SIZE),
        "pocket_ring_resnums": POCKET_RING_RESNUMS_FULL,
        "per_snapshot_raw": dock_results["per_snapshot"],
        "aggregated_per_ligand": aggregated,
        "wt_vs_mutant_comparison": wt_mut_compare,
        "ranking": [{"name": lig,
                     "mean_dG": s["mean_dG_kcal_mol"],
                     "Kd_uM": s["Kd_uM"],
                     "f_PC_at_10uM": s["f_PC_at_10uM_conservative"]}
                    for lig, s in ranked],
        "notes": [
            "Box centre derived per snapshot from Cα centroid of K1141 + pocket "
            "ring residues (full-length numbering).",
            "Same Vina settings as Phase 5b for direct comparability.",
            "Key comparison: f_PC(mut) vs f_PC(WT) per ligand — if mutant "
            "Kd < 2× WT Kd, pocket geometry is robust to E1659A and WT-based "
            "docking was a valid proxy. If Kd > 3× WT → retarget on mutant.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    print()
    print("=" * 84)
    print("Phase 5e — Mutant Ensemble Re-Docking Summary")
    print("=" * 84)
    print(f"\n{'Ligand':<28}{'mean ΔG':<12}{'±std':<10}"
          f"{'Kd (μM)':<10}{'f_PC@10μM':<12}{'ΔΔG vs WT':<12}")
    print("─" * 84)
    for lig, stats in ranked:
        cmp = wt_mut_compare.get(lig, {})
        ddg = cmp.get("ddG_mut_minus_wt", "—") if cmp.get("wt_available") else "—"
        print(f"{lig:<28}{stats['mean_dG_kcal_mol']:<12.3f}"
              f"{stats['std_dG_kcal_mol']:<10.3f}"
              f"{stats['Kd_uM']:<10.3f}"
              f"{stats['f_PC_at_10uM_conservative']:<12.3f}"
              f"{ddg}")
    print(f"\nJSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
