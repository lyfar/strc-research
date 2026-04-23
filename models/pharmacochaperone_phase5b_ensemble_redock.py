#!/usr/bin/env python3
"""
Phase 5b — Ensemble re-docking of Phase 4b shortlist on Phase 5a MD snapshots.

Rationale: single-structure Vina docking (Phase 4b) has ±2 kcal/mol absolute
affinity noise from lucky/unlucky poses and single-conformation bias. Ensemble
docking against MD-sampled receptor conformations reduces this noise and
produces a mean ± stdev affinity that is more reliable for Kd prediction.

Pipeline:
  1. Load Phase 5a snapshot PDB paths from JSON.
  2. For each snapshot: extract chain A (strip water + ions) → PDB → PDBQT
     via obabel (rigid receptor, Gasteiger charges, pH 7.4).
  3. For each lead + diflunisal positive control:
       Vina dock against EACH snapshot using SAME BOX as Phase 4b
       (centre [-22.027, -18.547, 2.215], 18×18×18 Å).
  4. Collect per-(ligand, snapshot) best ΔG_bind.
  5. Aggregate: mean, std, median per ligand.
  6. Convert to Kd (μM) via ΔG = -RT ln(Kd); compute bound fraction at
     therapeutic intracochlear [L] = 1 μM and 10 μM.
  7. Estimate f_PC = bound_fraction × 0.5  (assumption: 50% of bound
     molecules deliver functional rescue — Phase 6 refinement needed).
  8. Feed f_PC estimate into [[Misha Compound-Het Therapy Stack Model]]
     reachability check.

Assumptions & caveats:
  • Vina score is empirical, trained on PDBbind; systematic bias +2 to -1
    kcal/mol per ligand class (carboxylates in our case); use as RANKING
    + Kd-order-of-magnitude, not exact.
  • Ensemble docking reduces noise but doesn't correct systematic bias.
  • Phase 4b Kd is HIGHLY uncertain; Phase 5b narrows the random noise but
    for absolute Kd we still want MM-GBSA (Phase 5c, conda-env-gated).

Runtime: ~ 30 snapshots × 5 ligands × 8 CPU Vina @ 5 s per dock =
  ~ 20 min wall time. Scales linearly in snapshots and ligands.
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

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
PHASE5A_JSON = WORK / "pharmacochaperone_phase5a_apo_md_smoke.json"
PHASE4B_JSON = WORK / "pharmacochaperone_phase4b_vina_gnina_screen.json"
DOCKING_DIR = WORK / "docking_runs" / "5b"
LIGAND_DIR = WORK / "docking_runs" / "4b" / "ligands"
OUT_JSON = WORK / "pharmacochaperone_phase5b_ensemble_redock.json"

# Phase 4b box
BOX_CENTRE = (-22.027, -18.547, 2.215)
BOX_SIZE = (18.0, 18.0, 18.0)

# Leads to re-dock (top-5 from Phase 4b + diflunisal positive)
LEADS_WITH_POSITIVE = [
    "indole-3-acetic-acid",
    "naphthalene-2-carboxylic",
    "cyclopropane-phenyl-COOH",
    "salicylic-acid",
    "nicotinic-acid",
    "diflunisal",
]

# Thermodynamics
RT_KCAL = 0.5921  # kcal/mol at 298 K
# Therapeutic intracochlear drug concentration from Phase 4e PKPD window
L_THERAPEUTIC_uM_LOW = 1.0
L_THERAPEUTIC_uM_HIGH = 10.0
# f_PC = bound_fraction × rescue efficiency per bound molecule (eta)
RESCUE_ETA_CONSERVATIVE = 0.5


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


def strip_water_ions_to_pdb(pdb_in: Path, pdb_out: Path) -> int:
    """Keep only standard protein AAs (ATOM + chain A). Return num atoms kept."""
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.Polypeptide import is_aa

    class OnlyProtChainA(Select):
        def accept_chain(self, c):
            return c.id == "A"

        def accept_residue(self, r):
            return is_aa(r, standard=True)

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", str(pdb_in))
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(pdb_out), select=OnlyProtChainA())
    # count atoms in new file
    n = sum(1 for line in pdb_out.read_text().splitlines()
            if line.startswith(("ATOM", "HETATM")))
    return n


def pdb_to_pdbqt(pdb_in: Path, pdbqt_out: Path) -> None:
    cmd = ["obabel", str(pdb_in), "-O", str(pdbqt_out),
           "-xr", "--partialcharge", "gasteiger"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if r.returncode != 0 or not pdbqt_out.exists():
        raise RuntimeError(f"obabel failed: {r.stderr}")


def vina_dock(receptor_pdbqt: Path, ligand_pdbqt: Path,
              out_pdbqt: Path, log_txt: Path) -> float | None:
    """Run Vina, return best affinity kcal/mol."""
    cx, cy, cz = BOX_CENTRE
    sx, sy, sz = BOX_SIZE
    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", "16",    # reduced from 32 for ensemble speed
        "--num_modes", "5",
        "--cpu", "8",
        "--out", str(out_pdbqt),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    log_txt.write_text(r.stdout + "\n----STDERR----\n" + r.stderr)
    # Parse best affinity from stdout table
    for line in r.stdout.splitlines():
        m = re.match(r"^\s*1\s+(-?\d+\.\d+)", line)
        if m:
            return float(m.group(1))
    return None


def ensemble_redock(snapshots: list[Path]) -> dict:
    DOCKING_DIR.mkdir(parents=True, exist_ok=True)
    per_ligand: dict[str, list[float]] = {lig: [] for lig in LEADS_WITH_POSITIVE}
    per_snapshot: dict[str, dict[str, float | None]] = {}

    for i, snap in enumerate(snapshots):
        snap_id = f"snap_{i:03d}"
        log(f"=== {snap_id} / {len(snapshots)} ===")
        # Strip water/ions
        with tempfile.TemporaryDirectory() as tmpd:
            tmp = Path(tmpd)
            protein_pdb = tmp / f"{snap_id}_prot.pdb"
            n_atoms = strip_water_ions_to_pdb(snap, protein_pdb)
            log(f"  chain A extracted: {n_atoms} atoms")
            protein_pdbqt = DOCKING_DIR / f"{snap_id}_receptor.pdbqt"
            pdb_to_pdbqt(protein_pdb, protein_pdbqt)
            log(f"  receptor PDBQT ready")

            per_snapshot[snap_id] = {}
            for lig in LEADS_WITH_POSITIVE:
                lig_pdbqt = LIGAND_DIR / f"{lig}.pdbqt"
                if not lig_pdbqt.exists():
                    log(f"  MISSING ligand {lig_pdbqt}")
                    per_snapshot[snap_id][lig] = None
                    continue
                pose_out = DOCKING_DIR / f"{lig}__{snap_id}.pdbqt"
                log_out = DOCKING_DIR / f"{lig}__{snap_id}.log"
                aff = vina_dock(protein_pdbqt, lig_pdbqt, pose_out, log_out)
                per_snapshot[snap_id][lig] = aff
                if aff is not None:
                    per_ligand[lig].append(aff)
                log(f"  {lig:30s} ΔG={aff}")

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

        # Kd estimation from mean ΔG: Kd = exp(ΔG / RT) in M units
        # Vina ΔG is already in kcal/mol; Kd = exp(ΔG / RT_kcal)
        kd_M = math.exp(mean / RT_KCAL)
        kd_uM = kd_M * 1e6

        # Bound fraction at therapeutic [L]
        theta_low = L_THERAPEUTIC_uM_LOW / (kd_uM + L_THERAPEUTIC_uM_LOW)
        theta_high = L_THERAPEUTIC_uM_HIGH / (kd_uM + L_THERAPEUTIC_uM_HIGH)

        # Conservative f_PC estimate
        f_PC_at_low = theta_low * RESCUE_ETA_CONSERVATIVE
        f_PC_at_high = theta_high * RESCUE_ETA_CONSERVATIVE

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
            "f_PC_estimate_at_1uM_conservative": round(f_PC_at_low, 3),
            "f_PC_estimate_at_10uM_conservative": round(f_PC_at_high, 3),
        }
    return out


def misha_implications(aggregated: dict) -> dict:
    """
    Apply aggregated f_PC per lead to Misha compound-het stack model.

    NORMAL rescue threshold (per [[Misha Compound-Het Therapy Stack Model]]):
      mild E1659A (f_mat=0.40): f_PC ≥ 0.50 gives NORMAL monotherapy
      moderate (f_mat=0.25):    f_PC ≥ 0.60
      severe (f_mat=0.10):      f_PC ≥ 0.67
    """
    THRESHOLDS = {"mild_E1659A": 0.50,
                  "moderate_E1659A": 0.60,
                  "severe_E1659A": 0.67}

    impl = {}
    for lig, stats in aggregated.items():
        if not stats.get("valid"):
            impl[lig] = {"applicable": False}
            continue
        f_PC_10 = stats["f_PC_estimate_at_10uM_conservative"]
        meets = {scen: bool(f_PC_10 >= thr) for scen, thr in THRESHOLDS.items()}
        impl[lig] = {
            "f_PC_at_therapeutic_10uM": f_PC_10,
            "meets_NORMAL_threshold_per_scenario": meets,
            "supports_Misha_cure_monotherapy_if": [
                s for s, m in meets.items() if m],
        }
    return impl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase5a", default=str(PHASE5A_JSON))
    args = ap.parse_args()

    log("Phase 5b ensemble re-docking")
    phase5a = json.loads(Path(args.phase5a).read_text())
    snapshots = [Path(p) for p in phase5a["result"]["snapshot_pdb_paths"]]
    log(f"Loaded {len(snapshots)} Phase 5a snapshots")

    dock_results = ensemble_redock(snapshots)
    aggregated = aggregate(dock_results["per_ligand"])
    misha = misha_implications(aggregated)

    # Phase 4b comparison for context
    phase4b = json.loads(PHASE4B_JSON.read_text())
    phase4b_single = {r["name"]: r["best_affinity_kcal_mol"]
                      for r in phase4b["results"]}

    # Ranking decision support
    ranked = sorted(
        [(lig, stats) for lig, stats in aggregated.items() if stats.get("valid")],
        key=lambda kv: kv[1]["mean_dG_kcal_mol"])

    summary = {
        "batch": "pharmacochaperone_phase5b_ensemble_redock",
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_snapshots": len(snapshots),
        "phase5a_source": str(args.phase5a),
        "leads_docked": LEADS_WITH_POSITIVE,
        "box_centre": list(BOX_CENTRE),
        "box_size": list(BOX_SIZE),
        "per_snapshot_raw": dock_results["per_snapshot"],
        "aggregated_per_ligand": aggregated,
        "phase4b_single_structure_comparison": phase4b_single,
        "misha_stack_implications": misha,
        "ranking": [{"name": lig,
                     "mean_dG": s["mean_dG_kcal_mol"],
                     "Kd_uM": s["Kd_uM"],
                     "f_PC_at_10uM": s["f_PC_estimate_at_10uM_conservative"]}
                    for lig, s in ranked],
        "notes": [
            "Vina ensemble ΔG reduces single-structure noise but has systematic "
            "bias that absolute Kd may be off by 10× each direction",
            "f_PC estimate uses RESCUE_ETA=0.5 assumption — 50% of bound drug "
            "molecules provide functional rescue. Phase 6 would refine via "
            "mechanism-specific QSAR.",
            "Kd from Vina is IC50-like, not true equilibrium Kd; interpret as "
            "order-of-magnitude estimate only.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    # Human summary
    print()
    print("=" * 80)
    print("Phase 5b — Ensemble Re-Docking Summary")
    print("=" * 80)
    print(f"\n{'Ligand':<28}{'mean ΔG':<12}{'±std':<10}"
          f"{'Kd (μM)':<10}{'f_PC @ 10 μM':<14}{'Phase4b':<10}")
    print("─" * 84)
    for lig, stats in ranked:
        p4b_aff = phase4b_single.get(lig, "?")
        print(f"{lig:<28}{stats['mean_dG_kcal_mol']:<12.3f}"
              f"{stats['std_dG_kcal_mol']:<10.3f}"
              f"{stats['Kd_uM']:<10.3f}"
              f"{stats['f_PC_estimate_at_10uM_conservative']:<14.3f}"
              f"{p4b_aff}")
    print()
    print("Misha cure threshold check (f_PC at [L]=10 μM vs NORMAL threshold):")
    for lig, impl in misha.items():
        if not impl.get("applicable", True):
            continue
        scenarios = impl["supports_Misha_cure_monotherapy_if"]
        tag = "; ".join(scenarios) if scenarios else "none"
        print(f"  {lig:<28} → {tag}")
    print(f"\nJSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
