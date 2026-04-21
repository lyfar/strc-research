#!/usr/bin/env python3
"""
Phase 4c — WT vs E1659A decoy docking.

Re-dock the same 9-compound roster (Phase 4b) against both WT STRC
(job4-wildtype.cif) and E1659A mutant (job3-mutant.cif). Reuses
ligand PDBQTs from Phase 4b.

The pharmacochaperone hypothesis predicts: rescue ligands bind the
E1659A pocket BETTER than the WT pocket, because in WT the K1141
sidechain is already occupied by its native E1659 salt-bridge
partner, whereas in E1659A the Ala1659 leaves K1141 free to accept
a carboxylate ligand.

Expected signal:
  - Leads: ΔΔG(WT − E1659A) >= +1 kcal/mol (leads prefer mutant)
  - Diflunisal (non-STRC-specific positive control): ΔΔG ≈ 0
  - Negatives (glucose/urea/acetamide): ΔΔG ≈ 0 (noise)

Gate: >=3 of 5 leads show ΔΔG(WT − E1659A) >= 1 kcal/mol.

If leads bind WT equally or better than mutant, the pocket is not
E1659A-rescue-specific — it's just "a pocket on STRC" — and the
pharmacochaperone mechanism is falsified.

Threshold rationale: Vina baseline error ~1.5 kcal/mol per run but
~0.5-0.8 kcal/mol for paired relative scores on the same system.
1 kcal/mol gate is ~2x noise; tight but defensible.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pharmacochaperone_phase4_common import (  # noqa: E402
    gate_ready, print_env_report, TARGETS, FIXED_ROSTER,
    MODELS_DIR, WORK_DIR,
)
from pharmacochaperone_phase4b_vina_gnina_screen import (  # noqa: E402
    cif_to_chain_pdb, prep_receptor_pdbqt, box_centre_for_target,
    run_vina, prep_ligand_pdbqt,
)

PHASE = "4c"
OUT_JSON = WORK_DIR / f"pharmacochaperone_phase{PHASE}_wt_decoy.json"
RUN_DIR = WORK_DIR / "docking_runs" / PHASE
PHASE4B_LIG_DIR = WORK_DIR / "docking_runs" / "4b" / "ligands"
DELTA_GATE_KCAL = 1.0
MIN_LEADS_PASSING = 3


def dock_roster(target_key: str, work_dir: Path) -> list[dict]:
    centre, size = box_centre_for_target(target_key)
    t = TARGETS[target_key]
    work_dir.mkdir(parents=True, exist_ok=True)
    rec_pdb = work_dir / f"{target_key}_chainA.pdb"
    rec_pdbqt = work_dir / f"{target_key}_chainA.pdbqt"
    if not rec_pdbqt.exists():
        print(f"  preparing receptor PDBQT for {target_key}...", flush=True)
        cif_to_chain_pdb(MODELS_DIR / t["cif"], t["chain"], rec_pdb)
        prep_receptor_pdbqt(rec_pdb, rec_pdbqt)
    print(f"  target {target_key:<16s}  box_c={[round(float(x),2) for x in centre]}")

    results = []
    for cmpd in FIXED_ROSTER:
        lig_pdbqt = PHASE4B_LIG_DIR / f"{cmpd['name']}.pdbqt"
        if not lig_pdbqt.exists():
            print(f"  prepping ligand (fallback — 4b prep missing) {cmpd['name']}")
            lig_pdbqt = work_dir / "ligands" / f"{cmpd['name']}.pdbqt"
            lig_pdbqt.parent.mkdir(parents=True, exist_ok=True)
            try:
                prep_ligand_pdbqt(cmpd["name"], cmpd["smiles"], lig_pdbqt)
            except Exception as e:
                results.append({**cmpd, "best_affinity_kcal_mol": None, "error": str(e)})
                continue
        pose_pdbqt = work_dir / "poses" / f"{cmpd['name']}__{target_key}.pdbqt"
        pose_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        log_path = work_dir / "logs" / f"{cmpd['name']}__{target_key}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            print(f"  docking  {cmpd['name']:30s} on {target_key}", flush=True)
            score = run_vina(rec_pdbqt, lig_pdbqt, centre, size, pose_pdbqt, log_path)
            print(f"    ΔG = {score:+.2f} kcal/mol")
            results.append({
                **cmpd,
                "best_affinity_kcal_mol": round(float(score), 3),
                "pose_pdbqt": str(pose_pdbqt),
            })
        except Exception as e:
            print(f"  VINA FAILED for {cmpd['name']} on {target_key}: {e}")
            results.append({**cmpd, "best_affinity_kcal_mol": None, "error": str(e)})
    return results


def main():
    print_env_report(PHASE)
    ok, missing = gate_ready(PHASE)
    if not ok:
        print(f"BLOCKED: {missing}")
        return 2
    if not PHASE4B_LIG_DIR.exists():
        print(f"Phase 4b ligand prep directory missing: {PHASE4B_LIG_DIR}")
        return 2

    print("Phase 4c — WT vs E1659A dock, same 9-compound roster")
    wt_results  = dock_roster("wt_full",       RUN_DIR)
    mut_results = dock_roster("e1659a_mutant", RUN_DIR)

    # Merge and compute ΔΔG per compound.
    by_name_wt  = {r["name"]: r for r in wt_results}
    by_name_mut = {r["name"]: r for r in mut_results}
    paired = []
    for cmpd in FIXED_ROSTER:
        wt  = by_name_wt.get(cmpd["name"])
        mut = by_name_mut.get(cmpd["name"])
        wt_s  = wt.get("best_affinity_kcal_mol")  if wt  else None
        mut_s = mut.get("best_affinity_kcal_mol") if mut else None
        ddg = (wt_s - mut_s) if (wt_s is not None and mut_s is not None) else None
        paired.append({
            "name": cmpd["name"],
            "role": cmpd["role"],
            "smiles": cmpd["smiles"],
            "wt_dG_kcal_mol":     wt_s,
            "e1659a_dG_kcal_mol": mut_s,
            "delta_wt_minus_mut": (round(float(ddg), 3) if ddg is not None else None),
        })

    leads = [p for p in paired if p["role"] == "lead" and p["delta_wt_minus_mut"] is not None]
    pos   = [p for p in paired if p["role"] == "positive" and p["delta_wt_minus_mut"] is not None]
    negs  = [p for p in paired if p["role"] == "negative" and p["delta_wt_minus_mut"] is not None]

    leads_passing = sum(1 for l in leads if l["delta_wt_minus_mut"] >= DELTA_GATE_KCAL)
    gate_verdict = "PASS" if leads_passing >= MIN_LEADS_PASSING else "FAIL"

    summary = {
        "phase": PHASE,
        "targets_compared": ["wt_full", "e1659a_mutant"],
        "delta_gate_kcal_mol": DELTA_GATE_KCAL,
        "min_leads_passing": MIN_LEADS_PASSING,
        "leads_passing": leads_passing,
        "verdict": gate_verdict,
        "mean_delta_leads":    (round(float(np.mean([l["delta_wt_minus_mut"] for l in leads])), 3) if leads else None),
        "mean_delta_positive": (round(float(np.mean([l["delta_wt_minus_mut"] for l in pos])),   3) if pos   else None),
        "mean_delta_negative": (round(float(np.mean([l["delta_wt_minus_mut"] for l in negs])),  3) if negs  else None),
    }
    payload = {"summary": summary, "paired": paired,
               "wt_results": wt_results, "e1659a_results": mut_results}
    OUT_JSON.write_text(json.dumps(payload, indent=2))

    print()
    print(f"  compound                         ΔG_WT    ΔG_E1659A   Δ(WT−mut)")
    for p in paired:
        wt_s = p["wt_dG_kcal_mol"]; mut_s = p["e1659a_dG_kcal_mol"]
        d = p["delta_wt_minus_mut"]
        print(f"  {p['role']:8s} {p['name']:28s} "
              f"{wt_s if wt_s is not None else 'NA':>7} "
              f"{mut_s if mut_s is not None else 'NA':>10} "
              f"{(f'{d:+.3f}' if d is not None else 'NA'):>12}")
    print()
    print(f"  leads passing ΔΔG ≥ {DELTA_GATE_KCAL} kcal/mol: {leads_passing}/{len(leads)}")
    print(f"  mean ΔΔG — leads {summary['mean_delta_leads']}, "
          f"positive {summary['mean_delta_positive']}, "
          f"negatives {summary['mean_delta_negative']}")
    print()
    print(f"verdict: {gate_verdict}")
    print(f"wrote {OUT_JSON}")
    return 0 if gate_verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
