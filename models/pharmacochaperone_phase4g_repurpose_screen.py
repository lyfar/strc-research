#!/usr/bin/env python3
"""
Phase 4g — Repurposable chaperone docking screen.

Motivation: two recent papers (2026-04-21 hayase-4pba-dscaml1-chaperone and
kaur-ip045-synuclein-chaperone — indexed by the STRC Research Monitor) argue
that generic / FDA-approved chemical chaperones can rescue missense-trapped
membrane proteins in vivo. If ANY of these score competitively against our
Phase 3C top-5 leads on the E1659A pocket, the clinical path shortens
dramatically — 4PBA is already FDA-approved for urea cycle disorder.

Scope of this gate:
- Compounds: 4PBA (FDA), IP-045 (clinical chaperone, indole class),
  TMAO (generic osmolyte chaperone). The first two came from routine
  papers; TMAO is included as an osmotic-stabiliser baseline.
- Targets: ultra_x_tmem145 (clinical construct), wt_full (WT decoy),
  e1659a_mutant (mutant decoy).
- Metric: Vina ΔG best of 9 modes at exhaustiveness 32, per receptor.
  Selectivity = ΔG(E1659A) - ΔG(WT). Negative = prefers mutant.

Outputs:
- JSON: pharmacochaperone_phase4g_repurpose_screen.json
- Poses + logs: docking_runs/4g/poses/*.pdbqt, logs/*.log
- Summary table printed to stdout.

Runs with strc-mmgbsa conda env.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pharmacochaperone_phase4_common import (  # noqa: E402
    TARGETS, WORK_DIR, gate_ready, print_env_report,
)
from pharmacochaperone_phase4b_vina_gnina_screen import (  # noqa: E402
    prep_ligand_pdbqt, box_centre_for_target, run_vina,
    cif_to_chain_pdb, prep_receptor_pdbqt,
)
from pharmacochaperone_phase4_common import MODELS_DIR  # noqa: E402

PHASE = "4g"
OUT_JSON = WORK_DIR / f"pharmacochaperone_phase{PHASE}_repurpose_screen.json"
RUN_DIR = WORK_DIR / "docking_runs" / PHASE

REPURPOSE_ROSTER = [
    {
        "name": "4-phenylbutyrate",
        "smiles": "OC(=O)CCCc1ccccc1",
        "role": "fda-approved-chaperone",
        "source": "hayase 2026; DSCAML1 A2105T in-vivo rescue",
        "clinical_status": "FDA approved (Buphenyl — urea cycle disorder)",
    },
    {
        "name": "ip-045",
        "smiles": "O=C(OC1=CC=CC=C1F)CCc2c[nH]c3ccccc23",
        "role": "indole-chaperone",
        "source": "kaur 2026; α-synuclein rotenone-PD rat in-vivo",
        "clinical_status": "preclinical; indole chemotype matches our lead",
    },
    {
        "name": "tmao",
        "smiles": "C[N+](C)(C)[O-]",
        "role": "osmolyte-baseline",
        "source": "reference; trimethylamine-N-oxide",
        "clinical_status": "dietary; generic protein stabiliser",
    },
]


def dock_compound_on_target(cmpd: dict, target_key: str,
                            work_dir: Path) -> dict:
    """Run ligand prep (if needed) and Vina dock on one target. Returns
    dict with score + paths, or error."""
    t = TARGETS[target_key]

    # Receptor prep (skipped if already done in 4b/4c).
    rec_pdb = work_dir / f"{target_key}_chainA.pdb"
    rec_pdbqt = work_dir / f"{target_key}_chainA.pdbqt"
    if not rec_pdbqt.exists():
        # Fall back to 4b/4c prepared receptor if present (saves time + ensures
        # exact same receptor state as earlier phases for comparability).
        for earlier in ("4b", "4c"):
            alt = WORK_DIR / "docking_runs" / earlier / f"{target_key}_chainA.pdbqt"
            if alt.exists():
                rec_pdbqt = alt
                break
        else:
            work_dir.mkdir(parents=True, exist_ok=True)
            cif_to_chain_pdb(MODELS_DIR / t["cif"], t["chain"], rec_pdb)
            prep_receptor_pdbqt(rec_pdb, rec_pdbqt)

    # Ligand prep (cached across targets — same ligand, different receptor).
    lig_dir = work_dir / "ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)
    lig_pdbqt = lig_dir / f"{cmpd['name']}.pdbqt"
    if not lig_pdbqt.exists():
        print(f"  prepping ligand {cmpd['name']}", flush=True)
        try:
            prep_ligand_pdbqt(cmpd["name"], cmpd["smiles"], lig_pdbqt)
        except Exception as e:
            return {"error": f"ligand prep failed: {e}"}

    # Box derived from this target's K1141_CA + loop centroid.
    try:
        centre, size = box_centre_for_target(target_key)
    except Exception as e:
        return {"error": f"box derivation failed: {e}"}

    # Dock.
    pose_dir = work_dir / "poses"
    pose_dir.mkdir(parents=True, exist_ok=True)
    log_dir = work_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    pose = pose_dir / f"{cmpd['name']}__{target_key}.pdbqt"
    log = log_dir / f"{cmpd['name']}__{target_key}.log"
    try:
        print(f"  docking  {cmpd['name']} on {target_key}", flush=True)
        score = run_vina(rec_pdbqt, lig_pdbqt, centre, size, pose, log)
    except Exception as e:
        return {"error": f"vina failed: {e}", "log": str(log)}

    return {
        "best_affinity_kcal_mol": round(float(score), 3),
        "pose_pdbqt": str(pose),
        "log": str(log),
        "receptor_pdbqt": str(rec_pdbqt),
        "box_centre": [round(float(x), 3) for x in centre],
        "box_size": list(size),
    }


def main() -> int:
    ok, missing = gate_ready("4b")  # same env requirements as 4b
    if not ok:
        print_env_report("4b")
        print(f"\nPhase {PHASE} blocked; install: {missing}")
        return 2

    results: dict[str, dict[str, dict]] = {}
    for cmpd in REPURPOSE_ROSTER:
        print(f"\n=== {cmpd['name']} ({cmpd['role']}) ===")
        per_target = {}
        for target_key in ("ultra_x_tmem145", "wt_full", "e1659a_mutant"):
            per_target[target_key] = dock_compound_on_target(cmpd, target_key, RUN_DIR)
        results[cmpd["name"]] = {"meta": cmpd, "targets": per_target}

    # Selectivity: MUT - WT (negative means prefers mutant)
    def score(cmpd_res: dict, target: str):
        d = cmpd_res["targets"].get(target, {})
        return d.get("best_affinity_kcal_mol")

    summary = []
    for name, r in results.items():
        ultra = score(r, "ultra_x_tmem145")
        wt = score(r, "wt_full")
        mut = score(r, "e1659a_mutant")
        sel = (mut - wt) if (mut is not None and wt is not None) else None
        summary.append({
            "name": name,
            "role": r["meta"]["role"],
            "clinical_status": r["meta"]["clinical_status"],
            "dG_ultra": ultra,
            "dG_wt": wt,
            "dG_e1659a": mut,
            "selectivity_mut_minus_wt": round(sel, 3) if sel is not None else None,
        })

    payload = {
        "phase": PHASE,
        "run_date": "2026-04-22",
        "vina_exhaustiveness": 32,
        "vina_num_modes": 9,
        "box_size_A": 18.0,
        "roster": REPURPOSE_ROSTER,
        "targets": list(TARGETS.keys()),
        "results": results,
        "summary": summary,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"\nWrote {OUT_JSON}")

    print("\n=== SUMMARY ===")
    print(f"{'compound':25s}  {'ΔG_ultra':>9s}  {'ΔG_wt':>9s}  {'ΔG_mut':>9s}  {'mut-wt':>7s}  role")
    for s in summary:
        def fmt(x): return f"{x:+.2f}" if isinstance(x, (int, float)) else "  --  "
        print(f"{s['name']:25s}  {fmt(s['dG_ultra']):>9s}  {fmt(s['dG_wt']):>9s}  "
              f"{fmt(s['dG_e1659a']):>9s}  {fmt(s['selectivity_mut_minus_wt']):>7s}  {s['role']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
