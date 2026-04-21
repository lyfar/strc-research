#!/usr/bin/env python3
"""
Phase 5 — GROMACS MD stability check on Phase 4b top-5 leads
pre-bound in the K1141 pocket on Ultra-Mini × TMEM145.

Protocol:
  Force field: AmberFF19SB + GAFF2, TIP3P water, 0.15 M NaCl.
  System: Ultra-Mini × TMEM145 + ligand (receptor pre-minimised from
  Phase 4f), box = cubic, dodecahedron, 1 nm padding.
  Minimisation: 10000 steps steepest descent.
  NVT 100 ps, NPT 100 ps, production 50 ns × 3 replicates per ligand.

Per-ligand gates (averaged over 3 replicates, last 40 ns):
  - ligand heavy-atom RMSD <3 A
  - K1141-NZ <-> ligand acid distance <3.5 A in >60% of frames
  - loop 1642-1651 RMSD to WT Cα trajectory <2 A
  - MM-PBSA ΔG_bind <= -6 kcal/mol

Phase advance if >=1 of top-5 leads passes all four sub-gates.

Runtime: ~48 h on A100 for 24 trajectories (8 systems × 3 replicates
× 50 ns). Local GPU fallback: 20 ns × 1 replicate on top-3 only.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pharmacochaperone_phase4_common import (  # noqa: E402
    gate_ready, print_env_report, WORK_DIR,
)

PHASE = "5"
OUT_JSON = WORK_DIR / "pharmacochaperone_phase5_md.json"
DEP_JSON = WORK_DIR / "pharmacochaperone_phase4b_vina_gnina_screen.json"


def main():
    print_env_report(PHASE)
    ok, missing = gate_ready(PHASE)
    if not ok:
        print(f"Phase {PHASE} BLOCKED: missing {missing} "
              f"(local gmx OR A100 rental ~$100-500).")
        return 2
    if not DEP_JSON.exists():
        print(f"Phase {PHASE} waiting for Phase 4b output at {DEP_JSON}")
        return 3

    # --- TODO (after 4b/4f complete) ------------------------------
    # 1) Prep top-5 ligand-receptor systems with 4f minimised geom.
    # 2) Solvate + ionise with gmx pdb2gmx + solvate + genion.
    # 3) Minimise + NVT + NPT + 50 ns production × 3 replicates.
    # 4) gmx rms + gmx distance + MM-PBSA (gmx_MMPBSA).
    # 5) Write per-ligand gate verdicts to OUT_JSON.
    # --------------------------------------------------------------
    print(f"Phase {PHASE} scaffold present; awaiting Phase 4b + gmx install.")
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
