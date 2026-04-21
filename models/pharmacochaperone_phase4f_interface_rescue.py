#!/usr/bin/env python3
"""
Phase 4f — Interface-rescue MM-GBSA on Ultra-Mini × TMEM145 complex.

This is the single gate that tests the THERAPEUTIC claim (vs. mere
pocket binding): does a pharmacochaperone pre-bound in the K1141 pocket
improve the Ultra-Mini × TMEM145 interface energy toward WT levels?

Method:
  1. Take Ultra-Mini × TMEM145 complex CIF (AF3 job, ipTM 0.43, 23/41
     residues in GOLD zone aa 1603-1749).
  2. For each of the Phase 4b top-3 leads:
       a. Dock the lead into the K1141 pocket on chain A of the complex
          (reuse Phase 4b pose with best Vina score).
       b. Build receptor-ligand-TMEM145 system (ligand as cofactor).
       c. Minimise briefly (MMFF94s or AmberFF14SB + GAFF2 if available).
       d. Compute MM-GBSA on the A↔B chain interface with ligand bound.
       e. Compare to apo MM-GBSA on the same complex (no ligand).
  3. Rescue recovery = (ΔG_interface_apo - ΔG_interface_with_ligand)
     / 8.4 kcal/mol (the E1659A-induced binding-energy loss, from
     [[STRC Electrostatic Analysis E1659A]]).

Gate: >=1 of top-3 leads recovers >=30% of the 8.4 kcal/mol gap.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pharmacochaperone_phase4_common import (  # noqa: E402
    gate_ready, print_env_report, TARGETS, WORK_DIR,
)

PHASE = "4f"
OUT_JSON = WORK_DIR / "pharmacochaperone_phase4f_interface_rescue.json"
DEP_JSON = WORK_DIR / "pharmacochaperone_phase4b_vina_gnina_screen.json"


def main():
    print_env_report(PHASE)
    ok, missing = gate_ready(PHASE)
    if not ok:
        print(f"Phase {PHASE} BLOCKED: missing {missing}")
        return 2
    if not DEP_JSON.exists():
        print(f"Phase {PHASE} waiting for Phase 4b output at {DEP_JSON}")
        return 3

    # --- TODO -----------------------------------------------------
    # 1) Load Phase 4b top-3 leads + Vina best poses.
    # 2) Splice ligand into Ultra-Mini × TMEM145 complex (chain A
    #    pocket, chain B untouched).
    # 3) MMFF94s minimise + AmberFF14SB/GAFF2 fallback if available.
    # 4) MM-GBSA on chain A↔chain B interface per system (apo + 3
    #    ligand-bound).
    # 5) Compute rescue_recovery = (apo - bound) / 8.4 kcal/mol.
    # 6) Gate: >=1 of top-3 leads shows rescue_recovery >= 0.30.
    # --------------------------------------------------------------
    print(f"Phase {PHASE} scaffold present; awaiting Phase 4b output.")
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
