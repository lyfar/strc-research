#!/usr/bin/env python3
"""
Phase 4d — K1141A in-silico decoy.

Build STRC(E1659A, K1141A) by swapping Lys1141 sidechain to Ala on the
clinical Ultra-Mini × TMEM145 CIF. Re-dock the Phase 4b top-20 +
FIXED_ROSTER against the K1141A mutant pocket.

Gate: carboxylate leads (salicylic, nicotinic, indole-3-acetic,
cyclopropane-phenyl-COOH, naphthalene-2-COOH) lose >=2 kcal/mol vs.
the K1141-intact docking. If they don't, K1141 is not the load-bearing
salt-bridge partner the Phase 3A pharmacophore claims, and the hypothesis
is falsified.

Mutagenesis: PyMOL `mutate 1141, ALA` + brief local minimisation with
RDKit+MMFF or openmm. No repacking of distal residues (would alter the
pocket independently of the K1141 substitution).
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pharmacochaperone_phase4_common import (  # noqa: E402
    gate_ready, print_env_report, TARGETS, WORK_DIR,
)

PHASE = "4d"
OUT_JSON = WORK_DIR / "pharmacochaperone_phase4d_k1141a_decoy.json"
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

    # --- TODO (after 4b completes) --------------------------------
    # 1) Build K1141A mutant PDB from Ultra-Mini × TMEM145 CIF
    #    (chain A, CIF residue 1141-1074=67 -> mutate to ALA).
    # 2) Prep receptor PDBQT; reuse Phase 4b box centre.
    # 3) Dock Phase 4b top-20 + FIXED_ROSTER with same Vina settings.
    # 4) Per-lead score loss = ΔG_K1141A - ΔG_K1141_intact.
    # 5) Gate: >=3 of top-5 carboxylate leads show score loss >=2 kcal/mol.
    # --------------------------------------------------------------
    print(f"Phase {PHASE} scaffold present; awaiting Phase 4b output.")
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
