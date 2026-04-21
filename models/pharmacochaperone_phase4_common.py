#!/usr/bin/env python3
"""
Shared utilities for Phase 4b-4f pharmacochaperone docking runs.

- Environment probe (Vina / Smina / GNINA / Meeko / Open Babel / RDKit).
- Docking target registry with per-CIF residue-number offsets.
- Fixed compound roster (top-5 Phase 3C leads + diflunisal positive +
  three polar non-binder negatives). Every Phase 4b-4f run MUST dock
  this roster so each gate is comparable to every other.

Not runnable on its own — imported by phase4b/c/d/e/f driver scripts.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

MODELS_DIR = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models")
WORK_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
DOCKING_DIR = WORK_DIR / "docking_runs"

# --- Fixed compound roster --------------------------------------------------
# Every Phase 4 docking pass MUST dock these. Additional library compounds
# (DrugBank FDA / DSi-Poised / ZINC22 carboxylate tranche) are appended by
# each phase's driver. Keeping this roster fixed makes the gates comparable.
FIXED_ROSTER = [
    # --- top-5 leads from Phase 3C shape-fit ranking ---
    {"name": "salicylic-acid",             "smiles": "OC(=O)c1ccccc1O",                "role": "lead"},
    {"name": "nicotinic-acid",             "smiles": "OC(=O)c1cccnc1",                 "role": "lead"},
    {"name": "cyclopropane-phenyl-COOH",   "smiles": "OC(=O)C1CC1c1ccccc1",            "role": "lead"},
    {"name": "indole-3-acetic-acid",       "smiles": "OC(=O)Cc1c[nH]c2ccccc12",        "role": "lead"},
    {"name": "naphthalene-2-carboxylic",   "smiles": "OC(=O)c1ccc2ccccc2c1",           "role": "lead"},
    # --- positive control: TTR pharmacochaperone ---
    {"name": "diflunisal",                 "smiles": "OC(=O)c1cc(-c2ccc(F)cc2F)ccc1O", "role": "positive"},
    # --- negative controls: polar, no aromatic, no carboxylate ---
    {"name": "glucose",                    "smiles": "OCC1OC(O)C(O)C(O)C1O",           "role": "negative"},
    {"name": "urea",                       "smiles": "NC(=O)N",                        "role": "negative"},
    {"name": "acetamide",                  "smiles": "CC(=O)N",                        "role": "negative"},
]

# --- Docking target CIFs + residue-number offsets ---------------------------
# real_resnum = cif_resnum + offset ;  cif_resnum = real_resnum - offset
TARGETS = {
    "ultra_x_tmem145": {
        "cif":    "job-ultramini-x-tmem145-full.cif",
        "offset": 1074,
        "chain":  "A",
        "note":   "clinical construct; Phase 4b primary target",
    },
    "wt_full": {
        "cif":    "job4-wildtype.cif",
        "offset": 0,
        "chain":  "A",
        "note":   "Phase 4c WT decoy target",
    },
    "e1659a_mutant": {
        "cif":    "job3-mutant.cif",
        "offset": 0,
        "chain":  "A",
        "note":   "Phase 4d backbone for in-silico K1141A mutant build",
    },
}

# Docking box — derived from Phase 2B subpocket scan. Each driver rebuilds
# the box in its CIF's own frame by measuring from local K1141 Cα; the
# size parameter below is universal.
BOX_SIZE_A = (18.0, 18.0, 18.0)

# Real STRC residues referenced by all phases.
REAL_K1141 = 1141
REAL_D1140 = 1140
REAL_D1173 = 1173
REAL_F1646 = 1646
REAL_E1659 = 1659
REAL_LOOP = (1642, 1651)


@dataclass
class EnvStatus:
    rdkit:   bool
    biopdb:  bool
    openbabel_cli: bool
    vina_cli:      bool
    vina_py:       bool
    smina_cli:     bool
    gnina_cli:     bool
    meeko:         bool
    fpocket_cli:   bool
    adfr_prepare_receptor: bool
    gmx_cli:       bool

    def missing_for(self, phase: str) -> list[str]:
        """Return list of missing tools for a given phase. Empty = ready."""
        need = {
            "4b": ["rdkit", "biopdb", "vina_cli or vina_py", "meeko or openbabel_cli"],
            "4c": ["rdkit", "biopdb", "vina_cli or vina_py", "meeko or openbabel_cli"],
            "4d": ["rdkit", "biopdb", "vina_cli or vina_py", "meeko or openbabel_cli"],
            "4e": ["rdkit", "biopdb"],  # pure geometry + shape fit — runnable today
            "4f": ["rdkit", "biopdb", "vina_cli or vina_py"],
            "5":  ["gmx_cli"],
        }
        want = need[phase]
        missing = []
        for w in want:
            options = [o.strip() for o in w.split("or")]
            if not any(getattr(self, o) for o in options if hasattr(self, o)):
                missing.append(w)
        return missing


def probe_env() -> EnvStatus:
    def mod(name: str) -> bool:
        try:
            __import__(name)
            return True
        except Exception:
            return False

    def cli(name: str) -> bool:
        return shutil.which(name) is not None

    return EnvStatus(
        rdkit=mod("rdkit"),
        biopdb=mod("Bio.PDB"),
        openbabel_cli=cli("obabel"),
        vina_cli=cli("vina"),
        vina_py=mod("vina"),
        smina_cli=cli("smina"),
        gnina_cli=cli("gnina"),
        meeko=mod("meeko"),
        fpocket_cli=cli("fpocket"),
        adfr_prepare_receptor=cli("prepare_receptor"),
        gmx_cli=cli("gmx"),
    )


INSTALL_HINTS = {
    "rdkit":        "pip install rdkit",
    "biopdb":       "pip install biopython",
    "openbabel_cli":"brew install open-babel",
    "vina_cli":     "brew install autodock-vina",
    "vina_py":      "pip install vina",
    "smina_cli":    "brew install smina",
    "gnina_cli":    "Download from https://github.com/gnina/gnina/releases (macOS arm64 build)",
    "meeko":        "pip install meeko",
    "fpocket_cli":  "brew install fpocket",
    "adfr_prepare_receptor": "Install AutoDock ADFR suite from https://ccsb.scripps.edu/adfr/downloads/",
    "gmx_cli":      "brew install gromacs   (local MD; or rent A100 and run cloud-side)",
}


def gate_ready(phase: str) -> tuple[bool, list[str]]:
    env = probe_env()
    missing = env.missing_for(phase)
    return (not missing, missing)


def prepare_work_dir(phase: str) -> Path:
    p = DOCKING_DIR / phase
    p.mkdir(parents=True, exist_ok=True)
    return p


def print_env_report(phase: str) -> None:
    env = probe_env()
    print(f"phase {phase} environment probe:")
    for k, v in asdict(env).items():
        mark = "OK" if v else "--"
        hint = "" if v else f"   [install: {INSTALL_HINTS.get(k,'?')}]"
        print(f"  [{mark}] {k}{hint}")
    missing = env.missing_for(phase)
    if missing:
        print(f"  PHASE {phase} BLOCKED on: {', '.join(missing)}")
    else:
        print(f"  PHASE {phase} READY")


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


if __name__ == "__main__":
    for phase in ("4b", "4c", "4d", "4e", "4f", "5"):
        print_env_report(phase)
        print()
