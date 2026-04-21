#!/usr/bin/env python3
"""
Phase 4b — Vina docking against the K1141 pocket on Ultra-Mini ×
TMEM145 (clinical construct).

This version runs the fixed 9-compound roster (Phase 3C top-5 leads +
diflunisal positive + 3 polar negative controls) as a smoke test.
Library expansion to DrugBank FDA / DSi-Poised / ZINC22 carboxylate
tranche is a separate run invoked via `--library <name>`.

Pipeline (single-binary, no GNINA yet — rescore is Phase 4b-pt2):
  1. Receptor prep: CIF → chain A only PDB → obabel -xr -h PDBQT
     (rigid receptor, Gasteiger charges, polar H).
  2. Box: 18×18×18 Å at K1141_CA + 3 Å toward loop-1642-1651 centroid
     (reference-frame-free, derived per-CIF).
  3. Ligand prep per SMILES:
       RDKit ETKDGv3 (30 conformers) + MMFF94s minimise → SDF,
       obabel -p 7.4 protonation,
       meeko.MoleculePreparation → PDBQT.
  4. Dock: vina --exhaustiveness 32 --num_modes 9 --cpu 8
     per ligand; record best score.
  5. Emit JSON with per-compound score, pose file paths, gate verdict
     against diflunisal (positive control must score ≤ -5 kcal/mol
     AND ≥3 leads score ≤ diflunisal score).
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pharmacochaperone_phase4_common import (  # noqa: E402
    gate_ready, print_env_report, TARGETS, FIXED_ROSTER,
    MODELS_DIR, WORK_DIR, BOX_SIZE_A,
    REAL_K1141, REAL_LOOP,
)
from pharmacochaperone_phase2b_subpockets import (  # noqa: E402
    load_chain, residue_ca, loop_centroid,
)
from Bio.PDB import MMCIFParser, PDBIO, Select  # noqa: E402
from Bio.PDB.Polypeptide import is_aa  # noqa: E402

PHASE = "4b"
OUT_JSON = WORK_DIR / f"pharmacochaperone_phase{PHASE}_vina_gnina_screen.json"
RUN_DIR = WORK_DIR / "docking_runs" / PHASE
VINA_EXHAUSTIVENESS = 32
VINA_NUM_MODES = 9
VINA_CPU = 8


class StandardAAOnly(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, res):
        # Standard AA, skip waters, heteroatoms, ligands
        return is_aa(res, standard=True)


def cif_to_chain_pdb(cif_path: Path, chain_id: str, out_pdb: Path) -> None:
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", cif_path)
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out_pdb), select=StandardAAOnly(chain_id))


def prep_receptor_pdbqt(pdb_path: Path, pdbqt_path: Path) -> None:
    # obabel handles protonation (-p 7.4), Gasteiger charges, and PDBQT
    # writing with --partialcharge gasteiger -xr (rigid receptor).
    cmd = [
        "obabel", str(pdb_path), "-O", str(pdbqt_path),
        "-xr",                       # rigid receptor
        "--partialcharge", "gasteiger",
        "-p", "7.4",                 # add H at pH 7.4
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not pdbqt_path.exists():
        raise RuntimeError(f"receptor prep failed: {r.stderr}")


def box_centre_for_target(target_key: str) -> tuple[np.ndarray, tuple[float, float, float]]:
    t = TARGETS[target_key]
    chain = load_chain(MODELS_DIR / t["cif"], t["chain"])
    off = t["offset"]
    k_ca = residue_ca(chain, REAL_K1141 - off)
    loop_c = loop_centroid(chain, REAL_LOOP[0] - off, REAL_LOOP[1] - off)
    vec = loop_c - k_ca
    vec = vec / np.linalg.norm(vec)
    centre = k_ca + 3.0 * vec
    return centre, BOX_SIZE_A


def prep_ligand_pdbqt(name: str, smiles: str, out_pdbqt: Path) -> None:
    """RDKit ETKDG + MMFF94s → SDF → obabel pH 7.4 → meeko PDBQT."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"bad SMILES for {name}: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    confs = AllChem.EmbedMultipleConfs(mol, numConfs=30, params=params)
    if len(confs) == 0:
        # Fallback to single-shot embed + 2D coord generation.
        if AllChem.EmbedMolecule(mol, params) != 0:
            AllChem.Compute2DCoords(mol)
    # Minimise every conformer, keep lowest E.
    energies = []
    for cid in mol.GetConformers():
        try:
            ff = AllChem.MMFFGetMoleculeForceField(
                mol, AllChem.MMFFGetMoleculeProperties(mol),
                confId=cid.GetId()
            )
            if ff is None:
                energies.append((cid.GetId(), float("inf")))
                continue
            ff.Minimize(maxIts=500)
            energies.append((cid.GetId(), ff.CalcEnergy()))
        except Exception:
            energies.append((cid.GetId(), float("inf")))
    energies.sort(key=lambda x: x[1])
    best_cid = energies[0][0] if energies else 0

    sdf = out_pdbqt.with_suffix(".rdkit.sdf")
    w = Chem.SDWriter(str(sdf))
    w.write(mol, confId=best_cid)
    w.close()

    # obabel for pH 7.4 reprotonation, then meeko for PDBQT.
    sdf_pH = out_pdbqt.with_suffix(".pH74.sdf")
    r = subprocess.run(
        ["obabel", str(sdf), "-O", str(sdf_pH), "-p", "7.4"],
        capture_output=True, text=True,
    )
    if r.returncode != 0 or not sdf_pH.exists():
        raise RuntimeError(f"obabel pH fail for {name}: {r.stderr}")

    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mols = list(Chem.SDMolSupplier(str(sdf_pH), removeHs=False))
    if not mols or mols[0] is None:
        raise RuntimeError(f"rdkit couldn't re-read pH-adjusted SDF for {name}")
    prep = MoleculePreparation()
    setups = prep.prepare(mols[0])
    pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        raise RuntimeError(f"meeko PDBQT fail for {name}: {err}")
    out_pdbqt.write_text(pdbqt_str)


def run_vina(receptor_pdbqt: Path, ligand_pdbqt: Path,
             centre: np.ndarray, size: tuple[float, float, float],
             out_pose: Path, log_path: Path) -> float:
    """Run Vina; return best (most-negative) binding-affinity kcal/mol."""
    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand",   str(ligand_pdbqt),
        "--center_x", f"{centre[0]:.3f}",
        "--center_y", f"{centre[1]:.3f}",
        "--center_z", f"{centre[2]:.3f}",
        "--size_x",   f"{size[0]:.1f}",
        "--size_y",   f"{size[1]:.1f}",
        "--size_z",   f"{size[2]:.1f}",
        "--exhaustiveness", str(VINA_EXHAUSTIVENESS),
        "--num_modes", str(VINA_NUM_MODES),
        "--cpu", str(VINA_CPU),
        "--out", str(out_pose),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    log_path.write_text(r.stdout + "\n---STDERR---\n" + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"vina exit {r.returncode}:\n{r.stderr[:800]}")
    # Parse best score from output PDBQT (first "REMARK VINA RESULT" line).
    best = float("inf")
    for line in out_pose.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT"):
            parts = line.split()
            best = float(parts[3])
            break
    return best


def run_roster(target_key: str, work_dir: Path) -> dict:
    centre, size = box_centre_for_target(target_key)
    t = TARGETS[target_key]
    print(f"\nPhase 4b target: {target_key} ({t['cif']})")
    print(f"  box centre = {[round(float(x),2) for x in centre]}  size = {size}")
    work_dir.mkdir(parents=True, exist_ok=True)

    # Receptor prep (once).
    rec_pdb = work_dir / f"{target_key}_chainA.pdb"
    rec_pdbqt = work_dir / f"{target_key}_chainA.pdbqt"
    if not rec_pdbqt.exists():
        print("  preparing receptor PDBQT...", flush=True)
        cif_to_chain_pdb(MODELS_DIR / t["cif"], t["chain"], rec_pdb)
        prep_receptor_pdbqt(rec_pdb, rec_pdbqt)
        print(f"    wrote {rec_pdbqt.name} ({rec_pdbqt.stat().st_size} B)")
    else:
        print(f"  receptor already prepared: {rec_pdbqt.name}")

    # Dock each roster compound.
    results = []
    for cmpd in FIXED_ROSTER:
        lig_dir = work_dir / "ligands"
        lig_dir.mkdir(exist_ok=True)
        lig_pdbqt = lig_dir / f"{cmpd['name']}.pdbqt"
        try:
            if not lig_pdbqt.exists():
                print(f"  prepping ligand {cmpd['name']:32s}", flush=True)
                prep_ligand_pdbqt(cmpd["name"], cmpd["smiles"], lig_pdbqt)
        except Exception as e:
            print(f"  LIGAND PREP FAILED for {cmpd['name']}: {e}")
            results.append({
                **cmpd, "best_affinity_kcal_mol": None, "error": str(e),
            })
            continue

        pose_pdbqt = work_dir / "poses" / f"{cmpd['name']}__{target_key}.pdbqt"
        pose_pdbqt.parent.mkdir(exist_ok=True)
        log_path = work_dir / "logs" / f"{cmpd['name']}__{target_key}.log"
        log_path.parent.mkdir(exist_ok=True)
        try:
            print(f"  docking  {cmpd['name']:32s}", flush=True)
            score = run_vina(rec_pdbqt, lig_pdbqt, centre, size, pose_pdbqt, log_path)
            print(f"    ΔG_best = {score:+.2f} kcal/mol")
            results.append({
                **cmpd, "best_affinity_kcal_mol": round(float(score), 3),
                "pose_pdbqt": str(pose_pdbqt), "log": str(log_path),
            })
        except Exception as e:
            print(f"  VINA FAILED for {cmpd['name']}: {e}")
            results.append({
                **cmpd, "best_affinity_kcal_mol": None, "error": str(e),
                "log": str(log_path) if log_path.exists() else None,
            })

    # Gate: diflunisal score ≤ -5 kcal/mol AND ≥3 leads ≤ diflunisal.
    diflu = next(
        (r for r in results if r["name"] == "diflunisal"
         and r.get("best_affinity_kcal_mol") is not None),
        None,
    )
    leads = [r for r in results if r["role"] == "lead"
             and r.get("best_affinity_kcal_mol") is not None]
    negs  = [r for r in results if r["role"] == "negative"
             and r.get("best_affinity_kcal_mol") is not None]

    gate = {
        "positive_control_present": diflu is not None,
        "positive_score": diflu["best_affinity_kcal_mol"] if diflu else None,
        "positive_threshold_kcal_mol": -5.0,
        "positive_passes": bool(diflu and diflu["best_affinity_kcal_mol"] <= -5.0),
        "leads_beating_positive": sum(
            1 for l in leads
            if diflu and l["best_affinity_kcal_mol"] <= diflu["best_affinity_kcal_mol"]
        ),
        "min_leads_to_advance": 3,
    }
    gate["verdict"] = (
        "PASS" if gate["positive_passes"] and gate["leads_beating_positive"] >= 3
        else "FAIL"
    )

    payload = {
        "phase": PHASE,
        "target": target_key,
        "box_centre": [round(float(x), 3) for x in centre],
        "box_size_A": list(size),
        "vina_exhaustiveness": VINA_EXHAUSTIVENESS,
        "vina_num_modes": VINA_NUM_MODES,
        "gate": gate,
        "results": results,
        "ranked": sorted(
            [r for r in results if r.get("best_affinity_kcal_mol") is not None],
            key=lambda r: r["best_affinity_kcal_mol"],
        ),
    }
    return payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="ultra_x_tmem145",
                    choices=list(TARGETS.keys()))
    args = ap.parse_args()

    print_env_report(PHASE)
    ok, missing = gate_ready(PHASE)
    # openbabel_cli + meeko cover the meeko-or-obabel logic; rdkit + biopdb
    # + vina_cli are required.
    if not ok:
        print(f"BLOCKED: {missing}")
        return 2

    payload = run_roster(args.target, RUN_DIR)
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print()
    print(f"verdict: {payload['gate']['verdict']}")
    print(f"ranked:")
    for r in payload["ranked"]:
        print(f"  {r['role']:8s}  {r['name']:32s}  {r['best_affinity_kcal_mol']:+6.2f} kcal/mol")
    print(f"wrote {OUT_JSON}")
    return 0 if payload["gate"]["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
