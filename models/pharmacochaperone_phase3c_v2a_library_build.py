#!/usr/bin/env python3
"""
Phase 3c v2a — Expanded screening library for K1141 pocket.

Rationale: Phase 4b shortlist was 5 aromatic carboxylates + diflunisal; Phase 5b
RED-LIGHT (f_PC < 0.10) means the chemical space around simple aryl-COOH is
insufficient. This phase expands chemical space along three axes:

  1. Acid-group bioisosterism: COOH → tetrazole, acylsulfonamide, hydroxamic,
     phosphonate, CH2COOH (extended linker)
  2. Aromatic/heteroaromatic scaffold diversity: 25 cores (benzene,
     naphthalene, biphenyl, quinoline, indole, benzofuran, benzimidazole,
     thiophene-benzene fused, etc.)
  3. FDA-approved carboxylate drugs: 30 curated (including Phase 4b leads +
     established chaperone chemotypes like tafamidis analogs, TTR stabilizers)

Library is filtered for drug-like constraints (Lipinski with extended MW for
pharmacochaperone use case) + must contain ≥1 aromatic ring + ≥1 acidic group.

Output: `phase3c_v2a_library/*.pdbqt` + `phase3c_v2a_library.json` manifest.

Deps: rdkit, obabel. Runtime ~1 min on local Mac.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
OUT_DIR = WORK / "docking_runs" / "3c_v2" / "ligands"
OUT_JSON = WORK / "pharmacochaperone_phase3c_v2a_library_build.json"

# Drug-like window tuned for pharmacochaperone (pocket Kd target ≤ 5 μM)
MW_MIN, MW_MAX = 150.0, 450.0
LOGP_MIN, LOGP_MAX = -1.5, 4.5
HBD_MAX = 4
HBA_MAX = 8
ROTB_MAX = 7

# ---------------------------------------------------------------------------
# Scaffold SMILES (aromatic cores, leaving one H on a site where acid attaches)
# Use [*H] marker carbons where we'll attach
SCAFFOLDS_SMARTS = [
    # Simple aromatics
    ("benzene",                 "c1ccccc1"),
    ("toluene",                 "Cc1ccccc1"),
    ("phenol",                  "Oc1ccccc1"),
    ("anisole",                 "COc1ccccc1"),
    ("aniline",                 "Nc1ccccc1"),
    ("fluorobenzene",           "Fc1ccccc1"),
    ("chlorobenzene",           "Clc1ccccc1"),
    ("naphthalene",             "c1ccc2ccccc2c1"),
    ("biphenyl",                "c1ccc(-c2ccccc2)cc1"),
    # Heteroaromatics
    ("pyridine",                "c1ccncc1"),
    ("pyrimidine",              "c1cncnc1"),
    ("pyrazine",                "c1cnccn1"),
    ("thiophene",               "c1ccsc1"),
    ("furan",                   "c1ccoc1"),
    ("oxazole",                 "c1ocnc1"),
    ("imidazole",               "c1[nH]cnc1"),
    # Fused bicyclics
    ("quinoline",               "c1ccc2ncccc2c1"),
    ("isoquinoline",            "c1ccc2cnccc2c1"),
    ("indole",                  "c1ccc2[nH]ccc2c1"),
    ("benzimidazole",           "c1ccc2[nH]cnc2c1"),
    ("benzoxazole",             "c1ccc2ocnc2c1"),
    ("benzothiazole",           "c1ccc2scnc2c1"),
    ("benzofuran",              "c1ccc2occc2c1"),
    ("benzothiophene",          "c1ccc2sccc2c1"),
    ("indoline-2-one",          "O=C1Nc2ccccc2C1"),
]

# Acidic / acid-bioisostere appendage SMILES, attached via the leading atom
# (the acid group's pharmacophoric anchor)
ACID_GROUPS = [
    ("COOH",           "C(=O)O"),
    ("CH2COOH",        "CC(=O)O"),
    ("OCH2COOH",       "OCC(=O)O"),
    ("CONHOH",         "C(=O)NO"),        # hydroxamic
    ("SO2NH2",         "S(=O)(=O)N"),     # sulfonamide (weak acidic)
    ("CONHSO2Me",      "C(=O)NS(=O)(=O)C"),  # acylsulfonamide
    ("PO3H2",          "P(=O)(O)O"),      # phosphonate
    ("tetrazole",      "c1nnn[nH]1"),    # 1H-tetrazole, pKa ~4.9
    ("CH=CHCOOH",      "/C=C/C(=O)O"),   # acrylic
    ("NHCOCOOH",       "NC(=O)C(=O)O"),  # oxamic
]

# FDA-approved / literature-precedent carboxylate drugs and known chaperones.
# Includes Phase 4b leads as positive controls to confirm library-build parity.
CURATED_DRUGS = [
    # Phase 4b leads (positive controls)
    ("indole-3-acetic-acid",        "O=C(O)Cc1c[nH]c2ccccc12"),
    ("naphthalene-2-carboxylic",    "O=C(O)c1ccc2ccccc2c1"),
    ("cyclopropane-phenyl-COOH",    "O=C(O)C1(CC1)c1ccccc1"),
    ("salicylic-acid",              "O=C(O)c1ccccc1O"),
    ("nicotinic-acid",              "O=C(O)c1cccnc1"),
    ("diflunisal",                  "O=C(O)c1cc(-c2ccc(F)cc2F)ccc1O"),
    # TTR stabilizers and chaperones
    ("tafamidis-analog",            "O=C(O)c1cnc2ccc(Oc3ccccc3)cc2c1"),
    ("tolcapone",                   "O=C(c1cccc([N+](=O)[O-])c1)c1ccc(C)c(O)c1O"),
    ("entacapone-like",             "O=C(O)CC=C(C#N)c1ccc(O)c(O)c1"),
    # NSAIDs with COOH
    ("ibuprofen",                   "CC(C(=O)O)c1ccc(CC(C)C)cc1"),
    ("naproxen",                    "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
    ("diclofenac",                  "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    ("indomethacin",                "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1"),
    ("sulindac",                    "Cc1c(CC(=O)O)c2cc(F)ccc2/C1=C\\c1ccc(S(C)=O)cc1"),
    ("etodolac",                    "CCc1[nH]c2c(c1CCC(C)(O)C1=CC(=O)O)cccc2"),
    ("ketoprofen",                  "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1"),
    ("flurbiprofen",                "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1"),
    ("fenoprofen",                  "CC(C(=O)O)c1cccc(Oc2ccccc2)c1"),
    ("mefenamic-acid",              "Cc1cccc(Nc2ccccc2C(=O)O)c1C"),
    ("flufenamic-acid",             "O=C(O)c1ccccc1Nc1cccc(C(F)(F)F)c1"),
    ("niflumic-acid",               "O=C(O)c1cccnc1Nc1cccc(C(F)(F)F)c1"),
    ("tolfenamic-acid",             "Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1"),
    ("meclofenamic-acid",           "Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl"),
    # Other acid-bearing drugs
    ("probenecid",                  "CCCN(CCC)S(=O)(=O)c1ccc(C(=O)O)cc1"),
    ("furosemide",                  "O=C(O)c1cc(S(=O)(=O)N)cc(Cl)c1NCc1ccco1"),
    ("ethacrynic-acid",             "CCC(=C)C(=O)c1ccc(OCC(=O)O)c(Cl)c1Cl"),
    ("bumetanide",                  "CCCCNc1cc(C(=O)O)c(S(=O)(=O)N)cc1Oc1ccccc1"),
    ("valproic-acid",               "CCCC(CCC)C(=O)O"),   # control, non-aromatic
    ("cromolyn-half",               "O=C(O)c1cc(=O)c2cc(O)ccc2o1"),
    ("chlorambucil",                "OC(=O)CCCc1ccc(N(CCCl)CCCl)cc1"),
    ("sulfasalazine",               "O=C(O)c1cc(/N=N/c2ccc(S(=O)(=O)Nc3ncccn3)cc2)ccc1O"),
    ("olsalazine",                  "O=C(O)c1cc(/N=N/c2cc(C(=O)O)c(O)cc2)ccc1O"),
]


def enumerate_combinatorial() -> list[tuple[str, str]]:
    """For each (scaffold, acid_group), generate ortho-meta-para substitution
    patterns by attaching acid_group to each H-bearing aromatic carbon.
    Unique canonical SMILES kept."""
    seen: set[str] = set()
    out: list[tuple[str, str]] = []
    for sc_name, sc_smiles in SCAFFOLDS_SMARTS:
        sc_mol = Chem.MolFromSmiles(sc_smiles)
        if sc_mol is None:
            continue
        for ac_name, ac_smiles in ACID_GROUPS:
            # Generate substitution by hydrogens on aromatic carbons
            # Pick each aromatic H-bearing atom and build a new molecule
            for atom_idx, atom in enumerate(sc_mol.GetAtoms()):
                if not atom.GetIsAromatic():
                    continue
                if atom.GetTotalNumHs() < 1:
                    continue
                if atom.GetSymbol() not in ("C", "N"):
                    continue
                if atom.GetSymbol() == "N" and atom.GetFormalCharge() != 0:
                    continue
                try:
                    # Build SMILES by atom-map insertion
                    new_smiles = _attach(sc_smiles, atom_idx, ac_smiles)
                    if new_smiles is None:
                        continue
                    mol = Chem.MolFromSmiles(new_smiles)
                    if mol is None:
                        continue
                    canon = Chem.MolToSmiles(mol)
                    if canon in seen:
                        continue
                    seen.add(canon)
                    out.append((f"{sc_name}__{ac_name}__{atom_idx}", canon))
                except Exception:
                    continue
    return out


def _attach(scaffold_smiles: str, atom_idx: int, acid_smiles: str) -> str | None:
    """Build a combined molecule by attaching acid_smiles as a branch at
    atom_idx of scaffold.
    """
    scaffold = Chem.MolFromSmiles(scaffold_smiles)
    acid = Chem.MolFromSmiles(acid_smiles)
    if scaffold is None or acid is None:
        return None
    rw = Chem.RWMol(Chem.CombineMols(scaffold, acid))
    # Atoms of acid are appended after scaffold atoms
    acid_start = scaffold.GetNumAtoms()
    # Find first heavy atom of acid (the pharmacophoric anchor)
    acid_anchor = acid_start + 0
    # Bond scaffold atom_idx → acid_anchor
    try:
        rw.AddBond(atom_idx, acid_anchor, Chem.BondType.SINGLE)
    except Exception:
        return None
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol)


def drug_like(mol: Chem.Mol) -> bool:
    if mol is None:
        return False
    mw = Descriptors.MolWt(mol)
    if not (MW_MIN <= mw <= MW_MAX):
        return False
    logp = Descriptors.MolLogP(mol)
    if not (LOGP_MIN <= logp <= LOGP_MAX):
        return False
    hbd = Lipinski.NumHDonors(mol)
    if hbd > HBD_MAX:
        return False
    hba = Lipinski.NumHAcceptors(mol)
    if hba > HBA_MAX:
        return False
    rotb = Lipinski.NumRotatableBonds(mol)
    if rotb > ROTB_MAX:
        return False
    if rdMolDescriptors.CalcNumAromaticRings(mol) < 1:
        return False
    # Has acidic group
    acid_smarts = [
        "C(=O)[OH]",             # COOH
        "S(=O)(=O)[NH]",         # acylsulfonamide / sulfonamide
        "C(=O)N[OH]",            # hydroxamic
        "P(=O)([OH])[OH]",       # phosphonate
        "c1nnn[nH]1",            # tetrazole
    ]
    has_acid = any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in acid_smarts)
    return has_acid


def to_pdbqt(smiles: str, name: str, out_dir: Path) -> Path | None:
    """Embed 3D conformer, minimize briefly, write SDF, convert to PDBQT."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    try:
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if embed_result == -1:
            return None
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        return None
    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        sdf = tmp / f"{name}.sdf"
        w = Chem.SDWriter(str(sdf))
        w.write(mol)
        w.close()
        pdbqt = out_dir / f"{name}.pdbqt"
        r = subprocess.run(
            ["obabel", str(sdf), "-O", str(pdbqt),
             "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=60,
        )
        if r.returncode != 0 or not pdbqt.exists():
            return None
        return pdbqt


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


def sanitize_name(name: str) -> str:
    # Safe for filesystem + Vina log parsing
    return re.sub(r"[^A-Za-z0-9_\-]", "_", name)[:60]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    log("=== Phase 3c v2a library build ===")
    # 1. Combinatorial enumeration
    log("Enumerating scaffold × acid combinatorial...")
    combinatorial = enumerate_combinatorial()
    log(f"  raw combinatorial: {len(combinatorial)}")

    # 2. Filter for drug-likeness + acid presence
    filtered: list[tuple[str, str]] = []
    seen_canon = set()
    for name, smiles in combinatorial:
        mol = Chem.MolFromSmiles(smiles)
        if not drug_like(mol):
            continue
        canon = Chem.MolToSmiles(mol)
        if canon in seen_canon:
            continue
        seen_canon.add(canon)
        filtered.append((sanitize_name(name), canon))

    # 3. Add curated drugs (override dedup — they're canonical anchors)
    for name, smiles in CURATED_DRUGS:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            log(f"  skip invalid SMILES: {name}")
            continue
        canon = Chem.MolToSmiles(mol)
        if canon in seen_canon:
            # remove from filtered so we use the curated label
            filtered = [(n, s) for (n, s) in filtered if s != canon]
        seen_canon.add(canon)
        filtered.append((sanitize_name(name), canon))

    log(f"  drug-like post-filter + curated: {len(filtered)}")

    # 4. Write PDBQTs
    log("Converting to 3D PDBQT...")
    manifest = []
    n_ok = 0
    for i, (name, smi) in enumerate(filtered):
        pdbqt = to_pdbqt(smi, name, OUT_DIR)
        if pdbqt:
            mol = Chem.MolFromSmiles(smi)
            manifest.append({
                "name": name,
                "smiles": smi,
                "pdbqt_path": str(pdbqt),
                "MW": round(Descriptors.MolWt(mol), 1),
                "logP": round(Descriptors.MolLogP(mol), 2),
                "HBD": Lipinski.NumHDonors(mol),
                "HBA": Lipinski.NumHAcceptors(mol),
                "TPSA": round(Descriptors.TPSA(mol), 1),
                "rot_bonds": Lipinski.NumRotatableBonds(mol),
                "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            })
            n_ok += 1
        if (i + 1) % 50 == 0:
            log(f"  {i+1}/{len(filtered)} processed, {n_ok} OK")
    log(f"Final library: {n_ok} ligands written as PDBQT")

    OUT_JSON.write_text(json.dumps({
        "phase": "3c-v2a",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_scaffolds": len(SCAFFOLDS_SMARTS),
        "n_acid_groups": len(ACID_GROUPS),
        "n_curated": len(CURATED_DRUGS),
        "n_combinatorial_raw": len(combinatorial),
        "n_post_filter": len(filtered),
        "n_written": n_ok,
        "filter": {
            "MW": [MW_MIN, MW_MAX], "logP": [LOGP_MIN, LOGP_MAX],
            "HBD_max": HBD_MAX, "HBA_max": HBA_MAX, "rotB_max": ROTB_MAX,
        },
        "ligands": manifest,
    }, indent=2))
    log(f"Written {OUT_JSON}")


if __name__ == "__main__":
    main()
