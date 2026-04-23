#!/usr/bin/env python3
"""
Phase 3c v3 + Phase 6b — Fenamic-focused expansion library with optional
reversible covalent warheads targeting K1141 Lys ε-NH₂.

Phase 3c v2 identified 2-(arylamino)benzoic acid (fenamic) as the winning
scaffold family (niflumic / flufenamic / mefenamic cluster at top). Phase 3c
v3 enumerates deep around this core: core × N-aryl substituent × ring
substituent × acid bioisostere. Phase 6b adds reversible covalent warheads
(acrylamide, α-cyanoacrylamide, β-ketoamide, α-keto amide, salicylaldehyde-
like aldehyde) attached at ortho/meta/para positions of the N-aryl moiety
pointing toward K1141's lysine.

Library composition target:
  • Non-covalent fenamic expansion    ~1500
  • Covalent warhead variants           ~300
  • Total                                ~1800 ligands

Output:
  `docking_runs/3c_v3/ligands/*.pdbqt`
  `pharmacochaperone_phase3c_v3_fenamic_covalent_library.json`

Deps: rdkit, obabel. Runtime ~3 min.
"""

from __future__ import annotations

import json
import re
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
OUT_DIR = WORK / "docking_runs" / "3c_v3" / "ligands"
OUT_JSON = WORK / "pharmacochaperone_phase3c_v3_fenamic_covalent_library.json"

# Lipinski (loosened for pharmacochaperone; MW up to 500 for covalent warheads)
MW_MIN, MW_MAX = 180.0, 520.0
LOGP_MIN, LOGP_MAX = -1.5, 5.5
HBD_MAX = 5
HBA_MAX = 10
ROTB_MAX = 10

# ---------------------------------------------------------------------------
# Core scaffolds: anthranilic family (2-aminobenzoic acid + analogs)
# [*:1] = N-amine attachment (where arylamine will attach)
# [*:2] = acid group position (usually COOH, can be bioisostere)
# [*:3] = C-ring substituent position (R3 on the benzoic ring, varied by position)
CORES = [
    # Fenamic: anthranilic (2-aminobenzoic)
    ("anthranilic-C2",   "Nc1ccccc1C(=O)O"),
    # 2-amino-nicotinic (niflumic backbone)
    ("2-amino-nicotinic", "Nc1ncccc1C(=O)O"),
    # 2-amino-pyrazine-carboxylic
    ("2-amino-pyrazine",  "Nc1nccnc1C(=O)O"),
    # 2-amino-quinoline-3-carboxylic
    ("2-amino-quinoline-3", "Nc1nc2ccccc2cc1C(=O)O"),
    # 3-amino-benzoic (meta analog, for comparison)
    ("3-amino-benzoic",  "Nc1cccc(C(=O)O)c1"),
    # 3-amino-isonicotinic
    ("3-amino-isonicotinic",  "Nc1cnccc1C(=O)O"),
    # Benzofuran-2-carboxylic 3-amino
    ("3-amino-benzofuran-2-COOH", "Nc1c(C(=O)O)oc2ccccc12"),
    # Indole-3-acetic + NH2 substituent
    ("indole-4-amino-3-acetic", "Nc1cccc2c1c(CC(=O)O)c[nH]2"),
]

# Aryl amine substituents (60 diverse R-NH- groups for N-H replacement)
ARYL_SUBS = [
    # Simple phenyl variants
    ("phenyl",           "c1ccccc1"),
    ("4-MeO-phenyl",     "c1ccc(OC)cc1"),
    ("4-OH-phenyl",      "c1ccc(O)cc1"),
    ("4-F-phenyl",       "c1ccc(F)cc1"),
    ("4-Cl-phenyl",      "c1ccc(Cl)cc1"),
    ("4-Br-phenyl",      "c1ccc(Br)cc1"),
    ("4-CF3-phenyl",     "c1ccc(C(F)(F)F)cc1"),
    ("4-CN-phenyl",      "c1ccc(C#N)cc1"),
    ("4-NO2-phenyl",     "c1ccc([N+](=O)[O-])cc1"),
    ("4-NH2-phenyl",     "c1ccc(N)cc1"),
    ("4-SO2NH2-phenyl",  "c1ccc(S(=O)(=O)N)cc1"),
    ("4-COOH-phenyl",    "c1ccc(C(=O)O)cc1"),
    ("4-COMe-phenyl",    "c1ccc(C(C)=O)cc1"),
    ("4-CONH2-phenyl",   "c1ccc(C(N)=O)cc1"),
    ("3-CF3-phenyl",     "c1cccc(C(F)(F)F)c1"),
    ("3-Cl-phenyl",      "c1cccc(Cl)c1"),
    ("3-F-phenyl",       "c1cccc(F)c1"),
    ("3-OMe-phenyl",     "c1cccc(OC)c1"),
    ("3-OH-phenyl",      "c1cccc(O)c1"),
    ("2-Cl-phenyl",      "c1ccccc1Cl"),
    ("2-F-phenyl",       "c1ccccc1F"),
    ("2-Me-phenyl",      "c1ccccc1C"),
    ("2-CF3-phenyl",     "c1ccccc1C(F)(F)F"),
    ("2,6-diCl-phenyl",  "Clc1ccccc1Cl"),
    ("3,5-diCl-phenyl",  "Clc1cccc(Cl)c1"),
    ("3,5-diMe-phenyl",  "Cc1cccc(C)c1"),
    ("2,6-diMe-phenyl",  "Cc1ccccc1C"),
    ("3,4-diCl-phenyl",  "Clc1ccc(Cl)cc1"),
    ("3,4-OCH2O-phenyl", "c1cc2OCOc2cc1"),
    # Biphenyl
    ("biphenyl",         "c1ccc(-c2ccccc2)cc1"),
    ("4-F-biphenyl",     "c1ccc(-c2ccc(F)cc2)cc1"),
    # Naphthalene
    ("1-naphthyl",       "c1ccc2ccccc2c1"),
    ("2-naphthyl",       "c1ccc2ccccc2c1"),
    # Heteroaryl
    ("2-pyridyl",        "c1ccncc1"),
    ("3-pyridyl",        "c1cccnc1"),
    ("4-pyridyl",        "c1ccncc1"),
    ("2-pyrimidyl",      "c1cncnc1"),
    ("3-thienyl",        "c1ccsc1"),
    ("2-thienyl",        "c1ccsc1"),
    ("2-furyl",          "c1ccoc1"),
    ("3-furyl",          "c1ccoc1"),
    ("2-thiazolyl",      "c1cnsc1"),
    ("2-oxazolyl",       "c1cnoc1"),
    ("2-imidazolyl",     "c1[nH]cnc1"),
    ("2-benzimidazolyl", "c1ccc2[nH]cnc2c1"),
    ("2-benzoxazolyl",   "c1ccc2ocnc2c1"),
    ("2-benzothiazolyl", "c1ccc2scnc2c1"),
    ("8-quinolinyl",     "c1ccc2ncccc2c1"),
    ("5-indolyl",        "c1ccc2[nH]ccc2c1"),
    ("3-indazolyl",      "c1ccc2[nH]ncc2c1"),
    # Cyclic + fused
    ("cyclohexyl",       "C1CCCCC1"),
    ("cyclopropyl",      "C1CC1"),
    ("4-pyranyl",        "O1CCCCC1"),
    ("4-THF",            "C1CCOC1"),
    ("N-Me-piperidine-4-yl", "CN1CCCCC1"),
    ("morpholine-4-yl",  "N1CCOCC1"),
    ("pyrrolidine-1-yl", "N1CCCC1"),
]

# Ring substituents R3 (apply to core free aromatic H, various positions)
RING_SUBS = [
    ("-H",       ""),            # unsubstituted
    ("-Cl",      "Cl"),
    ("-F",       "F"),
    ("-CF3",     "C(F)(F)F"),
    ("-OMe",     "OC"),
    ("-OH",      "O"),
    ("-Me",      "C"),
    ("-NO2",     "[N+](=O)[O-]"),
]

# Acid bioisosteres (replace COOH on core)
ACID_BIOISOSTERES = [
    ("COOH",         "C(=O)O"),
    ("tetrazole",    "c1nnn[nH]1"),
    ("CONHOH",       "C(=O)NO"),
    ("CONHSO2Me",    "C(=O)NS(=O)(=O)C"),
    ("SO2NH2",       "S(=O)(=O)N"),
]

# Covalent warheads (Phase 6b). Reversible lysine-targeting.
# Placed as additional substituent on the N-aryl ring at ortho position
# (to reach K1141 Lys ε-NH₂). Reversibility comes from hemiaminal/β-ketoamide
# chemistry with fast off-rates, not permanent.
COVALENT_WARHEADS = [
    ("acrylamide",        "C(=O)/C=C"),
    ("α-cyanoacrylate",   "C(=O)/C(=C)/C#N"),
    ("β-ketoamide",       "C(=O)CC(C)=O"),
    ("α-ketoamide",       "C(=O)C(=O)C"),
    ("salicylaldehyde-like",  "C=O"),   # benzaldehyde equiv
    ("alpha-chloroacetamide", "C(=O)CCl"),  # reactive but less reversible
]


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_\-]", "_", s)[:80]


def build_fenamic(core_smiles: str, aryl_smiles: str,
                   acid_smiles: str, ring_sub: str,
                   warhead_smiles: str | None = None) -> str | None:
    """Construct full fenamic analog SMILES:

    core: N-[amine]-[ring-scaffold]-[acid] where [amine] is replaced by
    -NH-[aryl]; [acid] is replaced by acid_smiles; one aromatic H replaced
    by ring_sub; optionally one warhead on the aryl.

    Implementation: start from core, substitute the NH₂ with NH-Aryl
    (optional warhead on Aryl), substitute C(=O)O with acid bioisostere,
    substitute one H with ring_sub.
    """
    # Replace NH2 → NH-Aryl (build aryl with or without warhead)
    aryl_mol = Chem.MolFromSmiles(aryl_smiles)
    if aryl_mol is None:
        return None

    if warhead_smiles:
        # Add warhead at a free aromatic position of aryl
        aryl_with_wh = _add_substituent_to_aromatic(aryl_smiles, warhead_smiles)
        if aryl_with_wh is None:
            return None
        aryl_smiles_final = aryl_with_wh
    else:
        aryl_smiles_final = aryl_smiles

    # Build core with N-aryl
    # Find NH2 in core, replace with N(H)-Aryl
    core = Chem.MolFromSmiles(core_smiles)
    if core is None:
        return None

    # Find aromatic amine N
    n_idx = None
    for atom in core.GetAtoms():
        if atom.GetSymbol() == "N" and atom.GetTotalNumHs() == 2:
            n_idx = atom.GetIdx()
            break
    if n_idx is None:
        return None

    # Use SMILES manipulation: combine via RWMol
    try:
        aryl = Chem.MolFromSmiles(aryl_smiles_final)
        if aryl is None:
            return None
        combined = Chem.CombineMols(core, aryl)
        rw = Chem.RWMol(combined)
        aryl_offset = core.GetNumAtoms()
        # Find anchor atom on aryl (first aromatic carbon with H)
        anchor = None
        for a in aryl.GetAtoms():
            if a.GetIsAromatic() and a.GetTotalNumHs() >= 1 and a.GetSymbol() == "C":
                anchor = a.GetIdx()
                break
        if anchor is None:
            return None
        rw.AddBond(n_idx, aryl_offset + anchor, Chem.BondType.SINGLE)
        new_mol = rw.GetMol()
        Chem.SanitizeMol(new_mol)
        new_smiles = Chem.MolToSmiles(new_mol)
    except Exception:
        return None

    # Replace COOH with acid bioisostere if different
    if acid_smiles != "C(=O)O":
        m = Chem.MolFromSmiles(new_smiles)
        if m is None:
            return None
        cooh_pat = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        matches = m.GetSubstructMatches(cooh_pat)
        if matches:
            # Remove COOH, attach acid_smiles at the carbon that held COOH
            match = matches[0]
            cooh_carbon = match[0]
            # Find neighbor of cooh_carbon that is in the ring (the anchor)
            rw = Chem.RWMol(m)
            to_remove = list(match)  # remove all 3 atoms of COOH
            anchor_atom_idx = None
            for nb in rw.GetAtomWithIdx(cooh_carbon).GetNeighbors():
                if nb.GetIdx() not in match:
                    anchor_atom_idx = nb.GetIdx()
                    break
            if anchor_atom_idx is None:
                return None
            # Remove COOH atoms in reverse order to preserve indices
            for idx in sorted(to_remove, reverse=True):
                rw.RemoveAtom(idx)
            # Anchor idx may have shifted
            if anchor_atom_idx > max(to_remove):
                anchor_atom_idx -= len(to_remove)
            # Add acid bioisostere
            acid_mol = Chem.MolFromSmiles(acid_smiles)
            if acid_mol is None:
                return None
            combined = Chem.CombineMols(rw.GetMol(), acid_mol)
            rw2 = Chem.RWMol(combined)
            main_n = rw.GetMol().GetNumAtoms()
            # Find anchor of acid (first atom)
            acid_anchor = main_n + 0
            try:
                rw2.AddBond(anchor_atom_idx, acid_anchor, Chem.BondType.SINGLE)
                final = rw2.GetMol()
                Chem.SanitizeMol(final)
                new_smiles = Chem.MolToSmiles(final)
            except Exception:
                return None
        else:
            # No COOH found (maybe already different); skip bioisostere
            pass

    # Apply ring_sub if specified
    if ring_sub:
        new_smiles = _add_substituent_to_aromatic(new_smiles, ring_sub)
        if new_smiles is None:
            return None

    # Canonicalize
    final = Chem.MolFromSmiles(new_smiles)
    if final is None:
        return None
    return Chem.MolToSmiles(final)


def _add_substituent_to_aromatic(smiles: str, sub_smiles: str) -> str | None:
    """Attach sub_smiles at the first aromatic C with a free H."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for a in mol.GetAtoms():
        if a.GetIsAromatic() and a.GetTotalNumHs() >= 1 and a.GetSymbol() == "C":
            anchor = a.GetIdx()
            break
    else:
        return None
    sub_mol = Chem.MolFromSmiles(sub_smiles)
    if sub_mol is None:
        return None
    combined = Chem.CombineMols(mol, sub_mol)
    rw = Chem.RWMol(combined)
    sub_start = mol.GetNumAtoms()
    try:
        rw.AddBond(anchor, sub_start, Chem.BondType.SINGLE)
        final = rw.GetMol()
        Chem.SanitizeMol(final)
        return Chem.MolToSmiles(final)
    except Exception:
        return None


def drug_like(mol: Chem.Mol) -> bool:
    if mol is None:
        return False
    mw = Descriptors.MolWt(mol)
    if not (MW_MIN <= mw <= MW_MAX):
        return False
    logp = Descriptors.MolLogP(mol)
    if not (LOGP_MIN <= logp <= LOGP_MAX):
        return False
    if Lipinski.NumHDonors(mol) > HBD_MAX:
        return False
    if Lipinski.NumHAcceptors(mol) > HBA_MAX:
        return False
    if Lipinski.NumRotatableBonds(mol) > ROTB_MAX:
        return False
    if rdMolDescriptors.CalcNumAromaticRings(mol) < 1:
        return False
    acid_smarts = [
        "C(=O)[OH]",
        "S(=O)(=O)[NH]",
        "C(=O)N[OH]",
        "P(=O)([OH])[OH]",
        "c1nnn[nH]1",
    ]
    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in acid_smarts):
        return False
    return True


def to_pdbqt(smiles: str, name: str, out_dir: Path) -> Path | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    try:
        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
            return None
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        return None
    with tempfile.TemporaryDirectory() as tmpd:
        sdf = Path(tmpd) / f"{name}.sdf"
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


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log("=== Phase 3c v3 + 6b fenamic-focused + covalent library ===")

    seen: set[str] = set()
    entries: list[tuple[str, str, bool]] = []  # (name, smiles, is_covalent)

    # --- Non-covalent expansion (Phase 3c v3) ---
    log("Phase 3c v3 non-covalent enumeration...")
    n_nc = 0
    for core_name, core_smiles in CORES:
        for aryl_name, aryl_smiles in ARYL_SUBS:
            for acid_name, acid_smiles in ACID_BIOISOSTERES:
                for rs_name, ring_sub in RING_SUBS:
                    smi = build_fenamic(core_smiles, aryl_smiles,
                                         acid_smiles, ring_sub, warhead_smiles=None)
                    if smi is None or smi in seen:
                        continue
                    mol = Chem.MolFromSmiles(smi)
                    if not drug_like(mol):
                        continue
                    seen.add(smi)
                    name = safe_name(
                        f"nc__{core_name}__{aryl_name}__{acid_name}__{rs_name}"
                    )
                    entries.append((name, smi, False))
                    n_nc += 1
    log(f"  non-covalent unique: {n_nc}")

    # --- Covalent warhead analogs (Phase 6b) ---
    # Add warheads only on selected high-performing N-aryl substituents
    # (3-CF3-phenyl, biphenyl, 2,6-diCl-phenyl, 3,5-diCl-phenyl, 3-Cl-phenyl,
    #  4-Cl-phenyl, 3,4-diCl-phenyl, 3-F-phenyl, naphthyl — 2-aryls matching
    #  top Phase 3c v2 hits)
    log("Phase 6b covalent warhead enumeration...")
    HOT_ARYLS = [
        ("3-CF3-phenyl", "c1cccc(C(F)(F)F)c1"),
        ("4-Cl-phenyl", "c1ccc(Cl)cc1"),
        ("3-Cl-phenyl", "c1cccc(Cl)c1"),
        ("2,6-diCl-phenyl", "Clc1ccccc1Cl"),
        ("3,5-diCl-phenyl", "Clc1cccc(Cl)c1"),
        ("3,4-diCl-phenyl", "Clc1ccc(Cl)cc1"),
        ("biphenyl", "c1ccc(-c2ccccc2)cc1"),
        ("2-naphthyl", "c1ccc2ccccc2c1"),
        ("4-F-phenyl", "c1ccc(F)cc1"),
        ("3-F-phenyl", "c1cccc(F)c1"),
    ]
    n_cov = 0
    for core_name, core_smiles in CORES:
        for aryl_name, aryl_smiles in HOT_ARYLS:
            for wh_name, wh_smiles in COVALENT_WARHEADS:
                for acid_name, acid_smiles in [("COOH", "C(=O)O"),
                                                 ("tetrazole", "c1nnn[nH]1")]:
                    smi = build_fenamic(core_smiles, aryl_smiles,
                                         acid_smiles, "",
                                         warhead_smiles=wh_smiles)
                    if smi is None or smi in seen:
                        continue
                    mol = Chem.MolFromSmiles(smi)
                    if not drug_like(mol):
                        continue
                    seen.add(smi)
                    name = safe_name(
                        f"cov__{core_name}__{aryl_name}__{wh_name}__{acid_name}"
                    )
                    entries.append((name, smi, True))
                    n_cov += 1
    log(f"  covalent unique: {n_cov}")

    log(f"Total library: {len(entries)} ligands, converting to PDBQT...")

    manifest = []
    n_ok = 0
    for i, (name, smi, is_cov) in enumerate(entries):
        pdbqt = to_pdbqt(smi, name, OUT_DIR)
        if pdbqt is None:
            continue
        mol = Chem.MolFromSmiles(smi)
        manifest.append({
            "name": name,
            "smiles": smi,
            "pdbqt_path": str(pdbqt),
            "is_covalent": is_cov,
            "MW": round(Descriptors.MolWt(mol), 1),
            "logP": round(Descriptors.MolLogP(mol), 2),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "rotB": Lipinski.NumRotatableBonds(mol),
        })
        n_ok += 1
        if (i + 1) % 200 == 0:
            log(f"  {i+1}/{len(entries)} processed, {n_ok} OK")
    log(f"Final: {n_ok} ligands written.")

    OUT_JSON.write_text(json.dumps({
        "phase": "3c-v3 + 6b",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_cores": len(CORES),
        "n_aryl_subs": len(ARYL_SUBS),
        "n_hot_aryls_for_covalent": len(HOT_ARYLS),
        "n_ring_subs": len(RING_SUBS),
        "n_acid_isosteres": len(ACID_BIOISOSTERES),
        "n_warheads": len(COVALENT_WARHEADS),
        "n_non_covalent": n_nc,
        "n_covalent": n_cov,
        "n_unique": len(entries),
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
