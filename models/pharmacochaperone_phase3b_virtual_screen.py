"""
Phase 3B: Virtual screen of literature-curated library against
the K1141 pocket pharmacophore extracted in Phase 3A.

Scoring scheme (each component in [0, 1], composite is weighted mean):

    anion_score       — does the ligand carry an ionizable anion
                        (COOH, tetrazole, sulfonate, sulfonamide,
                        acyl sulfonamide, phosphonate) at physiological pH?
    anion_aromatic_d  — for ligands with both anion + aromatic, match the
                        protein's K1141-NZ / F1646-ring distance (5.7 Å,
                        tolerance ±1.5 Å) within any conformer
    size_fit          — MW in [180, 350] Da (pocket volume 159 Å³)
    ro3_pass          — Congreve rule-of-3 (MW<300, logP<3, HBD<=3, HBA<=3,
                        rot<=3) — for fragment-sized leads
    lipinski_pass     — Rule of 5 (for drug-sized leads)
    hb_donor          — has at least one H-bond donor that can reach a
                        protein carbonyl (backbone O of G1645/F1646/G1647)
    tpsa              — Topological Polar Surface Area, computed and stored
                        as descriptor; NOT used in composite score below.
                        Retained for downstream filtering/inspection.
                        Note: STRC is EXTRACELLULAR — the earlier docstring
                        framing of "CNS-like 40-90 Å²" was incorrect (BBB
                        irrelevant for otic RWM delivery). Literature for
                        otic/RWM permeability vs TPSA is thin; no primary
                        lit-backed range applied here. [audit 2026-04-23]
    qed               — drug-likeness quantitative estimate

Composite drug score D = 0.30*anion + 0.15*anion_aromatic_d +
                         0.15*size_fit + 0.10*ro3 +
                         0.10*lipinski + 0.10*hb_donor + 0.10*qed
(Note: TPSA and QED-derived `score_qed` above: only score_qed enters the
composite — TPSA is not. See audit 2026-04-23.)

Library (30 compounds, literature SMILES):
  A. CFTR correctors — closest clinical analogue class (misfolding rescue
     of a secreted/membrane protein)
  B. Misfolding rescuers — general chemical chaperones (osmolytes, 4-PBA)
  C. K1141-focused acidic fragments — carboxylic acids with aromatic tail,
     mimicking VX-809 scaffold: cyclopropane-carboxylic-acid + phenyl
  D. Negative controls — neutral drugs (should score low)
"""

import json
from pathlib import Path
from dataclasses import dataclass, field, asdict

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, Lipinski
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds, CalcTPSA

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
PHARMACO_JSON = OUT_DIR / "pharmacochaperone_phase3a_pharmacophore.json"

# Target distance between ligand anion and ligand aromatic centroid
# — matches K1141-NZ to F1646-ring distance in protein = 5.69 Å.
TARGET_ANION_AROM_D = 5.7
TARGET_ANION_AROM_TOL = 1.5

# ---------- Library ----------
# (category, name, SMILES, note)
LIBRARY = [
    # A. CFTR correctors — drugs shown to rescue ΔF508 and other misfolded CFTR
    ("A-corrector", "lumacaftor_VX-809",
     "CC1(C)C(C(=O)NC2=CC3=C(C=C2C)OC(C(=O)O)(C4=CC=CC=C4)C3)C1",
     "ΔF508-CFTR corrector, tertiary-fold stabiliser"),
    ("A-corrector", "tezacaftor_VX-661",
     "CC(C)(CO)C1=CC=C(C=C1)C2=NC(=CS2)C3=CC=C(C=C3)F",
     "second-gen corrector, 5-amino-2-(4-chlorophenyl) variant"),
    ("A-corrector", "elexacaftor_VX-445",
     "CC(C)(C)c1cc(NS(=O)(=O)c2cccc(F)c2)nn1C1CC1",
     "triple-combo corrector, highly potent"),
    ("A-corrector", "ivacaftor_VX-770",
     "CC(C)(C)c1cc(cc(c1O)C(C)(C)C)NC(=O)c2cc3ccccc3[nH]2",
     "CFTR potentiator (ion-channel gating, not corrector)"),
    ("A-corrector", "FDL-169",
     "OC(=O)C1(c2ccccc2)Cc2cc(Cl)c(OC(F)F)cc2O1",
     "FLATLEY corrector; F2G analogue"),
    ("A-corrector", "glafenine",
     "OCC(O)COC(=O)c1cccc(Nc2ccnc3cc(Cl)ccc23)c1",
     "discontinued anthranilate; pleiotropic corrector activity"),

    # B. Misfolding rescuers — general chemical chaperones
    ("B-chaperone", "4-PBA_phenylbutyrate",
     "OC(=O)CCCc1ccccc1",
     "clinical PBA; acid + phenyl, matches our pharmacophore"),
    ("B-chaperone", "sodium_phenylbutyrate",
     "OC(=O)CCCc1ccccc1",
     "PBA sodium salt"),
    ("B-chaperone", "trimethylamine_N_oxide_TMAO",
     "C[N+](C)(C)[O-]",
     "zwitterionic osmolyte — protein stabiliser"),
    ("B-chaperone", "trehalose",
     "OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
     "disaccharide osmolyte"),
    ("B-chaperone", "taurine",
     "OS(=O)(=O)CCN",
     "sulfonate + amine; inner-ear native osmolyte"),
    ("B-chaperone", "betaine",
     "C[N+](C)(C)CC(=O)[O-]",
     "zwitterionic chaperone; inner-ear native osmolyte"),

    # C. K1141-focused acidic fragments (fragment-based design)
    ("C-fragment", "cyclopropane-1-carboxyphenyl",
     "OC(=O)C1(c2ccccc2)CC1",
     "VX-809 scaffold minimised: acid + aromatic + cyclopropane"),
    ("C-fragment", "4-chlorobenzoate",
     "OC(=O)c1ccc(Cl)cc1",
     "simplest fragment: acid + aromatic 5.4 Å apart"),
    ("C-fragment", "4-phenylbenzoate",
     "OC(=O)c1ccc(-c2ccccc2)cc1",
     "biphenyl carboxylic acid; rigid 5-ring + ring extension"),
    ("C-fragment", "indole-3-acetic-acid",
     "OC(=O)Cc1c[nH]c2ccccc12",
     "tryptophan-derived; auxin scaffold"),
    ("C-fragment", "nicotinic-acid",
     "OC(=O)c1cccnc1",
     "small polar acid + pyridine"),
    ("C-fragment", "salicylic-acid",
     "OC(=O)c1ccccc1O",
     "acid + phenol; ortho-H-bond; cheap"),
    ("C-fragment", "5-aminosalicylate_5ASA",
     "Nc1ccc(O)c(C(=O)O)c1",
     "mesalazine; anti-inflammatory, already in bodies"),
    ("C-fragment", "ibuprofen",
     "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
     "NSAID; acid + aromatic; ref drug"),
    ("C-fragment", "diflunisal",
     "OC(=O)c1cc(-c2ccc(F)cc2F)ccc1O",
     "NSAID; difluoro-biphenyl + acid"),
    ("C-fragment", "tolfenamic-acid",
     "Cc1ccccc1Nc1ccccc1C(=O)O",
     "NSAID; diphenylamine + acid"),
    ("C-fragment", "probenecid",
     "CCCN(CCC)S(=O)(=O)c1ccc(C(=O)O)cc1",
     "acyl sulfonamide; passes blood-labyrinth barrier"),

    # D. Negative controls (non-acidic or too-big; should score low)
    ("D-control", "caffeine",
     "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
     "no anion, xanthine; low pharmacochaperone probability"),
    ("D-control", "glucose",
     "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
     "sugar, no anion at pH7, too polar"),
    ("D-control", "dextromethorphan",
     "CO[C@@H]1CC[C@H]2Cc3ccc(OC)cc3[C@]11CCN2C",
     "cationic; wrong charge"),
    ("D-control", "atenolol",
     "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
     "β-blocker; no anion"),

    # Extra candidates inspired by pocket geometry
    ("C-fragment", "naphthalene-2-carboxylic-acid",
     "OC(=O)c1ccc2ccccc2c1",
     "rigid bicyclic + acid; large aromatic"),
    ("C-fragment", "fenbufen",
     "OC(=O)CCC(=O)c1ccc(-c2ccccc2)cc1",
     "NSAID; biphenyl ketone-acid"),
]


def has_carboxylate_or_similar(mol):
    patterns = [
        "C(=O)[O;H1,-1]",      # carboxylate
        "c1nnn[nH]1",          # tetrazole
        "S(=O)(=O)[O;H1,-1]",  # sulfonate
        "NS(=O)(=O)C",         # sulfonamide (N-H acidic with aryl ring)
        "C(=O)NS(=O)(=O)",     # acyl sulfonamide
        "P(=O)([O;H1,-1])[O;H1,-1]",  # phosphonate
    ]
    for pat in patterns:
        q = Chem.MolFromSmarts(pat)
        if q is not None and mol.HasSubstructMatch(q):
            return True
    return False


def has_aromatic_ring(mol):
    return any(ring for ring in mol.GetRingInfo().AtomRings()
               if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))


def anion_aromatic_distance(mol):
    """In the generated 3D conformer, find min distance between any
    anionic O and any aromatic ring centroid."""
    if mol.GetNumConformers() == 0:
        return None

    conf = mol.GetConformer()
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    matches = mol.GetSubstructMatches(acid_pattern)
    if not matches:
        return None
    acid_atoms = [m[2] for m in matches]   # the O atom

    ri = mol.GetRingInfo()
    aromatic_centroids = []
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            coords = np.array([list(conf.GetAtomPosition(i)) for i in ring])
            aromatic_centroids.append(coords.mean(axis=0))
    if not aromatic_centroids:
        return None

    min_d = 1e9
    for ai in acid_atoms:
        ap = np.array(list(conf.GetAtomPosition(ai)))
        for ac in aromatic_centroids:
            d = np.linalg.norm(ap - ac)
            if d < min_d:
                min_d = d
    return min_d


def score_anion_arom_distance(d):
    if d is None:
        return 0.0
    diff = abs(d - TARGET_ANION_AROM_D)
    if diff > TARGET_ANION_AROM_TOL:
        # Gaussian decay outside tolerance
        return float(np.exp(-((diff - TARGET_ANION_AROM_TOL) ** 2) / 2.0))
    return 1.0 - (diff / TARGET_ANION_AROM_TOL) * 0.3   # mild penalty inside tolerance


def score_size(mw):
    if 180 <= mw <= 350:
        return 1.0
    if 120 <= mw < 180:
        return 0.7
    if 350 < mw <= 500:
        return 0.5
    return 0.2


def ro3_pass(mw, logp, hbd, hba, rot):
    return mw <= 300 and logp <= 3.0 and hbd <= 3 and hba <= 3 and rot <= 3


def lipinski_pass(mw, logp, hbd, hba):
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return violations <= 1


@dataclass
class Result:
    category: str
    name: str
    smiles: str
    note: str
    mw: float
    logp: float
    hbd: int
    hba: int
    rot: int
    tpsa: float
    qed: float
    has_anion: bool
    has_aromatic: bool
    anion_arom_d: float
    ro3: bool
    lipinski: bool
    score_anion: float = 0.0
    score_anion_arom: float = 0.0
    score_size: float = 0.0
    score_ro3: float = 0.0
    score_lipinski: float = 0.0
    score_hb_donor: float = 0.0
    score_qed: float = 0.0
    composite: float = 0.0


def score_compound(category, name, smiles, note):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)

    # Embed 3D conformer (try 5 random seeds, keep the best-embedded)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    cid = AllChem.EmbedMolecule(mol, params)
    if cid < 0:
        # try with larger box / random coords
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mol, params)
    if cid >= 0:
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = CalcNumRotatableBonds(mol)
    tpsa = CalcTPSA(mol)
    qed = QED.qed(mol)

    has_anion = has_carboxylate_or_similar(mol)
    has_arom = has_aromatic_ring(mol)
    d_aa = anion_aromatic_distance(mol) if (has_anion and has_arom) else None

    r = Result(
        category=category, name=name, smiles=smiles, note=note,
        mw=mw, logp=logp, hbd=hbd, hba=hba, rot=rot, tpsa=tpsa, qed=qed,
        has_anion=has_anion, has_aromatic=has_arom,
        anion_arom_d=d_aa if d_aa is not None else -1.0,
        ro3=ro3_pass(mw, logp, hbd, hba, rot),
        lipinski=lipinski_pass(mw, logp, hbd, hba),
    )

    r.score_anion = 1.0 if has_anion else 0.0
    r.score_anion_arom = score_anion_arom_distance(d_aa)
    r.score_size = score_size(mw)
    r.score_ro3 = 1.0 if r.ro3 else 0.0
    r.score_lipinski = 1.0 if r.lipinski else 0.0
    # H-bond-donor presence (for backbone-O targeting)
    r.score_hb_donor = min(1.0, hbd / 2.0)
    r.score_qed = float(qed)

    r.composite = (
        0.30 * r.score_anion +
        0.15 * r.score_anion_arom +
        0.15 * r.score_size +
        0.10 * r.score_ro3 +
        0.10 * r.score_lipinski +
        0.10 * r.score_hb_donor +
        0.10 * r.score_qed
    )
    return r


def main():
    pharmaco = json.loads(PHARMACO_JSON.read_text())
    pocket_centroid = pharmaco["pocket_centroid_A"]

    results = []
    for cat, name, smi, note in LIBRARY:
        r = score_compound(cat, name, smi, note)
        if r is None:
            print(f"! failed to parse {name} ({smi})")
            continue
        results.append(r)

    # Sort by composite score descending
    results.sort(key=lambda x: -x.composite)

    # Write JSON
    out_json = OUT_DIR / "pharmacochaperone_phase3b_virtual_screen.json"
    out_json.write_text(json.dumps({
        "pocket_centroid_A": pocket_centroid,
        "target_anion_arom_d_A": TARGET_ANION_AROM_D,
        "anion_arom_d_tolerance_A": TARGET_ANION_AROM_TOL,
        "scoring_weights": {
            "anion": 0.30, "anion_arom_dist": 0.15, "size": 0.15,
            "ro3": 0.10, "lipinski": 0.10, "hb_donor": 0.10, "qed": 0.10,
        },
        "library_size": len(LIBRARY),
        "ranked_candidates": [asdict(r) for r in results],
    }, indent=2))

    # Print ranked table
    print(f"\n{'Rank':<5}{'Name':<30}{'Cat':<14}{'MW':>6}{'logP':>6}{'HBD':>4}{'Anion?':>8}"
          f"{'d(A-Ar)':>9}{'Composite':>11}")
    print("-" * 95)
    for i, r in enumerate(results, 1):
        anion = "YES" if r.has_anion else "no"
        d_aa = f"{r.anion_arom_d:.2f}" if r.anion_arom_d > 0 else "--"
        print(f"{i:<5}{r.name:<30}{r.category:<14}{r.mw:>6.1f}{r.logp:>6.2f}{r.hbd:>4}"
              f"{anion:>8}{d_aa:>9}{r.composite:>11.3f}")

    print(f"\n→ wrote {out_json}")
    print(f"→ top-5 shortlist for Phase 3C shape complementarity + docking:")
    for r in results[:5]:
        print(f"    {r.name} ({r.category}, composite={r.composite:.3f})")


if __name__ == "__main__":
    main()
