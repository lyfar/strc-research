"""
Phase 3A: Extract 3D pharmacophore from the Phase 2B K1141 pocket.

Goal: convert the geometric pocket into a *negative image* — a
pharmacophore specification that a small-molecule can be scored
against for fit.

Pocket-defining anchor residues (from Phase 2B subpocket 1 + 2):
    D1140, K1141, R1169, D1173, G1645, F1646, G1647

For each pharmacophoric feature we record:
    kind  — one of {anion_target, cation_target, hb_donor_target,
                     hb_acceptor_target, aromatic_target, hydrophobic_target,
                     disease_residue}
    xyz   — 3D position in Job-4 WT frame
    residue, atom — source of the feature

The resulting JSON is the ligand design target. A ligand is a good
pharmacochaperone candidate if its 3D feature pattern COMPLEMENTS
this pharmacophore:

    protein cation (K1141 NZ) ⟶ ligand needs ANION within 3-4 Å
    protein aromatic (F1646)  ⟶ ligand needs AROMATIC within 4-6 Å (π-stack)
    protein hydrophobic pocket⟶ ligand needs HYDROPHOBIC fill
    protein backbone-NH       ⟶ ligand needs H-BOND ACCEPTOR
    protein backbone-O=C      ⟶ ligand needs H-BOND DONOR
    mutation site E1659→A     ⟶ ligand should NOT clash here (MUT has
                                 a cavity gain; ligand should compensate
                                 from the adjacent pocket)
"""

import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

MODELS_DIR = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models")
WT_CIF = MODELS_DIR / "job4-wildtype.cif"
MUT_CIF = MODELS_DIR / "job3-mutant.cif"
OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")

POCKET_RESIDUES = [
    1140, 1141, 1142, 1143, 1144,
    1168, 1169, 1170, 1172, 1173,
    1612,
    1642, 1643, 1644, 1645, 1646, 1647, 1648,
    1651, 1652, 1653, 1655, 1656,
    1659,
]

# Residue categorization
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS", "GLY"}
POLAR       = {"SER", "THR", "ASN", "GLN", "TYR", "HIS"}
CATIONIC    = {"LYS", "ARG"}
ANIONIC     = {"ASP", "GLU"}
AROMATIC    = {"PHE", "TYR", "TRP", "HIS"}

DONORS_SC = {
    "LYS": ["NZ"],
    "ARG": ["NE", "NH1", "NH2"],
    "TRP": ["NE1"],
    "ASN": ["ND2"],
    "GLN": ["NE2"],
    "HIS": ["ND1", "NE2"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
}
ACCEPTORS_SC = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "ASN": ["OD1"],
    "GLN": ["OE1"],
    "HIS": ["ND1", "NE2"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
}


def load_chain_A(cif_path: Path):
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", cif_path)
    model = next(struct.get_models())
    return model["A"]


def aromatic_ring_centroid(residue):
    if residue.resname == "PHE":
        atom_names = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    elif residue.resname == "TYR":
        atom_names = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    elif residue.resname == "TRP":
        atom_names = ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
    elif residue.resname == "HIS":
        atom_names = ["CG", "ND1", "CE1", "NE2", "CD2"]
    else:
        return None
    coords = []
    for name in atom_names:
        if name in residue:
            coords.append(residue[name].coord)
    if len(coords) < 4:
        return None
    return np.mean(coords, axis=0)


def extract_features(chain):
    """Return a list of pharmacophore feature dicts."""
    features = []

    for res in chain:
        if not is_aa(res, standard=True):
            continue
        resnum = res.id[1]
        if resnum not in POCKET_RESIDUES:
            continue
        rn = res.resname

        # Backbone N-H donor (except proline)
        if rn != "PRO" and "N" in res:
            features.append({
                "kind": "hb_acceptor_target",   # protein N-H donates → needs ligand acceptor
                "xyz": [float(x) for x in res["N"].coord],
                "residue": resnum,
                "resname": rn,
                "atom": "N",
                "source": "backbone",
            })

        # Backbone C=O acceptor
        if "O" in res:
            features.append({
                "kind": "hb_donor_target",      # protein C=O accepts → needs ligand donor
                "xyz": [float(x) for x in res["O"].coord],
                "residue": resnum,
                "resname": rn,
                "atom": "O",
                "source": "backbone",
            })

        # Cationic side chain → ligand needs anion
        if rn in CATIONIC:
            for name in DONORS_SC.get(rn, []):
                if name in res:
                    features.append({
                        "kind": "anion_target",
                        "xyz": [float(x) for x in res[name].coord],
                        "residue": resnum,
                        "resname": rn,
                        "atom": name,
                        "source": "sidechain",
                    })

        # Anionic side chain → ligand needs cation/donor
        if rn in ANIONIC:
            for name in ACCEPTORS_SC.get(rn, []):
                if name in res:
                    features.append({
                        "kind": "cation_target",
                        "xyz": [float(x) for x in res[name].coord],
                        "residue": resnum,
                        "resname": rn,
                        "atom": name,
                        "source": "sidechain",
                    })

        # Aromatic ring → ligand needs aromatic
        if rn in AROMATIC:
            ring = aromatic_ring_centroid(res)
            if ring is not None:
                features.append({
                    "kind": "aromatic_target",
                    "xyz": [float(x) for x in ring],
                    "residue": resnum,
                    "resname": rn,
                    "atom": "RING_CENTROID",
                    "source": "sidechain",
                })

        # Hydrophobic side chain centroid → ligand needs hydrophobic
        if rn in HYDROPHOBIC and rn not in ("GLY", "ALA"):
            side_atoms = [a for a in res if a.name not in ("N", "CA", "C", "O", "H")]
            if len(side_atoms) >= 2:
                centroid = np.mean([a.coord for a in side_atoms], axis=0)
                features.append({
                    "kind": "hydrophobic_target",
                    "xyz": [float(x) for x in centroid],
                    "residue": resnum,
                    "resname": rn,
                    "atom": "SC_CENTROID",
                    "source": "sidechain",
                })

    return features


def compute_disease_delta(wt_chain, mut_chain, resnum=1659):
    """E1659A mutation → cavity gain. Record ligand keep-out (or fill)."""
    if resnum in [r.id[1] for r in wt_chain if is_aa(r, standard=True)]:
        wt_res = [r for r in wt_chain if r.id[1] == resnum][0]
    else:
        return None
    mut_res = [r for r in mut_chain if r.id[1] == resnum][0]

    wt_sc = [a for a in wt_res if a.name in ("CG", "CD", "OE1", "OE2", "OD1", "OD2")]
    mut_sc = [a for a in mut_res if a.name in ("CB",)]

    if not wt_sc or not mut_sc:
        return None

    wt_sc_centroid = np.mean([a.coord for a in wt_sc], axis=0)
    mut_sc_centroid = mut_res["CB"].coord
    gap = np.linalg.norm(wt_sc_centroid - mut_sc_centroid)

    return {
        "wt_side_chain_centroid": [float(x) for x in wt_sc_centroid],
        "mut_cb": [float(x) for x in mut_sc_centroid],
        "distance_A": float(gap),
        "note": "Atoms lost in MUT: CG, CD, OE1, OE2 (4 heavy atoms; ~45 Å³ volume)",
    }


def compute_centroids_and_distances(features):
    """Distance matrix of feature centroids — used for pharmacophore triangles."""
    # Define the CORE anchor features used for docking-orientation constraints.
    # In rank order of importance:
    #   1. K1141 NZ (anion target) — THE salt-bridge anchor
    #   2. F1646 aromatic ring — π-stack
    #   3. G1645/G1647 backbone carbonyl — ligand H-bond donor target
    #   4. R1169 NE/NH2 — secondary anion target
    #   5. D1140/D1173 OD — keep-out for ligand acids (already-occupied cation)
    k1141_nz = next((f for f in features if f["residue"] == 1141 and f["atom"] == "NZ"), None)
    f1646_ring = next((f for f in features if f["residue"] == 1646 and f["atom"] == "RING_CENTROID"), None)
    g1645_o = next((f for f in features if f["residue"] == 1645 and f["atom"] == "O"), None)
    g1647_o = next((f for f in features if f["residue"] == 1647 and f["atom"] == "O"), None)

    def pair_dist(a, b):
        if a is None or b is None:
            return None
        return float(np.linalg.norm(np.array(a["xyz"]) - np.array(b["xyz"])))

    return {
        "k1141_nz_to_f1646_ring_A": pair_dist(k1141_nz, f1646_ring),
        "k1141_nz_to_g1645_o_A": pair_dist(k1141_nz, g1645_o),
        "k1141_nz_to_g1647_o_A": pair_dist(k1141_nz, g1647_o),
        "f1646_ring_to_g1645_o_A": pair_dist(f1646_ring, g1645_o),
        "anchor_triangle_vertices": {
            "K1141_NZ": k1141_nz["xyz"] if k1141_nz else None,
            "F1646_RING": f1646_ring["xyz"] if f1646_ring else None,
            "G1645_O": g1645_o["xyz"] if g1645_o else None,
        },
    }


def pocket_centroid_from_phase2b():
    """Pocket centroid from Phase 2B (combined subpocket 1+2)."""
    subpocket_1_centroid = np.array([9.52, -6.14, -39.83])
    subpocket_2_centroid = np.array([5.89, -4.68, -43.13])
    combined = np.mean([subpocket_1_centroid, subpocket_2_centroid], axis=0)
    return [float(x) for x in combined]


def summarize(features):
    by_kind = {}
    for f in features:
        by_kind.setdefault(f["kind"], []).append(f)
    return {k: len(v) for k, v in by_kind.items()}


def main():
    wt_chain = load_chain_A(WT_CIF)
    mut_chain = load_chain_A(MUT_CIF)

    features = extract_features(wt_chain)
    anchors = compute_centroids_and_distances(features)
    pocket_centroid = pocket_centroid_from_phase2b()

    # Distance of each feature to the pocket centroid (relevance weight)
    pc = np.array(pocket_centroid)
    for f in features:
        f["dist_to_pocket_A"] = float(np.linalg.norm(np.array(f["xyz"]) - pc))

    # Keep only features within 12 Å of pocket centroid — these reach the ligand
    near = [f for f in features if f["dist_to_pocket_A"] <= 12.0]

    disease_delta = compute_disease_delta(wt_chain, mut_chain, resnum=1659)
    summary = summarize(near)

    result = {
        "source": "Phase 3A pharmacophore extraction",
        "wt_cif": str(WT_CIF.name),
        "pocket_residues": POCKET_RESIDUES,
        "pocket_centroid_A": pocket_centroid,
        "feature_cutoff_to_pocket_A": 12.0,
        "n_features_total": len(features),
        "n_features_near_pocket": len(near),
        "feature_counts_by_kind": summary,
        "anchor_triangle": anchors,
        "disease_residue_delta_E1659A": disease_delta,
        "ligand_design_rules": {
            "required_anionic_group": True,
            "required_anion_xyz_target": [f["xyz"] for f in near
                                          if f["residue"] == 1141 and f["atom"] == "NZ"],
            "preferred_aromatic_xyz_target": [f["xyz"] for f in near
                                              if f["residue"] == 1646
                                              and f["atom"] == "RING_CENTROID"],
            "mw_range_Da": [150, 350],
            "logp_range": [-1.0, 4.0],
            "rule_of_3_preferred": True,
            "h_bond_donor_ligand_for_backbone_acceptor": True,
            "keep_out_of_e1659_mut_cavity": False,
            "notes": (
                "E1659A mutation loses the side-chain OE1/OE2 that in WT "
                "formed an intra-chain salt bridge with K1141. Pharmacochaperone "
                "ligand must donate a carboxylate within 3.5 Å of K1141 NZ to "
                "RESTORE the lost salt bridge (compensatory pharmacochaperone)."
            ),
        },
        "near_features": near,
    }

    out_json = OUT_DIR / "pharmacochaperone_phase3a_pharmacophore.json"
    out_json.write_text(json.dumps(result, indent=2))

    # Summary output
    print(f"Pocket centroid: {pocket_centroid}")
    print(f"Features near pocket (<=12 Å): {len(near)}")
    print(f"  by kind: {summary}")
    print(f"K1141-NZ to F1646 ring: {anchors['k1141_nz_to_f1646_ring_A']:.2f} Å")
    print(f"K1141-NZ to G1645 O:    {anchors['k1141_nz_to_g1645_o_A']:.2f} Å")
    print(f"K1141-NZ to G1647 O:    {anchors['k1141_nz_to_g1647_o_A']:.2f} Å")
    print(f"E1659 WT sidechain to MUT-CB: {disease_delta['distance_A']:.2f} Å")
    print(f"→ wrote {out_json}")


if __name__ == "__main__":
    main()
