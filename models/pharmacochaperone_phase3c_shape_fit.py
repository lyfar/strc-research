"""
Phase 3C: Shape complementarity scoring for the top-ranked Phase 3B
candidates against the Phase 2B K1141 pocket.

Algorithm (light-weight, no Vina required):
  1. Load Phase 2B subpocket centroid (7.7, -5.4, -41.5 Å) and K1141 NZ
     position from Phase 3A pharmacophore JSON.
  2. Load the WT protein heavy atoms within 15 Å of pocket centroid
     (from job4-wildtype.cif).
  3. For each Phase 3B top-10 candidate:
     a. Generate N conformers (ETKDG + MMFF optimisation).
     b. For each conformer:
          - Find the acid carbon (carboxylate-C) and aromatic centroid.
          - Translate the ligand so acid-C → (K1141-NZ - 2.8 Å along
            the NZ-to-pocket-centroid axis). This places the anion
            in salt-bridge distance from K1141 NZ.
          - For 24 rotations about the anion-protein-NZ axis,
            compute:
              pocket_fill = heavy atoms within pocket Bsphere (5.5 Å
                            of centroid)
              clash       = heavy atoms within 2.0 Å of any protein
                            heavy atom
              contact     = heavy atoms within 4.0 Å of any protein
                            heavy atom (productive contact)
            shape_score = (pocket_fill + contact) / (clash + 1) / total_heavy
     c. Keep the best (pose, score).

Output: ranked JSON with top pose per candidate + predicted binding mode
(acid-C coords, aromatic-centroid coords, rotation angle).
"""

import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from rdkit import Chem
from rdkit.Chem import AllChem

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
WT_CIF = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job4-wildtype.cif")
SCREEN_JSON = OUT_DIR / "pharmacochaperone_phase3b_virtual_screen.json"
PHARMACO_JSON = OUT_DIR / "pharmacochaperone_phase3a_pharmacophore.json"

POCKET_SPHERE_RADIUS = 5.5       # pocket centered ball
PROTEIN_CONTEXT_RADIUS = 15.0    # protein atoms within this from centroid
SALT_BRIDGE_DIST = 2.8           # K1141-NZ ... ligand-O distance (Å)
CLASH_CUTOFF = 2.0
CONTACT_CUTOFF = 4.0
N_CONFORMERS = 30
N_ROTATIONS = 24


def load_protein_atoms_near(cif_path, centroid_xyz, radius=PROTEIN_CONTEXT_RADIUS):
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", cif_path)
    model = next(struct.get_models())
    chain = model["A"]

    centroid = np.array(centroid_xyz)
    coords = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        for atom in res:
            if atom.element in ("H", "D"):
                continue
            c = atom.coord
            if np.linalg.norm(c - centroid) < radius:
                coords.append(c)
    return np.array(coords)


def get_acid_carbon_and_oxygen(mol):
    """Return (acid_C_idx, carboxylate_O_idx) for the first COOH/COO- match."""
    pat = Chem.MolFromSmarts("[CX3](=O)[O;H1,-1]")
    matches = mol.GetSubstructMatches(pat)
    if not matches:
        return None
    # (C, O_double, O_H_or_anion)
    c, o1, o2 = matches[0]
    # prefer the acidic O (single bond)
    o_acidic = o2
    return c, o_acidic


def get_aromatic_centroid_idx(mol, conf):
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            coords = np.array([list(conf.GetAtomPosition(i)) for i in ring])
            return coords.mean(axis=0), ring
    return None, None


def rotation_matrix_axis(axis, theta):
    a = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    C = 1.0 - c
    x, y, z = a
    return np.array([
        [c + x*x*C, x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C],
    ])


def score_pose(ligand_coords, protein_coords, pocket_centroid,
               aromatic_target_xyz=None):
    """Return (score, pocket_fill, contact, clash, aromatic_d)."""
    # distance to pocket centroid
    d_pc = np.linalg.norm(ligand_coords - pocket_centroid, axis=1)
    pocket_fill = int(np.sum(d_pc < POCKET_SPHERE_RADIUS))

    # distances to protein atoms
    # broadcasting: (N_lig, 1, 3) - (1, N_prot, 3)
    diff = ligand_coords[:, None, :] - protein_coords[None, :, :]
    d_prot = np.linalg.norm(diff, axis=-1)
    min_d_per_lig_atom = d_prot.min(axis=1)

    clash = int(np.sum(min_d_per_lig_atom < CLASH_CUTOFF))
    contact = int(np.sum((min_d_per_lig_atom >= CLASH_CUTOFF) &
                         (min_d_per_lig_atom < CONTACT_CUTOFF)))

    total = len(ligand_coords)
    score = (pocket_fill + contact) / max(1, clash + 1) / total

    aromatic_d = None
    if aromatic_target_xyz is not None:
        # find ligand aromatic centroid (approximate by mean of atoms within
        # 1.5 Å of any protein-adjacent ring centroid → use just ligand mean)
        # Actually, we're re-passing aromatic target by pose — skip for now.
        pass

    return score, pocket_fill, contact, clash, aromatic_d


def dock_candidate(name, smiles, k1141_nz, pocket_centroid, protein_coords,
                   aromatic_target):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=N_CONFORMERS, params=params)
    if len(ids) == 0:
        return None

    try:
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=200)
    except Exception:
        pass

    ac_atoms = get_acid_carbon_and_oxygen(mol)
    if ac_atoms is None:
        return None
    c_idx, o_idx = ac_atoms

    # Target position for the acidic oxygen:
    #   salt-bridge distance from K1141 NZ, on the axis NZ → pocket centroid.
    axis = pocket_centroid - k1141_nz
    axis /= np.linalg.norm(axis)
    target_o_pos = k1141_nz + SALT_BRIDGE_DIST * axis

    best = None

    for cid in ids:
        conf = mol.GetConformer(cid)
        n = mol.GetNumAtoms()
        coords = np.array([list(conf.GetAtomPosition(i)) for i in range(n)])

        # only score heavy atoms
        heavy_idx = [i for i in range(n) if mol.GetAtomWithIdx(i).GetAtomicNum() != 1]
        heavy = coords[heavy_idx]

        # translate so acid O is at target
        o_pos = coords[o_idx]
        shift = target_o_pos - o_pos
        heavy_trans = heavy + shift

        # axis for rotation: around the (target_o - k1141_nz) direction
        # i.e. the salt-bridge bond — rotating keeps O in salt-bridge
        # but swings the rest of the ligand around that axis
        rot_axis = target_o_pos - k1141_nz
        rot_axis /= np.linalg.norm(rot_axis)

        for k in range(N_ROTATIONS):
            theta = 2 * np.pi * k / N_ROTATIONS
            R = rotation_matrix_axis(rot_axis, theta)
            # rotate heavy atoms about target_o_pos (pivot)
            heavy_rot = (heavy_trans - target_o_pos) @ R.T + target_o_pos

            score, pf, ct, cl, _ = score_pose(
                heavy_rot, protein_coords, pocket_centroid, aromatic_target
            )
            if (best is None) or (score > best["score"]):
                # aromatic distance from pose centroid — approximate by ring centroid
                aromatic_d = None
                ri = mol.GetRingInfo()
                for ring in ri.AtomRings():
                    if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                        ring_heavy = [heavy_idx.index(i) for i in ring if i in heavy_idx]
                        ring_c = heavy_rot[ring_heavy].mean(axis=0)
                        aromatic_d = float(np.linalg.norm(ring_c - aromatic_target))
                        break
                best = {
                    "conformer": cid,
                    "rotation": k,
                    "theta_deg": float(np.degrees(theta)),
                    "score": float(score),
                    "pocket_fill": pf,
                    "contact": ct,
                    "clash": cl,
                    "heavy_atoms": len(heavy_idx),
                    "aromatic_d_to_F1646_A": aromatic_d,
                    "acid_c_xyz": [float(x) for x in (heavy_rot[heavy_idx.index(c_idx)] if c_idx in heavy_idx else heavy_rot[0])],
                }

    return best


def main():
    pharmaco = json.loads(PHARMACO_JSON.read_text())
    screen = json.loads(SCREEN_JSON.read_text())

    pocket_centroid = np.array(pharmaco["pocket_centroid_A"])
    k1141_nz_list = pharmaco["ligand_design_rules"]["required_anion_xyz_target"]
    if not k1141_nz_list:
        raise ValueError("K1141 NZ not found in pharmacophore")
    k1141_nz = np.array(k1141_nz_list[0])

    aromatic_list = pharmaco["ligand_design_rules"]["preferred_aromatic_xyz_target"]
    aromatic_target = np.array(aromatic_list[0]) if aromatic_list else None

    protein_coords = load_protein_atoms_near(WT_CIF, pocket_centroid)
    print(f"Protein context atoms: {len(protein_coords)}")
    print(f"Pocket centroid:  {pocket_centroid}")
    print(f"K1141 NZ xyz:     {k1141_nz}")
    print(f"F1646 ring xyz:   {aromatic_target}")

    # Take top-10 candidates with anion
    candidates = [c for c in screen["ranked_candidates"] if c["has_anion"]][:10]

    results = []
    for c in candidates:
        print(f"\n→ docking {c['name']} ...", end=" ", flush=True)
        best = dock_candidate(
            c["name"], c["smiles"], k1141_nz, pocket_centroid,
            protein_coords, aromatic_target,
        )
        if best is None:
            print("FAILED")
            continue
        out = {
            "category": c["category"],
            "name": c["name"],
            "smiles": c["smiles"],
            "composite_3B": c["composite"],
            "anion_arom_d_ligand_A": c["anion_arom_d"],
            "best_pose": best,
        }
        results.append(out)
        print(f"score={best['score']:.3f} fill={best['pocket_fill']} "
              f"contact={best['contact']} clash={best['clash']} "
              f"d_arom={best['aromatic_d_to_F1646_A']:.2f}"
              if best['aromatic_d_to_F1646_A'] is not None
              else f"score={best['score']:.3f} fill={best['pocket_fill']} "
                   f"contact={best['contact']} clash={best['clash']}")

    # Final composite: Phase 3B composite × shape score
    for r in results:
        r["final_score"] = r["composite_3B"] * r["best_pose"]["score"]
    results.sort(key=lambda x: -x["final_score"])

    out_json = OUT_DIR / "pharmacochaperone_phase3c_shape_fit.json"
    out_json.write_text(json.dumps({
        "pocket_centroid_A": pocket_centroid.tolist(),
        "k1141_nz_A": k1141_nz.tolist(),
        "f1646_ring_A": aromatic_target.tolist() if aromatic_target is not None else None,
        "protein_atoms_context": len(protein_coords),
        "n_conformers_per_lig": N_CONFORMERS,
        "n_rotations_per_conf": N_ROTATIONS,
        "ranked_poses": results,
    }, indent=2))

    print("\n=== Final shape-fit ranking ===")
    print(f"{'Rank':<5}{'Name':<32}{'3B':>7}{'Shape':>8}{'Final':>8}"
          f"{'Fill':>5}{'Clash':>6}{'d_arom':>8}")
    print("-" * 80)
    for i, r in enumerate(results, 1):
        p = r["best_pose"]
        d_arom = f"{p['aromatic_d_to_F1646_A']:.2f}" if p["aromatic_d_to_F1646_A"] is not None else "--"
        print(f"{i:<5}{r['name']:<32}{r['composite_3B']:>7.3f}"
              f"{p['score']:>8.3f}{r['final_score']:>8.3f}"
              f"{p['pocket_fill']:>5}{p['clash']:>6}{d_arom:>8}")

    print(f"\n→ wrote {out_json}")


if __name__ == "__main__":
    main()
