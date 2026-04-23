"""h26 Phase 1c — spatial re-cluster of Ultra-Mini homodimer inter-chain contacts.

Phase 1 mutagenesis targeted only the 1579-1581 ARM sub-cluster and failed
(all 4 mutants destabilized dimer further per [[STRC Engineered Homodimer Phase 1 Results]]).

Prior consensus analysis ([[ultramini_homodimer_consensus.json]]) revealed TWO
independent consensus clusters at 5 A cutoff:
  A) "stump"    1077-1131  — 13 strict-consensus residues (N-term of Ultra-Mini)
  B) "deep ARM" 1493-1590  — 16 strict-consensus residues (C-term ARM region)

Phase 1c does geometric spatial clustering on the full contact set (not just
consensus) and ranks mutation candidates per cluster, so Phase 1b can retry with
targets outside 1579-1581.

Method
  1. Parse job-ultramini-homodimer.cif (offset 1074 -> STRC canonical).
  2. NeighborSearch 5 A on chain A vs chain B, residue-pair min distances.
  3. Extract Ca coords of A-side contact residues.
  4. DBSCAN spatial cluster at eps 10 A, min 3 residues.
  5. Per cluster: residues, COM, C2-symmetry (does the partner residue on chain
     B sit in a symmetric cluster?), residue-type composition.
  6. Mutation candidate ranking:
     - bulky aromatic replacement of small hydrophobic (steric disruption)
     - charge swap (R/K <-> D/E) if neighboring oppositely charged partner
     - Cys introduction for disulfide engineering if A-side and B-side Cbeta
       within 4.5-7.5 A (C2-symmetric pair -> inter-chain disulfide)

Output: engineered_homodimer_phase1c_contact_cluster.json
"""
import json
import warnings
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch

warnings.filterwarnings("ignore")

CIF = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job-ultramini-homodimer.cif")
OUT = Path("/Users/egorlyfar/Brain/research/strc/models/engineered_homodimer_phase1c_contact_cluster.json")

CONTACT_CUTOFF = 5.0         # A, heavy-atom contact definition
CLUSTER_EPS = 6.5            # A, DBSCAN neighborhood on Ca coords (tuned to split stump vs deep-ARM)
CLUSTER_MIN_PTS = 3          # DBSCAN min samples
DISULFIDE_MIN = 4.5          # A, Cb-Cb minimum for Cys engineering
DISULFIDE_MAX = 7.5          # A, Cb-Cb maximum (native S-S typical ~5-6 A Cb-Cb)

# Ultra-Mini construct: chain residues 1..701 == STRC canonical 1075..1775
OFFSET = 1074                # canonical = local + OFFSET

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}
CHARGED_POS = {"R", "K", "H"}
CHARGED_NEG = {"D", "E"}
POLAR = {"S", "T", "N", "Q", "Y"}
HYDROPHOBIC = {"A", "V", "L", "I", "M", "F", "W", "P", "C", "G"}
AROMATIC = {"F", "W", "Y", "H"}
SMALL_HYDROPHOBIC = {"A", "V", "L", "I"}


def dbscan(points: np.ndarray, eps: float, min_pts: int) -> np.ndarray:
    """Minimal DBSCAN implementation.  Labels: -1 = noise, 0..k = clusters."""
    n = len(points)
    labels = -np.ones(n, dtype=int)
    visited = np.zeros(n, dtype=bool)
    dist_sq = np.sum((points[:, None, :] - points[None, :, :]) ** 2, axis=-1)
    neighborhoods = [np.where(dist_sq[i] <= eps * eps)[0] for i in range(n)]
    cluster_id = -1
    for i in range(n):
        if visited[i]:
            continue
        visited[i] = True
        if len(neighborhoods[i]) < min_pts:
            continue  # noise for now
        cluster_id += 1
        labels[i] = cluster_id
        queue = list(neighborhoods[i])
        while queue:
            j = queue.pop(0)
            if not visited[j]:
                visited[j] = True
                if len(neighborhoods[j]) >= min_pts:
                    queue.extend([k for k in neighborhoods[j] if not visited[k]])
            if labels[j] == -1:
                labels[j] = cluster_id
    return labels


def residue_charge_class(one_letter: str) -> str:
    if one_letter in CHARGED_POS:
        return "pos"
    if one_letter in CHARGED_NEG:
        return "neg"
    if one_letter in POLAR:
        return "polar"
    if one_letter in AROMATIC:
        return "aromatic"
    if one_letter in HYDROPHOBIC:
        return "hydrophobic"
    return "other"


def parse_homodimer(cif_path: Path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("x", cif_path)
    model = next(structure.get_models())
    chain_A = model["A"]
    chain_B = model["B"]
    atoms_A = [a for a in chain_A.get_atoms() if a.element != "H"]
    atoms_B = [a for a in chain_B.get_atoms() if a.element != "H"]
    return model, chain_A, chain_B, atoms_A, atoms_B


def residue_pair_distances(atoms_A, atoms_B):
    """Return {(canon_A, canon_B): (min_dist, A-restype, B-restype)}."""
    ns_B = NeighborSearch(atoms_B)
    pairs = {}
    for a in atoms_A:
        for b in ns_B.search(a.coord, CONTACT_CUTOFF):
            rA = a.get_parent()
            rB = b.get_parent()
            key = (rA.id[1] + OFFSET, rB.id[1] + OFFSET)
            d = a - b
            if key not in pairs or d < pairs[key][0]:
                pairs[key] = (d, rA.get_resname(), rB.get_resname())
    return pairs


def ca_coord(chain, local_id: int):
    r = chain[local_id]
    return np.asarray(r["CA"].coord, dtype=float)


def cb_or_ca_coord(chain, local_id: int):
    r = chain[local_id]
    if "CB" in r:
        return np.asarray(r["CB"].coord, dtype=float)
    return np.asarray(r["CA"].coord, dtype=float)


def main():
    model, chain_A, chain_B, atoms_A, atoms_B = parse_homodimer(CIF)
    pairs = residue_pair_distances(atoms_A, atoms_B)
    print(f"inter-chain residue pairs at <{CONTACT_CUTOFF} A: {len(pairs)}")

    # A-side contact residues (unique)
    a_contacts = {}
    for (rA, rB), (d, resname_A, resname_B) in pairs.items():
        if rA not in a_contacts:
            a_contacts[rA] = resname_A
    a_residues = sorted(a_contacts)
    print(f"unique A-side contact residues: {len(a_residues)}")

    # Cluster A-side residues on Ca coords
    coords = np.stack([ca_coord(chain_A, r - OFFSET) for r in a_residues])
    labels = dbscan(coords, CLUSTER_EPS, CLUSTER_MIN_PTS)
    n_clusters = labels.max() + 1 if labels.max() >= 0 else 0
    print(f"DBSCAN (eps={CLUSTER_EPS} A, min_pts={CLUSTER_MIN_PTS}): {n_clusters} clusters, {(labels==-1).sum()} noise")

    # Per-cluster summary
    clusters_summary = []
    for c in range(n_clusters):
        mask = labels == c
        res_canon = [a_residues[i] for i in range(len(a_residues)) if mask[i]]
        res_types = [a_contacts[r] for r in res_canon]
        res_one = [THREE_TO_ONE.get(t, "X") for t in res_types]
        com = coords[mask].mean(axis=0)
        composition = defaultdict(int)
        for o in res_one:
            composition[residue_charge_class(o)] += 1

        # C2-symmetric partner test: for each A-side residue r, does B-side r
        # also sit in a contact (homotypic) and in the same cluster box?
        c2_partner_count = 0
        homotypic_pair_count = 0
        for r in res_canon:
            b_partner_local = r - OFFSET
            # Find all B-side residues that pair with this A-side residue
            b_partners_canon = sorted(rB for (rA, rB) in pairs if rA == r)
            for rB in b_partners_canon:
                # Check whether rB (B-side) maps back to a position inside the
                # same cluster's canonical resid range -> C2-symmetric
                if min(res_canon) - 2 <= rB <= max(res_canon) + 2:
                    c2_partner_count += 1
                    if rB == r:
                        homotypic_pair_count += 1
                    break

        # Disulfide engineering candidates: C2-symmetric Cb-Cb within 4.5-7.5 A
        disulfide_candidates = []
        for r in res_canon:
            for rB in sorted(rBB for (rA, rBB) in pairs if rA == r):
                if not (min(res_canon) - 2 <= rB <= max(res_canon) + 2):
                    continue
                cb_a = cb_or_ca_coord(chain_A, r - OFFSET)
                cb_b = cb_or_ca_coord(chain_B, rB - OFFSET)
                d_cb = float(np.linalg.norm(cb_a - cb_b))
                if DISULFIDE_MIN <= d_cb <= DISULFIDE_MAX:
                    disulfide_candidates.append({
                        "A_res": r, "A_type": a_contacts[r],
                        "B_res": rB, "B_type": pairs[(r, rB)][2],
                        "cb_cb_A": round(d_cb, 2),
                    })

        # Mutation candidate ranking (steric / electrostatic / covalent)
        mutation_candidates = []
        for i, r in enumerate(res_canon):
            one = res_one[i]
            rationales = []
            if one in SMALL_HYDROPHOBIC:
                rationales.append(f"{one}->W or {one}->F: steric bulk disrupts packing")
            if one in CHARGED_POS:
                # look for adjacent pos partner on B -> charge-swap mutation
                for rB in sorted(rBB for (rA, rBB) in pairs if rA == r):
                    if not (min(res_canon) - 2 <= rB <= max(res_canon) + 2):
                        continue
                    rB_one = THREE_TO_ONE.get(pairs[(r, rB)][2], "X")
                    if rB_one in CHARGED_POS:
                        rationales.append(f"{one}->D/E: relieve electrostatic repulsion vs {rB_one}{rB}")
                    if rB_one in CHARGED_NEG:
                        rationales.append(f"{one}->{rB_one}: break inter-chain salt bridge (destabilize)")
            if one in CHARGED_NEG:
                for rB in sorted(rBB for (rA, rBB) in pairs if rA == r):
                    if not (min(res_canon) - 2 <= rB <= max(res_canon) + 2):
                        continue
                    rB_one = THREE_TO_ONE.get(pairs[(r, rB)][2], "X")
                    if rB_one in CHARGED_NEG:
                        rationales.append(f"{one}->R/K: relieve repulsion vs {rB_one}{rB}")
            if one in {"S", "T", "A", "V"}:
                # Cys engineering candidate if C2-symmetric
                for rB in sorted(rBB for (rA, rBB) in pairs if rA == r):
                    if rB == r:
                        cb_a = cb_or_ca_coord(chain_A, r - OFFSET)
                        cb_b = cb_or_ca_coord(chain_B, rB - OFFSET)
                        d_cb = float(np.linalg.norm(cb_a - cb_b))
                        if DISULFIDE_MIN <= d_cb <= DISULFIDE_MAX:
                            rationales.append(f"{one}->C / {one}->C: inter-chain disulfide (Cb-Cb {d_cb:.2f} A)")
            if rationales:
                mutation_candidates.append({
                    "res": r, "type": one,
                    "rationales": rationales,
                })

        cluster = {
            "id": c,
            "n_residues": int(mask.sum()),
            "residues_canonical": res_canon,
            "residue_types": {r: a_contacts[r] for r in res_canon},
            "center_of_mass_A": [round(float(com[i]), 2) for i in range(3)],
            "composition": dict(composition),
            "c2_symmetry_partner_count": c2_partner_count,
            "homotypic_pair_count": homotypic_pair_count,
            "disulfide_candidates": disulfide_candidates,
            "mutation_candidates": mutation_candidates,
            "zone_label": label_zone(res_canon),
        }
        clusters_summary.append(cluster)

    # Sort clusters by size descending
    clusters_summary.sort(key=lambda c: -c["n_residues"])
    for c in clusters_summary:
        c["rank"] = clusters_summary.index(c) + 1

    # Cluster -> literature fit
    for c in clusters_summary:
        residues = c["residues_canonical"]
        lo, hi = min(residues), max(residues)
        c["phase1_overlap"] = (1579 in residues) or (1580 in residues) or (1581 in residues)

    verdict = derive_verdict(clusters_summary)

    result = {
        "cif": str(CIF),
        "contact_cutoff_A": CONTACT_CUTOFF,
        "dbscan": {"eps_A": CLUSTER_EPS, "min_pts": CLUSTER_MIN_PTS},
        "total_inter_chain_residue_pairs": len(pairs),
        "total_a_side_contact_residues": len(a_residues),
        "n_clusters": n_clusters,
        "n_noise_residues": int((labels == -1).sum()),
        "clusters": clusters_summary,
        "verdict": verdict,
    }
    OUT.write_text(json.dumps(result, indent=2, default=float))

    # Pretty print
    print()
    for c in clusters_summary:
        residues = c["residues_canonical"]
        print(f"\ncluster {c['rank']} [{c['zone_label']}] — {c['n_residues']} residues "
              f"{min(residues)}-{max(residues)}")
        print(f"  composition: {c['composition']}")
        print(f"  C2-symmetric pair count: {c['c2_symmetry_partner_count']}, "
              f"homotypic (rA=rB): {c['homotypic_pair_count']}")
        if c["disulfide_candidates"]:
            print(f"  disulfide candidates: {len(c['disulfide_candidates'])}")
            for d in c["disulfide_candidates"][:3]:
                print(f"    A:{d['A_res']}{d['A_type']} <-> B:{d['B_res']}{d['B_type']}  Cb-Cb {d['cb_cb_A']} A")
        print(f"  overlap with Phase 1 (1579-1581): {c['phase1_overlap']}")
        if c["mutation_candidates"]:
            print(f"  top mutation candidates:")
            for m in c["mutation_candidates"][:5]:
                for r in m["rationales"]:
                    print(f"    {r} [res {m['res']} {m['type']}]")
    print(f"\nverdict: {verdict}")
    print(f"saved: {OUT}")


def label_zone(residues):
    lo, hi = min(residues), max(residues)
    if hi <= 1200:
        return "stump_1075_1200"
    if lo >= 1401:
        return "deep_arm_1401_1775"
    if lo >= 1201 and hi <= 1400:
        return "mid_1201_1400"
    return f"mixed_{lo}_{hi}"


def derive_verdict(clusters):
    # Gate: cluster outside 1579-1581 that is larger than the Phase 1 cluster
    outside = [c for c in clusters if not c["phase1_overlap"] and c["n_residues"] >= 5]
    with_disulfide = [c for c in clusters if c["disulfide_candidates"] and not c["phase1_overlap"]]
    if with_disulfide:
        return (f"phase_1c_green — {len(with_disulfide)} cluster(s) outside 1579-1581 "
                f"with disulfide engineering candidates; queue Phase 1d Cys-engineering AF3")
    if outside:
        return (f"phase_1c_yellow — {len(outside)} secondary cluster(s) outside 1579-1581 "
                f"(size >=5) without disulfide geometry; queue Phase 1b double-mutant AF3 "
                f"on largest")
    return "phase_1c_red — no productive secondary cluster; h26 moves to C-tier"


if __name__ == "__main__":
    main()
