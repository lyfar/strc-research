"""
Pharmacochaperone Phase 2: druggable-pocket scan on WT Job 4 CIF
around the mutation-responsive loop L1642-P-G-G-F-G-P-G-N1651
(identified in Phase 1B as the true structural signal of E1659A).

Strategy (LIGSITE-inspired, standalone numpy implementation):
  1. Load full-length WT chain A (Job 4).
  2. Focus region: 22 Å sphere around loop-1639-1651 C-alpha centroid.
  3. Grid (1 Å) → for each grid cell within region:
       - reject if inside vdW of any heavy atom (d ≤ 2.8 Å)
       - enclosure: cast 14 rays, count rays that hit protein ≤ 10 Å
       - keep if enclosure rays ≥ 8
  4. Cluster surviving grid cells (scipy.ndimage.label, 26-connectivity).
  5. For each cluster ≥ 80 Å³:
       - volume, lining residues (≤5 Å), hydrophobic fraction,
         H-bond donors/acceptors, distance to loop centroid,
         distance to E1659, druggability score.
  6. Export top pocket as PDB for visualisation.

Rationale: Phase 1B showed E1659A rearranges a proline-rich loop 17 Å from
TMEM145. Class-VX-809 tertiary stabilisers (lumacaftor etc.) bind allosteric
pockets that cage a misfolding loop. We need a real binding hole near the
loop — not the mutation site itself.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from scipy.ndimage import label as nd_label
from scipy.spatial import cKDTree

# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
MODELS_DIR = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models")
WT_CIF = MODELS_DIR / "job4-wildtype.cif"
OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")

LOOP_RANGE = (1639, 1651)        # mutation-responsive region from Phase 1B
E1659_RESNUM = 1659
SCAN_RADIUS_A = 22.0

GRID_SPACING = 1.0               # Å per voxel
VDW_CUTOFF = 2.8                 # reject grid cells inside heavy-atom vdW
RAY_LENGTH_A = 10.0              # probe distance for enclosure
RAYS_HIT_THRESHOLD = 8           # of 14 directions
MIN_CLUSTER_VOXELS = 80          # ≥ 80 Å³

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "HIS"}
CHARGED = {"LYS", "ARG", "ASP", "GLU"}

DONORS = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"],
          "TRP": ["NE1"], "ASN": ["ND2"], "GLN": ["NE2"],
          "HIS": ["ND1", "NE2"], "SER": ["OG"], "THR": ["OG1"],
          "TYR": ["OH"]}
ACCEPTORS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"],
             "ASN": ["OD1"], "GLN": ["OE1"], "HIS": ["ND1", "NE2"],
             "SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"]}

# 14 ray directions: 6 cardinal + 8 diagonals (corners of unit cube)
_RAYS = []
for dx in (-1, 0, 1):
    for dy in (-1, 0, 1):
        for dz in (-1, 0, 1):
            if (dx, dy, dz) == (0, 0, 0):
                continue
            # use only axis-aligned + body-diagonal (LIGSITE-style 7 axis)
            non_zero = (dx != 0) + (dy != 0) + (dz != 0)
            if non_zero in (1, 3):
                v = np.array([dx, dy, dz], dtype=float)
                v /= np.linalg.norm(v)
                _RAYS.append(v)
RAY_DIRS = np.stack(_RAYS, axis=0)  # shape (14, 3)


# --------------------------------------------------------------------------- #
# Structure parsing
# --------------------------------------------------------------------------- #
def load_chain(cif_path: Path, chain_id: str = "A"):
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", cif_path)
    model = next(struct.get_models())
    return model[chain_id]


def heavy_atom_coords(chain):
    coords, meta = [], []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        resnum = res.id[1]
        resname = res.get_resname()
        for atom in res:
            if atom.element == "H":
                continue
            coords.append(atom.coord)
            meta.append((resnum, resname, atom.get_name()))
    return np.array(coords), meta


def loop_centroid(chain, start: int, end: int):
    cas = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        if start <= res.id[1] <= end and "CA" in res:
            cas.append(res["CA"].coord)
    if not cas:
        raise RuntimeError("no loop CAs found")
    return np.mean(cas, axis=0)


def residue_ca(chain, resnum: int):
    for res in chain:
        if is_aa(res, standard=True) and res.id[1] == resnum and "CA" in res:
            return res["CA"].coord
    raise KeyError(resnum)


# --------------------------------------------------------------------------- #
# Pocket scan
# --------------------------------------------------------------------------- #
def build_grid(center: np.ndarray, radius: float, spacing: float):
    r = np.arange(-radius, radius + spacing, spacing)
    gx, gy, gz = np.meshgrid(r, r, r, indexing="ij")
    grid = np.stack([gx, gy, gz], axis=-1) + center  # shape (nx,ny,nz,3)
    return grid


def scan_pockets(grid_xyz: np.ndarray, atom_tree: cKDTree):
    """Return boolean mask over grid of cavity voxels."""
    shape = grid_xyz.shape[:3]
    flat = grid_xyz.reshape(-1, 3)

    # exclude voxels inside vdW of any atom
    vdw_neighbors = atom_tree.query_ball_point(flat, VDW_CUTOFF)
    outside_vdw = np.array([len(n) == 0 for n in vdw_neighbors], dtype=bool)

    mask = np.zeros(len(flat), dtype=bool)
    candidate_idx = np.where(outside_vdw)[0]
    step = 0.5
    steps = np.arange(step, RAY_LENGTH_A + step, step)  # sample points along ray

    for i in candidate_idx:
        origin = flat[i]
        rays_hit = 0
        for direction in RAY_DIRS:
            probe_points = origin + np.outer(steps, direction)  # (n_steps, 3)
            # fast: query single nearest
            d, _ = atom_tree.query(probe_points, k=1)
            if np.any(d <= 1.8):  # ~vdW of heavy atom
                rays_hit += 1
        if rays_hit >= RAYS_HIT_THRESHOLD:
            mask[i] = True
    return mask.reshape(shape)


def cluster_pocket(mask: np.ndarray, grid_xyz: np.ndarray, min_voxels: int):
    # 26-connectivity
    structure = np.ones((3, 3, 3), dtype=int)
    labels, n = nd_label(mask, structure=structure)
    clusters = []
    for k in range(1, n + 1):
        idx = np.where(labels == k)
        voxels = len(idx[0])
        if voxels < min_voxels:
            continue
        pts = grid_xyz[idx]
        clusters.append({"label": k, "voxels": voxels, "points": pts})
    return clusters


def lining_residues(pocket_pts: np.ndarray, chain, cutoff: float = 5.0):
    """Return set of (resnum, resname) within cutoff Å of any pocket point."""
    residues = {}
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        for atom in res:
            if atom.element == "H":
                continue
            d = np.linalg.norm(pocket_pts - atom.coord, axis=1)
            if np.any(d <= cutoff):
                residues[res.id[1]] = res.get_resname()
                break
    return residues


def druggability(volume_A3: float, hydrophobic_frac: float, nres: int,
                 donors: int, acceptors: int):
    """DoGSiteScorer-inspired heuristic (simplified, 0..1)."""
    v_opt = 300.0  # peak at 300 Å³
    v_score = np.exp(-((volume_A3 - v_opt) ** 2) / (2 * 200.0 ** 2))
    # hydrophobic 0.4–0.7 optimal
    if 0.4 <= hydrophobic_frac <= 0.7:
        h_score = 1.0
    else:
        h_score = max(0.0, 1.0 - 3 * abs(hydrophobic_frac - 0.55))
    # lining residues 12–25 optimal
    if 12 <= nres <= 25:
        r_score = 1.0
    else:
        r_score = max(0.0, 1.0 - abs(nres - 18) / 18)
    # at least 2 H-bond donors + 2 acceptors
    hb_score = min(1.0, (donors + acceptors) / 6.0)
    return float(0.40 * v_score + 0.25 * h_score + 0.20 * r_score + 0.15 * hb_score)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    print(f"[load] {WT_CIF.name}")
    chain = load_chain(WT_CIF, "A")
    all_atoms, atom_meta = heavy_atom_coords(chain)
    atom_tree = cKDTree(all_atoms)

    loop_ctr = loop_centroid(chain, *LOOP_RANGE)
    print(f"[loop] CA centroid 1639–1651 @ {loop_ctr}")
    e1659_ca = residue_ca(chain, E1659_RESNUM)

    # build grid centred on loop; slightly larger than scan radius
    grid = build_grid(loop_ctr, SCAN_RADIUS_A, GRID_SPACING)

    # restrict scan to points within SCAN_RADIUS_A from loop centroid
    rel = grid - loop_ctr
    dist_to_loop = np.linalg.norm(rel, axis=-1)
    in_sphere = dist_to_loop <= SCAN_RADIUS_A
    print(f"[grid] {in_sphere.sum()} voxels inside {SCAN_RADIUS_A:.1f} Å sphere "
          f"(spacing {GRID_SPACING} Å)")

    # zero-out voxels outside sphere by pre-filtering flat
    # Simpler: run scan on whole grid, then AND with sphere
    print("[scan] grid cavity detection (this may take ~1-2 min)...")
    pocket_mask = scan_pockets(grid, atom_tree)
    pocket_mask &= in_sphere
    print(f"[scan] {pocket_mask.sum()} cavity voxels (V≈{pocket_mask.sum()*GRID_SPACING**3:.0f} Å³ raw)")

    clusters = cluster_pocket(pocket_mask, grid, MIN_CLUSTER_VOXELS)
    print(f"[cluster] {len(clusters)} clusters ≥ {MIN_CLUSTER_VOXELS} voxels")

    # --- score each cluster --------------------------------------------- #
    scored = []
    for c in clusters:
        pts = c["points"]
        centroid = pts.mean(axis=0)
        d_loop = float(np.linalg.norm(centroid - loop_ctr))
        d_e1659 = float(np.linalg.norm(centroid - e1659_ca))

        lining = lining_residues(pts, chain, cutoff=5.0)
        n_res = len(lining)

        hydro = sum(1 for r in lining.values() if r in HYDROPHOBIC)
        pol = sum(1 for r in lining.values() if r in POLAR)
        charg = sum(1 for r in lining.values() if r in CHARGED)
        hydro_frac = hydro / max(1, n_res)

        donors = 0
        acceptors = 0
        for resnum, resname in lining.items():
            if resname in DONORS:
                donors += len(DONORS[resname])
            if resname in ACCEPTORS:
                acceptors += len(ACCEPTORS[resname])

        vol = c["voxels"] * GRID_SPACING ** 3
        drug = druggability(vol, hydro_frac, n_res, donors, acceptors)

        scored.append({
            "voxels": c["voxels"],
            "volume_A3": float(vol),
            "centroid": [float(x) for x in centroid],
            "dist_to_loop_centroid_A": round(d_loop, 2),
            "dist_to_E1659_CA_A": round(d_e1659, 2),
            "n_lining_residues": n_res,
            "hydrophobic_frac": round(hydro_frac, 3),
            "charged_count": charg,
            "polar_count": pol,
            "h_bond_donors": donors,
            "h_bond_acceptors": acceptors,
            "druggability": round(drug, 3),
            "lining_residues": sorted(lining.keys()),
            "_points": pts,
        })

    scored.sort(key=lambda s: -s["druggability"])

    # --- top cluster -------------------------------------------------- #
    if scored:
        top = scored[0]
        print(f"\n[top pocket] V={top['volume_A3']:.0f} Å³ | "
              f"druggability={top['druggability']:.2f} | "
              f"d_loop={top['dist_to_loop_centroid_A']:.1f} Å | "
              f"d_E1659={top['dist_to_E1659_CA_A']:.1f} Å | "
              f"lining={top['n_lining_residues']} res")
        # export as fake PDB: pocket points as HETATM DUM
        top_pdb = OUT_DIR / "pharmacochaperone_phase2_top_pocket.pdb"
        lines = []
        for i, p in enumerate(top["_points"]):
            lines.append(
                f"HETATM{i+1:5d}  C   DUM A{i%9999:4d}    "
                f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00 20.00           C"
            )
        top_pdb.write_text("\n".join(lines) + "\nEND\n")
        print(f"[write] {top_pdb.name}")

    # strip "_points" for JSON
    for s in scored:
        s.pop("_points", None)

    # --- report ------------------------------------------------------- #
    summary = {
        "wt_cif": WT_CIF.name,
        "loop_range": list(LOOP_RANGE),
        "scan_radius_A": SCAN_RADIUS_A,
        "grid_spacing_A": GRID_SPACING,
        "loop_centroid_A": [float(x) for x in loop_ctr],
        "e1659_ca_A": [float(x) for x in e1659_ca],
        "n_clusters": len(scored),
        "clusters": scored[:8],
        "verdict": "see_cluster_ranking",
    }

    if scored:
        top = scored[0]
        if top["druggability"] >= 0.5 and top["dist_to_loop_centroid_A"] <= 12:
            verdict = "DRUGGABLE ALLOSTERIC POCKET WITHIN 12 Å OF LOOP — advance to docking"
        elif top["druggability"] >= 0.4:
            verdict = "MARGINAL POCKET — fragment screen feasible, avoid early HTS"
        else:
            verdict = "NO DRUGGABLE POCKET NEAR LOOP — reconsider pharmacochaperone strategy"
        summary["verdict"] = verdict
        print(f"\n[verdict] {verdict}")
    else:
        summary["verdict"] = "NO CLUSTERS FOUND — loop too solvent-exposed"
        print("\n[verdict] NO CLUSTERS FOUND")

    out_json = OUT_DIR / "pharmacochaperone_phase2_results.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"[write] {out_json.name}")

    # --- figure ------------------------------------------------------- #
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))

    if scored:
        top_pts = top["_points"] if "_points" in top else None  # may be popped
        # re-find top points by centroid match
        pts_all = np.concatenate([
            grid[pocket_mask]
        ], axis=0) if pocket_mask.any() else np.zeros((0, 3))

        # residue CA cloud (within 25 Å of loop)
        res_ca = []
        res_nums = []
        for res in chain:
            if not is_aa(res, standard=True):
                continue
            if "CA" not in res:
                continue
            d = np.linalg.norm(res["CA"].coord - loop_ctr)
            if d <= 25:
                res_ca.append(res["CA"].coord)
                res_nums.append(res.id[1])
        res_ca = np.array(res_ca)

        ax[0].scatter(res_ca[:, 0], res_ca[:, 2], s=8, c="lightgray",
                      label="chain A CA (<25 Å of loop)")
        if pts_all.size:
            ax[0].scatter(pts_all[:, 0], pts_all[:, 2], s=2, c="tab:red",
                          alpha=0.45, label="cavity voxels")
        # loop CAs
        loop_cas = np.array([res["CA"].coord for res in chain
                             if is_aa(res, standard=True) and
                             LOOP_RANGE[0] <= res.id[1] <= LOOP_RANGE[1]
                             and "CA" in res])
        ax[0].scatter(loop_cas[:, 0], loop_cas[:, 2], s=40,
                      c="tab:blue", marker="D", label="loop 1639–1651")
        ax[0].scatter([e1659_ca[0]], [e1659_ca[2]], s=80, c="black",
                      marker="*", label="E1659 CA")
        ax[0].set_xlabel("x / Å"); ax[0].set_ylabel("z / Å")
        ax[0].set_title("Loop neighbourhood (WT Job 4) — XZ slab")
        ax[0].legend(loc="upper right", fontsize=8)
        ax[0].set_aspect("equal", "datalim")

        # bar chart of top-8 cluster scores
        labels = [f"#{i+1}" for i in range(len(scored[:8]))]
        drug_scores = [s["druggability"] for s in scored[:8]]
        vols = [s["volume_A3"] for s in scored[:8]]
        dlocs = [s["dist_to_loop_centroid_A"] for s in scored[:8]]
        ax2 = ax[1]
        x = np.arange(len(labels))
        bars = ax2.bar(x, drug_scores, color="tab:orange",
                       edgecolor="black", alpha=0.8)
        ax2.axhline(0.5, color="tab:red", ls="--", lw=1, label="druggable ≥0.5")
        ax2.axhline(0.4, color="tab:gray", ls=":", lw=1, label="marginal ≥0.4")
        for i, b in enumerate(bars):
            ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.01,
                     f"{vols[i]:.0f} Å³\nd={dlocs[i]:.1f}Å",
                     ha="center", fontsize=7)
        ax2.set_xticks(x)
        ax2.set_xticklabels(labels)
        ax2.set_ylim(0, 1.0)
        ax2.set_ylabel("druggability score")
        ax2.set_title("Top cavities near loop 1639–1651")
        ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "pharmacochaperone_phase2.png", dpi=160)
    print(f"[write] pharmacochaperone_phase2.png")


if __name__ == "__main__":
    main()
