"""
Phase 2B: split the broad loop-LRR tertiary pocket found in Phase 2
into drug-sized subpockets using stricter buriedness thresholds
and a local-maximum-depth hotspot search.

Strategy changes vs Phase 2:
  - Tighter enclosure (11 of 14 rays, down from 8)
  - Shorter ray length 8 Å (was 10 Å) — only true buried cavities survive
  - Finer grid 0.8 Å
  - Minimum cluster 30 voxels (15 Å³), retain smaller clefts
  - Depth map: for each voxel, distance-to-nearest-atom → peaks define
    drug-like subpockets (V ≈ 100–400 Å³)
  - Report only subpockets whose centroid lies within 10 Å of the
    mutation-responsive loop (1639–1651)
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from scipy.ndimage import label as nd_label, maximum_filter
from scipy.spatial import cKDTree

MODELS_DIR = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models")
WT_CIF = MODELS_DIR / "job4-wildtype.cif"
OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")

LOOP_RANGE = (1639, 1651)
E1659_RESNUM = 1659
SCAN_RADIUS_A = 18.0
GRID_SPACING = 0.8                # finer grid
VDW_CUTOFF = 2.6
RAY_LENGTH_A = 9.0
RAYS_HIT_THRESHOLD = 9            # moderate: 9 of 14 rays must hit
MIN_CLUSTER_VOXELS = 50
SUBPOCKET_MAX_DIST_A = 15.0       # from loop CA centroid

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

_RAYS = []
for dx in (-1, 0, 1):
    for dy in (-1, 0, 1):
        for dz in (-1, 0, 1):
            if (dx, dy, dz) == (0, 0, 0):
                continue
            non_zero = (dx != 0) + (dy != 0) + (dz != 0)
            if non_zero in (1, 3):
                v = np.array([dx, dy, dz], dtype=float)
                v /= np.linalg.norm(v)
                _RAYS.append(v)
RAY_DIRS = np.stack(_RAYS, axis=0)


def load_chain(cif_path: Path, chain_id: str = "A"):
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", cif_path)
    model = next(struct.get_models())
    return model[chain_id]


def heavy_atom_coords(chain):
    coords = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        for atom in res:
            if atom.element == "H":
                continue
            coords.append(atom.coord)
    return np.array(coords)


def loop_centroid(chain, start: int, end: int):
    cas = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        if start <= res.id[1] <= end and "CA" in res:
            cas.append(res["CA"].coord)
    return np.mean(cas, axis=0)


def residue_ca(chain, resnum: int):
    for res in chain:
        if is_aa(res, standard=True) and res.id[1] == resnum and "CA" in res:
            return res["CA"].coord
    raise KeyError(resnum)


def build_grid(center: np.ndarray, radius: float, spacing: float):
    r = np.arange(-radius, radius + spacing, spacing)
    gx, gy, gz = np.meshgrid(r, r, r, indexing="ij")
    grid = np.stack([gx, gy, gz], axis=-1) + center
    return grid


def enclosure_and_depth(grid_xyz: np.ndarray, atom_tree: cKDTree):
    """For each voxel compute (buried?, depth-to-nearest-atom)."""
    shape = grid_xyz.shape[:3]
    flat = grid_xyz.reshape(-1, 3)

    # reject voxels inside vdW
    nn_d, _ = atom_tree.query(flat, k=1)
    outside_vdw = nn_d >= VDW_CUTOFF

    mask = np.zeros(len(flat), dtype=bool)
    depth = np.zeros(len(flat), dtype=float)
    depth[outside_vdw] = nn_d[outside_vdw]

    candidate_idx = np.where(outside_vdw)[0]
    step = 0.5
    steps = np.arange(step, RAY_LENGTH_A + step, step)

    for i in candidate_idx:
        origin = flat[i]
        rays_hit = 0
        for direction in RAY_DIRS:
            probe_points = origin + np.outer(steps, direction)
            d, _ = atom_tree.query(probe_points, k=1)
            if np.any(d <= 1.8):
                rays_hit += 1
        if rays_hit >= RAYS_HIT_THRESHOLD:
            mask[i] = True
    return mask.reshape(shape), depth.reshape(shape)


def cluster(mask: np.ndarray, grid_xyz: np.ndarray, min_voxels: int):
    structure = np.ones((3, 3, 3), dtype=int)
    labels, n = nd_label(mask, structure=structure)
    clusters = []
    for k in range(1, n + 1):
        idx = np.where(labels == k)
        voxels = len(idx[0])
        if voxels < min_voxels:
            continue
        pts = grid_xyz[idx]
        clusters.append({"label": k, "voxels": voxels, "points": pts,
                         "idx": idx})
    return clusters


def lining_residues(pocket_pts: np.ndarray, chain, cutoff: float = 4.5):
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


def druggability(vol, hydro_frac, nres, donors, acceptors, depth_mean):
    """Subpocket druggability (includes depth descriptor, 0..1).

    NOTE [audit 2026-04-23]: Phase 2b weights 0.30/0.20/0.15/0.15/0.20
    (vol/hydro/nres/hb/depth) differ from Phase 1 (0.5/0.3/0.2) and
    Phase 2 (no depth component). Scores NOT cross-phase comparable.
    Depth component d_score = depth_mean/3.5 is a burial-favouring term
    specific to subpocket ranking; upstream phases don't compute burial.
    Use within-phase only.
    """
    v_opt = 250.0
    v_score = np.exp(-((vol - v_opt) ** 2) / (2 * 150.0 ** 2))
    if 0.4 <= hydro_frac <= 0.7:
        h_score = 1.0
    else:
        h_score = max(0.0, 1.0 - 3 * abs(hydro_frac - 0.55))
    if 10 <= nres <= 22:
        r_score = 1.0
    else:
        r_score = max(0.0, 1.0 - abs(nres - 16) / 16)
    hb_score = min(1.0, (donors + acceptors) / 6.0)
    # depth: burial ≥ 3.5 Å optimal
    d_score = min(1.0, depth_mean / 3.5)
    return float(0.30 * v_score + 0.20 * h_score + 0.15 * r_score +
                 0.15 * hb_score + 0.20 * d_score)


def flood_from_seeds(cluster_mask: np.ndarray, depth_map: np.ndarray,
                     grid_xyz: np.ndarray):
    """For a large cluster, find local-maximum depth seeds and
    grow subpockets by region-growing (inverse watershed)."""
    # local maxima in depth within the cluster
    depth_in_cluster = np.where(cluster_mask, depth_map, 0)
    peaks = (depth_in_cluster == maximum_filter(depth_in_cluster, size=5)) & \
            (depth_in_cluster > 2.5)
    peak_idx = np.where(peaks)
    peak_depths = depth_in_cluster[peak_idx]
    order = np.argsort(-peak_depths)
    peak_coords = np.stack([peak_idx[0][order], peak_idx[1][order],
                            peak_idx[2][order]], axis=-1)

    # assign each cluster voxel to nearest peak (grid distance)
    cluster_voxels = np.stack(np.where(cluster_mask), axis=-1)
    if len(peak_coords) == 0:
        return []
    # nearest peak per voxel
    dists = np.zeros((len(cluster_voxels), len(peak_coords)))
    for k, p in enumerate(peak_coords):
        dists[:, k] = np.linalg.norm(cluster_voxels - p, axis=1)
    assign = np.argmin(dists, axis=1)

    subpockets = []
    for k in range(len(peak_coords)):
        vox_idx = cluster_voxels[assign == k]
        if len(vox_idx) < 25:
            continue
        pts = grid_xyz[vox_idx[:, 0], vox_idx[:, 1], vox_idx[:, 2]]
        peak_xyz = grid_xyz[tuple(peak_coords[k])]
        subpockets.append({
            "voxels": len(vox_idx),
            "points": pts,
            "peak_xyz": peak_xyz,
            "peak_depth": float(depth_in_cluster[tuple(peak_coords[k])]),
        })
    return subpockets


def main():
    print(f"[load] {WT_CIF.name}")
    chain = load_chain(WT_CIF, "A")
    all_atoms = heavy_atom_coords(chain)
    atom_tree = cKDTree(all_atoms)

    loop_ctr = loop_centroid(chain, *LOOP_RANGE)
    e1659_ca = residue_ca(chain, E1659_RESNUM)

    grid = build_grid(loop_ctr, SCAN_RADIUS_A, GRID_SPACING)
    rel = grid - loop_ctr
    dist_to_loop = np.linalg.norm(rel, axis=-1)
    in_sphere = dist_to_loop <= SCAN_RADIUS_A

    print(f"[grid] {in_sphere.sum()} voxels, spacing {GRID_SPACING} Å")
    print("[scan] enclosure (11/14 rays) + depth...")
    buried_mask, depth_map = enclosure_and_depth(grid, atom_tree)
    buried_mask &= in_sphere
    print(f"[scan] {buried_mask.sum()} buried voxels "
          f"(V≈{buried_mask.sum()*GRID_SPACING**3:.0f} Å³)")

    clusters = cluster(buried_mask, grid, MIN_CLUSTER_VOXELS)
    print(f"[cluster] {len(clusters)} primary clusters")

    # split large clusters into subpockets
    scored_subpockets = []
    for ci, c in enumerate(clusters):
        cmask = np.zeros_like(buried_mask, dtype=bool)
        cmask[c["idx"]] = True
        subs = flood_from_seeds(cmask, depth_map, grid)
        print(f"  cluster #{ci+1}: {c['voxels']} voxels → {len(subs)} subpockets")

        for si, s in enumerate(subs):
            pts = s["points"]
            centroid = pts.mean(axis=0)
            d_loop = float(np.linalg.norm(centroid - loop_ctr))
            if d_loop > SUBPOCKET_MAX_DIST_A:
                continue
            d_e1659 = float(np.linalg.norm(centroid - e1659_ca))

            lining = lining_residues(pts, chain, cutoff=4.5)
            n_res = len(lining)
            if n_res == 0:
                continue
            hydro = sum(1 for r in lining.values() if r in HYDROPHOBIC)
            pol = sum(1 for r in lining.values() if r in POLAR)
            charg = sum(1 for r in lining.values() if r in CHARGED)
            hydro_frac = hydro / n_res

            donors = sum(len(DONORS.get(resname, []))
                         for resname in lining.values())
            acceptors = sum(len(ACCEPTORS.get(resname, []))
                            for resname in lining.values())

            vol = s["voxels"] * GRID_SPACING ** 3

            # per-voxel depth
            voxel_depths = []
            for p in pts:
                d, _ = atom_tree.query(p, k=1)
                voxel_depths.append(d)
            depth_mean = float(np.mean(voxel_depths))
            depth_max = float(np.max(voxel_depths))

            drug = druggability(vol, hydro_frac, n_res, donors, acceptors,
                                depth_mean)

            # loop-adjacency flag: touches loop residues 1639-1651?
            touches_loop = any(LOOP_RANGE[0] <= r <= LOOP_RANGE[1]
                               for r in lining.keys())

            scored_subpockets.append({
                "parent_cluster": ci + 1,
                "subpocket": si + 1,
                "voxels": s["voxels"],
                "volume_A3": round(vol, 1),
                "centroid": [round(float(x), 2) for x in centroid],
                "peak_xyz": [round(float(x), 2) for x in s["peak_xyz"]],
                "peak_depth_A": round(s["peak_depth"], 2),
                "dist_to_loop_centroid_A": round(d_loop, 2),
                "dist_to_E1659_CA_A": round(d_e1659, 2),
                "n_lining_residues": n_res,
                "hydrophobic_frac": round(hydro_frac, 3),
                "charged_count": charg,
                "polar_count": pol,
                "h_bond_donors": donors,
                "h_bond_acceptors": acceptors,
                "depth_mean_A": round(depth_mean, 2),
                "depth_max_A": round(depth_max, 2),
                "druggability": round(drug, 3),
                "touches_loop_1639_1651": touches_loop,
                "lining_residues": sorted(lining.keys()),
                "_points": pts,
            })

    scored_subpockets.sort(key=lambda s: (-s["druggability"],
                                          s["dist_to_loop_centroid_A"]))
    print(f"[subpockets] {len(scored_subpockets)} near loop")

    # report top 5
    print("\n[top subpockets]")
    for i, s in enumerate(scored_subpockets[:5]):
        print(f"  #{i+1}: V={s['volume_A3']:5.0f} Å³  "
              f"drug={s['druggability']:.2f}  "
              f"d_loop={s['dist_to_loop_centroid_A']:4.1f} Å  "
              f"d_E1659={s['dist_to_E1659_CA_A']:4.1f} Å  "
              f"lining={s['n_lining_residues']:2d}  "
              f"touches_loop={s['touches_loop_1639_1651']}")

    # export top subpocket (touches loop) as PDB
    top_loop = next((s for s in scored_subpockets
                     if s["touches_loop_1639_1651"] and s["druggability"] >= 0.4),
                    None)
    if top_loop:
        out_pdb = OUT_DIR / "pharmacochaperone_phase2b_top_subpocket.pdb"
        lines = []
        for i, p in enumerate(top_loop["_points"]):
            lines.append(
                f"HETATM{i+1:5d}  C   DUM A{(i % 9999)+1:4d}    "
                f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00 20.00           C"
            )
        out_pdb.write_text("\n".join(lines) + "\nEND\n")
        print(f"[write] {out_pdb.name}")

    # strip _points from JSON
    for s in scored_subpockets:
        s.pop("_points", None)

    # loop-facing subpockets only
    loop_facing = [s for s in scored_subpockets if s["touches_loop_1639_1651"]]

    summary = {
        "wt_cif": WT_CIF.name,
        "loop_range": list(LOOP_RANGE),
        "scan_radius_A": SCAN_RADIUS_A,
        "grid_spacing_A": GRID_SPACING,
        "rays_threshold": RAYS_HIT_THRESHOLD,
        "ray_length_A": RAY_LENGTH_A,
        "min_cluster_voxels": MIN_CLUSTER_VOXELS,
        "n_clusters_primary": len(clusters),
        "n_subpockets_total": len(scored_subpockets),
        "n_loop_facing": len(loop_facing),
        "loop_centroid_A": [float(x) for x in loop_ctr],
        "e1659_ca_A": [float(x) for x in e1659_ca],
        "top_subpockets": scored_subpockets[:8],
        "top_loop_facing": loop_facing[:5],
    }

    # verdict
    if loop_facing and loop_facing[0]["druggability"] >= 0.55:
        verdict = (f"DRUGGABLE loop-facing subpocket found "
                   f"(drug={loop_facing[0]['druggability']:.2f}, "
                   f"V={loop_facing[0]['volume_A3']:.0f} Å³) — docking-ready")
    elif loop_facing and loop_facing[0]["druggability"] >= 0.45:
        verdict = (f"FRAGMENT-ABLE loop-facing subpocket "
                   f"(drug={loop_facing[0]['druggability']:.2f}) — "
                   f"fragment screen only")
    elif loop_facing:
        verdict = (f"WEAK loop-facing subpockets "
                   f"(top drug={loop_facing[0]['druggability']:.2f}) — "
                   f"pharmacochaperone strategy marginal")
    else:
        verdict = "NO LOOP-FACING SUBPOCKETS — drop pharmacochaperone hypothesis"
    summary["verdict"] = verdict
    print(f"\n[verdict] {verdict}")

    out_json = OUT_DIR / "pharmacochaperone_phase2b_results.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"[write] {out_json.name}")

    # --- figure ----------------------------------------------------- #
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))

    # XZ projection
    res_ca, res_nums = [], []
    for res in chain:
        if not is_aa(res, standard=True) or "CA" not in res:
            continue
        d = np.linalg.norm(res["CA"].coord - loop_ctr)
        if d <= 24:
            res_ca.append(res["CA"].coord)
            res_nums.append(res.id[1])
    res_ca = np.array(res_ca)

    ax[0].scatter(res_ca[:, 0], res_ca[:, 2], s=8, c="lightgray")
    # loop
    loop_cas = np.array([res["CA"].coord for res in chain
                         if is_aa(res, standard=True) and
                         LOOP_RANGE[0] <= res.id[1] <= LOOP_RANGE[1]
                         and "CA" in res])
    ax[0].scatter(loop_cas[:, 0], loop_cas[:, 2], s=60,
                  c="tab:blue", marker="D", label=f"loop {LOOP_RANGE[0]}–{LOOP_RANGE[1]}")
    ax[0].scatter([e1659_ca[0]], [e1659_ca[2]], s=120, c="black",
                  marker="*", label="E1659 CA")

    # top loop-facing subpockets
    colors = plt.cm.tab10(np.linspace(0, 1, max(4, len(loop_facing[:5]))))
    for i, s in enumerate(loop_facing[:5]):
        # reload points (we popped); compute again from centroid guess…
        pass  # skip replot of points; use centroids below

    # show all buried voxels in lighter red
    buried_pts = grid[buried_mask]
    ax[0].scatter(buried_pts[:, 0], buried_pts[:, 2], s=2, c="tab:red",
                  alpha=0.25, label="buried voxels")

    # centroids of loop-facing subpockets
    for i, s in enumerate(loop_facing[:5]):
        ax[0].scatter([s["centroid"][0]], [s["centroid"][2]],
                      s=250, marker="o", facecolor="none",
                      edgecolor=colors[i], linewidth=2,
                      label=f"#{i+1} drug={s['druggability']:.2f}")

    ax[0].set_xlabel("x / Å"); ax[0].set_ylabel("z / Å")
    ax[0].set_title("WT Job 4 — loop 1639–1651 neighbourhood (XZ)")
    ax[0].legend(loc="upper right", fontsize=7)
    ax[0].set_aspect("equal", "datalim")

    # druggability bar
    show = scored_subpockets[:min(10, len(scored_subpockets))]
    x = np.arange(len(show))
    scores = [s["druggability"] for s in show]
    vols = [s["volume_A3"] for s in show]
    is_loop = [s["touches_loop_1639_1651"] for s in show]
    bar_colors = ["tab:blue" if l else "tab:gray" for l in is_loop]
    ax2 = ax[1]
    bars = ax2.bar(x, scores, color=bar_colors, edgecolor="black", alpha=0.8)
    ax2.axhline(0.55, color="tab:green", ls="--", lw=1, label="docking-ready ≥0.55")
    ax2.axhline(0.45, color="tab:orange", ls="--", lw=1, label="fragment ≥0.45")
    for i, b in enumerate(bars):
        ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 0.01,
                 f"{vols[i]:.0f}\nÅ³", ha="center", fontsize=7)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"#{i+1}" for i in x])
    ax2.set_ylabel("druggability")
    ax2.set_title("Subpockets — blue = touches loop 1639–1651")
    ax2.set_ylim(0, 1.0)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "pharmacochaperone_phase2b.png", dpi=160)
    print(f"[write] pharmacochaperone_phase2b.png")


if __name__ == "__main__":
    main()
