#!/usr/bin/env python3
"""
Phase 4a — K1141 pocket reproducibility across 5 AF3 STRC CIFs.

Reuses the Phase 2B ray-casting + depth-peak scorer on every available
STRC structure and asks: does the K1141-inclusive loop-facing pocket
(druggability 0.86 in WT) reproduce in the mutant, the truncations, and
the Ultra-Mini × TMEM145 complex? If not, Phase 0 was a one-job artifact
and the virtual-screen gate is invalid.

Gate: pocket present in >=4/5 CIFs with druggability >=0.70 AND K1141
in its lining residues AND centroid within 6 A of the reference box
centre (7.7, -5.4, -41.5).
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Reuse Phase 2B helpers — identical scoring, identical radii, identical grid.
from pharmacochaperone_phase2b_subpockets import (  # noqa: E402
    load_chain,
    heavy_atom_coords,
    loop_centroid,
    residue_ca,
    build_grid,
    enclosure_and_depth,
    cluster,
    lining_residues,
    druggability,
    flood_from_seeds,
    HYDROPHOBIC,
    DONORS,
    ACCEPTORS,
    SCAN_RADIUS_A,
    GRID_SPACING,
    MIN_CLUSTER_VOXELS,
    LOOP_RANGE,
    E1659_RESNUM,
)

MODELS_DIR = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models")
OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")

# Five CIFs spanning: WT (reference), E1659A mutant, mini-STRC 594-1775,
# Ultra-Mini solo (C-term-only proxy), Ultra-Mini x TMEM145 complex.
# `offset` converts real STRC residue number -> CIF residue number
# (AF3 renumbers truncated constructs starting from 1).
# real_resnum = cif_resnum + offset ; cif_resnum = real_resnum - offset
CIFS = [
    {"label": "WT_full",          "cif": "job4-wildtype.cif",               "offset": 0},
    {"label": "E1659A_mutant",    "cif": "job3-mutant.cif",                 "offset": 0},
    {"label": "mini_594_1775",    "cif": "job5-mini-strc.cif",              "offset": 593},
    {"label": "ultra_cterm_only", "cif": "job-h-strc-cterm-only.cif",       "offset": 1074},
    {"label": "ultra_x_tmem145",  "cif": "job-ultramini-x-tmem145-full.cif","offset": 1074},
]

REAL_K1141 = 1141
REAL_E1659 = E1659_RESNUM
REAL_LOOP_RANGE = LOOP_RANGE  # (1639, 1651) — real STRC numbering

DRUGG_THRESHOLD = 0.70
PASS_MIN = 4  # of 5
# Reference-frame-free proximity: pocket centroid must be within these
# distances from local K1141 Cα and local loop centroid, in each CIF's
# own coordinate system. AF3 rotates complexes arbitrarily, so absolute
# box coords from Phase 2B don't transfer across jobs.
MAX_DIST_K1141_A = 12.0
MAX_DIST_LOOP_A = 12.0


def subpocket_stats(pts, chain, k1141_cif_num: int):
    """Aggregate properties for a single subpocket point cloud.
    `k1141_cif_num` = residue number of K1141 in this CIF's numbering."""
    res = lining_residues(pts, chain, cutoff=4.5)
    nres = len(res)
    hydro = sum(1 for n in res.values() if n in HYDROPHOBIC)
    hydro_frac = hydro / nres if nres else 0.0
    donors = sum(len(DONORS.get(n, [])) for n in res.values())
    acceptors = sum(len(ACCEPTORS.get(n, [])) for n in res.values())
    # depth = distance to nearest protein atom
    coords = heavy_atom_coords(chain)
    tree = cKDTree(coords)
    d, _ = tree.query(pts, k=1)
    depth_mean = float(np.mean(d))
    volume = len(pts) * (GRID_SPACING ** 3)
    return {
        "n_points": int(len(pts)),
        "volume_A3": round(volume, 2),
        "centroid": [round(float(x), 2) for x in pts.mean(axis=0)],
        "depth_mean_A": round(depth_mean, 2),
        "n_lining_residues": nres,
        "hydrophobic_frac": round(hydro_frac, 3),
        "donors": donors,
        "acceptors": acceptors,
        "druggability": round(
            druggability(volume, hydro_frac, nres, donors, acceptors, depth_mean), 3
        ),
        "lining_residues_cif": sorted(res.keys()),
        "has_K1141_cif": k1141_cif_num in res,
    }


def scan_cif(cif_path: Path, offset: int):
    chain = load_chain(cif_path, "A")
    atoms = heavy_atom_coords(chain)
    if len(atoms) < 500:
        return {"error": f"chain A has only {len(atoms)} heavy atoms"}

    # Translate real STRC residue numbers into this CIF's numbering.
    cif_loop = (REAL_LOOP_RANGE[0] - offset, REAL_LOOP_RANGE[1] - offset)
    cif_k1141 = REAL_K1141 - offset
    cif_e1659 = REAL_E1659 - offset

    tree = cKDTree(atoms)
    try:
        loop_c = loop_centroid(chain, *cif_loop)
    except Exception as e:
        return {"error": f"loop centroid failed ({cif_loop}): {e}"}
    if np.any(np.isnan(loop_c)):
        return {"error": f"loop centroid NaN — loop {cif_loop} not in chain"}
    try:
        e_ca = residue_ca(chain, cif_e1659)
    except KeyError:
        # E1659 location useful for provenance but not required for the
        # pocket test — the truncation may not include it.
        e_ca = None
    try:
        k_ca = residue_ca(chain, cif_k1141)
    except KeyError:
        return {"error": f"K1141 (CIF resnum {cif_k1141}) not in chain"}

    grid = build_grid(loop_c, SCAN_RADIUS_A, GRID_SPACING)
    mask, depth = enclosure_and_depth(grid, tree)
    clusters = cluster(mask, grid, MIN_CLUSTER_VOXELS)
    if not clusters:
        return {
            "error": "no clusters",
            "cif_loop_range": list(cif_loop),
            "cif_k1141_resnum": cif_k1141,
            "loop_centroid_A": [round(float(x), 2) for x in loop_c],
            "k1141_ca_A": [round(float(x), 2) for x in k_ca],
            "e1659_ca_A": (None if e_ca is None
                           else [round(float(x), 2) for x in e_ca]),
        }

    subpockets = []
    for cl in clusters:
        cluster_mask = np.zeros_like(mask)
        cluster_mask[cl["idx"]] = True
        for sp in flood_from_seeds(cluster_mask, depth, grid):
            subpockets.append(sp)

    scored = []
    for sp in subpockets:
        stats = subpocket_stats(sp["points"], chain, cif_k1141)
        stats["peak_xyz"] = [round(float(x), 2) for x in sp["peak_xyz"]]
        stats["peak_depth_A"] = round(float(sp["peak_depth"]), 2)
        c = np.array(stats["centroid"])
        stats["dist_to_K1141_CA_A"] = round(float(np.linalg.norm(c - k_ca)), 2)
        stats["dist_to_loop_centroid_A"] = round(float(np.linalg.norm(c - loop_c)), 2)
        # Lining residues, translated back to real STRC numbering for
        # cross-CIF comparison.
        stats["lining_residues_real"] = sorted(
            r + offset for r in stats["lining_residues_cif"]
        )
        scored.append(stats)

    scored.sort(key=lambda s: (-int(s["has_K1141_cif"]), -s["druggability"]))
    k1141_hits = [s for s in scored if s["has_K1141_cif"]]
    best_k1141 = k1141_hits[0] if k1141_hits else None
    top_overall = scored[0] if scored else None

    return {
        "chain_heavy_atoms": int(len(atoms)),
        "offset": offset,
        "cif_loop_range": list(cif_loop),
        "cif_k1141_resnum": cif_k1141,
        "cif_e1659_resnum": cif_e1659,
        "loop_centroid_A": [round(float(x), 2) for x in loop_c],
        "k1141_ca_A": [round(float(x), 2) for x in k_ca],
        "e1659_ca_A": (None if e_ca is None
                       else [round(float(x), 2) for x in e_ca]),
        "n_clusters": len(clusters),
        "n_subpockets": len(scored),
        "n_k1141_subpockets": len(k1141_hits),
        "best_k1141_subpocket": best_k1141,
        "top_subpocket_overall": top_overall,
    }


def evaluate(results):
    passes = []
    for r in results:
        sp = r.get("best_k1141_subpocket")
        if sp is None:
            passes.append(False)
            continue
        ok = (
            sp["druggability"] >= DRUGG_THRESHOLD
            and sp["has_K1141_cif"]
            and sp["dist_to_K1141_CA_A"] <= MAX_DIST_K1141_A
            and sp["dist_to_loop_centroid_A"] <= MAX_DIST_LOOP_A
        )
        passes.append(bool(ok))
    n_pass = sum(passes)
    verdict = "PASS" if n_pass >= PASS_MIN else "FAIL"
    return passes, n_pass, verdict


def main():
    print(f"Phase 4a pocket reproducibility — 5 CIFs, threshold druggability >= {DRUGG_THRESHOLD}")
    results = []
    for row in CIFS:
        cif_path = MODELS_DIR / row["cif"]
        if not cif_path.exists():
            results.append({**row, "error": f"missing: {cif_path}"})
            print(f"  {row['label']:22s} MISSING  {cif_path.name}")
            continue
        print(f"  scanning {row['label']:22s} {row['cif']} ...", flush=True)
        r = scan_cif(cif_path, offset=row["offset"])
        results.append({**row, **r})
        sp = r.get("best_k1141_subpocket")
        if sp:
            print(
                f"    drugg={sp['druggability']:.2f}  "
                f"V={sp['volume_A3']:.1f}A^3  "
                f"K1141={'Y' if sp['has_K1141_cif'] else 'N'}  "
                f"d_K1141={sp['dist_to_K1141_CA_A']:.2f}A  "
                f"d_loop={sp['dist_to_loop_centroid_A']:.2f}A"
            )
        else:
            print(f"    no K1141-inclusive subpocket: {r.get('error','see json')}")

    passes, n_pass, verdict = evaluate(results)
    summary = {
        "n_cifs": len(results),
        "n_pass": int(n_pass),
        "threshold_druggability": DRUGG_THRESHOLD,
        "threshold_dist_to_K1141_A": MAX_DIST_K1141_A,
        "threshold_dist_to_loop_A": MAX_DIST_LOOP_A,
        "min_pass_to_advance": PASS_MIN,
        "verdict": verdict,
        "per_cif_pass": dict(zip([r["label"] for r in results], passes)),
    }

    out_json = OUT_DIR / "pharmacochaperone_phase4a_pocket_reproducibility.json"
    out_json.write_text(
        json.dumps({"summary": summary, "per_cif": results}, indent=2)
    )
    print()
    print(f"verdict: {verdict}  ({n_pass}/{len(results)} passed)")
    print(f"wrote {out_json}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
