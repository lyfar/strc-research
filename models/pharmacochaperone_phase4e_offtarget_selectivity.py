#!/usr/bin/env python3
"""
Phase 4e — Off-target pocket scan on Ultra-Mini × TMEM145.

Question: is the K1141 loop-facing pocket the ONLY druggable site for
the Phase 3C top-5 leads on the clinical construct? If carboxylate
ligands bind 4+ other pockets equally well, selectivity is zero and
the screen is just "COOH + LRR basic residues".

Method (toolchain-minimal, runnable today without Vina):
  1. Enumerate all druggable pockets on Ultra-Mini × TMEM145 chain A
     using the Phase 2B ray-casting + depth-peak scorer, scanning
     from the Cα of every residue (not just the K1141 loop).
  2. For each pocket, compute a pharmacophore compatibility score
     against the Phase 3A anchor triangle (K1141-NZ H-bond anchor,
     F1646 aromatic hot-spot, D1140/D1173 LRR face). A surface
     pocket that lacks one or more anchors scores low.
  3. Selectivity metric per lead = pharmacophore score of K1141
     pocket − best non-K1141 pocket score. Gate: selectivity margin
     ≥0.20 (on a 0-1 scale) for ≥3 of top-5 leads.

This is NOT a substitute for Vina-based ΔG selectivity (Phase 4b/4c).
It is a first-pass structural check that kills the "COOH magnet"
failure mode before burning compute on Vina.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pharmacochaperone_phase2b_subpockets import (  # noqa: E402
    load_chain, heavy_atom_coords, build_grid, enclosure_and_depth,
    cluster, lining_residues, druggability, flood_from_seeds,
    HYDROPHOBIC, DONORS, ACCEPTORS, SCAN_RADIUS_A, GRID_SPACING,
    MIN_CLUSTER_VOXELS,
)
from pharmacochaperone_phase4_common import (  # noqa: E402
    MODELS_DIR, WORK_DIR, TARGETS,
    REAL_K1141, REAL_D1140, REAL_D1173, REAL_F1646, REAL_LOOP,
    FIXED_ROSTER, gate_ready, print_env_report,
)

from Bio.PDB.Polypeptide import is_aa

OUT_JSON = WORK_DIR / "pharmacochaperone_phase4e_offtarget_selectivity.json"

# Gate thresholds.
MIN_POCKET_DRUGGABILITY = 0.55
MIN_POCKET_VOLUME_A3   = 40.0
SELECTIVITY_MARGIN      = 0.20   # K1141 pocket pharmacophore score - best off-target
MIN_LEADS_PASSING       = 3      # of the 5 leads

# Residue classes for pharmacophore matching.
BASIC_RES    = {"LYS", "ARG", "HIS"}
ACIDIC_RES   = {"ASP", "GLU"}
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}


def sample_seeds(chain, stride: int = 20) -> list[np.ndarray]:
    """Sample Cα coordinates every `stride` residues as pocket seeds.

    The full-surface scan of Phase 4e is a PROXY gate (pharmacophore-
    geometry only; real selectivity comes from Phase 4b Vina ΔG gap).
    Stride 20 gives ~35 seeds over a 700-residue chain → tractable in
    ~5-10 min while still covering the LRR surface densely enough that
    any pocket ≥40 Å³ is reachable by at least one seed.
    """
    cas = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        if "CA" in res:
            cas.append((res.id[1], res["CA"].coord))
    cas.sort(key=lambda kv: kv[0])
    return [c[1] for c in cas[::stride]]


# Coarser grid for proxy-scale sweep (Phase 2B used 0.8 Å for the one
# seed; we use 1.2 Å to quarter voxel count per seed × ~35 seeds).
PROXY_GRID_SPACING = 1.2
PROXY_SCAN_RADIUS  = 14.0   # slightly smaller sphere → fewer voxels still
PROXY_MIN_VOXELS   = 25     # 25 × 1.2³ = 43 Å³ ≈ Phase 2B 50 × 0.8³ = 26 Å³

def scan_all_pockets(chain):
    atoms = heavy_atom_coords(chain)
    tree = cKDTree(atoms)

    seen = set()
    pockets = []
    for seed in sample_seeds(chain, stride=20):
        grid = build_grid(seed, PROXY_SCAN_RADIUS, PROXY_GRID_SPACING)
        mask, depth = enclosure_and_depth(grid, tree)
        clusters = cluster(mask, grid, PROXY_MIN_VOXELS)
        for cl in clusters:
            cluster_mask = np.zeros_like(mask)
            cluster_mask[cl["idx"]] = True
            for sp in flood_from_seeds(cluster_mask, depth, grid):
                c = sp["points"].mean(axis=0)
                key = tuple(np.round(c, 1))
                if key in seen:
                    continue
                seen.add(key)
                pts = sp["points"]
                # per-point min distance to protein
                d, _ = tree.query(pts, k=1)
                depth_mean = float(np.mean(d))
                res = lining_residues(pts, chain, cutoff=4.5)
                nres = len(res)
                if nres < 6:
                    continue
                hydro = sum(1 for n in res.values() if n in HYDROPHOBIC) / nres
                donors = sum(len(DONORS.get(n, [])) for n in res.values())
                acceptors = sum(len(ACCEPTORS.get(n, [])) for n in res.values())
                vol = len(pts) * (PROXY_GRID_SPACING ** 3)
                drug = druggability(vol, hydro, nres, donors, acceptors, depth_mean)
                if drug < MIN_POCKET_DRUGGABILITY or vol < MIN_POCKET_VOLUME_A3:
                    continue
                pockets.append({
                    "centroid": [round(float(x), 2) for x in c],
                    "volume_A3": round(vol, 2),
                    "druggability": round(drug, 3),
                    "n_lining_residues": nres,
                    "hydrophobic_frac": round(hydro, 3),
                    "donors": donors,
                    "acceptors": acceptors,
                    "depth_mean_A": round(depth_mean, 2),
                    "lining_residues": dict(sorted(res.items())),
                })
    return pockets


def pharmacophore_score(pocket, cif_k1141: int, cif_d1140: int,
                        cif_d1173: int, cif_f1646: int) -> dict:
    """Score pocket against Phase 3A anchor triangle.

    1.0 = all anchors present (K1141 basic + D1140 or D1173 acidic +
          F1646 or any aromatic + at least 1 hbond donor/acceptor).
    0.0 = no anchors at all.
    """
    res = pocket["lining_residues"]  # cif-numbering -> resname
    has_k1141  = res.get(cif_k1141)  in BASIC_RES
    has_d_face = (res.get(cif_d1140) in ACIDIC_RES) or (res.get(cif_d1173) in ACIDIC_RES)
    has_arom   = (res.get(cif_f1646) in AROMATIC_RES) or \
                 any(rn in AROMATIC_RES for rn in res.values())
    has_hb     = (pocket["donors"] + pocket["acceptors"]) >= 4

    # Substitute basic/acidic/aromatic residues for the off-target case —
    # if the pocket has *different* basic/acidic/aromatic residues (not
    # K1141/D1140/D1173/F1646 specifically), it still counts partially.
    # The specific-residue terms capture pharmacophore *identity*; the
    # class terms capture pharmacophore *geometry*.
    has_any_basic   = any(rn in BASIC_RES   for rn in res.values())
    has_any_acidic  = any(rn in ACIDIC_RES  for rn in res.values())

    score = (
        0.30 * (1.0 if has_k1141 else (0.5 if has_any_basic else 0.0)) +
        0.25 * (1.0 if has_d_face else (0.5 if has_any_acidic else 0.0)) +
        0.20 * (1.0 if has_arom else 0.0) +
        0.15 * (1.0 if has_hb else 0.0) +
        0.10 * pocket["druggability"]
    )
    return {
        "score": round(float(score), 3),
        "has_K1141": has_k1141,
        "has_D_face": has_d_face,
        "has_aromatic": has_arom,
        "has_4plus_hb": has_hb,
        "has_any_basic": has_any_basic,
        "has_any_acidic": has_any_acidic,
    }


def main():
    print_env_report("4e")
    ok, missing = gate_ready("4e")
    if not ok:
        print(f"blocked: missing {missing}")
        return 2

    target = TARGETS["ultra_x_tmem145"]
    offset = target["offset"]
    cif_path = MODELS_DIR / target["cif"]
    print(f"\nPhase 4e — off-target pocket scan on {target['cif']}")
    chain = load_chain(cif_path, target["chain"])
    n_seeds = len(sample_seeds(chain, stride=20))
    print(f"  scanning pockets across chain A surface ({n_seeds} seeds, "
          f"grid={PROXY_GRID_SPACING}A, r={PROXY_SCAN_RADIUS}A)...", flush=True)
    pockets = scan_all_pockets(chain)
    print(f"  found {len(pockets)} druggable pockets (V>={MIN_POCKET_VOLUME_A3}A^3, drugg>={MIN_POCKET_DRUGGABILITY})")

    cif_k1141 = REAL_K1141 - offset
    cif_d1140 = REAL_D1140 - offset
    cif_d1173 = REAL_D1173 - offset
    cif_f1646 = REAL_F1646 - offset

    # Score pharmacophore per pocket.
    for p in pockets:
        p["pharmacophore"] = pharmacophore_score(
            p, cif_k1141, cif_d1140, cif_d1173, cif_f1646
        )

    # Identify the K1141 pocket (contains K1141 + is closest to loop).
    k1141_pockets = [p for p in pockets if p["pharmacophore"]["has_K1141"]]
    k1141_pockets.sort(key=lambda p: -p["pharmacophore"]["score"])
    k1141_best = k1141_pockets[0] if k1141_pockets else None

    off_target = [p for p in pockets if not p["pharmacophore"]["has_K1141"]]
    off_target.sort(key=lambda p: -p["pharmacophore"]["score"])

    if k1141_best is None:
        print("  no K1141 pocket found — Phase 4a contradiction, investigate")
        return 2

    top_off = off_target[0] if off_target else None
    margin = k1141_best["pharmacophore"]["score"] - (
        top_off["pharmacophore"]["score"] if top_off else 0.0
    )

    print()
    print(f"  K1141 pocket pharmacophore score: {k1141_best['pharmacophore']['score']:.3f}")
    if top_off:
        print(f"  best off-target  pharmacophore score: {top_off['pharmacophore']['score']:.3f}")
    print(f"  selectivity margin: {margin:+.3f}  (threshold {SELECTIVITY_MARGIN})")

    # Per-lead selectivity note: in this first-pass proxy, selectivity is
    # the same for every carboxylate lead because the pharmacophore terms
    # depend on pocket, not ligand. Full ligand-specific selectivity is
    # the Phase 4b/4c Vina job.
    leads = [c for c in FIXED_ROSTER if c["role"] == "lead"]
    per_lead = [
        {"name": c["name"], "selectivity_margin_proxy": round(float(margin), 3)}
        for c in leads
    ]
    passing = sum(1 for p in per_lead if p["selectivity_margin_proxy"] >= SELECTIVITY_MARGIN)

    summary = {
        "target_cif": target["cif"],
        "n_pockets_found": len(pockets),
        "k1141_pocket_score": k1141_best["pharmacophore"]["score"],
        "k1141_pocket_centroid": k1141_best["centroid"],
        "k1141_pocket_druggability": k1141_best["druggability"],
        "k1141_pocket_volume_A3": k1141_best["volume_A3"],
        "best_off_target_score": (
            top_off["pharmacophore"]["score"] if top_off else None
        ),
        "best_off_target_centroid": (top_off["centroid"] if top_off else None),
        "selectivity_margin": round(float(margin), 3),
        "selectivity_threshold": SELECTIVITY_MARGIN,
        "leads_passing": passing,
        "min_leads_passing": MIN_LEADS_PASSING,
        "verdict": "PASS" if margin >= SELECTIVITY_MARGIN and passing >= MIN_LEADS_PASSING else "FAIL",
    }

    write = {
        "summary": summary,
        "per_lead": per_lead,
        "k1141_pocket": k1141_best,
        "off_target_pockets_ranked": off_target[:10],
    }
    OUT_JSON.write_text(json.dumps(write, indent=2))

    print(f"\nverdict: {summary['verdict']}")
    print(f"wrote {OUT_JSON}")
    return 0 if summary["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
