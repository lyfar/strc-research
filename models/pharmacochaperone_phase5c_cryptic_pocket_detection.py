#!/usr/bin/env python3
"""
Phase 5c — K1141 pocket stability + grid-based cavity scan on Phase 5a MD ensemble.

Goal: decide whether the Phase 5b RED-LIGHT (all leads f_PC < 0.10) reflects
  (a) a fundamentally undruggable K1141 pocket → retarget Phase 3c v2
      to an alternative site found in the MD dynamics; or
  (b) a stable K1141 pocket → RED-LIGHT is chemistry-limited → Phase 3c v2
      runs the expanded chem-space screen on K1141.

Approach (fpocket/qhull broken on brew 4.0 build, so custom grid + analysis):
  1. K1141 stability: per-residue RMSF over 20 snapshots; local RMSD of
     pocket-lining residues; K1141 pocket volume time-series via grid void.
  2. Global cavity scan (on snap_010 midpoint): enumerate all cavity
     clusters ≥ 50 Å³, rank by volume, report their distance to K1141.
  3. Verdict: K1141_STABLE vs K1141_COLLAPSED; rank alternative cavities.

Deps: numpy, scipy, OpenMM (for PDB I/O with hybrid36 support).
Runtime: ~2 min.
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Reuse hybrid36-aware stripper from Phase 5b
sys.path.insert(0, "/Users/egorlyfar/Brain/research/strc/models")
from pharmacochaperone_phase5b_ensemble_redock import (  # type: ignore
    strip_water_ions_to_pdb,
)

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
SNAPSHOT_DIR = WORK / "artifacts" / "phase5a_snapshots"
OUT_JSON = WORK / "pharmacochaperone_phase5c_cryptic_pocket_detection.json"

# Phase 4b K1141 pocket centre (target of Phase 3/4/5b virtual screen)
K1141_CENTRE = np.array([-22.027, -18.547, 2.215])
POCKET_RESIDUE_CUTOFF_A = 6.0   # residues within this radius define "pocket lining"

# Grid cavity parameters
GRID_DX_A = 1.5
PAD_A = 4.0
# A grid point is "cavity" if:
#   - closest protein heavy atom is > PROBE_MIN from it (not inside protein),
#   - and < PROBE_MAX (buried, not bulk solvent),
#   - and it is enclosed by protein along ≥ DIRECTIONS_BURIED of 14 cone dirs
PROBE_MIN_A = 2.6
PROBE_MAX_A = 4.5
DIRECTIONS_BURIED = 16  # of 26 probe cones — ~60% burial


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


def read_heavy_atom_coords(pdb_path: Path) -> tuple[np.ndarray, list[tuple[int, str, str]]]:
    """Return coords (N×3) and parallel list of (res_id, res_name, atom_name)
    after H-strip. Uses OpenMM so hybrid36 numbering survives read, but we
    emit sequential integer residue IDs using residue object order."""
    from openmm.app import PDBFile
    pdb = PDBFile(str(pdb_path))
    coords = []
    meta = []
    res_counter = 0
    for chain in pdb.topology.chains():
        if chain.id != "A":
            continue
        for res in chain.residues():
            res_counter += 1
            for atom in res.atoms():
                if atom.element is None:
                    continue
                if atom.element.symbol == "H":
                    continue
                p = pdb.positions[atom.index].value_in_unit_system(
                    __import__("openmm.unit", fromlist=["md_unit_system"]).md_unit_system)
                # md_unit_system returns nm; convert to Å
                coords.append([p[0] * 10, p[1] * 10, p[2] * 10])
                meta.append((res_counter, res.name, atom.name))
    return np.array(coords), meta


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return rotation R and translation t that align P onto Q (both N×3, same N)."""
    pc = P.mean(axis=0)
    qc = Q.mean(axis=0)
    X = P - pc
    Y = Q - qc
    H = X.T @ Y
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = qc - R @ pc
    return R, t


def per_atom_rmsf(frames: list[np.ndarray]) -> np.ndarray:
    """RMSF per atom across pre-aligned frames (list of N×3), returned length-N."""
    F = np.stack(frames, axis=0)  # (F, N, 3)
    mean = F.mean(axis=0)          # (N, 3)
    sq = ((F - mean) ** 2).sum(axis=-1)  # (F, N)
    return np.sqrt(sq.mean(axis=0))      # (N,)


def k1141_pocket_residues(coords: np.ndarray, meta: list) -> set[int]:
    """Residues with at least one heavy atom within POCKET_RESIDUE_CUTOFF_A
    of K1141_CENTRE in a given snapshot."""
    d = np.linalg.norm(coords - K1141_CENTRE, axis=1)
    near = np.where(d <= POCKET_RESIDUE_CUTOFF_A)[0]
    return {meta[i][0] for i in near}


def voxel_cavity_scan(coords: np.ndarray, k1141: np.ndarray) -> dict:
    """Grid-based cavity detection over protein bounding box.
    Returns list of cavity clusters with volume, centre, distance to K1141,
    and per-cluster #points.
    """
    from scipy.spatial import cKDTree
    from scipy.ndimage import label

    lo = coords.min(axis=0) - PAD_A
    hi = coords.max(axis=0) + PAD_A
    nx = int(np.ceil((hi[0] - lo[0]) / GRID_DX_A))
    ny = int(np.ceil((hi[1] - lo[1]) / GRID_DX_A))
    nz = int(np.ceil((hi[2] - lo[2]) / GRID_DX_A))
    log(f"  voxel grid {nx}×{ny}×{nz} = {nx*ny*nz:,} voxels")

    xs = lo[0] + GRID_DX_A * np.arange(nx)
    ys = lo[1] + GRID_DX_A * np.arange(ny)
    zs = lo[2] + GRID_DX_A * np.arange(nz)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    grid_pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    tree = cKDTree(coords)
    d_nearest, _ = tree.query(grid_pts, k=1)
    # Cavity candidates: in a shell PROBE_MIN < d < PROBE_MAX from nearest heavy atom
    candidate = (d_nearest > PROBE_MIN_A) & (d_nearest < PROBE_MAX_A)
    log(f"  candidate voxels (by shell): {candidate.sum():,}")

    # Enclosedness: for each candidate, cast 14 rays in cube+face-diagonal
    # directions; count rays that hit protein within 8 Å. Dense but fast.
    dirs = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == dy == dz == 0:
                    continue
                v = np.array([dx, dy, dz], dtype=float)
                dirs.append(v / np.linalg.norm(v))
    dirs = np.array(dirs)  # 26 directions; we count ≥ DIRECTIONS_BURIED of 26 below

    cand_pts = grid_pts[candidate]
    # Probe each direction by querying KD-tree for atoms within a narrow tube:
    # step along dir and check distance to tree until hit or 8 Å exceeded.
    # To keep it fast: step 7 times at 1.2 Å, mark direction as "buried" if
    # any step finds a heavy atom within 2.0 Å.
    STEP = 1.0
    NSTEPS = 6
    HIT_R = 1.8
    buried_count = np.zeros(len(cand_pts), dtype=int)
    for d_vec in dirs:
        hit = np.zeros(len(cand_pts), dtype=bool)
        for k in range(1, NSTEPS + 1):
            probe_pts = cand_pts + (k * STEP) * d_vec
            nd, _ = tree.query(probe_pts, k=1)
            hit |= (nd < HIT_R) & (~hit)
        buried_count += hit.astype(int)

    enclosed = buried_count >= DIRECTIONS_BURIED
    log(f"  enclosed voxels: {enclosed.sum():,} / {len(cand_pts):,} candidates")

    # Build a 3D boolean cube of enclosed voxels to cluster connected components
    cube = np.zeros((nx, ny, nz), dtype=bool)
    cand_indices = np.where(candidate)[0]
    enclosed_indices = cand_indices[enclosed]
    kx, rem = np.divmod(enclosed_indices, ny * nz)
    ky, kz = np.divmod(rem, nz)
    cube[kx, ky, kz] = True

    labels, n_clust = label(cube)
    log(f"  cavity clusters: {n_clust}")

    clusters = []
    for cid in range(1, n_clust + 1):
        pts_mask = labels == cid
        n_pts = int(pts_mask.sum())
        if n_pts < 20:   # < ~70 Å³, too small for a drug pocket
            continue
        idx = np.argwhere(pts_mask)
        centres_xyz = (lo + idx * GRID_DX_A + GRID_DX_A / 2)
        cent = centres_xyz.mean(axis=0)
        volume_A3 = n_pts * (GRID_DX_A ** 3)
        d_k1141 = float(np.linalg.norm(cent - k1141))
        # Residues lining this cavity: protein residues with any heavy atom
        # within 4 Å of ANY cluster voxel.
        cluster_tree = cKDTree(centres_xyz)
        nearest_prot_d, nearest_prot_i = cluster_tree.query(coords, k=1)
        lining_atom_mask = nearest_prot_d < 4.0
        clusters.append({
            "volume_A3": round(volume_A3, 1),
            "n_voxels": n_pts,
            "centre": cent.tolist(),
            "dist_to_K1141_A": round(d_k1141, 2),
            "lining_atom_indices": np.where(lining_atom_mask)[0].tolist(),
        })
    clusters.sort(key=lambda c: -c["volume_A3"])
    return clusters


def k1141_pocket_volume_per_snapshot(coords_series: list[np.ndarray]) -> list[float]:
    """Quick sphere-void estimate: count grid voxels within 8 Å of K1141_CENTRE
    that have no heavy atom within 2.0 Å. Proportional to local free volume,
    a proxy for pocket openness."""
    from scipy.spatial import cKDTree
    # Precompute a local grid around K1141
    R = 8.0
    lo = K1141_CENTRE - R
    hi = K1141_CENTRE + R
    steps = int(np.ceil((hi[0] - lo[0]) / GRID_DX_A))
    ax = lo[0] + GRID_DX_A * np.arange(steps)
    ay = lo[1] + GRID_DX_A * np.arange(steps)
    az = lo[2] + GRID_DX_A * np.arange(steps)
    X, Y, Z = np.meshgrid(ax, ay, az, indexing="ij")
    gp = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    # Keep grid points within 8 Å of K1141
    within = np.linalg.norm(gp - K1141_CENTRE, axis=1) <= R
    gp = gp[within]

    vols = []
    for coords in coords_series:
        tree = cKDTree(coords)
        d, _ = tree.query(gp, k=1)
        void_mask = d > 2.0
        vols.append(round(void_mask.sum() * (GRID_DX_A ** 3), 1))
    return vols


def main():
    snapshots = sorted(SNAPSHOT_DIR.glob("snap_???.pdb"))
    log(f"=== Phase 5c MD stability + cavity scan ===")
    log(f"Loading {len(snapshots)} snapshots...")

    all_coords: list[np.ndarray] = []
    all_meta: list[list] = []
    for i, snap in enumerate(snapshots):
        with tempfile.TemporaryDirectory() as tmpd:
            prot = Path(tmpd) / f"{snap.stem}_prot.pdb"
            strip_water_ions_to_pdb(snap, prot)
            coords, meta = read_heavy_atom_coords(prot)
        all_coords.append(coords)
        all_meta.append(meta)
        if i == 0:
            log(f"  snap_000: {len(coords)} heavy atoms, "
                f"{len(set(m[0] for m in meta))} residues")

    # Align all snapshots to snap_000 on Cα
    ref_coords = all_coords[0]
    ref_meta = all_meta[0]
    ca_mask = np.array([m[2] == "CA" for m in ref_meta])
    ref_ca = ref_coords[ca_mask]
    log(f"  Cα atoms: {len(ref_ca)}")

    aligned = [ref_coords.copy()]
    for i in range(1, len(all_coords)):
        co = all_coords[i]
        me = all_meta[i]
        ca_mask_i = np.array([m[2] == "CA" for m in me])
        if ca_mask_i.sum() != len(ref_ca):
            log(f"  WARN frame {i} Cα count mismatch — skipping alignment")
            aligned.append(co)
            continue
        R, t = kabsch(co[ca_mask_i], ref_ca)
        aligned.append(co @ R.T + t)

    # Per-residue RMSF (on Cα only)
    ca_frames = [A[ca_mask] for A in aligned]
    rmsf_ca = per_atom_rmsf(ca_frames)
    mean_rmsf = float(rmsf_ca.mean())
    log(f"Global mean Cα RMSF: {mean_rmsf:.2f} Å")

    # K1141 pocket residues in snap_000
    pocket_res = k1141_pocket_residues(ref_coords, ref_meta)
    log(f"K1141 pocket lining residues (static snap_000): {sorted(pocket_res)}")
    pocket_ca_mask = np.array([
        (m[0] in pocket_res and m[2] == "CA") for m in ref_meta
    ])
    pocket_rmsf = rmsf_ca[pocket_ca_mask[ca_mask]] if pocket_ca_mask.any() else np.array([mean_rmsf])
    pocket_mean_rmsf = float(pocket_rmsf.mean()) if len(pocket_rmsf) else mean_rmsf
    log(f"K1141 pocket mean Cα RMSF: {pocket_mean_rmsf:.2f} Å "
        f"(vs global {mean_rmsf:.2f} Å)")

    # K1141 pocket volume time series (voxel void count around K1141)
    k1141_vol_series = k1141_pocket_volume_per_snapshot(aligned)
    log(f"K1141 local void volume series (Å³): "
        f"mean={np.mean(k1141_vol_series):.1f} "
        f"min={np.min(k1141_vol_series):.1f} "
        f"max={np.max(k1141_vol_series):.1f}")

    # Global cavity scan on mid-trajectory snapshot
    mid_i = len(aligned) // 2
    log(f"=== Global cavity scan on snap_{mid_i:03d} ===")
    cavities = voxel_cavity_scan(aligned[mid_i], K1141_CENTRE)
    log(f"Top 5 cavities by volume:")
    log(f"{'#':>3} {'vol_A3':>8} {'d_K1141':>8} {'centre':>30}")
    for i, c in enumerate(cavities[:5]):
        log(f"{i:>3} {c['volume_A3']:>8.1f} {c['dist_to_K1141_A']:>8.2f} "
            f"[{c['centre'][0]:>8.2f}, {c['centre'][1]:>8.2f}, {c['centre'][2]:>8.2f}]")

    # Find cavity containing K1141 (smallest dist_to_K1141_A)
    k1141_cavity = cavities[0] if cavities else None
    for c in cavities:
        if c["dist_to_K1141_A"] <= POCKET_RESIDUE_CUTOFF_A:
            k1141_cavity = c
            break

    # Residues lining K1141's cavity
    k1141_residues = set()
    if k1141_cavity:
        for atom_idx in k1141_cavity["lining_atom_indices"]:
            r = all_meta[mid_i][atom_idx][0]
            k1141_residues.add(r)

    # Alternative cavities: rank those NOT overlapping K1141 by volume
    alt_cavities = [c for c in cavities
                    if c["dist_to_K1141_A"] > POCKET_RESIDUE_CUTOFF_A][:5]

    # Verdict
    log("=== VERDICT ===")
    pocket_stable = pocket_mean_rmsf < 1.5 * mean_rmsf and min(k1141_vol_series) > 20
    if pocket_stable:
        log(f"K1141 pocket STABLE (RMSF {pocket_mean_rmsf:.2f} Å within "
            f"{1.5*mean_rmsf:.2f} Å bound; min void {min(k1141_vol_series):.0f} Å³)")
        log("→ Phase 5b RED-LIGHT is chemistry-limited; Phase 3c v2 should "
            "screen expanded chem-space against K1141 with ensemble targeting.")
    else:
        log(f"K1141 pocket UNSTABLE (RMSF {pocket_mean_rmsf:.2f} Å vs global "
            f"{mean_rmsf:.2f} Å; min void {min(k1141_vol_series):.0f} Å³)")
        log("→ Consider Phase 3c v2 retargeting to a conformer-specific "
            "K1141 pose (ensemble docking already caught this in 5b).")

    log(f"Alternative cavities (dist > {POCKET_RESIDUE_CUTOFF_A} Å from K1141): "
        f"{len(alt_cavities)} found, top vol {alt_cavities[0]['volume_A3'] if alt_cavities else 0:.1f} Å³")
    if alt_cavities and alt_cavities[0]["volume_A3"] > (k1141_cavity["volume_A3"] if k1141_cavity else 0) * 1.5:
        log(f"  → LARGE ALTERNATIVE CAVITY at {alt_cavities[0]['centre']} — "
            f"worth exploring in Phase 3c v2 parallel screen.")

    # Trim heavy lining_atom_indices from output for size
    for c in cavities:
        c["n_lining_atoms"] = len(c["lining_atom_indices"])
        del c["lining_atom_indices"]

    payload = {
        "phase": "5c",
        "approach": "grid-based cavity scan + K1141 stability (fpocket/qhull broken on brew 4.0)",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_snapshots": len(snapshots),
        "global_mean_rmsf_A": round(mean_rmsf, 3),
        "k1141_pocket_residues_static": sorted(pocket_res),
        "k1141_pocket_mean_rmsf_A": round(pocket_mean_rmsf, 3),
        "k1141_local_void_volume_A3_series": k1141_vol_series,
        "k1141_local_void_volume_mean_A3": round(float(np.mean(k1141_vol_series)), 1),
        "k1141_local_void_volume_min_A3": round(float(np.min(k1141_vol_series)), 1),
        "k1141_cavity_mid_frame": k1141_cavity,
        "top_cavities": cavities[:10],
        "alt_cavities_not_overlap_K1141": alt_cavities,
        "verdict": {
            "k1141_pocket_stable": bool(pocket_stable),
            "rmsf_ratio_pocket_to_global": round(pocket_mean_rmsf / mean_rmsf, 3),
            "has_larger_alt_cavity": bool(
                alt_cavities and k1141_cavity
                and alt_cavities[0]["volume_A3"] > k1141_cavity["volume_A3"] * 1.5
            ),
        },
        "parameters": {
            "grid_dx_A": GRID_DX_A,
            "probe_min_A": PROBE_MIN_A,
            "probe_max_A": PROBE_MAX_A,
            "directions_buried_threshold": DIRECTIONS_BURIED,
            "pocket_residue_cutoff_A": POCKET_RESIDUE_CUTOFF_A,
        },
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str))
    log(f"Written {OUT_JSON}")


if __name__ == "__main__":
    main()
