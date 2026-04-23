#!/usr/bin/env python3
"""
Phase 5c — Cryptic pocket detection on Phase 5a MD ensemble.

Rationale: Phase 4b/5b targeted the K1141 pocket identified on a single static
Ultra-Mini × TMEM145 conformation. Phase 5b RED-LIGHT (all leads f_PC < 0.10)
may reflect either (a) K1141 pocket is intrinsically undruggable, (b) pocket
K1141 opens differently in solution dynamics than in the static model, or
(c) a better cryptic pocket exists elsewhere that only opens during MD.

Run fpocket on each of 20 Phase 5a snapshots + starting apo structure, cluster
detected pockets spatially, compute per-cluster persistence + mean drug_score,
and flag:

  • K1141 persistence: does the Phase 4b target pocket survive the MD?
  • Cryptic pockets: clusters with persistence ≥ 10/20 AND NOT present in
    the static starting structure (distance > 10 Å from any static pocket
    or drug_score > 2× the best overlapping static pocket).
  • Druggability ranking: top pockets by mean drug_score across ensemble.

If a druggable cryptic pocket emerges → Phase 3c v2 virtual screen retargeted
against that pocket instead of K1141. If no cryptic pocket, K1141 confirmed
as the right target and red-light is chemistry-limited, not site-limited.

Dependencies: fpocket 4.0 (Homebrew), OpenMM (for hybrid36 PDB stripping).

Runtime: 20 snapshots × ~30 s fpocket + aggregation ~ 15 min wall.
"""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path

# Reuse the OpenMM-based hybrid36 stripper from Phase 5b
sys.path.insert(0, "/Users/egorlyfar/Brain/research/strc/models")
from pharmacochaperone_phase5b_ensemble_redock import (  # type: ignore
    strip_water_ions_to_pdb,
)

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
SNAPSHOT_DIR = WORK / "artifacts" / "phase5a_snapshots"
STATIC_APO = SNAPSHOT_DIR / "ultramini_chainA_prepped.pdb"
OUT_DIR = WORK / "artifacts" / "phase5c_fpocket"
OUT_JSON = WORK / "pharmacochaperone_phase5c_cryptic_pocket_detection.json"

# Phase 4b K1141 pocket centre (target of existing screen)
K1141_CENTRE = (-22.027, -18.547, 2.215)

# Clustering parameters
CLUSTER_RADIUS_A = 6.0           # Angstroms centre-to-centre for same pocket
STATIC_DISTANCE_A = 10.0         # > this distance from static = cryptic candidate
PERSISTENCE_THRESHOLD = 10       # ≥ 10/20 snapshots for "persistent"
DRUG_SCORE_CRYPTIC_RATIO = 2.0   # cryptic if ≥ 2× best static drug_score


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


@dataclass
class Pocket:
    snap_id: str
    rank: int
    score: float
    drug_score: float
    volume: float
    n_alpha_spheres: int
    centre: tuple[float, float, float]
    residue_ids: list[int]      # residue numbers lining the pocket


def strip_hydrogens(pdb_in: Path, pdb_out: Path) -> int:
    """fpocket qhull fails on H-heavy PDBs; keep only heavy atoms."""
    n = 0
    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                # Element in cols 77-78, atom name cols 13-16
                elem = line[76:78].strip()
                aname = line[12:16].strip()
                is_h = (elem == "H") or (elem == "" and aname.startswith("H"))
                if is_h:
                    continue
                fout.write(line)
                n += 1
            elif line.startswith(("TER", "END", "CRYST1", "REMARK")):
                fout.write(line)
    return n


def run_fpocket(pdb: Path, work_dir: Path) -> list[Pocket]:
    """Run fpocket, parse info.txt + alpha sphere centres, return list of
    Pocket objects tagged with snap_id from pdb.stem."""
    snap_id = pdb.stem
    # fpocket writes into {input_dir}/{input_stem}_out/ by default; copy PDB
    # into work_dir so we don't pollute snapshot_dir. Also strip hydrogens —
    # fpocket's qhull Delaunay triangulation fails on H-dense all-atom PDBs.
    # Use _nh suffix so we never read-and-write the same path.
    local_pdb = work_dir / f"{pdb.stem}_nh.pdb"
    n_heavy = strip_hydrogens(pdb, local_pdb)
    if n_heavy < 100:
        log(f"  too few heavy atoms ({n_heavy}) — fpocket skip")
        return []
    r = subprocess.run(
        ["fpocket", "-f", str(local_pdb)],
        capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        log(f"  fpocket FAILED: {r.stderr[:400]}")
        return []
    out_dir = work_dir / f"{local_pdb.stem}_out"
    info_txt = out_dir / f"{local_pdb.stem}_info.txt"
    if not info_txt.exists():
        log(f"  no info.txt at {info_txt}")
        return []

    # Parse info.txt — one block per pocket
    text = info_txt.read_text()
    blocks = re.split(r"Pocket\s+(\d+)\s*:\s*", text)
    pockets: list[Pocket] = []
    # re.split with capture gives: ['', '1', 'body1', '2', 'body2', ...]
    for i in range(1, len(blocks), 2):
        try:
            rank = int(blocks[i])
        except ValueError:
            continue
        body = blocks[i + 1]
        score = _parse_float(body, r"Score\s*:\s*(-?\d+\.\d+)")
        drug = _parse_float(body, r"Druggability Score\s*:\s*(-?\d+\.\d+)")
        vol = _parse_float(body, r"Volume\s*:\s*(-?\d+\.\d+)")
        n_alpha = _parse_int(body, r"Number of Alpha Spheres\s*:\s*(\d+)")

        # Parse alpha sphere centres from pocket{rank}_vert.pqr
        vert_pqr = out_dir / "pockets" / f"pocket{rank}_vert.pqr"
        centre, _n_spheres = _parse_alpha_centre(vert_pqr)

        # Residues lining the pocket from pocket{rank}_atm.pdb
        atm_pdb = out_dir / "pockets" / f"pocket{rank}_atm.pdb"
        residues = _parse_residues(atm_pdb)

        if centre is None:
            continue
        pockets.append(Pocket(
            snap_id=snap_id,
            rank=rank,
            score=score if score is not None else 0.0,
            drug_score=drug if drug is not None else 0.0,
            volume=vol if vol is not None else 0.0,
            n_alpha_spheres=n_alpha if n_alpha is not None else 0,
            centre=centre,
            residue_ids=residues,
        ))
    return pockets


def _parse_float(text: str, pattern: str) -> float | None:
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def _parse_int(text: str, pattern: str) -> int | None:
    m = re.search(pattern, text)
    return int(m.group(1)) if m else None


def _parse_alpha_centre(vert_pqr: Path) -> tuple[tuple[float, float, float] | None, int]:
    """Mean xyz across all ATOM/HETATM lines of pocket vertex PQR."""
    if not vert_pqr.exists():
        return None, 0
    xs, ys, zs = [], [], []
    for line in vert_pqr.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                xs.append(x); ys.append(y); zs.append(z)
            except ValueError:
                continue
    if not xs:
        return None, 0
    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)), len(xs)


def _parse_residues(atm_pdb: Path) -> list[int]:
    if not atm_pdb.exists():
        return []
    seen = set()
    for line in atm_pdb.read_text().splitlines():
        if line.startswith("ATOM"):
            try:
                r = int(line[22:26])
                seen.add(r)
            except ValueError:
                continue
    return sorted(seen)


def dist(a: tuple, b: tuple) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def cluster_pockets(all_pockets: list[Pocket]) -> list[dict]:
    """Spatial clustering by centre proximity.
    One pass greedy: each new pocket joins the nearest existing cluster if
    centre distance ≤ CLUSTER_RADIUS_A; otherwise spawns a new cluster."""
    clusters: list[list[Pocket]] = []
    for p in all_pockets:
        # Find nearest cluster centroid
        best_idx, best_d = -1, float("inf")
        for idx, c in enumerate(clusters):
            centroid = tuple(
                sum(q.centre[ax] for q in c) / len(c) for ax in range(3)
            )
            d = dist(p.centre, centroid)
            if d < best_d:
                best_d, best_idx = d, idx
        if best_d <= CLUSTER_RADIUS_A:
            clusters[best_idx].append(p)
        else:
            clusters.append([p])

    # Summarise each cluster
    summaries = []
    for i, c in enumerate(clusters):
        centroid = tuple(
            sum(p.centre[ax] for p in c) / len(c) for ax in range(3)
        )
        mean_drug = sum(p.drug_score for p in c) / len(c)
        mean_score = sum(p.score for p in c) / len(c)
        mean_vol = sum(p.volume for p in c) / len(c)
        persistence = len({p.snap_id for p in c})
        # Consensus residues: those appearing in ≥ 50% of cluster members
        residue_counts: dict[int, int] = {}
        for p in c:
            for r in p.residue_ids:
                residue_counts[r] = residue_counts.get(r, 0) + 1
        consensus = sorted([
            r for r, n in residue_counts.items() if n >= len(c) * 0.5
        ])
        summaries.append({
            "cluster_id": i,
            "centroid": centroid,
            "persistence": persistence,
            "n_total_members": len(c),
            "mean_drug_score": round(mean_drug, 3),
            "mean_pocket_score": round(mean_score, 3),
            "mean_volume_A3": round(mean_vol, 1),
            "consensus_residues": consensus,
            "dist_to_K1141_A": round(dist(centroid, K1141_CENTRE), 2),
        })
    summaries.sort(key=lambda s: (-s["persistence"], -s["mean_drug_score"]))
    return summaries


def classify_clusters(md_clusters: list[dict],
                      static_clusters: list[dict]) -> list[dict]:
    """Tag each MD cluster:
      • 'K1141_PRIMARY' if centroid within 6 Å of K1141_CENTRE
      • 'STATIC' if a static pocket exists within STATIC_DISTANCE_A
      • 'CRYPTIC' if persistence ≥ PERSISTENCE_THRESHOLD AND
                  no nearby static OR drug_score ≥ 2× best nearby static
      • 'TRANSIENT' if persistence below threshold
    """
    tagged = []
    for mc in md_clusters:
        nearest_static_d = float("inf")
        nearest_static = None
        for sc in static_clusters:
            d = dist(mc["centroid"], sc["centroid"])
            if d < nearest_static_d:
                nearest_static_d = d
                nearest_static = sc

        is_k1141 = dist(mc["centroid"], K1141_CENTRE) <= CLUSTER_RADIUS_A
        is_persistent = mc["persistence"] >= PERSISTENCE_THRESHOLD

        if is_k1141:
            tag = "K1141_PRIMARY"
        elif not is_persistent:
            tag = "TRANSIENT"
        else:
            # Persistent. Cryptic vs STATIC test.
            if (nearest_static is None or
                    nearest_static_d > STATIC_DISTANCE_A):
                tag = "CRYPTIC_NEW"
            elif mc["mean_drug_score"] >= DRUG_SCORE_CRYPTIC_RATIO * (
                    nearest_static["mean_drug_score"] + 0.001):
                tag = "CRYPTIC_UPGRADED"
            else:
                tag = "STATIC_ECHO"

        tagged.append({
            **mc,
            "tag": tag,
            "nearest_static_dist_A": round(nearest_static_d, 2) if nearest_static else None,
            "nearest_static_drug": nearest_static["mean_drug_score"] if nearest_static else None,
        })
    return tagged


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ----- Static baseline: apo starting structure -----
    log("=== STATIC apo fpocket ===")
    with tempfile.TemporaryDirectory() as tmpd:
        tmp = Path(tmpd)
        # Static apo is already clean; skip stripper
        static_pockets = run_fpocket(STATIC_APO, tmp)
    log(f"  {len(static_pockets)} static pockets detected")
    # Single-structure static: each pocket is its own "cluster"
    static_clusters = cluster_pockets(static_pockets)
    for sc in static_clusters[:10]:
        log(f"   static #{sc['cluster_id']:2d} d_score={sc['mean_drug_score']:.2f} "
            f"vol={sc['mean_volume_A3']:.0f} dK1141={sc['dist_to_K1141_A']:.1f}")

    # ----- MD ensemble -----
    snapshots = sorted(SNAPSHOT_DIR.glob("snap_???.pdb"))
    log(f"=== MD ensemble: {len(snapshots)} snapshots ===")
    all_md_pockets: list[Pocket] = []
    for snap in snapshots:
        with tempfile.TemporaryDirectory() as tmpd:
            tmp = Path(tmpd)
            prot = tmp / f"{snap.stem}_prot.pdb"
            n_atoms = strip_water_ions_to_pdb(snap, prot)
            log(f"{snap.stem}: stripped {n_atoms} atoms")
            pockets = run_fpocket(prot, tmp)
        log(f"  → {len(pockets)} pockets")
        all_md_pockets.extend(pockets)

    log(f"Total MD pockets across ensemble: {len(all_md_pockets)}")

    md_clusters = cluster_pockets(all_md_pockets)
    log(f"Clustered into {len(md_clusters)} unique MD pockets")

    tagged = classify_clusters(md_clusters, static_clusters)

    # ----- Report -----
    log("=== TOP-10 MD pockets by persistence × druggability ===")
    log(f"{'tag':17s} {'persist':>7s} {'drug':>6s} {'vol':>6s} "
        f"{'dK1141':>7s} {'dStatic':>8s} {'residues':s}")
    for c in tagged[:15]:
        ds = c["nearest_static_dist_A"]
        ds_s = f"{ds:.1f}" if ds is not None else "—"
        log(f"{c['tag']:17s} {c['persistence']:>7d} "
            f"{c['mean_drug_score']:>6.2f} {c['mean_volume_A3']:>6.0f} "
            f"{c['dist_to_K1141_A']:>7.1f} {ds_s:>8s} "
            f"{c['consensus_residues'][:8]}")

    # Verdict
    cryptic_new = [c for c in tagged if c["tag"] == "CRYPTIC_NEW"]
    cryptic_upg = [c for c in tagged if c["tag"] == "CRYPTIC_UPGRADED"]
    k1141 = next((c for c in tagged if c["tag"] == "K1141_PRIMARY"), None)

    log("=== VERDICT ===")
    if k1141:
        log(f"K1141 primary pocket: persistence {k1141['persistence']}/20, "
            f"drug_score {k1141['mean_drug_score']:.2f}, "
            f"volume {k1141['mean_volume_A3']:.0f} Å³")
    else:
        log("K1141 primary pocket NOT RECOVERED in MD — serious flag")

    log(f"Cryptic NEW pockets (persistent, not in static): {len(cryptic_new)}")
    for c in cryptic_new:
        log(f"   → d_score {c['mean_drug_score']:.2f} at {c['centroid']} "
            f"persist {c['persistence']}/20 vol {c['mean_volume_A3']:.0f}")

    log(f"Cryptic UPGRADED pockets (static echo × 2+ drug_score): {len(cryptic_upg)}")
    for c in cryptic_upg:
        log(f"   → d_score {c['mean_drug_score']:.2f} at {c['centroid']} "
            f"persist {c['persistence']}/20 vol {c['mean_volume_A3']:.0f}")

    payload = {
        "phase": "5c",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_snapshots": len(snapshots),
        "static_clusters": static_clusters,
        "md_clusters": tagged,
        "verdict": {
            "k1141_recovered": k1141 is not None,
            "k1141_summary": k1141,
            "n_cryptic_new": len(cryptic_new),
            "n_cryptic_upgraded": len(cryptic_upg),
            "top_cryptic_new": cryptic_new[:3],
            "top_cryptic_upgraded": cryptic_upg[:3],
        },
        "parameters": {
            "cluster_radius_A": CLUSTER_RADIUS_A,
            "static_distance_threshold_A": STATIC_DISTANCE_A,
            "persistence_threshold": PERSISTENCE_THRESHOLD,
            "drug_score_cryptic_ratio": DRUG_SCORE_CRYPTIC_RATIO,
            "k1141_centre": K1141_CENTRE,
        },
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str))
    log(f"Written {OUT_JSON}")


if __name__ == "__main__":
    main()
