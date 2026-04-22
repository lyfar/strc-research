#!/usr/bin/env python3
"""
ASO Phase 2 — STRCP1 cross-hybridization scan.

Question: do any of the 16 Phase 1 ASO candidates also bind the STRCP1
pseudogene transcript (chr15 ~43.69-43.72 Mb hg38, 97% identity to STRC)?

Same risk that just downgraded PE from A→B (see [[STRC PE Phase4 STRCP1
Paralog Off-Target]]). If an ASO targets a window with ≤2 mismatches in
STRCP1 mRNA, RNase-H1 will cleave both transcripts → loss of specificity.

Method:
  1. Load 16 ASO candidates from aso_phase1_design.json. Each has
     `target_sense` = pre-mRNA sequence the ASO hybridises to (5'→3' sense).
  2. Load hg38 chr15 FASTA into uint8 array, build reverse complement.
  3. For every ASO target_sense, slide along + and - strands of chr15,
     count Hamming-distance hits at <=2 mismatches.
  4. Bin hits by region:
       - STRC locus (~43.58-43.64 Mb hg38) — expected, on-target
       - STRCP1 locus (~43.68-43.72 Mb hg38) — paralog risk
       - elsewhere chr15 — random off-targets (low priority)
  5. Gate per ASO:
       - 0 STRCP1 hits at <=2 mm AND >=1 STRC hit at 0 mm → PASS
       - any STRCP1 hit at <=2 mm → FAIL (paralog risk)
       - 0 STRC hits at 0 mm → FAIL (lost on-target)

Note: only chr15 is scanned. Whole-genome BLAST is the next step if we
want to clear ASO for IND-enabling work, but STRCP1 is the dominant risk
because of the 97% paralog identity.

Run env: any conda env with numpy. No external bioinformatics tool needed.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

WORK_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
PHASE1_JSON = WORK_DIR / "aso_phase1_design.json"
OUT_JSON = WORK_DIR / "aso_phase2_strcp1_specificity.json"
GENOME_PATH = Path("/Users/egorlyfar/Brain/research/strc/genomes/hg38_chr15.fa")

MAX_MISMATCH = 2

# hg38 coordinates (verified via PE Phase 4 on-target locus chr15:43600548
# and STRCP1 paralog cluster chr15:43700346):
STRC_START   = 43_580_000
STRC_END     = 43_640_000
STRCP1_START = 43_680_000
STRCP1_END   = 43_720_000


def log(msg):
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}", flush=True)


def load_genome_uint8(path: Path):
    log(f"loading genome from {path}")
    parts = []
    chrom_name = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                chrom_name = line[1:].strip().split()[0]
            else:
                parts.append(line.strip().upper())
    seq = "".join(parts)
    log(f"  genome chrom={chrom_name} length={len(seq):,}")
    table = np.full(256, 4, dtype=np.uint8)
    table[ord("A")] = 0
    table[ord("C")] = 1
    table[ord("G")] = 2
    table[ord("T")] = 3
    arr_b = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return chrom_name, table[arr_b]


def revcomp_uint8(arr):
    comp = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
    return comp[arr][::-1]


def encode_query(seq: str):
    table = {"A": 0, "C": 1, "G": 2, "T": 3}
    return np.array([table.get(c, 4) for c in seq], dtype=np.uint8)


def search_strand(strand_arr, query_arr, max_mm: int):
    """Return list of (pos, mismatches) for hits with <=max_mm Hamming distance.

    Uses sliding-window vectorised comparison: build a 2D mismatch matrix
    only at the end via column-wise inequality counts.
    """
    L = len(strand_arr)
    qlen = len(query_arr)
    nwin = L - qlen + 1
    if nwin <= 0:
        return []
    # Sum mismatches column by column to keep memory bounded.
    mm = np.zeros(nwin, dtype=np.int32)
    for k in range(qlen):
        col = strand_arr[k : k + nwin]
        mm += (col != query_arr[k]).astype(np.int32)
        # Skip wildcard (N=4) — count as mismatch; conservative.
    hit_idx = np.where(mm <= max_mm)[0]
    return [(int(i), int(mm[i])) for i in hit_idx]


def classify_pos(pos_fwd: int) -> str:
    """Bin a + strand 0-indexed start position into a region label."""
    if STRC_START <= pos_fwd <= STRC_END:
        return "STRC"
    if STRCP1_START <= pos_fwd <= STRCP1_END:
        return "STRCP1"
    return "other_chr15"


def fwd_pos_from_rev(rev_pos: int, qlen: int, L: int) -> int:
    """Translate a position on the reverse-complement strand back to the +
    strand 0-indexed coordinate of the same window."""
    return L - rev_pos - qlen


def main():
    log("ASO Phase 2 — STRCP1 cross-hybridization scan starting")
    if not PHASE1_JSON.exists():
        log(f"FATAL — Phase 1 input missing: {PHASE1_JSON}")
        sys.exit(2)
    if not GENOME_PATH.exists():
        log(f"FATAL — genome FASTA missing: {GENOME_PATH}")
        sys.exit(2)

    with open(PHASE1_JSON) as f:
        phase1 = json.load(f)

    # Flatten candidates across all events.
    candidates = []
    for ev in phase1["events"]:
        for cand in ev["candidates"]:
            candidates.append({
                "event": ev["label"],
                "priority": ev["priority"],
                "window_tag": cand["window_tag"],
                "aso_5to3": cand["aso_5to3"],
                "target_sense": cand["target_sense"],
                "length": cand["length"],
                "tm_c": cand["tm_c"],
                "gc_pct": cand["gc_pct"],
            })
    log(f"loaded {len(candidates)} ASO candidates from {len(phase1['events'])} events")

    chrom, fwd = load_genome_uint8(GENOME_PATH)
    rev = revcomp_uint8(fwd)
    L = len(fwd)

    results = []
    pass_count = 0
    fail_paralog = 0
    fail_lost_target = 0

    for i, cand in enumerate(candidates, 1):
        q_arr = encode_query(cand["target_sense"])
        qlen = len(q_arr)

        fwd_hits = search_strand(fwd, q_arr, MAX_MISMATCH)
        rev_hits = search_strand(rev, q_arr, MAX_MISMATCH)

        bins = {"STRC": {0: 0, 1: 0, 2: 0},
                "STRCP1": {0: 0, 1: 0, 2: 0},
                "other_chr15": {0: 0, 1: 0, 2: 0}}
        sample_hits = []

        for pos, mm in fwd_hits:
            region = classify_pos(pos)
            bins[region][mm] += 1
            if len(sample_hits) < 6:
                sample_hits.append({"strand": "+", "pos_chr15_0idx": pos, "mismatches": mm, "region": region})
        for rpos, mm in rev_hits:
            pos_fwd = fwd_pos_from_rev(rpos, qlen, L)
            region = classify_pos(pos_fwd)
            bins[region][mm] += 1
            if len(sample_hits) < 6:
                sample_hits.append({"strand": "-", "pos_chr15_0idx": pos_fwd, "mismatches": mm, "region": region})

        strcp1_total = sum(bins["STRCP1"].values())
        strc_perfect = bins["STRC"][0]
        if strcp1_total > 0:
            verdict = "FAIL_paralog"
            fail_paralog += 1
        elif strc_perfect == 0:
            verdict = "FAIL_lost_target"
            fail_lost_target += 1
        else:
            verdict = "PASS"
            pass_count += 1

        results.append({
            "rank": i,
            "event": cand["event"],
            "window_tag": cand["window_tag"],
            "aso_5to3": cand["aso_5to3"],
            "length": cand["length"],
            "target_sense": cand["target_sense"],
            "verdict": verdict,
            "hits": bins,
            "sample_hits": sample_hits,
        })
        log(f"  [{i:02d}] {cand['event']}/{cand['window_tag']:14s} {cand['aso_5to3']:24s} "
            f"STRC={bins['STRC']} STRCP1={bins['STRCP1']} other={bins['other_chr15']} -> {verdict}")

    summary = {
        "meta": {
            "date": datetime.now(timezone.utc).date().isoformat(),
            "genome": str(GENOME_PATH),
            "chrom": chrom,
            "max_mismatch": MAX_MISMATCH,
            "strc_window_hg38": [STRC_START, STRC_END],
            "strcp1_window_hg38": [STRCP1_START, STRCP1_END],
        },
        "totals": {
            "n_candidates": len(candidates),
            "PASS": pass_count,
            "FAIL_paralog": fail_paralog,
            "FAIL_lost_target": fail_lost_target,
        },
        "results": results,
    }

    OUT_JSON.write_text(json.dumps(summary, indent=2))
    log(f"wrote {OUT_JSON}")
    log(f"verdict — PASS={pass_count} FAIL_paralog={fail_paralog} FAIL_lost_target={fail_lost_target}")


if __name__ == "__main__":
    main()
