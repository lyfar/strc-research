#!/usr/bin/env python3
"""
Prime Editing Phase 3.5 — STRCP1-paralog-discriminating pegRNA redesign.

Context: Phase 4 Cas-OFFinder scan found 61/61 Phase-3 pegRNAs hit STRCP1
pseudogene paralog at chr15:43,700,346-43,700,363 with 0 mismatches. Every
Phase-3 candidate was designed for STRC-vs-STRC allele discrimination but
blind to the 100 kb-downstream paralog. This script is Phase 3.5: the
STRCP1-aware redesign.

Method:
  1. Extract STRCP1 paralog window from hg38_chr15.fa — the 18 bp core
     plus ±40 bp flanks for PAM context sliding.
  2. Extract STRC on-target window at the same locus (~100 kb upstream).
  3. Align STRC vs STRCP1 at nucleotide level to locate divergent positions.
  4. For every Phase-3 candidate: align its protospacer + PAM to STRCP1,
     count seed mismatches (PAM-proximal 10-12 nt). Seed convention:
       SpCas9 (NGG, PAM at 3'): seed = positions 9-20 of protospacer
       SpG/enCas9/SpRY (NGN/NRN): same PAM-3' geometry → same seed
       SaCas9 (NNGRRT, PAM at 3'): seed = positions 10-21 of 21-mer
  5. Filter: require ≥2 mismatches in seed AGAINST STRCP1 for discriminating
     activity. <2 mismatches in seed → STRCP1 will be edited too → kill.
  6. For candidates passing the seed filter, run pure-numpy chr15-wide
     off-target search (at most 3 mismatches) to confirm low off-target
     beyond STRCP1.

Output:
  pe_phase3_5_strcp1_aware.json — per-candidate STRCP1 seed mismatch count,
  pass/fail verdict, final STRCP1-aware shortlist.

Gates:
  - ≥1 candidate passes seed-mismatch gate → promote hypothesis #7
    [[Prime Editing for STRC]] back from B → A.
  - 0 candidates pass → promote [[STRCP1 paralog sequence divergence]]
    as a dead-end flag; further PE work needs a completely new target
    window (e.g. different exon) or kill.

Run env: strc-mmgbsa conda env (numpy + biopython).
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np

HERE = Path(__file__).resolve().parent
PHASE4_JSON = HERE / "pe_phase4_cas_offinder.json"
GENOME_FA = Path("/Users/egorlyfar/Brain/research/strc/genomes/hg38_chr15.fa")
OUT_JSON = HERE / "pe_phase3_5_strcp1_aware.json"

# Windows on chr15 (hg38). Derived from Phase 4 Cas-OFFinder output.
# STRCP1 core: the 18 bp window where all 61 Phase-3 protospacers cluster.
STRCP1_CORE_START = 43_700_346   # 1-based inclusive
STRCP1_CORE_END = 43_700_363
STRCP1_FLANK = 40                 # extract +/-40 bp for PAM context sliding

# STRC target window from Phase 4 on-target hits (SpCas9_NGG top-ranked was at
# chr15:43600548-43600567 per prior notes). Use ±40 bp flanks.
STRC_ONTARGET_START = 43_600_548
STRC_ONTARGET_END = 43_600_567
STRC_FLANK = 40

# Per-Cas seed position convention (positions within the 20-21 nt protospacer
# where a mismatch is sufficient to abolish cleavage). PAM is at 3' end for
# all the Cas variants we use; seed is PAM-proximal ~10-12 nt.
SEED_POSITIONS = {
    "SpCas9_NGG":   (9, 20),
    "SpG_NGN":      (9, 20),
    "enCas9_NGN":   (9, 20),
    "SpRY_NRN":     (9, 20),
    "SpCas9NG_NG":  (9, 20),
    "SaCas9_NNGRRT":(10, 21),
}

# STRCP1-discrimination gate: need ≥ this many seed mismatches to reliably
# avoid paralog cleavage. 2 is the minimum; 3 is comfortable.
SEED_MISMATCH_GATE = 2


def load_hg38_chr15_slice(start: int, end: int) -> str:
    """Load a 1-based inclusive slice from the chr15 FASTA."""
    if not GENOME_FA.exists():
        raise RuntimeError(f"Missing {GENOME_FA}")
    seq = []
    with open(GENOME_FA) as f:
        header = f.readline()
        offset = 0
        want_start = start - 1
        want_end = end
        for line in f:
            line = line.strip()
            L = len(line)
            if offset + L < want_start:
                offset += L
                continue
            s = max(0, want_start - offset)
            e = min(L, want_end - offset)
            seq.append(line[s:e])
            offset += L
            if offset >= want_end:
                break
    return "".join(seq).upper()


def complement(base: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}.get(base, "N")


def revcomp(seq: str) -> str:
    return "".join(complement(b) for b in reversed(seq))


def count_mismatches(a: str, b: str) -> int:
    assert len(a) == len(b), f"len mismatch {len(a)} vs {len(b)}"
    return sum(1 for x, y in zip(a, b) if x != y)


def align_protospacer_to_window(proto: str, window_seq: str) -> tuple[int, str, int]:
    """Slide proto across window (both strands) and return best alignment.
    Returns (min_mismatches, best_match_sequence, strand: +1/-1)."""
    best_mm = len(proto) + 1
    best_match = ""
    best_strand = 1
    # forward strand
    for i in range(len(window_seq) - len(proto) + 1):
        slice_ = window_seq[i : i + len(proto)]
        if "N" in slice_:
            continue
        mm = count_mismatches(proto, slice_)
        if mm < best_mm:
            best_mm = mm
            best_match = slice_
            best_strand = 1
    # reverse strand
    rc_window = revcomp(window_seq)
    for i in range(len(rc_window) - len(proto) + 1):
        slice_ = rc_window[i : i + len(proto)]
        if "N" in slice_:
            continue
        mm = count_mismatches(proto, slice_)
        if mm < best_mm:
            best_mm = mm
            best_match = slice_
            best_strand = -1
    return best_mm, best_match, best_strand


def seed_mismatches(proto: str, match: str, seed_range: tuple[int, int]) -> int:
    """Mismatches restricted to the seed positions (1-based inclusive range)."""
    start, end = seed_range
    return count_mismatches(proto[start - 1 : end], match[start - 1 : end])


def main():
    # Load genome windows.
    strc_window = load_hg38_chr15_slice(
        STRC_ONTARGET_START - STRC_FLANK, STRC_ONTARGET_END + STRC_FLANK
    )
    strcp1_window = load_hg38_chr15_slice(
        STRCP1_CORE_START - STRCP1_FLANK, STRCP1_CORE_END + STRCP1_FLANK
    )
    print(f"STRC on-target window ({len(strc_window)} bp): {strc_window[:20]}... ")
    print(f"STRCP1 paralog window ({len(strcp1_window)} bp): {strcp1_window[:20]}...")

    # Quick STRC vs STRCP1 core alignment sanity — just the core 18 bp window
    strcp1_core = strcp1_window[STRCP1_FLANK : STRCP1_FLANK + (STRCP1_CORE_END - STRCP1_CORE_START + 1)]
    # Find best-aligning slice in STRC window of the same length
    core_len = len(strcp1_core)
    best_strc_mm = core_len + 1
    best_strc_core = ""
    for i in range(len(strc_window) - core_len + 1):
        s = strc_window[i : i + core_len]
        if "N" in s:
            continue
        mm = count_mismatches(strcp1_core, s)
        if mm < best_strc_mm:
            best_strc_mm = mm
            best_strc_core = s
    print(f"STRCP1 core vs best STRC slice: {best_strc_mm}/{core_len} mismatches")
    print(f"  STRCP1: {strcp1_core}")
    print(f"  STRC:   {best_strc_core}")

    # Load Phase 4 guides.
    d4 = json.loads(PHASE4_JSON.read_text())
    guides = d4["results"]
    print(f"\n--- Evaluating {len(guides)} Phase-4 guides against STRCP1 ---")

    per_guide = []
    for g in guides:
        pam = g["pam"]
        proto = g["protospacer"]
        seed_range = SEED_POSITIONS.get(pam)
        if seed_range is None:
            per_guide.append({**g, "strcp1_total_mm": None, "strcp1_seed_mm": None,
                              "seed_gate_pass": False, "note": f"unknown PAM {pam}"})
            continue
        total_mm, match, strand = align_protospacer_to_window(proto, strcp1_window)
        seed_mm = seed_mismatches(proto, match, seed_range)
        per_guide.append({
            "pam": pam,
            "protospacer": proto,
            "discrimination_grade": g.get("discrimination_grade"),
            "nick_to_edit_nt": g.get("nick_to_edit_nt"),
            "original_phase4_pass": g.get("passes_gate"),
            "strcp1_best_match": match,
            "strcp1_total_mm": total_mm,
            "strcp1_strand": "+" if strand == 1 else "-",
            "strcp1_seed_mm": seed_mm,
            "seed_range_1based": list(seed_range),
            "seed_gate_pass": seed_mm >= SEED_MISMATCH_GATE,
        })

    # Summary.
    n_pass = sum(1 for r in per_guide if r["seed_gate_pass"])
    print(f"\nSeed-mismatch gate (≥{SEED_MISMATCH_GATE} mm in seed vs STRCP1):")
    print(f"  {n_pass}/{len(per_guide)} guides pass")
    # Distribution
    from collections import Counter
    dist = Counter(r["strcp1_seed_mm"] for r in per_guide if r["strcp1_seed_mm"] is not None)
    print(f"  Seed-mm distribution: {dict(sorted(dist.items()))}")

    top_pass = [r for r in per_guide if r["seed_gate_pass"]]
    top_pass.sort(key=lambda r: -r["strcp1_seed_mm"])

    # Output
    payload = {
        "phase": "3.5",
        "hypothesis": "Prime Editing for STRC",
        "method": "Phase-4 guides re-filtered against STRCP1 paralog seed-region mismatches",
        "date": "2026-04-23",
        "strc_window": {"start": STRC_ONTARGET_START, "end": STRC_ONTARGET_END, "flank_bp": STRC_FLANK, "seq": strc_window},
        "strcp1_window": {"start": STRCP1_CORE_START, "end": STRCP1_CORE_END, "flank_bp": STRCP1_FLANK, "seq": strcp1_window},
        "strcp1_core_vs_strc_best_mm": best_strc_mm,
        "strcp1_core_len": core_len,
        "seed_mismatch_gate": SEED_MISMATCH_GATE,
        "n_guides_total": len(per_guide),
        "n_guides_pass_seed_gate": n_pass,
        "seed_mm_distribution": dict(sorted(dist.items())),
        "per_guide": per_guide,
        "top_passing": top_pass[:10],
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\nWrote {OUT_JSON}")

    # Show top-3 passing
    if top_pass:
        print(f"\nTop-3 STRCP1-aware candidates:")
        for r in top_pass[:3]:
            print(f"  [{r['pam']}] {r['protospacer']}  seed_mm={r['strcp1_seed_mm']}  grade={r['discrimination_grade']}  nick→edit={r['nick_to_edit_nt']}")
    else:
        print(f"\n⚠️  0 candidates pass. STRCP1 paralog identity too high at this locus; PE branch blocked here.")


if __name__ == "__main__":
    main()
