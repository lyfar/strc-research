#!/usr/bin/env python3
"""
OTOA Phase 1A — STRC x OTOA paralog sequence audit.

Question: are STRC and OTOA close enough paralogs that OTOA could be engineered
to take over STRC's hair-bundle function (OTOA paralog cross-rescue hypothesis)?

Cheapest filter: if global identity < ~15%, the paralog hypothesis is dead.
If >=20%, worth pursuing structural comparison (Phase 1B).

Inputs (fetched live via UniProt REST):
  - Q7RTU9 (human stereocilin, STRC)
  - Q7RTW9 (human otoancorin, OTOA)

Outputs:
  - otoa_phase1_sequence_alignment.json
      * global identity / similarity / gaps (Needleman-Wunsch, BLOSUM62)
      * mini-STRC (700-1775) local alignment into OTOA
      * ARM repeats (1603-1770) local alignment into OTOA
      * sliding-window identity profile (50 aa, step 20)
      * UniProt domain features for both proteins
      * mapped fate of anchor residues (K1141, E1659)

Deterministic; no randomness. Gap scores: open -11, extend -1 (BLAST defaults).

Replication:
    /opt/miniconda3/bin/python3 otoa_phase1_sequence_alignment.py
"""
import json
from io import StringIO
from pathlib import Path
from urllib.request import urlopen

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

OUT_DIR = Path(__file__).resolve().parent
OUT_JSON = OUT_DIR / "otoa_phase1_sequence_alignment.json"

STRC_ACC = "Q7RTU9"
OTOA_ACC = "Q7RTW8"  # canonical human OTOA (1153 aa). Q7RTW9 is a different protein.

MINI_STRC_START, MINI_STRC_END = 700, 1775
ARM_START, ARM_END = 1603, 1770
POCKET_RES = [1141, 1659, 1646, 1645]  # pharmacochaperone anchor triangle
WINDOW, STEP = 50, 20

BL62 = substitution_matrices.load("BLOSUM62")


def fetch_uniprot_fasta(acc: str):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urlopen(url) as r:
        data = r.read().decode("utf-8")
    return next(SeqIO.parse(StringIO(data), "fasta"))


def fetch_uniprot_json(acc: str):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    with urlopen(url) as r:
        return json.load(r)


def make_aligner(mode: str = "global"):
    a = Align.PairwiseAligner()
    a.substitution_matrix = BL62
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = mode
    return a


def metrics_from_alignment(aln, s1_orig: str, s2_orig: str) -> dict:
    """Walk aln.aligned block coords (indexes into original ungapped inputs) to
    count matches / similarities. s1_orig/s2_orig are the sequences passed to
    aligner.align()."""
    blocks1 = aln.aligned[0]
    blocks2 = aln.aligned[1]
    matches = 0
    sim = 0
    pair_n = 0
    for (st1, en1), (st2, en2) in zip(blocks1, blocks2):
        block1 = s1_orig[int(st1) : int(en1)]
        block2 = s2_orig[int(st2) : int(en2)]
        for a, b in zip(block1, block2):
            pair_n += 1
            if a == b:
                matches += 1
            try:
                if BL62[a, b] > 0:
                    sim += 1
            except Exception:
                pass
    # gap-included alignment length = total columns covered on either side
    covered_1 = sum(int(en) - int(st) for (st, en) in blocks1)
    covered_2 = sum(int(en) - int(st) for (st, en) in blocks2)
    gaps_in_aligned = max(covered_1, covered_2) - pair_n  # rough lower bound
    return {
        "score": float(aln.score),
        "ungapped_length": int(pair_n),
        "matches": int(matches),
        "identity_pct": round(100 * matches / pair_n, 2) if pair_n else 0.0,
        "similar": int(sim),
        "similarity_pct": round(100 * sim / pair_n, 2) if pair_n else 0.0,
        "block_count": int(len(blocks1)),
        "covered_s1": int(covered_1),
        "covered_s2": int(covered_2),
    }


def local_fragment_alignment(s1_full: str, s2_full: str, start1: int, end1: int):
    frag = s1_full[start1 - 1 : end1]
    aligner = make_aligner("local")
    aln = aligner.align(frag, s2_full)[0]
    m = metrics_from_alignment(aln, frag, s2_full)
    blocks2 = aln.aligned[1]
    otoa_blocks = [[int(b[0]) + 1, int(b[1])] for b in blocks2]  # 1-indexed, inclusive end
    return {
        "strc_start": int(start1),
        "strc_end": int(end1),
        "strc_fragment_len": int(end1 - start1 + 1),
        "otoa_matched_blocks_1idx": otoa_blocks,
        **m,
    }


def map_residue_to_otoa(aln_global, pos_strc: int):
    """For a STRC residue index (1-based), find the OTOA residue it aligned to (if any)."""
    pos0 = pos_strc - 1
    for (st1, en1), (st2, en2) in zip(aln_global.aligned[0], aln_global.aligned[1]):
        st1, en1, st2, en2 = int(st1), int(en1), int(st2), int(en2)
        if st1 <= pos0 < en1:
            offset = pos0 - st1
            return {"otoa_pos": int(st2 + offset + 1), "aligned": True}
    return {"otoa_pos": None, "aligned": False, "note": "falls in a gap"}


def sliding_profile(s1: str, s2: str, window: int, step: int):
    aligner = make_aligner("local")
    profile = []
    for start in range(1, len(s1) - window + 2, step):
        end = start + window - 1
        frag = s1[start - 1 : end]
        try:
            aln = aligner.align(frag, s2)[0]
            blocks1 = aln.aligned[0]
            blocks2 = aln.aligned[1]
            s1_block = frag
            matches = 0
            pair_n = 0
            for (st1, en1), (st2, en2) in zip(blocks1, blocks2):
                for a, b in zip(s1_block[st1:en1], s2[st2:en2]):
                    pair_n += 1
                    if a == b:
                        matches += 1
            idp = 100 * matches / pair_n if pair_n else 0.0
            score = float(aln.score)
        except Exception:
            idp, score = 0.0, 0.0
        profile.append({"start": start, "end": end, "identity_pct": round(idp, 1), "score": round(score, 1)})
    return profile


def main():
    print("Fetching UniProt sequences ...")
    strc_rec = fetch_uniprot_fasta(STRC_ACC)
    otoa_rec = fetch_uniprot_fasta(OTOA_ACC)
    print(f"  STRC {STRC_ACC}: {len(strc_rec.seq)} aa")
    print(f"  OTOA {OTOA_ACC}: {len(otoa_rec.seq)} aa")

    print("Fetching UniProt annotations ...")
    strc_json = fetch_uniprot_json(STRC_ACC)
    otoa_json = fetch_uniprot_json(OTOA_ACC)

    def extract_features(uj):
        return [
            {
                "type": f.get("type"),
                "description": f.get("description"),
                "start": f.get("location", {}).get("start", {}).get("value"),
                "end": f.get("location", {}).get("end", {}).get("value"),
            }
            for f in uj.get("features", [])
            if f.get("type") in ("Domain", "Region", "Repeat", "Signal", "Lipidation", "Glycosylation")
        ]

    strc_feat = extract_features(strc_json)
    otoa_feat = extract_features(otoa_json)

    print("Global pairwise alignment (Needleman-Wunsch, BLOSUM62) ...")
    aligner_g = make_aligner("global")
    s_strc = str(strc_rec.seq)
    s_otoa = str(otoa_rec.seq)
    aln_g = aligner_g.align(s_strc, s_otoa)[0]
    global_metrics = metrics_from_alignment(aln_g, s_strc, s_otoa)
    print(
        f"  Global score {global_metrics['score']:.1f} | "
        f"identity {global_metrics['identity_pct']}% | "
        f"similarity {global_metrics['similarity_pct']}% | "
        f"ungapped {global_metrics['ungapped_length']} residues"
    )

    print("Local: mini-STRC (700-1775) vs OTOA ...")
    mini_metrics = local_fragment_alignment(str(strc_rec.seq), str(otoa_rec.seq), MINI_STRC_START, MINI_STRC_END)
    print(
        f"  mini-STRC local identity {mini_metrics['identity_pct']}% "
        f"over {mini_metrics['ungapped_length']} matched residues"
    )

    print("Local: ARM repeats (1603-1770) vs OTOA ...")
    arm_metrics = local_fragment_alignment(str(strc_rec.seq), str(otoa_rec.seq), ARM_START, ARM_END)
    print(
        f"  ARM local identity {arm_metrics['identity_pct']}% "
        f"over {arm_metrics['ungapped_length']} matched residues; "
        f"hits OTOA blocks {arm_metrics['otoa_matched_blocks_1idx']}"
    )

    print(f"Sliding profile ({WINDOW} aa / step {STEP}) ...")
    prof = sliding_profile(str(strc_rec.seq), str(otoa_rec.seq), WINDOW, STEP)
    peak = max(prof, key=lambda p: p["identity_pct"])
    print(
        f"  windows={len(prof)} peak identity {peak['identity_pct']}% at STRC {peak['start']}-{peak['end']}"
    )

    print("Anchor residues mapped into OTOA (global alignment) ...")
    anchor_map = {str(r): map_residue_to_otoa(aln_g, r) for r in POCKET_RES}
    for r, v in anchor_map.items():
        print(f"  STRC {r} -> OTOA {v.get('otoa_pos')} ({'aligned' if v['aligned'] else 'gap'})")

    result = {
        "meta": {
            "date": "2026-04-20",
            "question": "Are STRC and OTOA close enough paralogs for cross-rescue?",
            "strc_accession": STRC_ACC,
            "otoa_accession": OTOA_ACC,
            "strc_length": len(strc_rec.seq),
            "otoa_length": len(otoa_rec.seq),
            "aligner": "Biopython PairwiseAligner",
            "matrix": "BLOSUM62",
            "gap_open": -11,
            "gap_extend": -1,
            "window": WINDOW,
            "step": STEP,
        },
        "global_alignment": global_metrics,
        "mini_strc_vs_otoa_local": mini_metrics,
        "arm_repeats_vs_otoa_local": arm_metrics,
        "sliding_profile": prof,
        "anchor_residues_map": anchor_map,
        "strc_domains": strc_feat,
        "otoa_domains": otoa_feat,
    }

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {OUT_JSON}")
    verdict = "PURSUE" if global_metrics["identity_pct"] >= 15 else "KILL"
    print(f"\nVerdict (global identity >=15% gate): {verdict}")


if __name__ == "__main__":
    main()
