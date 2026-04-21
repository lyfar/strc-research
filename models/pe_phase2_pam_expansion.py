#!/usr/bin/env python3
"""
PE Phase 2 — engineered-Cas9 PAM expansion sweep for STRC E1659A.

Phase 1 ([[STRC PE Phase1 pegRNA E1659A]]) identified two unavoidable
constraints with SpCas9 NGG:
  1. Nick-to-edit = 14 nt (Anzalone optimum 3-10 nt) → ~50% efficiency penalty
  2. No PE3b nicking sgRNA exists (no -strand NGG spans variant within 100 nt)

This script asks: do engineered Cas9 variants with relaxed PAM requirements
place a nick 3-10 nt from chr15:43600551 AND admit a PE3b-eligible nicking
sgRNA? If yes, those variants unblock efficient PE for Misha's allele.

Variants surveyed:
  - SpCas9 NGG (baseline — what we already have)
  - SpG        NGN  (Walton 2020)
  - SpRY       NRN  (Walton 2020; near-PAMless, NYN tolerated but weaker)
  - SpCas9-NG  NG   (Nishimasu 2018)
  - enCas9     NGN  (Hu 2018 enhanced NGG; equiv. to SpG for this purpose)
  - xCas9 3.7  NG/GAT/GAA  (Hu 2018; a subset of NG)
  - SaCas9     NNGRRT (Ran 2015, 21-24 nt protospacer)

Only SpCas9-family variants use 20-nt protospacer + 3 nt PAM + nick at pos 17/18.
SaCas9 uses 21-nt protospacer + 6 nt PAM + nick at pos 17 (consult Ran 2015 for
exact geometry). Here we use a 21-nt protospacer + nick between 17 and 18
for SaCas9 as a first-order model; the 22-24 nt protospacer range exists but
most published designs use 21.

All analyses on + strand genomic coords (GRCh38). Variant = chr15:43600551.

Replication:
    /opt/miniconda3/bin/python3 pe_phase2_pam_expansion.py
"""
import json
import re
from pathlib import Path
from urllib.request import Request, urlopen

OUT_DIR = Path(__file__).resolve().parent
OUT_JSON = OUT_DIR / "pe_phase2_pam_expansion.json"
ENSEMBL = "https://rest.ensembl.org"

VARIANT_POS = 43600551  # chr15 + strand
REGION_START = 43600400
REGION_END = 43600700

# Anzalone-optimum nick-to-edit window
NICK_TO_EDIT_OPT = (3, 10)
# PE3b protospacer-spans-variant + nick distance window (PE2/PE3 nicker zone)
NICKER_DIST_RANGE = (30, 100)

COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "R": "Y", "Y": "R"}


def revcomp(s):
    return "".join(COMP.get(b, "N") for b in s[::-1])


def fetch_region(chrom, start, end, strand=1):
    url = f"{ENSEMBL}/sequence/region/human/{chrom}:{start}..{end}:{strand}?coord_system_version=GRCh38"
    with urlopen(Request(url, headers={"Accept": "application/json"})) as r:
        return json.load(r)["seq"].upper()


def gc_pct(s):
    if not s:
        return 0.0
    return round(100 * (s.count("G") + s.count("C")) / len(s), 1)


# PAM match functions. Each returns True if candidate triplet (or 6mer for Sa)
# matches the variant's PAM pattern. Operate on + strand sequence read 5'->3'
# (for + strand PAMs) or reverse-complemented (for - strand PAMs).
def matches_ngg(pam):
    return len(pam) == 3 and pam[1:] == "GG"


def matches_ngn(pam):
    return len(pam) == 3 and pam[1] == "G"


def matches_nrn(pam):
    # SpRY prefers NRN (R = A/G); N/Y/N weaker but tolerated
    return len(pam) == 3 and pam[1] in "AG"


def matches_ng(pam):
    return len(pam) >= 2 and pam[1] == "G"


def matches_nngrrt(pam):
    # SaCas9 PAM = NNGRRT (6 nt)
    return (
        len(pam) == 6
        and pam[2] == "G"
        and pam[3] in "AG"
        and pam[4] in "AG"
        and pam[5] == "T"
    )


VARIANTS = {
    "SpCas9_NGG": {"pam_len": 3, "match": matches_ngg, "protospacer_len": 20, "nick_from_3end": 3, "baseline": True},
    "SpG_NGN": {"pam_len": 3, "match": matches_ngn, "protospacer_len": 20, "nick_from_3end": 3, "baseline": False},
    "enCas9_NGN": {"pam_len": 3, "match": matches_ngn, "protospacer_len": 20, "nick_from_3end": 3, "baseline": False},
    "SpRY_NRN": {"pam_len": 3, "match": matches_nrn, "protospacer_len": 20, "nick_from_3end": 3, "baseline": False},
    "SpCas9NG_NG": {"pam_len": 2, "match": matches_ng, "protospacer_len": 20, "nick_from_3end": 3, "baseline": False},
    "SaCas9_NNGRRT": {"pam_len": 6, "match": matches_nngrrt, "protospacer_len": 21, "nick_from_3end": 3, "baseline": False},
}


def scan_plus_strand_pams(plus_seq, region_start_1idx, variant: dict):
    """PAMs on + strand → pegRNA targets + strand, nick on + strand.
    Protospacer 5' of PAM (upstream in + strand coords)."""
    pam_len = variant["pam_len"]
    proto_len = variant["protospacer_len"]
    nick_from_3end = variant["nick_from_3end"]
    match_fn = variant["match"]
    out = []
    for i in range(proto_len, len(plus_seq) - pam_len + 1):
        pam = plus_seq[i : i + pam_len]
        if not match_fn(pam):
            continue
        proto = plus_seq[i - proto_len : i]
        pam_start = region_start_1idx + i  # 1-idx
        proto_start = region_start_1idx + i - proto_len
        proto_end = region_start_1idx + i - 1
        # Nick between protospacer positions (proto_len - nick_from_3end) and
        # (proto_len - nick_from_3end + 1) counted from 5'; i.e., 3 nt upstream
        # of PAM.
        nick_plus_5 = proto_end - nick_from_3end + 1 - 1  # nt just before nick
        nick_plus_3 = nick_plus_5 + 1  # nt just after nick
        # nick_plus_3 — nick_plus_5 == 1. The edit is downstream of nick on +.
        dist = VARIANT_POS - nick_plus_3  # can be negative if PAM past variant
        if abs(dist) > 200:
            continue
        out.append({
            "strand": "+",
            "pam_seq_plus": pam,
            "pam_coords_plus": [pam_start, pam_start + pam_len - 1],
            "protospacer_plus_5to3": proto,
            "protospacer_coords_plus": [proto_start, proto_end],
            "protospacer_gc_pct": gc_pct(proto),
            "nick_between_plus": [nick_plus_5, nick_plus_3],
            "nick_to_edit_nt": dist,
        })
    return out


def scan_minus_strand_pams(plus_seq, region_start_1idx, variant: dict):
    """PAMs on - strand → pegRNA targets - strand, nick on - strand.
    On + strand the PAM appears as revcomp. Protospacer is 3' of PAM in +strand
    coords. Nick on - strand: between positions 17 and 18 of the - strand
    protospacer counting from its 5' = + strand 3' edge of protospacer minus
    nick_from_3end."""
    pam_len = variant["pam_len"]
    proto_len = variant["protospacer_len"]
    nick_from_3end = variant["nick_from_3end"]
    match_fn = variant["match"]
    out = []
    for i in range(0, len(plus_seq) - pam_len - proto_len + 1):
        pam_plus = plus_seq[i : i + pam_len]
        pam_minus = revcomp(pam_plus)
        if not match_fn(pam_minus):
            continue
        # - strand protospacer is 3' of PAM on - strand = upstream of PAM on +
        # On + strand, protospacer region = [i + pam_len, i + pam_len + proto_len - 1]
        proto_plus_start = i + pam_len
        proto_plus_end = proto_plus_start + proto_len - 1
        if proto_plus_end >= len(plus_seq):
            continue
        proto_plus = plus_seq[proto_plus_start : proto_plus_end + 1]
        proto_minus = revcomp(proto_plus)
        pam_start = region_start_1idx + i
        proto_start_1idx = region_start_1idx + proto_plus_start
        proto_end_1idx = region_start_1idx + proto_plus_end

        # Nick on - strand between protospacer pos (proto_len - nick_from_3end)
        # and (proto_len - nick_from_3end + 1) counted from 5' of - strand.
        # 5' of - strand protospacer = + strand position proto_end_1idx.
        # Position k from 5' of - strand (1-indexed) = + strand position
        # proto_end_1idx - (k - 1).
        # Nick between positions (proto_len - nick_from_3end) and (+1).
        nick_minus_pos_k = proto_len - nick_from_3end  # 1-indexed from 5'
        nick_3_plus = proto_end_1idx - (nick_minus_pos_k - 1)
        nick_5_plus = nick_3_plus - 1
        # Edit vs nick: measured on - strand, edit at position 43600551.
        # dist = signed distance (nt) from nick to edit, counted on - strand
        # reading 5' -> 3' (i.e., + strand high to low coordinate).
        # On - strand, downstream of nick = lower + strand coord.
        # Edit on - strand is "downstream" if + strand position < nick_3_plus? No:
        # - strand reading direction is + strand high coord to low. So
        # downstream of nick on - strand means + strand coord LOWER than
        # nick_3_plus? Wait, nick_3_plus is the + strand coord of the
        # nucleotide IMMEDIATELY 3' OF NICK on - strand — so it is the LOWER
        # + strand coord of the two flanking nt. Let me re-derive.
        # Actually: - strand reading 5'->3' proceeds in the direction of
        # decreasing + strand coordinate. So if nick_5_plus=x and nick_3_plus=x+1,
        # on - strand the 5' nucleotide of the nick flank is at + strand x+1
        # and the 3' is at + strand x. We had nick_3_plus = higher coord
        # above... let me fix the naming. Easier: dist = distance in nt from
        # nick on - strand. Measured as abs diff in + strand coords between
        # nick flank and variant, since + strand and - strand distances are
        # identical in absolute terms.
        dist = abs(VARIANT_POS - ((nick_5_plus + nick_3_plus) / 2))
        # Protospacer-spans-variant test: variant position within + strand
        # protospacer range [proto_plus_start .. proto_plus_end]
        spans = proto_start_1idx <= VARIANT_POS <= proto_end_1idx
        out.append({
            "strand": "-",
            "pam_seq_minus_5to3": pam_minus,
            "pam_coords_plus": [pam_start, pam_start + pam_len - 1],
            "protospacer_minus_5to3": proto_minus,
            "protospacer_coords_plus": [proto_start_1idx, proto_end_1idx],
            "protospacer_gc_pct": gc_pct(proto_minus),
            "nick_between_plus": sorted([nick_5_plus, nick_3_plus]),
            "nick_to_edit_nt": round(dist, 1),
            "protospacer_spans_variant_for_PE3b": spans,
        })
    return out


def main():
    print(f"Fetching chr15:{REGION_START}-{REGION_END} ...")
    plus_seq = fetch_region("15", REGION_START, REGION_END, 1)
    print(f"  region {len(plus_seq)} nt; variant at offset {VARIANT_POS - REGION_START} = '{plus_seq[VARIANT_POS - REGION_START]}' (expect T)")
    assert plus_seq[VARIANT_POS - REGION_START] == "T"

    results = {
        "meta": {
            "date": "2026-04-21",
            "variant_chr15_plus_strand_pos": VARIANT_POS,
            "anzalone_nick_to_edit_optimum_nt": list(NICK_TO_EDIT_OPT),
            "region_scanned": [REGION_START, REGION_END],
            "cas9_variants_scanned": list(VARIANTS.keys()),
        },
        "variants": {},
    }

    summary_rows = []

    for var_name, var in VARIANTS.items():
        print(f"\n=== {var_name} ===")
        plus_hits = scan_plus_strand_pams(plus_seq, REGION_START, var)
        minus_hits = scan_minus_strand_pams(plus_seq, REGION_START, var)

        # For PE pegRNA: we want a + strand OR - strand PAM whose nick is
        # close to the edit (3-10 nt). Filter.
        opt_pegrna = []
        for h in plus_hits:
            d = h["nick_to_edit_nt"]
            if NICK_TO_EDIT_OPT[0] <= d <= NICK_TO_EDIT_OPT[1]:
                h2 = {**h, "pegrna_fit": "optimum", "edit_downstream_of_nick": True}
                opt_pegrna.append(h2)
        for h in minus_hits:
            # For - strand pegRNA, nick on - strand, edit at VARIANT_POS on +
            # = position on - strand too. Distance symmetric.
            d = h["nick_to_edit_nt"]
            if NICK_TO_EDIT_OPT[0] <= d <= NICK_TO_EDIT_OPT[1]:
                h2 = {**h, "pegrna_fit": "optimum"}
                opt_pegrna.append(h2)

        # PE3b nicker: protospacer spans variant (distance from nick to edit is
        # naturally <= 20 nt since protospacer is 20 nt). No separate distance
        # filter needed — spanning implies proximity.
        pe3b_on_minus = [
            h for h in minus_hits if h["protospacer_spans_variant_for_PE3b"]
        ]
        pe3b_on_plus = []
        for h in plus_hits:
            proto_span = h["protospacer_coords_plus"]
            spans = proto_span[0] <= VARIANT_POS <= proto_span[1]
            if spans:
                pe3b_on_plus.append({**h, "protospacer_spans_variant_for_PE3b": True})

        n_plus = len(plus_hits)
        n_minus = len(minus_hits)
        n_opt = len(opt_pegrna)
        n_pe3b_minus = len(pe3b_on_minus)
        n_pe3b_plus = len(pe3b_on_plus)
        print(f"  total PAMs: + strand {n_plus}, - strand {n_minus}")
        print(f"  pegRNA-OPTIMAL (nick-to-edit ∈ {NICK_TO_EDIT_OPT}): {n_opt}")
        if n_opt:
            for h in sorted(opt_pegrna, key=lambda x: x["nick_to_edit_nt"])[:4]:
                proto = h.get("protospacer_plus_5to3") or h.get("protospacer_minus_5to3")
                pam_show = h.get("pam_seq_plus") or h.get("pam_seq_minus_5to3")
                print(f"    ({h['strand']}) PAM {pam_show} dist={h['nick_to_edit_nt']} proto={proto} GC={h['protospacer_gc_pct']}%")
        print(f"  PE3b-eligible nickers (span variant + 30-100 nt): - strand {n_pe3b_minus}, + strand {n_pe3b_plus}")
        for h in pe3b_on_minus[:3]:
            print(f"    [- PE3b] PAM {h['pam_seq_minus_5to3']} dist={h['nick_to_edit_nt']} proto={h['protospacer_minus_5to3']} GC={h['protospacer_gc_pct']}%")
        for h in pe3b_on_plus[:3]:
            print(f"    [+ PE3b] PAM {h['pam_seq_plus']} dist={abs(h['nick_to_edit_nt'])} proto={h['protospacer_plus_5to3']} GC={h['protospacer_gc_pct']}%")

        results["variants"][var_name] = {
            "total_plus_pams": n_plus,
            "total_minus_pams": n_minus,
            "pegrna_optimal": opt_pegrna[:10],
            "pe3b_nickers_minus_strand": pe3b_on_minus[:10],
            "pe3b_nickers_plus_strand": pe3b_on_plus[:10],
        }

        summary_rows.append({
            "variant": var_name,
            "baseline": var["baseline"],
            "n_total_pams": n_plus + n_minus,
            "n_pegrna_in_3to10nt_window": n_opt,
            "n_pe3b_nickers": n_pe3b_minus + n_pe3b_plus,
            "best_pegrna_dist_nt": (
                min((h["nick_to_edit_nt"] if h["nick_to_edit_nt"] > 0 else abs(h["nick_to_edit_nt"])) for h in opt_pegrna) if opt_pegrna else None
            ),
        })

    results["summary"] = summary_rows
    # Comparative table
    print("\n=== SUMMARY ===")
    print(f"{'Variant':<18} {'Total PAMs':>10} {'pegRNA 3-10':>12} {'PE3b nickers':>13}")
    for s in summary_rows:
        print(f"{s['variant']:<18} {s['n_total_pams']:>10d} {s['n_pegrna_in_3to10nt_window']:>12d} {s['n_pe3b_nickers']:>13d}")

    # Verdict: which variant first unblocks both constraints?
    unblockers = [s for s in summary_rows if s["n_pegrna_in_3to10nt_window"] > 0 and s["n_pe3b_nickers"] > 0]
    if unblockers:
        unblockers.sort(key=lambda x: (-x["n_pegrna_in_3to10nt_window"], -x["n_pe3b_nickers"]))
        top = unblockers[0]
        results["verdict"] = f"BOTH CONSTRAINTS UNBLOCKED by {top['variant']}: {top['n_pegrna_in_3to10nt_window']} optimal pegRNAs + {top['n_pe3b_nickers']} PE3b nickers"
    else:
        pegrna_only = [s for s in summary_rows if s["n_pegrna_in_3to10nt_window"] > 0]
        pe3b_only = [s for s in summary_rows if s["n_pe3b_nickers"] > 0]
        results["verdict"] = (
            f"Partial: {len(pegrna_only)} variant(s) unblock pegRNA geometry; "
            f"{len(pe3b_only)} unblock PE3b; none unblock both"
        )
    print(f"\nVerdict: {results['verdict']}")

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
