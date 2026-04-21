"""
PE Phase 3 — PE3b allele discrimination audit.

For every PE3b-spanning nicker surfaced by Phase 2 (across all Cas9 variants),
verify that the nicker's sgRNA (= protospacer sequence 5'→3') matches the
EDITED genome and MISMATCHES the mutant genome at the variant position.
Classify the mismatch position within the protospacer (seed/mid/distal) to
score discrimination strength.

Geometry recap for Misha's c.4976A>C (p.E1659A):
- Variant on + strand (chr15:43600551): WT = T, MUT = G, corrected/edited = T
- Variant on - strand (complement):      WT = A, MUT = C, corrected/edited = A
- pegRNA spacer targets + strand; H840A nickase nicks the + strand (protospacer/
  non-target strand); RT installs T at chr15:43600551 on + strand.
- PE3b nicker should nick the OPPOSITE strand (- strand) after the edit, so
  MMR uses the edited + strand as template.
- PE3b discrimination: protospacer letter at the variant position must equal
  the EDITED base on the protospacer strand (A on - strand; T on + strand).

Seed/mid/distal convention for SpCas9:
- Protospacer read 5'→3', 20 nt, PAM is 3' of protospacer.
- SEED = positions 13-20 (PAM-proximal); mismatches here strongly block Cas9.
- MID  = positions 8-12.
- DISTAL = positions 1-7 (PAM-distal); mismatches here are tolerated.
PE3b discrimination grades:
- STRONG  if variant falls in positions 13-20
- MODERATE if positions 8-12
- WEAK    if positions 1-7
Reference: Anzalone 2019; Doench 2016 (mismatch tolerance by position).
"""

import json
from pathlib import Path

PHASE2 = Path("/Users/egorlyfar/Brain/research/strc/models/pe_phase2_pam_expansion.json")
OUT = Path("/Users/egorlyfar/Brain/research/strc/models/pe_phase3_allele_discrimination.json")

VARIANT_POS_PLUS = 43600551
PLUS_WT = "T"
PLUS_MUT = "G"
PLUS_EDITED = "T"  # correction restores WT
MINUS_WT = "A"
MINUS_MUT = "C"
MINUS_EDITED = "A"


def region_class(pos_1idx: int) -> tuple[str, str]:
    """Return (region, discrimination_grade) for a mismatch at `pos_1idx` in a 20-nt protospacer."""
    if pos_1idx >= 13:
        return "seed", "strong"
    if pos_1idx >= 8:
        return "mid", "moderate"
    return "distal", "weak"


def variant_position_in_protospacer(strand: str, plus_coords: list[int]) -> int | None:
    """1-indexed variant position within a 20-nt protospacer read 5'→3'."""
    start_plus, end_plus = plus_coords  # inclusive, always start<=end in + coords
    if not (start_plus <= VARIANT_POS_PLUS <= end_plus):
        return None
    if strand == "+":
        return (VARIANT_POS_PLUS - start_plus) + 1
    # - strand: 5' end = + coord end_plus; read 5'->3' decreasing + coord
    return (end_plus - VARIANT_POS_PLUS) + 1


def audit_candidate(cand: dict, pegrna_edited_strand: str) -> dict:
    strand = cand["strand"]
    proto = cand["protospacer_minus_5to3"] if strand == "-" else cand["protospacer_plus_5to3"]
    plus_coords = cand["protospacer_coords_plus"]
    pos = variant_position_in_protospacer(strand, plus_coords)
    if pos is None:
        return {**cand, "audit_verdict": "FAILED_SPAN", "reason": "variant not within protospacer"}

    letter_at_variant = proto[pos - 1].upper()

    # What letter the EDITED genome would have on this strand at the variant position:
    if strand == "+":
        expected_edited = PLUS_EDITED
        expected_mut = PLUS_MUT
    else:
        expected_edited = MINUS_EDITED
        expected_mut = MINUS_MUT

    matches_edited = letter_at_variant == expected_edited
    mismatches_mut = letter_at_variant != expected_mut
    allele_discriminating = matches_edited and mismatches_mut

    region, grade = region_class(pos)

    # Utility: PE3b should nick the OPPOSITE strand from the pegRNA-edited strand.
    if strand == ("-" if pegrna_edited_strand == "+" else "+"):
        strand_utility = "USEFUL"
    else:
        strand_utility = "REVERSE_PE3B"  # would nick the edited strand, counterproductive

    return {
        **cand,
        "variant_pos_in_protospacer_1idx": pos,
        "protospacer_letter_at_variant": letter_at_variant,
        "expected_letter_edited_strand": expected_edited,
        "expected_letter_mut_strand": expected_mut,
        "matches_edited_allele": matches_edited,
        "mismatches_mut_allele": mismatches_mut,
        "allele_discriminating": allele_discriminating,
        "mismatch_region": region,
        "discrimination_grade": grade,
        "strand_utility": strand_utility,
    }


def rank_candidates(audited: list[dict]) -> list[dict]:
    """Rank PE3b candidates primarily by discrimination (seed>mid>distal), then by distance."""
    grade_rank = {"strong": 0, "moderate": 1, "weak": 2}

    def key(c):
        return (
            0 if c.get("strand_utility") == "USEFUL" else 1,
            0 if c.get("allele_discriminating") else 1,
            grade_rank.get(c.get("discrimination_grade"), 9),
            abs(c.get("nick_to_edit_nt", 999)),
        )

    return sorted(audited, key=key)


def rank_conservative(audited: list[dict]) -> list[dict]:
    """
    Alternative ranking that weights safety: prefer nick-to-edit 10–80 nt
    (avoids concurrent-DSB risk from very-close PE3b nicks, and keeps efficiency
    reasonable). Among distance-safe candidates, higher discrimination wins.
    """
    grade_rank = {"strong": 0, "moderate": 1, "weak": 2}

    def distance_bucket(d: float) -> int:
        if 10 <= d <= 80:
            return 0  # safe
        if d < 10:
            return 1  # too close — concurrent DSB risk
        return 2  # too far

    def key(c):
        return (
            0 if c.get("strand_utility") == "USEFUL" else 1,
            0 if c.get("allele_discriminating") else 1,
            distance_bucket(c.get("nick_to_edit_nt", 999)),
            grade_rank.get(c.get("discrimination_grade"), 9),
            abs(c.get("nick_to_edit_nt", 999)),
        )

    return sorted(audited, key=key)


def main():
    data = json.loads(PHASE2.read_text())
    variants = data["variants"]
    pegrna_edited_strand = "+"  # pegRNA designs across Phase 2 all target + strand

    audit_by_variant = {}
    for var_name, vrec in variants.items():
        audited = []
        for c in vrec.get("pe3b_nickers_minus_strand", []):
            audited.append(audit_candidate(c, pegrna_edited_strand))
        for c in vrec.get("pe3b_nickers_plus_strand", []):
            audited.append(audit_candidate(c, pegrna_edited_strand))
        ranked_aggressive = rank_candidates(audited)
        ranked_conservative = rank_conservative(audited)
        audit_by_variant[var_name] = {
            "n_pe3b_total": len(ranked_aggressive),
            "n_useful": sum(1 for c in ranked_aggressive if c["strand_utility"] == "USEFUL"),
            "n_allele_discriminating": sum(1 for c in ranked_aggressive if c.get("allele_discriminating")),
            "n_strong_discrimination": sum(1 for c in ranked_aggressive if c.get("discrimination_grade") == "strong"),
            "n_moderate_discrimination": sum(1 for c in ranked_aggressive if c.get("discrimination_grade") == "moderate"),
            "n_weak_discrimination": sum(1 for c in ranked_aggressive if c.get("discrimination_grade") == "weak"),
            "top_candidate_aggressive": ranked_aggressive[0] if ranked_aggressive else None,
            "top_candidate_conservative": ranked_conservative[0] if ranked_conservative else None,
            "all_ranked": ranked_aggressive,
        }

    # Build condensed summary across variants.
    summary = {}
    for var_name, rec in audit_by_variant.items():
        agg = rec["top_candidate_aggressive"]
        con = rec["top_candidate_conservative"]

        def unpack(c):
            if not c:
                return {k: None for k in ["protospacer", "grade", "pos", "nick_to_edit", "discriminating"]}
            return {
                "protospacer": c.get("protospacer_minus_5to3") or c.get("protospacer_plus_5to3"),
                "grade": c.get("discrimination_grade"),
                "pos": c.get("variant_pos_in_protospacer_1idx"),
                "nick_to_edit": c.get("nick_to_edit_nt"),
                "discriminating": c.get("allele_discriminating"),
            }

        summary[var_name] = {
            "n_useful_discriminating": sum(
                1
                for c in rec["all_ranked"]
                if c["strand_utility"] == "USEFUL" and c.get("allele_discriminating")
            ),
            "n_strong_seed": sum(
                1
                for c in rec["all_ranked"]
                if c["strand_utility"] == "USEFUL"
                and c.get("allele_discriminating")
                and c.get("discrimination_grade") == "strong"
            ),
            "aggressive_top": unpack(agg),
            "conservative_top": unpack(con),
        }

    out = {
        "meta": {
            "phase": "3",
            "title": "PE3b allele-discrimination audit for E1659A",
            "variant_plus_coord": VARIANT_POS_PLUS,
            "edited_plus_base": PLUS_EDITED,
            "mut_plus_base": PLUS_MUT,
            "edited_minus_base": MINUS_EDITED,
            "mut_minus_base": MINUS_MUT,
            "pegrna_edited_strand": pegrna_edited_strand,
            "method": "Score each PE3b-spanning candidate: letter at variant position, region classification (seed/mid/distal), utility (nicks opposite strand?).",
        },
        "summary": summary,
        "audit_by_variant": audit_by_variant,
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(f"Wrote: {OUT}")
    for v, s in summary.items():
        agg = s["aggressive_top"]
        con = s["conservative_top"]
        print(
            f"  {v}: useful_disc={s['n_useful_discriminating']}, strong_seed={s['n_strong_seed']}"
        )
        print(
            f"    aggressive: {agg['protospacer']} (pos {agg['pos']}, {agg['grade']}, {agg['nick_to_edit']} nt)"
        )
        print(
            f"    conservative: {con['protospacer']} (pos {con['pos']}, {con['grade']}, {con['nick_to_edit']} nt)"
        )


if __name__ == "__main__":
    main()
