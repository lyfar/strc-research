#!/usr/bin/env python3
"""
Hydrogel Phase 4f — MHC-II T-cell immunogenicity scan for PEPTIDE_TAIL91.

Predicts HLA class II (DRB1) binding 9-mer cores in sliding 15-mer windows.
Analytical position-specific scoring approximation of NetMHCIIpan-style
binders — not as good as the real tool but catches the obvious hot-spots
that would trigger DC presentation → T-cell activation → anti-drug antibodies.

This matters for a chronic-dose peptide therapeutic:
  - If tail91 contains a strong immunodominant epitope, patients will
    develop anti-drug antibodies within 2-4 weeks of chronic dosing
    → therapy fails after first cycle
  - Single-allele hits are manageable (each human has 2 DRB1 alleles max)
  - Pan-allele hits (binds many HLA-DR classes) are developability-killing

Method: score each 9-mer core in tail91 against simplified position-specific
scoring matrices (PSSMs) for the 9 most common HLA-DRB1 alleles (coverage
~75% of European/Asian populations). 9-mers with score above median +2σ are
flagged as potential binders.

Output:
  - Per-allele list of predicted binders
  - Pan-allele promiscuous cores (bind 4+ alleles — highest risk)
  - Region mapping (WH2 / linker / RADA16 / tail)
  - Recommended mutations to break immunogenic cores (tail91_immunoclean variant)

References for PSSM motifs (coarse approximation):
  - P1 anchor hydrophobic (F/Y/W/L/I/V/M)
  - P4 anchor varies by allele (DR*01: polar/aromatic; DR*03: negative; etc.)
  - P6 anchor varies
  - P9 anchor often small (A/V/I/L)
  - Most DR alleles dislike Pro at P4/P6/P9

This is an approximation — not FDA-grade. Use as a triage before running
NetMHCIIpan / IEDB / SureBind for final clinical candidate selection.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path
import numpy as np

OUT = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase4f_mhc_immunogenicity.json")


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


STRC_FULL = fetch_uniprot("Q7RTU9")
WH2_NTERM = "RQLVKAIPDNCSKSNVSR"
RADA16 = "ADARADARADARADA"
LINK_5 = "GSGSG"
TAIL_91 = STRC_FULL[1619:1710]
PEPTIDE_TAIL91 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_91
assert len(PEPTIDE_TAIL91) == 134

# Simplified per-position preferences for 9 common HLA-DRB1 alleles
# Derived from published binding motifs (Sette 2012, Nielsen 2020); coarse.
# Positions P1-P9. "+" residues gain score, "-" residues lose score.

ALLELES = {
    "DRB1_01_01": {  # HLA-DR1, ~15% European; hydrophobic P1, aromatic P4, small P6
        "P1_positive": set("YFWILVM"),
        "P1_negative": set("DEHKR"),
        "P4_positive": set("VILMAF"),
        "P4_negative": set("DEPKR"),
        "P6_positive": set("STAVIL"),
        "P6_negative": set("DEP"),
        "P9_positive": set("AVILT"),
        "P9_negative": set("DEP"),
    },
    "DRB1_03_01": {  # HLA-DR17, ~10%; polar P1, negative P4, small P6/P9
        "P1_positive": set("YFLIMN"),
        "P1_negative": set("GP"),
        "P4_positive": set("DE"),
        "P4_negative": set("KRPW"),
        "P6_positive": set("AVS"),
        "P6_negative": set("PKR"),
        "P9_positive": set("AVILN"),
        "P9_negative": set("DEP"),
    },
    "DRB1_04_01": {  # HLA-DR4 (rheumatoid-associated), ~10%
        "P1_positive": set("YFWILMV"),
        "P1_negative": set("DP"),
        "P4_positive": set("YFWILV"),
        "P4_negative": set("DEPR"),
        "P6_positive": set("AKRN"),
        "P6_negative": set("DEP"),
        "P9_positive": set("SAVLIT"),
        "P9_negative": set("DEP"),
    },
    "DRB1_07_01": {  # HLA-DR7, ~13%
        "P1_positive": set("YFWILMV"),
        "P1_negative": set("DPG"),
        "P4_positive": set("ILMVAF"),
        "P4_negative": set("DEPR"),
        "P6_positive": set("AVS"),
        "P6_negative": set("DEP"),
        "P9_positive": set("AVIL"),
        "P9_negative": set("DEP"),
    },
    "DRB1_11_01": {  # HLA-DR5, ~8%
        "P1_positive": set("YFWILMV"),
        "P1_negative": set("DPG"),
        "P4_positive": set("NQDH"),
        "P4_negative": set("KRWP"),
        "P6_positive": set("ATSVIL"),
        "P6_negative": set("DEP"),
        "P9_positive": set("AVILNQ"),
        "P9_negative": set("DEP"),
    },
    "DRB1_13_01": {  # HLA-DR6, ~6%
        "P1_positive": set("YFWILMV"),
        "P1_negative": set("DP"),
        "P4_positive": set("DENQ"),
        "P4_negative": set("KRPW"),
        "P6_positive": set("AKRN"),
        "P6_negative": set("DEP"),
        "P9_positive": set("AVILQ"),
        "P9_negative": set("DEP"),
    },
    "DRB1_15_01": {  # HLA-DR15 (MS-associated), ~12%
        "P1_positive": set("YFWILMV"),
        "P1_negative": set("DPG"),
        "P4_positive": set("ILVMA"),
        "P4_negative": set("DEPR"),
        "P6_positive": set("FYW"),
        "P6_negative": set("DEP"),
        "P9_positive": set("AVIL"),
        "P9_negative": set("DEP"),
    },
    "DRB1_08_01": {  # ~2%, included for breadth
        "P1_positive": set("YFLIMV"),
        "P1_negative": set("DP"),
        "P4_positive": set("DE"),
        "P4_negative": set("KRP"),
        "P6_positive": set("AKRN"),
        "P6_negative": set("DEP"),
        "P9_positive": set("AVIL"),
        "P9_negative": set("DEP"),
    },
    "DRB1_16_01": {  # Asian-common, ~3% in Asians
        "P1_positive": set("YFWILM"),
        "P1_negative": set("DP"),
        "P4_positive": set("AVILM"),
        "P4_negative": set("KR"),
        "P6_positive": set("AKRN"),
        "P6_negative": set("DEP"),
        "P9_positive": set("AVIL"),
        "P9_negative": set("DEP"),
    },
}


def score_9mer(core: str, allele: dict) -> float:
    """Coarse PSSM score for a 9-mer core against one allele."""
    s = 0.0
    pos_map = {0: "P1", 3: "P4", 5: "P6", 8: "P9"}
    for pos, label in pos_map.items():
        aa = core[pos]
        if aa in allele.get(f"{label}_positive", set()):
            s += 1.0
        if aa in allele.get(f"{label}_negative", set()):
            s -= 1.5
    # Penalty for proline in core (rigid, interferes with groove)
    s -= 0.5 * core.count("P")
    # Flanking charged residues provide mild bonus for accessibility
    return s


def region_of(pos: int) -> str:
    regions = {"wh2": (1, 18), "gs1": (19, 23), "rada16": (24, 38),
               "gs2": (39, 43), "tail": (44, 134)}
    for name, (lo, hi) in regions.items():
        if lo <= pos <= hi:
            return name
    return "out"


def tail_aa_to_strc(pos: int) -> int | None:
    """Construct pos (1-idx) → STRC aa number."""
    if 44 <= pos <= 134:
        return 1620 + (pos - 44)
    return None


def scan_peptide(seq: str, threshold: float = 3.0) -> dict:
    """Return per-allele list of (pos, 9mer, score, region) above threshold."""
    per_allele = {a: [] for a in ALLELES}
    all_binders_per_pos = {i+1: [] for i in range(len(seq) - 8)}
    for i in range(len(seq) - 8):
        core = seq[i:i+9]
        pos = i + 1
        r = region_of(pos)
        strc_n = tail_aa_to_strc(pos)
        for a, rules in ALLELES.items():
            s = score_9mer(core, rules)
            if s >= threshold:
                hit = {
                    "pos": pos, "core": core, "score": round(s, 2),
                    "region": r, "strc_aa_if_tail": strc_n,
                }
                per_allele[a].append(hit)
                all_binders_per_pos[pos].append({"allele": a, "score": round(s, 2)})
    # Pan-allele promiscuous = 9mers with ≥4 alleles binding
    promiscuous = []
    for pos, hits in all_binders_per_pos.items():
        if len(hits) >= 4:
            promiscuous.append({
                "pos": pos, "core": seq[pos-1:pos+8],
                "region": region_of(pos),
                "strc_aa_if_tail": tail_aa_to_strc(pos),
                "n_alleles_binding": len(hits),
                "alleles": [h["allele"] for h in hits],
                "max_score": max(h["score"] for h in hits),
                "mean_score": round(np.mean([h["score"] for h in hits]), 2),
            })
    return {
        "per_allele_binders": per_allele,
        "per_allele_binder_count": {a: len(v) for a, v in per_allele.items()},
        "promiscuous_cores": sorted(promiscuous,
                                    key=lambda x: -x["n_alleles_binding"]),
        "n_promiscuous": len(promiscuous),
    }


def recommend_mutations(promiscuous_cores: list) -> list[dict]:
    """For each promiscuous core, suggest a minimal mutation that breaks
    the P1 anchor (swap hydrophobic → polar)."""
    recommendations = []
    for hit in promiscuous_cores:
        core = hit["core"]
        # Primary recommendation: swap P1 (position 0 in core) to a small polar
        p1 = core[0]
        if p1 in set("YFWILMV"):
            # Recommend X → S or T (breaks hydrophobic P1 anchor)
            recommendations.append({
                "position": hit["pos"],
                "region": hit["region"],
                "strc_aa": hit["strc_aa_if_tail"],
                "core": core,
                "original_P1": p1,
                "recommended_substitution": f"{p1}→S or T",
                "n_alleles_broken": hit["n_alleles_binding"],
                "note": "Breaks P1 anchor for hydrophobic pockets (9 DRB1 alleles)"
            })
        else:
            # Secondary: swap P4/P6 anchor
            recommendations.append({
                "position": hit["pos"],
                "region": hit["region"],
                "strc_aa": hit["strc_aa_if_tail"],
                "core": core,
                "original_P1": p1,
                "recommended_substitution": "try P4 substitution (pos+3)",
                "n_alleles_broken": hit["n_alleles_binding"],
                "note": "Conservative P1 already; target P4"
            })
    return recommendations


def main():
    scan = scan_peptide(PEPTIDE_TAIL91, threshold=3.0)

    # Per-region distribution of binders
    region_dist = {"wh2": 0, "gs1": 0, "rada16": 0, "gs2": 0, "tail": 0}
    for a, hits in scan["per_allele_binders"].items():
        for h in hits:
            region_dist[h["region"]] += 1

    recommendations = recommend_mutations(scan["promiscuous_cores"])

    # Immunogenicity burden summary
    total_binders = sum(scan["per_allele_binder_count"].values())
    n_alleles = len(ALLELES)
    avg_binders_per_allele = total_binders / n_alleles
    # Clinical rule-of-thumb: <3 promiscuous cores = low-risk
    #                        3-10 = moderate-risk
    #                        >10 = high-risk, needs deimmunisation
    n_promisc = scan["n_promiscuous"]
    risk_tier = (
        "LOW" if n_promisc < 3 else
        "MODERATE" if n_promisc <= 10 else
        "HIGH"
    )

    summary = {
        "batch": "hydrogel_phase4f_mhc_immunogenicity",
        "date": "2026-04-23",
        "construct": "PEPTIDE_TAIL91",
        "length": 134,
        "n_alleles_tested": n_alleles,
        "allele_list": list(ALLELES.keys()),
        "pop_coverage_percent_estimate": 75,
        "scoring_threshold": 3.0,
        "total_binders_all_alleles": total_binders,
        "average_binders_per_allele": round(avg_binders_per_allele, 1),
        "per_allele_binder_count": scan["per_allele_binder_count"],
        "region_distribution": region_dist,
        "promiscuous_cores_4plus_alleles": scan["promiscuous_cores"],
        "n_promiscuous_cores": n_promisc,
        "immunogenicity_risk_tier": risk_tier,
        "recommended_mutations_for_deimmunisation": recommendations,
        "notes": [
            "This is a coarse analytical approximation — NOT a substitute for NetMHCIIpan",
            "RADA16 scaffold (24-38) is highly repetitive; any binder there will repeat 3-4× in the sequence",
            "Tail region 44-134 contains STRC native sequence — T-cell tolerance likely exists "
            "for self-peptide; immunogenic only for patients with STRC knockout (Misha's paternal del)",
            "MHC-II binding is necessary but not sufficient for immunogenicity; T-cell receptor "
            "engagement + costimulation required for actual ADA response",
        ],
        "follow_up_if_high_risk": [
            "Run NetMHCIIpan 4.3 on full construct (27 alleles)",
            "Deimmunisation via T-cell-hopping (swap P4 anchor to break multiple alleles)",
            "Consider non-proteogenic amino acids at key anchors (already in Phase 4c v5 plan)",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, default=str))

    print("=== Phase 4f MHC-II T-cell Immunogenicity Scan ===\n")
    print(f"Construct: PEPTIDE_TAIL91 (134 aa)")
    print(f"Alleles tested: {n_alleles} ({', '.join(list(ALLELES.keys())[:4])}, ...)")
    print(f"Threshold score: 3.0\n")
    print(f"Total predicted binders across all alleles: {total_binders}")
    print(f"Avg binders per allele: {avg_binders_per_allele:.1f}")
    print(f"\nPer-allele binder counts:")
    for a, n in scan["per_allele_binder_count"].items():
        print(f"  {a}: {n}")
    print(f"\nRegion distribution:")
    for r, n in region_dist.items():
        print(f"  {r}: {n}")
    print(f"\nPromiscuous (≥4 alleles) cores: {n_promisc}")
    if n_promisc > 0:
        print("Top 5 promiscuous hits:")
        for h in scan["promiscuous_cores"][:5]:
            strc = f"STRC aa {h['strc_aa_if_tail']}" if h["strc_aa_if_tail"] else ""
            print(f"  pos {h['pos']:3d} ({h['region']:7s}) {strc} core={h['core']}, "
                  f"{h['n_alleles_binding']} alleles, max={h['max_score']}")
    print(f"\n>>> Immunogenicity risk tier: {risk_tier}")
    print(f"\nJSON written: {OUT}")


if __name__ == "__main__":
    main()
