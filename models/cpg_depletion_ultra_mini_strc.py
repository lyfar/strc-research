#!/usr/bin/env python3
"""Ultra-Mini STRC (aa 1075-1775) CpG depletion for AAV vector de-immunization.

Follow-on to [[STRC CpG Depletion Mini-STRC]] (which depleted aa 700-1775 CDS to
0 CpG at 3.5% CAI cost). After [[STRC Mini-STRC Truncation Interface Validation]]
proved Ultra-Mini (aa 1075-1775) preserves the TMEM145 binding pocket at sub-A
RMSD, the next clinical prerequisite is a CpG-free CDS for the shorter
construct. Expected CpG count scales with CDS length (~0.65x), and CAI cost
should also drop because there is less sequence space for unavoidable CpGs.

Reuses logic from cpg_depletion_mini_strc.py but runs on the 701-aa Ultra-Mini
protein. Protein sequence is extracted from the existing verified FASTA (which
was itself verified vs job-h-strc-cterm-only.cif chain A extraction).

Deterministic: seed 42. Replication: /opt/miniconda3/bin/python3 cpg_depletion_ultra_mini_strc.py
"""
import json
import random
from pathlib import Path

# Reuse depletion logic from sibling script
from cpg_depletion_mini_strc import (
    HUMAN_CODON_FREQ,
    CODON_TABLE,
    AA_TO_MAX,
    codon_optimize_max,
    count_cpg,
    cai,
    deplete_cpg,
    gc_percent,
)

random.seed(42)
OUT_DIR = Path(__file__).resolve().parent
OUT_JSON = OUT_DIR / "cpg_depletion_ultra_mini_strc.json"
OUT_FASTA_V0 = OUT_DIR / "cpg_depletion_ultra_mini_strc_v0.fasta"
OUT_FASTA_MAX = OUT_DIR / "cpg_depletion_ultra_mini_strc_max.fasta"

ULTRA_START, ULTRA_END = 1075, 1775
# Verified Ultra-Mini STRC protein sequence (701 aa, sha256=d1e49a41a686...)
# Extracted from two independent sources: (1) mini_STRC_700_1775 FASTA slice
# [375:1076]; (2) MMCIFParser extraction from job-h-strc-cterm-only.cif chain A.
# Both matched. See af3_jobs_2026-04-21/MANIFEST.json for provenance.
ULTRA_MINI_PROT = (
    "LVGACSCLAPELSRLSACQTAALLQTFRVKDGVKNMGTTGAGPAVCIPGQPIPTTWPDCLLPLLPLKLLQLDSL"
    "ALLANRRRYWELPWSEQQAQFLWKKMQVPTNLTLRNLQALGTLAGGMSCEFLQQINSMVDFLEVVHMIYQLPTR"
    "VRGSLRACIWAELQRRMAMPEPEWTTVGPELNGLDSKLLLDLPIQLMDRLSNESIMLVVELVQRAPEQLLALTP"
    "LHQAALAERALQNLAPKETPVSGEVLETLGPLVGFLGTESTRQIPLQILLSHLSQLQGFCLGETFATELGWLLL"
    "QESVLGKPELWSQDEVEQAGRLVFTLSTEAISLIPREALGPETLERLLEKQQSWEQSRVGQLCREPQLAAKKAA"
    "LVAGVVRPAAEDLPEPVPNCADVRGTFPAAWSATQIAEMELSDFEDCLTLFAGDPGLGPEELRAAMGKAKQLWG"
    "PPRGFRPEQILQLGRLLIGLGDRELQELILVDWGVLSTLGQIDGWSTTQLRIVVSSFLRQSGRHVSHLDFVHLT"
    "ALGYTLCGLRPEELQHISSWEFSQAALFLGTLHLQCSEEQLEVLAHLLVLPGGFGPISNWGPEIFTEIGTIAAG"
    "IPDLALSALLRGQIQGVTPLAISVIPPPKFAVVFSPIQLSSLTSAQAVAVTPEQMAFLSPEQRRAVAWAQHEGK"
    "ESPEQQGRSTAWGLQDWSRPSWSLVLTISFLGHLL"
)
assert len(ULTRA_MINI_PROT) == 701, f"expected 701 aa, got {len(ULTRA_MINI_PROT)}"

HUMAN_GENOME_CPG_PER_KB = 9.7


def main():
    print(f"Ultra-Mini STRC: aa {ULTRA_START}-{ULTRA_END} ({len(ULTRA_MINI_PROT)} aa)")

    print("\nBaseline codon optimization (max-freq synonym per aa) ...")
    v0 = codon_optimize_max(ULTRA_MINI_PROT)
    v0_cpg = count_cpg(v0)
    v0_gc = gc_percent(v0)
    v0_cai = cai(v0)
    v0_len = len(v0)
    v0_density = 1000 * v0_cpg / v0_len
    fold_vs_genome = v0_density / HUMAN_GENOME_CPG_PER_KB
    print(f"  length: {v0_len} bp ({v0_len // 3} codons incl stop)")
    print(f"  CpG count: {v0_cpg} ({v0_density:.2f}/kb, {fold_vs_genome:.2f}x human genome avg)")
    print(f"  GC%: {v0_gc}")
    print(f"  CAI: {v0_cai}")

    print("\nCpG depletion sweep ...")
    thresholds = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60]
    sweep = []
    for t in thresholds:
        depleted, swaps = deplete_cpg(v0, t)
        cpg = count_cpg(depleted)
        cai_val = cai(depleted)
        gc = gc_percent(depleted)
        reduction = 100 * (v0_cpg - cpg) / v0_cpg if v0_cpg else 0
        cai_delta_pct = 100 * (v0_cai - cai_val) / v0_cai
        sweep.append({
            "threshold": t,
            "swaps": swaps,
            "cpg_count": cpg,
            "cpg_per_kb": round(1000 * cpg / len(depleted), 2),
            "cpg_reduction_pct": round(reduction, 2),
            "cai": cai_val,
            "cai_delta_pct": round(cai_delta_pct, 2),
            "gc_pct": gc,
            "length_bp": len(depleted),
        })
        print(f"  threshold={t:.2f}  swaps={swaps:4d}  CpG={cpg:4d} (-{reduction:5.1f}%)  CAI={cai_val}  (-{cai_delta_pct:4.1f}%)  GC={gc}%")

    # Best = max CpG reduction with CAI cost <= 5%
    best = None
    for s in sweep:
        if s["cai_delta_pct"] <= 5.0:
            if best is None or s["cpg_reduction_pct"] > best["cpg_reduction_pct"]:
                best = s
    if best is None:
        best = sweep[0]
    print(f"\nBest (<=5% CAI cost): threshold={best['threshold']:.2f}  CpG-{best['cpg_reduction_pct']}%  CAI-{best['cai_delta_pct']}%")

    def write_fasta(path, header, seq):
        with open(path, "w") as f:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    max_depleted, _ = deplete_cpg(v0, thresholds[-1])
    write_fasta(
        OUT_FASTA_V0,
        f"ultra_mini_STRC_aa{ULTRA_START}-{ULTRA_END}_v0_max_CAI len={len(v0)} cpg={v0_cpg} cai={v0_cai}",
        v0,
    )
    write_fasta(
        OUT_FASTA_MAX,
        f"ultra_mini_STRC_aa{ULTRA_START}-{ULTRA_END}_CpG_depleted_t0.60 cpg={count_cpg(max_depleted)} cai={cai(max_depleted)}",
        max_depleted,
    )

    # Translation round-trip
    codon_to_aa = {c: aa for aa, cs in CODON_TABLE.items() for c in cs}
    def translate(cds):
        aas = []
        for i in range(0, len(cds) - 2, 3):
            c = cds[i:i + 3]
            if c in ("TAA", "TAG", "TGA"):
                break
            aas.append(codon_to_aa.get(c, "X"))
        return "".join(aas)

    v0_prot = translate(v0)
    max_prot = translate(max_depleted)
    prot_ok_v0 = v0_prot == ULTRA_MINI_PROT
    prot_ok_max = max_prot == ULTRA_MINI_PROT
    print(f"\nTranslation round-trip: v0 {'OK' if prot_ok_v0 else 'FAIL'}, max-depleted {'OK' if prot_ok_max else 'FAIL'}")

    # Comparison vs prior mini-STRC (700-1775) CpG depletion
    prior_path = OUT_DIR / "cpg_depletion_mini_strc.json"
    prior = json.loads(prior_path.read_text()) if prior_path.exists() else None

    comparison = None
    if prior:
        prior_v0 = prior["baseline_v0"]
        # Pick prior's final-depletion line (last threshold in sweep)
        prior_final = prior["depletion_sweep"][-1]
        comparison = {
            "mini_700_1775_baseline_cpg": prior_v0["cpg_count"],
            "mini_700_1775_baseline_bp": prior_v0["cds_length_bp"],
            "ultra_1075_1775_baseline_cpg": v0_cpg,
            "ultra_1075_1775_baseline_bp": v0_len,
            "cpg_reduction_due_to_shorter_cds": round(
                100 * (prior_v0["cpg_count"] - v0_cpg) / prior_v0["cpg_count"], 2
            ),
            "bp_reduction_due_to_shorter_cds": round(
                100 * (prior_v0["cds_length_bp"] - v0_len) / prior_v0["cds_length_bp"], 2
            ),
            "mini_700_1775_final_cpg": prior_final["cpg_count"],
            "mini_700_1775_final_cai_cost_pct": prior_final["cai_delta_pct"],
            "ultra_1075_1775_final_cpg": sweep[-1]["cpg_count"],
            "ultra_1075_1775_final_cai_cost_pct": sweep[-1]["cai_delta_pct"],
        }
        print("\nComparison vs mini-STRC (700-1775):")
        for k, v in comparison.items():
            print(f"  {k}: {v}")

    result = {
        "meta": {
            "date": "2026-04-21",
            "construct": "ultra_mini_STRC_aa1075-1775",
            "protein_length_aa": len(ULTRA_MINI_PROT),
            "parent_note": "STRC Mini-STRC Truncation Interface Validation",
            "prior_clinical_construct": "mini_STRC_aa700-1775 (cpg_depletion_mini_strc.json)",
            "codon_usage_source": "Kazusa Homo sapiens",
            "cai_method": "Sharp & Li 1987",
            "human_genome_cpg_per_kb_reference": HUMAN_GENOME_CPG_PER_KB,
        },
        "baseline_v0": {
            "cds_length_bp": v0_len,
            "cpg_count": v0_cpg,
            "cpg_per_kb": round(v0_density, 2),
            "cpg_fold_vs_human_genome": round(fold_vs_genome, 2),
            "cai": v0_cai,
            "gc_pct": v0_gc,
            "translation_roundtrip_ok": prot_ok_v0,
        },
        "depletion_sweep": sweep,
        "best_under_5pct_cai_cost": best,
        "max_depleted_roundtrip_ok": prot_ok_max,
        "comparison_vs_mini_700_1775": comparison,
        "outputs": {
            "v0_fasta": str(OUT_FASTA_V0),
            "max_depleted_fasta": str(OUT_FASTA_MAX),
            "results_json": str(OUT_JSON),
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {OUT_JSON}")
    print(f"Wrote {OUT_FASTA_V0}")
    print(f"Wrote {OUT_FASTA_MAX}")


if __name__ == "__main__":
    main()
