#!/usr/bin/env python3
"""
Mini-STRC CpG depletion for AAV vector de-immunization.

Context: TLR9-sensing of unmethylated CpG in AAV genomes is a major driver of
anti-capsid immunity and transgene silencing (Shao 2018, Faust 2013, Chan 2021).
The parent note [[STRC Mini-STRC Single-Vector Hypothesis]] flags:
"CpG depletion from current 156 (4.9x genome avg) — 87% reducible at 3% CAI
cost. **Must do before lab.**" This script backs that claim with numbers.

Pipeline:
  1. Fetch human STRC (UniProt Q7RTU9) FASTA.
  2. Extract mini-STRC protein window (residues 700..1775).
  3. Baseline: codon-optimize for Homo sapiens using the Kazusa codon-usage
     table (highest-frequency synonym per residue) -> "v0".
  4. Count CpG dinucleotides + CpG density (per kb) vs human genome average
     (~9.7 CpG/kb = 0.0097 fraction per nt, since CpG is suppressed genome-wide).
  5. Iterative depletion sweep: at each pass, scan for CpG dinucleotides and
     swap one codon of the offending pair with a synonym that breaks the CG
     dinucleotide, as long as CAI cost per swap <= threshold. Sweep thresholds
     {0.00, 0.02, 0.05, 0.10, 0.20} — produces the CpG-reduction vs CAI-cost
     curve.
  6. Output depleted FASTAs + JSON table + summary figure-data.

Deterministic: seed 42.

Replication:
    /opt/miniconda3/bin/python3 cpg_depletion_mini_strc.py
"""
import json
import random
import re
from io import StringIO
from pathlib import Path
from urllib.request import urlopen

from Bio import SeqIO

random.seed(42)
OUT_DIR = Path(__file__).resolve().parent
OUT_JSON = OUT_DIR / "cpg_depletion_mini_strc.json"
OUT_FASTA_V0 = OUT_DIR / "cpg_depletion_mini_strc_v0.fasta"
OUT_FASTA_MAX = OUT_DIR / "cpg_depletion_mini_strc_max.fasta"

STRC_ACC = "Q7RTU9"
MINI_START, MINI_END = 700, 1775

# Kazusa Homo sapiens codon usage (frequency per thousand codons).
# Source: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
# Per-codon frequencies; higher = more abundant tRNA pool.
HUMAN_CODON_FREQ = {
    # Phe
    "TTT": 17.6, "TTC": 20.3,
    # Leu
    "TTA": 7.7, "TTG": 12.9, "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
    # Ile
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5,
    # Met
    "ATG": 22.0,
    # Val
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
    # Ser
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4, "AGT": 12.1, "AGC": 19.5,
    # Pro
    "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    # Thr
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
    # Ala
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    # Tyr
    "TAT": 12.2, "TAC": 15.3,
    # Stop
    "TAA": 1.0, "TAG": 0.8, "TGA": 1.6,
    # His
    "CAT": 10.9, "CAC": 15.1,
    # Gln
    "CAA": 12.3, "CAG": 34.2,
    # Asn
    "AAT": 17.0, "AAC": 19.1,
    # Lys
    "AAA": 24.4, "AAG": 31.9,
    # Asp
    "GAT": 21.8, "GAC": 25.1,
    # Glu
    "GAA": 29.0, "GAG": 39.6,
    # Cys
    "TGT": 10.6, "TGC": 12.6,
    # Trp
    "TGG": 13.2,
    # Arg
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4, "AGA": 12.2, "AGG": 12.0,
    # Gly
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}

CODON_TABLE = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "N": ["AAT", "AAC"],
    "D": ["GAT", "GAC"],
    "C": ["TGT", "TGC"],
    "E": ["GAA", "GAG"],
    "Q": ["CAA", "CAG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "K": ["AAA", "AAG"],
    "M": ["ATG"],
    "F": ["TTT", "TTC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
}

AA_TO_MAX = {aa: max(codons, key=lambda c: HUMAN_CODON_FREQ[c]) for aa, codons in CODON_TABLE.items()}
AA_TO_MAXFREQ = {aa: HUMAN_CODON_FREQ[AA_TO_MAX[aa]] for aa in CODON_TABLE}


def fetch_protein(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urlopen(url) as r:
        return str(next(SeqIO.parse(StringIO(r.read().decode()), "fasta")).seq)


def codon_optimize_max(protein: str) -> str:
    """Baseline: highest-frequency synonym per residue. Max CAI = 1.0 by
    construction."""
    return "".join(AA_TO_MAX[aa] for aa in protein) + "TAA"  # add stop


def count_cpg(cds: str) -> int:
    # CpG = CG dinucleotide in DNA (uppercase)
    return len(re.findall("CG", cds))


def cai(cds: str) -> float:
    """Sharp 1987-style CAI: geometric mean of per-codon relative adaptiveness.
    Relative adaptiveness w_c = freq(c) / max_freq(aa(c))."""
    if len(cds) < 3:
        return 0.0
    codons = [cds[i : i + 3] for i in range(0, len(cds) - 2, 3)]
    codons = [c for c in codons if c in HUMAN_CODON_FREQ and c not in ("TAA", "TAG", "TGA")]
    if not codons:
        return 0.0
    import math
    # Reverse lookup codon -> aa
    codon_to_aa = {c: aa for aa, cs in CODON_TABLE.items() for c in cs}
    log_sum = 0.0
    for c in codons:
        aa = codon_to_aa[c]
        w = HUMAN_CODON_FREQ[c] / AA_TO_MAXFREQ[aa]
        log_sum += math.log(max(w, 1e-9))
    return round(math.exp(log_sum / len(codons)), 4)


def deplete_cpg(cds: str, cai_cost_per_swap_max: float) -> tuple[str, int]:
    """Scan CDS for CpG dinucleotides; try to eliminate each by swapping one of
    the two codons involved to a synonym that (a) still encodes the same aa
    and (b) incurs CAI penalty <= cai_cost_per_swap_max. Repeat until no more
    swaps improve or all residual CpGs lack a cheap synonym.

    Returns (depleted_cds, swaps_made)."""
    cds_list = list(cds)
    codon_to_aa = {c: aa for aa, cs in CODON_TABLE.items() for c in cs}
    swaps = 0

    def codon_index(nt_pos):
        return nt_pos - (nt_pos % 3)

    def try_swap(c_start: int) -> bool:
        """Try to swap codon starting at c_start to eliminate any CpG it
        participates in (internal or at boundary). Return True if swapped."""
        codon = "".join(cds_list[c_start : c_start + 3])
        if codon not in codon_to_aa:
            return False
        aa = codon_to_aa[codon]
        synonyms = CODON_TABLE[aa]
        original_w = HUMAN_CODON_FREQ[codon] / AA_TO_MAXFREQ[aa]
        # We want swap that:
        #  - removes any internal CpG in the codon
        #  - breaks any boundary CpG (codon starts with G after preceding C;
        #    codon ends with C before following G)
        prev_nt = cds_list[c_start - 1] if c_start > 0 else ""
        next_nt = cds_list[c_start + 3] if c_start + 3 < len(cds_list) else ""

        def cpg_count_with(codon_seq: str) -> int:
            # count CpGs in: prev_nt + codon + next_nt triple context
            ctx = prev_nt + codon_seq + next_nt
            return ctx.count("CG")

        baseline_cpgs = cpg_count_with(codon)
        if baseline_cpgs == 0:
            return False

        best_syn = None
        best_cost = 10
        for syn in synonyms:
            if syn == codon:
                continue
            new_cpgs = cpg_count_with(syn)
            if new_cpgs >= baseline_cpgs:
                continue
            new_w = HUMAN_CODON_FREQ[syn] / AA_TO_MAXFREQ[aa]
            cost = max(0, original_w - new_w)
            if cost > cai_cost_per_swap_max:
                continue
            # Prefer max CpG reduction, then smallest cost
            score = (baseline_cpgs - new_cpgs) * 100 - cost
            if score > -best_cost:
                best_syn = syn
                best_cost = -score
        if best_syn is None:
            return False
        for i, ch in enumerate(best_syn):
            cds_list[c_start + i] = ch
        return True

    # Iterate: keep sweeping until a full pass makes no swaps
    while True:
        swaps_this_pass = 0
        for i in range(0, len(cds_list) - 2, 3):
            if try_swap(i):
                swaps_this_pass += 1
        swaps += swaps_this_pass
        if swaps_this_pass == 0:
            break
    return "".join(cds_list), swaps


def gc_percent(cds: str) -> float:
    cds = cds.upper()
    return round(100 * (cds.count("G") + cds.count("C")) / len(cds), 2)


def main():
    print(f"Fetching STRC ({STRC_ACC}) ...")
    strc = fetch_protein(STRC_ACC)
    print(f"  full STRC length: {len(strc)} aa")

    mini = strc[MINI_START - 1 : MINI_END]
    print(f"  mini-STRC window: aa {MINI_START}-{MINI_END} ({len(mini)} aa)")

    print("\nBaseline codon optimization (max-freq synonym per aa) ...")
    v0 = codon_optimize_max(mini)
    v0_cpg = count_cpg(v0)
    v0_gc = gc_percent(v0)
    v0_cai = cai(v0)
    v0_len = len(v0)
    v0_density = 1000 * v0_cpg / v0_len
    print(f"  length: {v0_len} bp ({v0_len // 3} codons incl stop)")
    print(f"  CpG count: {v0_cpg} ({v0_density:.2f}/kb)")
    print(f"  GC%: {v0_gc}")
    print(f"  CAI: {v0_cai}")

    # Human genome CpG density ~9.7/kb (genome-wide avg; CpG is suppressed)
    # Housekeeping CDS average ~21-23/kb.
    HUMAN_GENOME_CPG_PER_KB = 9.7
    fold_vs_genome = v0_density / HUMAN_GENOME_CPG_PER_KB
    print(f"  CpG density fold vs human-genome avg ({HUMAN_GENOME_CPG_PER_KB}/kb): {fold_vs_genome:.2f}x")

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

    # Pick "best" = max CpG reduction with CAI cost <= 5%
    best = None
    for s in sweep:
        if s["cai_delta_pct"] <= 5.0:
            if best is None or s["cpg_reduction_pct"] > best["cpg_reduction_pct"]:
                best = s
    if best is None:
        best = sweep[0]
    print(f"\nBest (<=5% CAI cost): threshold={best['threshold']:.2f}  CpG-{best['cpg_reduction_pct']}%  CAI-{best['cai_delta_pct']}%")

    # Dump FASTAs for v0 and max-aggressive
    def write_fasta(path, header, seq):
        with open(path, "w") as f:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")

    max_depleted, _ = deplete_cpg(v0, thresholds[-1])
    write_fasta(OUT_FASTA_V0, f"mini_STRC_aa{MINI_START}-{MINI_END}_v0_max_CAI len={len(v0)} cpg={v0_cpg} cai={v0_cai}", v0)
    write_fasta(
        OUT_FASTA_MAX,
        f"mini_STRC_aa{MINI_START}-{MINI_END}_CpG_depleted_t0.40 cpg={count_cpg(max_depleted)} cai={cai(max_depleted)}",
        max_depleted,
    )

    # Sanity: translation round-trip check
    codon_to_aa = {c: aa for aa, cs in CODON_TABLE.items() for c in cs}

    def translate(cds):
        aas = []
        for i in range(0, len(cds) - 2, 3):
            c = cds[i : i + 3]
            if c in ("TAA", "TAG", "TGA"):
                break
            aas.append(codon_to_aa.get(c, "X"))
        return "".join(aas)

    v0_prot = translate(v0)
    max_prot = translate(max_depleted)
    prot_ok_v0 = v0_prot == mini
    prot_ok_max = max_prot == mini
    print(f"\nTranslation round-trip: v0 {'OK' if prot_ok_v0 else 'FAIL'}, max-depleted {'OK' if prot_ok_max else 'FAIL'}")

    result = {
        "meta": {
            "date": "2026-04-21",
            "source_accession": STRC_ACC,
            "protein_length": len(strc),
            "mini_strc_window_aa": [MINI_START, MINI_END],
            "mini_strc_protein_length_aa": len(mini),
            "human_genome_cpg_per_kb_reference": HUMAN_GENOME_CPG_PER_KB,
            "codon_usage_source": "Kazusa Homo sapiens codon usage",
            "cai_method": "Sharp & Li 1987 geometric mean of relative adaptiveness",
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
