#!/usr/bin/env python3
"""
Prime editing pegRNA design for STRC c.4976A>C (p.E1659A), Misha's maternal
allele.

Context: [[Prime Editing for STRC]] confirmed the only viable precision editor
is PE3/PE3b because (a) the mutation site has no usable A on either strand for
ABE/ACBE and (b) CBE gives wrong outcome. PAM survey done: optimal SpCas9 NGG
PAM = TGG at chr15:43600540 (+ strand), 14 nt from the variant.

This script designs:
  1. pegRNA spacer (20 nt, protospacer upstream of PAM on + strand).
  2. RT template — installs the correction (C -> A on coding = G -> T on
     + strand). Standard RT lengths swept {10, 12, 14, 16} nt.
  3. PBS (primer binding site) annealing to the nicked strand, lengths {11, 12, 13, 14, 15}.
  4. PE3b nicking sgRNA on the opposite strand to bias repair toward the edit,
     within 40-90 nt of the edit.
  5. Composes the full pegRNA scaffold: spacer + sgRNA scaffold + RT + PBS.

Conventions (Anzalone 2019, Chen 2021 Nat Protoc):
  - Nick occurs between positions 17 and 18 of the 20-nt protospacer (3 nt
    5' of PAM).
  - RT template is reverse-complement of the desired edited sequence, written
    5' -> 3' on the pegRNA, reading FROM the nick INTO the edit.
  - PBS is reverse-complement of the genomic sequence IMMEDIATELY 5' of the
    nick site, 11-15 nt long.

GRCh38 coords (confirmed in [[Prime Editing for STRC]]):
  - STRC strand: minus (-)
  - Variant chr15:43600551  (+ strand: T=WT, G=mut ; coding = A/C)
  - PAM TGG: chr15:43600540-43600542 (+ strand)
  - Protospacer: chr15:43600520-43600539 (+ strand, 20 nt upstream of PAM)
  - Nick between chr15:43600536 and 43600537 (between protospacer pos 17 & 18)

Replication:
    /opt/miniconda3/bin/python3 pe_phase1_pegrna_design.py
"""
import json
import re
from io import StringIO
from pathlib import Path
from urllib.request import Request, urlopen

OUT_DIR = Path(__file__).resolve().parent
OUT_JSON = OUT_DIR / "pe_phase1_pegrna_design.json"

ENSEMBL = "https://rest.ensembl.org"

# GRCh38, + strand (reference strand). STRC gene is on - strand.
VARIANT_POS = 43600551  # chr15, + strand
PAM_START = 43600540    # chr15, + strand; PAM = TGG
PAM_END = 43600542
NICK_5NT = 43600536     # + strand nucleotide 5' of the nick
NICK_3NT = 43600537     # + strand nucleotide 3' of the nick
PROTOSPACER_START = 43600520
PROTOSPACER_END = 43600539

REGION_START = 43600400
REGION_END = 43600700

COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def revcomp(s):
    return "".join(COMP.get(b, "N") for b in s[::-1])


def fetch_region(chrom, start, end, strand=1):
    url = f"{ENSEMBL}/sequence/region/human/{chrom}:{start}..{end}:{strand}?coord_system_version=GRCh38"
    with urlopen(Request(url, headers={"Accept": "application/json"})) as r:
        return json.load(r)["seq"].upper()


def gc_pct(s):
    return round(100 * (s.count("G") + s.count("C")) / len(s), 1) if s else 0


def find_ngg_pams(seq, region_start_1idx):
    """Return list of (pam_start_1idx, pam_seq, strand) for all NGG/CCN in seq."""
    hits = []
    for i in range(len(seq) - 2):
        triplet = seq[i : i + 3]
        if triplet[1:] == "GG":
            hits.append({"pam_start_1idx": region_start_1idx + i, "pam_seq": triplet, "strand": "+"})
        if triplet[:2] == "CC":
            # On - strand, this is NGG. PAM start in + strand coords = i+2.
            hits.append({"pam_start_1idx": region_start_1idx + i + 2, "pam_seq": revcomp(triplet), "strand": "-"})
    return hits


def design_pegrna(plus_seq, region_start_1idx):
    """Build pegRNA for + strand PAM TGG @ 43600540."""
    # Protospacer on + strand, 20 nt immediately 5' of PAM
    proto_offset = PROTOSPACER_START - region_start_1idx
    proto = plus_seq[proto_offset : proto_offset + 20]
    pam_offset = PAM_START - region_start_1idx
    pam = plus_seq[pam_offset : pam_offset + 3]
    assert pam == "TGG", f"PAM mismatch: expected TGG got {pam}"
    # Nick falls between protospacer pos 17 and 18 (1-indexed within 20-nt)
    # On + strand that's between chr15:43600536 and 43600537
    nick_5_offset = NICK_5NT - region_start_1idx  # last nt before nick (+ strand)

    # For nick on the + strand (non-target strand of pegRNA), the pegRNA
    # spacer = protospacer (same sequence as + strand). Nicking is on the
    # PAM-containing strand = + strand.

    # Wait — clarification: with a + strand PAM, Cas9 targets the + strand and
    # nicks the PROTOSPACER-CONTAINING strand = + strand. The pegRNA spacer is
    # the + strand sequence (read 5'->3'). Nick is 3 nt upstream of PAM on +
    # strand, i.e. between + strand positions 43600536 and 43600537.

    # Edit site: chr15:43600551 (+ strand). We need + strand G -> T at pos
    # 43600551 (i.e., reference T is the WT; mutant has G; correction turns
    # the MUT G back into WT T on + strand).
    # RT template must encode the WT sequence from the nick point onward,
    # then some flanking WT sequence (homology arm), reverse-complemented on
    # the pegRNA.

    # Reading WT + strand from nick_3 (43600537) through 43600565 (28 nt RT window):
    rt_window_end = 43600565
    wt_rt_region_plus = plus_seq[NICK_3NT - region_start_1idx : rt_window_end - region_start_1idx + 1]
    # ^ this is what the + strand looks like in WT (after correction).
    # For Misha's MUT allele, position 43600551 reads G instead of T.
    # So the pegRNA must write T at that position.
    # RT templates (pegRNA 3' extension) are revcomp of this WT sequence.

    rt_designs = []
    for rt_len in (10, 12, 14, 16, 20, 24):
        wt_stretch = wt_rt_region_plus[:rt_len]
        # must contain the edit position (43600551) within rt_len
        if (NICK_3NT + rt_len - 1) < VARIANT_POS:
            continue
        rt_template = revcomp(wt_stretch)
        pos_in_rt = VARIANT_POS - NICK_3NT  # 0-indexed in wt_stretch
        rt_designs.append({
            "rt_length_nt": rt_len,
            "wt_plus_stretch": wt_stretch,
            "rt_template_5to3": rt_template,
            "edit_nt_position_in_rt_template": rt_len - 1 - pos_in_rt,
            "edit_in_rt_context": rt_template[rt_len - 1 - pos_in_rt : rt_len - pos_in_rt],
        })

    # PBS: RC of 11-15 nt immediately 5' of the nick on + strand
    pbs_designs = []
    for pbs_len in (11, 12, 13, 14, 15):
        pbs_region = plus_seq[NICK_5NT - region_start_1idx - pbs_len + 1 : NICK_5NT - region_start_1idx + 1]
        pbs = revcomp(pbs_region)
        pbs_designs.append({
            "pbs_length_nt": pbs_len,
            "plus_stretch_5of_nick": pbs_region,
            "pbs_5to3": pbs,
            "gc_pct": gc_pct(pbs),
        })

    # Default recommended combination: RT_len 13-16 with edit 3-6 nt from 3' of RT
    # (Anzalone), PBS_len 13 (Chen 2021 default)
    return {
        "protospacer_plus_strand_5to3": proto,
        "pam_plus_strand": pam,
        "nick_position_plus_strand_chr15": [NICK_5NT, NICK_3NT],
        "spacer_gc_pct": gc_pct(proto),
        "rt_templates": rt_designs,
        "pbs_templates": pbs_designs,
    }


def find_pe3b_nicking_sites(plus_seq, region_start_1idx):
    """Look for NGG PAMs on the - strand such that the nick falls 40-90 nt
    from the edit site on the - strand. The - strand nick biases repair.

    For PE3b (preferentially), the nicking sgRNA should target a sequence
    that is PRESENT only in the EDITED allele -> it matches the post-edit
    sequence but not the pre-edit. Implemented here as a flag: any nicking
    sgRNA whose 20-nt protospacer spans the variant position gets flagged
    as PE3b-eligible."""
    hits = []
    for i in range(len(plus_seq) - 22):
        # - strand PAM as CCN on + strand at pos i
        triplet = plus_seq[i : i + 3]
        if triplet[:2] != "CC":
            continue
        pam_start_plus = region_start_1idx + i + 2
        # Protospacer on - strand = revcomp of +strand[i+3 : i+23]
        proto_plus = plus_seq[i + 3 : i + 23]
        proto_minus = revcomp(proto_plus)
        # Nick on - strand: between positions 17 and 18 of the 20-nt
        # protospacer counting from its 5' end (which is on the - strand).
        # On + strand, the protospacer region spans i+3 .. i+22; the nick on
        # - strand corresponds to + strand position between i+5 and i+6.
        nick_minus_between_plus = (region_start_1idx + i + 5, region_start_1idx + i + 6)
        nick_site_plus = nick_minus_between_plus[1]
        dist_to_edit = abs(nick_site_plus - VARIANT_POS)
        if dist_to_edit < 30 or dist_to_edit > 100:
            continue
        # PE3b test: does protospacer cover the edited nucleotide?
        proto_span_plus = (region_start_1idx + i + 3, region_start_1idx + i + 22)
        pe3b_eligible = proto_span_plus[0] <= VARIANT_POS <= proto_span_plus[1]
        hits.append({
            "strand": "-",
            "pam_plus_coords": [pam_start_plus, pam_start_plus + 2],
            "protospacer_minus_5to3": proto_minus,
            "gc_pct": gc_pct(proto_minus),
            "nick_between_plus_coords": list(nick_minus_between_plus),
            "dist_nt_from_edit": dist_to_edit,
            "pe3b_eligible_spans_variant": pe3b_eligible,
        })

    # Also + strand PAMs (for PE3 classic). PE3b prefers - strand since our
    # pegRNA is + strand.
    return hits


def main():
    print(f"Fetching region chr15:{REGION_START}-{REGION_END} ...")
    plus_seq = fetch_region("15", REGION_START, REGION_END, 1)
    print(f"  region length: {len(plus_seq)} nt")
    # Sanity: variant at 43600551 is T on + strand in reference
    var_offset = VARIANT_POS - REGION_START
    ref_nt_plus = plus_seq[var_offset]
    print(f"  + strand nt at chr15:{VARIANT_POS}: {ref_nt_plus} (ref expects T)")
    assert ref_nt_plus == "T", "Reference genome mismatch at variant position"

    # Confirm PAM
    pam = plus_seq[PAM_START - REGION_START : PAM_END - REGION_START + 1]
    print(f"  PAM at chr15:{PAM_START}-{PAM_END} (+ strand): {pam} (expects TGG)")
    assert pam == "TGG"

    print("\nDesigning pegRNA (+ strand PAM TGG) ...")
    peg = design_pegrna(plus_seq, REGION_START)
    print(f"  spacer: {peg['protospacer_plus_strand_5to3']} (GC {peg['spacer_gc_pct']}%)")
    print(f"  RT templates computed: {len(peg['rt_templates'])}")
    for r in peg["rt_templates"]:
        print(f"    len={r['rt_length_nt']:2d} RT={r['rt_template_5to3']} edit_nt_in_RT='{r['edit_in_rt_context']}' at pos {r['edit_nt_position_in_rt_template']}")
    print(f"  PBS candidates:")
    for p in peg["pbs_templates"]:
        print(f"    len={p['pbs_length_nt']:2d} PBS={p['pbs_5to3']} (GC {p['gc_pct']}%)")

    print("\nScanning PE3b nicking sgRNA candidates (- strand PAMs, 30-100 nt from edit) ...")
    nickers = find_pe3b_nicking_sites(plus_seq, REGION_START)
    print(f"  total candidates: {len(nickers)}")
    pe3b = [n for n in nickers if n["pe3b_eligible_spans_variant"]]
    print(f"  PE3b-eligible (protospacer spans variant): {len(pe3b)}")
    for n in sorted(nickers, key=lambda x: (not x["pe3b_eligible_spans_variant"], x["dist_nt_from_edit"]))[:6]:
        flag = "PE3b" if n["pe3b_eligible_spans_variant"] else "PE3 "
        print(f"    {flag} dist={n['dist_nt_from_edit']:3d} proto- (5'->3'): {n['protospacer_minus_5to3']} GC {n['gc_pct']}%")

    # Recommend a single best combination
    rt_recommended = next((r for r in peg["rt_templates"] if r["rt_length_nt"] == 14), peg["rt_templates"][1])
    pbs_recommended = next((p for p in peg["pbs_templates"] if p["pbs_length_nt"] == 13), peg["pbs_templates"][2])
    best_nicker = None
    if pe3b:
        best_nicker = sorted(pe3b, key=lambda x: abs(x["dist_nt_from_edit"] - 60))[0]
    elif nickers:
        best_nicker = sorted(nickers, key=lambda x: abs(x["dist_nt_from_edit"] - 60))[0]

    pegrna_sequence = f"[spacer 5'->3']: {peg['protospacer_plus_strand_5to3']} + [sgRNA scaffold] + [3' ext: RT({rt_recommended['rt_length_nt']}) + PBS({pbs_recommended['pbs_length_nt']})] = {rt_recommended['rt_template_5to3']}{pbs_recommended['pbs_5to3']}"
    print(f"\nRecommended pegRNA: {pegrna_sequence}")
    if best_nicker:
        flag = "PE3b (spans variant — discriminates edited from unedited)" if best_nicker["pe3b_eligible_spans_variant"] else "PE3 (does not span variant)"
        print(f"Recommended nicking sgRNA ({flag}):")
        print(f"  protospacer 5'->3' ({best_nicker['strand']} strand): {best_nicker['protospacer_minus_5to3']} (GC {best_nicker['gc_pct']}%)")
        print(f"  distance from edit: {best_nicker['dist_nt_from_edit']} nt (optimal range 40-90)")

    result = {
        "meta": {
            "date": "2026-04-21",
            "variant": "STRC c.4976A>C p.E1659A (Misha maternal)",
            "genome_build": "GRCh38",
            "chrom": "15",
            "variant_pos_plus_strand": VARIANT_POS,
            "pam": "TGG (+ strand)",
            "pam_coords_plus": [PAM_START, PAM_END],
            "nick_coords_plus_between": [NICK_5NT, NICK_3NT],
            "region_fetched": [REGION_START, REGION_END],
        },
        "pegrna": peg,
        "nicking_sgrnas": nickers,
        "recommendation": {
            "pegrna_spacer_5to3": peg["protospacer_plus_strand_5to3"],
            "rt_template_5to3": rt_recommended["rt_template_5to3"],
            "rt_length_nt": rt_recommended["rt_length_nt"],
            "pbs_5to3": pbs_recommended["pbs_5to3"],
            "pbs_length_nt": pbs_recommended["pbs_length_nt"],
            "nicking_sgrna": best_nicker,
        },
    }

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
