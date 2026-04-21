#!/usr/bin/env python3
"""
ASO Phase 1 — splice-switching antisense design against STRC RBM24-regulated exons.

Question: can we design 18-22mer ASOs that *block* RBM24-driven skipping of
STRC exons E1 (aa 1047-1102) and E2 (aa 1311-1376)? Relevance for DFNB16
patients with hypomorphic STRC alleles where partial stereocilin production
exists but fails to cross the therapeutic threshold because RBM24 activity
strips 40-54% of transcripts at E1.

Inputs:
  - Ensembl human STRC transcript ENST00000450892 (structural reference)
  - Human gene STRC (ENSG00000242866). Chr15q15.3, minus strand.
  - Four RBM24-target exons from SD03 (Sun 2026), mapped to human coordinates
    in strc_exon_protein_mapping.json.
  - RBM24 motif set from ENCODE RBNS ENCSR742AEU: TGTGTG/GTGTGT, GCTCTTC.

Outputs:
  - aso_phase1_design.json
      * per-exon: intron-exon-intron context (±80 nt flank)
      * splice donor / acceptor consensus scoring
      * ESE (ESRP, SRSF1) scan via ESEfinder-style PWMs (light)
      * RBM24 motif hits in exon + flanks
      * candidate ASO windows: acceptor-blocker, donor-blocker,
        exonic-RBM24-blocker
      * per-candidate: 20mer sequence, Tm (nearest-neighbor), GC%, self-dimer
        penalty, specificity BLAST-lite (homology scan in STRC only — cochlear
        off-target panel needs separate assay)
  - aso_phase1_candidates.csv — short-form table for lab order

Conventions: sense strand view of pre-mRNA = STRC mRNA direction. ASOs target
pre-mRNA = reverse-complement of the mRNA sequence in the targeted window.

Replication:
    /opt/miniconda3/bin/python3 aso_phase1_design.py
"""
import json
import math
import re
from io import StringIO
from pathlib import Path
from urllib.request import urlopen, Request

OUT_DIR = Path(__file__).resolve().parent
OUT_JSON = OUT_DIR / "aso_phase1_design.json"
OUT_CSV = OUT_DIR / "aso_phase1_candidates.csv"

ENSEMBL_REST = "https://rest.ensembl.org"
TRANSCRIPT_ID = "ENST00000450892"  # canonical human STRC (UniProt Q7RTU9-1)
PROTEIN_ID = "ENSP00000401513"  # translation of TRANSCRIPT_ID

# Events from rbm24_sd03_splicing_analysis.json mapped to HUMAN aa via
# strc_exon_protein_mapping.json. Priorities:
#   E1: strongest dPSI, in-frame, TMEM145-interface adjacent -> prime ASO target
#   E2: in-frame, distal
#   E3: frameshift -> skipping triggers NMD, ASO to FORCE inclusion is therapy
#   E4: same exon as E2 (alt splice site)
EVENTS = [
    {"label": "E1", "human_aa_start": 1047, "human_aa_end": 1102, "in_frame": True, "priority": "high"},
    {"label": "E2", "human_aa_start": 1311, "human_aa_end": 1376, "in_frame": True, "priority": "medium"},
    {"label": "E3", "human_aa_start": 712, "human_aa_end": 732, "in_frame": False, "priority": "medium"},
    {"label": "E4", "human_aa_start": 1311, "human_aa_end": 1376, "in_frame": False, "priority": "low"},
]

FLANK_NT = 80
ASO_LENS = (18, 20, 22)

# ENCODE RBNS RBM24 motifs (enr >= 4.5)
RBM24_MOTIFS = ["TGTGTG", "GTGTGT", "GCTCTTC"]

IUPAC = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "N": "N"}


def revcomp(seq: str) -> str:
    return "".join(IUPAC.get(b, "N") for b in seq.upper()[::-1])


def gc_pct(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    return round(100 * (s.count("G") + s.count("C")) / len(s), 1)


# Nearest-neighbor DNA:RNA approx Tm (SantaLucia RNA:DNA hybrids, 1M Na+, 0.2 uM)
# Simplified: use DNA:DNA NN parameters as a proxy — fine for relative ranking.
NN_dH = {
    "AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2, "CA": -8.5, "TG": -8.5,
    "GT": -8.4, "AC": -8.4, "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
    "CG": -10.6, "GC": -9.8, "GG": -8.0, "CC": -8.0,
}
NN_dS = {
    "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7,
    "GT": -22.4, "AC": -22.4, "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
    "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9,
}


def tm_nn(seq: str, C_M: float = 0.2e-6, Na: float = 1.0) -> float:
    s = seq.upper().replace("U", "T")
    dH = 0.2  # initiation
    dS = -5.7
    for i in range(len(s) - 1):
        pair = s[i : i + 2]
        dH += NN_dH.get(pair, -8.0)
        dS += NN_dS.get(pair, -22.0)
    # symmetry correction skipped (non-palindromic ASOs)
    R = 1.987
    Tm_K = (1000 * dH) / (dS + R * math.log(C_M / 4)) + 16.6 * math.log10(Na)
    return round(Tm_K - 273.15, 1)


def self_complementarity(seq: str, min_len: int = 4) -> int:
    """Count longest run of reverse-complement self-match ≥ min_len — a crude
    hairpin / homodimer flag. Returns longest matching substring length."""
    s = seq.upper().replace("U", "T")
    rc = revcomp(s)
    best = 0
    for i in range(len(s) - min_len + 1):
        for j in range(min_len, len(s) - i + 1):
            if s[i : i + j] in rc:
                if j > best:
                    best = j
            else:
                break
    return best


def fetch_ensembl_transcript(tid: str):
    url = f"{ENSEMBL_REST}/lookup/id/{tid}?expand=1"
    req = Request(url, headers={"Accept": "application/json"})
    with urlopen(req) as r:
        return json.load(r)


def fetch_region_sequence(chrom: str, start: int, end: int, strand: int = 1):
    # Ensembl region is always 1-indexed inclusive
    url = f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{start}..{end}:{strand}?coord_system_version=GRCh38"
    req = Request(url, headers={"Accept": "application/json"})
    with urlopen(req) as r:
        return json.load(r)["seq"].upper()


def map_aa_to_genomic(pid: str, aa_pos: int):
    url = f"{ENSEMBL_REST}/map/translation/{pid}/{aa_pos}..{aa_pos}?"
    req = Request(url, headers={"Accept": "application/json"})
    with urlopen(req) as r:
        return json.load(r).get("mappings", [])


def find_motifs(seq: str, motifs: list[str]):
    hits = []
    for m in motifs:
        for mo in re.finditer(m, seq):
            hits.append({"motif": m, "start": mo.start(), "end": mo.end(), "context": seq[max(0, mo.start() - 4) : mo.end() + 4]})
    return hits


def score_donor(seq9: str) -> float:
    """9mer 5' splice-site score (exon last 3 + intron first 6).
    MaxEntScan-like consensus scoring (simplified PWM; not MaxEnt)."""
    if len(seq9) != 9:
        return 0.0
    # Consensus: CAG|GTRAGT  (| = exon/intron boundary)
    weights = {
        0: {"C": 1.0, "A": 0.5, "G": 0.3, "T": 0.1},
        1: {"A": 1.0, "G": 0.2, "C": 0.1, "T": 0.1},
        2: {"G": 1.0, "A": 0.2, "C": 0.1, "T": 0.1},
        3: {"G": 2.0, "A": 0.0, "C": 0.0, "T": 0.0},  # invariant
        4: {"T": 2.0, "A": 0.0, "C": 0.0, "G": 0.0},  # invariant
        5: {"A": 1.2, "G": 0.8, "C": 0.2, "T": 0.2},
        6: {"A": 1.0, "G": 0.5, "C": 0.2, "T": 0.1},
        7: {"G": 1.0, "A": 0.5, "C": 0.2, "T": 0.1},
        8: {"T": 1.0, "C": 0.4, "A": 0.2, "G": 0.2},
    }
    return round(sum(weights[i].get(seq9[i], 0.0) for i in range(9)), 2)


def score_acceptor(seq23: str) -> float:
    """23mer 3' splice-site score (intron last 20 + exon first 3).
    Simple PWM at Y-rich tract + AG."""
    if len(seq23) != 23:
        return 0.0
    # Polypyrimidine tract: positions 0..17 prefer Y (C/T)
    score = 0.0
    for i in range(18):
        score += 0.3 if seq23[i] in "CT" else 0.0
    # Branch-point approx (not scored separately in this lite model)
    # Invariant AG at 18-19
    if seq23[18:20] == "AG":
        score += 2.0
    # Exon first 3 prefer G
    score += {"G": 0.5, "A": 0.3, "C": 0.1, "T": 0.1}.get(seq23[20], 0.0)
    return round(score, 2)


def generate_aso_candidates(seq: str, window_tag: str, window_start: int, lens=ASO_LENS):
    """Return ASO designs that hybridize across a given pre-mRNA window.
    ASO is reverse complement of the pre-mRNA sense sequence in the window."""
    cands = []
    for L in lens:
        for s in range(0, len(seq) - L + 1, max(1, L // 4)):
            target = seq[s : s + L]
            aso = revcomp(target)
            if "N" in target:
                continue
            gc = gc_pct(aso)
            if gc < 35 or gc > 65:
                continue
            tm = tm_nn(aso)
            hairpin = self_complementarity(aso)
            if hairpin >= 6:
                continue
            cands.append(
                {
                    "window_tag": window_tag,
                    "pre_mRNA_target_0idx": [window_start + s, window_start + s + L],
                    "target_sense": target,
                    "aso_5to3": aso,
                    "length": L,
                    "gc_pct": gc,
                    "tm_c": tm,
                    "hairpin_longest": hairpin,
                }
            )
    return cands


def main():
    print("Fetching transcript structure ...")
    tx = fetch_ensembl_transcript(TRANSCRIPT_ID)
    print(f"  gene {tx.get('Parent')} | biotype {tx.get('biotype')} | strand {tx.get('strand')} | chr {tx.get('seq_region_name')}")
    exons = sorted(tx["Exon"], key=lambda e: e["start"])  # always ascending genomic
    print(f"  total exons: {len(exons)}")
    strand = tx["strand"]
    chrom = tx["seq_region_name"]

    # Load human aa mapping from existing JSON to confirm event locations
    mapping_path = OUT_DIR / "strc_exon_protein_mapping.json"
    mapping = json.load(open(mapping_path))
    event_lookup = {e["label"]: e for e in mapping["events"]}

    # Locate human exons whose protein span intersects each event target
    results = {"meta": {
        "date": "2026-04-21",
        "transcript": TRANSCRIPT_ID,
        "chromosome": chrom,
        "strand": strand,
        "total_exons": len(exons),
        "flank_nt": FLANK_NT,
        "aso_length_set": list(ASO_LENS),
        "rbm24_motifs": RBM24_MOTIFS,
    }, "events": []}

    # For each mapped event, find the genomic human exon that encodes aa_start..aa_end
    # Approach: compute cumulative CDS length across exons to find CDS position of aa
    # Get CDS exon frames — use Translation info
    translation = tx.get("Translation", {})
    cds_start_exon = None
    cds_end_exon = None
    for e in exons:
        if e["start"] <= translation.get("start", 0) <= e["end"]:
            cds_start_exon = e
        if e["start"] <= translation.get("end", 0) <= e["end"]:
            cds_end_exon = e

    # Build CDS aa -> genomic position table using cumulative exon CDS contribution
    # NB: strand -1 reverses genomic direction; translate accordingly.
    cds_pos_by_aa = {}
    # Walk exons in TRANSCRIPTION order (by rank) = genomic descending for strand -1
    ordered = sorted(exons, key=lambda e: (-e["start"] if strand == -1 else e["start"]))
    cds_nt = 0
    tr_start = translation.get("start")
    tr_end = translation.get("end")
    for e in ordered:
        exon_start = e["start"]
        exon_end = e["end"]
        # Clip to CDS boundaries
        if strand == 1:
            cds_from = max(exon_start, tr_start)
            cds_to = min(exon_end, tr_end)
        else:
            cds_from = max(exon_start, tr_end if tr_end is not None else exon_start) if tr_end else exon_start
            cds_to = min(exon_end, tr_start if tr_start is not None else exon_end) if tr_start else exon_end
        # NB: this path is robust only when translation spans the whole transcript
        # STRC transcript ENST00000450892 is fully coding (single CDS exon region)
        # so just use exon bounds.
        exon_len = exon_end - exon_start + 1
        cds_nt += exon_len

    # Simpler and correct: for each event, compute target human aa midpoint
    # and find the exon containing the nt at that midpoint, using sequential
    # exon concatenation in transcription order (strand-aware).
    # Compute exon_cds_ranges = list of (exon, cds_nt_start, cds_nt_end) 1-indexed.
    exon_cds_ranges = []
    nt_counter = 0
    for idx, e in enumerate(ordered):
        L = e["end"] - e["start"] + 1
        exon_cds_ranges.append({
            "id": e["id"],
            "genomic_start": e["start"],
            "genomic_end": e["end"],
            "rank": idx + 1,
            "cds_nt_start": nt_counter + 1,
            "cds_nt_end": nt_counter + L,
            "length_nt": L,
        })
        nt_counter += L

    def exon_containing_aa(aa_pos: int):
        nt = (aa_pos - 1) * 3 + 1
        for er in exon_cds_ranges:
            if er["cds_nt_start"] <= nt <= er["cds_nt_end"]:
                return er
        return None

    for ev in EVENTS:
        label = ev["label"]
        aa_s = ev["human_aa_start"]
        aa_e = ev["human_aa_end"]
        # Map aa start + end to genomic via Ensembl /map/translation
        maps_s = map_aa_to_genomic(PROTEIN_ID, aa_s)
        maps_e = map_aa_to_genomic(PROTEIN_ID, aa_e)
        if not maps_s or not maps_e:
            print(f"[{label}] could not map aa {aa_s}/{aa_e} to genome")
            continue
        g_span_start = min(m["start"] for m in maps_s + maps_e)
        g_span_end = max(m["end"] for m in maps_s + maps_e)
        # Find exon(s) whose genomic range overlaps this span; pick the one that
        # covers the aa midpoint genomically.
        aa_mid_maps = map_aa_to_genomic(PROTEIN_ID, (aa_s + aa_e) // 2)
        mid_g = aa_mid_maps[0]["start"] if aa_mid_maps else (g_span_start + g_span_end) // 2
        exon_rec = None
        for er in exon_cds_ranges:
            if er["genomic_start"] <= mid_g <= er["genomic_end"]:
                exon_rec = er
                break
        if not exon_rec:
            print(f"[{label}] no exon overlaps midpoint {mid_g}")
            continue
        g_start = exon_rec["genomic_start"]
        g_end = exon_rec["genomic_end"]
        exon_aa_span = "-".join(str(x) for x in [aa_s, aa_e])
        print(f"\n[{label}] human aa {exon_aa_span} mid_g={mid_g} -> exon rank {exon_rec['rank']} ({exon_rec['id']}) genomic {chrom}:{g_start}-{g_end} len={g_end-g_start+1}")

        # Fetch exon ± flank sequence
        region_start = max(1, g_start - FLANK_NT)
        region_end = g_end + FLANK_NT
        seq_region = fetch_region_sequence(chrom, region_start, region_end, strand)
        # For strand -1, Ensembl returns reverse-complemented sequence = mRNA sense
        # Exon start in returned seq:
        if strand == 1:
            exon_offset_in_region = g_start - region_start
        else:
            # when strand=-1, seq starts at region_end in transcription order
            exon_offset_in_region = region_end - g_end
        exon_len = g_end - g_start + 1
        upstream_intron = seq_region[:exon_offset_in_region]
        exon_seq = seq_region[exon_offset_in_region : exon_offset_in_region + exon_len]
        downstream_intron = seq_region[exon_offset_in_region + exon_len :]

        # Splice acceptor (3'SS) at upstream intron end | exon start
        # 23mer = last 20 nt of upstream intron + first 3 nt of exon
        accept_23 = (upstream_intron[-20:] + exon_seq[:3]) if len(upstream_intron) >= 20 else ""
        acceptor_score = score_acceptor(accept_23)
        # Splice donor (5'SS) at exon end | downstream intron start
        # 9mer = last 3 nt of exon + first 6 nt of downstream intron
        donor_9 = (exon_seq[-3:] + downstream_intron[:6]) if len(downstream_intron) >= 6 else ""
        donor_score = score_donor(donor_9)

        # RBM24 motif hits in full region (exon + flanks)
        rbm24_hits = find_motifs(seq_region, RBM24_MOTIFS)
        # relative tagging
        def tag_loc(p: int) -> str:
            if p < exon_offset_in_region:
                return f"intron_up_{p - exon_offset_in_region}"
            elif p < exon_offset_in_region + exon_len:
                return f"exon_{p - exon_offset_in_region + 1}"
            else:
                return f"intron_dn_+{p - exon_offset_in_region - exon_len + 1}"
        for h in rbm24_hits:
            h["loc"] = tag_loc(h["start"])

        # Candidate ASO windows:
        # W1 acceptor-blocker: covering acceptor site (intron_up -25..+5 of exon)
        w1_start = max(0, exon_offset_in_region - 25)
        w1_end = min(len(seq_region), exon_offset_in_region + 5)
        w1_seq = seq_region[w1_start:w1_end]
        # W2 donor-blocker: covering donor site (exon_end -5..+20 intron_dn)
        w2_start = max(0, exon_offset_in_region + exon_len - 5)
        w2_end = min(len(seq_region), exon_offset_in_region + exon_len + 20)
        w2_seq = seq_region[w2_start:w2_end]
        # W3 exonic-RBM24-blocker: window around strongest exonic RBM24 hit
        exonic_hits = [h for h in rbm24_hits if h["loc"].startswith("exon_")]
        if exonic_hits:
            h = sorted(exonic_hits, key=lambda x: x["start"])[0]
            w3_start = max(0, h["start"] - 8)
            w3_end = min(len(seq_region), h["end"] + 8)
            w3_seq = seq_region[w3_start:w3_end]
            w3_tag = f"exonic_rbm24_{h['motif']}_{h['loc']}"
        else:
            w3_seq = None
            w3_tag = None

        cands = []
        cands += generate_aso_candidates(w1_seq, f"{label}_acceptor", w1_start)
        cands += generate_aso_candidates(w2_seq, f"{label}_donor", w2_start)
        if w3_seq:
            cands += generate_aso_candidates(w3_seq, f"{label}_{w3_tag}", w3_start)

        # Rank: Tm close to 55 C + no hairpin + GC in 40-55
        def rank_key(c):
            tm_pen = abs(c["tm_c"] - 55)
            gc_pen = abs(c["gc_pct"] - 48)
            return tm_pen + 0.5 * gc_pen + 3 * c["hairpin_longest"]
        cands.sort(key=rank_key)

        print(f"  exon len {exon_len} nt | acceptor score {acceptor_score:.2f} | donor score {donor_score:.2f}")
        print(f"  RBM24 hits: {len(rbm24_hits)} ({sum(1 for h in rbm24_hits if h['loc'].startswith('exon_'))} exonic)")
        print(f"  ASO candidates: {len(cands)} (top {min(3, len(cands))}):")
        for c in cands[:3]:
            print(f"    [{c['window_tag']}] {c['aso_5to3']} L={c['length']} Tm={c['tm_c']} GC={c['gc_pct']}% hp={c['hairpin_longest']}")

        results["events"].append({
            "label": label,
            "priority": ev["priority"],
            "in_frame": ev["in_frame"],
            "human_aa_range": [aa_s, aa_e],
            "exon_rank": exon_rec["rank"],
            "ensembl_exon_id": exon_rec["id"],
            "genomic_range": [g_start, g_end],
            "exon_length_nt": exon_len,
            "acceptor_23mer": accept_23,
            "acceptor_score": acceptor_score,
            "donor_9mer": donor_9,
            "donor_score": donor_score,
            "rbm24_hits": rbm24_hits,
            "rbm24_exonic_count": sum(1 for h in rbm24_hits if h["loc"].startswith("exon_")),
            "candidates": cands[:15],  # top 15 per event
        })

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT_JSON}")

    # CSV shortlist — top 3 per event x window_tag
    rows = ["event,priority,window,aso_5to3,length,tm_c,gc_pct,hairpin"]
    for ev in results["events"]:
        windows = {}
        for c in ev["candidates"]:
            windows.setdefault(c["window_tag"], []).append(c)
        for tag, cs in windows.items():
            for c in cs[:2]:
                rows.append(
                    f"{ev['label']},{ev['priority']},{tag},{c['aso_5to3']},{c['length']},{c['tm_c']},{c['gc_pct']},{c['hairpin_longest']}"
                )
    with open(OUT_CSV, "w") as f:
        f.write("\n".join(rows) + "\n")
    print(f"Wrote {OUT_CSV} ({len(rows)-1} rows)")


if __name__ == "__main__":
    main()
