#!/usr/bin/env python3
"""
Prime Editing Phase 4 — Cas-OFFinder off-target scan for E1659A pegRNAs.

Reads the discriminating candidates from
[[STRC pe_phase3_allele_discrimination]] (~/Brain/research/strc/models/
pe_phase3_allele_discrimination.json), runs Cas-OFFinder against hg38
chr15 (where STRC and the STRCP1 paralog both live), and counts
off-target sites per guide at increasing mismatch tolerance.

Method:
  1. For each PAM variant (SpCas9-NGG, SpG-NGN, enCas9-NGN, SpRY-NRN,
     SpCas9NG-NG, SaCas9-NNGRRT) take the top discriminating
     `all_ranked` entries.
  2. Build a single Cas-OFFinder input file per PAM (CPU mode, max
     4 mismatches — covers seed-region tolerance plus a buffer).
  3. Run cas-offinder, parse the TSV output.
  4. For each guide, count off-targets binned by mismatch count
     (0, 1, 2, 3, 4).
  5. Gate per guide:
       - ≤5 off-targets at <=3 mismatches → PASS
       - 0 off-targets at 0 mismatches off the on-target locus → REQUIRED

The on-target locus (chr15:43600548-43600567 for the SpCas9-NGG strong
hit) is excluded when counting off-targets.

Run env: any conda env with cas-offinder installed (see 0-setup.sh).
Genome: ~/Brain/research/strc/genomes/hg38_chr15.fa (downloaded by setup).

Note: chr15-only is sufficient for the STRC + STRCP1 paralog risk
(both on chr15q15.3). For full off-target risk across the genome,
swap GENOME_PATH to a full hg38 multi-FASTA and re-run; ~30x slower.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

WORK_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
PHASE3_JSON = WORK_DIR / "pe_phase3_allele_discrimination.json"
OUT_JSON = WORK_DIR / "pe_phase4_cas_offinder.json"
GENOME_PATH = Path("/Users/egorlyfar/Brain/research/strc/genomes/hg38_chr15.fa")

MAX_MISMATCHES = 4
PASS_GATE_MAX_OFFTARGETS_AT_3MM = 5

# PAM patterns in Cas-OFFinder syntax:
PAM_PATTERNS = {
    "SpCas9_NGG":  ("NNNNNNNNNNNNNNNNNNNN", "NGG"),
    "SpG_NGN":     ("NNNNNNNNNNNNNNNNNNNN", "NGN"),
    "enCas9_NGN":  ("NNNNNNNNNNNNNNNNNNNN", "NGN"),
    "SpRY_NRN":    ("NNNNNNNNNNNNNNNNNNNN", "NRN"),
    "SpCas9NG_NG": ("NNNNNNNNNNNNNNNNNNNN", "NG"),
    "SaCas9_NNGRRT": ("NNNNNNNNNNNNNNNNNNNNN", "NNGRRT"),
}

DEVICE = "C"  # CPU. Use "G0" or "A0" for GPU/Apple Silicon if Cas-OFFinder build supports it.


def log(msg):
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}", flush=True)


def env_check():
    if shutil.which("cas-offinder") is None:
        log("FATAL — cas-offinder not on PATH. Run 0-setup.sh first.")
        sys.exit(2)
    if not GENOME_PATH.exists():
        log(f"FATAL — genome FASTA missing: {GENOME_PATH}. Run 0-setup.sh first.")
        sys.exit(2)
    if not PHASE3_JSON.exists():
        log(f"FATAL — Phase 3 input missing: {PHASE3_JSON}")
        sys.exit(3)
    log(f"env OK: cas-offinder, genome ({GENOME_PATH.stat().st_size/1024/1024:.1f} MiB)")


def load_guides():
    """Pull discriminating top candidates per PAM. Dedupe by protospacer."""
    data = json.loads(PHASE3_JSON.read_text())
    audit = data["audit_by_variant"]

    guides_by_pam = defaultdict(list)
    seen = set()
    for pam_name, pam_block in audit.items():
        if pam_name not in PAM_PATTERNS:
            log(f"  skipping unknown PAM: {pam_name}")
            continue
        ranked = pam_block.get("all_ranked", [])
        for cand in ranked:
            if not cand.get("allele_discriminating"):
                continue
            # Strand-dependent field — Phase 3 stores the protospacer in 5'→3'
            # of whichever strand the guide reads on. Cas-OFFinder also wants
            # the guide in 5'→3' on its own strand.
            if cand.get("strand") == "+":
                proto = cand.get("protospacer_plus_5to3")
            else:
                proto = cand.get("protospacer_minus_5to3")
            if not proto:
                log(f"  skip {pam_name} cand: no protospacer field for strand={cand.get('strand')}")
                continue
            key = (pam_name, proto)
            if key in seen:
                continue
            seen.add(key)
            guides_by_pam[pam_name].append({
                "protospacer": proto,
                "strand": cand["strand"],
                "pam_coords_plus": cand["pam_coords_plus"],
                "protospacer_coords_plus": cand["protospacer_coords_plus"],
                "discrimination_grade": cand["discrimination_grade"],
                "nick_to_edit_nt": cand["nick_to_edit_nt"],
                "variant_pos_in_protospacer_1idx": cand["variant_pos_in_protospacer_1idx"],
            })

    total = sum(len(g) for g in guides_by_pam.values())
    log(f"loaded {total} discriminating guides across {len(guides_by_pam)} PAM variants:")
    for pam, gs in guides_by_pam.items():
        log(f"  {pam}: {len(gs)}")
    return guides_by_pam


_GENOME_CACHE: dict | None = None
IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"},
    "N": {"A", "C", "G", "T"},
}


def _load_genome_uint8():
    """Load chr15 FASTA into a single uppercase uint8 array (A=0,C=1,G=2,T=3,N=4)."""
    global _GENOME_CACHE
    if _GENOME_CACHE is not None:
        return _GENOME_CACHE
    import numpy as np
    log(f"  loading genome from {GENOME_PATH}")
    chrom_name = None
    parts = []
    with open(GENOME_PATH) as f:
        for line in f:
            if line.startswith(">"):
                chrom_name = line[1:].strip().split()[0]
            else:
                parts.append(line.strip().upper())
    seq = "".join(parts)
    log(f"  genome loaded: chrom={chrom_name}, length={len(seq)}")
    table = np.full(256, 4, dtype=np.uint8)
    table[ord("A")] = 0
    table[ord("C")] = 1
    table[ord("G")] = 2
    table[ord("T")] = 3
    arr_b = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    arr = table[arr_b]
    _GENOME_CACHE = {"chr": chrom_name, "fwd": arr, "length": len(seq)}
    return _GENOME_CACHE


def _revcomp_uint8(arr):
    """Reverse-complement a uint8 (A=0,C=1,G=2,T=3,N=4) array."""
    import numpy as np
    comp_table = np.array([3, 2, 1, 0, 4], dtype=np.uint8)  # A↔T, C↔G, N→N
    return comp_table[arr][::-1]


def _python_offtarget_search(pam_name: str, used_guides: list[dict]) -> dict[str, list[dict]]:
    """Pure-numpy fallback when cas-offinder OpenCL fails on macOS Apple Silicon.

    Returns hits_by_guide mapping {query_seq -> [hit dicts]} matching the
    cas-offinder output format. Only chr15 (the loaded FASTA) is searched.
    """
    import numpy as np
    template, pam = PAM_PATTERNS[pam_name]
    proto_len = len(template)
    pam_len = len(pam)
    full_len = proto_len + pam_len

    g = _load_genome_uint8()
    chrom = g["chr"]
    fwd = g["fwd"]
    rev = _revcomp_uint8(fwd)
    L = len(fwd)

    # PAM mask: precompute which bases are allowed at each PAM position
    pam_allowed = []
    for ch in pam:
        allowed = IUPAC[ch]
        mask = np.zeros(5, dtype=bool)
        for b in allowed:
            mask["ACGT".index(b)] = True
        pam_allowed.append(mask)

    def find_pam_positions(strand_arr):
        """Return positions i such that strand_arr[i:i+full_len] has matching PAM at the end."""
        nbases = L - full_len + 1
        if nbases <= 0:
            return np.empty(0, dtype=np.int64)
        # Check PAM positions: indices [proto_len, proto_len+1, ..., proto_len+pam_len-1]
        keep = np.ones(nbases, dtype=bool)
        for k in range(pam_len):
            col = strand_arr[proto_len + k : proto_len + k + nbases]
            allowed = pam_allowed[k]
            keep &= allowed[col]
        return np.where(keep)[0]

    hits_by_guide: dict[str, list[dict]] = {}
    for g_meta in used_guides:
        proto = g_meta["protospacer"]
        query_seq = proto + ("N" * pam_len)
        proto_arr = np.array([{"A": 0, "C": 1, "G": 2, "T": 3}.get(c, 4) for c in proto], dtype=np.uint8)

        guide_hits = []
        for strand_label, strand_arr in (("+", fwd), ("-", rev)):
            pam_positions = find_pam_positions(strand_arr)
            if len(pam_positions) == 0:
                continue
            # Build (M, proto_len) matrix of protospacer windows at PAM positions
            # Use stride trick: contiguous slices already
            for chunk_start in range(0, len(pam_positions), 200_000):
                chunk = pam_positions[chunk_start:chunk_start + 200_000]
                M = len(chunk)
                # Extract windows
                windows = np.empty((M, proto_len), dtype=np.uint8)
                for k in range(proto_len):
                    windows[:, k] = strand_arr[chunk + k]
                mismatches = (windows != proto_arr).sum(axis=1)
                # Also count Ns in the genome window as mismatches (cas-offinder excludes those)
                mismatches += (windows == 4).sum(axis=1)
                ok = mismatches <= MAX_MISMATCHES
                if not ok.any():
                    continue
                ok_idx = np.where(ok)[0]
                for j in ok_idx:
                    pos_in_strand = int(chunk[j])
                    if strand_label == "+":
                        plus_pos = pos_in_strand
                    else:
                        plus_pos = L - pos_in_strand - full_len
                    match_seq_arr = windows[j].tolist() + [strand_arr[chunk[j] + proto_len + k] for k in range(pam_len)]
                    match_seq = "".join("ACGTN"[b] for b in match_seq_arr)
                    guide_hits.append({
                        "chr": chrom,
                        "pos": plus_pos,
                        "match_seq": match_seq,
                        "strand": strand_label,
                        "mismatches": int(mismatches[j]),
                    })
        hits_by_guide[query_seq] = guide_hits
    return hits_by_guide


def run_cas_offinder(pam_name: str, guides: list[dict], outdir: Path) -> list[dict]:
    template, pam = PAM_PATTERNS[pam_name]
    expected_proto_len = len(template)

    input_file = outdir / f"{pam_name}.in"
    output_file = outdir / f"{pam_name}.out"

    lines = [str(GENOME_PATH), template + pam]
    used_guides = []
    for g in guides:
        proto = g["protospacer"]
        if len(proto) != expected_proto_len:
            log(f"  WARN — {pam_name}: protospacer length {len(proto)} != template {expected_proto_len}, skipping {proto}")
            continue
        # Pad protospacer with PAM-length N to match template+pam, then mismatch count
        query_seq = proto + ("N" * len(pam))
        lines.append(f"{query_seq} {MAX_MISMATCHES}")
        used_guides.append(g)

    if not used_guides:
        log(f"  {pam_name}: no usable guides")
        return []

    input_file.write_text("\n".join(lines) + "\n")
    log(f"  {pam_name}: running cas-offinder ({len(used_guides)} guides, mm<={MAX_MISMATCHES})")

    cmd = ["cas-offinder", str(input_file), DEVICE, str(output_file)]
    t0 = datetime.now()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = (datetime.now() - t0).total_seconds()
    used_python_fallback = False
    if proc.returncode != 0 or "Failed:" in proc.stderr:
        log(f"  cas-offinder failed (likely OpenCL on macOS Apple Silicon):")
        log(f"    stderr: {proc.stderr[:200]}")
        log(f"  -> falling back to pure-numpy off-target search")
        used_python_fallback = True

    if used_python_fallback:
        t0 = datetime.now()
        hits_by_guide = _python_offtarget_search(pam_name, used_guides)
        dt = (datetime.now() - t0).total_seconds()
        log(f"  {pam_name}: python fallback done in {dt:.1f}s")
    else:
        log(f"  {pam_name}: cas-offinder done in {dt:.1f}s, parsing {output_file.stat().st_size} bytes")
        # Parse output: TSV "query  chr  pos  match_seq  strand  mismatches"
        hits_by_guide = defaultdict(list)
        for line in output_file.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            query, chrom, pos, match_seq, strand, mm = parts[:6]
            hits_by_guide[query].append({
                "chr": chrom,
                "pos": int(pos),
                "match_seq": match_seq,
                "strand": strand,
                "mismatches": int(mm),
            })

    # Bucket per guide and exclude on-target locus
    results = []
    for g in used_guides:
        proto = g["protospacer"]
        query_seq = proto + ("N" * len(pam))
        all_hits = hits_by_guide.get(query_seq, [])

        on_target_lo = min(g["protospacer_coords_plus"]) - 5
        on_target_hi = max(g["protospacer_coords_plus"]) + 5

        offtargets = [
            h for h in all_hits
            if not (h["chr"] in {"chr15", "15"}
                    and on_target_lo <= h["pos"] <= on_target_hi)
        ]
        on_target_hits = [h for h in all_hits if h not in offtargets]

        bins = defaultdict(int)
        for h in offtargets:
            bins[h["mismatches"]] += 1

        n_le3 = sum(c for mm, c in bins.items() if mm <= 3)
        n_perfect_off = bins.get(0, 0)
        passes = (n_le3 <= PASS_GATE_MAX_OFFTARGETS_AT_3MM) and (n_perfect_off == 0)

        results.append({
            "pam": pam_name,
            "protospacer": proto,
            "strand": g["strand"],
            "discrimination_grade": g["discrimination_grade"],
            "nick_to_edit_nt": g["nick_to_edit_nt"],
            "on_target_hits": len(on_target_hits),
            "offtarget_count_total": len(offtargets),
            "offtargets_by_mismatch": dict(bins),
            "n_offtargets_le_3mm": n_le3,
            "n_perfect_offtargets": n_perfect_off,
            "passes_gate": passes,
            "offtarget_loci_top10": [
                {"chr": h["chr"], "pos": h["pos"], "mm": h["mismatches"], "strand": h["strand"]}
                for h in sorted(offtargets, key=lambda x: x["mismatches"])[:10]
            ],
        })
        log(f"    {proto} mm<=3 off-targets: {n_le3}  perfect off: {n_perfect_off}  -> {'PASS' if passes else 'FAIL'}")

    return results


def main():
    env_check()
    guides_by_pam = load_guides()

    workdir = WORK_DIR / "pe_phase4_workdir"
    workdir.mkdir(exist_ok=True, parents=True)

    all_results = []
    for pam_name, gs in guides_by_pam.items():
        log("=" * 60)
        log(f"PAM: {pam_name}")
        log("=" * 60)
        all_results.extend(run_cas_offinder(pam_name, gs, workdir))

    n_pass = sum(1 for r in all_results if r["passes_gate"])
    n_total = len(all_results)
    verdict = "PASS" if n_pass >= 1 else "FAIL"

    log("=" * 60)
    log(f"PE PHASE 4 VERDICT: {verdict}  ({n_pass}/{n_total} guides cleared off-target gate)")
    log("=" * 60)

    out = {
        "phase": "PE-4",
        "method": "Cas-OFFinder on hg38 chr15, max 4 mismatches",
        "genome": str(GENOME_PATH),
        "max_mismatches": MAX_MISMATCHES,
        "gate_max_offtargets_at_3mm": PASS_GATE_MAX_OFFTARGETS_AT_3MM,
        "n_guides_total": n_total,
        "n_guides_pass": n_pass,
        "verdict": verdict,
        "results": all_results,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    log(f"Wrote {OUT_JSON}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
