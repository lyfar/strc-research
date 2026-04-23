#!/usr/bin/env python3
"""
Build AlphaFold Server submission batch 2026-04-23 — SpyCatcher Assembly Phase 1b.

Sibling builder to af3_jobs_2026-04-22_builder.py. The 2026-04-22 batch split
STRC at aa 1074/1075 (Ultra-Mini boundary) and returned marginal MISS on both
gates (fold pTM 0.59 vs 0.60, binding ipTM 0.37 vs 0.40). The SpyCatcher proof
note identified the **Mini-STRC canonical domain boundary at aa 700** as the
alternative split point: aa 700 is a natural LRR-to-ARM hinge region with no
known tertiary contacts crossing it, so the 130-aa SpyCatcher/SpyTag insertion
should be less disruptive there than at aa 1074/1075.

Split is at aa 700/701:
  Fragment 1: STRC 1-700 (includes signal peptide + LRR region) + SpyCatcher
  Fragment 2: SpyTag + STRC 701-1775 (= Mini-STRC, 1075 aa, the clinical candidate)

Total reassembled chain: 700 + 5 + 114 + 2 + 13 + 5 + 1075 = 1914 aa (same as
aa 1074 split by construction — STRC fragment sum is always 1775).

Output files:
  af3_jobs_2026-04-23/
    MANIFEST.json
    af3_jobs_sequences.fasta
    strc_spy_reassembled_aa700_fold.json
    strc_spy_reassembled_aa700_x_tmem145.json

Submission: AlphaFold Server (alphafoldserver.com) → Import JSON, one per job.
Daily limit 20 jobs, this batch uses 2. Model seed 42 for reproducibility vs
prior STRC runs (same seed as 2026-04-22 batch).
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Source sequences
# -----------------------------------------------------------------------------

def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


STRC_FULL = fetch_uniprot("Q7RTU9")  # human stereocilin, 1775 aa
assert len(STRC_FULL) == 1775

# Split at aa 700 — Mini-STRC canonical domain boundary
STRC_N_1_700 = STRC_FULL[0:700]          # N-terminal fragment incl signal peptide + LRR, 700 aa
STRC_MINI = STRC_FULL[700:1775]          # Mini-STRC (clinical candidate), 1075 aa
assert len(STRC_N_1_700) == 700
assert len(STRC_MINI) == 1075

# TMEM145 — re-use verified sequence from prior AF3 batches
TMEM145_FULL = (
    "MEPLRAPALRRLLPPLLLLLLSLPPRARAKYVRGNLSSKEDWVFLTRFCFLSDYGRLDFRFRYPEAKCCQNILLYFDDPSQWPAVYKAGDKDCLAKESVIRPE"
    "NNQVINLTTQYAWSGCQVVSEEGTRYLSCSSGRSFRSGDGLQLEYEMVLTNGKSFWTRHFSADEFGILETDVTFLLIFILIFFLSCYFGYLLKGRQLLHTTY"
    "KMFMAAAGVEVLSLLFFCIYWGQYATDGIGNESVKILAKLLFSSSFLIFLLMLILLGKGFTVTRGRISHAGSVKLSVYMTLYTLTHVVLLIYEAEFFDPGQV"
    "LYTYESPAGYGLIGLQVAAYVWFCYAVLVSLRHFPEKQPFYVPFFAAYTLWFFAVPVMALIANFGIPKWAREKIVNGIQLGIHLYAHGVFLIMTRPSAANKN"
    "FPYHVRTSQIASAGVPGPGGSQSADKAFPQHVYGNVTFISDSVPNFTELFSIPPPATSPLPRAAPDSGLPLFRDLRPPGPLRDL"
)
assert len(TMEM145_FULL) == 493

# -----------------------------------------------------------------------------
# SpyCatcher/SpyTag
# -----------------------------------------------------------------------------

SPYCATCHER = (
    "MVDTLSGLSSEQGQSGDMTIEEDSATHIKFSKRDEDGKELAGATMELRDSSGKTISTWISDGQ"
    "VKDFYLYPGKYTFVETAAPDGYEVATAITFTVNEQGQVTVNGKATKGDAHI"
)
SPYTAG = "AHIVMVDAYKPTK"
assert len(SPYCATCHER) == 114
assert len(SPYTAG) == 13

LINK_5 = "GSGSG"  # flexible 5 aa linker
LINK_2 = "GG"     # 2 aa bridge at the pseudo-isopeptide seam

# Reassembled construct with alternate split
STRC_SPY_REASSEMBLED_AA700 = (
    STRC_N_1_700
    + LINK_5  # between STRC-N and SpyCatcher
    + SPYCATCHER
    + LINK_2  # pseudo-isopeptide seam
    + SPYTAG
    + LINK_5  # between SpyTag and STRC-C (Mini-STRC)
    + STRC_MINI
)
assert len(STRC_SPY_REASSEMBLED_AA700) == 700 + 5 + 114 + 2 + 13 + 5 + 1075 == 1914

# -----------------------------------------------------------------------------
# AF3 JSON builder
# -----------------------------------------------------------------------------

def job(name: str, chains: list[tuple[str, int]]) -> list[dict]:
    return [{
        "name": name,
        "modelSeeds": ["42"],
        "sequences": [
            {"proteinChain": {"sequence": seq, "count": count}}
            for seq, count in chains
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }]


jobs = {
    "strc_spy_reassembled_aa700_fold": {
        "json": job("strc_spy_reassembled_aa700_fold", [(STRC_SPY_REASSEMBLED_AA700, 1)]),
        "description": (
            "Phase 1b alternate-split SpyCatcher reassembly: STRC 1-700 + "
            "GSGSG + SpyCatcher (114) + GG + SpyTag + GSGSG + STRC 701-1775 "
            f"(Mini-STRC). {len(STRC_SPY_REASSEMBLED_AA700)} aa total. "
            "Split at aa 700 = Mini-STRC canonical LRR-to-ARM domain boundary."
        ),
        "success_criteria": (
            "pTM >= 0.60 for overall chain. Compared to 2026-04-22 batch which "
            "missed at pTM 0.59, pass at 0.60+ means the aa 700 boundary "
            "tolerates the SpyCatcher insertion better than the aa 1074 "
            "boundary did. C-terminal Mini-STRC region (chain positions "
            "840-1914) should adopt the ARM-repeat fold."
        ),
        "hypothesis_impact": (
            "pass -> SpyCatcher reassembly viable at canonical domain boundary; "
            "[[STRC In Situ SpyCatcher Assembly]] B -> A (evidence depth +2, "
            "backup to Mini-STRC single-vector). "
            "fail (pTM < 0.55) -> SpyCatcher insertion fundamentally "
            "incompatible with STRC fold regardless of split point; "
            "hypothesis demoted to C."
        ),
    },
    "strc_spy_reassembled_aa700_x_tmem145": {
        "json": job("strc_spy_reassembled_aa700_x_tmem145", [
            (STRC_SPY_REASSEMBLED_AA700, 1),
            (TMEM145_FULL, 1),
        ]),
        "description": (
            "aa700-split reassembled construct + TMEM145 full. "
            f"{len(STRC_SPY_REASSEMBLED_AA700)} + {len(TMEM145_FULL)} = "
            f"{len(STRC_SPY_REASSEMBLED_AA700) + len(TMEM145_FULL)} aa. "
            "Phase 1b: does the reassembled Mini-STRC preserve TMEM145 binding?"
        ),
        "success_criteria": (
            "ipTM >= 0.40 (matches Ultra-Mini baseline). STRC chain interface "
            "residues concentrate around native aa 1603-1749 (= chain positions "
            "1742-1888 in the reassembled product). PAE min between chains <= 9 A."
        ),
        "hypothesis_impact": (
            "pass -> full hypothesis green; Phase 2 cell-based expression + "
            "SpyCatcher reaction assay (~$5-10k, 8-12 weeks). "
            "fail -> even the canonical domain boundary can't preserve "
            "binding with SpyCatcher architecture; kill hypothesis."
        ),
    },
}

# -----------------------------------------------------------------------------
# Write files
# -----------------------------------------------------------------------------

for name, entry in jobs.items():
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(entry["json"], indent=2))

manifest = {
    "generated": "2026-04-23",
    "builder": "af3_jobs_2026-04-23_builder.py",
    "purpose": (
        "SpyCatcher Assembly Phase 1b with alternate split at aa 700 "
        "(Mini-STRC canonical LRR-to-ARM domain boundary). "
        "Prior 2026-04-22 batch split at aa 1074/1075 (Ultra-Mini boundary) "
        "and returned marginal MISS on both fold (pTM 0.59 vs 0.60) and "
        "binding (ipTM 0.37 vs 0.40). The escape path documented in "
        "[[STRC SpyCatcher Assembly Phase 1 Geometry]] is the aa 700 "
        "canonical domain boundary: structural biology literature suggests "
        "this hinge region has no cross-domain tertiary contacts, so the "
        "130-aa SpyCatcher/SpyTag insertion should be better tolerated."
    ),
    "source_sequences": {
        "STRC_full": "UniProt Q7RTU9, 1775 aa (fetched 2026-04-23)",
        "TMEM145_full": "verified from prior AF3 batches",
        "SpyCatcher": "original SpyCatcher (Zakeri 2012), 114 aa",
        "SpyTag": "original SpyTag, 13 aa",
    },
    "construct_assembly": {
        "STRC_SPY_REASSEMBLED_AA700": {
            "length": len(STRC_SPY_REASSEMBLED_AA700),
            "composition": (
                "STRC_1-700 (700) + GSGSG + SpyCatcher (114) + GG + "
                "SpyTag (13) + GSGSG + STRC_701-1775 (1075, = Mini-STRC) "
                "= 1914 aa"
            ),
            "rationale": (
                "Mini-STRC (aa 700-1775) is the clinical candidate construct "
                "(see [[STRC Mini-STRC Truncation Interface Validation]]). "
                "The aa 700 boundary is a natural LRR-to-ARM domain hinge; "
                "splitting there isolates the SpyCatcher insertion from "
                "both folded subdomains. Compare to aa 1074/1075 split which "
                "bisects the ARM-repeat region internally."
            ),
        },
    },
    "priority_order": [
        "strc_spy_reassembled_aa700_fold",
        "strc_spy_reassembled_aa700_x_tmem145",
    ],
    "jobs": {name: {k: v for k, v in entry.items() if k != "json"}
             for name, entry in jobs.items()},
    "how_to_submit": (
        "AlphaFold Server (alphafoldserver.com): Import from JSON for each "
        "file. Daily limit 20 jobs, this batch uses 2. Seed 42 for "
        "reproducibility with 2026-04-22 batch (same seed allows direct "
        "pTM/ipTM delta comparison)."
    ),
    "decision_rules": {
        "fold_pass_binding_pass": (
            "SpyCatcher hypothesis B -> A (evidence depth +2); "
            "proceed to Phase 2 wet-lab."
        ),
        "fold_pass_binding_fail": (
            "Insertion is geometrically OK but interrupts TMEM145 recognition; "
            "SpyCatcher B -> C. No further computational escape paths."
        ),
        "fold_fail": (
            "SpyCatcher architecture fundamentally incompatible with STRC; "
            "hypothesis B -> D (killed). Focus on Mini-STRC single-vector."
        ),
        "both_marginal_again": (
            "If both gates miss by <=0.03 again, the method is at its confidence "
            "noise floor; advance to MD-based assessment or wet-lab."
        ),
    },
}
(OUT_DIR / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

fasta_lines = [
    f">STRC_SPY_REASSEMBLED_AA700 | {len(STRC_SPY_REASSEMBLED_AA700)} aa | hypothesis #10 Phase 1b alternate-split reassembly",
    STRC_SPY_REASSEMBLED_AA700,
    f">TMEM145_FULL | {len(TMEM145_FULL)} aa | binding partner",
    TMEM145_FULL,
]
(OUT_DIR / "af3_jobs_sequences.fasta").write_text("\n".join(fasta_lines) + "\n")

print(f"Generated AF3 Phase 1b batch in {OUT_DIR}")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size} B)")
