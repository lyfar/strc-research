#!/usr/bin/env python3
"""
Build AlphaFold Server submission batch 2026-04-22 for two new Phase 1 tests:
hypotheses #10 (SpyCatcher Assembly) and #11 (TECTA Chimera).

Both hypotheses sit in B-tier with ZERO compute done. Moving them with one
lightweight AF3 pass either (a) validates the construct geometry so they earn
Phase 2 compute, or (b) kills them cheaply.

Output files:
  af3_jobs_2026-04-22/
    MANIFEST.json
    af3_jobs_sequences.fasta
    strc_spy_reassembled_fold.json
    strc_spy_reassembled_x_tmem145.json
    strc_tecta_chimera_fold.json
    strc_tecta_chimera_x_tmem145.json

Submission: AlphaFold Server (alphafoldserver.com) → Import JSON, one per job.
Daily limit 20 jobs. Model seed 42 for reproducibility vs prior STRC runs.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-22")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Source sequences
# -----------------------------------------------------------------------------

def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


STRC_FULL = fetch_uniprot("Q7RTU9")   # human stereocilin, 1775 aa
TECTA_FULL = fetch_uniprot("O75443")  # human alpha-tectorin, 2155 aa
assert len(STRC_FULL) == 1775
assert len(TECTA_FULL) == 2155

STRC_N_1_1074 = STRC_FULL[0:1074]           # N-terminal fragment for SpyCatcher
STRC_ULTRAMINI = STRC_FULL[1074:1775]       # Ultra-Mini 1075-1775, 701 aa
TECTA_ZP_REGION = TECTA_FULL[1654:2155]     # C-terminal 500 aa, spans ZP-N + ZP-C + CFCS

# TMEM145 — re-use verified sequences from prior AF3 batch.
TMEM145_FULL = (
    "MEPLRAPALRRLLPPLLLLLLSLPPRARAKYVRGNLSSKEDWVFLTRFCFLSDYGRLDFRFRYPEAKCCQNILLYFDDPSQWPAVYKAGDKDCLAKESVIRPE"
    "NNQVINLTTQYAWSGCQVVSEEGTRYLSCSSGRSFRSGDGLQLEYEMVLTNGKSFWTRHFSADEFGILETDVTFLLIFILIFFLSCYFGYLLKGRQLLHTTY"
    "KMFMAAAGVEVLSLLFFCIYWGQYATDGIGNESVKILAKLLFSSSFLIFLLMLILLGKGFTVTRGRISHAGSVKLSVYMTLYTLTHVVLLIYEAEFFDPGQV"
    "LYTYESPAGYGLIGLQVAAYVWFCYAVLVSLRHFPEKQPFYVPFFAAYTLWFFAVPVMALIANFGIPKWAREKIVNGIQLGIHLYAHGVFLIMTRPSAANKN"
    "FPYHVRTSQIASAGVPGPGGSQSADKAFPQHVYGNVTFISDSVPNFTELFSIPPPATSPLPRAAPDSGLPLFRDLRPPGPLRDL"
)
assert len(TMEM145_FULL) == 493

# -----------------------------------------------------------------------------
# Construct assemblies
# -----------------------------------------------------------------------------

# SpyCatcher (original, 113 aa) + SpyTag (13 aa). After reaction they form a
# covalent isopeptide bond (SpyCatcher K31 — SpyTag D7). AF3 can't model the
# isopeptide explicitly, but a single-chain representation where SpyCatcher
# and SpyTag sit back-to-back with a 2 aa pseudo-bond gives AF3 the assembly
# topology and forces it to find a fold consistent with covalent linkage.
SPYCATCHER = (
    "MVDTLSGLSSEQGQSGDMTIEEDSATHIKFSKRDEDGKELAGATMELRDSSGKTISTWISDGQ"
    "VKDFYLYPGKYTFVETAAPDGYEVATAITFTVNEQGQVTVNGKATKGDAHI"
)
SPYTAG = "AHIVMVDAYKPTK"
assert len(SPYCATCHER) == 114
assert len(SPYTAG) == 13

LINK_5 = "GSGSG"  # flexible 5 aa linker
LINK_2 = "GG"     # 2 aa bridge at the pseudo-isopeptide seam

# Hypothesis #10 reassembled construct (single chain, models the post-reaction product)
STRC_SPY_REASSEMBLED = (
    STRC_N_1_1074
    + LINK_5  # between STRC-N and SpyCatcher
    + SPYCATCHER
    + LINK_2  # pseudo-isopeptide seam
    + SPYTAG
    + LINK_5  # between SpyTag and STRC-C
    + STRC_ULTRAMINI
)
assert len(STRC_SPY_REASSEMBLED) == 1074 + 5 + 113 + 2 + 13 + 5 + 701 == 1913

# Hypothesis #11 TECTA chimera: TECTA ZP region as C-terminal anchor +
# GGSGSG linker + STRC Ultra-Mini as the N-terminal functional domain.
TECTA_CHIMERA = TECTA_ZP_REGION + "GGSGSG" + STRC_ULTRAMINI
assert len(TECTA_CHIMERA) == 500 + 6 + 701 == 1207

# -----------------------------------------------------------------------------
# AF3 JSON builder
# -----------------------------------------------------------------------------

def job(name: str, chains: list[tuple[str, int]]) -> list[dict]:
    """Single-job list for AlphaFold Server format."""
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
    "strc_spy_reassembled_fold": {
        "json": job("strc_spy_reassembled_fold", [(STRC_SPY_REASSEMBLED, 1)]),
        "description": (
            "Single-chain simulation of post-assembly SpyCatcher product: "
            f"STRC N-term 1-1074 + linker + SpyCatcher + SpyTag + linker + "
            f"STRC 1075-1775. {len(STRC_SPY_REASSEMBLED)} aa total. "
            "Purpose: does the reassembled product adopt a native-like STRC fold?"
        ),
        "success_criteria": (
            "pTM >= 0.60 for overall chain AND the C-terminal STRC region "
            "(1075-1775 equivalent positions, aa 1194-1913 of this chain) "
            "adopts an ARM-repeat fold comparable to Ultra-Mini solo "
            "(previous job-h-strc-cterm-only.cif at pTM 0.87)"
        ),
        "hypothesis_impact": (
            "pass -> SpyCatcher reassembly is geometrically viable; advance "
            "to Phase 1b TMEM145 binding test. "
            "fail (pTM < 0.40 OR C-term ARM fold destroyed) -> SpyCatcher "
            "tag topology incompatible with STRC C-term tertiary structure; "
            "kill or try alternate split points"
        ),
    },
    "strc_spy_reassembled_x_tmem145": {
        "json": job("strc_spy_reassembled_x_tmem145", [
            (STRC_SPY_REASSEMBLED, 1),
            (TMEM145_FULL, 1),
        ]),
        "description": (
            "Reassembled SpyCatcher construct + TMEM145 full length. "
            f"{len(STRC_SPY_REASSEMBLED)} + {len(TMEM145_FULL)} = 2406 aa. "
            "Phase 1b: does the reassembled STRC preserve the TMEM145 interface?"
        ),
        "success_criteria": (
            "ipTM >= 0.40 (matches prior Ultra-Mini x TMEM145 baseline). "
            "STRC chain interface residues concentrate in positions equivalent "
            "to aa 1603-1749 of native STRC (= positions 1603-1749 within the "
            "C-terminal segment of the reassembled chain, indexed from SpyTag+5)"
        ),
        "hypothesis_impact": (
            "pass -> SpyCatcher fragment complementation preserves TMEM145 "
            "binding. Hypothesis #10 earns Phase 2 (fragment expression + "
            "in vitro SpyCatcher-SpyTag reaction test). "
            "fail -> TMEM145 interface disrupted; kill"
        ),
    },
    "strc_tecta_chimera_fold": {
        "json": job("strc_tecta_chimera_fold", [(TECTA_CHIMERA, 1)]),
        "description": (
            f"TECTA ZP region (aa 1655-2155, 500 aa) + linker + STRC Ultra-Mini "
            f"(1075-1775, 701 aa). {len(TECTA_CHIMERA)} aa. "
            "Purpose: does the chimera adopt a stable two-domain fold?"
        ),
        "success_criteria": (
            "pTM >= 0.55 overall. TECTA ZP domain region (aa 1-500 of chimera) "
            "adopts ZP-fold. STRC Ultra-Mini region (aa 507-1207 of chimera) "
            "adopts ARM-repeat fold comparable to Ultra-Mini solo."
        ),
        "hypothesis_impact": (
            "pass -> chimera geometry plausible; advance to Phase 1b TMEM145 "
            "test. fail (domains clash or linker forces non-native fold) -> "
            "TECTA-STRC fusion incompatible; try alternate TECTA region "
            "(e.g. VWF D domain) or kill"
        ),
    },
    "strc_tecta_chimera_x_tmem145": {
        "json": job("strc_tecta_chimera_x_tmem145", [
            (TECTA_CHIMERA, 1),
            (TMEM145_FULL, 1),
        ]),
        "description": (
            f"TECTA-STRC chimera + TMEM145 full. "
            f"{len(TECTA_CHIMERA)} + {len(TMEM145_FULL)} = 1700 aa. "
            "Phase 1b: does chimera preserve TMEM145 binding despite TECTA scaffold?"
        ),
        "success_criteria": (
            "ipTM >= 0.40. STRC region of chimera engages TMEM145 interface; "
            "TECTA ZP region does NOT dominate the interface "
            "(false-positive binder)"
        ),
        "hypothesis_impact": (
            "pass -> TECTA-STRC chimera preserves the STRC-TMEM145 recognition "
            "surface while gaining ZP anchoring; advance to Phase 2 "
            "(tectorial membrane self-assembly model). "
            "fail -> chimera either breaks TMEM145 binding OR TECTA ZP "
            "dominates interface artifactually; kill or redesign"
        ),
    },
}

# -----------------------------------------------------------------------------
# Write files
# -----------------------------------------------------------------------------

for name, entry in jobs.items():
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(entry["json"], indent=2))

# Manifest
manifest = {
    "generated": "2026-04-22",
    "builder": "af3_jobs_2026-04-22_builder.py",
    "purpose": (
        "Phase 1 geometry validation for two B-tier STRC hypotheses currently "
        "at 0 compute: #10 [[STRC In Situ SpyCatcher Assembly]] and #11 "
        "[[STRC Engineered TECTA Chimera]]. Both sit in B-tier partly because "
        "no structural evidence exists either way. A single AF3 pass moves "
        "them to Phase 2 (if fold + TMEM145 binding preserved) or kills them "
        "cheaply."
    ),
    "source_sequences": {
        "STRC_full": "UniProt Q7RTU9, 1775 aa (fetched 2026-04-22)",
        "TECTA_full": "UniProt O75443, 2155 aa (fetched 2026-04-22)",
        "TMEM145_full": "verified from prior AF3 batch 2026-04-21 FASTA",
        "SpyCatcher": "original SpyCatcher (Zakeri 2012), 113 aa",
        "SpyTag": "original SpyTag, 13 aa",
    },
    "construct_assemblies": {
        "STRC_SPY_REASSEMBLED": {
            "length": len(STRC_SPY_REASSEMBLED),
            "composition": (
                "STRC_1-1074 (1074) + GSGSG + SpyCatcher (113) + GG + "
                "SpyTag (13) + GSGSG + STRC_1075-1775 (701) = 1913 aa"
            ),
            "rationale": (
                "Models the post-reaction covalent product of Fragment 1 "
                "(STRC-N + SpyCatcher) meeting Fragment 2 (SpyTag + "
                "STRC-C) in the secretory pathway or on the OHC surface."
            ),
        },
        "TECTA_CHIMERA": {
            "length": len(TECTA_CHIMERA),
            "composition": (
                "TECTA_1655-2155 (500, ZP region) + GGSGSG + "
                "STRC_1075-1775 (701) = 1207 aa"
            ),
            "rationale": (
                "TECTA ZP C-terminal anchor scaffold + STRC Ultra-Mini "
                "functional domain. Idea: TECTA ZP auto-polymerisation "
                "in tectorial membrane provides the anchoring, STRC "
                "Ultra-Mini provides the TMEM145-binding horizontal "
                "connector function."
            ),
        },
    },
    "priority_order": [
        "strc_spy_reassembled_fold",
        "strc_tecta_chimera_fold",
        "strc_spy_reassembled_x_tmem145",
        "strc_tecta_chimera_x_tmem145",
    ],
    "jobs": {name: {k: v for k, v in entry.items() if k != "json"}
             for name, entry in jobs.items()},
    "how_to_submit": (
        "AlphaFold Server (alphafoldserver.com): Import from JSON for each "
        "file. Daily limit 20 jobs, we use 4. Seed 42 for reproducibility."
    ),
}
(OUT_DIR / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

# FASTA for human reading
fasta_lines = [
    f">STRC_SPY_REASSEMBLED | {len(STRC_SPY_REASSEMBLED)} aa | hypothesis #10 post-assembly simulation",
    STRC_SPY_REASSEMBLED,
    f">TECTA_CHIMERA | {len(TECTA_CHIMERA)} aa | hypothesis #11 TECTA-STRC fusion",
    TECTA_CHIMERA,
    f">TMEM145_FULL | {len(TMEM145_FULL)} aa | partner for both binding tests",
    TMEM145_FULL,
]
(OUT_DIR / "af3_jobs_sequences.fasta").write_text("\n".join(fasta_lines) + "\n")

print(f"Generated AF3 batch in {OUT_DIR}")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size} B)")
