#!/usr/bin/env python3
"""
Build AlphaFold Server submission batch 2026-04-23b — Hydrogel Phase 2 AF3.

Phase 1 analytical pre-screen ([[STRC Hydrogel HTC Phase 1 Self-Assembly]])
passed 5/6 SAP candidates against four geometric gates. Top-3 advance to
Phase 2 AF3 validation of the two dual-interface binding handles:

  1. TMEM145-binding tail (derived from STRC aa 1669-1680 contact cluster)
  2. Actin-binding N-terminus (WH2 motif, binds G-actin monomers)

AF3 Phase 2 gate: ipTM >= 0.50 on BOTH interfaces per candidate.

Top-3 peptides (from hydrogel_phase1_self_assembly.json, ranked by composite):

  1. RADA16-WH2-native (45 aa, score 70.6)
  2. RADA16-WH2-Cmut   (44 aa, score 68.9)
  3. EAK16-WH2-denovo  (44 aa, score 68.9)

Jobs (6 total):
  - 3 x peptide x TMEM145 full (493 aa)
  - 3 x peptide x G-actin trimer (ACTG1, 3 x 375 = 1125 aa)

G-actin trimer chosen over full F-actin filament (5-7 monomers) for AF3 context
economy; WH2 biology is G-actin sequestration, so a trimer captures both the
binding pose and the nearest-neighbor filament geometry that matters for
stereocilium-bundle bridging.

Output files:
  af3_jobs_2026-04-23b_hydrogel/
    MANIFEST.json
    af3_jobs_sequences.fasta
    hydrogel_rada16_wh2_native_x_tmem145.json
    hydrogel_rada16_wh2_cmut_x_tmem145.json
    hydrogel_eak16_wh2_denovo_x_tmem145.json
    hydrogel_rada16_wh2_native_x_actin.json
    hydrogel_rada16_wh2_cmut_x_actin.json
    hydrogel_eak16_wh2_denovo_x_actin.json

Submission: alphafoldserver.com Import-from-JSON. Daily limit 20; seed 42.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23b_hydrogel")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


# -----------------------------------------------------------------------------
# Binding partners
# -----------------------------------------------------------------------------

# TMEM145 — verified from prior AF3 batches
TMEM145_FULL = (
    "MEPLRAPALRRLLPPLLLLLLSLPPRARAKYVRGNLSSKEDWVFLTRFCFLSDYGRLDFRFRYPEAKCCQNILLYFDDPSQWPAVYKAGDKDCLAKESVIRPE"
    "NNQVINLTTQYAWSGCQVVSEEGTRYLSCSSGRSFRSGDGLQLEYEMVLTNGKSFWTRHFSADEFGILETDVTFLLIFILIFFLSCYFGYLLKGRQLLHTTY"
    "KMFMAAAGVEVLSLLFFCIYWGQYATDGIGNESVKILAKLLFSSSFLIFLLMLILLGKGFTVTRGRISHAGSVKLSVYMTLYTLTHVVLLIYEAEFFDPGQV"
    "LYTYESPAGYGLIGLQVAAYVWFCYAVLVSLRHFPEKQPFYVPFFAAYTLWFFAVPVMALIANFGIPKWAREKIVNGIQLGIHLYAHGVFLIMTRPSAANKN"
    "FPYHVRTSQIASAGVPGPGGSQSADKAFPQHVYGNVTFISDSVPNFTELFSIPPPATSPLPRAAPDSGLPLFRDLRPPGPLRDL"
)
assert len(TMEM145_FULL) == 493

# G-actin (cytoplasmic gamma-1, ACTG1, P63261) — main stereocilia core actin
ACTG1 = fetch_uniprot("P63261")
assert len(ACTG1) == 375

# -----------------------------------------------------------------------------
# Top-3 peptide candidates (from hydrogel_phase1_self_assembly.json)
# -----------------------------------------------------------------------------

PEPTIDES = {
    "rada16_wh2_native": "RQLVKAIPDNCSKSNVSRADARADARADARADAAEDLPEPVPNCA",  # 45 aa
    "rada16_wh2_cmut":   "RQLVKAIPDNCSKSNVSRADARADARADARADAAEDLPEPVPNA",   # 44 aa
    "eak16_wh2_denovo":  "RQLVKAIPDNCSKSNVSAEAKAEAKAEAKAEAKEELPEPVPNYK",   # 44 aa
}
assert len(PEPTIDES["rada16_wh2_native"]) == 45
assert len(PEPTIDES["rada16_wh2_cmut"]) == 44
assert len(PEPTIDES["eak16_wh2_denovo"]) == 44


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


# -----------------------------------------------------------------------------
# Jobs: 3 x (peptide x TMEM145) + 3 x (peptide x G-actin trimer)
# -----------------------------------------------------------------------------

jobs = {}

for pname, pseq in PEPTIDES.items():
    jobs[f"hydrogel_{pname}_x_tmem145"] = {
        "json": job(f"hydrogel_{pname}_x_tmem145", [
            (pseq, 1),
            (TMEM145_FULL, 1),
        ]),
        "description": (
            f"Peptide {pname} ({len(pseq)} aa) + TMEM145 full (493 aa) = "
            f"{len(pseq) + 493} aa. Tests whether the TMEM145-binding tail "
            "(derived from STRC aa 1669-1680 contact cluster) retains "
            "recognition when fused to the RADA16/EAK16 self-assembling "
            "scaffold plus WH2 actin-binder."
        ),
        "success_criteria": (
            "ipTM >= 0.50 (Phase 2 gate, tighter than Phase 1 0.40 baseline). "
            "Peptide C-terminal tail residues (last ~12 aa, positions "
            f"{len(pseq)-11}-{len(pseq)}) drive the interface. TMEM145 "
            "contact residues cluster in the same zone as native STRC contacts "
            "(cross-referenced against [[STRC Ultra-Mini Full-Length TMEM145 AF3]] "
            "GOLD zone)."
        ),
        "hypothesis_impact": (
            "pass -> TMEM145-binding handle is functional at real peptide length; "
            "one of two interfaces validated. "
            "fail -> short-peptide Kd assumption (100 nM in Phase 1 Gate 3) "
            "falsified; Gate 3 f=0.70 estimate collapses; revisit design "
            "(longer tail, structured helical mimetic, or alternative TMEM145 "
            "epitope from STRC GOLD zone)."
        ),
    }

    jobs[f"hydrogel_{pname}_x_actin"] = {
        "json": job(f"hydrogel_{pname}_x_actin", [
            (pseq, 1),
            (ACTG1, 3),
        ]),
        "description": (
            f"Peptide {pname} ({len(pseq)} aa) + G-actin trimer (3 x 375 aa ACTG1) = "
            f"{len(pseq) + 3 * 375} aa. Tests whether the WH2 motif N-terminus "
            "retains actin binding when fused to the self-assembling scaffold. "
            "G-actin trimer captures near-filament geometry without the full "
            "F-actin AF3 context cost."
        ),
        "success_criteria": (
            "ipTM >= 0.50 on at least one peptide-actin pair. WH2 motif "
            "(peptide positions 1-17) engages the actin barbed-end cleft "
            "(actin subdomains 1 and 3, residues around D1/T148/D286)."
        ),
        "hypothesis_impact": (
            "pass -> actin-binding handle functional; one of two interfaces "
            "validated. "
            "fail -> WH2 fusion to RADA16/EAK16 scaffold disrupts actin "
            "recognition; peptide cannot localize to stereocilia cores; "
            "hypothesis architecturally broken unless an alternative "
            "actin-binder (espin mini, fascin, filamin-derived peptide) "
            "is engineered in."
        ),
    }

for name, entry in jobs.items():
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(entry["json"], indent=2))

manifest = {
    "generated": "2026-04-23",
    "builder": "af3_jobs_2026-04-23b_hydrogel_builder.py",
    "purpose": (
        "Hydrogel Phase 2 AF3-Multimer validation of top-3 SAP candidates "
        "from Phase 1 analytical pre-screen. Tests both dual-interface "
        "binding handles (TMEM145 + actin) at realistic peptide length, "
        "converting Phase 1 Kd=100 nM assumption from analytical hand-wave "
        "to structural evidence."
    ),
    "source_sequences": {
        "TMEM145_full": "verified from prior AF3 batches, 493 aa",
        "ACTG1": "UniProt P63261, cytoplasmic gamma-1 actin, 375 aa, "
                 "main stereocilia isoform (fetched 2026-04-23)",
        "Peptides": "top-3 from hydrogel_phase1_self_assembly.json results_ranked",
    },
    "peptides_tested": [
        {"name": "rada16_wh2_native", "aa": 45, "score": 70.6,
         "sequence": PEPTIDES["rada16_wh2_native"]},
        {"name": "rada16_wh2_cmut",   "aa": 44, "score": 68.9,
         "sequence": PEPTIDES["rada16_wh2_cmut"]},
        {"name": "eak16_wh2_denovo",  "aa": 44, "score": 68.9,
         "sequence": PEPTIDES["eak16_wh2_denovo"]},
    ],
    "priority_order": [
        "hydrogel_rada16_wh2_native_x_tmem145",
        "hydrogel_rada16_wh2_cmut_x_tmem145",
        "hydrogel_eak16_wh2_denovo_x_tmem145",
        "hydrogel_rada16_wh2_native_x_actin",
        "hydrogel_rada16_wh2_cmut_x_actin",
        "hydrogel_eak16_wh2_denovo_x_actin",
    ],
    "jobs": {name: {k: v for k, v in entry.items() if k != "json"}
             for name, entry in jobs.items()},
    "how_to_submit": (
        "AlphaFold Server (alphafoldserver.com): Import from JSON for each "
        "file. Daily limit 20; this batch uses 6. Seed 42 for consistency "
        "with 2026-04-22 and 2026-04-23 STRC batches."
    ),
    "decision_rules": {
        "both_interfaces_pass_for_any_candidate": (
            "Hypothesis #9 [[STRC Synthetic Peptide Hydrogel HTC]]: B -> A "
            "(evidence depth +2). Advance to Phase 2b Martini3 MD dose-response "
            "sweep and Phase 2c wet-lab cell-level actin-bundling assay."
        ),
        "tmem145_pass_actin_fail": (
            "B stays B. TMEM145 handle works, but the actin-bundling "
            "architecture is broken. Redesign N-terminus before Phase 2b."
        ),
        "tmem145_fail_actin_pass": (
            "B stays B. Short peptide TMEM145 Kd is the real gate "
            "(Phase 1 Gate 3 sensitivity confirmed). Retool tail."
        ),
        "both_fail": (
            "B -> C. Dual-interface short-peptide architecture is not "
            "structurally competent at AF3 confidence; hypothesis needs "
            "structured-mimetic approach (miniprotein, cyclic peptide) "
            "not SAP."
        ),
    },
}
(OUT_DIR / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

fasta_lines = []
for pname, pseq in PEPTIDES.items():
    fasta_lines += [f">{pname} | {len(pseq)} aa | hydrogel Phase 2 candidate", pseq]
fasta_lines += [f">TMEM145_FULL | {len(TMEM145_FULL)} aa | binding partner", TMEM145_FULL]
fasta_lines += [f">ACTG1 | {len(ACTG1)} aa | stereocilia gamma-actin", ACTG1]
(OUT_DIR / "af3_jobs_sequences.fasta").write_text("\n".join(fasta_lines) + "\n")

print(f"Generated AF3 Hydrogel Phase 2 batch in {OUT_DIR}")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size} B)")
