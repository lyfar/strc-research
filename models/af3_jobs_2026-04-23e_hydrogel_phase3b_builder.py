#!/usr/bin/env python3
"""
Build AlphaFold Server submission batch 2026-04-23e — Hydrogel Phase 3b
full-construct validation.

Phase 3 tail retool (af3_jobs_2026-04-23c) returned:
  51 aa tail (1660-1710, clusters 4+5 only)    ipTM 0.23  FAIL
  71 aa tail (1640-1710, +cluster 3 loop cap)  ipTM 0.63  PASS
  91 aa tail (1620-1710, +cluster 2)           ipTM 0.68  PASS (best)

Phase 3b fuses the winning tail(s) back to the full RADA16-WH2 scaffold
architecture validated in Phase 2 (actin ipTM 0.58-0.59). Tests whether:
  1. TMEM145 binding holds at full-construct scale (~140-160 aa peptide
     vs isolated 91 aa tail)
  2. Actin binding preserved (WH2 motif still engages G-actin trimer)
  3. Self-assembly backbone (RADA16) doesn't disrupt either interface

Full construct = WH2 actin-binder (17 aa) + GSGSG (5 aa) + RADA16 (16 aa)
                 + GSGSG (5 aa) + extended tail (71 or 91 aa)

Two candidates submitted (both with RADA16 scaffold, winning tails):
  1. rada16_wh2_tail91  — 91 aa tail (ipTM 0.68 solo) — best-performer carry-forward
  2. rada16_wh2_tail71  — 71 aa tail (ipTM 0.63 solo) — manufacturability-optimum

Each tested × TMEM145 AND × G-actin trimer = 4 AF3 jobs.

Gate per pair: TMEM145 ipTM >= 0.50 AND actin ipTM >= 0.50. Both must
hold for #9 Hydrogel hypothesis to confirm A-tier and advance to Phase 2b
Martini3 MD + Phase 2c wet-lab actin-bundling assay.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23e_hydrogel_phase3b")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


STRC_FULL = fetch_uniprot("Q7RTU9")
ACTG1 = fetch_uniprot("P63261")
assert len(STRC_FULL) == 1775
assert len(ACTG1) == 375

TMEM145_FULL = (
    "MEPLRAPALRRLLPPLLLLLLSLPPRARAKYVRGNLSSKEDWVFLTRFCFLSDYGRLDFRFRYPEAKCCQNILLYFDDPSQWPAVYKAGDKDCLAKESVIRPE"
    "NNQVINLTTQYAWSGCQVVSEEGTRYLSCSSGRSFRSGDGLQLEYEMVLTNGKSFWTRHFSADEFGILETDVTFLLIFILIFFLSCYFGYLLKGRQLLHTTY"
    "KMFMAAAGVEVLSLLFFCIYWGQYATDGIGNESVKILAKLLFSSSFLIFLLMLILLGKGFTVTRGRISHAGSVKLSVYMTLYTLTHVVLLIYEAEFFDPGQV"
    "LYTYESPAGYGLIGLQVAAYVWFCYAVLVSLRHFPEKQPFYVPFFAAYTLWFFAVPVMALIANFGIPKWAREKIVNGIQLGIHLYAHGVFLIMTRPSAANKN"
    "FPYHVRTSQIASAGVPGPGGSQSADKAFPQHVYGNVTFISDSVPNFTELFSIPPPATSPLPRAAPDSGLPLFRDLRPPGPLRDL"
)
assert len(TMEM145_FULL) == 493

# Phase 2 construct components (from af3_jobs_2026-04-23b_hydrogel_builder.py):
# rada16_wh2_native = RQLVKAIPDNCSKSNVSR (N-term WH2, 18 aa) + ADARADARADARADA (RADA16, 15 aa) +
#                     AEDLPEPVPNCA (wrong-epitope tail, 12 aa) = 45 aa
# For Phase 3b we replace the 12 aa wrong-tail with the correct 71 or 91 aa GOLD-zone tail.
WH2_NTERM = "RQLVKAIPDNCSKSNVSR"
assert len(WH2_NTERM) == 18
RADA16 = "ADARADARADARADA"
assert len(RADA16) == 15
LINK_5 = "GSGSG"

# Correct GOLD-zone tails (validated in Phase 3)
TAIL_71 = STRC_FULL[1639:1710]  # 71 aa, ipTM 0.63 solo
TAIL_91 = STRC_FULL[1619:1710]  # 91 aa, ipTM 0.68 solo
assert len(TAIL_71) == 71
assert len(TAIL_91) == 91

# Full constructs: WH2 + LINK + RADA16 + LINK + TAIL
PEPTIDE_TAIL71 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_71
PEPTIDE_TAIL91 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_91
assert len(PEPTIDE_TAIL71) == 18 + 5 + 15 + 5 + 71 == 114
assert len(PEPTIDE_TAIL91) == 18 + 5 + 15 + 5 + 91 == 134


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


CONSTRUCTS = {
    "rada16_wh2_tail71": (PEPTIDE_TAIL71, "71 aa GOLD-zone tail (clusters 3+4+5, ipTM 0.63 solo)"),
    "rada16_wh2_tail91": (PEPTIDE_TAIL91, "91 aa GOLD-zone tail (clusters 2+3+4+5, ipTM 0.68 solo)"),
}

jobs = {}
for cname, (pseq, desc) in CONSTRUCTS.items():
    jobs[f"hydrogel_{cname}_x_tmem145"] = {
        "json": job(f"hydrogel_{cname}_x_tmem145", [
            (pseq, 1),
            (TMEM145_FULL, 1),
        ]),
        "description": (
            f"Full-construct peptide × TMEM145. "
            f"Peptide = WH2 (18) + GSGSG + RADA16 (15) + GSGSG + tail. "
            f"Length {len(pseq)} aa. Tail: {desc}. "
            f"Tests whether scaffold fusion preserves the 0.63-0.68 isolated-tail binding."
        ),
        "success_criteria": "ipTM >= 0.50 on TMEM145 interface",
        "hypothesis_impact": (
            "pass -> scaffold fusion preserves TMEM145 binding; half of Phase 3b gate clear. "
            "fail -> scaffold/linker contaminates the binding surface; try flexible linker "
            "or N-terminal tail placement."
        ),
    }
    jobs[f"hydrogel_{cname}_x_actin"] = {
        "json": job(f"hydrogel_{cname}_x_actin", [
            (pseq, 1),
            (ACTG1, 3),
        ]),
        "description": (
            f"Full-construct peptide × G-actin trimer. "
            f"Peptide {len(pseq)} aa. Tail: {desc}. "
            f"Tests whether extended tail (71/91 aa) disturbs WH2 actin-binding from N-terminus."
        ),
        "success_criteria": "ipTM >= 0.50 on any peptide-actin pair (matches Phase 2 0.58-0.59)",
        "hypothesis_impact": (
            "pass -> full dual-interface construct works; advance to Phase 2b Martini3 MD + "
            "Phase 2c wet-lab actin-bundling assay. "
            "fail -> extended tail disrupts WH2 allosterically; try larger linker between "
            "RADA16 and tail."
        ),
    }

for name, entry in jobs.items():
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(entry["json"], indent=2))

manifest = {
    "generated": "2026-04-23",
    "builder": "af3_jobs_2026-04-23e_hydrogel_phase3b_builder.py",
    "hypothesis": "#9 STRC Synthetic Peptide Hydrogel HTC — Phase 3b full-construct validation",
    "purpose": (
        "Phase 3b fuses the Phase 3 winning GOLD-zone tails back to the "
        "Phase 2 RADA16-WH2 scaffold architecture. Tests whether isolated-"
        "tail TMEM145 binding (ipTM 0.63-0.68) is preserved when the tail "
        "is attached to the WH2 actin-binder N-terminus + RADA16 self-"
        "assembling scaffold. Also tests whether extended tail disturbs "
        "WH2 actin engagement from the far N-terminus."
    ),
    "baselines_for_comparison": {
        "phase_2_native_full_construct_x_tmem145_WRONG_EPITOPE": 0.39,
        "phase_2_full_construct_x_actin": "0.58-0.59 (consistent across 3 candidates)",
        "phase_3_tail_71_solo_x_tmem145": 0.63,
        "phase_3_tail_91_solo_x_tmem145": 0.68,
        "ultramini_solo_x_tmem145_full": 0.43,
        "ultramini_solo_x_tmem145_GOLD_pruned": 0.68,
    },
    "construct_designs": {
        "rada16_wh2_tail71": {
            "length": len(PEPTIDE_TAIL71),
            "composition": "WH2 (18) + GSGSG + RADA16 (15) + GSGSG + tail aa 1640-1710 (71) = 114 aa",
        },
        "rada16_wh2_tail91": {
            "length": len(PEPTIDE_TAIL91),
            "composition": "WH2 (18) + GSGSG + RADA16 (15) + GSGSG + tail aa 1620-1710 (91) = 134 aa",
        },
    },
    "priority_order": list(jobs.keys()),
    "jobs": {name: {k: v for k, v in entry.items() if k != "json"}
             for name, entry in jobs.items()},
    "how_to_submit": (
        "AlphaFold Server Import-from-JSON for each file. 4 jobs this batch."
    ),
    "decision_rules": {
        "both_pass_both_interfaces": (
            "Hypothesis #9 A-tier CONFIRMED. Full construct validated at "
            "AF3 level. Advance to Phase 2b Martini3 CG MD dose-response "
            "sweep and Phase 2c wet-lab cell-based actin-bundling assay. "
            "Select tail91 for max confidence or tail71 for manufacturability "
            "(0.05 ipTM cost vs $200/ear cheaper)."
        ),
        "tmem145_pass_actin_fail": (
            "Extended tail allosterically disturbs WH2 from N-term. "
            "Try longer inter-domain linker (GSGSGSGSGS, 10 aa) or "
            "alternative WH2 placement (C-terminal instead of N-terminal)."
        ),
        "tmem145_fail_actin_pass": (
            "Scaffold fusion contaminates TMEM145 binding. Try flexible "
            "linker between RADA16 and tail (longer GSGSG variants) "
            "or C-terminal scaffold placement."
        ),
        "both_fail": (
            "Scaffold architecture doesn't support full dual-interface peptide. "
            "Move to Option B (helical mimetic / miniprotein scaffold replacing RADA16)."
        ),
    },
}
(OUT_DIR / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

fasta_lines = []
for cname, (pseq, desc) in CONSTRUCTS.items():
    fasta_lines += [f">{cname} | {len(pseq)} aa | {desc}", pseq]
fasta_lines += [f">TMEM145_FULL | {len(TMEM145_FULL)} aa | binding partner", TMEM145_FULL]
fasta_lines += [f">ACTG1 | {len(ACTG1)} aa | stereocilia gamma-actin", ACTG1]
(OUT_DIR / "af3_jobs_sequences.fasta").write_text("\n".join(fasta_lines) + "\n")

print(f"Generated AF3 Hydrogel Phase 3b full-construct batch in {OUT_DIR}")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size} B)")
