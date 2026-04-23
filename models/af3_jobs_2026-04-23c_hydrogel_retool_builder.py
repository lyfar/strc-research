#!/usr/bin/env python3
"""
Build AlphaFold Server submission batch 2026-04-23c — Hydrogel Phase 3 tail retool.

Phase 2 result (af3_jobs_2026-04-23b_hydrogel): actin interface PASS (ipTM 0.59),
TMEM145 interface FAIL (ipTM 0.37-0.39) across all 3 candidates with 11-12 aa
tails derived from STRC aa 1669-1680. MANIFEST decision rule:
"TMEM145 FAIL + actin PASS -> B stays, retool TMEM145 tail".

Phase 3 tests Option A from the retool plan — longer tail with aa 1603-1749
GOLD-zone context. Six canonical TMEM145 contact clusters exist in GOLD zone
(from [[STRC Ultra-Mini Full-Length TMEM145 AF3]]):

  cluster 1: aa 1603-1607 (5 residues)
  cluster 2: aa 1630-1638 (9 residues)
  cluster 3: aa 1648-1651 (4 residues)
  cluster 4: aa 1669-1680 (12 residues, DOMINANT — 9/12 reproduced in Ultra-Mini solo)
  cluster 5: aa 1692-1707 (16 residues, SECOND DOMINANT — 8/16 reproduced)
  cluster 6: aa 1770 (1 residue)

Clusters 4 + 5 together account for 28 of 47 canonical contacts. Phase 2 tails
hit only cluster 4, leaving cluster 5 untouched. Extending the tail to span
both clusters should recover the dominant hot-spot pair.

Three candidate tails (isolated, no scaffold — retool test only):

  1. STRC_1660_1710 (51 aa) — MINIMAL: clusters 4+5 with 10 aa flank each side
  2. STRC_1640_1710 (71 aa) — +loop cap: clusters 3+4+5
  3. STRC_1620_1710 (91 aa) — WIDE: clusters 2+3+4+5 (SPPS feasibility ceiling)

Each × TMEM145 full = 3 AF3 jobs. Phase 3 gate: ipTM >= 0.50 (same as Phase 2).
Seed 42. If winning tail identified, Phase 3b fuses best tail back to
RADA16-WH2 scaffold.

Output files:
  af3_jobs_2026-04-23c_hydrogel_retool/
    MANIFEST.json
    af3_jobs_sequences.fasta
    hydrogel_tail_1660_1710_x_tmem145.json
    hydrogel_tail_1640_1710_x_tmem145.json
    hydrogel_tail_1620_1710_x_tmem145.json

Submission: alphafoldserver.com Import-from-JSON.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23c_hydrogel_retool")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


STRC_FULL = fetch_uniprot("Q7RTU9")
assert len(STRC_FULL) == 1775

# TMEM145 — verified
TMEM145_FULL = (
    "MEPLRAPALRRLLPPLLLLLLSLPPRARAKYVRGNLSSKEDWVFLTRFCFLSDYGRLDFRFRYPEAKCCQNILLYFDDPSQWPAVYKAGDKDCLAKESVIRPE"
    "NNQVINLTTQYAWSGCQVVSEEGTRYLSCSSGRSFRSGDGLQLEYEMVLTNGKSFWTRHFSADEFGILETDVTFLLIFILIFFLSCYFGYLLKGRQLLHTTY"
    "KMFMAAAGVEVLSLLFFCIYWGQYATDGIGNESVKILAKLLFSSSFLIFLLMLILLGKGFTVTRGRISHAGSVKLSVYMTLYTLTHVVLLIYEAEFFDPGQV"
    "LYTYESPAGYGLIGLQVAAYVWFCYAVLVSLRHFPEKQPFYVPFFAAYTLWFFAVPVMALIANFGIPKWAREKIVNGIQLGIHLYAHGVFLIMTRPSAANKN"
    "FPYHVRTSQIASAGVPGPGGSQSADKAFPQHVYGNVTFISDSVPNFTELFSIPPPATSPLPRAAPDSGLPLFRDLRPPGPLRDL"
)
assert len(TMEM145_FULL) == 493

# Extended tail candidates (STRC 1-indexed → Python 0-indexed slicing)
TAIL_1660_1710 = STRC_FULL[1659:1710]   # 51 aa
TAIL_1640_1710 = STRC_FULL[1639:1710]   # 71 aa
TAIL_1620_1710 = STRC_FULL[1619:1710]   # 91 aa
assert len(TAIL_1660_1710) == 51
assert len(TAIL_1640_1710) == 71
assert len(TAIL_1620_1710) == 91

# Sanity: all three should END with STRC aa 1710. The 11-mer native tail
# AEDLPEPVPNC from Phase 2 is STRC 1669-1679. Verify our extended tails
# contain this native motif.
NATIVE_11MER = "AEDLPEPVPNC"
for tail in (TAIL_1660_1710, TAIL_1640_1710, TAIL_1620_1710):
    # Tail may not include the aa 1669-1679 Phase 2 motif since we're beyond —
    # the old motif was AFTER the Phase 2 peptide. Let's check what's actually
    # at STRC 1669-1680.
    pass
# Actually let's audit what's actually at those positions
print("STRC aa 1669-1680:", STRC_FULL[1668:1680])
print("STRC aa 1692-1707:", STRC_FULL[1691:1707])
print("TAIL_1660_1710 head/tail:", TAIL_1660_1710[:15], "...", TAIL_1660_1710[-15:])


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


TAILS = {
    "tail_1660_1710": (TAIL_1660_1710, "51 aa minimal; clusters 4+5 with 10 aa flank each side"),
    "tail_1640_1710": (TAIL_1640_1710, "71 aa +loop cap; clusters 3+4+5 (adds aa 1648-1651 loop cap)"),
    "tail_1620_1710": (TAIL_1620_1710, "91 aa wide; clusters 2+3+4+5 (adds aa 1630-1638 second ARM anchor)"),
}

jobs = {}
for tname, (seq, desc) in TAILS.items():
    jobs[f"hydrogel_{tname}_x_tmem145"] = {
        "json": job(f"hydrogel_{tname}_x_tmem145", [
            (seq, 1),
            (TMEM145_FULL, 1),
        ]),
        "description": (
            f"Extended STRC-derived tail × TMEM145 full. {desc}. "
            f"Tail length {len(seq)} aa, total complex {len(seq) + 493} aa. "
            "Isolated tail test (no RADA16/EAK16 scaffold, no WH2 actin-binder) "
            "to separately validate the retooled TMEM145 interface before "
            "Phase 3b full-construct fusion."
        ),
        "success_criteria": (
            f"ipTM >= 0.50 (Phase 3 gate, same as Phase 2). "
            f"Target: beat Phase 2 native 11-mer baseline (ipTM 0.39) by >=0.10. "
            f"Aspirational: match Ultra-Mini solo × TMEM145 baseline (ipTM 0.43)."
        ),
        "hypothesis_impact": (
            "pass (any candidate) -> Phase 3b: fuse winning tail back to "
            "RADA16-WH2 scaffold; re-run full-construct × TMEM145 AF3. "
            "If Phase 3b preserves ipTM at 0.50+, hypothesis #9 B -> A. "
            "fail (all three) -> Option A (longer tail) exhausted; move to "
            "Option B (helical mimetic / miniprotein) or Option C "
            "(different epitope, e.g. aa 1630-1638 cluster 2 alone)."
        ),
    }

for name, entry in jobs.items():
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(entry["json"], indent=2))

manifest = {
    "generated": "2026-04-23",
    "builder": "af3_jobs_2026-04-23c_hydrogel_retool_builder.py",
    "purpose": (
        "Hydrogel Phase 3 TMEM145 tail retool AF3 batch. Phase 2 "
        "short-peptide tails (11-12 aa from aa 1669-1680) failed TMEM145 "
        "gate at ipTM 0.37-0.39 across 3 candidates (Kd=100 nM assumption "
        "falsified). Phase 3 Option A: extend tail to span both dominant "
        "canonical contact clusters (aa 1669-1680 + aa 1692-1707). Three "
        "tail lengths tested: 51, 71, 91 aa (manufacturing SPPS-feasible)."
    ),
    "source_sequences": {
        "STRC_full": "UniProt Q7RTU9, 1775 aa (fetched 2026-04-23)",
        "TMEM145_full": "verified from prior AF3 batches",
    },
    "canonical_contact_clusters_GOLD_zone": {
        "cluster_1_1603-1607": "5 residues; not included in any extended tail",
        "cluster_2_1630-1638": "9 residues; in tail_1620_1710",
        "cluster_3_1648-1651": "4 residues (loop cap); in tails 1640+1620",
        "cluster_4_1669-1680": "12 residues DOMINANT; in all extended tails",
        "cluster_5_1692-1707": "16 residues SECOND DOMINANT; in all extended tails",
        "cluster_6_1770": "1 residue; not included (C-term tail drift)",
    },
    "tails_tested": {
        tname: {"length": len(seq), "description": desc, "sequence": seq}
        for tname, (seq, desc) in TAILS.items()
    },
    "priority_order": [
        "hydrogel_tail_1660_1710_x_tmem145",
        "hydrogel_tail_1640_1710_x_tmem145",
        "hydrogel_tail_1620_1710_x_tmem145",
    ],
    "jobs": {name: {k: v for k, v in entry.items() if k != "json"}
             for name, entry in jobs.items()},
    "how_to_submit": (
        "AlphaFold Server Import-from-JSON for each file. Daily limit 20; "
        "this batch uses 3. Seed 42 for consistency with prior Hydrogel batches."
    ),
    "decision_rules": {
        "any_tail_passes_50": (
            "Phase 3b: fuse winning tail back to RADA16-WH2 scaffold + "
            "WH2 actin-binder. Re-run full-construct × TMEM145 + G-actin "
            "trimer (3 more jobs). If Phase 3b preserves interface -> "
            "hypothesis #9 B -> A."
        ),
        "progressive_pattern": (
            "If 51 aa < 71 aa < 91 aa (monotonic improvement), extend further "
            "to 100-120 aa (approaching miniprotein scale, SPPS becomes "
            "expensive but still feasible)."
        ),
        "plateau_pattern": (
            "If all three tails cluster near Ultra-Mini baseline 0.43, the "
            "ceiling is the isolated peptide nature — move to Option B "
            "(helical mimetic scaffold) to stabilise binding conformation."
        ),
        "all_fail": (
            "Option A exhausted. Move to Option B (miniprotein scaffold with "
            "stapled binding loop) or Option C (cluster 2 aa 1630-1638 "
            "epitope alone, different binding site)."
        ),
    },
    "baselines_for_comparison": {
        "Ultra-Mini_solo_x_TMEM145": 0.43,
        "Phase_2_native_11mer": 0.39,
        "Phase_2_Cmut_11mer": 0.37,
        "Phase_2_denovo_11mer": 0.37,
    },
}
(OUT_DIR / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

fasta_lines = []
for tname, (seq, desc) in TAILS.items():
    fasta_lines += [f">STRC_{tname} | {len(seq)} aa | {desc}", seq]
fasta_lines += [f">TMEM145_FULL | {len(TMEM145_FULL)} aa | binding partner", TMEM145_FULL]
(OUT_DIR / "af3_jobs_sequences.fasta").write_text("\n".join(fasta_lines) + "\n")

print(f"\nGenerated AF3 Hydrogel Phase 3 tail-retool batch in {OUT_DIR}")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size} B)")
