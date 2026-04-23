#!/usr/bin/env python3
"""
Build AF3 submission batch 2026-04-23f — Hydrogel Phase 4b.

Addresses 4 concerns raised by Phase 4a CIF forensics:
  A. Multi-seed robustness — Phase 3b used seed 42 only; ipTM 0.57/0.51
     could be a lucky-seed fluke. Submit tail91 × TMEM145 and × actin
     with 5 seeds each → 25 structures per job, σ-estimate on ipTM.
  B. Triple-complex test — Phase 4a revealed single-partner artifact
     (WH2 binding TMEM145 in x_tmem145 job because actin absent;
     WH2 not binding actin in x_actin job because TMEM145 absent,
     tail steals actin). Submit tail91 + TMEM145 + G-actin trimer
     (5-chain complex) forcing AF3 to partition correctly.
  C. WH2-ablation control — replace WH2 with a polyGS spacer of same
     length. If tail91_noWH2 × TMEM145 retains ipTM ~0.5, the WH2→TMEM145
     "contact mass" in Phase 3b was just a hydrophobic artifact and
     the tail alone carries the signal. If ipTM drops to <0.4, WH2
     contributes real binding energy.
  D. Tail-length robustness — include tail91 with single-seed (for 4b
     seed 42 baseline) AND a new tail_84 variant (aa 1627-1710, drops
     N-term 6 aa of cluster 2) to find the minimal tail length that
     clears scaffold invasion.

JOBS (prioritised under 20/day budget — this batch: 8 jobs):
  1. tail91_x_tmem145_5seeds         (robustness, pair)
  2. tail91_x_actin_5seeds           (robustness, pair)
  3. tail91_x_tmem145_x_actin_trimer (triple complex, seed 42)
  4. tail71_x_tmem145_x_actin_trimer (triple complex, seed 42)
  5. tail91_noWH2_x_tmem145          (ablation, seed 42)
  6. tail91_noWH2_x_actin            (ablation, seed 42)
  7. tail_84_x_tmem145               (tail-length sweep, seed 42)
  8. tail_84_x_actin                 (tail-length sweep, seed 42)
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23f_hydrogel_phase4b")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


STRC_FULL = fetch_uniprot("Q7RTU9")
ACTG1 = fetch_uniprot("P63261")
assert len(STRC_FULL) == 1775 and len(ACTG1) == 375

TMEM145_FULL = (
    "MEPLRAPALRRLLPPLLLLLLSLPPRARAKYVRGNLSSKEDWVFLTRFCFLSDYGRLDFRFRYPEAKCCQNILLYFDDPSQWPAVYKAGDKDCLAKESVIRPE"
    "NNQVINLTTQYAWSGCQVVSEEGTRYLSCSSGRSFRSGDGLQLEYEMVLTNGKSFWTRHFSADEFGILETDVTFLLIFILIFFLSCYFGYLLKGRQLLHTTY"
    "KMFMAAAGVEVLSLLFFCIYWGQYATDGIGNESVKILAKLLFSSSFLIFLLMLILLGKGFTVTRGRISHAGSVKLSVYMTLYTLTHVVLLIYEAEFFDPGQV"
    "LYTYESPAGYGLIGLQVAAYVWFCYAVLVSLRHFPEKQPFYVPFFAAYTLWFFAVPVMALIANFGIPKWAREKIVNGIQLGIHLYAHGVFLIMTRPSAANKN"
    "FPYHVRTSQIASAGVPGPGGSQSADKAFPQHVYGNVTFISDSVPNFTELFSIPPPATSPLPRAAPDSGLPLFRDLRPPGPLRDL"
)
assert len(TMEM145_FULL) == 493

WH2_NTERM = "RQLVKAIPDNCSKSNVSR"   # 18 aa
RADA16 = "ADARADARADARADA"          # 15 aa
LINK_5 = "GSGSG"

TAIL_71 = STRC_FULL[1639:1710]  # 71 aa, aa 1640-1710
TAIL_91 = STRC_FULL[1619:1710]  # 91 aa, aa 1620-1710
TAIL_84 = STRC_FULL[1626:1710]  # 84 aa, aa 1627-1710 (drops 6 aa N-term of cluster 2)
assert len(TAIL_71) == 71 and len(TAIL_91) == 91 and len(TAIL_84) == 84

PEPTIDE_TAIL91 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_91       # 134 aa
PEPTIDE_TAIL71 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_71       # 114 aa
PEPTIDE_TAIL84 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_84       # 127 aa

# WH2 ablation: replace 18 aa WH2 with 18 aa polyGS spacer (zero hydrophobic)
SPACER_18 = "GSGSGSGSGSGSGSGSGS"
assert len(SPACER_18) == 18
PEPTIDE_TAIL91_NOWH2 = SPACER_18 + LINK_5 + RADA16 + LINK_5 + TAIL_91  # 134 aa total (identical length)
assert len(PEPTIDE_TAIL91_NOWH2) == 134

FIVE_SEEDS = ["7", "13", "42", "99", "777"]


def job(name: str, chains: list[tuple[str, int]], seeds: list[str] = ("42",)) -> list[dict]:
    return [{
        "name": name,
        "modelSeeds": list(seeds),
        "sequences": [
            {"proteinChain": {"sequence": seq, "count": count}}
            for seq, count in chains
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }]


jobs = {}

# 1. Multi-seed robustness — tail91 × TMEM145
jobs["hydrogel_tail91_x_tmem145_5seeds"] = {
    "json": job("hydrogel_tail91_x_tmem145_5seeds",
                [(PEPTIDE_TAIL91, 1), (TMEM145_FULL, 1)],
                seeds=FIVE_SEEDS),
    "description": (
        "Multi-seed robustness for Phase 3b tail91 × TMEM145 result (ipTM 0.57, seed 42). "
        "5 seeds × 5 models = 25 structures. Measures σ(ipTM) to verify 0.57 is not a "
        "lucky-seed outcome. PASS = median ipTM >= 0.50 AND σ < 0.05."
    ),
    "success_criteria": "median ipTM >= 0.50 AND σ(ipTM) across 5 seeds < 0.05",
    "hypothesis_impact": (
        "Confirms/rejects Phase 3b positive result. Critical for Harvard-grade claim."
    ),
}

# 2. Multi-seed robustness — tail91 × actin
jobs["hydrogel_tail91_x_actin_5seeds"] = {
    "json": job("hydrogel_tail91_x_actin_5seeds",
                [(PEPTIDE_TAIL91, 1), (ACTG1, 3)],
                seeds=FIVE_SEEDS),
    "description": (
        "Multi-seed robustness for Phase 3b tail91 × actin trimer (ipTM 0.51, seed 42). "
        "5 seeds × 5 models. Same PASS criterion."
    ),
    "success_criteria": "median ipTM >= 0.50 AND σ(ipTM) < 0.05",
    "hypothesis_impact": "Same as above — actin interface is the weaker of the two (margin 0.01).",
}

# 3. Triple complex — tail91 + TMEM145 + actin trimer (FORCES correct partitioning)
jobs["hydrogel_tail91_x_tmem145_x_actin"] = {
    "json": job("hydrogel_tail91_x_tmem145_x_actin",
                [(PEPTIDE_TAIL91, 1), (TMEM145_FULL, 1), (ACTG1, 3)],
                seeds=["42"]),
    "description": (
        "5-chain complex: peptide + TMEM145 + G-actin trimer. Addresses Phase 4a finding "
        "that single-partner jobs produce cross-binding artifact (WH2 binds TMEM145 when "
        "actin is absent). This job forces AF3 to partition: WH2 must pick either TMEM145 "
        "or actin, not both. Expected outcome if mechanism is correct: WH2 → actin, tail → "
        "TMEM145. ipTM read on each pair (chain_pair_iptm) should hold."
    ),
    "success_criteria": "peptide-TMEM145 chain_pair_iptm >= 0.40 AND peptide-actin >= 0.40",
    "hypothesis_impact": (
        "The definitive mechanism test. If it passes, Phase 3b binding is real and "
        "architecturally correct — ready for S-tier promotion. If it fails (one partner "
        "dominates, the other collapses), we have cross-talk — demotion warranted."
    ),
}

# 4. Triple complex — tail71 variant for comparison
jobs["hydrogel_tail71_x_tmem145_x_actin"] = {
    "json": job("hydrogel_tail71_x_tmem145_x_actin",
                [(PEPTIDE_TAIL71, 1), (TMEM145_FULL, 1), (ACTG1, 3)],
                seeds=["42"]),
    "description": (
        "Same triple complex with tail71 construct. Phase 4a showed tail71 × actin has "
        "CORRECT WH2 binding (46% share) but tail91 × actin has NO WH2 engagement. "
        "The triple complex may vindicate tail71 as the cleaner architecture despite "
        "lower single-partner ipTM."
    ),
    "success_criteria": "peptide-TMEM145 >= 0.40 AND peptide-actin >= 0.40",
    "hypothesis_impact": (
        "If tail71 triple complex passes while tail91 fails → flip the carry-forward "
        "design to tail71 (shorter, cheaper, cleaner)."
    ),
}

# 5. WH2 ablation control × TMEM145
jobs["hydrogel_tail91_noWH2_x_tmem145"] = {
    "json": job("hydrogel_tail91_noWH2_x_tmem145",
                [(PEPTIDE_TAIL91_NOWH2, 1), (TMEM145_FULL, 1)],
                seeds=["42"]),
    "description": (
        "Negative control: replace WH2 (18 aa) with polyGS spacer of identical length. "
        "If tail91 × TMEM145 binding (ipTM 0.57) was driven by real tail engagement, "
        "this should retain ipTM ~0.5. If driven by WH2 hydrophobic artifact, ipTM should "
        "drop to tail-solo level (0.68 of isolated 91 aa tail) — but with the scaffold "
        "still present. Deconvolutes WH2 vs tail contribution."
    ),
    "success_criteria": "tail-region binding signal isolated",
    "hypothesis_impact": (
        "Quantifies the 'WH2 contact mass' artifact seen in Phase 4a (78.8% of TMEM145 "
        "contacts were WH2). If binding survives ablation, tail is the real driver."
    ),
}

# 6. WH2 ablation control × actin
jobs["hydrogel_tail91_noWH2_x_actin"] = {
    "json": job("hydrogel_tail91_noWH2_x_actin",
                [(PEPTIDE_TAIL91_NOWH2, 1), (ACTG1, 3)],
                seeds=["42"]),
    "description": (
        "Ablation × actin: WH2 is canonical G-actin binder. Removing it should collapse "
        "actin binding from 0.51 to near-zero. Confirms that without WH2, no actin binding."
    ),
    "success_criteria": "ipTM < 0.35 (collapse)",
    "hypothesis_impact": (
        "If actin binding survives WH2 ablation (ipTM stays high), Phase 4a finding that "
        "'tail91 tail binds actin instead of WH2' is confirmed and we need to redesign — "
        "tail has an unintended actin-binding motif. If actin binding collapses, WH2 is "
        "the right actin handle and Phase 4a was a single-partner artifact."
    ),
}

# 7. Tail-length sweep — 84 aa intermediate
jobs["hydrogel_tail84_x_tmem145"] = {
    "json": job("hydrogel_tail84_x_tmem145",
                [(PEPTIDE_TAIL84, 1), (TMEM145_FULL, 1)],
                seeds=["42"]),
    "description": (
        "Intermediate tail length (84 aa, aa 1627-1710 — drops 7 aa N-term of cluster 2) "
        "to find the minimal tail that buffers the scaffold. Phase 3b showed tail71 fails "
        "(0.35) and tail91 passes (0.57). 84 aa binary-searches the threshold. "
        "If 84 aa passes → confirms cluster 2 core (aa 1627+) is what matters; "
        "if 84 aa fails → full cluster 2 needed."
    ),
    "success_criteria": "ipTM >= 0.50",
    "hypothesis_impact": (
        "Optimises the manufacturability-vs-binding tradeoff: 7 aa shorter = ~5% cheaper "
        "SPPS if it still works. Refines Phase 2c peptide design."
    ),
}

# 8. Tail-length sweep — 84 aa × actin
jobs["hydrogel_tail84_x_actin"] = {
    "json": job("hydrogel_tail84_x_actin",
                [(PEPTIDE_TAIL84, 1), (ACTG1, 3)],
                seeds=["42"]),
    "description": (
        "84 aa tail × actin partner — confirms actin binding preserved at intermediate length."
    ),
    "success_criteria": "ipTM >= 0.50",
    "hypothesis_impact": "Completes the tail-length sweep on the actin axis.",
}


for name, entry in jobs.items():
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(entry["json"], indent=2))


manifest = {
    "generated": "2026-04-23",
    "builder": "af3_jobs_2026-04-23f_hydrogel_phase4b_builder.py",
    "hypothesis": "#9 STRC Synthetic Peptide Hydrogel HTC — Phase 4b comprehensive validation",
    "motivation": (
        "Phase 4a CIF interface forensics revealed single-partner binding artifacts: "
        "WH2 binds TMEM145 (not actin) when actin is absent; tail binds actin (not TMEM145) "
        "when TMEM145 is absent. This batch runs 4 orthogonal controls to disentangle "
        "and de-risk the Phase 3b positive result before committing to wet-lab."
    ),
    "n_jobs": len(jobs),
    "ipTM_gate": 0.50,
    "robustness_criterion_sigma": 0.05,
    "priority_order": list(jobs.keys()),
    "jobs": {name: {k: v for k, v in entry.items() if k != "json"}
             for name, entry in jobs.items()},
    "how_to_submit": (
        "AlphaFold Server Import-from-JSON for each file. 8 jobs this batch. "
        "Seed-list jobs (5seeds) count as 1 AF3 submission but produce 5x models. "
        "Triple-complex jobs exceed 1800 aa (134 + 493 + 3×375 = 1752 — within 5000 limit)."
    ),
    "expected_learnings": [
        "Statistical significance of 0.57/0.51 ipTM (robustness)",
        "Correct partitioning in triple complex (mechanism)",
        "WH2 vs tail contribution (ablation)",
        "Minimal tail length for scaffold insulation (optimization)",
    ],
    "go_no_go_matrix": {
        "all_pass": "Hypothesis #9 clears Phase 4 — advance to Martini3 CG MD (Phase 5) + Phase 2c wet-lab. Promotion candidate from A-tier to tentative-S-tier (pending wet-lab).",
        "robustness_pass_triple_fail": "Single-partner ipTM is reproducible but mechanism wrong (partners cross-bind). Redesign needed — likely shorter WH2 variant or C-terminal scaffold placement.",
        "robustness_fail": "Phase 3b was a seed artifact. Demote back to B-tier pending deeper search.",
        "ablation_shows_wh2_drives_tmem145": "Tail binding energy was overestimated in Phase 3 solo test. Real tail-only ipTM may be <0.4. Demote to B.",
        "tail84_pass": "Carry-forward manufacturability-optimum shifts to tail84 (134→127 aa, ~5% cheaper).",
    },
}

(OUT_DIR / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

fasta = []
for name, seq in [
    ("PEPTIDE_TAIL71", PEPTIDE_TAIL71),
    ("PEPTIDE_TAIL84", PEPTIDE_TAIL84),
    ("PEPTIDE_TAIL91", PEPTIDE_TAIL91),
    ("PEPTIDE_TAIL91_NOWH2", PEPTIDE_TAIL91_NOWH2),
    ("TMEM145_FULL", TMEM145_FULL),
    ("ACTG1", ACTG1),
]:
    fasta += [f">{name} | {len(seq)} aa", seq]
(OUT_DIR / "af3_jobs_sequences.fasta").write_text("\n".join(fasta) + "\n")

print(f"Generated AF3 Hydrogel Phase 4b batch in {OUT_DIR}")
print(f"  {len(jobs)} jobs, covering: multi-seed × 2, triple-complex × 2, ablation × 2, tail-length × 2")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size} B)")
