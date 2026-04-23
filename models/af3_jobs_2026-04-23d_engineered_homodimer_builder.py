#!/usr/bin/env python3
"""
Build AlphaFold Server submission batch 2026-04-23d — Engineered Mini-STRC
Homodimer Avidity (hypothesis #26) Phase 1 AF3 test.

Four single-point mutations tested at the AF3-identified weak homodimer
interface (ARM 1579-1581) of Ultra-Mini STRC (aa 1075-1775, 701 aa). Each
mutant is tested as homodimer × 2 TMEM145 (full 493 aa) — tests both
dimer stabilisation AND each subunit's preserved TMEM145 binding.

Mutation candidates (from [[STRC Engineered Homodimer Avidity]]):
  1. S1579W — tryptophan pi-stack across C2 axis (leucine-zipper precedent)
  2. L1581Y — tyrosine pi-pi interaction with partner L1581
  3. S1579F — phenylalanine symmetric stack (shorter range than W, less steric risk)
  4. L1581F — analogous F-F stack at partner position

Gate per candidate (both required):
  - homodimer ipTM >= 0.50 (up from wild-type 0.28-0.30)
  - each Ultra-Mini × TMEM145 ipTM >= 0.40 (preserves monomer baseline)

5 AF3 jobs total (4 mutants + 1 WT reference complex for baseline).
Seed 42 for consistency with prior STRC batches.

Output files:
  af3_jobs_2026-04-23d_engineered_homodimer/
    MANIFEST.json
    af3_jobs_sequences.fasta
    engineered_homodimer_WT_ref.json
    engineered_homodimer_S1579W.json
    engineered_homodimer_L1581Y.json
    engineered_homodimer_S1579F.json
    engineered_homodimer_L1581F.json
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23d_engineered_homodimer")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])


STRC_FULL = fetch_uniprot("Q7RTU9")
assert len(STRC_FULL) == 1775

# Ultra-Mini = aa 1075-1775 (701 aa). Python 0-indexed slice [1074:1775]
ULTRAMINI = STRC_FULL[1074:1775]
assert len(ULTRAMINI) == 701

# Verify positions
# Ultra-Mini position = STRC_position - 1074
# STRC aa 1579 = Ultra-Mini position 505 (0-indexed 504)
# STRC aa 1581 = Ultra-Mini position 507 (0-indexed 506)
STRC_1579_0idx = 1578  # STRC 1-indexed aa 1579 → 0-idx 1578
STRC_1581_0idx = 1580
ULTRAMINI_1579_0idx = STRC_1579_0idx - 1074  # = 504
ULTRAMINI_1581_0idx = STRC_1581_0idx - 1074  # = 506

print(f"STRC aa 1579 (from full): {STRC_FULL[STRC_1579_0idx]} (expect S or similar)")
print(f"STRC aa 1581 (from full): {STRC_FULL[STRC_1581_0idx]} (expect L or similar)")
print(f"Ultra-Mini position 505 (= STRC 1579): {ULTRAMINI[ULTRAMINI_1579_0idx]}")
print(f"Ultra-Mini position 507 (= STRC 1581): {ULTRAMINI[ULTRAMINI_1581_0idx]}")
print(f"Context around 1579-1581 in Ultra-Mini: ...{ULTRAMINI[ULTRAMINI_1579_0idx-5:ULTRAMINI_1581_0idx+6]}...")


def mutate(seq: str, pos_0idx: int, new_aa: str) -> str:
    return seq[:pos_0idx] + new_aa + seq[pos_0idx + 1:]


# Actual residues at 1579 and 1581 (verify below, may not be S/L as assumed)
AA_1579 = STRC_FULL[STRC_1579_0idx]
AA_1581 = STRC_FULL[STRC_1581_0idx]

# Build mutant Ultra-Mini sequences — named by SUBSTITUTION applied
MUTANTS = {
    "WT_ref": ULTRAMINI,
    f"{AA_1579}1579W": mutate(ULTRAMINI, ULTRAMINI_1579_0idx, "W"),
    f"{AA_1581}1581Y": mutate(ULTRAMINI, ULTRAMINI_1581_0idx, "Y"),
    f"{AA_1579}1579F": mutate(ULTRAMINI, ULTRAMINI_1579_0idx, "F"),
    f"{AA_1581}1581F": mutate(ULTRAMINI, ULTRAMINI_1581_0idx, "F"),
}

for name, seq in MUTANTS.items():
    assert len(seq) == 701, f"{name} length mismatch: {len(seq)}"

# TMEM145 — verified
TMEM145_FULL = (
    "MEPLRAPALRRLLPPLLLLLLSLPPRARAKYVRGNLSSKEDWVFLTRFCFLSDYGRLDFRFRYPEAKCCQNILLYFDDPSQWPAVYKAGDKDCLAKESVIRPE"
    "NNQVINLTTQYAWSGCQVVSEEGTRYLSCSSGRSFRSGDGLQLEYEMVLTNGKSFWTRHFSADEFGILETDVTFLLIFILIFFLSCYFGYLLKGRQLLHTTY"
    "KMFMAAAGVEVLSLLFFCIYWGQYATDGIGNESVKILAKLLFSSSFLIFLLMLILLGKGFTVTRGRISHAGSVKLSVYMTLYTLTHVVLLIYEAEFFDPGQV"
    "LYTYESPAGYGLIGLQVAAYVWFCYAVLVSLRHFPEKQPFYVPFFAAYTLWFFAVPVMALIANFGIPKWAREKIVNGIQLGIHLYAHGVFLIMTRPSAANKN"
    "FPYHVRTSQIASAGVPGPGGSQSADKAFPQHVYGNVTFISDSVPNFTELFSIPPPATSPLPRAAPDSGLPLFRDLRPPGPLRDL"
)
assert len(TMEM145_FULL) == 493


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


jobs = {}
for mut_name, mut_seq in MUTANTS.items():
    job_name = f"engineered_homodimer_{mut_name}"
    # Homodimer × 2 TMEM145: two copies of Ultra-Mini (mutant), two copies of TMEM145
    jobs[job_name] = {
        "json": job(job_name, [
            (mut_seq, 2),       # Ultra-Mini homodimer
            (TMEM145_FULL, 2),  # each subunit × its own TMEM145
        ]),
        "description": (
            f"Ultra-Mini {mut_name} homodimer × 2 TMEM145. "
            f"{len(mut_seq)} × 2 + {len(TMEM145_FULL)} × 2 = "
            f"{2*len(mut_seq) + 2*len(TMEM145_FULL)} aa total complex. "
            f"Tests whether single-point mutation at homodimer interface "
            f"(ARM 1579-1581) stabilises the weak wild-type self-association "
            f"(AF3 ipTM 0.28-0.30 baseline) while preserving each subunit's "
            f"TMEM145 binding (ipTM 0.43 monomer baseline)."
        ),
        "success_criteria": (
            "Both criteria required: "
            "(1) homodimer self-ipTM (chain-pair A↔B) >= 0.50 (up from 0.28-0.30 WT baseline); "
            "(2) each Ultra-Mini × TMEM145 cross-ipTM >= 0.40 (preserves monomer baseline of 0.43). "
            "If chain-pair shows A↔C and B↔D cross-ipTM >= 0.40 and A↔B homodimer >= 0.50, PASS."
        ),
        "hypothesis_impact": (
            "pass -> #26 Engineered Homodimer Avidity hypothesis elevates to A-tier (min=4). "
            "Phase 2: wet-lab SEC-AUC on winning mutant to confirm solution dimer at physiological conc. "
            "fail (homodimer still < 0.40) -> single-point mutation insufficient to stabilise; move to "
            "double mutant (S1579W + L1581Y) or alternative engineered dimer scaffolds. "
            "fail (TMEM145 cross lost) -> mutation allosterically disrupts binding surface; kill specific mutant, try next."
        ),
    }

for name, entry in jobs.items():
    path = OUT_DIR / f"{name}.json"
    path.write_text(json.dumps(entry["json"], indent=2))

manifest = {
    "generated": "2026-04-23",
    "builder": "af3_jobs_2026-04-23d_engineered_homodimer_builder.py",
    "hypothesis": "#26 STRC Engineered Homodimer Avidity — [[STRC Engineered Homodimer Avidity]]",
    "purpose": (
        "Phase 1 AF3 test of 4 single-point mutations at the Ultra-Mini "
        "weak homodimer interface (ARM 1579-1581). Each mutant tested as "
        "homodimer × 2 TMEM145 — simultaneously gates (a) dimer stabilisation "
        "and (b) preserved TMEM145 binding per subunit. WT reference included "
        "for direct delta comparison."
    ),
    "baselines_for_comparison": {
        "Ultra-Mini_WT_homodimer": "ipTM 0.28-0.30 (5 models consensus, AF3 weak-positive)",
        "Ultra-Mini_WT_x_TMEM145_monomer": "ipTM 0.43 (full-length TMEM145 multimer context)",
        "Ultra-Mini_x_TMEM145_GOLD_pruned": "ipTM 0.68 (Derstroff-style pruning)",
    },
    "mutation_candidates": {
        "AA_at_1579_in_WT": AA_1579,
        "AA_at_1581_in_WT": AA_1581,
        "substitutions_tested": [
            f"{AA_1579}1579W (tryptophan π-stack across C2)",
            f"{AA_1581}1581Y (tyrosine π-π with partner chain L1581)",
            f"{AA_1579}1579F (phenylalanine, shorter-range than W)",
            f"{AA_1581}1581F (F-F stack, symmetric with position 1579)",
        ],
    },
    "jobs": {name: {k: v for k, v in entry.items() if k != "json"}
             for name, entry in jobs.items()},
    "priority_order": list(jobs.keys()),
    "how_to_submit": (
        "AlphaFold Server Import-from-JSON for each file. Daily limit 20; "
        "this batch uses 5. Seed 42 for reproducibility vs prior STRC homodimer runs."
    ),
    "decision_rules": {
        "any_mutant_passes_both_gates": (
            "Hypothesis #26 B → A → phase 2 wet-lab SEC-AUC + HEK293 co-IP "
            "for winner. Mutation carried forward as single additional line "
            "on top of Mini-STRC clinical candidate (no new AAV payload cost)."
        ),
        "all_mutants_fail_homodimer_gate": (
            "Single-point engineering insufficient. Two options: (a) double "
            "mutant S1579W + L1581Y / S1579F + L1581F, (b) pivot to "
            "de novo dimer scaffold (e.g. leucine zipper fusion at Ultra-Mini C-term)."
        ),
        "homodimer_passes_but_TMEM145_lost": (
            "Mutation at 1579-1581 allosterically disrupts GOLD-zone binding. "
            "Mutation site changes to a more peripheral dimer contact residue "
            "identified in [[STRC Homodimer Interface From CIF]] — require "
            "re-analysis of that data with broader contact tolerance."
        ),
    },
}
(OUT_DIR / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))

fasta_lines = []
for mut_name, mut_seq in MUTANTS.items():
    fasta_lines += [f">UltraMini_{mut_name} | 701 aa | hypothesis #26 candidate", mut_seq]
fasta_lines += [f">TMEM145_FULL | {len(TMEM145_FULL)} aa | binding partner", TMEM145_FULL]
(OUT_DIR / "af3_jobs_sequences.fasta").write_text("\n".join(fasta_lines) + "\n")

print(f"\nGenerated AF3 Engineered Homodimer Avidity batch in {OUT_DIR}")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size} B)")
