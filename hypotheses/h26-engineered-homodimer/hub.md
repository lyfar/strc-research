---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h26]
status: B-tier
---

# h26 — Engineered Homodimer Avidity

**mech.** Single-point mutation at Ultra-Mini homodimer interface (ARM 1579-1581) tested for dimer stabilization; Phase 1 AF3 FAILED — all 4 mutants destabilize dimer further, R-R repulsion hypothesis falsified.
**delivery.** Would be drop-in to h03 Mini-STRC AAV (Anc80L65 + B8-IgK-Ultra-Mini-WPRE3-bGH) if any mutant passed; currently none does.
**patient-fit.** Misha-compatible by inheritance from h03; not independently therapeutic.

## status

B-tier (was A, 2026-04-23). next: Phase 1c structural review — re-run contact analysis at 5 Å cutoff to identify secondary contact residues outside 1579-1581. If new cluster → Phase 1b double-mutant AF3. If no cluster → C.

## evidence

- [[STRC Engineered Homodimer Avidity]] — main hypothesis (R-R repulsion prediction, now falsified)
- [[STRC Engineered Homodimer Phase 1 Results]] — Phase 1 AF3 results: 4/4 mutants fail homodimer gate
- [[STRC Homodimer Interface From CIF]] — source: weak homodimer at ARM 1579-1581 identified

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § h26)
- `af3_jobs_2026-04-23d_engineered_homodimer_builder.py` — 5 jobs: WT ref + R1581Y + R1581F + S1579W + S1579F
- `af3_jobs_2026-04-23d_engineered_homodimer_forensics.py` — Phase 1 results parser (4 mutants × 1-2 seeds, gate evaluation)

## artifacts

research/strc/hypotheses/h26-engineered-homodimer/artifacts/phase1_af3_2026-04-23/ — 6 AF3 runs (4 mutants + 2 restart seeds) + analysis_summary.json

## log

[[h26 log]]

## cross

- [[STRC Mini-STRC Single-Vector Hypothesis]] — h26 is a direct enhancement of h03 clinical candidate
- [[STRC Homodimer Interface From CIF]] — source analysis that identified the interface

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[about]` [[Misha]]
