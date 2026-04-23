---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h01]
status: A-tier
---

# h01 — Pharmacochaperone E1659A

**mech.** Small-molecule fold-stabilizer binds K1141 pocket on E1659A STRC, rescuing misfolded stereocilin.
**delivery.** Oral/topical small molecule; round-window permeability precedent exists.
**patient-fit.** Maternal allele c.4976A>C E1659A directly targeted; extracellular concern flags mechanism.

## status

A-tier held. **Phase 5 MD + ensemble rescoring DELIVERED 2026-04-23** (pipeline on local Mac, 62 ns/day, 2 ns MD + 20 snapshots + Vina re-docking of Phase 4b top-5 leads). Result: **RED-LIGHT — all leads f_PC < 0.10 at therapeutic [L]=10 μM** (best = diflunisal positive 0.083); Phase 4b single-structure was over-optimistic by 0.36-0.92 kcal/mol. See [[STRC h01 Phase 5 MD Ensemble Rescoring 2026-04-23]].

Hypothesis intact — only the current shortlist is insufficient. **Next: Phase 3c v2 expanded virtual screen** — DrugBank FDA library (~2500 approved) + ZINC22 carboxylate tranche (~50 k) + fragment-based screen at K1141 pocket on 5 ensemble-sampled receptor conformers; filter at ensemble ΔG ≤ −7.5 kcal/mol (target Kd ≤ 5 μM → f_PC ≥ 0.67 at [L]=10 μM).

## evidence

- [[STRC Pharmacochaperone Virtual Screen E1659A]] — main hypothesis note; Phases 0-3 complete
- [[STRC Pharmacochaperone K1141 Fragment Pocket]] — pocket characterization
- [[STRC Pharmacochaperone Phase 4 Plan]] — 7-gate validation ladder
- [[STRC Pharmacochaperone Phase 4a Pocket Reproducibility]] — PASS 4/5 CIFs
- [[STRC Pharmacochaperone Phase 4b Smoke Test]] — PASS on LE; leads beat diflunisal 29-57%
- [[STRC Pharmacochaperone Phase 4c WT Decoy]] — FAIL: leads prefer WT; aspect (a) flagged
- [[STRC Pharmacochaperone Phase 4e Off-Target Selectivity]] — soft-FAIL (margin 0.15 vs 0.20)
- [[STRC Pharmacochaperone Phase 4f Interface Rescue SMOKE]] — method inadequate; Phase 4f canceled
- [[STRC Pharmacochaperone Phase 4g Repurpose Screen]] — 4PBA mut-preferring (−0.32); branch later killed
- [[STRC Repurposed FDA Chaperone Branch]] — KILLED: STRC extracellular, no ER target for 4PBA
- [[STRC Pharmacochaperone Virtual Screen Ranked Leads]] — top-5 shortlist
- [[STRC Pharmacophore Model K1141 Pocket]] — Phase 3a pharmacophore

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § Hypothesis 1)
- `pharmacochaperone_target_prep.py` through `pharmacochaperone_phase4g_repurpose_screen.py` (16 scripts)
- `pharmacochaperone_phase5_md.py` — scaffold (GROMACS-based, superseded by 5a/5b OpenMM)
- `pharmacochaperone_phase5a_apo_md_smoke.py` — OpenMM + AMBER14SB apo MD, SMOKE + 2 ns production delivered 2026-04-23
- `pharmacochaperone_phase5b_ensemble_redock.py` — Vina ensemble re-docking on Phase 5a snapshots; f_PC estimates per lead

## log

[[h01 log]]

## cross

- [[STRC Electrostatic Analysis E1659A]] — electrostatic input to pocket design
- [[STRC Mini-STRC Single-Vector Hypothesis]] — h03 is active fallback if h01 mechanism fails

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[about]` [[Misha]]
