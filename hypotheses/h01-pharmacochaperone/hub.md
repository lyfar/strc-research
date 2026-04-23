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

A-tier. next: Phase 5 MD + ensemble MM-GBSA on 3 de-novo leads.

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
- `pharmacochaperone_phase5_md.py` — scaffold, deferred

## log

[[h01 log]]

## cross

- [[STRC Electrostatic Analysis E1659A]] — electrostatic input to pocket design
- [[STRC Mini-STRC Single-Vector Hypothesis]] — h03 is active fallback if h01 mechanism fails

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[about]` [[Misha]]
