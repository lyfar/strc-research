---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h26]
status: A-tier
---

# h26 — Engineered Homodimer Avidity

**mech.** Single-point mutation at Ultra-Mini homodimer interface (ARM 1579-1581, R-R charge repulsion) converts weak dimer (ipTM 0.28-0.30) to stable dimer → avidity drops effective Kd 100-1000×.
**delivery.** Same AAV construct as h03 Mini-STRC (Anc80L65 + B8-IgK-Ultra-Mini-WPRE3-bGH); zero payload cost.
**patient-fit.** Direct upgrade to S-tier h03 clinical candidate; Misha-compatible by inheritance.

## status

A-tier. next: submit Phase 1 AF3 (5 jobs, homodimer × 2 TMEM145). If any mutant homodimer ≥ 0.50 AND TMEM145 ≥ 0.40 → Phase 2 SEC-AUC wet-lab.

## evidence

- [[STRC Engineered Homodimer Avidity]] — main hypothesis; Phase 1 AF3 batch built
- [[STRC Homodimer Interface From CIF]] — source: weak homodimer at ARM 1579-1581 identified

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § Hypothesis 1)
- `af3_jobs_2026-04-23d_engineered_homodimer_builder.py` — 5 jobs: WT ref + R1581Y + R1581F + S1579W + S1579F

## log

[[h26 log]]

## cross

- [[STRC Mini-STRC Single-Vector Hypothesis]] — h26 is a direct enhancement of h03 clinical candidate
- [[STRC Homodimer Interface From CIF]] — source analysis that identified the interface

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[about]` [[Misha]]
