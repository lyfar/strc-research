---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h05]
hypothesis_num: 5
hypothesis_title: "STRC Calcium Oscillation Acoustic Therapy"
tier: "A"
mech: 3
deliv: 4
misha_fit: 2
next_step: "SPR/BLI E1659A affinity penalty"
lit_audit: fixed
lit_audit_date: 2026-04-23
---

# h05 — Calcium Oscillation Acoustic Therapy

**mech.** Chronic acoustic stimulation at 45-60 dB LAeq drives AC1-CREB Ca2+ oscillations, upregulating endogenous STRC expression via maternal E1659A allele.
**delivery.** Hearing aid or Touch Grass app; non-invasive, zero delivery risk.
**patient-fit.** Maternal allele only (E1659A contributes 0.5 dose); affinity penalty likely 10-100x puts rescue outside feasibility window.

## status

A-tier. next: SPR/BLI measurement of E1659A TMEM145 affinity penalty (wet-lab gate).

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## evidence

- [[STRC Calcium Oscillation Acoustic Therapy]] — main hypothesis; Phases 1-3 complete
- [[STRC AC1-CREB Alternative Hypothesis]] — AC1-CREB pivot from RBM24 path
- [[STRC Ca Oscillation Maternal-Only PD]] — Misha-specific PD: rescue only at penalty ≤3× and fold ≥1.3×

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § Hypothesis 5)
- `ca_oscillation_rbm24_ode.py`, `ca_oscillation_ac1_creb_pivot.py`
- `ca_oscillation_phase2_sensitivity.py`, `ca_oscillation_phase3_bifurcation.py`
- `ca_osc_maternal_allele_only.py`

## log

[[h05 log]]

## cross

- [[STRC mRNA Therapy Hypothesis]] — h06 Strategy A shares RBM24 regulatory context
- [[STRC Mini-STRC Single-Vector Hypothesis]] — h03 is primary; acoustic therapy is free adjunct regardless

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[about]` [[Misha]]
