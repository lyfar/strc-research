---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h02]
hypothesis_num: 2
hypothesis_title: "STRC Piezoelectric TM Bioelectronic Amplifier"
tier: "B"
mech: 2
deliv: 1
misha_fit: 5
next_step: "Find OHC ligand + FEM strain-share"
status: S-tier
---

# h02 — Piezoelectric TM Bioelectronic Amplifier

**mech.** PVDF-TrFE piezoelectric polymer deposited on tectorial membrane amplifies sound-driven displacement, bypassing STRC protein entirely.
**delivery.** NP delivery to TM surface; A666 prestin-peptide targeting precedent.
**patient-fit.** Bypasses STRC null entirely; highest Misha-fit (works on paternal 98 kb Δ equally).

## status

S-tier. next: ex-vivo PVDF-TrFE deposition assay feasibility.

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## evidence

- [[STRC Piezoelectric TM Bioelectronic Amplifier]] — main hypothesis; Phases 1-3 complete
- [[STRC Piezo Delivery Feasibility OHC Targeting]] — Phase 3 NP delivery model; 92% audiogram coverage

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § Hypothesis 2)
- `piezo_voltage_budget.py`, `piezo_phase2_frequency_bundle.py`, `piezo_phase3_delivery_feasibility.py`

## log

[[h02 log]]

## cross

- [[STRC Stereocilia Bundle Mechanics Model]] — mechanical model underpinning voltage budget

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[about]` [[Misha]]
