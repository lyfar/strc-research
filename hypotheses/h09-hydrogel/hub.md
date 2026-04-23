---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h09]
hypothesis_num: 9
hypothesis_title: "STRC Synthetic Peptide Hydrogel HTC"
tier: "A"
mech: 4
deliv: 3
misha_fit: 3
next_step: "Phase 2c WH2 bundling + 4b AF3"
lit_audit: partial
lit_audit_date: 2026-04-23
---

# h09 — Synthetic Peptide Hydrogel HTC

**mech.** Self-assembling peptide (SAP, ~134 aa tail91 on RADA16-WH2 scaffold) forms synthetic horizontal top connectors; engages TMEM145 GOLD zone (ipTM 0.68) + G-actin (ipTM 0.51).
**delivery.** Ototopical ear drop + LIFU round-window permeation; PKPD window 0.13-1.32 mg; $42/ear GMP via NCL.
**patient-fit.** Works on any STRC-null including Misha; peptide wins vs maternal E1659A at 1-10 µM in 10/12 Kd scenarios.

## status

A-tier. next priority reordered per [[STRC h09 Phase 4i Kd Sensitivity Sweep 2026-04-23]]:
1. **Phase 2c WH2 bundling** (HEK293 GFP-actin confocal) — gate #1, dominates therapeutic window
2. Phase 4b AF3 submission (8 jobs built, awaiting external AF3 server) — triple-complex + σ<0.05 needed
3. TMEM145 SPR/BLI — gate #2, already in operable ipSAE-narrowed band

S-tier promotion path: WH2 Kd wet-lab narrowing into PASS zone (< 50 μM) + 4b triple-complex σ<0.05 → S candidate.

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## evidence

- [[STRC Synthetic Peptide Hydrogel HTC]] — main hypothesis
- [[STRC Hydrogel HTC Phase 1 Self-Assembly]] — Phase 1-3b: analytical pre-screen → AF3 retool → full-construct validation
- [[STRC Hydrogel Phase 4 Computational Campaign]] — 8-axis campaign: 4a-4h all complete except 4b AF3 batch

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § Hypothesis 1, h09-specific rows)
- `hydrogel_phase1_self_assembly.py`
- `af3_jobs_2026-04-23b_hydrogel_builder.py` through `af3_jobs_2026-04-23f_hydrogel_phase4b_builder.py` (4 builders)
- `hydrogel_phase4a_cif_interface_forensics.py` through `hydrogel_phase4h_endogenous_strc_competition.py` (8 scripts)
- `hydrogel_phase4c_sequence_liabilities.py` — tail91_v2 Cys→Ser variant design

## log

[[h09 log]]

## cross

- [[STRC Mini-STRC Single-Vector Hypothesis]] — h03 is primary AAV path; h09 is independent materials path
- [[STRC Engineered Homodimer Avidity]] — h26 strengthens h03; h09 is parallel non-AAV track

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[about]` [[Misha]]
