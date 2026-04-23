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
status: A-tier
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

<!-- RECENT:START auto-synced from log.md by sync_strc_indices.py — do not hand-edit -->

## recent

- Phase 4i Kd × Kd sensitivity sweep (5×5, KD_TMEM145 10 nM-10 μM × KD_WH2_FACTIN 1 μM-1 mM) on therapeutic window: **f≥0.3 reachable on 80% of plausible grid** (PASS 56% single-dose + MARGINAL 24%); FAIL 20% concentrated at WH2 Kd = 1 mM column (Tβ4 analog floor, toxic ceiling beats f=0.3). Strategic finding: **WH2 Kd dominates, TMEM145 Kd is secondary**. Reorders wet-lab priority → Phase 2c WH2 bundling gate #1; TMEM145 SPR/BLI → #2 (ipSAE already operable-band-placed). Verdict `A_hold_multi_dose_ok`, not S-robust (70% PASS threshold). A held, Mech 4 held, Deliv 3 held. → [[STRC h09 Phase 4i Kd Sensitivity Sweep 2026-04-23]]
- Dunbrack 2025 ipSAE reassessment on Ultra-Mini × TMEM145 CIFs: **GOLD pruned ipSAE 0.591** sits in known-binder-zone (NFAT-CnB 0.55 … NFAT-CnA 0.78 physiological binders) → gate 3 (STRC×TMEM145 Kd) narrowed: placeholder 100 nM now lit-band-consistent (10 nM-10 µM Calcineurin-family band) instead of "arbitrary"; absolute Kd still wet-lab-gated. Full-length ipSAE 0.014 reaffirms known TM-collapse AF3 limit. Path B (AF3-only absolute Kd calibration) closed as lit-infeasible. A held, Mech 4 held. → [[STRC ipSAE Cross-Complex Reassessment 2026-04-23]]
- Phase 4e_v2 blend scaffold model: closes **1/3 S-tier gate** (118 aa tail > 12 aa empirical RADA16 limit). Recommended design point φ_tail=0.05 (5% RADA16-tail91, 95% plain), 2% w/v gel: PASS flag (within Gelain hybrid lit range), 6.4×10⁷ tails per OHC bundle shell, 2000× valency excess vs putative 3×10⁴ TMEM145 sites/bundle low est, 3.8× dose cost vs hypothetical 100% tail91. Gates 2 (WH2×F-actin Kd) + 3 (STRC×TMEM145 Kd) still open (wet-lab or Path B AF3 calibration). A held, Mech 4 held, Deliv 3 held. → [[STRC h09 Phase 4e_v2 Blend Scaffold]]
- Hub-note status banner added: both load-bearing unmeasured Kds (STRC×TMEM145 100 nM placeholder; WH2×F-actin 5 μM optimistic vs Tβ4 5-10 mM analog) + tail-91 > 12 aa RADA16 modification limit logged as promotion-to-S gates. Unsourced "48 h actin t½" in note replaced with Zhang 2012 Nature primary data (shaft stable months; tip β-actin t½ hours). A held. → [[STRC Synthetic Peptide Hydrogel HTC]]
- Parameter provenance audit: 3 critical gaps (118 aa tail, WH2×F-actin Kd, STRC×TMEM145 Kd); 7 phantom/wrong values; Mech 5→4 Deliv 4→3 (A held); triple promotion gate to S → [[STRC h09 Parameter Provenance Audit 2026-04-23]]

<!-- RECENT:END -->

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
