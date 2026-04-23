---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h05]
status: A-tier
---

# h05 — Calcium Oscillation Acoustic Therapy

**mech.** Chronic acoustic stimulation at 45-60 dB LAeq drives AC1-CREB Ca2+ oscillations, upregulating endogenous STRC expression via maternal E1659A allele.
**delivery.** Hearing aid or Touch Grass app; non-invasive, zero delivery risk.
**patient-fit.** Maternal allele only (E1659A contributes 0.5 dose); affinity penalty likely 10-100x puts rescue outside feasibility window.

## status

A-tier. next: SPR/BLI measurement of E1659A TMEM145 affinity penalty (wet-lab gate).

<!-- RECENT:START auto-synced from log.md by sync_strc_indices.py — do not hand-edit -->

## recent

- Misha-specific PD: realistic affinity penalty 10-100× outside feasibility window; acoustic therapy = free adjunct → [[STRC Ca Oscillation Maternal-Only PD]]
- No tier change; stays A. Next: SPR/BLI E1659A affinity penalty measurement.
- Parameter provenance audit complete + post-audit fix applied: 4 phantom cites removed (Wu 2011, Sharma 2018, Cha 2010, Krey 2015), Chao 2010→Chao 2011 (PMID 21458670), PKA K_cAMP cite corrected to Surdo 2017 (PMID 29074866). CREB-P dephos 2-4× faster than lit flagged, STRC mRNA/protein t½ cross-script mismatches (4× and 20×) flagged explicitly. Phase 3 topological verdict defensible. Tier A held. → [[STRC h05 Parameter Provenance Audit 2026-04-23]]
- Phase 3 globally stable; max Re(λ) = −5×10⁻⁶/s; true Hill n=4.33 confirmed → [[STRC Calcium Oscillation Acoustic Therapy]]
- Phase 2 robustness index 0.91; Phase 1 RBM24 ODE baseline established → [[STRC Calcium Oscillation Acoustic Therapy]]

<!-- RECENT:END -->

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
