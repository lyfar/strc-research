---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h02]
status: S-tier
---

# h02 — Piezoelectric TM Bioelectronic Amplifier

**mech.** PVDF-TrFE piezoelectric polymer deposited on tectorial membrane amplifies sound-driven displacement, bypassing STRC protein entirely.
**delivery.** NP delivery to TM surface; A666 prestin-peptide targeting precedent.
**patient-fit.** Bypasses STRC null entirely; highest Misha-fit (works on paternal 98 kb Δ equally).

## status

S-tier. next: ex-vivo PVDF-TrFE deposition assay feasibility.

<!-- RECENT:START auto-synced from log.md by sync_strc_indices.py — do not hand-edit -->

## recent

- AAV-LNP stack model confirms Piezo S-tier independent of AAV path → [[STRC AAV-LNP Stack PKPD]]
- Parameter provenance audit complete: 3 phantom cites caught (K_D A666-prestin = 10 nM; ETA_POLY = 0.7; d31 Phase1/Phase2 inconsistency). 9 additional ⚠ constants. New topic files: [[piezoelectric-materials]], [[tectorial-membrane]]. Full audit at `/tmp/strc-lit-params-h02.md`. New finding: TM stiffness (24–210 kPa) vs. PVDF-TrFE (3 GPa) mismatch not currently modeled — potential hypothesis-level mechanical flaw.
- Post-audit fix applied: scripts pinned phantoms (A666, ETA_POLY, SELECT_S), d31 reconciled to −12 pC/N (was silent 2× split). Phase 1 baseline no longer passes 60 dB. TM mismatch 10⁵× flagged in output JSON. Tier S → B. Mech 3→2, Deliv 2→1. → [[STRC h02 Parameter Provenance Audit 2026-04-23]]
- Phases 1-3 complete; delivery model 92% audiogram coverage → [[STRC Piezo Delivery Feasibility OHC Targeting]]
- Initial assignment: S-tier, Mech 3, Deliv 2, Misha 5 → [[STRC Hypothesis Ranking]]

<!-- RECENT:END -->

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
