---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h03]
status: S-tier
---

# h03 — Mini-STRC Single-Vector AAV

> 🔒 **LITERATURE AUDIT DEFERRED (2026-04-23, user directive).** Do NOT run a parameter-provenance audit on this hypothesis without explicit re-authorization. Reason: fear of destabilizing the active S-tier path while Holt lab is independently working on the same architecture. If parameters need retrieval for a specific proof, ask Egor first.

**mech.** Ultra-Mini STRC (aa 1075-1775, 2.1 kb CDS) preserves TMEM145-binding pocket and fits single AAV with OHC-specific promoter (B8+WPRE3).
**delivery.** Anc80L65 AAV; 60-100% cochlear transduction in OTOF trials precedent; single-vector.
**patient-fit.** Null-compatible (works on paternal 98 kb Δ); Shanghai Shu Yilai c.4976 knock-in mouse active.

## status

S-tier. next: order Ultra-Mini gBlock; clone pAAV B8-IgK-Ultra-Mini-WPRE3-bGH; Phase 4 HEK coIP.

<!-- RECENT:START auto-synced from log.md by sync_strc_indices.py — do not hand-edit -->

## recent

- Misha compound-het stack integration: **S held, scores unchanged**, endpoint framing clarified. h03 for Misha's specific genotype is physically ceiling-capped at MILD (~30-40 dB ABR) — at realistic transduction ε ≤ 0.5, ≥50% of OHCs remain non-transduced and cannot be rescued by AAV (maternal E1659A sub-threshold + paternal null → f_OHC < θ). "Meaningful rescue" endpoint intact (20-30 dB ABR improvement vs 64 dB baseline is clinically huge and matches OTOF trial responder criteria). "Cure" endpoint (NORMAL ≤ 25 dB) requires h01 co-therapy. No compute work needed on h03; this is endpoint semantics. → [[Misha Compound-Het Therapy Stack Model]]
- S-tier reinforced: AAV-LNP stack model formalises AAV as primary path → [[STRC AAV-LNP Stack PKPD]]
- S-tier reinforced: Iranfar 2026 CTM + Holt 2021 Science Advances + Regeneron AAV.104 → [[STRC Gene Therapy Landscape 2026]]
- Delivery 4→5 confirmed: all 3 AF3 gates closed (GOLD 0.68, homodimer, TMEM145 0.43) → [[STRC Mini-STRC Truncation Interface Validation]]
- Mech 4→5: Ultra-Mini sub-Å RMSD at TMEM145; unlocks 2.6 kb AAV headroom → [[STRC Mini-STRC Truncation Interface Validation]]

<!-- RECENT:END -->

## evidence

- [[STRC Mini-STRC Single-Vector Hypothesis]] — main hypothesis; full-stack validated 2026-04-21
- [[STRC Mini-STRC Truncation Interface Validation]] — sub-Å RMSD at TMEM145 pocket; Ultra-Mini = clinical 700-1775
- [[STRC Ultra-Mini CpG Depletion]] — 0 CpG CDS at 3.65% CAI cost; clinical CDS ready
- [[STRC Homodimer Interface From CIF]] — weak homodimer at ARM 1579-1581; real weak dimerization surface
- [[STRC Ultra-Mini Promoter Shortlist]] — B8+WPRE3-compact winner; OHC-exclusive
- [[STRC Ultra-Mini Full-Length TMEM145 AF3]] — ipTM 0.43, 23/41 contacts in GOLD zone
- [[STRC AF3 Ultra-Mini Multimer Validation]] — GOLD ipTM 0.68; homodimer weak-positive
- [[STRC Gene Therapy Landscape 2026]] — Iranfar 2026, Holt 2021, Regeneron AAV.104
- [[STRC Dual-Vector vs Single-Vector Transduction]] — engineering choice rationale

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § Hypothesis 3)
- `mini_strc_interface_preservation.py`, `cpg_depletion_mini_strc.py`, `cpg_depletion_ultra_mini_strc.py`
- `strc_homodimer_interface_from_cif.py`, `ultramini_homodimer_consensus.py`
- `ultra_mini_promoter_shortlist.py`, `ultramini_vector_cpg_audit.py`

## log

[[h03 log]]

## cross

- [[STRC Engineered Homodimer Avidity]] — h26 upgrades this candidate's Kd 100-1000× at zero payload cost
- [[STRC Synthetic Peptide Hydrogel HTC]] — h09 is independent path; if both pass → stack

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[about]` [[Misha]]
