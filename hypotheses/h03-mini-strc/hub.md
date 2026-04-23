---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h03]
hypothesis_num: 3
hypothesis_title: "STRC Mini-STRC Single-Vector Hypothesis"
tier: "S"
mech: 5
deliv: 5
misha_fit: 4
next_step: "Order gBlock, clone pAAV, coIP"
lit_audit: deferred
lit_audit_date: 2026-04-23
---

# h03 — Mini-STRC Single-Vector AAV

> 🔒 **LITERATURE AUDIT DEFERRED (2026-04-23, user directive).** Do NOT run a parameter-provenance audit on this hypothesis without explicit re-authorization. Reason: fear of destabilizing the active S-tier path while Holt lab is independently working on the same architecture. If parameters need retrieval for a specific proof, ask Egor first.

**mech.** Ultra-Mini STRC (aa 1075-1775, 2.1 kb CDS) preserves TMEM145-binding pocket and fits single AAV with OHC-specific promoter (B8+WPRE3).
**delivery.** Anc80L65 AAV; 60-100% cochlear transduction in OTOF trials precedent; single-vector.
**patient-fit.** Null-compatible (works on paternal 98 kb Δ); Shanghai Shu Yilai c.4976 knock-in mouse active.

## status

S-tier. next: order Ultra-Mini gBlock; clone pAAV B8-IgK-Ultra-Mini-WPRE3-bGH; Phase 4 HEK coIP.

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

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
