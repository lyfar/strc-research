---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h27]
hypothesis_num: 27
hypothesis_title: "STRC STRCP1 Activation Rescue"
tier: "C"
mech: 2
deliv: 2
misha_fit: 4
next_step: "Email Holt lab ribosome-profiling"
lit_audit: fixed
lit_audit_date: 2026-04-23
---

# h27 — STRCP1 Activation Rescue

**mech.** Activate endogenous STRCP1 paralog (100% identical to STRC at critical regions) via CRISPRa / dCas9-TET1 / steric-block ASO; if STRCP1 translates, it is a native STRC backup.
**delivery.** CRISPRa or steric-block ASO; delivery same class as active paths.
**patient-fit.** Works on any STRC-null including Misha; Gate 1 (translation evidence) is wet-lab-only block.

## status

C-tier. next: email Holt lab re STRCP1 ribosome-profiling data from STRC mouse work. If Gate 1 passes → A-tier.

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## evidence

- [[STRC STRCP1 Activation Rescue]] — main hypothesis
- [[STRC STRCP1 GTEx Expression Check]] — STRCP1 transcribed all 54 GTEx tissues (0.01-2.6 TPM)
- [[STRC PE Phase 3_5 STRCP1-Aware Redesign]] — STRCP1 core 100% identical at locus; reframed here

## scripts

None yet. Wet-lab gate 1 blocks compute.

## log

[[h27 log]]

## cross

- [[Prime Editing for STRC]] — h07 C-tier because of same STRCP1 identity problem
- [[STRC ASO Exon Skipping]] — h08 C-tier; h27 reframes their "problem" as an asset

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[about]` [[Misha]]
