---
date: 2026-04-23
type: hypothesis-hub
hypothesis_num: 8
hypothesis_title: "STRC ASO Exon Skipping"
tier: "C"
mech: 2
deliv: 3
misha_fit: 2
next_step: "Phase 3a morpholino + 3b gapmer"
tags: [strc, h08]
status: C-tier
---

# h08 — ASO Exon Skipping

Stub hub created 2026-04-23 during vault-stack migration. Detailed hypothesis content lives in [[STRC ASO Exon Skipping]]; this file exists so the Dataview ranking covers every numbered hypothesis and logs can accumulate here.

## status

C-tier. See [[STRC ASO Exon Skipping]] for mech/delivery/patient-fit reasoning and current state.

## active compute

_(none — add PIDs + ETAs here when jobs run)_

## next-step tree

Phase 3a morpholino + 3b gapmer

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC ASO Exon Skipping]]
