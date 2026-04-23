---
date: 2026-04-23
type: hypothesis-hub
hypothesis_num: 6
hypothesis_title: "STRC mRNA Therapy Hypothesis"
tier: "B"
mech: 2
deliv: 2
misha_fit: 2
next_step: "Retrieve OHC RBM24 titration lit"
tags: [strc, h06]
status: B-tier
---

# h06 — mRNA Therapy

Stub hub created 2026-04-23 during vault-stack migration. Detailed hypothesis content lives in [[STRC mRNA Therapy Hypothesis]]; this file exists so the Dataview ranking covers every numbered hypothesis and logs can accumulate here.

## status

B-tier. See [[STRC mRNA Therapy Hypothesis]] for mech/delivery/patient-fit reasoning and current state.

## active compute

_(none — add PIDs + ETAs here when jobs run)_

## next-step tree

Retrieve OHC RBM24 titration lit

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC mRNA Therapy Hypothesis]]
