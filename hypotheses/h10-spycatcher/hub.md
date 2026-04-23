---
date: 2026-04-23
type: hypothesis-hub
hypothesis_num: 10
hypothesis_title: "STRC In Situ SpyCatcher Assembly"
tier: "C"
mech: 2
deliv: 2
misha_fit: 2
next_step: "Wet-lab only if #3 fails"
tags: [strc, h10]
status: C-tier
---

# h10 — SpyCatcher Assembly

Stub hub created 2026-04-23 during vault-stack migration. Detailed hypothesis content lives in [[STRC In Situ SpyCatcher Assembly]]; this file exists so the Dataview ranking covers every numbered hypothesis and logs can accumulate here.

## status

C-tier. See [[STRC In Situ SpyCatcher Assembly]] for mech/delivery/patient-fit reasoning and current state.

## active compute

_(none — add PIDs + ETAs here when jobs run)_

## next-step tree

Wet-lab only if #3 fails

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC In Situ SpyCatcher Assembly]]
