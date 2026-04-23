---
date: 2026-04-23
type: hypothesis-hub
aliases:
- h26 hub
tags:
- strc
- h26
hypothesis_num: 26
hypothesis_title: STRC Engineered Homodimer Avidity
tier: B
mech: 2
deliv: 5
misha_fit: 4
next_step: Phase 1d AF3 A1078C/S1080C/S1579C
lit_audit: fixed
lit_audit_date: 2026-04-23
---

# h26 — Engineered Homodimer Avidity

**mech.** Single-point mutation at Ultra-Mini homodimer interface (ARM 1579-1581) tested for dimer stabilization; Phase 1 AF3 FAILED — all 4 mutants destabilize dimer further, R-R repulsion hypothesis falsified.
**delivery.** Would be drop-in to h03 Mini-STRC AAV (Anc80L65 + B8-IgK-Ultra-Mini-WPRE3-bGH) if any mutant passed; currently none does.
**patient-fit.** Misha-compatible by inheritance from h03; not independently therapeutic.

## status

B-tier. Phase 1 ARM R-R hypothesis falsified; Phase 1c GREEN on secondary stump-cluster 1077-1114 with native C1081 disulfide-engineering geometry. **Next: Phase 1d AF3 triple batch** (A1078C / S1080C / S1579C homotypic + A1078W negative control) — user-invoked external AF3.

## next-step tree

- Phase 1d AF3 PASS on any mutant (ipTM ≥0.50 + homodimer contacts ≥ WT) → **A-tier → drop-in to h03 Ultra-Mini as engineered-avidity variant** (zero payload cost)
- Phase 1d AF3 FAIL on all → C-tier (disulfide branch killed; no further compute)

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## evidence

- [[STRC Engineered Homodimer Avidity]] — main hypothesis (R-R repulsion prediction, now falsified)
- [[STRC Engineered Homodimer Phase 1 Results]] — Phase 1 AF3 results: 4/4 mutants fail homodimer gate
- [[STRC Homodimer Interface From CIF]] — source: weak homodimer at ARM 1579-1581 identified

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § h26)
- `af3_jobs_2026-04-23d_engineered_homodimer_builder.py` — 5 jobs: WT ref + R1581Y + R1581F + S1579W + S1579F
- `af3_jobs_2026-04-23d_engineered_homodimer_forensics.py` — Phase 1 results parser (4 mutants × 1-2 seeds, gate evaluation)
- `engineered_homodimer_phase1c_contact_cluster.py` — DBSCAN 6.5 Å spatial re-cluster on Ultra-Mini homodimer inter-chain contacts; per-cluster disulfide Cb-Cb scan 4.5-7.5 Å + mutation candidate ranking (steric / electrostatic / covalent)

## artifacts

research/strc/hypotheses/h26-engineered-homodimer/artifacts/phase1_af3_2026-04-23/ — 6 AF3 runs (4 mutants + 2 restart seeds) + analysis_summary.json

## log

[[h26 log]]

## cross

- [[STRC Mini-STRC Single-Vector Hypothesis]] — h26 is a direct enhancement of h03 clinical candidate
- [[STRC Homodimer Interface From CIF]] — source analysis that identified the interface

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[about]` [[Misha]]
