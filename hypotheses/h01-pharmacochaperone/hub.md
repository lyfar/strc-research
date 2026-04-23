---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h01]
hypothesis_num: 1
hypothesis_title: "STRC Pharmacochaperone Virtual Screen E1659A"
tier: "A"
mech: 3
deliv: 4
misha_fit: 4
next_step: "Phase 3c v3b dock + 5d E1659A MD running"
lit_audit: fixed
lit_audit_date: 2026-04-23
active_runs:
  - phase: "3c-v3b"
    name: "fenamic+covalent 12k-ligand dock"
    pid: 73601
    eta: "2026-04-24 ~05:10 local"
  - phase: "5d"
    name: "E1659A full-length mutant MD"
    pid: 28875
    eta: "2026-04-23 ~23:00 local"
  - phase: "5e"
    name: "mutant-ensemble re-dock (scaffold ready)"
    eta: "post-5d + post-v3b-Stage-2"
---

# h01 — Pharmacochaperone E1659A

**mech.** Small-molecule fold-stabilizer binds K1141 pocket on E1659A STRC, rescuing misfolded stereocilin.
**delivery.** Oral/topical small molecule; round-window permeability precedent exists.
**patient-fit.** Maternal allele c.4976A>C E1659A directly targeted; extracellular concern flags mechanism.

## status

A-tier. Reframed from "NORMAL monotherapy" → **"MILD-MODERATE adjunct lever with conditional NORMAL"** per [[Misha Compound-Het Therapy Stack Model]]. K1141 pocket site-druggable (Phase 5c GREEN), chemistry-limited on current shortlist (Phase 5b RED); fenamic parents not developable for pediatric target per [[STRC h01 Fenamic Scaffold Tox Audit 2026-04-23]], tafamidis-style bioisosteric optimization is the forward path.

## active compute

```dataviewjs
const runs = dv.current().active_runs
if (!runs || runs.length === 0) dv.paragraph("_(no active runs)_")
else dv.table(
  ["Phase", "Job", "PID", "ETA"],
  runs.map(r => [r.phase, r.name, r.pid ?? "—", r.eta ?? "—"])
)
```

## next-step tree

After v3b delivery:
1. **GREEN non-covalent** → Phase 4h MD-scored validation on top-5 + wet-lab triage [[STRC h01 Phase 8 Wet-Lab Triage SOP]] on tafamidis-style optimized analog (parents are MD probes only per tox audit)
2. **GREEN covalent** → Phase 6c selectivity audit vs proteome-wide Lys pockets + cochlear ion-channel panel (TRPM4 / Cx50 / BK / KCNQ4 / TMEM16A)
3. **YELLOW** → Phase 3c v4 fragment-growing on best cluster OR Phase 6d different warhead class
4. **RED** → Phase 3c v5 de novo RFdiffusion pocket design

After Phase 5d + 5e delivery:
- Mutant K1141 stable (ΔΔG < 0.3 kcal/mol, Kd ratio < 2×) → WT-based docking validated, continue medchem on v3b hits
- Mutant pocket shifts → **Phase 3c v6**: re-screen on mutant ensemble, or pivot to h11 dimer-interface rescue
- Mutant reveals new cryptic pocket → Phase 3c v7 new-pocket screen

## Recent activity

```dataviewjs
const p = dv.page(dv.current().file.folder + "/log.md")
if (p && p.file.lists.length) dv.list(p.file.lists.slice(0, 5).map(l => l.text))
else dv.paragraph("_(no log entries)_")
```

## evidence

- [[STRC Pharmacochaperone Virtual Screen E1659A]] — main hypothesis; Phases 0-3
- [[STRC Pharmacochaperone K1141 Fragment Pocket]] — pocket characterization
- [[STRC Pharmacochaperone Phase 4 Plan]] — 7-gate validation ladder
- [[STRC h01 Phase 4h Tafamidis Playbook Library 2026-04-23]] — 30-compound bioisostere seed
- [[STRC h01 Phase 8 Wet-Lab Triage SOP]] — 3-gate wet-lab protocol
- [[STRC h01 Fenamic Scaffold Tox Audit 2026-04-23]] — scaffold developability
- [[STRC h01 Phase 5 MD Ensemble Rescoring 2026-04-23]] — Phase 5a/5b
- [[STRC h01 Phase 5c Cryptic Pocket Analysis 2026-04-23]] — site stability
- [[STRC h01 Phase 3c v2 Expanded Screen 2026-04-23]] — expanded dock

## scripts

See [[STRC Computational Scripts Inventory]] § Hypothesis 1.

## log

[[h01 log]]

## cross

- [[STRC Electrostatic Analysis E1659A]] — electrostatic input to pocket design
- [[STRC Mini-STRC Single-Vector Hypothesis]] — h03 is parallel S-tier, not fallback
- [[Misha Compound-Het Therapy Stack Model]] — h01 is Misha's only monotherapy-to-NORMAL route

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[about]` [[Misha]]
