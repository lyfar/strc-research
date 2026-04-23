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
status: A-tier
---

# h01 — Pharmacochaperone E1659A

**mech.** Small-molecule fold-stabilizer binds K1141 pocket on E1659A STRC, rescuing misfolded stereocilin.
**delivery.** Oral/topical small molecule; round-window permeability precedent exists.
**patient-fit.** Maternal allele c.4976A>C E1659A directly targeted; extracellular concern flags mechanism.

## status

A-tier. Reframed from "NORMAL monotherapy" → **"MILD-MODERATE adjunct lever with conditional NORMAL"** per [[Misha Compound-Het Therapy Stack Model]]. K1141 pocket site-druggable (Phase 5c GREEN), chemistry-limited on current shortlist (Phase 5b RED); fenamic parents not developable for pediatric target per [[STRC h01 Fenamic Scaffold Tox Audit 2026-04-23]], tafamidis-style bioisosteric optimization is the forward path.

## active compute

- **Phase 3c v3b + 6b** fenamic+covalent 12 253-ligand dock (PID 73601 nohup, ETA ~05:10 local 2026-04-24)
- **Phase 5d** E1659A full-length mutant MD (PID 28875 nohup, ETA ~23:00 local 2026-04-23)
- **Phase 5e** mutant-ensemble re-dock script ready; runs post-5d + post-v3b-Stage-2

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

<!-- RECENT:START auto-synced from log.md by sync_strc_indices.py — do not hand-edit -->

## recent

- **Phase 4h Tafamidis-Playbook Library drafted**: 30-compound prioritized target list, 5 cores × 6 polar distal subs × 2 acid bioisosteres = 60 combinatorial pruned to top 30. Tier 1 commercial + Tier 2 2-step in-house = 8 compounds orderable this week (tafamidis #5, iododiflunisal #11/#12, benzoxazole-3CN #1, benzoxazole-3OMe #2). 3 Phase 6b covalent warheads (salicylaldehyde imine, α-cyanoacrylate Michael, acyl-hydrazone). Predicted ΔG improvement -0.5 to -2.0 kcal/mol over parent niflumic −6.18. 6-target cochlear channel off-target panel pre-scored (TRPM4/Cx50/BK/KCNQ4/TMEM16A/COX) queued as Phase 6c 30×6 = 180 docks, ~1 h local wall. A held. No empirical data yet; this seeds Phase 3c v4 virtual screen + Phase 8 wet-lab synthesis Stage 1. → [[STRC h01 Phase 4h Tafamidis Playbook Library 2026-04-23]]
- **Phase 5d E1659A MD LAUNCHED** (PID 28875 nohup, 19:34 local 2026-04-23, ETA ~23:00 local = 3.5 h wall): full-length AF3 mutant `job3-mutant.cif` chain A (1775 residues), verified K1141=LYS + E1659=ALA pre-run. 645,244-atom system (4× Phase 5a) solvated with TIP3P + 0.15 M NaCl, AMBER14SB, 2 ns production × 20 snapshots. SMOKE (50 ps × 2) passed at 13.86 ns/day. **Closes critical gap**: all prior Phase 5a/5b/5c were on WT Ultra-Mini (594-1294 full-length) which excludes E1659. Phase 5d is the first MD on the actual disease target. Phase 5e script written (mutant-ensemble Vina re-dock of LEGACY_LEADS + V2_HITS + tafamidis-analog, pocket box centre derived per-snapshot from K1141 + ring residue Cα centroid, 18×18×18 Å box, cpu=8 parallel). Phase 5e scheduled post-Phase 5d delivery + post-v3b Stage 2 (CPU schedule). Verdict deferred. → [[STRC h01 Phase 5d E1659A MD 2026-04-23]] expected
- **Fenamic Scaffold Tox Audit**: A held, scaffold-developability risk sharpened. Parent fenamates (niflumic / flufenamic / meclofenamic / mefenamic / tolfenamic) **NOT developable as-is for pediatric cochlear target** — ion-channel promiscuity (TRPM4 2.8 µM, Cx50 3 µM gap junctions, BK 25 µM, KCNQ-family, TMEM16A 12 µM) overlaps Kd ~30 µM window; mefenamic seizure risk ≥2.5 g; COX inhibition reduces cochlear blood flow (therapeutically antagonistic). Scaffold developable via tafamidis-style optimization (benzoxazole bioisostere + polar distal ring subs); direct precedent flufenamic → tafamidis → FDA 2019. **Phase 4h: parents as MD PROBES only**, not wet-lab triage candidates. → [[STRC h01 Fenamic Scaffold Tox Audit 2026-04-23]]
- **Phase 8 Wet-Lab Triage SOP drafted**: 3-gate protocol ready for GREEN-scenario from v3b. G1 ThermoFluor ΔTm (iododiflunisal positive, 3-day ~$3k), G2 MST Kd (2-wk ~$8k), G3 cryoEM + HEK293-FL-E1659A-GFP surface trafficking rescue (4-wk ~$25k) + parallel cochlear ion-channel tox pre-screen (TRPM4, Cx50, BK, KCNQ4, TMEM16A). Decision tree through Phase 9 mouse OHC ex-vivo + CRO. → [[STRC h01 Phase 8 Wet-Lab Triage SOP]]
- **Phase 3c v3b + 6b LAUNCHED** (in progress, PID 73601 nohup detached, ETA ~5.7h): fenamic-focused 12,253-ligand library (619-lig combinatorial → expanded to 12k via 8 cores × 57 N-aryl subs × 5 acid bioisosteres × 8 ring subs + 6 covalent warheads on 8 cores × 10 hot N-aryls × 2 acid bioisosteres). Two-stage ensemble dock identical to Phase 3c v2b: Stage 1 on snap_008 exh 8 × 3 modes (0.65 lig/s, ETA 313 min), Stage 2 on top-50 × 5 k-means-selected conformers exh 16 × 5 modes. Checkpoint JSON every 100 lig. Expected delivery ~04:30 local. Output → `pharmacochaperone_phase3c_v3b_ensemble_dock.json`. Will produce proof note [[STRC h01 Phase 3c v3b Fenamic + Covalent Screen 2026-04-23]] on completion. → [[STRC h01 Phase 3c v2 Expanded Screen 2026-04-23]] (precedent)

<!-- RECENT:END -->

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
