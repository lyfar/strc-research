---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h01]
status: A-tier
---

# h01 — Pharmacochaperone E1659A

**mech.** Small-molecule fold-stabilizer binds K1141 pocket on E1659A STRC, rescuing misfolded stereocilin.
**delivery.** Oral/topical small molecule; round-window permeability precedent exists.
**patient-fit.** Maternal allele c.4976A>C E1659A directly targeted; extracellular concern flags mechanism.

## status

A-tier held. **Phase 5 MD + ensemble rescoring DELIVERED 2026-04-23** (pipeline on local Mac, 62 ns/day, 2 ns MD + 20 snapshots + Vina re-docking of Phase 4b top-5 leads). Result: **RED-LIGHT — all leads f_PC < 0.10 at therapeutic [L]=10 μM** (best = diflunisal positive 0.083); Phase 4b single-structure was over-optimistic by 0.36-0.92 kcal/mol. See [[STRC h01 Phase 5 MD Ensemble Rescoring 2026-04-23]].

**Phase 5c cryptic pocket analysis DELIVERED 2026-04-23**: K1141 pocket GREEN-LIGHT on structural stability — Cα RMSF 0.62 Å vs global 1.23 Å (2× more rigid), local void volume 719-850 Å³ across 20 frames (does not collapse), no alt cavity > 152 Å³ on the protein. **Phase 5b RED-LIGHT is chemistry-limited not site-limited.** See [[STRC h01 Phase 5c Cryptic Pocket Analysis 2026-04-23]].

**Phase 3c v2 expanded virtual screen DELIVERED 2026-04-23**: 667-ligand two-stage ensemble dock (619 combinatorial scaffold × acid bioisostere + 29 curated FDA carboxylates; Stage 1 on snap_008 exh 8 + Stage 2 on top-30 × 5 k-means conformers exh 16; 13 min wall). **Top hits: niflumic-acid / flufenamic / sulfasalazine at Kd ≈ 30 μM, f_PC ≈ 0.125 at [L]=10 μM** — 1.7× tighter Kd than Phase 5b baseline (diflunisal 50 μM) but still RED vs NORMAL threshold. Fenamic-acid family (2-arylaminobenzoic) revealed as new scaffold direction. 0 GREEN / 0 YELLOW / 30 RED. To cross MILD-MODERATE f_PC ≥ 0.30 need ΔG ≈ -7.06 kcal/mol (0.88 below current best). See [[STRC h01 Phase 3c v2 Expanded Screen 2026-04-23]].

Hypothesis intact; h01 reframed from "NORMAL monotherapy path" to **"MILD-MODERATE adjunct lever with conditional NORMAL"**.

**Phase 3c v3b + 6b IN PROGRESS (launched 2026-04-23 10:49 UTC, PID 73601 nohup, ETA ~5.7 h)**: 12,253-ligand fenamic-focused expansion library (8 anthranilic/fenamic cores × 57 N-aryl subs × 5 acid bioisosteres × 8 ring subs = ~11,500 non-covalent + 8 cores × 10 hot aryls × 6 covalent warheads × 2 acids = ~700 covalent warhead analogs). Two-stage ensemble dock identical to Phase 3c v2b. Checkpoint JSON every 100 lig. Expected proof note [[STRC h01 Phase 3c v3b Fenamic + Covalent Screen 2026-04-23]].

**Next after v3b delivery**:
1. If GREEN non-covalent hit → Phase 4h MD-scored validation on top-5 + wet-lab triage design.
2. If GREEN covalent hit → Phase 6c selectivity audit vs proteome-wide Lys pockets (off-target cytotoxicity test).
3. If YELLOW → Phase 3c v4 fragment-growing on best cluster or Phase 6d different warhead class (thiol-Michael acceptor, isothiocyanate).
4. If RED → Phase 3c v5 de novo RFdiffusion-pocket ligand design (qualitative pivot).

## evidence

- [[STRC Pharmacochaperone Virtual Screen E1659A]] — main hypothesis note; Phases 0-3 complete
- [[STRC Pharmacochaperone K1141 Fragment Pocket]] — pocket characterization
- [[STRC Pharmacochaperone Phase 4 Plan]] — 7-gate validation ladder
- [[STRC Pharmacochaperone Phase 4a Pocket Reproducibility]] — PASS 4/5 CIFs
- [[STRC Pharmacochaperone Phase 4b Smoke Test]] — PASS on LE; leads beat diflunisal 29-57%
- [[STRC Pharmacochaperone Phase 4c WT Decoy]] — FAIL: leads prefer WT; aspect (a) flagged
- [[STRC Pharmacochaperone Phase 4e Off-Target Selectivity]] — soft-FAIL (margin 0.15 vs 0.20)
- [[STRC Pharmacochaperone Phase 4f Interface Rescue SMOKE]] — method inadequate; Phase 4f canceled
- [[STRC Pharmacochaperone Phase 4g Repurpose Screen]] — 4PBA mut-preferring (−0.32); branch later killed
- [[STRC Repurposed FDA Chaperone Branch]] — KILLED: STRC extracellular, no ER target for 4PBA
- [[STRC Pharmacochaperone Virtual Screen Ranked Leads]] — top-5 shortlist
- [[STRC Pharmacophore Model K1141 Pocket]] — Phase 3a pharmacophore

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § Hypothesis 1)
- `pharmacochaperone_target_prep.py` through `pharmacochaperone_phase4g_repurpose_screen.py` (16 scripts)
- `pharmacochaperone_phase5_md.py` — scaffold (GROMACS-based, superseded by 5a/5b OpenMM)
- `pharmacochaperone_phase5a_apo_md_smoke.py` — OpenMM + AMBER14SB apo MD, SMOKE + 2 ns production delivered 2026-04-23
- `pharmacochaperone_phase5b_ensemble_redock.py` — Vina ensemble re-docking on Phase 5a snapshots; f_PC estimates per lead
- `pharmacochaperone_phase5c_cryptic_pocket_detection.py` — K1141 stability (Cα RMSF, local void volume) + grid-based alt cavity scan on Phase 5a trajectory
- `pharmacochaperone_phase3c_v2a_library_build.py` — RDKit combinatorial library (25 scaffolds × 10 acid bioisosteres + 29 curated FDA drugs = 667 ligands in PDBQT)
- `pharmacochaperone_phase3c_v2b_ensemble_dock.py` — 2-stage Vina ensemble dock: Stage 1 breadth (snap_008 exh 8) + Stage 2 ensemble (top-30 × 5 k-means receptor conformers exh 16)

## log

[[h01 log]]

## cross

- [[STRC Electrostatic Analysis E1659A]] — electrostatic input to pocket design
- [[STRC Mini-STRC Single-Vector Hypothesis]] — h03 is active fallback if h01 mechanism fails

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[about]` [[Misha]]
