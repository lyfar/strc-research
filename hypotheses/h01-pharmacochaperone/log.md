---
type: hypothesis-log
hypothesis: h01
tags: [strc, h01]
---

# h01 log

## 2026-04-23

- **Phase 3c v3b + 6b LAUNCHED** (in progress, PID 73601 nohup detached, ETA ~5.7h): fenamic-focused 12,253-ligand library (619-lig combinatorial → expanded to 12k via 8 cores × 57 N-aryl subs × 5 acid bioisosteres × 8 ring subs + 6 covalent warheads on 8 cores × 10 hot N-aryls × 2 acid bioisosteres). Two-stage ensemble dock identical to Phase 3c v2b: Stage 1 on snap_008 exh 8 × 3 modes (0.65 lig/s, ETA 313 min), Stage 2 on top-50 × 5 k-means-selected conformers exh 16 × 5 modes. Checkpoint JSON every 100 lig. Expected delivery ~04:30 local. Output → `pharmacochaperone_phase3c_v3b_ensemble_dock.json`. Will produce proof note [[STRC h01 Phase 3c v3b Fenamic + Covalent Screen 2026-04-23]] on completion. → [[STRC h01 Phase 3c v2 Expanded Screen 2026-04-23]] (precedent)
- Phase 3c v2 expanded virtual screen: **RED overall, partial progress**. 667-ligand library (619 RDKit combinatorial scaffold × acid bioisostere + 29 curated FDA/literature drugs) × 2-stage ensemble docking (Stage 1 snap_008 exh 8 on all 667 → 8.5 min; Stage 2 top-30 × 5 k-means conformers exh 16 → 4 min, total 13 min wall). Top hits niflumic-acid / flufenamic-acid / sulfasalazine at **Kd ≈ 30 μM, f_PC ≈ 0.125 at [L]=10 μM** — 1.7× better Kd than Phase 5b diflunisal baseline. **Fenamic-acid family (2-arylaminobenzoic) = new scaffold direction**; tetrazole bioisostere competitive with COOH. 0 GREEN / 0 YELLOW / 30 RED. Ceiling analysis: to cross f_PC ≥ 0.30 (MILD-MODERATE) need ΔG ≈ -7.06 kcal/mol (0.88 below current best). A held, all scores unchanged. **Next-step reordered**: (1) Phase 3c v3 ZINC22 bioactives ~20k library, (2) Phase 3c v4 fragment-growing niflumic core, (3) Phase 6b reversible covalent Lys warhead. → [[STRC h01 Phase 3c v2 Expanded Screen 2026-04-23]]
- Phase 5c cryptic pocket analysis: **K1141 site GREEN-LIGHT**. Custom grid-cavity + Kabsch-aligned RMSF analysis on Phase 5a 20-snapshot trajectory (fpocket 4.0 brew build broken on qhull). K1141 pocket Cα RMSF 0.62 Å vs global 1.23 Å (2× more rigid than average); local void volume 719-850 Å³ across 20 frames (~15% CV, pocket does not collapse or shift geometry). Global cavity scan (snap_010, 26-direction burial ≥16/26): no alt cavity > 152 Å³, nearest at 18.7 Å (91 Å³). **Phase 5b RED-LIGHT confirmed chemistry-limited not site-limited.** Phase 3c v2 expanded screen targeted at K1141 with ensemble receptor docking (use all 20 Phase 5a snapshots) validated as next move. A held, all scores unchanged. → [[STRC h01 Phase 5c Cryptic Pocket Analysis 2026-04-23]]
- Phase 5a MD pipeline validated on local Mac: 62 ns/day OpenMM OpenCL (Metal backend) on 164 k-atom solvated Ultra-Mini × TMEM145 chain A. 2 ns production delivered 20 snapshots, 45 min wall time. Replaces deferred GROMACS/AmberTools scaffold which required A100 rental. Phase 5b Vina ensemble re-docking of Phase 4b top-5 leads + diflunisal positive control against 20 snapshots → **all f_PC at [L]=10 μM < 0.10** (best = diflunisal positive 0.083, best lead = naphthalene-2-COOH 0.047). Phase 4b single-structure was over-optimistic by +0.36 to +0.92 kcal/mol. RED-LIGHT: current shortlist insufficient for NORMAL-rescue monotherapy. **A held, scores unchanged** (mechanism intact, Misha fit intact; only shortlist insufficient). **Next step**: Phase 3c v2 expanded virtual screen (DrugBank FDA + ZINC22 carboxylate tranche + fragment-based; ensemble filter ΔG ≤ −7.5 kcal/mol). → [[STRC h01 Phase 5 MD Ensemble Rescoring 2026-04-23]]
- Misha compound-het stack integration: A held, scores unchanged, **framing strengthened**. Per [[Misha Compound-Het Therapy Stack Model]] binary-functional-OHC model, h01 PC is the **only monotherapy route to NORMAL (≤25 dB ABR) for Misha** — because PC reaches every OHC (including non-transduced ones where h03 AAV by definition does nothing), rescuing maternal E1659A to 0.5 × f_mat_treated ≥ θ. Required f_PC ≥ 0.50 (mild E1659A) to ≥ 0.75 (severe). Phase 5 MD ensemble MM-GBSA priority strengthened — f_PC is the critical-path deliverable. → [[Misha Compound-Het Therapy Stack Model]]
- Post-audit housekeeping: A held. TPSA docstring in `phase3b_virtual_screen.py` corrected (descriptor-only, not composite scoring; CNS TPSA bracket removed — STRC extracellular). Druggability cross-phase incomparability flagged in Phase 1/2/2b. No fabricated ranges added. → [[STRC h01 Parameter Provenance Audit 2026-04-23]]

## 2026-04-22

- S → A: mechanism 4→3; STRC extracellular protein, 4PBA ER mechanism mismatch → [[STRC Repurposed FDA Chaperone Branch]]
- Phase 4g: 4PBA first mut-preferring compound (−0.32 kcal/mol); repurposing branch opened then killed same day → [[STRC Pharmacochaperone Phase 4g Repurpose Screen]]
- Phase 4f SMOKE: ΔΔG 20,439 kcal/mol = method error; production canceled → [[STRC Pharmacochaperone Phase 4f Interface Rescue SMOKE]]

## 2026-04-21

- Phase 4a: PASS 4/5, K1141 pocket reproduces across CIFs → [[STRC Pharmacochaperone Phase 4a Pocket Reproducibility]]
- Phase 4b: PASS on LE; leads beat diflunisal 29-57% → [[STRC Pharmacochaperone Phase 4b Smoke Test]]
- Phase 4c: FAIL; all 5 leads prefer WT (mean −0.455 kcal/mol); aspect (a) flagged → [[STRC Pharmacochaperone Phase 4c WT Decoy]]
- Phase 4e: soft-FAIL margin 0.15 vs 0.20; re-opens as Phase 4b sub-test → [[STRC Pharmacochaperone Phase 4e Off-Target Selectivity]]
- Phase 4 plan written, 7-gate ladder → [[STRC Pharmacochaperone Phase 4 Plan]]
- Phases 0-3 complete; top-5 leads shortlist → [[STRC Pharmacochaperone Virtual Screen Ranked Leads]]

## Connections

- `[part-of]` [[h01 hub]]
