---
date: 2026-04-21
type: note
tags: [strc, misha, meta, ranking, priority, decision-rubric]
sites: [strc]
---

# STRC Hypothesis Ranking

> MANDATORY before any STRC work. Read tier table. Update after every proof. Missing update = incomplete work.

Change log: [[STRC Hypothesis Ranking Log]]

## Update protocol (for agents)

1. Every proof's atomic note closes with `## Ranking delta`.
2. Append 1-line entry to [[STRC Hypothesis Ranking Log]] (date | delta | link). Do NOT add to this file.
3. Edit table row in place (current state only).
4. Update per-hypothesis `h{N}/log.md` with 1-line entry.

## Rubric

Axes 1-5: **Mech** (biology plausible), **Deliv** (cochlear delivery 3-5 yr), **Misha-fit** (compound het, age, OHC window).

Tier heuristic: `min(Mech, Deliv, Misha-fit)`. S = top 5 active. A = backburner. B = watch. C = paused, needs external catalyst. D = killed.

## Active register (2026-04-23)

| #   | Hypothesis                                         | Mech | Deliv | Misha | Tier    | Next step (≤6 words)                    | Hub         |
| --- | -------------------------------------------------- | ---- | ----- | ----- | ------- | --------------------------------------- | ----------- |
| 1   | [[STRC Pharmacochaperone Virtual Screen E1659A]]   | 3    | 4     | 4     | **A**   | Phase 5 MD ensemble MM-GBSA             | [[h01 hub]] |
| 2   | [[STRC Piezoelectric TM Bioelectronic Amplifier]]  | 2    | 1     | 5     | **B**   | Find real OHC ligand + FEM strain-share | [[h02 hub]] |
| 3   | [[STRC Mini-STRC Single-Vector Hypothesis]]        | 5    | 5     | 4     | **S**   | Order gBlock, clone pAAV, coIP          | [[h03 hub]] |
| 4   | [[STRC mRNA-LNP Strategy B Full-Length]]           | 3    | 2     | 2     | **B**   | Academic/interim if AAV fails           | [[h04 hub]] |
| 5   | [[STRC Calcium Oscillation Acoustic Therapy]]      | 3    | 4     | 2     | **A**   | SPR/BLI E1659A affinity penalty         | [[h05 hub]] |
| 6   | [[STRC mRNA Therapy Hypothesis]]                   | 3    | 2     | 2     | **A**   | Adjunctive to Strategy B                | [[h06 hub]] |
| 7   | [[Prime Editing for STRC]]                         | 3    | 2     | 1     | **C**   | Wet-lab: ribosome-profiling STRCP1      | [[h07 hub]] |
| 8   | [[STRC ASO Exon Skipping]]                         | 2    | 3     | 2     | **C**   | Phase 3a morpholino + 3b gapmer         | [[h08 hub]] |
| 9   | [[STRC Synthetic Peptide Hydrogel HTC]]            | 4    | 3     | 3     | **A**   | 4b + blend model + lit-backed re-run    | [[h09 hub]] |
| 10  | [[STRC In Situ SpyCatcher Assembly]]               | 2    | 2     | 2     | **C**   | Wet-lab only if #3 fails                | [[h10 hub]] |
| 11  | [[STRC Engineered TECTA Chimera]]                  | 2    | 2     | 2     | **C**   | Alt scaffold deferred, S-tier bandwidth | [[h11 hub]] |
| 12  | [[Sonogenetic STRC Computational Proof]]           | 2    | 2     | 2     | **C**   | Speculative, no near-term path          | [[h12 hub]] |
| 13  | [[STRC Programmable Recombinases]]                 | 2    | 1     | 3     | **C**   | Technology-watch only                   | [[h13 hub]] |
| 14  | [[STRC Protein Replacement Therapy]]               | 2    | 1     | 3     | **C**   | No delivery route                       | [[h14 hub]] |
| 15  | [[STRC OTOA Paralog Cross-Rescue]]                 | 2    | 2     | 2     | **D/C** | Chimera dead; upregulation weak         | [[h15 hub]] |
| 16  | [[STRC ZP Domain Prion-Like Seeding]]              | 3    | 2     | 1     | **D**   | Wrong patient (paternal 98 kb Δ)        | [[h16 hub]] |
| 17  | [[STRC RBM24 Regulatory Hypothesis]]               | —    | —     | —     | ref     | Rolled into Strategy A                  | —           |
| 18  | [[STRC Stereocilia Bundle Mechanics Model]]        | —    | —     | —     | ref     | Supporting model                        | —           |
| 19  | [[STRC AAV Vector Design]]                         | —    | —     | —     | ref     | Mini-STRC implementation                | —           |
| 20  | [[STRC B8 Enhancer Selection]]                     | —    | —     | —     | ref     | AAV pipeline                            | —           |
| 21  | [[STRC Dual-Vector vs Single-Vector Transduction]] | —    | —     | —     | ref     | Engineering layer                       | —           |
| 22  | [[STRC Anti-AAV Immune Response Model]]            | —    | —     | —     | ref     | Re-dosing constraint                    | —           |
| 23  | [[STRC Electrostatic Analysis E1659A]]             | —    | —     | —     | ref     | Pharmacochaperone input                 | —           |
| 24  | [[Adult Treatment Window STRC]]                    | —    | —     | —     | ref     | Constraint all mechanisms               | —           |
| 25  | [[Alternative STRC Delivery Hypotheses]]           | —    | —     | —     | ref     | Delivery-layer survey                   | —           |
| 26  | [[STRC Engineered Homodimer Avidity]]              | 3    | 5     | 4     | **B**   | Phase 1c 5 Å contact re-cluster         | [[h26 hub]] |
| 27  | [[STRC STRCP1 Activation Rescue]]                  | 2    | 2     | 4     | **C**   | Email Holt lab ribosome-profiling       | [[h27 hub]] |

## S-tier (active compute now)

- **#3 Mini-STRC AAV** — single-dose, null-compatible. Iranfar 2026 + Holt 2021 validated. Regeneron AAV.104 commercial track. Shanghai Shu Yilai knock-in mouse active.
- (h02 Piezo demoted S→B after 2026-04-23 post-fix re-run: baseline no longer passes 60 dB after d31 reconcile; A666 ligand phantom; TM mismatch 10⁵× not modelled.)

## Kill list (D-tier)

- OTOA chimera — Cα RMSD 13.8 Å.
- ZP prion-like seeding — requires WT substrate; Misha paternal 98 kb Δ has none.
- Direct protein replacement — no delivery route.

## Decision gates

- **→ D**: mechanism falsified OR Misha-fit=1 OR delivery impossible 10 yr.
- **→ C**: no near-term computational step; log external catalyst.
- **→ B**: evidence neutral, shallow engagement.
- **→ A**: real progress possible, not highest-leverage today.
- **→ S**: top 5 max, advances Misha directly, actionable next step.

## Literature-audit status (2026-04-23)

Per [[feedback_literature_first]] rule + [[AGENTS §0c]]: before any computational proof, parameter-provenance must be verified against primary literature. Audit status per hypothesis:

| # | Hypothesis | Tier | Lit audit | Notes |
|---|---|---|---|---|
| 1 | [[STRC Pharmacochaperone Virtual Screen E1659A]] | A | ✅ **FIXED 2026-04-23** | 0 phantoms. Housekeeping closed: TPSA docstring corrected (descriptor-only, STRC extracellular so CNS bracket irrelevant); druggability cross-phase incomparability flagged in 3 functions. See [[STRC h01 Parameter Provenance Audit 2026-04-23]]. Defensible. |
| 2 | [[STRC Piezoelectric TM Bioelectronic Amplifier]] | **B** (was S) | ✅ **FIXED 2026-04-23** | Phantoms pinned as PHANTOM in-code; d31 reconciled to −12 pC/N; TM mismatch flag in JSON; Mech 3→2, Deliv 2→1; see [[STRC h02 Parameter Provenance Audit 2026-04-23]] |
| 3 | [[STRC Mini-STRC Single-Vector Hypothesis]] | S | 🔒 **DEFERRED by Egor directive 2026-04-23** | Do not audit without re-authorization. Likely inherits STRC×TMEM145 Kd gap. |
| 5 | [[STRC Calcium Oscillation Acoustic Therapy]] | A | ✅ **FIXED 2026-04-23** | Phantoms removed/flagged in-code (Wu 2011, Sharma 2018, Cha 2010, Krey 2015); PKA K_cAMP cite corrected (Surdo 2017); CREB-P dephos 2-4× faster than lit flagged; STRC mRNA/protein t½ cross-script mismatches flagged. Phase 3 topological result robust. Tier A held. See [[STRC h05 Parameter Provenance Audit 2026-04-23]] |
| 9 | [[STRC Synthetic Peptide Hydrogel HTC]] | A | ⚠ **3 critical lit gaps 2026-04-23** | See [[STRC h09 Parameter Provenance Audit 2026-04-23]] |
| 26 | [[STRC Engineered Homodimer Avidity]] | B | ✅ **FIXED 2026-04-23** | Scripts already clean; prose phantoms (100 nM Kd from ipTM, 1 nM derived, 100-1000× absolute) reframed as ratio-on-unmeasured-baseline; status banner added; blocked at SPR/BLI + Phase-1-AF3-fail |
| 4 | [[STRC mRNA-LNP Strategy B Full-Length]] | B | ✅ **AUDITED 2026-04-23 Batch 2** | 0 phantoms in scripts. Inherits h06 Hill/t½ risk. All suspicious constants explicitly labeled hypothetical or estimate. |
| 6 | [[STRC mRNA Therapy Hypothesis]] | A | ✅ **PARTIAL FIX 2026-04-23** | STRC_HL_D + RBM24_HL_D lit-anchored (Zhang 2012 PMID 22246323 + Mauriac 2024 PMID 39320919 + Schwanhäusser 2011 PMID 21593866). Gao 2020 PMID 32493791 PHANTOM documented. Hill K_M/n/max_boost remain CIRCULAR FIT — no OHC splicing dose-response lit exists. See [[STRC h06 Parameter Provenance Audit 2026-04-23]]. |
| 7 | [[Prime Editing for STRC]] | C | ✅ **7/8 CONFIRMED 2026-04-23 Batch 2** | Same Gao 2020 PMID phantom as h06 (now corrected in prose); other 7 cites verified real (Chen 2024, Chemla 2025, Fang 2021, Villiger 2021, Anzalone 2019, Kim 2023, Zhang 2025). Best-cited note in portfolio. |
| 8 | [[STRC ASO Exon Skipping]] | C | ✅ **CLEAN 2026-04-23 Batch 2** | Cleanest scripts in portfolio. Ensembl REST + SantaLucia 1998 NN (verified exact), 4/4 FDA-PMO precedents confirmed. |
| 27 | [[STRC STRCP1 Activation Rescue]] | C | ✅ **CLEAN 2026-04-23 Batch 2** | Live GTEx API query, zero embedded bio-constants. Correct GENCODE IDs. |
| 10–17 | SpyCatcher, TECTA, Sono, Recombinases, Protein Replacement, OTOA, ZP-prion, ref rows | C/D/ref | needs audit if promoted | Not scanned; flag before any new compute |

All 4 scanning agents target `research/strc/models/*.py` scripts belonging to their hypothesis, produce parameter-provenance tables, retrieve missing papers via Anna's Archive + open sources, MinerU-parse into `~/BookLibrary/mineru-output/`, write paper notes into `research/strc/papers/`, and append to `research/strc/literature-params/*.md` topic files.

## Connections

- `[part-of]` [[STRC Research Portal]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[see-also]` [[STRC Mini-STRC Single-Vector Hypothesis]]
- `[see-also]` [[STRC Pharmacochaperone Virtual Screen E1659A]]
- `[see-also]` [[STRC Piezoelectric TM Bioelectronic Amplifier]]
- `[see-also]` [[STRC Calcium Oscillation Acoustic Therapy]]
- `[see-also]` [[STRC Synthetic Peptide Hydrogel HTC]]
- `[about]` [[Misha]]
