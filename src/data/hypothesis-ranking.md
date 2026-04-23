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

| #   | Hypothesis                                         | Mech | Deliv | Misha | Tier    | Next step (≤6 words)                     | Hub         |
| --- | -------------------------------------------------- | ---- | ----- | ----- | ------- | ---------------------------------------- | ----------- |
| 1   | [[STRC Pharmacochaperone Virtual Screen E1659A]]   | 3    | 4     | 4     | **A**   | Phase 5 MD ensemble MM-GBSA              | [[h01 hub]] |
| 2   | [[STRC Piezoelectric TM Bioelectronic Amplifier]]  | 2    | 1     | 5     | **B**   | Find real OHC ligand + FEM strain-share  | [[h02 hub]] |
| 3   | [[STRC Mini-STRC Single-Vector Hypothesis]]        | 5    | 5     | 4     | **S**   | Order gBlock, clone pAAV, coIP           | [[h03 hub]] |
| 4   | [[STRC mRNA-LNP Strategy B Full-Length]]           | 3    | 2     | 2     | **B**   | Academic/interim if AAV fails            | [[h04 hub]] |
| 5   | [[STRC Calcium Oscillation Acoustic Therapy]]      | 3    | 4     | 2     | **A**   | SPR/BLI E1659A affinity penalty          | [[h05 hub]] |
| 6   | [[STRC mRNA Therapy Hypothesis]]                   | 2    | 2     | 2     | **B**   | Retrieve OHC RBM24 titration lit         | [[h06 hub]] |
| 7   | [[Prime Editing for STRC]]                         | 3    | 2     | 1     | **C**   | Wet-lab: ribosome-profiling STRCP1       | [[h07 hub]] |
| 8   | [[STRC ASO Exon Skipping]]                         | 2    | 3     | 2     | **C**   | Phase 3a morpholino + 3b gapmer          | [[h08 hub]] |
| 9   | [[STRC Synthetic Peptide Hydrogel HTC]]            | 4    | 3     | 3     | **A**   | Phase 2c WH2 bundling (gate #1) + 4b AF3 | [[h09 hub]] |
| 10  | [[STRC In Situ SpyCatcher Assembly]]               | 2    | 2     | 2     | **C**   | Wet-lab only if #3 fails                 | [[h10 hub]] |
| 11  | [[STRC Engineered TECTA Chimera]]                  | 2    | 2     | 2     | **C**   | Alt scaffold deferred, S-tier bandwidth  | [[h11 hub]] |
| 12  | [[Sonogenetic STRC Computational Proof]]           | 2    | 2     | 2     | **C**   | Speculative, no near-term path           | [[h12 hub]] |
| 13  | [[STRC Programmable Recombinases]]                 | 2    | 1     | 3     | **C**   | Technology-watch only                    | [[h13 hub]] |
| 14  | [[STRC Protein Replacement Therapy]]               | 2    | 1     | 3     | **C**   | No delivery route                        | [[h14 hub]] |
| 15  | [[STRC OTOA Paralog Cross-Rescue]]                 | 2    | 2     | 2     | **D/C** | Chimera dead; upregulation weak          | [[h15 hub]] |
| 16  | [[STRC ZP Domain Prion-Like Seeding]]              | 3    | 2     | 1     | **D**   | Wrong patient (paternal 98 kb Δ)         | [[h16 hub]] |
| 17  | [[STRC RBM24 Regulatory Hypothesis]]               | —    | —     | —     | ref     | Rolled into Strategy A                   | —           |
| 18  | [[STRC Stereocilia Bundle Mechanics Model]]        | —    | —     | —     | ref     | Supporting model                         | —           |
| 19  | [[STRC AAV Vector Design]]                         | —    | —     | —     | ref     | Mini-STRC implementation                 | —           |
| 20  | [[STRC B8 Enhancer Selection]]                     | —    | —     | —     | ref     | AAV pipeline                             | —           |
| 21  | [[STRC Dual-Vector vs Single-Vector Transduction]] | —    | —     | —     | ref     | Engineering layer                        | —           |
| 22  | [[STRC Anti-AAV Immune Response Model]]            | —    | —     | —     | ref     | Re-dosing constraint                     | —           |
| 23  | [[STRC Electrostatic Analysis E1659A]]             | —    | —     | —     | ref     | Pharmacochaperone input                  | —           |
| 24  | [[Adult Treatment Window STRC]]                    | —    | —     | —     | ref     | Constraint all mechanisms                | —           |
| 25  | [[Alternative STRC Delivery Hypotheses]]           | —    | —     | —     | ref     | Delivery-layer survey                    | —           |
| 26  | [[STRC Engineered Homodimer Avidity]]              | 2    | 5     | 4     | **B**   | Phase 1c 5 Å contact re-cluster          | [[h26 hub]] |
| 27  | [[STRC STRCP1 Activation Rescue]]                  | 2    | 2     | 4     | **C**   | Email Holt lab ribosome-profiling        | [[h27 hub]] |
| 28  | [[Misha Compound-Het Therapy Stack Model]]         | —    | —     | —     | ref     | Clinical plan: parallel h01+h03 stack    | —           |

## S-tier (active compute now)

- **#3 Mini-STRC AAV** — single-dose, null-compatible. Iranfar 2026 + Holt 2021 validated. Regeneron AAV.104 commercial track. Shanghai Shu Yilai knock-in mouse active. **NEW 2026-04-23**: per [[Misha Compound-Het Therapy Stack Model]], ceiling-capped at MILD (30-40 dB ABR) for Misha's compound-het genotype — cannot reach NORMAL alone because non-transduced OHCs retain maternal E1659A + null paternal. "Meaningful rescue" endpoint, not "cure" endpoint.
- (h02 Piezo demoted S→B after 2026-04-23 post-fix re-run: baseline no longer passes 60 dB after d31 reconcile; A666 ligand phantom; TM mismatch 10⁵× not modelled.)

**Clinical plan for Misha** per [[Misha Compound-Het Therapy Stack Model]]: pursue #1 (PC for maternal E1659A) + #3 (AAV for paternal null) in parallel. h01 is the only monotherapy route to NORMAL (≤ 25 dB ABR); stack provides redundancy and lower drug/AAV burden. Phase 5 MD on #1 is the critical-path compute that produces the f_PC parameter.

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
| 9 | [[STRC Synthetic Peptide Hydrogel HTC]] | A | ⚠ **Gate 3 NARROWED, 2 gates remain 2026-04-23** | Gate 1 (118 aa tail > 12 aa limit) CLOSED via 5% blend scaffold model (φ_tail=0.05 PASS, 2000× TMEM145 valency excess). Gate 3 (STRC×TMEM145 Kd) NARROWED via Dunbrack 2025 ipSAE: UM × TMEM145 GOLD ipSAE 0.591 in lit binder-band (NFAT-CnB 0.55 … NFAT-CnA 0.78), placeholder 100 nM now band-consistent not arbitrary; absolute Kd still wet-lab-gated. Gate 2 (WH2×F-actin Kd) open. Path B (AF3 abs Kd calibration) lit-closed. See [[STRC ipSAE Cross-Complex Reassessment 2026-04-23]] + [[STRC h09 Phase 4e_v2 Blend Scaffold]] + [[STRC h09 Parameter Provenance Audit 2026-04-23]]. |
| 26 | [[STRC Engineered Homodimer Avidity]] | B | ✅ **FIXED 2026-04-23** + ⚠ **ipSAE 0.000 2026-04-23** | Scripts already clean; prose phantoms reframed as ratio-on-unmeasured-baseline; status banner added; blocked at SPR/BLI + Phase-1-AF3-fail. **Dunbrack 2025 ipSAE reassessment: UM homodimer ipSAE = 0 (zero interface residues pass PAE<10 Å) — ipTM 0.28-0.30 uplift now identified as disordered-region false positive**; Mech 3→2. Structural symmetry + aa 1579-1581 self-contacts legs survive. See [[STRC ipSAE Cross-Complex Reassessment 2026-04-23]]. |
| 4 | [[STRC mRNA-LNP Strategy B Full-Length]] | B | ✅ **AUDITED 2026-04-23 Batch 2** | 0 phantoms in scripts. Inherits h06 Hill/t½ risk. All suspicious constants explicitly labeled hypothetical or estimate. |
| 6 | [[STRC mRNA Therapy Hypothesis]] | **B** (was A) | ✅ **BOUNDS 2026-04-23** | STRC_HL_D + RBM24_HL_D lit-anchored (Zhang 2012 + Mauriac 2024 + Schwanhäusser 2011). Hill K_M/n/max_boost CIRCULAR FIT. **Hill-sensitivity sweep (125-pt grid) shows per-OHC ROBUST but cochlea-mean therapeutic ARCHITECTURALLY IMPOSSIBLE at LNP eff ≤ 20% regardless of Hill params** — Mech 3→2; tier A→B. See [[STRC h06 Hill Sensitivity Sweep 2026-04-23]] + [[STRC h06 Parameter Provenance Audit 2026-04-23]]. |
| 7 | [[Prime Editing for STRC]] | C | ✅ **7/8 CONFIRMED 2026-04-23 Batch 2** | Same Gao 2020 PMID phantom as h06 (now corrected in prose); other 7 cites verified real (Chen 2024, Chemla 2025, Fang 2021, Villiger 2021, Anzalone 2019, Kim 2023, Zhang 2025). Best-cited note in portfolio. |
| 8 | [[STRC ASO Exon Skipping]] | C | ✅ **CLEAN 2026-04-23 Batch 2** | Cleanest scripts in portfolio. Ensembl REST + SantaLucia 1998 NN (verified exact), 4/4 FDA-PMO precedents confirmed. |
| 27 | [[STRC STRCP1 Activation Rescue]] | C | ✅ **CLEAN 2026-04-23 Batch 2** | Live GTEx API query, zero embedded bio-constants. Correct GENCODE IDs. |
| 10 | [[STRC In Situ SpyCatcher Assembly]] | C | ✅ **CLEAN 2026-04-23** | No scripts; prose verified (Zakeri 2012 SpyCatcher real; k_on 10³ M⁻¹s⁻¹ consistent). See [[STRC h10 h11 Parameter Provenance Audit 2026-04-23]] |
| 11 | [[STRC Engineered TECTA Chimera]] | C | ✅ **FIXED 2026-04-23** | No scripts; 1 phantom ("Song 2021" prestin N-glycans) replaced with Matsuda 2004 PMID 15140192 + Zheng 2009 *JARO*. Other cites real (Verpy 2008, Anc80L65, AAV-ie, AAV-PHP.B). See [[STRC h10 h11 Parameter Provenance Audit 2026-04-23]] |
| 12–17 | Sono, Recombinases, Protein Replacement, OTOA, ZP-prion, ref rows | C/D/ref | needs audit if promoted | Not scanned; flag before any new compute |

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
