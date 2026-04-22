# STRC / DFNB16 — Problem Context for Research Routines

> **Agent context file.** Every scheduled routine working in this repo reads this first to understand the current research state before acting. Regenerated from `~/Brain/research/strc/_context.md` and `~/Brain/notes/STRC Hypothesis Ranking.md`.
>
> Last synced: 2026-04-22

---

## Core problem

STRC encodes stereocilin — an extracellular protein that maintains the structure of the hair bundle on outer hair cells (OHCs) of the cochlea. Stereocilin forms the horizontal top connectors between stereocilia and attaches the tallest stereocilia to the tectorial membrane. Without it, hair bundles collapse, OHCs lose mechanosensitivity, and hearing fails.

Mutations in STRC cause DFNB16 — autosomal recessive non-syndromic hearing loss. Target patient: **Misha** (compound heterozygote — paternal 98 kb deletion / maternal c.4976A>C E1659A).

## Central engineering constraint

STRC CDS is ~5,320 bp. AAV payload limit ~4,700 bp. STRC doesn't fit a single AAV. Every mechanism hypothesis addresses this: truncate, split, edit in situ, use a non-AAV chassis, or replace function via materials.

---

## Active hypothesis register

Source of truth: `~/Brain/notes/STRC Hypothesis Ranking.md`. Full table below.

Scoring: M = Mechanism, D = Delivery, F = Misha-fit (each 1-5). Tier ≈ min(M, D, F).
Tiers: **S** (primary, active compute, top-5 max) / **A** (active backburner) / **B** (watch) / **C** (paused) / **D** (killed) / **reference** (supporting models).


| #   | Hypothesis                                            | Mech | Deliv | Misha | Evidence                                               | Status                                  | Tier    | Next step                                                        |
| --- | ----------------------------------------------------- | ---- | ----- | ----- | ------------------------------------------------------ | --------------------------------------- | ------- | ---------------------------------------------------------------- |
| 1   | [[STRC Pharmacochaperone Virtual Screen E1659A]]      | 4    | 4     | 4     | Phases 0–3 done; top-5 leads shortlist; [[STRC Pharmacochaperone Phase 4 Plan]] queued; **Phase 4a PASS 4/5** pocket reproduces across CIFs ([[STRC Pharmacochaperone Phase 4a Pocket Reproducibility]]); **Phase 4e proxy soft-FAIL** (margin 0.15 vs 0.20); **Phase 4b smoke PASS on LE** — leads beat diflunisal 29–57% ([[STRC Pharmacochaperone Phase 4b Smoke Test]]); **Phase 4c FAIL** — all 5 leads prefer WT over E1659A by 0.40–0.55 kcal/mol (mean −0.455) vs diflunisal flat (−0.04) and negatives near-flat (−0.14) → aspect (a) "carboxylate replaces E1659" sub-mechanism flagged; aspect (b) "loop-capping fold-stabilizer" (VX-809 class, induced-fit) untested by rigid Vina ([[STRC Pharmacochaperone Phase 4c WT Decoy]]). Mechanism axis flagged pending Phase 4f interface-rescue MM-GBSA + Phase 5 MD | active                                  | **S**   | run Phase 4f interface-rescue MM-GBSA on Ultra-Mini × TMEM145 (skip Phase 4d — same rigid-receptor limitation) |
| 2   | [[STRC Piezoelectric TM Bioelectronic Amplifier]]     | 3    | 2     | 5     | Phases 1–3 done; delivery model 92% audiogram coverage | active                                  | **S**   | ex-vivo PVDF-TrFE deposition assay feasibility                   |
| 3   | [[STRC Mini-STRC Single-Vector Hypothesis]]           | 4→**5** | 4→**5**     | 4     | Full-stack 2026-04-21 Ultra-Mini validation: sub-Å TMEM145 interface preservation, 0 CpG CDS at 3.65% CAI, B8+WPRE3 vector fits, AF3 GOLD ipTM 0.68 (PAE 2.26 Å, 21/21 residues in zone), AF3 homodimer real-weak-interface (C2 94%, SASA 98%, ARM 1579-1581), full-vector CpG at AAV floor (13.9/kb, 67% from ITRs), **AF3 Ultra-Mini × TMEM145 full-length ipTM 0.43 with 23/41 contacts in GOLD zone and four of six canonical ARM clusters reproduced — matches full-STRC precedent, no regression from truncation**. **Delivery score officially upgraded 4 → 5; all three AF3 gates (GOLD pruned, homodimer, full TMEM145) cleared 2026-04-21** | active                                  | **S**   | order Ultra-Mini gBlock; clone pAAV B8-IgK-Ultra-Mini-WPRE3-bGH; Phase 4 HEK coIP wet-lab gate |
| 4   | [[STRC mRNA-LNP Strategy B Full-Length]]              | 3    | 2     | 5     | PK/PD + audiogram done 2026-04-21; cochlear LNP tropism literature scan 2026-04-22 ([[STRC LNP Cochlear Tropism Literature Scan]]) — Delivery 2/5 confirmed; 5% scenario plausible-but-unbuilt (prestin-A666-LNP-mRNA hybrid not yet published), 20% aspirational, 50% sci-fi; one direct LNP-mRNA cochlear precedent (Nazari 2020 ssPalm-BDNF guinea pig, OHC% not quantified) | active                                  | **S**   | (a) extend Strategy B PK/PD model with 1-3% OHC tropism row + ER/UPR ~5× soft ceiling; (b) ask Holt for unpublished cochlear LNP groups |
| 5   | [[STRC Calcium Oscillation Acoustic Therapy]]         | 3    | 4     | 2     | Phases 1–3 done; n=4.3 Hill, globally stable           | active                                  | **A**   | maternal-allele-only; Touch Grass integration path               |
| 6   | [[STRC mRNA Therapy Hypothesis]] (Strategy A / RBM24) | 3    | 2     | 2     | PK/PD + audiogram done                                 | backburner                              | **A**   | adjunctive to Strategy B; alone subtherapeutic for Misha         |
| 7   | [[Prime Editing for STRC]]                            | 4    | 2     | 2     | Phases 1–3 pegRNA design done; Phase 4 Cas-OFFinder scan 2026-04-22 ([[STRC PE Phase4 STRCP1 Paralog Off-Target]]) — **0 / 61 candidates pass**, every PE3b discriminating pegRNA hits STRCP1 paralog at chr15:43,700,346-363 (single 18 bp window, 100 kb downstream of STRC) with 0 mismatches. Phase 3 design was STRC-vs-STRC discriminating but blind to STRCP1 paralog. Hypothesis intact, candidate pool = 0 pending Phase 3.5 STRCP1-aware redesign. | active                                  | **B**   | Phase 3.5 STRCP1-discriminating pegRNA redesign: pull STRCP1 sequence at variant-equivalent locus, filter Phase 3 candidates by ≥2 mismatches in seed region against STRCP1, re-run Phase 4 |
| 8   | [[STRC ASO Exon Skipping]]                            | 3    | 3     | 2     | Phase 1 design (54 candidates); local self-fold check **PASS 54/54** ([[STRC ASO Phase1 Fold Check]]); heteroduplex ΔG + RNAplfold accessibility still pending | backburner                              | **B**   | maternal-only; heteroduplex ΔG (NUPACK) + RNAplfold sliding-window accessibility on full intron/exon context |
| 9   | [[STRC Synthetic Peptide Hydrogel HTC]]               | 3    | 3     | 3     | none                                                   | active                                  | **B**   | Phase 1 self-assembly geometry model                             |
| 10  | [[STRC In Situ SpyCatcher Assembly]]                  | 3    | 2     | 3     | none                                                   | active                                  | **B**   | Phase 1 fragment complementation geometry                        |
| 11  | [[STRC Engineered TECTA Chimera]]                     | 2    | 2     | 3     | none                                                   | active                                  | **B**   | Phase 1 TECTA-STRC fusion structural check                       |
| 12  | [[Sonogenetic STRC Computational Proof]]              | 2    | 2     | 2     | bifurcation analysis done                              | paused                                  | **C**   | speculative mechanism; no near-term path                         |
| 13  | [[STRC Programmable Recombinases]]                    | 2    | 1     | 3     | none                                                   | paused                                  | **C**   | technology-watch only                                            |
| 14  | [[STRC Protein Replacement Therapy]]                  | 2    | 1     | 3     | none                                                   | paused                                  | **C**   | no concrete delivery route beyond mRNA                           |
| 15  | [[STRC OTOA Paralog Cross-Rescue]]                    | 2    | 2     | 2     | Phases 1A–1B done; chimera killed                      | killed (chimera); paused (upregulation) | **D/C** | chimera dead (RMSD 13.8 Å); upregulation branch open but weak    |
| 16  | [[STRC ZP Domain Prion-Like Seeding]]                 | 3    | 2     | **1** | none                                                   | killed (Misha)                          | **D**   | requires WT substrate; Misha paternal 98 kb Δ has no WT to seed  |
| 17  | [[STRC RBM24 Regulatory Hypothesis]]                  | —    | —     | —     | parent of Strategy A                                   | reference                               | —       | rolled into Strategy A                                           |
| 18  | [[STRC Stereocilia Bundle Mechanics Model]]           | —    | —     | —     | supporting model                                       | reference                               | —       | used by other hypotheses                                         |
| 19  | [[STRC AAV Vector Design]]                            | —    | —     | —     | implementation helper                                  | reference                               | —       | used by Mini-STRC + Prime Editing                                |
| 20  | [[STRC B8 Enhancer Selection]]                        | —    | —     | —     | implementation helper                                  | reference                               | —       | used by AAV pipeline                                             |
| 21  | [[STRC Dual-Vector vs Single-Vector Transduction]]    | —    | —     | —     | implementation choice                                  | reference                               | —       | engineering layer of Mini-STRC                                   |
| 22  | [[STRC Anti-AAV Immune Response Model]]               | —    | —     | —     | supporting model                                       | reference                               | —       | constraint on re-dosing                                          |
| 23  | [[STRC Electrostatic Analysis E1659A]]                | —    | —     | —     | supporting model                                       | reference                               | —       | input to pharmacochaperone design                                |
| 24  | [[Adult Treatment Window STRC]]                       | —    | —     | —     | reference note                                         | reference                               | —       | constraint across all mechanisms                                 |
| 25  | [[Alternative STRC Delivery Hypotheses]]              | —    | —     | —     | sibling to mRNA strategies                             | reference                               | —       | delivery-layer survey                                            |


---

## Hypothesis → synthesis file slug

When Part 2 (hypothesis synthesis) needs to write a Recent Papers entry, use these slugs. `null` = no synthesis file yet; routine should skip synthesis for that hypothesis.

| Hypothesis | Synthesis file |
|---|---|
| STRC Pharmacochaperone Virtual Screen E1659A | `null` |
| STRC Piezoelectric TM Bioelectronic Amplifier | `null` |
| STRC Mini-STRC Single-Vector Hypothesis | `hypotheses/mini-strc.md` |
| STRC mRNA-LNP Strategy B Full-Length | `null` |
| STRC Calcium Oscillation Acoustic Therapy | `null` |
| STRC mRNA Therapy Hypothesis | `null` |
| Prime Editing for STRC | `hypotheses/prime-editing.md` |
| STRC ASO Exon Skipping | `hypotheses/exon-skipping.md` |
| STRC Synthetic Peptide Hydrogel HTC | `null` |
| STRC In Situ SpyCatcher Assembly | `null` |
| STRC Engineered TECTA Chimera | `null` |
| Sonogenetic STRC Computational Proof | `hypotheses/sonogenetics.md` |
| STRC Programmable Recombinases | `hypotheses/recombinases.md` |
| STRC Protein Replacement Therapy | `null` |
| STRC OTOA Paralog Cross-Rescue | `null` |
| STRC ZP Domain Prion-Like Seeding | `null` |
| STRC RBM24 Regulatory Hypothesis | `null` |
| STRC Stereocilia Bundle Mechanics Model | `hypotheses/bundle-mechanics.md` |
| STRC AAV Vector Design | `null` |
| STRC B8 Enhancer Selection | `hypotheses/enhancer.md` |
| STRC Dual-Vector vs Single-Vector Transduction | `hypotheses/dual-vector.md` |
| STRC Anti-AAV Immune Response Model | `hypotheses/immune.md` |
| STRC Electrostatic Analysis E1659A | `hypotheses/electrostatics.md` |
| Adult Treatment Window STRC | `null` |
| Alternative STRC Delivery Hypotheses | `hypotheses/delivery.md` |

---

## Kill list (D-tier — do not pursue)

- **OTOA chimera** — Phase 1B killed it (Cα RMSD 13.8 Å; paralog structure too divergent).
- **ZP prion-like seeding** — requires WT STRC substrate; Misha paternal 98 kb deletion has no WT to seed. Wrong patient.
- **Direct protein replacement** (without a delivery route) — no concrete pathway; supplanted by Strategy B mRNA.

---

## How to use this file (agent directive)

1. **Before scanning**: read this entire file. Understand S-tier (top compute), A-tier (backburner), killed. Skew relevance scoring toward S/A — papers advancing S-tier hypotheses are HIGH value. Papers on D-tier should only be indexed if they credibly flip the kill.
2. **For paper dedup**: see routine Step 1 (builds DOI/PMID/bioRxiv lookup sets).
3. **For hypothesis synthesis** (Part 2): match papers to hypotheses by tag + relevance + mechanism overlap. Append "Recent Papers" entries only to files listed in the slug table. If `null`: log `SKIP-NO-SYNTH-FILE: <hypothesis>` and move on.
4. **If this file feels stale** (S-tier hypothesis with completed next-step, or paper mentions hypothesis not in the table): operate from this file anyway, flag drift in commit message.
