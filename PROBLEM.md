# STRC / DFNB16 — Problem Context for Research Routines

> **Agent context file.** Every scheduled routine working in this repo reads this first to understand the current research state before acting. Regenerated from `~/Brain/research/strc/_context.md` and `~/Brain/notes/STRC Hypothesis Ranking.md`.
>
> Last synced: 2026-04-21

---

## Core problem

STRC encodes stereocilin — an extracellular protein that maintains the structure of the hair bundle on outer hair cells (OHCs) of the cochlea. Stereocilin forms the horizontal top connectors between stereocilia and attaches the tallest stereocilia to the tectorial membrane. Without it, hair bundles collapse, OHCs lose mechanosensitivity, and hearing fails.

Mutations in STRC cause DFNB16 — autosomal recessive non-syndromic hearing loss. ~1 in 1,000 in consanguineous populations, one of the most common causes of moderate-to-severe genetic hearing loss globally.

Target patient: **Misha** (compound heterozygote — paternal 98 kb deletion / maternal c.4976A>C E1659A). Any hypothesis is evaluated for Misha-fit first.

## Central engineering constraint

STRC CDS is ~5,320 bp. Standard AAV payload limit: ~4,700 bp including promoter and regulatory elements. STRC doesn't fit a single AAV. Every mechanism hypothesis has to address this: truncate, split, edit in situ, use a non-AAV chassis, or replace function via materials.

## Current engineering gates

- **Mini-STRC AAV**: Ultra-Mini 1075-1775 CDS validated; 0 CpG at 3.65% CAI; B8+WPRE3 fits 4,543 bp vector. Next gate: final AF3 `strc_ultramini_x_tmem145_full` job (in progress), then gBlock order.
- **mRNA-LNP (Strategy B, full-length)**: Hill ceiling broken (2.18× cochlea-mean at 5% LNP). Works on Misha's paternal null allele. Engineering gate: LNP cochlear tropism ≥ 20% with OHC targeting.
- **Pharmacochaperone E1659A**: Virtual screen complete, top-5 leads shortlisted. Wet-lab gate: bind assay on salicylic acid / indole-3-acetic lead pair.
- **Piezoelectric TM amplifier**: 92% audiogram coverage modeled. Wet-lab gate: ex-vivo PVDF-TrFE deposition feasibility.

---

## Active hypothesis register (25 entries, 2026-04-21 snapshot)

Scoring: M = Mechanism (does biology work?), D = Delivery (can we get it into Misha's cochlea in 3-5 yrs?), F = Misha-fit (works for his specific compound het + age + OHC window?). Each 1-5.

Tier ≈ min(M, D, F). Tiers: **S** (primary, active compute now, top-5 max) / **A** (active backburner) / **B** (watch) / **C** (paused, external catalyst needed) / **D** (killed) / **reference** (supporting models, not standalone hypotheses).

| # | Hypothesis | Tier | M/D/F | Status | Evidence | Next step |
|---|---|---|---|---|---|---|
| 1 | STRC Pharmacochaperone Virtual Screen E1659A | **S** | 4/4/4 | active | Phases 0-3; top-5 shortlist | wet-lab bind assay |
| 2 | STRC Piezoelectric TM Bioelectronic Amplifier | **S** | 3/2/5 | active | Phases 1-3; 92% audiogram coverage | PVDF-TrFE deposition feasibility |
| 3 | STRC Mini-STRC Single-Vector Hypothesis | **S** | 4→5/4/4 | active | Ultra-Mini validated 2026-04-21 | await AF3 x_tmem145_full, order gBlock |
| 4 | STRC mRNA-LNP Strategy B Full-Length | **S** | 3/2/5 | active | PK/PD + audiogram done 2026-04-21 | LNP tropism literature scan |
| 5 | STRC Calcium Oscillation Acoustic Therapy | **A** | 3/4/2 | active | Phases 1-3; Hill n=4.3, globally stable | maternal-only; Touch Grass integration |
| 6 | STRC mRNA Therapy Hypothesis (Strategy A / RBM24) | **A** | 3/2/2 | backburner | PK/PD + audiogram done | adjunctive to Strategy B only |
| 7 | Prime Editing for STRC | **A** | 4/2/2 | active | Phase 1-3 pegRNA design | maternal-only; Cas-OFFinder scan |
| 8 | STRC ASO Exon Skipping | **B** | 3/3/2 | backburner | Phase 1 design; lead ASO | maternal-only; ViennaRNA fold check |
| 9 | STRC Synthetic Peptide Hydrogel HTC | **B** | 3/3/3 | active | none | Phase 1 self-assembly geometry |
| 10 | STRC In Situ SpyCatcher Assembly | **B** | 3/2/3 | active | none | Phase 1 fragment complementation |
| 11 | STRC Engineered TECTA Chimera | **B** | 2/2/3 | active | none | Phase 1 TECTA-STRC fusion check |
| 12 | Sonogenetic STRC Computational Proof | **C** | 2/2/2 | paused | bifurcation done | speculative; no near-term path |
| 13 | STRC Programmable Recombinases | **C** | 2/1/3 | paused | none | technology-watch |
| 14 | STRC Protein Replacement Therapy | **C** | 2/1/3 | paused | none | no delivery beyond mRNA |
| 15 | STRC OTOA Paralog Cross-Rescue | **D/C** | 2/2/2 | partially killed | Phases 1A-1B | chimera dead (RMSD 13.8 Å); upregulation weak |
| 16 | STRC ZP Domain Prion-Like Seeding | **D** | 3/2/1 | killed | none | no WT to seed for Misha paternal Δ |
| 17 | STRC RBM24 Regulatory Hypothesis | ref | — | reference | parent of Strategy A | rolled into Strategy A |
| 18 | STRC Stereocilia Bundle Mechanics Model | ref | — | reference | supporting model | used by others |
| 19 | STRC AAV Vector Design | ref | — | reference | implementation helper | used by Mini-STRC + PE |
| 20 | STRC B8 Enhancer Selection | ref | — | reference | implementation helper | used by AAV pipeline |
| 21 | STRC Dual-Vector vs Single-Vector Transduction | ref | — | reference | implementation choice | engineering layer |
| 22 | STRC Anti-AAV Immune Response Model | ref | — | reference | supporting model | constraint on re-dosing |
| 23 | STRC Electrostatic Analysis E1659A | ref | — | reference | supporting model | input to pharmacochaperone |
| 24 | Adult Treatment Window STRC | ref | — | reference | reference note | constraint across mechanisms |
| 25 | Alternative STRC Delivery Hypotheses | ref | — | reference | sibling to mRNA | delivery-layer survey |

---

## Hypothesis → synthesis file slug

When Part 2 (hypothesis synthesis) needs to write a "Recent Papers" note entry, use these slugs. `null` = no synthesis file yet; the routine should skip synthesis for that hypothesis (no file to update, no new file to create automatically — the human decides when to open a synthesis file).

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
| Sonogenetic STRC Computational Proof | `hypotheses/sonogenetics.md` |
| STRC Programmable Recombinases | `hypotheses/recombinases.md` |
| STRC Stereocilia Bundle Mechanics Model | `hypotheses/bundle-mechanics.md` |
| STRC B8 Enhancer Selection | `hypotheses/enhancer.md` |
| STRC Dual-Vector vs Single-Vector Transduction | `hypotheses/dual-vector.md` |
| STRC Anti-AAV Immune Response Model | `hypotheses/immune.md` |
| STRC Electrostatic Analysis E1659A | `hypotheses/electrostatics.md` |
| Alternative STRC Delivery Hypotheses | `hypotheses/delivery.md` |
| (everything else) | `null` |

---

## Kill list (D-tier — do not pursue)

- **OTOA chimera** — Phase 1B killed it (Cα RMSD 13.8 Å; paralog too divergent). Chimera branch dead. Upregulation branch paused.
- **ZP prion-like seeding** — requires WT STRC substrate; Misha paternal 98 kb deletion has no WT to seed. Wrong-patient.
- **Direct protein replacement (without a delivery route)** — no concrete pathway; supplanted by Strategy B mRNA.

Do not advance these without explicit authorization. Papers relevant to these are LOW priority unless they open a specifically unblocked sub-branch (e.g. upregulation for OTOA is paused, not killed).

---

## How to use this file (agent directive)

1. **Before scanning**: read this entire file. Understand which 4 hypotheses are S-tier (top compute), which are A (active backburner), which are killed. Skew relevance scoring toward S/A tier — papers that advance S-tier hypotheses are HIGH value. Papers on D-tier mechanisms should only be indexed if they open a path that flips the tier.
2. **For paper dedup**: see the routine's own Step 1 instructions (builds DOI/PMID/bioRxiv lookup sets from `papers/*.md`).
3. **For hypothesis synthesis** (Part 2): match new papers to hypotheses by tag + relevance_type + mechanism overlap. Append "Recent Papers" entries only to files listed in the slug table above. If a hypothesis has no synthesis file (`null` in the table), log `SKIP-NO-SYNTH-FILE: <hypothesis>` and move on — do not create new files.
4. **If this file feels stale** (e.g. a hypothesis you see mentioned in a recent paper isn't in the table, or an S-tier hypothesis has a next-step that was completed): still operate from this file. Flag the drift in your commit message so the human knows to resync.
