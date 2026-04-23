---
title: "Pharmacokinetic principles in the inner ear: Influence of drug properties on intratympanic applications"
type: source
authors: ["Salt AN", "Plontke SK"]
journal: "Hearing Research"
date: 2018-03-01
pubmed_id: "29551306"
pmc: "PMC6133771"
doi: "10.1016/j.heares.2018.01.017"
relevance_type: critical
relevance_score: 9/10
tags: [strc, pkpd, inner-ear, cochlea, pharmacokinetics, round-window, perilymph, intratympanic]
hypotheses: [h09]
---

## What they found

Comprehensive review of pharmacokinetic principles for intratympanic drug delivery. Uses dexamethasone-phosphate and triamcinolone-acetonide as examples to illustrate how molecular properties (logP, TPSA) govern RWM crossing and perilymph distribution.

Key argument: dexamethasone-phosphate is poorly suited for local ear therapy because its polar, water-soluble form crosses lipid membranes poorly and is eliminated rapidly from basal perilymph without reaching apical regions. Triamcinolone-acetonide has more favorable pharmacokinetics.

## Key findings — Middle ear clearance

Figure 4 shows time-course of drug concentration decline in the middle ear niche after a single drop application:
- Gentamicin: fell to **46% of initial at 83 min** (t½ ≈ 83 × ln2 / ln(100/46) ≈ 70–75 min)
- Dexamethasone-phosphate: fell to **10% at 93 min** (t½ ≈ 40 min, extremely fast)

No single mucociliary half-life is given for all substances — it is drug-dependent and dominated by Eustachian tube drainage, not just ciliary clearance. The model's K_CLEAR_MIDDLE_EAR = 0.35/h (t½ = 2 h) is **slower than all measured values** in this paper.

## Key findings — RWM permeability

No numerical permeability coefficients vs. MW table. Paper discusses lipophilicity (logP) and TPSA as drivers, not MW directly. Substances with high TPSA cross lipid membranes poorly. Does NOT provide a permeability value for 14 kDa peptides.

## Key findings — Cochlear aqueduct / perilymph clearance

Mentions "cochlear aqueduct provides the outlet for fluid" and CSF-perilymph exchange is important in rodents. No quantitative human aqueduct clearance rate provided.

## Critical gap for model

No explicit t½ for middle ear clearance of unprotected peptide. No RWM permeability value for 14 kDa peptide. Perilymph clearance half-time not given for humans.

## Numbers that matter

| Parameter | Value | Context |
|---|---|---|
| Middle ear t½ (gentamicin) | ~70–75 min (≈0.55/h) | Guinea pig, single drop |
| Middle ear t½ (dex-phosphate) | ~40 min (≈1.0/h) | Guinea pig, single drop |
| RWM permeability (coefficients) | Not tabulated by MW | Discussed qualitatively via logP/TPSA |

## Supplementary materials

None identified.

## Access

PMC open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC6133771/

## Connections

- `[source]` [[hydrogel_phase4e_cochlear_pkpd]] — K_CLEAR_MIDDLE_EAR parameter; this paper does NOT support 2h half-life
- [[STRC Hypothesis Ranking]] — h09 cochlear PK basis
