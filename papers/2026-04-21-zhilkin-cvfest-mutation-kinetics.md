---
title: "Guiding Peptide Kinetics via Collective-Variable Tuning of Free-Energy Barriers"
authors: ["Zhilkin, A", "Medaparambath, M", "Mendels, D"]
journal: "Journal of Chemical Theory and Computation"
date: 2026-04-20
pubmed_id: "42007551"
doi: "10.1021/acs.jctc.6c00418"
relevance_type: lateral
relevance_score: 4/10
tags: [molecular-dynamics, mutation-kinetics, protein-folding, computational, free-energy, collective-variables]
status: unread
date_added: 2026-04-21
---

## What they found
The CV-FEST framework uses Harmonic Linear Discriminant Analysis (HLDA) collective variables from short MD trajectories to predict how point mutations affect protein folding/unfolding kinetics — without exhaustive sampling. Validated on Chignolin miniprotein mutants, the HLDA-derived residue-level scores predict whether mutations accelerate or slow conformational transitions. The leading HLDA eigenvalue correlates significantly with transition rates across mutations. This enables kinetic prediction from minimal local sampling.

## Lateral connection
The STRC E1659A mutation likely alters the folding kinetics of stereocilin's local domain, potentially trapping it in a misfolded intermediate that prevents trafficking. CV-FEST could predict whether E1659A accelerates or slows the folding transition compared to wild-type, and quantify the kinetic barrier that a pharmacochaperone must overcome. This is complementary to our existing electrostatic analysis of E1659A and could inform the pharmacochaperone binding energy threshold needed for therapeutic rescue.

## Hypothesis suggested
CV-FEST analysis of the STRC region around position 1659 could quantify the folding kinetic penalty of E→A substitution and predict whether our pharmacochaperone leads have sufficient binding energy to compensate. Maps to **STRC Pharmacochaperone Virtual Screen E1659A** (S-tier) and **STRC Electrostatic Analysis E1659A** (reference). Testable: apply HLDA to short MD trajectories of WT vs. E1659A STRC domain to predict kinetic barrier difference.

## What could be computed
HLDA collective variables from WT and E1659A STRC domain MD. Free energy barrier estimation for local folding transition. Correlation with pharmacochaperone binding free energies from our virtual screen. Residue-level mutation sensitivity scores for positions adjacent to E1659 (could identify synergistic stabilizing mutations for protein engineering).

## Links
- PubMed: https://pubmed.ncbi.nlm.nih.gov/42007551/
- DOI: https://doi.org/10.1021/acs.jctc.6c00418
