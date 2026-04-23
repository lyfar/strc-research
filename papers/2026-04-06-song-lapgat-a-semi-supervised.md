---
title: "LapGAT: A Semi-Supervised Learning Framework for Drug-Target Interaction Prediction."
authors: ["Song, Lianjun", "Yuan, Wei", "Pei, Xinyu"]
journal: "Molecular informatics"
date: 2026-04-06
pubmed_id: "42011802"
doi: "10.1002/minf.70028"
relevance_type: lateral
relevance_score: 3/10
tags: [strc]
status: unread
date_added: 2026-04-23
publication_status: published
---

## What they found
Drug-target interaction (DTI) prediction is a fundamental task in the field of drug discovery, with direct implications for the identification of novel therapeutic candidates and the repositioning of existing drugs. However, the practical application of DTI prediction remains hindered by several persistent challenges, including data scarcity, the absence of reliable negative samples, and limited model generalization across diverse biological contexts. Addressing these limitations is crucial for developing robust and generalizable predictive frameworks. In this study, we present LapGAT, a semi-supervised framework combining graph-enhanced Laplacian regularized least squares (LapRLS) with a graph attention network (GAT) to address these issues. In the upstream stage, LapRLS fuses multiple drug-target similarity matrices, applies Laplacian regularization, and selects top- and bottom-scoring pairs as high-confidence positive and negative samples. In the downstream stage, a multilayer GAT learns from these pseudo-labeled interactions, capturing both local graph structures and nonlinear dependencies. We validate LapGAT on four target categories (enzymes, G-protein-coupled receptors (GPCRs), ion channels, nuclear receptors), demonstrating robust performance in computational and experimental validation. Molecular docking (AutoDock Vina) confirms the physical plausibility of top-ranked predictions, with binding affinities ranging from -4.5 to -7.3&#x2009;kcal/mol. Literature-based validation achieves accuracies of 81.3% (enzymes), 100% (GPCRs), 71.4% (ion channels), and 88.9% (nuclear receptors). This work offers a scalable and flexible computational tool for accelerating drug discovery efforts, with the potential for broad applicability across various therapeutic domains.

## How this applies to our program
Needs manual review. Relevance type auto-inferred as `lateral` (score 3/10).

## Key numbers
- TBD — auto-scan caught this paper. Numbers extraction deferred to MinerU full-text parse if PDF retrieved.

## Links
- pubmed_id: https://pubmed.ncbi.nlm.nih.gov/42011802/
- DOI: https://doi.org/10.1002/minf.70028

## Connections

- `[source]` auto-indexed 2026-04-23 by [[strc-lit-watch]]
