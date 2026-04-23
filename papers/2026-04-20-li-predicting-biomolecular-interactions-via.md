---
title: "Predicting biomolecular interactions via a dual-stream graph neural network with motif constraint and diffusion-based regularization."
authors: ["Li, Danyu", "Huang, Rubing", "Zhou, Ling", "Tian, Jinyu", "Guo, Shikai", "Zou, Bin"]
journal: "Computational biology and chemistry"
date: 2026-04-20
pubmed_id: "42019100"
doi: "10.1016/j.compbiolchem.2026.109056"
relevance_type: lateral
relevance_score: 3/10
tags: [strc]
status: unread
date_added: 2026-04-23
publication_status: ahead-of-print
---

## What they found
Biomolecular interactions, such as RNA-Protein Interactions (RPIs) and Protein-Protein Interactions (PPIs), are fundamental to life; however, accurately predicting them remains a central challenge in computational biology. Current deep learning methods, despite their promise, often suffer from limited interpretability, poor generalization to novel biomolecules, and sensitivity to sparse and noisy data. To overcome these limitations, we propose a novel framework, DSG-BIP, for Biomolecular Interaction Prediction (BIP). DSG-BIP employs a Dual-Stream Graph (DSG) neural network that separately models the topological structure and node features of the biomolecular interaction network, and adaptively balances structural and functional features via a collaborative gated attention mechanism. The core innovation of DSG-BIP is the integration of learnable motif constraints and an improved diffusion-based regularization. Unlike prior methods that use motifs as static, precomputed features, DSG-BIP introduces learnable motif constraints that are dynamically optimized based on both sequence conservation and network context, enhancing biological interpretability and mitigating data sparsity. Concurrently, the improved diffusion-based regularization term stabilizes representation learning by modeling the distributional uncertainty in sparse and noisy data. Additionally, an adaptive masking mechanism with placeholder design enhances robustness against both data sparsity and class imbalance. Extensive experiments on multiple benchmark datasets for RPI and PPI demonstrate that DSG-BIP achieves prediction performance comparable to that of state-of-the-art methods, while offering improved interpretability and generalization.

## How this applies to our program
Needs manual review. Relevance type auto-inferred as `lateral` (score 3/10).

## Key numbers
- TBD — auto-scan caught this paper. Numbers extraction deferred to MinerU full-text parse if PDF retrieved.

## Links
- pubmed_id: https://pubmed.ncbi.nlm.nih.gov/42019100/
- DOI: https://doi.org/10.1016/j.compbiolchem.2026.109056

## Connections

- `[source]` auto-indexed 2026-04-23 by [[strc-lit-watch]]
