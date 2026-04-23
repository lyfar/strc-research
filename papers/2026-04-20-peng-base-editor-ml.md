---
title: "Rapid customization of base editors via machine learning-powered combinatorial mutagenesis"
type: source
authors: ["Peng, Jiaxing", "Chan, Don C T", "Chu, Hoi Yee", "Fong, John H C", "Lam, Yin Kau", "Cheung, Maggie S H", "Leung, Wing", "Choi, Gigi C G", "Wong, Alan S L"]
journal: "Mol Cell"
date: 2026-04-16
pubmed_id: "41997156"
doi: "10.1016/j.molcel.2026.03.030"
relevance_type: lateral
relevance_score: 5/10
tags: [base-editing, machine-learning, protein-engineering, precision-editing, bystander-edits, prime-editing]
status: unread
date_added: 2026-04-20
---

## What they found
Combined combinatorial mutagenesis with ML to engineer base editors with motif-specific activity. Profiled 160,000 evoAPOBEC1 and 64 million TadA variants in human cells, using only 0.004% of the mutational landscape for predictions. Identified variants that eliminated residual adenine editing in cytosine base editors. Achieved undetectable bystander edits in 50% of >800 disease-associated mutations tested. A structure-based deep learning model predicted functional TadA variants with 63% success across 20^26 variants spanning 26 amino acid sites without experimental data.

## Lateral connection
STRC mutations causing DFNB16 include point mutations that could theoretically be corrected by base editing. The challenge has been bystander edits in the surrounding sequence context. This paper's ML-powered approach to customizing base editors for specific motifs means that for any given STRC point mutation, a custom editor with minimal off-target activity could be computationally designed. The 0.004% sampling efficiency makes this experimentally tractable. Additionally, the principle of ML-guided enzyme engineering applies directly to our prime editing hypothesis — the same approach could optimize prime editors for STRC correction efficiency.

## Hypothesis suggested
For prevalent STRC point mutations (e.g., missense mutations in conserved domains), ML-customized base editors designed for the specific sequence context could achieve precise correction with undetectable bystander editing, potentially delivered via AAV to outer hair cells.

## What could be computed
Catalog all known pathogenic STRC point mutations, determine which fall within base editor windows (C>T or A>G), and computationally predict optimal editor variants using the structure-based deep learning approach described. Estimate what fraction of DFNB16 patients could be addressed by customized base editors vs. requiring the mini-gene approach.

## Links
- PubMed: https://pubmed.ncbi.nlm.nih.gov/41997156/
- DOI: https://doi.org/10.1016/j.molcel.2026.03.030

## Connections

- `[source]` auto-indexed 2026-04-20 by [[strc-lit-watch]]
