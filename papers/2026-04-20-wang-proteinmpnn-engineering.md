---
title: "Deep Learning-Enabled Engineering of a Hyper-Stable and Soluble MPB70-83 Antigen for Sensitive Bovine Tuberculosis Surveillance"
authors: ["Wang, Wen-Hao", "Ge, Jia-Zhen", "Xie, Ying-Ying"]
journal: "Microb Biotechnol"
date: 2026-04
pubmed_id: "41988819"
doi: "10.1111/1751-7915.70348"
relevance_type: lateral
relevance_score: 4/10
tags: [protein-engineering, proteinmpnn, deep-learning, protein-stability, molecular-dynamics, scaffold-design]
status: unread
date_added: 2026-04-20
---

## What they found
Used ProteinMPNN deep learning framework to redesign non-epitope scaffold regions of a fusion protein for enhanced thermodynamic stability. Molecular dynamics simulations confirmed the redesigned construct achieved a rigid, compact native state. The optimized protein achieved high-yield soluble expression in E. coli without inclusion bodies, demonstrating that computational scaffold redesign can dramatically improve protein stability and expression while preserving functional regions.

## Lateral connection
Mini-STRC truncation creates new domain junctions and exposed surfaces that may destabilize the protein. The ProteinMPNN approach of redesigning non-functional scaffold regions while preserving functional epitopes is directly applicable: the truncation junction in mini-STRC (where amino acids 699 and 1776 are joined) could be computationally redesigned using ProteinMPNN to optimize stability of the new interface without altering the functional N-terminal and C-terminal domains. The MD validation step provides a framework for computationally vetting designs before expensive in vivo testing.

## Hypothesis suggested
ProteinMPNN-guided redesign of the truncation junction and exposed surfaces in mini-STRC could produce a more thermodynamically stable protein than simple truncation, improving folding efficiency, secretion, and ultimately stereocilia localization.

## What could be computed
(1) Generate AlphaFold structure of mini-STRC with the 700-1775 truncation. (2) Use ProteinMPNN to redesign 5-10 residues flanking the junction while constraining the N-terminal and C-terminal functional domains. (3) Run MD simulations comparing stability of wild-type junction vs. ProteinMPNN-optimized junction. (4) Predict aggregation propensity using tools like CamSol.

## Links
- PubMed: https://pubmed.ncbi.nlm.nih.gov/41988819/
- DOI: https://doi.org/10.1111/1751-7915.70348
