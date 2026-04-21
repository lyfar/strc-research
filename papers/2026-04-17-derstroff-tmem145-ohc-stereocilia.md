---
title: "TMEM145 is a principal component of outer hair cell stereocilia"
authors: ["Derstroff, D", "Flook, M", "Löhnes, A", "Kreye, P", "Newton, S", "Renigunta, V", "Hanemaaijer, S", "Aguilar, C", "Holt, JR", "Bowl, MR", "Oliver, D", "Reimann, K"]
journal: "Neuron"
date: 2026-03-31
pubmed_id: "41923617"
doi: "10.1016/j.neuron.2026.03.007"
relevance_type: direct
relevance_score: 9/10
tags: [strc, stereocilin, tmem145, tm-acs, tectorial-membrane, outer-hair-cell, stereocilia, protein-complex]
status: unread
date_added: 2026-04-17
---

## What they found

TMEM145 is a newly identified transmembrane protein that forms the core of TM-ACs (tectorial membrane attachment crowns) on outer hair cell stereocilia. Without TMEM145, stereocilin is lost from OHC stereocilia, the hair bundle physically disconnects from the tectorial membrane, and profound hearing loss results. TMEM145 appears to act as an anchor or scaffold that retains stereocilin at the stereocilia tip — not just a passenger in the complex.

## How this applies to mini-STRC

This changes the structural model for mini-STRC design. Stereocilin doesn't sit at the TM-AC alone — it depends on TMEM145 for localization. Mini-STRC (residues 700–1775) must preserve the TMEM145-binding interface, or the truncated protein will express but fail to localize correctly. This is the critical test: does the truncated region contain the TMEM145 interaction domain? Computational priority: AlphaFold-Multimer model of STRC-TMEM145 complex, then map which STRC residues form the interface, then verify those residues fall within our 700–1775 window.

## Key numbers

- Without TMEM145: stereocilin lost from OHC stereocilia entirely
- Result: hair bundle disconnects from tectorial membrane → profound HL
- TMEM145 is transmembrane (OHC-expressed), stereocilin is extracellular

## Links

- PubMed: https://pubmed.ncbi.nlm.nih.gov/41923617/
- DOI: https://doi.org/10.1016/j.neuron.2026.03.007

## Connections

- `[informs]` [[STRC Mini-STRC Single-Vector Hypothesis]] — TMEM145 as anchor; mini-STRC must preserve binding interface
- `[informs]` [[STRC Stereocilia Bundle Mechanics Model]] — TM attachment depends on TMEM145–stereocilin link
- `[see-also]` [[Derstroff et al 2026 TMEM145 Paper]] — full analysis note
