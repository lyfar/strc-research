---
title: "LinkCraft: An interactive tool for the design of flexible linkers"
type: source
authors: ["Pajkos, M", "Barrera, P", "Clerc, I", "Zanon, C", "Bernardo, P", "Cortes, J"]
journal: "Journal of Molecular Biology"
date: 2026-04-17
pubmed_id: "42001978"
doi: "10.1016/j.jmb.2026.169814"
relevance_type: lateral
relevance_score: 4/10
tags: [protein-engineering, linker-design, multi-domain, computational, intrinsically-disordered]
status: unread
date_added: 2026-04-21
---

## What they found
LinkCraft is a computational tool for rational design of intrinsically disordered linkers (IDLs) in multi-domain proteins. It suggests optimal IDL lengths as a function of inter-domain distance, supports custom or auto-generated linker sequences based on desired physicochemical properties, and enables ensemble-based structural modeling of the complete multi-modular protein. The tool treats linkers as active, tunable determinants of molecular function rather than passive connectors.

## Lateral connection
Mini-STRC (S-tier) involves truncating stereocilin domains while preserving function. The remaining domains must maintain proper spatial relationships, and the truncation junctions are effectively new linker regions. LinkCraft could help design optimal junction sequences that maintain the inter-domain geometry of full-length STRC after domain deletion. Additionally, the STRC In Situ SpyCatcher Assembly hypothesis (B-tier) requires two STRC fragments to assemble with proper inter-domain spacing — linker design is critical for both the SpyCatcher/SpyTag junction and the overall protein architecture.

## Hypothesis suggested
Rational linker design at Mini-STRC truncation junctions could improve folding efficiency and functional preservation of the minimized protein. Maps to **STRC Mini-STRC Single-Vector Hypothesis** (S-tier) and **STRC In Situ SpyCatcher Assembly** (B-tier). Testable: use LinkCraft to model the truncation junction in Ultra-Mini-STRC and compare predicted ensemble dynamics against full-length domain-domain distance distributions from MD.

## What could be computed
LinkCraft modeling of the Ultra-Mini-STRC truncation junctions. Ensemble conformational sampling of Mini-STRC with optimized vs. native junction sequences. Inter-domain distance distributions compared to full-length STRC. For SpyCatcher assembly: optimal linker length between SpyCatcher/SpyTag and the STRC fragment termini.

## Links
- PubMed: https://pubmed.ncbi.nlm.nih.gov/42001978/
- DOI: https://doi.org/10.1016/j.jmb.2026.169814

## Connections

- `[source]` auto-indexed 2026-04-21 by [[strc-lit-watch]]
