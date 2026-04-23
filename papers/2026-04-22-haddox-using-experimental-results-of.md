---
title: "Using experimental results of protein design to guide biomolecular energy-function development."
authors: ["Haddox, Hugh K", "Rocklin, Gabriel J", "Motta, Francis C", "Strickland, Devin", "Halabiya, Samer F", "Cordray, Cameron", "Park, Hahnbeom", "Klavins, Eric", "Baker, David", "DiMaio, Frank"]
journal: "PLoS computational biology"
date: 2026-04-22
pubmed_id: "42018576"
doi: "10.1371/journal.pcbi.1014215"
relevance_type: lateral
relevance_score: 4/10
tags: [de-novo-design]
status: unread
date_added: 2026-04-23
publication_status: ahead-of-print
---

## What they found
Computational models of macromolecules have many applications in biochemistry, but physical inaccuracies limit their utility. One class of models uses energy functions rooted in classical mechanics. The standard datasets used to train these models are limited in diversity, pointing to a need for new training data. Here, we sought to explore a new paradigm for training an energy function, where the Rosetta energy function was used to design de novo proteins. Experimental results on these designs were then used to identify failure modes of design, which were subsequently used as a "guiding principle" to retrain the energy function. Specifically, we examined a diverse set of de novo protein designs experimentally tested for their ability to stably fold, identifying unstable designs that were predicted to be stable by the Rosetta energy function. Using deep mutational scanning, we identified single amino-acid mutations that rescued the stability of these designs, providing insight into common failure modes of the energy function. We identified one key failure mode, involving steric clashing in protein cores. We identified similar overpacking when using Rosetta to refine high-resolution protein crystal structures, quantified the degree of overpacking, and refit a small set of energy-function parameters to better recapitulate native-like packing. Following fitting, we largely eliminated the failure mode in the refinement task, while retaining performance on other benchmarks, resulting in an updated version of the Rosetta energy function. This work shows how learning from protein designs can guide energy-function development.

## How this applies to our program
Needs manual review. Relevance type auto-inferred as `lateral` (score 4/10).

## Key numbers
- TBD — auto-scan caught this paper. Numbers extraction deferred to MinerU full-text parse if PDF retrieved.

## Links
- pubmed_id: https://pubmed.ncbi.nlm.nih.gov/42018576/
- DOI: https://doi.org/10.1371/journal.pcbi.1014215

## Connections

- `[source]` auto-indexed 2026-04-23 by [[strc-lit-watch]]
