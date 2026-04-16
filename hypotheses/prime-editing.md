---
title: Prime Editing — Fix the Single Mutation In Situ
status: exploratory
stage: feasibility-analysis
priority: alternative
tags: [prime-editing, base-editing, gene-correction, CRISPR, PE3]
---

## Core claim

Instead of replacing the entire STRC gene (5,325 bp), correct just the single mutated base. Misha's mutation is c.4976A>C → p.E1659A: a C>A transversion. Only prime editing (PE3/PE3b) can make this correction.

## Why only prime editing works

- CBE (cytosine base editor): can only do C→T. Can't make C→A.
- ABE (adenine base editor): wrong direction entirely.
- PE3/PE3b: C→A included in pegRNA design. ✓

## PAM site analysis

Prime editing requires NGG PAM within ~15 bp of target. Found: PAM site 4 bp from variant (CC on reverse complement = NGG). pegRNA design is feasible.

## Reality check

- Prime editing not tested in inner ear hair cells in vivo
- Delivery problem: prime editor protein (~180 kDa) + pegRNA → large payload, hard to deliver
- Efficiency in post-mitotic cells (OHCs don't divide) unknown
- Off-target edits possible

## Comparison to mini-STRC

Mini-STRC (gene replacement) is currently the more validated path. Prime editing is theoretically more elegant (fixes the actual mutation, leaves rest of genome untouched) but faces unsolved delivery and efficiency challenges in cochlear cells.

## Open questions

- What PE3 efficiency is achievable in non-dividing OHCs?
- Can the editor + pegRNA be co-delivered in a single vector?
- What are off-target rates at this genomic locus?
