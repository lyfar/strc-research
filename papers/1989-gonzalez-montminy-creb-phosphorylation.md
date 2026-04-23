---
date: 2026-04-23
type: source
tags: [strc, h05, paper, creb, phosphorylation, transcription, cAMP]
hypotheses: [h05]
title: "Cyclic AMP stimulates somatostatin gene transcription by phosphorylation of CREB at serine 133"
authors: ["Gonzalez GA", "Montminy MR"]
journal: Cell
year: 1989
volume: 59
pages: "675-680"
pubmed_id: "2573431"
doi: "10.1016/0092-8674(89)90013-5"
---

## What they found

Foundational paper identifying Ser133 as the PKA phosphorylation site on CREB and demonstrating that pCREB activates the somatostatin CRE. Established the PKA→CREB→CRE transcriptional cascade. Does not directly measure dephosphorylation kinetics of pCREB in cochlear hair cells.

## Numbers that matter

- PKA phosphorylates CREB at Ser133: confirmed
- Dephosphorylation half-life of pCREB: **t½ ~5–10 min in HeLa/COS cells** (later work by the Montminy lab and others; not directly measured in this 1989 paper)
- t½ 5–10 min → k_dephos = 0.0012–0.0023/s
- **The h05 model uses K_CREB_DEPHOS_S = 0.005/s (t½ 2.3 min) and cites this paper — a 2–4× discrepancy.** The model dephosphorylates CREB faster than Gonzalez & Montminy's data implies.
- Hair cell-specific CREB-P kinetics not measured anywhere; the 2–4× discrepancy may be acceptable given tissue differences, but should be flagged.

## Fit to h05

The citation is real and correct for the CREB phosphorylation mechanism. However, the specific rate constant K_CREB_DEPHOS_S = 0.005/s is faster than the cited t½ range implies. This should be explicitly noted as a deviation. A sensitivity analysis varying this parameter by 4× is warranted.

## Connections

- [[STRC Calcium Oscillation Acoustic Therapy]] — CREB dephosphorylation rate
- `[part-of]` [[calcium-oscillation]] (literature-params topic)
