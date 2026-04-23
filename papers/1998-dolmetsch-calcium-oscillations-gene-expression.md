---
date: 2026-04-23
type: source
tags: [strc, h05, paper, calcium-oscillation, camkii, calcineurin, gene-expression, frequency-decoding]
hypotheses: [h05]
title: "Calcium oscillations increase the efficiency and specificity of gene expression"
authors: ["Dolmetsch RE", "Xu K", "Lewis RS"]
journal: Nature
year: 1998
volume: 392
pages: "933-936"
pubmed_id: "9582075"
doi: "10.1038/31960"
---

## What they found

In Jurkat T cells, Ca²⁺ oscillations (rather than sustained Ca²⁺ elevations) activate NF-κB, NF-AT, and Oct/OAP transcription factors with high efficiency and specificity. The key finding: different transcription factors have different threshold frequencies — NF-AT is activated at low frequencies; NF-κB requires higher frequency. CaMKII autophosphorylation accumulates preferentially with high-frequency oscillations; calcineurin (which dephosphorylates and activates NF-AT) is activated by lower-frequency sustained signals. This establishes the concept of frequency decoding of Ca²⁺ signals.

## Numbers that matter

- Optimal oscillation frequency for NF-κB: ~0.2–0.5 Hz
- Calcineurin (NF-AT pathway) activated by low-frequency/sustained Ca²⁺: EC50 ~ 200–500 nM
- CaMKII switches to sustained activity at higher frequencies due to autophosphorylation memory
- **Rate constants `k_on_CaN = 0.3/s` and `k_off_CaN = 0.05/s` in the h05 model are not directly from this paper** — they are fitted values consistent with the conceptual framework

## Fit to h05

The RBM24 ODE Phase 1 model is built on the Dolmetsch frequency-decoding concept. Dolmetsch 1998 confirms the biological plausibility of the CaMKII:CaN ratio as a Ca²⁺ frequency readout. However, the specific rate constants used in the model are estimates, not direct measurements from this paper. The paper does not study cochlear hair cells.

## Connections

- [[STRC Calcium Oscillation Acoustic Therapy]] — conceptual basis for RBM24 ODE frequency decoder
- `[see-also]` [[1994-stemmer-klee-calcineurin-dual-calcium]] — calcineurin kinetics
- `[part-of]` [[calcium-oscillation]] (literature-params topic)
