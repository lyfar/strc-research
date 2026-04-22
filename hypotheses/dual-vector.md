---
title: Dual-Vector vs Single-Vector Transduction Analysis
status: active
stage: computational-modeling
priority: supporting
tags: [aav, dual-vector, single-vector, transduction, statistics, gamma-poisson, strc, synthesis, hypothesis]
date: 2026-04-22
type: synthesis
---

## Core claim

Mini-STRC (single vector) is mathematically superior to dual-vector approaches for STRC delivery. At standard clinical titer (3.75×10¹² GC/mL), single-vector delivers functional protein to 2.8x more hair cells than dual-vector.

## The dual-vector problem

STRC CDS = 5,325 bp. AAV limit = 4,700 bp. Dual-vector splits CDS across two AAVs. Problem: both viruses must enter the same cell AND the two halves must recombine intracellularly. Each failure point compounds.

## Model

Viral particle uptake follows gamma-Poisson (negative binomial) distribution, accounting for extreme cell-to-cell heterogeneity in AAV uptake. Calibrated against Omichi et al. (2020): single AAV2 = 83.9% OHC transduction, dual AAV2 = 65.6%.

## Results

| Titer (GC/mL) | Single-vector | Dual-vector (R=50%) | Gap |
|---|---|---|---|
| 10¹⁰ | 1.3% | 0.0% | 1.3% |
| 3×10¹⁰ | 3.8% | 0.0% | 3.8% |
| 10¹¹ | 11.5% | 0.4% | 11.1% |
| 3×10¹¹ | 27.2% | 2.5% | 24.7% |
| 3.75×10¹² (clinical) | ~67% | ~24% | 2.8x advantage |

Gap widens at clinically realistic titers. Only at extreme titers (>10¹³, difficult to manufacture) does dual-vector become comparable.

## Conclusion

Mini-STRC single-vector approach is not incremental — it's a 2.8x improvement in therapeutic coverage at clinical titers. This is the quantitative argument for prioritizing mini-STRC over dual-vector full-STRC.

## Recent Papers

### 2026-04-16 — [Gene Therapy for Infants and Children With Otoferlin-Related Auditory Neuropathy Spectrum Disorder](https://pubmed.ncbi.nlm.nih.gov/41979424/)
*Ear and Hearing · relevance: 7/10 · fda_path*

OTOF (~6kb CDS) uses dual-AAV to deliver the full gene, making it the most prominent real-world test of dual-vector cochlear gene therapy. Clinical efficacy data from OTOF trials will provide the first human evidence on dual-vector recombination efficiency in hair cells — a key parameter in our transduction model. If OTOF dual-vector achieves high efficacy, it may narrow our modeled 2.8x single-vector advantage; if efficacy is lower than expected, it strengthens the case for mini-STRC.

## Connections

- `[part-of]` [[STRC]]
- `[synthesizes]` [[STRC Dual-Vector vs Single-Vector Transduction]]
- `[see-also]` [[STRC Hypothesis Ranking]]
