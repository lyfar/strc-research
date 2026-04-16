---
title: Dual-Vector vs Single-Vector Transduction Analysis
status: active
stage: computational-modeling
priority: supporting
tags: [aav, dual-vector, single-vector, transduction, statistics, gamma-poisson]
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
