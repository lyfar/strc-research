---
title: "Cochlear outer hair cell horizontal top connectors mediate mature stereocilia bundle mechanics"
authors: ["Dulon D", "Papal S", "Castilla-Eder C", "Daulat A", "Petit C", "El-Amraoui A"]
journal: "Science Advances"
date: 2019-02-22
pubmed_id: "30801007"
pmc: "PMC6382404"
doi: "10.1126/sciadv.aat9934"
relevance_type: model_calibration
relevance_score: 10/10
tags: [htc, bundle-mechanics, stiffness, afm, stereocilia, strc, direct-measurement]
---

## What they found

First direct measurement of HTC spring contribution to OHC hair bundle stiffness using
non-contact acoustic FM-AFM. Measured at P9–P15 in Strc-/-, Strc+/-, and Strc-/-/Tecta-/-
mice, isolating the specific contribution of horizontal top connectors.

## Key numbers (directly from paper)

| Condition | OHC bundle stiffness | IHC bundle stiffness |
|---|---|---|
| With HTCs (Strc+/−/Tecta−/−) | **5.12 ± 0.46 pN/nm** | 2.34 ± 0.64 pN/nm |
| Without HTCs (Strc−/−/Tecta−/−) | **2.05 ± 0.15 pN/nm** | 2.74 ± 0.50 pN/nm |

- **Reduction: ~60%** in OHC bundle stiffness without HTCs
- IHC stiffness is UNCHANGED — HTCs are OHC-specific, IHC are unaffected
- Measurements at P13–P15 (post HTC maturation)

## Unit conversion for model

- 5.12 pN/nm = **5.12 mN/m** (per bundle, ~20 stereocilia)
- HTC contribution: 5.12 − 2.05 = **3.07 mN/m** per bundle
- Per stereocilium: ~0.154 mN/m
- Per HTC link (19 links in 20-stereocilia chain): ~0.16 mN/m per link

This is **~5.7× larger** than my original k_HTC estimate (0.028 mN/m).
Model needs recalibration with these values.

## Developmental timeline

Top connector maturation coincides with hearing onset (P12–P15).
Significant increase in OHC (not IHC) stiffness from P12→P15 = HTC maturation window.
No adult measurements — limitation of the study.

## Relevance to model

Provides ground-truth calibration for k_HTC in stereocilia_bundle_mechanics.py.
The 60% reduction figure is confirmed with direct measurement, not inferred.
Critical: IHC stiffness unchanged = tip-link parameters remain valid for HTC-free state.

## Links

- PubMed: https://pubmed.ncbi.nlm.nih.gov/30801007/
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6382404/
- DOI: https://doi.org/10.1126/sciadv.aat9934

## Connections

- `[validates]` [[STRC Stereocilia Bundle Mechanics Model]] — direct k_HTC measurement confirms 60% stiffness loss
