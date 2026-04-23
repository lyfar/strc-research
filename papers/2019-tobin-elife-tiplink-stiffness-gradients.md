---
title: "Stiffness and tension gradients of the hair cell's tip-link complex in the mammalian cochlea"
authors: ["Tobin M", "Chaiyasitdhi A", "Michel V", "Michalski N", "Martin P"]
journal: "eLife"
date: 2019-04-01
pubmed_id: "30932811"
pmc: "PMC6464607"
doi: "10.7554/eLife.43473"
relevance_type: model_calibration
relevance_score: 7/10
tags: [stereocilia, bundle-mechanics, tip-link, stiffness, gating-spring, tonotopy, rat, ohc, ihc]
---

## What they found

Tonotopic gradients of hair-bundle mechanics in rat cochlea (apical region, 1–4 kHz, P7–P10).
Measured total bundle stiffness K_HB, gating-spring contribution K_GS, and pivot contribution K_SP.
This is a tip-link/gating-spring paper — it does NOT measure horizontal top connector (HTC) stiffness.

## Key numbers (directly from paper)

**OHC total bundle stiffness K_HB (rat, apical half, P7–P10, fluid-jet stimulation):**

| CF location | K_HB |
|---|---|
| 1 kHz (apex) | 2.5 ± 0.2 mN/m (= 2.5 pN/nm), n=19 |
| ~2 kHz | ~4–5 mN/m interpolated |
| 4 kHz | 8.6 ± 0.5 mN/m (= 8.6 pN/nm), n=21 |

**IHC total bundle stiffness K_HB (same preparation):**

| CF location | K_HB |
|---|---|
| 1 kHz | 1.7 ± 0.2 mN/m, n=19 |
| 4 kHz | 3.8 ± 0.4 mN/m, n=19 |
| 15 kHz | 5.5 ± 0.4 mN/m, n=14 |

**Gating-spring stiffness K_GS (single spring per tip-link complex):**
- OHC at 1 kHz: 1.3 ± 0.4 mN/m
- OHC at 4 kHz: 3.7 ± 0.7 mN/m
- Bundle stiffness increased ~240% over 2 octaves (1–4 kHz) for OHC

**Tip-link tension at rest:**
- 1 kHz: ~5 pN
- 4 kHz OHC: 34 ± 8 pN

**Species/age:** Sprague-Dawley rat, P7–P10 (pre-hearing, juvenile)
**Method:** calibrated fluid-jet stimulation

## Critical correction for model

**The model cites `K_HTC_PN_PER_NM = 7.5 (approx, Tobin 2019)` — this citation is WRONG.**

Tobin 2019 measures tip-link complex stiffness (K_GS + K_SP). It contains no data on horizontal top connectors (HTCs). The paper's formula is `K_HB = K_GS + K_SP` — no HTC term appears at all. HTC stiffness is measured in Dulon 2019 (Science Advances) and modeled in Kozlov 2011 (Nature).

**Correct HTC stiffness sources:**
- Dulon 2019 net contribution: ~3.07 pN/nm per bundle, ~0.16 pN/nm per link (mouse apical, P13–P15)
- Kozlov 2011 (Nature, bullfrog sacculus model): 20 mN/m = 20 pN/nm aggregate for whole bundle

**The model's `WT_BUNDLE_STIFFNESS_PN_PER_NM = 1400` has no basis in this paper.** Tobin 2019 reports 2.5–8.6 pN/nm for juvenile rat OHC. Even adult basal-turn OHC would not reach 1400 pN/nm (literature max is ~50–100 mN/m = 50–100 pN/nm for highest-frequency locations). The 1400 value is also dead code — never referenced in gate3_bundle_stiffness().

## Links

- PubMed: https://pubmed.ncbi.nlm.nih.gov/30932811/
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6464607/
- DOI: https://doi.org/10.7554/eLife.43473
- PDF: ~/BookLibrary/incoming/tobin-2019-elife-43473-tiplink-stiffness.pdf

## Connections

- `[corrects]` [[STRC Stereocilia Bundle Mechanics Model]] — K_HTC citation wrong; this paper measures tip-link/gating-spring stiffness, not HTC stiffness
- `[informs]` [[STRC Stereocilia Bundle Mechanics Model]] — real OHC bundle stiffness range: 2.5–8.6 pN/nm (juvenile rat 1–4 kHz), NOT 1400 pN/nm
