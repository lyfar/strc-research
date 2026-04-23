---
date: 2026-04-23
type: source
tags: [strc, h05, paper, met-channel, mechanotransduction, hair-cell, review]
hypotheses: [h05]
title: "The physiology of mechanoelectrical transduction channels in hearing"
authors: ["Fettiplace R", "Kim KX"]
journal: "Physiological Reviews"
year: 2014
volume: 94
pages: "951-986"
pubmed_id: "24987009"
doi: "10.1152/physrev.00038.2013"
---

## What they found

Comprehensive review of MET channel biology: biophysics, gating, Ca²⁺ permeation, adaptation, molecular identity (TMC1/TMC2), and tonotopy. Synthesizes data on channel numbers per stereocilium (~2), bundle conductance (~1–6 nA peak), single-channel conductance gradients (apical ~100 pS to basal ~300 pS in mouse), and the role of Ca²⁺ in fast and slow adaptation.

## Numbers that matter

- MET channels per stereocilium: ~1–2 (i.e., ~50–100 per OHC for a bundle of 50–70 stereocilia)
- **n_channels = 134 in the h05 RBM24 model cites "Fettiplace 2017" — the correct citation is Fettiplace & Kim 2014 (this paper). The year "2017" is wrong.** 134 channels is plausible for a mid-frequency OHC (67 stereocilia × 2 channels each).
- Ca²⁺ permeation fraction f_Ca: ~15% at physiological ionic conditions (consistent with Beurg 2006)
- Single-channel conductance: 145–210 pS apical-to-basal in rat (matches Beurg 2006)

## Fit to h05

This review is the correct source for `n_channels = 134` (replacing the phantom "Fettiplace 2017"). The value itself is plausible per the review's channel density data. Note: the review itself notes variability in channel number across cochlear position and species; 134 should be treated as an order-of-magnitude estimate, not a precise measurement.

## Connections

- `[informs]` [[STRC Calcium Oscillation Acoustic Therapy]] — MET channel count (n_channels)
- `[see-also]` [[2006-beurg-met-channel-conductance]] — g_MET value
- `[part-of]` [[calcium-oscillation]] (literature-params topic)
