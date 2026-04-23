---
date: 2026-04-23
type: source
tags: [strc, h05, paper, adenylyl-cyclase, ac1, calmodulin, calcium, kinetics]
hypotheses: [h05]
title: "Distinct mechanisms of calmodulin binding and regulation of adenylyl cyclases 1 and 8"
authors: ["Masada N", "Schaks S", "Jackson SE", "Sinz A", "Cooper DM"]
journal: Biochemistry
year: 2012
pubmed_id: "22971080"
doi: "10.1021/bi300943d"
---

## What they found

Masada et al. used GST pull-down, mass spectrometry, and stopped-flow fluorescence to dissect how calmodulin (CaM) binds and activates adenylyl cyclase 1 (AC1) vs AC8. Key finding: AC1 is activated by an initial encounter complex between CaM and the IQ-like domain in the N-terminus, forming a 1:1 AC1:CaM complex at physiological Ca²⁺ concentrations. Activation is cooperative with Ca²⁺ occupancy of CaM's C-lobe driving the stimulatory interaction. This establishes the mechanistic basis for AC1's Ca²⁺/CaM sensitivity.

## Numbers that matter

- AC1 activation by Ca²⁺/CaM: cooperative, driven by C-lobe Ca²⁺ binding
- Calmodulin binding stoichiometry: 1:1 (AC1:CaM)
- Mechanism: IQ-like domain in AC1 N-terminus is the primary CaM docking site
- AC8 uses a different (IQ-motif independent) mechanism — not applicable to h05
- **K_Ca for AC1 activation not given as a single Kd number** — determined by CaM affinity constants (~nM–µM range depending on Ca²⁺ saturation state)
- Note: the model's `K_CA_AC1_NM = 150 nM` and `AC1_VMAX_NM_S = 2000 nM/s` are not directly tabulated here; these likely derive from earlier biochemical work (Willoughby & Cooper 2007 review synthesis).

## Fit to h05

This is the AC1 kinetics reference for the pivot model. The "Wu 2011" citation in the scripts is a phantom — Masada 2012 is the real paper. However, Masada 2012 gives mechanistic insight and relative rates, not the exact nM/s Vmax values used in the ODE. Those values require supplemental justification from Willoughby & Cooper 2007 or explicit labeling as estimates.

## Connections

- `[informs]` [[STRC Calcium Oscillation Acoustic Therapy]] — AC1 kinetics underpins pivot model
- `[see-also]` [[2007-willoughby-cooper-adenylyl-cyclase-review]] — companion canonical review
- `[part-of]` [[calcium-oscillation]] (literature-params topic)
