---
date: 2026-04-23
type: source
tags: [strc, h05, paper, adenylyl-cyclase, ac1, calcium, cAMP, review]
hypotheses: [h05]
title: "Organization and Ca2+ regulation of adenylyl cyclases in cAMP microdomains"
authors: ["Willoughby D", "Cooper DM"]
journal: "Physiological Reviews"
year: 2007
volume: 87
pages: "965-1010"
pubmed_id: "17615394"
doi: "10.1152/physrev.00049.2006"
---

## What they found

Comprehensive review of all mammalian adenylyl cyclase isoforms, with emphasis on Ca²⁺-regulated isoforms (AC1, AC3, AC8 as Ca²⁺/CaM-stimulated; AC5, AC6 as Ca²⁺-inhibited). For AC1 specifically: Ca²⁺/CaM stimulation is cooperative, EC50 for Ca²⁺ is in the 100–500 nM range (depending on assay conditions), and Vmax in reconstituted systems is in the µM/min range. Review synthesizes data from multiple labs on cAMP microdomain signaling and the spatial organization of AC isoforms with PKA scaffolds (AKAPs).

## Numbers that matter

- AC1 EC50 for Ca²⁺ (via CaM): ~100–500 nM (assay-dependent; K_Ca = 150 nM used in h05 model is within this range)
- AC1 Vmax: ~2–5 µM cAMP/min in membrane preparations; Vmax = 2000 nM/s ≈ 120 µM/min is likely an overestimate unless at saturation in a concentrated cell compartment
- AC1 Hill coefficient for Ca²⁺ activation: 1.5–2.5 (n=2 in model is consistent)
- AC1 mRNA/protein expression: highest in neurons, cochlea, olfactory bulb; expressed in cochlear IHC and OHC (cited separately as Visel 1997/Vorobiova 1997 series)
- **This is the canonical reference that should replace the phantom "Wu 2011" citation**

## Fit to h05

Willoughby & Cooper 2007 is the correct source for the AC1 kinetic rationale in the h05 pivot model. The model's K_CA_AC1_NM = 150 nM falls within the reported EC50 range; AC1_VMAX_NM_S = 2000 nM/s needs to be checked against cell-level measurements (membrane prep values converted to cytoplasmic concentrations depend on AC1 copy number, which is not specified in the model).

## Connections

- `[informs]` [[STRC Calcium Oscillation Acoustic Therapy]] — replaces phantom "Wu 2011"
- `[see-also]` [[2012-masada-ac1-ac8-calmodulin-kinetics]] — mechanistic follow-up
- `[part-of]` [[calcium-oscillation]] (literature-params topic)
