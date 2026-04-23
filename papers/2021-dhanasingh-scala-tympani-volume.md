---
title: "Cochlear Size Assessment Predicts Scala Tympani Volume and Electrode Insertion Force — Implications in Robotic Assisted Cochlear Implant Surgery"
type: source
authors: ["Dhanasingh A", "Swords C", "Bance M", "Van Rompaey V", "Van de Heyning P"]
journal: "Frontiers in Surgery"
date: 2021-09-27
pubmed_id: "34660703"
pmc: "PMC8514755"
doi: "10.3389/fsurg.2021.723897"
relevance_type: critical
relevance_score: 8/10
tags: [strc, pkpd, inner-ear, cochlea, anatomy, perilymph, scala-tympani, volume]
hypotheses: [h09]
---

## What they found

Measured human scala tympani (ST) volume from 30 micro-CT cadaveric temporal bone datasets using 3D segmentation (Slicer 4.10.2 or Stradwin 6.1, 24–30 μm isotropic voxels). Primary purpose was cochlear implant surgical planning, but yields the best direct human ST volume measurement.

## Numbers that matter

| Parameter | Value | Units | Notes |
|---|---|---|---|
| Human ST volume, mean | 34.2 ± 7 | μL | n=30 cadaveric temporal bones |
| Human ST volume, range | 23–50 | μL | Full range observed |

**Implication for model:** The model uses `PERILYMPH_VOL_UL = 70 μL` for "human cochlea." This is inconsistent with direct anatomical measurement:
- ST alone: 34.2 μL mean
- Total perilymph (ST + SV), from Ekdale 2016 μCT: ~93 μL
- 70 μL is not directly supported by either measurement; it may be a rough midpoint or from an older anatomical estimate not traced in the model

The drug-accessible compartment for RWM-delivered agents is primarily ST (drug enters basally, distributes apically by diffusion). Using 93 μL would underestimate concentration; 34 μL would overestimate. 70 μL has no traceable primary citation.

## Access

PMC open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC8514755/

## Connections

- `[source]` [[hydrogel_phase4e_cochlear_pkpd]] — PERILYMPH_VOL_UL discrepancy; 34 μL (ST only) vs 93 μL (total) vs 70 μL (model, unsourced)
- [[STRC Hypothesis Ranking]] — h09 PK volume parameter
