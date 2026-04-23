---
title: "Perilymph Kinetics of FITC-Dextran Reveals Homeostasis Dominated by the Cochlear Aqueduct and Cerebrospinal Fluid"
type: source
authors: ["Salt AN", "Hartsock JJ", "Plontke SK", "DeMott JE", "Riveros A"]
journal: "Journal of the Association for Research in Otolaryngology"
date: 2015-08-01
pubmed_id: "26001373"
pmc: "PMC4417094"
doi: "10.1007/s10162-015-0512-1"
relevance_type: important
relevance_score: 8/10
tags: [strc, pkpd, inner-ear, cochlea, perilymph, clearance, cochlear-aqueduct]
hypotheses: [h09]
---

## What they found

Used FITC-dextran (large MW tracer — slower diffusion, less tissue uptake) to characterize perilymph homeostasis mechanisms. Found perilymph clearance is dominated by CSF exchange via cochlear aqueduct. All studies in guinea pig.

## Numbers that matter

| Parameter | Value | Units | Notes |
|---|---|---|---|
| CSF entry rate via cochlear aqueduct | ~30 | nL/min | Guinea pig sealed cochlea |
| ST→tissue communication half-time | 6 | min | Fast exchange with spiral ligament etc. |
| SV/vestibule/SCCs elimination half-time | 230 | min | Slow; dextran "lost more slowly than other substances" |
| CSF-driven dilution of drug | 53% | of gentamicin entering inner ear lost to CSF | Modeled outcome |

## Implications for model

The K_PERILYMPH_CLEAR = 0.35/h (t½ = 2h) in the model is in the right ballpark for slow-clearing large molecules (230 min = 3.8h → 0.18/h) but faster than what this paper shows for high-MW tracers. For a 14 kDa peptide (smaller than FITC-dextran but larger than steroids), the actual clearance would be between these extremes.

CSF inflow rate of 30 nL/min corresponds to ~1.8 μL/h. In 70 μL perilymph, this gives a dilution half-life of ~27 h from CSF turnover alone — suggesting CSF exchange is NOT the dominant clearance route for most drugs within a 48h window; tissue uptake and aqueduct bulk flow dominate.

## Access

PMC open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC4417094/

## Connections

- `[primary source]` [[hydrogel_phase4e_cochlear_pkpd]] — K_PERILYMPH_CLEAR; 2h model value is between guinea pig fast-MW (60 min) and dextran-slow (230 min) — plausible but not directly measured for 14 kDa peptide
- `[informs]` [[2001-salt-ma-quantification-rwm-permeability]] — context for ST clearance parameter
