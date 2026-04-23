---
title: "Quantification of solute entry into cochlear perilymph through the round window membrane"
type: source
authors: ["Salt AN", "Ma Y"]
journal: "Hearing Research"
date: 2001-04-01
pubmed_id: "11423219"
pmc: null
doi: "10.1016/S0378-5955(01)00223-4"
relevance_type: critical
relevance_score: 10/10
tags: [strc, pkpd, inner-ear, cochlea, round-window, permeability, pharmacokinetics]
hypotheses: [h09]
---

## What they found

Primary quantitative source for RWM permeability and scala tympani clearance. Used TMPA (trimethylphenylammonium — a low-MW ion) as marker in guinea pig. Applied to intact RWM for 90 min; measured distribution with ion-selective microelectrodes at turns 1 and 2.

After 90 min:
- Turn 1 (1.4 mm from base): 330 ± 147 μM (n=8)
- Turn 2 (7.5 mm from base): 15 ± 33 μM (n=5)

Simulation fitting yielded the canonical parameter set for the WUSTL Cochlear Fluids Simulator default (guinea pig).

## Numbers that matter

| Parameter | Value | Units | Notes |
|---|---|---|---|
| RWM permeability (TMPA, low MW) | 1.9 × 10⁻⁸ | cm/s | Guinea pig; low-MW ion |
| ST clearance half-time | 60 | min | = 0.693/h → K_clear ≈ 0.69/h |
| Longitudinal perilymph flow | 4.4 | nL/min (base→apex) | Guinea pig |

**Critical limitation:** TMPA is a small ion (MW ~166 Da). The paper does NOT provide permeability for 14 kDa peptides. The model's K_RWM = 0.02/h citing "Salt 2011" cannot be traced to this paper — and the 2001 value would be for a very different (small) molecule.

## What this paper does NOT provide

- No MW scaling table for RWM permeability
- No human-specific parameters (all guinea pig)
- No middle ear clearance rate
- Guinea pig ST clearance t½ = 60 min ≠ human

## Access status

**PAYWALLED** — ScienceDirect, Hearing Research Vol.154(1-2):88-97. Not on PMC. Needs institutional access or Sci-Hub.

## Connections

- `[primary source]` [[hydrogel_phase4e_cochlear_pkpd]] — K_PERILYMPH_CLEAR baseline (guinea pig 60 min t½); K_RWM is not from this paper
- `[cited by]` [[2018-salt-plontke-pharmacokinetic-principles-inner-ear]]
