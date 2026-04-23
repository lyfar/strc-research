---
date: 2026-04-23
type: source
tags: [strc, h05, paper, calcineurin, calcium, calmodulin, kinetics]
hypotheses: [h05]
title: "Dual calcium ion regulation of calcineurin by calmodulin and calcineurin B"
authors: ["Stemmer PM", "Klee CB"]
journal: Biochemistry
year: 1994
volume: 33
pages: "6859-6866"
pubmed_id: "8204620"
doi: "10.1021/bi00188a015"
---

## What they found

Characterization of calcineurin's dual Ca²⁺-sensing mechanism: calmodulin binding to the CaM-binding domain increases catalytic rate and modulates Ca²⁺ response; calcineurin B (regulatory subunit) sets baseline Ca²⁺ sensitivity. Together they produce highly cooperative Ca²⁺ activation (Hill n = 2.8–3). Half-maximum Ca²⁺ concentration for activation: 0.6–1.3 µM depending on CaM concentration.

## Numbers that matter

- Ca²⁺ for half-max calcineurin activation: **0.6–1.3 µM** (model Kd_CaN = 500 nM is slightly below range — a minor underestimate)
- Hill coefficient: **2.8–3** (model n_CaN_Hill = 4.0 is slightly steeper than measured)
- CaM increases catalytic turnover rate ~10-fold at saturating Ca²⁺
- Calcineurin B shifts Ca²⁺ sensitivity to lower concentrations

## Fit to h05

The RBM24 ODE model uses `Kd_CaN = 500 nM` and `n_CaN_Hill = 4.0`. Stemmer & Klee 1994 gives 600–1300 nM and Hill n = 2.8–3. Both values deviate modestly: Kd is ~20% low (conservative — makes CaN activate earlier), Hill n is ~35% steep (makes CaN more switch-like). Neither is a fatal error for the qualitative CaMKII:CaN decoding hypothesis.

## Connections

- `[informs]` [[STRC Calcium Oscillation Acoustic Therapy]] — CaN kinetics in RBM24 ODE
- `[see-also]` [[1998-dolmetsch-calcium-oscillations-gene-expression]] — downstream of CaN
- `[part-of]` [[calcium-oscillation]] (literature-params topic)
