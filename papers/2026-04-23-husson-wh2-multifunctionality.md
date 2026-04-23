---
title: "Multifunctionality of the β-thymosin/WH2 module: G-actin sequestration, actin filament growth, nucleation, and severing"
authors: ["Husson C", "Cantrelle FX", "Roblin P", "Didry D", "Le KH", "Perez J", "Guittet E", "Van Heijenoort C", "Renault L", "Carlier MF"]
journal: "Annals of the New York Academy of Sciences"
year: 2010
volume: 1194
pages: "44–52"
doi: "10.1111/j.1749-6632.2010.05473.x"
date: 2026-04-23
type: source
tags: [strc, paper, actin, wh2, beta-thymosin, kinetics, h09]
hypotheses: [h09]
---

## What they found

Comprehensive review/research article from the Carlier lab (MF Carlier is the authoritative actin
dynamics biochemist) on the β-thymosin/WH2 module. Used chimeric protein engineering, actin
polymerization assays, ITC, NMR, and SAXS to dissect how a single βT/WH2 domain can function
as G-actin sequesterer, filament barbed-end deliverer, nucleator, or severer depending on
C-terminal sequence context.

## Numbers that matter

**Thymosin-β4 (Tβ4) × G-actin:**
- Kd = 1 μM for ATP-actin (moderate affinity, the reference baseline)
- Kd = 80–100 μM for ADP-actin (50–100× weaker — nucleotide-state dependent)

**Tβ4 × F-actin (side-binding):**
- Kd = 5–10 mM (weak cooperative binding at [Tβ4] > 20 μM)
- This is 5,000–10,000× weaker than G-actin binding
- Mechanism: C-terminal α-helix of Tβ4 interferes with filament contacts

**Implication for WH2 × F-actin:**
- WH2 lacks the C-terminal α-helix that interferes with F-actin in Tβ4
- But WH2 canonical binding site (barbed-end groove) is also involved in longitudinal
  actin–actin contacts in filament — still not a clean side-binding site
- No direct measurement of WH2 × F-actin Kd reported in this paper
- The paper does not support the existence of a low-μM WH2 × F-actin interaction

**Context for model parameter WH2_KD_FACTIN_M = 5 μM:**
- Based on all available literature, 5 μM would be VERY optimistic for a WH2 × F-actin
  side-binding event. It is possible if constructs are specifically engineered, but not for
  a canonical isolated WH2 peptide. True value is likely ≥ 1 mM or unmeasurable.

## Connections

- `[primary source for]` [[STRC H09 WH2 F-actin Bundling Hypothesis]] — Tβ4/WH2 Kd values
  and F-actin binding data; shows F-actin side-binding is mM-range or absent
- `[informs]` [[Hydrogel Phase 4d F-actin Bundling Model]] — key constraint on
  WH2_KD_FACTIN_M: the 5 μM value in model is 1000× more optimistic than Tβ4 literature
- `[kinetics reference]` [[Actin Treadmilling Stereocilia]] — Carlier lab is primary
  authority on actin dynamics constants
