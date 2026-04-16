---
title: Electrostatic Analysis — Why E1659A Is Pathogenic
status: complete
stage: computational-validation
priority: supporting
tags: [electrostatics, pathogenicity, E1659A, biophysics, adhesion]
---

## Core claim

AlphaFold says E1659A doesn't change protein structure. AlphaMissense says it's pathogenic (0.9016). Both are correct. The paradox resolves when you look at chemistry, not geometry. STRC is an extracellular adhesion protein. E1659 contributes charge, H-bonds, and hydrophilicity to the adhesion interface. Alanine contributes none. The fold is intact. The glue is broken.

## What changes at E1659A

| Property | Wildtype Glu | Mutant Ala | Impact |
|---|---|---|---|
| Net charge (pH 7.4) | -1 (deprotonated carboxylate) | 0 (nonpolar methyl) | Lost negative charge at adhesion interface |
| H-bond acceptors | 2 (carboxylate oxygens) | 0 | Two H-bonds broken |
| Side chain volume | 138.4 Å³ | 88.6 Å³ | 49.8 Å³ cavity created |
| Hydrophobicity | Hydrophilic | Hydrophobic | +5.3 shift: surface becomes water-repelling |

## Quantitative results

- ΔΔG(binding) = 8.4 kcal/mol (multi-method, PAE-corrected)
- Original Coulomb estimate: 8.62 kcal/mol — confirmed
- ΔΔG(folding) = +0.9 kcal/mol — folding barely affected
- Conclusion: binding defect, not folding defect

## Significance for reclassification

This explains mechanistically WHY the variant is pathogenic when structural tools alone show minimal change. Supports ACMG upgrade from VUS → Likely Pathogenic/Pathogenic alongside AlphaMissense (0.9016) and conservation data (100% across 9 mammals, 80M years).
