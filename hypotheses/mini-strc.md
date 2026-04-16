---
title: Mini-STRC Single-Vector Hypothesis
status: active
stage: computational-validation
priority: primary
tags: [aav, gene-therapy, truncation, mini-gene]
---

## Core claim

STRC (5,325 bp CDS) is too big for single AAV (4,700 bp limit). Truncating to residues 700–1775 produces a 3,228 bp mini-gene that fits in single AAV with 1,472 bp headroom, retains both LRR domain and C-terminal functional core, and maintains structural integrity.

## Current evidence (all computational)

- pTM 0.86 (700–1775 construct), 4% disordered — best fold-to-coverage ratio
- N-terminal (1–699) is intrinsically disordered: pTM 0.27, 38% unstructured
- Removing N-terminal improves folding 0.63 → 0.86 (explained by hierarchical constraint propagation, Li et al. 2026)
- TMEM145 GOLD domain binding preserved: ipTM 0.43 vs 0.47 full-length
- GPI-anchor signal intact: NetGPI 1.1 confirms omega site S1749
- Signal peptide: IgK SP scores 99.97% on SignalP 6.0 (vs native 93.4%)
- Localization confirmed: DeepLoc 2.1 → Extracellular + Lipid anchor (72.1%)
- Full construct fits AAV: SP + B8 enhancer + CDS + polyA = 4,103 bp (597 bp headroom)
- 5 of 14 N-glycosylation sites retained, all high-confidence (NetNGlyc 1.0: 0.52–0.72)
- Normal Mode Analysis: N-terminal dynamically decoupled (modes 1-20 localized >92%), subspace overlap 0.905
- Population variants: 87% of pathogenic variants in retained region (ClinVar, 427 variants)
- Comparative genomics: birds lack STRC yet hear; Pfam LRR boundary = our cut point
- CpG hazard: 156 CpGs (4.9x genome avg) — 87% reducible by synonymous swaps at 3% CAI cost
- ΔΔG confirmed: binding destabilization 8.4 kcal/mol, folding stable +0.9 kcal/mol
- Precedent: micro-dystrophin (FDA approved) validates mini-gene approach

## Backup construct

Residues 1075–1775: pTM 0.87, 2,597 bp headroom. Better fold but removes 60% of protein — needs functional validation before using.

## Open questions (lab required)

- Does mini-STRC localize correctly to stereocilia tips?
- Does it form horizontal top connectors with tectorial membrane?
- Does hair cell function recover (ABR, DPOAE in STRC knockout mice)?
- Does loss of 9 N-glycosylation sites affect trafficking in vivo?

## Computational next steps

- Molecular dynamics (GROMACS/OpenMM) at 37°C to verify fold stability
- Commercial codon optimization — current CAI 0.711, GC 54.5% (can improve 2-10x expression)
- AF3 with membrane context and glycan modifications
- Systematic exon-deletion sweep (29 AF3 jobs) to find structurally dispensable exons
- HADDOCK/ClusPro docking as orthogonal validation of AF3 interaction predictions
- CpG depletion — **must do before lab**
