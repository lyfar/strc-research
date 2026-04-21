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

## Recent Papers

### 2026-04-21 — [LinkCraft: An interactive tool for the design of flexible linkers](https://pubmed.ncbi.nlm.nih.gov/42001978/)
*Journal of Molecular Biology · relevance: 4/10 · lateral*

Computational tool for rational design of intrinsically disordered linkers in multi-domain proteins, with ensemble-based structural modeling. Directly applicable to optimizing the truncation junction (aa 699↔700) in Mini-STRC and the Ultra-Mini construct. LinkCraft could model whether the junction region behaves as an unstructured linker and suggest sequence modifications that maintain proper inter-domain spacing and dynamics while minimizing misfolding risk at the truncation boundary.

### 2026-04-16 — [Gene Therapy for Infants and Children With Otoferlin-Related Auditory Neuropathy Spectrum Disorder](https://pubmed.ncbi.nlm.nih.gov/41979424/)
*Ear and Hearing · relevance: 7/10 · fda_path*

OTOF gene therapy is the first successful inner ear gene therapy and establishes the regulatory and clinical precedent for all subsequent programs including mini-STRC. The dual-AAV approach used for otoferlin (~6kb CDS) directly validates that oversized hearing loss genes can be delivered to the cochlea, though our mini-gene strategy offers a single-vector alternative. Clinical endpoints and pediatric patient selection criteria from OTOF trials will likely shape FDA expectations for DFNB16.

### 2026-04-16 — [Postnatal Slc26a4 gene therapy improves hearing and structural integrity in a hereditary hearing loss model](https://pubmed.ncbi.nlm.nih.gov/41701544/)
*The Journal of Clinical Investigation · relevance: 6/10 · applicable*

Demonstrates postnatal AAV gene therapy can restore hearing in another hereditary HL model (DFNB4), validating the general approach. The identification of a critical therapeutic window is particularly important — we need to determine if DFNB16 has a similar window before irreversible OHC damage occurs. Their use of single-AAV for the ~2.3kb SLC26A4 CDS provides a regulatory reference for single-vector cochlear gene therapy.

### 2026-04-20 — [Generation of human iPSC lines from a patient with OTOF-related deafness](https://pubmed.ncbi.nlm.nih.gov/42000544/)
*Stem Cell Res · relevance: 5/10 · applicable*

iPSC disease modeling platform for DFNB9 (OTOF compound heterozygous mutations). The carrier line concept — comparing homozygous patient vs. heterozygous carrier — is directly transferable to STRC/DFNB16 for testing mini-STRC dose-response in human-derived cells. Could generate equivalent lines for STRC mutations.

### 2026-04-20 — [Localized active transport shapes nanoscopic features at mechanosensory cilia tips](https://pubmed.ncbi.nlm.nih.gov/41941305/)
*J Cell Biol · relevance: 7/10 · lateral*

Kif19A coordinates molecular enrichment at mechanosensory cilia tips via localized active transport counteracting diffusion. Critical implication for mini-STRC: truncation reduces protein size, increasing diffusion rate, which could impair tip concentration. The transport-binding-diffusion balance model should be applied to predict whether mini-stereocilin achieves sufficient tip localization.

### 2026-04-20 — [Seed-based design + ribosome display for nanobody-like scaffolds](https://pubmed.ncbi.nlm.nih.gov/41995290/)
*Biosci Biotechnol Biochem · relevance: 4/10 · lateral*

Conserved "seed framework" residues (just 10 residues) can nucleate correct protein folding in randomized contexts. Inverse of our truncation approach: instead of removing internal domains, they kept minimal seeds. MSA-based identification of conserved framework residues in stereocilin orthologs could define the irreducible core that must be preserved, potentially enabling a smaller, more stable mini-STRC.

### 2026-04-20 — [Tri-modal contrastive learning for protein representation](https://pubmed.ncbi.nlm.nih.gov/41990738/)
*Cell Rep Methods · relevance: 3/10 · lateral*

ProteinAligner integrates sequence + structure + literature for protein function prediction, outperforming sequence-only models. Stereocilin is poorly characterized experimentally — a multimodal model incorporating literature context could predict functional impact of truncations more accurately than AlphaFold alone.

### 2026-04-20 — [ProteinMPNN engineering for hyper-stable protein redesign](https://pubmed.ncbi.nlm.nih.gov/41988819/)
*Microb Biotechnol · relevance: 4/10 · lateral*

ProteinMPNN redesigned non-functional scaffold regions for enhanced stability while preserving functional epitopes, achieving soluble expression without inclusion bodies. Directly applicable: ProteinMPNN could optimize the truncation junction (aa 699↔1776) in mini-STRC, redesigning exposed surfaces for improved thermodynamic stability without altering functional domains.

### 2026-04-20 — [Mechanical confinement regulates ECM remodeling via MAPK and Hedgehog](https://pubmed.ncbi.nlm.nih.gov/41985400/)
*Biomaterials · relevance: 3/10 · lateral*

Mechanical confinement level directly regulates ECM production through cilia-dependent MAPK/Hedgehog signaling. OHCs exist in a mechanically confined space defined partly by stereocilin-mediated TM attachment. Loss of stereocilin may trigger aberrant confinement signaling beyond simple mechanical decoupling — mini-STRC must restore sufficient TM attachment to normalize this.

### 2026-04-20 — [Cilia biophysics in PCD: quantitative characterization of RSPH4A variant](https://pubmed.ncbi.nlm.nih.gov/41972697/)
*Cells · relevance: 3/10 · lateral*

High-speed video microscopy quantified biophysical consequences of a specific genetic variant on ciliary mechanics (angular excursion, beat amplitude). Methodology adaptable for stereocilia bundle mechanics: measuring deflection amplitude, stiffness, and recovery dynamics to functionally compare mini-STRC vs. WT bundles beyond morphology alone.

### 2026-04-20 — [Shear stress mechanosensing in endothelial cells](https://pubmed.ncbi.nlm.nih.gov/41969033/)
*J Cell Sci · relevance: 4/10 · lateral*

Endothelial cells use distributed mechanosensor networks (glycocalyx, ion channels, cilia) rather than single receptors. Stereocilin may function like the glycocalyx — not transducing force directly but shaping frequency-dependent mechanical coupling to the MET apparatus. Mini-STRC must preserve this filtering function, not just physical attachment.
