---
date: 2026-04-24
type: synthesis
title: Mini-STRC Single-Vector Hypothesis
tags:
- aav
- gene-therapy
- truncation
- mini-gene
- strc
- synthesis
- hypothesis
- legacy-hypothesis-sheet
status: active
priority: primary
tier: S
lit_audit: fixed
lit_audit_date: 2026-04-24
stage: clinical-vector-ready
---

> Canonical hypothesis state now lives in `~/STRC/hypotheses/h03-mini-strc/index.md` (post vault-split 2026-04-24). Full research wiki at [wiki.strc.egor.lol](https://wiki.strc.egor.lol/hypotheses/h03-mini-strc). This page is the legacy synthesis sheet kept in sync for marketing.

## Core claim

STRC (5,325 bp CDS) is too big for single AAV (4,700 bp limit). The **Ultra-Mini construct (aa 1075–1775, 2,106 bp CDS)** preserves the TMEM145-binding C-terminal core and fits into a single AAV with comfortable regulatory headroom. The original Mini construct (aa 700–1775, 3,228 bp CDS) retains more of the protein but leaves much less room for regulatory elements. Ultra-Mini is the current lead architecture; Mini 700–1775 stays as backup.

## Lead construct — Ultra-Mini 1075–1775

- CDS 2,106 bp (701 aa), pTM 0.87, 2% disordered
- TMEM145 GOLD-domain binding preserved: AF3 ipTM 0.68 (vs 0.47 full-length)
- GPI-anchor signal intact: NetGPI 1.1, omega site S1749
- Signal peptide: IgK leader scores 99.97% on SignalP 6.0
- Localization: DeepLoc 2.1 → Extracellular + Lipid anchor (72.1%)
- CpG depletion: 0 CpG achieved at 3.65% CAI cost (`cpg_depletion_ultra_mini_strc.py`)
- 5 of 14 N-glycosylation sites retained, all high-confidence (NetNGlyc 1.0: 0.52–0.72)
- ΔΔG of binding destabilization 8.4 kcal/mol; folding stable +0.9 kcal/mol

## Vector architecture (B8 + Ultra-Mini + WPRE3)

Full assembled cassette, 5'-ITR to 3'-ITR:

| Element | Size (bp) | Source |
|---|---|---|
| 5' ITR (AAV2) | 145 | Samulski 1987 |
| B8 OHC-specific enhancer | **706** | Zhao et al. 2025, *Neuron* 113(10):1579–1596 · PMID 40262614 |
| Kozak + IgK signal peptide | 63 | canonical / Choi 2014 |
| Ultra-Mini CDS (CpG-depleted) | 2,106 | this project |
| Stop codon | 3 | — |
| WPRE3-compact | 219 | Choi et al. 2014, *Mol Brain* 7:17 · PMID 24618276 |
| bGH polyA | 208 | pAAV-MCS canonical |
| 3' ITR (AAV2) | 145 | Samulski 1987 |
| **Total** | **3,565** | 75.9% of 4,700 bp AAV capacity — 1,135 bp spare |

B8 composition (Zhao 2025 Fig 6A / Fig S6A-B): E1P3×2 + E2P2×2 + E2P3×2, with modules E1P3 = 93 bp, E2P2 = 132 bp, E2P3 = 128 bp — contiguous back-calculation 706 bp. Exact cloned sequence including KpnI/XbaI linkers lives in Zhao 2025 Table S2 (SI Excel); treated as a cloning-QC gate before GMP submission.

## Literature audit — closed 2026-04-24

A vault-wide audit on 2026-04-23/25 flagged phantom citations across most hypotheses. For h03 the fixes were:

- B8 enhancer source: "Yoshimura 2018" (phantom) → **Zhao et al. 2025 Neuron, PMID 40262614**
- Myo15 truncated promoter: "Zhao 2024" (phantom) → **Hu et al. 2024 Research (AAAS), PMID 38665848**
- WPRE3: "Choi 2014 Cell 157" (wrong journal) → **Choi et al. 2014 Mol Brain 7:17, PMID 24618276**
- B8 size: script value 587 bp (no literature backing) retracted; set to 706 bp per Zhao 2025 module back-calc

Vector scripts re-run clean with corrected parameters. `lit_audit: partial → fixed`; h03 re-promoted A → S.

## Backup construct — Mini 700–1775

- CDS 3,228 bp (1,076 aa), pTM 0.86, 4% disordered — retains LRR domain
- N-terminal 1–699 is intrinsically disordered (pTM 0.27, 38% unstructured); removing it improves folding 0.63 → 0.86 (hierarchical constraint propagation, Li et al. 2026)
- Comparative genomics: birds lack STRC yet hear; Pfam LRR boundary coincides with this cut
- 87% of pathogenic variants land in the retained region (ClinVar, 427 variants)
- Fits AAV only with minimal regulatory elements — no room for WPRE3 + B8 together
- Kept as fallback if Ultra-Mini fails functional validation

## Open questions (lab required)

- Does Ultra-Mini localize correctly to stereocilia tips?
- Does it form horizontal top connectors with the tectorial membrane?
- Does hair cell function recover (ABR, DPOAE in STRC-knockout mice)?
- Does loss of 9 N-glycosylation sites affect trafficking in vivo?

## Next steps

1. Order Ultra-Mini gBlock + B8+WPRE3 cassette
2. Verify B8 exact cloned sequence against Zhao 2025 Table S2 (KpnI/XbaI linkers)
3. Clone pAAV B8-IgK-Ultra-Mini-WPRE3-bGH
4. Phase 4 HEK coIP — does Ultra-Mini pull down TMEM145?
5. Shanghai Shu Yilai c.4976 knock-in mouse — transduction + ABR rescue

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

## Connections

- `[part-of]` [[STRC]]
- [[STRC Mini-STRC Single-Vector Hypothesis]]
- `[see-also]` [[STRC Hypothesis Ranking]]
