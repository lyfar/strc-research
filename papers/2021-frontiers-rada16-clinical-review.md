---
title: "Clinical Use of the Self-Assembling Peptide RADA16: A Review of Current and Future Trends in Biomedicine"
type: source
authors: ["Various"]
journal: "Frontiers in Bioengineering and Biotechnology"
date: 2021-06-01
pubmed_id: "34164387"
doi: "10.3389/fbioe.2021.679525"
relevance_type: direct
relevance_score: 8/10
tags: [strc, paper, peptide, rada16, sap, clinical, in-vivo, degradation, pk, puramatrix, h09]
hypotheses: [h09]
status: read-partial
date_added: 2026-04-23
---

## What they found

Comprehensive clinical review of RADA16 (PuraMatrix, PuraStat, PuraSinus) scaffold in vivo performance across tissue types. Best available PK/degradation data: bulk gel dissolution observed within 3–7 days in rat liver injury model, complete by 2 weeks; rat brain TBI model shows scaffold present at 1 week, absent at 3 weeks. Modified/crosslinked variants: 61% degraded at 35 days, 82% non-crosslinked. Degradation products are L-amino acids — no inflammation. Fibril dimensions: individual monomer ~6 nm; nanofibers smaller than wavelength of light; pores 50–200 nm. Modified RADA16 (RADA16-RGD, RADA16-IKVAV etc.): adding ≤12 residues preserves gelation; longer sequences progressively reduce beta-sheet signal and viscoelasticity.

## How this applies to h09 — CRITICAL PK FLAG

This is the primary source for in vivo RADA16 degradation data. The scaffold persists 1–3 weeks in soft CNS/liver tissue. This is INCONSISTENT with the model's K_PROTEOLYSIS = 1.4/h (t½ = 30 min). That rate constant is likely applicable to free monomeric peptide in serum, not the assembled hydrogel scaffold. The assembled gel is sterically protected: literature shows days-to-weeks persistence, implying effective gel t½ on the order of ~1–7 days (k ~0.004–0.029/h for gel, not 1.4/h).

The 118 aa tail (h09 construct) exceeds the 12-residue functional appendage limit established here. Literature consensus: the longer the appended motif, the lower the beta-sheet signal and viscoelasticity. A 118 aa WH2 domain appended to RADA16 is untested — no paper has characterized assembly for such long fusions. This is a critical open question.

## Key numbers — CRITICAL DEGRADATION DATA

- Rat liver: scaffold visible 3–7 days, complete dissolution by 14 days
- Rat brain TBI: scaffold present at 7 days, absent at ~21 days
- Crosslinked RADA16: 61% degraded at day 35; uncrosslinked: 82% degraded at day 35
- Maximum appended sequence preserving assembly: ~12 residues (empirical rule)
- Nanofiber individual monomer: ~6 nm length
- Pore size: 50–200 nm

## Links

- PubMed: https://pubmed.ncbi.nlm.nih.gov/34164387/
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8216384/
- DOI: https://doi.org/10.3389/fbioe.2021.679525

## Connections

- [[STRC Phase 4d F-actin Bundling Model]] — K_PROTEOLYSIS flag; gel t½ days not minutes
- [[STRC Horizontal Top Connector Hydrogel Hypothesis]] — clinical PK; 118 aa tail assembly concern
- `[see-also]` [[Yokoi 2005 RADA16 Reassembly]] — fibril geometry
- `[see-also]` [[Modification Strategies RADA16 PMC9739689]] — tail length limits
