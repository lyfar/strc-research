---
title: Delivery — Sonoporation + LNP Through Round Window Membrane
status: active
stage: computational-modeling
priority: secondary
tags: [delivery, sonoporation, LNP, round-window, non-invasive, ODE]
---

## Core claim

Standard cochlear gene therapy requires surgery (round window membrane injection under general anesthesia — risky for 4-year-olds, one-shot only due to AAV immune response). Alternative: ultrasound + microbubbles create transient pores in RWM, LNP-packaged mRNA passes through non-invasively. Repeatable.

## Why current delivery has problems

- General anesthesia required (risk for young children)
- One-shot only: anti-AAV antibodies form after first dose
- Surgical risk: cochlear damage, further hearing loss
- Dual-vector: full STRC needs 2 viruses entering same cell (low probability)

## Proposed mechanism

1. Gel with lipid microbubbles (SonoVue, 2–5 µm) + LNP-packaged mRNA placed in ear canal against RWM
2. 1 MHz ultrasound (3 W/cm², MI 0.254) → microbubbles cavitate → transient pores ~110 nm in RWM
3. LNPs (80 nm) pass through pores via hindered diffusion (Renkin equation, ~7.5x reduction vs free diffusion)
4. OHCs endocytose LNPs → ionizable lipids destabilize endosome → mRNA released into cytoplasm
5. RWM recovers fully within 24 hours (Shih et al. 2019: zero ABR shift)

## ODE model results

- Scenario 1 (standard LNPs): **0.39% transduction — FAILED.** 258 sessions needed. Bottleneck: 2% endosomal escape.
- Scenario 2 (optimized LNPs): **78% transduction in 2 sessions.** Changes: 10x LNP concentration, ionizable lipids (SM-102 class, 30% escape), optimized timing.
- Scenario 3 (hybrid): Year 0 standard AAV surgery (Anc80L65, 80% OHC transduction). Year 5+ if expression drops: single optimized LNP sonoporation session as top-up.

## Sensitivity analysis

Biggest lever: endosomal escape rate (not ultrasound parameters). Path to viability = better LNP chemistry, not stronger ultrasound.

## Status

Computational model only. Sonoporation through RWM demonstrated in guinea pigs (Shih 2019). LNP + sonoporation combination not tested in cochlea. This is a theoretical framework.

## Open questions

- Can optimized LNP formulations achieve 30% endosomal escape in OHCs?
- What is actual pore diameter and permeability in human RWM vs animal models?
- Does repeated sonoporation cause cumulative RWM damage?
- Can mRNA expression be sustained long enough for functional recovery?

## Recent Papers

### 2026-04-21 — [Noninvasive Focal Gene Delivery into the Cerebellum of Non-Human Primates using Focused Ultrasound](https://pubmed.ncbi.nlm.nih.gov/42003738/)
*Advanced Science · relevance: 5/10 · lateral*

FUS + IV microbubbles achieved non-invasive, spatially precise AAV9 delivery to macaque deep cerebellar nuclei with robust neuronal transduction. The blood-labyrinth barrier is structurally analogous to the BBB, and this NHP demonstration of FUS-mediated AAV delivery strengthens the feasibility case for non-invasive cochlear gene delivery. Combined with our existing ODE model (Scenario 2/3), FUS-mediated barrier opening could bypass the surgical round window approach entirely, enabling bilateral treatment and repeat dosing of AAV-mini-STRC or mRNA-LNPs.

### 2026-04-16 — [Gene Therapy for Infants and Children With Otoferlin-Related Auditory Neuropathy Spectrum Disorder](https://pubmed.ncbi.nlm.nih.gov/41979424/)
*Ear and Hearing · relevance: 7/10 · fda_path*

OTOF clinical trials establish real-world data on cochlear gene therapy delivery in pediatric patients, including surgical approach, dose volumes, and safety monitoring. This clinical framework directly informs regulatory expectations for cochlear delivery — both surgical (standard approach) and any non-surgical alternatives like our sonoporation concept. Understanding delivery-related adverse events from OTOF trials will help benchmark the safety case for alternative delivery routes.

### 2026-04-20 — [FUS + perfluorocarbon nanodroplets for barrier disruption and enhanced delivery](https://pubmed.ncbi.nlm.nih.gov/41993780/)
*Int J Nanomedicine · relevance: 3/10 · lateral*

Achieved 28-fold reduction in systemic clearance and 33% long-term survival in glioblastoma using 62 nm perfluorocarbon nanodroplets combined with FUS-mediated BBB disruption. The blood-labyrinth barrier is structurally analogous to the BBB. This dual approach (barrier opening + nanocarrier) could enable non-invasive cochlear delivery of AAV-mini-STRC or mRNA-LNPs without surgical round window access, potentially enabling bilateral treatment and repeat dosing. The 62 nm particle size is within range for cochlear tissue penetration.

### 2026-04-20 — [FUS-mediated gene therapy: preclinical systematic review and meta-analysis](https://pubmed.ncbi.nlm.nih.gov/41955864/)
*Biomed Pharmacother · relevance: 4/10 · lateral*

Meta-analysis of 9 studies: FUS enhanced gene expression 6.34x (95% CI 2.21-18.18) across viral, non-viral, nanoparticle, and exosome vectors. The vector-agnostic nature of FUS enhancement is key — it would likely work with AAV-mini-STRC. At the 6.34x enhancement factor, the required AAV dose for therapeutic OHC transduction could drop proportionally, directly addressing our sensitivity analysis showing endosomal escape as the primary bottleneck. Combined with our Scenario 3 hybrid approach (AAV + LNP top-up), FUS could make the initial AAV dose more efficient.
