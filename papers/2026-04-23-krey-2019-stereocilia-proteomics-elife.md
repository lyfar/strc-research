---
title: "Single-cell proteomics reveals changes in expression during hair-cell development"
authors: ["Krey JF", "Chatterjee P", "Halford J", "Cunningham CL", "Perrin BJ", "Bharat TAM", "Bharat R", "Bharat A", "Barr-Gillespie PG"]
journal: "eLife"
year: 2019
volume: 8
pages: "e50777"
pubmed_id: "31724949"
pmc: "PMC6855842"
doi: "10.7554/eLife.50777"
date: 2026-04-23
type: source
tags: [strc, paper, stereocilia, proteomics, actin, hair-cell, barr-gillespie, h09, h26]
hypotheses: [h09, h26]
---

## What they found

Highly sensitive single-cell proteomics (SCoPE-MS) applied to E15 chick utricle hair cells.
Quantified relative molar fractions of ~1500 proteins per cell type. Provides the most precise
molecular-level census of hair cell protein composition, including actin copy numbers that
can be used to estimate local [F-actin] in stereocilia.

## Numbers that matter

**Actin copy numbers:**
- E15 chick utricle hair cell: ~15,000,000 total actin molecules per cell
- Per stereocilium: ~200,000 actin molecules (72 stereocilia counted per utricle hair cell)
- E20 (later development): ~400,000 actin molecules per stereocilium (estimate)
- ACTG1 group (total actin) = 0.043 ± 0.001 molar fraction in FM1-43high cells

**Derived local [F-actin] in stereocilia (calculation, not in paper):**
- OHC stereocilium: ~200 nm diameter, ~3 μm long → volume ≈ 9.4 × 10⁻¹⁷ L
- 200,000 molecules / (6.022×10²³ × 9.4×10⁻¹⁷) ≈ 3,500 μM total actin
- Assuming 90% filamentous → [F-actin] ≈ 3,150 μM ≈ 3.15 mM
- Cross-check with ~100 nm diameter: [actin] ≈ 14,000 μM (extreme upper bound, likely for
  developing stereocilia with denser core)
- Practical working estimate: [F-actin] in mature stereocilium core = 1–5 mM

**Critical implication for WH2 × F-actin bundling model:**
- [F-actin] is 1–5 mM — this is far above any plausible WH2_KD_FACTIN_M
- If WH2_KD_FACTIN_M = 5 μM (as in model), fractional occupancy at 1 mM [F-actin] would
  approach 100% — but this is irrelevant if the binding event doesn't exist in the first place
- The high local [F-actin] is not a bottleneck; the question is purely whether WH2 can
  side-bind at all

**STRC/stereocilin:** Not quantified in this study (chick utricle dataset; STRC is mammalian
cochlear-specific in its critical role).

## Connections

- [[Hydrogel Phase 4d F-actin Bundling Model]] — provides
  200,000 actin/stereocilium baseline for [F-actin] calculation
- [[STRC Normal OHC Concentration Parameter]] — STRC_NORMAL_OHC_M needs own
  source; this paper doesn't cover STRC
- [[Barr-Gillespie Lab Stereocilia Atlas]] — anchor paper for
  quantitative hair cell proteomics
