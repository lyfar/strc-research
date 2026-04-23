---
date: '2026-04-24'
type: hypothesis-log
aliases:
- h26 log
tags:
- strc
- h26
hypothesis: h26
---

# h26 log

## 2026-04-23

- Phase 1c contact re-cluster: **GREEN**. DBSCAN eps=6.5 Å on Ultra-Mini homodimer CIF found 3 clusters. **Cluster 1 "stump" 1077-1114 (25 residues, 21 Cb-Cb disulfide geometries) ENTIRELY outside Phase 1 target 1579-1581.** Native C1079+C1081 on both chains sit 6.87-7.09 Å from A1078/S1080 on opposing chain → single A1078C or S1080C engineers inter-chain disulfide without changing cystine count. Cluster 2 ARM has homotypic S1579C/S1579C at Cb-Cb 6.94 Å (re-opens Phase 1 target with disulfide chemistry). B held pending Phase 1d AF3 triple batch (A1078C + S1080C + S1579C homotypic + A1078W negative control). Promotion path: Phase 1d ipTM ≥0.50 + homodimer contacts ≥ WT → A-tier → drop into h03 Ultra-Mini as engineered-avidity variant (zero payload cost). → [[STRC h26 Phase 1c Contact Re-Cluster 2026-04-23]]
- Added to register: A-tier, Mech 4, Deliv 5, Misha 4. R-R repulsion at ARM 1579-1581 identified as dimer destabilizer → [[STRC Engineered Homodimer Avidity]]
- Phase 1 AF3 batch built: 5 jobs (WT + 4 single-point mutants × 2 TMEM145), awaiting submission → [[STRC Engineered Homodimer Avidity]]
- Phase 1 AF3 returned 4 mutants (WT_ref not downloaded): all FAIL homodimer gate (mean 0.15-0.19 vs WT 0.28-0.30); 3/4 also destroy TMEM145 binding → [[STRC Engineered Homodimer Phase 1 Results]]
- Tier A → B. Mech 4→3. R-R repulsion falsified by R1581F (arg removed, homodimer still weak) → [[STRC Engineered Homodimer Phase 1 Results]]
- Dunbrack 2025 ipSAE reassessment: Ultra-Mini homodimer **ipSAE 0.000** (zero interface residues pass PAE<10 Å cutoff at default AF3 settings) — ipTM 0.28-0.30 uplift from baseline 0.20-0.24 now identified as disordered-region false-positive pattern ipSAE was designed to catch. Mech 3→2, B held. Structural symmetry (94% C2) + aa 1579-1581 self-contacts legs survive the reassessment (they are geometry observations independent of the confidence metric). Phase 1c 5 Å re-cluster still the next viable move. → [[STRC ipSAE Cross-Complex Reassessment 2026-04-23]]

## Connections

- `[part-of]` [[h26 hub]]
