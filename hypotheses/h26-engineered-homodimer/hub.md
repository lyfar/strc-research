---
date: 2026-04-23
type: hypothesis-hub
tags: [strc, h26]
status: B-tier
---

# h26 — Engineered Homodimer Avidity

**mech.** Single-point mutation at Ultra-Mini homodimer interface (ARM 1579-1581) tested for dimer stabilization; Phase 1 AF3 FAILED — all 4 mutants destabilize dimer further, R-R repulsion hypothesis falsified.
**delivery.** Would be drop-in to h03 Mini-STRC AAV (Anc80L65 + B8-IgK-Ultra-Mini-WPRE3-bGH) if any mutant passed; currently none does.
**patient-fit.** Misha-compatible by inheritance from h03; not independently therapeutic.

## status

B-tier (was A, 2026-04-23). **Phase 1c DELIVERED 2026-04-23 GREEN**: DBSCAN eps=6.5 Å re-cluster found 3 clusters; Cluster 1 "stump" 1077-1114 (25 residues, 21 Cb-Cb disulfide geometries 4.5-7.5 Å) is entirely outside Phase 1's 1579-1581 target. Native C1081 forms inter-chain Cb-Cb 6.87-7.09 Å with A1078/C1079/S1080 on opposing chain → single A1078C or S1080C engineers inter-chain disulfide. Cluster 2 ARM also has homotypic S1579→C Cb-Cb 6.94 Å (re-opens Phase 1 target with disulfide chemistry). **Phase 1d AF3 triple batch queued** (A1078C / S1080C / S1579C homotypic + A1078W negative control, user-invoked external AF3). Promotion path: Phase 1d ipTM ≥0.50 + homodimer contacts ≥ WT → A-tier → drop-in to h03 Ultra-Mini. See [[STRC h26 Phase 1c Contact Re-Cluster 2026-04-23]].

<!-- RECENT:START auto-synced from log.md by sync_strc_indices.py — do not hand-edit -->

## recent

- Phase 1c contact re-cluster: **GREEN**. DBSCAN eps=6.5 Å on Ultra-Mini homodimer CIF found 3 clusters. **Cluster 1 "stump" 1077-1114 (25 residues, 21 Cb-Cb disulfide geometries) ENTIRELY outside Phase 1 target 1579-1581.** Native C1079+C1081 on both chains sit 6.87-7.09 Å from A1078/S1080 on opposing chain → single A1078C or S1080C engineers inter-chain disulfide without changing cystine count. Cluster 2 ARM has homotypic S1579C/S1579C at Cb-Cb 6.94 Å (re-opens Phase 1 target with disulfide chemistry). B held pending Phase 1d AF3 triple batch (A1078C + S1080C + S1579C homotypic + A1078W negative control). Promotion path: Phase 1d ipTM ≥0.50 + homodimer contacts ≥ WT → A-tier → drop into h03 Ultra-Mini as engineered-avidity variant (zero payload cost). → [[STRC h26 Phase 1c Contact Re-Cluster 2026-04-23]]
- Added to register: A-tier, Mech 4, Deliv 5, Misha 4. R-R repulsion at ARM 1579-1581 identified as dimer destabilizer → [[STRC Engineered Homodimer Avidity]]
- Phase 1 AF3 batch built: 5 jobs (WT + 4 single-point mutants × 2 TMEM145), awaiting submission → [[STRC Engineered Homodimer Avidity]]
- Phase 1 AF3 returned 4 mutants (WT_ref not downloaded): all FAIL homodimer gate (mean 0.15-0.19 vs WT 0.28-0.30); 3/4 also destroy TMEM145 binding → [[STRC Engineered Homodimer Phase 1 Results]]
- Tier A → B. Mech 4→3. R-R repulsion falsified by R1581F (arg removed, homodimer still weak) → [[STRC Engineered Homodimer Phase 1 Results]]

<!-- RECENT:END -->

## evidence

- [[STRC Engineered Homodimer Avidity]] — main hypothesis (R-R repulsion prediction, now falsified)
- [[STRC Engineered Homodimer Phase 1 Results]] — Phase 1 AF3 results: 4/4 mutants fail homodimer gate
- [[STRC Homodimer Interface From CIF]] — source: weak homodimer at ARM 1579-1581 identified

## scripts

research/strc/models/ (legacy pool — see [[STRC Computational Scripts Inventory]] § h26)
- `af3_jobs_2026-04-23d_engineered_homodimer_builder.py` — 5 jobs: WT ref + R1581Y + R1581F + S1579W + S1579F
- `af3_jobs_2026-04-23d_engineered_homodimer_forensics.py` — Phase 1 results parser (4 mutants × 1-2 seeds, gate evaluation)
- `engineered_homodimer_phase1c_contact_cluster.py` — DBSCAN 6.5 Å spatial re-cluster on Ultra-Mini homodimer inter-chain contacts; per-cluster disulfide Cb-Cb scan 4.5-7.5 Å + mutation candidate ranking (steric / electrostatic / covalent)

## artifacts

research/strc/hypotheses/h26-engineered-homodimer/artifacts/phase1_af3_2026-04-23/ — 6 AF3 runs (4 mutants + 2 restart seeds) + analysis_summary.json

## log

[[h26 log]]

## cross

- [[STRC Mini-STRC Single-Vector Hypothesis]] — h26 is a direct enhancement of h03 clinical candidate
- [[STRC Homodimer Interface From CIF]] — source analysis that identified the interface

## Connections

- `[part-of]` [[STRC Hypothesis Ranking]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[about]` [[Misha]]
