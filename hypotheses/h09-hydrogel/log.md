---
type: hypothesis-log
hypothesis: h09
tags: [strc, h09]
---

# h09 log

## 2026-04-23

- Phase 4e_v2 blend scaffold model: closes **1/3 S-tier gate** (118 aa tail > 12 aa empirical RADA16 limit). Recommended design point φ_tail=0.05 (5% RADA16-tail91, 95% plain), 2% w/v gel: PASS flag (within Gelain hybrid lit range), 6.4×10⁷ tails per OHC bundle shell, 2000× valency excess vs putative 3×10⁴ TMEM145 sites/bundle low est, 3.8× dose cost vs hypothetical 100% tail91. Gates 2 (WH2×F-actin Kd) + 3 (STRC×TMEM145 Kd) still open (wet-lab or Path B AF3 calibration). A held, Mech 4 held, Deliv 3 held. → [[STRC h09 Phase 4e_v2 Blend Scaffold]]
- Hub-note status banner added: both load-bearing unmeasured Kds (STRC×TMEM145 100 nM placeholder; WH2×F-actin 5 μM optimistic vs Tβ4 5-10 mM analog) + tail-91 > 12 aa RADA16 modification limit logged as promotion-to-S gates. Unsourced "48 h actin t½" in note replaced with Zhang 2012 Nature primary data (shaft stable months; tip β-actin t½ hours). A held. → [[STRC Synthetic Peptide Hydrogel HTC]]
- Parameter provenance audit: 3 critical gaps (118 aa tail, WH2×F-actin Kd, STRC×TMEM145 Kd); 7 phantom/wrong values; Mech 5→4 Deliv 4→3 (A held); triple promotion gate to S → [[STRC h09 Parameter Provenance Audit 2026-04-23]]
- Phase 4d/4e/4h re-run with lit-backed parameters: PERILYMPH_VOL 70→34 μL (Dhanasingh 2021); K_RWM 0.02→0.003 /h (Salt & Ma 2001, prior PHANTOM); K_CLEAR_ME 0.35→0.7 (Salt & Plontke 2018); K_PERILYMPH 0.35→0.18 (Salt & Hartsock 2015); K_PROTEOLYSIS 1.4→0.05 (⚠ still unmeasured); STEREOCILIA_SPACING 12→9 nm (Krey 2016 plastin, was espin-specific). Therapeutic window (1-10 μM, 2 log units) survives swap. WH2×F-actin Kd + STRC×TMEM145 Kd remain load-bearing placeholders. A held. → [[STRC h09 Parameter Provenance Audit 2026-04-23]]
- Phase 3b full-construct: tail91 PASSES TMEM145 0.57 + actin 0.51; tail71 FAILS TMEM145 catastrophically (0.35) → [[STRC Hydrogel HTC Phase 1 Self-Assembly]]
- B → A: mechanism 3→5; Phase 3 tail retool ipTM 0.68 = Ultra-Mini GOLD baseline → [[STRC Hydrogel HTC Phase 1 Self-Assembly]]
- Phase 4 8-axis campaign: Delivery 3→4; PKPD window 0.13-1.32 mg; $42/ear GMP; MODERATE immunogenicity → [[STRC Hydrogel Phase 4 Computational Campaign]]
- Phase 4b AF3 batch built (8 jobs); awaiting submission → [[STRC Hydrogel Phase 4 Computational Campaign]]
- Phase 3 tail retool: design bug found (wrong epitope Phase 1); 91 aa correct-epitope tail ipTM 0.68 → [[STRC Hydrogel HTC Phase 1 Self-Assembly]]
- Phase 3 retool batch built (3 jobs): 51/71/91 aa tails; gate ipTM ≥ 0.50 → [[STRC Hydrogel HTC Phase 1 Self-Assembly]]
- Phase 2 AF3: actin PASS (0.59), TMEM145 FAIL (0.37-0.39, wrong-epitope tail confirmed) → [[STRC Hydrogel HTC Phase 1 Self-Assembly]]
- Phase 2 AF3 batch built (6 jobs: 3 candidates × TMEM145 + actin trimer) → [[STRC Hydrogel HTC Phase 1 Self-Assembly]]
- Phase 1 analytical: 5/6 pass all gates; top-3 shortlisted → [[STRC Hydrogel HTC Phase 1 Self-Assembly]]

## Connections

- `[part-of]` [[h09 hub]]
