---
date: '2026-04-24'
type: hypothesis-log
aliases:
- h02 log
tags:
- strc
- h02
hypothesis: h02
---

# h02 log

## 2026-04-23

- AAV-LNP stack model confirms Piezo S-tier independent of AAV path → [[STRC AAV-LNP Stack PKPD]]
- Parameter provenance audit complete: 3 phantom cites caught (K_D A666-prestin = 10 nM; ETA_POLY = 0.7; d31 Phase1/Phase2 inconsistency). 9 additional ⚠ constants. New topic files: [[piezoelectric-materials]], [[tectorial-membrane]]. Full audit at `/tmp/strc-lit-params-h02.md`. New finding: TM stiffness (24–210 kPa) vs. PVDF-TrFE (3 GPa) mismatch not currently modeled — potential hypothesis-level mechanical flaw.
- Post-audit fix applied: scripts pinned phantoms (A666, ETA_POLY, SELECT_S), d31 reconciled to −12 pC/N (was silent 2× split). Phase 1 baseline no longer passes 60 dB. TM mismatch 10⁵× flagged in output JSON. Tier S → B. Mech 3→2, Deliv 2→1. → [[STRC h02 Parameter Provenance Audit 2026-04-23]]

## 2026-04-20

- Phases 1-3 complete; delivery model 92% audiogram coverage → [[STRC Piezo Delivery Feasibility OHC Targeting]]
- Initial assignment: S-tier, Mech 3, Deliv 2, Misha 5 → [[STRC Hypothesis Ranking]]

## Connections

- `[part-of]` [[h02 hub]]
