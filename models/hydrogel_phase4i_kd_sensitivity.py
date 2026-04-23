#!/usr/bin/env python3
"""
Hydrogel Phase 4i — Kd × Kd therapeutic-window sensitivity sweep.

Gate question: does the h09 therapeutic claim survive across the full plausible
Kd range for both load-bearing unmeasured parameters?

  • KD_TMEM145 (h09 gate 3): ipSAE reassessment 2026-04-23 places
    UM × TMEM145 GOLD ipSAE 0.591 in Calcineurin-family binder band
    (NFAT-CnB 0.55 @ ~0.5-2 μM … NFAT-CnA 0.78 @ ~0.1-1 μM).
    Sweep: 10 nM → 10 μM (3 log units).

  • WH2_KD_FACTIN (h09 gate 2): no primary measurement published.
    Optimistic current assumption 5 μM (25× weakening vs G-actin).
    Tβ4 F-actin analog 5-10 mM floor.
    Sweep: 1 μM → 1 mM (3 log units).

For each (Kd_T, Kd_W) grid point, compute:
  (a) Minimum local perilymph [L] required for functional HTC fraction
        f = θ_T × θ_W ≥ 0.3  (22 dB ABR rescue target per Tobin 2019 Hopf coupling)
      θ_T = [L]/(Kd_T + [L])   — TMEM145 anchor occupancy
      θ_W = [L]/(Kd_W + [L])   — F-actin hook occupancy
  (b) Intratympanic dose required to achieve that peak [L] in perilymph
      via 2-compartment PKPD (reuses Phase 4e model — same rate constants).
  (c) PASS (dose ≤ 5 mg single intratympanic, tolerable per AM-101 precedent) /
      MARGINAL (5-20 mg, feasible with split dose or multi-depot) /
      FAIL (> 20 mg or [L] exceeds 100 μM toxicity ceiling for G-actin
      depletion).

Promotion-to-S decision rule:
  • ≥ 70% of 5×5 grid PASS → robust therapeutic window → evidence for
    S-tier promotion without wet-lab Kd narrowing.
  • 30-70% PASS → window exists but requires wet-lab Kd narrowing to
    identify the operable corner. Stay A.
  • < 30% PASS → therapy is parametric on unmeasured Kds. Stay A.
    Wet-lab is gate.

═══════════════════════════════════════════════════════════════════════════
CAVEAT: This is a SENSITIVITY sweep, not a closure. Wet-lab SPR/BLI for
TMEM145 and actin-bundling assay for WH2 remain the authoritative gates
for absolute Kd. This script only asks: "given the parametric uncertainty
we have today, is the design robust enough for S-tier action?"
═══════════════════════════════════════════════════════════════════════════
"""

from __future__ import annotations

import json
import numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp

OUT = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase4i_kd_sensitivity.json")

# PKPD params — Phase 4e post-audit, all lit-anchored except K_PROTEOLYSIS
PEPTIDE_MW_KDA = 14.2
PERILYMPH_VOL_UL = 34.0          # Dhanasingh 2021
K_CLEAR_MIDDLE_EAR = 0.7         # Salt & Plontke 2018
K_RWM = 0.003                    # Stokes-Einstein from Salt & Ma 2001
K_PERILYMPH_CLEAR = 0.18         # Salt & Hartsock 2015
K_PROTEOLYSIS = 0.05             # still ⚠ unmeasured, Phase 4e flag
TOXIC_CEILING_uM = 100.0         # Phase 4d G-actin depletion threshold
DOSE_TOLERABLE_MG = 5.0          # AM-101 intratympanic precedent (0.5-2 mg typical)
DOSE_MARGINAL_MG = 20.0          # above this = not a single clinic visit

# HTC functional thresholds (from [[STRC Stereocilia Bundle Mechanics Model]])
F_TARGET_CONSERVATIVE = 0.3      # 22 dB ABR rescue
F_TARGET_STRETCH = 0.6           # 30 dB ABR rescue

# Kd sweep ranges
KD_TMEM145_GRID_M = np.array([10e-9, 100e-9, 500e-9, 2e-6, 10e-6])
KD_WH2_FACTIN_GRID_M = np.array([1e-6, 5e-6, 50e-6, 200e-6, 1000e-6])


def pkpd_derivative(t, y):
    c1, c2, _ = y
    flux_12 = c1 * K_RWM
    flux_1_out = c1 * K_CLEAR_MIDDLE_EAR
    flux_2_clear = c2 * K_PERILYMPH_CLEAR
    flux_2_deg = c2 * K_PROTEOLYSIS
    dc1 = -flux_12 - flux_1_out
    dc2 = flux_12 - flux_2_clear - flux_2_deg
    ddeg = flux_2_deg + flux_1_out + flux_2_clear
    return [dc1, dc2, ddeg]


def amount_to_conc_uM(amount_mg: float) -> float:
    moles = (amount_mg * 1e-3) / (PEPTIDE_MW_KDA * 1000)
    vol_L = PERILYMPH_VOL_UL * 1e-6
    return (moles / vol_L) * 1e6


def peak_perilymph_uM(dose_mg: float) -> float:
    y0 = [dose_mg, 0.0, 0.0]
    sol = solve_ivp(pkpd_derivative, [0, 48], y0,
                    t_eval=np.linspace(0, 48, 400),
                    method="LSODA", rtol=1e-6, atol=1e-12)
    peri_uM = np.array([amount_to_conc_uM(x) for x in sol.y[1]])
    return float(peri_uM.max())


def dose_for_peak(target_uM: float, lo=1e-4, hi=200.0) -> float:
    """Bisection: find dose (mg) giving peak perilymph == target_uM."""
    for _ in range(50):
        mid = (lo + hi) / 2
        if peak_perilymph_uM(mid) < target_uM:
            lo = mid
        else:
            hi = mid
    return mid


def min_L_for_f(f_target: float, kd_t_uM: float, kd_w_uM: float,
                search_L_range_uM=(0.01, 1e5)) -> float | None:
    """
    Solve [L] such that f = ([L]/(Kd_T+[L])) * ([L]/(Kd_W+[L])) = f_target.
    Returns smallest [L] that hits target, or None if impossible.
    """
    def f_at_L(L):
        return (L / (kd_t_uM + L)) * (L / (kd_w_uM + L))

    lo, hi = search_L_range_uM
    if f_at_L(hi) < f_target:
        return None    # even at max [L], can't hit target
    for _ in range(80):
        mid = np.sqrt(lo * hi)   # log-scale bisection
        if f_at_L(mid) < f_target:
            lo = mid
        else:
            hi = mid
    return mid


def flag_outcome(L_required_uM: float | None, dose_mg: float | None) -> str:
    if L_required_uM is None:
        return "FAIL_unreachable"
    if L_required_uM > TOXIC_CEILING_uM:
        return "FAIL_toxic"
    if dose_mg is None or dose_mg > DOSE_MARGINAL_MG:
        return "FAIL_dose"
    if dose_mg > DOSE_TOLERABLE_MG:
        return "MARGINAL"
    return "PASS"


def main():
    grid_results = []
    for kd_t in KD_TMEM145_GRID_M:
        row = []
        for kd_w in KD_WH2_FACTIN_GRID_M:
            kd_t_uM, kd_w_uM = kd_t * 1e6, kd_w * 1e6
            L_conserv = min_L_for_f(F_TARGET_CONSERVATIVE, kd_t_uM, kd_w_uM)
            L_stretch = min_L_for_f(F_TARGET_STRETCH, kd_t_uM, kd_w_uM)
            dose_conserv = dose_for_peak(L_conserv) if (
                L_conserv and L_conserv <= TOXIC_CEILING_uM) else None
            dose_stretch = dose_for_peak(L_stretch) if (
                L_stretch and L_stretch <= TOXIC_CEILING_uM) else None
            flag_c = flag_outcome(L_conserv, dose_conserv)
            flag_s = flag_outcome(L_stretch, dose_stretch)
            row.append({
                "kd_tmem145_uM": kd_t_uM,
                "kd_wh2_factin_uM": kd_w_uM,
                "min_L_for_f0.3_uM": (round(L_conserv, 3)
                                      if L_conserv else None),
                "dose_for_f0.3_mg": (round(dose_conserv, 3)
                                     if dose_conserv else None),
                "flag_f0.3": flag_c,
                "min_L_for_f0.6_uM": (round(L_stretch, 3)
                                      if L_stretch else None),
                "dose_for_f0.6_mg": (round(dose_stretch, 3)
                                     if dose_stretch else None),
                "flag_f0.6": flag_s,
            })
        grid_results.append(row)

    # Robustness: % of grid cells that PASS or MARGINAL (achievable) vs FAIL
    flat = [c for row in grid_results for c in row]
    n_total = len(flat)
    n_pass_c = sum(1 for c in flat if c["flag_f0.3"] == "PASS")
    n_marginal_c = sum(1 for c in flat if c["flag_f0.3"] == "MARGINAL")
    n_fail_c = sum(1 for c in flat if c["flag_f0.3"].startswith("FAIL"))
    n_pass_s = sum(1 for c in flat if c["flag_f0.6"] == "PASS")
    n_marginal_s = sum(1 for c in flat if c["flag_f0.6"] == "MARGINAL")
    n_fail_s = sum(1 for c in flat if c["flag_f0.6"].startswith("FAIL"))

    pct_pass_c = n_pass_c / n_total * 100
    pct_reachable_c = (n_pass_c + n_marginal_c) / n_total * 100
    pct_pass_s = n_pass_s / n_total * 100
    pct_reachable_s = (n_pass_s + n_marginal_s) / n_total * 100

    if pct_pass_c >= 70:
        promotion_verdict = ("S_candidate_robust: ≥70% of plausible Kd space "
                             "is single-dose tolerable at f=0.3")
    elif pct_reachable_c >= 70:
        promotion_verdict = ("A_hold_multi_dose_ok: therapy reachable on ≥70% "
                             "grid but requires marginal doses in parts — "
                             "wet-lab Kd narrowing to find operable corner")
    elif pct_reachable_c >= 30:
        promotion_verdict = ("A_hold_wetlab_gated: partial robustness; wet-lab "
                             "must narrow at least one Kd before promotion")
    else:
        promotion_verdict = ("A_at_risk: therapy parametric on unmeasured Kds; "
                             "broad plausible range is FAIL — wet-lab is gate")

    summary = {
        "batch": "hydrogel_phase4i_kd_sensitivity",
        "date": "2026-04-23",
        "model_assumptions": {
            "pkpd_params": "Phase 4e post-audit, lit-anchored except K_PROTEOLYSIS",
            "f_target_conservative": F_TARGET_CONSERVATIVE,
            "f_target_stretch": F_TARGET_STRETCH,
            "toxic_ceiling_perilymph_uM": TOXIC_CEILING_uM,
            "dose_tolerable_mg": DOSE_TOLERABLE_MG,
            "dose_marginal_ceiling_mg": DOSE_MARGINAL_MG,
            "kd_tmem145_grid_M": KD_TMEM145_GRID_M.tolist(),
            "kd_wh2_factin_grid_M": KD_WH2_FACTIN_GRID_M.tolist(),
            "notes": [
                "f = theta_T * theta_W assumes independent binding sites at HTC",
                "KD_TMEM145 range based on ipSAE placement in Calcineurin "
                "family + 10x margin; absolute Kd still wet-lab gated",
                "KD_WH2_FACTIN range spans optimistic 25x weakening from G-actin "
                "to Tbeta4 analog floor 1 mM (10x margin below Tbeta4 low end)",
                "Dose reflects single intratympanic injection; multi-depot or "
                "daily dosing possible but adds compliance burden",
            ],
        },
        "grid_results": grid_results,
        "robustness_stats": {
            "n_grid_points": n_total,
            "f_0.3_conservative": {
                "pct_PASS": round(pct_pass_c, 1),
                "pct_MARGINAL": round(n_marginal_c / n_total * 100, 1),
                "pct_FAIL": round(n_fail_c / n_total * 100, 1),
                "pct_reachable_PASS_or_MARGINAL": round(pct_reachable_c, 1),
            },
            "f_0.6_stretch": {
                "pct_PASS": round(pct_pass_s, 1),
                "pct_MARGINAL": round(n_marginal_s / n_total * 100, 1),
                "pct_FAIL": round(n_fail_s / n_total * 100, 1),
                "pct_reachable_PASS_or_MARGINAL": round(pct_reachable_s, 1),
            },
        },
        "promotion_verdict": promotion_verdict,
    }

    OUT.write_text(json.dumps(summary, indent=2))

    # Human-readable
    print("=== Phase 4i Kd × Kd Therapeutic Window Sensitivity ===\n")
    print("Grid: KD_TMEM145 (rows) × KD_WH2_FACTIN (cols) → flag for f≥0.3\n")
    header = "Kd_T\\Kd_W    " + "  ".join(
        f"{k*1e6:7.2g}μM" for k in KD_WH2_FACTIN_GRID_M)
    print(header)
    for i, kd_t in enumerate(KD_TMEM145_GRID_M):
        row_str = f"{kd_t*1e6:7.2g}μM     "
        for j in range(len(KD_WH2_FACTIN_GRID_M)):
            cell = grid_results[i][j]
            row_str += f"{cell['flag_f0.3']:>9s} "
        print(row_str)

    print(f"\nf=0.3 conservative (22 dB ABR): PASS {pct_pass_c:.0f}%, "
          f"MARGINAL {n_marginal_c/n_total*100:.0f}%, "
          f"FAIL {n_fail_c/n_total*100:.0f}% "
          f"(reachable {pct_reachable_c:.0f}%)")
    print(f"f=0.6 stretch      (30 dB ABR): PASS {pct_pass_s:.0f}%, "
          f"MARGINAL {n_marginal_s/n_total*100:.0f}%, "
          f"FAIL {n_fail_s/n_total*100:.0f}% "
          f"(reachable {pct_reachable_s:.0f}%)")
    print(f"\nVerdict: {promotion_verdict}")
    print(f"\nJSON written: {OUT}")


if __name__ == "__main__":
    main()
