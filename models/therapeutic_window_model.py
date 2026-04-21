#!/usr/bin/env python3
"""
Therapeutic Window Model for STRC/DFNB16 Gene Therapy
======================================================
ODE model of OHC stereocilia bundle degradation and hair cell death
in STRC-null inner ears. Identifies postnatal intervention deadline
for mini-STRC gene therapy.

Biology:
- STRC forms horizontal top connectors (HTCs) between stereocilia rows
- HTCs required for proper MET channel gating geometry
- Without HTCs: bundles disorganize → MET fails → OHC functional death
- OHC structural death follows (Ca²⁺ dysregulation, loss of trophic signals)

Mouse calibration data (Verpy et al. 2011, Nat Neurosci; Grillet et al. 2009):
- Strc−/− mice: no HTCs at P12 (hearing onset), bundles already disorganized
- P12: bundles disrupted, ABR ~60 dB (normal ~20 dB) — 40 dB shift
- P30: ABR ~70-75 dB — progressive
- P90: ABR ~80 dB or absent — near-complete
- OHC count: ~70% survive to P30, ~40% survive to P90 (estimated from DPOAE)

Mouse → Human scaling (Rubel & Fritzsch 2002; Müller & Bhatt 2016):
- Mouse inner ear development compressed ~10-fold vs human
- Mouse P12 ≈ human newborn (gestational hearing onset ~28 wk)
- Mouse P30 ≈ human ~3 months postnatal
- Mouse P90 ≈ human ~9-12 months postnatal

Author: Egor Lyfar (computational)
Date: 2026-04-16
Motivation: Tsai et al. 2026 (JCI) — postnatal Slc26a4 gene therapy
identified critical therapeutic window for hereditary HL. STRC/DFNB16
requires equivalent characterization for IND planning.
"""

import numpy as np
from scipy.integrate import solve_ivp
import json

# ─── Parameters ───────────────────────────────────────────────────────────────

# Bundle degradation rate constant [1/day], mouse
# Calibrated: integrity drops to ~0.2 by P30 (40 dB shift observed)
# OHC bundles intact at P0, fully disrupted by ~P60-P90 in Strc-/-
K_DEGRADE_MOUSE = 0.045      # [1/day]

# Below this bundle integrity threshold, OHC begins functional death
BUNDLE_THRESHOLD = 0.35      # normalized [0-1]

# OHC death rate constant once below threshold [1/day]
K_DEATH = 0.025              # calibrated: ~40% OHC loss by P90

# Rescue efficiency: fraction of functional OHCs restored by mini-STRC
# (if delivered before irreversible death)
# Based on Anc80L65 OHC transduction: 60-80% (Landegger 2017)
RESCUE_EFFICIENCY = 0.70

# Human/mouse developmental scaling factor (Rubel & Fritzsch 2002)
HUMAN_MOUSE_SCALE = 10.0     # P-day × 10 ≈ postnatal days in human equivalent

# Minimum OHC survival for meaningful hearing rescue (ABR threshold <50 dB)
# Calibrated from Bredberg 1968 OHC-count vs threshold correlation
MIN_OHC_FOR_RESCUE = 0.30    # 30% survival minimum for partial rescue
MIN_OHC_FOR_GOOD_OUTCOME = 0.55  # 55% for <40 dB threshold (functional hearing)

# ─── ODE Model ────────────────────────────────────────────────────────────────

def ohc_degradation(t, y, strc_present=False, rescue_time=None):
    """
    ODE system for STRC-null OHC stereocilia degradation.

    State vector y:
      y[0] = bundle_integrity  (0=fully degraded, 1=intact)
      y[1] = ohc_survival      (0=all dead, 1=all alive)

    With strc_present=True: bundle integrity maintained (therapy worked)
    rescue_time: P-day when mini-STRC was delivered
    """
    bundle_integrity, ohc_survival = y

    # STRC restores HTC formation — stops further degradation
    if strc_present and rescue_time is not None and t >= rescue_time:
        k_deg = 0.0
    elif strc_present and rescue_time is None:
        k_deg = 0.0
    else:
        k_deg = K_DEGRADE_MOUSE

    d_bundle = -k_deg * bundle_integrity

    # OHC death only kicks in below bundle integrity threshold
    deficit = max(0.0, BUNDLE_THRESHOLD - bundle_integrity)
    d_ohc = -K_DEATH * deficit * ohc_survival

    return [d_bundle, d_ohc]


def run_scenario(rescue_day_mouse=None, t_max=120, dt=0.5):
    """
    Simulate OHC survival from P0 to t_max (P-days, mouse).

    rescue_day_mouse: P-day when mini-STRC is delivered (None = untreated)
    Returns dict with time array and outcome metrics.
    """
    t_span = (0, t_max)
    t_eval = np.arange(0, t_max + dt, dt)
    y0 = [1.0, 1.0]  # intact bundles, all OHCs alive at P0

    strc = rescue_day_mouse is not None

    sol = solve_ivp(
        lambda t, y: ohc_degradation(t, y, strc_present=strc, rescue_time=rescue_day_mouse),
        t_span, y0, t_eval=t_eval, method='RK45', rtol=1e-6
    )

    bundle = sol.y[0]
    survival = sol.y[1]
    final_survival = survival[-1]

    # Apply rescue efficiency: transduced OHCs are rescued, non-transduced continue to die
    if rescue_day_mouse is not None:
        idx_rescue = np.searchsorted(sol.t, rescue_day_mouse)
        survival_at_rescue = survival[idx_rescue] if idx_rescue < len(survival) else survival[-1]
        # Rescued fraction: transduced OHCs stop dying
        # Non-transduced: continue degradation to t_max
        rescued = survival_at_rescue * RESCUE_EFFICIENCY
        sol2 = solve_ivp(
            lambda t, y: ohc_degradation(t, y, strc_present=False),
            (rescue_day_mouse, t_max), [bundle[idx_rescue], survival_at_rescue * (1 - RESCUE_EFFICIENCY)],
            t_eval=np.arange(rescue_day_mouse, t_max + dt, dt), method='RK45', rtol=1e-6
        )
        final_survival = rescued + sol2.y[1][-1]

    human_deadline_days = (rescue_day_mouse or t_max) * HUMAN_MOUSE_SCALE

    return {
        "rescue_day_mouse": rescue_day_mouse,
        "human_equivalent_days": human_deadline_days,
        "human_equivalent_months": round(human_deadline_days / 30.4, 1),
        "final_ohc_survival_fraction": round(float(final_survival), 3),
        "viable_for_partial_rescue": bool(final_survival >= MIN_OHC_FOR_RESCUE),
        "viable_for_good_outcome": bool(final_survival >= MIN_OHC_FOR_GOOD_OUTCOME),
        "t": sol.t.tolist(),
        "bundle_integrity": bundle.tolist(),
        "ohc_survival_raw": survival.tolist(),
    }


# ─── Run Analysis ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Sweep intervention timing from P5 to P90 (mouse)
    rescue_days = [None, 5, 10, 15, 20, 30, 45, 60, 90]

    results = []
    print(f"\n{'Rescue day':>12} {'Human equiv':>14} {'Final OHC%':>12} {'Partial rescue':>16} {'Good outcome':>14}")
    print("-" * 72)

    for day in rescue_days:
        r = run_scenario(rescue_day_mouse=day)
        results.append(r)
        label = f"P{day}" if day else "untreated"
        human_mo = r["human_equivalent_months"] if day else "—"
        ohc_pct = f"{r['final_ohc_survival_fraction']*100:.1f}%"
        partial = "YES" if r["viable_for_partial_rescue"] else "NO"
        good = "YES" if r["viable_for_good_outcome"] else "NO"
        print(f"{label:>12} {str(human_mo)+' mo':>14} {ohc_pct:>12} {partial:>16} {good:>14}")

    # Find intervention deadline (last day with good outcome)
    good_days = [r["rescue_day_mouse"] for r in results
                 if r["rescue_day_mouse"] and r["viable_for_good_outcome"]]
    partial_days = [r["rescue_day_mouse"] for r in results
                    if r["rescue_day_mouse"] and r["viable_for_partial_rescue"]]

    deadline_good_mouse = max(good_days) if good_days else None
    deadline_partial_mouse = max(partial_days) if partial_days else None

    print(f"\n{'='*72}")
    if deadline_good_mouse:
        hm = round(deadline_good_mouse * HUMAN_MOUSE_SCALE / 30.4, 1)
        print(f"Good outcome deadline: mouse P{deadline_good_mouse} → human ~{hm} months postnatal")
    if deadline_partial_mouse:
        hm = round(deadline_partial_mouse * HUMAN_MOUSE_SCALE / 30.4, 1)
        print(f"Partial rescue deadline: mouse P{deadline_partial_mouse} → human ~{hm} months postnatal")

    print(f"\nClinical implication:")
    print(f"  - Newborn hearing screening (UNHS) detects DFNB16 at birth")
    print(f"  - Genetic confirmation: ~1-3 months")
    print(f"  - Treatment must begin before ~{round((deadline_good_mouse or 30)*HUMAN_MOUSE_SCALE/30.4, 0):.0f} months")
    print(f"  - This is compatible with standard newborn screening → diagnosis → IND pathway")

    # Save results (without time series for brevity)
    summary = []
    for r in results:
        s = {k: v for k, v in r.items() if k not in ("t", "bundle_integrity", "ohc_survival_raw")}
        summary.append(s)

    with open("therapeutic_window_results.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nFull results → therapeutic_window_results.json")
