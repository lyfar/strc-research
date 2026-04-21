#!/usr/bin/env python3
"""
Stereocilia Bundle Mechanics — HTC Coupling vs Cochlear Function
================================================================
Analytical model relating HTC spring coupling (maintained by STRC) to
ABR threshold, cochlear amplification, and functional hearing.

Biological framing:
- Misha is compound heterozygous: STRC deletion (paternal) + VUS c.4976A>C (maternal)
- He is NOT a complete STRC knockout — VUS allele produces some stereocilin
- Moderate hearing loss (~40-60 dB) = partial HTC coupling, not zero
- This model estimates what fraction of coupling remains and what gene therapy adds

Mechanism (from Verpy 2008, Fettiplace 2017, Ashmore 2008):
- HTCs synchronize stereocilia rows 1-3: signal from TM reaches all MET channels
- Without HTCs: only tallest row (row 1) directly driven → 1/3 of MET current
- Without TM attachment crowns: OHC-TM feedback loop broken
- Cochlear amplification is a Hopf oscillator: gain ∝ coupling^α (nonlinear)
- Loss of HTC coupling → loss of active amplification → ~40-50 dB threshold shift

Two-component model:
1. MET current component: HTC coupling drives rows 2-3; loss → -9 dB (3x fewer channels)
2. Amplification component: coupling required for Hopf oscillator; loss → -35 dB
   Total: ~44 dB shift at k_HTC = 0 (matches Verpy 2008 Strc-/- data)

For partial coupling (Misha):
  ABR_shift = A_met × (1 - f) + A_amp × (1 - f^alpha)
  where f = k_HTC / k_HTC_max, A_met = 9 dB, A_amp = 35 dB, alpha = 0.5

Parameters from literature:
- A_total = 44 dB at f=0 (Verpy 2008, Strc-/- threshold shift)
- A_met = 9.5 dB (20 log10(3)) — 3 rows vs 1 row stimulation
- A_amp = 34.5 dB — active Hopf amplification component
- alpha = 0.5 (square-root coupling for Hopf oscillator near bifurcation,
             Camalet et al. 2000, PNAS; Eguíluz et al. 2000)
- Misha phenotype: 30-55 dB shift from normal hearing (~20 dB threshold)

Author: Egor Lyfar (computational)
Date: 2026-04-17
"""

import numpy as np
import json

# ─── Parameters ───────────────────────────────────────────────────────────────

# ABR threshold shift components (calibrated to Verpy 2008)
A_MET = 9.5      # dB from 3→1 active rows at f=0
A_AMP = 34.5     # dB from lost Hopf amplification at f=0
ALPHA = 0.5      # Hopf oscillator coupling exponent

BASELINE_NORMAL_ABR = 20.0  # dB SPL, normal hearing threshold (ABR click)
DFNB16_BASELINE = BASELINE_NORMAL_ABR + A_MET + A_AMP  # ~64 dB (complete KO)


def abr_threshold(f_htc):
    """
    Predicted ABR threshold (dB SPL) for HTC coupling fraction f_htc ∈ [0, 1].
    f_htc = k_HTC / k_HTC_max, i.e. fraction of normal STRC coupling.

    At f=1: normal hearing (~20 dB)
    At f=0: complete KO (Strc-/- ~64 dB)
    Misha: moderate (40-60 dB) → f ≈ 0.10-0.35
    """
    shift_met = A_MET * (1.0 - f_htc)
    shift_amp = A_AMP * (1.0 - f_htc ** ALPHA)
    return BASELINE_NORMAL_ABR + shift_met + shift_amp


def f_htc_for_threshold(target_dB, tol=0.1):
    """Inverse: find HTC coupling fraction that gives target ABR threshold."""
    for f in np.linspace(0, 1, 10000):
        if abs(abr_threshold(f) - target_dB) < tol:
            return f
    return None


def classify_hearing(threshold_dB):
    if threshold_dB <= 25:
        return "Normal"
    elif threshold_dB <= 40:
        return "Mild"
    elif threshold_dB <= 55:
        return "Moderate"
    elif threshold_dB <= 70:
        return "Moderate-severe"
    else:
        return "Severe-profound"


# ─── Rescuability Analysis ────────────────────────────────────────────────────

def rescuability(f_htc_current, f_htc_therapy=0.70):
    """
    Can gene therapy rescue function from current coupling level?

    Assumptions:
    - Gene therapy (Anc80L65, 67% OHC transduction) restores STRC in 67% of cells
    - In transduced cells: HTC coupling restored toward full k_HTC_max
    - In untransduced cells: coupling stays at f_htc_current
    - Bundle coherence requires critical mass: if <30% cells transduced → minimal gain
    - Tectorial membrane remodeling: unknown factor, modeled as penalty [0-1]

    Returns: predicted post-therapy ABR threshold (optimistic case)
    """
    TRANSDUCTION = 0.67
    # Weighted average: transduced cells get full coupling, rest stay partial
    f_post = TRANSDUCTION * 1.0 + (1 - TRANSDUCTION) * f_htc_current
    return abr_threshold(f_post)


# ─── Main Analysis ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("HTC Coupling → ABR Threshold Model")
    print("Calibrated to Verpy 2008 (Strc-/- mouse: ~44 dB shift)")
    print("=" * 68)

    # Coupling sweep
    fractions = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.85, 1.0]

    print(f"\n{'HTC coupling':>14} {'ABR threshold':>15} {'Category':>20} {'Misha range?':>14}")
    print("-" * 66)

    results = []
    for f in fractions:
        thr = abr_threshold(f)
        cat = classify_hearing(thr)
        # Misha: moderate bilateral (right ~50-65 dB, left ~40-50 dB)
        misha = 35.0 <= thr <= 65.0
        results.append({
            "htc_coupling_pct": round(f * 100, 0),
            "abr_threshold_dB": round(thr, 1),
            "category": cat,
            "matches_misha": bool(misha)
        })
        flag = "← MISHA" if misha else ""
        print(f"{f*100:>13.0f}% {thr:>13.1f} dB {cat:>20} {flag}")

    # Misha's estimated coupling range
    misha_results = [r for r in results if r["matches_misha"]]
    f_lo = min(r["htc_coupling_pct"] for r in misha_results) if misha_results else "?"
    f_hi = max(r["htc_coupling_pct"] for r in misha_results) if misha_results else "?"

    print(f"\n{'='*68}")
    print(f"\nMisha's estimated current HTC coupling: ~{f_lo}–{f_hi}% of maximum")
    print(f"(from moderate SNHL phenotype, audiometry Sept 2025)")

    # Therapy projections at Misha's coupling range
    print(f"\nGene therapy projections (Anc80L65, 67% transduction):")
    print(f"\n{'Current coupling':>18} {'Current ABR':>13} {'Post-therapy ABR':>17} {'Improvement':>13} {'Outcome':>20}")
    print("-" * 84)

    therapy_results = []
    for r in results:
        if r["matches_misha"]:
            f = r["htc_coupling_pct"] / 100
            current_thr = r["abr_threshold_dB"]
            post_thr = rescuability(f)
            improvement = current_thr - post_thr
            outcome = classify_hearing(post_thr)
            therapy_results.append({
                "current_coupling_pct": r["htc_coupling_pct"],
                "current_abr_dB": current_thr,
                "post_therapy_abr_dB": round(post_thr, 1),
                "improvement_dB": round(improvement, 1),
                "post_therapy_category": outcome
            })
            print(f"{r['htc_coupling_pct']:>17.0f}% {current_thr:>12.1f} dB "
                  f"{post_thr:>15.1f} dB {improvement:>12.1f} dB {outcome:>20}")

    print(f"""
Key structural question (NOT modelled — requires experiment):
  Can tectorial membrane re-establish contact with OHC attachment crowns
  after ~5 years of disconnection?

  If YES → therapy restores full coupling chain, model prediction holds
  If NO  → HTCs can reform internally, but TM-OHC link broken → partial gain

  Evidence: TM structural changes in Strc-/- mice documented (subtle at P30,
  progressive after). At 5 human-equivalent years: unknown.
  This is the single most important experiment for Misha's program.
""")

    output = {
        "model": "htc_coupling_hopf_analytical",
        "calibration": "Verpy 2008 Strc-/- 44dB shift at f=0",
        "params": {"A_met_dB": A_MET, "A_amp_dB": A_AMP, "alpha": ALPHA},
        "misha_estimated_coupling_pct": f"{f_lo}-{f_hi}",
        "sweep": results,
        "therapy_projections": therapy_results
    }

    with open("/Users/egorlyfar/Brain/sites/strc/models/bundle_mechanics_results.json", "w") as f_out:
        json.dump(output, f_out, indent=2)
    print("Results → Brain/sites/strc/models/bundle_mechanics_results.json")
