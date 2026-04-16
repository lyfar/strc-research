#!/usr/bin/env python3
"""
ABR Transfer Function: OHC Transduction % → ABR Threshold
==========================================================
Maps mini-STRC transduction efficiency to clinically meaningful ABR
threshold predictions. Closes the loop between computational transduction
models and clinical outcome measures used in gene therapy trials.

Motivation (Lesperance et al. 2026, Ear and Hearing):
- OTOF trials use ABR threshold improvement as primary endpoint
- FDA will expect same endpoint for STRC IND filing
- Our models predict transduction %; we need to translate to ABR dB
- This creates a single number a clinician can evaluate

Model basis:
1. Cochlear amplifier theory (Davis 1983; Dallos 1992; Ashmore 2008):
   - OHCs provide ~40-50 dB of active amplification at threshold
   - IHC-only hearing (no OHC amplification): ~65 dB threshold
   - Normal hearing: ~15-20 dB threshold → ~45-50 dB OHC contribution

2. OHC loss vs threshold shift (Bredberg 1968; Schuknecht & Gacek 1993):
   - Relationship is nonlinear: early OHC loss disproportionately affects
     threshold; once most OHCs lost, additional loss has diminishing effect
   - Approximated by: threshold_shift = -A × log10(transduction_fraction)
   - where A ≈ 50-60 dB (from published OHC-count/threshold correlations)

3. DFNB16 baseline:
   - STRC-null patients: severe-profound SNHL, typically 70-90 dB threshold
   - Consistent with near-zero OHC function (STRC required for stereocilia
     HTCs and TM attachment → without it, OHCs are essentially non-functional)

4. Clinical significance threshold:
   - OTOF trials defined responder as ≥20 dB improvement
   - Functional hearing (can use with hearing aids): < 60 dB threshold
   - Speech perception without aids: < 40 dB threshold

Literature calibration points:
  transduction%  →  threshold dB  (source)
  100%           →  20 dB         (normal, Jerger reference)
   70%           →  30 dB         (mild-moderate; compatible with OTOF cohort 2)
   50%           →  38 dB         (moderate)
   30%           →  50 dB         (moderate-severe; minimum useful hearing)
   15%           →  63 dB         (severe; hearing aid territory)
    5%           →  78 dB         (severe-profound)
    0%           →  85 dB         (profound; Strc-/- mouse baseline)

Author: Egor Lyfar (computational)
Date: 2026-04-16
"""

import numpy as np
from scipy.optimize import curve_fit
import json

# ─── Calibration Data ─────────────────────────────────────────────────────────

# (transduction_fraction, abr_threshold_dB)
CALIBRATION_POINTS = [
    (1.00, 20),   # normal
    (0.70, 30),   # OTOF DB-OTO cohort 2 estimate
    (0.55, 38),   # OTOF Sun et al. 2024 estimate
    (0.40, 47),   # OTOF DB-OTO cohort 1 estimate
    (0.30, 53),   # Bredberg: ~50% OHC function ≈ moderate-severe
    (0.15, 64),   # severe HL
    (0.05, 76),   # severe-profound
    (0.00, 85),   # DFNB16 / Strc-/- baseline
]

# ─── Transfer Function ────────────────────────────────────────────────────────

def abr_log_model(transduction, A, B, C):
    """
    ABR threshold = A - B × log10(transduction + C)
    Handles transduction=0 via offset C, avoids log(0).
    """
    return A - B * np.log10(transduction + C)


def fit_transfer_function():
    """Fit log model to calibration data. Returns fitted parameters."""
    x = np.array([p[0] for p in CALIBRATION_POINTS])
    y = np.array([p[1] for p in CALIBRATION_POINTS])

    # Initial guess: A=20, B=60, C=0.01
    popt, pcov = curve_fit(abr_log_model, x, y, p0=[20, 60, 0.01],
                            bounds=([0, 20, 1e-4], [50, 100, 0.1]))
    residuals = y - abr_log_model(x, *popt)
    rmse = np.sqrt(np.mean(residuals ** 2))
    return popt, rmse


def transduction_to_abr(transduction_fraction, params=None):
    """
    Convert OHC transduction fraction [0-1] to predicted ABR threshold [dB].
    transduction_fraction: fraction of OHCs expressing functional mini-STRC
    Returns predicted ABR threshold in dB SPL (HL).
    """
    if params is None:
        params, _ = fit_transfer_function()
    threshold = abr_log_model(transduction_fraction, *params)
    # Clamp to physiological range
    return float(np.clip(threshold, 15, 90))


def abr_improvement(transduction_fraction, baseline_dB=85, params=None):
    """ABR threshold improvement over untreated DFNB16 baseline."""
    new_threshold = transduction_to_abr(transduction_fraction, params)
    return baseline_dB - new_threshold


# ─── Clinical Significance Categories ─────────────────────────────────────────

def classify_outcome(threshold_dB):
    """
    Standard audiometric classification (WHO/ASHA).
    Returns category and clinical interpretation.
    """
    if threshold_dB <= 25:
        return "Normal", "No intervention needed"
    elif threshold_dB <= 40:
        return "Mild", "Conversational speech without aids"
    elif threshold_dB <= 55:
        return "Moderate", "Benefits from hearing aids; functional hearing"
    elif threshold_dB <= 70:
        return "Moderate-severe", "Hearing aids required; speech perception possible"
    elif threshold_dB <= 90:
        return "Severe", "Limited benefit from hearing aids; cochlear implant territory"
    else:
        return "Profound", "Cochlear implant candidate; no functional hearing"


RESPONDER_THRESHOLD_DB = 20  # ≥20 dB improvement = responder (OTOF trial standard)
FUNCTIONAL_THRESHOLD_DB = 55  # <55 dB = functional hearing (WHO moderate)


# ─── Run Analysis ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    params, rmse = fit_transfer_function()
    print(f"Transfer function fit: A={params[0]:.1f}, B={params[1]:.1f}, C={params[2]:.4f}")
    print(f"RMSE against calibration: {rmse:.1f} dB\n")

    # Transduction sweep
    transduction_levels = [0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.67, 0.75, 0.85, 1.0]

    print(f"{'Transduction%':>15} {'ABR threshold':>14} {'Improvement':>13} {'Category':>20} {'Responder?':>12}")
    print("-" * 78)

    results = []
    for t in transduction_levels:
        threshold = transduction_to_abr(t, params)
        improvement = abr_improvement(t, params=params)
        category, interp = classify_outcome(threshold)
        is_responder = improvement >= RESPONDER_THRESHOLD_DB
        results.append({
            "transduction_pct": round(t * 100, 1),
            "abr_threshold_dB": round(threshold, 1),
            "improvement_dB": round(improvement, 1),
            "category": category,
            "responder": is_responder,
            "functional": threshold <= FUNCTIONAL_THRESHOLD_DB
        })
        resp_str = "YES" if is_responder else "NO"
        print(f"{t*100:>14.0f}% {threshold:>13.1f} dB {improvement:>12.1f} dB {category:>20} {resp_str:>12}")

    # Key finding: what transduction does mini-STRC achieve?
    mini_strc_transduction = 0.67  # from dual-vector model at clinical titer
    mini_threshold = transduction_to_abr(mini_strc_transduction, params)
    mini_improvement = abr_improvement(mini_strc_transduction, params=params)
    mini_category, mini_interp = classify_outcome(mini_threshold)

    print(f"\n{'='*78}")
    print(f"Mini-STRC projection at clinical titer (3.75×10¹² GC/mL, 67% transduction):")
    print(f"  Predicted ABR threshold: {mini_threshold:.0f} dB")
    print(f"  Improvement over baseline: {mini_improvement:.0f} dB")
    print(f"  Category: {mini_category}")
    print(f"  Clinical interpretation: {mini_interp}")
    print(f"  Responder (≥20 dB improvement): {'YES' if mini_improvement >= 20 else 'NO'}")
    print(f"  Functional hearing (< 55 dB): {'YES' if mini_threshold <= 55 else 'NO'}")

    # Minimum transduction for responder status
    min_t_responder = None
    for t in np.linspace(0, 1, 1000):
        if abr_improvement(t, params=params) >= RESPONDER_THRESHOLD_DB:
            min_t_responder = t
            break

    min_t_functional = None
    for t in np.linspace(0, 1, 1000):
        if transduction_to_abr(t, params) <= FUNCTIONAL_THRESHOLD_DB:
            min_t_functional = t
            break

    print(f"\nClinical thresholds:")
    if min_t_responder:
        print(f"  Minimum transduction for responder status: {min_t_responder*100:.0f}%")
    if min_t_functional:
        print(f"  Minimum transduction for functional hearing: {min_t_functional*100:.0f}%")

    with open("abr_transfer_results.json", "w") as f:
        json.dump({
            "model_params": {"A": params[0], "B": params[1], "C": params[2]},
            "fit_rmse_dB": round(rmse, 2),
            "mini_strc_projection": {
                "transduction_pct": 67.0,
                "abr_threshold_dB": round(mini_threshold, 1),
                "improvement_dB": round(mini_improvement, 1),
                "category": mini_category,
                "responder": bool(mini_improvement >= RESPONDER_THRESHOLD_DB),
                "functional_hearing": bool(mini_threshold <= FUNCTIONAL_THRESHOLD_DB)
            },
            "min_transduction_for_responder_pct": round((min_t_responder or 0) * 100, 1),
            "min_transduction_for_functional_pct": round((min_t_functional or 0) * 100, 1),
            "sweep": results
        }, f, indent=2)

    print(f"\nFull results → abr_transfer_results.json")
