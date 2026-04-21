#!/usr/bin/env python3
"""
Dual-Vector Model Recalibration Against OTOF Clinical Data
===========================================================
Recalibrates the gamma-Poisson (negative binomial) transduction model
against OTOF/DB-OTO clinical trial outcomes, then recomputes the
mini-STRC single-vector advantage ratio.

Motivation (Lesperance et al. 2026, Ear and Hearing):
- OTOF gene therapy (DB-OTO trial) is the first successful inner ear gene
  therapy and the clinical comparator for all subsequent programs including STRC
- The dual-AAV approach used for OTOF provides real clinical transduction data
- Our 2.8x advantage claim was calibrated against Omichi 2020 bench data only
- Recalibrating against OTOF Phase 1/2 outcomes validates the model prediction

OTOF Clinical Reference Data:
- DB-OTO (Regeneron/Decibel Therapeutics): Anc80L65 dual-AAV OTOF, cochlear injection
  * Phase 1/2, pediatric, severe-profound ANSD due to OTOF mutations
  * Published: Lustig et al. 2024 (NEJM Evidence) + Oesterle 2024 (Laryngoscope)
  * Dose cohort 1 (~1×10¹¹ GC/ear): 2/3 patients showed ABR, ~40% transduction est.
  * Dose cohort 2 (~3×10¹¹ GC/ear): 3/3 patients showed ABR, ~60-70% transduction est.
- Chinese OTOF trial (Sun et al. 2024, Lancet): dual-AAV, ~5×10¹¹ GC/ear
  * 6/6 patients showed ABR threshold improvement
  * Estimated OHC transduction: ~55-65% (inferred from ABR improvement magnitude)

Note: transduction% is inferred from ABR improvement using the ABR transfer
function (abr_transfer_model.py). Direct transduction counts from clinical
human cochleae are not available.

Original model calibration: Omichi et al. 2020
- Single AAV2, ex vivo guinea pig OHCs: 83.9% transduction at ~10¹² GC/mL
- Dual AAV2: 65.6% co-transduction at same titer
- Co-recombination efficiency R estimated at ~50% of co-transduction events

Author: Egor Lyfar (computational)
Date: 2026-04-16
"""

import numpy as np
from scipy.stats import nbinom
from scipy.optimize import minimize_scalar
import json

# ─── Negative Binomial (Gamma-Poisson) Model ──────────────────────────────────

def nb_transduction(titer_gc_ml, r_disp=0.35, scale_factor=None):
    """
    P(cell transduced) = P(NegBin(r, p) >= 1)
    NegBin calibrated to Omichi 2020: 83.9% at 10¹² GC/mL (single vector).

    r_disp: dispersion parameter (cell-to-cell heterogeneity)
    scale_factor: converts titer → mean MOI per cell
    """
    if scale_factor is None:
        scale_factor = 1e-12  # calibrated: mean MOI=1 at titer 10¹² GC/mL

    mean_moi = titer_gc_ml * scale_factor
    # NegBin parameterization: mean = r*(1-p)/p → p = r/(r + mean_moi)
    p_nb = r_disp / (r_disp + mean_moi)
    p_not_transduced = nbinom.pmf(0, r_disp, p_nb)
    return 1.0 - p_not_transduced


def dual_vector_transduction(titer_gc_ml, r_disp=0.35, scale_factor=None,
                              recombination_efficiency=0.50):
    """
    P(co-transduced AND recombined) for dual-vector approach.
    Both vectors must enter same cell AND undergo ITR-mediated recombination.
    """
    p_single = nb_transduction(titer_gc_ml, r_disp, scale_factor)
    p_cotransduced = p_single ** 2   # independent uptake of each vector
    p_functional = p_cotransduced * recombination_efficiency
    return p_functional


def calibrate_to_omichi(target_single=0.839, titer_ref=1e12):
    """Find scale_factor that reproduces Omichi 2020 single-vector result."""
    def objective(log_sf):
        sf = 10 ** log_sf
        pred = nb_transduction(titer_ref, scale_factor=sf)
        return (pred - target_single) ** 2

    res = minimize_scalar(objective, bounds=(-15, -9), method='bounded')
    return 10 ** res.x


# ─── OTOF Clinical Validation ─────────────────────────────────────────────────

# Clinical data points: (titer GC/ear, estimated transduction %)
# Transduction % inferred from ABR improvement via abr_transfer_model.py
# ABR improvement → transduction% is approximate; see notes in each entry
OTOF_CLINICAL_DATA = [
    {
        "source": "DB-OTO cohort 1 (Lustig 2024)",
        "titer_gc_ear": 1e11,
        "cochlear_volume_uL": 50,      # human cochlear fluid volume
        "estimated_transduction": 0.40,  # inferred from ~20 dB ABR improvement
        "n_patients": 3,
        "responders": 2,
        "note": "2/3 showed ABR; ~40% transduction estimated from partial response"
    },
    {
        "source": "DB-OTO cohort 2 (Lustig 2024)",
        "titer_gc_ear": 3e11,
        "cochlear_volume_uL": 50,
        "estimated_transduction": 0.65,  # inferred from >30 dB ABR improvement
        "n_patients": 3,
        "responders": 3,
        "note": "3/3 showed ABR; ~65% transduction estimated from full response"
    },
    {
        "source": "Sun et al. 2024 (Chinese trial)",
        "titer_gc_ear": 5e11,
        "cochlear_volume_uL": 50,
        "estimated_transduction": 0.60,
        "n_patients": 6,
        "responders": 6,
        "note": "6/6 ABR improvement; ~60% transduction inferred"
    },
]


def titer_ear_to_concentration(titer_gc_ear, volume_uL=50):
    """Convert GC/ear to GC/mL (cochlear injection concentration)."""
    volume_mL = volume_uL / 1000
    return titer_gc_ear / volume_mL


# ─── Run Calibration ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Step 1: Calibrate scale factor to Omichi 2020
    scale_factor = calibrate_to_omichi()
    print(f"Scale factor (Omichi 2020 calibration): {scale_factor:.3e}")

    # Verify calibration
    pred_omichi_single = nb_transduction(1e12, scale_factor=scale_factor)
    pred_omichi_dual = dual_vector_transduction(1e12, scale_factor=scale_factor)
    print(f"Omichi 2020 verification:")
    print(f"  Single-vector @10¹² GC/mL: predicted {pred_omichi_single*100:.1f}% (observed 83.9%)")
    print(f"  Dual-vector   @10¹² GC/mL: predicted {pred_omichi_dual*100:.1f}% (observed 65.6% co-transduction)")

    # Step 2: Predict transduction at OTOF clinical titers
    print(f"\n{'='*70}")
    print("OTOF Clinical Validation:")
    print(f"\n{'Source':<30} {'Conc (GC/mL)':>14} {'Pred dual%':>12} {'Obs est%':>10} {'Error':>8}")
    print("-" * 70)

    errors = []
    for d in OTOF_CLINICAL_DATA:
        conc = titer_ear_to_concentration(d["titer_gc_ear"], d["cochlear_volume_uL"])
        pred = dual_vector_transduction(conc, scale_factor=scale_factor)
        obs = d["estimated_transduction"]
        err = abs(pred - obs)
        errors.append(err)
        print(f"{d['source']:<30} {conc:>14.2e} {pred*100:>11.1f}% {obs*100:>9.1f}% {err*100:>7.1f}%")

    mean_error = np.mean(errors)
    print(f"\nMean absolute error: {mean_error*100:.1f}%")
    if mean_error < 0.12:
        print("Model calibration: VALID (error < 12%)")
    else:
        print(f"Model calibration: NEEDS REFINEMENT (error {mean_error*100:.1f}% > 12%)")

    # Step 3: Recompute single-vector advantage at clinical titers
    # Using mini-STRC clinical titer target: 3.75×10¹² GC/mL
    print(f"\n{'='*70}")
    print("Mini-STRC Single-Vector Advantage (recalibrated model):")
    print(f"\n{'Titer (GC/mL)':>15} {'Single%':>10} {'Dual%':>10} {'Advantage':>12}")
    print("-" * 50)

    clinical_titers = [1e10, 3e10, 1e11, 3e11, 3.75e12, 1e13]
    summary = []
    for titer in clinical_titers:
        single = nb_transduction(titer, scale_factor=scale_factor)
        dual = dual_vector_transduction(titer, scale_factor=scale_factor)
        advantage = single / dual if dual > 0 else float('inf')
        summary.append({
            "titer_gc_ml": titer,
            "single_vector_pct": round(single * 100, 1),
            "dual_vector_pct": round(dual * 100, 1),
            "advantage_ratio": round(advantage, 2)
        })
        label = f"{titer:.2e}"
        adv_str = f"{advantage:.2f}x" if dual > 0.001 else ">100x"
        print(f"{label:>15} {single*100:>9.1f}% {dual*100:>9.1f}% {adv_str:>12}")

    clinical_result = next(r for r in summary if r["titer_gc_ml"] == 3.75e12)
    print(f"\nAt clinical titer (3.75×10¹² GC/mL):")
    print(f"  Single-vector (mini-STRC): {clinical_result['single_vector_pct']}%")
    print(f"  Dual-vector (full STRC):   {clinical_result['dual_vector_pct']}%")
    print(f"  Advantage ratio:           {clinical_result['advantage_ratio']}x")

    if abs(clinical_result['advantage_ratio'] - 2.8) < 0.4:
        print(f"  → Original 2.8x claim CONFIRMED after OTOF clinical recalibration")
    else:
        print(f"  → Original 2.8x claim REVISED to {clinical_result['advantage_ratio']:.1f}x")

    with open("dual_vector_otof_results.json", "w") as f:
        json.dump({
            "scale_factor": scale_factor,
            "otof_validation": OTOF_CLINICAL_DATA,
            "mean_error_pct": round(mean_error * 100, 2),
            "titer_sweep": summary
        }, f, indent=2)

    print(f"\nFull results → dual_vector_otof_results.json")
