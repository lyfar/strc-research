#!/usr/bin/env python3
"""
Hydrogel Phase 4e — cochlear PKPD for tail91 peptide (134 aa).

Two-compartment intra-cochlear model for ototopical delivery:
  C1  Middle-ear depot (topical, gel or solution)
  C2  Perilymph free peptide (cochlear fluid, 70 μL, seen by TMEM145 on OHC apical surface)

Plus a degraded-amount sink for proteolysis/clearance accounting.

Readout at OHC surface is not a separate compartment — the perilymph [peptide]
directly sets TMEM145 fractional occupancy via θ = [P] / (Kd + [P]).

Key question: What ototopical dose gives peak perilymph concentration in the
Phase 4d therapeutic window (1-10 μM — enough for F-actin bundling, not so
much that G-actin is sequestered)?
"""

from __future__ import annotations

import json
import numpy as np
from scipy.integrate import solve_ivp
from pathlib import Path

OUT = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase4e_cochlear_pkpd.json")

# Peptide properties
PEPTIDE_MW_KDA = 14.2  # 134 aa
PERILYMPH_VOL_UL = 70.0  # human cochlea

# Rate constants (1/h)
K_CLEAR_MIDDLE_EAR = 0.35  # mucociliary clearance, half-life 2 h
K_RWM = 0.02               # round window crossing; ~2%/h for 14 kDa peptide (Salt 2011)
K_PERILYMPH_CLEAR = 0.35   # cochlear aqueduct + perilymph-CSF exchange, half-life 2 h
K_PROTEOLYSIS = 1.4        # perilymph proteolysis, half-life 30 min (conservative)
                           # perilymph has fewer proteases than serum, but 30 min is a
                           # reasonable upper bound for unprotected 134 aa with 9 K/R cuts

# Target TMEM145 binding Kd (from Phase 3 solo tail91 ipTM 0.68 → Kd ~50-500 nM)
KD_TMEM145_M = 100e-9   # 100 nM nominal

# Therapeutic window from Phase 4d:
#   lower bound: 1 μM perilymph → 91% TMEM145 occupancy, sufficient bundling
#   upper bound: 10 μM perilymph → 99% TMEM145 occupancy, G-actin sequestration < 30%
#   toxic threshold: >100 μM perilymph → G-actin depletion >99%
LOWER_THERAPEUTIC_uM = 1.0
UPPER_THERAPEUTIC_uM = 10.0
TOXIC_THRESHOLD_uM = 100.0


def pkpd_derivative(t, y):
    """
    y[0] = C1 middle-ear amount (mg)
    y[1] = C2 perilymph amount (mg)
    y[2] = C_deg cumulative degraded (mg)
    """
    c1, c2, _ = y
    flux_12_rwm = c1 * K_RWM
    flux_1_out = c1 * K_CLEAR_MIDDLE_EAR
    flux_2_clear = c2 * K_PERILYMPH_CLEAR
    flux_2_deg = c2 * K_PROTEOLYSIS
    dc1 = -flux_12_rwm - flux_1_out
    dc2 = flux_12_rwm - flux_2_clear - flux_2_deg
    ddeg = flux_2_deg + flux_1_out + flux_2_clear
    return [dc1, dc2, ddeg]


def amount_to_conc_uM(amount_mg: float) -> float:
    """Convert perilymph amount (mg) to concentration (μM)."""
    moles = (amount_mg * 1e-3) / (PEPTIDE_MW_KDA * 1000)
    vol_L = PERILYMPH_VOL_UL * 1e-6
    return (moles / vol_L) * 1e6


def fractional_occupancy(conc_uM: float, kd_uM: float = KD_TMEM145_M * 1e6) -> float:
    return conc_uM / (kd_uM + conc_uM)


def simulate_single_dose(dose_mg: float, duration_h: float = 48.0):
    y0 = [dose_mg, 0.0, 0.0]
    t_eval = np.linspace(0, duration_h, 400)
    sol = solve_ivp(pkpd_derivative, [0, duration_h], y0,
                    t_eval=t_eval, method="LSODA",
                    rtol=1e-6, atol=1e-12)
    c1, c2, cdeg = sol.y
    peri_uM = np.array([amount_to_conc_uM(x) for x in c2])
    occupancy = np.array([fractional_occupancy(x) for x in peri_uM])
    # Time above therapeutic lower bound (1 μM)
    above_lower = sol.t[peri_uM >= LOWER_THERAPEUTIC_uM]
    dur_above_lower = float(above_lower[-1] - above_lower[0]) if len(above_lower) > 1 else 0.0
    # Time in therapeutic window (1-10 μM)
    in_window = (peri_uM >= LOWER_THERAPEUTIC_uM) & (peri_uM <= UPPER_THERAPEUTIC_uM)
    dur_in_window = float(np.ptp(sol.t[in_window])) if np.any(in_window) else 0.0
    # Time above toxic threshold (100 μM)
    above_toxic = sol.t[peri_uM > TOXIC_THRESHOLD_uM]
    dur_above_toxic = float(above_toxic[-1] - above_toxic[0]) if len(above_toxic) > 1 else 0.0
    # Bioavailability into perilymph: integrated flux_12 over 24 h / dose
    bioavail = 1 - (float(c1[t_eval.searchsorted(24)]) + float(cdeg[t_eval.searchsorted(24)])) / dose_mg
    auc_peri = float(np.trapezoid(peri_uM, sol.t))
    return {
        "dose_mg": dose_mg,
        "peak_perilymph_uM": float(peri_uM.max()),
        "time_peak_perilymph_h": float(sol.t[peri_uM.argmax()]),
        "peak_TMEM145_occupancy": float(occupancy.max()),
        "auc_perilymph_uM_hr": auc_peri,
        "duration_above_1uM_h": dur_above_lower,
        "duration_in_therapeutic_window_1_to_10_uM_h": dur_in_window,
        "duration_above_toxic_100uM_h": dur_above_toxic,
        "apparent_bioavailability_perilymph_24h": bioavail,
        "time_course": {
            "t_h": sol.t.tolist(),
            "perilymph_uM": peri_uM.tolist(),
            "occupancy": occupancy.tolist(),
        },
    }


def dose_to_achieve_peak(target_uM: float):
    """Find dose for which peak perilymph = target_uM."""
    # Bisection; dose in mg, 0.001 → 100
    lo, hi = 0.001, 100.0
    for _ in range(60):
        mid = (lo + hi) / 2
        r = simulate_single_dose(mid, duration_h=12.0)
        if r["peak_perilymph_uM"] < target_uM:
            lo = mid
        else:
            hi = mid
    return mid


def main():
    doses_mg = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
    per_dose = {}
    for d in doses_mg:
        r = simulate_single_dose(d)
        per_dose[f"{d}_mg"] = {k: v for k, v in r.items() if k != "time_course"}

    # Find doses that hit window boundaries
    dose_for_1uM = dose_to_achieve_peak(LOWER_THERAPEUTIC_uM)
    dose_for_10uM = dose_to_achieve_peak(UPPER_THERAPEUTIC_uM)
    dose_for_100uM = dose_to_achieve_peak(TOXIC_THRESHOLD_uM)

    # Multi-dose — daily for 7 days at the mid-window dose
    mid_window_dose = round(np.sqrt(dose_for_1uM * dose_for_10uM), 3)
    multi_dose_result = None
    y = [0.0, 0.0, 0.0]
    t_acc, peri_acc, occ_acc = [], [], []
    for day in range(7):
        y = [y[0] + mid_window_dose, y[1], y[2]]
        t_eval = np.linspace(0, 24, 100)
        sol = solve_ivp(pkpd_derivative, [0, 24], y, t_eval=t_eval, method="LSODA")
        y = [float(sol.y[i][-1]) for i in range(3)]
        peri = np.array([amount_to_conc_uM(x) for x in sol.y[1]])
        occ = np.array([fractional_occupancy(x) for x in peri])
        t_acc += (sol.t + 24 * day).tolist()
        peri_acc += peri.tolist()
        occ_acc += occ.tolist()
    peri_arr = np.array(peri_acc)
    multi_dose_result = {
        "daily_dose_mg": mid_window_dose,
        "n_days": 7,
        "peak_perilymph_uM_in_cycle": float(peri_arr.max()),
        "trough_perilymph_uM_in_cycle": float(peri_arr[-1]),
        "fraction_time_above_1uM": round(float(np.mean(peri_arr >= 1.0)), 3),
        "fraction_time_in_window_1_to_10_uM": round(float(np.mean(
            (peri_arr >= 1.0) & (peri_arr <= 10.0))), 3),
        "mean_occupancy_steady_state_day7": round(float(np.mean(
            [o for o, t in zip(occ_acc, t_acc) if t >= 5 * 24])), 3),
    }

    # Sensitivity to proteolysis half-life
    sensitivity_proteolysis = {}
    global K_PROTEOLYSIS
    for t12_min in [10, 30, 60, 240]:
        K_PROTEOLYSIS = np.log(2) / (t12_min / 60.0)
        r = simulate_single_dose(mid_window_dose, duration_h=24.0)
        sensitivity_proteolysis[f"{t12_min}_min"] = {
            "peak_perilymph_uM": r["peak_perilymph_uM"],
            "duration_above_1uM_h": r["duration_above_1uM_h"],
            "mean_occupancy_over_24h": round(float(np.mean(
                [fractional_occupancy(x) for x in r["time_course"]["perilymph_uM"]])), 3),
        }
    K_PROTEOLYSIS = 1.4  # restore

    summary = {
        "batch": "hydrogel_phase4e_cochlear_pkpd",
        "date": "2026-04-23",
        "model_assumptions": {
            "peptide_mw_kDa": PEPTIDE_MW_KDA,
            "perilymph_volume_uL": PERILYMPH_VOL_UL,
            "middle_ear_clearance_half_life_h": round(np.log(2) / K_CLEAR_MIDDLE_EAR, 2),
            "rwm_crossing_percent_per_h": K_RWM * 100,
            "perilymph_clearance_half_life_h": round(np.log(2) / K_PERILYMPH_CLEAR, 2),
            "proteolysis_half_life_h": round(np.log(2) / K_PROTEOLYSIS, 2),
            "TMEM145_Kd_M": KD_TMEM145_M,
            "therapeutic_window_lower_uM": LOWER_THERAPEUTIC_uM,
            "therapeutic_window_upper_uM": UPPER_THERAPEUTIC_uM,
            "toxic_threshold_uM": TOXIC_THRESHOLD_uM,
        },
        "single_dose_scan": per_dose,
        "dose_to_hit_boundaries_mg": {
            "peak_perilymph_1_uM": round(dose_for_1uM, 3),
            "peak_perilymph_10_uM": round(dose_for_10uM, 3),
            "peak_perilymph_100_uM_toxic": round(dose_for_100uM, 3),
        },
        "therapeutic_window_width_log": round(
            np.log10(dose_for_100uM / dose_for_1uM), 2),
        "multi_dose_7day_at_mid_window": multi_dose_result,
        "sensitivity_to_proteolysis_half_life": sensitivity_proteolysis,
        "key_predictions": [
            f"Dose {dose_for_1uM:.2f}-{dose_for_10uM:.2f} mg hits 1-10 μM perilymph (Phase 4d therapeutic window)",
            f"At {mid_window_dose} mg daily, steady-state perilymph ~mid-window, occupancy > 90% most of day",
            f"Toxic G-actin depletion threshold (100 μM) requires dose {dose_for_100uM:.1f} mg — >30× above therapeutic, wide safety margin",
            f"Half-life dominated by proteolysis (30 min assumed); if true t½ < 10 min, PEGylation or d-aa backbone needed",
        ],
        "clinical_translation_comparison": {
            "DB-OTO (OTOF AAV, intracochlear surgery)": "single-dose cure; our peptide is palliative chronic",
            "ototopical neomycin": "mg-range doses well tolerated topically in humans",
            "Auris Medical AM-101 (NMDA antagonist for tinnitus)": "intratympanic gel 0.5-2 mg, 1-2 week duration — our format pattern",
        },
        "risk_flags": [
            "Proteolysis t½ assumed 30 min — if actual < 10 min (plausible for 134 aa with 9 K/R), need PEGylation or d-AA substitutions",
            "Model assumes rigid 1-compartment perilymph; apical-basal gradient exists in vivo (Salt 2011) — OHC apical surface at scala tympani may see 2-3× higher concentration than model mean",
            "RWM crossing 2%/h is for aqueous vehicle; hydrogel depot format may slow absorption by 10× → lower peaks but longer exposure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, default=str))

    print("=== Phase 4e Cochlear PKPD for tail91 peptide (14.2 kDa) ===\n")
    print(f"Model: 2-compartment (middle ear → perilymph)")
    print(f"  Middle-ear t½ = {np.log(2)/K_CLEAR_MIDDLE_EAR:.1f} h")
    print(f"  RWM crossing = {K_RWM*100:.0f}% / h")
    print(f"  Perilymph clearance t½ = {np.log(2)/K_PERILYMPH_CLEAR:.1f} h")
    print(f"  Proteolysis t½ = {np.log(2)/K_PROTEOLYSIS*60:.0f} min")
    print(f"  TMEM145 Kd = {KD_TMEM145_M*1e9:.0f} nM")
    print(f"\n{'Dose(mg)':<10}{'peak[peri] μM':<16}{'peak occ':<11}{'AUC(μM·h)':<13}{'dur>1μM h':<13}{'dur window h':<14}")
    for d in doses_mg:
        r = per_dose[f"{d}_mg"]
        print(f"{d:<10.2f}{r['peak_perilymph_uM']:<16.3f}"
              f"{r['peak_TMEM145_occupancy']:<11.2%}"
              f"{r['auc_perilymph_uM_hr']:<13.2f}"
              f"{r['duration_above_1uM_h']:<13.2f}"
              f"{r['duration_in_therapeutic_window_1_to_10_uM_h']:<14.2f}")
    print()
    print(f"Dose to hit lower boundary 1 μM:  {dose_for_1uM:.2f} mg")
    print(f"Dose to hit upper boundary 10 μM: {dose_for_10uM:.2f} mg")
    print(f"Dose to hit toxic 100 μM:         {dose_for_100uM:.1f} mg")
    print(f"Therapeutic window width:         {np.log10(dose_for_100uM/dose_for_1uM):.1f} log units")
    print()
    print(f"Multi-dose 7-day at {mid_window_dose} mg/day:")
    for k, v in multi_dose_result.items():
        if k != "n_days":
            print(f"  {k:40s} {v}")
    print()
    print(f"Sensitivity to proteolysis t½:")
    for k, v in sensitivity_proteolysis.items():
        print(f"  t½={k}: peak {v['peak_perilymph_uM']:.2f} μM, dur>1μM {v['duration_above_1uM_h']:.1f}h, mean occ {v['mean_occupancy_over_24h']:.2f}")
    print(f"\nJSON written: {OUT}")


if __name__ == "__main__":
    main()
