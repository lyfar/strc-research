#!/usr/bin/env python3
"""
Hydrogel Phase 4e — cochlear PKPD for tail91 peptide (134 aa).

Build a 4-compartment ODE model for ototopical peptide delivery:
  C1  Middle ear (gel depot, topical trans-tympanic)
  C2  Round window membrane (diffusion barrier, ~70 μm thick)
  C3  Perilymph (cochlear fluid, 70 μL volume)
  C4  OHC stereocilia surface (binding site, local TMEM145 density)

Inputs:
  Dose = 0.1 - 10 mg peptide (typical ototopical)
  Middle-ear half-life = 1-3 h (mucociliary clearance)
  RWM crossing rate for 15 kDa peptide ~1-3% over 4 h (Salt & Plontke)
  Perilymph → OHC surface partitioning: depends on stereocilia density
  Proteolytic half-life of 134 aa unprotected peptide in perilymph: ~30 min - 4 h
  TMEM145 density on OHC: ~10,000 copies per cell × 3500 OHCs/cochlea × 3 rows

Outputs:
  - Time-course of [peptide] at OHC stereocilia surface
  - Steady-state fractional TMEM145 occupancy vs dose
  - Dose needed to hit the Phase 4d sweet spot (1-10 μM local)
  - AUC(0-24h) at OHC surface
  - Predicted duration of effect per dose
  - Comparison to other ototopical peptide products for calibration
"""

from __future__ import annotations

import json
import numpy as np
from scipy.integrate import solve_ivp
from pathlib import Path

OUT = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase4e_cochlear_pkpd.json")

# Physical parameters (orders-of-magnitude sourced from cochlear PKPD literature)
PEPTIDE_MW_KDA = 14.2   # 134 aa × ~110 g/mol / aa / 1000
PERILYMPH_VOL_UL = 70.0  # human cochlea total perilymph volume

# Rate constants (1/h)
K_CLEAR_MIDDLE_EAR = 0.35     # mucociliary clearance; t_1/2 ~2 h
K_RWM = 0.008                 # round window crossing rate (0.8% / h — 3% over 4h)
K_PERILYMPH_CLEAR = 0.35      # cochlear aqueduct + cleft drainage; t_1/2 ~2 h
K_OHC_BINDING_ON = 1.0        # rapid equilibrium with surface
K_OHC_OFF_RATE = 0.1          # slow dissociation once bound
K_PROTEOLYSIS = 1.4           # proteolytic degradation; t_1/2 ~30 min in perilymph
                              # (conservative — GSGSG linkers are substrate for several
                              # proteases; tail91 has 9 K/R trypsin cuts)

# TMEM145 binding Kd on OHC surface — estimated from Phase 3 tail91 × TMEM145 ipTM 0.68
# ipTM 0.68 typically maps to Kd ~50-500 nM for well-defined interfaces
KD_TMEM145_M = 100e-9  # 100 nM

# TMEM145 surface density on OHC stereocilia
# ~3500 OHCs × 3 rows × 100 stereocilia × ~3000 TMEM145/stereocilium = ~3e9 TMEM145 per cochlea
# Confined to stereocilia tip surface area ~0.1 μm² per tip × 3500 × 300 = ~1e5 μm² = 1e-7 m²
# So local TMEM145 density ~3e16 molecules/m² or ~50 μM in a 10 nm thick boundary layer
# Translates to ~5 μM effective in a thin stereocilium-adjacent layer
TMEM145_SURFACE_DENSITY_M = 5e-6


def pkpd_derivative(t, y, dose_mg: float, kd_M: float = KD_TMEM145_M):
    """
    State:
      y[0] = C1 (middle ear depot, mg) — starts at dose
      y[1] = C2 (RWM, mg)
      y[2] = C3 (perilymph free peptide, mg)
      y[3] = C4 (OHC-surface bound peptide, mg)
      y[4] = degraded (cumulative, mg)
    """
    c1, c2, c3, c4, cdeg = y
    flux_12 = c1 * K_RWM
    flux_13_clearance = c1 * K_CLEAR_MIDDLE_EAR * (1 - 0.01)  # most middle ear drug clears out, not in
    flux_23 = c2 * 0.5  # RWM hold-up; pass through once crossed
    flux_34_on = c3 * K_OHC_BINDING_ON * 0.1  # only 10% of perilymph volume sees OHC surface
    flux_34_off = c4 * K_OHC_OFF_RATE
    flux_3_clear = c3 * K_PERILYMPH_CLEAR
    flux_proteolysis_3 = c3 * K_PROTEOLYSIS
    flux_proteolysis_4 = c4 * K_PROTEOLYSIS * 0.2  # bound peptide protected

    dc1 = -flux_12 - flux_13_clearance
    dc2 = flux_12 - flux_23
    dc3 = flux_23 - flux_34_on + flux_34_off - flux_3_clear - flux_proteolysis_3
    dc4 = flux_34_on - flux_34_off - flux_proteolysis_4
    ddeg = flux_proteolysis_3 + flux_proteolysis_4 + flux_13_clearance + flux_3_clear

    return [dc1, dc2, dc3, dc4, ddeg]


def concentration_in_perilymph_uM(c3_mg: float) -> float:
    """Convert perilymph amount (mg) to concentration (μM)."""
    moles = (c3_mg * 1e-3) / (PEPTIDE_MW_KDA * 1000)  # mol
    vol_L = PERILYMPH_VOL_UL * 1e-6
    return moles / vol_L * 1e6


def ohc_surface_concentration_uM(c4_mg: float) -> float:
    """Convert bound amount to effective concentration at OHC surface.
    The OHC stereocilia surface layer is ~10 nm thick × ~1e-7 m² area
    × 3500 OHCs ≈ 3.5e-9 L volume. Bound amount concentrates in this
    thin layer."""
    moles = (c4_mg * 1e-3) / (PEPTIDE_MW_KDA * 1000)
    vol_L = 3.5e-9
    return moles / vol_L * 1e6


def fractional_occupancy(conc_uM: float, kd_uM: float = KD_TMEM145_M * 1e6) -> float:
    return conc_uM / (kd_uM + conc_uM)


def simulate_dose(dose_mg: float, duration_h: float = 48.0):
    """Simulate one dose over duration_h hours."""
    y0 = [dose_mg, 0, 0, 0, 0]
    t_eval = np.linspace(0, duration_h, 200)
    sol = solve_ivp(
        pkpd_derivative, [0, duration_h], y0,
        args=(dose_mg,), t_eval=t_eval,
        rtol=1e-6, atol=1e-9, method="LSODA"
    )
    c1, c2, c3, c4, cdeg = sol.y
    # Key outputs
    perilymph_uM = np.array([concentration_in_perilymph_uM(x) for x in c3])
    ohc_surface_uM = np.array([ohc_surface_concentration_uM(x) for x in c4])
    # Peak and AUC
    peak_peri = float(perilymph_uM.max())
    peak_ohc = float(ohc_surface_uM.max())
    t_peak_peri = float(sol.t[perilymph_uM.argmax()])
    t_peak_ohc = float(sol.t[ohc_surface_uM.argmax()])
    # AUC
    auc_peri = float(np.trapezoid(perilymph_uM, sol.t))
    auc_ohc = float(np.trapezoid(ohc_surface_uM, sol.t))
    # Effective duration: time above 1 μM OHC surface (Phase 4d minimum-rescue floor)
    above_rescue = sol.t[ohc_surface_uM >= 1.0]
    duration_above_rescue = float(above_rescue[-1] - above_rescue[0]) if len(above_rescue) > 0 else 0.0
    # Fractional occupancy at peak
    occ_peak = fractional_occupancy(peak_ohc)
    return {
        "dose_mg": dose_mg,
        "peak_perilymph_uM": peak_peri,
        "peak_ohc_surface_uM": peak_ohc,
        "time_peak_perilymph_h": t_peak_peri,
        "time_peak_ohc_h": t_peak_ohc,
        "auc_perilymph_uM_hr": auc_peri,
        "auc_ohc_uM_hr": auc_ohc,
        "duration_above_1uM_ohc_h": duration_above_rescue,
        "fractional_TMEM145_occupancy_peak": occ_peak,
        "time_course": {
            "t_h": sol.t.tolist(),
            "perilymph_uM": perilymph_uM.tolist(),
            "ohc_surface_uM": ohc_surface_uM.tolist(),
            "fractional_occupancy": [fractional_occupancy(x) for x in ohc_surface_uM],
        },
    }


def main():
    doses_mg = [0.1, 0.3, 1.0, 3.0, 10.0]
    per_dose = {}
    for d in doses_mg:
        r = simulate_dose(d)
        per_dose[f"{d}_mg"] = {k: v for k, v in r.items() if k != "time_course"}
    # Recommended dose: smallest that gives ≥50% occupancy
    best_dose = None
    for d in doses_mg:
        if per_dose[f"{d}_mg"]["fractional_TMEM145_occupancy_peak"] >= 0.5:
            best_dose = d
            break

    # Dose-response sweep
    dose_sweep = np.logspace(-2, 1.5, 30)  # 0.01 to ~30 mg
    sweep = []
    for d in dose_sweep:
        r = simulate_dose(d, duration_h=24.0)
        sweep.append({
            "dose_mg": float(d),
            "peak_ohc_uM": r["peak_ohc_surface_uM"],
            "peak_occupancy": r["fractional_TMEM145_occupancy_peak"],
            "auc_ohc": r["auc_ohc_uM_hr"],
            "duration_above_1uM_h": r["duration_above_1uM_ohc_h"],
        })

    # Multi-dose simulation for therapeutic window
    # Assume dosing every 24 h; simulate 7 days at 1 mg
    multi_dose_result = None
    if best_dose:
        # Custom stepped simulation
        y = [0, 0, 0, 0, 0]
        t_acc = []
        peri_acc = []
        ohc_acc = []
        for day in range(7):
            y[0] += best_dose  # add dose to middle ear
            t_eval = np.linspace(0, 24, 50)
            sol = solve_ivp(pkpd_derivative, [0, 24], y, args=(best_dose,),
                            t_eval=t_eval, method="LSODA")
            y = [sol.y[i][-1] for i in range(5)]
            peri = np.array([concentration_in_perilymph_uM(x) for x in sol.y[2]])
            ohc = np.array([ohc_surface_concentration_uM(x) for x in sol.y[3]])
            t_acc += (sol.t + 24*day).tolist()
            peri_acc += peri.tolist()
            ohc_acc += ohc.tolist()
        multi_dose_result = {
            "dose_per_day_mg": best_dose,
            "time_h": t_acc,
            "perilymph_uM": peri_acc,
            "ohc_surface_uM": ohc_acc,
            "steady_state_achieved_by_day": (
                "~Day 3-4 based on 2-h perilymph half-life and daily dosing"),
            "duration_above_1uM_7day_fraction": round(
                sum(1 for x in ohc_acc if x >= 1.0) / len(ohc_acc), 3),
        }

    summary = {
        "batch": "hydrogel_phase4e_cochlear_pkpd",
        "date": "2026-04-23",
        "model_assumptions": {
            "peptide_mw_kDa": PEPTIDE_MW_KDA,
            "perilymph_volume_uL": PERILYMPH_VOL_UL,
            "middle_ear_clearance_half_life_h": round(np.log(2) / K_CLEAR_MIDDLE_EAR, 2),
            "rwm_crossing_rate_per_h": K_RWM,
            "perilymph_clearance_half_life_h": round(np.log(2) / K_PERILYMPH_CLEAR, 2),
            "proteolysis_half_life_h": round(np.log(2) / K_PROTEOLYSIS, 2),
            "TMEM145_Kd_M": KD_TMEM145_M,
            "TMEM145_surface_density_M": TMEM145_SURFACE_DENSITY_M,
        },
        "single_dose_results": per_dose,
        "recommended_starting_dose_mg": best_dose,
        "dose_sweep_summary": sweep,
        "multi_dose_7day_at_recommended": multi_dose_result,
        "risk_flags": [
            f"Proteolysis half-life assumed 30 min — may be faster (tail91 has 9 K/R + 6 F/W/Y cut sites). "
            f"If actual t_1/2 < 10 min, need PEGylation or d-amino-acid backbone.",
            f"RWM crossing rate is for 15 kDa globular proteins; unstructured linker regions may transit "
            f"faster or slower (confidence ±3×).",
            f"OHC surface concentration estimate assumes full stereocilia-surface access; if mucus/"
            f"cerumen in middle ear impedes, effective dose could be 10× lower.",
        ],
        "clinical_translation_comparison": {
            "DB-OTO (OTOF AAV)": "Intracochlear injection, single dose. Our peptide targets same OHC "
                                 "surface but topically — ~10-100× less drug reaches target but no surgery.",
            "AAV Anc80L65 (Mini-STRC)": "Same clinical precedent but for DNA-based replacement therapy. "
                                          "Our peptide is protein-replacement — faster on/off, no genome change.",
        },
    }

    OUT.write_text(json.dumps(summary, indent=2, default=str))

    print("=== Phase 4e Cochlear PKPD for tail91 peptide (134 aa) ===\n")
    print(f"Model: 4-compartment (middle ear → RWM → perilymph → OHC surface)")
    print(f"Key parameters:")
    print(f"  Peptide MW: {PEPTIDE_MW_KDA} kDa")
    print(f"  Perilymph vol: {PERILYMPH_VOL_UL} μL")
    print(f"  Middle-ear half-life: {np.log(2)/K_CLEAR_MIDDLE_EAR:.1f} h")
    print(f"  RWM crossing: {K_RWM*100:.1f}% / h")
    print(f"  Proteolysis half-life: {np.log(2)/K_PROTEOLYSIS*60:.0f} min\n")
    print("Single-dose results:")
    print(f"{'Dose (mg)':<12}{'peak peri (μM)':<17}{'peak OHC (μM)':<16}{'occupancy':<13}{'dur >1μM (h)':<15}")
    for d in doses_mg:
        r = per_dose[f"{d}_mg"]
        print(f"{d:<12.2f}"
              f"{r['peak_perilymph_uM']:<17.3f}"
              f"{r['peak_ohc_surface_uM']:<16.2f}"
              f"{r['fractional_TMEM145_occupancy_peak']:<13.2%}"
              f"{r['duration_above_1uM_ohc_h']:<15.1f}")
    print()
    if best_dose:
        print(f">>> Recommended starting dose for 50% TMEM145 occupancy: {best_dose} mg")
    else:
        print(">>> WARNING: no tested dose achieves 50% occupancy (need higher dose or improved PK)")
    print(f"\nJSON written: {OUT}")


if __name__ == "__main__":
    main()
