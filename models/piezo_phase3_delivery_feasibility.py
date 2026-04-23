"""
Piezo Phase 3 — Delivery feasibility for OHC-apical conformal PVDF-TrFE film.

⚠ PARAMETER PROVENANCE WARNING — read before using numeric outputs.
Several constants below are NOT backed by primary literature and were
carried over from earlier model drafts. Treat this entire script as a
**geometry + PK/PD scaffold**, not a quantitative prediction:

  • A666 prestin-binding peptide — PHANTOM. No "Zou 2015" paper exists
    describing A666. Targeted delivery ligand for OHC apical membrane
    is a research gap; prestin-specific extracellular peptides have
    not been published to our knowledge (2026-04-23 audit).
    → K_OFF_A666, K_D ≈ 10 nM values are fabricated.

  • SELECT_S = 80 (fold selectivity OHC vs IHC/support) — PHANTOM.
    No primary source. Kept as model parameter to illustrate the
    *form* of the selectivity equation, not a measured ratio.

  • ETA_POLY = 0.7 (in-situ VDF/TrFE polymerisation efficiency) —
    PHANTOM. No primary paper reports cochlear in-situ polymerisation
    of VDF/TrFE. Range "0.6-0.8" in earlier comment was invented.

  • K_BIND_BASELINE derives Berg-Purcell diffusion-to-sink (real
    physics) but multiplies by an arbitrary "A666 contact probability"
    of 10% — this multiplier has no source.

What *is* defensible in this script:
  • Cochlear geometry (N_OHC, A_apical, NP monolayer budget)
  • Stokes-Einstein D_NP for 15 nm particle in water
  • Berg-Purcell rate constant form 4πDR
  • Perilymph volume 1 µL (mouse); RWM bolus t½ ~30 min
    (Salt & Plontke cochlear-fluid simulator)
  • Hydrogel IT t½ ~24 h (Salt review on sustained RWM release)

CONCLUSION: Phase 3 in its current form CANNOT support a delivery-feasibility
claim. Re-run after (a) an actual OHC-targeting ligand is identified
with measured Kd, and (b) VDF/TrFE in-situ polymerisation is demonstrated
ex vivo on cochlear-cell culture.

──────────────────────────────────────────────────────────────────────

Model layers (kept for reference; outputs are UPPER BOUNDS / SCAFFOLD):

  1. Capture-efficiency model (per dose, per scenario):
         f_captured = k_bind / (k_bind + k_clear)

     Three scenarios (k_clear literature-anchored; k_bind unvalidated):
        A. Bolus intratympanic (IT)          t½_clear = 30 min
        B. Hydrogel IT (chitosan-GP)         t½_clear = 24 h
        C. Direct intracochlear injection    t½_clear = 24 h, k_bind ×10

  2. Selectivity gate — FORM ONLY, coefficient PHANTOM.
  3. Polymerisation efficiency — PHANTOM.

  4. Voltage gate (from Phase 2 wall-curvature model):
         V_wall(SPL, coverage) = V_wall_100(SPL) × coverage
         Pass if V_wall ≥ 10 mV (prestin threshold) at ≥50 dB SPL
         across the clinical audiogram.

  5. Dose-response + monthly repeat-dose accumulation:
         Steady-state multi-dose coverage:
            θ_ss = 1 - (1 - θ_single_dose)^N_doses
         (each dose captures an additional 1-θ fraction of remaining
          free sites)

Parameters (literature order-of-magnitude):
    OHC count (mouse)                3,000
    OHC apical area (single cell)   25 µm²
    NP diameter                     15 nm (fits stria slit, ~30 nm max)
    NP footprint                    ~π(7.5 nm)² = 177 nm²
    Monolayer sites total           75,000 µm² / 177 nm² = 4.24e8
    Perilymph volume (mouse)        1 µL
    NP diff coefficient             3e-8 cm²/s (Stokes–Einstein, 15 nm)
    Capture radius (OHC half-area)  100 µm characteristic
    k_clear (bolus IT)              ln2/30 min = 3.85e-4 s^-1
    k_clear (hydrogel IT)           ln2/24 h   = 8.0e-6 s^-1 (48× slower)
    k_bind (baseline, diffusion)    estimated 1e-5 s^-1 (see derivation)
    k_bind (direct cochlear)        1e-4 s^-1 (10× local concentration)
    A666 selectivity                80×
    η_poly                          0.5
    Monthly dose interval           30 d

Outputs: ranked scenarios by achievable coverage, minimum effective
dose, and cumulative coverage under repeated dosing.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")

# ===== CONSTANTS =====
N_OHC = 3000
A_APICAL_UM2 = 25.0
A_APICAL_TOTAL_CM2 = N_OHC * A_APICAL_UM2 * 1e-8
NP_DIAM_NM = 15.0
NP_FOOTPRINT_NM2 = np.pi * (NP_DIAM_NM / 2) ** 2
N_MONOLAYER_TOTAL = int(A_APICAL_TOTAL_CM2 * 1e14 / NP_FOOTPRINT_NM2)   # ~4.2e8
V_PERILYMPH_ML = 1e-3

# Binding rate — derived from diffusion-to-sink (Berg-Purcell)
# 4πDR for N absorbing targets of radius R in volume V:
#   k_bind = 4πDR × N / V      (first-order in [NP], as s^-1)
# D = 3e-8 cm²/s, R = 7.5e-7 cm (NP), but effective absorber is the
# OHC array, with each cell being the "receptor" of radius ~5 µm.
D_NP = 3e-8
R_OHC_CM = 5e-4   # 5 µm OHC radius
V_CM3 = V_PERILYMPH_ML   # 1e-3 cm³
K_BIND_DIFF = 4 * np.pi * D_NP * R_OHC_CM * N_OHC / V_CM3
# ≈ 4·π·3e-8·5e-4·3000 / 1e-3 = 5.66e-4 s^-1
# Plus A666 contact probability (~10% at each encounter) = 5.7e-5
# BUT once A666 binds, k_off is slow (~1e-4), so effective net = ~5e-5

# ⚠ PHANTOM PARAMETERS — see module docstring for full provenance audit.
# Kept as named constants so downstream code reflects "model scaffold"
# status rather than "measured". Re-run Phase 3 only after real ligand + Kd.
K_BIND_BASELINE = 5e-5    # s^-1 — PHANTOM (Berg-Purcell × arbitrary 10% contact)
K_OFF_A666 = 1e-4         # s^-1 — PHANTOM (A666 peptide does not exist in lit)

ETA_POLY = 0.7      # PHANTOM — no primary source for in-situ VDF/TrFE polymerisation
SELECT_S = 80       # PHANTOM — no primary source; prestin-specific ligand unpublished


def scenario_params():
    return {
        "A_bolus_IT": {
            "description": "Single bolus intratympanic (round-window membrane diffusion)",
            "k_clear_sinv": np.log(2) / (30 * 60),
            "k_bind_sinv": K_BIND_BASELINE,
            "color": "#e63946",
        },
        "B_hydrogel_IT": {
            "description": "Chitosan-glycerophosphate hydrogel IT sustained release over ~24 h",
            "k_clear_sinv": np.log(2) / (24 * 3600),
            "k_bind_sinv": K_BIND_BASELINE,
            "color": "#f4a261",
        },
        "C_direct_cochlear": {
            "description": "Direct intracochlear injection via cochleostomy (Meniere's style)",
            "k_clear_sinv": np.log(2) / (24 * 3600),
            "k_bind_sinv": K_BIND_BASELINE * 10,   # 10× higher local conc
            "color": "#2a9d8f",
        },
    }


def capture_efficiency(k_bind, k_clear):
    return k_bind / (k_bind + k_clear)


def theta_single_dose(dose_particles, f_capture):
    n_bound = min(dose_particles * f_capture, N_MONOLAYER_TOTAL)
    return n_bound / N_MONOLAYER_TOTAL


def theta_cumulative(dose_particles, f_capture, n_doses):
    """After N equal doses, each capturing f_capture * (1 - θ_prev) × dose/N_total."""
    single = theta_single_dose(dose_particles, f_capture)
    # Each dose fills (1-θ) fraction of remaining sites
    theta = 0.0
    for _ in range(n_doses):
        theta = theta + single * (1 - theta)
    return theta


def selectivity_gate():
    N_IHC = 1000
    A_IHC_UM2 = 25
    N_SUPPORT = 12000
    A_SUPPORT_UM2 = 60
    N_off_sites_um2 = N_IHC * A_IHC_UM2 + N_SUPPORT * A_SUPPORT_UM2
    N_off_sites = int(N_off_sites_um2 * 1e6 / NP_FOOTPRINT_NM2)
    N_OHC_sites = N_MONOLAYER_TOTAL
    f_OHC = (SELECT_S * N_OHC_sites) / (SELECT_S * N_OHC_sites + N_off_sites)
    return {
        "N_OHC_sites": int(N_OHC_sites),
        "N_off_target_sites": int(N_off_sites),
        "A666_selectivity_S": SELECT_S,
        "fraction_dose_to_OHC": float(f_OHC),
        "fraction_dose_off_target": float(1 - f_OHC),
    }


def audiogram_check(coverage_eff):
    """Phase 2 wall-curvature model full audiogram (freq × SPL).
    Returns number of cells passing 10 mV prestin threshold, out of 30."""
    # From STRC Piezo Frequency Response Bundle Mechanics, 100% coverage
    v_100 = {
        250:  {40: 2.2,  50: 6.8,   60: 21.6, 70: 68.0,  80: 216.0},
        500:  {40: 3.7,  50: 11.6,  60: 36.6, 70: 116.0, 80: 366.0},
        1000: {40: 6.0,  50: 19.0,  60: 60.1, 70: 190.0, 80: 601.0},
        2000: {40: 8.9,  50: 28.3,  60: 89.4, 70: 283.0, 80: 894.0},
        4000: {40: 11.5, 50: 36.4,  60: 115.0, 70: 364.0, 80: 1151.0},
        8000: {40: 10.2, 50: 32.3,  60: 102.0, 70: 323.0, 80: 1021.0},
    }
    cells = 0
    passes = 0
    passes_50db_plus = 0
    cells_50db_plus = 0
    by_cell = {}
    for freq, row in v_100.items():
        for db, v in row.items():
            v_at = v * coverage_eff
            ok = v_at >= 10.0
            cells += 1
            if ok:
                passes += 1
            if db >= 50:
                cells_50db_plus += 1
                if ok:
                    passes_50db_plus += 1
            by_cell[f"{freq}Hz_{db}dB"] = {"V_wall_mV": round(v_at, 2), "pass": ok}
    return {
        "pass_count": passes,
        "total_cells": cells,
        "pass_fraction": passes / cells,
        "pass_count_50dB_plus": passes_50db_plus,
        "total_cells_50dB_plus": cells_50db_plus,
        "pass_fraction_50dB_plus": passes_50db_plus / cells_50db_plus,
        "by_cell": by_cell,
    }


def run_scenarios():
    scenarios = scenario_params()
    selectivity = selectivity_gate()

    doses = np.logspace(8, 13, 26)
    results = {}

    for name, p in scenarios.items():
        f_cap = capture_efficiency(p["k_bind_sinv"], p["k_clear_sinv"])
        f_cap_selected = f_cap * selectivity["fraction_dose_to_OHC"]

        rows = []
        for dose in doses:
            theta = theta_single_dose(dose, f_cap_selected)
            cov_eff_single = theta * ETA_POLY

            # repeat dosing 12 x monthly
            theta_repeat = theta_cumulative(dose, f_cap_selected, 12)
            cov_eff_repeat = theta_repeat * ETA_POLY

            ag_single = audiogram_check(cov_eff_single)
            ag_repeat = audiogram_check(cov_eff_repeat)
            rows.append({
                "dose_NPs": dose,
                "dose_mg": dose * (4/3) * np.pi * (NP_DIAM_NM * 1e-7 / 2) ** 3 * 1.8,
                "f_capture_raw": f_cap,
                "f_capture_OHC": f_cap_selected,
                "theta_single": theta,
                "coverage_single": cov_eff_single,
                "audiogram_single_pass_50dB_plus": ag_single["pass_fraction_50dB_plus"],
                "theta_12mo_repeat": theta_repeat,
                "coverage_12mo_repeat": cov_eff_repeat,
                "audiogram_12mo_pass_50dB_plus": ag_repeat["pass_fraction_50dB_plus"],
            })

        # Find minimum dose that covers at least 80% of the ≥50 dB audiogram
        min_single = next((r["dose_NPs"] for r in rows
                           if r["audiogram_single_pass_50dB_plus"] >= 0.8), None)
        min_repeat = next((r["dose_NPs"] for r in rows
                           if r["audiogram_12mo_pass_50dB_plus"] >= 0.8), None)

        results[name] = {
            "description": p["description"],
            "k_clear_sinv": p["k_clear_sinv"],
            "k_bind_sinv": p["k_bind_sinv"],
            "f_capture_raw": f_cap,
            "f_capture_selected": f_cap_selected,
            "min_dose_single_NPs": min_single,
            "min_dose_12mo_repeat_NPs": min_repeat,
            "dose_response_table": rows,
        }

    return results, selectivity


def plot(results, selectivity):
    scenarios = scenario_params()
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # A. Dose-response (single dose)
    ax = axes[0, 0]
    for name, data in results.items():
        doses = [r["dose_NPs"] for r in data["dose_response_table"]]
        cov = [r["coverage_single"] for r in data["dose_response_table"]]
        ax.semilogx(doses, cov, "-o", color=scenarios[name]["color"], markersize=3,
                    label=name.replace("_", " "))
    ax.axhline(0.6, ls="--", color="gray", label="60% voltage gate")
    ax.set_xlabel("Dose (NPs in 1 µL perilymph)")
    ax.set_ylabel("Effective coverage after single dose")
    ax.set_title("Single-dose dose-response")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # B. Dose-response (12 months monthly)
    ax = axes[0, 1]
    for name, data in results.items():
        doses = [r["dose_NPs"] for r in data["dose_response_table"]]
        cov = [r["coverage_12mo_repeat"] for r in data["dose_response_table"]]
        ax.semilogx(doses, cov, "-o", color=scenarios[name]["color"], markersize=3,
                    label=name.replace("_", " "))
    ax.axhline(0.6, ls="--", color="gray", label="60% voltage gate")
    ax.set_xlabel("Per-dose NPs (× 12 monthly)")
    ax.set_ylabel("Cumulative 12-month coverage")
    ax.set_title("Monthly-repeat dose-response (12 doses)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # C. Capture efficiency bar chart
    ax = axes[1, 0]
    names = list(results.keys())
    fcap = [results[n]["f_capture_raw"] for n in names]
    fcap_sel = [results[n]["f_capture_selected"] for n in names]
    x = np.arange(len(names))
    ax.bar(x - 0.2, fcap, width=0.4, label="raw capture", color="#264653")
    ax.bar(x + 0.2, fcap_sel, width=0.4, label="× A666 selectivity",
           color="#2a9d8f")
    ax.set_xticks(x)
    ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=8)
    ax.set_ylabel("Fraction of dose landing on OHC apical")
    ax.set_yscale("log")
    ax.set_title("Capture efficiency per scenario")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    # D. Voltage gate at 60% and actual coverage for best scenario
    ax = axes[1, 1]
    db_labels = ["40 dB", "50 dB", "60 dB", "70 dB", "80 dB"]
    best_sc = "C_direct_cochlear"
    best_dose = results[best_sc]["min_dose_12mo_repeat_NPs"]
    if best_dose:
        cov_at_min = 0.6
    else:
        cov_at_min = results[best_sc]["dose_response_table"][-1]["coverage_12mo_repeat"]

    v_at = [11.5 * cov_at_min, 36.4 * cov_at_min, 115 * cov_at_min,
            364 * cov_at_min, 1151 * cov_at_min]
    colors = ["#2a9d8f" if v >= 10 else "#e63946" for v in v_at]
    ax.bar(db_labels, v_at, color=colors, alpha=0.8)
    ax.axhline(10, ls="--", color="black", label="prestin 10 mV threshold")
    ax.set_xlabel("SPL (dB) — audiogram sweep at 4 kHz")
    ax.set_ylabel("V_wall (mV)")
    ax.set_yscale("log")
    ax.set_title(f"Voltage gate @ {best_sc} at 60% coverage")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    fig.savefig(OUT_DIR / "piezo_phase3_delivery.png", dpi=120)
    plt.close(fig)


def main():
    results, selectivity = run_scenarios()
    plot(results, selectivity)

    result = {
        "parameters": {
            "N_OHC": N_OHC,
            "N_monolayer_total": int(N_MONOLAYER_TOTAL),
            "NP_diam_nm": NP_DIAM_NM,
            "NP_footprint_nm2": float(NP_FOOTPRINT_NM2),
            "V_perilymph_mL": V_PERILYMPH_ML,
            "k_bind_baseline_sinv": K_BIND_BASELINE,
            "eta_polymerisation": ETA_POLY,
            "A666_selectivity": SELECT_S,
        },
        "selectivity": selectivity,
        "scenarios": results,
    }

    out_json = OUT_DIR / "piezo_phase3_delivery.json"
    out_json.write_text(json.dumps(result, indent=2, default=str))

    print("=== Piezo Phase 3 Delivery Feasibility ===\n")
    print(f"Monolayer target (mouse):        {N_MONOLAYER_TOTAL:.2e} NPs")
    print(f"A666 selectivity to OHC:          {selectivity['fraction_dose_to_OHC']*100:.1f}%")
    print()
    print(f"{'Scenario':<24}{'f_cap':>10}{'f_cap×S':>12}{'min-single':>14}{'min-12mo':>14}{'ag@max':>10}")
    print("-" * 90)
    for name, r in results.items():
        ms = f"{r['min_dose_single_NPs']:.1e}" if r['min_dose_single_NPs'] else "FAIL"
        mr = f"{r['min_dose_12mo_repeat_NPs']:.1e}" if r['min_dose_12mo_repeat_NPs'] else "FAIL"
        # audiogram at max dose 12mo
        max_row = r["dose_response_table"][-1]
        ag = max_row["audiogram_12mo_pass_50dB_plus"]
        print(f"{name:<24}{r['f_capture_raw']:>10.4f}{r['f_capture_selected']:>12.4f}"
              f"{ms:>14}{mr:>14}{ag*100:>8.0f}%")

    print(f"\n→ wrote {out_json}")
    print(f"→ wrote {OUT_DIR / 'piezo_phase3_delivery.png'}")


if __name__ == "__main__":
    main()
