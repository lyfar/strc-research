"""
STRC mRNA-LNP Strategy B (full-length) audiogram rescue model.

Parallel to `mrna_lnp_audiogram_rescue.py` (Strategy A: RBM24 mRNA) but uses
Strategy B (full-length STRC mRNA) per-OHC pharmacology from
`strc_mrna_strategy_b_pkpd.json`. Three-layer composition:

1. Per-OHC absolute-STRC criterion (not fold-over-baseline):
     total_STRC_per_OHC = (s_endo_frac + s_exog_tx / s_WT)
     threshold_1x: total >= 1.0 (functional restoration)
     threshold_2x: total >= 2.0 (robust bundle coupling)

2. Tonotopic LNP distribution: basal 2x, mid 1x, apical 0.5x of eff_mean,
   clamped at 1.0. Same geometry as Strategy A audiogram model.

3. ABR transfer function (same calibration as Strategy A):
   ABR_dB = A - B * log10(functional_fraction + C).

Key difference from Strategy A audiogram: scenarios distinguish by
endogenous STRC fraction (s_endo_frac):

     WT_ref       s_endo_frac = 1.0  (already at WT; therapy is a boost)
     DFNB16_null  s_endo_frac = 0.0  (biallelic null; Strategy A fails here)
     misha        s_endo_frac = 0.15 (paternal null + maternal E1659A residual)

For DFNB16 null, Strategy A cannot functionally rescue any OHC (no pre-mRNA
substrate for RBM24 splicing). Strategy B can. This is the core rescue
asymmetry: Strategy A audiogram for null alleles is 0 dB; Strategy B
produces a real audiogram.

Reference regimen: cheapest Strategy B regimen that passes per-OHC threshold
per scenario, per threshold tier. Tabulated below.

Misha baseline audiogram: 65/55/40 dB basal/mid/apical (moderate sloping).
DFNB16 null reference: 85 dB all zones. Floor effect: post = min(baseline,
rescue_projection).

Runtime: seconds.
"""

from __future__ import annotations

import json
import pathlib
import sys
from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

ROOT = pathlib.Path("/Users/egorlyfar/Brain/research/strc/models")

sys.path.insert(0, str(ROOT))
from abr_transfer_model import (  # type: ignore
    CALIBRATION_POINTS,
    abr_log_model,
    classify_outcome,
    FUNCTIONAL_THRESHOLD_DB,
    RESPONDER_THRESHOLD_DB,
)

STRAT_B = json.loads((ROOT / "strc_mrna_strategy_b_pkpd.json").read_text())
REGIMENS = STRAT_B["regimens"]


def fit_abr():
    x = np.array([p[0] for p in CALIBRATION_POINTS])
    y = np.array([p[1] for p in CALIBRATION_POINTS])
    popt, _ = curve_fit(
        abr_log_model, x, y,
        p0=[20, 60, 0.01],
        bounds=([0, 20, 1e-4], [50, 100, 0.1]),
    )
    return popt


ABR_PARAMS = fit_abr()


def abr_dB(functional_fraction: float) -> float:
    f = float(np.clip(functional_fraction, 1e-4, 1.0))
    return float(np.clip(abr_log_model(f, *ABR_PARAMS), 15, 90))


TONOTOPIC_ZONES = {
    "basal_4_8kHz":    {"gradient": 2.0, "weight": 1 / 3},
    "mid_1_4kHz":      {"gradient": 1.0, "weight": 1 / 3},
    "apical_250_1kHz": {"gradient": 0.5, "weight": 1 / 3},
}

SCENARIOS = {
    "WT_ref":      {"s_endo_frac": 1.0, "baseline_dB": 20.0, "note": "homozygous WT (no hearing loss)"},
    "DFNB16_null": {"s_endo_frac": 0.0, "baseline_dB": 85.0, "note": "biallelic STRC null; Strategy A fails (no substrate)"},
    "misha":       {"s_endo_frac": 0.15, "baseline_dB_zones": {"basal_4_8kHz": 65.0, "mid_1_4kHz": 55.0, "apical_250_1kHz": 40.0}, "note": "paternal 98 kb Δ + maternal E1659A"},
}


def tonotopic_eff(eff_mean: float, gradient: float) -> float:
    return float(np.clip(eff_mean * gradient, 0.0, 1.0))


def find_min_regimen(scenario_key: str, threshold_total: float) -> dict | None:
    """Cheapest m1ψ regimen (min annual_extra_total_per_OHC) that passes
    s_endo_frac + s_exog_tx_trough >= threshold_total per transfected OHC."""
    s_endo = SCENARIOS[scenario_key]["s_endo_frac"]
    viable = [
        r for r in REGIMENS
        if r["modification"] == "m1psi"
        and r["targeting"] == "untargeted"
        and s_endo + r["tx_s_exog_trough"] >= threshold_total
    ]
    if not viable:
        return None
    return min(viable, key=lambda r: r["annual_extra_total_per_OHC"])


def baseline_for_zone(scenario_key: str, zone: str) -> float:
    scen = SCENARIOS[scenario_key]
    if "baseline_dB_zones" in scen:
        return scen["baseline_dB_zones"][zone]
    return scen["baseline_dB"]


def simulate_scenario_lnp(
    scenario_key: str,
    eff_mean: float,
    per_ohc_passed: bool,
) -> dict:
    per_zone: dict = {}
    weighted_post = 0.0
    weighted_improvement = 0.0
    for zone, spec in TONOTOPIC_ZONES.items():
        baseline = baseline_for_zone(scenario_key, zone)
        eff_z = tonotopic_eff(eff_mean, spec["gradient"])
        rescued_fraction = eff_z if per_ohc_passed else 0.0
        thr_rescue = abr_dB(rescued_fraction)
        thr_post = float(min(baseline, thr_rescue))
        improvement = baseline - thr_post
        category, interp = classify_outcome(thr_post)
        per_zone[zone] = {
            "baseline_dB": baseline,
            "eff_ohc": eff_z,
            "rescued_fraction": rescued_fraction,
            "abr_post_dB": round(thr_post, 1),
            "improvement_dB": round(improvement, 1),
            "category": category,
            "responder": improvement >= RESPONDER_THRESHOLD_DB,
            "functional": thr_post <= FUNCTIONAL_THRESHOLD_DB,
            "interpretation": interp,
        }
        weighted_post += thr_post * spec["weight"]
        weighted_improvement += improvement * spec["weight"]

    return {
        "per_zone": per_zone,
        "weighted_post_dB": round(weighted_post, 1),
        "weighted_improvement_dB": round(weighted_improvement, 1),
        "any_band_responder": any(z["responder"] for z in per_zone.values()),
        "all_bands_responder": all(z["responder"] for z in per_zone.values()),
    }


LNP_SCENARIOS = [
    ("untargeted_0.8pct",     0.008),
    ("cochlear_tropic_5pct",  0.05),
    ("OHC_targeted_20pct",    0.20),
    ("hypothetical_50pct",    0.50),
    ("Anc80L65_ref_67pct",    0.67),
]


def format_zones_row(name, eff, result):
    z = result["per_zone"]
    b = z["basal_4_8kHz"]
    m = z["mid_1_4kHz"]
    a = z["apical_250_1kHz"]
    return (
        f"{name:<24} {eff:>7.3f} "
        f"{b['baseline_dB']:>3.0f}->{b['abr_post_dB']:>3.0f} "
        f"{m['baseline_dB']:>3.0f}->{m['abr_post_dB']:>3.0f} "
        f"{a['baseline_dB']:>3.0f}->{a['abr_post_dB']:>3.0f} "
        f"wavg_imp={result['weighted_improvement_dB']:>5.1f}dB "
        f"{'RESP' if result['any_band_responder'] else '    '}"
    )


def main() -> None:
    print(f"Strategy B audiogram rescue (uses `strc_mrna_strategy_b_pkpd.json`)\n")

    report: dict = {
        "model": "STRC mRNA-LNP Strategy B audiogram rescue",
        "abr_transfer_fit": {"A": float(ABR_PARAMS[0]), "B": float(ABR_PARAMS[1]), "C": float(ABR_PARAMS[2])},
        "tonotopic_gradient": {z: s["gradient"] for z, s in TONOTOPIC_ZONES.items()},
        "scenarios_clinical": SCENARIOS,
        "scenarios_lnp": dict(LNP_SCENARIOS),
        "per_ohc_reference_regimens": {},
        "rescue_tables": {},
    }

    for threshold_name, threshold_total in [("1x_WT", 1.0), ("2x_WT", 2.0)]:
        print(f"\n======================== per-OHC threshold: {threshold_name} ========================")
        report["per_ohc_reference_regimens"][threshold_name] = {}
        report["rescue_tables"][threshold_name] = {}

        for scenario_key in SCENARIOS:
            ref = find_min_regimen(scenario_key, threshold_total)
            per_ohc_passed = ref is not None
            if ref is None:
                print(f"\n  {scenario_key:<13} per-OHC {threshold_name}: INFEASIBLE under any Strategy B dose")
                report["per_ohc_reference_regimens"][threshold_name][scenario_key] = None
            else:
                total_trough = SCENARIOS[scenario_key]["s_endo_frac"] + ref["tx_s_exog_trough"]
                print(
                    f"\n  {scenario_key:<13} per-OHC {threshold_name}: "
                    f"min m1ψ T={ref['interval_d']:.0f}d D={ref['dose_intra']:.0f} "
                    f"(annual {ref['annual_extra_total_per_OHC']:.2e} extra/OHC) "
                    f"total_trough={total_trough:.2f}"
                )
                report["per_ohc_reference_regimens"][threshold_name][scenario_key] = {
                    "interval_d": ref["interval_d"],
                    "dose_intra": ref["dose_intra"],
                    "annual_doses": ref["annual_doses"],
                    "annual_intra_total_per_OHC": ref["annual_intra_total_per_OHC"],
                    "annual_extra_total_per_OHC": ref["annual_extra_total_per_OHC"],
                    "tx_s_exog_trough": ref["tx_s_exog_trough"],
                    "total_per_OHC_trough": total_trough,
                }

            rescue_by_lnp = {}
            for name, eff in LNP_SCENARIOS:
                r = simulate_scenario_lnp(scenario_key, eff, per_ohc_passed)
                rescue_by_lnp[name] = {"eff_mean": eff, **r}
                print(f"    {format_zones_row(name, eff, r)}")
            report["rescue_tables"][threshold_name][scenario_key] = rescue_by_lnp

    out_path = ROOT / "strategy_b_audiogram_rescue.json"
    out_path.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
