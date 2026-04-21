"""
STRC mRNA-LNP audiogram rescue model.

Goal: translate the PK/PD-plus-LNP-targeting numbers from
`mrna_lnp_pkpd_integration.json` into a predicted ABR audiogram for Misha
under different LNP tropism scenarios. Closes the loop between
pharmacology (solved) and clinical endpoint (per-tonotopic-band ABR dB).

Three layers:

1. Per-OHC PK/PD -- from PK/PD integration model, if the regimen keeps
   transfected-OHC trough >= 2x, that OHC is "functionally rescued".
   We take the Q6W m1psi 200 mol/OHC regimen as the reference (solved).

2. LNP tonotopic distribution -- round-window injection produces a
   base-to-apex gradient. Parametric model:

       eff_basal  = eff_mean * 2.0
       eff_mid    = eff_mean * 1.0
       eff_apical = eff_mean * 0.5

   Preserves cochlea-mean at eff_mean but introduces realistic basal
   bias consistent with published Anc80L65/ExoAAV/LNP cochlear delivery
   (Geng 2017, Landegger 2017, Leclere 2024). Where eff_zone exceeds
   1.0 (saturation) it is clamped at 1.0.

3. Functional-fraction -> ABR transfer function (per zone) -- reuses
   the calibrated transfer function in `abr_transfer_model.py`
   (A - B * log10(transduction + C)). Applied per frequency zone.

Scenarios swept:

    untargeted       eff_mean = 0.008  (standard LNP, 96/12000 OHC)
    cochlear_tropic  eff_mean = 0.05   (hypothetical next-gen LNP)
    OHC_targeted     eff_mean = 0.20   (hypothetical ligand LNP)
    hypothetical_50  eff_mean = 0.50   (required for cochlea-mean >=2x)
    Anc80L65_ref     eff_mean = 0.67   (reference AAV for comparison)

For Misha we use DFNB16 baseline 85 dB untreated threshold (from the
abr_transfer_model calibration), then predict per-band post-therapy
threshold and improvement in dB.

Runtime: seconds.
"""

from __future__ import annotations

import json
import pathlib
import sys

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

PKPD = json.loads((ROOT / "mrna_lnp_pkpd_integration.json").read_text())


def fit_abr():
    x = np.array([p[0] for p in CALIBRATION_POINTS])
    y = np.array([p[1] for p in CALIBRATION_POINTS])
    popt, _ = curve_fit(abr_log_model, x, y, p0=[20, 60, 0.01], bounds=([0, 20, 1e-4], [50, 100, 0.1]))
    return popt


ABR_PARAMS = fit_abr()


def abr_dB(functional_fraction: float) -> float:
    """Functional OHC fraction per zone -> ABR threshold dB."""
    f = float(np.clip(functional_fraction, 1e-4, 1.0))
    return float(np.clip(abr_log_model(f, *ABR_PARAMS), 15, 90))


TONOTOPIC_ZONES = {
    "basal_4_8kHz": {"gradient": 2.0, "weight": 1 / 3},
    "mid_1_4kHz": {"gradient": 1.0, "weight": 1 / 3},
    "apical_250_1kHz": {"gradient": 0.5, "weight": 1 / 3},
}


def tonotopic_eff(eff_mean: float, gradient: float) -> float:
    return float(np.clip(eff_mean * gradient, 0.0, 1.0))


def therapeutic_fraction(eff_zone: float, per_ohc_therapeutic: bool) -> float:
    """If regimen fails per-OHC, no transfected OHC is rescued -> 0.
    If regimen works per-OHC, transfected OHCs are rescued -> eff_zone."""
    return eff_zone if per_ohc_therapeutic else 0.0


def simulate_scenario(
    name: str,
    eff_mean: float,
    per_ohc_therapeutic: bool = True,
    baseline_dB: float = 85.0,
) -> dict:
    per_zone: dict = {}
    weighted_threshold = 0.0
    weighted_improvement = 0.0
    for zone, spec in TONOTOPIC_ZONES.items():
        eff_z = tonotopic_eff(eff_mean, spec["gradient"])
        rescued = therapeutic_fraction(eff_z, per_ohc_therapeutic)
        thr = abr_dB(rescued)
        improvement = baseline_dB - thr
        category, interp = classify_outcome(thr)
        per_zone[zone] = {
            "eff_ohc": eff_z,
            "rescued_fraction": rescued,
            "abr_threshold_dB": round(thr, 1),
            "improvement_dB": round(improvement, 1),
            "category": category,
            "responder": improvement >= RESPONDER_THRESHOLD_DB,
            "functional": thr <= FUNCTIONAL_THRESHOLD_DB,
            "interpretation": interp,
        }
        weighted_threshold += thr * spec["weight"]
        weighted_improvement += improvement * spec["weight"]

    any_responder = any(z["responder"] for z in per_zone.values())
    all_responder = all(z["responder"] for z in per_zone.values())
    any_functional = any(z["functional"] for z in per_zone.values())

    return {
        "name": name,
        "eff_mean": eff_mean,
        "per_ohc_therapeutic": per_ohc_therapeutic,
        "baseline_dB": baseline_dB,
        "per_zone": per_zone,
        "weighted_threshold_dB": round(weighted_threshold, 1),
        "weighted_improvement_dB": round(weighted_improvement, 1),
        "any_band_responder": any_responder,
        "all_bands_responder": all_responder,
        "any_band_functional": any_functional,
    }


MISHA = {
    "baseline_dB": {"basal_4_8kHz": 65.0, "mid_1_4kHz": 55.0, "apical_250_1kHz": 40.0},
    "note": "moderate-severe configuration, sloping; approximated from audiogram.",
}


def simulate_misha(name: str, eff_mean: float, per_ohc_therapeutic: bool = True) -> dict:
    per_zone: dict = {}
    weighted_threshold = 0.0
    weighted_improvement = 0.0
    for zone, spec in TONOTOPIC_ZONES.items():
        baseline = MISHA["baseline_dB"][zone]
        eff_z = tonotopic_eff(eff_mean, spec["gradient"])
        rescued = therapeutic_fraction(eff_z, per_ohc_therapeutic)
        thr_rescue = abr_dB(rescued)
        thr = float(min(baseline, thr_rescue))
        improvement = baseline - thr
        category, interp = classify_outcome(thr)
        per_zone[zone] = {
            "baseline_dB": baseline,
            "eff_ohc": eff_z,
            "rescued_fraction": rescued,
            "abr_post_dB": round(thr, 1),
            "improvement_dB": round(improvement, 1),
            "category": category,
            "responder": improvement >= RESPONDER_THRESHOLD_DB,
            "interpretation": interp,
        }
        weighted_threshold += thr * spec["weight"]
        weighted_improvement += improvement * spec["weight"]

    return {
        "name": name,
        "eff_mean": eff_mean,
        "per_ohc_therapeutic": per_ohc_therapeutic,
        "per_zone": per_zone,
        "weighted_post_dB": round(weighted_threshold, 1),
        "weighted_improvement_dB": round(weighted_improvement, 1),
        "any_band_responder": any(z["responder"] for z in per_zone.values()),
        "all_bands_responder": all(z["responder"] for z in per_zone.values()),
    }


SCENARIOS = [
    ("untargeted_0.8pct", 0.008),
    ("cochlear_tropic_5pct", 0.05),
    ("OHC_targeted_20pct", 0.20),
    ("hypothetical_50pct", 0.50),
    ("Anc80L65_ref_67pct", 0.67),
]


def main() -> None:
    regimens = PKPD["regimens"]
    ref = next(
        r
        for r in regimens
        if r["modification"] == "m1psi" and r["interval_d"] == 42.0 and r["dose_intra"] == 200.0 and r["targeting"] == "untargeted"
    )
    per_ohc_ok = ref["therapeutic_per_OHC"]
    tx_trough = ref["tx_trough_fold"]

    print(f"reference regimen: m1psi Q6W 200 mol/OHC, tx_trough={tx_trough:.2f}, ok={per_ohc_ok}\n")

    dfnb16_scen = []
    print("=== DFNB16 baseline (untreated 85 dB) ===")
    print(f"{'scenario':<25} {'e_mean':>8} {'basal':>7} {'mid':>7} {'apical':>7} {'wavg_post':>10} {'wavg_imp':>9} {'any_resp':>9}")
    for name, eff in SCENARIOS:
        s = simulate_scenario(name, eff, per_ohc_ok)
        dfnb16_scen.append(s)
        b = s["per_zone"]["basal_4_8kHz"]
        m = s["per_zone"]["mid_1_4kHz"]
        a = s["per_zone"]["apical_250_1kHz"]
        print(
            f"{name:<25} {eff:>8.3f} "
            f"{b['abr_threshold_dB']:>6.0f}dB {m['abr_threshold_dB']:>6.0f}dB {a['abr_threshold_dB']:>6.0f}dB "
            f"{s['weighted_threshold_dB']:>9.1f}dB {s['weighted_improvement_dB']:>8.1f}dB {str(s['any_band_responder']):>9}"
        )

    print()
    print("=== Misha per-zone audiogram (moderate-severe sloping) ===")
    print(f"{'scenario':<25} {'e_mean':>8} {'basal_post':>11} {'mid_post':>9} {'apical_post':>12} {'wavg_imp':>9}")
    misha_scen = []
    for name, eff in SCENARIOS:
        s = simulate_misha(name, eff, per_ohc_ok)
        misha_scen.append(s)
        b = s["per_zone"]["basal_4_8kHz"]
        m = s["per_zone"]["mid_1_4kHz"]
        a = s["per_zone"]["apical_250_1kHz"]
        print(
            f"{name:<25} {eff:>8.3f} "
            f"{b['baseline_dB']:>3.0f}->{b['abr_post_dB']:>3.0f} {m['baseline_dB']:>3.0f}->{m['abr_post_dB']:>3.0f} "
            f"{a['baseline_dB']:>3.0f}->{a['abr_post_dB']:>3.0f} {s['weighted_improvement_dB']:>8.1f}dB"
        )

    out = {
        "model": "STRC mRNA-LNP audiogram rescue (tonotopic x transfer function)",
        "abr_transfer_fit": {"A": float(ABR_PARAMS[0]), "B": float(ABR_PARAMS[1]), "C": float(ABR_PARAMS[2])},
        "tonotopic_gradient": {z: s["gradient"] for z, s in TONOTOPIC_ZONES.items()},
        "pkpd_reference_regimen": {
            "modification": ref["modification"],
            "interval_d": ref["interval_d"],
            "dose_intra": ref["dose_intra"],
            "tx_trough_fold": ref["tx_trough_fold"],
            "per_ohc_therapeutic": ref["therapeutic_per_OHC"],
        },
        "dfnb16_baseline_85dB": dfnb16_scen,
        "misha_per_zone": misha_scen,
        "misha_baseline": MISHA,
    }
    out_path = ROOT / "mrna_lnp_audiogram_rescue.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
