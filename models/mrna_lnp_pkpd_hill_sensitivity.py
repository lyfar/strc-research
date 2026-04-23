"""
Hill-coupling sensitivity sweep for Strategy A (h06 mRNA-LNP).

MOTIVATION
----------
The base Strategy A model uses HILL_K=200, HILL_N=2, MAX_BOOST=3 for the
RBM24→STRC splicing coupling. These are ⚠ CIRCULAR FIT: no primary measurement
of OHC-RBM24→STRC splicing dose-response exists (post-audit 2026-04-23).
They were chosen so the ODE hits "therapeutic" at a reasonable dose; by
construction the model cannot fail.

This script bounds the conclusion across the biologically plausible Hill
parameter space. Output tells us whether Strategy A's claim of "therapeutic
per-OHC at moderate mRNA dose" survives Hill-parameter uncertainty, or
collapses outside the narrow point-estimate.

METHOD
------
Fixed reference regimen (middle of earlier sweep grid, chosen for clarity):
    modification = m1psi (mRNA t½ = 12 h in OHC)
    interval_d   = 28 d
    dose_intra   = 1000 mol/OHC per dose
    horizon      = 365 d

Grid:
    HILL_K    ∈ {50, 100, 200, 300, 500}   — affinity (μ-units of RBM24 prot)
    HILL_N    ∈ {1, 1.5, 2, 3, 4}           — cooperativity
    MAX_BOOST ∈ {1.5, 2, 3, 4, 5}           — STRC fold ceiling

Total: 125 per-OHC trajectories.

Per-OHC metric: tail-30% STRC fold trough, mean, peak.
Population metric: cochlea_mean = eff·fold_tx + (1-eff)·1 at LNP eff ∈
    {0.008 untargeted, 0.05 cochlear_tropic, 0.20 OHC_targeted}.
Therapeutic gate: tx_trough ≥ 2.0 (per-OHC definition inherited from h06 main).

OUTPUT
------
  1. Full 125-row regimen table (K, n, max_boost, tx_trough, tx_mean, cochlea_*).
  2. Min MAX_BOOST required for therapeutic per-OHC, grouped by (K, n).
  3. Ceiling analysis: max cochlea_mean at each targeting across whole grid.
  4. Robustness verdict: fraction of grid where therapeutic_per_OHC holds.

This is post-hoc bounds analysis, NOT a cure for the circular fit. The only
fix for the fit itself is primary RBM24-titration data in OHCs.

RUNTIME: ~1-2 min on single CPU.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
from scipy.integrate import solve_ivp

ROOT = pathlib.Path("/Users/egorlyfar/Brain/research/strc/models")
STABILITY = json.loads((ROOT / "mrna_stability_cochlear_results.json").read_text())
DOSE_RESP = json.loads((ROOT / "rbm24_mrna_dose_response_results.json").read_text())

MRNA_HL_M1PSI_H = STABILITY["stability"]["m1psi_OHC"]["hl_h"]
RBM24_HL_D = DOSE_RESP["params"]["rbm24_hl_d"]
STRC_HL_D = DOSE_RESP["params"]["strc_hl_d"]
THRESHOLD = DOSE_RESP["params"]["threshold_fold"]

LNP_UNTARGETED = 0.008
LNP_COCHLEAR_TROPIC = 0.05
LNP_OHC_TARGETED = 0.20

K_TRANSLATE = 1.0
K_STRC_BASAL = np.log(2.0) / (STRC_HL_D * 24.0 * 3600.0)

# ─────────────────── Reference regimen (fixed) ────────────────────
REF_INTERVAL_D = 28.0
REF_DOSE = 1000.0
REF_MRNA_HL_H = MRNA_HL_M1PSI_H
HORIZON_D = 365.0

# ─────────────────── Sweep grid ───────────────────────────────────
HILL_K_GRID = [50.0, 100.0, 200.0, 300.0, 500.0]
HILL_N_GRID = [1.0, 1.5, 2.0, 3.0, 4.0]
MAX_BOOST_GRID = [1.5, 2.0, 3.0, 4.0, 5.0]


def ode_pkpd_parametric(
    t_eval_s: np.ndarray,
    dose_times_s: list[float],
    dose_intra: float,
    mrna_hl_h: float,
    hill_k: float,
    hill_n: float,
    max_boost: float,
) -> np.ndarray:
    """3-compartment ODE with parametric Hill. Returns fold STRC over baseline."""
    k_m = np.log(2.0) / (mrna_hl_h * 3600.0)
    k_p = np.log(2.0) / (RBM24_HL_D * 24.0 * 3600.0)
    k_f = K_STRC_BASAL

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        m, p, f = y
        hill = (p**hill_n) / (hill_k**hill_n + p**hill_n)
        f_target = 1.0 + (max_boost - 1.0) * hill
        return np.array([
            -k_m * m,
            K_TRANSLATE * m - k_p * p,
            k_f * (f_target - f),
        ])

    y = np.array([0.0, 0.0, 1.0])
    segments_f: list[np.ndarray] = []
    segments_t: list[np.ndarray] = []

    t_breaks = sorted(set([float(t_eval_s[0])] + dose_times_s + [float(t_eval_s[-1])]))
    t_breaks = [t for t in t_breaks if t_eval_s[0] <= t <= t_eval_s[-1]]

    for i in range(len(t_breaks) - 1):
        t0, t1 = t_breaks[i], t_breaks[i + 1]
        if t0 in dose_times_s:
            y = y.copy()
            y[0] += dose_intra
        mask = (t_eval_s >= t0) & (t_eval_s <= t1)
        t_seg = t_eval_s[mask]
        if len(t_seg) == 0:
            continue
        sol = solve_ivp(rhs, (t0, t1), y, t_eval=t_seg, method="LSODA", rtol=1e-6, atol=1e-9)
        segments_t.append(sol.t)
        segments_f.append(sol.y[2])
        y = sol.y[:, -1]

    return np.concatenate(segments_f) if segments_f else np.array([])


def simulate_point(hill_k: float, hill_n: float, max_boost: float) -> dict:
    horizon_s = HORIZON_D * 24.0 * 3600.0
    t_eval = np.linspace(0.0, horizon_s, 4000)
    n_doses = int(round(365.0 / REF_INTERVAL_D))
    dose_times = [i * REF_INTERVAL_D * 24.0 * 3600.0 for i in range(n_doses + 1)]
    dose_times = [t for t in dose_times if t < horizon_s]

    f_tx = ode_pkpd_parametric(
        t_eval, dose_times, REF_DOSE, REF_MRNA_HL_H,
        hill_k, hill_n, max_boost,
    )

    tail_mask = t_eval >= 0.7 * horizon_s
    tail = f_tx[tail_mask]
    tx_trough = float(np.min(tail))
    tx_peak = float(np.max(tail))
    tx_mean = float(np.mean(tail))

    def coch(eff: float) -> tuple[float, float, float]:
        return (
            eff * tx_trough + (1.0 - eff),
            eff * tx_peak + (1.0 - eff),
            eff * tx_mean + (1.0 - eff),
        )

    ct, cp, cm = coch(LNP_UNTARGETED)
    ct2, cp2, cm2 = coch(LNP_COCHLEAR_TROPIC)
    ct3, cp3, cm3 = coch(LNP_OHC_TARGETED)

    return {
        "hill_K": hill_k,
        "hill_n": hill_n,
        "max_boost": max_boost,
        "tx_trough_fold": tx_trough,
        "tx_peak_fold": tx_peak,
        "tx_mean_fold": tx_mean,
        "therapeutic_per_OHC": tx_trough >= THRESHOLD,
        "cochlea_untargeted_trough": ct,
        "cochlea_untargeted_mean": cm,
        "cochlea_tropic_trough": ct2,
        "cochlea_tropic_mean": cm2,
        "cochlea_OHC_targeted_trough": ct3,
        "cochlea_OHC_targeted_mean": cm3,
        "therapeutic_cochlea_OHC_targeted": ct3 >= THRESHOLD,
    }


def min_max_boost_for_therapeutic(rows: list[dict]) -> dict:
    """For each (HILL_K, HILL_N), find the minimum MAX_BOOST where
    tx_trough ≥ 2.0."""
    out: dict = {}
    for k in HILL_K_GRID:
        for n in HILL_N_GRID:
            candidates = sorted(
                [r for r in rows if r["hill_K"] == k and r["hill_n"] == n
                 and r["therapeutic_per_OHC"]],
                key=lambda r: r["max_boost"],
            )
            key = f"K={int(k):<3d}_n={n}"
            if not candidates:
                out[key] = {
                    "min_max_boost_needed": None,
                    "note": f"therapeutic NEVER reached up to MAX_BOOST={MAX_BOOST_GRID[-1]}",
                }
            else:
                out[key] = {
                    "min_max_boost_needed": candidates[0]["max_boost"],
                    "tx_trough_at_min": candidates[0]["tx_trough_fold"],
                }
    return out


def ceiling_analysis(rows: list[dict]) -> dict:
    """Max cochlea_mean achievable across the entire Hill grid, per targeting.
    Compares against the pure max_boost ceiling (no hill_K saturation)."""
    out: dict = {}
    for tgt, eff in [
        ("untargeted", LNP_UNTARGETED),
        ("cochlear_tropic", LNP_COCHLEAR_TROPIC),
        ("OHC_targeted", LNP_OHC_TARGETED),
    ]:
        col = f"cochlea_{tgt.replace('cochlear_','').replace('tropic','tropic')}_mean" \
            if tgt == "cochlear_tropic" else \
            f"cochlea_{tgt}_mean"
        key_map = {
            "untargeted": "cochlea_untargeted_mean",
            "cochlear_tropic": "cochlea_tropic_mean",
            "OHC_targeted": "cochlea_OHC_targeted_mean",
        }
        col = key_map[tgt]
        max_observed = max(r[col] for r in rows)
        # theoretical ceiling = eff * max_boost + (1 - eff), at max_boost=5
        ceiling_if_saturated = eff * MAX_BOOST_GRID[-1] + (1.0 - eff)
        out[tgt] = {
            "eff_ohc": eff,
            "max_cochlea_mean_over_grid": max_observed,
            "theoretical_ceiling_at_max_boost_5": ceiling_if_saturated,
            "gap_to_2x": 2.0 - max_observed,
            "therapeutic_possible_anywhere_on_grid": max_observed >= 2.0,
        }
    return out


def robustness(rows: list[dict]) -> dict:
    n_total = len(rows)
    n_ther = sum(1 for r in rows if r["therapeutic_per_OHC"])
    n_ther_coch_tropic = sum(1 for r in rows if r["cochlea_tropic_trough"] >= 2.0)
    n_ther_coch_ohc = sum(1 for r in rows if r["cochlea_OHC_targeted_trough"] >= 2.0)
    return {
        "n_grid_points": n_total,
        "frac_therapeutic_per_OHC": n_ther / n_total,
        "frac_therapeutic_cochlea_tropic_5pct": n_ther_coch_tropic / n_total,
        "frac_therapeutic_cochlea_OHC_targeted_20pct": n_ther_coch_ohc / n_total,
        "verdict_per_OHC": (
            "ROBUST" if n_ther / n_total >= 0.5
            else "FRAGILE" if n_ther / n_total >= 0.2
            else "CIRCULAR-ONLY (only base point works)"
        ),
    }


def main() -> None:
    print(f"Hill sensitivity sweep — ref regimen: m1psi T={REF_INTERVAL_D}d D={REF_DOSE}")
    print(f"Grid: {len(HILL_K_GRID)}×{len(HILL_N_GRID)}×{len(MAX_BOOST_GRID)} = "
          f"{len(HILL_K_GRID)*len(HILL_N_GRID)*len(MAX_BOOST_GRID)} points")
    print()

    rows: list[dict] = []
    i = 0
    for k in HILL_K_GRID:
        for n in HILL_N_GRID:
            for mb in MAX_BOOST_GRID:
                i += 1
                r = simulate_point(k, n, mb)
                rows.append(r)
                if i % 10 == 0 or i == 1:
                    print(
                        f"  [{i:3d}/{len(HILL_K_GRID)*len(HILL_N_GRID)*len(MAX_BOOST_GRID)}] "
                        f"K={k:<5.0f} n={n:<4.1f} max_boost={mb:<4.1f}  "
                        f"tx_trough={r['tx_trough_fold']:.3f}  "
                        f"ther_per_OHC={r['therapeutic_per_OHC']}",
                        flush=True,
                    )

    min_mb = min_max_boost_for_therapeutic(rows)
    ceiling = ceiling_analysis(rows)
    robust = robustness(rows)

    out = {
        "model": "Strategy A Hill-coupling sensitivity sweep",
        "reference_regimen": {
            "modification": "m1psi",
            "interval_d": REF_INTERVAL_D,
            "dose_intra_per_OHC": REF_DOSE,
            "horizon_d": HORIZON_D,
            "threshold_fold": THRESHOLD,
        },
        "grid": {
            "hill_K": HILL_K_GRID,
            "hill_n": HILL_N_GRID,
            "max_boost": MAX_BOOST_GRID,
        },
        "base_params": {
            "hill_K": DOSE_RESP["params"]["hill_km"],
            "hill_n": DOSE_RESP["params"]["hill_n"],
            "max_boost": DOSE_RESP["params"]["max_boost"],
            "note": "⚠ CIRCULAR FIT — no OHC RBM24→STRC dose-response primary lit",
        },
        "points": rows,
        "min_max_boost_for_therapeutic_by_Kn": min_mb,
        "ceiling_analysis": ceiling,
        "robustness": robust,
        "interpretation_key": {
            "ROBUST": "≥50% of Hill grid yields therapeutic per-OHC — conclusion tolerates circular fit",
            "FRAGILE": "20-50% — conclusion depends on Hill assumptions",
            "CIRCULAR-ONLY": "<20% — therapeutic flag is an artifact of the fit",
        },
    }

    out_path = ROOT / "mrna_lnp_pkpd_hill_sensitivity.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {out_path}")
    print(f"\nrobustness: {json.dumps(robust, indent=2)}")
    print(f"\nceiling analysis:\n{json.dumps(ceiling, indent=2)}")


if __name__ == "__main__":
    main()
