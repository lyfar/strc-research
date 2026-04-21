#!/usr/bin/env python3
"""
STRC Sonogenetics Phase 3 — AC1-CREB bifurcation analysis.

Goal: classify the AC1→cAMP→PKA→CREB→STRC steady-state structure as a
function of SPL (the bifurcation parameter), and answer:
  1. Is there a single stable fixed point across the physiological SPL range?
  2. Where is the switching threshold (SPL at half-max protein)?
  3. What is the Hill slope of the SPL → protein dose-response?
  4. Is there bistability / hysteresis under quasi-static SPL ramps?
  5. How fast do transients settle after a step SPL change?

This extends [[STRC AC1-CREB Parameter Robustness]] (Phase 2) with
linear-stability theory and quasi-static continuation rather than point-wise
perturbation. The ODE system is imported from `ca_oscillation_ac1_creb_pivot.py`.

Method summary:
- Fine SPL sweep (20-110 dB, 1 dB resolution).
- Steady state via long integration (LSODA, 72 h).
- Jacobian via central finite difference (h = 1e-4 in relative units).
- Eigenvalues via numpy.linalg.eig; report max(Re(λ)).
- Hill fit: P(SPL) = P_min + (P_max-P_min) * SPL^n / (K^n + SPL^n).
- Hysteresis: quasi-static ramp 20→110→20 dB with 48 h dwell, compare.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

MODEL_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
MODEL_FILE = MODEL_DIR / "ca_oscillation_ac1_creb_pivot.py"
OUT = MODEL_DIR / "ca_oscillation_phase3_bifurcation.json"


def load_model():
    spec = importlib.util.spec_from_file_location("ac1_creb", MODEL_FILE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ac1_creb"] = mod
    spec.loader.exec_module(mod)
    return mod


def steady_state_long(mod, SPL_dB: float, y0: np.ndarray | None = None, t_hours: float = 72.0) -> np.ndarray:
    """Integrate to steady state using LSODA. Return final state vector."""
    if y0 is None:
        y0 = np.array([
            mod.CA_REST_NM,
            0.05,
            100.0,
            0.05,
            0.02,
            mod.K_STRC_TXN_BASAL_S / mod.STRC_MRNA_DECAY_S,
            (mod.K_STRC_TXN_BASAL_S / mod.STRC_MRNA_DECAY_S) * mod.K_TRANSLATION / mod.K_PROT_DECAY_S,
        ])
    t_end = t_hours * 3600.0
    sol = solve_ivp(
        mod.rhs, [0, t_end], y0, args=(SPL_dB, 0.0),
        method="LSODA", rtol=1e-5, atol=1e-8, max_step=60.0,
    )
    if not sol.success:
        raise RuntimeError(f"ODE solver failed at SPL={SPL_dB}")
    return sol.y[:, -1]


def jacobian_fd(mod, y_ss: np.ndarray, SPL_dB: float, h_rel: float = 1e-4) -> np.ndarray:
    """Numerical Jacobian of rhs at y_ss via central differences."""
    n = len(y_ss)
    J = np.zeros((n, n))
    for j in range(n):
        h = max(h_rel * abs(y_ss[j]), 1e-8)
        y_p = y_ss.copy(); y_p[j] += h
        y_m = y_ss.copy(); y_m[j] -= h
        f_p = mod.rhs(0.0, y_p, SPL_dB, 0.0)
        f_m = mod.rhs(0.0, y_m, SPL_dB, 0.0)
        J[:, j] = (f_p - f_m) / (2 * h)
    return J


def hill(spl, P_min, P_max, K, n):
    spl = np.asarray(spl, dtype=float)
    return P_min + (P_max - P_min) * (spl ** n) / (K ** n + spl ** n)


def fit_hill(spl_arr, P_arr):
    P_min0 = float(P_arr.min())
    P_max0 = float(P_arr.max())
    K0 = float(spl_arr[np.argmin(np.abs(P_arr - (P_min0 + P_max0) / 2))])
    try:
        popt, _ = curve_fit(
            hill, spl_arr, P_arr,
            p0=[P_min0, P_max0, K0, 4.0],
            bounds=([0, 0, 1, 0.1], [np.inf, np.inf, 200, 60]),
            maxfev=20000,
        )
        P_fit = hill(spl_arr, *popt)
        ss_res = float(np.sum((P_arr - P_fit) ** 2))
        ss_tot = float(np.sum((P_arr - P_arr.mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        return {
            "P_min": float(popt[0]),
            "P_max": float(popt[1]),
            "K_dB": float(popt[2]),
            "hill_n": float(popt[3]),
            "r2": r2,
        }
    except Exception as e:
        return {"error": str(e)}


def hysteresis_ramp(mod, spl_seq: np.ndarray, y_start: np.ndarray | None, dwell_h: float = 96.0) -> tuple[np.ndarray, np.ndarray]:
    """Quasi-static SPL ramp: integrate at each SPL for `dwell_h`, carry state forward."""
    y = y_start
    proteins = []
    for spl in spl_seq:
        y_next = steady_state_long(mod, float(spl), y0=y, t_hours=dwell_h)
        proteins.append(float(y_next[-1]))
        y = y_next
    return np.array(proteins), y


def main():
    mod = load_model()

    # --- 1. Fine SPL sweep (independent initial conditions, long settle) ---
    spl_arr = np.arange(20.0, 110.0 + 1e-9, 2.0)  # 2 dB resolution
    protein_arr = np.zeros_like(spl_arr)
    max_real_eig = np.zeros_like(spl_arr)
    eig_all = []

    print("SPL sweep with per-point steady state + Jacobian eigenvalues:")
    print(f"{'SPL (dB)':>9s} {'protein':>12s} {'max Re(λ)':>12s}  {'stable?':>8s}")
    for i, spl in enumerate(spl_arr):
        y_ss = steady_state_long(mod, float(spl), y0=None, t_hours=72.0)
        protein_arr[i] = y_ss[-1]
        J = jacobian_fd(mod, y_ss, float(spl))
        eigs = np.linalg.eigvals(J)
        max_re = float(np.max(eigs.real))
        max_real_eig[i] = max_re
        eig_all.append({"SPL": float(spl), "eigs_real": eigs.real.tolist(), "eigs_imag": eigs.imag.tolist()})
        print(f"{spl:>9.1f} {protein_arr[i]:>12.1f} {max_re:>12.2e}  {'yes' if max_re < 0 else 'no':>8s}")

    stable_everywhere = bool(np.all(max_real_eig < 0))

    # --- 2. Hill fit ---
    hill_fit = fit_hill(spl_arr, protein_arr)

    # --- 3. Hysteresis: forward vs reverse quasi-static ramp ---
    # Pre-converge initial states: rest for forward ramp, high-SPL SS for reverse ramp.
    print("\nPre-converging steady states (400 h each)...")
    y_rest = steady_state_long(mod, 20.0, y0=None, t_hours=400.0)
    y_top = steady_state_long(mod, 110.0, y0=None, t_hours=400.0)
    print(f"  rest (20 dB) → protein {y_rest[-1]:.1f}")
    print(f"  top  (110 dB) → protein {y_top[-1]:.1f}")

    ramp_spl = np.arange(20.0, 110.0 + 1e-9, 10.0)
    print("\nForward quasi-static ramp (20→110, starting from rest SS):")
    p_up, _ = hysteresis_ramp(mod, ramp_spl, y_start=y_rest, dwell_h=96.0)
    for s, p in zip(ramp_spl, p_up):
        print(f"  {s:>5.1f} dB  →  protein {p:.1f}")

    print("\nReverse quasi-static ramp (110→20, starting from high-SPL SS):")
    p_down_rev, _ = hysteresis_ramp(mod, ramp_spl[::-1], y_start=y_top, dwell_h=96.0)
    # Align back to ascending order for comparison
    p_down = p_down_rev[::-1]
    for s, p in zip(ramp_spl, p_down):
        print(f"  {s:>5.1f} dB  →  protein {p:.1f}")

    diff_abs = np.abs(p_up - p_down)
    diff_rel = diff_abs / np.maximum(np.abs(p_up), 1.0)
    max_rel_gap = float(diff_rel.max())
    hysteresis_detected = max_rel_gap > 0.05

    # --- 3b. Hysteresis-closure test: hold 30 dB for 1000 h from both start states ---
    # If protein converges to the same value from both ICs, there's no true bistability
    # and the apparent gap is a slow-decay / dwell artifact.
    print("\nHysteresis closure test at 30 dB (1000 h dwell, both initial states):")
    y_30_from_rest = steady_state_long(mod, 30.0, y0=y_rest, t_hours=1000.0)
    y_30_from_top = steady_state_long(mod, 30.0, y0=y_top, t_hours=1000.0)
    p30_rest = float(y_30_from_rest[-1])
    p30_top = float(y_30_from_top[-1])
    closure_gap_rel = abs(p30_rest - p30_top) / max(abs(p30_rest), 1.0)
    closes = closure_gap_rel < 0.01
    print(f"  from rest SS → {p30_rest:.2f}")
    print(f"  from top  SS → {p30_top:.2f}")
    print(f"  closure gap  → {closure_gap_rel*100:.3f}%  ({'CLOSED — no bistability' if closes else 'OPEN — possible bistability'})")

    # --- 3c. Smooth-Ca control: remove the SPL<=30 step, re-sweep ---
    # The published ca_forcing has an if/else at 30 dB. We monkey-patch a smooth
    # sigmoid over all SPL to check whether the apparent switch at 31 dB is a
    # model artifact rather than a true bifurcation.
    CA_REST = mod.CA_REST_NM
    CA_CEIL = mod.CA_CEILING_NM
    SPL_HALF = mod.SPL_HALF_SATURATION

    def smooth_ca_forcing(t, SPL_dB, am_freq_hz):
        x = (SPL_dB - SPL_HALF) / 15.0
        base = CA_REST + (CA_CEIL - CA_REST) * (1.0 / (1.0 + np.exp(-x)))
        if am_freq_hz <= 0:
            return base
        phase = (t * am_freq_hz) % 1.0
        gate = 1.0 if phase < 0.5 else 0.1
        return CA_REST + (base - CA_REST) * gate

    orig_ca = mod.ca_forcing
    mod.ca_forcing = smooth_ca_forcing
    print("\nSmooth-Ca control sweep (monkey-patched, no 30 dB step):")
    spl_smooth = np.arange(20.0, 110.0 + 1e-9, 2.0)
    prot_smooth = np.zeros_like(spl_smooth)
    for i, s in enumerate(spl_smooth):
        y_ss_s = steady_state_long(mod, float(s), y0=None, t_hours=72.0)
        prot_smooth[i] = y_ss_s[-1]
    mod.ca_forcing = orig_ca  # restore
    hill_fit_smooth = fit_hill(spl_smooth, prot_smooth)
    print(f"  Hill fit (smooth): {hill_fit_smooth}")

    # --- 4. Summary verdict ---
    if not stable_everywhere:
        verdict = "UNSTABLE: at least one SPL has a positive-real eigenvalue; Hopf bifurcation possible."
    elif hysteresis_detected and not closes:
        verdict = f"BISTABLE: forward/reverse ramps differ up to {max_rel_gap*100:.1f}% AND closure test fails at 30 dB ({closure_gap_rel*100:.2f}% gap)."
    elif hysteresis_detected and closes:
        verdict = (
            f"MONO-STABLE with slow relaxation: apparent ramp gap up to {max_rel_gap*100:.1f}% "
            f"closes to {closure_gap_rel*100:.2f}% under 1000 h dwell → no true bistability; path-dependence is protein-decay artifact."
        )
    else:
        verdict = f"MONO-STABLE: single stable FP per SPL; Hill n={hill_fit.get('hill_n'):.2f}, K={hill_fit.get('K_dB'):.1f} dB."

    # --- 5. Pack results ---
    out = {
        "meta": {
            "phase": "3",
            "model_file": str(MODEL_FILE),
            "variables": ["Ca_apex_nM", "AC1", "cAMP_nM", "PKAc", "CREB_P", "mRNA", "protein"],
            "bifurcation_parameter": "SPL (dB)",
            "spl_range_dB": [float(spl_arr[0]), float(spl_arr[-1])],
            "spl_step_dB": float(spl_arr[1] - spl_arr[0]),
            "method": "LSODA long-time SS + central-difference Jacobian + eigen decomposition; quasi-static ramp for hysteresis.",
        },
        "sweep": {
            "SPL_dB": spl_arr.tolist(),
            "protein_ss": protein_arr.tolist(),
            "max_real_eig": max_real_eig.tolist(),
        },
        "stability": {
            "stable_everywhere": stable_everywhere,
            "max_real_eig_overall": float(max_real_eig.max()),
        },
        "hill_fit": hill_fit,
        "hysteresis": {
            "ramp_SPL_dB": ramp_spl.tolist(),
            "protein_forward_ramp": p_up.tolist(),
            "protein_reverse_ramp": p_down.tolist(),
            "max_relative_gap": max_rel_gap,
            "hysteresis_detected": hysteresis_detected,
            "closure_test_30dB": {
                "from_rest": p30_rest,
                "from_top": p30_top,
                "closure_gap_rel": closure_gap_rel,
                "closes": closes,
            },
        },
        "smooth_ca_control": {
            "SPL_dB": spl_smooth.tolist(),
            "protein_ss": prot_smooth.tolist(),
            "hill_fit": hill_fit_smooth,
        },
        "verdict": verdict,
        "eigenvalues_per_SPL": eig_all,
    }

    OUT.write_text(json.dumps(out, indent=2))
    print(f"\nWrote: {OUT}")
    print(f"Stable everywhere? {stable_everywhere} | Hysteresis? {hysteresis_detected}")
    print(f"Hill fit: {hill_fit}")
    print(f"Verdict: {verdict}")


if __name__ == "__main__":
    main()
