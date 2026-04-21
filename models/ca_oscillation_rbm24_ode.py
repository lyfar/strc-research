#!/usr/bin/env python3
"""
STRC Ca²⁺ Oscillation Acoustic Therapy — Phase 1 ODE Proof

8-variable ODE extending the sonogenetic 5-variable system with:
  - CaMKII autophosphorylation (fast Ca²⁺ responder)
  - Calcineurin (slow Ca²⁺ integrator) — kept from Sonogenetic model
  - RBM24 phospho state (CaMKII substrate hypothesis)
  - RBM24 nuclear fraction (phospho-gated shuttling)
  - STRC exon inclusion fraction (RBM24 nuclear → splicing outcome,
    calibrated against Sun et al. 2026 PNAS SD03: dominant exon dPSI = 0.542)

State:
  Ca_apex     [nM]   apical free calcium
  CaMKII_auto [0-1]  autonomously active CaMKII (post-T286)
  CaN_active  [0-1]  active calcineurin
  RBM24_phos  [0-1]  phosphorylated RBM24 fraction
  RBM24_nuc   [0-1]  nuclear RBM24 fraction
  STRC_incl   [0-1]  functional-isoform fraction (exon inclusion)
  STRC_mRNA   [#]    total STRC mRNA per OHC
  STRC_prot   [#]    stereocilin protein per OHC

Goal: sweep AM modulation frequency (0.001–5 Hz) at 70 dB SPL carrier,
      find frequency that maximizes STRC_prot steady-state.

Output:
  ca_oscillation_rbm24_results.json
  ca_oscillation_rbm24.png
"""

from __future__ import annotations
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# Parameters (literature-grounded)
# ============================================================

@dataclass
class Params:
    # --- OHC geometry & Ca²⁺ influx (from Sonogenetic base model) ---
    n_channels: int = 134           # Fettiplace 2017 MET channels per OHC
    g_MET: float = 150e-12          # S (single channel, Beurg 2006)
    V_drive: float = 0.1            # V driving force
    f_Ca: float = 0.15              # Ca²⁺ fraction of MET current (Lumpkin 2001)
    buffer_ratio: float = 1000.0    # apical Ca²⁺ buffer (calbindin + oncomodulin high in OHC)
    V_apex: float = 5e-14           # L apical compartment volume (Lumpkin 2001: 0.05 pL)
    Ca_rest: float = 30e-9          # M resting apical [Ca²⁺]

    # --- Ca²⁺ clearance (Mammano 1999: PMCA is fast in OHCs) ---
    k_extrusion: float = 30.0       # /s PMCA+NCX clearance (fast in OHC stereocilia)

    # --- CaMKII (Chao et al. 2010) ---
    Kd_CaMKII: float = 1e-6         # M, half-max Ca²⁺/CaM binding
    n_CaMKII_Hill: float = 2.0
    k_auto: float = 0.5             # /s autophosphorylation rate at saturating Ca²⁺
    k_dephos_CaMKII: float = 0.02   # /s PP1-mediated dephosphorylation (slow)

    # --- Calcineurin (Stemmer & Klee 1994; Dolmetsch 1998) ---
    # NOTE: k_on/k_off ratio tuned so CaN is unsaturated at physiologic Ca²⁺
    # (otherwise CaMKII:CaN ratio becomes Ca²⁺-level-insensitive).
    Kd_CaN: float = 500e-9          # M
    n_CaN_Hill: float = 4.0
    k_on_CaN: float = 0.3           # /s activation
    k_off_CaN: float = 0.05         # /s deactivation (t ~ 20 s)

    # --- RBM24 phospho dynamics (hypothesized CaMKII substrate) ---
    k_phos_RBM24: float = 0.1       # /s rate of CaMKII-mediated phosphorylation
    k_dephos_RBM24: float = 0.05    # /s rate of CaN-mediated dephosphorylation

    # --- RBM24 nuclear shuttling ---
    k_import_RBM24: float = 0.05    # /s
    k_export_RBM24: float = 0.05    # /s

    # --- STRC splicing (Sun et al. 2026 SD03) ---
    # Dominant exon event: dPSI = 0.542 between RBM24 present vs absent
    # i.e. inclusion = 0.40 (no RBM24 nuclear) -> 0.94 (full RBM24 nuclear)
    strc_inclusion_min: float = 0.40
    strc_inclusion_max: float = 0.94
    k_splice: float = 0.005         # /s splicing update rate (reflects mRNA half-life minutes)

    # --- STRC transcription / translation ---
    # Calibrated so that at STRC_inclusion = 0.40 (no RBM24 nuclear), prot_ss ≈ 3000
    # and at STRC_inclusion = 0.94 (full RBM24 nuclear), prot_ss ≈ target 15000.
    # Steady state (linear regime): prot_ss = k_translation × k_txn × incl / (k_mRNA_deg × k_prot_deg)
    k_txn_basal: float = 4.05e-5    # mRNA molecules/s per full-inclusion unit
    k_mRNA_deg: float = 0.0004      # /s mRNA degradation (~30 min t½)
    k_translation: float = 0.02     # proteins per mRNA per second
    k_prot_deg: float = 2.7e-7      # /s protein degradation (~30 days)

    # Target
    target_protein: int = 15000      # molecules/OHC (functional HTC density, Krey 2015)
    baseline_steady_prot: int = 3000 # predicted silent/baseline output


# ============================================================
# Acoustic driver: carrier at CF, amplitude-modulated at f_mod
# ============================================================

def P_open(SPL_dB: float) -> float:
    """MET channel open probability at instantaneous SPL (Fettiplace & Kim 2014)."""
    x0, dx = 70.0, 20.0
    return 0.05 + 0.80 / (1.0 + np.exp(-(SPL_dB - x0) / dx))


def sound_SPL(t: float, SPL_base: float, f_mod: float, m: float,
              mode: str = "pulse") -> float:
    """
    Two AM modes:
    - "sine": pressure modulated sinusoidally (small dynamic range)
    - "pulse": square-wave AM, 50% duty cycle — pressure toggles between
              SPL_base and SPL_base + m×20 dB (much larger Ca²⁺ swing)
    """
    if mode == "sine":
        pressure_env = 1.0 + m * np.sin(2 * np.pi * f_mod * t)
        pressure_env = max(0.1, pressure_env)
        return SPL_base + 20.0 * np.log10(pressure_env)
    else:  # pulse
        phase = (t * f_mod) % 1.0
        on = phase < 0.5
        # On: SPL_base + 20m dB; Off: SPL_base - 20m dB (symmetric)
        return SPL_base + (20.0 * m if on else -20.0 * m)


# ============================================================
# ODE system
# ============================================================

def ode_system(t, y, p: Params, SPL_base: float, f_mod: float, m: float,
               mode: str = "pulse"):
    Ca, CaMKII_a, CaN_a, RBM24_p, RBM24_n, STRC_incl, mRNA, prot = y

    # ---- Ca²⁺ dynamics ----
    instant_SPL = sound_SPL(t, SPL_base, f_mod, m, mode)
    Po = P_open(instant_SPL)
    I_MET = p.n_channels * p.g_MET * Po * p.V_drive
    I_Ca = p.f_Ca * abs(I_MET)
    J_Ca = I_Ca / (2.0 * 96485.0)  # mol/s
    dCa = J_Ca / (p.buffer_ratio * p.V_apex) - p.k_extrusion * (Ca - p.Ca_rest)

    # ---- CaMKII autophosphorylation ----
    CaMKII_bound = Ca ** p.n_CaMKII_Hill / (p.Kd_CaMKII ** p.n_CaMKII_Hill + Ca ** p.n_CaMKII_Hill)
    dCaMKII_a = (p.k_auto * CaMKII_bound * (1.0 - CaMKII_a)
                 - p.k_dephos_CaMKII * CaMKII_a)

    # ---- Calcineurin (Hill-activated on Ca²⁺) ----
    CaN_signal = Ca ** p.n_CaN_Hill / (p.Kd_CaN ** p.n_CaN_Hill + Ca ** p.n_CaN_Hill)
    dCaN_a = p.k_on_CaN * CaN_signal * (1.0 - CaN_a) - p.k_off_CaN * CaN_a

    # ---- RBM24 phosphorylation (CaMKII → phos, CaN → dephos) ----
    dRBM24_p = (p.k_phos_RBM24 * CaMKII_a * (1.0 - RBM24_p)
                - p.k_dephos_RBM24 * CaN_a * RBM24_p)

    # ---- RBM24 nuclear shuttling (phospho-gated) ----
    dRBM24_n = (p.k_import_RBM24 * RBM24_p * (1.0 - RBM24_n)
                - p.k_export_RBM24 * (1.0 - RBM24_p) * RBM24_n)

    # ---- STRC exon inclusion (RBM24_nuc determines inclusion fraction) ----
    incl_target = p.strc_inclusion_min + (p.strc_inclusion_max - p.strc_inclusion_min) * RBM24_n
    dSTRC_incl = p.k_splice * (incl_target - STRC_incl)

    # ---- STRC mRNA ----
    # Correctly-spliced mRNA accumulates only from the included fraction
    dmRNA = p.k_txn_basal * STRC_incl - p.k_mRNA_deg * mRNA

    # ---- STRC protein ----
    max_prot = p.target_protein * 2.0
    sat = max(0.0, 1.0 - prot / max_prot)
    dprot = p.k_translation * mRNA * sat - p.k_prot_deg * prot

    return [dCa, dCaMKII_a, dCaN_a, dRBM24_p, dRBM24_n, dSTRC_incl, dmRNA, dprot]


# ============================================================
# Simulation runner
# ============================================================

def run(p: Params, SPL_base: float, f_mod: float, m: float, T_hours: float = 72.0,
        mode: str = "pulse"):
    T = T_hours * 3600.0
    y0 = [p.Ca_rest, 0.0, 0.0, 0.0, 0.0, p.strc_inclusion_min, 0.0, 0.0]
    # Timestep capped at 2 s (downstream dynamics are slow); for f_mod > 1 Hz
    # we accept sub-period resolution since Ca²⁺ buffer filters >>0.5 Hz anyway.
    max_step = min(2.0, 0.2 / max(f_mod, 0.01))
    sol = solve_ivp(
        ode_system,
        (0, T),
        y0,
        args=(p, SPL_base, f_mod, m, mode),
        t_eval=np.linspace(0, T, 400),
        method="LSODA",
        rtol=1e-4,
        atol=1e-8,
        max_step=max_step,
    )
    return sol


def steady_from_tail(series, last_frac=0.25):
    n = int(len(series) * last_frac)
    return float(np.mean(series[-n:]))


# ============================================================
# Main analysis
# ============================================================

def analyze():
    p = Params()
    T_HOURS = 72.0  # 3 days — sufficient for mRNA equilibration; protein still accumulating but ranking preserved

    # -- Silence baseline --
    sol_silence = run(p, SPL_base=20.0, f_mod=0.1, m=0.0, T_hours=T_HOURS, mode="pulse")
    silence_prot = steady_from_tail(sol_silence.y[7])

    # -- AM-frequency sweep: mean 60 dB with pulse AM (m=0.5 → swings between 50 and 70 dB) --
    # Skip >1 Hz: Ca²⁺ buffering (τ~33 ms) filters fast oscillations before they can drive Dolmetsch decoder.
    f_mods = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0])
    sweep = []
    trace_at = {}
    for fm in f_mods:
        sol = run(p, SPL_base=60.0, f_mod=float(fm), m=0.5, T_hours=T_HOURS, mode="pulse")
        ca_mean = steady_from_tail(sol.y[0]) * 1e9
        ca_peak = float(np.max(sol.y[0][-200:]) * 1e9)
        CaMKII_mean = steady_from_tail(sol.y[1])
        CaN_mean = steady_from_tail(sol.y[2])
        ratio = CaMKII_mean / (CaN_mean + 1e-6)
        RBM24_n_mean = steady_from_tail(sol.y[4])
        STRC_incl = steady_from_tail(sol.y[5])
        mRNA = steady_from_tail(sol.y[6])
        prot = steady_from_tail(sol.y[7])
        sweep.append({
            "f_mod_Hz": float(fm),
            "Ca_mean_nM": ca_mean,
            "Ca_peak_nM": ca_peak,
            "CaMKII_a_mean": CaMKII_mean,
            "CaN_a_mean": CaN_mean,
            "CaMKII_over_CaN": float(ratio),
            "RBM24_nuclear_mean": RBM24_n_mean,
            "STRC_inclusion": STRC_incl,
            "STRC_mRNA_mean": mRNA,
            "STRC_protein_mean": prot,
            "fold_vs_silence": prot / max(silence_prot, 1.0),
        })
        # Keep full time traces for 3 representative frequencies
        if fm in (0.01, 0.1, 1.0):
            trace_at[fm] = {
                "t_hours": (sol.t / 3600).tolist(),
                "Ca_nM": (sol.y[0] * 1e9).tolist(),
                "CaMKII_a": sol.y[1].tolist(),
                "CaN_a": sol.y[2].tolist(),
                "RBM24_n": sol.y[4].tolist(),
                "STRC_prot": sol.y[7].tolist(),
            }

    # -- Find best frequency --
    best_idx = int(np.argmax([r["STRC_protein_mean"] for r in sweep]))
    best = sweep[best_idx]

    # -- Compare to un-modulated 60 dB (m=0) --
    sol_unmod = run(p, SPL_base=60.0, f_mod=0.1, m=0.0, T_hours=T_HOURS, mode="pulse")
    unmod_prot = steady_from_tail(sol_unmod.y[7])

    # -- Compare to un-modulated 70 dB (constant sound) --
    sol_70dB_const = run(p, SPL_base=70.0, f_mod=0.1, m=0.0, T_hours=T_HOURS, mode="pulse")
    prot_70dB_const = steady_from_tail(sol_70dB_const.y[7])

    # -- AM-depth sweep at best f_mod --
    depths = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
    depth_sweep = []
    for m in depths:
        sol = run(p, SPL_base=60.0, f_mod=float(best["f_mod_Hz"]), m=float(m), T_hours=T_HOURS, mode="pulse")
        prot = steady_from_tail(sol.y[7])
        depth_sweep.append({
            "m": float(m),
            "STRC_protein": prot,
            "fold_vs_unmod": prot / max(unmod_prot, 1.0),
        })

    results = {
        "silence_baseline_protein": silence_prot,
        "unmodulated_60dB_protein": unmod_prot,
        "unmodulated_70dB_protein": prot_70dB_const,
        "target_protein": p.target_protein,
        "frequency_sweep_60dB_m0.5_pulse": sweep,
        "best_frequency": best,
        "depth_sweep_at_best_f": depth_sweep,
        "representative_traces": trace_at,
        "simulation_hours": T_HOURS,
        "params": asdict(p),
    }
    return results


# ============================================================
# Plot
# ============================================================

def plot(results: dict, outpath: Path):
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(3, 3, hspace=0.5, wspace=0.4)
    sweep = results["frequency_sweep_60dB_m0.5_pulse"]

    # Panel 1: CaMKII:CaN ratio vs frequency (Dolmetsch decoder)
    ax = fig.add_subplot(gs[0, 0])
    fs = [r["f_mod_Hz"] for r in sweep]
    camk = [r["CaMKII_a_mean"] for r in sweep]
    can = [r["CaN_a_mean"] for r in sweep]
    ax.plot(fs, camk, "o-", label="CaMKII_auto", color="C2")
    ax.plot(fs, can, "s-", label="CaN_active", color="C3")
    ax.set_xscale("log")
    ax.set_xlabel("AM frequency (Hz)")
    ax.set_ylabel("Active fraction")
    ax.set_title("Dolmetsch frequency decoder")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 2: CaMKII:CaN ratio
    ax = fig.add_subplot(gs[0, 1])
    ratio = [r["CaMKII_over_CaN"] for r in sweep]
    ax.plot(fs, ratio, "o-", color="C4")
    ax.set_xscale("log")
    ax.set_xlabel("AM frequency (Hz)")
    ax.set_ylabel("CaMKII / CaN")
    ax.set_title("Frequency-decoded ratio")
    ax.grid(True, alpha=0.3)

    # Panel 3: RBM24 nuclear fraction
    ax = fig.add_subplot(gs[0, 2])
    rbm = [r["RBM24_nuclear_mean"] for r in sweep]
    ax.plot(fs, rbm, "o-", color="C5")
    ax.set_xscale("log")
    ax.set_xlabel("AM frequency (Hz)")
    ax.set_ylabel("Nuclear RBM24 fraction")
    ax.set_title("RBM24 nuclear retention")
    ax.grid(True, alpha=0.3)

    # Panel 4: STRC inclusion and protein
    ax = fig.add_subplot(gs[1, :2])
    prot = [r["STRC_protein_mean"] for r in sweep]
    ax2 = ax.twinx()
    ax.plot(fs, prot, "o-", color="C1", label="STRC protein (molecules/OHC)")
    incl = [r["STRC_inclusion"] for r in sweep]
    ax2.plot(fs, incl, "^-", color="C6", label="STRC exon inclusion")
    ax.axhline(results["target_protein"], linestyle="--", color="red", alpha=0.5, label="target")
    ax.axhline(results["silence_baseline_protein"], linestyle="--", color="gray", alpha=0.5, label="silence")
    ax.axhline(results["unmodulated_60dB_protein"], linestyle=":", color="black", alpha=0.5, label="unmod 60 dB")
    ax.set_xscale("log")
    ax.set_xlabel("AM frequency (Hz)")
    ax.set_ylabel("STRC protein (#)")
    ax2.set_ylabel("STRC inclusion fraction")
    ax.set_title("STRC output vs AM frequency (60 dB carrier, pulse AM m=0.5)")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], fontsize=8, loc="center right")
    ax.grid(True, alpha=0.3)

    # Panel 5: fold vs silence
    ax = fig.add_subplot(gs[1, 2])
    fold = [r["fold_vs_silence"] for r in sweep]
    ax.bar(range(len(fs)), fold, color="C1")
    ax.set_xticks(range(len(fs)))
    ax.set_xticklabels([f"{f:.3g}" for f in fs], rotation=45, fontsize=7)
    ax.set_xlabel("AM frequency (Hz)")
    ax.set_ylabel("Fold vs silence")
    ax.set_title("Protein fold-change vs silence")
    ax.axhline(1.0, linestyle="--", color="gray")
    ax.grid(True, alpha=0.3, axis="y")

    # Panel 6-8: traces at three frequencies
    trace = results["representative_traces"]
    for i, fm in enumerate(sorted(trace.keys())):
        ax = fig.add_subplot(gs[2, i])
        tr = trace[fm]
        t = np.array(tr["t_hours"])
        ax.plot(t, np.array(tr["Ca_nM"]), color="C0", alpha=0.5, label="Ca nM")
        ax.set_ylabel("Ca²⁺ (nM)", color="C0")
        ax.tick_params(axis="y", labelcolor="C0")
        ax_prot = ax.twinx()
        ax_prot.plot(t, tr["STRC_prot"], color="C1", label="STRC #")
        ax_prot.set_ylabel("STRC protein", color="C1")
        ax_prot.tick_params(axis="y", labelcolor="C1")
        ax.set_xlabel("Time (h)")
        ax.set_title(f"AM f = {fm} Hz")
        ax.grid(True, alpha=0.3)

    plt.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    out_dir = Path(__file__).parent
    print("Running 8-variable Ca²⁺ oscillation ODE — this takes ~1-3 minutes...")
    results = analyze()

    json_path = out_dir / "ca_oscillation_rbm24_results.json"
    png_path = out_dir / "ca_oscillation_rbm24.png"

    def default(o):
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        return str(o)

    with json_path.open("w") as f:
        json.dump(results, f, indent=2, default=default)
    plot(results, png_path)

    print("=" * 78)
    print("Ca²⁺ Oscillation Acoustic Therapy — Phase 1 ODE Proof")
    print("=" * 78)
    print(f"Simulation duration: {results['simulation_hours']:.0f} hours")
    print(f"Silence baseline STRC protein: {results['silence_baseline_protein']:.0f} /OHC")
    print(f"Unmodulated 60 dB STRC protein: {results['unmodulated_60dB_protein']:.0f} /OHC")
    print(f"Unmodulated 70 dB STRC protein: {results['unmodulated_70dB_protein']:.0f} /OHC")
    print(f"Target (Krey 2015): {results['target_protein']} /OHC")
    print()
    print("Frequency sweep @ 60 dB mean, pulse AM m=0.5:")
    print(f"  {'f (Hz)':>8} | {'Ca peak':>8} | {'CaMKII':>6} | {'CaN':>6} | "
          f"{'ratio':>6} | {'RBM24_n':>7} | {'Incl':>5} | {'Prot':>6} | {'×silence':>8}")
    for r in results["frequency_sweep_60dB_m0.5_pulse"]:
        print(f"  {r['f_mod_Hz']:>8.3g} | {r['Ca_peak_nM']:>8.1f} | "
              f"{r['CaMKII_a_mean']:>6.3f} | {r['CaN_a_mean']:>6.3f} | "
              f"{r['CaMKII_over_CaN']:>6.2f} | {r['RBM24_nuclear_mean']:>7.3f} | "
              f"{r['STRC_inclusion']:>5.3f} | {r['STRC_protein_mean']:>6.0f} | "
              f"{r['fold_vs_silence']:>8.2f}")
    print()
    b = results["best_frequency"]
    print(f"BEST AM FREQUENCY: {b['f_mod_Hz']} Hz → "
          f"{b['STRC_protein_mean']:.0f} proteins/OHC "
          f"({b['fold_vs_silence']:.2f}× silence, "
          f"{b['STRC_protein_mean'] / max(results['unmodulated_60dB_protein'], 1):.2f}× unmodulated 60 dB)")
    print()
    print("Depth sweep at best f:")
    for d in results["depth_sweep_at_best_f"]:
        print(f"  m={d['m']:.2f} → prot={d['STRC_protein']:.0f}, "
              f"{d['fold_vs_unmod']:.2f}× unmodulated")
    print()
    print(f"Written: {json_path.name}, {png_path.name}")


if __name__ == "__main__":
    main()
