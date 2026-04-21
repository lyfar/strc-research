#!/usr/bin/env python3
"""
STRC Piezoelectric TM — Phase 2: Frequency Response + Bundle Mechanics

Extends Phase 1 (piezo_voltage_budget.py) in three dimensions:

  1. FREQUENCY: includes film-membrane impedance divider H(ω) — at low
     frequency the membrane leaks charge through Rm, so V_wall rolls off.
  2. BUNDLE MECHANICS: hinge model (rigid stereocilium, rotational pivot)
     vs beam model (Euler-Bernoulli bending). Strain location matters.
  3. COCHLEAR PLACE: basilar-membrane transfer varies by frequency;
     compute bundle deflection at best-frequency for 4 places (200, 1k, 4k, 8k Hz).

Goal: is the hypothesis viable across a clinical audiogram?
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS0 = 8.854e-12  # F/m

# ----------------------------- Piezo materials -----------------------------
MATERIALS = {
    "PVDF-TrFE (baseline)":     {"d31": -25e-12, "eps_r": 10,  "E_film": 3.0e9},
    "PVDF-TrFE (thin)":         {"d31": -25e-12, "eps_r": 10,  "E_film": 3.0e9},
    "PVDF-TrFE terpolymer":     {"d31": -40e-12, "eps_r": 40,  "E_film": 1.2e9},
    "PLLA":                     {"d14": 10e-12,  "eps_r": 4,   "E_film": 3.5e9},  # d14 shear, treat as eq.
    "ZnO (biocompatible)":      {"d31": -5e-12,  "eps_r": 8,   "E_film": 140e9},
}

# OHC/membrane constants
C_OHC_MEM_PER_AREA = 9e-3   # F/m² (0.9 µF/cm², Ashmore 1987)
R_OHC_MEM_PER_AREA = 1e-3   # Ω·m² (typical specific resistance ~kΩ·cm²)
V_PRESTIN_THRESHOLD = 10e-3 # V, Santos-Sacchi 1991

# ----------------------------- Bundle mechanics -----------------------------
# OHC tallest stereocilium (base of cochlea)
L_STEREO = 1.5e-6            # m (shorter at base, taller at apex)
R_STEREO = 0.1e-6            # m (100 nm radius)
E_ACTIN = 2.0e9              # Pa (actin Young modulus; overestimate for bundle)
I_STEREO = np.pi * R_STEREO**4 / 4  # second moment of area

# Hinge pivot stiffness (Karavitaki 2010, OHC base): ~1-10 pN·rad for small bundle
KAPPA_PIVOT = 5e-6          # N/rad (rotational stiffness)

# Reissner membrane and tectorial coupling
# At 60 dB SPL, TM shear displacement at best-frequency site ≈ 5-30 nm
# (Gueta 2006, Ren 2011; scales linearly with pressure)
TM_DISP_60dB = {200: 30e-9, 1000: 20e-9, 4000: 10e-9, 8000: 5e-9}   # m


def spl_to_disp(spl_db: float, freq_hz: float, cf_disp_table=TM_DISP_60dB) -> float:
    """Linear pressure scaling: disp = disp_60dB × 10^((SPL-60)/20)."""
    # interpolate CF displacement
    cfs = sorted(cf_disp_table.keys())
    disp_60 = np.interp(np.log10(freq_hz), np.log10(cfs), [cf_disp_table[c] for c in cfs])
    return disp_60 * 10 ** ((spl_db - 60) / 20)


# ----------------------------- Strain models -----------------------------
def strain_hinge(delta: float, L: float, r: float) -> float:
    """Rigid stereocilium rotating about base — maximum strain at the rootlet,
    zero along the shaft. For conformal film along shaft, effective strain is 0.
    Use rootlet strain: angular deflection × r / r_rootlet.
    Approximation: if rootlet length ~50 nm, curvature κ = (δ/L)/rootlet_L,
    and strain = r × κ.
    """
    rootlet_L = 50e-9  # m
    theta = delta / L  # radians
    kappa = theta / rootlet_L
    return r * kappa


def strain_beam(delta: float, L: float, r: float, z_frac: float = 0.0) -> float:
    """Euler-Bernoulli cantilever with tip displacement delta.
    Curvature κ(z) = (3δ/L³)(L-z); max at base (z=0).
    Strain at surface = r × κ(z).
    z_frac = 0 base, 1 tip.
    """
    z = z_frac * L
    kappa = (3 * delta / L**3) * (L - z)
    return r * kappa


def strain_wall_curvature(delta: float, R_wall: float) -> float:
    """Film on a surface with effective radius of curvature under deflection.
    Large flat structures with small deflection: κ ≈ δ/L_span; strain = t_film/R.
    This is what Phase 1 did. Include for reference — not bundle mechanics.
    """
    return delta / R_wall


# ----------------------------- Piezo voltage -----------------------------
def voc_strain(d_coeff: float, E_film: float, strain: float, t_film: float) -> float:
    """Open-circuit voltage from strain."""
    # V_oc = g × stress × t_film = (d/ε)·(E·ε_strain)·t_film = d·E·ε_strain·t_film/(eps0·eps_r)
    # But d_coeff already has sign and units; keep simple.
    stress = E_film * strain  # Pa
    return abs(d_coeff) * stress * t_film / (EPS0 * 1)  # units-check: V


def voc_from_properties(d31: float, eps_r: float, E_film: float, strain: float, t_film: float) -> float:
    """Same as above but includes eps_r explicitly."""
    stress = E_film * strain  # Pa (N/m²)
    return abs(d31) * stress * t_film / (EPS0 * eps_r)


# ----------------------------- Impedance divider -----------------------------
def h_transfer(omega: float, t_film: float, eps_r: float, area_m2: float = 1e-9) -> float:
    """Frequency-dependent transfer: V_mem / V_film.

    Film: series capacitor C_film = eps0·eps_r·area / t_film
    Membrane: parallel C_mem + R_mem (per area → per patch)

    H(ω) = (1/Z_mem) / (1/Z_film + 1/Z_mem) where Z = 1/(jωC + 1/R)
    For high ω, purely capacitive: H_inf = C_film / (C_film + C_mem)
    For low ω, membrane shorts out through R_mem: H → 0.

    Returns |H(ω)|.
    """
    C_film = EPS0 * eps_r * area_m2 / t_film
    C_mem = C_OHC_MEM_PER_AREA * area_m2
    R_mem = R_OHC_MEM_PER_AREA / area_m2  # patch resistance (Ω)

    # Y_film = jωC_film (series element)
    # Z_mem = R_mem || (1/jωC_mem) = R_mem / (1 + jωRC_mem)
    # V_mem / V_source = Z_mem / (Z_film + Z_mem)
    # = (R_mem/(1+jωRC_mem)) / (1/(jωC_film) + R_mem/(1+jωRC_mem))
    # Manipulate:
    # = (jωC_film·R_mem) / ((1+jωRC_mem) + jωC_film·R_mem)
    # = jωR_mem·C_film / (1 + jωR_mem·(C_mem + C_film))

    num = 1j * omega * R_mem * C_film
    den = 1 + 1j * omega * R_mem * (C_mem + C_film)
    return abs(num / den)


def h_lowfreq_corner(t_film: float, eps_r: float, area_m2: float = 1e-9) -> float:
    """-3dB frequency below which membrane leaks signal."""
    C_film = EPS0 * eps_r * area_m2 / t_film
    C_mem = C_OHC_MEM_PER_AREA * area_m2
    R_mem = R_OHC_MEM_PER_AREA / area_m2
    tau = R_mem * (C_mem + C_film)
    return 1 / (2 * np.pi * tau)


# ----------------------------- Composite voltage -----------------------------
def v_wall_frequency(freq_hz: float, spl_db: float, material: dict, t_film: float,
                      R_wall: float = 100e-9, strain_model: str = "beam") -> dict:
    d31 = material.get("d31") or material.get("d14")
    eps_r = material["eps_r"]
    E_film = material["E_film"]

    delta = spl_to_disp(spl_db, freq_hz)
    if strain_model == "wall_curvature":
        strain = strain_wall_curvature(delta, R_wall)
    elif strain_model == "beam":
        strain = strain_beam(delta, L_STEREO, R_STEREO, z_frac=0.1)
    elif strain_model == "hinge":
        strain = strain_hinge(delta, L_STEREO, R_STEREO)
    else:
        strain = 0

    V_oc = voc_from_properties(d31, eps_r, E_film, strain, t_film)
    H = h_transfer(2 * np.pi * freq_hz, t_film, eps_r)
    V_wall = V_oc * H

    return {
        "freq_hz": freq_hz,
        "spl_db": spl_db,
        "delta_nm": round(delta * 1e9, 3),
        "strain": round(strain, 6),
        "V_oc_mV": round(V_oc * 1e3, 4),
        "H_transfer": round(H, 4),
        "V_wall_mV": round(V_wall * 1e3, 4),
        "pass_threshold": bool(V_wall >= V_PRESTIN_THRESHOLD),
    }


# ----------------------------- Main sweep -----------------------------
def main():
    out_dir = Path("/Users/egorlyfar/Brain/research/strc/models")

    # Frequency sweep at 60 dB for each strain model
    freqs = np.logspace(np.log10(125), np.log10(16000), 30)
    spl = 60
    mat = MATERIALS["PVDF-TrFE terpolymer"]
    t_film = 100e-9

    print(f"=== Piezo Phase 2 — frequency sweep at {spl} dB SPL ===")
    print(f"Material: PVDF-TrFE terpolymer, film thickness {t_film*1e9:.0f} nm\n")

    all_results = {"frequency_sweep": {}, "audiogram": [], "spl_sweep_at_1kHz": []}

    for model in ("hinge", "beam", "wall_curvature"):
        rows = [v_wall_frequency(f, spl, mat, t_film, strain_model=model) for f in freqs]
        all_results["frequency_sweep"][model] = rows
        passes = sum(1 for r in rows if r["pass_threshold"])
        print(f"  {model:>18s}: V_wall(1kHz)={next(r['V_wall_mV'] for r in rows if abs(r['freq_hz']-1000)<200):.3f} mV, "
              f"passes {passes}/{len(freqs)} freqs")

    # Audiogram: 4 clinical frequencies, full SPL range
    print(f"\n=== Audiogram: V_wall at clinical frequencies ===")
    audiogram_freqs = [250, 500, 1000, 2000, 4000, 8000]
    spl_levels = [40, 50, 60, 70, 80, 90]
    print(f"{'freq Hz':>10s} {'model':>10s}  " + "  ".join(f"{s}dB" for s in spl_levels))
    for f in audiogram_freqs:
        for model in ("beam", "wall_curvature"):
            row = {"freq_hz": f, "model": model, "v_at_spl_mV": {}}
            line = f"{f:>10d} {model:>10s}  "
            for spl in spl_levels:
                r = v_wall_frequency(f, spl, mat, t_film, strain_model=model)
                row["v_at_spl_mV"][f"{spl}dB"] = r["V_wall_mV"]
                pass_mark = "✓" if r["V_wall_mV"] >= V_PRESTIN_THRESHOLD * 1e3 else " "
                line += f" {r['V_wall_mV']:>4.1f}{pass_mark}"
            print(line)
            all_results["audiogram"].append(row)

    # SPL sweep at 1 kHz for each model
    print(f"\n=== SPL sweep at 1 kHz ===")
    for spl in range(30, 110, 5):
        res = v_wall_frequency(1000, spl, mat, t_film, strain_model="beam")
        all_results["spl_sweep_at_1kHz"].append(res)

    # Material comparison at 1 kHz, 60 dB, beam model
    print(f"\n=== Material comparison at 1 kHz 60 dB (beam model) ===")
    all_results["material_comparison"] = []
    for name, mat_props in MATERIALS.items():
        r = v_wall_frequency(1000, 60, mat_props, t_film, strain_model="beam")
        r["material"] = name
        all_results["material_comparison"].append(r)
        print(f"  {name:>30s}: V_wall={r['V_wall_mV']:>6.2f} mV, "
              f"pass={r['pass_threshold']}")

    # Low-freq corner — audiogram viability floor
    corner = h_lowfreq_corner(t_film, mat["eps_r"])
    print(f"\n=== System corner frequency (signal rolloff) ===")
    print(f"-3dB corner: {corner:.1f} Hz  (below this, membrane leaks charge)")
    all_results["lowfreq_corner_hz"] = round(corner, 3)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: V_wall vs frequency for 3 strain models
    for model, rows in all_results["frequency_sweep"].items():
        freqs_plot = [r["freq_hz"] for r in rows]
        vs = [r["V_wall_mV"] for r in rows]
        axes[0,0].plot(freqs_plot, vs, marker="o", ms=3, label=model)
    axes[0,0].axhline(V_PRESTIN_THRESHOLD * 1e3, color="red", ls="--",
                     label="prestin threshold 10 mV")
    axes[0,0].set_xscale("log"); axes[0,0].set_yscale("log")
    axes[0,0].set_xlabel("Frequency (Hz)"); axes[0,0].set_ylabel("V_wall (mV)")
    axes[0,0].set_title(f"V_wall vs frequency @ 60 dB SPL (terpolymer 100 nm)")
    axes[0,0].legend(); axes[0,0].grid(alpha=0.3, which="both")

    # Panel B: Transfer function H(ω)
    omegas = 2 * np.pi * freqs
    H_vals = [h_transfer(w, t_film, mat["eps_r"]) for w in omegas]
    axes[0,1].semilogx(freqs, H_vals, "k-", lw=2)
    axes[0,1].axhline(0.707, color="red", ls="--", label="-3 dB (0.707)")
    axes[0,1].axvline(corner, color="blue", ls=":", label=f"corner {corner:.1f} Hz")
    axes[0,1].set_xlabel("Frequency (Hz)"); axes[0,1].set_ylabel("|H(ω)|")
    axes[0,1].set_title("Film→membrane transfer function")
    axes[0,1].legend(); axes[0,1].grid(alpha=0.3)

    # Panel C: SPL sweep at 1 kHz (beam model)
    rows = all_results["spl_sweep_at_1kHz"]
    axes[1,0].plot([r["spl_db"] for r in rows], [r["V_wall_mV"] for r in rows],
                   marker="o", color="steelblue")
    axes[1,0].axhline(V_PRESTIN_THRESHOLD * 1e3, color="red", ls="--",
                     label="prestin 10 mV")
    axes[1,0].set_xlabel("SPL (dB)"); axes[1,0].set_ylabel("V_wall (mV)")
    axes[1,0].set_title("SPL input-output @ 1 kHz, beam model")
    axes[1,0].set_yscale("log"); axes[1,0].legend(); axes[1,0].grid(alpha=0.3, which="both")

    # Panel D: material comparison
    mat_names = [r["material"] for r in all_results["material_comparison"]]
    vs = [r["V_wall_mV"] for r in all_results["material_comparison"]]
    colors = ["green" if v >= 10 else "firebrick" for v in vs]
    axes[1,1].barh(mat_names, vs, color=colors)
    axes[1,1].axvline(V_PRESTIN_THRESHOLD * 1e3, color="red", ls="--",
                     label="prestin 10 mV")
    axes[1,1].set_xlabel("V_wall (mV)")
    axes[1,1].set_title("Material comparison @ 1 kHz 60 dB")
    axes[1,1].set_xscale("log"); axes[1,1].legend(); axes[1,1].grid(alpha=0.3, axis="x")

    plt.tight_layout()
    plt.savefig(out_dir / "piezo_phase2_frequency_bundle.png", dpi=150)
    plt.close()

    # Composite verdict
    audiogram_beam = [r for r in all_results["audiogram"] if r["model"] == "beam"]
    passes_60db = sum(1 for r in audiogram_beam if r["v_at_spl_mV"]["60dB"] >= 10)
    total_freqs = len(audiogram_beam)
    verdict = (
        f"Beam model passes 60 dB threshold at {passes_60db}/{total_freqs} audiogram frequencies. "
        f"System -3dB corner at {corner:.1f} Hz — clinical speech band (200 Hz - 8 kHz) largely within pass band. "
        f"{'VIABLE' if passes_60db >= total_freqs // 2 else 'MARGINAL'} for bundle-scale conformal delivery."
    )
    all_results["verdict"] = verdict

    (out_dir / "piezo_phase2_results.json").write_text(json.dumps(all_results, indent=2))
    print(f"\nVerdict: {verdict}")


if __name__ == "__main__":
    main()
