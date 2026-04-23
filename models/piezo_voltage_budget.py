#!/usr/bin/env python3
"""
STRC Piezoelectric TM Bioelectronic Amplifier — Physics Proof (v2)

Phase 1 computational check for [[STRC Piezoelectric TM Bioelectronic Amplifier]].

Question: can a thin PVDF-TrFE film deposited near OHC stereocilia bundles
generate enough extracellular voltage to activate prestin at normal
conversational SPL?

Physics:
  1) Piezo film bends under TM displacement -> open-circuit voltage V_oc
        V_oc = g31 * E * strain * t_film
        strain = TM_disp / R_curv  (bending over curvature radius)
  2) Film couples capacitively to OHC membrane through perilymph gap
     (perilymph is conductive at audio frequencies -> quasi-static).
     Use a PER-UNIT-AREA specific-capacitance voltage divider:
        eta_transfer = C_film_area / (C_film_area + C_membrane_area)
        C_film_area      = eps_r_film * eps_0 / t_film           (F/m²)
        C_membrane_area  = C_spec_membrane                       (F/m²)
     V_wall = V_oc * eta_transfer
  3) Prestin activates when |V_wall| >= 10 mV (Santos-Sacchi 1991).

Gate: V_wall at 60 dB SPL >= 10 mV for realistic parameters.

Outputs:
    piezo_voltage_budget_results.json
    piezo_voltage_budget.png    (5-panel figure)
"""

from __future__ import annotations
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


EPS_0 = 8.854e-12  # F/m


# ============================================================
# Material library
# ============================================================

@dataclass
class Material:
    name: str
    d31: float          # C/N  (transverse piezo coefficient)
    eps_r: float        # relative permittivity
    E_mod: float = 3e9  # Young's modulus (Pa)
    biocompat: str = "good"

    @property
    def g31(self) -> float:
        # g = d / (eps_r * eps_0)
        return self.d31 / (self.eps_r * EPS_0)


# d31 / eps_r values: Arkema Piezotech PVDF-TrFE datasheet (manufacturer);
# BaTiO3 values: textbook (bulk ceramic). PLLA d14 (shear) approximated as d31 —
# MODEL APPROXIMATION (PLLA is shear-mode, not transverse). See
# literature-params/piezoelectric-materials.md for provenance.
MATERIALS = [
    Material("PVDF-TrFE (70/30, beta)", d31=-12e-12, eps_r=10.0, biocompat="proven"),
    # Terpolymer d31: manufacturer-reported d33 ≈ 40 pC/N; d31 typically lower.
    # Use conservative −20 pC/N; true d31 composition-dependent, not primary-lit backed.
    Material("P(VDF-TrFE-CFE) terpolymer", d31=-20e-12, eps_r=50.0, biocompat="less data"),
    Material("PLLA (biodegradable, d14 approx.)", d31=-8e-12, eps_r=3.0, biocompat="proven"),
    Material("BaTiO3 (reference only)", d31=-78e-12, eps_r=1200.0, biocompat="toxic Ba release"),
]


# ============================================================
# Cochlear parameters (literature-sourced)
# ============================================================

@dataclass
class Cochlea:
    # OHC specific membrane capacitance 0.9 µF/cm² — Gentet et al. 2000
    # Biophys J 79:314–320. (Ashmore 1987 J Physiol 388:323 established
    # OHC electromotility but did not quantify C_spec.)
    C_spec_membrane: float = 9e-3   # F/m²  == 0.9 µF/cm²
    # Prestin NLC voltage sensitivity — Santos-Sacchi 1991 J Neurosci 11:3096.
    # V_threshold_mV = signal-modulation amplitude assumed sufficient for
    # amplification. NLC operates −150 to +50 mV; this is a model simplification.
    V_threshold_mV: float = 10.0
    # V_saturation_mV: no direct primary citation found. Flag as model ceiling
    # until Santos-Sacchi saturation data are sourced (literature-params/piezoelectric-materials.md).
    V_saturation_mV: float = 70.0  # UNSOURCED — model ceiling only


# ============================================================
# TM displacement model (literature)
#   Ren 2002 / Olson 2001 / Gao 2014 middle-turn measurements
# ============================================================

def TM_displacement_nm(SPL_dB: float) -> float:
    """Peak TM displacement at the characteristic-frequency place."""
    # Anchor: 0.01 nm at 20 dB, linear-in-pressure to 80 dB, compressive above
    ref_SPL, ref_disp = 20.0, 0.01
    lin = ref_disp * 10 ** ((SPL_dB - ref_SPL) / 20.0)
    if SPL_dB > 80:
        lin *= 0.3 ** ((SPL_dB - 80) / 20.0)
    return lin


# ============================================================
# Piezo electrical budget
# ============================================================

@dataclass
class FilmDesign:
    material: Material
    t_film_nm: float = 100.0        # film thickness
    R_curv_um: float = 1.0          # bending-curvature radius

    @property
    def t_film_m(self) -> float:
        return self.t_film_nm * 1e-9

    @property
    def R_curv_m(self) -> float:
        return self.R_curv_um * 1e-6

    @property
    def C_film_area(self) -> float:
        """Per-unit-area capacitance (F/m²)."""
        return self.material.eps_r * EPS_0 / self.t_film_m


def V_open_circuit(design: FilmDesign, TM_disp_nm: float) -> float:
    """Open-circuit voltage at film surface (V)."""
    strain = (TM_disp_nm * 1e-9) / design.R_curv_m   # dimensionless
    stress = design.material.E_mod * strain          # Pa
    V_oc = abs(design.material.g31) * stress * design.t_film_m
    return V_oc


def transfer_fraction(design: FilmDesign, cochlea: Cochlea) -> float:
    """Per-area capacitive divider between film and OHC membrane."""
    Cf = design.C_film_area
    Cm = cochlea.C_spec_membrane
    return Cf / (Cf + Cm)


def V_lateral_wall(design: FilmDesign, cochlea: Cochlea, SPL_dB: float) -> Dict:
    disp = TM_displacement_nm(SPL_dB)
    V_oc = V_open_circuit(design, disp)
    eta = transfer_fraction(design, cochlea)
    V_wall = V_oc * eta
    return {
        "SPL_dB": SPL_dB,
        "TM_disp_nm": disp,
        "V_oc_mV": V_oc * 1000,
        "transfer_eta": eta,
        "V_wall_mV": V_wall * 1000,
        "activates_prestin": (V_wall * 1000) >= cochlea.V_threshold_mV,
    }


# ============================================================
# Analysis: baseline + sweeps
# ============================================================

def analysis() -> dict:
    cochlea = Cochlea()

    # Baseline: PVDF-TrFE 100 nm film, 1-µm curvature (macroscopic TM bending)
    baseline = FilmDesign(material=MATERIALS[0], t_film_nm=100, R_curv_um=1.0)
    aggressive = FilmDesign(material=MATERIALS[0], t_film_nm=100, R_curv_um=0.1)  # bundle-scale curvature
    terpolymer = FilmDesign(material=MATERIALS[1], t_film_nm=100, R_curv_um=1.0)  # higher-eps material
    plla = FilmDesign(material=MATERIALS[2], t_film_nm=100, R_curv_um=1.0)        # biodegradable

    scenarios = {
        "baseline (PVDF-TrFE, 100 nm, R=1 µm)": baseline,
        "bundle-scale (PVDF-TrFE, 100 nm, R=100 nm)": aggressive,
        "terpolymer (P(VDF-TrFE-CFE), 100 nm, R=1 µm)": terpolymer,
        "biodegradable PLLA (100 nm, R=1 µm)": plla,
    }

    SPL_sweep = np.arange(30, 101, 5)
    results_by_scenario = {}
    for name, design in scenarios.items():
        rows = [V_lateral_wall(design, cochlea, float(SPL)) for SPL in SPL_sweep]
        # Find minimum SPL that activates prestin
        passing = [r for r in rows if r["activates_prestin"]]
        min_SPL_pass = passing[0]["SPL_dB"] if passing else None
        results_by_scenario[name] = {
            "design": {
                "material": design.material.name,
                "d31": design.material.d31,
                "eps_r": design.material.eps_r,
                "t_film_nm": design.t_film_nm,
                "R_curv_um": design.R_curv_um,
                "C_film_area_uF_per_cm2": design.C_film_area * 1e-4 * 1e6,
                "transfer_eta": transfer_fraction(design, cochlea),
            },
            "SPL_sweep": rows,
            "min_SPL_pass": min_SPL_pass,
            "V_wall_60dB_mV": next(r["V_wall_mV"] for r in rows if r["SPL_dB"] == 60),
            "V_wall_70dB_mV": next(r["V_wall_mV"] for r in rows if r["SPL_dB"] == 70),
            "V_wall_80dB_mV": next(r["V_wall_mV"] for r in rows if r["SPL_dB"] == 80),
        }

    # Parameter sweep: film thickness × curvature, PVDF-TrFE at 60 dB
    thicknesses = [10, 30, 50, 100, 200, 500]  # nm
    curvatures = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]  # µm
    grid = []
    for t in thicknesses:
        row = []
        for R in curvatures:
            design = FilmDesign(material=MATERIALS[0], t_film_nm=t, R_curv_um=R)
            out = V_lateral_wall(design, cochlea, 60.0)
            row.append(out["V_wall_mV"])
        grid.append(row)

    # Gate verdict
    baseline_pass = results_by_scenario["baseline (PVDF-TrFE, 100 nm, R=1 µm)"]["V_wall_60dB_mV"] >= cochlea.V_threshold_mV
    any_path_passes = any(
        s["V_wall_60dB_mV"] >= cochlea.V_threshold_mV
        for s in results_by_scenario.values()
    )

    verdict = {
        "threshold_mV": cochlea.V_threshold_mV,
        "baseline_passes_60dB": baseline_pass,
        "any_engineering_path_passes_60dB": any_path_passes,
        "engineering_levers_required_for_60dB": not baseline_pass and any_path_passes,
    }

    return {
        "cochlea": asdict(cochlea),
        "scenarios": results_by_scenario,
        "sweep_thickness_x_curvature_at_60dB_PVDF": {
            "thicknesses_nm": thicknesses,
            "curvatures_um": curvatures,
            "V_wall_mV_grid": grid,
        },
        "verdict": verdict,
    }


# ============================================================
# Plot
# ============================================================

def plot(results: dict, outpath: Path):
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)
    scenarios = results["scenarios"]
    thr = results["cochlea"]["V_threshold_mV"]

    # Panel 1: V_wall vs SPL for each scenario
    ax1 = fig.add_subplot(gs[0, :2])
    colors = ["C0", "C1", "C2", "C3"]
    for (name, s), c in zip(scenarios.items(), colors):
        SPLs = [r["SPL_dB"] for r in s["SPL_sweep"]]
        V = [r["V_wall_mV"] for r in s["SPL_sweep"]]
        ax1.plot(SPLs, V, "o-", label=name, color=c, markersize=5)
    ax1.axhline(thr, linestyle="--", color="red", alpha=0.6, label=f"Prestin threshold {thr:.0f} mV")
    ax1.axvline(60, linestyle=":", color="gray", alpha=0.6, label="60 dB (conversation)")
    ax1.axvline(70, linestyle=":", color="gray", alpha=0.4)
    ax1.set_yscale("log")
    ax1.set_xlabel("SPL (dB)")
    ax1.set_ylabel("V_wall (mV)")
    ax1.set_title("Extracellular voltage at OHC lateral wall vs SPL")
    ax1.legend(loc="upper left", fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel 2: min SPL for activation
    ax2 = fig.add_subplot(gs[0, 2])
    names = list(scenarios.keys())
    min_SPL = [scenarios[n]["min_SPL_pass"] or 110 for n in names]
    colors_bar = [("green" if s <= 60 else "orange" if s <= 75 else "red") for s in min_SPL]
    y_pos = np.arange(len(names))
    ax2.barh(y_pos, min_SPL, color=colors_bar)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([n.split(" (")[0] for n in names], fontsize=8)
    ax2.axvline(60, linestyle="--", color="red", label="60 dB target")
    ax2.set_xlabel("Min SPL to activate prestin (dB)")
    ax2.set_title("Activation threshold by scenario")
    ax2.legend(fontsize=7)
    ax2.grid(True, axis="x", alpha=0.3)

    # Panel 3: thickness × curvature grid
    ax3 = fig.add_subplot(gs[1, :2])
    sweep = results["sweep_thickness_x_curvature_at_60dB_PVDF"]
    grid = np.array(sweep["V_wall_mV_grid"])
    im = ax3.imshow(
        np.log10(grid + 1e-4),
        aspect="auto",
        origin="lower",
        cmap="RdYlGn",
        vmin=-1,
        vmax=2,
    )
    ax3.set_xticks(range(len(sweep["curvatures_um"])))
    ax3.set_xticklabels(sweep["curvatures_um"])
    ax3.set_yticks(range(len(sweep["thicknesses_nm"])))
    ax3.set_yticklabels(sweep["thicknesses_nm"])
    ax3.set_xlabel("Curvature radius R (µm)")
    ax3.set_ylabel("Film thickness (nm)")
    ax3.set_title("V_wall (mV, log₁₀) at 60 dB — PVDF-TrFE sweep")
    # overlay the 10-mV contour
    X, Y = np.meshgrid(range(len(sweep["curvatures_um"])), range(len(sweep["thicknesses_nm"])))
    cs = ax3.contour(X, Y, grid, levels=[10.0], colors="black", linewidths=2)
    ax3.clabel(cs, fmt="10 mV gate")
    fig.colorbar(im, ax=ax3, label="log₁₀(V_wall / mV)")

    # Panel 4: transfer fraction vs film thickness for each material
    ax4 = fig.add_subplot(gs[1, 2])
    ts = np.linspace(5, 500, 50)
    cochlea = Cochlea()
    for mat, c in zip(MATERIALS[:3], colors):  # skip BaTiO3 (toxic)
        etas = []
        for t in ts:
            design = FilmDesign(material=mat, t_film_nm=float(t), R_curv_um=1.0)
            etas.append(transfer_fraction(design, cochlea))
        ax4.plot(ts, etas, "-", label=mat.name.split(" (")[0], color=c)
    ax4.set_xlabel("Film thickness (nm)")
    ax4.set_ylabel("Transfer fraction η")
    ax4.set_title("Coupling efficiency vs thickness")
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    plt.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    out_dir = Path(__file__).parent
    results = analysis()

    json_path = out_dir / "piezo_voltage_budget_results.json"
    png_path = out_dir / "piezo_voltage_budget.png"

    # Convert numpy types for JSON
    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    with json_path.open("w") as f:
        json.dump(results, f, indent=2, default=default)
    plot(results, png_path)

    print("=" * 78)
    print("STRC Piezoelectric TM Bioelectronic Amplifier — Phase 1 Physics Check v2")
    print("=" * 78)
    cochlea = results["cochlea"]
    v = results["verdict"]
    for name, s in results["scenarios"].items():
        gate = "✅" if s["V_wall_60dB_mV"] >= cochlea["V_threshold_mV"] else "❌"
        print(f"  {gate} {name}")
        print(f"     V_oc@60dB: {s['SPL_sweep'][6]['V_oc_mV']:7.2f} mV | "
              f"η={s['design']['transfer_eta']:.3f} | "
              f"V_wall@60dB={s['V_wall_60dB_mV']:6.2f} mV | "
              f"min pass SPL = {s['min_SPL_pass']} dB")

    print()
    print(f"Prestin activation threshold: {cochlea['V_threshold_mV']} mV")
    print(f"Baseline passes at 60 dB?            {v['baseline_passes_60dB']}")
    print(f"Any engineered path passes at 60 dB? {v['any_engineering_path_passes_60dB']}")
    print(f"Engineering levers required?         {v['engineering_levers_required_for_60dB']}")
    print()
    print(f"Written: {json_path.name}, {png_path.name}")


if __name__ == "__main__":
    main()
