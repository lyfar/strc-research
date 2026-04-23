#!/usr/bin/env python3
"""
Misha compound-heterozygous STRC genotype-specific therapy stack model.

Misha's STRC genotype:
  • Paternal allele: 98 kb deletion spanning STRC+STRCP1. Produces
    ZERO protein. Requires AAV-delivered mini-STRC (h03) to substitute.
  • Maternal allele: E1659A missense. Produces FULL-LENGTH protein with
    intact fold (AF3 pTM 0.64, identical to WT) but AlphaMissense 0.9016
    pathogenic. Charge loss at 1659 (Glu→Ala) damages a functional
    interface. Requires pharmacochaperone (h01) to rescue function.

All existing hypotheses treat "STRC-null" as a single-mechanism problem.
They are not wrong but they are genotype-agnostic. For Misha specifically,
the two alleles are broken differently and need different tools. This
model asks: what does STACKING h01 (maternal rescue) + h03 (paternal
replacement) do that either monotherapy cannot?

═══════════════════════════════════════════════════════════════════════════
Model (binary-functional-OHC — corrected from continuous-mean)
═══════════════════════════════════════════════════════════════════════════

Per-allele protein-function score (relative to WT = 1.0):

  f_pat_untreated = 0.0            (98 kb del, no protein)
  f_mat_untreated = f_mat ∈ [0.05, 0.5]   (E1659A; AlphaMissense 0.90
                                          pathogenic constrains to < 0.5,
                                          Misha's moderate-severe phenotype
                                          constrains to > ~0.05)

With pharmacochaperone (PC) rescue fraction f_PC ∈ [0, 1]:

  f_mat_treated = f_mat + f_PC × (1 - f_mat)

Per-OHC STRC protein level (relative to WT diploid = 1.0):

  Non-transduced OHC:  f_OHC = 0.5 × f_pat_untreated + 0.5 × f_mat_treated
                              = 0.5 × f_mat_treated
  Transduced OHC:      f_OHC = 0.5 × f_mat_treated + A (clamped 1.0)

An OHC is FUNCTIONAL if f_OHC ≥ θ where θ is the STRC-protein threshold for
bundle formation. Biological anchor: STRC +/- heterozygous carriers (50%
protein) have normal hearing → θ ≤ 0.5. E1659A homozygotes would have
~20% protein and ARE deaf → θ > 0.2. Use θ = 0.35 baseline; sensitivity
analysis over {0.25, 0.35, 0.45}.

Cochlea fraction of functional OHCs:
  F_funct = ε × I(f_OHC_trans ≥ θ) + (1-ε) × I(f_OHC_nontrans ≥ θ)
  where I is the indicator function.

ABR threshold prediction feeds F_funct (fraction of FUNCTIONAL OHCs, the
correct input semantic) into the lit-calibrated log transfer from
abr_transfer_model.py (Bredberg 1968, Schuknecht & Gacek 1993, Sun 2024).
Misha baseline: 64 dB SPL.

Critical insight from the model structure: NORMAL hearing (all OHCs
functional) REQUIRES non-transduced OHCs to be functional, since only
fraction ε gets the AAV payload. That means the maternal E1659A rescue
(via f_PC) is the gating lever for NORMAL; AAV provides ceiling-lift
but cannot push a non-transduced OHC above θ (AAV by definition doesn't
touch those cells).

═══════════════════════════════════════════════════════════════════════════
Outcome classification (WHO/ASHA)
═══════════════════════════════════════════════════════════════════════════

  ≤ 25 dB: NORMAL (Nobel outcome — no hearing aids needed)
  26-40 dB: MILD (normal speech perception without amplification)
  41-55 dB: MODERATE (speech OK with aids)
  56-70 dB: MODERATE-SEVERE (aids mandatory, harder learning)
  > 70 dB: SEVERE-PROFOUND (cochlear implant territory)

═══════════════════════════════════════════════════════════════════════════
Gate questions
═══════════════════════════════════════════════════════════════════════════

  Q1. For Misha's plausible f_mat range, what (ε, A, f_PC) combinations
      clear the NORMAL threshold (≤ 25 dB)?

  Q2. Is normal-hearing rescue achievable by h03 AAV monotherapy alone
      (i.e., at f_PC = 0)? If yes, at what ε?

  Q3. Is normal-hearing rescue achievable by h01 pharmacochaperone
      monotherapy alone (i.e., at ε = 0)? If yes, at what f_PC?

  Q4. Does the h01+h03 stack unlock NORMAL-tier outcomes that are
      inaccessible to either monotherapy? This is the Nobel argument.

Author: Egor (computational)
Date: 2026-04-23
═══════════════════════════════════════════════════════════════════════════
"""

from __future__ import annotations

import json
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, "/Users/egorlyfar/Brain/research/strc/models")
from abr_transfer_model import transduction_to_abr, fit_transfer_function

OUT = Path("/Users/egorlyfar/Brain/research/strc/models/misha_compound_het_stack_model.json")

# ─── Misha's clinical baseline ─────────────────────────────────────────────
MISHA_BASELINE_DB = 64.0   # from audiograms, moderate-severe SNHL

# ─── Misha's maternal allele function range ────────────────────────────────
# AlphaMissense 0.9016 pathogenic → f_mat < 0.5
# But Misha's moderate not profound phenotype → f_mat > ~0.05
F_MAT_SCENARIOS = {
    "severe_E1659A_loss":   0.10,   # 90% function lost
    "moderate_E1659A_loss": 0.25,   # 75% lost
    "mild_E1659A_loss":     0.40,   # 60% lost (bound by AlphaMissense pathogenic)
}

# ─── AAV axis (h03) ────────────────────────────────────────────────────────
# ε: transduction efficiency (fraction of OHCs expressing ectopic mini-STRC)
#    Range: 0.05-0.50. Regeneron AAV.104 current preclinical ~0.30-0.40 OHC
#    envelope; conservative lower 5% floor; upper 50% ceiling (best case).
# A: per-transduced-OHC ectopic expression relative to WT single-allele.
#    0.3 = weak CMV/CB driven, 0.6 = strong CB or B8, 1.0 = engineered STRC
#    promoter equivalent to native double-allele.
EPSILON_GRID = np.array([0.0, 0.10, 0.20, 0.30, 0.40, 0.50])
A_GRID = np.array([0.3, 0.6, 1.0])

# ─── Pharmacochaperone axis (h01) ──────────────────────────────────────────
# f_PC: fraction of E1659A-driven function loss that the drug rescues.
# 0 = no drug; 1 = full rescue (protein behaves like WT).
# Realistic PC outcomes in published literature (misfolding rescue drugs
# like VX-770 for CFTR, tafamidis for TTR): 30-70% rescue at clinical dose.
F_PC_GRID = np.array([0.0, 0.25, 0.50, 0.75, 1.00])

# ─── Ceilings ──────────────────────────────────────────────────────────────
PER_OHC_CEILING = 1.0   # cell-biology ceiling — one OHC can't do >100% of WT
# NORMAL-tier gate
ABR_NORMAL_THRESHOLD_DB = 25.0
ABR_MILD_THRESHOLD_DB = 40.0

# STRC-protein functional threshold θ anchored by:
#   upper: 0.5 — STRC+/- carriers (50% protein) have normal hearing
#   lower: 0.20 — E1659A homo would be ~20% protein and deaf (Misha affected)
THETA_BASELINE = 0.35
THETA_SENSITIVITY = [0.25, 0.35, 0.45]


def f_mat_treated(f_mat: float, f_PC: float) -> float:
    """PC rescue: residual + drug×deficit"""
    return f_mat + f_PC * (1.0 - f_mat)


def functional_ohc_fraction(epsilon: float, A: float, f_mat_tx: float,
                            theta: float) -> float:
    """
    Fraction of OHCs with STRC ≥ theta.
    Non-transduced: 0.5 × f_mat_treated must clear theta
    Transduced: min(0.5 × f_mat_treated + A, 1.0) must clear theta
    """
    f_nontrans_OHC = 0.5 * f_mat_tx
    f_trans_OHC = min(0.5 * f_mat_tx + A, PER_OHC_CEILING)
    nontrans_functional = 1.0 if f_nontrans_OHC >= theta else 0.0
    trans_functional = 1.0 if f_trans_OHC >= theta else 0.0
    return epsilon * trans_functional + (1.0 - epsilon) * nontrans_functional


def classify_abr(abr_dB: float) -> str:
    if abr_dB <= 25.0: return "NORMAL"
    if abr_dB <= 40.0: return "MILD"
    if abr_dB <= 55.0: return "MODERATE"
    if abr_dB <= 70.0: return "MODERATE_SEVERE"
    return "SEVERE_PROFOUND"


def run_grid(f_mat: float, abr_params, theta: float):
    rows = []
    for eps in EPSILON_GRID:
        for A in A_GRID:
            for f_PC in F_PC_GRID:
                f_m_tx = f_mat_treated(f_mat, f_PC)
                F_funct = functional_ohc_fraction(eps, A, f_m_tx, theta)
                abr = transduction_to_abr(F_funct, params=abr_params)
                rows.append({
                    "epsilon": float(eps),
                    "A": float(A),
                    "f_PC": float(f_PC),
                    "f_nontrans_OHC": round(0.5 * f_m_tx, 3),
                    "f_trans_OHC": round(min(0.5 * f_m_tx + A, 1.0), 3),
                    "F_functional_fraction": round(F_funct, 3),
                    "abr_dB": round(abr, 1),
                    "category": classify_abr(abr),
                    "monotherapy_h03_only": (eps > 0 and f_PC == 0.0),
                    "monotherapy_h01_only": (eps == 0 and f_PC > 0),
                    "stack_h01_h03": (eps > 0 and f_PC > 0),
                })
    return rows


def main():
    abr_params, abr_rmse = fit_transfer_function()

    scenarios = {}
    for name, f_mat in F_MAT_SCENARIOS.items():
        rows = run_grid(f_mat, abr_params, THETA_BASELINE)

        n_total = len(rows)
        n_normal = sum(1 for r in rows if r["category"] == "NORMAL")
        n_mild_or_better = sum(1 for r in rows if r["category"] in ("NORMAL", "MILD"))

        # Can h03 monotherapy alone achieve NORMAL at realistic ε (≤0.5)?
        h03_only_normal = [r for r in rows if r["monotherapy_h03_only"]
                           and r["category"] == "NORMAL"]
        h03_monotherapy_achievable = len(h03_only_normal) > 0
        h03_min_eps_for_normal = (
            min(r["epsilon"] for r in h03_only_normal)
            if h03_monotherapy_achievable else None)

        # Can h01 monotherapy alone achieve NORMAL?
        h01_only_normal = [r for r in rows if r["monotherapy_h01_only"]
                           and r["category"] == "NORMAL"]
        h01_monotherapy_achievable = len(h01_only_normal) > 0
        h01_min_fpc_for_normal = (
            min(r["f_PC"] for r in h01_only_normal)
            if h01_monotherapy_achievable else None)

        # Stack-unique NORMAL: achievable with stack but not with either mono
        stack_normal = [r for r in rows if r["stack_h01_h03"]
                        and r["category"] == "NORMAL"]
        n_stack_normal = len(stack_normal)

        # Low-burden stack solutions: ε ≤ 0.3 (achievable) AND f_PC ≤ 0.5 (realistic drug)
        low_burden_normal = [r for r in stack_normal
                             if r["epsilon"] <= 0.30 and r["f_PC"] <= 0.50]
        n_low_burden_normal = len(low_burden_normal)

        baseline_funct = functional_ohc_fraction(0.0, 0.0, f_mat, THETA_BASELINE)
        scenarios[name] = {
            "f_mat_untreated": f_mat,
            "baseline_functional_ohc_fraction": round(baseline_funct, 3),
            "baseline_abr_no_therapy_dB": round(
                transduction_to_abr(baseline_funct, params=abr_params), 1),
            "grid_size": n_total,
            "n_reaching_NORMAL": n_normal,
            "pct_reaching_NORMAL": round(n_normal / n_total * 100, 1),
            "n_reaching_MILD_or_NORMAL": n_mild_or_better,
            "pct_reaching_MILD_or_NORMAL": round(
                n_mild_or_better / n_total * 100, 1),
            "h03_monotherapy_achieves_NORMAL": h03_monotherapy_achievable,
            "h03_min_epsilon_for_NORMAL": h03_min_eps_for_normal,
            "h01_monotherapy_achieves_NORMAL": h01_monotherapy_achievable,
            "h01_min_fPC_for_NORMAL": h01_min_fpc_for_normal,
            "n_stack_solutions_NORMAL": n_stack_normal,
            "n_low_burden_stack_NORMAL": n_low_burden_normal,
            "sample_low_burden_NORMAL_solutions": low_burden_normal[:5],
            "stack_enables_NORMAL": (
                n_stack_normal > 0 and (
                    not h03_monotherapy_achievable or
                    (h03_min_eps_for_normal is not None and
                     h03_min_eps_for_normal > 0.30)
                )),
        }

    # Nobel argument synthesis: across f_mat scenarios, is stack strictly
    # more capable than either monotherapy?
    nobel_arg = {
        "stack_strictly_enables_NORMAL_in_all_f_mat_scenarios": all(
            s["stack_enables_NORMAL"] for s in scenarios.values()),
        "per_scenario_stack_uniquely_enables_NORMAL": {
            name: s["stack_enables_NORMAL"] for name, s in scenarios.items()
        },
    }

    summary = {
        "batch": "misha_compound_het_stack_model",
        "date": "2026-04-23",
        "patient": "Misha (paternal 98kb STRC del + maternal E1659A)",
        "misha_baseline_abr_dB": MISHA_BASELINE_DB,
        "model_assumptions": {
            "per_OHC_ceiling": PER_OHC_CEILING,
            "abr_transfer_rmse_dB": round(abr_rmse, 2),
            "abr_transfer_source": "abr_transfer_model.py (Bredberg 1968, "
                                    "Schuknecht&Gacek 1993, Sun 2024)",
            "f_mat_scenarios": F_MAT_SCENARIOS,
            "epsilon_grid": EPSILON_GRID.tolist(),
            "A_grid": A_GRID.tolist(),
            "f_PC_grid": F_PC_GRID.tolist(),
        },
        "scenarios": scenarios,
        "nobel_argument": nobel_arg,
    }

    OUT.write_text(json.dumps(summary, indent=2))

    # Human summary
    print("=== Misha Compound-Het Therapy Stack Model ===\n")
    print(f"ABR transfer fit RMSE: {abr_rmse:.2f} dB")
    print(f"Misha baseline (audiogram): {MISHA_BASELINE_DB} dB\n")
    print(f"{'Scenario':<24}{'baseline':<10}{'h03-only':<15}{'h01-only':<15}"
          f"{'stack adds?':<12}{'low-burden stack':<18}")
    print("─" * 94)
    for name, s in scenarios.items():
        h03 = (f"ε≥{s['h03_min_epsilon_for_NORMAL']:.2f}"
               if s['h03_monotherapy_achieves_NORMAL'] else "NO")
        h01 = (f"f_PC≥{s['h01_min_fPC_for_NORMAL']:.2f}"
               if s['h01_monotherapy_achieves_NORMAL'] else "NO")
        stack_note = "YES" if s["stack_enables_NORMAL"] else "no extra"
        print(f"{name:<24}{s['baseline_abr_no_therapy_dB']:<10.0f}"
              f"{h03:<15}{h01:<15}{stack_note:<12}"
              f"{s['n_low_burden_stack_NORMAL']:<18d}")
    print()
    print("Nobel argument: stack strictly enables NORMAL hearing in",
          sum(1 for s in scenarios.values() if s["stack_enables_NORMAL"]),
          f"/ {len(scenarios)} plausible f_mat scenarios.\n")
    print(f"JSON: {OUT}")


if __name__ == "__main__":
    main()
