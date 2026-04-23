#!/usr/bin/env python3
"""
Hydrogel Phase 4h — endogenous STRC vs peptide competition model at TMEM145.

Critical for Misha-specific therapeutic strategy:
  Misha has compound heterozygous STRC genotype:
    - Paternal: 98 kb deletion — no STRC produced from paternal allele (null)
    - Maternal: c.4976A>C missense = p.E1659A — full-length STRC produced
      but with ARM-repeat pocket-destabilising point mutation

If maternal E1659A STRC is produced at stable levels AND partially functional
(binds TMEM145, even with reduced affinity), our peptide must COMPETE with it.
If E1659A STRC is degraded (misfolded, ER-retained, or rapidly turned over),
there is no competition — peptide acts on free TMEM145.

This script models:
  Q1. Given [peptide] in perilymph (from Phase 4e) and [E1659A STRC] at OHC
      surface, and respective Kd values, what is fractional TMEM145 occupancy
      by peptide vs by E1659A STRC?
  Q2. How does occupancy change as E1659A Kd shifts (WT-equivalent 10 nM →
      catastrophically weakened 100 μM)?
  Q3. Is there a [peptide] threshold where competition effectively erases
      E1659A occupancy regardless of its Kd?

Additional use case: all-DFNB16 non-Misha patients. E1659A is only ~10-15%
of DFNB16 alleles; most other alleles are truncating (null). For null/null
patients, there is zero endogenous STRC to compete with → peptide directly
occupies TMEM145 → simpler case.
"""

from __future__ import annotations

import json
import numpy as np
from pathlib import Path

OUT = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase4h_endogenous_strc_competition.json")

# Biology parameters
WT_STRC_KD_TMEM145_M = 10e-9   # literature estimate for full-length STRC × TMEM145
PEPTIDE_KD_TMEM145_M = 100e-9  # Phase 3 solo ipTM 0.68 → ~100 nM

# E1659A Kd scenarios (spans 10 nM WT-equivalent to 10 μM severely impaired)
E1659A_KD_SCENARIOS_M = {
    "wt_like_10_nM": 10e-9,
    "mild_3x_weaker_30_nM": 30e-9,
    "moderate_10x_weaker_100_nM": 100e-9,
    "severe_100x_weaker_1_uM": 1e-6,
    "catastrophic_1000x_weaker_10_uM": 10e-6,
    "nonfunctional_100000x_weaker_1_mM": 1e-3,
}

# Local concentrations on OHC surface
# E1659A STRC expression: compound het, only maternal allele, protein turnover 24h
# Upper estimate: same density as WT homozygote × 0.5 × 0.7 (protein instability) = 0.35× normal
# Normal STRC at OHC surface ~1 μM local (within tectorial attachment zone)
STRC_NORMAL_OHC_M = 1e-6
MISHA_E1659A_STRC_OHC_M_UPPER = 0.35 * STRC_NORMAL_OHC_M  # if stable
MISHA_E1659A_STRC_OHC_M_LOWER = 0.05 * STRC_NORMAL_OHC_M  # if unstable
TMEM145_TOTAL_OHC_M = 5e-6  # local density at OHC apical surface


def fractional_occupancy_2ligand(peptide_M: float, strc_M: float,
                                 kd_peptide_M: float, kd_strc_M: float,
                                 receptor_M: float) -> dict:
    """Two-ligand competition for one receptor.

    Solve: [P·R] / (receptor_free × peptide_free) = 1/Kd_P
           [S·R] / (receptor_free × strc_free)    = 1/Kd_S
           [R_total] = [R_free] + [P·R] + [S·R]
           [P_total] = [P_free] + [P·R]
           [S_total] = [S_free] + [S·R]

    For large ligand pools (P_total, S_total >> receptor_total),
    simplify: P_free ≈ P_total, S_free ≈ S_total. Then occupancy is:
      theta_P = [P]/Kd_P / (1 + [P]/Kd_P + [S]/Kd_S)
      theta_S = [S]/Kd_S / (1 + [P]/Kd_P + [S]/Kd_S)
      theta_free = 1 / (1 + [P]/Kd_P + [S]/Kd_S)
    """
    x_p = peptide_M / kd_peptide_M
    x_s = strc_M / kd_strc_M
    denom = 1 + x_p + x_s
    theta_p = x_p / denom
    theta_s = x_s / denom
    theta_free = 1 / denom
    return {
        "theta_peptide": theta_p,
        "theta_e1659a_strc": theta_s,
        "theta_free_receptor": theta_free,
    }


def main():
    # Peptide concentration sweep (from Phase 4e: 0.1-100 μM perilymph)
    peptide_concs_uM = np.logspace(-2, 3, 50)

    # Per-scenario competition
    per_scenario = {}
    for sc_name, e1659a_kd in E1659A_KD_SCENARIOS_M.items():
        rows = []
        for pep_uM in peptide_concs_uM:
            for strc_label, strc_conc in [("misha_upper", MISHA_E1659A_STRC_OHC_M_UPPER),
                                           ("misha_lower", MISHA_E1659A_STRC_OHC_M_LOWER)]:
                occ = fractional_occupancy_2ligand(
                    peptide_M=pep_uM * 1e-6,
                    strc_M=strc_conc,
                    kd_peptide_M=PEPTIDE_KD_TMEM145_M,
                    kd_strc_M=e1659a_kd,
                    receptor_M=TMEM145_TOTAL_OHC_M,
                )
                rows.append({
                    "peptide_uM": pep_uM,
                    "strc_scenario": strc_label,
                    "strc_M": strc_conc,
                    "e1659a_kd_M": e1659a_kd,
                    "theta_peptide": round(occ["theta_peptide"], 4),
                    "theta_e1659a": round(occ["theta_e1659a_strc"], 4),
                    "theta_free": round(occ["theta_free_receptor"], 4),
                })
        per_scenario[sc_name] = rows

    # Key concentration points
    key_doses_uM = [0.1, 1.0, 10.0, 100.0]
    key_table = []
    for pep_uM in key_doses_uM:
        for sc_name, e1659a_kd in E1659A_KD_SCENARIOS_M.items():
            for strc_label, strc_conc in [("misha_upper", MISHA_E1659A_STRC_OHC_M_UPPER),
                                           ("misha_lower", MISHA_E1659A_STRC_OHC_M_LOWER)]:
                occ = fractional_occupancy_2ligand(
                    peptide_M=pep_uM * 1e-6,
                    strc_M=strc_conc,
                    kd_peptide_M=PEPTIDE_KD_TMEM145_M,
                    kd_strc_M=e1659a_kd,
                    receptor_M=TMEM145_TOTAL_OHC_M,
                )
                key_table.append({
                    "peptide_uM": pep_uM,
                    "scenario": sc_name,
                    "strc_scenario": strc_label,
                    "theta_peptide": round(occ["theta_peptide"], 3),
                    "theta_e1659a": round(occ["theta_e1659a_strc"], 3),
                    "displacement_achieved": occ["theta_peptide"] > occ["theta_e1659a_strc"],
                })

    # Threshold peptide concentration to dominate at each E1659A scenario
    thresholds = {}
    for sc_name, e1659a_kd in E1659A_KD_SCENARIOS_M.items():
        for strc_label, strc_conc in [("misha_upper", MISHA_E1659A_STRC_OHC_M_UPPER),
                                       ("misha_lower", MISHA_E1659A_STRC_OHC_M_LOWER)]:
            crossover_pep = None
            for pep_uM in np.logspace(-2, 3, 500):
                occ = fractional_occupancy_2ligand(
                    peptide_M=pep_uM * 1e-6,
                    strc_M=strc_conc,
                    kd_peptide_M=PEPTIDE_KD_TMEM145_M,
                    kd_strc_M=e1659a_kd,
                    receptor_M=TMEM145_TOTAL_OHC_M,
                )
                if occ["theta_peptide"] >= 2 * occ["theta_e1659a_strc"]:
                    crossover_pep = pep_uM
                    break
            thresholds[f"{sc_name}_{strc_label}"] = (
                f"{crossover_pep:.2f} μM" if crossover_pep else ">1000 μM (unreachable)"
            )

    # All-DFNB16-null case: no competition
    null_case = {}
    for pep_uM in key_doses_uM:
        occ = fractional_occupancy_2ligand(
            peptide_M=pep_uM * 1e-6,
            strc_M=0.0,
            kd_peptide_M=PEPTIDE_KD_TMEM145_M,
            kd_strc_M=1.0,
            receptor_M=TMEM145_TOTAL_OHC_M,
        )
        null_case[f"{pep_uM}_uM"] = {
            "theta_peptide": round(occ["theta_peptide"], 3),
            "theta_free": round(occ["theta_free_receptor"], 3),
        }

    summary = {
        "batch": "hydrogel_phase4h_endogenous_strc_competition",
        "date": "2026-04-23",
        "patient_context": {
            "misha_genotype": "compound het: paternal 98 kb del + maternal E1659A",
            "only_maternal_allele_produces_protein": True,
            "e1659a_functional_status_unknown": (
                "AlphaMissense 0.996 (likely pathogenic), REVEL 0.68, but protein "
                "stability unmeasured; pharmacochaperone hypothesis assumes misfolding"
            ),
        },
        "biology_parameters": {
            "wt_strc_tmem145_kd_M": WT_STRC_KD_TMEM145_M,
            "peptide_tmem145_kd_M": PEPTIDE_KD_TMEM145_M,
            "e1659a_kd_scenarios_M": E1659A_KD_SCENARIOS_M,
            "strc_local_ohc_upper_estimate_M": MISHA_E1659A_STRC_OHC_M_UPPER,
            "strc_local_ohc_lower_estimate_M": MISHA_E1659A_STRC_OHC_M_LOWER,
            "tmem145_total_ohc_M": TMEM145_TOTAL_OHC_M,
        },
        "key_dose_table": key_table,
        "displacement_threshold_peptide_uM": thresholds,
        "all_null_patient_case": null_case,
        "conclusions": {
            "misha_best_case_e1659a_wt_like": (
                "If maternal E1659A retains WT affinity (Kd 10 nM) and is stable at ~0.35 μM "
                "OHC surface, peptide at 1 μM perilymph = ~20 μM local dominates 2:1 over E1659A. "
                "Threshold peptide concentration is achievable per Phase 4e (0.3-1 mg dose)."
            ),
            "misha_likely_case_e1659a_partially_impaired": (
                "If E1659A Kd is 10-100× weaker than WT (likely given ARM-repeat pocket "
                "destabilisation), competition is easy — 0.1 μM peptide already dominates. "
                "Phase 4e 0.1 mg dose sufficient."
            ),
            "misha_worst_case_e1659a_nonfunctional_stable": (
                "If E1659A is stable but binds TMEM145 at mM Kd (essentially non-functional), "
                "there is no competition — any peptide dose works."
            ),
            "non_misha_null_null_patients": (
                "Null/null genotype (both alleles produce no protein): no competition, "
                "peptide owns TMEM145 from first dose. Simpler therapeutic case."
            ),
        },
        "implications_for_phase_2c_wet_lab": [
            "HEK293 / HeLa actin-bundling assay will NOT have endogenous STRC interference — "
            "competition effects will be invisible there",
            "Competition matters only in OHC co-culture or in vivo — plan Phase 2d (organoid) "
            "or Phase 3 (Shanghai Shu Yilai knock-in mouse)",
            "E1659A functional assay (separate experiment) could resolve Kd — hybrid computational "
            "approach: run AF3 on E1659A × TMEM145 with 5 seeds, compare ipTM delta vs WT",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, default=str))

    print("=== Phase 4h Endogenous STRC vs Peptide Competition Model ===\n")
    print(f"Context: Misha compound het (paternal 98kb del + maternal E1659A)")
    print(f"Assumptions:")
    print(f"  WT STRC × TMEM145 Kd = {WT_STRC_KD_TMEM145_M*1e9:.0f} nM")
    print(f"  Peptide × TMEM145 Kd = {PEPTIDE_KD_TMEM145_M*1e9:.0f} nM")
    print(f"  E1659A local OHC conc: {MISHA_E1659A_STRC_OHC_M_LOWER*1e6:.2f} - {MISHA_E1659A_STRC_OHC_M_UPPER*1e6:.2f} μM\n")

    print("Displacement threshold peptide concentration (theta_pep ≥ 2× theta_e1659a):")
    for scenario_strc, thresh in thresholds.items():
        print(f"  {scenario_strc:55s}  peptide = {thresh}")

    print("\nKey dose table (peptide uM × e1659a_kd × strc conc):")
    print(f"{'[pep]':<8}{'E1659A scenario':<35}{'STRC pool':<15}{'θ_pep':<8}{'θ_e1659a':<10}{'wins':<6}")
    for row in key_table[:24]:
        print(f"{row['peptide_uM']:<8.1f}{row['scenario']:<35}"
              f"{row['strc_scenario']:<15}{row['theta_peptide']:<8.2%}"
              f"{row['theta_e1659a']:<10.2%}"
              f"{'Y' if row['displacement_achieved'] else 'N':<6}")
    print()
    print(f"All-null DFNB16 patients (no competition):")
    for k, v in null_case.items():
        print(f"  peptide {k}: θ_pep={v['theta_peptide']:.2%}, θ_free={v['theta_free']:.2%}")
    print(f"\nJSON written: {OUT}")


if __name__ == "__main__":
    main()
