"""
h09 Phase 4e_v2 — Blend-strategy scaffold geometry + TMEM145 valency.

MOTIVATION
----------
h09 Phase 3b validated a 118 aa RADA16-tail91 construct in AF3 (ipTM 0.68).
BUT: empirical literature (Frontiers 2021 RADA16 review, PMC8216384) caps
preserved β-sheet assembly at ≤12 residues of appended sequence; 12–50
residues show monotonic viscoelasticity decline; >50 aa never tested. Our
tail91 (118 aa) sits ~10× beyond the characterized limit — assembly at 100%
tail91 composition is NOT a safe assumption.

Standard workaround (Gelain osteogenic hybrids, Zhang group): **molar blend**
of plain RADA16 : RADA16-tail-fusion at 8:2 to 9.5:0.5 ratio. Plain RADA16
carries β-sheet scaffold; tail displays binding motifs on surface. Drawback:
every dose/valency calculation in current Phase 4d/4e assumes 100% tail91
contributes — at 9:1 blend, only 10% of mass is functional.

This script bounds the blend-strategy design: at each molar fraction
phi_tail ∈ [0.05, 0.50], compute:
  1. Fibril-level tail density (tails per nm, tails per fibril)
  2. Bundle-level tail density (bundles per stereocilia unit area)
  3. Cochlear-volume-level tail count at standard gel concentrations
  4. Comparison to putative TMEM145 site density at stereocilia
  5. Assembly feasibility flag (based on Frontiers 2021 empirical limits)
  6. Dose implication: how much blend peptide is needed to match the
     functional tail91 mass assumed in current Phase 4d/4e

LITERATURE-ANCHORED PARAMETERS (from literature-params/rada16-geometry.md)
  RADA16 monomers/nm fibril axis:  4.26  (Paravastu 2014 + Cormier 2012)
  Mean fibril length:              615 ± 104 nm  (Yokoi 2005)
  Fibril cross-section (tape):     3 nm × 2.5 nm  (Paravastu 2014)
  Bundle diameter:                 10-20 nm  (Zhang 1993)
  Fibers per bundle (derived):     3-7
  Empirical tail limit (assembled):≤12 residues  (Frontiers 2021)
  Monotonic decline:               12–50 residues
  Untested regime:                 >50 residues  (our tail91 = 118 aa)
  Standard gel concentration:      1-2% w/v → 6-12 mM monomer
  Stereocilia spacing:             9 nm  (Krey 2016 plastin)
  Perilymph volume:                34 μL  (Dhanasingh 2021)
  Stereocilia bundle geometry:     ~30 per OHC, ~5 μm × 0.5 μm each

UNMEASURED / PLACEHOLDER (flagged in output):
  TMEM145 sites per stereocilium — never quantified in primary lit. Order
    estimate used: 10^3 per stereocilium (geometric capacity at 9 nm spacing
    over ~5 μm² membrane area). This is a CAPACITY upper bound, not a
    measured density.
  Assembly feasibility at phi_tail > 0.05 — extrapolation; Frontiers 2021
    only characterizes ≤50 aa tails in pure (not blended) constructs.

OUTPUT
------
  JSON table + summary. No PNG.

RUNTIME: <1 s.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np

ROOT = pathlib.Path("/Users/egorlyfar/Brain/research/strc/models")

# ─────────────── Literature-anchored constants ────────────────────
RADA16_MONOMERS_PER_NM = 4.26                 # Paravastu 2014 + Cormier 2012
FIBRIL_LENGTH_MEAN_NM = 615.0                 # Yokoi 2005
FIBRIL_LENGTH_SD_NM = 104.0                   # Yokoi 2005
FIBRIL_TAPE_WIDTH_NM = 3.0                    # Paravastu 2014
FIBRIL_TAPE_HEIGHT_NM = 2.5                   # Paravastu 2014
BUNDLE_DIAMETER_NM_MEAN = 15.0                # Zhang 1993, 10-20 nm range
FIBERS_PER_BUNDLE_MEAN = 5.0                  # 10-20 nm / 3 nm, range 3-7
RADA16_MW_DA = 1712.0                         # Ac-(RADA)4-CONH2
TAIL91_MW_DA = 1712.0 + 91 * 110.0            # + 91 aa @ avg 110 Da = ~11.7 kDa total

# Empirical assembly limits (Frontiers 2021 review, PMC8216384)
TAIL_LIMIT_PRESERVED_AA = 12                  # ≤12 aa: β-sheet preserved
TAIL_LIMIT_MONOTONIC_AA = 50                  # 12–50 aa: monotonic decline
TAIL91_LENGTH_AA = 91                         # our construct

# Stereocilia geometry (Krey 2016 + standard hair-cell anatomy)
STEREOCILIA_SPACING_NM = 9.0                  # Krey 2016 plastin crosslink
STEREOCILIA_LENGTH_UM = 5.0                   # OHC tallest row (species-avg)
STEREOCILIA_DIAM_UM = 0.5                     # tallest row diameter
STEREOCILIA_PER_OHC = 30.0                    # tallest row count
OHC_TOTAL = 12000.0                           # human cochlea

# Cochlear delivery (literature-backed post-audit 2026-04-23)
PERILYMPH_VOL_UL = 34.0                       # Dhanasingh 2021

# ─────────────── Blend grid ───────────────────────────────────────
PHI_TAIL_GRID = [0.01, 0.05, 0.10, 0.20, 0.33, 0.50]
GEL_CONC_PCT_WV_GRID = [1.0, 2.0]             # standard RADA16 range


# ─────────────── Derived quantities ───────────────────────────────
def monomers_per_fibril(fibril_length_nm: float = FIBRIL_LENGTH_MEAN_NM) -> float:
    return RADA16_MONOMERS_PER_NM * fibril_length_nm


def gel_concentration_mM(gel_pct_wv: float, phi_tail: float) -> dict:
    """Convert % w/v to mM, split into plain + tail91 fractions."""
    # Total mass-fraction: 1% w/v = 10 mg/mL. Mass-average MW depends on phi.
    mass_avg_mw = (1.0 - phi_tail) * RADA16_MW_DA + phi_tail * TAIL91_MW_DA
    tot_mg_per_ml = gel_pct_wv * 10.0
    total_monomer_mM = tot_mg_per_ml / mass_avg_mw * 1000.0
    return {
        "total_monomer_mM": total_monomer_mM,
        "plain_monomer_mM": (1.0 - phi_tail) * total_monomer_mM,
        "tail91_monomer_mM": phi_tail * total_monomer_mM,
    }


def assembly_feasibility(phi_tail: float) -> dict:
    """Heuristic feasibility flag based on Frontiers 2021 empirical trends.

    Note: Frontiers 2021 characterizes pure constructs with variable tail
    length, NOT blends with constant tail but variable molar ratio. We
    extrapolate: the disruption scales with the *fraction of β-sheet
    interfaces occupied by a >50 aa tail*, i.e., ≈ phi_tail.
    """
    if phi_tail <= 0.05:
        flag = "PASS"
        note = "≤5% tail-bearing monomer; β-sheet scaffold dominated by plain RADA16. Standard Gelain-range blend."
    elif phi_tail <= 0.15:
        flag = "WARN"
        note = "5-15%: within reported hybrid range for osteogenic/RGD RADA16 (Gelain group) but tail >50 aa is untested — viscoelastic decline likely."
    elif phi_tail <= 0.33:
        flag = "RISK"
        note = "15-33%: no literature precedent for >50 aa tail at this loading; expect significant viscoelastic decline and possible phase separation."
    else:
        flag = "FAIL"
        note = ">33% tail91: extrapolation unsupported; assembly likely collapses per Frontiers 2021 ALK-tag analog."
    return {"flag": flag, "note": note}


def tails_per_fibril(phi_tail: float, fibril_length_nm: float = FIBRIL_LENGTH_MEAN_NM) -> dict:
    mon = monomers_per_fibril(fibril_length_nm)
    tails = mon * phi_tail
    # If tail91 is randomly distributed, linear density:
    tails_per_nm = RADA16_MONOMERS_PER_NM * phi_tail
    # Mean gap between tails along fibril axis
    gap_nm = 1.0 / tails_per_nm if tails_per_nm > 0 else np.inf
    return {
        "tails_per_fibril_mean": tails,
        "tails_per_nm_axis": tails_per_nm,
        "mean_gap_along_axis_nm": gap_nm,
    }


def tails_per_bundle_length(phi_tail: float) -> dict:
    """Tails per nm along a bundle (contains FIBERS_PER_BUNDLE fibrils in parallel)."""
    per_nm = RADA16_MONOMERS_PER_NM * phi_tail * FIBERS_PER_BUNDLE_MEAN
    return {
        "tails_per_nm_bundle_axis": per_nm,
        "tails_per_um_bundle_axis": per_nm * 1000.0,
    }


def tails_per_cm3_gel(phi_tail: float, gel_pct_wv: float) -> float:
    """Total tail91 count per cm³ at given gel concentration."""
    conc = gel_concentration_mM(gel_pct_wv, phi_tail)
    tails_mM = conc["tail91_monomer_mM"]
    # mM → molec/cm³:   mM = 10^-3 mol/L = 10^-6 mol/cm³
    #                   × 6.022e23 mol^-1 = 6.022e17 / cm³ per mM
    return tails_mM * 6.022e17


def valency_at_stereocilia_bundle(phi_tail: float, gel_pct_wv: float,
                                   gel_thickness_at_bundle_um: float = 1.0) -> dict:
    """How many tail91 sites are within the gel shell surrounding one OHC bundle?

    OHC bundle approximation: 30 stereocilia, each cylinder 5 μm × 0.5 μm.
    Outer surface area ≈ 30 × π × 0.5 × 5 = 235 μm² per bundle.
    Gel shell volume at 1 μm thickness ≈ 235 × 1 μm³ = 2.35 × 10^-10 cm³.
    """
    shell_area_um2 = STEREOCILIA_PER_OHC * np.pi * STEREOCILIA_DIAM_UM * STEREOCILIA_LENGTH_UM
    shell_vol_um3 = shell_area_um2 * gel_thickness_at_bundle_um
    shell_vol_cm3 = shell_vol_um3 * 1e-12
    tails_density = tails_per_cm3_gel(phi_tail, gel_pct_wv)
    tails_in_shell = tails_density * shell_vol_cm3
    # Putative TMEM145 site count per bundle (CAPACITY upper bound — not measured):
    #   30 stereocilia × (5 μm length / 9 nm spacing)² × geometric packing ≈
    #   30 × (5000/9)^2  ≈ 30 × 3.1e5  ≈ 9.3e6 sites/bundle MAX
    # More conservative: ~10^3-10^4 per stereocilium (membrane-receptor order)
    tmem145_sites_estimate_per_bundle_low = STEREOCILIA_PER_OHC * 1e3
    tmem145_sites_estimate_per_bundle_high = STEREOCILIA_PER_OHC * 1e4
    return {
        "shell_area_per_bundle_um2": shell_area_um2,
        "shell_vol_per_bundle_um3": shell_vol_um3,
        "tails_in_shell_per_bundle": tails_in_shell,
        "tmem145_sites_per_bundle_low_estimate": tmem145_sites_estimate_per_bundle_low,
        "tmem145_sites_per_bundle_high_estimate": tmem145_sites_estimate_per_bundle_high,
        "valency_ratio_tails_per_TMEM145_low": tails_in_shell / tmem145_sites_estimate_per_bundle_low,
        "valency_ratio_tails_per_TMEM145_high": tails_in_shell / tmem145_sites_estimate_per_bundle_high,
        "note": ("TMEM145 sites per bundle are order-estimate (10^3-10^4/ster). "
                 "Primary lit measurement missing — gate to S-tier."),
    }


def dose_implication_vs_100pct_tail(phi_tail: float, target_tail_mass_ug: float) -> dict:
    """At blend ratio phi_tail, how much total peptide is needed to deliver
    the same tail91 mass as a hypothetical 100% tail91 dose?"""
    # Mass fraction of tail91 in blend:
    mass_avg_mw = (1.0 - phi_tail) * RADA16_MW_DA + phi_tail * TAIL91_MW_DA
    tail_mass_fraction = (phi_tail * TAIL91_MW_DA) / mass_avg_mw
    total_peptide_mass_ug = target_tail_mass_ug / tail_mass_fraction
    return {
        "mass_fraction_tail91": tail_mass_fraction,
        "total_peptide_mass_ug_needed": total_peptide_mass_ug,
        "ratio_vs_100pct": total_peptide_mass_ug / target_tail_mass_ug,
    }


def cochlea_level_coverage(phi_tail: float, gel_pct_wv: float,
                           gel_volume_ul: float = 10.0) -> dict:
    """Given gel volume injected into perilymph, how does scaffold distribute?

    10 μL gel in 34 μL perilymph → ~22% of perilymph volume occupied by gel
    initially. Scaffold will localize near injection + migrate via diffusion
    + gel dissolution (t½ 7-14 days per Frontiers 2021).

    Upper bound: all 12000 OHCs contacted by gel scaffold.
    """
    total_tails_delivered = tails_per_cm3_gel(phi_tail, gel_pct_wv) * gel_volume_ul * 1e-3  # μL → cm³ ×1e-3
    per_OHC_if_uniform = total_tails_delivered / OHC_TOTAL
    return {
        "gel_volume_ul": gel_volume_ul,
        "total_tail91_molecules_delivered": total_tails_delivered,
        "per_OHC_if_uniform_distribution": per_OHC_if_uniform,
        "per_stereocilium_if_uniform": per_OHC_if_uniform / STEREOCILIA_PER_OHC,
        "note": "Uniform distribution is optimistic; real delivery has apico-basal gradient (Salt & Plontke 2018).",
    }


def main() -> None:
    rows: list[dict] = []
    for phi_tail in PHI_TAIL_GRID:
        for gel_c in GEL_CONC_PCT_WV_GRID:
            feas = assembly_feasibility(phi_tail)
            conc = gel_concentration_mM(gel_c, phi_tail)
            fibril = tails_per_fibril(phi_tail)
            bundle = tails_per_bundle_length(phi_tail)
            sterio = valency_at_stereocilia_bundle(phi_tail, gel_c)
            cover = cochlea_level_coverage(phi_tail, gel_c)
            # Reference comparison: 100% tail91 mass needed (hypothetical) = 100 μg
            # ("standard" Phase 4e injection). Quote ratio.
            dose_100ug = dose_implication_vs_100pct_tail(phi_tail, 100.0)
            rows.append({
                "phi_tail": phi_tail,
                "gel_pct_wv": gel_c,
                "assembly_flag": feas["flag"],
                "assembly_note": feas["note"],
                "concentration": conc,
                "fibril_level": fibril,
                "bundle_level": bundle,
                "stereocilia_valency": sterio,
                "cochlea_coverage": cover,
                "dose_to_match_100ug_tail91": dose_100ug,
            })

    # Headline summary: recommended blend
    recommended = [r for r in rows
                   if r["assembly_flag"] in ("PASS", "WARN") and r["gel_pct_wv"] == 2.0]
    best = max(recommended,
               key=lambda r: r["stereocilia_valency"]["tails_in_shell_per_bundle"]
               if r["assembly_flag"] == "PASS" else -1)

    out = {
        "model": "h09 Phase 4e_v2 blend-strategy scaffold geometry + valency",
        "literature_anchors": {
            "RADA16_monomers_per_nm": RADA16_MONOMERS_PER_NM,
            "fibril_length_mean_nm": FIBRIL_LENGTH_MEAN_NM,
            "empirical_tail_limit_preserved_aa": TAIL_LIMIT_PRESERVED_AA,
            "empirical_tail_limit_monotonic_aa": TAIL_LIMIT_MONOTONIC_AA,
            "our_tail_length_aa": TAIL91_LENGTH_AA,
            "perilymph_vol_ul_Dhanasingh_2021": PERILYMPH_VOL_UL,
            "stereocilia_spacing_nm_Krey_2016": STEREOCILIA_SPACING_NM,
        },
        "flags_and_caveats": {
            "tail91_118aa_vs_empirical_limit_12aa": (
                "Our construct is ~10× beyond characterized safe limit. "
                "Blend at phi_tail ≤ 0.15 is only defensible design."),
            "TMEM145_sites_per_bundle_unmeasured": (
                "All valency ratios depend on 10^3-10^4/stereocilium estimate. "
                "Primary lit measurement is the S-tier promotion gate."),
            "assembly_feasibility_extrapolated": (
                "Frontiers 2021 does not characterize blends with >50 aa tails; "
                "feasibility flags are extrapolations."),
        },
        "grid_rows": rows,
        "recommended_blend": {
            "phi_tail": best["phi_tail"],
            "gel_pct_wv": best["gel_pct_wv"],
            "assembly_flag": best["assembly_flag"],
            "tails_in_shell_per_bundle": best["stereocilia_valency"]["tails_in_shell_per_bundle"],
            "valency_ratio_low_est": best["stereocilia_valency"]["valency_ratio_tails_per_TMEM145_low"],
            "valency_ratio_high_est": best["stereocilia_valency"]["valency_ratio_tails_per_TMEM145_high"],
            "dose_to_match_100ug_pure_tail91": best["dose_to_match_100ug_tail91"]["total_peptide_mass_ug_needed"],
        },
    }
    out_path = ROOT / "hydrogel_phase4e_v2_blend_scaffold.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"wrote {out_path}")
    print()
    print("Blend grid (φ_tail × gel%):")
    print(f"{'φ_tail':>8} {'gel%':>5} {'flag':>6} {'tails/fibril':>14} "
          f"{'tails/bundle':>14} {'valency-low':>14} {'dose-vs-100%':>14}")
    for r in rows:
        sv = r["stereocilia_valency"]
        print(f"{r['phi_tail']:>8.2f} {r['gel_pct_wv']:>4.1f}% {r['assembly_flag']:>6s} "
              f"{r['fibril_level']['tails_per_fibril_mean']:>14.1f} "
              f"{sv['tails_in_shell_per_bundle']:>14.2e} "
              f"{sv['valency_ratio_tails_per_TMEM145_low']:>14.2f} "
              f"{r['dose_to_match_100ug_tail91']['ratio_vs_100pct']:>14.1f}×")
    print()
    print(f"Recommended: φ_tail={best['phi_tail']}, gel {best['gel_pct_wv']}% w/v "
          f"({best['assembly_flag']}).")
    print(f"  Tails per OHC bundle shell: {best['stereocilia_valency']['tails_in_shell_per_bundle']:.2e}")
    print(f"  Valency ratio vs TMEM145 (low est): "
          f"{best['stereocilia_valency']['valency_ratio_tails_per_TMEM145_low']:.2f}")
    print(f"  Dose cost vs 100% tail91: "
          f"{best['dose_to_match_100ug_tail91']['ratio_vs_100pct']:.1f}×")


if __name__ == "__main__":
    main()
