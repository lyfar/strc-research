#!/usr/bin/env python3
"""
Hydrogel Phase 4d — F-actin bundling geometric + stoichiometric model.

CRITICAL GATE: WH2 is canonically a G-actin sequestering motif. For our
hypothesis #9 to work, the tail91 peptide-hydrogel must BUNDLE F-actin
filaments (side-by-side crosslinking), not sequester G-actin monomers.

The AF3 Phase 3b × G-actin trimer test showed ipTM 0.51 but is a poor
proxy for F-actin filament binding because:
  (a) G-actin trimer is not a filament — no pitch, no subdomain 1/3 cleft
      on the filament side surface.
  (b) WH2 binding to G-actin is high-affinity (Kd ~100 nM) but WH2 on
      F-actin side-binding is lower (Kd ~1-10 μM, if at all).
  (c) Canonical F-actin side-binders (cofilin, tropomyosin, ABPs like
      α-actinin, fimbrin/plastin, fascin, espin) have very different
      folds than WH2.

This script builds an analytical model to ask:

  Q1. Given WH2-G-actin Kd ~0.1-1 μM and peptide-hydrogel local
      concentration (when RADA16 self-assembles into a β-sheet ribbon),
      does F-actin bundling emerge?

  Q2. Does the geometric arrangement of multiple peptides on a RADA16
      fibril surface allow two F-actin filaments to be contacted
      simultaneously (the minimum requirement for bundling)?

  Q3. What is the minimum peptide concentration and RADA16 fibril density
      needed to achieve stereocilia-relevant F-actin organisation
      (inter-filament spacing ~12 nm, observed in normal stereocilia)?

Outputs:
  - Bundling propensity vs peptide concentration (mM) and RADA16 pitch
  - Predicted F-actin crosslink density (crosslinks / μm²)
  - Prediction for Phase 2c wet-lab readout (HEK293 GFP-actin, confocal
    image should show punctate co-localisation at peptide concentration
    > 10 μM if bundling mechanism operates)
  - Clear pass/fail gate for the hypothesis

Method: simple chemical physics — no Martini3 MD (that is Phase 5).
We model:
  (i)   WH2-G-actin binding kinetics (Kd = 100 nM literature for WASP WH2)
  (ii)  F-actin binding affinity adjusted for filament context (Kd = 5 μM
        estimated — weaker because of limited accessible surface)
  (iii) RADA16 β-sheet fibril geometry (2 peptides per 4.6 Å, so 4.35 peptides/nm
        along fibril axis, and ~4-6 fibrils cross-section for ~10 nm diameter)
  (iv)  F-actin filament geometry (13 subunits per 36 nm pitch, so 1 subunit
        every 2.77 nm along axis)
  (v)   Bundling requires ≥2 peptides contacting two different filaments
        within a persistence-length interval (~3 nm for stiff actin)
"""

from __future__ import annotations

import json
import numpy as np
from pathlib import Path

OUT = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase4d_factin_bundling_model.json")

# Physical constants
AVOGADRO = 6.022e23

# ═══════════════════════════════════════════════════════════════════════
# POST-AUDIT 2026-04-23 parameter updates — see literature-params/actin-kinetics.md
# ═══════════════════════════════════════════════════════════════════════

# WH2 × G-actin Kd — range 50 nM (WAVE) to 3 μM (WASP V-domain).
# Chereau 2005 PNAS (SI Table 2 was inaccessible; value inferred from
# narrative). Dominguez 2016 review confirms 50 nM-3 μM range for WH2 family.
# Padrick/Kim 2011 reports WASP V-domain × G-actin at 3.1 ± 0.5 μM (fluor.
# anisotropy) — primary. Current 200 nM is in the mid-range of WAVE-like.
WH2_KD_GACTIN_M = 200e-9   # Chereau 2005 PNAS + Dominguez 2016 review range

# WH2 × F-actin SIDE-binding: ⚠ NOT MEASURED in any published paper
# (literature search exhaustive, literature-params/actin-kinetics.md).
# Closest analog: Tβ4 × F-actin = 5-10 mM (1000-2000× weaker than WH2×G).
# The 5 μM here assumes a 25× weakening vs G-actin; this is OPTIMISTIC
# compared to Tβ4 analog. LOAD-BEARING RISK for h09: central therapeutic
# mechanism (F-actin side-bundling) requires this assumption + avidity
# from multi-WH2 RADA16 scaffold. Phase 2c wet-lab bundling assay is the
# only way to resolve.
WH2_KD_FACTIN_M = 5e-6     # ⚠ LOAD-BEARING unmeasured; Tβ4 analog = 5-10 mM

# RADA16 fibril geometry — Yokoi/Kinoshita 2005 PNAS + Zhang lab seminal
# β-sheet rise. 4.35 peptides/nm consistent with ~4.6 Å rise per pair.
RADA16_MONOMERS_PER_NM_FIBRIL = 4.35  # Yokoi/Kinoshita/Zhang 2005 PNAS
RADA16_FIBRIL_DIAMETER_NM = 10.0      # Yokoi 2005 AFM
RADA16_FIBRILS_PER_CROSSSECTION = 5   # Yokoi 2005 bundling estimate

# F-actin filament geometry — Oda & Namba 2009 Nature 457:441 cryo-EM
FACTIN_SUBUNITS_PER_NM = 13.0 / 36.0  # 0.361 subunits/nm
FACTIN_DIAMETER_NM = 7.0              # Oda & Namba 2009
FACTIN_PERSISTENCE_NM = 10000.0       # Isambert 1995 — standard textbook

# Stereocilia inter-filament spacing:
# Native OHC uses plastin-1/fimbrin → 7.9-9.7 nm (Krey 2016 J Cell Biol).
# Prior value 12 nm corresponds to ESPIN-crosslinked stereocilia — wrong
# paralog for OHC (OHC expresses plastin-1 + fimbrin, NOT espin dominantly).
# Updated to 9.0 nm (mid of Krey 2016 range).
STEREOCILIA_INTER_FILAMENT_NM = 9.0  # Krey 2016 J Cell Biol (was 12 espin)
STEREOCILIA_FILAMENT_DENSITY_PER_UM2 = 1 / (np.pi * (STEREOCILIA_INTER_FILAMENT_NM * 1e-3)**2) * 1.1
# With 9 nm spacing, density is ~4300 filaments/μm² (was ~2400 at 12 nm)

# Peptide concentration ranges (M)
peptide_conc_range = np.logspace(-7, -3, 50)  # 100 nM to 1 mM
rada16_fibril_pitch_range = np.linspace(2.0, 10.0, 9)  # peptides per nm (density)


def fraction_bound_wh2(peptide_conc_M: float, kd_M: float) -> float:
    """Fraction of actin sites occupied by WH2 at equilibrium.
    Simplified 1:1 binding: θ = [P] / (Kd + [P])."""
    return peptide_conc_M / (kd_M + peptide_conc_M)


def bundling_propensity(peptide_conc_M: float,
                        peptide_density_per_nm_fibril: float,
                        kd_facetin_M: float = WH2_KD_FACTIN_M) -> dict:
    """
    Bundling requires a RADA16 fibril with multiple peptides exposed on
    its surface, where each peptide has its WH2 hooked to a different
    F-actin filament on either side of the fibril.

    Probabilistic model:
      - Let p = fraction of peptides on fibril whose WH2 is actin-bound
        = fraction_bound_wh2(conc, Kd)
      - Let N = peptides in a 'bundling segment' = persistence length of
        F-actin × peptide_density_per_nm_fibril. Since actin is extremely
        stiff, we limit by RADA16 fibril length scale (~1 μm = 1000 nm).
        Effective N for a ~10 nm decision volume (2× actin diameter):
          N_eff = 10 nm × peptide_density_per_nm
      - Two filaments on opposite sides = 2 independent peptides, each
        with probability p.
      - Bundling probability per 10 nm segment: p² (two peptides, each
        independently bound).
      - Crosslinks per μm of fibril = p² × peptide_density_per_nm × 1000
    """
    p_bound = fraction_bound_wh2(peptide_conc_M, kd_facetin_M)
    # Effective peptides per decision volume
    decision_volume_nm = 10.0  # 2× actin diameter
    N_eff = decision_volume_nm * peptide_density_per_nm_fibril
    # Probability at least 2 simultaneously bound (binomial with N_eff trials)
    # P(X >= 2) = 1 - P(0) - P(1)
    if N_eff < 2:
        prob_2plus = 0.0
    else:
        p0 = (1 - p_bound) ** N_eff
        p1 = N_eff * p_bound * (1 - p_bound) ** (N_eff - 1)
        prob_2plus = 1 - p0 - p1
    # Crosslinks per μm of fibril (fibril length × crosslink events)
    fibril_length_nm = 1000.0  # 1 μm typical
    expected_crosslinks_per_um = prob_2plus * (fibril_length_nm / decision_volume_nm)
    return {
        "p_bound": round(p_bound, 4),
        "N_eff_peptides_per_10nm": round(N_eff, 2),
        "prob_two_plus_bound": round(prob_2plus, 4),
        "expected_crosslinks_per_um_fibril": round(expected_crosslinks_per_um, 2),
    }


def bundling_heatmap():
    """Generate heatmap of bundling propensity across conditions."""
    rows = []
    for conc in peptide_conc_range:
        row = []
        for pitch in rada16_fibril_pitch_range:
            bp = bundling_propensity(conc, pitch)
            row.append(bp["expected_crosslinks_per_um_fibril"])
        rows.append(row)
    return rows


def minimum_concentration_for_rescue(kd_factin_M: float = WH2_KD_FACTIN_M,
                                     target_crosslink_density: float = 10.0) -> float:
    """Find peptide concentration at which 10 crosslinks/μm are expected
    (empirical threshold for F-actin bundle formation in vitro)."""
    for conc in peptide_conc_range:
        bp = bundling_propensity(conc, RADA16_MONOMERS_PER_NM_FIBRIL)
        if bp["expected_crosslinks_per_um_fibril"] >= target_crosslink_density:
            return conc
    return None


def actin_monomer_sequestration_fraction(peptide_conc_M: float,
                                         gactin_pool_M: float = 50e-6) -> float:
    """What fraction of free G-actin pool is sequestered by WH2?
    OHC cytosolic G-actin pool ~50 μM. If we sequester too much, we
    starve F-actin polymerisation — a BAD outcome, we want filament
    bundling, not monomer depletion."""
    # Simple competition: at equilibrium, [P·G] / [G_free] = [P_free] / Kd
    Kd = WH2_KD_GACTIN_M
    # Quadratic: [PG]² - ([P0]+[G0]+Kd)[PG] + [P0][G0] = 0
    P0 = peptide_conc_M; G0 = gactin_pool_M
    b = P0 + G0 + Kd
    disc = b*b - 4*P0*G0
    PG = (b - np.sqrt(max(disc, 0))) / 2
    return PG / G0 if G0 > 0 else 0.0


def main():
    # Sweep conditions
    heatmap = bundling_heatmap()

    # Key dose points (ototopical delivery realism)
    key_concs = {
        "100 nM (literature Kd)": 100e-9,
        "1 μM (clinical-minimum floor)": 1e-6,
        "10 μM (Phase 2c wet-lab default)": 10e-6,
        "100 μM (topical max before solubility limit)": 100e-6,
        "1 mM (supersaturation)": 1e-3,
    }
    key_results = {}
    for label, conc in key_concs.items():
        key_results[label] = {
            "concentration_M": conc,
            "bundling_at_average_density": bundling_propensity(
                conc, RADA16_MONOMERS_PER_NM_FIBRIL),
            "bundling_at_max_density": bundling_propensity(
                conc, RADA16_MONOMERS_PER_NM_FIBRIL * RADA16_FIBRILS_PER_CROSSSECTION),
            "gactin_sequestration_fraction": round(
                actin_monomer_sequestration_fraction(conc), 4),
        }

    min_conc = minimum_concentration_for_rescue()

    # Competing-affinity analysis: what if our WH2-F-actin Kd is 10x
    # weaker than literature estimate?
    sensitivity = {}
    for kd_mult, label in [(0.1, "optimistic (stronger binder)"),
                           (1.0, "nominal"),
                           (10.0, "pessimistic (weaker binder)")]:
        kd = WH2_KD_FACTIN_M * kd_mult
        sensitivity[label] = {
            "kd_factin_M": kd,
            "bundling_at_10uM": bundling_propensity(10e-6, RADA16_MONOMERS_PER_NM_FIBRIL,
                                                    kd_facetin_M=kd),
        }

    # Minimum peptide conc needed for stereocilia-relevant bundling density
    # Target: 10 crosslinks per μm of RADA16 fibril → produces ~1-3 F-actin
    # bundles per μm² at reasonable fibril density

    # Decision logic
    bundling_pass = any(
        k["bundling_at_average_density"]["expected_crosslinks_per_um_fibril"] >= 10
        for k in key_results.values()
    )
    clinical_pass = any(
        conc <= 100e-6 and  # achievable topically
        key_results[label]["bundling_at_average_density"]["expected_crosslinks_per_um_fibril"] >= 10
        for label, conc in key_concs.items()
    )
    sequestration_safe = all(
        v["gactin_sequestration_fraction"] < 0.3  # <30% of G-actin pool
        for v in key_results.values() if v["concentration_M"] <= 100e-6
    )

    if clinical_pass and sequestration_safe:
        verdict = "PASS — F-actin bundling achievable at ≤100 μM with safe G-actin depletion."
    elif clinical_pass and not sequestration_safe:
        verdict = "PARTIAL — bundling works but G-actin sequestration exceeds 30% at needed dose. Reduce WH2 affinity or use F-actin-only binder."
    elif not clinical_pass:
        verdict = "FAIL — clinically achievable concentrations don't produce sufficient bundling. Consider F-actin side-binder (espin/fascin-like) instead of WH2."

    summary = {
        "batch": "hydrogel_phase4d_factin_bundling_model",
        "date": "2026-04-23",
        "parameters": {
            "WH2_Kd_Gactin_M": WH2_KD_GACTIN_M,
            "WH2_Kd_Factin_M_estimated": WH2_KD_FACTIN_M,
            "RADA16_peptides_per_nm_fibril": RADA16_MONOMERS_PER_NM_FIBRIL,
            "stereocilia_target_crosslinks_per_um": 10,
            "gactin_pool_ohc_M": 50e-6,
        },
        "key_concentration_results": key_results,
        "minimum_concentration_for_10_crosslinks_per_um_M": min_conc,
        "sensitivity_to_Kd_estimate": sensitivity,
        "heatmap_crosslinks_per_um": {
            "conc_M_axis": peptide_conc_range.tolist(),
            "peptide_density_per_nm_axis": rada16_fibril_pitch_range.tolist(),
            "values": heatmap,
        },
        "bundling_gate_pass": bundling_pass,
        "clinically_achievable_gate_pass": clinical_pass,
        "sequestration_safe_gate_pass": sequestration_safe,
        "verdict": verdict,
        "wet_lab_predictions_phase_2c": {
            "HEK293_GFP_actin_10uM_peptide": (
                "Punctate co-localisation expected in fluorescence microscopy within 30 min. "
                f"Bundles should form at crosslink density ~{key_results['10 μM (Phase 2c wet-lab default)']['bundling_at_average_density']['expected_crosslinks_per_um_fibril']}/μm. "
                "Negative control: WH2-ablated variant (Phase 4b job 6) should show diffuse actin."
            ),
            "TIRF_single_filament_10uM_peptide": (
                "F-actin filaments should bundle into thick cables within 5-15 min. "
                "Measurable by TIRF imaging of Alexa647-actin."
            ),
            "critical_concentration_for_bundling": f"{min_conc * 1e6:.1f} μM" if min_conc else "not achievable below 1 mM",
        },
        "risk_flags": [
            "WH2 is canonically a G-ACTIN sequestering motif, not an F-actin side-binder. Our model assumes WH2 retains weak F-actin affinity (Kd 5 μM) — this must be verified in Phase 2c.",
            "If WH2 truly cannot bind F-actin side, hypothesis fails regardless of ipTM — switch to fascin/espin/plastin-derived F-actin binder.",
            "G-actin monomer sequestration at 100 μM peptide dose: ~{:.1%} of free G-actin pool — could stress OHC actin turnover. Monitor via phalloidin staining in Phase 2c.".format(
                actin_monomer_sequestration_fraction(100e-6)),
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, default=str))

    print("=== Phase 4d F-actin Bundling Geometric Model ===\n")
    print(f"Parameters: WH2-G-actin Kd = {WH2_KD_GACTIN_M*1e9:.0f} nM")
    print(f"            WH2-F-actin Kd (est) = {WH2_KD_FACTIN_M*1e6:.0f} μM")
    print(f"            RADA16 density = {RADA16_MONOMERS_PER_NM_FIBRIL:.2f} peptides/nm fibril\n")
    print(f"Key concentration results:")
    for label, v in key_results.items():
        cl = v["bundling_at_average_density"]["expected_crosslinks_per_um_fibril"]
        gs = v["gactin_sequestration_fraction"]
        print(f"  {label:45s}  crosslinks/μm={cl:5.2f}, G-actin seq={gs:.2%}")
    print()
    print(f"Minimum concentration for rescue-relevant bundling: {min_conc*1e6:.1f} μM" if min_conc else "NOT ACHIEVABLE below 1 mM")
    print()
    print(f"Sensitivity to Kd estimate (at 10 μM peptide):")
    for label, v in sensitivity.items():
        print(f"  {label:30s}  Kd={v['kd_factin_M']*1e6:.1f} μM, crosslinks/μm={v['bundling_at_10uM']['expected_crosslinks_per_um_fibril']:.2f}")
    print()
    print(f">>> VERDICT: {verdict}")
    print(f"\nRisk flags:")
    for r in summary["risk_flags"]:
        print(f"  - {r}")


if __name__ == "__main__":
    main()
