#!/usr/bin/env python3
"""
Sonoporation-Mediated Delivery Model for STRC Gene Therapy
===========================================================
Models the physics of non-invasive delivery through the round window membrane (RWM)
using ultrasound + microbubbles (sonoporation) to deliver LNP-packaged mini-STRC mRNA
to outer hair cells (OHCs).

Three sub-models:
1. Sonoporation pore formation in RWM (pore size, density, lifetime)
2. Nanoparticle diffusion through transient pores into perilymph  
3. LNP uptake by OHCs and comparison with standard AAV injection

Literature sources:
- Zhou et al. 2009 (PMC2752487): pore radius ~110 nm ± 40 nm
- Shih et al. 2019 (PMC7823126): RWM sonoporation in guinea pigs
  * US: 1 MHz, 3 W/cm², MI=0.254, 50% duty cycle, 1-min courses
  * MBs: SonoVue, 2-5 × 10⁸ bubbles/mL, 200 µL volume
  * Enhanced permeability: 3.1× (3 courses) to 5.2× (5 courses)
  * Recovery: full RWM integrity restored within 24h
  * Safety: no ABR threshold shift, no hair cell damage
- Landegger et al. 2017 (PMC5340646): AAV-Anc80L65 OHC transduction 60-100%
- Iversen et al. 2022 (Nat Commun): NHP AAV cochlear delivery, Anc80L65 up to 90% IHC

Author: Egor Lyfar (computational, not experimentally validated)
Date: 2026-03-17
"""

import numpy as np
import json
from dataclasses import dataclass, asdict

# ============================================================
# PHYSICAL CONSTANTS
# ============================================================
k_B = 1.381e-23      # Boltzmann constant [J/K]
T = 310.15            # Body temperature [K] (37°C)
eta_perilymph = 0.72e-3  # Perilymph viscosity [Pa·s] (close to water at 37°C)
N_A = 6.022e23        # Avogadro's number

# ============================================================
# MODEL 1: SONOPORATION PORE FORMATION IN RWM
# ============================================================
@dataclass
class SonoporationParams:
    """Parameters for ultrasound + microbubble sonoporation of RWM"""
    # Ultrasound parameters (from Shih et al. 2019)
    frequency_MHz: float = 1.0          # [MHz]
    intensity_W_cm2: float = 3.0        # [W/cm²] spatial average
    mechanical_index: float = 0.254     # MI = P_neg / sqrt(f_MHz)
    duty_cycle: float = 0.5             # 50%
    pulse_duration_ms: float = 2.0      # burst duration [ms]
    pulse_rate_Hz: float = 250          # burst rate [Hz]
    exposure_courses: int = 5           # number of 1-min courses
    exposure_time_per_course_s: float = 60.0  # [s]
    
    # Microbubble parameters (SonoVue)
    mb_concentration: float = 3.5e8     # [bubbles/mL] (midpoint 2-5 × 10⁸)
    mb_diameter_um: float = 2.5         # mean diameter [µm] (SonoVue 1.5-4 µm)
    mb_volume_uL: float = 200.0        # volume applied to bulla [µL]
    
    # RWM geometry (from guinea pig, scaled to human)
    rwm_thickness_um: float = 70.0      # human RWM thickness [µm] (GP: 10-40, human: 55-100)
    rwm_area_mm2: float = 2.0           # human RWM area [mm²] (GP: 1.1, human: 2.0-2.3)
    
    # Pore parameters (from Zhou et al. 2009 voltage clamp measurements)
    pore_radius_nm: float = 110.0       # mean pore radius [nm]
    pore_radius_std_nm: float = 40.0    # std dev [nm]
    pore_lifetime_s: float = 5.0        # single-cell pores reseal in ~seconds
    # For tissue (RWM): pores remain open longer due to multi-layer structure
    rwm_pore_lifetime_min: float = 30.0  # enhanced permeability persists ~30 min
    
    # Permeability enhancement (from Shih et al. 2019, Figure 5)
    perm_enhancement_3course: float = 3.1  # fold increase after 3 courses
    perm_enhancement_5course: float = 5.2  # fold increase after 5 courses


@dataclass
class LNPParams:
    """Lipid nanoparticle parameters"""
    diameter_nm: float = 80.0           # LNP diameter [nm] (typical 60-100)
    mRNA_copies_per_LNP: int = 4        # mRNA molecules per LNP (typical 1-7)
    concentration_particles_mL: float = 1e12  # typical LNP concentration
    volume_applied_uL: float = 50.0     # volume applied to RWM
    endosomal_escape_pct: float = 2.0   # % of LNPs that escape endosome (1-4%)
    mRNA_translation_proteins: int = 100 # proteins per mRNA before degradation


@dataclass
class CochlearParams:
    """Inner ear geometry and cell parameters"""
    # Perilymph volume
    scala_tympani_volume_uL: float = 30.0  # human scala tympani [µL]
    
    # Hair cell counts (human)
    n_OHC_total: int = 12000           # total outer hair cells
    n_IHC_total: int = 3500            # total inner hair cells
    
    # OHC parameters
    ohc_surface_area_um2: float = 500   # apical surface area [µm²]
    ohc_volume_pL: float = 1.0          # cell volume [pL]
    
    # Stereocilin requirement
    stereocilin_molecules_per_OHC: int = 15000  # from our NFAT model
    stereocilin_mw_kDa: float = 131.0   # mini-STRC molecular weight

    # AAV comparison (Landegger 2017, Iversen 2022)
    aav_ohc_transduction_pct: float = 80.0  # Anc80L65 OHC transduction (60-100%)
    aav_ihc_transduction_pct: float = 95.0  # Anc80L65 IHC transduction (near 100%)


# ============================================================
# MODEL 2: DIFFUSION THROUGH SONOPORATION PORES
# ============================================================

def stokes_einstein_diffusion(particle_diameter_nm, viscosity=eta_perilymph, temp=T):
    """Stokes-Einstein diffusion coefficient for spherical particle"""
    r = particle_diameter_nm * 1e-9 / 2  # radius in meters
    D = k_B * temp / (6 * np.pi * viscosity * r)
    return D  # [m²/s]

def hindered_diffusion_factor(particle_radius_nm, pore_radius_nm):
    """
    Hindered diffusion through a cylindrical pore (Renkin equation).
    Returns ratio D_pore / D_bulk.
    Valid for particle_radius / pore_radius < 0.4 (good accuracy)
    """
    lam = particle_radius_nm / pore_radius_nm
    if lam >= 1.0:
        return 0.0  # particle too large for pore
    # Renkin equation (1954)
    H = (1 - lam)**2 * (1 - 2.104*lam + 2.089*lam**3 - 0.948*lam**5)
    return max(H, 0.0)

def flux_through_rwm(sono_params, lnp_params):
    """
    Calculate nanoparticle flux through sonoporated RWM.
    Returns number of LNPs entering perilymph per second.
    """
    # Bulk diffusion coefficient of LNP
    D_bulk = stokes_einstein_diffusion(lnp_params.diameter_nm)
    
    # Hindered diffusion through pores
    H = hindered_diffusion_factor(
        lnp_params.diameter_nm / 2, 
        sono_params.pore_radius_nm
    )
    D_pore = D_bulk * H
    
    # Pore density estimate from permeability enhancement
    # Enhanced permeability = baseline + pore contribution
    # Shih et al.: 5.2× enhancement with 5 courses
    # This means pore area ≈ 4.2× baseline membrane area for transport
    # Estimate: ~1000-5000 pores per mm² based on enhancement factor
    
    # From Zhou et al.: single pore radius ~110 nm
    # Effective pore area per pore: π × (110e-9)² = 3.8e-14 m²
    pore_area_single = np.pi * (sono_params.pore_radius_nm * 1e-9)**2
    
    # RWM total area
    rwm_area = sono_params.rwm_area_mm2 * 1e-6  # [m²]
    
    # Estimate pore density from enhancement factor
    # Enhancement = 1 + (n_pores × pore_area × D_pore) / (rwm_area × D_baseline_rwm)
    # D_baseline_rwm ≈ D_bulk × baseline_hindrance (tight junctions, ~0.001)
    baseline_hindrance = 0.001  # RWM is quite impermeable
    enhancement = sono_params.perm_enhancement_5course
    
    # Solve for n_pores × pore_area
    total_pore_area = rwm_area * (enhancement - 1) * baseline_hindrance / H
    n_pores = total_pore_area / pore_area_single
    n_pores_per_mm2 = n_pores / sono_params.rwm_area_mm2
    
    # Concentration gradient across RWM
    # Outside (middle ear): C_0 = LNP concentration
    C_0 = lnp_params.concentration_particles_mL * 1e6  # [particles/m³]
    # Inside (perilymph): initially 0
    C_inside = 0
    
    # Fick's first law: J = D_pore × (C_0 - C_inside) / L × total_pore_area
    L = sono_params.rwm_thickness_um * 1e-6  # membrane thickness [m]
    
    J = D_pore * (C_0 - C_inside) / L * total_pore_area  # [particles/s]
    
    return {
        'D_bulk_m2_s': D_bulk,
        'D_pore_m2_s': D_pore,
        'hindrance_factor': H,
        'pore_area_single_m2': pore_area_single,
        'n_pores_estimated': int(n_pores),
        'n_pores_per_mm2': n_pores_per_mm2,
        'total_pore_area_um2': total_pore_area * 1e12,
        'flux_particles_per_s': J,
        'flux_particles_per_min': J * 60,
    }

# ============================================================
# MODEL 3: LNP UPTAKE AND PROTEIN PRODUCTION
# ============================================================

def lnp_delivery_simulation(sono_params, lnp_params, cochlear_params, 
                             exposure_min=10.0, dt_s=1.0):
    """
    Time-dependent simulation of LNP delivery through sonoporated RWM.
    
    Phases:
    1. Sonoporation exposure (pores open, LNPs diffuse through)
    2. Post-exposure (pores gradually close, diffusion slows)
    3. LNP distribution in perilymph and uptake by OHCs
    """
    # Get flux parameters
    flux_data = flux_through_rwm(sono_params, lnp_params)
    
    # Time array
    total_time_s = exposure_min * 60 + sono_params.rwm_pore_lifetime_min * 60
    n_steps = int(total_time_s / dt_s)
    t = np.arange(n_steps) * dt_s
    
    # Track LNPs in perilymph
    LNP_in_perilymph = np.zeros(n_steps)
    LNP_uptaken_by_cells = np.zeros(n_steps)
    
    # Pore permeability over time (exponential decay after exposure ends)
    exposure_end_s = exposure_min * 60
    tau_pore_close = sono_params.rwm_pore_lifetime_min * 60 / 3  # time constant
    
    # LNP uptake rate by cells (endocytosis, ~0.1-1% per minute)
    k_uptake = 0.005 / 60  # 0.5% per minute [1/s]
    
    # Perilymph volume
    V_perilymph = cochlear_params.scala_tympani_volume_uL * 1e-9  # [L]
    
    initial_flux = flux_data['flux_particles_per_s']
    
    for i in range(1, n_steps):
        # Pore permeability factor (1 during exposure, decays after)
        if t[i] <= exposure_end_s:
            perm_factor = 1.0
        else:
            perm_factor = np.exp(-(t[i] - exposure_end_s) / tau_pore_close)
        
        # Flux into perilymph (decreasing as concentration gradient drops)
        C_outside = lnp_params.concentration_particles_mL * 1e6  # constant reservoir
        C_inside = LNP_in_perilymph[i-1] / (V_perilymph * 1e3)  # particles/m³
        gradient_factor = max(0, 1 - C_inside / C_outside)
        
        flux_in = initial_flux * perm_factor * gradient_factor * dt_s
        
        # Cellular uptake
        uptake = k_uptake * LNP_in_perilymph[i-1] * dt_s
        
        LNP_in_perilymph[i] = LNP_in_perilymph[i-1] + flux_in - uptake
        LNP_uptaken_by_cells[i] = LNP_uptaken_by_cells[i-1] + uptake
    
    # Final accounting
    total_LNPs_delivered = LNP_in_perilymph[-1] + LNP_uptaken_by_cells[-1]
    total_LNPs_uptaken = LNP_uptaken_by_cells[-1]
    
    # Distribute uptaken LNPs among hair cells
    # OHCs are primary target (most exposed to perilymph in scala tympani)
    # ~60% of uptake goes to OHCs (they outnumber IHCs ~3.5:1 and are closer)
    ohc_fraction = 0.60
    lnps_per_ohc = (total_LNPs_uptaken * ohc_fraction) / cochlear_params.n_OHC_total
    
    # Endosomal escape
    escaped_lnps_per_ohc = lnps_per_ohc * lnp_params.endosomal_escape_pct / 100
    
    # mRNA copies reaching cytoplasm
    mrna_per_ohc = escaped_lnps_per_ohc * lnp_params.mRNA_copies_per_LNP
    
    # Protein production (each mRNA → ~100 proteins before degradation)
    proteins_per_ohc = mrna_per_ohc * lnp_params.mRNA_translation_proteins
    
    # Fraction of therapeutic requirement met
    therapeutic_fraction = proteins_per_ohc / cochlear_params.stereocilin_molecules_per_OHC
    
    return {
        'total_LNPs_entered_perilymph': int(total_LNPs_delivered),
        'total_LNPs_uptaken_by_cells': int(total_LNPs_uptaken),
        'LNPs_per_OHC': round(lnps_per_ohc, 1),
        'escaped_LNPs_per_OHC': round(escaped_lnps_per_ohc, 2),
        'mRNA_copies_per_OHC': round(mrna_per_ohc, 1),
        'proteins_per_OHC_per_dose': int(proteins_per_ohc),
        'therapeutic_requirement': cochlear_params.stereocilin_molecules_per_OHC,
        'fraction_of_therapeutic_dose': round(therapeutic_fraction, 4),
        'percent_of_therapeutic_dose': round(therapeutic_fraction * 100, 2),
        'doses_needed_for_full': max(1, int(np.ceil(1.0 / therapeutic_fraction))) if therapeutic_fraction > 0 else float('inf'),
        'time_profile': {
            'total_time_min': round(total_time_s / 60, 1),
            'exposure_time_min': exposure_min,
            'pore_open_time_min': sono_params.rwm_pore_lifetime_min,
        },
        'flux_data': flux_data,
    }


def compare_delivery_methods(sono_params, lnp_params, cochlear_params):
    """Compare sonoporation+LNP vs standard AAV injection"""
    
    # Sonoporation + LNP results
    sono_result = lnp_delivery_simulation(sono_params, lnp_params, cochlear_params)
    
    # AAV intracochlear injection (standard of care)
    aav_result = {
        'method': 'AAV-Anc80L65 intracochlear injection',
        'invasiveness': 'Surgical (general anesthesia, round window injection)',
        'OHC_transduction_pct': cochlear_params.aav_ohc_transduction_pct,
        'IHC_transduction_pct': cochlear_params.aav_ihc_transduction_pct,
        'OHCs_transduced': int(cochlear_params.n_OHC_total * cochlear_params.aav_ohc_transduction_pct / 100),
        'repeatable': 'No (anti-AAV antibodies after first dose)',
        'expression_duration': 'Years (episomal DNA in post-mitotic cells)',
        'payload_limit_kb': 4.7,
        'time_to_expression_days': '7-14',
        'risks': 'Surgical damage, hearing loss risk, immune response',
        'source': 'Landegger et al. 2017 (PMC5340646), Iversen et al. 2022 (Nat Commun)',
    }
    
    # Sonoporation + LNP summary
    sono_summary = {
        'method': 'Sonoporation + LNP (non-invasive)',
        'invasiveness': 'Non-invasive (gel + ultrasound probe, 10 min, no anesthesia)',
        'proteins_per_OHC_per_dose': sono_result['proteins_per_OHC_per_dose'],
        'percent_therapeutic_per_dose': sono_result['percent_of_therapeutic_dose'],
        'doses_for_full_therapy': sono_result['doses_needed_for_full'],
        'repeatable': 'Yes (no immune memory against LNPs)',
        'expression_duration': 'Transient (mRNA: days to weeks per dose)',
        'payload_limit_kb': 'Unlimited (mRNA can be any length)',
        'time_to_expression_hours': '4-12 (mRNA translation is fast)',
        'risks': 'Minimal (RWM fully recovers in 24h, no ABR shift)',
        'source': 'Shih et al. 2019 (PMC7823126), Zhou et al. 2009 (PMC2752487)',
    }
    
    return {
        'sonoporation_LNP': sono_summary,
        'aav_injection': aav_result,
        'sonoporation_detailed': sono_result,
        'key_insight': (
            f"Single sonoporation+LNP session delivers ~{sono_result['percent_of_therapeutic_dose']}% "
            f"of therapeutic dose per OHC. {sono_result['doses_needed_for_full']} sessions needed. "
            f"But: repeatable (monthly), non-invasive, zero surgical risk. "
            f"AAV: one surgery, {cochlear_params.aav_ohc_transduction_pct}% OHC transduction, "
            f"but irreversible, non-repeatable, surgical risk."
        ),
    }


# ============================================================
# SENSITIVITY ANALYSIS
# ============================================================

def sensitivity_analysis():
    """Vary key parameters to understand model robustness"""
    base_sono = SonoporationParams()
    base_lnp = LNPParams()
    base_cochlear = CochlearParams()
    
    results = {}
    
    # Vary pore radius (most uncertain parameter)
    pore_radii = [70, 90, 110, 130, 150, 200]
    pore_results = []
    for r in pore_radii:
        params = SonoporationParams(pore_radius_nm=r)
        res = lnp_delivery_simulation(params, base_lnp, base_cochlear)
        pore_results.append({
            'pore_radius_nm': r,
            'percent_therapeutic': res['percent_of_therapeutic_dose'],
            'doses_needed': res['doses_needed_for_full'],
            'hindrance': res['flux_data']['hindrance_factor'],
        })
    results['vary_pore_radius'] = pore_results
    
    # Vary LNP concentration
    concentrations = [1e11, 5e11, 1e12, 5e12, 1e13]
    conc_results = []
    for c in concentrations:
        params = LNPParams(concentration_particles_mL=c)
        res = lnp_delivery_simulation(base_sono, params, base_cochlear)
        conc_results.append({
            'concentration_per_mL': f"{c:.0e}",
            'percent_therapeutic': res['percent_of_therapeutic_dose'],
            'doses_needed': res['doses_needed_for_full'],
        })
    results['vary_LNP_concentration'] = conc_results
    
    # Vary exposure time
    times = [5, 10, 15, 20, 30]
    time_results = []
    for t in times:
        res = lnp_delivery_simulation(base_sono, base_lnp, base_cochlear, exposure_min=t)
        time_results.append({
            'exposure_min': t,
            'percent_therapeutic': res['percent_of_therapeutic_dose'],
            'doses_needed': res['doses_needed_for_full'],
        })
    results['vary_exposure_time'] = time_results
    
    # Vary endosomal escape (biggest bottleneck)
    escapes = [1.0, 2.0, 5.0, 10.0, 20.0]
    escape_results = []
    for e in escapes:
        params = LNPParams(endosomal_escape_pct=e)
        res = lnp_delivery_simulation(base_sono, params, base_cochlear)
        escape_results.append({
            'endosomal_escape_pct': e,
            'percent_therapeutic': res['percent_of_therapeutic_dose'],
            'doses_needed': res['doses_needed_for_full'],
        })
    results['vary_endosomal_escape'] = escape_results
    
    return results


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("SONOPORATION-MEDIATED DELIVERY MODEL FOR STRC GENE THERAPY")
    print("=" * 70)
    
    sono = SonoporationParams()
    lnp = LNPParams()
    cochlear = CochlearParams()
    
    # Run comparison
    comparison = compare_delivery_methods(sono, lnp, cochlear)
    
    print("\n--- DELIVERY METHOD COMPARISON ---\n")
    
    print("SONOPORATION + LNP (non-invasive):")
    for k, v in comparison['sonoporation_LNP'].items():
        print(f"  {k}: {v}")
    
    print(f"\nAAV INTRACOCHLEAR INJECTION (surgical):")
    for k, v in comparison['aav_injection'].items():
        print(f"  {k}: {v}")
    
    print(f"\n--- KEY INSIGHT ---")
    print(comparison['key_insight'])
    
    # Detailed sonoporation results
    detail = comparison['sonoporation_detailed']
    print(f"\n--- DETAILED SONOPORATION MODEL ---")
    print(f"  Flux data:")
    for k, v in detail['flux_data'].items():
        if isinstance(v, float):
            print(f"    {k}: {v:.4e}")
        else:
            print(f"    {k}: {v}")
    print(f"\n  Delivery chain:")
    print(f"    Total LNPs entered perilymph: {detail['total_LNPs_entered_perilymph']:,}")
    print(f"    Total LNPs uptaken by cells:  {detail['total_LNPs_uptaken_by_cells']:,}")
    print(f"    LNPs per OHC:                 {detail['LNPs_per_OHC']}")
    print(f"    Escaped endosome per OHC:     {detail['escaped_LNPs_per_OHC']}")
    print(f"    mRNA copies per OHC:          {detail['mRNA_copies_per_OHC']}")
    print(f"    Proteins per OHC (per dose):  {detail['proteins_per_OHC_per_dose']:,}")
    print(f"    Therapeutic requirement:       {detail['therapeutic_requirement']:,}")
    print(f"    % of therapeutic dose:         {detail['percent_of_therapeutic_dose']}%")
    print(f"    Doses needed for full:         {detail['doses_needed_for_full']}")
    
    # Sensitivity analysis
    print(f"\n--- SENSITIVITY ANALYSIS ---")
    sensitivity = sensitivity_analysis()
    
    for param_name, param_results in sensitivity.items():
        print(f"\n  {param_name}:")
        for r in param_results:
            parts = [f"{k}={v}" for k, v in r.items()]
            print(f"    {', '.join(parts)}")
    
    # Save results
    all_results = {
        'comparison': {
            'sonoporation_LNP': comparison['sonoporation_LNP'],
            'aav_injection': comparison['aav_injection'],
            'key_insight': comparison['key_insight'],
        },
        'detailed': {
            k: v for k, v in detail.items() 
            if k != 'flux_data'
        },
        'flux': detail['flux_data'],
        'sensitivity': sensitivity,
        'parameters': {
            'sonoporation': asdict(sono),
            'lnp': asdict(lnp),
            'cochlear': asdict(cochlear),
        },
        'sources': [
            'Zhou et al. 2009. Pore size ~110 nm. PMC2752487',
            'Shih et al. 2019. RWM sonoporation in guinea pig. PMC7823126',
            'Landegger et al. 2017. AAV-Anc80L65. PMC5340646',
            'Iversen et al. 2022. NHP cochlear AAV. Nat Commun',
        ],
    }
    
    with open('/Users/egorlyfar/DeepResearch/strc/sonoporation_results.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\n✅ Results saved to sonoporation_results.json")
    
    # ============================================================
    # OPTIMIZED SCENARIO
    # ============================================================
    print(f"\n{'=' * 70}")
    print("OPTIMIZED SCENARIO: High-concentration LNP + enhanced escape + longer exposure")
    print("=" * 70)
    
    # Optimized LNP (high concentration, ionizable lipid for better escape)
    opt_lnp = LNPParams(
        concentration_particles_mL=1e13,  # 10× higher (achievable with concentration)
        endosomal_escape_pct=15.0,        # Ionizable lipids (SM-102, ALC-0315) achieve 10-20%
        mRNA_copies_per_LNP=6,            # Optimized loading
    )
    
    opt_result = lnp_delivery_simulation(sono, opt_lnp, cochlear, exposure_min=20)
    
    print(f"\n  Optimized parameters:")
    print(f"    LNP concentration:       1×10¹³/mL (concentrated)")
    print(f"    Endosomal escape:         15% (ionizable lipid, SM-102 class)")
    print(f"    mRNA per LNP:             6 (optimized loading)")
    print(f"    Exposure time:            20 min")
    print(f"\n  Results:")
    print(f"    LNPs per OHC:             {opt_result['LNPs_per_OHC']}")
    print(f"    Escaped per OHC:          {opt_result['escaped_LNPs_per_OHC']}")
    print(f"    mRNA per OHC:             {opt_result['mRNA_copies_per_OHC']}")
    print(f"    Proteins per OHC/dose:    {opt_result['proteins_per_OHC_per_dose']:,}")
    print(f"    % therapeutic per dose:   {opt_result['percent_of_therapeutic_dose']}%")
    print(f"    Doses needed:             {opt_result['doses_needed_for_full']}")
    print(f"    → Monthly sessions:       {opt_result['doses_needed_for_full']} months")
    
    # Hybrid scenario: AAV (one-time) + LNP maintenance
    print(f"\n{'=' * 70}")
    print("HYBRID SCENARIO: AAV initial + LNP maintenance top-up")
    print("=" * 70)
    
    aav_initial_pct = cochlear.aav_ohc_transduction_pct
    # If AAV gives 80% transduction but expression declines over years:
    # Year 1-3: full expression
    # Year 3-5: ~60% (some episome loss)  
    # Year 5+: ~40%
    # LNP top-up could maintain the remaining deficit
    
    deficit_year5_pct = 40  # remaining deficit
    opt_lnp_topup = opt_result['percent_of_therapeutic_dose']
    sessions_to_cover_deficit = max(1, int(np.ceil(deficit_year5_pct / opt_lnp_topup)))
    
    print(f"\n  Year 1-3: AAV alone provides {aav_initial_pct}% OHC transduction")
    print(f"  Year 5+:  AAV expression drops to ~{100-deficit_year5_pct}% (episome loss)")
    print(f"  Deficit:  ~{deficit_year5_pct}% needs supplementation")
    print(f"  LNP top-up: {opt_lnp_topup}% per session")
    print(f"  Sessions to cover deficit: {sessions_to_cover_deficit}")
    print(f"  → One session every {max(1, 30 // sessions_to_cover_deficit)} days")
    print(f"\n  This is the most realistic path:")
    print(f"  1. AAV surgery (one-time, under anesthesia, high transduction)")
    print(f"  2. Years later if expression fades: LNP top-ups (non-invasive)")
    print(f"  3. Combines strength of both: AAV durability + LNP repeatability")
    
    # Save optimized results
    all_results['optimized'] = {
        'params': asdict(opt_lnp),
        'results': opt_result,
    }
    all_results['hybrid'] = {
        'aav_initial_pct': aav_initial_pct,
        'deficit_year5_pct': deficit_year5_pct,
        'lnp_per_session_pct': opt_lnp_topup,
        'sessions_to_cover_deficit': sessions_to_cover_deficit,
    }
    
    with open('/Users/egorlyfar/DeepResearch/strc/sonoporation_results.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\n✅ All results saved to sonoporation_results.json")
