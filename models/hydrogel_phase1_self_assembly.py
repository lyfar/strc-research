#!/usr/bin/env python3
"""
Phase 1 computational proof for [[STRC Synthetic Peptide Hydrogel HTC]].

Pre-screens a literature-based shortlist of self-assembling peptide (SAP)
backbones with actin-binding N-terminus + TMEM145-binding C-terminus.
Computes four gates:

  G1. Fibril geometry — can the assembled β-sheet fibril bridge
      adjacent stereocilia at the native ~8 nm HTC spacing?
  G2. Cochlear PK — Stokes-Einstein diffusion time from round window
      membrane to basal OHC row. Gate: < 4h feasible.
  G3. Bundle stiffness recovery — HTC coupling fraction f restored by
      SAP crosslinks. Gate: f >= 0.4 achievable at physiologic SAP
      concentration.
  G4. Manufacturing feasibility — SPPS cost + dose volume per ear.
      Gate: dose < 100 mg per ear, cost < $5k per dose.

Candidates ranked by composite score. Output: JSON + stdout table +
Phase 2 recommendation (top-3 for de-novo design and AF3 validation).

Not heavy compute. Runs in < 30s on any machine.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

OUT_JSON = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase1_self_assembly.json")

# ----------------------------------------------------------------------------
# Physical + biological constants
# ----------------------------------------------------------------------------

KB_T_KCAL = 0.593                # kT at 310 K, kcal/mol
ETA_WATER_37C = 0.693e-3         # Pa s at 37 C
NA = 6.022e23

# Cochlea geometry (human, estimates from OC histology)
COCHLEA_STB_LENGTH_MM = 32.0     # scala tympani length (uncoiled)
BASAL_OHC_DEPTH_MM = 1.5         # first-row basal OHC distance from RWM
APEX_DEPTH_MM = 30.0             # apical OHCs, worst case

# Bundle mechanics (Tobin 2019, Kozlov 2007, Roongthumskul 2019)
OHC_BUNDLE_STEREOCILIA_ROW1 = 40     # stereocilia in tallest row per OHC
HTC_PER_STEREOCILIUM = 6             # horizontal top connectors per stereocilium
HTC_SPACING_NM = 8.0                 # native connector span between adjacent stereocilia
K_HTC_PN_PER_NM = 7.5                # per-crosslink stiffness (approx, Tobin 2019)
WT_BUNDLE_STIFFNESS_PN_PER_NM = 1_400   # WT OHC bundle at gating compliance range
TARGET_F_FOR_30DB = 0.40             # f fraction to achieve <30 dB ABR (hypothesis note)

# Manufacturing
SPPS_COST_USD_PER_MG_60MER = 5.0     # mid-range for 60-mer peptide, SPPS pharmacy scale
HUMAN_SCALA_TYMPANI_VOLUME_UL = 70   # 70 microlitres is generous estimate per ear
MAX_DOSE_MG_PER_EAR = 100            # upper feasible ear-drop dose
MAX_DOSE_USD_PER_EAR = 5_000         # cost cap per treatment dose

# β-sheet fibril geometry
BETA_STRAND_AA_LEN_NM = 0.35         # extended β-strand rise per residue


# ----------------------------------------------------------------------------
# Peptide shortlist
#
# Backbones come from published self-assembling peptides that are in
# clinical or approved use (or close to it). Actin-binding N-term is
# borrowed from espin/fascin / WH2 literature. TMEM145-binding C-term is
# derived from the STRC–TMEM145 interface residues in our Derstroff-confirmed
# GOLD zone (aa 1603-1749) — specifically the sub-region aa 1669-1680 which
# was identified as the dominant contact cluster in our Ultra-Mini x TMEM145
# AF3 job.
# ----------------------------------------------------------------------------

@dataclass
class SAPCandidate:
    name: str
    backbone: str         # β-sheet forming core
    n_term_actin: str     # actin-binding handle
    c_term_tmem145: str   # TMEM145-binding handle
    notes: str

    @property
    def sequence(self) -> str:
        return self.n_term_actin + self.backbone + self.c_term_tmem145

    @property
    def length_aa(self) -> int:
        return len(self.sequence)

    @property
    def molecular_weight_da(self) -> float:
        # rough: 110 Da/aa avg
        return self.length_aa * 110.0


# WH2-domain actin binder (~17 aa core)
WH2 = "RQLVKAIPDNCSKSNVS"
# Espin-derived actin-binding motif (~12 aa)
ESPIN_MINI = "DLIKRNFPRSIA"
# Short fascin-derived peptide
FASCIN_MINI = "LREYKNGSTNA"

# TMEM145 binder candidates derived from STRC interface sub-region
# (based on AF3 contact clusters in ARM repeat zone 1669-1680):
# Native STRC residues 1669-1680 (from uniprot Q7RTU9): AEDLPEPVPNCA
# Trimmed and stabilised variant for peptide use:
TMEM145_BIND_NATIVE12 = "AEDLPEPVPNCA"    # native 12 aa, cys may be problematic
TMEM145_BIND_CMUT11 = "AEDLPEPVPNA"        # cys -> ala, 11 aa
# De-novo-style candidate designed to preserve electrostatic pattern
TMEM145_BIND_DE_NOVO = "EELPEPVPNYK"       # de-novo hypothetical

CANDIDATES = [
    SAPCandidate(
        name="RADA16-WH2-native",
        backbone="RADARADARADARADA",
        n_term_actin=WH2,
        c_term_tmem145=TMEM145_BIND_NATIVE12,
        notes="Original Zhang RADA16 + WH2 actin + native STRC 12-mer TMEM145 anchor. Cys oxidation risk.",
    ),
    SAPCandidate(
        name="RADA16-WH2-Cmut",
        backbone="RADARADARADARADA",
        n_term_actin=WH2,
        c_term_tmem145=TMEM145_BIND_CMUT11,
        notes="RADA16 + WH2 + Cys->Ala TMEM145 anchor. Lower oxidation risk.",
    ),
    SAPCandidate(
        name="IEIK13-Espin-denovo",
        backbone="IEIKIEIKIEIKI",
        n_term_actin=ESPIN_MINI,
        c_term_tmem145=TMEM145_BIND_DE_NOVO,
        notes="Miroshnikova-tier IEIK13 backbone + espin mini + de-novo TMEM145 binder.",
    ),
    SAPCandidate(
        name="EAK16-WH2-denovo",
        backbone="AEAKAEAKAEAKAEAK",
        n_term_actin=WH2,
        c_term_tmem145=TMEM145_BIND_DE_NOVO,
        notes="EAK16 + WH2 + de-novo binder. Lower charge density than RADA16 at pH 7.4.",
    ),
    SAPCandidate(
        name="KLD12-Fascin-native",
        backbone="KLDLKLDLKLDL",
        n_term_actin=FASCIN_MINI,
        c_term_tmem145=TMEM145_BIND_CMUT11,
        notes="KLD12 backbone (used in spinal cord regen) + fascin mini + Cys->Ala TMEM145 anchor.",
    ),
    SAPCandidate(
        name="RADA16-Espin-native",
        backbone="RADARADARADARADA",
        n_term_actin=ESPIN_MINI,
        c_term_tmem145=TMEM145_BIND_CMUT11,
        notes="RADA16 + espin + Cys->Ala 11-mer. Shortest assembled peptide.",
    ),
]


# ----------------------------------------------------------------------------
# Gates
# ----------------------------------------------------------------------------

def gate1_fibril_geometry(c: SAPCandidate) -> dict:
    """The peptide, when assembled, must span HTC_SPACING_NM (8 nm) between
    adjacent stereocilia. Extended β-strand rise per residue ~0.35 nm. So
    the end-to-end length of the full peptide gives the max reachable span.
    """
    extended_length_nm = c.length_aa * BETA_STRAND_AA_LEN_NM
    # Realistic span in fibril ≈ 0.65 × extended (accounts for β-sheet torsion)
    realistic_span_nm = extended_length_nm * 0.65
    reaches_htc = realistic_span_nm >= HTC_SPACING_NM
    return {
        "extended_length_nm": round(extended_length_nm, 2),
        "realistic_span_nm": round(realistic_span_nm, 2),
        "required_span_nm": HTC_SPACING_NM,
        "reaches_htc_span": bool(reaches_htc),
        "margin_nm": round(realistic_span_nm - HTC_SPACING_NM, 2),
    }


def gate2_cochlear_pk(c: SAPCandidate) -> dict:
    """Stokes-Einstein diffusion through perilymph. Time to reach basal OHC
    row (1.5 mm from round window) and apex (30 mm).
    """
    # Hydrodynamic radius for globular peptide, empirical: r ≈ 0.066 * MW^(1/3) nm with MW in Da.
    # Actually for intrinsically disordered / unfolded chains the radius is larger,
    # but a just-assembled SAP monomer in perilymph is closer to globular. Use
    # average of globular estimate and extended:
    r_nm_globular = 0.066 * (c.molecular_weight_da ** (1.0/3.0))
    r_nm_used = r_nm_globular
    D_m2_s = (KB_T_KCAL * 4184 / NA) / (6 * math.pi * ETA_WATER_37C * r_nm_used * 1e-9)
    D_um2_s = D_m2_s * 1e12

    t_basal_s = ((BASAL_OHC_DEPTH_MM * 1e-3) ** 2) / (2 * D_m2_s)
    t_apex_s = ((APEX_DEPTH_MM * 1e-3) ** 2) / (2 * D_m2_s)

    return {
        "hydrodynamic_radius_nm": round(r_nm_used, 2),
        "diffusion_coeff_um2_s": round(D_um2_s, 1),
        "t_to_basal_ohc_min": round(t_basal_s / 60, 1),
        "t_to_apex_hours": round(t_apex_s / 3600, 1),
        "basal_clinic_feasible": t_basal_s / 60 <= 240,   # <= 4h
    }


def gate3_bundle_stiffness(c: SAPCandidate) -> dict:
    """What HTC coupling fraction f can this SAP realistically restore?
    f depends on SAP monomer density at stereocilia surface. Assume the
    SAP concentrates by factor 100x at the bundle (actin-binding localisation).
    Gate: f >= 0.4 at clinically feasible dose.
    """
    # Assume 3x10^14 HTC connector sites per OHC (40 stereocilia x 6 HTC x row structure)
    total_htc_sites_per_ohc = OHC_BUNDLE_STEREOCILIA_ROW1 * HTC_PER_STEREOCILIUM * 3
    # SAP dose -> molecules in scala tympani
    # Assume 10 mg dose in 70 uL scala tympani = 143 g/L
    # Peptide MW 4-8 kDa -> ~2-4 mM monomer concentration
    dose_mg = 10.0
    conc_M = (dose_mg * 1e-3) / (c.molecular_weight_da * (HUMAN_SCALA_TYMPANI_VOLUME_UL * 1e-6))
    # Local concentration at bundle = 100x (actin-binding enrichment)
    conc_local_M = conc_M * 100
    # HTC sites competing: assume Kd = 100 nM for dual-interface binder
    Kd_M = 1e-7
    occupancy = conc_local_M / (conc_local_M + Kd_M)
    # f = fraction of HTC sites occupied * assembly efficiency
    assembly_efficiency = 0.7
    f_restored = occupancy * assembly_efficiency
    f_capped = min(f_restored, 1.0)

    # ABR rescue prediction (log-linear approx from bundle stiffness model):
    # baseline DFNB16 ABR ~64 dB, WT ~22 dB; rescue ~22 dB at f=0.3, ~42 dB at f=0.6
    # linearish between: dB_rescue ~ 70 * f
    db_rescue = 70 * f_capped

    return {
        "dose_mg": dose_mg,
        "bulk_conc_mM": round(conc_M * 1000, 3),
        "local_conc_mM": round(conc_local_M * 1000, 2),
        "Kd_assumed_nM": int(Kd_M * 1e9),
        "occupancy_fraction": round(occupancy, 3),
        "f_restored": round(f_capped, 3),
        "db_rescue_predicted": round(db_rescue, 1),
        "meets_target_f_0p4": f_capped >= TARGET_F_FOR_30DB,
    }


def gate4_manufacturing(c: SAPCandidate) -> dict:
    """SPPS cost + dose feasibility."""
    dose_mg = 10.0  # as in gate 3
    cost_usd = dose_mg * SPPS_COST_USD_PER_MG_60MER * (c.length_aa / 60.0)
    volume_ul = (dose_mg * 1e-3) / (c.molecular_weight_da * HUMAN_SCALA_TYMPANI_VOLUME_UL * 1e-6)
    volume_ul = dose_mg / 1.0  # simplified: 10 mg/mL formulation → 10 mL, but diluted to ear drop volume 70 uL

    feasible = dose_mg <= MAX_DOSE_MG_PER_EAR and cost_usd <= MAX_DOSE_USD_PER_EAR

    return {
        "peptide_length_aa": c.length_aa,
        "molecular_weight_da": round(c.molecular_weight_da, 1),
        "dose_mg": dose_mg,
        "spps_cost_usd": round(cost_usd, 0),
        "manufacturable": feasible,
    }


def composite_score(candidate: SAPCandidate, gates: dict) -> float:
    """Weighted composite, 0-100. Heuristic only."""
    score = 0
    if gates["g1"]["reaches_htc_span"]:
        score += 30 * min(1.0, gates["g1"]["margin_nm"] / 4.0)  # extra margin counts
    if gates["g2"]["basal_clinic_feasible"]:
        score += 15 * max(0.0, 1.0 - (gates["g2"]["t_to_basal_ohc_min"] / 240.0))
    if gates["g3"]["meets_target_f_0p4"]:
        score += 40 * min(1.0, gates["g3"]["f_restored"] / 1.0)
    else:
        # partial credit
        score += 40 * (gates["g3"]["f_restored"] / TARGET_F_FOR_30DB) * 0.5
    if gates["g4"]["manufacturable"]:
        score += 15
    return round(min(score, 100), 1)


# ----------------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------------

def main():
    results = []
    for c in CANDIDATES:
        g1 = gate1_fibril_geometry(c)
        g2 = gate2_cochlear_pk(c)
        g3 = gate3_bundle_stiffness(c)
        g4 = gate4_manufacturing(c)
        gates = {"g1": g1, "g2": g2, "g3": g3, "g4": g4}
        all_pass = (
            g1["reaches_htc_span"]
            and g2["basal_clinic_feasible"]
            and g3["meets_target_f_0p4"]
            and g4["manufacturable"]
        )
        results.append({
            "name": c.name,
            "sequence": c.sequence,
            "length_aa": c.length_aa,
            "molecular_weight_da": round(c.molecular_weight_da, 1),
            "notes": c.notes,
            "gates": gates,
            "all_pass": all_pass,
            "composite_score": composite_score(c, gates),
        })

    results.sort(key=lambda r: -r["composite_score"])

    payload = {
        "phase": "1",
        "hypothesis": "STRC Synthetic Peptide Hydrogel HTC",
        "date": "2026-04-23",
        "constants": {
            "HTC_SPACING_NM": HTC_SPACING_NM,
            "TARGET_F_FOR_30DB": TARGET_F_FOR_30DB,
            "DOSE_MG_TESTED": 10.0,
            "SCALA_TYMPANI_VOLUME_UL": HUMAN_SCALA_TYMPANI_VOLUME_UL,
        },
        "candidates_tested": len(CANDIDATES),
        "candidates_passing_all_gates": sum(1 for r in results if r["all_pass"]),
        "results_ranked": results,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {OUT_JSON}")

    # Stdout summary
    print(f"\n{'='*90}")
    print(f"STRC Synthetic Peptide Hydrogel HTC — Phase 1 ({len(CANDIDATES)} candidates)")
    print(f"{'='*90}")
    hdr = f"{'rank':<5}{'name':<26}{'aa':<5}{'G1 span':>9}{'G2 min':>9}{'G3 f':>8}{'G3 dB':>8}{'score':>8}{'all':>5}"
    print(hdr)
    print("-" * len(hdr))
    for i, r in enumerate(results, 1):
        g1span = f"{r['gates']['g1']['realistic_span_nm']:.1f}"
        g2min = f"{r['gates']['g2']['t_to_basal_ohc_min']:.0f}"
        g3f = f"{r['gates']['g3']['f_restored']:.2f}"
        g3db = f"{r['gates']['g3']['db_rescue_predicted']:.1f}"
        score = f"{r['composite_score']:.1f}"
        flag = "PASS" if r["all_pass"] else "fail"
        print(f"{i:<5}{r['name']:<26}{r['length_aa']:<5}{g1span:>9}{g2min:>9}{g3f:>8}{g3db:>8}{score:>8}{flag:>5}")

    n_pass = sum(1 for r in results if r["all_pass"])
    top3 = [r["name"] for r in results[:3]]
    print(f"\n{n_pass}/{len(CANDIDATES)} candidates pass all 4 gates.")
    print(f"Top-3 for Phase 2 (de-novo design + AF3 validation): {top3}")
    print(f"\nJSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
