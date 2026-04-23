#!/usr/bin/env python3
"""
Calcium Oscillation Acoustic Therapy — maternal-allele-only pharmacodynamic
model for Misha (#5 A-tier next step).

Question: given Misha's compound-het genotype (paternal 98 kb deletion null +
maternal c.4976A>C E1659A missense), can acoustic-driven 1.3-2.6× fold
induction of the E1659A mutant protein compensate its reduced binding
affinity enough to reach therapeutic occupancy?

Model:

  Total functional STRC occupancy =
      allele_dosage × acoustic_induction × (1 / affinity_penalty)

where:
  allele_dosage = 0.5  (only maternal, diploid scaling)
  acoustic_induction ∈ [1.0, 2.6]  (from Phases 1-3 AC1-CREB model)
  affinity_penalty ∈ [1x, 1000x]  (E1659A reduction; 8.4 kcal/mol
      electrostatic estimate = ~1_000_000x but that's destabilisation,
      not binding; realistic affinity penalty from structural biology
      of ARM-repeat missense variants is typically 3-30x)

Therapeutic gates:
  - Partial rescue: occupancy >= 0.3 (hypothesis note's f=0.3 ABR shift
    22 dB rescue)
  - Full rescue: occupancy >= 1.0 (WT function restored)

Output: sweep JSON + feasibility map + decision:
  at what affinity penalty can acoustic therapy deliver partial / full
  rescue, and is it within the biologically plausible range for ARM-repeat
  missense mutations.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

OUT_JSON = Path("/Users/egorlyfar/Brain/research/strc/models/ca_osc_maternal_allele_only.json")

ALLELE_DOSAGE_MATERNAL_ONLY = 0.5   # Misha expresses 1 of 2 alleles, and only maternal
ACOUSTIC_FOLD_RANGE = [1.0, 1.3, 1.6, 1.9, 2.2, 2.6]   # from AC1-CREB Phases 1-3
AFFINITY_PENALTY_RANGE = [1.0, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 300.0, 1000.0]

PARTIAL_RESCUE_THRESHOLD = 0.30   # f=0.3 → 22 dB ABR rescue per bundle mechanics
FULL_RESCUE_THRESHOLD = 1.0       # restore WT function

# Reference: 8.4 kcal/mol electrostatic penalty → ΔΔG for MM binding affinity
# is modulated by (delta ΔG / RT). For salt bridge loss (ionic): ΔΔG ≈ 2-4 kcal
# contributes to binding, others to destabilisation (folding, packing). Upper
# bound of pure binding penalty at 4 kcal/mol → 600x affinity reduction at 310K.
# Realistic arm-repeat extracellular missense: 3-30x (literature average).
REALISTIC_AFFINITY_PENALTY_BAND = (3.0, 30.0)   # for interpretation
KCAL_PER_MOL_TO_FOLD = lambda dG_kcal: math.exp(dG_kcal * 4184 / (8.314 * 310)) / 1.0


def occupancy(acoustic_fold: float, affinity_penalty: float) -> float:
    return ALLELE_DOSAGE_MATERNAL_ONLY * acoustic_fold * (1.0 / affinity_penalty)


def main():
    sweep = []
    for fold, penalty in itertools.product(ACOUSTIC_FOLD_RANGE, AFFINITY_PENALTY_RANGE):
        occ = occupancy(fold, penalty)
        sweep.append({
            "acoustic_fold": fold,
            "affinity_penalty_x": penalty,
            "effective_occupancy": round(occ, 4),
            "partial_rescue_pass": occ >= PARTIAL_RESCUE_THRESHOLD,
            "full_rescue_pass": occ >= FULL_RESCUE_THRESHOLD,
        })

    # Feasibility boundaries
    partial_boundary = {}
    for fold in ACOUSTIC_FOLD_RANGE:
        max_penalty_pass = None
        for penalty in sorted(AFFINITY_PENALTY_RANGE):
            if occupancy(fold, penalty) >= PARTIAL_RESCUE_THRESHOLD:
                max_penalty_pass = penalty
            else:
                break
        partial_boundary[fold] = max_penalty_pass

    full_boundary = {}
    for fold in ACOUSTIC_FOLD_RANGE:
        max_penalty_pass = None
        for penalty in sorted(AFFINITY_PENALTY_RANGE):
            if occupancy(fold, penalty) >= FULL_RESCUE_THRESHOLD:
                max_penalty_pass = penalty
            else:
                break
        full_boundary[fold] = max_penalty_pass

    # In realistic penalty band (3-30x), does acoustic therapy pass partial rescue?
    realistic_pass_partial = sum(
        1 for r in sweep
        if REALISTIC_AFFINITY_PENALTY_BAND[0] <= r["affinity_penalty_x"] <= REALISTIC_AFFINITY_PENALTY_BAND[1]
        and r["partial_rescue_pass"]
    )
    realistic_total = sum(
        1 for r in sweep
        if REALISTIC_AFFINITY_PENALTY_BAND[0] <= r["affinity_penalty_x"] <= REALISTIC_AFFINITY_PENALTY_BAND[1]
    )

    payload = {
        "model": "Ca-oscillation acoustic therapy — maternal-allele-only PD for Misha",
        "date": "2026-04-23",
        "assumptions": {
            "allele_dosage_maternal_only": ALLELE_DOSAGE_MATERNAL_ONLY,
            "partial_rescue_threshold_f": PARTIAL_RESCUE_THRESHOLD,
            "full_rescue_threshold_f": FULL_RESCUE_THRESHOLD,
            "realistic_affinity_penalty_band_x": REALISTIC_AFFINITY_PENALTY_BAND,
        },
        "ranges": {
            "acoustic_fold": ACOUSTIC_FOLD_RANGE,
            "affinity_penalty_x": AFFINITY_PENALTY_RANGE,
        },
        "sweep": sweep,
        "partial_rescue_feasibility_boundary_max_penalty": partial_boundary,
        "full_rescue_feasibility_boundary_max_penalty": full_boundary,
        "realistic_penalty_band_partial_pass_rate": f"{realistic_pass_partial}/{realistic_total}",
        "kcal_reference": {
            "kcal_4_to_fold_penalty": round(KCAL_PER_MOL_TO_FOLD(4.0), 1),
            "kcal_2_to_fold_penalty": round(KCAL_PER_MOL_TO_FOLD(2.0), 1),
            "kcal_1_to_fold_penalty": round(KCAL_PER_MOL_TO_FOLD(1.0), 1),
            "note": "in-silico 8.4 kcal is destabilisation energy (folding), not pure binding; binding penalty alone is typically 3-30x for ARM-repeat missense",
        },
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {OUT_JSON}")

    # Pretty summary
    print(f"\n=== Calcium Oscillation — Maternal-Allele-Only PD Sweep ===")
    print(f"Misha: paternal 98 kb del (null) + maternal E1659A. Allele dosage 0.5.")
    print(f"Partial rescue threshold: f = {PARTIAL_RESCUE_THRESHOLD} (≈22 dB ABR)")
    print(f"Full rescue threshold: f = {FULL_RESCUE_THRESHOLD} (WT restore)")
    print()
    print(f"{'':>15}", end="")
    for p in AFFINITY_PENALTY_RANGE:
        print(f"{p:>8.0f}x", end="")
    print()
    print("-" * (15 + 9 * len(AFFINITY_PENALTY_RANGE)))
    for fold in ACOUSTIC_FOLD_RANGE:
        print(f"fold {fold:>4.1f}x     ", end="")
        for penalty in AFFINITY_PENALTY_RANGE:
            occ = occupancy(fold, penalty)
            if occ >= FULL_RESCUE_THRESHOLD:
                mark = "FULL"
            elif occ >= PARTIAL_RESCUE_THRESHOLD:
                mark = "PART"
            else:
                mark = "  - "
            print(f" {occ:>5.2f}{mark[0]:>1}", end="" if len(mark)==4 else "")
            print(f"{mark[1:]:>3}", end="") if len(mark)==4 else None
        print()

    # Text rows for partial boundary
    print(f"\nPartial-rescue feasibility boundary (max penalty that still passes):")
    for fold, max_p in partial_boundary.items():
        print(f"  acoustic {fold:.1f}x → max penalty {max_p}x" if max_p else f"  acoustic {fold:.1f}x → cannot reach partial rescue at any tested penalty")

    print(f"\nFull-rescue feasibility boundary:")
    for fold, max_p in full_boundary.items():
        print(f"  acoustic {fold:.1f}x → max penalty {max_p}x" if max_p else f"  acoustic {fold:.1f}x → cannot reach full rescue (even at 1x penalty, max occupancy = {ALLELE_DOSAGE_MATERNAL_ONLY*fold:.2f})")

    print(f"\nIn realistic affinity-penalty band (3-30x): {realistic_pass_partial}/{realistic_total} fold×penalty combos pass partial rescue.")


if __name__ == "__main__":
    main()
