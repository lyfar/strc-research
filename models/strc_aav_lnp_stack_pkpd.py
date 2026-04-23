#!/usr/bin/env python3
"""
Integrated AAV Mini-STRC + Strategy B mRNA-LNP stack model.

Answers the question left open in [[STRC mRNA-LNP PKPD v2 OHC Tropism]]:
does Strategy B mRNA-LNP as a maintenance layer on AAV-untransduced OHC
meaningfully lift cochlea-mean fold above threshold, or is AAV doing
all the real work and Strategy B dilutes out?

Model partitions the OHC population into four compartments:

  1. AAV-transduced (fraction: f_aav)          → fold f_aav_fold (from Iranfar 2026 / Holt 2021)
  2. LNP-tropic AAV-untransduced               → fold f_lnp_fold (from v2 PKPD per-OHC result)
  3. LNP-untropic AAV-untransduced             → fold 1.0 (baseline)
  4. (Dead / no rescue possible)               → fold 0 but we assume OHC bodies survive per
                                                  GeneReviews STRC (no mass apoptosis in DFNB16)

Cochlea-mean fold = weighted sum. Therapeutic threshold: ≥ 2.0 cochlea-mean.

Sweep:
  - AAV transduction efficiency: 0, 10, 20, 30, 40, 50, 60, 70, 80%
  - AAV per-OHC fold (when transduced): 3, 5, 7, 10 (Iranfar 2026 range)
  - LNP tropism: 0.8% (untargeted), 1%, 3% (realistic), 5%, 10%, 20% (aspirational)
  - LNP per-OHC fold (when tropic): 2.85 (v2 best) or 2.5 (conservative)

Output:
  strc_aav_lnp_stack_pkpd.json — full sweep grid + best regimen per AAV level
  + answer to "does Strategy B add meaningful lift?"
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

OUT_JSON = Path("/Users/egorlyfar/Brain/research/strc/models/strc_aav_lnp_stack_pkpd.json")

# Parameter ranges.
AAV_TRANSDUCTION = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
AAV_FOLD_WHEN_TRANSDUCED = [3.0, 5.0, 7.0, 10.0]
LNP_TROPISM = [0.008, 0.01, 0.03, 0.05, 0.10, 0.20]
LNP_FOLD_WHEN_TROPIC = [2.85, 2.5]

THRESHOLD = 2.0  # cochlea-mean fold needed for therapeutic effect
BASELINE_FOLD = 1.0  # DFNB16 OHC without rescue (structural connectors lost but cells survive)


def stack_cochlea_mean(f_aav: float, aav_fold: float,
                       lnp_trop: float, lnp_fold: float) -> dict:
    """Partition OHC population and compute weighted cochlea-mean fold.

    AAV transduction is independent of LNP tropism (different delivery
    routes, different cellular uptake mechanisms — can be treated as
    Bernoulli-independent for a first-order model).
    """
    # Four compartments:
    # AAV+ (size f_aav, any LNP status — AAV rescue dominates)
    # AAV- LNP+ (size (1-f_aav) * lnp_trop)
    # AAV- LNP- (size (1-f_aav) * (1-lnp_trop))

    frac_aav_pos = f_aav
    frac_aav_neg_lnp_pos = (1 - f_aav) * lnp_trop
    frac_aav_neg_lnp_neg = (1 - f_aav) * (1 - lnp_trop)

    # Fold per compartment. Assume AAV rescue dominates when both (AAV+LNP co-treats additive?
    # For simplicity and conservatively: in AAV+ cells, LNP doesn't help further
    # since functional STRC is already present; this is CI (competitive inhibition
    # on UPR load).
    cochlea_mean = (
        frac_aav_pos * aav_fold
        + frac_aav_neg_lnp_pos * lnp_fold
        + frac_aav_neg_lnp_neg * BASELINE_FOLD
    )

    # Strategy B "added value" = stack cochlea-mean MINUS AAV-alone cochlea-mean.
    aav_alone = frac_aav_pos * aav_fold + (1 - f_aav) * BASELINE_FOLD
    added_value = cochlea_mean - aav_alone
    relative_added_pct = (added_value / aav_alone * 100) if aav_alone > 0 else 0.0

    return {
        "frac_aav_pos": round(frac_aav_pos, 3),
        "frac_aav_neg_lnp_pos": round(frac_aav_neg_lnp_pos, 4),
        "frac_aav_neg_lnp_neg": round(frac_aav_neg_lnp_neg, 3),
        "cochlea_mean_fold": round(cochlea_mean, 3),
        "aav_alone_cochlea_mean_fold": round(aav_alone, 3),
        "strategy_b_added_fold": round(added_value, 3),
        "strategy_b_added_pct_relative": round(relative_added_pct, 2),
        "passes_threshold_2x": cochlea_mean >= THRESHOLD,
        "aav_alone_passes_2x": aav_alone >= THRESHOLD,
    }


def main():
    sweep = []
    for f_aav, aav_fold, lnp_trop, lnp_fold in itertools.product(
        AAV_TRANSDUCTION, AAV_FOLD_WHEN_TRANSDUCED, LNP_TROPISM, LNP_FOLD_WHEN_TROPIC
    ):
        row = {
            "aav_transduction": f_aav,
            "aav_fold_transduced": aav_fold,
            "lnp_tropism": lnp_trop,
            "lnp_fold_tropic": lnp_fold,
            **stack_cochlea_mean(f_aav, aav_fold, lnp_trop, lnp_fold),
        }
        sweep.append(row)

    # Summaries of interest:
    # 1. Does Strategy B ever change pass/fail verdict (AAV alone fails but stack passes)?
    rescue_pivots = [r for r in sweep if r["passes_threshold_2x"] and not r["aav_alone_passes_2x"]]

    # 2. What's the max relative added value from Strategy B across sweep?
    max_added_pct = max(sweep, key=lambda r: r["strategy_b_added_pct_relative"])

    # 3. At realistic AAV (50-70% transduction at 5× fold), what's Strategy B contribution?
    realistic = [r for r in sweep if r["aav_transduction"] in (0.5, 0.6, 0.7)
                 and r["aav_fold_transduced"] == 5.0
                 and r["lnp_tropism"] in (0.008, 0.01, 0.03)
                 and r["lnp_fold_tropic"] == 2.85]

    # 4. Worst-case AAV (20% transduction at 3×): does Strategy B save?
    worst_aav = [r for r in sweep if r["aav_transduction"] == 0.2
                 and r["aav_fold_transduced"] == 3.0
                 and r["lnp_fold_tropic"] == 2.85]

    payload = {
        "model": "Integrated AAV Mini-STRC + Strategy B mRNA-LNP stack (cochlea-mean fold)",
        "date": "2026-04-23",
        "threshold_fold": THRESHOLD,
        "baseline_fold_no_rescue": BASELINE_FOLD,
        "n_sweep_rows": len(sweep),
        "ranges": {
            "aav_transduction": AAV_TRANSDUCTION,
            "aav_fold_when_transduced": AAV_FOLD_WHEN_TRANSDUCED,
            "lnp_tropism": LNP_TROPISM,
            "lnp_fold_when_tropic": LNP_FOLD_WHEN_TROPIC,
        },
        "rescue_pivots_count": len(rescue_pivots),
        "rescue_pivots": rescue_pivots[:20],
        "max_strategy_b_added_pct": max_added_pct,
        "realistic_aav_scenarios": realistic,
        "worst_aav_scenarios": worst_aav,
        "full_sweep": sweep,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {OUT_JSON}")

    print(f"\n=== STRC AAV Mini-STRC + Strategy B mRNA-LNP Stack Model ===")
    print(f"Threshold: {THRESHOLD}x cochlea-mean fold")
    print(f"Sweep rows: {len(sweep)}")

    print(f"\n--- Does Strategy B EVER change pass/fail verdict (AAV alone fails, stack passes)?")
    print(f"  Rescue pivots: {len(rescue_pivots)} / {len(sweep)}")
    if rescue_pivots:
        for r in rescue_pivots[:5]:
            print(f"    AAV {r['aav_transduction']*100:.0f}% @{r['aav_fold_transduced']:.1f}x + LNP {r['lnp_tropism']*100:.1f}% @{r['lnp_fold_tropic']:.1f}x : stack={r['cochlea_mean_fold']:.3f}, alone={r['aav_alone_cochlea_mean_fold']:.3f}")

    print(f"\n--- Max Strategy B relative added value across entire sweep:")
    m = max_added_pct
    print(f"  {m['strategy_b_added_pct_relative']:.2f}% added, at AAV {m['aav_transduction']*100:.0f}% @{m['aav_fold_transduced']:.1f}x, LNP {m['lnp_tropism']*100:.1f}% @{m['lnp_fold_tropic']:.1f}x")

    print(f"\n--- Realistic AAV (50-70%) + realistic LNP (0.8-3%):")
    for r in realistic:
        print(f"  AAV {r['aav_transduction']*100:.0f}% @5x + LNP {r['lnp_tropism']*100:.1f}% @2.85x: "
              f"stack {r['cochlea_mean_fold']:.3f} (AAV alone {r['aav_alone_cochlea_mean_fold']:.3f}) "
              f"— Strategy B adds {r['strategy_b_added_pct_relative']:.2f}%")

    print(f"\n--- Worst-case AAV (20% @3x) + various LNP:")
    for r in worst_aav:
        print(f"  LNP {r['lnp_tropism']*100:.1f}% @2.85x: stack {r['cochlea_mean_fold']:.3f} "
              f"(AAV alone {r['aav_alone_cochlea_mean_fold']:.3f}) "
              f"— threshold {r['passes_threshold_2x']} (alone: {r['aav_alone_passes_2x']})")


if __name__ == "__main__":
    main()
