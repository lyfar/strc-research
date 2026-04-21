"""
Ca²⁺ Phase 2: parameter sensitivity sweep of the AC1-CREB pivot.

Goal: the Phase 1B pivot gave 1.82× silence→45 dB STRC protein fold and
monotonic Ca→protein up-regulation. Is that robust to parameter choice,
or did we accidentally calibrate into a narrow sweet spot?

Method (one-at-a-time local sensitivity):
  - baseline (Phase 1B parameters)
  - for each of 8 pharmacology parameters:
      run 0.5× and 2× perturbations
      evaluate at SPL = 30 (silence), 45 dB, 75 dB, 100 dB
  - report:
      * fold-change silence → 45 dB (the clinical regime)
      * fold-change silence → 100 dB (saturation)
      * monotonicity with 0.5% tol
      * CREB-P dynamic range
  - flag variants that: (a) lose monotonicity, (b) give <1.2× fold at 45 dB,
    or (c) saturate below 30 dB (silence already maxed).

Expected outcome: a "robustness map" showing which parameters are load-bearing
and which are safe. Results drive both (a) confidence in Touch Grass rationale
and (b) wet-lab priority — parameters with high sensitivity are where we need
literature/experimental pinning.
"""

from __future__ import annotations
import copy
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")

# --------------------------------------------------------------------------- #
# Baseline parameters (from Phase 1B pivot)
# --------------------------------------------------------------------------- #
BASELINE = {
    "CA_REST_NM": 25.0,
    "CA_CEILING_NM": 3000.0,
    "SPL_HALF_SATURATION": 70.0,
    "K_EXTRUSION_S": 30.0,
    "K_CA_AC1_NM": 150.0,
    "N_CA_AC1": 2.0,
    "K_AC1_ACT_S": 0.1,
    "K_AC1_INACT_S": 0.2,
    "AC1_VMAX_NM_S": 2000.0,
    "PDE4_VMAX_NM_S": 8000.0,
    "PDE4_KM_NM": 4000.0,
    "CAMP_BASAL_PROD_NM_S": 10.0,
    "K_CAMP_PKA_NM": 300.0,
    "N_PKA": 2.0,
    "K_PKA_ACT_S": 0.5,
    "K_CREB_PHOS_S": 0.02,
    "K_CREB_DEPHOS_S": 0.005,
    "K_STRC_TXN_BASAL_S": 1e-4,
    "K_CREB_CRE_HALF": 0.2,
    "N_CRE": 1.0,
    "K_TXN_MAX_FOLD": 6.0,
    "STRC_MRNA_DECAY_S": 1e-4,
    "K_TRANSLATION": 5e-4,
    "K_PROT_DECAY_S": 5e-6,
}

# Parameters to perturb (load-bearing pharmacology)
SWEEP_PARAMS = [
    "K_CA_AC1_NM",           # AC1 Ca sensitivity
    "AC1_VMAX_NM_S",         # AC1 max rate
    "PDE4_VMAX_NM_S",        # cAMP degradation rate
    "K_CAMP_PKA_NM",         # PKA activation threshold
    "K_CREB_DEPHOS_S",       # CREB dephosphorylation rate
    "K_TXN_MAX_FOLD",        # maximal CRE induction
    "CA_CEILING_NM",         # Ca drive saturation
    "SPL_HALF_SATURATION",   # dB at which Ca=half-max
]

SPL_LEVELS = [30, 45, 75, 100]  # dB


# --------------------------------------------------------------------------- #
# ODE
# --------------------------------------------------------------------------- #
def ca_forcing(t, SPL_dB, am_freq_hz, P):
    if SPL_dB <= 30:
        base = P["CA_REST_NM"]
    else:
        x = (SPL_dB - P["SPL_HALF_SATURATION"]) / 15.0
        ca_peak = P["CA_REST_NM"] + (P["CA_CEILING_NM"] - P["CA_REST_NM"]) / (
            1.0 + np.exp(-x))
        base = ca_peak
    if am_freq_hz <= 0:
        return base
    phase = (t * am_freq_hz) % 1.0
    gate = 1.0 if phase < 0.5 else 0.1
    return P["CA_REST_NM"] + (base - P["CA_REST_NM"]) * gate


def rhs(t, y, SPL_dB, am_freq_hz, P):
    Ca_drive, AC1, cAMP, PKAc, CREB_P, mRNA, prot = y

    Ca_target = ca_forcing(t, SPL_dB, am_freq_hz, P)
    dCa = P["K_EXTRUSION_S"] * (Ca_target - Ca_drive)

    ca4cam_frac = (Ca_drive ** P["N_CA_AC1"]) / (
        Ca_drive ** P["N_CA_AC1"] + P["K_CA_AC1_NM"] ** P["N_CA_AC1"])
    dAC1 = P["K_AC1_ACT_S"] * ca4cam_frac * (1 - AC1) - P["K_AC1_INACT_S"] * AC1

    prod = P["CAMP_BASAL_PROD_NM_S"] + P["AC1_VMAX_NM_S"] * AC1
    deg = P["PDE4_VMAX_NM_S"] * cAMP / (cAMP + P["PDE4_KM_NM"])
    dcAMP = prod - deg

    pka_target = (cAMP ** P["N_PKA"]) / (
        cAMP ** P["N_PKA"] + P["K_CAMP_PKA_NM"] ** P["N_PKA"])
    dPKAc = P["K_PKA_ACT_S"] * (pka_target - PKAc)

    dCREB_P = (P["K_CREB_PHOS_S"] * PKAc * (1 - CREB_P)
               - P["K_CREB_DEPHOS_S"] * CREB_P)

    cre_drive = (CREB_P ** P["N_CRE"]) / (
        CREB_P ** P["N_CRE"] + P["K_CREB_CRE_HALF"] ** P["N_CRE"])
    txn_rate = P["K_STRC_TXN_BASAL_S"] * (
        1 + (P["K_TXN_MAX_FOLD"] - 1) * cre_drive)
    dmRNA = txn_rate - P["STRC_MRNA_DECAY_S"] * mRNA

    dprot = P["K_TRANSLATION"] * mRNA - P["K_PROT_DECAY_S"] * prot

    return np.array([dCa, dAC1, dcAMP, dPKAc, dCREB_P, dmRNA, dprot])


def steady_state(SPL_dB, P, am_freq_hz=0.0, t_hours=20.0):
    y0 = np.array([
        P["CA_REST_NM"], 0.05, 100.0, 0.05, 0.02,
        P["K_STRC_TXN_BASAL_S"] / P["STRC_MRNA_DECAY_S"],
        (P["K_STRC_TXN_BASAL_S"] / P["STRC_MRNA_DECAY_S"])
        * P["K_TRANSLATION"] / P["K_PROT_DECAY_S"],
    ])
    t_end = t_hours * 3600
    max_step = 60.0 if am_freq_hz == 0 else min(20.0, 0.1 / am_freq_hz)
    sol = solve_ivp(
        rhs, [0, t_end], y0, args=(SPL_dB, am_freq_hz, P),
        method="LSODA", rtol=1e-3, atol=1e-6, max_step=max_step,
        dense_output=False,
    )
    n_tail = max(10, sol.t.size // 10)
    y_tail = sol.y[:, -n_tail:].mean(axis=1)
    return {
        "Ca_mean_nM": float(y_tail[0]),
        "AC1_mean": float(y_tail[1]),
        "cAMP_mean_nM": float(y_tail[2]),
        "PKAc_mean": float(y_tail[3]),
        "CREB_P_mean": float(y_tail[4]),
        "mRNA_mean": float(y_tail[5]),
        "protein_mean": float(y_tail[6]),
        "success": bool(sol.success),
    }


# --------------------------------------------------------------------------- #
# Sweep
# --------------------------------------------------------------------------- #
def evaluate(P, tag=""):
    """Return dict with protein fold-changes, CREB range, monotonicity."""
    results = {}
    for spl in SPL_LEVELS:
        r = steady_state(spl, P)
        results[spl] = r
    proteins = [results[s]["protein_mean"] for s in SPL_LEVELS]
    creb = [results[s]["CREB_P_mean"] for s in SPL_LEVELS]

    tol = 0.005 * max(proteins)
    monotonic = all(proteins[i] <= proteins[i+1] + tol
                    for i in range(len(proteins)-1))
    fold_45 = proteins[1] / proteins[0] if proteins[0] > 0 else np.inf
    fold_100 = proteins[-1] / proteins[0] if proteins[0] > 0 else np.inf
    creb_fold = (max(creb) / max(min(creb), 1e-4))
    return {
        "tag": tag,
        "proteins": [round(p, 2) for p in proteins],
        "creb_p": [round(c, 4) for c in creb],
        "ca_mean_nM": [round(results[s]["Ca_mean_nM"], 1)
                       for s in SPL_LEVELS],
        "ac1_mean": [round(results[s]["AC1_mean"], 3)
                     for s in SPL_LEVELS],
        "fold_silence_to_45dB": round(fold_45, 2),
        "fold_silence_to_100dB": round(fold_100, 2),
        "creb_p_fold": round(creb_fold, 1),
        "monotonic_5_tol": monotonic,
    }


def main():
    print("=" * 80)
    print("Ca²⁺ Phase 2 sensitivity sweep")
    print("=" * 80)
    print(f"SPL levels: {SPL_LEVELS} dB")
    print(f"Sweeping {len(SWEEP_PARAMS)} parameters × (0.5×, 2×) = "
          f"{2*len(SWEEP_PARAMS)+1} runs\n")

    # Baseline
    print("[baseline]")
    base_res = evaluate(BASELINE, tag="baseline")
    print(f"  proteins @ SPL {SPL_LEVELS}: {base_res['proteins']}")
    print(f"  fold silence→45dB = {base_res['fold_silence_to_45dB']}×  "
          f"silence→100dB = {base_res['fold_silence_to_100dB']}×")
    print(f"  CREB-P fold = {base_res['creb_p_fold']}×  "
          f"monotonic = {base_res['monotonic_5_tol']}\n")

    sweep_data = {"baseline": base_res}

    for p in SWEEP_PARAMS:
        for factor_label, factor in [("0.5x", 0.5), ("2x", 2.0)]:
            tag = f"{p}_{factor_label}"
            P = copy.deepcopy(BASELINE)
            P[p] = BASELINE[p] * factor
            print(f"[{tag}] {p} {BASELINE[p]:.3g} → {P[p]:.3g}")
            r = evaluate(P, tag=tag)
            sweep_data[tag] = r
            fold_ratio_vs_base = r["fold_silence_to_45dB"] / max(
                base_res["fold_silence_to_45dB"], 1e-3)
            flag = ""
            if not r["monotonic_5_tol"]:
                flag = "  ⚠ NON-MONOTONIC"
            elif r["fold_silence_to_45dB"] < 1.2:
                flag = "  ⚠ WEAK (<1.2× at 45 dB)"
            print(f"  proteins: {r['proteins']}  "
                  f"fold45={r['fold_silence_to_45dB']:.2f}×  "
                  f"fold100={r['fold_silence_to_100dB']:.2f}×  "
                  f"(Δvs_base={fold_ratio_vs_base-1:+.1%}){flag}")

    # Robustness summary: which params break the hypothesis?
    broken = [tag for tag, r in sweep_data.items()
              if tag != "baseline" and (
                  not r["monotonic_5_tol"]
                  or r["fold_silence_to_45dB"] < 1.2)]
    fragile = []
    for tag, r in sweep_data.items():
        if tag == "baseline":
            continue
        ratio = r["fold_silence_to_45dB"] / base_res["fold_silence_to_45dB"]
        if ratio < 0.7 or ratio > 1.5:
            fragile.append((tag, round((ratio - 1) * 100, 1)))

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Baseline fold silence→45 dB: {base_res['fold_silence_to_45dB']}×")
    print(f"Variants that BREAK hypothesis "
          f"(non-monotonic or fold<1.2×): {len(broken)}")
    for t in broken:
        print(f"   - {t}")
    print(f"\nVariants with >30% deviation in fold change: {len(fragile)}")
    for t, pct in fragile:
        print(f"   - {t}: {pct:+.1f}%")

    robustness_score = 1.0 - (len(broken) + 0.5 * len(fragile)) / (
        2 * len(SWEEP_PARAMS))
    print(f"\nRobustness index (1.0 = perfect): {robustness_score:.2f}")

    verdict = "ROBUST" if robustness_score >= 0.7 else (
        "FRAGILE" if robustness_score >= 0.4 else "UNSTABLE")
    print(f"Verdict: {verdict}")

    # Write JSON
    out = {
        "baseline_params": BASELINE,
        "sweep_params": SWEEP_PARAMS,
        "SPL_levels_dB": SPL_LEVELS,
        "sweep_results": sweep_data,
        "broken_variants": broken,
        "fragile_variants": fragile,
        "robustness_index": round(robustness_score, 3),
        "verdict": (
            f"{verdict} — {len(broken)} of {2*len(SWEEP_PARAMS)} "
            f"perturbations break hypothesis; "
            f"{len(fragile)} have >30% fold-change deviation."
        ),
    }
    out_json = OUT_DIR / "ca_oscillation_phase2_sensitivity.json"
    out_json.write_text(json.dumps(out, indent=2))
    print(f"\n[write] {out_json.name}")

    # Plot: tornado chart of 45dB-fold change for all perturbations
    labels = ["baseline"]
    fold45 = [base_res["fold_silence_to_45dB"]]
    for p in SWEEP_PARAMS:
        for fac in ("0.5x", "2x"):
            tag = f"{p}_{fac}"
            labels.append(tag.replace("_", " "))
            fold45.append(sweep_data[tag]["fold_silence_to_45dB"])

    fig, ax = plt.subplots(2, 1, figsize=(13, 10))

    # Tornado: fold silence → 45 dB
    y = np.arange(len(labels))
    colors = []
    for lbl, f in zip(labels, fold45):
        if lbl == "baseline":
            colors.append("tab:green")
        elif f < 1.2:
            colors.append("tab:red")
        elif f > 2.5:
            colors.append("tab:orange")
        else:
            colors.append("tab:blue")
    ax[0].barh(y, fold45, color=colors, edgecolor="black", alpha=0.8)
    ax[0].axvline(base_res["fold_silence_to_45dB"], color="tab:green",
                  ls="--", lw=1, label=f"baseline {base_res['fold_silence_to_45dB']}×")
    ax[0].axvline(1.0, color="black", ls=":", lw=1, label="no change (1×)")
    ax[0].axvline(1.2, color="tab:red", ls=":", lw=1, label="weak threshold")
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(labels, fontsize=7)
    ax[0].set_xlabel("fold-change silence → 45 dB (STRC protein)")
    ax[0].set_title("Sensitivity: which parameters break the 45-dB response?")
    ax[0].legend(loc="lower right", fontsize=8)
    ax[0].grid(alpha=0.3, axis="x")

    # Per-SPL curves for baseline vs worst perturbation
    worst_tag = None
    worst_45 = base_res["fold_silence_to_45dB"]
    for tag, r in sweep_data.items():
        if tag != "baseline" and r["fold_silence_to_45dB"] < worst_45:
            worst_tag = tag
            worst_45 = r["fold_silence_to_45dB"]
    best_tag = None
    best_45 = base_res["fold_silence_to_45dB"]
    for tag, r in sweep_data.items():
        if tag != "baseline" and r["fold_silence_to_45dB"] > best_45:
            best_tag = tag
            best_45 = r["fold_silence_to_45dB"]

    ax[1].plot(SPL_LEVELS, base_res["proteins"], "o-",
               label=f"baseline (fold45={base_res['fold_silence_to_45dB']}×)",
               lw=2, color="tab:green")
    if worst_tag:
        ax[1].plot(SPL_LEVELS, sweep_data[worst_tag]["proteins"], "s--",
                   label=f"{worst_tag} (fold45={sweep_data[worst_tag]['fold_silence_to_45dB']}×)",
                   color="tab:red")
    if best_tag:
        ax[1].plot(SPL_LEVELS, sweep_data[best_tag]["proteins"], "^--",
                   label=f"{best_tag} (fold45={sweep_data[best_tag]['fold_silence_to_45dB']}×)",
                   color="tab:orange")
    ax[1].set_xlabel("SPL (dB)")
    ax[1].set_ylabel("STRC protein (a.u.)")
    ax[1].set_title("Dose-response: baseline vs worst / best perturbations")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / "ca_oscillation_phase2_sensitivity.png", dpi=160)
    print(f"[write] ca_oscillation_phase2_sensitivity.png")


if __name__ == "__main__":
    main()
