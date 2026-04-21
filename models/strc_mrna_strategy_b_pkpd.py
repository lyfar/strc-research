"""
Strategy B (full-length STRC mRNA-LNP) PK/PD integration model.

Parallel to `mrna_lnp_pkpd_integration.py` (Strategy A) but:
  - 2-compartment ODE (no RBM24 intermediate): mRNA -> STRC protein directly
  - No Hill saturation: linear dose-response in STRC
  - Calibrated against WT endogenous STRC SS (s_WT = 3021 units from Strategy A)
  - Scenarios: WT reference / DFNB16 null / Misha compound het

Key differences from Strategy A:
  1. Works on NULL alleles (Misha paternal 98 kb deletion). Strategy A
     (RBM24 splicing) requires residual pre-mRNA substrate; Strategy B does not.
  2. No Hill cap on STRC fold. Dose can compensate for low LNP efficiency,
     unlike Strategy A where cochlea-mean is bounded by eff * 3 + (1-eff).
  3. Per-OHC therapeutic dose is SUBSTANTIALLY higher than Strategy A because
     direct protein replacement lacks the catalytic amplification of RBM24
     splicing (1 RBM24 protein corrects many STRC transcripts).
  4. Payload size ~4x larger (~6 kb STRC CDS vs ~1.5 kb RBM24).

Calibration:
  - Endogenous WT STRC mRNA SS = 50 copies/OHC (biological estimate)
  - Endogenous WT STRC protein SS = 3021 units (from Strategy A baseline)
  - K_TRANSLATE_STRC = s_WT * k_s / ENDO_MRNA_SS
  - STRC protein t1/2 = 14 d (inherited)

Model per OHC:
    dm/dt = pulse(t; D, T) - k_m * m
    ds/dt = K_TRANSLATE_STRC * m - k_s * s
    (where s is EXOGENOUS STRC from mRNA therapy;
     endogenous STRC contributes s_endo_frac * s_WT baseline,
     scenario-dependent)

Population layer (cochlea-mean fold over s_WT):
    transfected OHC:    total = s_endo_frac + s_exog/s_WT
    non-transfected:    total = s_endo_frac
    cochlea_mean = eff * (s_endo_frac + s_exog_tx/s_WT) + (1-eff) * s_endo_frac
                 = s_endo_frac + eff * s_exog_tx/s_WT

Notably: cochlea-mean ceiling is UNBOUNDED in dose (no Hill cap).
Therapeutic threshold per-OHC: s_exog_trough/s_WT >= 1.0 (1x WT) or >= 2.0 (2x WT).

Sweeps: 2 mods x 3 targetings x 6 intervals x 8 doses = 288 regimens, 365d horizon.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass, field, asdict

import numpy as np
from scipy.integrate import solve_ivp

ROOT = pathlib.Path("/Users/egorlyfar/Brain/research/strc/models")
STABILITY = json.loads((ROOT / "mrna_stability_cochlear_results.json").read_text())
STRATEGY_A = json.loads((ROOT / "rbm24_mrna_dose_response_results.json").read_text())

S_WT = STRATEGY_A["baseline_strc_prot_ss"]           # 3021.1 units (arb.)
STRC_HL_D = STRATEGY_A["params"]["strc_hl_d"]        # 14 d
MRNA_HL_UNMOD_H = STABILITY["stability"]["unmod_OHC"]["hl_h"]   # 4 h
MRNA_HL_M1PSI_H = STABILITY["stability"]["m1psi_OHC"]["hl_h"]   # 12 h
LNP_ENDO_ESCAPE = 0.02

ENDO_MRNA_SS = 50.0
K_STRC = np.log(2.0) / (STRC_HL_D * 24.0 * 3600.0)
K_TRANSLATE_STRC = S_WT * K_STRC / ENDO_MRNA_SS

THRESHOLD_1X = 1.0
THRESHOLD_2X = 2.0

LNP_UNTARGETED = 0.008
LNP_COCHLEAR_TROPIC = 0.05
LNP_OHC_TARGETED = 0.20

SCENARIOS = {
    "WT_ref":      {"s_endo_frac": 1.0, "note": "homozygous WT, therapy adds on top"},
    "DFNB16_null": {"s_endo_frac": 0.0, "note": "biallelic null (paternal 98kb del/del)"},
    "misha":       {"s_endo_frac": 0.15, "note": "paternal null + maternal E1659A residual"},
}


@dataclass
class Regimen:
    modification: str
    mrna_hl_h: float
    targeting: str
    eff_ohc: float
    interval_d: float
    dose_intra: float

    annual_doses: int = field(init=False)
    annual_intra_total: float = field(init=False)
    annual_extra_total: float = field(init=False)
    tx_s_exog_over_s_wt_trough: float = field(init=False, default=np.nan)
    tx_s_exog_over_s_wt_peak: float = field(init=False, default=np.nan)
    tx_s_exog_over_s_wt_mean: float = field(init=False, default=np.nan)
    therapeutic_per_OHC_1x: bool = field(init=False, default=False)
    therapeutic_per_OHC_2x: bool = field(init=False, default=False)

    def __post_init__(self) -> None:
        self.annual_doses = int(round(365.0 / self.interval_d))
        self.annual_intra_total = self.annual_doses * self.dose_intra
        self.annual_extra_total = self.annual_intra_total / LNP_ENDO_ESCAPE


def ode_strategy_b(
    t_eval_s: np.ndarray,
    dose_times_s: list[float],
    dose_intra: float,
    mrna_hl_h: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate 2-compartment ODE: mRNA -> STRC protein with delta-pulse doses."""
    k_m = np.log(2.0) / (mrna_hl_h * 3600.0)

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        m, s = y
        return np.array([-k_m * m, K_TRANSLATE_STRC * m - K_STRC * s])

    y = np.array([0.0, 0.0])
    segments_t: list[np.ndarray] = []
    segments_s: list[np.ndarray] = []

    t_breaks = sorted(set([float(t_eval_s[0])] + list(dose_times_s) + [float(t_eval_s[-1])]))
    t_breaks = [t for t in t_breaks if t_eval_s[0] <= t <= t_eval_s[-1]]

    for i in range(len(t_breaks) - 1):
        t0 = t_breaks[i]
        t1 = t_breaks[i + 1]
        if t0 in dose_times_s:
            y = y.copy()
            y[0] += dose_intra
        mask = (t_eval_s >= t0) & (t_eval_s <= t1)
        t_seg = t_eval_s[mask]
        if len(t_seg) == 0:
            continue
        sol = solve_ivp(rhs, (t0, t1), y, t_eval=t_seg, method="LSODA", rtol=1e-6, atol=1e-9)
        segments_t.append(sol.t)
        segments_s.append(sol.y[1])
        y = sol.y[:, -1]

    t_out = np.concatenate(segments_t) if segments_t else np.array([])
    s_out = np.concatenate(segments_s) if segments_s else np.array([])
    return t_out, s_out


def simulate_regimen(reg: Regimen, horizon_d: float = 365.0) -> dict:
    horizon_s = horizon_d * 24.0 * 3600.0
    t_eval = np.linspace(0.0, horizon_s, 4000)
    dose_times = [i * reg.interval_d * 24.0 * 3600.0 for i in range(reg.annual_doses + 1)]
    dose_times = [t for t in dose_times if t < horizon_s]
    _, s_exog = ode_strategy_b(t_eval, dose_times, reg.dose_intra, reg.mrna_hl_h)
    s_exog_over_s_wt = s_exog / S_WT

    tail_mask = t_eval >= 0.7 * horizon_s
    tail = s_exog_over_s_wt[tail_mask]
    reg.tx_s_exog_over_s_wt_trough = float(np.min(tail))
    reg.tx_s_exog_over_s_wt_peak = float(np.max(tail))
    reg.tx_s_exog_over_s_wt_mean = float(np.mean(tail))
    reg.therapeutic_per_OHC_1x = reg.tx_s_exog_over_s_wt_trough >= THRESHOLD_1X
    reg.therapeutic_per_OHC_2x = reg.tx_s_exog_over_s_wt_trough >= THRESHOLD_2X

    return {
        "t_d": (t_eval / 86400.0).tolist(),
        "s_exog_over_s_wt": s_exog_over_s_wt.tolist(),
    }


def sweep_schedules() -> list[Regimen]:
    intervals = [14.0, 21.0, 28.0, 42.0, 56.0, 84.0]
    doses = [200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0, 30000.0, 100000.0]
    mods = {"unmod": MRNA_HL_UNMOD_H, "m1psi": MRNA_HL_M1PSI_H}
    targets = {
        "untargeted": LNP_UNTARGETED,
        "cochlear_tropic": LNP_COCHLEAR_TROPIC,
        "OHC_targeted": LNP_OHC_TARGETED,
    }

    regs: list[Regimen] = []
    for mod_name, mod_hl in mods.items():
        for tgt_name, tgt_eff in targets.items():
            for T in intervals:
                for D in doses:
                    regs.append(Regimen(
                        modification=mod_name,
                        mrna_hl_h=mod_hl,
                        targeting=tgt_name,
                        eff_ohc=tgt_eff,
                        interval_d=T,
                        dose_intra=D,
                    ))
    return regs


def minimum_therapeutic_per_OHC(regs: list[Regimen], flag_attr: str) -> dict:
    out: dict = {}
    for mod in ["unmod", "m1psi"]:
        out[mod] = {}
        for T in sorted({r.interval_d for r in regs}):
            viable = [
                r for r in regs
                if r.modification == mod
                and r.interval_d == T
                and r.targeting == "untargeted"
                and getattr(r, flag_attr)
            ]
            if not viable:
                out[mod][f"T_{int(T)}d"] = None
                continue
            r_min = min(viable, key=lambda x: x.dose_intra)
            out[mod][f"T_{int(T)}d"] = {
                "min_dose_intra": r_min.dose_intra,
                "annual_intra_total_per_OHC": r_min.annual_intra_total,
                "annual_extra_total_per_OHC": r_min.annual_extra_total,
                "tx_trough_s_exog_over_s_wt": r_min.tx_s_exog_over_s_wt_trough,
                "tx_peak_s_exog_over_s_wt": r_min.tx_s_exog_over_s_wt_peak,
                "annual_doses": r_min.annual_doses,
            }
    return out


def cochlea_mean_by_scenario(regs: list[Regimen]) -> dict:
    """For each (targeting, scenario) find the best achievable cochlea-mean
    fold-over-WT across all regimens. Cochlea-mean = s_endo_frac + eff * s_exog_tx/s_WT."""
    out: dict = {}
    for tgt in ["untargeted", "cochlear_tropic", "OHC_targeted"]:
        eff = {"untargeted": LNP_UNTARGETED,
               "cochlear_tropic": LNP_COCHLEAR_TROPIC,
               "OHC_targeted": LNP_OHC_TARGETED}[tgt]
        out[tgt] = {"eff_ohc": eff}
        subset = [r for r in regs if r.targeting == tgt]
        if not subset:
            continue
        # Best regimen by trough s_exog
        best = max(subset, key=lambda x: x.tx_s_exog_over_s_wt_trough)
        s_exog_tx = best.tx_s_exog_over_s_wt_trough
        for scen, scen_conf in SCENARIOS.items():
            s_endo_frac = scen_conf["s_endo_frac"]
            cochlea_mean = s_endo_frac + eff * s_exog_tx
            out[tgt][scen] = {
                "s_endo_frac": s_endo_frac,
                "best_dose_intra": best.dose_intra,
                "best_interval_d": best.interval_d,
                "best_modification": best.modification,
                "s_exog_tx_trough": s_exog_tx,
                "cochlea_mean_fold_over_s_wt": cochlea_mean,
                "therapeutic_cochlea_mean_1x": cochlea_mean >= 1.0,
                "therapeutic_cochlea_mean_2x": cochlea_mean >= 2.0,
            }
    return out


def compare_to_strategy_a() -> dict:
    """Head-to-head summary: minimum annual extra burden per-OHC at m1psi Q6W.

    Strategy A (from rbm24_mrna_dose_response_results.json + PKPD result):
        m1psi Q6W D=200 mol/OHC intra -> 1800 mol/yr intra = 90000 mol/yr extra
        per TRANSFECTED OHC, gives trough 2.17x STRC fold (RBM24 -> Hill sat)

    Strategy B (from this run): find minimum dose at same interval/modification
    that achieves per-OHC trough s_exog_over_s_wt >= 2.0 (equivalent to
    doubling WT STRC in the transfected OHC).
    """
    return {
        "strategy_A_m1psi_Q6W": {
            "min_dose_intra_per_OHC": 200.0,
            "annual_intra_total_per_OHC": 1800.0,
            "annual_extra_total_per_OHC": 90000.0,
            "trough_strc_fold": 2.17,
            "mechanism": "RBM24 -> Hill-saturated STRC splicing boost (catalytic)",
            "works_on_null_allele": False,
            "cochlea_mean_ceiling": "eff * 3 + (1-eff) (hard cap from Hill max_boost=3)",
        },
    }


def main() -> None:
    regs = sweep_schedules()
    print(f"simulating {len(regs)} Strategy B regimens (K_TRANSLATE_STRC={K_TRANSLATE_STRC:.3e})")

    traces_m1psi_OHC_T42: dict = {}
    for i, reg in enumerate(regs):
        trace = simulate_regimen(reg)
        if reg.modification == "m1psi" and reg.targeting == "OHC_targeted" and reg.interval_d == 42.0:
            traces_m1psi_OHC_T42[f"D_{int(reg.dose_intra)}"] = {
                "t_d": trace["t_d"][::20],
                "s_exog_over_s_wt": trace["s_exog_over_s_wt"][::20],
            }
        if i % 40 == 0:
            print(
                f"  [{i+1}/{len(regs)}] {reg.modification:6s} {reg.targeting:17s} "
                f"T={reg.interval_d:4.0f}d D={reg.dose_intra:8.0f} "
                f"s_exog_trough={reg.tx_s_exog_over_s_wt_trough:.3f} "
                f"1x={reg.therapeutic_per_OHC_1x} 2x={reg.therapeutic_per_OHC_2x}",
                flush=True,
            )

    min_ther_1x = minimum_therapeutic_per_OHC(regs, "therapeutic_per_OHC_1x")
    min_ther_2x = minimum_therapeutic_per_OHC(regs, "therapeutic_per_OHC_2x")
    coverage = cochlea_mean_by_scenario(regs)
    head_to_head = compare_to_strategy_a()

    out = {
        "model": "STRC mRNA-LNP Strategy B (full-length, 2-compartment, no Hill)",
        "params": {
            "s_WT_units": S_WT,
            "endo_mRNA_SS_copies_per_OHC": ENDO_MRNA_SS,
            "K_TRANSLATE_STRC_units_per_mRNA_copy_per_s": K_TRANSLATE_STRC,
            "strc_hl_d": STRC_HL_D,
            "mrna_hl_unmod_OHC_h": MRNA_HL_UNMOD_H,
            "mrna_hl_m1psi_OHC_h": MRNA_HL_M1PSI_H,
            "lnp_endo_escape": LNP_ENDO_ESCAPE,
            "lnp_targeting": {
                "untargeted": LNP_UNTARGETED,
                "cochlear_tropic": LNP_COCHLEAR_TROPIC,
                "OHC_targeted": LNP_OHC_TARGETED,
            },
            "scenarios": SCENARIOS,
            "threshold_1x": THRESHOLD_1X,
            "threshold_2x": THRESHOLD_2X,
        },
        "regimens": [
            {
                "modification": r.modification,
                "targeting": r.targeting,
                "interval_d": r.interval_d,
                "dose_intra": r.dose_intra,
                "annual_doses": r.annual_doses,
                "annual_intra_total_per_OHC": r.annual_intra_total,
                "annual_extra_total_per_OHC": r.annual_extra_total,
                "tx_s_exog_trough": r.tx_s_exog_over_s_wt_trough,
                "tx_s_exog_peak": r.tx_s_exog_over_s_wt_peak,
                "tx_s_exog_mean": r.tx_s_exog_over_s_wt_mean,
                "therapeutic_per_OHC_1x": r.therapeutic_per_OHC_1x,
                "therapeutic_per_OHC_2x": r.therapeutic_per_OHC_2x,
            }
            for r in regs
        ],
        "minimum_therapeutic_per_OHC_1x_WT": min_ther_1x,
        "minimum_therapeutic_per_OHC_2x_WT": min_ther_2x,
        "cochlea_mean_by_scenario": coverage,
        "strategy_A_comparison": head_to_head,
        "exemplar_traces_m1psi_OHC_T42d": traces_m1psi_OHC_T42,
    }

    out_path = ROOT / "strc_mrna_strategy_b_pkpd.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {out_path}")
    print("\nmin therapeutic per-OHC (1x WT, m1psi):")
    print(json.dumps(min_ther_1x["m1psi"], indent=2))
    print("\nmin therapeutic per-OHC (2x WT, m1psi):")
    print(json.dumps(min_ther_2x["m1psi"], indent=2))
    print("\ncochlea-mean by scenario:")
    print(json.dumps(coverage, indent=2))


if __name__ == "__main__":
    main()
