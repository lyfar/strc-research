"""
STRC mRNA-LNP PK/PD integration v2 — OHC tropism realism + ER/UPR ceiling.

Refines v1 ([[STRC mRNA-LNP Strategy B Full-Length]]) per
[[STRC LNP Cochlear Tropism Literature Scan]] 2026-04-22 next step:

  (a) 1-3% OHC tropism row: between v1's "untargeted 0.8%" and the aspirational
      "cochlear-tropic 5%". Per the literature scan, the 1-3% range is the
      realistic plausible-but-unbuilt window: peptide-prestin-A666 LNP-mRNA
      hybrids could reach this in 12-18 mo of engineering work; A666 alone
      with drug payload (Wang 2018) demonstrates ~1% OHC targeting in vivo.
  (b) ER/UPR ~5× soft ceiling on STRC fold per OHC: chronic mRNA dosing that
      pushes per-cell protein production beyond ~5× endogenous triggers
      unfolded-protein-response (UPR) and ER stress, capping the achievable
      benefit independent of dose. Implemented as soft tanh cap in the ODE.

Outputs everything v1 did, plus per-regimen `er_stress_risk` flag and the
new realistic-low (1%) + realistic-high (3%) tropism rows in coverage.

Run env: any conda env with numpy + scipy. ~2 min on CPU.
"""

from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

ROOT = pathlib.Path("/Users/egorlyfar/Brain/research/strc/models")
STABILITY = json.loads((ROOT / "mrna_stability_cochlear_results.json").read_text())
DOSE_RESP = json.loads((ROOT / "rbm24_mrna_dose_response_results.json").read_text())

params = DOSE_RESP["params"]
HILL_K = params["hill_km"]
HILL_N = params["hill_n"]
MAX_BOOST = params["max_boost"]
THRESHOLD = params["threshold_fold"]

MRNA_HL_UNMOD_H = STABILITY["stability"]["unmod_OHC"]["hl_h"]
MRNA_HL_M1PSI_H = STABILITY["stability"]["m1psi_OHC"]["hl_h"]

RBM24_HL_D = params["rbm24_hl_d"]
STRC_HL_D = params["strc_hl_d"]

LNP_ENDO_ESCAPE = 0.02
LNP_UNTARGETED = 0.008
LNP_REALISTIC_LOW = 0.01            # NEW v2 — 1% OHC tropism, realistic floor
LNP_REALISTIC_HIGH = 0.03           # NEW v2 — 3% OHC tropism, realistic ceiling
LNP_COCHLEAR_TROPIC = 0.05          # v1 — aspirational
LNP_OHC_TARGETED = 0.20             # v1 — sci-fi for current chemistry

# NEW v2 — ER stress soft ceiling on STRC fold per OHC.
# 5× endogenous is the upper bound before unfolded-protein-response kicks in
# for chronic exogenous-mRNA-driven protein production. Past that, BiP/Grp78
# binds nascent stereocilin, ATF6/IRE1/PERK arms activate, translation rate
# drops globally (PERK/eIF2α), and the cell prioritizes folding over output.
# Source: chronic-mRNA UPR literature in muscle and liver (Sahin 2014,
# Hadas 2019); cochlea-specific data is missing but OHC ER capacity is
# almost certainly LOWER than hepatocytes given OHC's small cytoplasm and
# heavy mechanical load.
ER_CEILING_FOLD = 5.0

# Soft-cap activation threshold: regimens with peak per-OHC fold above this
# fraction of ER_CEILING_FOLD are flagged as "approaching ER stress".
ER_STRESS_FLAG_THRESHOLD = 0.8 * ER_CEILING_FOLD  # 4.0×

K_TRANSLATE = 1.0
K_STRC_BASAL = np.log(2) / (STRC_HL_D * 24.0 * 3600.0)


def _er_capped_target(f_target_uncapped: float) -> float:
    """Soft tanh cap on STRC fold target.

    Maps:
      f_target_uncapped = 1 -> 1
      f_target_uncapped = 5 (ceiling) -> ~4.05
      f_target_uncapped -> infinity -> ER_CEILING_FOLD
    """
    excess = f_target_uncapped - 1.0
    if excess <= 0.0:
        return f_target_uncapped
    span = ER_CEILING_FOLD - 1.0
    return 1.0 + span * np.tanh(excess / span)


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
    tx_peak_fold: float = field(init=False, default=np.nan)
    tx_trough_fold: float = field(init=False, default=np.nan)
    tx_mean_fold: float = field(init=False, default=np.nan)
    tx_time_above_threshold_frac: float = field(init=False, default=np.nan)
    cochlea_peak_fold: float = field(init=False, default=np.nan)
    cochlea_trough_fold: float = field(init=False, default=np.nan)
    cochlea_mean_fold: float = field(init=False, default=np.nan)
    therapeutic_per_OHC: bool = field(init=False, default=False)
    therapeutic_cochlea: bool = field(init=False, default=False)
    er_stress_risk: bool = field(init=False, default=False)

    def __post_init__(self) -> None:
        self.annual_doses = int(round(365.0 / self.interval_d))
        self.annual_intra_total = self.annual_doses * self.dose_intra
        self.annual_extra_total = self.annual_intra_total / LNP_ENDO_ESCAPE


def ode_pkpd(
    t_eval_s: np.ndarray,
    dose_times_s: list[float],
    dose_intra: float,
    mrna_hl_h: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    k_m = np.log(2.0) / (mrna_hl_h * 3600.0)
    k_p = np.log(2.0) / (RBM24_HL_D * 24.0 * 3600.0)
    k_f = K_STRC_BASAL

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        m, p, f = y
        hill = (p**HILL_N) / (HILL_K**HILL_N + p**HILL_N)
        f_target_uncapped = 1.0 + (MAX_BOOST - 1.0) * hill
        f_target = _er_capped_target(f_target_uncapped)
        dm = -k_m * m
        dp = K_TRANSLATE * m - k_p * p
        df = k_f * (f_target - f)
        return np.array([dm, dp, df])

    y = np.array([0.0, 0.0, 1.0])
    segments_t: list[np.ndarray] = []
    segments_f: list[np.ndarray] = []
    t_breaks = sorted(set([float(t_eval_s[0])] + dose_times_s + [float(t_eval_s[-1])]))
    t_breaks = [t for t in t_breaks if t_eval_s[0] <= t <= t_eval_s[-1]]
    for i in range(len(t_breaks) - 1):
        t0 = t_breaks[i]
        t1 = t_breaks[i + 1]
        if t0 in dose_times_s:
            y = y.copy()
            y[0] += dose_intra
        t_span = (t0, t1)
        mask = (t_eval_s >= t0) & (t_eval_s <= t1)
        t_seg = t_eval_s[mask]
        if len(t_seg) == 0:
            continue
        sol = solve_ivp(rhs, t_span, y, t_eval=t_seg, method="LSODA", rtol=1e-6, atol=1e-9)
        segments_t.append(sol.t)
        segments_f.append(sol.y[2])
        y = sol.y[:, -1]

    t_out = np.concatenate(segments_t) if segments_t else np.array([])
    f_out = np.concatenate(segments_f) if segments_f else np.array([])
    return t_out, np.zeros_like(f_out), f_out


def cochlea_fold(f_tx: np.ndarray, eff_ohc: float) -> np.ndarray:
    return eff_ohc * f_tx + (1.0 - eff_ohc) * 1.0


def simulate_regimen(reg: Regimen, horizon_d: float = 365.0) -> dict:
    horizon_s = horizon_d * 24.0 * 3600.0
    t_eval = np.linspace(0.0, horizon_s, 4000)
    dose_times = [i * reg.interval_d * 24.0 * 3600.0 for i in range(reg.annual_doses + 1)]
    dose_times = [t for t in dose_times if t < horizon_s]
    _, _, f_tx = ode_pkpd(t_eval, dose_times, reg.dose_intra, reg.mrna_hl_h)
    f_cochlea = cochlea_fold(f_tx, reg.eff_ohc)

    tail_mask = t_eval >= 0.7 * horizon_s
    tx_tail = f_tx[tail_mask]
    co_tail = f_cochlea[tail_mask]
    reg.tx_peak_fold = float(np.max(tx_tail))
    reg.tx_trough_fold = float(np.min(tx_tail))
    reg.tx_mean_fold = float(np.mean(tx_tail))
    reg.tx_time_above_threshold_frac = float(np.mean(tx_tail >= THRESHOLD))
    reg.cochlea_peak_fold = float(np.max(co_tail))
    reg.cochlea_trough_fold = float(np.min(co_tail))
    reg.cochlea_mean_fold = float(np.mean(co_tail))
    reg.therapeutic_per_OHC = reg.tx_trough_fold >= THRESHOLD
    reg.therapeutic_cochlea = reg.cochlea_trough_fold >= THRESHOLD
    reg.er_stress_risk = reg.tx_peak_fold >= ER_STRESS_FLAG_THRESHOLD

    return {
        "t_d": (t_eval / 86400.0).tolist(),
        "f_tx": f_tx.tolist(),
        "f_cochlea": f_cochlea.tolist(),
    }


TARGETS = {
    "untargeted": LNP_UNTARGETED,
    "realistic_1pct": LNP_REALISTIC_LOW,
    "realistic_3pct": LNP_REALISTIC_HIGH,
    "cochlear_tropic": LNP_COCHLEAR_TROPIC,
    "OHC_targeted": LNP_OHC_TARGETED,
}


def sweep_schedules() -> list[Regimen]:
    intervals = [14.0, 21.0, 28.0, 42.0, 56.0, 84.0]
    doses = [200.0, 500.0, 1000.0, 2000.0, 5000.0]
    mods = {"unmod": MRNA_HL_UNMOD_H, "m1psi": MRNA_HL_M1PSI_H}

    regimens: list[Regimen] = []
    for mod_name, mod_hl in mods.items():
        for tgt_name, tgt_eff in TARGETS.items():
            for T in intervals:
                for D in doses:
                    regimens.append(Regimen(
                        modification=mod_name, mrna_hl_h=mod_hl,
                        targeting=tgt_name, eff_ohc=tgt_eff,
                        interval_d=T, dose_intra=D,
                    ))
    return regimens


def minimum_therapeutic_per_OHC(regs: list[Regimen]) -> dict:
    out: dict = {}
    for mod in ["unmod", "m1psi"]:
        out[mod] = {}
        for T in sorted({r.interval_d for r in regs}):
            viable = [r for r in regs
                      if r.modification == mod and r.interval_d == T
                      and r.targeting == "untargeted" and r.therapeutic_per_OHC]
            if not viable:
                out[mod][f"T_{int(T)}d"] = None
                continue
            r_min = min(viable, key=lambda x: x.dose_intra)
            out[mod][f"T_{int(T)}d"] = {
                "min_dose_intra": r_min.dose_intra,
                "annual_intra_total_per_OHC": r_min.annual_intra_total,
                "annual_extra_total_per_OHC": r_min.annual_extra_total,
                "tx_trough_fold": r_min.tx_trough_fold,
                "tx_peak_fold": r_min.tx_peak_fold,
                "tx_tat": r_min.tx_time_above_threshold_frac,
                "annual_doses": r_min.annual_doses,
                "er_stress_risk": r_min.er_stress_risk,
            }
    return out


def coverage_summary(regs: list[Regimen]) -> dict:
    out: dict = {}
    for tgt, eff in TARGETS.items():
        subset = [r for r in regs if r.targeting == tgt and r.therapeutic_per_OHC]
        if not subset:
            out[tgt] = {
                "lnp_eff_ohc": eff,
                "max_cochlea_mean_fold": None,
                "max_cochlea_trough_fold": None,
                "theoretical_ceiling_cochlea_mean": eff * MAX_BOOST + (1.0 - eff),
                "n_therapeutic_regimens": 0,
                "n_er_stress_flagged": 0,
            }
            continue
        out[tgt] = {
            "lnp_eff_ohc": eff,
            "max_cochlea_mean_fold": max(r.cochlea_mean_fold for r in subset),
            "max_cochlea_trough_fold": max(r.cochlea_trough_fold for r in subset),
            "theoretical_ceiling_cochlea_mean": eff * MAX_BOOST + (1.0 - eff),
            "n_therapeutic_regimens": len(subset),
            "n_er_stress_flagged": sum(1 for r in subset if r.er_stress_risk),
        }
    return out


def main() -> None:
    regs = sweep_schedules()
    print(f"v2 sim: {len(regs)} regimens "
          f"({len(TARGETS)} tropisms × 2 mods × 6 intervals × 5 doses)", flush=True)

    for i, reg in enumerate(regs):
        simulate_regimen(reg)
        if i % 60 == 0:
            er = "⚠" if reg.er_stress_risk else " "
            print(
                f"  [{i+1:3d}/{len(regs)}] {reg.modification:6s} {reg.targeting:15s} "
                f"T={reg.interval_d:4.0f}d D={reg.dose_intra:6.0f} "
                f"tx_trough={reg.tx_trough_fold:.2f} peak={reg.tx_peak_fold:.2f} "
                f"ther={reg.therapeutic_per_OHC} {er}ER",
                flush=True,
            )

    out = {
        "model": "STRC mRNA-LNP PK/PD integration v2 (1-3% tropism + ER/UPR cap)",
        "params": {
            "hill_K": HILL_K, "hill_n": HILL_N, "max_boost": MAX_BOOST,
            "threshold_fold": THRESHOLD,
            "mrna_hl_unmod_OHC_h": MRNA_HL_UNMOD_H,
            "mrna_hl_m1psi_OHC_h": MRNA_HL_M1PSI_H,
            "rbm24_hl_d": RBM24_HL_D, "strc_hl_d": STRC_HL_D,
            "lnp_endo_escape": LNP_ENDO_ESCAPE,
            "er_ceiling_fold": ER_CEILING_FOLD,
            "er_stress_flag_threshold": ER_STRESS_FLAG_THRESHOLD,
            "targeting": TARGETS,
        },
        "regimens": [
            {
                "modification": r.modification, "targeting": r.targeting,
                "interval_d": r.interval_d, "dose_intra": r.dose_intra,
                "annual_doses": r.annual_doses,
                "annual_intra_total_per_OHC": r.annual_intra_total,
                "annual_extra_total_per_OHC": r.annual_extra_total,
                "tx_trough_fold": r.tx_trough_fold,
                "tx_peak_fold": r.tx_peak_fold,
                "tx_mean_fold": r.tx_mean_fold,
                "tx_tat": r.tx_time_above_threshold_frac,
                "cochlea_trough_fold": r.cochlea_trough_fold,
                "cochlea_peak_fold": r.cochlea_peak_fold,
                "cochlea_mean_fold": r.cochlea_mean_fold,
                "therapeutic_per_OHC": r.therapeutic_per_OHC,
                "therapeutic_cochlea": r.therapeutic_cochlea,
                "er_stress_risk": r.er_stress_risk,
            }
            for r in regs
        ],
        "minimum_therapeutic_per_OHC": minimum_therapeutic_per_OHC(regs),
        "coverage_summary": coverage_summary(regs),
    }

    out_path = ROOT / "mrna_lnp_pkpd_integration_v2_ohc_tropism.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {out_path}")
    print(f"coverage_summary: {json.dumps(out['coverage_summary'], indent=2)}")


if __name__ == "__main__":
    main()
