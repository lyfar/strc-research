"""
STRC mRNA-LNP PK/PD integration model.

Combines single-dose data from:
  - mrna_stability_cochlear_results.json (half-lives per compartment × modification)
  - rbm24_mrna_dose_response_results.json (peak STRC fold at single-dose sweep)

...into a chronic multi-dose schedule model. Goal: find the minimum annual
mRNA burden that keeps steady-state trough STRC fold >= 2.0x at the cochlea
level.

Model (per-OHC, deterministic 3-compartment):
    dm/dt = pulse(t) - k_m * m
    dp/dt = k_tr * m - k_p * p
    df/dt = k_f * (1 + (max_boost - 1) * Hill(p)) - k_f * f

    m  = RBM24 mRNA intracellular (molecules/OHC)
    p  = RBM24 protein (relative units)
    f  = STRC protein fold over baseline
    Hill(p) = p^n / (K^n + p^n)

Parameters come from rbm24_mrna_dose_response_results.json:
    mrna_hl_h = 12 h (m1psi in OHC)
    rbm24_hl_d = 2 d
    strc_hl_d = 14 d (dominant eigenvalue - sets persistence)
    hill_Km = 200, hill_n = 2, max_boost = 3

Pulse delivery: intracellular dose D (mol/OHC) is injected as a delta in m at
t = 0, T, 2T, ..., where T is dose interval. Extracellular dose E is inferred
from D / endosomal_escape (LNP escape ~ 2 %).

Population layer: LNP transfection efficiency eff_ohc (fraction of OHCs that
receive any dose). Fold for non-transfected OHCs = 1. Cochlea-mean fold =
eff * fold_tx + (1 - eff) * 1.

Sweeps:
  - dose interval T in {14, 21, 28, 42, 56, 84} days
  - per-dose intracellular D in {200, 500, 1000, 2000, 5000} mol/OHC
  - modification in {unmod, m1psi}
  - LNP targeting in {untargeted 0.8 %, cochlear-tropic 5 %, OHC-targeted 20 %}

For each combination: simulate 365 days, record steady-state trough/peak
fold_cochlea, compute annual extracellular dose per OHC, flag "therapeutic"
when trough >= 2.0x.

Runtime: fast (few minutes on CPU).
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

# ═══════════════════════════════════════════════════════════════════════
# POST-AUDIT 2026-04-23 — parameter provenance summary.
#
# HILL_K=200, HILL_N=2, MAX_BOOST=3 are ⚠ CIRCULAR FITS — values set so the
# ODE model hits "therapeutic" at a reasonable dose; no published RBM24→STRC
# splicing dose-response experiment in OHCs exists. Cannot fail by construction.
#
# STRC_HL_D=14, RBM24_HL_D=2 are ⚠ UNSOURCED from primary literature.
# Plausible magnitudes (ECM proteins and RBPs) but no cochlea-specific paper.
# Traced to rbm24_mrna_dose_response_results.json; generating script not present.
#
# LNP_ENDO_ESCAPE=0.02 is ⚠ unsourced specifically but consistent with
# published <10% endosomal escape range (Patel 2019 PNAS would anchor).
#
# LNP_UNTARGETED=0.008 = 96/12,000 OHC ⚠ unsourced specifically; no primary
# paper pins LNP-to-OHC efficiency. Prior "Gao 2020 PMID 32493791" citation
# in related notes was PHANTOM (COVID-19 paper). Candidate real sources:
# Gao 2018 Nature Beethoven (lipid-RNP), Yeh 2018 Nat Biomed Eng (CBE).
#
# LNP_COCHLEAR_TROPIC=0.05, LNP_OHC_TARGETED=0.20 are EXPLICITLY HYPOTHETICAL
# (next-gen / ligand-targeted LNPs). Correctly flagged.
# ═══════════════════════════════════════════════════════════════════════
params = DOSE_RESP["params"]
HILL_K = params["hill_km"]          # ⚠ CIRCULAR FIT
HILL_N = params["hill_n"]            # ⚠ CIRCULAR FIT
MAX_BOOST = params["max_boost"]      # ⚠ CIRCULAR FIT — no exp. STRC splicing ceiling
THRESHOLD = params["threshold_fold"]

MRNA_HL_UNMOD_H = STABILITY["stability"]["unmod_OHC"]["hl_h"]
MRNA_HL_M1PSI_H = STABILITY["stability"]["m1psi_OHC"]["hl_h"]

RBM24_HL_D = params["rbm24_hl_d"]    # ⚠ UNSOURCED — no OHC-RBM24 t½ measurement
STRC_HL_D = params["strc_hl_d"]      # ⚠ UNSOURCED — no stereocilin t½ measurement

LNP_ENDO_ESCAPE = 0.02    # ⚠ within published <10% range, specific cite missing
LNP_UNTARGETED = 0.008    # ⚠ unsourced specifically (prior "Gao 2020" was PHANTOM)
LNP_COCHLEAR_TROPIC = 0.05  # hypothetical next-gen
LNP_OHC_TARGETED = 0.20     # hypothetical ligand-LNP

K_TRANSLATE = 1.0
K_STRC_BASAL = np.log(2) / (STRC_HL_D * 24.0 * 3600.0)


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
    """Integrate the 3-compartment PKPD ODE with delta-pulse doses in m."""
    k_m = np.log(2.0) / (mrna_hl_h * 3600.0)
    k_p = np.log(2.0) / (RBM24_HL_D * 24.0 * 3600.0)
    k_f = K_STRC_BASAL

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        m, p, f = y
        hill = (p**HILL_N) / (HILL_K**HILL_N + p**HILL_N)
        f_target = 1.0 + (MAX_BOOST - 1.0) * hill
        dm = -k_m * m
        dp = K_TRANSLATE * m - k_p * p
        df = k_f * (f_target - f)
        return np.array([dm, dp, df])

    y = np.array([0.0, 0.0, 1.0])
    segments_m: list[np.ndarray] = []
    segments_p: list[np.ndarray] = []
    segments_f: list[np.ndarray] = []
    segments_t: list[np.ndarray] = []

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
        segments_m.append(sol.y[0])
        segments_p.append(sol.y[1])
        segments_f.append(sol.y[2])
        y = sol.y[:, -1]

    t_out = np.concatenate(segments_t) if segments_t else np.array([])
    m_out = np.concatenate(segments_m) if segments_m else np.array([])
    p_out = np.concatenate(segments_p) if segments_p else np.array([])
    f_out = np.concatenate(segments_f) if segments_f else np.array([])
    return t_out, p_out, f_out


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

    return {
        "t_d": (t_eval / 86400.0).tolist(),
        "f_tx": f_tx.tolist(),
        "f_cochlea": f_cochlea.tolist(),
    }


def sweep_schedules() -> list[Regimen]:
    intervals = [14.0, 21.0, 28.0, 42.0, 56.0, 84.0]
    doses = [200.0, 500.0, 1000.0, 2000.0, 5000.0]
    mods = {
        "unmod": MRNA_HL_UNMOD_H,
        "m1psi": MRNA_HL_M1PSI_H,
    }
    targets = {
        "untargeted": LNP_UNTARGETED,
        "cochlear_tropic": LNP_COCHLEAR_TROPIC,
        "OHC_targeted": LNP_OHC_TARGETED,
    }

    regimens: list[Regimen] = []
    for mod_name, mod_hl in mods.items():
        for tgt_name, tgt_eff in targets.items():
            for T in intervals:
                for D in doses:
                    regimens.append(
                        Regimen(
                            modification=mod_name,
                            mrna_hl_h=mod_hl,
                            targeting=tgt_name,
                            eff_ohc=tgt_eff,
                            interval_d=T,
                            dose_intra=D,
                        )
                    )
    return regimens


def minimum_therapeutic_per_OHC(regs: list[Regimen]) -> dict:
    """For each (modification, interval) -- targeting does not affect per-OHC
    PK/PD -- find the smallest dose where transfected-OHC trough >= 2.0x."""
    out: dict = {}
    for mod in ["unmod", "m1psi"]:
        key = mod
        out[key] = {}
        for T in sorted({r.interval_d for r in regs}):
            viable = [
                r
                for r in regs
                if r.modification == mod
                and r.interval_d == T
                and r.targeting == "untargeted"
                and r.therapeutic_per_OHC
            ]
            if not viable:
                out[key][f"T_{int(T)}d"] = None
                continue
            r_min = min(viable, key=lambda x: x.dose_intra)
            out[key][f"T_{int(T)}d"] = {
                "min_dose_intra": r_min.dose_intra,
                "annual_intra_total_per_OHC": r_min.annual_intra_total,
                "annual_extra_total_per_OHC": r_min.annual_extra_total,
                "tx_trough_fold": r_min.tx_trough_fold,
                "tx_peak_fold": r_min.tx_peak_fold,
                "tx_tat": r_min.tx_time_above_threshold_frac,
                "annual_doses": r_min.annual_doses,
            }
    return out


def best_per_OHC(regs: list[Regimen]) -> dict:
    therapeutic = [r for r in regs if r.therapeutic_per_OHC and r.targeting == "untargeted"]
    if not therapeutic:
        return {"note": "no therapeutic per-OHC regimen found"}
    best = min(therapeutic, key=lambda r: r.annual_extra_total)
    return {
        "modification": best.modification,
        "interval_d": best.interval_d,
        "dose_intra_mol_per_OHC": best.dose_intra,
        "annual_doses": best.annual_doses,
        "annual_intra_total_per_OHC": best.annual_intra_total,
        "annual_extra_total_per_OHC": best.annual_extra_total,
        "tx_trough_fold": best.tx_trough_fold,
        "tx_peak_fold": best.tx_peak_fold,
        "tx_tat": best.tx_time_above_threshold_frac,
    }


def coverage_summary(regs: list[Regimen]) -> dict:
    """Cochlea-mean fold as a function of LNP targeting, assuming per-OHC
    regimen is already therapeutic. Separates delivery from PK/PD."""
    out: dict = {}
    for tgt in ["untargeted", "cochlear_tropic", "OHC_targeted"]:
        subset = [r for r in regs if r.targeting == tgt and r.therapeutic_per_OHC]
        if not subset:
            out[tgt] = None
            continue
        max_mean = max(r.cochlea_mean_fold for r in subset)
        max_trough = max(r.cochlea_trough_fold for r in subset)
        eff = {
            "untargeted": LNP_UNTARGETED,
            "cochlear_tropic": LNP_COCHLEAR_TROPIC,
            "OHC_targeted": LNP_OHC_TARGETED,
        }[tgt]
        out[tgt] = {
            "lnp_eff_ohc": eff,
            "max_cochlea_mean_fold": max_mean,
            "max_cochlea_trough_fold": max_trough,
            "theoretical_ceiling_cochlea_mean": eff * MAX_BOOST + (1.0 - eff),
        }
    return out


def main() -> None:
    regs = sweep_schedules()
    print(f"simulating {len(regs)} regimens", flush=True)

    exemplar_traces: dict = {}
    for i, reg in enumerate(regs):
        trace = simulate_regimen(reg)
        if reg.targeting == "cochlear_tropic" and reg.dose_intra == 1000.0 and reg.modification == "m1psi":
            exemplar_traces[f"T_{int(reg.interval_d)}d"] = {
                "t_d": trace["t_d"][::20],
                "f_tx": trace["f_tx"][::20],
                "f_cochlea": trace["f_cochlea"][::20],
            }
        if i % 30 == 0:
            print(
                f"  [{i+1}/{len(regs)}] {reg.modification:6s} {reg.targeting:17s} "
                f"T={reg.interval_d:4.0f}d D={reg.dose_intra:6.0f} "
                f"tx_trough={reg.tx_trough_fold:.2f} ther_per_OHC={reg.therapeutic_per_OHC}",
                flush=True,
            )

    min_ther = minimum_therapeutic_per_OHC(regs)
    best = best_per_OHC(regs)
    coverage = coverage_summary(regs)

    out = {
        "model": "STRC mRNA-LNP PK/PD integration (3-compartment, multi-dose)",
        "params": {
            "hill_K": HILL_K,
            "hill_n": HILL_N,
            "max_boost": MAX_BOOST,
            "threshold_fold": THRESHOLD,
            "mrna_hl_unmod_OHC_h": MRNA_HL_UNMOD_H,
            "mrna_hl_m1psi_OHC_h": MRNA_HL_M1PSI_H,
            "rbm24_hl_d": RBM24_HL_D,
            "strc_hl_d": STRC_HL_D,
            "lnp_endo_escape": LNP_ENDO_ESCAPE,
            "targeting": {
                "untargeted": LNP_UNTARGETED,
                "cochlear_tropic": LNP_COCHLEAR_TROPIC,
                "OHC_targeted": LNP_OHC_TARGETED,
            },
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
                "tx_trough_fold": r.tx_trough_fold,
                "tx_peak_fold": r.tx_peak_fold,
                "tx_mean_fold": r.tx_mean_fold,
                "tx_tat": r.tx_time_above_threshold_frac,
                "cochlea_trough_fold": r.cochlea_trough_fold,
                "cochlea_peak_fold": r.cochlea_peak_fold,
                "cochlea_mean_fold": r.cochlea_mean_fold,
                "therapeutic_per_OHC": r.therapeutic_per_OHC,
                "therapeutic_cochlea": r.therapeutic_cochlea,
            }
            for r in regs
        ],
        "minimum_therapeutic_per_OHC": min_ther,
        "best_per_OHC": best,
        "coverage_summary": coverage,
        "exemplar_traces_m1psi_cochlear_tropic_D1000": exemplar_traces,
    }

    out_path = ROOT / "mrna_lnp_pkpd_integration.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {out_path}")
    print(f"best per-OHC: {json.dumps(best, indent=2)}")
    print(f"coverage:     {json.dumps(coverage, indent=2)}")


if __name__ == "__main__":
    main()
