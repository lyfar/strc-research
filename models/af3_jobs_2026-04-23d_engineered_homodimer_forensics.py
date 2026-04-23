#!/usr/bin/env python3
"""h26 Engineered Homodimer Phase 1 forensics.

Parses AF3 summary_confidences_*.json for each mutant run (± restart seed).
Chain layout (job_request.json): A,B = Ultra-Mini mutant (count=2); C,D = TMEM145 (count=2).

Gates (from MANIFEST):
  - homodimer_iptm  = chain_pair_iptm[A][B]                 >= 0.50 (up from WT 0.28-0.30)
  - tmem145_pair_1  = chain_pair_iptm[A][C]                 >= 0.40 (monomer baseline 0.43)
  - tmem145_pair_2  = chain_pair_iptm[B][D]                 >= 0.40

Aggregates: report max across models per run; also report mean ± sd.
PASS per run = best-model satisfies both gates.
"""
from __future__ import annotations

import json
import statistics
from pathlib import Path

UNPACK = Path("/tmp/h26_unpack")
OUT = Path(__file__).with_suffix("")  # dir for json output
OUT.mkdir(exist_ok=True)
SUMMARY_JSON = OUT.parent / "af3_jobs_2026-04-23d_engineered_homodimer_analysis.json"

JOBS = {
    "R1581F":           "fold_engineered_homodimer_r1581f",
    "R1581Y":           "fold_engineered_homodimer_r1581y",
    "S1579F":           "fold_engineered_homodimer_s1579f",
    "S1579F_restart":   "fold_engineered_homodimer_s1579f_restart",
    "S1579W":           "fold_engineered_homodimer_s1579w",
    "S1579W_restart":   "fold_engineered_homodimer_s1579w_restart",
}

# Pairings to extract
PAIRS = {
    "homodimer_AB": (0, 1),
    "umini_tmem_AC": (0, 2),
    "umini_tmem_BD": (1, 3),
    "umini_tmem_AD_cross": (0, 3),
    "umini_tmem_BC_cross": (1, 2),
    "tmem_tmem_CD": (2, 3),
}

GATE_HOMODIMER = 0.50
GATE_TMEM145 = 0.40

WT_BASELINE_HOMODIMER = (0.28, 0.30)
WT_MONOMER_TMEM145 = 0.43


def parse_run(run_dir: Path) -> dict:
    prefix = run_dir.name  # fold_engineered_homodimer_*
    models = []
    for i in range(5):
        sc = run_dir / f"{prefix}_summary_confidences_{i}.json"
        if not sc.exists():
            continue
        data = json.loads(sc.read_text())
        m = data["chain_pair_iptm"]
        rec = {"model": i, "iptm_global": data["iptm"], "ptm": data["ptm"], "ranking_score": data["ranking_score"]}
        for name, (a, b) in PAIRS.items():
            rec[name] = round((m[a][b] + m[b][a]) / 2, 4)  # symmetrize
        models.append(rec)
    return {"n_models": len(models), "models": models}


def aggregate(models: list[dict], key: str) -> dict:
    vals = [m[key] for m in models]
    return {"max": max(vals), "mean": round(statistics.mean(vals), 4),
            "sd": round(statistics.stdev(vals), 4) if len(vals) > 1 else 0.0,
            "values": vals}


def gate_pass(models: list[dict]) -> dict:
    """A mutant PASSES if any single model satisfies all 3 gates simultaneously,
    AND the mean homodimer ipTM across models is also >= gate (not a one-off fluke)."""
    per_model = []
    for m in models:
        p = (m["homodimer_AB"] >= GATE_HOMODIMER and
             m["umini_tmem_AC"] >= GATE_TMEM145 and
             m["umini_tmem_BD"] >= GATE_TMEM145)
        per_model.append(p)
    any_model = any(per_model)
    mean_homo = statistics.mean(m["homodimer_AB"] for m in models)
    mean_tmem_pair1 = statistics.mean(m["umini_tmem_AC"] for m in models)
    mean_tmem_pair2 = statistics.mean(m["umini_tmem_BD"] for m in models)
    mean_pass = (mean_homo >= GATE_HOMODIMER and mean_tmem_pair1 >= GATE_TMEM145 and mean_tmem_pair2 >= GATE_TMEM145)
    return {
        "any_model_passes_all_gates": any_model,
        "mean_passes_all_gates": mean_pass,
        "PASS_strict": any_model and mean_pass,
        "PASS_lenient_any_model": any_model,
        "per_model_pass_flags": per_model,
    }


def main():
    report = {"gates": {"homodimer_AB": GATE_HOMODIMER, "tmem145_per_subunit": GATE_TMEM145},
              "baselines": {"WT_homodimer_iptm": WT_BASELINE_HOMODIMER,
                            "WT_monomer_x_TMEM145_iptm": WT_MONOMER_TMEM145},
              "runs": {}}

    for label, run_name in JOBS.items():
        run_dir = UNPACK / run_name
        if not run_dir.exists():
            report["runs"][label] = {"error": f"dir missing: {run_dir}"}
            continue
        parsed = parse_run(run_dir)
        if parsed["n_models"] == 0:
            report["runs"][label] = {"error": "no summary_confidences files"}
            continue
        agg = {k: aggregate(parsed["models"], k) for k in PAIRS.keys()}
        report["runs"][label] = {
            "n_models": parsed["n_models"],
            "per_model": parsed["models"],
            "aggregate": agg,
            "gate_decision": gate_pass(parsed["models"]),
        }

    # Cross-run combined view per mutant (pool original + restart)
    mutants = ["R1581F", "R1581Y", "S1579F", "S1579W"]
    combined = {}
    for mut in mutants:
        pool = []
        for label in (mut, f"{mut}_restart"):
            r = report["runs"].get(label, {})
            if "per_model" in r:
                pool.extend(r["per_model"])
        if not pool:
            continue
        agg = {k: aggregate(pool, k) for k in PAIRS.keys()}
        combined[mut] = {
            "n_models_combined": len(pool),
            "aggregate": agg,
            "gate_decision": gate_pass(pool),
        }
    report["combined_by_mutant"] = combined

    SUMMARY_JSON.write_text(json.dumps(report, indent=2))

    # Human-readable stdout
    print("=" * 76)
    print("h26 Phase 1 — Engineered Homodimer AF3 forensics")
    print(f"Gates: homodimer A<->B >= {GATE_HOMODIMER};  Ultra-Mini x TMEM145 pair >= {GATE_TMEM145}")
    print(f"WT baseline homodimer: {WT_BASELINE_HOMODIMER};  WT x TMEM145 monomer: {WT_MONOMER_TMEM145}")
    print("=" * 76)
    for label, r in report["runs"].items():
        if "error" in r:
            print(f"\n[{label}] ERROR: {r['error']}")
            continue
        agg = r["aggregate"]
        g = r["gate_decision"]
        print(f"\n[{label}]  n={r['n_models']}")
        print(f"  homodimer_AB      max={agg['homodimer_AB']['max']:.3f}  mean={agg['homodimer_AB']['mean']:.3f} sd={agg['homodimer_AB']['sd']:.3f}")
        print(f"  Ultra-Mini x TMEM A<->C max={agg['umini_tmem_AC']['max']:.3f}  mean={agg['umini_tmem_AC']['mean']:.3f} sd={agg['umini_tmem_AC']['sd']:.3f}")
        print(f"  Ultra-Mini x TMEM B<->D max={agg['umini_tmem_BD']['max']:.3f}  mean={agg['umini_tmem_BD']['mean']:.3f} sd={agg['umini_tmem_BD']['sd']:.3f}")
        print(f"  TMEM x TMEM C<->D max={agg['tmem_tmem_CD']['max']:.3f}  mean={agg['tmem_tmem_CD']['mean']:.3f}")
        print(f"  gate PASS_strict={g['PASS_strict']}  any_model={g['any_model_passes_all_gates']}  mean={g['mean_passes_all_gates']}")

    print("\n" + "=" * 76)
    print("COMBINED by mutant (pool original + restart, 10 models each)")
    print("=" * 76)
    for mut, r in combined.items():
        agg = r["aggregate"]
        g = r["gate_decision"]
        print(f"\n[{mut} combined]  n={r['n_models_combined']}")
        print(f"  homodimer_AB      max={agg['homodimer_AB']['max']:.3f}  mean={agg['homodimer_AB']['mean']:.3f} sd={agg['homodimer_AB']['sd']:.3f}")
        print(f"  Ultra-Mini x TMEM A<->C max={agg['umini_tmem_AC']['max']:.3f}  mean={agg['umini_tmem_AC']['mean']:.3f}")
        print(f"  Ultra-Mini x TMEM B<->D max={agg['umini_tmem_BD']['max']:.3f}  mean={agg['umini_tmem_BD']['mean']:.3f}")
        print(f"  gate PASS_strict={g['PASS_strict']}  any_model={g['any_model_passes_all_gates']}  mean={g['mean_passes_all_gates']}")

    # Verdict
    print("\n" + "=" * 76)
    print("VERDICT")
    print("=" * 76)
    winners = [mut for mut, r in combined.items() if r["gate_decision"]["PASS_strict"]]
    lenient = [mut for mut, r in combined.items() if r["gate_decision"]["PASS_lenient_any_model"] and mut not in winners]
    if winners:
        print(f"  STRICT PASS: {winners}  -> h26 -> A-tier confirmed, queue Phase 2 SEC-AUC")
    if lenient:
        print(f"  LENIENT (any model only, mean fails): {lenient}  -> weak signal, consider double-mutant")
    if not winners and not lenient:
        print("  ALL FAIL: single-point mutation insufficient. Next: double-mutant (S1579W + L1581Y) or de novo dimer scaffold.")
    print(f"\nFull report written to: {SUMMARY_JSON}")


if __name__ == "__main__":
    main()
