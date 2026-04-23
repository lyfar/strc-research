#!/usr/bin/env python3
"""Parse AF3 2026-04-23 batch results (SpyCatcher Phase 1b, aa 700 split); emit consolidated summary JSON."""

from __future__ import annotations

import json
import glob
from pathlib import Path

BATCH_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23")
RESULTS = BATCH_DIR / "results"

JOBS = {
    "fold_strc_spy_reassembled_aa700_fold": {
        "hypothesis": "#10 STRC In Situ SpyCatcher Assembly (Phase 1b aa700 split)",
        "gate": {"metric": "pTM", "threshold": 0.60, "target": "overall fold"},
    },
    "fold_strc_spy_reassembled_aa700_x_tmem145": {
        "hypothesis": "#10 STRC In Situ SpyCatcher Assembly (Phase 1b TMEM145 binding)",
        "gate": {"metric": "ipTM", "threshold": 0.40, "target": "STRC-TMEM145 interface"},
    },
}

# Phase 1a baseline (aa 1074/1075 split) for direct delta comparison
PRIOR_BASELINE = {
    "fold": {"best_pTM": 0.59, "gate": 0.60, "verdict": "MISS", "margin": -0.01},
    "binding": {"best_ipTM": 0.37, "gate": 0.40, "verdict": "MISS", "margin": -0.03,
                "chain_pair_iptm": [0.51, 0.37], "pae_min": 9.5},
}


def extract_job(name: str) -> dict:
    summaries = sorted(glob.glob(str(RESULTS / name / f"{name}_summary_confidences_*.json")))
    if not summaries:
        return {"error": "no summaries"}
    models = []
    for s in summaries:
        d = json.load(open(s))
        m = int(s.rsplit("_", 1)[1].split(".")[0])
        entry = {
            "model": m,
            "pTM": d.get("ptm"),
            "ipTM": d.get("iptm"),
            "ranking_score": d.get("ranking_score"),
            "fraction_disordered": d.get("fraction_disordered"),
            "has_clash": d.get("has_clash"),
            "chain_pair_iptm": d.get("chain_pair_iptm"),
            "chain_pair_pae_min": d.get("chain_pair_pae_min"),
        }
        models.append(entry)
    models.sort(key=lambda m: -(m["ranking_score"] or 0))
    best = models[0]
    gate = JOBS[name]["gate"]
    best_val = best.get(gate["metric"])
    return {
        "hypothesis": JOBS[name]["hypothesis"],
        "n_models": len(models),
        "best_model_ranked": best["model"],
        "best_pTM": best["pTM"],
        "best_ipTM": best["ipTM"],
        "best_ranking_score": best["ranking_score"],
        "best_chain_pair_iptm": best.get("chain_pair_iptm"),
        "best_chain_pair_pae_min": best.get("chain_pair_pae_min"),
        "best_fraction_disordered": best.get("fraction_disordered"),
        "best_has_clash": best.get("has_clash"),
        "gate_metric": gate["metric"],
        "gate_threshold": gate["threshold"],
        "gate_best_value": best_val,
        "gate_margin": (best_val - gate["threshold"]) if best_val is not None else None,
        "gate_verdict": "PASS" if (best_val or 0) >= gate["threshold"] else "FAIL",
        "all_models": models,
    }


def main():
    all_results = {name: extract_job(name) for name in JOBS}

    f = all_results["fold_strc_spy_reassembled_aa700_fold"]
    b = all_results["fold_strc_spy_reassembled_aa700_x_tmem145"]
    f_pass = f["gate_verdict"] == "PASS"
    b_pass = b["gate_verdict"] == "PASS"
    if f_pass and b_pass:
        overall = "BOTH_PASS"
    elif f_pass and not b_pass:
        overall = "FOLD_OK_BINDING_FAIL"
    elif not f_pass and b_pass:
        overall = "FOLD_FAIL_BINDING_OK"
    else:
        overall = "BOTH_FAIL"

    # Deltas vs Phase 1a baseline (aa 1074/1075 split) — seed 42 in both batches,
    # so deltas are directly attributable to split-point choice
    fold_delta = f["best_pTM"] - PRIOR_BASELINE["fold"]["best_pTM"]
    binding_delta = b["best_ipTM"] - PRIOR_BASELINE["binding"]["best_ipTM"]

    summary = {
        "batch": "af3_jobs_2026-04-23",
        "description": "SpyCatcher Assembly Phase 1b — alternate split at aa 700 (Mini-STRC canonical LRR-to-ARM domain boundary)",
        "date_run": "2026-04-23",
        "n_jobs": len(JOBS),
        "jobs": all_results,
        "hypothesis_level": {
            "#10 SpyCatcher Assembly (Phase 1b)": {
                "fold_verdict": f["gate_verdict"],
                "binding_verdict": b["gate_verdict"],
                "fold_best_pTM": f["best_pTM"],
                "binding_best_ipTM": b["best_ipTM"],
                "overall": overall,
            },
        },
        "phase_1a_vs_1b_delta": {
            "fold_pTM": {
                "phase_1a_aa1074": PRIOR_BASELINE["fold"]["best_pTM"],
                "phase_1b_aa700": f["best_pTM"],
                "delta": round(fold_delta, 3),
            },
            "binding_ipTM": {
                "phase_1a_aa1074": PRIOR_BASELINE["binding"]["best_ipTM"],
                "phase_1b_aa700": b["best_ipTM"],
                "delta": round(binding_delta, 3),
            },
        },
    }

    out = BATCH_DIR / "analysis_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {out}")

    print("\n=== AF3 2026-04-23 Batch Results (SpyCatcher Phase 1b, aa700 split) ===")
    print(f"\nFold ({f['best_model_ranked']}):")
    print(f"  pTM {f['best_pTM']} vs gate {f['gate_threshold']} -> {f['gate_verdict']} (margin {f['gate_margin']:+.3f})")
    print(f"  disorder {f['best_fraction_disordered']}, clash {f['best_has_clash']}")
    print(f"\nBinding (model {b['best_model_ranked']}):")
    print(f"  ipTM {b['best_ipTM']} vs gate {b['gate_threshold']} -> {b['gate_verdict']} (margin {b['gate_margin']:+.3f})")
    print(f"  chain_pair_iptm {b['best_chain_pair_iptm']}")
    print(f"  chain_pair_pae_min {b['best_chain_pair_pae_min']}")
    print(f"\nOverall: {overall}")
    print(f"\nPhase 1a (aa 1074) vs 1b (aa 700) deltas:")
    print(f"  fold pTM: {PRIOR_BASELINE['fold']['best_pTM']} -> {f['best_pTM']} ({fold_delta:+.3f})")
    print(f"  binding ipTM: {PRIOR_BASELINE['binding']['best_ipTM']} -> {b['best_ipTM']} ({binding_delta:+.3f})")

    print(f"\nAll fold models:")
    for m in f["all_models"]:
        print(f"  model {m['model']}: pTM={m['pTM']}, rs={m['ranking_score']}")
    print(f"\nAll binding models:")
    for m in b["all_models"]:
        print(f"  model {m['model']}: ipTM={m['ipTM']}, pTM={m['pTM']}, rs={m['ranking_score']}, cp_iptm={m['chain_pair_iptm']}")


if __name__ == "__main__":
    main()
