#!/usr/bin/env python3
"""Parse AF3 2026-04-22 batch results; emit consolidated summary JSON."""

from __future__ import annotations

import json
import glob
from pathlib import Path

BATCH_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-22")
RESULTS = BATCH_DIR / "results"

JOBS = {
    "fold_strc_spy_reassembled_fold": {
        "hypothesis": "#10 STRC In Situ SpyCatcher Assembly",
        "gate": {"metric": "pTM", "threshold": 0.60, "target": "overall fold"},
    },
    "fold_strc_spy_reassembled_x_tmem145": {
        "hypothesis": "#10 STRC In Situ SpyCatcher Assembly (TMEM145 binding)",
        "gate": {"metric": "ipTM", "threshold": 0.40, "target": "STRC-TMEM145 interface"},
    },
    "fold_strc_tecta_chimera_fold": {
        "hypothesis": "#11 STRC Engineered TECTA Chimera",
        "gate": {"metric": "pTM", "threshold": 0.55, "target": "chimera fold"},
    },
    "fold_strc_tecta_chimera_x_tmem145": {
        "hypothesis": "#11 STRC Engineered TECTA Chimera (TMEM145 binding)",
        "gate": {"metric": "ipTM", "threshold": 0.40, "target": "chimera-TMEM145 interface"},
    },
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
        "gate_metric": gate["metric"],
        "gate_threshold": gate["threshold"],
        "gate_best_value": best_val,
        "gate_margin": (best_val - gate["threshold"]) if best_val is not None else None,
        "gate_verdict": "PASS" if (best_val or 0) >= gate["threshold"] else "FAIL",
        "all_models": models,
    }


def main():
    all_results = {}
    for job_name in JOBS:
        all_results[job_name] = extract_job(job_name)

    # Hypothesis-level verdicts
    def hypothesis_verdict(fold_job, binding_job):
        f = all_results[fold_job]
        b = all_results[binding_job]
        f_pass = f["gate_verdict"] == "PASS"
        b_pass = b["gate_verdict"] == "PASS"
        if f_pass and b_pass:
            return "BOTH_PASS"
        if f_pass and not b_pass:
            return "FOLD_OK_BINDING_FAIL"
        if not f_pass and b_pass:
            return "FOLD_FAIL_BINDING_OK"
        return "BOTH_FAIL"

    summary = {
        "batch": "af3_jobs_2026-04-22",
        "date_run": "2026-04-23",
        "n_jobs": len(JOBS),
        "jobs": all_results,
        "hypothesis_level": {
            "#10 SpyCatcher Assembly": {
                "fold_verdict": all_results["fold_strc_spy_reassembled_fold"]["gate_verdict"],
                "binding_verdict": all_results["fold_strc_spy_reassembled_x_tmem145"]["gate_verdict"],
                "fold_best_pTM": all_results["fold_strc_spy_reassembled_fold"]["best_pTM"],
                "binding_best_ipTM": all_results["fold_strc_spy_reassembled_x_tmem145"]["best_ipTM"],
                "overall": hypothesis_verdict("fold_strc_spy_reassembled_fold", "fold_strc_spy_reassembled_x_tmem145"),
            },
            "#11 TECTA Chimera": {
                "fold_verdict": all_results["fold_strc_tecta_chimera_fold"]["gate_verdict"],
                "binding_verdict": all_results["fold_strc_tecta_chimera_x_tmem145"]["gate_verdict"],
                "fold_best_pTM": all_results["fold_strc_tecta_chimera_fold"]["best_pTM"],
                "binding_best_ipTM": all_results["fold_strc_tecta_chimera_x_tmem145"]["best_ipTM"],
                "overall": hypothesis_verdict("fold_strc_tecta_chimera_fold", "fold_strc_tecta_chimera_x_tmem145"),
            },
        },
    }

    out = BATCH_DIR / "analysis_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {out}")

    print("\n=== AF3 2026-04-22 Batch Results ===")
    for hyp, h in summary["hypothesis_level"].items():
        print(f"\n{hyp}")
        print(f"  fold:    pTM {h['fold_best_pTM']} vs gate {all_results['fold_strc_spy_reassembled_fold' if 'SpyCatcher' in hyp else 'fold_strc_tecta_chimera_fold']['gate_threshold']} -> {h['fold_verdict']}")
        print(f"  binding: ipTM {h['binding_best_ipTM']} vs gate 0.40 -> {h['binding_verdict']}")
        print(f"  overall: {h['overall']}")


if __name__ == "__main__":
    main()
