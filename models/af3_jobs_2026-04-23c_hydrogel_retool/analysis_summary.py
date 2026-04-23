#!/usr/bin/env python3
"""Parse AF3 2026-04-23c Hydrogel Phase 3 tail-retool results."""

from __future__ import annotations

import json
import glob
from pathlib import Path

BATCH_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23c_hydrogel_retool")
RESULTS = BATCH_DIR / "results"

JOBS = [
    ("fold_hydrogel_tail_1660_1710_x_tmem145", "51 aa minimal; clusters 4+5 with 10 aa flanks"),
    ("fold_hydrogel_tail_1640_1710_x_tmem145", "71 aa +loop cap; clusters 3+4+5"),
    ("fold_hydrogel_tail_1620_1710_x_tmem145", "91 aa wide; clusters 2+3+4+5 (SPPS ceiling)"),
]

GATE_IPTM = 0.50

# Phase 2 baselines (wrong-epitope, invalidated by design bug discovery)
PHASE2_BASELINES = {
    "rada16_wh2_native (aa 1454-1465 wrong epitope)": 0.39,
    "rada16_wh2_cmut": 0.37,
    "eak16_wh2_denovo": 0.37,
}
# Correct reference
ULTRAMINI_SOLO_X_TMEM145 = 0.43


def extract_job(name: str) -> dict:
    summaries = sorted(glob.glob(str(RESULTS / name / f"{name}_summary_confidences_*.json")))
    if not summaries:
        return {"status": "MISSING"}
    models = []
    for s in summaries:
        d = json.load(open(s))
        m = int(s.rsplit("_", 1)[1].split(".")[0])
        models.append({
            "model": m,
            "pTM": d.get("ptm"),
            "ipTM": d.get("iptm"),
            "ranking_score": d.get("ranking_score"),
            "fraction_disordered": d.get("fraction_disordered"),
            "has_clash": d.get("has_clash"),
            "chain_pair_iptm": d.get("chain_pair_iptm"),
            "chain_pair_pae_min": d.get("chain_pair_pae_min"),
        })
    models.sort(key=lambda m: -(m["ranking_score"] or 0))
    best = models[0]
    return {
        "status": "DONE",
        "n_models": len(models),
        "best_model": best["model"],
        "best_pTM": best["pTM"],
        "best_ipTM": best["ipTM"],
        "best_ranking_score": best["ranking_score"],
        "best_chain_pair_iptm": best["chain_pair_iptm"],
        "best_chain_pair_pae_min": best["chain_pair_pae_min"],
        "best_fraction_disordered": best["fraction_disordered"],
        "best_has_clash": best["has_clash"],
        "gate_verdict": "PASS" if (best["ipTM"] or 0) >= GATE_IPTM else "FAIL",
        "gate_margin": round((best["ipTM"] or 0) - GATE_IPTM, 3),
        "all_models_iptm": [(m["model"], m["ipTM"]) for m in models],
    }


def main():
    all_results = {}
    for name, desc in JOBS:
        all_results[name] = {"description": desc, **extract_job(name)}

    print("=== Hydrogel Phase 3 Tail Retool AF3 Results (2026-04-23) ===\n")
    print(f"Gate: ipTM >= {GATE_IPTM}")
    print(f"Reference: Ultra-Mini solo × TMEM145 = {ULTRAMINI_SOLO_X_TMEM145}")
    print(f"Phase 2 baselines (wrong epitope, invalid): {PHASE2_BASELINES}\n")

    summary_rows = []
    for name, desc in JOBS:
        r = all_results[name]
        if r.get("status") != "DONE":
            continue
        tail_aa = int(name.split("_")[3])  # 1660/1640/1620
        end_aa = int(name.split("_")[4])
        tail_len = end_aa - tail_aa + 1
        print(f"{name}")
        print(f"  {desc}")
        print(f"  best ipTM={r['best_ipTM']} ({r['gate_verdict']}, margin {r['gate_margin']:+.3f})")
        print(f"  best pTM={r['best_pTM']}, rank score={r['best_ranking_score']}")
        print(f"  chain_pair_iptm={r['best_chain_pair_iptm']}")
        print(f"  chain_pair_pae_min={r['best_chain_pair_pae_min']}")
        print(f"  all models: {r['all_models_iptm']}")
        print()
        summary_rows.append((tail_len, r['best_ipTM'], r['gate_verdict']))

    # Pattern detection: monotonic? ceiling?
    print("\n--- Progression pattern (tail length → ipTM) ---")
    for tl, ipt, v in sorted(summary_rows):
        print(f"  {tl} aa → ipTM {ipt} ({v})")

    # Hypothesis-level verdict
    any_pass = any(r.get("gate_verdict") == "PASS" for r in all_results.values())
    best_iptm = max((r.get("best_ipTM") or 0) for r in all_results.values())
    if any_pass:
        hyp = f"PASS — best tail crosses gate ({best_iptm}); advance to Phase 3b fuse winning tail to RADA16-WH2 scaffold"
    elif best_iptm >= ULTRAMINI_SOLO_X_TMEM145:
        hyp = f"MATCHES Ultra-Mini baseline ({best_iptm} vs {ULTRAMINI_SOLO_X_TMEM145}); Option A longer tail hits the isolated-peptide ceiling. Move to Option B (helical mimetic/miniprotein)"
    elif best_iptm > max(PHASE2_BASELINES.values()):
        hyp = f"IMPROVED over Phase 2 wrong-epitope baseline (best {best_iptm} vs Phase 2 0.39) but doesn't cross gate 0.50. Epitope correctness matters; longer extensions OR helical scaffolds next"
    else:
        hyp = f"No improvement. Floor effect suspected — AF3 assigns ~0.35-0.40 to any short peptide × TMEM145 regardless of epitope. Move to full miniprotein scaffold"

    summary = {
        "batch": "af3_jobs_2026-04-23c_hydrogel_retool",
        "description": "Hydrogel Phase 3 tail retool — correct GOLD-zone epitopes (fixes Phase 1 design bug)",
        "date_run": "2026-04-23",
        "gate_ipTM": GATE_IPTM,
        "reference_baseline_ultramini_solo": ULTRAMINI_SOLO_X_TMEM145,
        "phase2_baselines_wrong_epitope_invalid": PHASE2_BASELINES,
        "jobs": all_results,
        "hypothesis_verdict": hyp,
    }
    (BATCH_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n>>> Hypothesis-level verdict:\n{hyp}\n")


if __name__ == "__main__":
    main()
