#!/usr/bin/env python3
"""Parse AF3 2026-04-23b Hydrogel Phase 2 results; emit consolidated summary JSON.

Expects 6 job directories; gracefully skips any missing. Phase 2 gate: ipTM >= 0.50
on BOTH interfaces (TMEM145 + G-actin trimer) per candidate peptide.
"""

from __future__ import annotations

import json
import glob
from pathlib import Path

BATCH_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23b_hydrogel")
RESULTS = BATCH_DIR / "results"

CANDIDATES = ["rada16_wh2_native", "rada16_wh2_cmut", "eak16_wh2_denovo"]
PARTNERS = ["tmem145", "actin"]

GATE_IPTM = 0.50


def extract_job(job_name: str) -> dict:
    # Check canonical name first, then re-run variants (AF3 sometimes requires
    # name change on resubmit). For rada16_wh2_native_x_actin the first submit
    # hit a transient server issue; re-submitted as *_2nndrun_1st_loading.
    candidates = [job_name, f"{job_name}_2nndrun_1st_loading"]
    job_dir = None
    resolved_name = job_name
    for cand in candidates:
        p = RESULTS / cand
        if p.exists():
            job_dir = p
            resolved_name = cand
            break
    if job_dir is None:
        return {"status": "MISSING", "job": job_name}
    summaries = sorted(glob.glob(str(job_dir / f"{resolved_name}_summary_confidences_*.json")))
    if not summaries:
        return {"status": "EMPTY", "job": job_name}
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
    return {
        "status": "DONE",
        "job": job_name,
        "n_models": len(models),
        "best_model_ranked": best["model"],
        "best_pTM": best["pTM"],
        "best_ipTM": best["ipTM"],
        "best_ranking_score": best["ranking_score"],
        "best_chain_pair_iptm": best["chain_pair_iptm"],
        "best_chain_pair_pae_min": best["chain_pair_pae_min"],
        "best_fraction_disordered": best["fraction_disordered"],
        "best_has_clash": best["has_clash"],
        "gate_verdict": (
            "PASS" if (best["ipTM"] or 0) >= GATE_IPTM
            else "FAIL"
        ),
        "gate_margin": round((best["ipTM"] or 0) - GATE_IPTM, 3),
        "all_models_iptm": [(m["model"], m["ipTM"]) for m in models],
    }


def main():
    all_results = {}
    for cand in CANDIDATES:
        for partner in PARTNERS:
            job_name = f"fold_hydrogel_{cand}_x_{partner}"
            all_results[job_name] = extract_job(job_name)

    # Per-candidate verdicts: BOTH interfaces must pass (gate ipTM >= 0.50)
    candidate_verdicts = {}
    for cand in CANDIDATES:
        tmem = all_results[f"fold_hydrogel_{cand}_x_tmem145"]
        actin = all_results[f"fold_hydrogel_{cand}_x_actin"]
        t_status = tmem.get("gate_verdict") or tmem.get("status")
        a_status = actin.get("gate_verdict") or actin.get("status")
        both_done = tmem.get("status") == "DONE" and actin.get("status") == "DONE"
        if both_done:
            if t_status == "PASS" and a_status == "PASS":
                overall = "BOTH_PASS"
            elif t_status == "PASS" and a_status == "FAIL":
                overall = "TMEM145_PASS_ACTIN_FAIL"
            elif t_status == "FAIL" and a_status == "PASS":
                overall = "TMEM145_FAIL_ACTIN_PASS"
            else:
                overall = "BOTH_FAIL"
        else:
            missing = [k for k, d in [("tmem145", tmem), ("actin", actin)]
                       if d.get("status") != "DONE"]
            overall = f"INCOMPLETE (missing: {','.join(missing)})"
        candidate_verdicts[cand] = {
            "tmem145_ipTM": tmem.get("best_ipTM"),
            "tmem145_verdict": t_status,
            "actin_ipTM": actin.get("best_ipTM"),
            "actin_verdict": a_status,
            "overall": overall,
        }

    # Hypothesis-level worst-case (best-performing candidate across both interfaces)
    any_both_pass = any(v["overall"] == "BOTH_PASS" for v in candidate_verdicts.values())
    any_tmem_pass = any(v.get("tmem145_verdict") == "PASS" for v in candidate_verdicts.values())
    any_actin_pass = any(v.get("actin_verdict") == "PASS" for v in candidate_verdicts.values())

    if any_both_pass:
        hyp_verdict = "ADVANCE (B -> A; Phase 2b Martini3 + wet-lab)"
    elif any_tmem_pass and any_actin_pass:
        hyp_verdict = "SPLIT (different candidates pass each interface; retool for any single candidate)"
    elif any_tmem_pass and not any_actin_pass:
        hyp_verdict = "B stays (TMEM145 handle works; actin architecture broken, redesign N-terminus)"
    elif not any_tmem_pass and any_actin_pass:
        hyp_verdict = "B stays (actin handle works; TMEM145 tail falsified, retool tail)"
    else:
        hyp_verdict = "B -> C (dual-interface short-peptide architecture not structurally competent)"

    summary = {
        "batch": "af3_jobs_2026-04-23b_hydrogel",
        "description": "Hydrogel Phase 2 AF3-Multimer: top-3 SAP candidates x TMEM145 + x G-actin trimer",
        "date_run": "2026-04-23",
        "n_jobs_expected": 6,
        "n_jobs_present": sum(1 for r in all_results.values() if r.get("status") == "DONE"),
        "gate_ipTM": GATE_IPTM,
        "jobs": all_results,
        "candidate_verdicts": candidate_verdicts,
        "hypothesis_verdict": hyp_verdict,
    }

    out = BATCH_DIR / "analysis_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {out}")

    print("\n=== Hydrogel Phase 2 AF3 Results (2026-04-23) ===")
    print(f"\nGate: ipTM >= {GATE_IPTM} on BOTH interfaces per candidate")
    print(f"Jobs present: {summary['n_jobs_present']} / 6\n")
    for cand, v in candidate_verdicts.items():
        print(f"{cand}:")
        print(f"  x TMEM145: ipTM={v['tmem145_ipTM']} -> {v['tmem145_verdict']}")
        print(f"  x actin:   ipTM={v['actin_ipTM']} -> {v['actin_verdict']}")
        print(f"  => {v['overall']}")
        print()
    print(f"Hypothesis-level: {hyp_verdict}")

    # Detailed breakdown for each TMEM145 job (most important, 3/3 present)
    print(f"\n--- Per-candidate × TMEM145 model spread ---")
    for cand in CANDIDATES:
        r = all_results[f"fold_hydrogel_{cand}_x_tmem145"]
        if r.get("status") != "DONE":
            continue
        print(f"{cand} x TMEM145: all models ipTM = {r['all_models_iptm']}")
        print(f"  chain_pair_iptm = {r['best_chain_pair_iptm']}")
        print(f"  chain_pair_pae_min = {r['best_chain_pair_pae_min']}")


if __name__ == "__main__":
    main()
