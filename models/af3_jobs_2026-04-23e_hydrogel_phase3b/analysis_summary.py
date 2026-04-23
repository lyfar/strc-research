#!/usr/bin/env python3
"""Parse AF3 2026-04-23e Hydrogel Phase 3b full-construct results."""

from __future__ import annotations

import json
import glob
from pathlib import Path

BATCH_DIR = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23e_hydrogel_phase3b")
RESULTS = BATCH_DIR / "results"

JOBS = [
    ("fold_hydrogel_rada16_wh2_tail71_x_tmem145", "tail71", "tmem145",
     "114 aa peptide × TMEM145"),
    ("fold_hydrogel_rada16_wh2_tail91_x_tmem145", "tail91", "tmem145",
     "134 aa peptide × TMEM145"),
    ("fold_hydrogel_rada16_wh2_tail71_x_actin",   "tail71", "actin",
     "114 aa peptide × G-actin trimer"),
    ("fold_hydrogel_rada16_wh2_tail91_x_actin",   "tail91", "actin",
     "134 aa peptide × G-actin trimer"),
]

GATE_IPTM = 0.50

# Baselines
PHASE2_NATIVE_WRONG_EPITOPE_X_TMEM145 = 0.39
PHASE2_NATIVE_X_ACTIN = 0.58  # consistent 0.58-0.59 across 3 Phase 2 candidates
PHASE3_TAIL71_SOLO_X_TMEM145 = 0.63
PHASE3_TAIL91_SOLO_X_TMEM145 = 0.68


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
    for name, construct, partner, desc in JOBS:
        all_results[name] = {
            "construct": construct,
            "partner": partner,
            "description": desc,
            **extract_job(name),
        }

    print("=== Hydrogel Phase 3b Full-Construct AF3 Results (2026-04-23) ===\n")
    print(f"Gate: ipTM >= {GATE_IPTM} on BOTH interfaces per candidate")
    print(f"Baselines:")
    print(f"  Phase 2 native wrong-epitope × TMEM145 = {PHASE2_NATIVE_WRONG_EPITOPE_X_TMEM145}")
    print(f"  Phase 2 full construct × actin        = {PHASE2_NATIVE_X_ACTIN} (consistent across 3)")
    print(f"  Phase 3 tail71 solo × TMEM145         = {PHASE3_TAIL71_SOLO_X_TMEM145}")
    print(f"  Phase 3 tail91 solo × TMEM145         = {PHASE3_TAIL91_SOLO_X_TMEM145}\n")

    for name, construct, partner, desc in JOBS:
        r = all_results[name]
        if r.get("status") != "DONE":
            print(f"{name}  MISSING")
            continue
        print(f"{name}")
        print(f"  {desc}")
        print(f"  best ipTM={r['best_ipTM']} ({r['gate_verdict']}, margin {r['gate_margin']:+.3f})")
        print(f"  best pTM={r['best_pTM']}, rank score={r['best_ranking_score']}")
        print(f"  chain_pair_iptm={r['best_chain_pair_iptm']}")
        print(f"  chain_pair_pae_min={r['best_chain_pair_pae_min']}")
        print(f"  all models: {r['all_models_iptm']}")
        print()

    # Per-candidate dual-interface gate
    print("\n--- Per-candidate dual-interface verdict ---")
    candidate_verdicts = {}
    for construct in ("tail71", "tail91"):
        tmem_r = all_results[f"fold_hydrogel_rada16_wh2_{construct}_x_tmem145"]
        actin_r = all_results[f"fold_hydrogel_rada16_wh2_{construct}_x_actin"]
        tmem_iptm = tmem_r.get("best_ipTM") or 0
        actin_iptm = actin_r.get("best_ipTM") or 0
        tmem_pass = tmem_iptm >= GATE_IPTM
        actin_pass = actin_iptm >= GATE_IPTM
        both_pass = tmem_pass and actin_pass
        verdict = "PASS" if both_pass else "FAIL"
        candidate_verdicts[construct] = {
            "tmem145_ipTM": tmem_iptm,
            "actin_ipTM": actin_iptm,
            "tmem145_pass": tmem_pass,
            "actin_pass": actin_pass,
            "verdict": verdict,
        }
        solo_baseline = (PHASE3_TAIL71_SOLO_X_TMEM145 if construct == "tail71"
                         else PHASE3_TAIL91_SOLO_X_TMEM145)
        print(f"  {construct}: TMEM145 {tmem_iptm} ({'PASS' if tmem_pass else 'FAIL'}) "
              f"[solo baseline {solo_baseline}, delta {tmem_iptm - solo_baseline:+.3f}], "
              f"actin {actin_iptm} ({'PASS' if actin_pass else 'FAIL'}) "
              f"[Phase 2 baseline {PHASE2_NATIVE_X_ACTIN}, delta {actin_iptm - PHASE2_NATIVE_X_ACTIN:+.3f}]  "
              f"=> {verdict}")

    # Hypothesis-level verdict
    n_pass = sum(1 for v in candidate_verdicts.values() if v["verdict"] == "PASS")
    if n_pass >= 1:
        winners = [c for c, v in candidate_verdicts.items() if v["verdict"] == "PASS"]
        hyp = (f"PASS — {n_pass}/2 candidates ({', '.join(winners)}) cross BOTH interface gates. "
               f"Hypothesis #9 Hydrogel A-tier CONFIRMED at full-construct AF3 level. "
               f"Advance to Phase 2b Martini3 CG MD + Phase 2c wet-lab cell-based actin-bundling assay.")
    else:
        # Partial / fail analysis
        tmem_any_pass = any(v["tmem145_pass"] for v in candidate_verdicts.values())
        actin_any_pass = any(v["actin_pass"] for v in candidate_verdicts.values())
        if tmem_any_pass and not actin_any_pass:
            hyp = ("PARTIAL — TMEM145 preserved but scaffold fusion disturbs WH2 actin binding. "
                   "Try longer inter-domain linker or alternative WH2 placement.")
        elif actin_any_pass and not tmem_any_pass:
            hyp = ("PARTIAL — Actin preserved but scaffold fusion contaminates TMEM145 epitope. "
                   "Try flexible linker between RADA16 and tail, or C-terminal scaffold placement.")
        else:
            hyp = ("FAIL — Scaffold architecture doesn't support full dual-interface peptide. "
                   "Move to Option B (helical mimetic / miniprotein scaffold replacing RADA16).")

    summary = {
        "batch": "af3_jobs_2026-04-23e_hydrogel_phase3b",
        "description": "Hydrogel Phase 3b full-construct validation — winning GOLD-zone tails fused to WH2+RADA16 scaffold",
        "date_run": "2026-04-23",
        "gate_ipTM": GATE_IPTM,
        "baselines": {
            "phase_2_native_wrong_epitope_x_tmem145": PHASE2_NATIVE_WRONG_EPITOPE_X_TMEM145,
            "phase_2_full_construct_x_actin": PHASE2_NATIVE_X_ACTIN,
            "phase_3_tail71_solo_x_tmem145": PHASE3_TAIL71_SOLO_X_TMEM145,
            "phase_3_tail91_solo_x_tmem145": PHASE3_TAIL91_SOLO_X_TMEM145,
        },
        "jobs": all_results,
        "candidate_verdicts": candidate_verdicts,
        "hypothesis_verdict": hyp,
    }
    (BATCH_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n>>> Hypothesis-level verdict:\n{hyp}\n")


if __name__ == "__main__":
    main()
