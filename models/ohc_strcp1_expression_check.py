#!/usr/bin/env python3
"""
OHC STRCP1 expression check — pull GTEx v8 + Ensembl/Protein Atlas metadata
to assess whether the STRC pseudogene paralog STRCP1 is transcriptionally
silent or co-expressed with STRC.

Cross-cutting unblock: PE Phase 4 ([[STRC PE Phase4 STRCP1 Paralog Off-Target]])
and ASO Phase 2 ([[STRC ASO Phase2 STRCP1 Paralog Cross-Hybridization]]) both
killed all candidate pools on STRCP1 cross-hybridization. The shared escape
path is "STRCP1 is transcriptionally silent in OHC, so the off-target has no
RNA substrate." This script tests the *necessary precondition*: if STRCP1
is silent in *every* tissue GTEx covers, the OHC-silent assumption is
plausible. If STRCP1 is co-expressed with STRC in any STRC-expressing
tissue, the escape path likely fails.

GTEx limitation: GTEx has 54 tissues but no cochlea / inner ear. We use
the highest-STRC tissues (cerebellum, testis, skin) as proxies.

Output: JSON with per-tissue median TPM for STRCP1 and STRC + a derived
STRC:STRCP1 ratio in tissues where STRC > 0.1 TPM.

Run env: any conda env with `requests` (or just stdlib `urllib`).
No heavy compute, no GPU, ~3 seconds total runtime.
"""

from __future__ import annotations

import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

WORK_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
OUT_JSON = WORK_DIR / "ohc_strcp1_expression_check.json"

# GENCODE v26 IDs as used by GTEx v8:
GENCODE_STRC = "ENSG00000242866.9"
GENCODE_STRCP1 = "ENSG00000166763.7"

GTEX_URL_TEMPLATE = (
    "https://gtexportal.org/api/v2/expression/medianGeneExpression"
    "?gencodeId={gid}&datasetId=gtex_v8"
)


def log(msg):
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}", flush=True)


def fetch_gtex(gencode_id: str) -> list[dict]:
    url = GTEX_URL_TEMPLATE.format(gid=gencode_id)
    log(f"  GET {url}")
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read())
    return data.get("data", [])


def main():
    log("OHC STRCP1 expression check — querying GTEx v8")

    strc_data = fetch_gtex(GENCODE_STRC)
    strcp1_data = fetch_gtex(GENCODE_STRCP1)
    log(f"  STRC tissues: {len(strc_data)}")
    log(f"  STRCP1 tissues: {len(strcp1_data)}")

    strc_by_tissue = {row["tissueSiteDetailId"]: row["median"] for row in strc_data}
    strcp1_by_tissue = {row["tissueSiteDetailId"]: row["median"] for row in strcp1_data}

    tissues = sorted(set(strc_by_tissue) | set(strcp1_by_tissue))
    rows = []
    for t in tissues:
        s = strc_by_tissue.get(t, 0.0)
        p = strcp1_by_tissue.get(t, 0.0)
        ratio = (s / p) if p > 0 else float("inf")
        rows.append({
            "tissue": t,
            "strc_tpm": s,
            "strcp1_tpm": p,
            "strc_to_strcp1_ratio": ratio if ratio != float("inf") else None,
        })
    rows.sort(key=lambda r: r["strc_tpm"], reverse=True)

    # Summary stats
    strc_expressed = [r for r in rows if r["strc_tpm"] > 0.1]
    strcp1_expressed = [r for r in rows if r["strcp1_tpm"] > 0.1]
    co_expressed = [r for r in rows if r["strc_tpm"] > 0.1 and r["strcp1_tpm"] > 0.1]

    summary = {
        "meta": {
            "date": datetime.now(timezone.utc).date().isoformat(),
            "source": "GTEx v8 medianGeneExpression API",
            "strc_gencode": GENCODE_STRC,
            "strcp1_gencode": GENCODE_STRCP1,
            "limitation": "GTEx has no cochlea/inner ear; highest-STRC tissues used as proxies",
        },
        "totals": {
            "n_tissues": len(tissues),
            "n_strc_expressed_tpm_gt_0.1": len(strc_expressed),
            "n_strcp1_expressed_tpm_gt_0.1": len(strcp1_expressed),
            "n_co_expressed_tpm_gt_0.1": len(co_expressed),
        },
        "verdict": _derive_verdict(rows),
        "rows_top10_by_strc": rows[:10],
        "rows_full": rows,
    }

    OUT_JSON.write_text(json.dumps(summary, indent=2))
    log(f"wrote {OUT_JSON}")
    log("")
    log(f"verdict: {summary['verdict']['label']}")
    log(f"  reason: {summary['verdict']['reason']}")
    log("")
    log("Top 10 STRC-expressing tissues (proxy for OHC where direct data missing):")
    log(f"  {'tissue':42s} {'STRC TPM':>10s} {'STRCP1 TPM':>12s} {'ratio':>8s}")
    for r in rows[:10]:
        ratio = f"{r['strc_to_strcp1_ratio']:8.2f}" if r['strc_to_strcp1_ratio'] is not None else "    inf"
        log(f"  {r['tissue']:42s} {r['strc_tpm']:10.4f} {r['strcp1_tpm']:12.4f} {ratio}")


def _derive_verdict(rows: list[dict]) -> dict:
    """Bin verdict on whether STRCP1 silent across STRC-expressing tissues."""
    # Look at top 5 STRC-expressing tissues
    top5 = rows[:5]
    ratios = [r["strc_to_strcp1_ratio"] for r in top5
              if r["strcp1_tpm"] > 0.05 and r["strc_to_strcp1_ratio"] is not None]
    if not ratios:
        return {
            "label": "STRCP1_SILENT_PROBABLE",
            "reason": "Top-5 STRC-expressing tissues all show STRCP1 < 0.05 TPM — paralog likely silent in OHC too. Both PE and ASO escape path PLAUSIBLE.",
        }
    median_ratio = sorted(ratios)[len(ratios) // 2]
    if median_ratio >= 10:
        return {
            "label": "STRCP1_LOW_FAVORABLE",
            "reason": f"STRC:STRCP1 median ratio {median_ratio:.1f}:1 across top STRC-expressing tissues. ASO loses ~10% of dose to STRCP1; tolerable. PE off-target risk still present per-cell.",
        }
    if median_ratio >= 3:
        return {
            "label": "STRCP1_CO_EXPRESSED_PARTIAL",
            "reason": f"STRC:STRCP1 median ratio {median_ratio:.1f}:1. ASO loses 25-50% of dose to STRCP1 — significant but not fatal. PE off-target rate stays unchanged per copy.",
        }
    return {
        "label": "STRCP1_CO_EXPRESSED_KILL",
        "reason": f"STRC:STRCP1 median ratio {median_ratio:.1f}:1 (essentially co-expressed). ASO loses >50% of dose to paralog cleavage — kills splice-modulation ROI. PE off-target risk unchanged. Escape path likely UNAVAILABLE.",
    }


if __name__ == "__main__":
    sys.exit(main())
