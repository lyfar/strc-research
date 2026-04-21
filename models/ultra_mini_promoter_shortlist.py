"""Ultra-Mini STRC vector regulatory architecture — promoter/enhancer shortlist.

Question: given Ultra-Mini (aa 1075-1775, 2,103 bp CDS) leaves ~1,994 bp AAV
headroom vs ~869 bp for the clinical 700-1775 mini, what regulatory upgrade
beats the current B8 enhancer alone? Score candidates on size fit, OHC
specificity (vs IHC and off-target cells), expression strength, and AAV track
record for DFNB-adjacent applications.

Candidates:
  1. B8 alone (current design, ARBITER synthetic enhancer, 587 bp) - BASELINE
  2. B8 + WPRE3 (posttranscriptional element, +~600 bp)
  3. Myo15 truncated promoter (Zhao 2024, 956 bp variant)
  4. Myo15 longer (Zhao 2024, 1,157 bp variant)
  5. Myo15 full (Liu 2024, 1,611 bp native promoter)
  6. B8 + Myo15-core hybrid (~800 bp, belt-and-suspenders specificity)
  7. Prestin endogenous promoter (native Slc26a5 regulatory region, ~1,800 bp est)

Scoring rubric (1-5 per axis, tier = min):
  - size_fit     : 1=doesn't fit, 5=fits with >500 bp spare
  - ohc_specificity : 1=ubiquitous, 3=HC-wide, 5=OHC-exclusive
  - expression_strength : 1=weak, 5=strong endogenous-equivalent
  - aav_track_record : 1=untested, 5=in-vivo hearing-rescue published
"""
import json
from pathlib import Path

AAV_ITR_CAPACITY = 4700       # industry standard limit (ITR to ITR)
ITR_OVERHEAD = 300            # 2 x ~150 bp ITRs
KOZAK_SP_POLYA = 450          # Kozak (~20) + IgK SP (~60) + stop (3) + bGH polyA (~225) + misc spacer (~140)
ULTRA_MINI_CDS = 2106
MINI_700_1775_CDS = 3231

available_ultra = AAV_ITR_CAPACITY - ITR_OVERHEAD - KOZAK_SP_POLYA - ULTRA_MINI_CDS
available_mini  = AAV_ITR_CAPACITY - ITR_OVERHEAD - KOZAK_SP_POLYA - MINI_700_1775_CDS

print(f"Ultra-Mini regulatory budget: {available_ultra} bp")
print(f"Mini 700-1775 regulatory budget: {available_mini} bp")
print(f"Upgrade gain: +{available_ultra - available_mini} bp for Ultra-Mini")

# Size of regulatory modules (bp)
MODULES = {
    "B8_enhancer":           {"size": 587,  "source": "ARBITER synthetic panel; [[STRC B8 Enhancer Selection]]"},
    "WPRE3_compact":         {"size": 247,  "source": "truncated WPRE (Choi 2014); mRNA stability + translation boost"},
    "WPRE_full":             {"size": 593,  "source": "full WPRE; larger but higher expression boost"},
    "Myo15_956":             {"size": 956,  "source": "Zhao 2024 Research; HC-exclusive in AAV-PHP.eB"},
    "Myo15_1157":            {"size": 1157, "source": "Zhao 2024 Research; HC-exclusive, stronger than 956"},
    "Myo15_1611_native":     {"size": 1611, "source": "Liu 2024 Mol Ther Nucleic Acids; OTOF rescue in DFNB9 mice"},
    "Prestin_native_est":    {"size": 1800, "source": "estimated from Slc26a5 regulatory region; used in some OHC constructs"},
    "Myo15_core":            {"size": 300,  "source": "Myo15 minimal proximal promoter estimate (TATA + BRE + 5'UTR)"},
    "bGH_polyA":             {"size": 225,  "source": "already counted in KOZAK_SP_POLYA"},
    "SV40_polyA":            {"size": 135,  "source": "alternative polyA; shorter than bGH"},
    "IgK_SP":                {"size": 63,   "source": "already counted in KOZAK_SP_POLYA"},
}

# Candidate architectures
CANDIDATES = [
    {
        "name": "01_B8_alone_baseline",
        "description": "Current design: B8 enhancer + Ultra-Mini CDS + IgK SP + polyA. No WPRE.",
        "modules": ["B8_enhancer"],
        "ohc_specificity": 5,           # zero ectopic per ARBITER
        "expression_strength": 3,       # strong but no WPRE boost
        "aav_track_record": 4,          # ARBITER tested in mouse cochlea, specific STRC use pending
    },
    {
        "name": "02_B8_plus_WPRE3",
        "description": "B8 enhancer + WPRE3-compact for posttranscriptional boost. Expression upgrade at modest size cost.",
        "modules": ["B8_enhancer", "WPRE3_compact"],
        "ohc_specificity": 5,
        "expression_strength": 4,       # WPRE3 reliably boosts steady-state mRNA
        "aav_track_record": 4,
    },
    {
        "name": "03_B8_plus_WPRE_full",
        "description": "B8 enhancer + full WPRE. Max expression boost, larger footprint.",
        "modules": ["B8_enhancer", "WPRE_full"],
        "ohc_specificity": 5,
        "expression_strength": 5,
        "aav_track_record": 4,
    },
    {
        "name": "04_Myo15_956_alone",
        "description": "Zhao 2024 truncated Myo15 promoter, HC-exclusive. Replaces B8 entirely with a hair-cell-restricted natural promoter.",
        "modules": ["Myo15_956"],
        "ohc_specificity": 3,           # HC-exclusive, but both IHCs and OHCs (off-target for OHC-specific STRC need)
        "expression_strength": 4,       # published rescue in DFNB9
        "aav_track_record": 5,          # in-vivo hearing rescue published
    },
    {
        "name": "05_Myo15_1611_full",
        "description": "Liu 2024 full native Myo15 promoter, clinically published for OTOF rescue.",
        "modules": ["Myo15_1611_native"],
        "ohc_specificity": 3,           # HC-exclusive, both IHC+OHC
        "expression_strength": 5,       # strongest published hair-cell promoter in AAV
        "aav_track_record": 5,          # OTOF/DFNB9 rescue in Otof-/- mice
    },
    {
        "name": "06_B8_plus_Myo15_core_plus_WPRE3",
        "description": "Belt-and-suspenders: B8 (OHC-restricted enhancer) + Myo15 minimal core (HC-restricted basal promoter) + WPRE3. Dual specificity layer + expression boost.",
        "modules": ["B8_enhancer", "Myo15_core", "WPRE3_compact"],
        "ohc_specificity": 5,           # B8 drives OHC restriction
        "expression_strength": 5,       # Myo15 core as basal + WPRE3 boost
        "aav_track_record": 3,          # novel combination, not directly published
    },
    {
        "name": "07_Prestin_native_plus_WPRE3",
        "description": "Native Slc26a5 (prestin) regulatory region + WPRE3. Uses endogenous OHC-specific promoter. Matches STRC's natural expression pattern (OHC-exclusive).",
        "modules": ["Prestin_native_est", "WPRE3_compact"],
        "ohc_specificity": 5,           # Slc26a5 is OHC-exclusive
        "expression_strength": 4,
        "aav_track_record": 3,          # used in prestin-CreERT2 mice, less AAV data
    },
]


def size_fit_score(total_bp, budget):
    if total_bp > budget:
        return 1, "does not fit"
    spare = budget - total_bp
    if spare < 100:
        return 2, f"fits, only {spare} bp spare"
    if spare < 300:
        return 3, f"fits, {spare} bp spare"
    if spare < 500:
        return 4, f"fits, {spare} bp spare"
    return 5, f"fits, {spare} bp spare"


def score_candidate(c, budget):
    total = sum(MODULES[m]["size"] for m in c["modules"])
    size_score, size_note = size_fit_score(total, budget)
    tier = min(size_score, c["ohc_specificity"], c["expression_strength"], c["aav_track_record"])
    return {
        "name": c["name"],
        "description": c["description"],
        "modules": c["modules"],
        "total_regulatory_bp": total,
        "fits_ultra_mini": total <= available_ultra,
        "spare_bp_ultra_mini": available_ultra - total,
        "fits_mini_700_1775": total <= available_mini,
        "spare_bp_mini_700_1775": available_mini - total,
        "size_fit_score": size_score,
        "size_note": size_note,
        "ohc_specificity_score": c["ohc_specificity"],
        "expression_strength_score": c["expression_strength"],
        "aav_track_record_score": c["aav_track_record"],
        "tier": tier,
    }


scored = [score_candidate(c, available_ultra) for c in CANDIDATES]

print(f"\n{'='*100}")
print(f"{'Candidate':40s}  {'total bp':>8s}  {'U-spare':>8s}  {'M-spare':>8s}  {'size':>4s}  {'OHC':>4s}  {'expr':>4s}  {'AAV':>4s}  {'tier':>4s}")
print("="*100)
for s in scored:
    print(
        f"{s['name']:40s}  "
        f"{s['total_regulatory_bp']:8d}  "
        f"{s['spare_bp_ultra_mini']:8d}  "
        f"{s['spare_bp_mini_700_1775']:8d}  "
        f"{s['size_fit_score']:4d}  "
        f"{s['ohc_specificity_score']:4d}  "
        f"{s['expression_strength_score']:4d}  "
        f"{s['aav_track_record_score']:4d}  "
        f"{s['tier']:4d}"
    )

# Top candidate = highest tier, tie-broken by total module complementarity
best = sorted(scored, key=lambda s: (-s["tier"], s["total_regulatory_bp"]))
print(f"\nTop candidates (Ultra-Mini, ranked by min-axis tier):")
for b in best[:3]:
    print(f"  [{b['tier']}] {b['name']}  ({b['total_regulatory_bp']} bp, spare={b['spare_bp_ultra_mini']} bp)")

# Which candidates do NOT fit in mini 700-1775 but DO fit in Ultra-Mini?
unique_ultra = [s for s in scored if s["fits_ultra_mini"] and not s["fits_mini_700_1775"]]
print(f"\nCandidates enabled exclusively by Ultra-Mini (cannot fit in clinical 700-1775):")
for s in unique_ultra:
    print(f"  - {s['name']}: {s['total_regulatory_bp']} bp (mini 700-1775 would have {s['spare_bp_mini_700_1775']} bp spare)")

out = {
    "aav_capacity_bp": AAV_ITR_CAPACITY,
    "itr_overhead_bp": ITR_OVERHEAD,
    "kozak_sp_polya_overhead_bp": KOZAK_SP_POLYA,
    "ultra_mini_cds_bp": ULTRA_MINI_CDS,
    "mini_700_1775_cds_bp": MINI_700_1775_CDS,
    "ultra_mini_regulatory_budget_bp": available_ultra,
    "mini_700_1775_regulatory_budget_bp": available_mini,
    "upgrade_gain_bp": available_ultra - available_mini,
    "modules": MODULES,
    "candidates": scored,
    "top_3_for_ultra_mini": [b["name"] for b in best[:3]],
    "exclusively_ultra_mini_enabled": [s["name"] for s in unique_ultra],
}
Path("/Users/egorlyfar/Brain/research/strc/models/ultra_mini_promoter_shortlist.json").write_text(json.dumps(out, indent=2))
print(f"\nSaved: ultra_mini_promoter_shortlist.json")
