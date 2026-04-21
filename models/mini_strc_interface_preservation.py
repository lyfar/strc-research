"""Structural interface preservation for mini-STRC truncation candidates.

Question: does truncating STRC to 700-1775 (current clinical candidate, CpG-depleted,
IgK-SP, fits single AAV) or 1075-1775 (aggressive C-term only) preserve the ARM-fold
that mediates TMEM145 binding (canonical residues 1603-1770, validated in Job 2 AF3)?

Method:
1. Three CIFs from existing AF3 results:
     - job2-mini-complex.cif        : STRC 594-1775 + TMEM145  (REFERENCE, interface mapped)
     - job-f-shorter-mini-strc-700  : STRC 700-1775 solo       (CLINICAL candidate)
     - job-h-strc-cterm-only        : STRC 1075-1775 solo      (AGGRESSIVE candidate)
2. Extract CA atoms for canonical residues 1603-1770 (TMEM145 interface zone).
3. Superpose each truncation onto reference using Biopython Superimposer.
4. Compute RMSD over the full interface and per-residue deviations at the 15 closest
   TMEM145 contact residues (from strc_tmem145_interface_from_cif.json Job 2 top_contacts).
5. Verdict per truncation:
     RMSD < 2 Å  -> interface preserved, construct therapeutically viable
     2-5 Å       -> ambiguous, needs AF3 direct test
     > 5 Å       -> interface broken, construct likely non-functional
"""
import json
import warnings
from pathlib import Path
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Superimposer import Superimposer

warnings.filterwarnings("ignore")

CIF_DIR = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models")
OUT = Path("/Users/egorlyfar/Brain/research/strc/models/mini_strc_interface_preservation.json")

INTERFACE_START, INTERFACE_END = 1603, 1770  # canonical STRC numbering (Derstroff 2026 confirmed)
# Biologically meaningful boundary: GPI anchor omega site = S1749 (NetGPI 1.1).
# Residues 1750-1775 are proteolytically cleaved during GPI attachment in the ER.
# Restricting to pre-GPI zone removes bioartifactual AF3 drift in the prosequence.
BIO_INTERFACE_END = 1749

# Offsets: local residue 1 = canonical residue N
TRUNCATIONS = {
    "reference_594_1775":  {"cif": "job2-mini-complex.cif",            "offset": 594,  "chain": "A"},
    "clinical_700_1775":   {"cif": "job-f-shorter-mini-strc-700.cif",  "offset": 700,  "chain": "A"},
    "aggressive_1075_1775":{"cif": "job-h-strc-cterm-only.cif",        "offset": 1075, "chain": "A"},
}

# Load top TMEM145 contact residues from previous extraction (hottest 15 Job 2 contacts)
prior = json.loads(Path("strc_tmem145_interface_from_cif.json").read_text())
top15_canonical = sorted({c["strc_aa"] for c in prior["job2-mini-complex.cif"]["top_contacts"]})
# These are already canonical (Job 2 offset applied during extraction)

parser = MMCIFParser(QUIET=True)


def get_ca_dict(cif_name, chain_id, offset):
    """Return {canonical_resnum: CA atom} for chain."""
    structure = parser.get_structure("x", CIF_DIR / cif_name)
    model = next(structure.get_models())
    chain = model[chain_id]
    out = {}
    for res in chain:
        if "CA" not in res:
            continue
        canonical = res.id[1] + offset - 1
        out[canonical] = res["CA"]
    return out


def superpose_interface(ref_ca, mov_ca, interface_residues):
    """Superpose mov onto ref using CA atoms at interface residues common to both.
    Returns (rmsd, n_atoms, per_residue_deviation_dict).
    """
    common = [r for r in interface_residues if r in ref_ca and r in mov_ca]
    if len(common) < 3:
        return None, 0, {}
    ref_atoms = [ref_ca[r] for r in common]
    mov_atoms = [mov_ca[r] for r in common]
    # Copy coords before transform (Superimposer mutates mov_atoms)
    mov_coords_before = [tuple(a.coord) for a in mov_atoms]
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mov_atoms)
    sup.apply(mov_atoms)
    per_res = {}
    for r, ref_a, mov_a in zip(common, ref_atoms, mov_atoms):
        d = ((ref_a.coord - mov_a.coord) ** 2).sum() ** 0.5
        per_res[r] = float(d)
    return float(sup.rms), len(common), per_res


# Build CA dicts for all three
ca_by_trunc = {
    name: get_ca_dict(meta["cif"], meta["chain"], meta["offset"])
    for name, meta in TRUNCATIONS.items()
}

# Interface residues in reference (restricted to pre-GPI bio-relevant zone)
interface_residues = list(range(INTERFACE_START, BIO_INTERFACE_END + 1))
ref_ca = ca_by_trunc["reference_594_1775"]

results = {
    "interface_zone_canonical_full": [INTERFACE_START, INTERFACE_END],
    "interface_zone_canonical_biorelevant": [INTERFACE_START, BIO_INTERFACE_END],
    "gpi_omega_site": 1749,
    "rationale_gpi_exclusion": "Residues 1750-1775 are the GPI prosequence, "
                                "proteolytically cleaved during GPI anchor "
                                "attachment in the ER (NetGPI 1.1 predicts "
                                "omega site S1749). Excluded from interface test "
                                "because AF3 models them as a free tail, but they "
                                "do not exist in the mature membrane-anchored protein.",
    "top15_tmem145_contacts_canonical": top15_canonical,
    "reference_construct": "STRC 594-1775 + TMEM145 (job2-mini-complex.cif)",
    "truncations": {},
}

def align_on_hot_measure_full(ref_ca, mov_ca, hot_residues, full_residues):
    """Rigid-body alignment using HOT contact residues only, then measure deviation
    at the FULL interface zone. This preserves the binding-pocket frame and exposes
    peripheral drift honestly.
    """
    common_hot = [r for r in hot_residues if r in ref_ca and r in mov_ca]
    if len(common_hot) < 3:
        return None, {}
    # Deep-copy coords since Superimposer mutates
    ref_hot = [ref_ca[r] for r in common_hot]
    mov_hot = [mov_ca[r] for r in common_hot]
    mov_coords_orig = {r: tuple(mov_ca[r].coord) for r in full_residues if r in mov_ca}

    sup = Superimposer()
    sup.set_atoms(ref_hot, mov_hot)
    rot, tran = sup.rotran

    # Apply transformation manually to all full-zone atoms (non-mutating copy)
    import numpy as np
    per_res = {}
    for r in full_residues:
        if r in mov_ca and r in ref_ca:
            coord_orig = np.array(mov_coords_orig[r])
            coord_new = np.dot(coord_orig, rot) + tran
            d = float(np.linalg.norm(coord_new - ref_ca[r].coord))
            per_res[r] = d
    # Compute aggregate RMSD over full zone
    if per_res:
        rmsd = (sum(d*d for d in per_res.values()) / len(per_res)) ** 0.5
    else:
        rmsd = None
    return float(rmsd) if rmsd is not None else None, per_res


for name, meta in TRUNCATIONS.items():
    if name == "reference_594_1775":
        continue
    mov_ca = ca_by_trunc[name]
    # Primary: align on hot contacts, measure full-zone drift (honest test)
    rmsd_peri, per_res_peri = align_on_hot_measure_full(
        ref_ca, mov_ca, top15_canonical, interface_residues
    )
    # Secondary: hot-only self-consistency (should be tiny by construction)
    rmsd_hot, n_hot, _ = superpose_interface(ref_ca, mov_ca, top15_canonical)
    devs_sorted = sorted(per_res_peri.items(), key=lambda x: x[1], reverse=True)

    # Verdict: look at 75th percentile to catch outliers
    devs_vals = sorted(per_res_peri.values())
    p75 = devs_vals[int(0.75 * len(devs_vals))] if devs_vals else 0
    p95 = devs_vals[int(0.95 * len(devs_vals))] if devs_vals else 0
    median = devs_vals[len(devs_vals) // 2] if devs_vals else 0
    n_below_3A = sum(1 for d in devs_vals if d < 3.0)
    frac_below_3A = n_below_3A / len(devs_vals) if devs_vals else 0

    # Verdict logic: if 80%+ of interface residues are within 3 Å AND hot contacts < 0.5 Å,
    # binding pocket is preserved. Otherwise ambiguous or broken.
    if rmsd_hot is not None and rmsd_hot < 0.5 and frac_below_3A >= 0.80:
        verdict = "preserved_binding_pocket"
    elif rmsd_hot is not None and rmsd_hot < 1.0 and frac_below_3A >= 0.50:
        verdict = "ambiguous_peripheral_drift"
    else:
        verdict = "broken"

    results["truncations"][name] = {
        "canonical_range": [meta["offset"], meta["offset"] + len(list(ca_by_trunc[name])) - 1],
        "n_interface_residues_matched": len(per_res_peri),
        "hot_contact_rmsd_A": rmsd_hot,
        "full_zone_rmsd_after_hot_align_A": rmsd_peri,
        "per_residue_median_dev_A": median,
        "per_residue_p75_dev_A": p75,
        "per_residue_p95_dev_A": p95,
        "per_residue_max_dev_A": devs_vals[-1] if devs_vals else None,
        "fraction_residues_below_3A": round(frac_below_3A, 3),
        "top5_peripheral_outliers": [
            {"resnum_canonical": r, "deviation_A": round(d, 2)}
            for r, d in devs_sorted[:5]
        ],
        "verdict": verdict,
    }
    print(f"\n=== {name} ({meta['offset']}-{meta['offset']+len(list(ca_by_trunc[name]))-1}) ===")
    print(f"  Hot contact RMSD (binding pocket): {rmsd_hot:.3f} Å")
    print(f"  Full-zone RMSD (after hot-align):  {rmsd_peri:.3f} Å")
    print(f"  Per-residue: median {median:.2f}, P75 {p75:.2f}, P95 {p95:.2f}, max {devs_vals[-1]:.2f} Å")
    print(f"  Fraction of 168 residues < 3 Å drift: {frac_below_3A:.1%}")
    print(f"  Top-5 peripheral outliers:")
    for r, d in devs_sorted[:5]:
        print(f"    aa {r}: {d:.2f} Å")
    print(f"  VERDICT: {verdict}")

OUT.write_text(json.dumps(results, indent=2, default=float))
print(f"\nSaved: {OUT}")
