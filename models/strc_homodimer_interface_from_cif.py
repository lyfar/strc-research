"""STRC homodimer inter-chain contact analysis from existing AF3 CIFs.

Motivation: [[2026-04-17-liang-pcdh15-cryo-em-tip-link]] showed PCDH15 (same
stereocilia mechanosensory system as STRC) forms obligate right-handed double-
helix dimers via a defined oligomerization domain. Parent Mini-STRC hypothesis
flagged: "AF3-Multimer STRC homodimer -- if STRC self-assembles like PCDH15,
oligomerization domain must be in 700-1775. Check first."

Prior AF3 multimer jobs reported low overall confidence (ipTM 0.24 for full
STRC homodimer, 0.20 for mini 594-1775 homodimer). Low ipTM != no contacts -- it
means AF3 is uncertain about overall packing, but any geometric contacts it did
place are still worth inspecting. This script extracts them.

Key question for Ultra-Mini: if ANY geometric contact between the two STRC
chains falls in the LRR stretch aa 700-1074 (the region Ultra-Mini deletes),
Ultra-Mini will break STRC self-assembly. If all contacts are in aa 1075-1775
or if there are none at all, Ultra-Mini's dimerization capacity is unaffected
by the aggressive truncation.

Method:
1. Parse job-b-full-strc-homodimer.cif (2x full STRC, offset=1 per chain) and
   job-c-mini-strc-homodimer.cif (2x mini 594-1775, offset=594 per chain).
2. Use Biopython NeighborSearch at 5 A cutoff between chain A and chain B.
3. Map local resnums -> canonical STRC numbering via offset.
4. Bucket contacts by zone: N-term disordered (1-699), LRR (700-1074),
   Ultra-Mini (1075-1775). Report contact counts per zone, hot contact table.
"""
import json
import warnings
from pathlib import Path
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch

warnings.filterwarnings("ignore")

CIF_DIR = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models")
OUT = Path("/Users/egorlyfar/Brain/research/strc/models/strc_homodimer_interface_from_cif.json")

CONTACT_CUTOFF = 5.0

JOBS = [
    {
        "label": "full_strc_homodimer",
        "cif": "job-b-full-strc-homodimer.cif",
        "offset_A": 1,   "offset_B": 1,
        "construct": "STRC 1-1775 x 2",
    },
    {
        "label": "mini_strc_homodimer",
        "cif": "job-c-mini-strc-homodimer.cif",
        "offset_A": 594, "offset_B": 594,
        "construct": "mini-STRC 594-1775 x 2",
    },
]

# Canonical zone boundaries (STRC full numbering)
ZONES = [
    ("n_term_disordered", 1,    699),
    ("lrr_stretch",       700,  1074),
    ("ultra_mini",        1075, 1775),
]


def zone_of(resnum):
    for name, start, end in ZONES:
        if start <= resnum <= end:
            return name
    return "out_of_range"


parser = MMCIFParser(QUIET=True)

results = {
    "cutoff_A": CONTACT_CUTOFF,
    "canonical_zones": {name: [start, end] for name, start, end in ZONES},
    "rationale": (
        "Ultra-Mini deletes aa 1-1074 (N-term disordered + LRR stretch) "
        "and retains aa 1075-1775. Any inter-chain homodimer contact falling "
        "in aa 1-1074 = contact lost by Ultra-Mini = dimerization compromised. "
        "Any contact in aa 1075-1775 = preserved in Ultra-Mini."
    ),
    "jobs": {},
}

for job in JOBS:
    print(f"\n{'='*72}\n{job['label']}  ({job['construct']})\n{'='*72}")
    structure = parser.get_structure("x", CIF_DIR / job["cif"])
    model = next(structure.get_models())
    chain_A = model["A"]
    chain_B = model["B"]
    atoms_A = [a for a in chain_A.get_atoms() if a.element != "H"]
    atoms_B = [a for a in chain_B.get_atoms() if a.element != "H"]
    ns_B = NeighborSearch(atoms_B)

    # For each atom in A, search B within cutoff
    pairs = {}   # (res_A_canonical, res_B_canonical) -> min atom-atom distance
    for a in atoms_A:
        near = ns_B.search(a.coord, CONTACT_CUTOFF)
        for b in near:
            res_A_local = a.get_parent().id[1]
            res_B_local = b.get_parent().id[1]
            res_A_canon = res_A_local + job["offset_A"] - 1
            res_B_canon = res_B_local + job["offset_B"] - 1
            key = (res_A_canon, res_B_canon)
            d = a - b
            if key not in pairs or d < pairs[key]:
                pairs[key] = d

    # Bucket A-side residues by zone
    zone_counts_A = {name: 0 for name, _, _ in ZONES}
    zone_counts_A["out_of_range"] = 0
    zone_residues_A = {name: set() for name, _, _ in ZONES}
    zone_residues_A["out_of_range"] = set()
    for (rA, rB), d in pairs.items():
        zA = zone_of(rA)
        zone_counts_A[zA] += 1
        zone_residues_A[zA].add(rA)

    # Same for B-side
    zone_counts_B = {name: 0 for name, _, _ in ZONES}
    zone_counts_B["out_of_range"] = 0
    zone_residues_B = {name: set() for name, _, _ in ZONES}
    zone_residues_B["out_of_range"] = set()
    for (rA, rB), d in pairs.items():
        zB = zone_of(rB)
        zone_counts_B[zB] += 1
        zone_residues_B[zB].add(rB)

    # Top contacts
    top = sorted(pairs.items(), key=lambda kv: kv[1])[:20]

    n_total = len(pairs)
    n_ultra_mini_compatible = sum(
        1 for (rA, rB), _ in pairs.items()
        if zone_of(rA) == "ultra_mini" and zone_of(rB) == "ultra_mini"
    )
    n_lrr_involved = sum(
        1 for (rA, rB), _ in pairs.items()
        if zone_of(rA) == "lrr_stretch" or zone_of(rB) == "lrr_stretch"
    )
    n_n_term_involved = sum(
        1 for (rA, rB), _ in pairs.items()
        if zone_of(rA) == "n_term_disordered" or zone_of(rB) == "n_term_disordered"
    )

    print(f"Total inter-chain atom pairs (<{CONTACT_CUTOFF} A): {n_total}")
    print(f"A-side contact residues by zone:")
    for name, _, _ in ZONES:
        print(f"  {name:20s}  {zone_counts_A[name]:4d} pairs  {len(zone_residues_A[name]):3d} unique residues")
    print(f"B-side contact residues by zone:")
    for name, _, _ in ZONES:
        print(f"  {name:20s}  {zone_counts_B[name]:4d} pairs  {len(zone_residues_B[name]):3d} unique residues")
    print(f"\nUltra-Mini-compatible pairs (both residues in aa 1075-1775): {n_ultra_mini_compatible} / {n_total}")
    print(f"LRR-involving pairs (any residue in aa 700-1074):            {n_lrr_involved} / {n_total}")
    print(f"N-term-involving pairs (any residue in aa 1-699):            {n_n_term_involved} / {n_total}")

    print(f"\nTop-20 closest contacts:")
    for (rA, rB), d in top:
        zA, zB = zone_of(rA), zone_of(rB)
        marker = " [ULTRA-MINI OK]" if zA == "ultra_mini" and zB == "ultra_mini" else ""
        if zA == "lrr_stretch" or zB == "lrr_stretch":
            marker = " [LRR - lost by Ultra-Mini]"
        if zA == "n_term_disordered" or zB == "n_term_disordered":
            marker = " [N-term - lost by Ultra-Mini]"
        print(f"  A:{rA:4d} ({zA:20s}) <-> B:{rB:4d} ({zB:20s})  {d:.2f} A{marker}")

    results["jobs"][job["label"]] = {
        "construct": job["construct"],
        "cif": job["cif"],
        "total_inter_chain_atom_pairs": n_total,
        "zone_counts_A_side": zone_counts_A,
        "zone_counts_B_side": zone_counts_B,
        "unique_contact_residues_A": {k: sorted(v) for k, v in zone_residues_A.items()},
        "unique_contact_residues_B": {k: sorted(v) for k, v in zone_residues_B.items()},
        "n_ultra_mini_compatible_pairs": n_ultra_mini_compatible,
        "n_lrr_involving_pairs": n_lrr_involved,
        "n_n_term_involving_pairs": n_n_term_involved,
        "top_20_contacts": [
            {"chain_A_aa": rA, "chain_B_aa": rB, "distance_A": round(d, 2),
             "chain_A_zone": zone_of(rA), "chain_B_zone": zone_of(rB)}
            for (rA, rB), d in top
        ],
    }

# Verdict for Ultra-Mini
print(f"\n{'='*72}\nVERDICT FOR ULTRA-MINI DIMERIZATION\n{'='*72}")
for job_label, job_data in results["jobs"].items():
    n_total = job_data["total_inter_chain_atom_pairs"]
    n_ok = job_data["n_ultra_mini_compatible_pairs"]
    n_lost = job_data["n_lrr_involving_pairs"] + job_data["n_n_term_involving_pairs"]
    frac_preserved = n_ok / n_total if n_total else None
    if n_total == 0:
        v = "no_inter_chain_contacts (AF3 did not place chains in contact)"
    elif frac_preserved is None:
        v = "undetermined"
    elif frac_preserved >= 0.90:
        v = "ultra_mini_preserves_dimerization (>=90% of contacts in 1075-1775)"
    elif frac_preserved >= 0.50:
        v = "partial_preservation (50-90% of contacts in 1075-1775)"
    else:
        v = "ultra_mini_breaks_dimerization (<50% of contacts in 1075-1775)"
    print(f"  {job_label}: {v}")
    job_data["ultra_mini_verdict"] = v
    if n_total:
        job_data["fraction_contacts_preserved_by_ultra_mini"] = round(frac_preserved, 3)

OUT.write_text(json.dumps(results, indent=2, default=float))
print(f"\nSaved: {OUT}")
