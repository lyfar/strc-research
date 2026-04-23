#!/usr/bin/env python3
"""
Hydrogel Phase 4a — CIF interface forensics.

Parses the 4 Phase 3b AF3 CIFs (tail91/tail71 × TMEM145/actin) and
extracts the ACTUAL binding mode:
  - contact residues at 4.5 Å heavy-atom cutoff
  - per-residue B-factor (pLDDT*100 in AF3 CIFs)
  - SASA by residue
  - which peptide residues dock where

Primary questions:
  Q1. Does tail91 × TMEM145 (ipTM 0.57) use the GOLD-zone contacts we
      designed for? (clusters 2/3/4/5 = aa 1620-1707)
  Q2. Does tail91 × actin (ipTM 0.51) use canonical WH2 actin contacts?
      (WH2 Y/F/I hydrophobic triad → G-actin subdomain 1/3 cleft)
  Q3. Why does tail71 × TMEM145 FAIL (0.35)? Which contacts are lost?
  Q4. Is the RADA16 scaffold inserting into either interface (would
      explain tail71 TMEM145 failure)?

No Biopython — parse CIF ATOM records directly. numpy + scipy for KD-tree.
"""

from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

BATCH = Path("/Users/egorlyfar/Brain/research/strc/models/af3_jobs_2026-04-23e_hydrogel_phase3b")
OUT = Path("/Users/egorlyfar/Brain/research/strc/models/hydrogel_phase4a_cif_interface_forensics.json")

CUTOFF = 4.5  # Å, standard interface heavy-atom cutoff

# Construct reference (from builder)
WH2_NTERM = "RQLVKAIPDNCSKSNVSR"   # 18 aa
RADA16 = "ADARADARADARADA"          # 15 aa
LINK_5 = "GSGSG"

# GOLD-zone tail from STRC Q7RTU9
def fetch_uniprot(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    return "".join(lines[1:])

STRC_FULL = fetch_uniprot("Q7RTU9")
TAIL_71 = STRC_FULL[1639:1710]  # 71 aa, aa 1640-1710
TAIL_91 = STRC_FULL[1619:1710]  # 91 aa, aa 1620-1710

PEPTIDE_TAIL71 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_71
PEPTIDE_TAIL91 = WH2_NTERM + LINK_5 + RADA16 + LINK_5 + TAIL_91
assert len(PEPTIDE_TAIL71) == 114
assert len(PEPTIDE_TAIL91) == 134

# Canonical GOLD-zone contact clusters from Derstroff et al. + our Phase 1 design
GOLD_CLUSTERS = {
    "cluster_2": (1620, 1638),   # N-term hump, aa 1620-1638 (19 aa) — only in tail91
    "cluster_3": (1648, 1651),   # loop cap, aa 1648-1651 (4 aa) — critical
    "cluster_4": (1669, 1680),   # dominant, aa 1669-1680 (12 aa)
    "cluster_5": (1692, 1707),   # second dominant, aa 1692-1707 (16 aa)
}

# WH2 canonical hydrophobic triad (typical YxxI/F motif in WH2 helix)
WH2_HYDROPHOBIC_RESIDUES = ["I", "L", "V", "F", "Y", "W", "M"]

# Peptide residue index boundaries (1-indexed in our construct)
def peptide_regions(length: int):
    """Residue ranges for WH2 / GS1 / RADA16 / GS2 / tail in the construct."""
    r = {}
    p = 1
    r["wh2"] = (p, p + 18 - 1);    p += 18
    r["gs1"] = (p, p + 5 - 1);     p += 5
    r["rada16"] = (p, p + 15 - 1); p += 15
    r["gs2"] = (p, p + 5 - 1);     p += 5
    r["tail"] = (p, p + (length - p) - 1 + 1)
    return r


def parse_cif_atoms(cif_path: Path):
    """Return list[dict] of atom records with chain, resnum, resname, name,
    element, xyz, B-factor. Parses only _atom_site loop."""
    lines = cif_path.read_text().splitlines()
    # Find _atom_site loop
    in_loop = False
    headers = []
    atoms = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "loop_":
            # Next lines are headers until we get a non-header
            j = i + 1
            headers = []
            while j < len(lines) and lines[j].startswith("_"):
                headers.append(lines[j].strip())
                j += 1
            if any(h.startswith("_atom_site.") for h in headers):
                in_loop = True
                # parse rows from j until '#' or next 'loop_' or blank
                col = {h: k for k, h in enumerate(headers)}
                while j < len(lines):
                    row = lines[j].strip()
                    if row == "" or row == "#" or row == "loop_":
                        break
                    parts = row.split()
                    if len(parts) < len(headers):
                        j += 1; continue
                    atoms.append({
                        "chain": parts[col.get("_atom_site.auth_asym_id",
                                               col.get("_atom_site.label_asym_id", 0))],
                        "resnum": int(parts[col.get("_atom_site.auth_seq_id",
                                                    col.get("_atom_site.label_seq_id", 8))]),
                        "resname": parts[col.get("_atom_site.auth_comp_id",
                                                 col.get("_atom_site.label_comp_id", 5))],
                        "atom_name": parts[col.get("_atom_site.label_atom_id", 3)].strip('"'),
                        "element": parts[col.get("_atom_site.type_symbol", 2)],
                        "xyz": np.array([
                            float(parts[col["_atom_site.Cartn_x"]]),
                            float(parts[col["_atom_site.Cartn_y"]]),
                            float(parts[col["_atom_site.Cartn_z"]]),
                        ]),
                        "bfac": float(parts[col.get("_atom_site.B_iso_or_equiv", 14)]),
                    })
                    j += 1
                i = j
                continue
            else:
                i = j
                continue
        i += 1
    return atoms


def contacts_between_chains(atoms, chain_a, chain_b, cutoff=CUTOFF):
    """Return list[(resA_num, resA_name, resB_num, resB_name, min_dist)]."""
    heavy_a = [a for a in atoms if a["chain"] == chain_a and a["element"] != "H"]
    heavy_b = [a for a in atoms if a["chain"] == chain_b and a["element"] != "H"]
    if not heavy_a or not heavy_b:
        return []
    xyz_a = np.array([a["xyz"] for a in heavy_a])
    xyz_b = np.array([a["xyz"] for a in heavy_b])
    tree_b = cKDTree(xyz_b)
    pairs = {}
    for i, a in enumerate(heavy_a):
        dists, idxs = tree_b.query(xyz_a[i], k=min(50, len(xyz_b)),
                                   distance_upper_bound=cutoff)
        dists = np.atleast_1d(dists); idxs = np.atleast_1d(idxs)
        for d, j in zip(dists, idxs):
            if np.isinf(d) or j >= len(heavy_b):
                continue
            b = heavy_b[j]
            key = (a["resnum"], a["resname"], b["resnum"], b["resname"])
            if key not in pairs or d < pairs[key]:
                pairs[key] = d
    return [(k[0], k[1], k[2], k[3], float(v)) for k, v in pairs.items()]


def residues_by_chain(atoms):
    """Return dict[chain] -> sorted list of (resnum, resname)."""
    out = {}
    seen = set()
    for a in atoms:
        k = (a["chain"], a["resnum"])
        if k in seen:
            continue
        seen.add(k)
        out.setdefault(a["chain"], []).append((a["resnum"], a["resname"]))
    for c in out:
        out[c].sort()
    return out


def best_model_dir(job_name: str) -> Path:
    """Pick model 0 (highest rank) CIF."""
    jd = BATCH / "results" / job_name
    return jd / f"{job_name}_model_0.cif"


def peptide_residue_region(resnum: int, construct_len: int) -> str:
    regs = peptide_regions(construct_len)
    for name, (lo, hi) in regs.items():
        if lo <= resnum <= hi:
            return name
    return "out_of_range"


def tail_resnum_to_strc(resnum: int, construct_len: int) -> int | None:
    """Map construct peptide residue to STRC position (aa 1620 or 1640 start)."""
    regs = peptide_regions(construct_len)
    tlo, thi = regs["tail"]
    if not (tlo <= resnum <= thi):
        return None
    offset = resnum - tlo  # 0-indexed into tail
    # tail91 starts at STRC 1620; tail71 starts at STRC 1640
    tail_start_strc = 1620 if construct_len == 134 else 1640
    return tail_start_strc + offset


def classify_strc_cluster(strc_resnum: int) -> str:
    for name, (lo, hi) in GOLD_CLUSTERS.items():
        if lo <= strc_resnum <= hi:
            return name
    if strc_resnum < 1620:
        return "pre_cluster_2"
    if 1638 < strc_resnum < 1648:
        return "linker_23"
    if 1651 < strc_resnum < 1669:
        return "linker_34"
    if 1680 < strc_resnum < 1692:
        return "linker_45"
    if strc_resnum > 1707:
        return "post_cluster_5"
    return "unclassified"


def analyze_tmem145_job(job_name: str, construct_len: int) -> dict:
    cif = best_model_dir(job_name)
    atoms = parse_cif_atoms(cif)
    by_chain = residues_by_chain(atoms)
    # Peptide is chain A (count 1 first in job def), TMEM145 is chain B
    chains = sorted(by_chain.keys())
    # Peptide chain = the shorter one
    chain_lens = {c: len(by_chain[c]) for c in chains}
    peptide_chain = min(chain_lens, key=chain_lens.get)
    tmem145_chain = max(chain_lens, key=chain_lens.get)

    contacts = contacts_between_chains(atoms, peptide_chain, tmem145_chain)
    contacts.sort(key=lambda x: x[4])

    # Classify peptide contacts by region
    region_counts = {"wh2": 0, "gs1": 0, "rada16": 0, "gs2": 0, "tail": 0, "out_of_range": 0}
    tail_cluster_counts = {k: 0 for k in GOLD_CLUSTERS}
    tail_linker_counts = {"pre_cluster_2": 0, "linker_23": 0, "linker_34": 0,
                          "linker_45": 0, "post_cluster_5": 0, "unclassified": 0}
    for (rnum, rname, _tnum, _tname, _d) in contacts:
        region = peptide_residue_region(rnum, construct_len)
        region_counts[region] += 1
        if region == "tail":
            strc_rn = tail_resnum_to_strc(rnum, construct_len)
            cluster = classify_strc_cluster(strc_rn)
            if cluster in tail_cluster_counts:
                tail_cluster_counts[cluster] += 1
            else:
                tail_linker_counts[cluster] = tail_linker_counts.get(cluster, 0) + 1

    # Per-residue pLDDT on peptide (avg B-factor of CA atoms)
    pep_plddt = {}
    for a in atoms:
        if a["chain"] == peptide_chain and a["atom_name"] == "CA":
            pep_plddt[a["resnum"]] = a["bfac"]

    # Which contacts involve RADA16 scaffold?
    scaffold_contacts = [
        (rnum, rname, tnum, tname, d)
        for (rnum, rname, tnum, tname, d) in contacts
        if peptide_residue_region(rnum, construct_len) in ("rada16", "gs1", "gs2")
    ]

    return {
        "job_name": job_name,
        "construct_len": construct_len,
        "peptide_chain": peptide_chain,
        "tmem145_chain": tmem145_chain,
        "n_contacts_total": len(contacts),
        "contacts_by_peptide_region": region_counts,
        "tail_contacts_by_gold_cluster": tail_cluster_counts,
        "tail_contacts_by_non_cluster": tail_linker_counts,
        "n_scaffold_contacts": len(scaffold_contacts),
        "scaffold_contact_details": [
            {"pep_residue": f"{r[1]}{r[0]}", "tmem145_residue": f"{r[3]}{r[2]}",
             "distance": round(r[4], 2)}
            for r in scaffold_contacts[:20]
        ],
        "top_10_contacts": [
            {"peptide": f"{r[1]}{r[0]}",
             "peptide_region": peptide_residue_region(r[0], construct_len),
             "tmem145": f"{r[3]}{r[2]}",
             "distance": round(r[4], 2),
             "peptide_plddt": round(pep_plddt.get(r[0], 0), 1)}
            for r in contacts[:10]
        ],
        "peptide_plddt_mean_at_tail": round(float(np.mean([
            pep_plddt[r] for r in pep_plddt
            if peptide_regions(construct_len)["tail"][0] <= r <= peptide_regions(construct_len)["tail"][1]
        ])), 1),
        "peptide_plddt_mean_at_scaffold": round(float(np.mean([
            pep_plddt[r] for r in pep_plddt
            if peptide_regions(construct_len)["rada16"][0] <= r <= peptide_regions(construct_len)["rada16"][1]
        ])), 1),
        "peptide_plddt_mean_at_wh2": round(float(np.mean([
            pep_plddt[r] for r in pep_plddt
            if peptide_regions(construct_len)["wh2"][0] <= r <= peptide_regions(construct_len)["wh2"][1]
        ])), 1),
    }


def analyze_actin_job(job_name: str, construct_len: int) -> dict:
    cif = best_model_dir(job_name)
    atoms = parse_cif_atoms(cif)
    by_chain = residues_by_chain(atoms)
    chains = sorted(by_chain.keys())
    chain_lens = {c: len(by_chain[c]) for c in chains}
    peptide_chain = min(chain_lens, key=chain_lens.get)
    actin_chains = [c for c in chains if c != peptide_chain]

    # Contacts to each actin protomer
    per_actin = {}
    total_contacts = 0
    wh2_contacts_total = 0
    for ac in actin_chains:
        cs = contacts_between_chains(atoms, peptide_chain, ac)
        wh2_cs = [c for c in cs if peptide_residue_region(c[0], construct_len) == "wh2"]
        per_actin[ac] = {
            "n_contacts": len(cs),
            "n_wh2_contacts": len(wh2_cs),
            "top_wh2_contacts": [
                {"peptide": f"{c[1]}{c[0]}", "actin": f"{c[3]}{c[2]}",
                 "distance": round(c[4], 2)}
                for c in sorted(wh2_cs, key=lambda x: x[4])[:6]
            ],
        }
        total_contacts += len(cs)
        wh2_contacts_total += len(wh2_cs)

    # Aggregate region distribution (from peptide side, to any actin)
    all_cs = []
    for ac in actin_chains:
        all_cs += contacts_between_chains(atoms, peptide_chain, ac)
    region_counts = {"wh2": 0, "gs1": 0, "rada16": 0, "gs2": 0, "tail": 0, "out_of_range": 0}
    for (rnum, _rname, _tnum, _tname, _d) in all_cs:
        region_counts[peptide_residue_region(rnum, construct_len)] += 1

    # Canonical WH2 residue analysis (our WH2: RQLVKAIPDNCSKSNVSR)
    # Hydrophobic triad typically at positions 2-4, 6-7 in WH2 helix
    # Our seq positions 3=L, 4=V, 6=A, 7=I → LV…AI likely the anchor triad
    wh2_hydrophobic_positions = []
    for i, aa in enumerate(WH2_NTERM, start=1):  # 1-indexed within WH2
        if aa in WH2_HYDROPHOBIC_RESIDUES:
            wh2_hydrophobic_positions.append(i)  # constant offset to construct 1-indexed

    return {
        "job_name": job_name,
        "construct_len": construct_len,
        "peptide_chain": peptide_chain,
        "actin_chains": actin_chains,
        "n_contacts_total": total_contacts,
        "n_contacts_wh2": wh2_contacts_total,
        "wh2_share": round(wh2_contacts_total / total_contacts, 3) if total_contacts else 0,
        "contacts_by_peptide_region": region_counts,
        "per_actin_protomer": per_actin,
        "wh2_hydrophobic_positions": wh2_hydrophobic_positions,
        "wh2_sequence": WH2_NTERM,
    }


def main():
    tmem_jobs = [
        ("fold_hydrogel_rada16_wh2_tail71_x_tmem145", 114),
        ("fold_hydrogel_rada16_wh2_tail91_x_tmem145", 134),
    ]
    actin_jobs = [
        ("fold_hydrogel_rada16_wh2_tail71_x_actin", 114),
        ("fold_hydrogel_rada16_wh2_tail91_x_actin", 134),
    ]

    results = {"tmem145": {}, "actin": {}}
    for name, clen in tmem_jobs:
        print(f"Parsing {name}…")
        results["tmem145"][name] = analyze_tmem145_job(name, clen)
    for name, clen in actin_jobs:
        print(f"Parsing {name}…")
        results["actin"][name] = analyze_actin_job(name, clen)

    # Interpretation
    tmem91 = results["tmem145"]["fold_hydrogel_rada16_wh2_tail91_x_tmem145"]
    tmem71 = results["tmem145"]["fold_hydrogel_rada16_wh2_tail71_x_tmem145"]

    # Fraction of peptide contacts that are from tail vs scaffold
    def tail_frac(r):
        tot = r["n_contacts_total"]
        return round(r["contacts_by_peptide_region"]["tail"] / tot, 3) if tot else 0

    def scaffold_frac(r):
        tot = r["n_contacts_total"]
        s = (r["contacts_by_peptide_region"]["rada16"]
             + r["contacts_by_peptide_region"]["gs1"]
             + r["contacts_by_peptide_region"]["gs2"])
        return round(s / tot, 3) if tot else 0

    interpretation = {
        "Q1_tail91_tmem145_uses_gold_zone": {
            "tail_contact_fraction": tail_frac(tmem91),
            "gold_cluster_counts": tmem91["tail_contacts_by_gold_cluster"],
            "scaffold_contact_fraction": scaffold_frac(tmem91),
            "verdict": (
                "DESIGNED_EPITOPE" if tail_frac(tmem91) > 0.75
                and tmem91["tail_contacts_by_gold_cluster"]["cluster_4"] > 0
                else "OFF_TARGET_BINDING_MODE"
            ),
        },
        "Q3_why_tail71_fails": {
            "tail71_tail_frac": tail_frac(tmem71),
            "tail71_scaffold_frac": scaffold_frac(tmem71),
            "tail71_gold_counts": tmem71["tail_contacts_by_gold_cluster"],
            "tail91_tail_frac": tail_frac(tmem91),
            "tail91_scaffold_frac": scaffold_frac(tmem91),
            "delta_scaffold_frac": round(scaffold_frac(tmem71) - scaffold_frac(tmem91), 3),
            "verdict": (
                "SCAFFOLD_INVASION" if scaffold_frac(tmem71) > scaffold_frac(tmem91) + 0.1
                else "WEAK_TAIL_ANCHOR"
            ),
        },
        "Q2_actin_uses_wh2": {
            name: {
                "n_wh2_contacts": r["n_contacts_wh2"],
                "wh2_share": r["wh2_share"],
                "verdict": "CANONICAL_WH2" if r["wh2_share"] > 0.3 else "NON_CANONICAL"
            }
            for name, r in results["actin"].items()
        },
    }

    # Headline summary
    summary = {
        "batch": "hydrogel_phase4a_cif_interface_forensics",
        "date": "2026-04-23",
        "cutoff_angstrom": CUTOFF,
        "construct_designs": {
            "tail71": {"length": 114, "peptide": PEPTIDE_TAIL71},
            "tail91": {"length": 134, "peptide": PEPTIDE_TAIL91},
        },
        "gold_clusters": GOLD_CLUSTERS,
        "raw_results": results,
        "interpretation": interpretation,
    }
    OUT.write_text(json.dumps(summary, indent=2))

    print("\n=== Phase 4a CIF Interface Forensics ===\n")
    print(f"Q1. Does tail91 × TMEM145 engage the GOLD zone we designed?")
    q1 = interpretation["Q1_tail91_tmem145_uses_gold_zone"]
    print(f"    tail contact fraction = {q1['tail_contact_fraction']}")
    print(f"    scaffold contact fraction = {q1['scaffold_contact_fraction']}")
    print(f"    GOLD cluster contact breakdown = {q1['gold_cluster_counts']}")
    print(f"    >>> {q1['verdict']}")
    print()
    print(f"Q3. Why does tail71 × TMEM145 FAIL?")
    q3 = interpretation["Q3_why_tail71_fails"]
    print(f"    tail71 scaffold frac = {q3['tail71_scaffold_frac']} vs tail91 {q3['tail91_scaffold_frac']}")
    print(f"    tail71 cluster counts = {q3['tail71_gold_counts']}")
    print(f"    >>> {q3['verdict']}")
    print()
    print(f"Q2. Does tail91 × actin use canonical WH2 contacts?")
    for name, v in interpretation["Q2_actin_uses_wh2"].items():
        print(f"    {name}: WH2 share = {v['wh2_share']}, verdict = {v['verdict']}")
    print()
    print(f"JSON written: {OUT}")


if __name__ == "__main__":
    main()
