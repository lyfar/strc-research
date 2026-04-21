#!/usr/bin/env python3
"""
STRC Pharmacochaperone Phase 1B — High-Confidence Realignment

Fix for Phase 1: the 14.9 Å global RMSD came from aligning on disordered
low-pLDDT regions. Realign WT onto MUT using only high-confidence residues
(pLDDT > 70), then recompute mutation impact.

Output:
  pharmacochaperone_phase1b_results.json
  pharmacochaperone_phase1b.png
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, Superimposer
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

JOB_MUT = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job3-mutant.cif")
JOB_WT  = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job4-wildtype.cif")

E1659 = 1659
K1141 = 1141
PLDDT_CUT = 70

HIGH_CONF_WINDOW_A = 10.0


def load(path: Path):
    p = MMCIFParser(QUIET=True)
    return p.get_structure(path.stem, str(path))


def chain_residues(s) -> Dict[int, object]:
    return {r.id[1]: r for r in s[0]["A"] if r.id[0] == " "}


def high_conf_residues(chain_dict: Dict, cut: float = PLDDT_CUT) -> set:
    out = set()
    for n, r in chain_dict.items():
        if "CA" in r and r["CA"].get_bfactor() >= cut:
            out.add(n)
    return out


def superimpose_on(wt_dict, mut_dict, keep_residues: set) -> Tuple[float, Superimposer]:
    wt_ca, mut_ca = [], []
    for n in sorted(keep_residues):
        if n in wt_dict and n in mut_dict and "CA" in wt_dict[n] and "CA" in mut_dict[n]:
            wt_ca.append(wt_dict[n]["CA"])
            mut_ca.append(mut_dict[n]["CA"])
    sup = Superimposer()
    sup.set_atoms(wt_ca, mut_ca)
    # Apply to MUT structure (transform all atoms)
    mut_model = next(iter(mut_dict.values())).get_parent().get_parent()
    sup.apply(list(mut_model.get_atoms()))
    return sup.rms, sup


def positive_tip(res):
    if res is None: return None
    atoms = {a.name: a for a in res.get_atoms()}
    if res.resname == "LYS" and "NZ" in atoms:
        return np.array(atoms["NZ"].coord)
    if res.resname == "ARG":
        t = [atoms[n].coord for n in ("NE","NH1","NH2") if n in atoms]
        if t: return np.mean(t, axis=0)
    return None


def negative_tip(res):
    if res is None: return None
    atoms = {a.name: a for a in res.get_atoms()}
    if res.resname == "GLU":
        t = [atoms[n].coord for n in ("OE1","OE2") if n in atoms]
        if t: return np.mean(t, axis=0)
    if res.resname == "ASP":
        t = [atoms[n].coord for n in ("OD1","OD2") if n in atoms]
        if t: return np.mean(t, axis=0)
    return None


def side_chain_centroid(res):
    side = [a for a in res.get_atoms()
            if a.name not in ("N","CA","C","O") and a.element != "H"]
    if not side: return np.array(res["CA"].coord)
    return np.mean([a.coord for a in side], axis=0)


def local_rmsd(wt_dict, mut_dict, center: np.ndarray, radius: float):
    pairs = []
    for n, r_wt in wt_dict.items():
        if "CA" not in r_wt: continue
        if np.linalg.norm(np.array(r_wt["CA"].coord) - center) > radius: continue
        r_mut = mut_dict.get(n)
        if r_mut is None or "CA" not in r_mut: continue
        pairs.append((np.array(r_wt["CA"].coord), np.array(r_mut["CA"].coord), n,
                      r_wt["CA"].get_bfactor(), r_mut["CA"].get_bfactor()))
    if not pairs:
        return 0.0, 0, []
    sq = [np.sum((a-b)**2) for a,b,_,_,_ in pairs]
    return float(np.sqrt(np.mean(sq))), len(pairs), pairs


def per_residue_displacement(wt_dict, mut_dict, resrange):
    rows = []
    for n in resrange:
        if n in wt_dict and n in mut_dict and "CA" in wt_dict[n] and "CA" in mut_dict[n]:
            d = np.linalg.norm(np.array(wt_dict[n]["CA"].coord) - np.array(mut_dict[n]["CA"].coord))
            rows.append({"resnum": n, "displacement_A": round(float(d), 3),
                         "WT_plddt": round(wt_dict[n]["CA"].get_bfactor(), 1),
                         "MUT_plddt": round(mut_dict[n]["CA"].get_bfactor(), 1),
                         "WT_resname": wt_dict[n].resname,
                         "MUT_resname": mut_dict[n].resname})
    return rows


def plot_phase1b(displacements, pairs_local, out: Path):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: CA displacement along sequence
    ax = axes[0, 0]
    nums = [d["resnum"] for d in displacements]
    vals = [d["displacement_A"] for d in displacements]
    plddt = [(d["WT_plddt"] + d["MUT_plddt"]) / 2 for d in displacements]
    sc = ax.scatter(nums, vals, c=plddt, cmap="viridis", s=12, vmin=50, vmax=100)
    ax.axvline(E1659, color="red", ls="--", label=f"E{E1659}A")
    ax.axvline(K1141, color="blue", ls=":", label=f"K{K1141}")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("Cα displacement WT→MUT (Å)")
    ax.set_title("Per-residue CA displacement (after high-conf alignment)")
    plt.colorbar(sc, ax=ax, label="mean pLDDT")
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel 2: histogram of all-residue displacements
    ax = axes[0, 1]
    vals_hc = [d["displacement_A"] for d in displacements if d["WT_plddt"] >= 70 and d["MUT_plddt"] >= 70]
    vals_lc = [d["displacement_A"] for d in displacements if d["WT_plddt"] < 70 or d["MUT_plddt"] < 70]
    if vals_hc:
        ax.hist(vals_hc, bins=50, alpha=0.7, label=f"pLDDT≥70 (n={len(vals_hc)})", color="steelblue")
    if vals_lc:
        ax.hist(vals_lc, bins=50, alpha=0.5, label=f"pLDDT<70 (n={len(vals_lc)})", color="lightcoral")
    ax.set_xlabel("Displacement (Å)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of CA displacements")
    ax.legend()
    ax.set_yscale("log")

    # Panel 3: zoomed E1659 region
    ax = axes[1, 0]
    zoom = [d for d in displacements if 1600 <= d["resnum"] <= 1720]
    zn = [d["resnum"] for d in zoom]
    zv = [d["displacement_A"] for d in zoom]
    zp = [(d["WT_plddt"] + d["MUT_plddt"]) / 2 for d in zoom]
    sc = ax.scatter(zn, zv, c=zp, cmap="viridis", s=40, vmin=50, vmax=100)
    ax.axvline(E1659, color="red", ls="--", label="E1659A")
    ax.plot(zn, zv, "k-", lw=0.5, alpha=0.5)
    ax.set_xlabel("Residue number")
    ax.set_ylabel("Displacement (Å)")
    ax.set_title("Zoom: E1659 neighborhood (1600-1720)")
    plt.colorbar(sc, ax=ax, label="pLDDT")
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel 4: zoomed K1141 region
    ax = axes[1, 1]
    zoom = [d for d in displacements if 1110 <= d["resnum"] <= 1180]
    if zoom:
        zn = [d["resnum"] for d in zoom]
        zv = [d["displacement_A"] for d in zoom]
        zp = [(d["WT_plddt"] + d["MUT_plddt"]) / 2 for d in zoom]
        sc = ax.scatter(zn, zv, c=zp, cmap="viridis", s=40, vmin=50, vmax=100)
        ax.axvline(K1141, color="blue", ls="--", label="K1141")
        ax.plot(zn, zv, "k-", lw=0.5, alpha=0.5)
        ax.set_xlabel("Residue number")
        ax.set_ylabel("Displacement (Å)")
        ax.set_title("Zoom: K1141 neighborhood (1110-1180)")
        plt.colorbar(sc, ax=ax, label="pLDDT")
        ax.legend()
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def main():
    wt_s = load(JOB_WT)
    mut_s = load(JOB_MUT)
    wt_dict = chain_residues(wt_s)
    mut_dict = chain_residues(mut_s)

    wt_hc = high_conf_residues(wt_dict, PLDDT_CUT)
    mut_hc = high_conf_residues(mut_dict, PLDDT_CUT)
    common_hc = wt_hc & mut_hc
    print(f"High-confidence (pLDDT≥{PLDDT_CUT}) common residues: {len(common_hc)} of {len(wt_dict)}")

    # Realign on high-confidence common residues
    rmsd_hc, _ = superimpose_on(wt_dict, mut_dict, common_hc)
    print(f"High-confidence RMSD: {rmsd_hc:.3f} Å")

    # Global RMSD over all shared residues (same transform now applied)
    all_pairs = []
    for n in wt_dict:
        if n in mut_dict and "CA" in wt_dict[n] and "CA" in mut_dict[n]:
            all_pairs.append((np.array(wt_dict[n]["CA"].coord),
                              np.array(mut_dict[n]["CA"].coord)))
    sq = [np.sum((a-b)**2) for a,b in all_pairs]
    rmsd_all = float(np.sqrt(np.mean(sq)))
    print(f"All-residue CA RMSD (same transform): {rmsd_all:.3f} Å")

    # Now meaningful local RMSDs
    wt_1659 = wt_dict[E1659]
    mut_1659 = mut_dict[E1659]
    wt_center = side_chain_centroid(wt_1659)
    loc_rmsd, loc_n, loc_pairs = local_rmsd(wt_dict, mut_dict, wt_center, HIGH_CONF_WINDOW_A)
    print(f"Local CA RMSD within {HIGH_CONF_WINDOW_A} Å of E{E1659}: {loc_rmsd:.3f} Å ({loc_n} res)")

    # K1141 displacement post-realignment
    wt_k = wt_dict[K1141]; mut_k = mut_dict[K1141]
    wt_k_nz = positive_tip(wt_k)
    mut_k_nz = positive_tip(mut_k)
    k_disp = float(np.linalg.norm(wt_k_nz - mut_k_nz)) if wt_k_nz is not None and mut_k_nz is not None else None
    print(f"K1141 NZ displacement (post-realign): {k_disp:.3f} Å")

    # Salt bridge in WT
    wt_e_tip = negative_tip(wt_1659)
    sb_wt = float(np.linalg.norm(wt_e_tip - wt_k_nz)) if wt_e_tip is not None and wt_k_nz is not None else None
    # And in MUT: distance K_NZ → CB of A1659 (stand-in for what USED TO BE E)
    mut_a_cb = np.array(mut_1659["CB"].coord) if "CB" in mut_1659 else None
    sb_mut = float(np.linalg.norm(mut_a_cb - mut_k_nz)) if mut_a_cb is not None and mut_k_nz is not None else None
    print(f"Salt-bridge-like distance E{E1659}(OE)→K{K1141}(NZ) WT: {sb_wt:.3f} Å")
    print(f"Distance A{E1659}(CB)→K{K1141}(NZ) MUT (post-realign): {sb_mut:.3f} Å")

    # Per-residue displacements
    print("\nComputing per-residue displacements...")
    all_disp = per_residue_displacement(wt_dict, mut_dict, range(1, 1776))

    # pLDDT-filtered local shell
    local_window = [d for d in all_disp if abs(d["resnum"] - E1659) <= 30]
    hot = sorted([d for d in local_window if d["WT_plddt"] >= 80 and d["MUT_plddt"] >= 80],
                 key=lambda x: -x["displacement_A"])[:10]
    print("\nTop residues displaced (pLDDT≥80, within ±30 of E1659):")
    for d in hot[:8]:
        print(f"  res{d['resnum']} ({d['WT_resname']}→{d['MUT_resname']}): {d['displacement_A']:.2f} Å "
              f"(pLDDT WT={d['WT_plddt']}, MUT={d['MUT_plddt']})")

    # Plot
    out_dir = Path("/Users/egorlyfar/Brain/research/strc/models")
    plot_phase1b(all_disp, loc_pairs, out_dir / "pharmacochaperone_phase1b.png")

    # Verdict
    if loc_rmsd < 1.0 and (k_disp or 0) < 1.5:
        verdict = "LOCAL STABILITY — mutation impact minimal → pharmacochaperone challenging (no hole to fill, no shift to correct)"
    elif loc_rmsd < 2.5 and (k_disp or 0) < 4:
        verdict = "MODERATE LOCAL REARRANGEMENT — classic pharmacochaperone niche (small-molecule reseats K1141 / rescues local fold)"
    else:
        verdict = "LARGE LOCAL MISFOLD — pharmacochaperone must engage WT fold intermediate; may need allosteric site"

    # JSON
    results = {
        "alignment_method": f"Kabsch on pLDDT≥{PLDDT_CUT} common residues (n={len(common_hc)})",
        "rmsd": {
            "high_confidence_A": round(rmsd_hc, 3),
            "all_residues_same_transform_A": round(rmsd_all, 3),
            "local_10A_shell_A": round(loc_rmsd, 3),
            "local_n_residues": loc_n,
        },
        "salt_bridge_geometry": {
            "wt_E1659_OE_to_K1141_NZ_A": round(sb_wt, 3) if sb_wt else None,
            "mut_A1659_CB_to_K1141_NZ_A": round(sb_mut, 3) if sb_mut else None,
            "K1141_NZ_displacement_WT_to_MUT_A": round(k_disp, 3) if k_disp else None,
        },
        "top_displaced_highconf": hot[:10],
        "verdict": verdict,
    }
    (out_dir / "pharmacochaperone_phase1b_results.json").write_text(json.dumps(results, indent=2))
    print("\nVerdict:", verdict)


if __name__ == "__main__":
    main()
