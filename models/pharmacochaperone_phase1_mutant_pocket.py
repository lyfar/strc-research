#!/usr/bin/env python3
"""
STRC Pharmacochaperone Phase 1 — Mutation Impact & Pocket Characterization

Compares AF3 wildtype (Job 4) vs E1659A mutant (Job 3) full-length STRC:

  1. Global structural alignment + RMSD
  2. Local RMSD within 10 Å shell of residue 1659
  3. K1141 side-chain displacement (salt-bridge partner)
  4. Acidic cluster analysis: E1655-E1659-E1164 triad
     (Ca²⁺ binding signature?)
  5. Pocket void volume in MUT — does E→A leave a drug-fillable cavity?
  6. Druggability metrics: volume, hydrophobic/polar ratio, H-bond donors/acceptors
  7. SiteMap-style scoring

Input:
  ~/Sites/site-strc-egor-lol/public/models/job3-mutant.cif
  ~/Sites/site-strc-egor-lol/public/models/job4-wildtype.cif

Output:
  pharmacochaperone_phase1_results.json
  pharmacochaperone_phase1.png
  pharmacochaperone_phase1_pocket_WT.pdb  (for PyMOL visualization)
  pharmacochaperone_phase1_pocket_MUT.pdb
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
from Bio.PDB import MMCIFParser, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------- Inputs -----------------------------
JOB_MUT = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job3-mutant.cif")
JOB_WT  = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job4-wildtype.cif")

E1659_RES = 1659
K1141_RES = 1141
E1655_RES = 1655
E1164_RES = 1164

LOCAL_SHELL_A = 10.0  # Å for local RMSD and pocket analysis
GRID_SPACING_A = 0.5   # Å for void volume calculation
PROBE_RADIUS_A = 1.4   # water probe
CARBON_VDW = 1.7       # simplified van der Waals

# Amino acid classes
POSITIVE_AA = {"LYS", "ARG", "HIS"}
NEGATIVE_AA = {"ASP", "GLU"}
HYDROPHOBIC_AA = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POLAR_AA = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
AROMATIC_AA = {"PHE", "TYR", "TRP", "HIS"}

# Backbone donor/acceptor atoms
HBOND_DONORS = {"N", "NZ", "NE", "NH1", "NH2", "NE1", "ND1", "ND2", "NE2",
                "OG", "OG1", "OH", "SG"}
HBOND_ACCEPTORS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OD", "OE"}


# ----------------------------- Loaders -----------------------------
def load(path: Path, label: str) -> Structure:
    parser = MMCIFParser(QUIET=True)
    s = parser.get_structure(label, str(path))
    return s


def get_chain_a(s: Structure):
    for model in s:
        return model["A"]
    raise ValueError("No chain A")


def residue_by_num(chain, num: int) -> Optional[Residue]:
    for res in chain:
        if res.id[1] == num and res.id[0] == " ":
            return res
    return None


def side_chain_centroid(res: Residue) -> np.ndarray:
    side = [a for a in res.get_atoms()
            if a.name not in ("N", "CA", "C", "O") and a.element != "H"]
    if not side:
        ca = res["CA"]
        return np.array(ca.coord)
    return np.mean([a.coord for a in side], axis=0)


def positive_tip(res: Residue) -> Optional[np.ndarray]:
    """NZ for LYS, centroid(NE,NH1,NH2) for ARG, centroid(ND1,NE2) for HIS."""
    atoms = {a.name: a for a in res.get_atoms()}
    if res.resname == "LYS" and "NZ" in atoms:
        return np.array(atoms["NZ"].coord)
    if res.resname == "ARG":
        tips = [atoms[n].coord for n in ("NE", "NH1", "NH2") if n in atoms]
        if tips: return np.mean(tips, axis=0)
    if res.resname == "HIS":
        tips = [atoms[n].coord for n in ("ND1", "NE2") if n in atoms]
        if tips: return np.mean(tips, axis=0)
    return None


def negative_tip(res: Residue) -> Optional[np.ndarray]:
    atoms = {a.name: a for a in res.get_atoms()}
    if res.resname == "GLU":
        tips = [atoms[n].coord for n in ("OE1", "OE2") if n in atoms]
        if tips: return np.mean(tips, axis=0)
    if res.resname == "ASP":
        tips = [atoms[n].coord for n in ("OD1", "OD2") if n in atoms]
        if tips: return np.mean(tips, axis=0)
    return None


# ----------------------------- Alignment + RMSD -----------------------------
def superimpose_on_distal(wt_chain, mut_chain, exclude_window=(1620, 1700)) -> Tuple[float, Superimposer]:
    """CA superimpose on residues OUTSIDE the exclude_window, then return global RMSD."""
    wt_ca, mut_ca = [], []
    for res_wt in wt_chain:
        if res_wt.id[0] != " ":
            continue
        n = res_wt.id[1]
        if exclude_window[0] <= n <= exclude_window[1]:
            continue
        res_mut = residue_by_num(mut_chain, n)
        if res_mut is None:
            continue
        if "CA" in res_wt and "CA" in res_mut:
            wt_ca.append(res_wt["CA"])
            mut_ca.append(res_mut["CA"])
    sup = Superimposer()
    sup.set_atoms(wt_ca, mut_ca)
    # Apply transform to MUT so it aligns onto WT
    all_atoms = [a for a in mut_chain.get_atoms()]
    sup.apply(all_atoms)
    return sup.rms, sup


def local_rmsd(wt_chain, mut_chain, center: np.ndarray, radius: float) -> Tuple[float, int]:
    """CA-CA RMSD within `radius` of center (center in WT frame after alignment)."""
    sq = []
    count = 0
    for res_wt in wt_chain:
        if res_wt.id[0] != " " or "CA" not in res_wt:
            continue
        ca_wt = np.array(res_wt["CA"].coord)
        if np.linalg.norm(ca_wt - center) > radius:
            continue
        res_mut = residue_by_num(mut_chain, res_wt.id[1])
        if res_mut is None or "CA" not in res_mut:
            continue
        ca_mut = np.array(res_mut["CA"].coord)
        sq.append(np.sum((ca_wt - ca_mut) ** 2))
        count += 1
    if not sq:
        return 0.0, 0
    return float(np.sqrt(np.mean(sq))), count


# ----------------------------- Pocket / Void -----------------------------
def all_heavy_atoms_within(chain, center: np.ndarray, radius: float) -> List[Atom]:
    atoms = []
    for res in chain:
        if res.id[0] != " ":
            continue
        for a in res.get_atoms():
            if a.element == "H":
                continue
            if np.linalg.norm(a.coord - center) <= radius:
                atoms.append(a)
    return atoms


def void_volume(chain, center: np.ndarray, radius: float,
                grid: float = GRID_SPACING_A, probe: float = PROBE_RADIUS_A) -> Dict:
    """Grid-based void detection — count grid points inside sphere(center, radius)
    that are >probe+vdw from any heavy atom (= accessible void)."""
    atoms = all_heavy_atoms_within(chain, center, radius + 3.0)  # pad for shell
    if not atoms:
        return {"void_A3": 0.0, "total_A3": 0.0, "fill_ratio": 0.0}
    coords = np.array([a.coord for a in atoms])

    # Build grid over sphere
    n = int(2 * radius / grid) + 1
    axis = np.linspace(-radius, radius, n)
    X, Y, Z = np.meshgrid(axis, axis, axis, indexing="ij")
    grid_pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1) + center
    in_sphere = np.linalg.norm(grid_pts - center, axis=1) <= radius
    total_pts = int(in_sphere.sum())
    interior = grid_pts[in_sphere]

    # Chunked distance to avoid huge matrices
    void_count = 0
    chunk = 5000
    threshold = probe + CARBON_VDW
    for i in range(0, len(interior), chunk):
        pts = interior[i:i+chunk]
        d = np.linalg.norm(pts[:, None, :] - coords[None, :, :], axis=2).min(axis=1)
        void_count += int((d > threshold).sum())

    vol_per_pt = grid ** 3
    return {
        "void_A3": round(void_count * vol_per_pt, 2),
        "total_A3": round(total_pts * vol_per_pt, 2),
        "fill_ratio": round(1 - void_count / total_pts, 3) if total_pts else 0.0,
    }


def pocket_composition(chain, center: np.ndarray, radius: float) -> Dict:
    """Classify residues and H-bond atoms in pocket."""
    residues_in_pocket = {}  # resnum -> (resname, min_dist)
    donor_count = acceptor_count = 0
    for res in chain:
        if res.id[0] != " ":
            continue
        min_d = 1e9
        for a in res.get_atoms():
            if a.element == "H":
                continue
            d = np.linalg.norm(a.coord - center)
            if d < min_d:
                min_d = d
            if d <= radius:
                if a.name in HBOND_DONORS:
                    donor_count += 1
                if a.name in HBOND_ACCEPTORS:
                    acceptor_count += 1
        if min_d <= radius:
            residues_in_pocket[res.id[1]] = (res.resname, float(min_d))

    hydrophobic = sum(1 for _, (rn, _) in residues_in_pocket.items() if rn in HYDROPHOBIC_AA)
    polar = sum(1 for _, (rn, _) in residues_in_pocket.items() if rn in POLAR_AA)
    positive = sum(1 for _, (rn, _) in residues_in_pocket.items() if rn in POSITIVE_AA)
    negative = sum(1 for _, (rn, _) in residues_in_pocket.items() if rn in NEGATIVE_AA)
    aromatic = sum(1 for _, (rn, _) in residues_in_pocket.items() if rn in AROMATIC_AA)
    total = len(residues_in_pocket)

    return {
        "total_residues": total,
        "hydrophobic": hydrophobic,
        "polar": polar,
        "positive": positive,
        "negative": negative,
        "aromatic": aromatic,
        "hydrophobic_ratio": round(hydrophobic / total, 3) if total else 0,
        "hbond_donors": donor_count,
        "hbond_acceptors": acceptor_count,
        "residues": sorted([(n, rn, round(d, 2)) for n, (rn, d) in residues_in_pocket.items()],
                           key=lambda x: x[2])[:30],
    }


def druggability_score(comp: Dict, void: Dict) -> Dict:
    """SiteMap-inspired heuristic: volume, hydrophobic ratio, H-bonds.
    Returns score in [0..1] and interpretation."""
    vol = void["void_A3"]
    # Ideal drug pocket: 200-700 Å³ (tiny fragments) or up to 1500 (lead-like)
    if vol < 50:
        vol_score = 0.1
    elif vol < 150:
        vol_score = 0.4
    elif vol < 600:
        vol_score = 1.0
    elif vol < 1200:
        vol_score = 0.8
    else:
        vol_score = 0.5

    hp = comp["hydrophobic_ratio"]
    # Ideal: 0.5-0.7 hydrophobic
    hp_score = 1.0 - abs(hp - 0.6) / 0.6

    # H-bond balance: want some of each
    hba = comp["hbond_acceptors"]
    hbd = comp["hbond_donors"]
    hb_score = min(hba, hbd) / max(hba, hbd, 1)

    # Composite (volume weighted heaviest)
    score = 0.5 * vol_score + 0.3 * hp_score + 0.2 * hb_score

    if score >= 0.7:
        verdict = "HIGHLY DRUGGABLE"
    elif score >= 0.5:
        verdict = "DRUGGABLE"
    elif score >= 0.35:
        verdict = "CHALLENGING"
    else:
        verdict = "UNDRUGGABLE"

    return {
        "score": round(score, 3),
        "volume_score": round(vol_score, 3),
        "hydrophobic_score": round(hp_score, 3),
        "hbond_balance_score": round(hb_score, 3),
        "verdict": verdict,
    }


# ----------------------------- PDB export -----------------------------
def export_pocket_pdb(chain, center: np.ndarray, radius: float, out: Path, title: str):
    """Minimal PDB writer for pocket residues — used by PyMOL/ChimeraX."""
    lines = [f"REMARK  {title}",
             f"REMARK  Center: {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}",
             f"REMARK  Radius: {radius} A"]
    atom_id = 1
    for res in chain:
        if res.id[0] != " ":
            continue
        min_d = min(np.linalg.norm(np.array(a.coord) - center)
                    for a in res.get_atoms() if a.element != "H")
        if min_d > radius:
            continue
        for a in res.get_atoms():
            if a.element == "H":
                continue
            x, y, z = a.coord
            lines.append(
                f"ATOM  {atom_id:5d} {a.name:<4s} {res.resname} A{res.id[1]:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {a.element:>2s}"
            )
            atom_id += 1
    lines.append("END")
    out.write_text("\n".join(lines))


# ----------------------------- Plotting -----------------------------
def plot_results(wt_chain, mut_chain, wt_center, mut_center, out: Path):
    fig = plt.figure(figsize=(14, 10))

    # Panel 1: 3D view of pocket, WT and MUT overlaid
    ax1 = fig.add_subplot(2, 2, 1, projection="3d")
    _scatter_pocket(ax1, wt_chain, wt_center, LOCAL_SHELL_A, color="steelblue", label="WT")
    _scatter_pocket(ax1, mut_chain, mut_center, LOCAL_SHELL_A, color="firebrick", label="MUT", alpha=0.5)
    ax1.set_title("Pocket around residue 1659 (WT blue vs MUT red)")
    ax1.legend()

    # Panel 2: residue-by-residue shift map
    ax2 = fig.add_subplot(2, 2, 2)
    shifts = []
    nums = []
    for res_wt in wt_chain:
        if res_wt.id[0] != " " or "CA" not in res_wt:
            continue
        n = res_wt.id[1]
        if not (1600 <= n <= 1720):
            continue
        res_mut = residue_by_num(mut_chain, n)
        if res_mut is None or "CA" not in res_mut:
            continue
        d = np.linalg.norm(np.array(res_wt["CA"].coord) - np.array(res_mut["CA"].coord))
        nums.append(n)
        shifts.append(d)
    ax2.plot(nums, shifts, "k-", lw=0.8)
    ax2.axvline(1659, color="red", ls="--", label="E1659A")
    ax2.axvline(1655, color="orange", ls=":", alpha=0.7, label="E1655")
    ax2.axvline(1660, color="gray", ls=":", alpha=0.5)
    ax2.set_xlabel("Residue number")
    ax2.set_ylabel("Cα displacement (Å)")
    ax2.set_title("Backbone shift around E1659A mutation")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    # Panel 3: contact map change
    ax3 = fig.add_subplot(2, 2, 3)
    win = range(1640, 1680)
    wt_map = _distance_map(wt_chain, win)
    mut_map = _distance_map(mut_chain, win)
    diff = mut_map - wt_map
    im = ax3.imshow(diff, cmap="RdBu_r", vmin=-4, vmax=4,
                    extent=[win[0], win[-1], win[-1], win[0]])
    ax3.set_xlabel("Residue")
    ax3.set_ylabel("Residue")
    ax3.set_title("Contact-map delta (MUT − WT, CA-CA Å)")
    plt.colorbar(im, ax=ax3, shrink=0.8)

    # Panel 4: K1141 distance to E1659-(would-be) position
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis("off")
    ax4.text(0.02, 0.98, "Salt-bridge partner (K1141) geometry",
             fontsize=12, fontweight="bold", va="top", transform=ax4.transAxes)
    ax4.text(0.02, 0.90, "— see pocket_summary in JSON —",
             fontsize=10, va="top", transform=ax4.transAxes, family="monospace")

    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def _scatter_pocket(ax, chain, center, radius, color, label, alpha=1.0):
    xs, ys, zs = [], [], []
    for res in chain:
        if res.id[0] != " " or "CA" not in res:
            continue
        c = np.array(res["CA"].coord)
        if np.linalg.norm(c - center) <= radius:
            xs.append(c[0]); ys.append(c[1]); zs.append(c[2])
    ax.scatter(xs, ys, zs, c=color, s=40, alpha=alpha, label=label)


def _distance_map(chain, residue_range) -> np.ndarray:
    nums = list(residue_range)
    n = len(nums)
    m = np.zeros((n, n))
    for i, ni in enumerate(nums):
        ri = residue_by_num(chain, ni)
        if ri is None or "CA" not in ri: continue
        ci = np.array(ri["CA"].coord)
        for j, nj in enumerate(nums):
            rj = residue_by_num(chain, nj)
            if rj is None or "CA" not in rj: continue
            cj = np.array(rj["CA"].coord)
            m[i, j] = np.linalg.norm(ci - cj)
    return m


# ----------------------------- Main -----------------------------
def main():
    print("Loading structures...")
    wt = load(JOB_WT, "WT")
    mut = load(JOB_MUT, "MUT")
    wt_chain = get_chain_a(wt)
    mut_chain = get_chain_a(mut)

    # Verify mutation
    wt_1659 = residue_by_num(wt_chain, E1659_RES)
    mut_1659 = residue_by_num(mut_chain, E1659_RES)
    print(f"WT residue 1659: {wt_1659.resname}")
    print(f"MUT residue 1659: {mut_1659.resname}")
    assert wt_1659.resname == "GLU", f"Expected GLU in WT, got {wt_1659.resname}"
    assert mut_1659.resname == "ALA", f"Expected ALA in MUT, got {mut_1659.resname}"

    # Global alignment on residues outside the mutation neighborhood
    print("\nSuperimposing on distal backbone...")
    global_rmsd, _ = superimpose_on_distal(wt_chain, mut_chain, exclude_window=(1620, 1700))
    print(f"Global CA RMSD (excluding 1620-1700): {global_rmsd:.3f} Å")

    # E1659 side-chain centroids after alignment
    wt_center = side_chain_centroid(wt_1659)
    # In mutant, use CB of Ala + projected side chain direction for center
    mut_cb = mut_1659["CB"] if "CB" in mut_1659 else None
    mut_center = np.array(mut_cb.coord) if mut_cb else side_chain_centroid(mut_1659)

    # Local RMSD
    print("\nLocal RMSD around residue 1659...")
    loc_rmsd_wt, loc_n_wt = local_rmsd(wt_chain, mut_chain, wt_center, LOCAL_SHELL_A)
    print(f"Local CA RMSD within {LOCAL_SHELL_A} Å: {loc_rmsd_wt:.3f} Å ({loc_n_wt} residues)")

    # K1141 position
    wt_k1141 = residue_by_num(wt_chain, K1141_RES)
    mut_k1141 = residue_by_num(mut_chain, K1141_RES)
    wt_k_tip = positive_tip(wt_k1141)
    mut_k_tip = positive_tip(mut_k1141)
    wt_e_tip = negative_tip(wt_1659)
    k_shift = float(np.linalg.norm(wt_k_tip - mut_k_tip)) if (wt_k_tip is not None and mut_k_tip is not None) else None
    salt_bridge_wt_dist = float(np.linalg.norm(wt_e_tip - wt_k_tip)) if (wt_e_tip is not None and wt_k_tip is not None) else None
    print(f"K1141 NZ displacement WT→MUT: {k_shift:.3f} Å" if k_shift else "K1141 tip missing")
    print(f"Salt bridge E1659-K1141 in WT: {salt_bridge_wt_dist:.3f} Å" if salt_bridge_wt_dist else "")

    # Acidic triad
    wt_e1655 = negative_tip(residue_by_num(wt_chain, E1655_RES))
    wt_e1164 = negative_tip(residue_by_num(wt_chain, E1164_RES))
    wt_e1659_tip = wt_e_tip
    triad = {}
    if wt_e1655 is not None and wt_e1659_tip is not None:
        triad["E1655_to_E1659_A"] = float(np.linalg.norm(wt_e1655 - wt_e1659_tip))
    if wt_e1164 is not None and wt_e1659_tip is not None:
        triad["E1164_to_E1659_A"] = float(np.linalg.norm(wt_e1164 - wt_e1659_tip))
    if wt_e1164 is not None and wt_e1655 is not None:
        triad["E1655_to_E1164_A"] = float(np.linalg.norm(wt_e1164 - wt_e1655))
    # Centroid of triad
    if wt_e1164 is not None and wt_e1655 is not None and wt_e1659_tip is not None:
        triad_cen = (wt_e1164 + wt_e1655 + wt_e1659_tip) / 3
        triad["centroid_A"] = triad_cen.tolist()
        triad["to_K1141_NZ"] = float(np.linalg.norm(triad_cen - wt_k_tip)) if wt_k_tip is not None else None
    print(f"\nAcidic triad (WT): {triad}")

    # Void volume in MUT at former E1659 position
    print("\nVoid volume analysis...")
    # Use the wt side-chain centroid as the reference point in the ALIGNED mut frame.
    # After superimpose_on_distal, MUT has been moved onto WT — so wt_center is valid reference.
    void_mut = void_volume(mut_chain, wt_center, radius=5.5)
    void_wt  = void_volume(wt_chain,  wt_center, radius=5.5)
    print(f"Void (MUT): {void_mut}")
    print(f"Void (WT):  {void_wt}")
    delta_void = void_mut["void_A3"] - void_wt["void_A3"]
    print(f"Δvoid (MUT − WT): {delta_void:.1f} Å³  (= hole from E→A)")

    # Pocket composition (docking-relevant shell)
    print("\nPocket composition (8 Å shell)...")
    comp_mut = pocket_composition(mut_chain, wt_center, radius=8.0)
    comp_wt  = pocket_composition(wt_chain, wt_center, radius=8.0)
    drug_mut = druggability_score(comp_mut, void_mut)
    drug_wt  = druggability_score(comp_wt, void_wt)
    print(f"MUT druggability: {drug_mut}")
    print(f"WT  druggability: {drug_wt}")

    # Export pocket PDBs for visualization
    print("\nExporting pocket PDBs...")
    out_dir = Path("/Users/egorlyfar/Brain/research/strc/models")
    export_pocket_pdb(wt_chain, wt_center, LOCAL_SHELL_A,
                      out_dir / "pharmacochaperone_phase1_pocket_WT.pdb", "WT pocket")
    export_pocket_pdb(mut_chain, wt_center, LOCAL_SHELL_A,
                      out_dir / "pharmacochaperone_phase1_pocket_MUT.pdb", "MUT pocket (aligned)")

    # Plot
    print("Plotting...")
    plot_results(wt_chain, mut_chain, wt_center, wt_center,
                 out_dir / "pharmacochaperone_phase1.png")

    # JSON output
    results = {
        "inputs": {
            "wildtype_cif": str(JOB_WT),
            "mutant_cif": str(JOB_MUT),
        },
        "mutation_verified": {"WT": wt_1659.resname, "MUT": mut_1659.resname},
        "alignment": {
            "global_CA_rmsd_excluding_1620_1700_A": round(global_rmsd, 4),
            "local_CA_rmsd_10A_shell_around_1659_A": round(loc_rmsd_wt, 4),
            "n_local_residues": loc_n_wt,
        },
        "salt_bridge": {
            "E1659_to_K1141_WT_tip_A": round(salt_bridge_wt_dist, 3) if salt_bridge_wt_dist else None,
            "K1141_NZ_displacement_WT_to_MUT_A": round(k_shift, 3) if k_shift else None,
            "interpretation": (
                "Intact WT salt bridge; mutation breaks it and K1141 shifts"
                if salt_bridge_wt_dist and salt_bridge_wt_dist < 4.0 and k_shift and k_shift > 1.0
                else "No dramatic rearrangement detected"
            ),
        },
        "acidic_triad": {
            "E1655-E1659-E1164_distances_A": {k: round(v, 3) if isinstance(v, float) else v
                                                for k, v in triad.items()},
            "interpretation": (
                "Ca²⁺-binding geometry (EEE triad ≤8 Å)"
                if all(isinstance(triad.get(k), float) and triad[k] < 8
                       for k in ("E1655_to_E1659_A", "E1164_to_E1659_A"))
                else "Not a tight Ca²⁺ site (distances too large)"
            ),
        },
        "void_analysis": {
            "mutant_void_A3": void_mut["void_A3"],
            "wildtype_void_A3": void_wt["void_A3"],
            "delta_void_A3_from_E_to_A": round(delta_void, 2),
            "shell_radius_A": 5.5,
            "grid_spacing_A": GRID_SPACING_A,
        },
        "pocket_composition_8A": {
            "mutant": comp_mut,
            "wildtype": comp_wt,
        },
        "druggability": {
            "mutant": drug_mut,
            "wildtype": drug_wt,
        },
        "verdict": _compose_verdict(drug_mut, drug_wt, delta_void, salt_bridge_wt_dist, k_shift, triad),
    }

    (out_dir / "pharmacochaperone_phase1_results.json").write_text(json.dumps(results, indent=2))
    print(f"\nSaved → pharmacochaperone_phase1_results.json")
    print("\nVerdict:", results["verdict"])


def _compose_verdict(drug_mut, drug_wt, delta_void, sb_dist, k_shift, triad) -> str:
    parts = []
    parts.append(f"Phase 1 pharmacochaperone assessment:")
    if drug_mut["score"] >= 0.5:
        parts.append(f"  MUT pocket {drug_mut['verdict']} (score {drug_mut['score']})")
    else:
        parts.append(f"  MUT pocket {drug_mut['verdict']} (score {drug_mut['score']}) — virtual screen may be low-yield")
    if delta_void > 20:
        parts.append(f"  Visible void from E→A: +{delta_void:.0f} Å³ — drug can fill this")
    elif delta_void > 5:
        parts.append(f"  Modest void: +{delta_void:.1f} Å³")
    else:
        parts.append(f"  No void expansion from mutation (Δ{delta_void:+.1f} Å³) — MUT may be compensated by sidechain flex")
    if sb_dist and sb_dist < 5:
        parts.append(f"  WT salt bridge E1659-K1141: {sb_dist:.2f} Å (YES, intact)")
        if k_shift and k_shift > 2:
            parts.append(f"  K1141 shifts {k_shift:.2f} Å in MUT — salt bridge broken")
    if "E1655_to_E1659_A" in triad and triad["E1655_to_E1659_A"] < 10 and \
       "E1164_to_E1659_A" in triad and triad["E1164_to_E1659_A"] < 12:
        parts.append(f"  Acidic cluster (E1655-E1659-E1164) within ~10 Å — POSSIBLE Ca²⁺ binding site")
        parts.append(f"  Suggests pharmacochaperone strategy: Ca²⁺-mimetic small molecule or positive-charge anchor")
    return " | ".join(parts)


if __name__ == "__main__":
    main()
