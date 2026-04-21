#!/usr/bin/env python3
"""
STRC Pharmacochaperone Virtual Screen E1659A — Phase 0 Target Prep

Extracts docking-ready target from AF3 Job 2 CIF:
  - Locate E1659 in STRC chain (within mini-STRC numbering)
  - Find TMEM145 counterion residue (positively charged, closest)
  - Compute grid box center & dimensions
  - Enumerate binding-site residues (4 and 6 Å shells)
  - Emit JSON spec for AutoDock Vina / Smina

Input:  ~/Sites/site-strc-egor-lol/public/models/job2-mini-complex.cif
Output: pharmacochaperone_target_prep.json
        pharmacochaperone_target_prep.png (cavity geometry figure)
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


CIF = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job2-mini-complex.cif")

# Mini-STRC in AF3 Job 2 spans residues 594-1775 of full STRC
# Chain A numbering: 1..1182 maps to full 594..1775
FULL_TO_CHAIN_OFFSET = 594 - 1   # full = chain_num + offset
FULL_E1659 = 1659                # the mutation of interest
CHAIN_A_E1659 = FULL_E1659 - FULL_TO_CHAIN_OFFSET  # = 1066

# Docking grid box
GRID_BOX_A = 18.0     # 18 × 18 × 18 Å (standard for small molecule)
POCKET_RADIUS_A = 6.0 # residues within this radius are the binding pocket

# Positively charged amino acids (salt-bridge partners)
POSITIVE_AA = {"LYS", "ARG", "HIS"}
NEGATIVE_AA = {"ASP", "GLU"}


def load_structure(path: Path) -> Structure:
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure("job2", str(path))


def residue_centroid(residue: Residue) -> np.ndarray:
    coords = [atom.coord for atom in residue.get_atoms()]
    return np.mean(coords, axis=0)


def residue_side_chain_centroid(residue: Residue) -> np.ndarray:
    """Centroid of side-chain heavy atoms (fallback to full residue for GLY)."""
    side_atoms = [a for a in residue.get_atoms()
                  if a.name not in ("N", "CA", "C", "O") and not a.name.startswith("H")]
    if not side_atoms:
        return residue_centroid(residue)
    return np.mean([a.coord for a in side_atoms], axis=0)


def find_E1659(structure: Structure) -> Residue:
    for model in structure:
        chain_a = model["A"]
        # Chain A residue = CHAIN_A_E1659 (= 1066) in mini-STRC numbering
        for res in chain_a:
            if res.id[1] == CHAIN_A_E1659:
                if res.resname not in ("GLU", "GLN"):
                    raise ValueError(f"Expected GLU at mini-STRC position {CHAIN_A_E1659} (full 1659), "
                                     f"got {res.resname} — numbering mismatch?")
                return res
    raise ValueError("E1659 not found in AF3 Job 2 CIF")


def interface_residues(structure: Structure, center: np.ndarray, radius: float) -> List[Dict]:
    """All residues with any heavy atom within `radius` of `center`."""
    model = structure[0]
    hits = []
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            full_num = res.id[1] + FULL_TO_CHAIN_OFFSET if chain.id == "A" else res.id[1]
            for atom in res.get_atoms():
                if atom.element == "H":
                    continue
                d = np.linalg.norm(atom.coord - center)
                if d <= radius:
                    hits.append({
                        "chain": chain.id,
                        "chain_resnum": int(res.id[1]),
                        "full_resnum": int(full_num) if chain.id == "A" else int(res.id[1]),
                        "resname": res.resname,
                        "min_atom_dist_A": float(d),
                        "atom": atom.name,
                    })
                    break  # only record closest atom per residue
    return sorted(hits, key=lambda r: r["min_atom_dist_A"])


def find_closest_positive(structure: Structure, e1659: Residue,
                          chain_filter: str | None = None) -> Dict | None:
    """
    Find closest positively-charged residue (any chain, or filter by chain)
    to the E1659 carboxylate tip (OE1/OE2 average).
    Measures side-chain tip distance (OE → NZ / NH / NE).
    """
    # E1659 carboxylate tip: OE1, OE2 centroid
    oe_atoms = [a for a in e1659 if a.name in ("OE1", "OE2")]
    if not oe_atoms:
        return None
    oe_centroid = np.mean([a.coord for a in oe_atoms], axis=0)

    model = structure[0]
    best = None
    best_dist = float("inf")
    for chain in model:
        if chain_filter and chain.id != chain_filter:
            continue
        for res in chain:
            if res.resname not in POSITIVE_AA:
                continue
            # Tip atoms for salt-bridge
            if res.resname == "LYS":
                tip_atoms = [a for a in res if a.name == "NZ"]
            elif res.resname == "ARG":
                tip_atoms = [a for a in res if a.name in ("NH1", "NH2", "NE")]
            else:  # HIS
                tip_atoms = [a for a in res if a.name in ("ND1", "NE2")]
            if not tip_atoms:
                continue
            tip_centroid = np.mean([a.coord for a in tip_atoms], axis=0)
            d = np.linalg.norm(tip_centroid - oe_centroid)
            if d < best_dist:
                best_dist = d
                full_num = res.id[1] + FULL_TO_CHAIN_OFFSET if chain.id == "A" else res.id[1]
                best = {
                    "chain": chain.id,
                    "chain_resnum": int(res.id[1]),
                    "full_resnum": int(full_num) if chain.id == "A" else int(res.id[1]),
                    "resname": res.resname,
                    "tip_distance_A": float(d),
                    "tip_centroid_coord": tip_centroid.tolist(),
                }
    return best


def compute_cavity_axis(e1659_cb: np.ndarray, counterion_centroid: np.ndarray) -> Dict:
    vec = counterion_centroid - e1659_cb
    dist = np.linalg.norm(vec)
    unit = vec / dist
    return {
        "distance_A": float(dist),
        "unit_vector": unit.tolist(),
    }


# ============================================================
# Plot: 3D view of cavity + counterion axis + pocket residues
# ============================================================

def plot_cavity(structure: Structure, e1659: Residue, counterion_any: Dict | None,
                counterion_tmem: Dict | None, pocket_4A: List[Dict],
                pocket_6A: List[Dict], outpath: Path):
    model = structure[0]
    # Gather coordinates
    e1659_cb = np.array([a.coord for a in e1659 if a.name == "CB"][0])
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Plot all STRC + TMEM145 CA atoms as background
    for chain_id, color, label in [("A", "C0", "STRC (mini)"), ("B", "C1", "TMEM145")]:
        ca_coords = []
        for res in model[chain_id]:
            if res.id[0] != " ":
                continue
            if "CA" in res:
                ca_coords.append(res["CA"].coord)
        ca_coords = np.array(ca_coords)
        ax.scatter(ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2],
                   s=2, alpha=0.2, color=color, label=label)

    # Mark E1659 Cβ
    ax.scatter(*e1659_cb, s=200, color="red", edgecolor="black", label="E1659 Cβ (mutation site)")

    # Mark primary counterion (intra or inter chain, whichever is closest)
    if counterion_any:
        sc = np.array(counterion_any["tip_centroid_coord"])
        location = "STRC" if counterion_any["chain"] == "A" else "TMEM145"
        label_num = counterion_any.get("full_resnum", counterion_any["chain_resnum"])
        ax.scatter(*sc, s=200, color="blue", edgecolor="black",
                   label=f"{location} {counterion_any['resname']}{label_num} "
                         f"(tip-tip {counterion_any['tip_distance_A']:.1f} Å)")
        ax.plot(*zip(e1659_cb, sc), color="black", linewidth=2, alpha=0.6)

    # Mark TMEM145 closest too if different
    if counterion_tmem and (not counterion_any or counterion_any["chain"] != "B"):
        sc_t = np.array(counterion_tmem["tip_centroid_coord"])
        ax.scatter(*sc_t, s=120, color="orange", edgecolor="black", marker="s",
                   label=f"TMEM145 {counterion_tmem['resname']}{counterion_tmem['chain_resnum']} "
                         f"(too far: {counterion_tmem['tip_distance_A']:.1f} Å)")

    # Shade pocket residues
    for res_info in pocket_6A[:25]:  # top 25 closest
        chain_id = res_info["chain"]
        rn = res_info["chain_resnum"]
        try:
            res = model[chain_id][rn]
            centroid = residue_centroid(res)
            d = res_info["min_atom_dist_A"]
            color = "purple" if d <= 4.0 else "orange"
            ax.scatter(*centroid, s=30, color=color, alpha=0.6)
        except KeyError:
            continue

    # Draw the docking grid box (18 × 18 × 18 Å centered at E1659 Cβ)
    cx, cy, cz = e1659_cb
    h = GRID_BOX_A / 2.0
    verts = np.array([[cx + dx, cy + dy, cz + dz]
                      for dx in (-h, h) for dy in (-h, h) for dz in (-h, h)])
    edges = [
        (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
        (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)
    ]
    for a, b in edges:
        ax.plot(*zip(verts[a], verts[b]), color="green", alpha=0.3, linewidth=1)
    ax.plot([], [], color="green", alpha=0.5, label=f"Docking box {GRID_BOX_A:.0f}×{GRID_BOX_A:.0f}×{GRID_BOX_A:.0f} Å")

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.set_title("STRC E1659 cavity for pharmacochaperone docking (AF3 Job 2)")
    ax.legend(loc="upper left", fontsize=7)

    # Zoom to cavity
    zoom = 15
    ax.set_xlim(cx - zoom, cx + zoom)
    ax.set_ylim(cy - zoom, cy + zoom)
    ax.set_zlim(cz - zoom, cz + zoom)

    plt.tight_layout()
    plt.savefig(outpath, dpi=130)
    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    out_dir = Path(__file__).parent
    print(f"Loading AF3 Job 2 CIF: {CIF.name}")
    structure = load_structure(CIF)
    e1659 = find_E1659(structure)
    e1659_cb = np.array([a.coord for a in e1659 if a.name == "CB"][0])
    print(f"Found E1659 at mini-STRC residue {CHAIN_A_E1659} ({e1659.resname})")
    print(f"E1659 Cβ coords: ({e1659_cb[0]:.2f}, {e1659_cb[1]:.2f}, {e1659_cb[2]:.2f})")

    # Find closest positive residue (anywhere)
    counterion_any = find_closest_positive(structure, e1659, chain_filter=None)
    counterion_strc = find_closest_positive(structure, e1659, chain_filter="A")
    counterion_tmem = find_closest_positive(structure, e1659, chain_filter="B")

    print()
    print("Closest salt-bridge partners (carboxylate tip → basic tip):")
    if counterion_any:
        src = "intra-STRC" if counterion_any["chain"] == "A" else "TMEM145"
        print(f"  OVERALL closest: {src} {counterion_any['resname']}"
              f"{counterion_any.get('full_resnum', counterion_any['chain_resnum'])} at {counterion_any['tip_distance_A']:.2f} Å")
    if counterion_strc:
        print(f"  Intra-STRC closest: {counterion_strc['resname']}{counterion_strc['full_resnum']} "
              f"(chain A {counterion_strc['chain_resnum']}) at {counterion_strc['tip_distance_A']:.2f} Å")
    if counterion_tmem:
        print(f"  TMEM145 closest: {counterion_tmem['resname']}{counterion_tmem['chain_resnum']} at "
              f"{counterion_tmem['tip_distance_A']:.2f} Å")

    pocket_4A = interface_residues(structure, e1659_cb, 4.0)
    pocket_6A = interface_residues(structure, e1659_cb, 6.0)
    pocket_10A = interface_residues(structure, e1659_cb, 10.0)
    print(f"Pocket residues (within 4 Å of E1659 Cβ): {len(pocket_4A)}")
    print(f"Pocket residues (within 6 Å of E1659 Cβ): {len(pocket_6A)}")
    print(f"Pocket residues (within 10 Å): {len(pocket_10A)}")

    # Focused list for docking sidebar
    strc_pocket = [r for r in pocket_6A if r["chain"] == "A"]
    tmem_pocket = [r for r in pocket_6A if r["chain"] == "B"]

    # Check for other charged residues nearby
    negative_in_pocket = [r for r in pocket_10A if r["resname"] in NEGATIVE_AA]
    positive_in_pocket = [r for r in pocket_10A if r["resname"] in POSITIVE_AA]

    result = {
        "source_cif": str(CIF),
        "mini_strc_chain": "A",
        "tmem145_chain": "B",
        "full_to_chain_offset": FULL_TO_CHAIN_OFFSET,
        "e1659": {
            "full_resnum": FULL_E1659,
            "chain_a_resnum": CHAIN_A_E1659,
            "resname_in_structure": e1659.resname,
            "Cbeta_coords": e1659_cb.tolist(),
            "sidechain_centroid": residue_side_chain_centroid(e1659).tolist(),
        },
        "closest_positive_any": counterion_any,
        "closest_intra_strc_positive": counterion_strc,
        "closest_tmem145_positive": counterion_tmem,
        "cavity_axis_to_primary_counterion": (
            compute_cavity_axis(e1659_cb, np.array(counterion_any["tip_centroid_coord"]))
            if counterion_any else None
        ),
        "grid_box": {
            "center_A": e1659_cb.tolist(),
            "dimensions_A": [GRID_BOX_A, GRID_BOX_A, GRID_BOX_A],
            "description": "Centered on E1659 Cβ; long axis aligned with TMEM145 counterion in post-processing."
        },
        "pocket_residues_4A": pocket_4A,
        "pocket_residues_6A": pocket_6A,
        "strc_pocket_residues": strc_pocket,
        "tmem145_pocket_residues": tmem_pocket,
        "nearby_charge_landscape_10A": {
            "negative_residues": negative_in_pocket,
            "positive_residues": positive_in_pocket,
        },
        "docking_spec": {
            "tool": "AutoDock Vina or Smina",
            "exhaustiveness": 32,
            "num_poses": 20,
            "grid_center_A": e1659_cb.tolist(),
            "grid_size_A": [GRID_BOX_A, GRID_BOX_A, GRID_BOX_A],
            "pharmacophore": {
                "required_features": [
                    "anionic anchor (carboxylate, phosphonate, or sulfonate) within 6 Å of TMEM145 counterion",
                    "cavity-filling body (40-50 Å³) occupying void left by E1659 Ala substitution",
                    "≥1 H-bond acceptor matching TMEM145 backbone amide donors",
                    "MW 250-400 Da, LogP 1.5-3.5, HBA 3-5, HBD 0-2 (RWM crossing envelope)",
                ]
            }
        },
    }

    json_path = out_dir / "pharmacochaperone_target_prep.json"
    png_path = out_dir / "pharmacochaperone_target_prep.png"

    def default(o):
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, (np.floating,)): return float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        return str(o)

    with json_path.open("w") as f:
        json.dump(result, f, indent=2, default=default)

    plot_cavity(structure, e1659, counterion_any, counterion_tmem, pocket_4A, pocket_6A, png_path)

    print()
    print("=" * 78)
    print("Pharmacochaperone Target Prep — AF3 Job 2 CIF Analysis")
    print("=" * 78)
    print(f"STRC E1659 (mini-STRC chain A residue {CHAIN_A_E1659})")
    print(f"Cβ coordinates: ({e1659_cb[0]:.2f}, {e1659_cb[1]:.2f}, {e1659_cb[2]:.2f}) Å")
    print(f"Docking grid box: {GRID_BOX_A} × {GRID_BOX_A} × {GRID_BOX_A} Å centered on Cβ")
    print()
    if counterion_any:
        src = "intra-STRC" if counterion_any["chain"] == "A" else "TMEM145"
        label_num = counterion_any.get("full_resnum", counterion_any["chain_resnum"])
        print(f"Primary counterion: {src} {counterion_any['resname']}{label_num} at {counterion_any['tip_distance_A']:.2f} Å (tip-tip)")
        print(f"  → hypothesis {'STABILIZED-FOLD (intra-STRC)' if counterion_any['chain'] == 'A' else 'INTERFACE (STRC-TMEM145)'} salt bridge")
    print()
    print(f"STRC pocket residues (within 6 Å of E1659 Cβ): {len(strc_pocket)}")
    for r in strc_pocket[:10]:
        print(f"  - {r['resname']}{r['full_resnum']} (chain A {r['chain_resnum']}): {r['min_atom_dist_A']:.2f} Å")
    print()
    print(f"TMEM145 pocket residues (within 6 Å of E1659 Cβ): {len(tmem_pocket)}")
    for r in tmem_pocket[:10]:
        print(f"  - {r['resname']}{r['tmem145_resnum']}: {r['min_atom_dist_A']:.2f} Å")
    print()
    print(f"Nearby positive residues (10 Å, potential salt-bridge partners): "
          f"{len(positive_in_pocket)}")
    for r in positive_in_pocket[:5]:
        print(f"  - {r['resname']} chain {r['chain']} res {r['full_resnum']}: {r['min_atom_dist_A']:.2f} Å")
    print()
    print(f"Written: {json_path.name}, {png_path.name}")


if __name__ == "__main__":
    main()
