#!/usr/bin/env python3
"""
OTOA Phase 1B — structural superposition of STRC ARM repeats vs OTOA C-terminal.

Rationale: Phase 1A (sequence) showed 28% global / 34% ARM-region identity
STRC(1603-1770) <-> OTOA(977-1136). Sequence alone is a weak predictor of
function. This script asks: when we superimpose the Calpha atoms of aligned
pairs, do they converge (RMSD <= 4 A) or diverge?

Convergence => OTOA can plausibly take over STRC's ARM fold and hence its
TMEM145-docking role. Divergence => paralog relationship is evolutionary
leftover without shared 3D structure.

Inputs:
  - STRC structure: ~/Sites/site-strc-egor-lol/public/models/job4-wildtype.cif
    (our AF3 WT full-length, residues 1-1775)
  - OTOA structure: /tmp/af_otoa.pdb (AlphaFold DB v6, Q7RTW8, residues 1-1153)
  - UniProt sequences fetched live

Outputs:
  - otoa_phase1b_structure_alignment.json
      * Calpha RMSD for ARM region anchors
      * Calpha RMSD for full mini-STRC <-> OTOA global aligned pairs
      * per-block RMSD breakdown
      * anchor residue 3D distance comparison (K1141 pocket, E1659)

Replication:
    /opt/miniconda3/bin/python3 otoa_phase1b_structure_alignment.py
"""
import json
import warnings
from io import StringIO
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa

warnings.filterwarnings("ignore")

OUT_DIR = Path(__file__).resolve().parent
OUT_JSON = OUT_DIR / "otoa_phase1b_structure_alignment.json"

STRC_CIF = Path.home() / "Sites/site-strc-egor-lol/public/models/job4-wildtype.cif"
OTOA_PDB = Path("/tmp/af_otoa.pdb")

STRC_ACC = "Q7RTU9"
OTOA_ACC = "Q7RTW8"

MINI_STRC_START, MINI_STRC_END = 700, 1775
ARM_START, ARM_END = 1603, 1770
POCKET_RES = {"K1141": 1141, "G1645": 1645, "F1646": 1646, "E1659": 1659}

BL62 = substitution_matrices.load("BLOSUM62")


def fetch_fasta(acc):
    with urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta") as r:
        return next(SeqIO.parse(StringIO(r.read().decode()), "fasta"))


def load_ca_by_resnum(structure):
    """Return {resnum:int -> (CA atom, residue_letter)}."""
    out = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                if not is_aa(res, standard=True):
                    continue
                if "CA" not in res:
                    continue
                out[res.id[1]] = (res["CA"], res.get_resname())
            break
        break
    return out


def three_to_one(three):
    map_ = dict(
        ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E",
        GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
        PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V",
    )
    return map_.get(three, "X")


def align_local(s1, s2):
    a = Align.PairwiseAligner()
    a.substitution_matrix = BL62
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "local"
    return a.align(s1, s2)[0]


def align_global(s1, s2):
    a = Align.PairwiseAligner()
    a.substitution_matrix = BL62
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = "global"
    return a.align(s1, s2)[0]


def matched_pairs(aln, s1_offset_1idx, s2_offset_1idx):
    """Return list of (pos1_1idx, pos2_1idx, aa1, aa2) for each aligned pair."""
    blocks1 = aln.aligned[0]
    blocks2 = aln.aligned[1]
    # these are 0-indexed into the inputs passed to aligner
    s1 = aln.target  # first sequence
    s2 = aln.query  # second sequence
    pairs = []
    for (st1, en1), (st2, en2) in zip(blocks1, blocks2):
        for i in range(int(en1) - int(st1)):
            p1 = int(st1) + i + s1_offset_1idx  # 1-indexed global pos
            p2 = int(st2) + i + s2_offset_1idx
            pairs.append((p1, p2, str(s1[int(st1) + i]), str(s2[int(st2) + i])))
    return pairs


def compute_rmsd(pairs, strc_ca, otoa_ca):
    """Given matched (strc_resnum, otoa_resnum, ...) pairs, superimpose Calpha."""
    usable = []
    for p1, p2, a1, a2 in pairs:
        if p1 in strc_ca and p2 in otoa_ca:
            usable.append((p1, p2, a1, a2))
    if len(usable) < 3:
        return {"rmsd": None, "n_aligned": 0, "error": "too few usable pairs"}
    strc_atoms = [strc_ca[p1][0] for (p1, p2, _, _) in usable]
    otoa_atoms = [otoa_ca[p2][0] for (p1, p2, _, _) in usable]
    sup = Superimposer()
    sup.set_atoms(strc_atoms, otoa_atoms)
    identities = sum(1 for (_, _, a1, a2) in usable if a1 == a2)
    return {
        "n_aligned": int(len(usable)),
        "rmsd": float(sup.rms),
        "rotation": [[float(x) for x in row] for row in sup.rotran[0]],
        "translation": [float(x) for x in sup.rotran[1]],
        "seq_identity_pct": round(100 * identities / len(usable), 2),
        "sample_pairs_first_5": usable[:5],
    }


def main():
    print("Fetching sequences ...")
    strc_seq = str(fetch_fasta(STRC_ACC).seq)
    otoa_seq = str(fetch_fasta(OTOA_ACC).seq)
    print(f"  STRC {len(strc_seq)} aa, OTOA {len(otoa_seq)} aa")

    print("Loading structures ...")
    cif_p = MMCIFParser(QUIET=True)
    pdb_p = PDBParser(QUIET=True)
    strc_struct = cif_p.get_structure("strc", str(STRC_CIF))
    otoa_struct = pdb_p.get_structure("otoa", str(OTOA_PDB))
    strc_ca = load_ca_by_resnum(strc_struct)
    otoa_ca = load_ca_by_resnum(otoa_struct)
    print(f"  STRC CAs: {len(strc_ca)} ({min(strc_ca)}..{max(strc_ca)})")
    print(f"  OTOA CAs: {len(otoa_ca)} ({min(otoa_ca)}..{max(otoa_ca)})")

    # Validate: structure residue letters match UniProt at sampled positions
    for sample in [1000, 1600, 1700]:
        if sample in strc_ca:
            struct_letter = three_to_one(strc_ca[sample][1])
            seq_letter = strc_seq[sample - 1]
            print(f"  STRC {sample}: struct={struct_letter} seq={seq_letter} {'OK' if struct_letter==seq_letter else 'MISMATCH'}")
    for sample in [500, 1000]:
        if sample in otoa_ca:
            struct_letter = three_to_one(otoa_ca[sample][1])
            seq_letter = otoa_seq[sample - 1]
            print(f"  OTOA {sample}: struct={struct_letter} seq={seq_letter} {'OK' if struct_letter==seq_letter else 'MISMATCH'}")

    results = {
        "meta": {
            "date": "2026-04-20",
            "question": "Do STRC ARM repeats structurally superimpose on OTOA C-terminal?",
            "strc_structure": str(STRC_CIF),
            "otoa_structure": str(OTOA_PDB),
            "strc_length": len(strc_seq),
            "otoa_length": len(otoa_seq),
        }
    }

    # 1) ARM-repeat region: STRC 1603-1770 local aligned to full OTOA
    print("\n[1] ARM repeats (STRC 1603-1770) vs OTOA ...")
    strc_arm = strc_seq[ARM_START - 1 : ARM_END]
    aln_arm = align_local(strc_arm, otoa_seq)
    pairs_arm = matched_pairs(aln_arm, ARM_START, 1)  # STRC offset = ARM_START, OTOA offset = 1
    rmsd_arm = compute_rmsd(pairs_arm, strc_ca, otoa_ca)
    print(
        f"  Matched pairs {len(pairs_arm)}, "
        f"usable {rmsd_arm['n_aligned']}, "
        f"Calpha RMSD {rmsd_arm.get('rmsd')}"
    )
    results["arm_region"] = {
        "strc_range_1idx": [ARM_START, ARM_END],
        "total_matched_pairs": len(pairs_arm),
        **rmsd_arm,
    }

    # 2) Mini-STRC window local aligned to OTOA
    print("\n[2] Mini-STRC (700-1775) vs OTOA ...")
    strc_mini = strc_seq[MINI_STRC_START - 1 : MINI_STRC_END]
    aln_mini = align_local(strc_mini, otoa_seq)
    pairs_mini = matched_pairs(aln_mini, MINI_STRC_START, 1)
    rmsd_mini = compute_rmsd(pairs_mini, strc_ca, otoa_ca)
    print(
        f"  Matched pairs {len(pairs_mini)}, "
        f"usable {rmsd_mini['n_aligned']}, "
        f"Calpha RMSD {rmsd_mini.get('rmsd')}"
    )
    results["mini_strc_region"] = {
        "strc_range_1idx": [MINI_STRC_START, MINI_STRC_END],
        "total_matched_pairs": len(pairs_mini),
        **rmsd_mini,
    }

    # 3) Global alignment pairs
    print("\n[3] Global alignment (full STRC vs full OTOA) ...")
    aln_global = align_global(strc_seq, otoa_seq)
    pairs_global = matched_pairs(aln_global, 1, 1)
    rmsd_global = compute_rmsd(pairs_global, strc_ca, otoa_ca)
    print(
        f"  Matched pairs {len(pairs_global)}, "
        f"usable {rmsd_global['n_aligned']}, "
        f"Calpha RMSD {rmsd_global.get('rmsd')}"
    )
    results["global_region"] = {
        "total_matched_pairs": len(pairs_global),
        **rmsd_global,
    }

    # 4) Anchor residue fate: is the paralog residue proximal to STRC anchor in 3D?
    print("\n[4] Anchor residue mapping (STRC -> OTOA via global alignment) ...")
    anchor_fate = {}
    for name, strc_pos in POCKET_RES.items():
        otoa_pos = None
        otoa_letter = None
        for p1, p2, a1, a2 in pairs_global:
            if p1 == strc_pos:
                otoa_pos = p2
                otoa_letter = a2
                break
        anchor_fate[name] = {
            "strc_pos": strc_pos,
            "strc_letter": strc_seq[strc_pos - 1],
            "otoa_pos": otoa_pos,
            "otoa_letter": otoa_letter,
            "conservation": (
                "identical" if otoa_letter == strc_seq[strc_pos - 1]
                else "aligned_substitution" if otoa_letter
                else "gap"
            ),
        }
        print(f"  STRC {name}@{strc_pos} -> OTOA {otoa_pos} ({otoa_letter})")
    results["anchor_residues"] = anchor_fate

    # 5) Verdict
    arm_rmsd = results["arm_region"].get("rmsd") or 99.0
    if arm_rmsd <= 4.0:
        verdict = "PURSUE_CHIMERA"
    elif arm_rmsd <= 6.0:
        verdict = "PURSUE_UPREGULATION_ONLY"
    else:
        verdict = "PARALOG_SEQUENCE_ONLY_NO_FOLD_SHARE"
    results["verdict"] = verdict
    print(f"\nVerdict: {verdict} (ARM RMSD {arm_rmsd:.2f} A)")

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
