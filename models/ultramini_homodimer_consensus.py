"""Ultra-Mini homodimer — 5-model consensus + pseudo-2-fold symmetry check.

The ranked-0 AF3 model of the Ultra-Mini homodimer showed 114 inter-chain
contacts clustered symmetrically in two zones: aa 1077-1131 (right after the
truncation cut point) and aa 1422-1590 (deep in the ARM repeats). Overall ipTM
0.28-0.30 is low but higher than prior homodimer AF3 jobs (0.20 for mini
594-1775, 0.24 for full STRC), with the improvement coming after we deleted
the disordered N-term.

Key question: is the aa 1422-1590 cluster a real weak-binding dimerization
surface, or an AF3 artifact from low-confidence packing?

Three falsification tests:

1. CROSS-MODEL CONSENSUS. If all 5 AF3 models converge on the same 1422-1590
   zone, that is real geometric signal. If models disagree on which residues
   contact, it is sampling noise.

2. PSEUDO-2-FOLD SYMMETRY. A real homodimer obeys C2 symmetry: residue X on
   chain A contacts residue X on chain B (homotypic self-contacts) or
   swapped-symmetric pairs (A:r1 <-> B:r2 AND A:r2 <-> B:r1). If contacts are
   asymmetric (A:r1 <-> B:r2 without the reverse), it is AF3 random packing.

3. SOLVENT EXPOSURE IN MONOMER. If the contact residues are buried in the
   monomer's own fold (SASA < 20%), they cannot be a real binding surface.
   Only solvent-exposed residues can mediate inter-chain interaction.
"""
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

warnings.filterwarnings("ignore")

RESULTS_DIR = Path("/Users/egorlyfar/DeepResearch/strc/af3-results/job-ultramini-homodimer")
MONO_CIF = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job-h-strc-cterm-only.cif")
OFFSET = 1075
CUTOFF = 5.0
OUT = Path("/Users/egorlyfar/Brain/research/strc/models/ultramini_homodimer_consensus.json")

parser = MMCIFParser(QUIET=True)


def extract_contacts(cif_path):
    s = parser.get_structure("x", cif_path)
    m = next(s.get_models())
    chains = list(m.get_chains())
    cA, cB = chains[0], chains[1]
    atoms_A = [a for a in cA.get_atoms() if a.element != "H"]
    atoms_B = [a for a in cB.get_atoms() if a.element != "H"]
    ns = NeighborSearch(atoms_B)
    pairs = {}
    for a in atoms_A:
        for b in ns.search(a.coord, CUTOFF):
            rA = a.get_parent().id[1] + OFFSET - 1
            rB = b.get_parent().id[1] + OFFSET - 1
            d = a - b
            k = (rA, rB)
            if k not in pairs or d < pairs[k]:
                pairs[k] = d
    return pairs, cA, cB


# ==== TEST 1: 5-model consensus ====
print("=" * 80)
print("TEST 1 — 5-model consensus on interface residues")
print("=" * 80)

model_pairs = {}
all_contact_residues = Counter()  # how many models show each (chain_A) residue
all_contact_pairs = Counter()     # how many models show each (A,B) pair

for i in range(5):
    pairs, _, _ = extract_contacts(RESULTS_DIR / f"fold_strc_ultramini_homodimer_model_{i}.cif")
    model_pairs[i] = pairs
    res_A_this_model = {k[0] for k in pairs}
    for r in res_A_this_model:
        all_contact_residues[r] += 1
    for k in pairs:
        all_contact_pairs[k] += 1

# Residues present in >=3 of 5 models
consensus_residues = sorted([r for r, n in all_contact_residues.items() if n >= 3])
strict_residues = sorted([r for r, n in all_contact_residues.items() if n == 5])

# Zone classification
zone_early = [r for r in consensus_residues if 1075 <= r <= 1200]  # stump zone
zone_mid = [r for r in consensus_residues if 1200 < r <= 1400]
zone_late = [r for r in consensus_residues if 1400 < r <= 1775]    # deep ARM

print(f"Residues contacting in ALL 5 models (strict): {len(strict_residues)}  {strict_residues[:20]}")
print(f"Residues contacting in >=3/5 models (consensus): {len(consensus_residues)}")
print(f"  Stump zone (1075-1200): {len(zone_early)}  {zone_early}")
print(f"  Mid (1201-1400):        {len(zone_mid)}  {zone_mid}")
print(f"  Deep ARM (1401-1775):   {len(zone_late)}  {zone_late}")

total_pairs_per_model = [len(model_pairs[i]) for i in range(5)]
print(f"\nPer-model pair counts: {total_pairs_per_model} (mean {np.mean(total_pairs_per_model):.0f})")

pairs_in_all_5 = sorted([k for k, n in all_contact_pairs.items() if n == 5], key=lambda x: x[0])
print(f"Pairs identical across ALL 5 models: {len(pairs_in_all_5)}")
for p in pairs_in_all_5[:10]:
    print(f"  A {p[0]:4d} <-> B {p[1]:4d}")

# ==== TEST 2: pseudo-2-fold symmetry check ====
print("\n" + "=" * 80)
print("TEST 2 — pseudo-2-fold C2 symmetry check")
print("=" * 80)

def symmetry_score(pairs):
    """For each (A:r1, B:r2), check if (A:r2, B:r1) also exists. Return fraction
    of symmetric pairs. True C2 dimer = fraction ~1.0; AF3 random = ~0.0.
    Homotypic self-contacts (r1==r2) count as automatic symmetry."""
    matched = 0
    total = len(pairs)
    for (rA, rB) in pairs.keys():
        if rA == rB:
            matched += 1
        elif (rB, rA) in pairs:
            matched += 1
    return matched, total


def homotypic_self_contacts(pairs):
    """Count A:r <-> B:r same-residue contacts."""
    return sorted([r for (rA, rB) in pairs.keys() if rA == rB for r in [rA]])


for i in range(5):
    sym, total = symmetry_score(model_pairs[i])
    self_contacts = homotypic_self_contacts(model_pairs[i])
    print(f"Model {i}: {sym}/{total} pairs are C2-symmetric ({100*sym/total:.0f}%), "
          f"{len(self_contacts)} homotypic self-contacts: {self_contacts[:8]}{'...' if len(self_contacts)>8 else ''}")

# Self-contacts consensus
self_contact_residues = Counter()
for i in range(5):
    for r in homotypic_self_contacts(model_pairs[i]):
        self_contact_residues[r] += 1
consensus_self = sorted([r for r, n in self_contact_residues.items() if n >= 3])
print(f"\nHomotypic self-contacts (same residue on both chains) consensus (>=3/5 models):")
print(f"  {consensus_self}")

# ==== TEST 3: solvent accessibility in monomer fold ====
print("\n" + "=" * 80)
print("TEST 3 — solvent accessibility of contact residues in monomer fold")
print("=" * 80)

# Use job-h-strc-cterm-only.cif as the Ultra-Mini monomer reference fold.
mono_struct = parser.get_structure("mono", MONO_CIF)
sr = ShrakeRupley()
sr.compute(mono_struct, level="R")
mono_chain = list(next(mono_struct.get_models()).get_chains())[0]
res_sasa = {}
for res in mono_chain:
    if res.id[0] != " ":
        continue
    local = res.id[1]
    canonical = local + OFFSET - 1
    res_sasa[canonical] = res.sasa

# Normalize against max-CA Shrake-Rupley estimate (crude proxy ~ 200 A^2 per residue)
# Use relative SASA as absolute_SASA / max_sasa_for_that_aa (approximate).
# Simpler: use absolute SASA threshold — residues with SASA < 20 A^2 are buried.
BURIED_THRESHOLD = 20.0

print(f"Consensus contact residues and their monomer SASA (A^2):")
print(f"  {'resnum':>6s}  {'SASA':>6s}  {'status':10s}  {'models':>6s}")
for r in consensus_residues[:30]:
    sasa = res_sasa.get(r, -1)
    status = "exposed" if sasa >= BURIED_THRESHOLD else "buried" if sasa >= 0 else "missing"
    print(f"  {r:6d}  {sasa:6.1f}  {status:10s}  {all_contact_residues[r]:>4d}/5")
if len(consensus_residues) > 30:
    print(f"  ... ({len(consensus_residues)-30} more)")

exposed_consensus = [r for r in consensus_residues if res_sasa.get(r, -1) >= BURIED_THRESHOLD]
buried_consensus = [r for r in consensus_residues if 0 <= res_sasa.get(r, -1) < BURIED_THRESHOLD]
print(f"\nConsensus contact residues:")
print(f"  SOLVENT-EXPOSED (SASA >= {BURIED_THRESHOLD} A^2): {len(exposed_consensus)}/{len(consensus_residues)}")
print(f"  BURIED (SASA <  {BURIED_THRESHOLD} A^2): {len(buried_consensus)}/{len(consensus_residues)}")

frac_exposed = len(exposed_consensus) / len(consensus_residues) if consensus_residues else 0.0

# ==== VERDICT ====
print("\n" + "=" * 80)
print("FALSIFICATION VERDICT")
print("=" * 80)

consensus_strength = len(strict_residues) / max(1, len(all_contact_residues))
sym_fractions = [symmetry_score(model_pairs[i])[0] / max(1, symmetry_score(model_pairs[i])[1]) for i in range(5)]
mean_sym = np.mean(sym_fractions)

print(f"1. Consensus strength:")
print(f"   {len(strict_residues)}/{len(all_contact_residues)} residues contact in ALL 5 models = {100*consensus_strength:.0f}%")
print(f"   Consensus zone: {'deep ARM (1422-1590)' if len(zone_late) > len(zone_early) else 'stump (1075-1200)' if zone_early else 'mixed'}")
print(f"2. Symmetry score (mean across models): {100*mean_sym:.0f}% pairs are C2-symmetric")
print(f"3. Solvent exposure: {100*frac_exposed:.0f}% of consensus residues are surface-exposed in monomer")

if consensus_strength >= 0.3 and mean_sym >= 0.5 and frac_exposed >= 0.5:
    verdict = "real_weak_dimerization_surface"
elif mean_sym < 0.3 or frac_exposed < 0.3:
    verdict = "af3_artifact"
else:
    verdict = "ambiguous"
print(f"\nFINAL VERDICT: {verdict}")

out = {
    "test_1_consensus": {
        "residues_in_all_5_models_strict": strict_residues,
        "residues_in_3_of_5_or_more": consensus_residues,
        "zone_breakdown": {
            "stump_1075_1200": zone_early,
            "mid_1201_1400": zone_mid,
            "deep_arm_1401_1775": zone_late,
        },
        "pairs_identical_in_all_5_models": [{"chain_A_aa": p[0], "chain_B_aa": p[1]} for p in pairs_in_all_5],
        "per_model_pair_counts": total_pairs_per_model,
        "consensus_strength_pct": round(100 * consensus_strength, 1),
    },
    "test_2_symmetry": {
        "per_model_c2_fraction": [round(100*s, 1) for s in sym_fractions],
        "mean_c2_fraction_pct": round(100 * mean_sym, 1),
        "homotypic_self_contacts_consensus": consensus_self,
    },
    "test_3_solvent_accessibility": {
        "buried_threshold_A2": BURIED_THRESHOLD,
        "exposed_consensus_residues": exposed_consensus,
        "buried_consensus_residues": buried_consensus,
        "fraction_exposed_pct": round(100 * frac_exposed, 1),
        "per_residue_sasa": {str(r): round(res_sasa.get(r, -1), 2) for r in consensus_residues},
    },
    "verdict": verdict,
    "decision_criteria": {
        "real_weak_dimerization_surface": "consensus >=30% AND symmetry >=50% AND exposure >=50%",
        "af3_artifact": "symmetry <30% OR exposure <30%",
        "ambiguous": "everything else",
    },
}
OUT.write_text(json.dumps(out, indent=2, default=float))
print(f"\nSaved: {OUT}")
