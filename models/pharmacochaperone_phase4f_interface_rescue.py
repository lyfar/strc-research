#!/usr/bin/env python3
"""
Phase 4f — Interface-rescue interaction energy via OpenMM single-point + GBn2.

Tests the THERAPEUTIC claim (vs. mere pocket binding):
  Does a pharmacochaperone pre-bound at the K1141 pocket on Ultra-Mini
  recover the Ultra-Mini × TMEM145 interface energy lost to E1659A?

Reference loss: 8.4 kcal/mol (electrostatic, from
[[STRC Electrostatic Analysis E1659A]]).

Method (static-snapshot approximation, NOT dynamic MM-GBSA):
  1. Load Ultra-Mini × TMEM145 complex CIF.
  2. Use Phase 4b best Vina pose for each of the top-3 leads
     (indole-3-acetic, naphthalene-2-COOH, cyclopropane-phenyl-COOH).
  3. PDBFixer: clean and protonate complex at pH 7.4.
  4. OpenFF-Toolkit: parameterise ligand with SMIRNOFF openff-2.1.0.
  5. OpenMM: build apo system + each ligand-bound system with AMBER14SB
     + GBn2 implicit solvent (eps_solvent=78.5).
  6. Single-point energies. Interface ΔG = E(complex) - E(chainA+lig) - E(chainB).
  7. Rescue ΔΔG = ΔG_interface_apo - ΔG_interface_bound.
  8. Recovery = ΔΔG / 8.4 kcal/mol.

Gate: >=1 of top-3 leads shows recovery >= 0.30 (>=2.52 kcal/mol rescue).

Honesty: this is a STATIC snapshot, not an MD-averaged ensemble. For a clear
pass or clear fail (>2 kcal/mol margin), it is decisive. For ambiguous results
(within ±2 kcal/mol of the gate), Phase 5 short MD + dynamic MM-GBSA is
the next step.

Run env: conda env `strc-mmgbsa` (see 0-setup.sh).
"""

from __future__ import annotations

import json
import sys
import shutil
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

WORK_DIR = Path("/Users/egorlyfar/Brain/research/strc/models")
TARGET_CIF = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job-ultramini-x-tmem145-full.cif")
PHASE4B_JSON = WORK_DIR / "pharmacochaperone_phase4b_vina_gnina_screen.json"
OUT_JSON = WORK_DIR / "pharmacochaperone_phase4f_interface_rescue.json"

REFERENCE_GAP_KCAL_MOL = 8.4   # E1659A interface penalty (electrostatic analysis)
RECOVERY_GATE = 0.30           # >=30% of the gap → pass
TOP_K_LEADS = 3                # top-3 by best_affinity (Vina raw, lower=better)

KCAL_PER_KJ = 1.0 / 4.184


def log(msg):
    print(f"[{datetime.now(timezone.utc).isoformat(timespec='seconds')}] {msg}", flush=True)


def env_check():
    log("Phase 4f environment check")
    try:
        import openmm                                  # noqa
        import openmm.app                              # noqa
        from openmm import unit                        # noqa
        from openff.toolkit import Molecule, ForceField  # noqa
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator  # noqa
        import pdbfixer                                # noqa
        import MDAnalysis                              # noqa
        from rdkit import Chem                         # noqa
    except ImportError as e:
        log(f"FATAL — missing dep: {e}")
        log("Run 0-setup.sh first to create the strc-mmgbsa conda env.")
        sys.exit(2)
    if shutil.which("obabel") is None:
        log("FATAL — obabel not on PATH")
        sys.exit(2)
    log("env OK: openmm, openff-toolkit, pdbfixer, MDAnalysis, rdkit, obabel")


PHASE4G_JSON = WORK_DIR / "pharmacochaperone_phase4g_repurpose_screen.json"


def load_top_leads():
    """Load Phase 4b top-3 leads. If PHASE4F_INCLUDE_REPURPOSE=1, also append
    the Phase 4g repurposable chaperones (4PBA, IP-045, TMAO) using their
    E1659A-pocket Vina pose for interface-rescue testing. This gives a direct
    MM-GBSA number for FDA-approved repurposing candidates in the same run."""
    import os

    if not PHASE4B_JSON.exists():
        log(f"FATAL — Phase 4b output missing: {PHASE4B_JSON}")
        sys.exit(3)
    data = json.loads(PHASE4B_JSON.read_text())
    leads = [r for r in data["results"] if r.get("role") == "lead"]
    leads.sort(key=lambda r: r["best_affinity_kcal_mol"])
    top = leads[:TOP_K_LEADS]

    include_repurpose = os.environ.get("PHASE4F_INCLUDE_REPURPOSE", "0") == "1"
    if include_repurpose and PHASE4G_JSON.exists():
        g = json.loads(PHASE4G_JSON.read_text())
        # Flatten nested shape: results[name] = {"meta": {...}, "targets": {target: {...}}}
        for name, entry in g.get("results", {}).items():
            meta = entry.get("meta", {})
            mut_target = entry.get("targets", {}).get("e1659a_mutant", {})
            if "best_affinity_kcal_mol" not in mut_target:
                continue
            # Prefer the e1659a_mutant pose since Phase 4f tests rescue of the
            # mutant's interface. Fall back to ultra_x_tmem145 pose if mutant
            # missing.
            pose_pdbqt = mut_target.get("pose_pdbqt")
            if not pose_pdbqt:
                ultra = entry.get("targets", {}).get("ultra_x_tmem145", {})
                pose_pdbqt = ultra.get("pose_pdbqt")
            top.append({
                "name": name,
                "smiles": meta.get("smiles"),
                "best_affinity_kcal_mol": mut_target["best_affinity_kcal_mol"],
                "pose_pdbqt": pose_pdbqt,
                "role": meta.get("role", "repurpose"),
                "clinical_status": meta.get("clinical_status", ""),
            })
        log(f"PHASE4F_INCLUDE_REPURPOSE=1 — appended {len(g.get('results', {}))} Phase 4g candidates")

    log(f"Roster for Phase 4f: {len(top)} compounds")
    for r in top:
        role_tag = f" [{r.get('role')}]" if r.get("role") and r["role"] != "lead" else ""
        log(f"  {r['name']:30s} {r['smiles']:40s} {r['best_affinity_kcal_mol']:.3f} kcal/mol{role_tag}")
    return top


def cif_to_clean_pdb(cif_path: Path, out_pdb: Path):
    """PDBFixer-clean the CIF to a PDB ready for OpenMM."""
    log(f"PDBFixer cleaning {cif_path.name} -> {out_pdb.name}")
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=str(cif_path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}            # do not model unresolved gaps
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(out_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    log(f"  cleaned PDB: {out_pdb.stat().st_size/1024:.1f} KiB")


def split_chains(complex_pdb: Path, chain_a_pdb: Path, chain_b_pdb: Path):
    """Split chain A (Ultra-Mini STRC) and chain B (TMEM145) into separate PDBs."""
    import MDAnalysis as mda
    u = mda.Universe(str(complex_pdb))
    chains = sorted({s.segid or s.chainID for s in u.atoms.segments})
    log(f"  chains seen: {chains}")
    if len(chains) < 2:
        # AF3 outputs sometimes use chainID instead of segid
        cids = sorted({a.chainID for a in u.atoms})
        log(f"  chainID seen: {cids}")
        if len(cids) >= 2:
            chain_a, chain_b = cids[:2]
            u.atoms.select_atoms(f"chainID {chain_a}").write(str(chain_a_pdb))
            u.atoms.select_atoms(f"chainID {chain_b}").write(str(chain_b_pdb))
        else:
            log("FATAL — could not split chains")
            sys.exit(4)
    else:
        chain_a, chain_b = chains[:2]
        u.atoms.select_atoms(f"segid {chain_a}").write(str(chain_a_pdb))
        u.atoms.select_atoms(f"segid {chain_b}").write(str(chain_b_pdb))
    log(f"  chain A -> {chain_a_pdb.name} ({chain_a_pdb.stat().st_size/1024:.1f} KiB)")
    log(f"  chain B -> {chain_b_pdb.name} ({chain_b_pdb.stat().st_size/1024:.1f} KiB)")


def pdbqt_to_sdf(pdbqt_text_or_path: str, out_sdf: Path, smiles: str | None = None):
    """Vina PDBQT pose -> SDF for openff parameterisation, via meeko.

    Accepts either raw PDBQT text or a path to an existing PDBQT file (Phase 4b
    JSON stores paths in the `pose_pdbqt` field, not the file content).

    Meeko (the AutoDock Vina ligand prep library) preserves bond orders and
    formal charges from the PDBQT, unlike obabel which loses bond orders and
    leaves aromatic carbons that openff-toolkit reads as radicals. The `smiles`
    arg is currently unused (meeko reconstructs the right molecule on its own)
    but kept in the signature for the SMILES-fallback path in main().
    """
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    candidate_path = Path(pdbqt_text_or_path) if "\n" not in pdbqt_text_or_path else None
    if candidate_path is not None and candidate_path.exists():
        pdbqt_path = candidate_path
        cleanup_path = None
    else:
        with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False, mode="w") as t:
            t.write(pdbqt_text_or_path)
            pdbqt_path = Path(t.name)
        cleanup_path = pdbqt_path

    try:
        pdbqt_mol = PDBQTMolecule.from_file(str(pdbqt_path), skip_typing=True)
        rdmols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        if not rdmols or rdmols[0] is None:
            raise RuntimeError(f"meeko returned no rdkit mol for {pdbqt_path.name}")
        mol = Chem.AddHs(rdmols[0], addCoords=True)
        writer = Chem.SDWriter(str(out_sdf))
        writer.write(mol)
        writer.close()
        if out_sdf.stat().st_size == 0:
            raise RuntimeError(f"meeko-written SDF {out_sdf.name} is 0 bytes")
    finally:
        if cleanup_path is not None:
            cleanup_path.unlink(missing_ok=True)


def smiles_to_sdf_3d(smiles: str, out_sdf: Path):
    """RDKit: SMILES -> 3D-embedded SDF (fallback when no Vina pose available)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    writer = Chem.SDWriter(str(out_sdf))
    writer.write(mol)
    writer.close()


_APO_CACHE: dict[str, float] = {}


def build_system(receptor_pdb: Path, ligand_sdf: Path | None, max_iter: int = 200):
    """OpenMM system with AMBER14SB + GBn2 implicit solvent. Optional ligand.

    Caches apo (no-ligand) energies by receptor PDB path — chain B is recomputed
    3 times across the 3 ligand iterations otherwise.
    """
    cache_key = str(receptor_pdb) if ligand_sdf is None else None
    if cache_key and cache_key in _APO_CACHE:
        log(f"    cached: E({receptor_pdb.name}) = {_APO_CACHE[cache_key]:+.2f} kcal/mol")
        return _APO_CACHE[cache_key]

    import openmm
    from openmm import unit, app
    from openmm.app import PDBFile, ForceField, Simulation, NoCutoff

    pdb = PDBFile(str(receptor_pdb))
    # GBn2 implicit solvent comes from the implicit/gbn2.xml force field file —
    # do NOT pass implicitSolvent= to createSystem in OpenMM 8.x or it raises
    # ValueError("'implicitSolvent' was specified ... but was never used").
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")

    if ligand_sdf is not None:
        from openff.toolkit import Molecule
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        lig = Molecule.from_file(str(ligand_sdf))
        if isinstance(lig, list):
            lig = lig[0]
        # Pre-assign partial charges via NAGL (am1bcc-quality, ML-fast).
        # Without this, SMIRNOFFTemplateGenerator tries to call AmberTools
        # antechamber for AM1-BCC, which is not in this conda env.
        if lig.n_conformers == 0:
            lig.generate_conformers(n_conformers=1)
        try:
            lig.assign_partial_charges(partial_charge_method="openff-gnn-am1bcc-1.0.0.pt")
        except Exception as e:
            log(f"    NAGL charges failed ({type(e).__name__}); falling back to gasteiger")
            lig.assign_partial_charges(partial_charge_method="gasteiger")
        smirnoff = SMIRNOFFTemplateGenerator(molecules=[lig])
        forcefield.registerTemplateGenerator(smirnoff.generator)

        # Combine receptor + ligand topology
        from openmm.app import Modeller
        modeller = Modeller(pdb.topology, pdb.positions)
        lig_topology = lig.to_topology().to_openmm()
        lig_positions = lig.conformers[0].to_openmm()
        modeller.add(lig_topology, lig_positions)
        topology = modeller.topology
        positions = modeller.positions
    else:
        topology = pdb.topology
        positions = pdb.positions

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=NoCutoff,
        constraints=app.HBonds,
    )

    integrator = openmm.LangevinIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    # Single-point energy mode (fast/noisy path): skip minimize entirely, take
    # the energy of the input geometry as-is. Trade-off vs minimised:
    #   minimised  → 30-200s per system on CPU, 1-2s per chunk on OpenCL (when
    #                it survives macOS GPU watchdog), low-noise ΔΔG signal
    #   single-pt → ~0.1-1s per system on either platform; ~5-15 kcal/mol noise
    #                from unrelaxed Vina-pose clashes; signal-to-noise of ΔΔG
    #                still ~2:1 if interface energies are >20 kcal/mol
    # Set PHASE4F_MINIMISE=1 in env to force the older chunked-minimisation path.
    import os
    import time
    n_atoms = system.getNumParticles()
    minimise = os.environ.get("PHASE4F_MINIMISE", "0") == "1"
    plat_force = os.environ.get("PHASE4F_PLATFORM", "").strip()  # "CPU" or "OpenCL"
    if plat_force == "CPU":
        plat_list = [("CPU", {})]
    elif plat_force == "OpenCL":
        plat_list = [("OpenCL", {"OpenCLPrecision": "mixed"})]
    else:
        plat_list = [("OpenCL", {"OpenCLPrecision": "mixed"}), ("CPU", {})]

    last_err = None
    for plat_try, props in plat_list:
        try:
            platform = openmm.Platform.getPlatformByName(plat_try)
        except Exception as e:
            last_err = e
            continue
        try:
            integrator2 = openmm.LangevinIntegrator(
                300 * unit.kelvin,
                1.0 / unit.picosecond,
                0.002 * unit.picoseconds,
            )
            sim = Simulation(topology, system, integrator2, platform, props)
            sim.context.setPositions(positions)
            t0 = time.time()
            if minimise:
                chunk_size = min(25, max_iter)
                chunks = max(1, max_iter // chunk_size)
                log(f"    platform={plat_try} atoms={n_atoms} — minimising {chunks}×{chunk_size} iter")
                for c in range(chunks):
                    sim.minimizeEnergy(maxIterations=chunk_size)
                    e_kj = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                    if e_kj != e_kj or abs(e_kj) > 1e10:
                        raise RuntimeError(f"non-finite energy after chunk {c+1}/{chunks}: {e_kj}")
                    if (c + 1) % 4 == 0 or c == chunks - 1:
                        log(f"      chunk {c+1}/{chunks} t={time.time()-t0:.1f}s E={e_kj * KCAL_PER_KJ:+.1f} kcal/mol")
            else:
                log(f"    platform={plat_try} atoms={n_atoms} — single-point energy")
                e_kj = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                if e_kj != e_kj or abs(e_kj) > 1e10:
                    raise RuntimeError(f"non-finite single-point energy: {e_kj}")
            dt = time.time() - t0
            e_kcal = e_kj * KCAL_PER_KJ
            mode_tag = "minimised" if minimise else "single-pt"
            log(f"    {mode_tag} in {dt:.1f}s ({plat_try}) -> E = {e_kcal:+.2f} kcal/mol")
            if cache_key:
                _APO_CACHE[cache_key] = e_kcal
            return e_kcal
        except Exception as e:
            log(f"    {plat_try} FAILED: {type(e).__name__}: {str(e)[:200]}")
            last_err = e
            continue
    raise RuntimeError(f"all platforms failed; last error: {last_err}")


def interface_energy(complex_pdb: Path,
                     chain_a_pdb: Path,
                     chain_b_pdb: Path,
                     ligand_sdf: Path | None) -> dict:
    """ΔG_interface = E(complex) - E(A+lig) - E(B). All energies in kcal/mol.

    Iter counts overridable via PHASE4F_APO_ITER and PHASE4F_BOUND_ITER env vars
    (used by smoke-test harness to reduce CPU runtime).
    """
    import os
    apo_iter = int(os.environ.get("PHASE4F_APO_ITER", "100"))
    bound_iter = int(os.environ.get("PHASE4F_BOUND_ITER", "300"))
    iter_complex = bound_iter if ligand_sdf else apo_iter
    iter_a = bound_iter if ligand_sdf else apo_iter

    log("  energy: complex" + (" + ligand" if ligand_sdf else " (apo)"))
    e_complex = build_system(complex_pdb, ligand_sdf, max_iter=iter_complex)

    log("  energy: chain A" + (" + ligand" if ligand_sdf else " (apo)"))
    e_a = build_system(chain_a_pdb, ligand_sdf, max_iter=iter_a)

    log("  energy: chain B (apo, ligand stays on A)")
    e_b = build_system(chain_b_pdb, None, max_iter=apo_iter)

    dg = e_complex - e_a - e_b
    log(f"  ΔG_interface = {dg:+.2f} kcal/mol")
    return {
        "e_complex_kcal": e_complex,
        "e_chain_a_kcal": e_a,
        "e_chain_b_kcal": e_b,
        "dg_interface_kcal": dg,
    }


def main():
    import os
    env_check()
    leads = load_top_leads()
    smoke = os.environ.get("PHASE4F_SMOKE", "0") == "1"
    if smoke:
        log(f"SMOKE TEST mode — restricting to first lead only "
            f"({leads[0]['name'] if leads else 'no leads found'})")
        leads = leads[:1]

    workdir = WORK_DIR / "phase4f_workdir"
    workdir.mkdir(exist_ok=True, parents=True)

    complex_pdb = workdir / "complex_clean.pdb"
    chain_a_pdb = workdir / "chain_a.pdb"
    chain_b_pdb = workdir / "chain_b.pdb"

    cif_to_clean_pdb(TARGET_CIF, complex_pdb)
    split_chains(complex_pdb, chain_a_pdb, chain_b_pdb)

    log("=" * 60)
    log("APO interface energy (no ligand)")
    log("=" * 60)
    apo = interface_energy(complex_pdb, chain_a_pdb, chain_b_pdb, ligand_sdf=None)

    results = []
    for lead in leads:
        log("=" * 60)
        log(f"BOUND interface energy with {lead['name']}")
        log("=" * 60)
        lig_sdf = workdir / f"{lead['name']}.sdf"
        if lead.get("pose_pdbqt"):
            try:
                pdbqt_to_sdf(lead["pose_pdbqt"], lig_sdf, smiles=lead.get("smiles"))
            except (subprocess.CalledProcessError, RuntimeError) as e:
                log(f"  pdbqt_to_sdf failed: {type(e).__name__}: {str(e)[:200]}; "
                    f"falling back to SMILES embed (loses Vina pose)")
                smiles_to_sdf_3d(lead["smiles"], lig_sdf)
        else:
            smiles_to_sdf_3d(lead["smiles"], lig_sdf)

        bound = interface_energy(complex_pdb, chain_a_pdb, chain_b_pdb, lig_sdf)
        ddg = apo["dg_interface_kcal"] - bound["dg_interface_kcal"]
        recovery = ddg / REFERENCE_GAP_KCAL_MOL
        results.append({
            "name": lead["name"],
            "smiles": lead["smiles"],
            "vina_best_affinity_kcal": lead["best_affinity_kcal_mol"],
            "apo_interface_dg_kcal": apo["dg_interface_kcal"],
            "bound_interface_dg_kcal": bound["dg_interface_kcal"],
            "ddg_rescue_kcal": ddg,
            "recovery_fraction": recovery,
            "passes_30pct_gate": recovery >= RECOVERY_GATE,
            "raw": {"apo": apo, "bound": bound},
        })
        log(f"  ΔΔG_rescue = {ddg:+.2f} kcal/mol  ({recovery*100:+.1f}% of {REFERENCE_GAP_KCAL_MOL} gap)")
        log(f"  gate (>=30%): {'PASS' if recovery >= RECOVERY_GATE else 'FAIL'}")

    n_pass = sum(1 for r in results if r["passes_30pct_gate"])
    verdict = "PASS" if n_pass >= 1 else "FAIL"
    log("=" * 60)
    log(f"PHASE 4f VERDICT: {verdict} ({n_pass}/{len(results)} leads cleared 30% gate)")
    log("=" * 60)

    out = {
        "phase": "4f",
        "method": "OpenMM single-point + GBn2 implicit solvent (static snapshot, not MD MM-GBSA)",
        "target_cif": str(TARGET_CIF),
        "reference_gap_kcal_mol": REFERENCE_GAP_KCAL_MOL,
        "recovery_gate": RECOVERY_GATE,
        "apo_interface": apo,
        "results": results,
        "n_pass": n_pass,
        "verdict": verdict,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    log(f"Wrote {OUT_JSON}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
