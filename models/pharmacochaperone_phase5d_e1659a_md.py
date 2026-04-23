#!/usr/bin/env python3
"""
Phase 5d — E1659A mutant MD.

WHY THIS
════════
Phase 5a ran MD on WT Ultra-Mini × TMEM145 complex chain A (residues 1-701,
sequential numbering for a 594-1294 full-length construct). That construct
does NOT include residue E1659 — the disease mutation sits 365 residues
outside the Ultra-Mini boundary. Phase 5a/5b/5c therefore validated the
K1141 pocket on a WT fold that never carried the mutation.

If E1659A distorts the K1141 pocket (shifts its residues, alters void
volume, or opens an alternative cavity), our entire Phase 3c/4b/5b docking
strategy has been aiming at the wrong target geometry.

This phase closes that gap: run MD on the AF3-predicted E1659A full-length
structure (job3-mutant.cif, residues 1-1775, chain A, with both K1141 and
E1659A) and compare K1141 pocket behavior to Phase 5a WT.

Phase 5e (follow-up) will re-dock Phase 3c v2 top hits (niflumic, diflunisal,
flufenamic) against the mutant ensemble and compare f_PC(mut) vs f_PC(WT).

═══════════════════════════════════════════════════════════════════════════
METHOD
═══════════════════════════════════════════════════════════════════════════

1. Source: job3-mutant.cif (AF3 prediction of full-length E1659A STRC),
   chain A only. Verify K1141=LYS, E1659=ALA before proceeding.
2. Prep with PDBFixer: add missing atoms/H, protonate at pH 7.4.
3. Force field: AMBER14SB + TIP3P water, NaCl 0.15 M neutralisation.
   (Identical FF/conditions to Phase 5a so trajectories are comparable.)
4. Box: cubic, 1 nm pad around solute. Expect ~400k atoms (2.5× Phase 5a).
5. Minimisation: 5000 steps.
6. NVT: 50 ps @ 310 K, Langevin thermostat, 2 fs timestep.
7. NPT: 50 ps @ 1 bar Monte Carlo barostat.
8. Production NPT: 2 ns (same as Phase 5a for direct comparability).
9. Save 20 snapshots, full system (water kept so we can strip consistently
   in 5e).

═══════════════════════════════════════════════════════════════════════════
RUN
═══════════════════════════════════════════════════════════════════════════

    python3.11 pharmacochaperone_phase5d_e1659a_md.py                 # SMOKE
    python3.11 pharmacochaperone_phase5d_e1659a_md.py --ns 2 --snapshots 20

═══════════════════════════════════════════════════════════════════════════
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

WORK = Path("/Users/egorlyfar/Brain/research/strc/models")
TARGET_CIF = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job3-mutant.cif")
OUT_JSON = WORK / "pharmacochaperone_phase5d_e1659a_md.json"
SNAPSHOTS_DIR = WORK / "artifacts" / "phase5d_snapshots"

# Sanity checks on AF3 mutant structure
EXPECTED_K1141_RESNAME = "LYS"
EXPECTED_E1659_RESNAME = "ALA"  # E1659A mutation → ALA in mutant CIF


def log(msg: str):
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    print(f"[{ts}] {msg}", flush=True)


def env_check():
    try:
        import openmm
        from openmm import app, unit, Platform, LangevinIntegrator, MonteCarloBarostat
        import pdbfixer
        from Bio.PDB import MMCIFParser, PDBIO, Select
        from Bio.PDB.Polypeptide import is_aa
    except ImportError as e:
        log(f"FATAL import: {e}")
        sys.exit(2)
    log(f"openmm {openmm.version.version} OK")


def verify_mutant_cif(cif_path: Path) -> dict:
    """Confirm the input CIF actually carries E1659A.

    Prevents silently running MD on the wrong structure (e.g. WT CIF
    accidentally swapped in)."""
    from Bio.PDB import MMCIFParser

    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("mut", str(cif_path))
    model = next(iter(struct))
    chain = model["A"]
    residues = {r.id[1]: r.get_resname() for r in chain if r.id[0] == " "}

    k1141 = residues.get(1141, "MISSING")
    e1659 = residues.get(1659, "MISSING")

    if k1141 != EXPECTED_K1141_RESNAME:
        raise RuntimeError(
            f"K1141 sanity check failed: expected {EXPECTED_K1141_RESNAME}, got {k1141}. "
            f"Input CIF {cif_path} is not the STRC E1659A mutant construct."
        )
    if e1659 != EXPECTED_E1659_RESNAME:
        raise RuntimeError(
            f"E1659 mutation check failed: expected {EXPECTED_E1659_RESNAME} (mutant), "
            f"got {e1659}. Input CIF {cif_path} appears to be WT STRC, not E1659A."
        )

    first = min(residues)
    last = max(residues)
    log(f"Mutant CIF verified: chain A res {first}-{last} ({len(residues)} residues)")
    log(f"  K1141 = {k1141} ✓   E1659 = {e1659} ✓ (ALA = mutant)")
    return {"first": first, "last": last, "n_residues": len(residues)}


def extract_chain_a(cif_path: Path, out_pdb: Path):
    """CIF → chain A only PDB, standard AA only (no het, no waters)."""
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from Bio.PDB.Polypeptide import is_aa

    class ChainAStandard(Select):
        def accept_chain(self, chain):
            return chain.id == "A"

        def accept_residue(self, res):
            return is_aa(res, standard=True)

    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("mut", str(cif_path))
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out_pdb), select=ChainAStandard())


def run_md(production_ps: float, n_snapshots: int,
           platform_name: str | None = None) -> dict:
    import openmm
    from openmm import app, unit, Platform
    from openmm import LangevinIntegrator, MonteCarloBarostat
    import pdbfixer

    SNAPSHOTS_DIR.mkdir(parents=True, exist_ok=True)
    chain_pdb = SNAPSHOTS_DIR / "mutant_chainA.pdb"
    log(f"Extracting chain A → {chain_pdb.name}")
    extract_chain_a(TARGET_CIF, chain_pdb)

    log("PDBFixer: adding missing atoms + protonating pH 7.4")
    fixer = pdbfixer.PDBFixer(filename=str(chain_pdb))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    prepped = SNAPSHOTS_DIR / "mutant_chainA_prepped.pdb"
    with open(prepped, "w") as fh:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    log(f"Prepped PDB: {prepped} ({len(fixer.positions)} atoms)")

    log("Building AMBER14SB + TIP3P system, 1 nm pad, 0.15 M NaCl")
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(forcefield,
                         padding=1.0 * unit.nanometer,
                         ionicStrength=0.15 * unit.molar,
                         neutralize=True,
                         model="tip3p")
    n_atoms = modeller.topology.getNumAtoms()
    log(f"Solvated system: {n_atoms} atoms")

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=5e-4,
    )

    temperature = 310 * unit.kelvin
    dt = 2 * unit.femtoseconds
    integrator = LangevinIntegrator(temperature, 1.0 / unit.picoseconds, dt)

    platform = None
    if platform_name:
        try:
            platform = Platform.getPlatformByName(platform_name)
        except Exception as e:
            log(f"Requested platform {platform_name} not available: {e}")
    if platform is None:
        for name in ("CUDA", "OpenCL", "CPU"):
            try:
                platform = Platform.getPlatformByName(name)
                log(f"Using OpenMM platform: {name}")
                break
            except Exception:
                continue

    simulation = app.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    log("Minimizing (5000 steps)")
    simulation.minimizeEnergy(maxIterations=5000)

    log("NVT 50 ps equilibration")
    simulation.context.setVelocitiesToTemperature(temperature)
    t0 = time.time()
    simulation.step(25_000)
    t_nvt = time.time() - t0

    log("NPT 50 ps equilibration")
    barostat = MonteCarloBarostat(1.0 * unit.bar, temperature, 25)
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)
    t0 = time.time()
    simulation.step(25_000)
    t_npt = time.time() - t0

    total_steps = int(production_ps / 2.0 * 1000)
    step_per_snapshot = max(1, total_steps // n_snapshots)
    log(f"Production: {production_ps} ps = {total_steps} steps, "
        f"snapshot every {step_per_snapshot} steps")

    snapshot_paths = []
    t0 = time.time()
    for i in range(n_snapshots):
        simulation.step(step_per_snapshot)
        state = simulation.context.getState(getPositions=True)
        positions = state.getPositions()
        snap_path = SNAPSHOTS_DIR / f"snap_{i:03d}.pdb"
        with open(snap_path, "w") as fh:
            app.PDBFile.writeFile(simulation.topology, positions, fh, keepIds=True)
        snapshot_paths.append(str(snap_path))
        log(f"  snapshot {i+1}/{n_snapshots} saved")
    t_prod = time.time() - t0

    ns_per_day = (production_ps / 1000.0) / (t_prod / 86400.0)

    return {
        "n_atoms_solvated": n_atoms,
        "production_ps": production_ps,
        "n_snapshots": n_snapshots,
        "timing_sec": {
            "nvt_equil_50ps": round(t_nvt, 1),
            "npt_equil_50ps": round(t_npt, 1),
            "production": round(t_prod, 1),
        },
        "ns_per_day_extrap": round(ns_per_day, 2),
        "production_10ns_eta_hours": round(10.0 / ns_per_day * 24, 1),
        "snapshot_pdb_paths": snapshot_paths,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=float, default=0.1,
                    help="Production duration in nanoseconds (SMOKE default 0.1)")
    ap.add_argument("--snapshots", type=int, default=5,
                    help="Number of snapshots to save during production")
    ap.add_argument("--platform", default=None,
                    help="Force OpenMM platform (CUDA, OpenCL, CPU)")
    args = ap.parse_args()

    env_check()
    mut_info = verify_mutant_cif(TARGET_CIF)

    production_ps = args.ns * 1000.0
    t0 = time.time()
    result = run_md(production_ps, args.snapshots, args.platform)
    total_elapsed = time.time() - t0

    summary = {
        "batch": "pharmacochaperone_phase5d_e1659a_md",
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "target_cif": str(TARGET_CIF),
        "chain_used": "A (full-length E1659A STRC mutant)",
        "mutation_verified": {"K1141": "LYS", "E1659": "ALA"},
        "residue_range": mut_info,
        "force_field": "AMBER14SB + TIP3P, PME, HBond constraints, 2 fs step "
                       "(identical to Phase 5a for comparability)",
        "temperature_K": 310,
        "pressure_bar": 1.0,
        "ionic_strength_M": 0.15,
        "run_args": vars(args),
        "result": result,
        "total_elapsed_sec": round(total_elapsed, 1),
        "interpretation": (
            f"ns/day extrapolation: {result['ns_per_day_extrap']:.2f}. "
            f"System is {result['n_atoms_solvated']} atoms vs Phase 5a ~164k. "
            "Phase 5e will run Vina ensemble re-dock on these snapshots to "
            "compare f_PC(mut) vs f_PC(WT) from Phase 5b."
        ),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    log(f"JSON: {OUT_JSON}")
    log(f"Snapshots: {SNAPSHOTS_DIR}")
    log(f"Total wall time: {total_elapsed:.1f} s")
    log(f"ns/day: {result['ns_per_day_extrap']:.2f}")


if __name__ == "__main__":
    main()
