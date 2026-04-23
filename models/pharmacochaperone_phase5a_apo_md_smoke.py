#!/usr/bin/env python3
"""
Phase 5a — Apo-receptor short MD to produce snapshot ensemble for re-docking.

WHY THIS, NOT FULL MM-GBSA
══════════════════════════
Phase 4f attempted single-point interface rescue via implicit GBn2 and produced
numerically unusable absolute energies (Phase 4f JSON shows +20,440 kcal/mol
"rescue" — method-inadequate, documented in scripts inventory).

Phase 5 proper MD + MM-GBSA on ligand-bound complexes needs AmberTools / OpenFF
for ligand parameterisation, a conda env with GAFF2, and overnight run on 5
ligands × 3 replicates × 10 ns. That is Phase 5b scope.

Phase 5a (THIS script) is a narrower scope: short MD on the APO receptor
only, extract N snapshots, save as PDB. Phase 5b will take those snapshots
and re-dock the Phase 4b shortlist to each, computing ensemble-averaged
Vina/Gnina ΔG_bind. Ensemble docking is the standard way to reduce
single-structure lucky-pose noise without needing full MM-GBSA.

Deliverable for Misha cure argument:
    Phase 5a (this)     — MD snapshots saved
    Phase 5b (follow-up) — ensemble re-docking → improved ΔG_bind per lead
    → Kd estimate with ±1 kcal/mol precision
    → bound-fraction at [L] = 10 μM therapeutic intracochlear conc
    → f_PC point estimate for each lead
    → feeds into [[Misha Compound-Het Therapy Stack Model]]

═══════════════════════════════════════════════════════════════════════════
METHOD
═══════════════════════════════════════════════════════════════════════════

1. Source: Ultra-Mini × TMEM145 CIF, chain A only (Ultra-Mini receptor).
2. Prep with PDBFixer: add missing atoms/H, protonate at pH 7.4.
3. Force field: AMBER14SB + TIP3P water, NaCl 0.15 M neutralisation.
4. Box: cubic, 1 nm pad around solute.
5. Minimisation: 5000 steps.
6. NVT: 50 ps @ 310 K, Langevin thermostat, 2 fs timestep.
7. NPT: 50 ps @ 1 bar Monte Carlo barostat.
8. Production NPT: configurable (SMOKE 100 ps; production 5-10 ns).
9. Save N snapshots (receptor chain only, no water/ions) as PDB.

═══════════════════════════════════════════════════════════════════════════
RUN MODES
═══════════════════════════════════════════════════════════════════════════

    python3.11 phase5a_apo_md_smoke.py                 # SMOKE: 100 ps / 5 snapshots
    python3.11 phase5a_apo_md_smoke.py --ns 5          # production: 5 ns / 50 snapshots
    python3.11 phase5a_apo_md_smoke.py --ns 10 --snapshots 100

SMOKE goal: prove the pipeline runs end-to-end on this machine and measure
ns/day, so we know whether production is feasible locally vs cloud.

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
TARGET_CIF = Path("/Users/egorlyfar/Sites/site-strc-egor-lol/public/models/job-ultramini-x-tmem145-full.cif")
OUT_JSON = WORK / "pharmacochaperone_phase5a_apo_md_smoke.json"
SNAPSHOTS_DIR = WORK / "artifacts" / "phase5a_snapshots"


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
        log("Missing dep. Install: pip3.11 install --user openmm pdbfixer")
        sys.exit(2)
    log(f"openmm {openmm.version.version} OK")


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
    struct = parser.get_structure("ultramini", str(cif_path))
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out_pdb), select=ChainAStandard())


def run_md(production_ps: float, n_snapshots: int,
           platform_name: str | None = None) -> dict:
    """
    End-to-end: prep → minimize → NVT → NPT → production → snapshots.
    Returns timing info and snapshot list.
    """
    import openmm
    from openmm import app, unit, Platform
    from openmm import LangevinIntegrator, MonteCarloBarostat
    import pdbfixer

    SNAPSHOTS_DIR.mkdir(parents=True, exist_ok=True)
    chain_pdb = SNAPSHOTS_DIR / "ultramini_chainA.pdb"
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

    prepped = SNAPSHOTS_DIR / "ultramini_chainA_prepped.pdb"
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
        # Prefer CUDA > OpenCL > CPU
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

    # NVT equilibration 50 ps
    log("NVT 50 ps equilibration")
    simulation.context.setVelocitiesToTemperature(temperature)
    t0 = time.time()
    simulation.step(25_000)   # 50 ps @ 2 fs
    t_nvt = time.time() - t0

    # NPT equilibration 50 ps
    log("NPT 50 ps equilibration")
    barostat = MonteCarloBarostat(1.0 * unit.bar, temperature, 25)
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)
    t0 = time.time()
    simulation.step(25_000)
    t_npt = time.time() - t0

    # Production
    total_steps = int(production_ps / 2.0 * 1000)   # fs → 2 fs step
    step_per_snapshot = max(1, total_steps // n_snapshots)
    log(f"Production: {production_ps} ps = {total_steps} steps, "
        f"snapshot every {step_per_snapshot} steps")

    snapshot_paths = []
    t0 = time.time()
    for i in range(n_snapshots):
        simulation.step(step_per_snapshot)
        state = simulation.context.getState(getPositions=True)
        positions = state.getPositions()
        # Save full system; downstream re-docking strips water/ions
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
    production_ps = args.ns * 1000.0
    t0 = time.time()
    result = run_md(production_ps, args.snapshots, args.platform)
    total_elapsed = time.time() - t0

    summary = {
        "batch": "pharmacochaperone_phase5a_apo_md_smoke",
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "target_cif": str(TARGET_CIF),
        "chain_used": "A (Ultra-Mini receptor)",
        "force_field": "AMBER14SB + TIP3P, PME, HBond constraints, 2 fs step",
        "temperature_K": 310,
        "pressure_bar": 1.0,
        "ionic_strength_M": 0.15,
        "run_args": vars(args),
        "result": result,
        "total_elapsed_sec": round(total_elapsed, 1),
        "interpretation": (
            f"ns/day extrapolation: {result['ns_per_day_extrap']:.2f}. "
            f"10 ns production would take ~{result['production_10ns_eta_hours']:.1f} h. "
            "If >24 h on local → cloud compute recommended for Phase 5b."
        ),
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    log(f"JSON: {OUT_JSON}")
    log(f"Snapshots: {SNAPSHOTS_DIR}")
    log(f"Total wall time: {total_elapsed:.1f} s")
    log(f"ns/day: {result['ns_per_day_extrap']:.2f} → "
        f"10 ns ETA ~{result['production_10ns_eta_hours']:.1f} h")


if __name__ == "__main__":
    main()
