#!/usr/bin/python3
import warnings
import mdtraj
import sys

import parmed as pmd

from restraints import set_restraints
from parmed.tools import change

from simtk.openmm import openmm
from simtk.openmm import app
from simtk import unit

# avoid warnings
warnings.filterwarnings('ignore')


def setup_simulation(top,
                     positions,
                     update,
                     box_vectors=None,
                     restraint_dict=None,
                     ligand_atom_list=None):
    '''Setup the openMM system with the current topology and
    the input coordinates or the current positions depending on
    the value of update.
    Standard conditions are assumed (298K, 1bar)
    Input:
    top : Topology object from OpenMM o ParmEd (Gromacs or Amber)
    positions: current positions of atoms
    update: integer of charge update cycle
    Returns:
    Simulation (OpenMM class)
    '''
    system = top.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=1 * unit.nanometer, constraints=app.HBonds)
    system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 298 * unit.kelvin))
    if restraint_dict is not None:
        set_restraints(top,
                       system,
                       positions,
                       restraint_dict,
                       ligand_atom_list)
    integrator = openmm.LangevinIntegrator(
        298 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(top.topology, system, integrator)
    simulation.reporters.append(app.StateDataReporter(
        sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True, density=True))
    simulation.reporters.append(app.DCDReporter(f'traj_{update}.dcd', 50000))
    simulation.context.setPositions(positions)
    if box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*box_vectors)
    simulation.minimizeEnergy()
    return simulation, system


def get_positions(simulation, system):
    positions = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=True).getPositions()
    _box_vectors = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(_box_vectors)
    app.PDBFile.writeFile(simulation.topology, positions,
                          open('output.pdb', 'w'))
    return positions


def image_molecule(pdbfile='output.pdb'):
    traj = mdtraj.load(pdbfile)
    traj.image_molecules()
    pos = mdtraj.utils.in_units_of(traj.xyz[0], traj._distance_unit, "angstroms")
    ucl = mdtraj.utils.in_units_of(traj.unitcell_lengths[0], traj._distance_unit, "angstroms")
    pdb_recentered = mdtraj.formats.PDBTrajectoryFile("output_recenter.pdb", "w")
    pdb_recentered.write(pos, traj.top, modelIndex=None, unitcell_lengths=ucl,
                         unitcell_angles=traj.unitcell_angles[0])


def make_new_top(top_file, box_vectors, charges_mean, ligand_selection):

    _top = pmd.load_file(top_file)
    _top.box_vectors = box_vectors
    for i, atom in enumerate(_top[ligand_selection].atoms):
        mask = f'{ligand_selection}&@{atom.name}'
        action = change(_top, mask, 'charge', round(charges_mean[i], 6))
        action.execute()
    _top.save(top_file, overwrite=True)
    return _top
