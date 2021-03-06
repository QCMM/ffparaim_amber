#!/usr/bin/python3

import parmed as pmd

from yank.restraints import Harmonic, FlatBottom, Boresch, RestraintParameterError
from openmmtools import states
from yank.yank import Topography
from simtk import unit


def get_ligand_atom_list(top_file,
                         coords_file,
                         ligand_selection):

    _top = pmd.load_file(top_file, xyz=coords_file)
    atom_list = _top.atoms
    lig = _top[ligand_selection].atoms
    ligand_atom_list = [idx for at in lig for idx, atom in enumerate(atom_list) if (
        at.residue.number, at.name) == (atom.residue.number, atom.name)]
    return ligand_atom_list


def set_restraints(top,
                   system,
                   positions,
                   restraint_dict=None,
                   ligand_atom_list=None):
    print('Setting restraints ...')
    if restraint_dict is not None:
        thermodynamic_state = states.ThermodynamicState(system, 298 * unit.kelvin)
        sampler_state = states.SamplerState(positions)
        if ligand_atom_list is not None:
            topography = Topography(top.topology, ligand_atoms=ligand_atom_list)
        if 'harmonic' in restraint_dict.keys():
            restraint = Harmonic(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                 restrained_receptor_atoms=topography.select(restraint_dict['harmonic']['receptor']),
                                 restrained_ligand_atoms=topography.select(restraint_dict['harmonic']['ligand']))
        if 'flatbottom' in restraint_dict.keys():
            restraint = FlatBottom(spring_constant=0.2 * unit.kilocalories_per_mole / unit.angstrom**2,
                                   well_radius=10.0 / unit.angstrom,
                                   restrained_receptor_atoms=topography.select(restraint_dict['flatbottom']['receptor']),
                                   restrained_ligand_atoms=topography.select(restraint_dict['flatbottom']['ligand']))

        if 'boresch' in restraint_dict.keys():
            restraint = Boresch(restrained_receptor_atoms=topography.select(restraint_dict['boresch']['receptor']),
                                restrained_ligand_atoms=topography.select(restraint_dict['boresch']['ligand']),
                                K_r=20.0 * unit.kilocalories_per_mole / unit.angstrom ** 2,
                                r_aA0=0.35 * unit.nanometer)
        try:
            restraint.restrain_state(thermodynamic_state)
        except RestraintParameterError:
            print('Choosing restraint parameters automatically.')
            restraint.determine_missing_parameters(thermodynamic_state, sampler_state, topography)
            restraint.restrain_state(thermodynamic_state)

    return
