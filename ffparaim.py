#!/usr/bin/python3
import time
import os
import shutil

import ffparaim_amber.mdtools as mdt
import ffparaim_amber.qmtools as qmt
import parmed as pmd

from ffparaim_amber.ffderiv import ForceFieldDerivation
from ffparaim_amber.orcaff import OrcaForceField
from ffparaim_amber import stats
from ffparaim_amber import io
from ffparaim_amber import utils
from simtk.openmm import app
from iodata import IOData


class FFparAIM(object):
    """docstring for ffparaim."""

    def __init__(self, top_file, coords_file, qm_charge=0, ligand_selection=':1', receptor_selection=None,
                 n_charge_updates=3, sampling_time=25, total_qm_calculations=100, method='B3LYP',
                 basis='def2-TZVP'):

        self.top_file = top_file
        self.coords_file = coords_file
        self.qm_charge = qm_charge  # total qm charge
        # residue index for molecule to calculate charges
        self.ligand_selection = ligand_selection
        if receptor_selection is not None:
            self.receptor_selection = receptor_selection
        self.n_charge_updates = n_charge_updates
        self.sampling_time = sampling_time  # sampling time in ns
        self.total_qm_calculations = total_qm_calculations
        self.method = method
        self.basis = basis
        self.data = dict()
        # format-dependent topology attributes
        if self.top_file.endswith('.top'):
            self.top_format = 'gromacs'
            self._gromacs_top_path = f'{shutil.which("gmx").split("bin")[0]}share/gromacs/top'
            pmd.gromacs.GROMACS_TOPDIR = self._gromacs_top_path
            self.coords = app.gromacsgrofile.GromacsGroFile(self.coords_file)
            self.box_vectors = self.coords.getPeriodicBoxVectors()
            self.top = app.gromacstopfile.GromacsTopFile(
                self.top_file, periodicBoxVectors=self.box_vectors, includeDir=self._gromacs_top_path)
        elif self.top_file.endswith('.prmtop'):
            self.top_format = 'amber'
            self.coords = app.amberinpcrdfile.AmberInpcrdFile(self.coords_file)
            self.box_vectors = self.coords.boxVectors
            self.top = app.amberprmtopfile.AmberPrmtopFile(self.top_file)
        else:
            raise Exception(
                'Topology not implemented ...')

    def run(self, compl=False, output=True, json=False):

        begin_time = time.time()
        mm_time = 0
        qm_time = 0
        deriv_time = 0
        molgrid_time = 0
        part_time = 0
        # create topology file backup
        shutil.copyfile(self.top_file, f'{self.top_file}.old')
        # apply restraints for complex simulations
        if compl:
            print('Setting restraints ...')
            restraint = True
            ligand_atom_list, receptor_atom_list = mdt.set_restrained_atoms(
                self.top_file, self.coords_file, self.ligand_selection, self.receptor_selection)
        else:
            restraint = False
            ligand_atom_list = None
            receptor_atom_list = None
        # get box info
        b_vectors = self.box_vectors if self.top_format is 'amber' else None
        for update in range(self.n_charge_updates):
            # list to store charges and polarization energies for each update
            self.data[update] = list()
            # creating OpenMM simulation class
            if update == 0:
                positions = self.coords.positions
            simulation, system = mdt.setup_simulation(
                self.top, positions, update, b_vectors, restraint, ligand_atom_list, receptor_atom_list)
            # write Orca forcefield file
            ff = OrcaForceField(self.top_file, self.coords_file)
            prms = ff.parse_prms()
            ff.write_prmsfile(prms)
            qm_calculations = int(
                self.total_qm_calculations / self.n_charge_updates) * (update + 1)
            # starting loop to calculate atomic charges from different conformations
            for i in range(qm_calculations):
                step = int(self.sampling_time * 500000 / qm_calculations)
                mm_begin = time.time()
                simulation.step(step)
                mm_end = time.time()
                # calculate charges and polarization energy for current configuration
                print('Calculating charges ...')
                positions = mdt.get_positions(simulation, system)
                mdt.image_molecule()
                qm_region = qmt.set_qm_atoms(self.ligand_selection)
                qmt.write_qmmm_pdb(qm_region)
                for inp in ('qmmm', 'pol_corr'):
                    lig = self.ligand_selection if inp is 'pol_corr' else None
                    qmt.write_orca_input(inp, ligand_selection=lig, method=self.method,
                                         basis=self.basis, qm_charge=self.qm_charge)
                qm_begin = time.time()
                qmt.exec_orca()
                qm_end = time.time()
                # parameter derivation
                deriv_begin = time.time()
                ffderiv = ForceFieldDerivation('orca_qmmm.molden.input')
                molgrid_begin = time.time()
                ffderiv.set_molgrid(75, 110)
                molgrid_end = time.time()
                part_begin = time.time()
                ffderiv.do_partitioning(method='mbis')
                part_end = time.time()
                # get data
                charges = ffderiv.get_charges()
                epol = ffderiv.get_epol()
                rcubed = ffderiv.get_rcubed()
                # store data
                self.data[update].append(IOData(atffparams={'charges': charges, 'rcubed': rcubed},
                                                extra={'epol': epol}))
                deriv_end = time.time()
                # add time exec
                mm_time += utils.get_time(mm_begin, mm_end)
                qm_time += utils.get_time(qm_begin, qm_end)
                molgrid_time += utils.get_time(molgrid_begin, molgrid_end)
                part_time += utils.get_time(part_begin, part_end)
                deriv_time += utils.get_time(deriv_begin, deriv_end)
            # update topology
            print('Creating new topology ...')
            new_charges = stats.charge_stats(self.data[update])[0]
            self.top = mdt.make_new_top(
                self.top_file, self.box_vectors, new_charges, self.ligand_selection)
        if output:
            io.write_output(self.data)
        if json:
            io.write_json(self.data)
        end_time = time.time()
        total_time = utils.get_time(begin_time, end_time)
        extras_time = total_time - (mm_time + qm_time + deriv_time)
        print('Summary:')
        print(f'Accumulated time in MM calculations: {round(mm_time, 2)} hours')
        print(f'Accumulated time in QM calculations: {round(qm_time, 2)} hours')
        print(f'Accumulated time in parameter derivation: {round(deriv_time, 2)} hours')
        print(f'\tAccumulated time in grid setting: {round(molgrid_time, 2)} hours')
        print(f'\tAccumulated time in partitioning: {round(part_time, 2)} hours')
        print(f'Accumulated time in extra functions: {round(extras_time, 2)} hours')
        print(f'Total time: {round(total_time, 2)} hours')
        return

    def validation(self, parm=None, parm_vals=[], overwrite=False, compl=False):

        parm_opt = ['sampling_time',
                    'n_charge_updates',
                    'total_qm_calculations',
                    'method',
                    'basis']
        if str(parm) in parm_opt and isinstance(parm_vals, list):
            parm_dict = dict({parm: parm_vals})
        else:
            raise Exception(
                '''Usage:
                        parm: str ('sampling_time',
                        'n_charge_updates',
                        'total_qm_calculations',
                        'method',
                        'basis')
                        parm_vals: list
                        ''')
        for val in parm_dict[parm]:
            if parm == 'sampling_time':
                self.sampling_time = float(val)
            elif parm == 'total_qm_calculations':
                self.total_qm_calculations = int(val)
            elif parm == 'method':
                self.method = str(val)
            elif parm == 'basis':
                self.basis = str(val)
            else:
                self.n_charge_updates = int(val)
            parm_dir = f'{parm[0]}_{val}'
            io.create_parm_dir(self.top_file, self.coords_file, parm_dir, overwrite)
            os.chdir(parm_dir)
            self.run(compl)
            os.chdir('..')
        return
