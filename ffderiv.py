#!/usr/bin/python3

from iodata import load_one
from denspart.adapters.horton3 import prepare_input
from denspart.mbis import partition
from denspart.properties import compute_rcubed
from ffparaim_amber import units


import numpy as np

# catch warnings
np.seterr(all='warn')


class ForceFieldDerivation(object):
    """docstring for ForceFieldDerivation."""

    def __init__(self, infile):
        self.data = load_one(infile)  # data from the orca mkl file
        self.grid = None
        self.rho = None
        self.pro_model = None
        self.localgrids = None

    # set molgrid
    def set_molgrid(self, nrad, nang):
        self.grid, self.rho = prepare_input(self.data, nrad, nang, 10000)  # molgrid
        return

    # do partitioning depending on method
    def do_partitioning(self, method='mbis'):
        if method == 'hi':
            print('Pending')
            return
        elif method == 'mbis':
            # do MBIS partitioning
            print('MBIS partitioning ...')
            self.pro_model, self.localgrids = partition(
                self.data.atnums, self.data.atcoords, self.grid, self.rho)
        else:
            print('Invalid method')
            return
        return

    def get_charges(self):
        return self.pro_model.charges.tolist()

    def get_epol(self):
        orcalog = load_one('orca_pol_corr.out', fmt="orcalog")
        epol = (
            orcalog.extra["scf_energies"][0] - orcalog.extra["scf_energies"][-1]) * units.E_AU
        return epol

    def get_rcubed(self):
        rc = compute_rcubed(self.pro_model, self.grid, self.rho, self.localgrids)
        return rc.tolist()
