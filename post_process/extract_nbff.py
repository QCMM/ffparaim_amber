#!/usr/bin/env python

import json
import argparse

import numpy as np
import parmed as pmd

from iodata.utils import kcalmol, angstrom


# Values in atomic units for r^3 at B3LYP/def2-TZVP level of theory.
# CAUTION: These values were derived using ffderiv functions, implemented for MBIS partitioning.
# Temporary implementation.

rcubed_table = {
    1: 7.902,
    6: 34.935,
    7: 26.351,
    8: 21.969,
    16: 74.164,
    17: 64.324
}

# Values in atomic units, taken from:
# Chu, X., & Dalgarno, A. (2004).
# Linear response time-dependent density functional theory for van der Waals coefficients.
# The Journal of Chemical Physics, 121(9), 4083--4088. http://doi.org/10.1063/1.1779576


c6_table = {
    1: 6.5,
    3: 1393.0,
    6: 46.6,
    7: 24.2,
    8: 15.6,
    9: 9.5,
    10: 6.38,
    14: 305.0,
    16: 134.0,
    17: 94.6,
    18: 64.3,
    35: 162.0,
    36: 130.0,
}

alpha_table = {
    1: 4.5,
    3: 164.0,
    6: 12.0,
    7: 7.4,
    8: 5.4,
    9: 3.8,
    10: 2.67,
    14: 37.0,
    16: 19.6,
    17: 15.0,
    18: 11.1,
    35: 20.0,
    36: 16.7,
}


header = """\
Units
    q: e
    scale: 1
    alpha: a.u.
    A: kcal * angstrom ** 12 / mol
    B: kcal * angstrom ** 6 / mol
    eps: kcal / mol
    rmin: angstrom
-----------------------------------------------------------------------------
  i   Z       q    scale   alpha           A           B         eps    rmin
-----------------------------------------------------------------------------"""


def main():
    top_file, lig_sel = parse_args()
    top = pmd.load_file(top_file)
    atoms = [atom for atom in top[lig_sel].atoms]
    # 1) Get the volumes and charges for the atoms in the molecule.
    rcubed_aim = dict()
    charges_aim = dict()
    with open('ffparaim.json', 'r') as f:
        datas = json.load(f)
        for key in datas.keys():
            rcubed_aim[key] = np.array([data['atffparams']['rcubed'] for data in datas[key]])
            charges_aim[key] = np.array([data['atffparams']['charges'] for data in datas[key]])
    # 2) Compute FF parameters
    print(header)
    nbff = dict()
    for idx, atom in enumerate(atoms):
        charge = charges_aim[list(datas.keys())[-1]][:, idx].mean().round(3)
        vol_aim = rcubed_aim[list(datas.keys())[-1]][:, idx].mean().round(3)
        vol_isolated = rcubed_table[atom.atomic_number]
        scaling = vol_aim / vol_isolated
        # scaling /= (populations[iatom]/number)
        alpha = alpha_table[atom.atomic_number] * scaling
        c6 = c6_table[atom.atomic_number] * scaling ** 2
        radius = 2.54 * alpha ** (1.0 / 7.0)
        rmin = 2 * radius
        epsilon = c6 / (2 * rmin ** 6.0)
        A = epsilon * rmin ** 12.0
        print("{:3}  {:2d}  {:6.3f}  {:7.4f}  {:6.3f}  {:10.3f}  {:10.3f}  {:10.6f}  {:6.3f}".format(
            atom.name, atom.atomic_number, charge,
            scaling, alpha,
            A / (kcalmol * angstrom ** 12),
            c6 / (kcalmol * angstrom ** 6),
            epsilon / kcalmol,
            rmin / angstrom,
        ))
        nbff[atom.name] = {'Z': atom.atomic_number,
                           'q': charge,
                           'scale': scaling,
                           'A': A / (kcalmol * angstrom ** 12),
                           'B': c6 / (kcalmol * angstrom ** 6),
                           'eps': epsilon / kcalmol,
                           'rmin': rmin / angstrom
                           }
    with open('nbff.json', 'w') as json_out:
        json.dump(nbff, json_out)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('top_file', help='topology file')
    parser.add_argument('lig_sel', help='ligand selection mask (amber)')
    args = parser.parse_args()
    return args.top_file, args.lig_sel


if __name__ == '__main__':
    main()
