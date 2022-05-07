import argparse

import parmed as pmd
import numpy as np

from parmed.tools import change, addLJType


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('top_file', help='topology file')
    parser.add_argument('lig_sel', help='ligand selection mask (amber)')
    parser.add_argument('--mol2_file', action='store', dest='mol2_file',
                        help='mol2 file with symmetrized charges')
    parser.add_argument('--lj_table_file', action='store', dest='lj_table_file',
                        help='LJ parameters table file (rmin and epsilon)')
    parser.add_argument('--output', action='store', dest='output',
                        help='output filename')
    args = parser.parse_args()
    return args.top_file, args.lig_sel, args.mol2_file, args.lj_table_file, args.output


def get_params(mol2_file=None, lj_table_file=None, charge=False, lj=False):
    params = dict()
    if mol2_file is not None:
        lig = pmd.load_file(mol2_file)
        params['charge'] = np.array([atom.charge for atom in lig.atoms])
    if lj_table_file is not None:
        lj_table = np.loadtxt(lj_table_file)
        rmin, eps = lj_table
        params['lj'] = rmin, eps
    print(params)
    return params


def write_params(top_file, lig_sel=':1', params=None, output=None):
    if params is not None:
        top = make_new_top(top_file, lig_sel, params, output)
        return top


def make_new_top(top_file, lig_sel, params, output=None):
    _top = pmd.load_file(top_file)
    for key in params.keys():
        for i, atom in enumerate(_top[lig_sel].atoms):
            mask = f'{lig_sel}&@{atom.name}'
            if key is 'charge':
                action = change(_top, mask, 'charge', round(params[key][i], 3))
            if key is 'lj':
                action = addLJType(_top, mask, 'radius', round(
                    params[key][0][i] / 2, 3), 'epsilon', round(params[key][1][i], 6))
            action.execute()
    out = f'mbis_{top_file}' if output is None else output
    _top.save(out, overwrite=True)
    return _top


if __name__ == '__main__':
    top_file, lig_sel, mol2_file, lj_table_file, output = parse_args()
    params = get_params(mol2_file, lj_table_file)
    if bool(params):
        write_params(top_file, lig_sel, params, output)
