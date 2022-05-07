#!/usr/bin/python

import numpy as np
import os as os
import argparse

from openeye.oechem import *
from openeye.oedepict import *


class charges(object):
    '''class that performs different things with charges
       and stores all the relevant information in a dictionary
    '''

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def read(cls, filename='charges.out'):
        '''reads the time step and the values of all the charges
           into a dictionary "charges" which has the charge indices as a keys
           '''
        charges = {}
        indices = []
        if os.path.exists(filename):
            f = open(filename, 'r')
        else:
            print(f"You have to provide a file with the charges")
            exit(1)
        for line in f:
            words = line.rsplit()
            if words[0][0] == "#":
                first_line = next(f)
                words = first_line.rsplit()
                for i in range(len(words)):
                    indices.append(str(i))
                    charges.update({str(i): [float(words[i])]})
            else:
                for i in range(len(words)):
                    charges[indices[i]].append(float(words[i]))
        f.close()
        return cls(**charges)

    def symmetrize(self, mol):
        '''Function to symmetrize charges
        '''
        charges = []
        data = vars(self).copy()
        for i in data.keys():
            # The last charge array has the converged values
            charges.append(data[i][-1])
        if not len(charges) == mol.NumAtoms():
            raise Exception('Error: Number of charges and atoms do not match')
        for (idx, atom) in enumerate(mol.GetAtoms()):
            atom.SetPartialCharge(charges[idx])

        OEPerceiveSymmetry(mol)
        symmetry_charges = {}
        for atom in mol.GetAtoms():
            sym_class = atom.GetSymmetryClass()
            if sym_class not in symmetry_charges:
                symmetry_charges[sym_class] = [atom.GetPartialCharge()]
            else:
                symmetry_charges[sym_class].append(atom.GetPartialCharge())
        for sym_class in symmetry_charges:
            symmetry_charges[sym_class] = np.mean(symmetry_charges[sym_class])
        # replace charges
        for atom in mol.GetAtoms():
            sym_class = atom.GetSymmetryClass()
            charge = symmetry_charges[sym_class]
            print(charge)
            atom.SetPartialCharge(charge)
        return mol


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('mol2_file', help='mol2 file')
    args = parser.parse_args()
    return args.mol2_file


def main():
    # Here starts the actual commands for processing the charges
    # Read in charges from the file "charges.out"
    mol2_file = parse_args()
    charge1 = charges.read()
    flavor = oechem.OEIFlavor_MOL2_Forcefield
    mol2_orig = mol2_file
    istream = oemolistream(mol2_orig)
    istream.SetFlavor(oechem.OEFormat_MOL2, flavor)
    mol = OEMol()
    OEReadMolecule(istream, mol)
    istream.close()
    mol2 = charge1.symmetrize(mol)
    tot_chg = 0.
    for atom in mol2.GetAtoms():
        tot_chg += atom.GetPartialCharge()
    print(f'Total charge: {tot_chg:1.3f}')
    ostream = oemolostream('mbis_' + mol2_orig)
    flavor2 = oechem.OEOFlavor_MOL2_Forcefield
    ostream.SetFlavor(oechem.OEFormat_MOL2, flavor2)
    OEWriteMolecule(ostream, mol2)


if __name__ == '__main__':
    main()
