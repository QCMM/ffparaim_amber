#!/usr/bin/python3
"""Helper variables used for Orca 4.2.1 enviroment execution."""
import os

# List of elements according to atomic number.
elements = [
    None,
    'H', 'He',
    'Li', 'Be',
    'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
    'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub'
]

# Default template for Orca forcefield.
orcaff_template = """\
$$atoms
${atoms_block}
$$bonds
${bonds_block}
$$angles
${angles_block}
$$dihedrals
${dihedrals_block}\
"""

# String format for Orca forcefield blocks.
orcaff_blocks = {
    'header': '{0:1d} {1:1d} {2:1d}\n',
    'atoms': '{0:6d}   {1:2s}   {2:10.6f}   {3:10.6f}   {4:10.6f}   {5:10.6f}   {6:10.6f}',
    'bonds': '{0:6d} {1:10d} {2:10.6f} {3:10.6f}',
    'angles': '{0:6d} {1:10d} {2:10.6f} {3:10.6f}',
    'dihedrals': '{0:6d}   {1:6d}   {2:6d}   {3:6d}   {4:10.6f}   {5:10.6f}   {6:6d}'
}

# Default template for Orca QM/MM input.
orca_qmmm_template = """\
! QMMM ${method} ${basis} Grid4 TightSCF NOFINALGRID KeepDens
%output PrintLevel Mini Print[ P_Mulliken ] 1 Print[P_AtCharges_M] 1 end
%pal nprocs ${nproc} end
%qmmm
    ORCAFFFilename "ORCAFF.prms"
    Use_QM_InfoFromPDB true     # get QM atoms from pdb file
end
*pdbfile ${qm_charge} ${qm_mult} output_qmmm.pdb
"""

# Default template for Orca polarization correction input.
orca_pol_corr_template = """\
! ${method} ${basis} Grid4 TightSCF NOFINALGRID KeepDens
%output PrintLevel Mini Print[ P_Mulliken ] 1 Print[P_AtCharges_M] 1 end
%pal nprocs ${nproc} end
%coords
    CTyp xyz
    Charge ${qm_charge}
    Mult ${qm_mult}
    Units Angs
    coords
${geometry}
    end
end"""


def get_nproc():
    """Get the number of processes for QM calculation."""
    if 'OMP_NUM_THREADS' in os.environ:
        nproc = int(os.environ['OMP_NUM_THREADS'])
    elif 'SLURM_NTASKS' in os.environ:
        nproc = int(os.environ['SLURM_NTASKS']) - 4
    else:
        nproc = 1
    return nproc


def get_time(begin_time, end_time):
    return (end_time - begin_time) / 3600
