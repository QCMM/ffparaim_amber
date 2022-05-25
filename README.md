# ffparaim_amber
Derivation of non-bonded non-polarizable force field parameters from Atom in Molecule properties. Compatibility with Amber and Gromacs topologies. Legacy code.

## Dependencies

- YANK 0.25.2
- ParmEd 3.4.3
- Horton3 packages (IOData, DensPart, Grid, GBasis) https://qcdevs.org/
- ORCA 4.2.1

## Usage

```python
from ffparaim import FFparAIM

top_file = 'system.prmtop'
coords_file = 'system.inpcrd'

nb_params = FFparAIM(top_file,
                     coords_file,
                     qm_charge=0,
                     ligand_selection=':1',
                     n_updates=2,
                     sampling_time=1,
                     total_qm_calculations=5)
nb_params.run(json=True) 
```
## Post-processing

Charge symmetrization script requieres an OpenEye license.

`python symmetrize_charges.py lig.mol2`

`python extract_nbff.py system.prmtop ':1'`

`python create_lj_table.py `

`python write_nbff.py system.prmtop ':1' --mol2_file mbis_lig.mol2 --lj_table_file lj_table.dat --output system_mbis_lj_charges.prmtop`
