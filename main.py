#!/usr/bin/python3
import sys

from ffparaim_amber import FFparAIM

top_file = sys.argv[1]
coords_file = sys.argv[2]


def main():
    charges_param = FFparAIM(top_file, coords_file, qm_charge=0, ligand_selection=f':1',
                             n_charge_updates=2, sampling_time=1, total_qm_calculations=5)
    charges_param.run(json=True)


if __name__ == '__main__':
    main()
