#!/usr/bin/python3

import sys

from ffparaim import FFparAIM

top_file = sys.argv[1]
coords_file = sys.argv[2]


def main():
    nb_params = FFparAIM(top_file,
                         coords_file,
                         qm_charge=0,
                         ligand_selection=':1',
                         n_updates=2,
                         sampling_time=1,
                         total_qm_calculations=5)
    nb_params.run(json=True)


if __name__ == '__main__':
    main()
