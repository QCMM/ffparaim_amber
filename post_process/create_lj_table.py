#!/usr/bin/python

import json
import numpy as np

if __name__ == '__main__':
    with open('nbff.json', 'r') as f:
        datas = json.load(f)
    rmin = [datas[atom]['rmin'] for atom in datas.keys()]
    eps = [datas[atom]['eps'] for atom in datas.keys()]
    lj_table = np.array([rmin, eps])
    np.savetxt('lj_table.dat', lj_table)

