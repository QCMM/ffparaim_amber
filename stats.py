#!/usr/bin/python3
import numpy as np

def charge_stats(datas):

    charges = np.array([data.atffparams['charges'] for data in datas])
    charges_mean, charges_std = charges.mean(
        axis=0, dtype=np.float64), charges.std(axis=0, dtype=np.float64)
    return charges_mean, charges_std


def epol_stats(datas):

    epol = np.array([data.extra['epol'] for data in datas])
    epol_mean, epol_std = epol.mean(), epol.std()
    return epol_mean, epol_std
