#-----------------------------------------------------------------------
#Copyright 2019 Centrum Wiskunde & Informatica, Amsterdam
#
#Author: Daniel M. Pelt
#Contact: D.M.Pelt@cwi.nl
#Website: http://dmpelt.github.io/foam_ct_phantom/
#License: MIT
#
#This file is part of foam_ct_phantom, a Python package for generating
#foam-like phantoms for CT.
#-----------------------------------------------------------------------

from .utils import FILE_VERSION
from . import ccode
import h5py
import numpy as np
import scipy.optimize as sco
import tqdm

def apply_poisson_noise(input_file, output_file, seed, flux, absorption_factor):
    ccode.setseed(seed)
    with h5py.File(input_file, 'r') as f:
        projs = f['projs']
        spatt = dict(f['projs'].attrs)
        with h5py.File(output_file, 'w') as fo:
            fo.attrs['FILE_VERSION'] = FILE_VERSION
            dset = fo.create_dataset('projs', projs.shape, dtype='f4')
            for i in tqdm.trange(projs.shape[0]):
                curproj = projs[i][:]
                ccode.applypoisson(curproj,flux,absorption_factor)
                dset[i] = curproj
            att = fo['projs'].attrs
            for key, val in spatt.items():
                att[key] = val
            att['seed'] = seed
            att['flux'] = flux
            att['absorption_factor'] = absorption_factor
            if 'times' in f:
                fo['times'] = f['times'][:]

def estimate_absorption_factor(input_file, average_absorption_ratio):
    with h5py.File(input_file, 'r') as f:
        projs = f['projs'][:]
    projs = projs[projs>0]

    def error_function(fac):
        rat = 1 - np.exp(-projs*fac[0]).mean()
        return (rat - average_absorption_ratio)**2
    
    return sco.fmin(error_function,1)[0]
    