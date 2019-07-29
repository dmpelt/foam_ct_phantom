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

import h5py

FILE_VERSION=1.0


def load_volume(fn):  
    with h5py.File(fn,'r') as f:
        vol = f['volume'][:]
    return vol

def load_projections(fn):
    with h5py.File(fn,'r') as f:
        projs = f['projs'][:]
    return projs