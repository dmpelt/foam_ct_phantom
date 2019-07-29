#-----------------------------------------------------------------------
#Copyright 2019 Centrum Wiskunde & Informatica, Amsterdam
#
#Author: Daniel M. Pelt
#Contact: D.M.Pelt@cwi.nl
#Website: http://cicwi.github.io/foam_ct_phantom/
#License: MIT
#
#This file is part of foam_ct_phantom, a Python package for generating
#foam-like phantoms for CT.
#-----------------------------------------------------------------------

import h5py
import numpy as np

from . import ccode

def single_par_projection(phantom, nx, ny, pixsize, angle, cx=0, cy=0, rotcx=0, rotcy=0, supersampling=1):
    if isinstance(phantom, str):
        with h5py.File(phantom, 'r') as f:
            spheres = f['spheres'][:]
    else:
        spheres = phantom
    proj = ccode.genparproj(spheres,supersampling*nx,supersampling*ny,pixsize/supersampling,angle,cx,cy,rotcx,rotcy)
    if supersampling>1:
        return ccode.average2d(proj, supersampling)
    return proj

def single_cone_projection(phantom, nx, ny, pixsize, angle, sod, sdd, zoff=0, supersampling=1, usecuda=True):
    if isinstance(phantom, str):
        with h5py.File(phantom, 'r') as f:
            spheres = f['spheres'][:]
    else:
        spheres = phantom
    if usecuda:
        proj = ccode.genconeproj_cuda(spheres, supersampling*nx, supersampling*ny, pixsize/supersampling, angle, sod, sdd, zoff=zoff)
    else:
        proj = ccode.genconeproj(spheres, supersampling*nx, supersampling*ny, pixsize/supersampling, angle, sod, sdd, zoff=zoff)
    if supersampling>1:
        return ccode.average2d(proj, supersampling)
    return proj