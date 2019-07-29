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

from . import project,generate,geometry
from .utils import FILE_VERSION

import numpy as np
import h5py
import tqdm

def generate_verticalmovement(outfile, phantom_file, seed, zmin, zmax, random_move=0.1, regularization=0.01, npoints=1024):
    with h5py.File(phantom_file, 'r') as f:
        spheres = f['spheres'][:]
        spatt = dict(f['spheres'].attrs)

    np.random.seed(seed)
    ypos = np.zeros(npoints, dtype=np.float32)
    v=1
    for i in range(1,npoints): 
        v = v + np.random.normal()*random_move + regularization*(1-v)
        if v<0:
            v = 0
        ypos[i] = ypos[i-1] + v
    ypos /= ypos[-1]
    ypos[:] = zmin + ypos*(zmax-zmin)
    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        f['spheres'] = spheres
        att = f['spheres'].attrs
        for key, val in spatt.items():
            att[key] = val
        f['ypos'] = ypos
        att = f['ypos'].attrs
        att['seed'] = seed
        att['zmin'] = zmin
        att['zmax'] = zmax
        att['random_move'] = random_move
        att['regularization'] = regularization
        att['npoints'] = npoints

def single_projection_verticalmovement(time, phantom, geom, angle):
    with h5py.File(phantom, 'r') as f:
        spheres = f['spheres'][:]
        ypos = f['ypos'][:]
    ypos_time = np.interp(time, np.linspace(0,1,ypos.size), ypos)
    spheres[2::5] -= ypos_time
    if type(geom) == geometry.ParallelGeometry:
        return project.single_par_projection(spheres.ravel(), geom.nx, geom.ny, geom.pixsize, angle, geom.cx, geom.cy, geom.rotcx, geom.rotcy, geom.supersampling)
    elif type(geom) == geometry.ConeGeometry:
        return project.single_cone_projection(spheres.ravel(), geom.nx, geom.ny, geom.pixsize, angle, geom.sod, geom.sod + geom.odd, zoff=geom.zoff, supersampling=geom.supersampling, usecuda=geom.usecuda)

def genvol_verticalmovement(time, outfile, phantom, geom):
    with h5py.File(phantom, 'r') as f:
        spheres = f['spheres'][:]
        ypos = f['ypos'][:]
    ypos_time = np.interp(time, np.linspace(0,1,ypos.size), ypos)
    spheres[2::5] -= ypos_time
    generate.genvol(outfile, spheres, geom)

def gen_dataset_verticalmovement(outfile, phantom, geom):
    angles = geom.angles
    ny = geom.ny
    nx = geom.nx
    times = np.linspace(0, 1, len(angles))
    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        dset = f.create_dataset('projs', (len(angles),ny,nx), dtype='f4')
        for i in tqdm.trange(len(angles)):
            dset[i] = single_projection_verticalmovement(times[i],phantom,geom,angles[i])
        att = dset.attrs
        for key, val in geom.to_dict().items():
            att[key] = val
        att['phantom'] = phantom
        att['times'] = times