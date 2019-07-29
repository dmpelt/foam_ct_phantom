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

from . import project,generate,geometry
from .utils import FILE_VERSION

import numpy as np
import h5py
import tqdm

def generate_expand(outfile, phantom_file, seed, start_size=0.25, random_move=0.1, regularization=0.01, static_after_fraction=0.1, npoints=1024):
    with h5py.File(phantom_file, 'r') as f:
        spheres = f['spheres'][:]
        spatt = dict(f['spheres'].attrs)

    np.random.seed(seed)
    sizes = np.zeros(npoints, dtype=np.float32)
    v=1
    for i in range(1,npoints): 
        v = v + np.random.normal()*random_move + regularization*(1-v)
        if v<0:
            v = 0
        sizes[i] = sizes[i-1] + v
    sizes /= sizes[-1]
    sizes[:] = start_size + (1-start_size)*sizes
    sizes = np.interp(np.linspace(0,1,npoints), np.linspace(0,1-static_after_fraction,npoints), sizes)
    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        f['spheres'] = spheres
        att = f['spheres'].attrs
        for key, val in spatt.items():
            att[key] = val
        f['sizes'] = sizes
        att = f['sizes'].attrs
        att['seed'] = seed
        att['start_size'] = start_size
        att['random_move'] = random_move
        att['regularization'] = regularization
        att['static_after_fraction'] = static_after_fraction
        att['npoints'] = npoints

def __get_material_volume(spheres):
    minz = (spheres[2::5] - spheres[3::5]).min()
    maxz = (spheres[2::5] + spheres[3::5]).max()
    zdiff = maxz-minz
    cyl_vol = np.pi*zdiff
    sphere_vol = (4*np.pi*spheres[3::5]**3/3).sum()
    return cyl_vol-sphere_vol

def single_projection_expand(time, phantom, geom, angle):
    with h5py.File(phantom, 'r') as f:
        spheres = f['spheres'][:]
        sizes = f['sizes'][:]
    start_vol = __get_material_volume(spheres)
    size_time = np.interp(time, np.linspace(0,1,sizes.size), sizes)
    spheres[3::5] *= size_time
    cur_vol = __get_material_volume(spheres)
    factor = (start_vol/cur_vol)**(1/3)
    if type(geom) == geometry.ParallelGeometry:
        return project.single_par_projection(spheres.ravel(), geom.nx, geom.ny, geom.pixsize/factor, angle, geom.cx/factor, geom.cy/factor, geom.rotcx/factor, geom.rotcy/factor, geom.supersampling)*factor
    elif type(geom) == geometry.ConeGeometry:
        return project.single_cone_projection(spheres.ravel(), geom.nx, geom.ny, geom.pixsize/factor, angle, geom.sod/factor, (geom.sod + geom.odd)/factor, zoff=geom.zoff/factor, supersampling=geom.supersampling, usecuda=geom.usecuda)*factor

def genvol_expand(time, outfile, phantom, geom):
    with h5py.File(phantom, 'r') as f:
        spheres = f['spheres'][:]
        sizes = f['sizes'][:]
    start_vol = __get_material_volume(spheres)
    size_time = np.interp(time, np.linspace(0,1,sizes.size), sizes)
    spheres[3::5] *= size_time
    cur_vol = __get_material_volume(spheres)
    factor = (start_vol/cur_vol)**(1/3)
    generate.genvol(outfile, spheres, geom, zoomfactor=1/factor)

def gen_dataset_expand(outfile, phantom, geom):
    angles = geom.angles
    ny = geom.ny
    nx = geom.nx
    times = np.linspace(0, 1, len(angles))
    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        dset = f.create_dataset('projs', (len(angles),ny,nx), dtype='f4')
        for i in tqdm.trange(len(angles)):
            dset[i] = single_projection_expand(times[i],phantom,geom,angles[i])
        att = dset.attrs
        for key, val in geom.to_dict().items():
            att[key] = val
        att['phantom'] = phantom
        att['times'] = times