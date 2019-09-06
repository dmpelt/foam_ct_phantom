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
import numpy as np
import random
import tqdm
from . import ccode,project,generate,geometry
from .utils import FILE_VERSION

def generate_inflitration(outfile, phantom_file, seed, fluid_value, startz=1.5, rand_width=1, cutoff=1e-5):
    random.seed(seed)
    with h5py.File(phantom_file, 'r') as f:
        spheres = f['spheres'][:]
        spatt = dict(f['spheres'].attrs)
    
    nsph = spheres.shape[0]//5

    
    all_touching = []
    for i in range(nsph):
        all_touching.append(ccode.gettouching(spheres,i,cutoff))

    arrival = np.zeros(nsph, dtype=np.float32)
    arrival[:] = np.inf

    to_be_processed = set()
    for i in range(nsph):
        if np.abs(spheres[5*i+2]-startz) < spheres[5*i+3]:
            arrival[i] = 0
            to_be_processed.add(i)
    while(len(to_be_processed)>0):
        proc = to_be_processed.pop()
        for i in all_touching[proc]:
            newtime = arrival[proc] + 1 + (random.random()-0.5)*rand_width
            if newtime<arrival[i]:
                arrival[i] = newtime
                to_be_processed.add(i)
    
    arrival /= arrival[np.isinf(arrival)==False].max()
    
    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        f['spheres'] = spheres
        att = f['spheres'].attrs
        for key, val in spatt.items():
            att[key] = val
        f['arrival'] = arrival
        att = f['arrival'].attrs
        att['seed'] = seed
        att['startz'] = startz
        att['rand_width'] = rand_width
        att['cutoff'] = cutoff
        att['fluid_value'] = fluid_value

def single_projection_infiltrate(time, phantom, geom, angle):
    with h5py.File(phantom, 'r') as f:
        spheres = f['spheres'][:].reshape((-1,5))
        arrival = f['arrival'][:]
        fluid_value = f['arrival'].attrs['fluid_value']
    
    spheres[arrival<=time,-1] = fluid_value
    spheres[arrival>time,-1] = 0

    if type(geom) == geometry.ParallelGeometry:
        return project.single_par_projection(spheres.ravel(), geom.nx, geom.ny, geom.pixsize, angle, geom.cx, geom.cy, geom.rotcx, geom.rotcy, geom.supersampling)
    elif type(geom) == geometry.ConeGeometry:
        return project.single_cone_projection(spheres.ravel(), geom.nx, geom.ny, geom.pixsize, angle, geom.sod, geom.sod + geom.odd, zoff=geom.zoff, supersampling=geom.supersampling, usecuda=geom.usecuda)

def genvol_infiltrate(time, outfile, phantom, geom):
    with h5py.File(phantom, 'r') as f:
        spheres = f['spheres'][:].reshape((-1,5))
        arrival = f['arrival'][:]
        fluid_value = f['arrival'].attrs['fluid_value']
    
    spheres[arrival<=time,-1] = fluid_value
    spheres[arrival>time,-1] = 0
    return generate.genvol(outfile, spheres.ravel(), geom)

def gen3d_infiltrate(time, phantom, nx, ny, pixsize, angle, tilt1, tilt2, maxz=1.5, cutout=0, cutoff=-np.inf):
    with h5py.File(phantom, 'r') as f:
        spheres = f['spheres'][:].reshape((-1,5))
        arrival = f['arrival'][:]
        fluid_value = f['arrival'].attrs['fluid_value']
    
    spheres[arrival<=time,-1] = fluid_value
    spheres[arrival>time,-1] = 0
    return ccode.gen3dproj(spheres.ravel(), nx, ny, pixsize, angle, tilt1, tilt2, maxz=maxz, cutout=cutout, cutoff=fluid_value)

def gen_dataset_infiltrate(outfile, phantom, geom):
    angles = geom.angles
    nx = geom.nx
    ny = geom.ny
    times = np.linspace(0, 1, len(angles))
    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        dset = f.create_dataset('projs', (len(angles),ny,nx), dtype='f4')
        for i in tqdm.trange(len(angles)):
            dset[i] = single_projection_infiltrate(times[i],phantom,geom,angles[i])
        att = dset.attrs
        for key, val in geom.to_dict().items():
            att[key] = val
        att['phantom'] = phantom
        att['times'] = times


