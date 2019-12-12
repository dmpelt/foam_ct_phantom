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

import numpy as np
import random
import sortedcollections
import tqdm
import h5py
import inspect

from . import ccode, project, geometry
from .utils import FILE_VERSION



def genphantom(outfile, seed, nspheres_per_unit=100000, ntrials_per_unit=1000000, maxsize=0.2, zrange=1.5):
    n = int(ntrials_per_unit*zrange)
    nsph = int(nspheres_per_unit*zrange)
    
    pos3 = np.zeros(n*3,dtype=np.float32)                                                                                                              
    random.seed(seed)
    ccode.setseed(random.randint(0,4294967295))
    
    ds = np.zeros(n,dtype=np.float32)
    
    ccode.drawnewpositions(pos3, ds, zrange)

    if callable(maxsize):
        maxsizes = maxsize(pos3[::3],pos3[1::3], pos3[2::3])
        msk = ds<-maxsizes
        ds[msk] = -maxsizes[msk]

    upd = np.zeros(n, dtype=np.uint32)
    spheres = np.zeros(nsph*5, dtype=np.float32)

    sd = sortedcollections.ValueSortedDict(zip(range(ds.size), ds))

    for i in tqdm.trange(nsph):
        itms = sd.items()
        if callable(maxsize)==False and itms[0][1]<-maxsize:
            allchoices = []
            for itm in itms:
                if itm[1] >= -maxsize:
                    break
                allchoices.append(itm[0])
            ch = random.choice(allchoices)
            spheres[5*i+3] = maxsize
        else:
            allchoices = [itms[0][0],]
            curmax = itms[0][1]
            for itm in range(1,len(itms)):
                if itms[itm][1] == curmax:
                    allchoices.append(itms[itm][0])
                else:
                    break
            ch = random.choice(allchoices)
            spheres[5*i+3] = -sd[ch]
        spheres[5*i] = pos3[3*ch]
        spheres[5*i+1] = pos3[3*ch+1]
        spheres[5*i+2] = pos3[3*ch+2]
        nupd = ccode.newsphere(pos3, ds, spheres[:5*(i+1)], zrange, upd)
        if callable(maxsize):
            maxsizes = maxsize(pos3[3*upd[:nupd]],pos3[3*upd[:nupd]+1],pos3[3*upd[:nupd]+2])
            msk = ds[upd[:nupd]] < -maxsizes
            ds[upd[:nupd][msk]] = -maxsizes[msk]
        for ky in upd[:nupd]:
            sd[ky] = ds[ky]
    
    with h5py.File(outfile,'w') as f:
        f['spheres'] = spheres
        f.attrs['FILE_VERSION'] = FILE_VERSION
        att = f['spheres'].attrs
        att['seed'] = seed
        att['nspheres_per_unit'] = nspheres_per_unit
        att['ntrials_per_unit'] = ntrials_per_unit
        if callable(maxsize):
            att['maxsize'] = inspect.getsource(maxsize)
        else:
            att['maxsize'] = maxsize
        att['zrange'] = zrange


def genvol(outfile, phantom, geom, zoomfactor=1):
    nx = geom.nx
    ny = geom.ny
    nz = geom.nz
    voxsize = geom.voxsize*zoomfactor
    supersampling = geom.supersampling
    if isinstance(phantom, str):
        with h5py.File(phantom, 'r') as f:
            spheres = f['spheres'][:]
    else:
        spheres = phantom
    mi = np.argmin([nx,ny,nz])
    vol = np.zeros((nz, ny, nx), dtype=np.float32)
    for i in tqdm.trange(nz):
        ccode.genvol(spheres, vol, nx, ny, nz, voxsize, i, cx=geom.cx*zoomfactor, cy=geom.cy*zoomfactor, cz=geom.cz*zoomfactor, supersampling=supersampling)

    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        f['volume'] = vol
        att = f['volume'].attrs
        for key, val in geom.to_dict().items():
            att[key] = val
        if isinstance(phantom, str):
            att['phantom'] = phantom

def gen3d(phantom, nx, ny, pixsize, angle, tilt1, tilt2, maxz=1.5, cutout=0, cutoff=-np.inf):
    if isinstance(phantom, str):
        with h5py.File(phantom, 'r') as f:
            spheres = f['spheres'][:]
    else:
        spheres = phantom
    return ccode.gen3dproj(spheres, nx, ny, pixsize, angle, tilt1, tilt2, maxz=maxz, cutout=cutout, cutoff=cutoff)

def gen_dataset(outfile, phantom, geom):
    angles = geom.angles
    nx = geom.nx
    ny = geom.ny
    pixsize = geom.pixsize
    supersampling = geom.supersampling
    with h5py.File(outfile, 'w') as f:
        f.attrs['FILE_VERSION'] = FILE_VERSION
        dset = f.create_dataset('projs', (len(angles),ny,nx), dtype='f4')
        for i in tqdm.trange(len(angles)):
            if type(geom) == geometry.ParallelGeometry:
                dset[i] = project.single_par_projection(phantom,nx,ny,pixsize,angles[i],cx=geom.cx,cy=geom.cy,rotcx=geom.rotcx,rotcy=geom.rotcy, supersampling=supersampling)
            elif type(geom) == geometry.ConeGeometry:
                dset[i] = project.single_cone_projection(phantom, nx, ny, pixsize, angles[i], geom.sod, geom.sod + geom.odd, zoff=geom.zoff, supersampling=supersampling, usecuda=geom.usecuda)
        att = dset.attrs
        for key, val in geom.to_dict().items():
            if key == "projgeom_angles":
                f['projgeom_angles'] = val
            else:
                att[key] = val
        att['phantom'] = phantom

