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
import tqdm
import h5py
import inspect
import queue
import threading

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

    ccode.init_skiplist(ds)

    for i in tqdm.trange(nsph):
        itr = ccode.skiplist_iter()
        smallidx = next(itr)
        if callable(maxsize)==False and ds[smallidx]<-maxsize:
            allchoices = [smallidx,]
            for itm in itr:
                if ds[itm] >= -maxsize:
                    break
                allchoices.append(itm)
            ch = random.choice(allchoices)
            spheres[5*i+3] = maxsize
        else:
            allchoices = [smallidx,]
            curmax = ds[smallidx]
            for itm in itr:
                if ds[itm] == curmax:
                    allchoices.append(itm)
                else:
                    break
            ch = random.choice(allchoices)
            spheres[5*i+3] = -ds[ch]
        spheres[5*i] = pos3[3*ch]
        spheres[5*i+1] = pos3[3*ch+1]
        spheres[5*i+2] = pos3[3*ch+2]
        nupd = ccode.newsphere(pos3, ds, spheres[:5*(i+1)], zrange, upd)
        if callable(maxsize):
            maxsizes = maxsize(pos3[3*upd[:nupd]],pos3[3*upd[:nupd]+1],pos3[3*upd[:nupd]+2])
            msk = ds[upd[:nupd]] < -maxsizes
            ds[upd[:nupd][msk]] = -maxsizes[msk]
        ccode.update_skiplist(ds,upd[:nupd])
    
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
        if hasattr(geom, 'usecuda') and hasattr(geom.usecuda, '__iter__'):
            # Use multiple GPUs
            ngpus = len(geom.usecuda)
            tmpdata = np.zeros((ngpus, ny, nx), dtype=np.float32)
            def worker(threadid):
                ccode.set_cuda_device(geom.usecuda[threadid])
                while True:
                    item = q.get()
                    if item is None:
                        break
                    tmpdata[threadid] = project.single_cone_projection(phantom, nx, ny, pixsize, angles[item], geom.sod, geom.sod + geom.odd, zoff=geom.zoff, supersampling=supersampling, usecuda=True)
                    q.task_done()
            q = queue.Queue()
            threads = [threading.Thread(target=worker, args=(i,)) for i in range(ngpus)]
            for t in threads:
                t.start()
            pbar = tqdm.tqdm(total=len(angles))
            for i in range(len(angles)/ngpus):
                for j in range(i*ngpus, (i+1)*ngpus):
                    q.put(j)
                q.join()
                for j in range(ngpus):
                    dset[i*ngpus+j] = tmpdata[j]
                pbar.update(4)
            for i in range(int(len(angles)/ngpus)*ngpus, len(angles)):
                dset[i] = project.single_cone_projection(phantom, nx, ny, pixsize, angles[i], geom.sod, geom.sod + geom.odd, zoff=geom.zoff, supersampling=supersampling, usecuda=True)
                pbar.update(1)
            pbar.close()
            for i in range(ngpus):
                q.put(None)
            for t in threads:
                t.join()
        else:
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

