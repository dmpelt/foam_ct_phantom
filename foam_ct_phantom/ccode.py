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

import ctypes
from pathlib2 import Path
import sys, os
import numpy as np

if sys.platform == 'darwin':
    libpath = Path(__file__).parent / 'libccode_c.dylib'
    lib = ctypes.CDLL(str(libpath))
elif os.name == 'nt':
    libpath = Path(__file__).parent.parent.parent.parent / 'bin' / 'ccode_c.dll'
    if not libpath.exists():
        libpath = Path(__file__).parent.parent / 'bin' / 'ccode_c.dll'
    lib = ctypes.WinDLL(str(libpath))
else:
    libpath = Path(__file__).parent / 'libccode_c.so'
    lib = ctypes.CDLL(str(libpath))


asuint = ctypes.c_uint32
asfloat = ctypes.c_float

cfloatp = ctypes.POINTER(ctypes.c_float)
def asfloatp(arr):
    return arr.ctypes.data_as(cfloatp)

cuintp = ctypes.POINTER(ctypes.c_uint32)
def asuintp(arr):
    return arr.ctypes.data_as(cuintp)

lib.newsphere.restype = ctypes.c_uint32
lib.gettouching.restype = ctypes.c_uint32
lib.iter_skiplist.restype = ctypes.c_int32

def drawnewpositions(pos3, ds, zrange):
    lib.drawnewpositions(asfloatp(pos3), asfloatp(ds), asuint(ds.size), asfloat(zrange))

def newsphere(pos3, ds, spheres, zrange, updated):
    return lib.newsphere(asfloatp(pos3), asfloatp(ds), asfloatp(spheres), asuint(ds.size), asuint(spheres.size//5), asfloat(zrange), asuintp(updated))

def gettouching(spheres, i, cutoff):
    touching = np.zeros(spheres.size//5, dtype=np.uint32)
    ntouch =  lib.gettouching(asfloatp(spheres), asuint(spheres.size//5), asuint(i), asfloat(cutoff), asuintp(touching))
    return touching[:ntouch].copy()

def setthreads(nthrds):
    lib.set_threads(asuint(nthrds))

def setseed(seed):
    lib.setseed(asuint(seed))

def genvol(spheres, vol, nx, ny, nz, voxsize, iz, cx=0, cy=0, cz=0, supersampling=1):
    n = np.array([nx,ny,nz],dtype=np.uint32)
    c = np.array([cx,cy,cz],dtype=np.float32)
    lib.genvol(asfloatp(spheres), asuint(spheres.size//5), asfloatp(vol.ravel()), asuintp(n), asfloat(voxsize), asfloatp(c), asuint(supersampling), asuint(iz))

def average2d(vol, supersampling):
    volnew = np.zeros((vol.shape[0]//supersampling, vol.shape[1]//supersampling), dtype=np.float32)
    lib.average2d(asfloatp(vol.ravel()), asfloatp(volnew.ravel()), asuint(volnew.shape[1]), asuint(volnew.shape[0]), asuint(supersampling))
    return volnew


def genparproj(spheres, nx, ny, pixsize, angle, cx=0, cy=0, rotcx=0, rotcy=0):
    proj = np.zeros((ny,nx),dtype=np.float32)
    n = np.array([nx,ny],dtype=np.uint32)
    c = np.array([cx,cy],dtype=np.float32)
    rotc = np.array([rotcx,rotcy],dtype=np.float32)
    lib.genparproj(asfloatp(spheres), asuint(spheres.size//5), asfloatp(proj.ravel()), asuintp(n), asfloat(pixsize), asfloatp(c), asfloat(angle), asfloatp(rotc))
    return proj

def gen3dproj(spheres, nx, ny, pixsize, angle, tilt1, tilt2, maxz=1.5, cutout=0, cutoff=-np.inf):
    sph_rot = spheres.copy()
    trot1 = np.array([[1,0,0],[0,np.cos(tilt1), -np.sin(tilt1)],[0,np.sin(tilt1),np.cos(tilt1)]])
    trot2 = np.array([[np.cos(tilt2),0,np.sin(tilt2)],[0,1,0],[-np.sin(tilt2), 0, np.cos(tilt2)]])
    arot = np.array([[np.cos(angle), -np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    rot = np.dot(trot2,np.dot(trot1,arot))
    pos = np.array([spheres[::5],spheres[1::5],spheres[2::5]])
    pos_rot = np.dot(rot,pos)
    trot1 = np.array([[1,0,0],[0,np.cos(-tilt1), -np.sin(-tilt1)],[0,np.sin(-tilt1),np.cos(-tilt1)]])
    trot2 = np.array([[np.cos(-tilt2),0,np.sin(-tilt2)],[0,1,0],[-np.sin(-tilt2), 0, np.cos(-tilt2)]])
    arot = np.array([[np.cos(-angle), -np.sin(-angle),0],[np.sin(-angle),np.cos(-angle),0],[0,0,1]])
    rot = np.dot(arot,np.dot(trot1,trot2)).astype(np.float32)
    sph_rot[::5] = pos_rot[0]
    sph_rot[1::5] = pos_rot[1]
    sph_rot[2::5] = pos_rot[2]
    proj = np.zeros((ny,nx),dtype=np.float32)
    n = np.array([nx,ny],dtype=np.uint32)
    lib.gen3dproj(asfloatp(sph_rot), asuint(spheres.size//5), asfloatp(proj.ravel()), asuintp(n), asfloat(pixsize), asfloat(maxz), asfloatp(rot), asuint(cutout), asfloat(cutoff))
    return proj

def genconeproj(spheres, nx, ny, pixsize, angle, sod, sdd, zoff=0):
    proj = np.zeros((ny,nx),dtype=np.float32)
    n = np.array([nx,ny],dtype=np.uint32)
    ca = np.cos(angle)
    sa = np.sin(angle)
    sph_rot = spheres.copy()
    rotx = ca*sph_rot[::5] - sa*sph_rot[1::5]
    sph_rot[1::5] = sa*sph_rot[::5] + ca*sph_rot[1::5]
    sph_rot[::5] = rotx
    lib.genconeproj(asfloatp(sph_rot), asuint(sph_rot.size//5), asfloatp(proj.ravel()), asuintp(n), asfloat(pixsize), asfloat(zoff), asfloat(sod), asfloat(sdd))
    return proj
    

def applypoisson(proj, flux, factor):
    lib.applypoisson(asfloatp(proj.ravel()), asuint(proj.size), asfloat(flux), asfloat(factor))

def init_skiplist(sizes):
    lib.init_skiplist(asfloatp(sizes), asuint(sizes.size))

def update_skiplist(sizes, idx):
    lib.update_skiplist(asfloatp(sizes), asuintp(idx), asuint(idx.size))

def skiplist_iter():
    lib.reset_iter_skiplist()
    val = lib.iter_skiplist()
    while val != -1:
        yield val
        val = lib.iter_skiplist()

# Try to set number of threads to number of physical cores
try:
    import psutil
    ncpu = psutil.cpu_count(logical=False)
    try:
        naff = len(psutil.Process().cpu_affinity())
        if naff < ncpu:
            ncpu = naff
    except AttributeError:
        pass
    setthreads(ncpu)
except ImportError:
    pass

try:
    from numba import cuda, float32, int32
    import math

    def genconeproj_cuda(spheres, nx, ny, pixsize, angle, sod, sdd, zoff=0):
        proj = cuda.to_device(np.zeros((ny,nx),dtype=np.float32))
        bpg0 = (ny + (32 - 1)) // 32
        bpg1 = (nx + (32 - 1)) // 32
        ca = np.cos(angle)
        sa = np.sin(angle)
        sph_rot = spheres.copy()
        rotx = ca*sph_rot[::5] - sa*sph_rot[1::5]
        sph_rot[1::5] = sa*sph_rot[::5] + ca*sph_rot[1::5]
        sph_rot[::5] = rotx
        sph_cuda = cuda.to_device(sph_rot)
        cudaconeproj_impl[(bpg0, bpg1), (32,32)](sph_cuda, proj, np.float32(pixsize), np.float32(sod), np.float32(sdd), np.float32(zoff))
        return proj.copy_to_host()

    @cuda.jit(fastmath=True)
    def cudaconeproj_impl(spheres, proj, pixsize, sod, sdd, zoff):
        iy, ix = cuda.grid(2)
        if iy<proj.shape[0] and ix<proj.shape[1]:
            x = (ix+0.5)*pixsize - proj.shape[1]*pixsize/2
            y = (iy+0.5)*pixsize - proj.shape[0]*pixsize/2
            
            # https://math.stackexchange.com/questions/2613781/line-cylinder-intersection
            x0 = -sod
            k = x/sdd
            l = y/sdd
            df = 1 - (x0*x0 - 1)*k*k
            if df>0:
                t1 = -(math.sqrt(df)+x0)/(k*k+1)
                t2 = (math.sqrt(df)-x0)/(k*k+1)
                dx = t2 - t1
                dy = k*dx
                dz = l*dx
                tmp = math.sqrt(dx*dx+dy*dy+dz*dz)
            else:
                proj[iy,ix]=0
                return

            sz = math.sqrt(1+k*k+l*l)
            tx = 1/sz
            ty = k/sz
            tz = l/sz

            # https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
            for j in range(spheres.shape[0]//5):
                sz = spheres[5*j+2]-zoff
                sx = spheres[5*j]
                sy = spheres[5*j+1]
                s2 = spheres[5*j+3]*spheres[5*j+3]
                psz = (x0-sx)*tx - sy*ty - sz*tz
                pdx = x0-sx - psz*tx
                pdy = -sy - psz*ty
                pdz = -sz - psz*tz
                dist = pdx*pdx+pdy*pdy+pdz*pdz
                if dist<s2:
                    tmp -= 2*(1-spheres[5*j+4])*math.sqrt(s2 - dist)
            proj[iy,ix] = tmp
except:
    pass