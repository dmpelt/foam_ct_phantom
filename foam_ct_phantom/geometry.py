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

class VolumeGeometry(object):
    def __init__(self, nx, ny, nz, voxsize, cx=0, cy=0, cz=0, supersampling=1):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.voxsize = voxsize
        self.cx = cx
        self.cy = cy
        self.cz = cz
        self.supersampling = supersampling
    
    def to_astra(self, single_slice=False):
        if self.cx != 0 or self.cy != 0 or self.cz != 0:
            raise ValueError('Shifted center not supported yet')
        import astra
        if single_slice:
            return astra.create_vol_geom((self.ny, self.nx))
        else:
            return astra.create_vol_geom((self.nz, self.ny, self.nx))
    
    def to_dict(self):
        dct = {}
        for att in ['nx', 'ny', 'nz', 'voxsize', 'cx', 'cy', 'cz', 'supersampling']:
            dct['volgeom_'+att] = getattr(self,att)
        return dct
    
    @classmethod
    def from_file(cls, fn):
        args = []
        kwargs = {}
        with h5py.File(fn, 'r') as f:
            att = f['volume'].attrs
            for key in ['nx', 'ny', 'nz', 'voxsize']:
                args.append(att['volgeom_'+key])
            for key in ['cx', 'cy', 'cz', 'supersampling']:
                kwargs[key] = att['volgeom_'+key]
        return cls(*args, **kwargs)

class ParallelGeometry(object):

    def __init__(self, nx, ny, angles, pixsize, cx=0, cy=0, rotcx=0, rotcy=0, supersampling=1):
        self.nx = nx
        self.ny = ny
        self.angles = angles
        self.pixsize = pixsize
        self.cx = cx
        self.cy = cy
        self.rotcx = rotcx
        self.rotcy = rotcy
        self.supersampling = supersampling
    
    def to_astra(self, single_slice=False):
        import astra
        if single_slice:
            return astra.create_proj_geom('parallel', 1, self.ny, self.angles)
        else:
            return astra.create_proj_geom('parallel3d', 1, 1, self.ny, self.nx, self.angles)
    
    def to_dict(self):
        dct = {}
        dct['projgeom_type'] = 'parallel'
        for att in ['nx', 'ny', 'angles', 'pixsize', 'cx', 'cy', 'rotcx', 'rotcy', 'supersampling']:
            dct['projgeom_'+att] = getattr(self,att)
        return dct
    
    @classmethod
    def from_file(cls, fn):
        args = []
        kwargs = {}
        with h5py.File(fn, 'r') as f:
            att = f['projs'].attrs
            for key in ['nx', 'ny', 'angles', 'pixsize']:
                if key == 'angles' and 'projgeom_angles' in f.keys():
                    args.append(f['projgeom_angles'][:])
                else:
                    args.append(att['projgeom_'+key])
            for key in ['cx', 'cy', 'rotcx', 'rotcy', 'supersampling']:
                kwargs[key] = att['projgeom_'+key]
        return cls(*args, **kwargs)

class ConeGeometry(object):

    def __init__(self, nx, ny, angles, pixsize, sod=None, sdd=None, odd=None, supersampling=1, zoff=0, usecuda=True):
        self.nx = nx
        self.ny = ny
        self.angles = angles
        self.pixsize = pixsize
        if sod is not None and sdd is not None and odd is not None:
            raise ValueError("Two of SOD, SDD, or ODD have to be given")
        elif sod is not None and sdd is not None:
            self.sod = sod
            self.odd = sdd - sod
        elif sod is not None and odd is not None:
            self.sod = sod
            self.odd = odd
        elif sdd is not None and odd is not None:
            self.sod = sdd - odd
            self.odd = odd
        else:
            raise ValueError("Two of SOD, SDD, or ODD have to be given")
        self.supersampling = supersampling
        self.zoff = zoff
        self.usecuda = usecuda
    
    def to_astra(self):
        import astra
        return astra.create_proj_geom('cone', 1, 1, self.ny, self.nx, self.angles, self.sod/self.pixsize, self.odd/self.pixsize)
    
    def to_dict(self):
        dct = {}
        dct['projgeom_type'] = 'cone'
        for att in ['nx', 'ny', 'angles', 'pixsize', 'sod', 'odd', 'supersampling', 'zoff', 'usecuda']:
            dct['projgeom_'+att] = getattr(self,att)
        return dct
    
    @classmethod
    def from_file(cls, fn, usecuda=True):
        args = []
        kwargs = {}
        with h5py.File(fn, 'r') as f:
            att = f['projs'].attrs
            for key in ['nx', 'ny', 'angles', 'pixsize']:
                if key == 'angles' and 'projgeom_angles' in f.keys():
                    args.append(f['projgeom_angles'][:])
                else:
                    args.append(att['projgeom_'+key])
            for key in ['sod', 'odd', 'supersampling', 'zoff']:
                kwargs[key] = att['projgeom_'+key]
            kwargs['usecuda'] = usecuda
        return cls(*args, **kwargs)