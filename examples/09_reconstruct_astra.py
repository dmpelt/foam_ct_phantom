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

"""
Example 09: Reconstruct data using ASTRA
========================================
"""

import foam_ct_phantom
import astra
import numpy as np
import h5py
import pylab as pl
pl.gray()

# Run 03_generate_projections_parallel.py and 04_generate_projections_cone first

projs = foam_ct_phantom.load_projections('test_projs_par.h5')

vol_geom = foam_ct_phantom.VolumeGeometry(256, 256, 256, 3/256)

proj_geom = foam_ct_phantom.ParallelGeometry.from_file('test_projs_par.h5')

pg = proj_geom.to_astra(single_slice=True)
vg = vol_geom.to_astra(single_slice=True)

pid = astra.create_projector('cuda', pg, vg)
w = astra.OpTomo(pid)

mid_slice = w.reconstruct('FBP_CUDA', projs[:,projs.shape[1]//2])

pl.imshow(mid_slice)
pl.show()

astra.projector.delete(pid)

projs = foam_ct_phantom.load_projections('test_projs_cone.h5')

proj_geom = foam_ct_phantom.ConeGeometry.from_file('test_projs_cone.h5', usecuda=False)

pg3d = proj_geom.to_astra()
vg3d = vol_geom.to_astra()

pid = astra.create_projector('cuda3d', pg3d, vg3d)
w = astra.OpTomo(pid)

cone_rec = w.reconstruct('FDK_CUDA', projs.swapaxes(0,1))

pl.imshow(cone_rec[cone_rec.shape[0]//2],vmin=0,vmax=3/256)
pl.show()





