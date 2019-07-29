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

"""
Example 03: Generate parallel-beam projections
==============================================
"""

import foam_ct_phantom
import numpy as np
import h5py
import pylab as pl
pl.gray()

phantom = foam_ct_phantom.FoamPhantom('test_phantom.h5')

geom = foam_ct_phantom.ParallelGeometry(256,256,np.linspace(0,np.pi,128,False),3/256)

phantom.generate_projections('test_projs_par.h5',geom)

projs = foam_ct_phantom.load_projections('test_projs_par.h5')

pl.imshow(projs[0])
pl.show()