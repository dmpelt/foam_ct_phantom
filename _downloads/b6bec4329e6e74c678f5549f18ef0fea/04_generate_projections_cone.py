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
Example 04: Generate cone-beam projections
==========================================
"""

import foam_ct_phantom
import numpy as np
import h5py
import pylab as pl
pl.gray()

phantom = foam_ct_phantom.FoamPhantom('test_phantom.h5')

geom = foam_ct_phantom.ConeGeometry(256,256,np.linspace(0,2*np.pi,256),3/256,sod=2,odd=0,usecuda=False)

phantom.generate_projections('test_projs_cone.h5',geom)

projs = foam_ct_phantom.load_projections('test_projs_cone.h5')

pl.imshow(projs[0])
pl.show()
