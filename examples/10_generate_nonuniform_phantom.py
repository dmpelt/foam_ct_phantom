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
Example 10: Generate a nonuniform phantom
=========================================
"""

import foam_ct_phantom
import numpy as np
import pylab as pl
pl.gray()

random_seed = 12345

def maxsize_func(x, y, z):
    return 0.2 - 0.1*np.abs(z)

# Note that nspheres_per_unit is set to a low value to reduce the computation time here.
# The default value is 100000.
foam_ct_phantom.FoamPhantom.generate('test_uniform.h5',random_seed,nspheres_per_unit=1000, maxsize=maxsize_func)

geom = foam_ct_phantom.VolumeGeometry(1,256,256,3/256)

phantom = foam_ct_phantom.FoamPhantom('test_uniform.h5')

phantom.generate_volume('test_uniform_midslice.h5', geom)

vol = foam_ct_phantom.load_volume('test_uniform_midslice.h5')

pl.imshow(vol[...,0])
pl.show()