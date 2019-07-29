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
Example 02: Plot a phantom
==========================
"""

import foam_ct_phantom
import h5py
import pylab as pl
pl.gray()

phantom = foam_ct_phantom.FoamPhantom('test_phantom.h5')

geom = foam_ct_phantom.VolumeGeometry(256,256,1,3/256)

phantom.generate_volume('test_midslice.h5', geom)

vol = foam_ct_phantom.load_volume('test_midslice.h5')

pl.imshow(vol[0])
pl.show()