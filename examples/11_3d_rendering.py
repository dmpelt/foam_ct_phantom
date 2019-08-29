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
Example 11: Generate 3D renderings of voids
===========================================
"""

import foam_ct_phantom
import numpy as np
import pylab as pl
pl.gray()

phantom = foam_ct_phantom.FoamPhantom('test_phantom.h5')
rendering = phantom.generate_3d(1024, 1024, 4/1024, 0, -np.pi/8, -np.pi/8)
pl.imshow(rendering)

expanding_phantom = foam_ct_phantom.ExpandingFoamPhantom('test_expand.h5')
rendering = expanding_phantom.generate_3d(1024, 1024, 4/1024, 0, -np.pi/8, -np.pi/8, time=0.5)
pl.figure()
pl.imshow(rendering)

moving_phantom = foam_ct_phantom.MovingFoamPhantom('test_moving.h5')
rendering = moving_phantom.generate_3d(1024, 1024, 4/1024, 0, -np.pi/8, -np.pi/8, time=0.5, maxz=0.75)
pl.figure()
pl.imshow(rendering)

infiltration_phantom = foam_ct_phantom.InfiltrationFoamPhantom('test_infiltration.h5')
rendering = infiltration_phantom.generate_3d(1024, 1024, 4/1024, 0, -np.pi/8, -np.pi/8, time=0.5)
pl.figure()
pl.imshow(rendering)

pl.show()