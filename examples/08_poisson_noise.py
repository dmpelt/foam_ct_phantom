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
Example 08: Apply Poission noise to projection data
===================================================
"""

import foam_ct_phantom
import numpy as np
import h5py
import pylab as pl
pl.gray()

# Run 03_generate_projections_parallel.py first

fac = foam_ct_phantom.estimate_absorption_factor('test_projs_par.h5',0.5)

foam_ct_phantom.apply_poisson_noise('test_projs_par.h5','test_projs_noisy.h5',1234,100,fac)

projs = foam_ct_phantom.load_projections('test_projs_noisy.h5')

pl.imshow(projs[0])
pl.show()