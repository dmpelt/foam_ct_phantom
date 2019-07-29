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

from skbuild import setup
from setuptools import find_packages
setup(
    name='foam_ct_phantom',
    packages=find_packages(),
    version=open('VERSION').read().strip(),
    include_package_data=True,
    cmake_languages=['C',],
    classifiers=[
        "License :: OSI Approved :: MIT License",
    ],
)
