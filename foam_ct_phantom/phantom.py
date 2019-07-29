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

from . import generate, infiltrate, project, verticalmovement, expand

import abc
import os.path

class Phantom(object):
    
    __metaclass__ = abc.ABCMeta

    def __init__(self, filename):
        self.filename = filename
        if not os.path.exists(self.filename):
            raise ValueError("{} does not exist.".format(self.filename))
    
    @abc.abstractmethod
    def generate_projections(self, outfile, geom):
        pass
    
    @abc.abstractmethod
    def generate_volume(self, outfile, geom, time=0):
        pass


class FoamPhantom(Phantom):

    @staticmethod
    def generate(filename, seed, nspheres_per_unit=100000, ntrials_per_unit=1000000, maxsize=0.2, zrange=1.5):
        generate.genphantom(filename, seed, nspheres_per_unit=nspheres_per_unit, ntrials_per_unit=ntrials_per_unit, maxsize=maxsize, zrange=zrange)
    
    def generate_projections(self, outfile, geom):
        generate.gen_dataset(outfile,self.filename, geom)
    
    def generate_volume(self, outfile, geom):
        generate.genvol(outfile, self.filename,geom)

class MovingFoamPhantom(Phantom):

    @staticmethod
    def generate(filename, phantom_file, seed, zmin, zmax, random_move=0.1, regularization=0.01, npoints=1024):
        verticalmovement.generate_verticalmovement(filename, phantom_file, seed, zmin, zmax, random_move=random_move, regularization=regularization, npoints=npoints)
    
    def generate_projections(self, outfile, geom):
        verticalmovement.gen_dataset_verticalmovement(outfile,self.filename, geom)
    
    def generate_volume(self, outfile, geom, time=0):
        verticalmovement.genvol_verticalmovement(time, outfile, self.filename, geom)

class ExpandingFoamPhantom(Phantom):

    @staticmethod
    def generate(filename, phantom_file, seed, start_size=0.25, random_move=0.1, regularization=0.01, static_after_fraction=0.1, npoints=1024):
        expand.generate_expand(filename, phantom_file, seed, start_size=start_size, random_move=random_move, regularization=regularization, static_after_fraction=0.1, npoints=npoints)
    
    def generate_projections(self, outfile, geom):
        expand.gen_dataset_expand(outfile,self.filename, geom)
    
    def generate_volume(self, outfile, geom, time=0):
        expand.genvol_expand(time, outfile, self.filename, geom)

class InfiltrationFoamPhantom(Phantom):

    @staticmethod
    def generate(filename, phantom_file, seed, fluid_value, startz=-1.5, rand_width=1, cutoff=1e-5):
        infiltrate.generate_inflitration(filename, phantom_file, seed, fluid_value, startz=startz, rand_width=rand_width, cutoff=cutoff)

    def generate_projections(self, outfile, geom):
        infiltrate.gen_dataset_infiltrate(outfile,self.filename,geom)
    
    def generate_volume(self, outfile, geom, time=0):
        infiltrate.genvol_infiltrate(time, outfile, self.filename,geom)