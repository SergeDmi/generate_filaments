#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright Serge Dmitrieff
# www.biophysics.fr

import numpy as np
import sys
import yaml
import pandas as pd
import time
from scipy.stats import special_ortho_group

__VERSION__ = "0.0.0"

"""
# SYNOPSIS
    ---

# DESCRIPTION
    ---

# SYNTAX
    $ ./filament_makery.py [CONFIG_FILE] [OPTIONS]

    CONFIG_FILE is a suitable yaml-type config file, e.g.

---
# OPTIONS
    See config file ; no other options available for now.

    The position of the filament is the position of the top left corner in the yz plane


# OUTPUT :



# @TODO : *


"""


def version():
    return __VERSION__

class Maker:
    """
        Class containing the filament maker

        Reads config file and creates the filaments
    """
    def __init__(self,*args,fname='config.yaml',get_radius_c=None,get_length=None,out="exported_filaments.txt"):
        self.config=read_config(fname)
        self.points=[]
        self.get_radius_c=get_radius_c
        self.get_length=get_length
        self.filaments=[]
        self.out=out

    def initialize(self,*args,**kwargs):
        """ Initializes the maker """
        config=self.config
        np.random.seed(int(time.time()))
        if self.get_radius_c is None:
            self.get_radius_c=lambda r: config['curvature_radius']['center']+r*(config['curvature_radius']['edge']-config['curvature_radius']['center'])
        if self.get_length is None:
            self.get_length   =lambda r: config['length']['center']+r*(config['length']['edge']-config['length']['center'])
        if "export" in config.keys():
            self.out=config['export']
        if "random_rotation" in config['setup'].keys():
            self.random_rotation=config['setup']['random_rotation']
        else:
            self.random_rotation=1
        self.Zmax=config['setup']['Zmax']
        self.Rmax=config['setup']['Rmax']
        self.step=config['setup']['step']
        self.n_filaments=0

    def make_all_filaments(self):
        for i in range(self.config["setup"]["number"]):
            self.add_filament()


    def add_filament(self,*args,position=None,**kwargs):
        """ Add a filament to self.filaments
        """
        # r,t,z are the cylindrical coordinates fo the filament center
        if position is None:
            position=self.get_filament_center()
        self.filaments.append(self.make_filament(position))
        self.n_filaments+=1


    def get_filament_center(self):
        r=np.sqrt(np.random.random())
        t=2.0*np.pi*np.random.random()
        z=np.random.random()
        return [r,t,z]

    def save_filaments(self):
        fid=open(self.out,"w")
        total=0
        for i,filament in enumerate(self.filaments):
            for j,segment in enumerate(filament):
                str_coords=",".join([str(pt) for pt in segment])
                fid.write("%s,%s,%s,%s\n" %(total,str_coords,i,j))
                total+=1
        fid.close()


    def make_filament(self,position):
        # cylindrical coordinates
        [r,t,z]=position[:3]
        # Bam !
        radius_c=self.get_radius_c(r)      # radius of curvature
        length=self.get_length(r)        # length of filament
        n=int(length/self.step)
        filament=np.zeros((n,3))
        thetas=np.arange(n)*(self.step/radius_c)

        # Creation and centering
        filament[:,0]=radius_c*(np.cos(thetas)-np.mean(np.cos(thetas)))
        filament[:,1]=radius_c*(np.sin(thetas)-np.mean(np.sin(thetas)))
        # rotation, using a scipy package to have a uniform 3D rotation
        if self.random_rotation:
            filament[:,:]=np.dot(filament,special_ortho_group.rvs(3))
        # Translation
        filament[:,0]+=self.Rmax*r*np.cos(t)
        filament[:,1]+=self.Rmax*r*np.sin(t)
        filament[:,2]+=z*self.Zmax
        return filament


def read_config(fname):
    # Reads a config file to a dictionary
    #   a wrapper to yaml.load
    file=open(fname,'r')
    config=yaml.load(file, Loader=yaml.SafeLoader)
    file.close()
    return config







if __name__ == "__main__":
    # The program as called from the command line
    #   Here we just gather arguments to create and run the Integrator

    # the config file, if any, should be the first argument
    nargs=len(sys.argv)
    args=[]
    if nargs<2:
        # if no config file name was given
        fname="config.cym"
    else:
        # a config file name was given !
        fname=sys.argv[1]
        if nargs>2:
            # if there are more arguments (there shouldn't for now)
            args=sys.argv[2:]

    # Now we create the integrator
    maker=Maker(*args,fname=fname)
    # We initialize and run the simulation
    maker.initialize()
    print("# Starting to make")
    maker.make_all_filaments()
    # and save
    print("# Starting to save")
    maker.save_filaments()
