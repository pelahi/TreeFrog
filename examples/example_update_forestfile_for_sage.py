# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:05:33 2017

Example of loading a forest file and updating for
use with sage. This sorts forests (prunes as well for all
forests with <2 halos, and adds meta information)
@author: Pascal Jahan Elahi

"""

import sys
import os
import glob

#load python routines
scriptpath=os.path.abspath(__file__)
basecodedir=scriptpath.split('examples/')[0]+'/tools/'
sys.path.append(basecodedir)
#load the cythonized code if compiled
if (len(glob.glob(basecodedir+'velociraptor_python_tools_cython.*.so'))==1):
    print('using cython VR+TF toolkit')
    import velociraptor_python_tools_cython as vpt
else:
    print('using python VR+TF toolkit')
    import velociraptor_python_tools as vpt

#base file name of the forest file
fname=sys.argv[1]

vpt.PruneForest(fname)
vpt.ForestSorter(fname)
vpt.ForceBiDirectionalTreeInForestFile(fname)
vpt.ForestFileAddMetaData(fname)
