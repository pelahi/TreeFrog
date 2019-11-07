# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:05:33 2017

Example of loading a raw tree, building the head/tail, root head/ root tail info writing a walkable tree.
@author: Pascal Jahan Elahi

"""

import sys
import os
import glob
import time
import numpy as np
import copy
import h5py

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

#base raw tree file name to load the raw tree if necessary
basetreefname=sys.argv[1]

#director that contains the halo catalogs
halocatalogdir=sys.argv[2]

#file name for the simplified tree file
outputfname=sys.argv[3]

#requested halo fields
requestedfields=[
    'ID', 'hostHaloID',
    ]

# define different types of input
ASCIIINPUT=0
HDFINPUT=2

#here we can easily get the version of TF
TFVERSION = np.loadtxt(basecodedir+'../VERSION')
#halo finder would need some updates
HFNAME = 'VELOCIraptor'
HFVERSION = 1.50

# could alter to have user indicate input type but currently assume all is HDF
RAWTREEFORMAT=HDFINPUT
RAWPROPFORMAT=HDFINPUT

# load the tree information stored in the file
# such as temporal halo id, number of snapshots searched when producing the tree
if (RAWTREEFORMAT == HDFINPUT):
    hfile = h5py.File(basetreefname+'.snapshot_000.VELOCIraptor.tree')
    numsnaps = hfile.attrs['Number_of_snapshots']
    TEMPORALHALOIDVAL = hfile.attrs['Temporal_halo_id_value']
    NSNAPSEARCH = hfile.attrs['Nsteps_search_new_links']
    TREEDIRECTION = hfile.attrs['Tree_direction']
    hfile.close()
#if ascii, output of tree frog can be parsed
else :
    #currently the ascii file does not store the temporal halo id but this will be updated
    treefile = open(basetreename, 'r')
    numsnaps = np.int32(treefile.readline())
    description = treefile.readline()
    treedirectionstring = description.split('Produce tree in direction  ').split(' |')[0]
    if (treedirectionstring == 'progenitors'):
        TREEDIRECTION = 0
    elif (treedirectionstring == 'descendants'):
        TREEDIRECTION = 1
    NSNAPSEARCH = np.int32(description.split('Tree built using ').split(' temporal steps')[0])
    TEMPORALHALOIDVAL = 1000000000000

# number of particles used in the halo catalog
# this can also be extracted from the halo catalog files directly, specifically
# the configuration files
NPARTTHRESHOLD=20

rawtreedata = None
#read raw descendant tree along with merits, don't reverse snap order,
if (RAWTREEFORMAT == HDFINPUT):
    #for hdf input produce file listing input files
    snaptreelist=open(basetreefname+'.snaptreelist.txt','w')
    for i in range(numsnaps):
        snaptreelist.write(basetreefname+'.snapshot_%03d.VELOCIraptor\n'%i)
    snaptreelist.close()
    fname = basetreefname+'.snaptreelist.txt'
else:
    fname = basetreename

if (TREEDIRECTION == 1):
    rawtreedata=vpt.ReadHaloMergerTreeDescendant(fname, False, RAWTREEFORMAT, 1, True)
elif (TREEDIRECTION == 0):
    print('Warning, progenitor based trees are less useful when building halo merger trees for SAMs.')
    rawtreedata=vpt.ReadHaloMergerTree(fname, RAWTREEFORMAT, 1, True)
else:
    print('Full graphs to walkable trees not implemented yet.')

print('Finished reading raw tree')
numhalos=np.zeros(numsnaps,dtype=np.uint64)
halodata=[dict() for i in range(numsnaps)]
atime=np.zeros(numsnaps)
for i in range(numsnaps):
    halodata[i],numhalos[i]=vpt.ReadPropertyFile(halocatalogdir+'/snapshot_%03d.VELOCIraptor'%i, 2, 0, 1, requestedfields)
    atime[i]=halodata[i]['SimulationInfo']['ScaleFactor']
    for key in halodata[i].keys():
        if (key == 'SimulationInfo' or key == 'UnitInfo' or key == "ConfigurationInfo"): continue
        if (halodata[i][key].dtype==np.float64):
            halodata[i][key] = np.array(halodata[i][key],dtype=np.float32)
print('Finished reading halo properties')

#produce head tail in ascending order
start=time.clock()
print("Building head/tail ")
vpt.BuildTemporalHeadTailDescendant(numsnaps, rawtreedata, numhalos, halodata,
    TEMPORALHALOIDVAL)
print("Finished head/tail ", time.clock()-start)

#store the description
DescriptionInfo={
        'Title':'Walkable Tree',
        'TreeBuilder' : {
            'Name' : 'TreeFrog',
            'Version' : TFVERSION,
            'Temporal_linking_length' : NSNAPSEARCH,
            'Temporal_halo_id_value' : TEMPORALHALOIDVAL,
            'Tree_direction' : TREEDIRECTION,
        },
        'HaloFinder' : {
            'Name' : HFNAME, 'Version' : HFVERSION,
            'Particle_num_threshold' : NPARTTHRESHOLD,
            },
        }

vpt.WriteWalkableHDFTree(outputfname, numsnaps, rawtreedata, numhalos, halodata,
                         atime, DescriptionInfo)
