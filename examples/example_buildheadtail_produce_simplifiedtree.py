# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:05:33 2017

Example of loading a raw tree, building the head/tail, root head/ root tail info and fixing the code for branch swapping events and writing a simplified HDF tree.
@author: Pascal Jahan Elahi


"""

import sys
import os
import glob
import time
import numpy as np
import copy

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

#base halo properties
basepropfname=sys.argv[2]

#number of snapshots
numsnaps=int(sys.argv[3])

#file name for the simplified tree file
outputfname=sys.argv[4]

#if post-process for branch swaps, more time consuming and requires more memory.
ibranchfix=bool(int(sys.argv[5]))
#define properties of interest
if (ibranchfix):
    requestedfields=[
        'ID', 'hostHaloID',
        'npart',
        'Xc', 'Yc', 'Zc',
        'VXc', 'VYc', 'VZc',
        'Rmax', 'Vmax', 'R_200crit',
        'Structuretype'
        ]
else:
    requestedfields=[
        'ID', 'hostHaloID',
        ]


#define different types of input
ASCIIINPUT=0
HDFINPUT=2

#could alter to have user indicate input type but currently assume all is HDF
RAWTREEFORMAT=HDFINPUT
RAWPROPFORMAT=HDFINPUT
TEMPORALHALOIDVAL=1000000000000
#number of snapshots searched when producing the tree
NSNAPSEARCH=4
#number of particles used in the halo catalog
NPARTTHRESHOLD=20

#read raw descendant tree along with merits, don't reverse snap order,
if (RAWTREEFORMAT == HDFINPUT):
    #for hdf input produce file listing input files
    snaptreelist=open(basetreefname+'.snaptreelist.txt','w')
    for i in range(numsnaps):
        snaptreelist.write(basetreefname+'.snapshot_%03d.VELOCIraptor\n'%i)
    snaptreelist.close()
    rawtreedata=vpt.ReadHaloMergerTreeDescendant(basetreefname+'.snaptreelist.txt', False, RAWTREEFORMAT, 1, True)
else:
    rawtreedata=vpt.ReadHaloMergerTreeDescendant(basetreefname, False, RAWTREEFORMAT, 0, True)
print('Finished reading raw tree')
numhalos=np.zeros(numsnaps,dtype=np.uint64)
halodata=[dict() for i in range(numsnaps)]
atime=np.zeros(numsnaps)
for i in range(numsnaps):
    halodata[i],numhalos[i]=vpt.ReadPropertyFile(basepropfname+'%03d.VELOCIraptor'%i, 2, 0, 1, requestedfields)
    atime[i]=halodata[i]['SimulationInfo']['ScaleFactor']
    for key in halodata[i].keys():
        if (key == 'SimulationInfo' or key == 'UnitInfo' or key == "ConfigurationInfo"): continue
        if (halodata[i][key].dtype==np.float64):
            halodata[i][key] = np.array(halodata[i][key],dtype=np.float32)
print('Finished reading halo properties')

#produce head tail in ascending order
start=time.clock()
print("Building head/tail ")
vpt.BuildTemporalHeadTailDescendant(numsnaps,rawtreedata,numhalos,halodata,TEMPORALHALOIDVAL)
print("Finished head/tail ", time.clock()-start)

if (ibranchfix):
    #set the limits to use when trying to fix branch swapping events.
    npartlim=NPARTTHRESHOLD*10
    #suggested values that reduce the number of outliers (well resolved objects with no progenitor)
    meritlim=0.01
    xdifflim=5.0
    vdifflim=2.5
    numsnapsearch = np.int32(np.ceil(NSNAPSEARCH*1.5))
    descendantsearchdepth = 2
    iswaphalosubhaloflag = True
    iverbose=0
    vpt.FixTruncationBranchSwapsInTreeDescendant(numsnaps, rawtreedata, halodata, numhalos,
                                                 npartlim, meritlim, xdifflim, vdifflim,
                                                 numsnapsearch,
                                                 descendantsearchdepth, iswaphalosubhaloflag,
                                                 TEMPORALHALOIDVAL,
                                                 iverbose)

#write the tree and quit
DescriptionInfo={
        'Title':'Walkable Tree',
        'TreeBuilder':'TreeFrog', 'TreeBuilder_version':1.20, 'Temporal_linking_length':NSNAPSEARCH, 'Temporal_halo_id_value':TEMPORALHALOIDVAL,
        'HaloFinder':'VELOCIraptor', 'HaloFinder_version':1.11, 'Particle_num_threshold':20,
        }

vpt.WriteWalkableHDFTree(outputfname, numsnaps, rawtreedata, numhalos, halodata,
                         atime, DescriptionInfo)
quit()
