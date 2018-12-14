# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:05:33 2017

Example of loading a tree, associated halo properties and applying
the branch fix code to deal with some tree pathologies. Such fixes
are typically unnecessary if a halo tracking code has been applied
but some final clean-up can be useful.
@author: Pascal Jahan Elahi

"""

import sys
import os
import glob
import psutil
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


#base raw tree file name to load the raw tree
basetreefname=sys.argv[1]

#number of snaphots
numsnaps=int(sys.argv[2])

#base halo properties file name
basepropfname=sys.argv[3]

#output name for cleaned walkable tree file
outputfname=sys.argv[4]

inputhdftreefname='null'
if (len(sys.argv)>5):
    #load the hdf tree file, means do not build
    #prelinary walkable tree
    inputhdftreefname=sys.argv[5]

#define properties of interest
requestedfields=[
    'ID', 'hostHaloID',
    'numSubStruct', 'npart',
    'R_200crit',
    'Xc', 'Yc', 'Zc',
    'VXc', 'VYc', 'VZc',
    'Rmax', 'Vmax',
    'Structuretype'
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

#load halo prop
numhalos=np.zeros(numsnaps,dtype=np.int64)
atime=np.zeros(numsnaps)
#load halo properties file
for i in range(numsnaps):
    fname=basepropfname+'_%03d.VELOCIraptor'%i
    halodata[i],numhalos[i] = vpt.ReadPropertyFile(fname,RAWPROPFORMAT,0,0,requestedfields)
    atime[i]=halodata[i]['SimulationInfo']['ScaleFactor']
    for key in halodata[i].keys():
        if (key == 'SimulationInfo' or key == 'UnitInfo'): continue
        if (halodata[i][key].dtype==np.float64):
            halodata[i][key] = np.array(halodata[i][key],dtype=np.float32)
#load raw tree info
snaptreelist=open(basetreefname+'.snaptreelist.txt','w')
for i in range(numsnaps):
    snaptreelist.write(basetreefname+'.snapshot_%03d.VELOCIraptor\n'%i)
snaptreelist.close()
rawtreedata=vpt.ReadHaloMergerTreeDescendant(basetreefname+'.snaptreelist.txt',False,HDFINPUT,1,True)

#load walkable tree information.
if (os.path.exists(inputhdftreefname)):
    treedata,numsnaps2=vpt.ReadWalkableHDFTree(inputhdftreefname)
    if (numsnaps != numsnaps2):
        print('ERROR, walkable tree number of snapshots does not match stored number of snapshots')
        print(numsnaps,numsnaps2)
        quit()
    #if loading walkable tree, assume tree processed for branch swapping events and need to update the halodata dictionary
    #link tree data to halo data
    for i in range(numsnaps):
        for key in treedata[i].keys():
            halodata[i][key] = treedata[i][key]
else:
    #build tree
    start=time.clock()
    #produce head tail in ascending order
    vpt.BuildTemporalHeadTailDescendant(numsnaps,rawtreedata,numhalos,halodata,TEMPORALHALOIDVAL)
    print("finished head tail ", time.clock()-start)

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
#and then clean secondary progenitors to make sure no secondary points to an object without
#a primary progenitor, useful for SAMs
npartlim=NPARTTHRESHOLD*2
vpt.CleanSecondaryProgenitorsFromNoPrimaryProgenObjectsTreeDescendant(numsnaps,
                                            halodata, halodata, numhalos,
                                            npartlim,
                                            TEMPORALHALOIDVAL,
                                            iverbose)


#write information
DescriptionInfo={
        'Title':'Tree and Halo', 'HaloFinder':'VELOCIraptor', 'TreeBuilder':'TreeFrog',
        'HaloFinder_version':1.25, 'TreeBuilder_version':1.2,
        'Particle_num_threshold':NPARTTHRESHOLD, 'Temporal_linking_length':NSNAPSEARCH, 'Temporal_halo_id_value':TEMPORALHALOIDVAL,
        }
vpt.WriteWalkableHDFTree(outputfname, numsnaps, rawtreedata, numhalos, halodata,
                         atime, DescriptionInfo)
