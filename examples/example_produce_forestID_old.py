# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:05:33 2017

Example of loading a tree, associated halo properties and building sublinks, progenitor links and forest IDs
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

#load the hdf tree file
inputhdftreefname=sys.argv[1]

#base raw tree file name to load the raw tree if necessary
basetreefname=sys.argv[2]

#base halo properties
basepropfname=sys.argv[3]

#file name for the unified file
outputfname=sys.argv[4]

#define properties of interest
requestedfields=[
    'ID', 'hostHaloID',
    'numSubStruct', 'npart',
    'Mass_tot', 'Mass_FOF', 'Mass_200mean', 'Mass_200crit',
    'R_size', 'R_HalfMass', 'R_200mean', 'R_200crit',
    'Xc', 'Yc', 'Zc',
    'VXc', 'VYc', 'VZc',
    'lambda_B',
    'Lx','Ly','Lz',
    'RVmax_Lx','RVmax_Ly','RVmax_Lz',
    'sigV', 'RVmax_sigV',
    'Rmax', 'Vmax',
    'cNFW',
    'Efrac','Structuretype'
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


#load walkable tree information.
if (os.path.exists(inputhdftreefname)):
    rawtreedata=[]
    halodata,numsnaps=vpt.ReadWalkableHDFTree(inputhdftreefname)
    #if loading walkable tree, assume tree processed for branch swapping events and need to update the halodata dictionary
    numhalos=np.zeros(numsnaps,dtype=np.int64)
    atime=np.zeros(numsnaps)
    print(numsnaps)
    #load halo properties file
    for i in range(numsnaps):
        fname=basepropfname+'%03d.VELOCIraptor'%i
        print(fname)
        halos,numhalos[i] = vpt.ReadPropertyFile(fname,RAWPROPFORMAT,0,0,requestedfields)
        atime[i]=halos['SimulationInfo']['ScaleFactor']
        halodata[i].update(halos)
        for key in halodata[i].keys():
            if (key == 'SimulationInfo' or key == 'UnitInfo'): continue
            if (halodata[i][key].dtype==np.float64):
                halodata[i][key] = np.array(halodata[i][key],dtype=np.float32)
else:
    #load raw tree information if necessary. if tree is HDF input then write a txt file that lists all the snapshots produced by treefrog
    snaptreelist=open(basetreefname+'.snaptreelist.txt','w')
    for i in range(numsnaps):
        snaptreelist.write(basetreefname+'.snapshot_%03d.VELOCIraptor\n'%i)
    snaptreelist.close()
    rawtreedata=vpt.ReadHaloMergerTreeDescendant(basetreefname+'.snaptreelist.txt',False,HDFINPUT,1,True)
    #if raw tree, then next load the halodata properties
    numhalos=np.zeros(numsnaps,dtype=np.int64)
    atime=np.zeros(numsnaps)
    #load halo properties file
    for i in range(numsnaps):
        fname=basepropfname+'_%03d.VELOCIraptor'%i
        halodata[i],numhalos[i] = vpt.ReadPropertyFile(fname,RAWPROPFORMAT,0,0,requestedfields)
        atime[i]=halodata[i]['SimulationInfo']['ScaleFactor']
        for key in halodata[i].keys():
            if (halodata[i][key].dtype==np.float64):
                halodata[i][key] = np.array(halodata[i][key],dtype=np.float32)
    #build tree
    start=time.clock()
    #produce head tail in ascending order
    vpt.BuildTemporalHeadTailDescendant(numsnaps,rawtreedata,numhalos,halodata,TEMPORALHALOIDVAL)
    print("finished head tail ", time.clock()-start)

#given walkable tree, determine the largest difference in snapshots between an object and its head
maxnsnapsearch=0
for i in range(numsnaps):
    if (numhalos[i] == 0): continue
    headsnap = np.int64(halodata[i]['Head']/TEMPORALHALOIDVAL)
    maxs = np.max(headsnap-i)
    maxnsnapsearch = max(maxnsnapsearch, maxs)
print('Walkable tree has maximum snaps search of ', maxnsnapsearch)
sys.stdout.flush()

#generate subhalo links
vpt.GenerateSubhaloLinks(numsnaps,numhalos,halodata)
#generate progenitor links
vpt.GenerateProgenitorLinks(numsnaps,numhalos,halodata)
#building forest
ireverseorder = False
iverbose = 1
iforestcheck = True
ForestStats=vpt.GenerateForest(numsnaps,numhalos,halodata,atime,maxnsnapsearch, ireverseorder, TEMPORALHALOIDVAL, iverbose, iforestcheck))

#strip out simulation and unit data
SimulationInfo=copy.deepcopy(halodata[0]['SimulationInfo'])
UnitInfo=copy.deepcopy(halodata[0]['UnitInfo'])
if SimulationInfo['Cosmological_Sim']:
    if not UnitInfo['Comoving_or_Physical']:
        #convert period to comoving little h
        SimulationInfo['Period']*=SimulationInfo['h_val']/SimulationInfo['ScaleFactor']
    del SimulationInfo['ScaleFactor']
    #lets update the names for genesis style
    UnitInfo['Comoving_unit_flag'] = UnitInfo['Comoving_or_Physical']

#remove the unit info from each snapshot as it is unnecessary
for i in range(numsnaps):
    del halodata[i]['SimulationInfo']
    del halodata[i]['UnitInfo']

#currently only dark matter runs
igas=istar=ibh=0
#write the unified file that contains forest ids, properties,
#description will have to be updated so as to use appropriate version numbers
DescriptionInfo={
        'Title':'Tree and Halo', 'HaloFinder':'VELOCIraptor', 'TreeBuilder':'TreeFrog',
        'HaloFinder_version':1.25, 'TreeBuilder_version':1.2,
        'Particle_num_threshold':20, 'Temporal_linking_length':4, 'Temporal_halo_id_value':TEMPORALHALOIDVAL,
        'Flag_gas':(igas==1), 'Flag_star':(istar==1), 'Flag_bh':(ibh==1),
        'Flag_subhalo_links':True, 'Flag_progenitor_links':True, 'Flag_forest_ids':True, 'Flag_sorted_forest':False
        }
vpt.WriteUnifiedTreeandHaloCatalog(outputfname, numsnaps, rawtreedata, numhalos, halodata, atime,
                                   DescriptionInfo, SimulationInfo, UnitInfo)
