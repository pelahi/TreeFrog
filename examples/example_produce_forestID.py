# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 11:05:33 2017

Example of loading a walkable tree, associated halo properties and
building sublinks, progenitor links and forest ID.
The information is then saved in a forest file.
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
walkabletreefile=sys.argv[1]

#base halo properties
basepropfname=sys.argv[2]

#file name for the unified file
outputfname=sys.argv[3]

#define properties of interest
#this list can be updated as desired
#for the information you would like
#accessible in the resulting forest file.
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

# define different types of input
ASCIIINPUT=0
HDFINPUT=2

#here we can easily get the version of TF
TFVERSION = np.loadtxt('../VERSION')
#halo finder would need some updates
HFNAME = 'VELOCIraptor'
HFVERSION = 1.50

# could alter to have user indicate input type but currently assume all is HDF
RAWTREEFORMAT=HDFINPUT
RAWPROPFORMAT=HDFINPUT

# load the tree information stored in the file
# such as temporal halo id, number of snapshots searched when producing the tree
treedata,numsnaps=vpt.ReadWalkableHDFTree(inputhdftreefname)
numsnaps = treedata['Header']['TreeBuilder']['Number_of_snapshots']
TEMPORALHALOIDVAL = treedata['Header']['TreeBuilder']['Temporal_halo_id_value']
NSNAPSEARCH = treedata['Header']['TreeBuilder']['Nsteps_search_new_links']
TREEDIRECTION = treedata['Header']['TreeBuilder']['Tree_direction']

#alis tree data to halo data as forest file will combine the data
halodata = treedata['Snapshots']

numhalos = np.zeros(numsnaps, dtype=np.int64)
scalefactors = np.zeros(numsnaps)
halodata = [None for i in range(numsnaps)]

#load halo properties file
for i in range(numsnaps):
    fname=basepropfname+'%03d.VELOCIraptor'%i
    halos, numhalos[i] = vpt.ReadPropertyFile(fname,RAWPROPFORMAT,0,0,requestedfields)
    scalefactor[i]=halos['SimulationInfo']['ScaleFactor']
    if (ireducememfootprintflag and numhalos[i] > 0):
        for key in requestedfields:
            if (halos[key].dtype==np.float64):
                halos[key] = np.array(halos[key],dtype=np.float32)
    halodata[i].update(halos)

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
iforestcheck = True #this uses extra compute and is generally unnecessary
forestdata = vpt.GenerateForest(numsnaps, numhalos, halodata, scalefactor,
    maxnsnapsearch, ireverseorder, TEMPORALHALOIDVAL, iverbose, iforestcheck))

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
        'Title':'Forest',
        'TreeBuilder' : copy.deepcopy(treedata['Header']['TreeBuilder']),
        'HaloFinder' : copy.deepcopy(treedata['Header']['HaloFinder']),
        'Flag_gas':(igas==1), 'Flag_star':(istar==1), 'Flag_bh':(ibh==1),
        'Flag_subhalo_links':True, 'Flag_progenitor_links':True, 'Flag_forest_ids':True, 'Flag_sorted_forest':False
        }

vpt.WriteForest(outputfname, numsnaps, numhalos, halodata, forestdata, scalefactor,
                                   DescriptionInfo, SimulationInfo, UnitInfo)
