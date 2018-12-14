
# coding: utf-8

# # Plots to test tree metrics

#load other useful packages
import numpy as np
import astropy as ap
import h5py, copy
#for general python stuff
import sys,os,string,time,re,struct,glob
import math,operator
#for useful scipy stuff
from scipy.stats.mstats import mquantiles
from scipy.misc import comb
import scipy.interpolate as scipyinterp
import scipy.integrate as scipyint
import scipy.optimize as scipyopt
import scipy.special as scipysp
from scipy.stats import gaussian_kde
from sklearn.neighbors.kde import KernelDensity
import itertools
#for useful mathematical tools
from sklearn.neighbors import NearestNeighbors
import scipy.spatial as spatial
import emcee, corner
#matplotlib items
import matplotlib
matplotlib.use('pdf') #this backend works fine on python (I think) but not ipython
from matplotlib.pylab import *
#to load specific functions defined in another python file

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

NumArg=5
if (len(sys.argv) < NumArg+1):
    print ('Number of Args provided ',len(sys.argv)-1, 'Number that must be provided', NumArg)
    print ('This script produces some tree metris')
    print ('Input : ')
    print ('1) Label of analysis')
    print ('2) Input walkable tree file name')
    print ('3) Input halo properties base file name')
    print ('4) Base plot directory ')
    print ('Extra Options :')
    print ('5) Extra Diagnostics flag (0, none; >=1 Merger rates; >=2 Merger radius statistics;')
    quit()

def GetHaloLiftime(treedata,haloindex,halosnap,haloid,atime,TEMPORALHALOIDVAL=1000000000000):
    """
    Get the length of a halo's progenitors
    """
    halolifetime=1
    curid=haloid
    curindex=int(curid%TEMPORALHALOIDVAL-1)
    cursnap=int(curid/TEMPORALHALOIDVAL)
    curhead=treedata[cursnap]["Head"][curindex]
    while (curid!=curhead):
        halolifetime+=1
        curid=curhead
        curindex=int(curid%TEMPORALHALOIDVAL-1)
        cursnap=int(curid/TEMPORALHALOIDVAL)
        curhead=treedata[cursnap]["Head"][curindex]
    return halolifetime



# ## Loading Halos


#define formats
ASCIIFORMAT=0
HDFFORMAT=2

#base filename
labelname=sys.argv[1]
walkabletreename=sys.argv[2]
basesnapname=sys.argv[3]
baseplotdir=sys.argv[4]
ExtraDiagnostics = False
if (len(sys.argv) == NumArg+1):
    ExtraDiagnostics = bool(int(sys.argv[5]))

#set default values
TEMPORALHALOIDVAL=1000000000000
NSNAPSEARCH=4

print('Reading tree ',walkabletreename, 'snapshots ',
      basesnapname,' and will plot and save data to ', baseplotdir,
      'under label', labelname)
treedata,numsnaps = vpt.ReadWalkableHDFTree(walkabletreename, True)
hfile=h5py.File(walkabletreename)
#if key present in header, use these
TEMPORALHALOIDVAL=hfile['Header/TreeBuilder'].attrs['Temporal_halo_id_value']
NSNAPSEARCH=hfile['Header/TreeBuilder'].attrs['Temporal_linking_length']
hfile.close()

numhalosintree = np.zeros(numsnaps,dtype=np.uint64)
for isnap in range(numsnaps):
    numhalosintree[isnap] = treedata[isnap]['ID'].size

#now we try loading halo property data
poskeys=['Xc','Yc','Zc','R_200crit','R_200mean','Rmax','R_size']
desiredkeys=['npart','Efrac', 'ID', 'hostHaloID','Structuretype','Mass_tot', 'Xc', 'Yc', 'Zc', 'R_200crit','Mass_200crit']
intkeys=['npart','ID','hostHaloID','Structuretype','SimulationInfo','UnitInfo']
halopropdata = [dict() for isnap in range(numsnaps)]
atime = np.zeros(numsnaps)
numhalos = np.zeros(numsnaps,dtype=np.uint64)
for isnap in range(numsnaps):
    inputfname=basesnapname+'%03d.VELOCIraptor'%isnap
    print(inputfname)
    halopropdata[isnap],numhalos[isnap] = vpt.ReadPropertyFile(inputfname,HDFFORMAT,0,0,desiredkeys)
    for key in halopropdata[isnap].keys():
        if key not in intkeys:
            halopropdata[isnap][key]=np.array(halopropdata[isnap][key],dtype=np.float32)
    atime[isnap] = halopropdata[isnap]['SimulationInfo']['ScaleFactor']
    numhalosintree[isnap] = treedata[isnap]['ID'].size
    if (numhalos[isnap] != numhalosintree[isnap]):
        print('ERORR: Number of halo catalogs disagrees with number of halos in tree! At snapshot',
            isnap, 'containing', numhalos[isnap],
            'with tree containing', numhalosintree[isnap])
        print('CHECK tree, halo catalogs or input files names')
        exit(9)

#store simulation info
SimulationInfo=dict()
for key in halopropdata[-1]['SimulationInfo'].keys():
    SimulationInfo[key] = halopropdata[-1]['SimulationInfo'][key]
#should generalize to whether comoving or not
SimulationInfo['Period'] = halopropdata[-1]['SimulationInfo']['Period'] / halopropdata[-1]['SimulationInfo']['ScaleFactor']

UnitInfo=dict()
for key in halopropdata[-1]['UnitInfo'].keys():
    UnitInfo[key]=halopropdata[-1]['UnitInfo'][key]

period = SimulationInfo['Period']
converttocomove = ['Xc', 'Yc', 'Zc', 'Rmax', 'R_200crit']
keys = halopropdata[0].keys()
for key in converttocomove:
    if key not in keys:
        converttocomove.remove(key)
# convert positions and sizes to comoving if necesary
if (UnitInfo['Comoving_or_Physical'] == 0 and SimulationInfo['Cosmological_Sim'] == 1):
    print('Converting to comoving')
    for i in range(numsnaps):
        if (numhalos[i] == 0):
            continue
        for key in converttocomove:
            halopropdata[i][key] /= atime[i]

print('Finished reading data, now analyse')

#set plotting parameters
rcParams['figure.figsize'] = (10,10)
rcParams['font.size'] = 24
rc('text', usetex=True)
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
rcParams['xtick.direction']='in'
rcParams['ytick.direction']='in'
rcParams['axes.linewidth']=2
rcParams['ytick.major.size']=7.5
rcParams['ytick.minor.size']=7.5
rcParams['ytick.major.width']=1.5
rcParams['ytick.minor.width']=1.5
rcParams['xtick.major.size']=7.5
rcParams['xtick.minor.size']=7.5
rcParams['xtick.major.width']=1.5
rcParams['xtick.minor.width']=1.5

#
fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Scale Factor'
y1title=r'(sub)Halo Number density'
xvals=atime[np.where(numhalos>0)]
yvals=numhalos[np.where(numhalos>0)]/(SimulationInfo['Period']*UnitInfo['Length_unit_to_kpc']*1000.)**3.0
ax.plot(xvals,yvals,color='black',ls='solid',lw=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,1)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-halonumberdesnity.pdf')
fig1.clf()

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
y1title=r'(sub)Halo Number density'
xvals=np.arange(numsnaps)[np.where(numhalos>0)]
yvals=numhalos[np.where(numhalos>0)]/(SimulationInfo['Period']*UnitInfo['Length_unit_to_kpc']*1000.)**3.0
ax.plot(xvals,yvals,color='black',ls='solid',lw=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-halonumberdesnity.pdf')
fig1.clf()
plt.close(fig1)

#
fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Scale Factor'
y1title=r'Mass in Halos density'
xvals=atime[np.where(numhalos>0)]
icount=0
for i in np.where(numhalos>0)[0]:
    yvals[icount]=np.sum(halopropdata[i]['Mass_tot'])*UnitInfo['Mass_unit_to_solarmass']/(SimulationInfo['Period']*UnitInfo['Length_unit_to_kpc']*1000.)**3.0
    icount+=1
    print
ax.plot(xvals,yvals,color='black',ls='solid',lw=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,1)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-massinhalosdesnity.pdf')
fig1.clf()
plt.close(fig1)

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
y1title=r'Mass in Halos density'
xvals=np.arange(numsnaps)[np.where(numhalos>0)]
icount=0
for i in np.where(numhalos>0)[0]:
    yvals[icount]=np.sum(halopropdata[i]['Mass_tot'])*UnitInfo['Mass_unit_to_solarmass']/(SimulationInfo['Period']*UnitInfo['Length_unit_to_kpc']*1000.)**3.0
    icount+=1
    print
ax.plot(xvals,yvals,color='black',ls='solid',lw=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-massinhalosdesnity.pdf')
fig1.clf()
plt.close(fig1)

#get no descendant information
numdata={'Low':np.zeros(numsnaps),'Intermediate':np.zeros(numsnaps),'High':np.zeros(numsnaps), 'All':np.zeros(numsnaps)}
numnodescen={'Low':np.zeros(numsnaps),'Intermediate':np.zeros(numsnaps),'High':np.zeros(numsnaps), 'All':np.zeros(numsnaps)}
numnoimmediatedescen={'Low':np.zeros(numsnaps),'Intermediate':np.zeros(numsnaps),'High':np.zeros(numsnaps), 'All':np.zeros(numsnaps)}
npartlim=20
nodescennumdistrib=np.zeros([numsnaps,7])
noimmediatedescennumdistrib=np.zeros([numsnaps,7])
for i in range(numsnaps):
    if (numhalos[i] == 0): continue
    numdata['All'][i] = numhalos[i]
    numdata['Low'][i] = np.where((halopropdata[i]['npart'] < 2*npartlim))[0].size
    numdata['Intermediate'][i] = np.where((halopropdata[i]['npart'] >= 2*npartlim)*(halopropdata[i]['npart'] < 4*npartlim))[0].size
    numdata['High'][i] = np.where((halopropdata[i]['npart'] >= 4*npartlim))[0].size

    numnodescen['All'][i] = np.where(treedata[i]['Head'] == treedata[i]['ID'])[0].size
    numnodescen['Low'][i] = np.where((treedata[i]['Head'] == treedata[i]['ID'])*(halopropdata[i]['npart'] < 2*npartlim))[0].size
    numnodescen['Intermediate'][i] = np.where((treedata[i]['Head'] == treedata[i]['ID'])*(halopropdata[i]['npart'] >= 2*npartlim)*(halopropdata[i]['npart'] < 4*npartlim))[0].size
    numnodescen['High'][i] = np.where((treedata[i]['Head'] == treedata[i]['ID'])*(halopropdata[i]['npart'] >= 4*npartlim))[0].size
    if (numnodescen['All'][i] == 0): continue
    data = halopropdata[i]['npart'][np.where(treedata[i]['Head'] == treedata[i]['ID'])]
    nodescennumdistrib[i][0] = np.min(data)
    nodescennumdistrib[i][6] = np.max(data)
    nodescennumdistrib[i][1:6] = np.percentile(data,[2.5, 16.0, 50.0, 84.0, 97.5])

    numnoimmediatedescen['All'][i] = np.where(treedata[i]['HeadSnap'] > i+1)[0].size
    numnoimmediatedescen['Low'][i] = np.where((treedata[i]['HeadSnap'] > i+1)*(halopropdata[i]['npart'] < 2*npartlim))[0].size
    numnoimmediatedescen['Intermediate'][i] = np.where((treedata[i]['HeadSnap'] > i+1)*(halopropdata[i]['npart'] >= 2*npartlim)*(halopropdata[i]['npart'] < 4*npartlim))[0].size
    numnoimmediatedescen['High'][i] = np.where((treedata[i]['HeadSnap'] > i+1)*(halopropdata[i]['npart'] >= 4*npartlim))[0].size
    if (numnoimmediatedescen['All'][i] == 0): continue
    data = halopropdata[i]['npart'][np.where(treedata[i]['HeadSnap'] > i+1)]
    noimmediatedescennumdistrib[i][0] = np.min(data)
    noimmediatedescennumdistrib[i][6] = np.max(data)
    noimmediatedescennumdistrib[i][1:6] = np.percentile(data,[2.5, 16.0, 50.0, 84.0, 97.5])

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Scale Factor'
y1title=r'$f$'
colorval={'All': 'black', 'Low':'blue','Intermediate':'DarkOrange','High':'Crimson'}
labelset=['']
xvals=atime[:-1]
for key in numdata.keys():
    yvals=(numnodescen[key]/numdata[key])[:-1]
    ax.plot(xvals,yvals,color=colorval[key],ls='solid',lw=2,label=key)
ax.legend(loc='lower right',numpoints=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,1)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-nodescenstats.pdf')
fig1.clf()
plt.close(fig1)

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshot'
y1title=r'$f$'
colorval={'All': 'black', 'Low':'blue','Intermediate':'DarkOrange','High':'Crimson'}
labelset=['']
xvals=np.arange(numsnaps)[:-1]
for key in numdata.keys():
    yvals=(numnodescen[key]/numdata[key])[:-1]
    ax.plot(xvals,yvals,color=colorval[key],ls='solid',lw=2,label=key)
ax.legend(loc='lower right',numpoints=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-nodescenstats-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Scale Factor'
y1title=r'$N_p(a_{\rm disrupt})$'
labelset=['']
xvals=atime[:-1]
yvals=nodescennumdistrib[:-1].transpose()
ymean=yvals[3]
ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4)
yel,yeh=yvals[2],yvals[4]
ax.fill_between(xvals,yel,yeh,
    facecolor='None',edgecolor='Blue',alpha=0.75,interpolate=True,zorder=1,linewidth=4,linestyle='dashed')
yell,yehh=yvals[1],yvals[5]
ax.fill_between(xvals,yell,yehh,
    facecolor='None',edgecolor='Teal',alpha=0.75,interpolate=True,zorder=1,linewidth=2,linestyle='dotted')
yelll,yehhh=yvals[0],yvals[6]
ax.plot(xvals,yehhh,color='Cyan',alpha=0.75,zorder=1,linewidth=1,linestyle='solid')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,1)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(10,110)
ax.set_yscale("linear")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-nodescen-numdistrib.pdf')
fig1.clf()

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
y1title=r'$N_p(a_{\rm disrupt})$'
labelset=['']
xvals=np.arange(numsnaps)[:-1]
yvals=nodescennumdistrib[:-1].transpose()
ymean=yvals[3]
ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4)
yel,yeh=yvals[2],yvals[4]
ax.fill_between(xvals,yel,yeh,
    facecolor='None',edgecolor='Blue',alpha=0.75,interpolate=True,zorder=1,linewidth=4,linestyle='dashed')
yell,yehh=yvals[1],yvals[5]
ax.fill_between(xvals,yell,yehh,
    facecolor='None',edgecolor='Teal',alpha=0.75,interpolate=True,zorder=1,linewidth=2,linestyle='dotted')
yelll,yehhh=yvals[0],yvals[6]
ax.plot(xvals,yehhh,color='Cyan',alpha=0.75,zorder=1,linewidth=1,linestyle='solid')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(10,110)
ax.set_yscale("linear")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-nodescen-numdistrib-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

#get no progenitor information
numdata={'Low':np.zeros(numsnaps),'Intermediate':np.zeros(numsnaps),'High':np.zeros(numsnaps), 'All':np.zeros(numsnaps)}
numnoprogen={'Low':np.zeros(numsnaps),'Intermediate':np.zeros(numsnaps),'High':np.zeros(numsnaps), 'All':np.zeros(numsnaps)}
numnoprogentype={'Halo':np.zeros(numsnaps),'Subhalo':np.zeros(numsnaps),'Merger Remnants':np.zeros(numsnaps)}
npartlim=20
noprogennumdistrib=np.zeros([numsnaps,7])
for i in range(numsnaps):
    if (numhalos[i] == 0): continue
    numdata['All'][i] = numhalos[i]
    numdata['Low'][i] = np.where((halopropdata[i]['npart'] < 2*npartlim))[0].size
    numdata['Intermediate'][i] = np.where((halopropdata[i]['npart'] >= 2*npartlim)*(halopropdata[i]['npart'] < 4*npartlim))[0].size
    numdata['High'][i] = np.where((halopropdata[i]['npart'] >= 4*npartlim))[0].size

    numnoprogen['All'][i] = np.where(treedata[i]['Tail'] == treedata[i]['ID'])[0].size
    numnoprogen['Low'][i] = np.where((treedata[i]['Tail'] == treedata[i]['ID'])*(halopropdata[i]['npart'] < 2*npartlim))[0].size
    numnoprogen['Intermediate'][i] = np.where((treedata[i]['Tail'] == treedata[i]['ID'])*(halopropdata[i]['npart'] >= 2*npartlim)*(halopropdata[i]['npart'] < 4*npartlim))[0].size
    numnoprogen['High'][i] = np.where((treedata[i]['Tail'] == treedata[i]['ID'])*(halopropdata[i]['npart'] >= 4*npartlim))[0].size
    if (numnoprogen['All'][i] == 0): continue
    data = halopropdata[i]['npart'][np.where(treedata[i]['Tail'] == treedata[i]['ID'])]
    noprogennumdistrib[i][0] = np.min(data)
    noprogennumdistrib[i][6] = np.max(data)
    noprogennumdistrib[i][1:6] = np.percentile(data,[2.5, 16.0, 50.0, 84.0, 97.5])

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Scale Factor'
y1title=r'$f$'
colorval={'All': 'black', 'Low':'blue','Intermediate':'DarkOrange','High':'Crimson'}
labelset=['']
xvals=atime
for key in numdata.keys():
    yvals=numnoprogen[key]/(numdata[key]+1e-32)
    ax.plot(xvals,yvals,color=colorval[key],ls='solid',lw=2,label=key)
ax.legend(loc='lower right',numpoints=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,1)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-noprogenitor.pdf')
fig1.clf()
plt.close(fig1)

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
y1title=r'$f$'
colorval={'All': 'black', 'Low':'blue','Intermediate':'DarkOrange','High':'Crimson'}
labelset=['']
xvals=np.arange(numsnaps)
for key in numdata.keys():
    yvals=numnoprogen[key]/(numdata[key]+1e-32)
    ax.plot(xvals,yvals,color=colorval[key],ls='solid',lw=2,label=key)
ax.legend(loc='lower right',numpoints=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_yscale("log")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-noprogenitor-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Scale Factor'
y1title=r'$N_p(a_{\rm form})$'
labelset=['']
xvals=atime
yvals=noprogennumdistrib.transpose()
ymean=yvals[3]
ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4)
yel,yeh=yvals[2],yvals[4]
ax.fill_between(xvals,yel,yeh,
    facecolor='None',edgecolor='Blue',alpha=0.75,interpolate=True,zorder=1,linewidth=4,linestyle='dashed')
yell,yehh=yvals[1],yvals[5]
ax.fill_between(xvals,yell,yehh,
    facecolor='None',edgecolor='Teal',alpha=0.75,interpolate=True,zorder=1,linewidth=2,linestyle='dotted')
yelll,yehhh=yvals[0],yvals[6]
ax.plot(xvals,yehhh,color='Cyan',alpha=0.75,zorder=1,linewidth=1,linestyle='solid')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,1)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(10,110)
ax.set_yscale("linear")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-noprogenitor-numdistrib.pdf')
fig1.clf()
plt.close(fig1)

fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
y1title=r'$N_p(a_{\rm form})$'
labelset=['']
xvals=np.arange(numsnaps)
yvals=noprogennumdistrib.transpose()
ymean=yvals[3]
ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4)
yel,yeh=yvals[2],yvals[4]
ax.fill_between(xvals,yel,yeh,
    facecolor='None',edgecolor='Blue',alpha=0.75,interpolate=True,zorder=1,linewidth=4,linestyle='dashed')
yell,yehh=yvals[1],yvals[5]
ax.fill_between(xvals,yell,yehh,
    facecolor='None',edgecolor='Teal',alpha=0.75,interpolate=True,zorder=1,linewidth=2,linestyle='dotted')
yelll,yehhh=yvals[0],yvals[6]
ax.plot(xvals,yehhh,color='Cyan',alpha=0.75,zorder=1,linewidth=1,linestyle='solid')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(10,110)
ax.set_yscale("linear")
ax.set_ylabel(y1title)
plt.savefig(baseplotdir+labelname+'-noprogenitor-numdistrib-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

#get unusual events (like objects with no progenitor that are large)
noprogendata={
            'Npart' : np.array([]),
            'Npart_subtohost' : np.array([]),
            'Npart_hosttosubs' : np.array([]),
            'Npart_stats' : np.zeros([numsnaps,7]),
            'Npart_subtohost_stats' : np.zeros([numsnaps,7]),
            'Npart_hosttosubs_stats' : np.zeros([numsnaps,8]),
            'Npart_structype_stats' : {'Halos': np.zeros([numsnaps,7]), 'Subhalos': np.zeros([numsnaps,7]), 'Merger Remnants': np.zeros([numsnaps,7])},
            'Zform' : np.array([]),
            'Zfinal' : np.array([]),
            'Npart_descen' : np.array([]),
            'Branch_type' : np.array([]),
            'Branch_lifetime' : np.array([]),
            'Sub_type' : np.array([]),
            'Sub_type_descen' : np.array([]),
            'Descen_type' : np.array([]),
            'Environ_type' : np.array([]),
            }
noprogenfrac={
    'All': np.zeros(numsnaps),
    'Masscut' : np.zeros(numsnaps),
    'Halos': np.zeros(numsnaps),
    'Subs': np.zeros(numsnaps),
    'MainBranch': np.zeros(numsnaps),
    'ShortBranch': np.zeros(numsnaps),
    'DescenType' : {
        'Same': np.zeros(numsnaps),
        'Change': np.zeros(numsnaps),
        'HaloToSub': np.zeros(numsnaps),
        'SubToHalo': np.zeros(numsnaps),
        }
    }
for i in range(numsnaps):
    if (numhalos[i] == 0):
        continue
    sizewdata = np.where((halopropdata[i]['npart'] > 10*npartlim))[0]
    if (sizewdata.size == 0):
        continue
    noprogenwdata = np.where(treedata[i]['Tail'][sizewdata] == treedata[i]['ID'][sizewdata])[0]
    if (noprogenwdata.size == 0):
        continue
    numval=noprogenwdata.size
    noprogenfrac['All'][i]=float(numval)/float(numhalos[i])
    noprogenfrac['Masscut'][i]=float(numval)/float(sizewdata.size)

    noprogenfrac['Halos'][i]=float(np.where(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]==10)[0].size)/float(numval)
    noprogenfrac['Subs'][i]=float(np.where(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]!=10)[0].size)/float(numval)

    noprogendata['Npart'] = np.concatenate([noprogendata['Npart'],
        halopropdata[i]['npart'][sizewdata][noprogenwdata]])
    noprogendata['Npart_stats'][i][1:6] = np.percentile(halopropdata[i]['npart'][sizewdata][noprogenwdata],[2.5,16.0,50.,84.0,97.5])
    noprogendata['Npart_stats'][i][0] = np.min(halopropdata[i]['npart'][sizewdata][noprogenwdata])
    noprogendata['Npart_stats'][i][-1] = np.max(halopropdata[i]['npart'][sizewdata][noprogenwdata])

    wdatatype = np.where(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]==10)[0]
    if (wdatatype.size > 0):
        noprogendata['Npart_structype_stats']['Halos'][i][1:6] = np.percentile(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype],[2.5,16.0,50.,84.0,97.5])
        noprogendata['Npart_structype_stats']['Halos'][i][0] = np.min(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype])
        noprogendata['Npart_structype_stats']['Halos'][i][-1] = np.max(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype])
    wdatatype = np.where(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]>=20)[0]
    if (wdatatype.size > 0):
        noprogendata['Npart_structype_stats']['Subhalos'][i][1:6] = np.percentile(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype],[2.5,16.0,50.,84.0,97.5])
        noprogendata['Npart_structype_stats']['Subhalos'][i][0] = np.min(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype])
        noprogendata['Npart_structype_stats']['Subhalos'][i][-1] = np.max(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype])
    wdatatype = np.where(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]==15)[0]
    if (wdatatype.size > 0):
        noprogendata['Npart_structype_stats']['Merger Remnants'][i][1:6] = np.percentile(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype],[2.5,16.0,50.,84.0,97.5])
        noprogendata['Npart_structype_stats']['Merger Remnants'][i][0] = np.min(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype])
        noprogendata['Npart_structype_stats']['Merger Remnants'][i][-1] = np.max(halopropdata[i]['npart'][sizewdata][noprogenwdata][wdatatype])
    #subs
    wdata = np.where(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]!=10)[0]
    if (wdata.size > 0):
        tempdata = np.zeros(numval)
        hosts = np.uint64(halopropdata[i]['hostHaloID'][sizewdata][noprogenwdata] % TEMPORALHALOIDVAL - 1)

        tempdata[wdata]=halopropdata[i]['npart'][sizewdata][noprogenwdata][wdata] / halopropdata[i]['npart'][hosts[wdata]]
        noprogendata['Npart_subtohost'] = np.concatenate([noprogendata['Npart_subtohost'],
            tempdata])
        noprogendata['Npart_subtohost_stats'][i][1:6] = np.percentile(tempdata[wdata],[2.5,16.0,50.,84.0,97.5])
        noprogendata['Npart_subtohost_stats'][i][0] = np.min(tempdata[wdata])
        noprogendata['Npart_subtohost_stats'][i][-1] = np.max(tempdata[wdata])

    noprogendata['Sub_type'] = np.concatenate([noprogendata['Sub_type'],
        halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]])
    noprogendata['Zform'] = np.concatenate([noprogendata['Zform'],
        np.ones(noprogenwdata.size)*(1.0/halopropdata[i]['SimulationInfo']['ScaleFactor']-1.0)])
    #hosts
    wdata = np.where(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]==10)[0]
    if (wdata.size > 0):
        tempdata = np.zeros(numval)
        subs = np.where(np.in1d(halopropdata[i]['hostHaloID'], halopropdata[i]['ID'][sizewdata][noprogenwdata][wdata]))[0]
        if (subs.size>0):
            hostmap = dict(zip(halopropdata[i]['ID'][sizewdata][noprogenwdata][wdata], wdata))
            hostmassmap = dict(zip(halopropdata[i]['ID'][sizewdata][noprogenwdata][wdata], halopropdata[i]['npart'][sizewdata][noprogenwdata][wdata]))
            for isub in subs:
                ihost=halopropdata[i]['hostHaloID'][isub]
                hostindex=hostmap[ihost]
                tempdata[hostindex]+=halopropdata[i]['npart'][isub]/hostmassmap[ihost]
            noprogendata['Npart_hosttosubs'] = np.concatenate([noprogendata['Npart_hosttosubs'],
                tempdata])
            #store number of halos with no substructure
            noprogendata['Npart_hosttosubs_stats'][i][7] = float(np.where(tempdata[wdata]==0)[0].size)/float(wdata.size)
            wdata2 = np.where(tempdata>0)[0]
            if (wdata2.size >0):
                tempdata=tempdata[wdata2]
                noprogendata['Npart_hosttosubs_stats'][i][1:6] = np.percentile(tempdata,[2.5,16.0,50.,84.0,97.5])
                noprogendata['Npart_hosttosubs_stats'][i][0] = np.min(tempdata)
                noprogendata['Npart_hosttosubs_stats'][i][6] = np.max(tempdata)

    temptreedata=np.zeros(numval, dtype=np.int64)
    temptreedata2=np.zeros(numval, dtype=np.int64)
    temptreedata3=np.zeros(numval, dtype=np.int64)

    #look at the immediate descendants of these objects
    heads = treedata[i]['Head'][sizewdata][noprogenwdata]
    headssnap = np.int32(np.floor(heads/TEMPORALHALOIDVAL))
    headsindex = np.uint64(heads%TEMPORALHALOIDVAL-1)
    for j in range(numval):
        isnap=headssnap[j]
        index=headsindex[j]
        temptreedata[j] = halopropdata[isnap]['npart'][index]
        temptreedata2[j] = halopropdata[isnap]['Structuretype'][index]
        temptreedata3[j] = treedata[isnap]['Tail'][index]
    temptreedata3=np.int32(temptreedata3==treedata[i]['ID'][sizewdata][noprogenwdata])

    noprogendata['Npart_descen'] = np.concatenate([noprogendata['Npart_descen'],
        temptreedata])
    noprogendata['Sub_type_descen'] = np.concatenate([noprogendata['Sub_type_descen'],
        temptreedata2])
    noprogendata['Descen_type'] = np.concatenate([noprogendata['Descen_type'],
        temptreedata3])
    noprogenfrac['DescenType']['Same'][i]=float(np.where(temptreedata2==halopropdata[i]['Structuretype'][sizewdata][noprogenwdata])[0].size)/float(numval)
    noprogenfrac['DescenType']['Change'][i]=float(np.where(temptreedata2!=halopropdata[i]['Structuretype'][sizewdata][noprogenwdata])[0].size)/float(numval)
    noprogenfrac['DescenType']['SubToHalo'][i]=float(np.where((temptreedata2==10)*(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]>10))[0].size)/float(numval)
    noprogenfrac['DescenType']['HaloToSub'][i]=float(np.where((temptreedata2>10)*(halopropdata[i]['Structuretype'][sizewdata][noprogenwdata]==10))[0].size)/float(numval)

    #check if objects are main branches and get their lifetime as a main branch
    rootheads = treedata[i]['RootHead'][sizewdata][noprogenwdata]
    roottails = treedata[i]['RootTail'][sizewdata][noprogenwdata]
    rootheadssnap = np.int32(np.floor(rootheads/TEMPORALHALOIDVAL))
    rootheadsindex = np.uint64(rootheads % TEMPORALHALOIDVAL-1)
    branchlife = np.zeros(numval,dtype=np.int32)
    for j in range(numval):
        temptreedata[j] = np.int32(treedata[rootheadssnap[j]]['RootTail'][rootheadsindex[j]] == roottails[j])
        curHalo = roottails[j]
        curSnap = np.uint64(curHalo / TEMPORALHALOIDVAL)
        curIndex = np.uint64(curHalo % TEMPORALHALOIDVAL - 1)
        curRootHead = treedata[curSnap]['RootHead'][curIndex]
        curRootTail = treedata[curSnap]['RootTail'][curIndex]
        branchlife[j] = -1
        while (curRootTail == roottails[j]):
            if (curHalo == curRootHead):
                break
            branchlife[j] += 1
            curHalo = treedata[curSnap]['Head'][curIndex]
            curSnap = np.uint64(curHalo / TEMPORALHALOIDVAL)
            curIndex = np.uint64(curHalo % TEMPORALHALOIDVAL - 1)
            curRootTail = treedata[curSnap]['RootTail'][curIndex]

    noprogendata['Branch_type'] = np.concatenate([noprogendata['Branch_type'],
        temptreedata])
    noprogendata['Branch_lifetime'] = np.concatenate([noprogendata['Branch_lifetime'],
        branchlife])
    noprogenfrac['MainBranch'][i]=float(np.where(temptreedata==1)[0].size)/float(numval)
    noprogenfrac['ShortBranch'][i]=float(np.where(branchlife<=3)[0].size)/float(numval)


#plot the information of no progenitors
#plot the fraction as a function of redshift
fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
ytitle=r'$N_{\rm outliers}/N_i$'
xvals=np.arange(numsnaps)
yvals=noprogenfrac['All']
ax.plot(xvals,yvals,color='Navy',ls='solid',lw=4, label='All')
yvals=noprogenfrac['Masscut']
ax.plot(xvals,yvals,color='Crimson',ls='solid',lw=4, label='Mass cut')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(1e-4,1)
ax.set_yscale("log")
ax.set_ylabel(ytitle)
ax.legend()
plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-fraction-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

fig1, ax = plt.subplots(figsize=(10,6))
ytitle=r'$N_i/N_{\rm outliers}$'
xvals=np.arange(numsnaps)
yvals=noprogenfrac['Halos']
ax.plot(xvals,yvals,color='Navy', ls='solid', lw=4, label='Halos')
yvals=noprogenfrac['Subs']
ax.plot(xvals,yvals,color='Crimson', ls='solid', lw=4, label='Subs')
yvals=noprogenfrac['MainBranch']
ax.plot(xvals,yvals,color='DarkGreen', ls='solid', lw=4, label='Main Branch')
yvals=noprogenfrac['ShortBranch']
ax.plot(xvals,yvals,color='LimeGreen', ls='solid', lw=4, label='Short Branch')
yvals=noprogenfrac['DescenType']['Same']
ax.plot(xvals,yvals,color='Magenta', ls='dashed', lw=2, label=r'$S_{\rm type}=S_{\rm type,d}$')
yvals=noprogenfrac['DescenType']['Change']
ax.plot(xvals,yvals,color='DarkOrange', ls='dashed', lw=2, label=r'$S_{\rm type}\neq S_{\rm type,d}$')
yvals=noprogenfrac['DescenType']['HaloToSub']
ax.plot(xvals,yvals,color='Gold', ls='dotted', lw=2, label=r'Halo$\rightarrow$Sub')
yvals=noprogenfrac['DescenType']['SubToHalo']
ax.plot(xvals,yvals,color='Olive', ls='dotted', lw=2, label=r'Sub$\rightarrow$Halo')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(0,1.2)
ax.set_yscale("linear")
ax.set_ylabel(ytitle)
ax.legend(fontsize=16, ncol=3, loc=0)
plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-subpopulations-fraction-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

#plot number of particles distribution as function of redshift
fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
ytitle=r'$N_p(z_{\rm form})$'
labelset=['']
xvals=np.arange(numsnaps)
yvals=noprogendata['Npart_stats'].transpose()
ymean=yvals[3]
ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4, zorder=4)
yel,yeh=yvals[2],yvals[4]
ax.fill_between(xvals,yel,yeh,
    facecolor='Blue',edgecolor='None',alpha=0.75,interpolate=True,zorder=3,linewidth=0,linestyle='None')
yell,yehh=yvals[1],yvals[5]
ax.fill_between(xvals,yell,yehh,
    facecolor='Teal',edgecolor='None',alpha=0.5,interpolate=True,zorder=2,linewidth=0,linestyle='None')
yelll,yehhh=yvals[0],yvals[6]
ax.fill_between(xvals,yelll,yehhh,
    facecolor='None',edgecolor='Cyan',alpha=0.25,interpolate=True,zorder=1,linewidth=2,linestyle='solid')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(100,1e6)
ax.set_yscale("log")
ax.set_ylabel(ytitle)
plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-numdistrib-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

#plot number of particles distribution as function of redshift split into types
selectionlabellist=['Halos','Subhalos','Merger Remnants']
selectionplotlabellist=['halos','subhalos','mergers']
numsel=3
for i in range(numsel):
    fig1, ax = plt.subplots(figsize=(10,6))
    xtitle=r'Snapshots'
    ytitle=r'$N_p(z_{\rm form})$'
    selectionlabel=selectionlabellist[i]
    selectionplotlabel=selectionplotlabellist[i]
    xvals=np.arange(numsnaps)
    yvals=noprogendata['Npart_structype_stats'][selectionlabel].transpose()
    ymean=yvals[3]
    ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4, zorder=4)
    yel,yeh=yvals[2],yvals[4]
    ax.fill_between(xvals,yel,yeh,
        facecolor='Blue',edgecolor='None',alpha=0.75,interpolate=True,zorder=3,linewidth=0,linestyle='None')
    yell,yehh=yvals[1],yvals[5]
    ax.fill_between(xvals,yell,yehh,
        facecolor='Teal',edgecolor='None',alpha=0.5,interpolate=True,zorder=2,linewidth=0,linestyle='None')
    yelll,yehhh=yvals[0],yvals[6]
    ax.fill_between(xvals,yelll,yehhh,
        facecolor='None',edgecolor='Cyan',alpha=0.25,interpolate=True,zorder=1,linewidth=2,linestyle='solid')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.set_xlim(0,numsnaps)
    ax.set_xscale("linear")
    ax.set_xlabel(xtitle)
    ax.set_ylim(100,1e6)
    ax.set_yscale("log")
    ax.set_ylabel(ytitle)
    plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-numdistrib-numsnaps-'+selectionplotlabel+'.pdf')
    fig1.clf()
    plt.close(fig1)

#plot subhalo to host halo fraction as a function of redshift
fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
ytitle=r'$N_{p,{\rm S}}/N_{p,{\rm H}}$'
labelset=['']
xvals=np.arange(numsnaps)
yvals=noprogendata['Npart_subtohost_stats'].transpose()
ymean=yvals[3]
ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4, zorder=4)
yel,yeh=yvals[2],yvals[4]
ax.fill_between(xvals,yel,yeh,
    facecolor='Blue',edgecolor='None',alpha=0.75,interpolate=True,zorder=3,linewidth=0,linestyle='None')
yell,yehh=yvals[1],yvals[5]
ax.fill_between(xvals,yell,yehh,
    facecolor='Teal',edgecolor='None',alpha=0.5,interpolate=True,zorder=2,linewidth=0,linestyle='None')
yelll,yehhh=yvals[0],yvals[6]
ax.fill_between(xvals,yelll,yehhh,
    facecolor='None',edgecolor='Cyan',alpha=0.25,interpolate=True,zorder=1,linewidth=2,linestyle='solid')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(1e-4,1)
ax.set_yscale("log")
ax.set_ylabel(ytitle)
plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-subhalos-numfracdistrib-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

#plot host halo's subhalo mass fraction as a function of redshift
fig1, ax = plt.subplots(figsize=(10,6))
xtitle=r'Snapshots'
ytitle=r'$\sum_i N_{p,{\rm S_i}}/N_{p,{\rm H}}$'
labelset=['']
xvals=np.arange(numsnaps)
yvals=noprogendata['Npart_hosttosubs_stats'].transpose()
ymean=yvals[3]
ax.plot(xvals,ymean,color='Navy',ls='solid',lw=4, zorder=4)
yel,yeh=yvals[2],yvals[4]
ax.fill_between(xvals,yel,yeh,
    facecolor='Blue',edgecolor='None',alpha=0.75,interpolate=True,zorder=3,linewidth=0,linestyle='None')
yell,yehh=yvals[1],yvals[5]
ax.fill_between(xvals,yell,yehh,
    facecolor='Teal',edgecolor='None',alpha=0.5,interpolate=True,zorder=2,linewidth=0,linestyle='None')
yelll,yehhh=yvals[0],yvals[6]
ax.fill_between(xvals,yelll,yehhh,
    facecolor='None',edgecolor='Cyan',alpha=0.25,interpolate=True,zorder=1,linewidth=2,linestyle='solid')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_xlim(0,numsnaps)
ax.set_xscale("linear")
ax.set_xlabel(xtitle)
ax.set_ylim(1e-4,2)
ax.set_yscale("log")
ax.set_ylabel(ytitle)
plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-host-numfracdistrib-numsnaps.pdf')
fig1.clf()
plt.close(fig1)

#plot first the time of occurence and number of particles split along type of subtructure
#coloured by the ratio of descendant mass
xtitle=r'$N_p$'
ytitle=r'$z_{\rm form}$'
ztitle=r'$N_p(z_{\rm descen})/N_p(z_{\rm form})$'
selectionlabellist=['Halos','Subhalos','Merger Remnants']
selectionplotlabellist=['halos','subhalos','mergers']
wdatalist=[np.where(noprogendata['Sub_type']==10)[0],np.where(noprogendata['Sub_type']>=20)[0], np.where(noprogendata['Sub_type']==15)[0]]
numsel=3
for i in range(numsel):
    fig1, ax = plt.subplots(figsize=(10,10))
    selectionlabel=selectionlabellist[i]
    selectionplotlabel=selectionplotlabellist[i]
    wdata=wdatalist[i]
    xvals=noprogendata['Npart'][wdata]
    yvals=noprogendata['Zform'][wdata]
    zvals=noprogendata['Npart_descen'][wdata]/noprogendata['Npart'][wdata]
    zzvals=noprogendata['Branch_lifetime'][wdata]

    #all points
    scatter1 = ax.scatter(xvals, yvals, c=zvals, s=10, marker='o', cmap='jet', norm = matplotlib.colors.LogNorm(vmin=1e-1,vmax=1e1), alpha=0.5, zorder=1, edgecolors='none')
    cbar = fig1.colorbar(scatter1, ax=ax)
    cbar.set_label(ztitle)

    #median and scatter
    xmean=np.median(xvals)
    ymean=np.median(yvals)
    zmean=np.median(zvals)
    zzmean=np.median(zzvals)
    ax.scatter([xmean],[ymean], c=[zmean], s=100, marker='s', cmap='jet', norm = matplotlib.colors.LogNorm(vmin=1e-1,vmax=1e1), alpha=0.9, edgecolors='k', linewidth=2, zorder=3)
    yerrval=np.array(np.fabs(np.percentile(yvals,[16.0,84.0])-ymean))
    xerrval=np.array(np.fabs(np.percentile(xvals,[16.0,84.0])-xmean))
    zerrval=np.array(np.fabs(np.percentile(zvals,[16.0,84.0])-zmean))
    zzerrval=np.array(np.fabs(np.percentile(zzvals,[16.0,84.0])-zzmean))
    ax.errorbar(xmean, ymean, yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])], xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])], color='k',ls='solid',lw=4, alpha=0.5, zorder=2)
    xerrval2=np.array(np.fabs(np.percentile(xvals,[2.5,97.5])-np.median(xvals)))
    yerrval2=np.array(np.fabs(np.percentile(yvals,[2.5,97.5])-np.median(yvals)))
    ax.errorbar(xmean, ymean, yerr=[np.array([yerrval2[0]]),np.array([yerrval2[1]])], xerr=[np.array([xerrval2[0]]),np.array([xerrval2[1]])], color='k',ls='solid',lw=2, alpha=0.25, zorder=2)
    t=ax.annotate(selectionlabel,
        xy=(0.95,0.90),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$N=%1.0f$'%xvals.size,
        xy=(0.95,0.85),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$N_p=%1.0f_{-%1.0f}^{+%1.0f}$'%
        (xmean,xerrval[0],xerrval[1]),
        xy=(0.95,0.80),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$z_{\rm form}=%1.2f_{-%1.2f}^{+%1.2f}$'%
        (ymean,yerrval[0],yerrval[1]),
        xy=(0.95,0.75),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$N_p(z_{\rm descen})/N_p(z_{\rm form})=%1.2f_{-%1.2f}^{+%1.2f}$'%
        (zmean,zerrval[0],zerrval[1]),
        xy=(0.95,0.70),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$\mathcal{B}_l=%1.0f_{-%1.0f}^{+%1.0f}$'%
        (zzmean,zzerrval[0],zzerrval[1]),
        xy=(0.95,0.65),xycoords='axes fraction', horizontalalignment='right')

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.set_xlim(100,1e6)
    ax.set_xscale("log")
    ax.set_xlabel(xtitle)
    ax.set_ylim(0,10)
    ax.set_yscale("linear")
    ax.set_ylabel(ytitle)
    plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-'+selectionplotlabel+'.pdf')
    fig1.clf()
    plt.close(fig1)

#similar to above but color code by branch length and also split by branch type (main or not main)
xtitle=r'$N_p$'
ytitle=r'$z_{\rm form}$'
ztitle=r'$\mathcal{B}_l$'
selectionlabellist=['Halos','Subhalos','Merger Remnants']
selectionplotlabellist=['halos','subhalos','mergers']
wdatalist=[np.where(noprogendata['Sub_type']==10)[0],np.where(noprogendata['Sub_type']>=20)[0], np.where(noprogendata['Sub_type']==15)[0]]
numsel=3
for i in range(numsel):
    fig1, ax = plt.subplots(figsize=(10,10))
    selectionlabel=selectionlabellist[i]
    selectionplotlabel=selectionplotlabellist[i]
    wdata=wdatalist[i]
    xvals=noprogendata['Npart'][wdata]
    yvals=noprogendata['Zform'][wdata]
    zvals=noprogendata['Branch_lifetime'][wdata]
    zzvals=noprogendata['Branch_type'][wdata]
    wdata2=np.where(noprogendata['Branch_type'][wdata]==1)[0]

    #all points
    scatter1 = ax.scatter(xvals, yvals, c=zvals, s=10, marker='o', cmap='jet', norm = matplotlib.colors.Normalize(vmin=0,vmax=10), alpha=0.5, zorder=1, edgecolors='none')
    cbar = fig1.colorbar(scatter1, ax=ax)
    cbar.set_label(ztitle)
    scatter2 = ax.scatter(xvals[wdata2], yvals[wdata2], facecolor='None', edgecolor='black', linewidth=2, s=10, marker='o', alpha=0.5, zorder=2)

    #median and scatter
    xmean=np.median(xvals)
    ymean=np.median(yvals)
    zmean=np.median(zvals)
    zzmean=np.median(zzvals)
    ax.scatter([xmean],[ymean], c=[zmean], s=100, marker='s', cmap='jet', norm = matplotlib.colors.Normalize(vmin=0,vmax=10), alpha=0.9, edgecolors='k', linewidth=2, zorder=3)
    yerrval=np.array(np.fabs(np.percentile(yvals,[16.0,84.0])-ymean))
    xerrval=np.array(np.fabs(np.percentile(xvals,[16.0,84.0])-xmean))
    zerrval=np.array(np.fabs(np.percentile(zvals,[16.0,84.0])-zmean))
    zzerrval=np.array(np.fabs(np.percentile(zzvals,[16.0,84.0])-zzmean))
    ax.errorbar(xmean, ymean, yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])], xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])], color='k',ls='solid',lw=4, alpha=0.5, zorder=2)
    xerrval2=np.array(np.fabs(np.percentile(xvals,[2.5,97.5])-np.median(xvals)))
    yerrval2=np.array(np.fabs(np.percentile(yvals,[2.5,97.5])-np.median(yvals)))
    ax.errorbar(xmean, ymean, yerr=[np.array([yerrval2[0]]),np.array([yerrval2[1]])], xerr=[np.array([xerrval2[0]]),np.array([xerrval2[1]])], color='k',ls='solid',lw=2, alpha=0.25, zorder=2)
    t=ax.annotate(selectionlabel,
        xy=(0.95,0.90),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$N=%1.0f$'%xvals.size,
        xy=(0.95,0.85),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$f_{\rm main branch}=%1.2f$'%
        (float(wdata2.size)/float(wdata.size)),
        xy=(0.95,0.80),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$N_p=%1.0f_{-%1.0f}^{+%1.0f}$'%
        (xmean,xerrval[0],xerrval[1]),
        xy=(0.95,0.75),xycoords='axes fraction', horizontalalignment='right')
    t=ax.annotate(r'$\mathcal{B}_l=%1.0f_{-%1.0f}^{+%1.0f}$'%
        (zmean,zerrval[0],zerrval[1]),
        xy=(0.95,0.70),xycoords='axes fraction', horizontalalignment='right')

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.set_xlim(100,1e6)
    ax.set_xscale("log")
    ax.set_xlabel(xtitle)
    ax.set_ylim(0,10)
    ax.set_yscale("linear")
    ax.set_ylabel(ytitle)
    plt.savefig(baseplotdir+labelname+'-noprogenitor-outliers-branchproperties-'+selectionplotlabel+'.pdf')
    fig1.clf()
    plt.close(fig1)

#Fakhouri style merger rates
#split z=0 halos/subhalos into several mass bins,
if (ExtraDiagnostics >= 1):
    dlogMhost,dz,dlogxi=0.25,0.5,0.25

    massbinedges=10.0**np.arange(10.5-0.5*dlogMhost,13.5+0.5*dlogMhost,dlogMhost)
    massbins=10**(0.5*(np.log10(massbinedges[:-1])+np.log10(massbinedges[1:])))
    dMhost=np.log(10.0)*massbins*dlogMhost
    nmassbins=len(massbins)
    scalefactorbinedges=1.0/(1.0+np.arange(0,1.0+0.5*dz,dz))[::-1]
    scalefactorbins=0.5*(scalefactorbinedges[:-1]+scalefactorbinedges[1:])
    nscalefactorbins=len(scalefactorbins)
    massratiobinedges=10**np.arange(-3.5-0.5*dlogxi,0.1+0.5*dlogxi,dlogxi)
    massratiobins=10**(0.5*(np.log10(massratiobinedges[:-1])+np.log10(massratiobinedges[1:])))
    dxi=np.log(10.0)*massratiobins*dlogxi
    nmassratiobins=len(massratiobins)
    #here, lets choose a hostHalo and a subhalo
    npartlim=40
    maxnsample=10000
    unitfac=UnitInfo['Mass_unit_to_solarmass']
    print('Merger rates')
    print('mass bins ',nmassbins,np.log10(massbins),dMhost)
    print('redshift bins',nscalefactorbins,scalefactorbins,1.0/scalefactorbins-1.0)
    print('mass ratio bins ',nmassratiobins, dxi)

    numhostsall=np.zeros([nmassbins,nscalefactorbins],dtype=np.int32)
    numdataall=np.zeros([nmassbins,nscalefactorbins,nmassratiobins],dtype=np.float64)
    mergerratedataall=np.zeros([nmassbins,nscalefactorbins,nmassratiobins])
    mergerratedatanormalizedall=np.zeros([nmassbins,nscalefactorbins,nmassratiobins])
    mergerratestatsall=np.zeros([nmassbins,nscalefactorbins,nmassratiobins,7])
    mergerratestatsnormalizedall=np.zeros([nmassbins,nscalefactorbins,nmassratiobins,7])

    start0=time.clock()
    for ibin in range(nmassbins):
        start=time.clock()
        for iscalebin in range(nscalefactorbins):
            numhostsatscalefactor=np.zeros(numsnaps,dtype=np.int32)
            snaps=np.where((atime>scalefactorbinedges[iscalebin])*(atime<=scalefactorbinedges[iscalebin+1]))[0]
            hostids=np.array([],dtype=np.int64)
            for curSnap in snaps:
                wdata=np.where((halopropdata[curSnap]['Mass_tot']*unitfac>=massbinedges[ibin])*(halopropdata[curSnap]['Mass_tot']*unitfac<massbinedges[ibin+1]))
                numhostsatscalefactor[curSnap]=len(wdata[0])
                hostids=np.concatenate([hostids,np.int64(halopropdata[curSnap]['ID'][wdata])])
            numhostsall[ibin][iscalebin]=len(hostids)
            if (numhostsall[ibin][iscalebin]==0):
                continue
            #if number of objects is very large subsample
            if (numhostsall[ibin][iscalebin]>=maxnsample):
                hostids=np.random.choice(hostids,maxnsample,replace=False)
                numhostsall[ibin][iscalebin]=maxnsample
            print('looking at bin ',ibin,massbins[ibin],'scale factor bin', iscalebin, 'initially containing',numhostsall[ibin][iscalebin])
            numactive=0
            start2=time.clock()
            for ihost in hostids:
                curHalo=ihost
                curSnap=np.uint64(curHalo/TEMPORALHALOIDVAL)
                curIndex=np.uint64(curHalo%TEMPORALHALOIDVAL-1)
                curTail=treedata[curSnap]['Tail'][curIndex]
                if (curHalo==curTail): continue
                curMass=halopropdata[curSnap]['Mass_tot'][curIndex]
                #if looking at true complete mergers then want objects that fully merge with object
                #find all progenitors of this (sub)halo
                progens=np.array([],dtype=np.int64)
                for snap in np.arange(curSnap-1,curSnap-NSNAPSEARCH-1,-1,dtype=np.int32):
                    wdata=np.where((halopropdata[snap]["Efrac"]>=0.1)*(halopropdata[snap]["npart"]>=npartlim)*(treedata[snap]["Head"]==curHalo)*(halopropdata[snap]["ID"]!=curTail))
                    progens=np.concatenate([progens,halopropdata[snap]['Mass_tot'][wdata]])
                if (progens.size==0): continue
                primaryprogenMass=(curMass-np.sum(progens))
                if (primaryprogenMass==0): continue
                #number of mergers within this redshift bin
                y,x=np.histogram(progens/primaryprogenMass,bins=massratiobinedges)
                mergerratedataall[ibin][iscalebin]+=np.float32(y)
                #number of mergers normalized by number of halos at given mass bin at current redshift
                mergerratedatanormalizedall[ibin][iscalebin]+=np.float32(y)/float(numhostsatscalefactor[curSnap])
                numactive+=1
            print('finished snap ', numactive, np.sum(numhostsatscalefactor),np.max(mergerratedatanormalizedall[ibin][iscalebin]),' bin in',time.clock()-start2)
            #once have done all snaps, get statistics
            tempdata=mergerratedataall[ibin][iscalebin]
            numdataall[ibin][iscalebin]=tempdata
            for iratio in range(nmassratiobins):
                normfac=1.0/(SimulationInfo['Period']**3.0*dz*dxi[iratio]*dMhost[ibin])
                #normalize by comoving volume, redshift bin size, mass ratio bin size and mass bin size
                mergerratestatsall[ibin][iscalebin][iratio][2:]=np.percentile(tempdata[iratio],[2.5,16.0,50.,84.0,97.5])*normfac
                mergerratestatsall[ibin][iscalebin][iratio][0]=np.average(tempdata[iratio])*normfac
                mergerratestatsall[ibin][iscalebin][iratio][1]=np.std(tempdata[iratio])*normfac
            tempdata=mergerratedatanormalizedall[ibin][iscalebin]
            for iratio in range(nmassratiobins):
                normfac=1.0/(dz*dxi[iratio])
                mergerratestatsnormalizedall[ibin][iscalebin][iratio][2:]=np.percentile(tempdata[iratio],[2.5,16.0,50.,84.0,97.5])*normfac
                mergerratestatsnormalizedall[ibin][iscalebin][iratio][0]=np.average(tempdata[iratio])*normfac
                mergerratestatsnormalizedall[ibin][iscalebin][iratio][1]=np.std(tempdata[iratio])*normfac
        print('finished mass bin in',time.clock()-start)
    print('finished processing all',time.clock()-start0)

    print('save data')
    fname=baseplotdir+'mergerrate-data.hdf5'
    hfile=h5py.File(fname,'w')
    hfile.create_dataset('massratiobins',data=massratiobins)
    hfile.create_dataset('numhostsall',data=numhostsall)
    hfile.create_dataset('numdataall',data=numdataall)
    hfile.create_dataset('mergerratedataall',data=mergerratedataall)
    hfile.create_dataset('mergerratedatanormalizedall',data=mergerratedatanormalizedall)
    hfile.create_dataset('mergerratestatsall',data=mergerratestatsall)
    hfile.create_dataset('mergerratestatsnormalizedall',data=mergerratestatsnormalizedall)
    hfile.close()



    def xifunc(xi,params):
        [xinorm,beta,gamma]=params
        return xi**beta*np.exp((xi/xinorm)**gamma)
    def meanmergerrate(xi,M,z,params):
        xiparams=params[3:]
        [A,alpha,eta]=params[:3]
        return A*(M/1e12)**alpha*(1+z)**eta*xifunc(xi,xiparams)
    FMparams=[0.0104,0.133,0.0993,0.00972,-1.995,0.263]
    Pooleparams=[0.104,0.122,0.069,0.483,-1.751,0.67]
    Poolesubparams=[0.039,0.254,0.872,0.483,-1.751,0.67]

    #fig1, ax = plt.subplots(figsize=(10,6))
    #xtitle=r'$\xi\equiv M_{P,s,{\rm tot}}/M_{P,p,{\rm tot}}$'
    #ytitle=r'$\frac{dN}{dzd\xi}\equiv \frac{B(M_D,\xi,z)}{n(M_D,z)}$'
    #ztitle=r'$M_D [h^{-1}{\rm M}_\odot]$'
    #cmap = cm.autumn
    #gs1=matplotlib.gridspec.GridSpec(1,20, hspace=0.0,wspace=0.0)
    #ax=plt.subplot(gs1[0, 0:18])
    #iscalebin=1
    #for ibin in range(nmassbins):
    #    yvals=mergerratestatsnormalizedall[ibin][iscalebin].transpose()[0]
    #    wdata=np.where(yvals>0)[0][1:]
    #    ax.plot(massratiobins[wdata],yvals[wdata],color=cmap(ibin/float(nmassbins)),lw=5,ls='solid',zorder=12,alpha=0.7)
    #modelmassratiobins=10**np.arange(-4,0.05,0.05)
    #ax.plot(modelmassratiobins,meanmergerrate(modelmassratiobins,10**12.5/SimulationInfo['h_val'],1.0/scalefactorbins[iscalebin]-1.0,FMparams),zorder=9,color='gray',alpha=0.9,lw=4,ls='dashed',label=r'FM10')
    #ax.plot(modelmassratiobins,meanmergerrate(modelmassratiobins,10**12.5/SimulationInfo['h_val'],1.0/scalefactorbins[iscalebin]-1.0,Pooleparams),zorder=8,color='g',alpha=0.9,lw=4,ls='solid',label='P17 FOF')
    #ax.plot(modelmassratiobins,meanmergerrate(modelmassratiobins,10**12.5/SimulationInfo['h_val'],1.0/scalefactorbins[iscalebin]-1.0,Poolesubparams),zorder=9,color='Blue',alpha=0.9,lw=4,ls='solid',label='P17 sub')
    #ax.plot(modelmassratiobins,1.0/1.3*meanmergerrate(modelmassratiobins,10**12.5/SimulationInfo['h_val']*0.87,1.0/scalefactorbins[iscalebin]-1.0,Poolesubparams),zorder=9,color='Cyan',alpha=0.5,lw=6,ls='dotted',label='modified P17 sub')
    #ax.xaxis.set_ticks_position('both')
    #ax.yaxis.set_ticks_position('both')
    #ax.set_xlim([5e-4,1.1])
    #ax.set_xscale("log")
    #ax.set_xlabel(xtitle)
    #ax.set_ylim([0.075,5e4])
    #ax.set_yscale("log")
    #ax.set_ylabel(ytitle)
    #leg1=legend(fontsize=16)
    #plt.gca().add_artist(leg1)
    #custom_lines = [
    #                Line2D([0], [0], color=cmap(0.5), lw=5,ls='solid'),
    #                Line2D([0], [0], color=cmap(0.5), lw=3,ls='dashed',path_effects=[pe.Stroke(linewidth=6, foreground='k'),pe.Normal()]),
    #               ]
    #custom_labels=[r'Mergers', r'Accretion (FOF Mergers)', ]
    #leg2=legend(custom_lines,custom_labels, fontsize=16, ncol=1,loc='lower left',numpoints=4)
    #
    #ax=plt.subplot(gs1[0:2, 18:19])
    #cb1=mpl.colorbar.ColorbarBase(ax=ax,cmap='autumn',norm=matplotlib.colors.LogNorm(vmin=3.16e10,vmax=3.16e13), orientation='vertical')
    #cb1.set_label(ztitle)
    #plt.savefig(baseplotdir+'mergerrate.pdf')

    #get mass accretion history
    #split z=0 halos into several mass bins,
    #massbinedges=10.0**np.array([9.5,10.5,11.5,12.5,13.5,14.5])
    #massbins=10**np.array([10,11,12,13,14])

    massbinedges=10.0**np.array([10.5,11.5,12.5,13.5])
    massbins=10**np.array([11,12,13])
    nmassbins=len(massbins)

    #here, lets choose a hostHalo and a subhalo
    numhosts=np.zeros(nmassbins,dtype=np.int32)
    numnomergers=np.zeros(nmassbins)
    numfewmergers=np.zeros(nmassbins)
    nummostlysubhalo=np.zeros(nmassbins)
    nummostlyhalo=np.zeros(nmassbins)
    accmassratiostats=np.zeros([nmassbins,numsnaps,7])
    accmajormergermassratiostats=np.zeros([nmassbins,numsnaps,7])
    accminormergermassratiostats=np.zeros([nmassbins,numsnaps,7])
    minormergermassratiostats=np.zeros([nmassbins,numsnaps,6])
    majormergermassratiostats=np.zeros([nmassbins,numsnaps,6])
    cumaccmassratiostats=np.zeros([nmassbins,6])
    cumaccmajormergermassratiostats=np.zeros([nmassbins,6])
    cumaccminormergermassratiostats=np.zeros([nmassbins,6])

    maxnsample=2000

    for ibin in range(nmassbins):
        start=time.clock()
        curSnap=-1
        hostindices=np.where((halopropdata[curSnap]['Mass_200crit']*UnitInfo['Mass_unit_to_solarmass']*SimulationInfo['h_val']>=massbinedges[ibin])*
            (halopropdata[curSnap]['Mass_200crit']*UnitInfo['Mass_unit_to_solarmass']*SimulationInfo['h_val']<massbinedges[ibin+1])*(halopropdata[curSnap]['hostHaloID']==-1))[0]
        numhosts[ibin]=len(hostindices)
        if (numhosts[ibin]==0):
            continue
        if (numhosts[ibin]>maxnsample):
            hostindices=np.random.choice(hostindices,size=maxnsample,replace=False)
            numhosts[ibin]=maxnsample
        print('looking at bin ',massbins[ibin],'initially containing',numhosts[ibin])
        numaccdata=np.zeros([numsnaps,numhosts[ibin]])
        numaccmajormergerdata=np.zeros([numsnaps,numhosts[ibin]])
        numaccminormergerdata=np.zeros([numsnaps,numhosts[ibin]])
        nummajormergerdata=np.zeros([numsnaps,numhosts[ibin]])
        numminormergerdata=np.zeros([numsnaps,numhosts[ibin]])
        accmassratiodata=np.zeros([numsnaps,numhosts[ibin]])
        accmajormergermassratiodata=np.zeros([numsnaps,numhosts[ibin]])
        accminormergermassratiodata=np.zeros([numsnaps,numhosts[ibin]])
        minormergermassratiodata=np.zeros([numsnaps,numhosts[ibin]])
        majormergermassratiodata=np.zeros([numsnaps,numhosts[ibin]])
        cummassratiodata=np.zeros([numsnaps,numhosts[ibin]])
        cummajormergermassratiodata=np.zeros([numsnaps,numhosts[ibin]])
        cumminormergermassratiodata=np.zeros([numsnaps,numhosts[ibin]])

        ihostcount=0
        for ihost in hostindices:
            curSnap=-1
            #for a given host, find all objects that fully merge with it.
            roothosthead=treedata[curSnap]['RootHead'][ihost]
            roothosttail=treedata[curSnap]['RootTail'][ihost]
            #now can trace this objects to see how they contribute to the total mass of the object
            #lets get the position of the host halo main branch
            halolifetime=GetHaloLiftime(treedata,np.uint64(roothosthead%TEMPORALHALOIDVAL-1),np.uint32(roothosthead/TEMPORALHALOIDVAL),roothosttail,atime,TEMPORALHALOIDVAL)
            #if object is not long lived ignore,
            if (halolifetime<10):
                continue
            hostatime=np.zeros(halolifetime)
            hostmasstot=np.zeros(halolifetime)
            #print(ihost,halolifetime)
            curHalo=roothosttail
            curSnap=np.uint64(curHalo/TEMPORALHALOIDVAL)
            curIndex=np.uint64(curHalo%TEMPORALHALOIDVAL-1)
            curHead=treedata[curSnap]["Head"][curIndex]
            halolifetime=0
            assubhalolifetime=0
            while (True):
                hostatime[halolifetime]=atime[curSnap]
                hostmasstot[halolifetime]=halopropdata[curSnap]['Mass_tot'][curIndex]
                curhostid=halopropdata[curSnap]['hostHaloID'][curIndex]
                if (curhostid==-1):
                    subs=np.where((halopropdata[curSnap]['hostHaloID']==halopropdata[curSnap]['ID'][curIndex]))
                    hostmasstot[halolifetime]+=np.sum(halopropdata[curSnap]['Mass_tot'][subs])
                halolifetime+=1
                assubhalolifetime+=1*(halopropdata[curSnap]['hostHaloID'][curIndex]!=-1)
                if (curHalo==curHead): break
                curHalo=curHead
                curIndex=int(curHalo%TEMPORALHALOIDVAL-1)
                curSnap=int(curHalo/TEMPORALHALOIDVAL)
                curHead=treedata[curSnap]["Head"][curIndex]
            #if object spends most of its time as a subhalo ignore
            if (assubhalolifetime>0.5*halolifetime):
                nummostlysubhalo[ibin]+=1.0
                continue
            nummostlyhalo[ibin]+=1.0
            hostmasstot=scipyinterp.interp1d(hostatime,np.log10(hostmasstot))

            #move to things that merge
            tails=np.array([],dtype=np.uint64)
            for snap in range(numsnaps-1):
                if numhalos[snap] == 0 : continue
                wdata=np.where((treedata[snap]['RootHead']==roothosthead)*(treedata[snap]['RootTail']!=roothosttail))
                if len(wdata[0]) == 0 : continue
                tails=np.concatenate([tails,np.uint64(treedata[snap]['RootTail'][wdata])])
            #ignore if there are no secondary progenitors
            if (tails.size == 0):
                numnomergers[ibin]+=1.0
                continue
            tails=np.unique(tails)
            #if the number of unique secondary progenitors is low, ignore
            if (tails.size < 5):
                numfewmergers[ibin]+=1.0
                continue

            #to store the merger channels
            hostRootTailSnap=np.int32(roothosttail/TEMPORALHALOIDVAL)
            totalhalolifetime=numsnaps-hostRootTailSnap
            #atimedata[ihostcount]=atime[hostRootTailSnap:]
            #accmassratiodata[ihostcount]=np.zeros(totalhalolifetime)
            #minormergermassratiodata[ihostcount]=np.zeros(totalhalolifetime)
            #majormergermassratiodata[ihostcount]=np.zeros(totalhalolifetime)

            for itail in np.uint64(tails):
                curHalo=itail
                curSnap=np.uint64(curHalo/TEMPORALHALOIDVAL)
                curIndex=np.uint64(curHalo%TEMPORALHALOIDVAL-1)
                curHead=treedata[curSnap]["Head"][curIndex]
                curRootTail=treedata[curSnap]["RootTail"][curIndex]
                curHost=halopropdata[curSnap]["hostHaloID"][curIndex]
                curMass=halopropdata[curSnap]["Mass_tot"][curIndex]
                refRootTail=curRootTail
                accMass=0
                mergeMass=0
                accSnap=curSnap
                mergeSnap=curSnap
                accScaleFactor=atime[curSnap]
                mergeScaleFactor=atime[curSnap]
                prevMass=curMass
                prevSnap=curSnap
                #ignore things that start off as subhalos
                #if (curHost!=-1): continue
                halolifetime=0
                while (True):
                    #if object is subhalo, check to see if its host is the root host, store the accretion mass
                    if (curHost!=-1):
                        curHostIndex=np.uint64(curHost%TEMPORALHALOIDVAL-1)
                        curHostRootTail=treedata[curSnap]["RootTail"][curHostIndex]
                        if (curHostRootTail==roothosttail and accMass==0):
                            accMass=prevMass
                            accSnap=curSnap

                    if (curHalo==curHead): break
                    #if object nolonger main branch then has merged with something
                    if (refRootTail!=curRootTail):
                        if (curRootTail==roothosttail):
                            mergeMass=prevMass
                            mergeSnap=curSnap
                        break
                    prevMass=curMass
                    prevSnap=curSnap
                    curHalo=curHead
                    curIndex=int(curHalo%TEMPORALHALOIDVAL-1)
                    curSnap=int(curHalo/TEMPORALHALOIDVAL)
                    curRootTail=treedata[curSnap]["RootTail"][curIndex]
                    curHost=halopropdata[curSnap]["hostHaloID"][curIndex]
                    curHead=treedata[curSnap]["Head"][curIndex]
                    curMass=halopropdata[curSnap]["Mass_tot"][curIndex]
                if (accMass==0): accMass,accSnap=mergeMass,mergeSnap
                sublossMass=max([0,accMass-mergeMass])
                if (mergeMass>0):
                    accMassRatio=accMass/np.power(10.0,hostmasstot(atime[accSnap]))
                    numaccdata[np.int32(accSnap)][ihostcount]+=1
                    accmassratiodata[np.int32(accSnap)][ihostcount]+=accMassRatio
                    mergeMassRatio=mergeMass/np.power(10.0,hostmasstot(atime[mergeSnap]))
                    cummassratiodata[np.int32(accSnap)][ihostcount]+=accMass/np.power(10.0,hostmasstot(1.0))
                    if (accMassRatio>=5e-2):
                        majormergermassratiodata[np.int32(mergeSnap)][ihostcount]+=mergeMassRatio
                        nummajormergerdata[np.int32(accSnap)][ihostcount]+=1
                        accmajormergermassratiodata[np.int32(mergeSnap)][ihostcount]+=accMassRatio
                        numaccmajormergerdata[np.int32(accSnap)][ihostcount]+=1
                        cummajormergermassratiodata[np.int32(accSnap)][ihostcount]+=accMass/np.power(10.0,hostmasstot(1.0))
                    else:
                        minormergermassratiodata[np.int32(mergeSnap)][ihostcount]+=mergeMassRatio
                        numminormergerdata[np.int32(accSnap)][ihostcount]+=1
                        accminormergermassratiodata[np.int32(mergeSnap)][ihostcount]+=accMassRatio
                        numaccminormergerdata[np.int32(accSnap)][ihostcount]+=1
                        cumminormergermassratiodata[np.int32(accSnap)][ihostcount]+=accMass/np.power(10.0,hostmasstot(1.0))
            #print('host ',ihost)
            #print(accmassratiodata)
            #print(majormergermassratiodata)
            #print(minormergermassratiodata)
            ihostcount+=1
        #finished bin, get average across cosmic time
        numaccdata=numaccdata.transpose()[:ihostcount].transpose()
        numaccmajormergerdata=numaccmajormergerdata.transpose()[:ihostcount].transpose()
        numaccminormergerdata=numaccminormergerdata.transpose()[:ihostcount].transpose()
        nummajormergerdata=nummajormergerdata.transpose()[:ihostcount].transpose()
        numminormergerdata=numminormergerdata.transpose()[:ihostcount].transpose()
        accmassratiodata=accmassratiodata.transpose()[:ihostcount].transpose()
        accminormergermassratiodata=accminormergermassratiodata.transpose()[:ihostcount].transpose()
        accmajormergermassratiodata=accmajormergermassratiodata.transpose()[:ihostcount].transpose()
        minormergermassratiodata=minormergermassratiodata.transpose()[:ihostcount].transpose()
        majormergermassratiodata=majormergermassratiodata.transpose()[:ihostcount].transpose()
        cummassratiodata=cummassratiodata.transpose()[:ihostcount].transpose()
        cumminormergermassratiodata=cumminormergermassratiodata.transpose()[:ihostcount].transpose()
        cummajormergermassratiodata=cummajormergermassratiodata.transpose()[:ihostcount].transpose()
        numhosts[ibin]=ihostcount
        print('finished with bin',ibin,'containing',numhosts[ibin],' and total number of events ',np.sum(numaccdata), 'in ',time.clock()-start)
        print('number of objects with few mergers are ',numfewmergers[ibin],'no mergers',numnomergers[ibin],'and spend their life mostly as a subhalo',nummostlysubhalo[ibin])
        for snap in range(numsnaps):
            if (np.sum(numaccdata[snap])>=2):
                accmassratiostats[ibin][snap][0]=np.sum(numaccdata[snap])
                accmassratiostats[ibin][snap][1:3]=np.percentile(accmassratiodata[snap],[16.0,84.0])
                accmassratiostats[ibin][snap][4:6]=np.percentile(accmassratiodata[snap][np.where(accmassratiodata[snap]>0)],[16.0,84.0])
                accmassratiostats[ibin][snap][3]=np.average(accmassratiodata[snap])
                accmassratiostats[ibin][snap][6]=np.std(accmassratiodata[snap])
            if (np.sum(numaccmajormergerdata[snap])>=1):
                accmajormergermassratiostats[ibin][snap][0]=np.sum(numaccmajormergerdata[snap])
                accmajormergermassratiostats[ibin][snap][1:3]=np.percentile(accmajormergermassratiodata[snap],[16.0,84.0])
                accmajormergermassratiostats[ibin][snap][4:6]=np.percentile(accmajormergermassratiodata[snap][np.where(accmassratiodata[snap]>0)],[16.0,84.0])
                accmajormergermassratiostats[ibin][snap][3]=np.average(accmajormergermassratiodata[snap])
                accmajormergermassratiostats[ibin][snap][6]=np.std(accmajormergermassratiodata[snap])
            if (np.sum(numaccminormergerdata[snap])>=2):
                accminormergermassratiostats[ibin][snap][0]=np.sum(numaccminormergerdata[snap])
                accminormergermassratiostats[ibin][snap][1:3]=np.percentile(accminormergermassratiodata[snap],[16.0,84.0])
                accminormergermassratiostats[ibin][snap][4:6]=np.percentile(accminormergermassratiodata[snap][np.where(accmassratiodata[snap]>0)],[16.0,84.0])
                accminormergermassratiostats[ibin][snap][3]=np.average(accminormergermassratiodata[snap])
                accminormergermassratiostats[ibin][snap][6]=np.std(accminormergermassratiodata[snap])
            if (np.sum(nummajormergerdata[snap])>=1):
                majormergermassratiostats[ibin][snap][0]=np.sum(nummajormergerdata[snap])
                majormergermassratiostats[ibin][snap][1:]=np.nanpercentile(majormergermassratiodata[snap],[2.5,16.0,50.0,84.0,97.5])
            if (np.sum(numminormergerdata[snap])>=2):
                minormergermassratiostats[ibin][snap][0]=np.sum(numminormergerdata[snap])
                minormergermassratiostats[ibin][snap][1:]=np.nanpercentile(minormergermassratiodata[snap],[2.5,16.0,50.0,84.0,97.5])

            cumstat=cummassratiodata.transpose()
            cumstat=np.array([np.sum(cumval) for cumval in cumstat])
            cumaccmassratiostats[ibin][1:]=np.percentile(cumstat,[2.5,16.0,50.,84.0,97.5])

            cumstat=cummajormergermassratiodata.transpose()
            cumstat=np.array([np.sum(cumval) for cumval in cumstat])
            cumaccmajormergermassratiostats[ibin][1:]=np.percentile(cumstat,[2.5,16.0,50.,84.0,97.5])

            cumstat=cumminormergermassratiodata.transpose()
            cumstat=np.array([np.sum(cumval) for cumval in cumstat])
            cumaccminormergermassratiostats[ibin][1:]=np.percentile(cumstat,[2.5,16.0,50.,84.0,97.5])
    #now save data
    print('save data')
    fname=baseplotdir+labelname+'-mergergrowth-data.hdf5'
    hfile=h5py.File(fname,'w')
    hfile.create_dataset('numaccdata',data=numaccdata)
    hfile.create_dataset('numaccmajormergerdata',data=numaccmajormergerdata)
    hfile.create_dataset('numaccminormergerdata',data=numaccminormergerdata)
    hfile.create_dataset('nummajormergerdata',data=nummajormergerdata)
    hfile.create_dataset('numminormergerdata',data=numminormergerdata)
    hfile.create_dataset('accmassratiodata',data=accmassratiodata)
    hfile.create_dataset('accminormergermassratiodata',data=accminormergermassratiodata)
    hfile.create_dataset('accmajormergermassratiodata',data=accmajormergermassratiodata)
    hfile.create_dataset('minormergermassratiodata',data=minormergermassratiodata)
    hfile.create_dataset('majormergermassratiodata',data=majormergermassratiodata)
    hfile.create_dataset('cummassratiodata',data=cummassratiodata)
    hfile.create_dataset('cumminormergermassratiodata',data=cumminormergermassratiodata)
    hfile.create_dataset('cummajormergermassratiodata',data=cummajormergermassratiodata)
    hfile.create_dataset('cumaccmassratiostats',data=cumaccmassratiostats)
    hfile.create_dataset('cumaccmajormergermassratiostats',data=cumaccmajormergermassratiostats)
    hfile.create_dataset('cumaccminormergermassratiostats',data=cumaccminormergermassratiostats)
    hfile.close()

if (ExtraDiagnostics >= 2):
    #merger statistics averaged over cosmic time
    snaplist=np.arange(numsnaps,dtype=np.int32)[::-1]
    #allocate enough memory for all halos ???
    relraddata=np.zeros(np.sum(numhalos))
    numpartsdata={'Merging':np.zeros(np.sum(numhalos)),'Host':np.zeros(np.sum(numhalos))}
    print('Getting radial and particle distribution of merging events')
    noffset = 0
    maxnsample=50000
    for i in snaplist:
        print('starting at snapshot ',i)
        start = time.clock()
        activehostindex=np.array([],dtype=np.int64)
        activehostsnap=np.array([],dtype=np.int32)
        activeprogenindex=np.array([],dtype=np.int64)
        activeprogensnap=np.array([],dtype=np.int32)
        nactive = 0
        for ihost in halopropdata[i]['ID']:
            curHalo = ihost
            curSnap = np.uint64(curHalo / TEMPORALHALOIDVAL)
            curIndex = np.uint64(curHalo % TEMPORALHALOIDVAL - 1)
            curTail = treedata[curSnap]['Tail'][curIndex]
            if (curHalo == curTail):
                continue
            #find all progenitors of this (sub)halo
            progenindex=np.array([], dtype=np.int64)
            progensnap=np.array([], dtype=np.int32)
            for snap in np.arange(curSnap-1, curSnap-NSNAPSEARCH,-1, dtype=np.int32):
                if (numhalos[snap] == 0):
                    continue
                wdata = np.where((halopropdata[snap]["Efrac"] >= 0.1)*
                    (treedata[snap]["Head"] == curHalo)*
                    (halopropdata[snap]["ID"] != curTail))[0]
                progenindex = np.concatenate([progenindex,wdata])
                progensnap = np.concatenate([progensnap,snap*np.ones(wdata.size,dtype=np.int32)])
            if (progenindex.size == 0):
                continue
            #get relative position, radii and npart for both primary and secondary
            relpos = np.zeros([progenindex.size,3])
            npartmerging = np.zeros(progenindex.size)
            for icount in range(progenindex.size):
                iprog,isnap = progenindex[icount],progensnap[icount]
                relpos[icount] = np.array([halopropdata[isnap]['Xc'][iprog] - halopropdata[curSnap]['Xc'][curIndex],
                        halopropdata[isnap]['Yc'][iprog] - halopropdata[curSnap]['Yc'][curIndex],
                        halopropdata[isnap]['Zc'][iprog] - halopropdata[curSnap]['Zc'][curIndex]])
                relpos[icount][np.where(relpos[icount]>0.5*SimulationInfo['Period'])] -= SimulationInfo['Period']
                relpos[icount][np.where(relpos[icount]<-0.5*SimulationInfo['Period'])] += SimulationInfo['Period']
                npartmerging[icount] = halopropdata[isnap]['npart'][iprog]
            #print(progenindex.size,nactive)
            #relraddata = np.concatenate([relraddata,np.linalg.norm(relpos,axis=1)/(halopropdata[curSnap]['R_200crit'][curIndex])])
            #numpartsdata['Merging'] = np.concatenate([numpartsdata['Merging'],npartmerging])
            #numpartsdata['Host'] = np.concatenate([numpartsdata['Host'],np.ones(progenindex.size)*halopropdata[curSnap]['npart'][curIndex]])
            relraddata[noffset:noffset+progenindex.size] = np.linalg.norm(relpos,axis=1)/(halopropdata[curSnap]['R_200crit'][curIndex])
            numpartsdata['Merging'][noffset:noffset+progenindex.size] = npartmerging
            numpartsdata['Host'][noffset:noffset+progenindex.size] = np.ones(progenindex.size)*halopropdata[curSnap]['npart'][curIndex]
            noffset += progenindex.size
            nactive += progenindex.size
            if (nactive > maxnsample):
                break
        print('Done processing ', i, 'with', nactive, 'mergers in', time.clock()-start)
    print('Finshed processing all snapshots, have ',noffset,'mergers')
    relraddata=relraddata[:noffset]
    numpartsdata['Merging'][:noffset]
    numpartsdata['Host'][:noffset]

    #now plot
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()
    colorval=['blue','DarkOrange','Crimson']
    cmapval=['Blues','Purples', 'Reds']
    itree=0

    # Setup the dimentons of the plot
    left, width = 0.12, 0.64
    bottom, height = 0.1, 0.7
    left_h = left + width + 0.0
    bottom_h = bottom + height +0.0

    rect_hex = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.12]
    rect_histy = [left_h, bottom, 0.12, height]

    fig=plt.figure( figsize=(10, 8))
    itree=0

    #Setup the axis
    axHex = fig.add_axes(rect_hex)
    axHistx = fig.add_axes(rect_histx)
    axHisty = fig.add_axes(rect_histy)

    #Log10 the data
    for itree in range(len(treelabel)):
        print(itree)
        mergernpart=treedata[treelabel[itree]]['mergerstats'][1]
        ratiomergerradius=treedata[treelabel[itree]]['mergerstats'][0]
        x =np.log10(mergernpart)
        y =np.log10(ratiomergerradius)

        ### Do the hexbin plot ###

        #Setup the hexbin plot
        #hb = axHex.hexbin(x,y,gridsize=20,extent=[1,5.5,-2.0,0.5],cmap=cmapval[itree],mincnt=5.0,vmax=-1,vmin=-5,bins="log",alpha=0.05,edgecolor=colorval[itree],facecolor='white')#,reduce_C_function=np.sum)
        #hb = axHex.hexbin(x,y,gridsize=25,extent=[1,5.5,-2.0,0.5],cmap='binary',mincnt=10.0,vmax=3,vmin=0.5,bins="log",alpha=0.5,edgecolor=colorval[itree],linewidths=2, facecolor='White')#,reduce_C_function=np.sum)
        #instead of generating hexabins plot contours
        #generate 2d histo
        counts,xbins,ybins=np.histogram2d(x,y,bins=[np.arange(1.5,5.1,(5-1.5)/25.0),np.arange(-2.25,0.5,(0.5+2.25)/25.0)])
        axHex.contour(counts.T,extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                      levels=[10.0,1000.0],colors=colorval[itree],linewidths=[2,4])

        #get mean,median,quantiles, etc
        """
        axHex.plot([np.median(x)],[np.median(y)],marker='s',mew=0,ms=5,alpha=0.7,zorder=10,color=colorval[itree])
        xerrval=[np.median(x)-np.percentile(x,[16.])[0],np.percentile(x,[84.])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[16.])[0],np.percentile(y,[84.])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=4,capsize=9)
        xerrval=[np.median(x)-np.percentile(x,[2.5])[0],np.percentile(x,[97.5])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[2.5])[0],np.percentile(y,[97.5])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=2,capsize=4)
        """
        #now get objects slight above the limit and perhaps later times
        wdata=np.where((treedata[treelabel[itree]]['mergerstats'][1]>=40.0)*(treedata[treelabel[itree]]['mergerstats'][2]>=40.0))
        mergernpart=treedata[treelabel[itree]]['mergerstats'][1][wdata]
        ratiomergerradius=treedata[treelabel[itree]]['mergerstats'][0][wdata]
        x =np.log10(mergernpart)
        y =np.log10(ratiomergerradius)
        axHex.plot([np.median(x)],[np.median(y)],marker='o',mew=0,ms=15,alpha=0.7,zorder=10,color=colorval[itree])
        xerrval=[np.median(x)-np.percentile(x,[16.])[0],np.percentile(x,[84.])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[16.])[0],np.percentile(y,[84.])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=4,capsize=9)
        xerrval=[np.median(x)-np.percentile(x,[2.5])[0],np.percentile(x,[97.5])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[2.5])[0],np.percentile(y,[97.5])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=2,capsize=4)
        """
        wdata=np.where((treedata[treelabel[itree]]['mergerstats'][1]>=100.0)*(treedata[treelabel[itree]]['mergerstats'][2]>=100.*treedata[treelabel[itree]]['mergerstats'][1]))
        mergernpart=treedata[treelabel[itree]]['mergerstats'][1][wdata]
        ratiomergerradius=treedata[treelabel[itree]]['mergerstats'][0][wdata]
        x =np.log10(mergernpart)
        y =np.log10(ratiomergerradius)
        axHex.plot([np.median(x)],[np.median(y)],marker='s',mew=0,ms=15,alpha=0.7,zorder=10,color=colorval[itree])
        xerrval=[np.median(x)-np.percentile(x,[16.])[0],np.percentile(x,[84.])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[16.])[0],np.percentile(y,[84.])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=4,capsize=9)
        xerrval=[np.median(x)-np.percentile(x,[2.5])[0],np.percentile(x,[97.5])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[2.5])[0],np.percentile(y,[97.5])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=2,capsize=4)

        wdata=np.where((treedata[treelabel[itree]]['mergerstats'][1]>=100.0)*(treedata[treelabel[itree]]['mergerstats'][2]<=100.*treedata[treelabel[itree]]['mergerstats'][1]))
        mergernpart=treedata[treelabel[itree]]['mergerstats'][1][wdata]
        ratiomergerradius=treedata[treelabel[itree]]['mergerstats'][0][wdata]
        x =np.log10(mergernpart)
        y =np.log10(ratiomergerradius)
        axHex.plot([np.median(x)],[np.median(y)],marker='D',mew=0,ms=15,alpha=0.7,zorder=10,color=colorval[itree])
        xerrval=[np.median(x)-np.percentile(x,[16.])[0],np.percentile(x,[84.])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[16.])[0],np.percentile(y,[84.])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=4,capsize=9)
        xerrval=[np.median(x)-np.percentile(x,[2.5])[0],np.percentile(x,[97.5])[0]-np.median(x)]
        yerrval=[np.median(y)-np.percentile(y,[2.5])[0],np.percentile(y,[97.5])[0]-np.median(y)]
        axHex.errorbar([np.median(x)],[np.median(y)],xerr=[np.array([xerrval[0]]),np.array([xerrval[1]])],yerr=[np.array([yerrval[0]]),np.array([yerrval[1]])],marker='None',alpha=0.7,zorder=9,color=colorval[itree],lw=2,capsize=4)
        """
        #add custom legend
        custom_lines = [
                    Line2D([0], [0], color='blue', lw=4,ls='solid'),
                    Line2D([0], [0], color='DarkOrange', lw=4,ls='solid'),
                    Line2D([0], [0], color='Crimson', lw=4,ls='solid'),
                   ]
        custom_labels=[r'3DFOF.t1', r'6DFOF.SUBS.t1',r'6DFOF.SUBS.t4', ]
        axHex.legend(custom_lines,custom_labels, fontsize=16, loc='lower right',numpoints=2)

        #Line for the npart limit
        axHex.axvline(np.log10(20),color="k",ls="--")
        axHex.axhline(np.log10(1),color="k",ls="--")

        #Set labels and ticks
        axHex.set_xlabel(r'$\log N_p$')
        axHex.set_ylabel(r'$\log r_{\mathrm{merge}}/R_{200\rho_c}$')
        axHex.axis([1,5.5,-2.0,0.5])
        axHex.set_xticks([1,2,3,4,5])
        axHex.set_yticks([-2,-1.5,-1,-0.5,0])
        axHex.set_yticks([-2,-1.5,-1,-0.5,0])
        axHex.xaxis.set_ticks_position('both')
        axHex.yaxis.set_ticks_position('both')


        #reset x,y
        mergernpart=treedata[treelabel[itree]]['mergerstats'][1]
        ratiomergerradius=treedata[treelabel[itree]]['mergerstats'][0]
        x =np.log10(mergernpart)
        y =np.log10(ratiomergerradius)

        #####  Projection of number (top subpanel)


        #Generate the histogram data
        normedCount = np.zeros(100)
        bins = np.linspace(0.5,6,101)
        count, bins = np.histogram(x,bins=bins)

        #Log10 the data
        normedCount[count>0] = np.log10(count[count>0])

        #Set any bin that is 0 to nan so doesn't change the axis
        normedCount[normedCount==0]=np.nan

        #Plot the histogram
        axHistx.plot(bins[:-1],normedCount, color=colorval[itree],lw=4,alpha=0.7)
        #axHistx.fill_between(bins[:-1],np.zeros(len(normedCount)),normedCount,alpha=0.2, color=colorval[itree])

        #Plot the median
        axHistx.axvline(np.median(x),lw=3,ls='-', color=colorval[itree],alpha=0.7)

        #Set the limits and labels
        axHistx.set_ylim(0,6)
        axHistx.set_yticks([0,2.5,5])
        axHistx.set_xlim(axHex.get_xlim())
        axHistx.set_xticks(axHex.get_xticks())
        axHistx.set_ylabel(r'$\log N$')
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHistx.xaxis.set_ticks_position('both')
        axHistx.yaxis.set_ticks_position('both')

        #####  Projection of Radius (right subpanel)

        #Generate the histogram data
        normedCount = np.zeros(100)
        bins = np.linspace(-2.0,0.5,101)
        count, bins = np.histogram(y,bins=bins)

        #Log10
        normedCount[count>0] = np.log10(count[count>0]) #- np.log10(VOL)

        #Set any bin that is 0 to nan so doesn't change the axis
        normedCount[normedCount==0]=np.nan

        #Plot the histogram
        axHisty.plot(normedCount,bins[:-1], color=colorval[itree], lw=4,alpha=0.7)
        #axHisty.fill_betweenx(bins[:-1],np.zeros(len(normedCount)),normedCount,alpha=0.2, color=colorval[itree])

        #Plot the median
        axHisty.axhline(np.median(y),lw=2,ls='--', color=colorval[itree],alpha=0.7)
        axHisty.axhline(np.median(y[np.where(x>=np.log10(40.0))]),lw=3,ls='-', color=colorval[itree],alpha=0.7)

        #Set the limits and labels
        axHisty.set_xlim(0,6)
        axHisty.set_xticks([0,2.5,5])
        axHisty.set_ylim(axHex.get_ylim())
        axHisty.set_yticks(axHex.get_yticks())
        axHisty.set_xlabel(r'$\log N$')
        axHisty.yaxis.set_major_formatter(nullfmt)
        axHisty.xaxis.set_ticks_position('both')
        axHisty.yaxis.set_ticks_position('both')

    # Put in the color bar
    #rect_cb = [left+width+0.11+0.01,bottom,0.02,height]
    #axCB = fig.add_axes(rect_cb)
    #cb = fig.colorbar(hb,cax=axCB)
    #cb.set_label(r"$\log N$",fontsize=14)

    plt.tight_layout()
    plt.savefig(baseplotdir+labelname+'-mergerstats.pdf')
