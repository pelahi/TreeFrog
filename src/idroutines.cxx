/*! \file idroutines.cxx
 *  \brief this file contains routines for manipulating/updating ids, particles and halos
 */

#include "TreeFrog.h"


/// \name if halo ids need to be adjusted
//@{
void UpdateHaloIDs(Options &opt, HaloTreeData *&pht) {
    Int_t i,j;
    if (opt.haloidval>0) {
        for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
        if (i>=StartSnap && i<EndSnap) {
#endif
            if(pht[i].numhalos>0) for (j=0;j<pht[i].numhalos;j++) pht[i].Halo[j].haloID+=opt.haloidval*(i+opt.snapshotvaloffset)+opt.haloidoffset;
#ifdef USEMPI
        }
#endif
        }
    }
}
//@}

/// \name if particle ids need to be mapped to indices
//@{

///see if map already exists and read it. otherwise generate it and alter data
void MemoryEfficientMap(Options &opt,HaloTreeData *&pht)
{
#ifndef USEMPI
    int ThisTask=0;
#endif
    map<IDTYPE, IDTYPE> idmap;
    //try reading information and if it does not suceed then produce map
    if (ReadPIDStoIndexMap(opt,idmap)==0) {
        if (ThisTask==0) cout<<"Generating unique memory efficent mapping for particle IDS to index"<<endl;
        idmap=ConstructMemoryEfficientPIDStoIndexMap(opt, pht);
        SavePIDStoIndexMap(opt,idmap);
    }
    MapPIDStoIndex(opt,pht, idmap);
    idmap.clear();
}

///builds a map by identifying the ids of particles in structure across snapshots
///which is memory efficient as only needs array size of maximum number of particles in structures
map<IDTYPE, IDTYPE> ConstructMemoryEfficientPIDStoIndexMap(Options &opt, HaloTreeData *&pht) {
    Int_t offsetsnap;
    double time1,time2;
#ifndef USEMPI
    int ThisTask=0,NProcs=1,NSnap=opt.numsnapshots,StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    //to handle overlapping snapshots between mpi domains
    if (ThisTask==0) offsetsnap=0;
    else offsetsnap=opt.numsteps;

    if (ThisTask==0) cout<<"Mapping PIDS to index "<<endl;
    //place ids in a set so have unique ordered set of ids
    vector<IDTYPE> idvec;
    vector<IDTYPE> idtempvec;
    map<IDTYPE, IDTYPE> idmap;
    unordered_set<IDTYPE> idset;
    time1=MyGetTime();
    time2=MyGetTime();
    Int_t reservesize=0;
    for (auto j=0;j<pht[EndSnap-1].numhalos;j++) reservesize+=pht[EndSnap-1].Halo[j].NumberofParticles;
    if (opt.iexclusiveids) {
        idvec.reserve(reservesize);
        cout<<ThisTask<<" initial reserve of memory for construction of id map is "<<reservesize*sizeof(IDTYPE)/1024./1024./1024.<<"GB"<<endl;
    }
    else {
        idset.reserve(reservesize);
        idvec.reserve(reservesize);
        cout<<ThisTask<<" Non exclusive ids, must generate memory heavy set, initial reserve of memory for construction of id map is "<<reservesize*(sizeof(IDTYPE)+32)/1024./1024./1024.<<"GB"<<endl;
    }
    for (auto i=EndSnap-1;i>=StartSnap+offsetsnap;i--) {
        idtempvec.clear();
        //initialize vector assuming uniqueness of particle ids! Otherwise need to use set's and these are VERY memory heavy relative to the vector
        //otherwise, need to insert into a set for the first round and then generate a vector from this
        //set that ensures uniqueness, which is then sorted
        if (i==EndSnap-1) {
            if (opt.iexclusiveids) {
                for (auto j=0;j<pht[i].numhalos;j++) {
                    for (auto k=0;k<pht[i].Halo[j].NumberofParticles;k++) idvec.push_back(pht[i].Halo[j].ParticleID[k]);
                }
                sort(idvec.begin(),idvec.end());
            }
            else {
                for (auto j=0;j<pht[i].numhalos;j++) {
                    for (auto k=0;k<pht[i].Halo[j].NumberofParticles;k++) idset.insert(pht[i].Halo[j].ParticleID[k]);
                }
            }
        }
        //once done and have sorted vector, generate new vector containing any additional new IDs
        else {
            if (opt.iexclusiveids) {
                for (auto j=0;j<pht[i].numhalos;j++) {
                    for (auto k=0;k<pht[i].Halo[j].NumberofParticles;k++) {
                        if (!(binary_search(idvec.begin(), idvec.end(), pht[i].Halo[j].ParticleID[k]))) {
                            idtempvec.push_back(pht[i].Halo[j].ParticleID[k]);
                        }
                    }
                }
                //append and then resort
                idvec.insert(idvec.end(), idtempvec.begin(), idtempvec.end());
                sort(idvec.begin(),idvec.end());
            }
            else {
                for (auto j=0;j<pht[i].numhalos;j++) {
                    for (auto k=0;k<pht[i].Halo[j].NumberofParticles;k++) idset.insert(pht[i].Halo[j].ParticleID[k]);
                }
            }
        }
        cout<<"Finished inserting "<<i<<endl;
    }
    if (opt.iexclusiveids==0) {
        for (auto setval: idset) idvec.push_back(setval);
        sort(idvec.begin(),idvec.end());
        idset.clear();
    }
    time1=MyGetTime()-time1;
    if (opt.iverbose) cout<<ThisTask<<" finished getting unique ids of "<<idvec.size()<<" in "<<time1<<endl;

#ifdef USEMPI
    idmap = MPIGatherIDs(opt, idvec);
#else
    opt.MaxIDValue = idvec.size();
    time1=MyGetTime();
    for (auto i=0; i<opt.MaxIDValue;i++) idmap.insert(pair<IDTYPE, IDTYPE>(idvec[i],(IDTYPE)i));
    time1=MyGetTime()-time1;
    if (opt.iverbose) cout<<ThisTask<<" constructed map in "<<time1<<endl;
#endif
    idvec.clear();

    time2=MyGetTime()-time2;
    if (opt.iverbose) {
        cout<<ThisTask<<" finished getting GLOBAL unique ids of "<<opt.MaxIDValue<<" in "<<time2<<endl;
    }
    return idmap;
}


void MapPIDStoIndex(Options &opt, HaloTreeData *&pht, map<IDTYPE, IDTYPE> &idmap) {
#ifndef USEMPI
        int ThisTask=0,StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    Int_t i,j,k;
    if (ThisTask==0) cout<<"Mapping PIDS to index "<<endl;
    for (i=StartSnap;i<EndSnap;i++) {
        #ifdef USEOPENMP
        #pragma omp parallel default(shared) \
        private(j,k)
        {
        #pragma omp for schedule(dynamic) nowait
        #endif
        for (j=0;j<pht[i].numhalos;j++) {
            for (k=0;k<pht[i].Halo[j].NumberofParticles;k++){
                pht[i].Halo[j].ParticleID[k]=idmap[pht[i].Halo[j].ParticleID[k]];
            }
        }
        #ifdef USEOPENMP
        }
        #endif
    }
}

void MapPIDStoIndex(Options &opt, HaloTreeData *&pht) {
#ifndef USEMPI
    int ThisTask=0,StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    Int_t i,j,k;
    if (ThisTask==0) cout<<"Mapping PIDS to index "<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=StartSnap;i<EndSnap;i++) {
        for (j=0;j<pht[i].numhalos;j++)
            for (k=0;k<pht[i].Halo[j].NumberofParticles;k++)
                opt.mappingfunc(pht[i].Halo[j].ParticleID[k]);
    }
#ifdef USEOPENMP
}
#endif
}

void simplemap(IDTYPE &i) {}

//@}

///check to see if ID data compatible with accessing index array allocated
void IDcheck(Options &opt, HaloTreeData *&pht){
    int ierrorflag=0,ierrorsumflag=0;
#ifndef USEMPI
    int ThisTask=0,NProcs=1,NSnap=opt.numsnapshots,StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    for (int i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        ierrorflag=0;
        for (Int_t j=0;j<pht[i].numhalos;j++) {
            for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++)
                if (pht[i].Halo[j].ParticleID[k]<0||pht[i].Halo[j].ParticleID[k]>opt.MaxIDValue) {
                    cout<<ThisTask<<" snapshot "<<i<<" particle id out of range "<<pht[i].Halo[j].ParticleID[k]<<" not in [0,"<<opt.MaxIDValue<<")"<<endl;
                    ierrorflag=1;
                    break;
                }
            if (ierrorflag==1) {ierrorsumflag+=1;break;}
        }
#ifdef USEMPI
    }
#endif
    }
#ifdef USEMPI
    int mpierrorflag;
    MPI_Allreduce(&ierrorsumflag,&mpierrorflag,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    ierrorsumflag=mpierrorflag;
#endif
    if (ierrorsumflag>0) {
#ifdef USEMPI
        if (ThisTask==0)
#endif
        cout<<"Error in particle ids, outside of range. Exiting"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
}


//then allocate simple array used for accessing halo ids of particles through their IDs
void AllocatePFOFMem(Options &opt, PFOFTYPE *&pfof, PFOFTYPE *&prank) {
    pfof=new PFOFTYPE[opt.MaxIDValue];
    for (auto i=0;i<opt.MaxIDValue;i++) pfof[i]=0;
    if (opt.imerittype==MERITRankWeightedBoth) {
        prank=new PFOFTYPE[opt.MaxIDValue];
        for (auto i=0;i<opt.MaxIDValue;i++) prank[i]=0;
    }
}

///set the pfof array based on halo ids
void SetPFOF(Options &opt, PFOFTYPE *&pfof, const Int_t nhalos, HaloData *&halos) {
    for (auto j=0;j<nhalos;j++) {
        for (auto k=0;k<halos[j].NumberofParticles;k++) {
            pfof[halos[j].ParticleID[k]]=j+1;
        }
    }
    //now if also doing core weighting then update the halo id associated with the particle so that
    //it is its current halo core ID + total number of halos
    if (opt.icorematchtype!=PARTLISTNOCORE && opt.particle_frac<1 && opt.particle_frac>0) {
        for (auto j=0;j<nhalos;j++) {
            unsigned long long newnp=max((unsigned long long)(halos[j].NumberofParticles*opt.particle_frac), (unsigned long long)opt.min_numpart);
            newnp=min((unsigned long long)halos[j].NumberofParticles, newnp);
            for (auto k=0;k<newnp;k++)
                pfof[halos[j].ParticleID[k]]=j+1+nhalos;
        }
    }
}

void SetPRank(Options &opt, PFOFTYPE *&prank, const Int_t nhalos, HaloData *&halos) {
    if (opt.imerittype==MERITRankWeightedBoth) {
        for (auto j=0;j<nhalos;j++) {
            for (auto k=0;k<halos[j].NumberofParticles;k++) {
                prank[halos[j].ParticleID[k]]=k+1;
            }
        }
    }
}

///reset the index halo array
void ResetPFOF(Options &opt, PFOFTYPE *&pfof, const Int_t nhalos, HaloData *&halos)
{
    for (auto j=0;j<nhalos;j++) {
        for (auto k=0;k<halos[j].NumberofParticles;k++) {
            pfof[halos[j].ParticleID[k]]=0;
        }
    }
}
