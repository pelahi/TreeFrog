/*! \file utilities.cxx
 *  \brief this file contains an assortment of utilities
 */

#include "TreeFrog.h"

/// Time functions
//@{
int GetMilliCount()
{
  // Something like GetTickCount but portable
  // It rolls over every ~ 12.1 days (0x100000/24/60/60)
  // Use GetMilliSpan to correct for rollover
  timeb tb;
  ftime( &tb );
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}

int GetMilliSpan( int nTimeStart )
{
  int nSpan = GetMilliCount() - nTimeStart;
  if ( nSpan < 0 )
    nSpan += 0x100000 * 1000;
  return nSpan;
}

//@}

double MyGetTime(){
#ifdef USEOPENMP
    return omp_get_wtime();
#else
    return (clock() /( (double)CLOCKS_PER_SEC));
#endif
}

//to free up some memory, no need to keep particle ids
void FreeParticleIDMemory(const Int_t numhalos, HaloData *&h){
    if (numhalos == 0) return;
    for (auto j=0;j<numhalos;j++) {
        if (h[j].ParticleID == NULL) continue;
        delete[] h[j].ParticleID;
        h[j].ParticleID=NULL;
    }
}


void FreeProgenitorMemory(Options &opt, ProgenitorData **&pprogen, DescendantDataProgenBased **&pprogendescen)
{
    if (opt.isearchdirection == SEARCHDESCEN) return;
    for (auto i=0;i<opt.numsnapshots;i++) if (pprogen[i]!=NULL) delete[] pprogen[i];
    delete[] pprogen;
    //if used multiple time steps
    if (opt.numsteps>1) {
        for (auto i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
        //check if data is load making sure i is in appropriate range
        if (i>=StartSnap && i<EndSnap) {
#endif
            if (pprogendescen[i]!=NULL) delete[] pprogendescen[i];
#ifdef USEMPI
        }
#endif
        }
        delete[] pprogendescen;
    }
}

void FreeDescendantMemory(Options &opt, DescendantData **&pdescen, ProgenitorDataDescenBased **&pdescenprogen)
{
    if (opt.isearchdirection == SEARCHPROGEN) return;
    for (auto i=0;i<opt.numsnapshots;i++) if (pdescen[i]!=NULL) delete[] pdescen[i];delete[] pdescen;
    //if used multiple time steps
    if (opt.numsteps>1) {
        for (auto i=0;i<opt.numsnapshots;i++) {
#ifdef USEMPI
        //check if data is load making sure i is in appropriate range
        if (i>=StartSnap && i<EndSnap) {
#endif
            if (pdescenprogen[i]!=NULL) delete[] pdescenprogen[i];
#ifdef USEMPI
        }
#endif
        }
    }
    delete[] pdescenprogen;
}
