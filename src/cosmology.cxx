/*! \file stfio.cxx
 *  \brief this file contains routines related to cosmology
 */

#include "TreeFrog.h"


///\name Simple cosmology related functions
//@{
void CalcOmegak(Options &opt) {
    opt.Omega_k=(1-opt.Omega_m-opt.Omega_Lambda-opt.Omega_r-opt.Omega_nu-opt.Omega_de);
}
void CalcCriticalDensity(Options &opt, Double_t a){
    Double_t Hubble=GetHubble(opt,a);
    opt.rhocrit=3.*Hubble*Hubble/(8.0*M_PI*opt.Gravity);
}
void CalcBackgroundDensity(Options &opt, Double_t a){
    Double_t Hubble=GetHubble(opt,1.0);
    opt.rhobg=3.*Hubble*Hubble/(8.0*M_PI*opt.Gravity)*opt.Omega_m/(a*a*a);
}
void CalcCosmoParams(Options &opt, Double_t a){
    CalcOmegak(opt);
    CalcCriticalDensity(opt,a);
    CalcBackgroundDensity(opt,a);
}

Double_t GetHubble(Options &opt, Double_t a){
    return opt.hval*opt.H*sqrt(opt.Omega_k*pow(a,-2.0)+opt.Omega_m*pow(a,-3.0)+opt.Omega_r*pow(a,-4.0)+opt.Omega_Lambda+opt.Omega_de*pow(a,-3.0*(1+opt.w_de)));
}

double GetInvaH(double a, void * params) {
    double Omega_m = ((double*)params)[0];
    double Omega_Lambda = ((double*)params)[1];
    double Omega_r = ((double*)params)[2];
    double Omega_nu = ((double*)params)[3];
    double Omega_k = ((double*)params)[4];
    double Omega_de = ((double*)params)[5];
    double w_de = ((double*)params)[6];

    double H=sqrt(Omega_k*pow(a,-2.0)+Omega_m*pow(a,-3.0)+Omega_r*pow(a,-3.0)+Omega_Lambda+Omega_de*pow(a,-3.0*(1+w_de)));
    return 1.0/(a*H);
}
//return cosmic time in years
Double_t CalcCosmicTime(Options &opt, Double_t a1, Double_t a2){
    Double_t cosmictime;
    double result, error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    double params[10];
    params[0]=opt.Omega_m;
    params[1]=opt.Omega_Lambda;
    params[2]=opt.Omega_k;
    params[3]=opt.Omega_r;
    params[4]=opt.Omega_nu;
    params[5]=opt.Omega_de;
    params[6]=opt.w_de;
    F.function = &GetInvaH;
    F.params = (void*)params;
    gsl_integration_qags (&F, a1, a2, 0, 1e-7, 1000, w, &result, &error);
    gsl_integration_workspace_free (w);
    cosmictime = 1./(opt.hval*opt.H*1.02269032e-9)*result;
    return cosmictime;
}

Double_t CalcFreeFallTimeFromOverdensity(Options &opt, Double_t a){
    Double_t tff;
    Double_t H = GetHubble(opt, a);
    tff = M_PI/(2.0*sqrt(opt.deltarho)*H);
    tff *= opt.HubbletoGyrs;
    return tff;
}

void FillNumStepsArray(Options &opt){
#ifndef USEMPI
    int ThisTask = 0;
#endif
    if (opt.delta_time > 0) FillNumStepsArrayBasedOnTime(opt);
    else if (opt.delta_scalefactor > 0) FillNumStepsArrayBasedOnScaleFactor(opt);
    else if (opt.delta_dynamical_time_fraction > 0) FillNumStepsArrayBasedOnDynamicalTime(opt);
    if (opt.iverbose>=1 && ThisTask == 0) {
        cout<<" Number of steps used at a given snapshot are :"<<endl;
        for (auto i=0; i<opt.numsnapshots;i++) {
            cout<<i<<" "<<opt.numstepsarray[i]<<endl;
        }
    }
}

void FillNumStepsArrayBasedOnTime(Options &opt) {
    opt.numstepsarray[opt.numsnapshots-1] = 0;
    for (auto i=0;i<opt.numsnapshots-1;i++){
        opt.numstepsarray[i] = -1;
        for (auto j=i+1;j<opt.numsnapshots;j++) {
            if (opt.snapshot_time[i]+opt.delta_time<opt.snapshot_time[j]){
                opt.numstepsarray[i] = j-i;
                break;
            }
        }
        if (opt.numstepsarray[i] == -1) {
           opt.numstepsarray[i] = opt.numsnapshots-1-i;
        }
    }
    opt.numsteps = 0;
    for (auto &x:opt.numstepsarray) if (opt.numsteps < x) opt.numsteps = x;
}

void FillNumStepsArrayBasedOnScaleFactor(Options &opt) {
    opt.numstepsarray[opt.numsnapshots-1] = 0;
    for (auto i=0;i<opt.numsnapshots-1;i++){
        opt.numstepsarray[i] = -1;
        for (auto j=i+1;j<opt.numsnapshots;j++) {
            if (opt.snapshot_scalefactor[i]+opt.delta_scalefactor<opt.snapshot_scalefactor[j]){
                opt.numstepsarray[i] = j-i;
                break;
            }
        }
        if (opt.numstepsarray[i] == -1) {
           opt.numstepsarray[i] = opt.numsnapshots-1-i;
        }
    }
    opt.numsteps = 0;
    for (auto &x:opt.numstepsarray) if (opt.numsteps < x) opt.numsteps = x;
}

void FillNumStepsArrayBasedOnDynamicalTime(Options &opt) {
    vector<Double_t> tend(opt.numsnapshots,0);
    for (auto i=0;i<opt.numsnapshots;i++) {
        tend[i] = CalcFreeFallTimeFromOverdensity(opt, opt.snapshot_scalefactor[i])*1e9*opt.delta_dynamical_time_fraction;
    }
    opt.numstepsarray[opt.numsnapshots-1] = 0;
    for (auto i=0;i<opt.numsnapshots-1;i++){
        opt.numstepsarray[i] = -1;
        for (auto j=i+1;j<opt.numsnapshots;j++) {
            if (opt.snapshot_time[i]+tend[i]<opt.snapshot_time[j]){
                opt.numstepsarray[i] = j-i;
                break;
            }
        }
        if (opt.numstepsarray[i] == -1) {
           opt.numstepsarray[i] = opt.numsnapshots-1-i;
        }
    }
    opt.numsteps = 0;
    for (auto &x:opt.numstepsarray) if (opt.numsteps < x) opt.numsteps = x;
}

//@}
