//#ifndef __boundaries_h_
//#define __boundaries_h_

#include <cstdlib>
#include <string>
#include <TMath.h>

#define NOBJECT_MAX 16384

//#define pi 3.14159265359
//const double pi=(double)TMath::Pi();
const double pi=3.14159265359;

// string arrays, etc.
const char *cdir[]  = {"PP"};
const char *ccent[] = {"PP"};
const int ncen=1;

const std::string algo ="ak";// strings of jet algo
//const std::string algo[] ={"ak","akPu"};// strings of jet algo
const std::string jetType="PF";// strings of jet type
//const std::string jetType[]={"PF","Calo","Cs"};// strings of jet type

const std::string srad[] ={"4"};// strings of jet radii
const int Nrad =sizeof(srad)/sizeof(std::string);//wasn't in original code, 10 arbitrarily chosen

const std::string PFCandType[] = {"unknown",
				  "chargedHadron",
				  "electron",  "muon", "photon",
				  "neutralHadron",
				  "HadEnergyinHF", "EMEnergyinHF"};  
const int PFType = sizeof(PFCandType)/sizeof(std::string);


// weights, entries, xsec
const int nentries_file[] = { 0, 333206, 250567, 395126, 368126, 366982, 392206, 181018, 50455};
const int N_nentries_file=sizeof(nentries_file)/sizeof(int);

const double xsecs[] = {5.269E-01, 3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05 , 1.223E-05, 
			3.031E-06 , 7.746E-07, 1.410E-07 , 3.216E-08, 1.001E-08 , 0.0};
const int N_xsecs=sizeof(xsecs)/sizeof(double);

const double weight_xsec[] = { 7.20357e-07, 4.51655e-08, 2.6964e-09, 2.77274e-10, 3.1878e-11, 
			       3.87126e-12, 1.62138e-12, 1.09471e-12, 4.40012e-13};
const int N_weight_xsec=sizeof(weight_xsec)/sizeof(double);


// momentum bins/boundaries
const double pthat[]={ 15, 30, 50, 80, 120, 170, 
		       220, 280, 370, 460, 540, 2000 };
const int Nbins_pthat=sizeof(pthat)/sizeof(double);

const double ptbins_jec[]={ 17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 
			    97, 114, 133, 153, 174, 196, 
			    220, 245, 272, 300, 350, 400, 550, 790, 1000};
const int nbins_pt_jec = sizeof(ptbins_jec)/sizeof(double) -1;

const int ptbins[] = {15 , 30 , 50 , 80 , 120, 170, 220, 300, 500 };
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;
// giving a TH1F an array for binning is fussy; demands root-type numbers or else it gives 
//Error in <TAxis::TAxis::Set>: bins must be in increasing order
const Double_t ptbins_debug[] = {15., 30., 50., 80., 120., 170., 220., 300., 500.};
const Int_t nbins_pt_debug = nbins_pt;
const double xmin=ptbins[0];
const double xmax=ptbins[nbins_pt];

const double ptbins_ana[] = {
  50, 60, 70, 80, 90,
  100, 110, 130, 150, 170, 190,
  210, 240, 270, 300};
const int nbins_ana=sizeof(ptbins_ana)/sizeof(double)-1;

// const double boundaries_pt[nbins_pt+1] = {
//  49, 56, 64, 74, 84, 97, 114, 133,
//  153, 174, 196, 220, 245, 272, 300
//};
// const int nbins_pt = sizeof(boundaries_pt)/sizeof(double);-1


// eta bins/boundaries
const double etabins[] = {
  -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, 
  -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
  -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
  -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
  -0.435, -0.348, -0.261, -0.174, -0.087, 
  +0.000, 
  +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
  +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
  +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, 
  +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, 
  +4.363, +4.538, +4.716, +4.889, +5.191 };
const int nbins_eta = sizeof(etabins)/sizeof(double) -1;

//#endif

