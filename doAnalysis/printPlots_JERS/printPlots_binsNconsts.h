//#ifndef __boundaries_h_
//#define __boundaries_h_

#include <cstdlib>
#include <string>
#include <TMath.h>

#define NOBJECT_MAX 16384

//#define pi 3.14159265359
//const double pi=(double)TMath::Pi();
//const double pi=3.14159265359;

// string arrays, etc.
const char *cdir[]  = {"PP"};
const char *ccent[] = {"PP"};
const int ncen=1;

const std::string algo ="ak";// strings of jet algo
//const std::string algo[] ={"ak","akPu"};// strings of jet algo
const std::string jetType="PF";// strings of jet type
//const std::string jetType[]={"PF","Calo","Cs"};// strings of jet type

//const std::string srad[] ={"4"};// strings of jet radii
//const int Nrad =sizeof(srad)/sizeof(std::string);//wasn't in original code, 10 arbitrarily chosen

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
//const double pthat[]={ 15, 30, 50, 80, 120, 170, 
//		       220, 280, 370, 460, 540, 2000 };
const double pthat[]={ 15, 30, 50, 80, 120, 170, 
		       220, 280, 370, 9999 };
const int Nbins_pthat=sizeof(pthat)/sizeof(double);

// giving a TH1F an array for binning is fussy; demands root-type numbers or else it gives 
// Error in <TAxis::TAxis::Set>: bins must be in increasing order
const Double_t ptbins_debug[]={
  //15., 18., 21., 24.,
  //  28.,

  32.,
  37.,
  43., //garbage bins
  49. , 56.,
  64., 74., 84., 97., 114.,
  133., 153., 174., 196.,
  220., 245., 272.,
  300., 330., 362., 395., 430., 468., 507.,
  548., 592., 638., 686.,
  737., 790., 846.,
  905.,   
  967.,
  1032.,//  //1000.,
  1101.,
  1172.,
  1248.,
  1327.,
  1410.
  
  // BINS USED FOR APS DNP BELOW
  //  //5.,
  //  //15., 18., 21., 24., 28.,
  //  32., 
  //  37., 
  //  43., 
  //  49., 
  //  56., 
  //  64., 
  //  74., 
  //  84., 97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362., 395., 430., 468.,
  //  507., 548., 592., 638., 686., // junk bins below
  //  1000.,
  //  1500.
  //  //  1050.
  //  //  737, 790, 
  //  //  846, 905, 
  //  //  967
}; 
const int nbins_pt_debug = sizeof(ptbins_debug)/sizeof(Double_t)-1;

// eta bins/boundaries


const Double_t absetabins[]={
  0.0,
  0.5, 
  1.0, 
  1.5,
  2.0//, 
  //2.5, 
  //3.0, 
  //3.2,
  //4.7
};
std::string absetabins_str[]={
  "0.0", 
  "0.5",  // 0.0, 0.5, 
  "1.0", 
  "1.5",  // 1.0, 1.5, 
  "2.0"//,
  //2.5,
  //3.0, 
  //3.2,
  //4.7
};
//const double absetabins_lo[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.2};
//const double absetabins_hi[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.7};
const int nbins_abseta = 8; //sizeof(absetabins)/sizeof(Double_t)-1;


































//const double ptbins_ana[] = {
//  50, 60, 70, 80, 90,
//  100, 110, 130, 150, 170, 190,
//  210, 240, 270, 300};
//const int nbins_ana=sizeof(ptbins_ana)/sizeof(double)-1;


//const double ptbins_jec[]={ 17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 
//			    97, 114, 133, 153, 174, 196, 
//			    220, 245, 272, 300, 350, 400, 550, 790, 1000};
//const int nbins_pt_jec = sizeof(ptbins_jec)/sizeof(double) -1;

//const int ptbins[] = {15 , 30 , 50 , 80 , 120, 170, 220, 300, 500 };
//const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;

// const double boundaries_pt[nbins_pt+1] = {
//  49, 56, 64, 74, 84, 97, 114, 133,
//  153, 174, 196, 220, 245, 272, 300
//};
// const int nbins_pt = sizeof(boundaries_pt)/sizeof(double);-1


//const double xmin=ptbins[0];
//const double xmax=ptbins[nbins_pt];
