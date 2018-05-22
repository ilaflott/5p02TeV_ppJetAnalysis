#include <cstdlib>
#include <string>
#include <TMath.h>

//#define NOBJECT_MAX 16384

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


//// momentum bins/boundaries
//const double pthat[]={ 15, 30, 50, 80, 120, 170, 
//		       220, 280, 370, 460, 540, 2000 };
////const double pthat[]={ 15, 30, 50, 80, 120, 170, 
////		       220, 280, 370, 9999 };
//const int Nbins_pthat=sizeof(pthat)/sizeof(double);

// giving a TH1F an array for binning is fussy; demands root-type numbers or else it gives 
// Error in <TAxis::TAxis::Set>: bins must be in increasing order
const Double_t ptbins_debug[]={
  
  //32.,
  //37.,  //40., //this junk bin added in for cosmetics
  43., 
  49., //junk bins from here an above
  56.,   64.,   74.,   84., 97., 114., 
  133., 153., 174., 196., 220., 
  245., 272., 300., 
  330., 362., 395., 430., 
  468., 507., 
  548., 592., 638., 
  686., 
  //737., 
  //790., 
  //846., 
  //905., 
  //967 ., // end of SMP Official bin limits 
  1000.,  // to catch overflow or other entries?
  1050. // junk bin, cosmetics only
}; 
const int nbins_pt_debug = sizeof(ptbins_debug)/sizeof(Double_t)-1;


// eta bins/boundaries
const Double_t absetabins[]={
  0.0, 0.5,  
  1.0, 1.5,  
  2.0, 2.5,  
  3.0, 3.2,  
  4.7
};
const int nbins_abseta = sizeof(absetabins)/sizeof(Double_t)-1;











//    5.,
//    15., 18., 21., 24., 28.,
//  32., 37., 43., 49., 56.,
//  64., 74., 84., 97., 114.,
//  133., 153., 174., 196.,
//  220., 245., 272., 300., 
//  330., 362., 395., 430.,
//  468., 507., 548., 592.,
//  638., 686., 1000.
// my tweaked version

//  7., 9., 12., /*3., 4., 5., */
//  15., 18., 21., 24., 28.,
//  32., 37., 43., 49., 56.,
//  64., 74., 84., 97., 114.,
//  133., 153., 174., 196.,
//  220., 245., 272., 300., 
//  330., 362., 395., 430.,
//  468., 507., 548., 592.,
//  638., 686., 1000.//, 1500
//}; //raghavs original sugg. JER genpt binning
