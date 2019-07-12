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
const Double_t ptbins_debug[]={
  //1., 5., 6., 8., 10., 12., 
  //15., 
  //18., 
  //21., 
  //24., 
  28.,
  32., 
  37.,
  43., //garbage bins
  49., //typical gen pt cut
  56., //typical recopt cut
  64.,
  74.,
  84.,
  97.,
  114.,
  133.,
  153.,
  174.,
  196.,
  220.,
  245.,
  272.,
  300.,
  330.,
  362.,
  395.,
  430.,
  468.,
  507.,
  548.,
  592.,
  638.,
  686.,
  737.,
  790.,
  846.,
  905.,
  967.,//,
  //1000.,
  1032., //generally, garbage bins
  1101., 
  1172., 
  1248., 
  1327.,
  1410.//,
  //1497., 1588., 1684., 1784., 1890., 2000., 2116.

}; 
const int nbins_pt_debug = sizeof(ptbins_debug)/sizeof(Double_t)-1;

// eta bins/boundaries


const Double_t absetabins[]={
  0.0,
  0.5, 
  1.0, 
  1.5,
  2.0, 
  2.5, 
  3.0//, 
  //3.2,
  //4.7
};
std::string absetabins_str[]={
  "0.0", 
  "0.5",  // 0.0, 0.5, 
  "1.0", 
  "1.5",  // 1.0, 1.5, 
  "2.0",
  "2.5",
  "3.0"//, 
  //3.2,
  //4.7
};
//const double absetabins_lo[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.2};
//const double absetabins_hi[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.7};
const int nbins_abseta = 6; //sizeof(absetabins)/sizeof(Double_t)-1;







const int fit_ptlo_bin_00y05=3, fit_pthi_bin_00y05=nbins_pt_debug-7;
const double fgaus_width_arr_00y05[]={
  0.27,// gen pt low edge =43., //garbage bins
  0.26,// gen pt low edge =49., //typical gen pt cut
  0.25,// gen pt low edge =56., //typical recopt cut
  0.24,// gen pt low edge =64.,
  0.23,// gen pt low edge =74.,
  0.22,// gen pt low edge =84.,
  0.21,// gen pt low edge =97.,
  0.21,// gen pt low edge =114.,
  0.20,// gen pt low edge =133.,
  0.18,// gen pt low edge =153.,
  0.18,// gen pt low edge =174.,
  0.18,// gen pt low edge =196.,
  0.17,// gen pt low edge =220.,
  0.17,// gen pt low edge =245.,
  0.17,// gen pt low edge =272.,
  0.15,// gen pt low edge =300.,
  0.15,// gen pt low edge =330.,
  0.15,// gen pt low edge =362.,
  0.13,// gen pt low edge =395.,
  0.13,// gen pt low edge =430.,
  0.13,// gen pt low edge =468.,
  0.12,// gen pt low edge =507.,
  0.12,// gen pt low edge =548.,
  0.12,// gen pt low edge =592.,
  0.10,// gen pt low edge =638.,//fit must reach at least here for this |y| range
  0.10,// gen pt low edge =686.,
  0.10,// gen pt low edge =737.,
  0.10,// gen pt low edge =790.,
  0.10// gen pt low edge =846.,
  //0.10// gen pt low edge =905.,
  
};
const int n_00y05_widths=sizeof(fgaus_width_arr_00y05)/sizeof(double);







const int fit_ptlo_bin_05y10=3, fit_pthi_bin_05y10=nbins_pt_debug-9;
const double fgaus_width_arr_05y10[]={
  0.27,// gen pt low edge =43., //garbage bins
  0.26,// gen pt low edge =49., //typical gen pt cut
  0.25,// gen pt low edge =56., //typical recopt cut
  0.24,// gen pt low edge =64.,
  0.23,// gen pt low edge =74.,
  0.22,// gen pt low edge =84.,
  0.21,// gen pt low edge =97.,
  0.21,// gen pt low edge =114.,
  0.20,// gen pt low edge =133.,
  0.18,// gen pt low edge =153.,
  0.18,// gen pt low edge =174.,
  0.18,// gen pt low edge =196.,
  0.17,// gen pt low edge =220.,
  0.17,// gen pt low edge =245.,
  0.17,// gen pt low edge =272.,
  0.15,// gen pt low edge =300.,
  0.15,// gen pt low edge =330.,
  0.15,// gen pt low edge =362.,
  0.13,// gen pt low edge =395.,
  0.13,// gen pt low edge =430.,
  0.13,// gen pt low edge =468.,
  0.12,// gen pt low edge =507.,
  0.12,// gen pt low edge =548.,
  0.12,// gen pt low edge =592.,
  0.10,// gen pt low edge =638.,//fit must reach at least here for this |y| range
  0.10,// gen pt low edge =686.,
  0.10// gen pt low edge =737.,
  //0.10,// gen pt low edge =790.,
  //0.10 // gen pt low edge =846.,
  //0.10// gen pt low edge =905.,
  
};
const int n_05y10_widths=sizeof(fgaus_width_arr_05y10)/sizeof(double);








const int fit_ptlo_bin_10y15=3, fit_pthi_bin_10y15=nbins_pt_debug-10;
const double fgaus_width_arr_10y15[]={
  0.20,//0.20,// gen pt low edge =43., //garbage bins
  0.20,//0.20,// gen pt low edge =49., //typical gen pt cut
  0.20,//0.20,// gen pt low edge =56., //typical recopt cut
  0.20,//0.20,// gen pt low edge =64.,
  0.20,//0.20,// gen pt low edge =74.,
  0.20,//0.20,// gen pt low edge =84.,
  0.20,//0.20,// gen pt low edge =97.,
  0.20,//0.20,// gen pt low edge =114.,
  0.20,//0.20,// gen pt low edge =133.,
  0.20,//0.20,// gen pt low edge =153.,
  0.20,//0.20,// gen pt low edge =174.,
  0.20,//0.20,// gen pt low edge =196.,
  0.17,//0.20,// gen pt low edge =220.,
  0.17,//0.20,// gen pt low edge =245.,
  0.17,//0.20,// gen pt low edge =272.,
  0.17,//0.20,// gen pt low edge =300.,
  0.16,//0.20,// gen pt low edge =330.,
  0.16,//0.20,// gen pt low edge =362.,
  0.16,//0.20,// gen pt low edge =395.,
  0.15,//0.20,// gen pt low edge =430.,
  0.15,//0.19,// gen pt low edge =468.,
  0.15,//0.18,// gen pt low edge =507.,
  0.14,//0.18,// gen pt low edge =548.,
  0.12,//0.18,// gen pt low edge =592.,//fit must reach at least here for this |y| range
  0.11,//0.14,// gen pt low edge =638.,
  0.11//,0.14 //// gen pt low edge =686.,
  //0.10// gen pt low edge =737.,
  //0.10,// gen pt low edge =790.,
  //0.10,// gen pt low edge =846.,
  //0.10// gen pt low edge =905.,
  
};
const int n_10y15_widths=sizeof(fgaus_width_arr_10y15)/sizeof(double);








const int fit_ptlo_bin_15y20=3, fit_pthi_bin_15y20=nbins_pt_debug-13;
const double fgaus_width_arr_15y20[]={
  0.26,// gen pt low edge =43., //garbage bins
  0.25,// gen pt low edge =49., //typical gen pt cut
  0.24,// gen pt low edge =56., //typical recopt cut
  0.23,// gen pt low edge =64.,
  0.22,// gen pt low edge =74.,
  0.22,// gen pt low edge =84.,
  0.21,// gen pt low edge =97.,
  0.21,// gen pt low edge =114.,
  0.20,// gen pt low edge =133.,
  0.20,// gen pt low edge =153.,
  0.20,// gen pt low edge =174.,
  0.18,//18,// gen pt low edge =196.,
  0.17,//15,// gen pt low edge =220.,
  0.16,//15,// gen pt low edge =245.,
  0.16,//15,// gen pt low edge =272.,
  0.16,//13,// gen pt low edge =300.,
  0.15,//13,// gen pt low edge =330.,
  0.15,//13,// gen pt low edge =362.,
  0.14,//11,// gen pt low edge =395.,
  0.13,//11,// gen pt low edge =430.,
  0.13,//11,// gen pt low edge =468.,
  0.12,//11,// gen pt low edge =507.,//fit must reach at least here for this |y| range
  0.11//11,// gen pt low edge =548.,
  //0.12,// gen pt low edge =592.,
  //0.10// gen pt low edge =638.,
  //0.10,// gen pt low edge =686.,
  //0.10,// gen pt low edge =737.,
  //0.10,// gen pt low edge =790.,
  //0.10,// gen pt low edge =846.,
  //0.10// gen pt low edge =905.,
  
};
const int n_15y20_widths=sizeof(fgaus_width_arr_15y20)/sizeof(double);

double* setWidthArr(int absetabin=0, int *nwidths=NULL){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in setWidthArr"<<std::endl;
  double* widtharr=NULL;
  if(!nwidths)assert(false);
  if(absetabin==0){
    widtharr=(double*) fgaus_width_arr_00y05;
    *nwidths= n_00y05_widths;
  }
  else if(absetabin==1){
    widtharr=(double*) fgaus_width_arr_05y10;
    *nwidths= n_05y10_widths;
  }
  else if(absetabin==2){
    widtharr=(double*) fgaus_width_arr_10y15;
    *nwidths= n_10y15_widths;
  }
  else if(absetabin==3){
    widtharr=(double*) fgaus_width_arr_05y10;
    *nwidths= n_15y20_widths;
  }
  else{
    std::cout<<"ERROR! setWidthArr absetabin not valid. exit."<<std::endl;
    assert(false);
  }
  
  if(funcDebug)std::cout<<"*nwidths="<<*nwidths<<std::endl;
  if(funcDebug)for(int i=0; i<(*nwidths);i++)std::cout<<"widtharr[i="<<i<<"]="<<widtharr[i]<<std::endl;
  if(funcDebug)std::cout<<std::endl<<"exiting setWidthArr"<<std::endl<<std::endl;
  return widtharr;
}





























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
