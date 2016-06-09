// Originally written by: Kurt Jung, December 2015
// Overhaul, Ian Laflotte, May 1st 2016
// Macro for deriving the dijet relative and absolute response from pp data/MC forests
// Uses JME-13-004 as a reference - "A" and "B" formulae are defined there


// C includes
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <vector>
// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1F.h>

//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1+1;//always one for ./ itself, then one more for derive or sum

//// deriveResponse
const int defStartFile=0;
const int defEndFile=2; //exclusive boundary, range of files run over doesnt include endfile
const std::string defInFilelist = "../filelists/5p02TeV_HighPtJet80_9Files_debug_forests.txt";
const std::string defDeriveOutputName = "deriveResponse_defOut.root";
const int defRadius=4;
const bool defIsMCderive=false;
const bool defDebugMode = true;
int deriveResponse(int startfile = defStartFile, int endfile = defEndFile,
                   std::string infile_Forest=defInFilelist, std::string outfile = defDeriveOutputName,
                   int radius = defRadius, bool isMC = defIsMCderive, bool debugMode = defDebugMode );
const int deriveRespArgCount=7+minArgs;

//// sumDijetResponse
const std::string defSumInput="deriveResponse_defOut.root";
const std::string defSumOutput="sumResponse_defOut.root";
const bool defIsMCsum=false;
const bool defDoDraw=true;
const bool defDebugModeSum=true;
int sumResponse(std::string filename=defSumInput, std::string outFileName=defSumOutput,
		bool isMC=defIsMCsum , bool doDraw=defDoDraw ,  bool debugMode=defDebugModeSum);
const int sumRespArgCount=5+minArgs;


//// CONSTANTS
// ---------------------------------------------------------------------------------------------------------------
// useful, code/function-wide constants

const bool atMIT  = true;
const std::string hadoopDir    = "/mnt/hadoop/cms";
const std::string xrootdDirect = "root://cmsxrootd.fnal.gov/";

const int colors[] = {1,2,4,8,20};
const int N_colors= sizeof(colors)/sizeof(int);

const double xbins_pt[] = {40, 60, 80, 110, 200, 1000};
const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;

const double xbins_eta[] = {
  -5.191, -4.716, -4.363, -4.013, -3.664, -3.314, 
  -2.964, -2.650, -2.322, -2.043, 
  -1.830, -1.653, -1.479, -1.305, -1.131, 
  -0.957, -0.783, -0.609, -0.435, -0.261, -0.087,
  0.000,
  0.087,  0.261,  0.435,  0.609,  0.783,  0.957, 
  1.131,  1.305,  1.479,  1.653,  1.830,
  2.043,  2.322,  2.650,  2.964,  
  3.314,   3.664, 4.013,  4.363,  4.716,  5.191};
const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;


//// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// mini functions not called by the frontend/main


// compares given value to boundaries of bins, returns bin number
int findBin(float value, int nbins, const double *xbins){
  for(int i=0; i<nbins; i++) if(value<xbins[i+1]) return i;
  return nbins-1;
}


// missEt = - Sum(Et you can see), i need a better description
TLorentzVector findMissEt(int nPFpart_,  								    
			  std::vector<int> *pfId_, std::vector<float> *pfPt_, 				    
			  std::vector<float> *pfEta_, std::vector<float> *pfPhi_, 			    
			  int nref, float *pt_F,  float *eta_F, float *phi_F, 			     
			  float *m_F, float *rawpt_F, float *eSum,  float *phoSum){                           

  bool missEtDebugMode=false;
  if(missEtDebugMode)std::cout<<std::endl<<"missEtDebugMode is ON"<<std::endl;

  //missing Et from particles  
  TLorentzVector missEt(0,0,0,0);
  for(int i=0; i<nPFpart_; i++){//PF Loop
    if(missEtDebugMode&&i==0)std::cout<<std::endl<<"in findMissEt's PF loop: "<<std::endl;
    TLorentzVector pfp;
    double pfMass=0;
    if(pfId_->at(i)==1 || pfId_->at(i)==4) pfMass = 0.1395702;//pi+/- meson, is this always justified?
    //pfp.SetPtEtaPhiM(pfPt[i],pfEta[i],pfPhi[i],pfMass);
    pfp.SetPtEtaPhiM(pfPt_->at(i),pfEta_->at(i),pfPhi_->at(i),pfMass);
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfPt   (i="<<i<<")  =  "<<pfPt_->at(i)<<std::endl;
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfEta  (i="<<i<<")  =  "<<pfEta_->at(i)<<std::endl;
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfPhi  (i="<<i<<")  =  "<<pfPhi_->at(i)<<std::endl;
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfMass (i="<<i<<")  =  "<<pfMass<<std::endl;
    pfp.SetPz(0.); //2d projection in x-y plane
    missEt += pfp;
  }//end PF loop
  
  //missing Et from jets
  for(int i=0; i<nref; i++){//begin jet loop
    if(missEtDebugMode&&i==0)std::cout<<std::endl<<"in findMissEt's jet loop"<<std::endl;
    if(pt_F[i]<=15) continue;
    TLorentzVector jtTmp(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
    double jetEnergy=jtTmp.E();
    if(missEtDebugMode&&(i%10==0))std::cout<<"eSum[i="<<i<<"] = "<<eSum[i]<<std::endl;
    if( (eSum[i]+phoSum[i])/jetEnergy<0.9 ){
      TLorentzVector jt;
      jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
      jt.SetPz(0.); //project to x-y plane
      missEt += jt;
    }
  }//end jet loop
  
  missEt*=-1; // rotate 180 deg, to get missing Et,
  missEt.SetE(0.);  // only care about Et, not E
  return missEt;

  // old, consider removing
  ////missing Et from jets
  //for(int i=0; i<nref; i++){//begin jet loop
  //  if(pf_F[i]<=15) continue;//quick cut before expensive computations
  //  TLorentzVector jt(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
  //  if( (eSum[i]+phoSum[i])/jt.E()<0.9 ){
  //    jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
  //    jt.SetPz(0.); //project to x-y plane
  //    missEt += jt;
  //  }
  //}//end jet loop
  

}

// missEt = - Sum(Et you can see), i need a better description
// this version uses array as input for PF particles
TLorentzVector findMissEt(int nPFpart_,  								    
			  int *pfId_, float *pfPt_, 				    
			  float *pfEta_, float *pfPhi_, 			    
			  int nref, float *pt_F,  float *eta_F, float *phi_F, 			     
			  float *m_F, float *rawpt_F, float *eSum,  float *phoSum){                           

  bool missEtDebugMode=false;
  if(missEtDebugMode)std::cout<<std::endl<<"missEtDebugMode is ON"<<std::endl;

  //missing Et from particles  
  TLorentzVector missEt(0,0,0,0);
  for(int i=0; i<nPFpart_; i++){//PF Loop
    if(missEtDebugMode&&i==0)std::cout<<std::endl<<"in findMissEt's PF loop: "<<std::endl;
    TLorentzVector pfp;
    double pfMass=0;
    if(pfId_[i]==1 || pfId_[i]==4) pfMass = 0.1395702;//pi+/- meson, is this always justified?
    //pfp.SetPtEtaPhiM(pfPt[i],pfEta[i],pfPhi[i],pfMass);
    pfp.SetPtEtaPhiM(pfPt_[i],pfEta_[i],pfPhi_[i],pfMass);
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfPt  [i="<<i<<"]  =  "<<pfPt_[i]<<std::endl;
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfEta [i="<<i<<"]  =  "<<pfEta_[i]<<std::endl;
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfPhi [i="<<i<<"]  =  "<<pfPhi_[i]<<std::endl;
    if(missEtDebugMode &&(i%100==0)) std::cout<<"pfMass[i="<<i<<"]  =  "<<pfMass<<std::endl;
    pfp.SetPz(0.); //2d projection in x-y plane
    missEt += pfp;
  }//end PF loop
  
  //missing Et from jets
  for(int i=0; i<nref; i++){//begin jet loop
    if(missEtDebugMode&&i==0)std::cout<<std::endl<<"in findMissEt's jet loop"<<std::endl;
    if(pt_F[i]<=15) continue;
    TLorentzVector jtTmp(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
    double jetEnergy=jtTmp.E();
    if(missEtDebugMode&&(i%10==0))std::cout<<"eSum[i="<<i<<"] = "<<eSum[i]<<std::endl;
    if( (eSum[i]+phoSum[i])/jetEnergy<0.9 ){
      TLorentzVector jt;
      jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
      jt.SetPz(0.); //project to x-y plane
      missEt += jt;
    }
  }//end jet loop
  
  missEt*=-1; // rotate 180 deg, to get missing Et,
  missEt.SetE(0.);  // only care about Et, not E
  return missEt;

  // old, consider removing
  ////missing Et from jets
  //for(int i=0; i<nref; i++){//begin jet loop
  //  if(pf_F[i]<=15) continue;//quick cut before expensive computations
  //  TLorentzVector jt(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
  //  if( (eSum[i]+phoSum[i])/jt.E()<0.9 ){
  //    jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
  //    jt.SetPz(0.); //project to x-y plane
  //    missEt += jt;
  //  }
  //}//end jet loop
  

}
