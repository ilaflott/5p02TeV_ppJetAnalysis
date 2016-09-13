// Author; Ian Laflotte
// Inspired by code inhereted fro,, Raghav K., Kurt J.
// 
// reads and writes hists from pp data/MC readFiles output to pdf files
// produces jet, trigger, and evt QA spectra
// compile with...
// g++ printQAPlots.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests_ppData.exe

////////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
// ROOTSYS
#include <TSystem.h>
#include <TProfile.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TEventList.h>
#include <TCut.h>
// utilities
#include <TMath.h>
#include <TRandom3.h>
#include <TStopwatch.h>
// plotting
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
// histos
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

// custom
#include "ppMC_etc.h"
#include "ppData_etc.h"

//CMSSW_BASE
const std::string CMSSW_BASE = 
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string condor_dir = "saved_outputCondor/";
//ppMC/ppData input directories, filenames, and fullFilenames

const std::string input_ppData_dir =
  "src/readFiles/readForests_ppData/";//+condor_dir;
const std::string input_ppMC_dir =
  "src/readFiles/readForests_ppMC/";//+condor_dir;

const std::string outputDir=CMSSW_BASE+"src/printQAPlots/output/";

// other useful strings
const int radius_int=4;
const std::string radius_int_str=std::to_string(radius_int);

const std::string jetType="PF";
const std::string etaWidth = "20_eta_20";

// combos
const std::string radius="R"+radius_int_str+"_";
const std::string pfRad_etaWidth=radius+etaWidth;
const std::string jetAna="ak"+radius+jetType+"JetAnalyzer";

// variable names, string array
const std::string var[] = {   
  "jtpt" ,  "rawpt",  
  "jteta", "jtphi", //normalize to data
  "trkMax", "trkSum", "trkHardSum", 
  "chMax",  "chSum",  "chHardSum", 
  "phMax",  "phSum",  "phHardSum", 
  "neMax",  "neSum", 
  "eMax",   "eSum", 
  "muMax",  "muSum", 
  "Aj",     "xj" , "dphi", "leadJetPt"
};
const int N_vars = sizeof(var)/sizeof(std::string);

// variable names, string array
const std::string var_xAx_Titles[] = {   
  "jtpt (GeV)" ,  "rawpt (GeV)",  
  "jteta", "jtphi (rad)", 
  "trkMax/rawpt", "trkMax/rawpt", "trkHardSum/rawpt", 
  "chMax/rawpt",  "chSum/rawpt",  "chHardSum/rawpt", 
  "phMax/rawpt",  "phSum/rawpt",  "phHardSum/rawpt", 
  "neMax/rawpt",  "neSum/rawpt", 
  "eMax/rawpt",   "eSum/rawpt", 
  "muMax/rawpt",  "muSum/rawpt", 
  "Aj",     "xj" , "dphi (radians)", "leadJetPt (GeV)"
};
const int N_vars_xAx_Titles = sizeof(var_xAx_Titles)/sizeof(std::string);

// variable names, string array
const std::string genVars[] = {   
  "pt","eta","phi" //,"drjt"
};
const int N_genVars = sizeof(genVars)/sizeof(std::string);
const std::string genVars_ptrat[] = {   
  "eta","phi" ,"drjt"
};
const int N_genVars_ptrat = sizeof(genVars_ptrat)/sizeof(std::string);

const std::string HLTName[] = {
  "HLT40","HLT60","HLT80","HLT100","HLTComb"// ,"TrgCombTest"
};
const int N_trigs=sizeof(HLTName)/sizeof(std::string);

//// MC+DATA COMMON BINNING
// ---------------------------------------------------------------------------------------------------------------

// binning arrays
const int ptbins[] = { 15, 30, 50, 80, 120, 170, 220, 300, 500 };
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one


const double JEC_ptbins[] = {
  17, 22, 27,    //15-30
  33, 39, 47,    //30-50
  55, 64, 74,    //50-80
  84, 97, 114,   //80-120
  133, 153,      //120-170
  174, 196,      //170-220
  220, 245, 272, //220-300
  300, 350, 400, //300-500
  550, 790, 1000 //500-inf
};
const int nbins_JEC_ptbins = sizeof(JEC_ptbins)/sizeof(double)-1;

const double etabins[] = {
  -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, 
  -3.839, -3.664, -3.489, -3.314, -3.139, 
  -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -2.043, 
  -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, 
  -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, 
  -0.087, +0.000, +0.087, 
  +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, +0.879, +0.957, 
  +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, +1.653, +1.740, +1.830, +1.930, 
  +2.043, +2.172, +2.322, +2.500, +2.650, +2.853, +2.964, 
  +3.139, +3.314, +3.489, +3.664, +3.839, 
  +4.013, +4.191, +4.363, +4.538, +4.716, +4.889, +5.191
};
const int nbins_eta = sizeof(etabins)/sizeof(double)-1;


// hist painting ------------------------
const int theDataOverlayMarkerColor=1, theMCOverlayMarkerColor=1,theRatioMarkerColor=1;//black
const int theDataOverlayLineColor=2, theMCOverlayLineColor=4,theRatioLineColor=9;//red, blue, olive

const int altOverlayLineColor=3; 
const int altRatioLineColor1=8, altRatioLineColor2=7;//green, violet or something

const int theTrigOverlayMarkerColor=1;
const int theTrigOverlayLineColor[]={9,8,7,2};//,3};


//// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------

void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){//binsX loop 
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//end nbinsX loop
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  return;
}
