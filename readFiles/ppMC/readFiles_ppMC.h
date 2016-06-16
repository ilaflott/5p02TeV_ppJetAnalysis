// based on RAA_read_mc_pp.C by Raghav Kunnawalkam Elayavalli
// current iteration written by Ian Laflotte
// Jun 13 th 2016
// reads and writes jets from pp MC forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readFiles_ppMC.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readFiles_ppMC.exe

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

//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//// readFiles_ppMC
const int defStartFile=0;
const int defEndFile=1; //inclusive boundary
const std::string defInFilelist="../../filelists/ppMC/5p02TeV_Py8_CUETP8M1_QCDjet15_20Files_debug_forests.txt";
const int defRadius=4;
const std::string defJetType="PF";
const std::string defOutputName="readFiles_ppMC_defOut.root";
const bool defDebugMode=true;
int readFiles_ppMC(int startfile=defStartFile , int endfile=defEndFile ,
                     std::string inFilelist=defInFilelist , std::string outfile=defOutputName ,
                     int radius=defRadius , std::string jetType=defJetType , bool debugMode=defDebugMode );
const int readFilesArgCount=7+minArgs;

//// helper functions
double trigComb(bool *trg, int *pscl, double pt);
void divideBinWidth(TH1 *h);
float deltaphi(float phi1, float phi2);

//// CONSTANTS
// ---------------------------------------------------------------------------------------------------------------
// useful, code/function-wide constants

// for convenience during testing, coding, etc.
const bool atMIT=true;
const std::string hadoopDir="/mnt/hadoop/cms";
const std::string xrootdDirect="root://cmsxrootd.fnal.gov/";

// pt binning
const int ptbins[]={ 15, 30, 50, 80, 120, 170, 220, 300, 500 };
const int nbins_pt=sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one

// jec
const double ptbins_jec[]={
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
const int nbins_pt_jec=sizeof(ptbins_jec)/sizeof(double)-1;

// eta binning
const double etabins[]={
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
const int nbins_eta=sizeof(etabins)/sizeof(double)-1;

//static const char *etaWidth=(char*)"20_eta_20";
const char *etaWidth=(char*)"20_eta_20";

// root file directories + tree names
const std::string treeNames[]={ 
  "GARBAGE ENTRY" , //make ak4PFJetAnalyzer string and replace later
  //"GARBAGE ENTRY" , //make ak4CaloJetAnalyzer string and replace later?
  "skimanalysis/HltTree"  ,
  "hltanalysis/HltTree"   ,
  "hiEvtAnalyzer/HiTree"  ,
  //"ppTrack"       ,
  //"hltobject/HLT_AK4CaloJet40_Eta5p1_v"     , //hlt40
  //"hltobject/HLT_AK4CaloJet60_Eta5p1_v"     , //hlt60
  //"hltobject/HLT_AK4CaloJet80_Eta5p1_v"     , //hlt80
  //"hltobject/HLT_AK4CaloJet100_Eta5p1_v"     //hlt100
}; 
const int N_trees=sizeof(treeNames)/sizeof(std::string);

// Jet variable names
const std::string var[]={   
  "jtpt" ,  "rawpt",  "jteta", "jtphi", 
  "trkMax", "trkSum", "trkHardSum", 
  "chMax",  "chSum",  "chHardSum", 
  "phMax",  "phSum",  "phHardSum", 
  "neMax",  "neSum", 
  "eMax",   "eSum", 
  "muMax",  "muSum", 
  "Aj",     "xj" 
};
const int N_vars=sizeof(var)/sizeof(std::string);

const std::string trigNames[]={"HLT40","HLT60","HLT80","HLT100","Combined"};
const int trigValue=sizeof(trigNames)/sizeof(std::string);

//// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// mini functions not called by the frontend/main

// divide by bin width //
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){//binsX loop 
    Float_t val=h->GetBinContent(i);
    Float_t valErr=h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//end nbinsX loop
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  return;
}

// trigger combination //
double trigComb(bool *trg, int *pscl, double pt){
  double weight=0;
  if(trg[3] && pt>=100 )          weight=pscl[3];
  if(trg[2] && pt>=80  && pt<100) weight=pscl[2];
  if(trg[1] && pt>=60  && pt<80 ) weight=pscl[1];
  if(trg[0] && pt>=40  && pt<60 ) weight=pscl[0];
  return weight;
}

// delta phi //
float deltaphi(float phi1, float phi2){
  float pi=TMath::Pi(); 
  float dphi=TMath::Abs(phi1-phi2);
  if(dphi > pi)dphi -= 2*pi;
  return TMath::Abs(dphi);
}
