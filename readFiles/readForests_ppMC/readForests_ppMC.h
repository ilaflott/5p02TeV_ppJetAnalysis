// original author: Raghav Kunnawalkam Elayavalli
// Nov 24th 2015
// Overhaul, Ian Laflotte
// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests_ppMC.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests_ppMC.exe

////////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
// global ROOT setups and includes
#include <TSystem.h>
#include <TProfile.h>
#include <TStyle.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TEventList.h>
#include <TCut.h>
// plotting/drawing, tools and functions
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
// graphs/errors
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
// hists
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
// fits
#include <TF1.h>
// misc utilities
#include <TMath.h>
#include <TRandom3.h>
#include <TStopwatch.h>

//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//// readForests_ppMC
const int defStartFile=0;
const int defEndFile=1; //inclusive boundary
const std::string defInFilelist = "filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
const int defRadius=4;
const std::string defJetType="PF";
const std::string defOutputName = "readForests_ppMC_defOut.root";
const bool defDebugMode = true;
int readForests_ppMC(int startfile = defStartFile , int endfile = defEndFile ,
                     std::string inFilelist = defInFilelist , std::string outfile = defOutputName ,
                     int radius = defRadius , std::string jetType=defJetType , bool debugMode = defDebugMode );
const int readFilesArgCount=7+minArgs;

/// other switches
const bool fillUnfoldingHists=false, fillMCEffHists=true;//, fillJECJERHists=false;

//// helper functions
float trigComb(bool *trg, int *pscl, float pt);
void divideBinWidth(TH1 *h);
float deltaphi(float phi1, float phi2);

//// CONSTANTS
// ---------------------------------------------------------------------------------------------------------------

// misc 
const char *etaWidth = (char*)"20_eta_20";

// binning arrays
const float ptbins[] = { 15., 30., 50., 80., 120., 170., 220., 300., 500. };
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one


const float pthatbins[]={ 15., 30., 50., 80., 120., 170., 220., 280., 370., 460., 540., 9999. };
const int nbins_pthat=sizeof(pthatbins)/sizeof(int)-1;

const float pthatWeights[]={//from evtPthatWeights/evtPthatWeights_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets.txt
/*(pthat>=15)&&(pthat<30), n[0]=932778, xsDiff=0.49235 mb^-3, weight=*/        5.27832e-07,
/*(pthat>=30)&&(pthat<50), n[1]=903567, xsDiff=0.030482 mb^-3, weight=*/       3.37352e-08,
/*(pthat>=50)&&(pthat<80), n[2]=983531, xsDiff=0.0035721 mb^-3, weight=*/      3.63191e-09,
/*(pthat>=80)&&(pthat<120), n[3]=1820782, xsDiff=0.00042494 mb^-3, weight=*/   2.33383e-10,
/*(pthat>=120)&&(pthat<170), n[4]=1080554, xsDiff=5.873e-05 mb^-3, weight=*/   5.43517e-11,
/*(pthat>=170)&&(pthat<220), n[5]=836152, xsDiff=9.199e-06 mb^-3, weight=*/    1.10016e-11,
/*(pthat>=220)&&(pthat<280), n[6]=954396, xsDiff=2.2564e-06 mb^-3, weight=*/   2.36422e-12,
/*(pthat>=280)&&(pthat<370), n[7]=1083994, xsDiff=6.336e-07 mb^-3, weight=*/   5.84505e-13,
/*(pthat>=370)&&(pthat<460), n[8]=948240, xsDiff=1.0884e-07 mb^-3, weight=*/   1.14781e-13,
/*(pthat>=460)&&(pthat<540), n[9]=1558268, xsDiff=2.215e-08 mb^-3, weight=*/   1.42145e-14,
/*(pthat>=540)&&(pthat<9999), n[10]=2597338, xsDiff=1.001e-08 mb^-3, weight=*/ 3.85395e-15
};

//these apply when comparing the MC to HLTAk4Calo-triggered jets in data
const float vzWeights[]={//from vzEvtWeights 9.14.16, //
/*for i=1, -15<vz<=-14.5, vzWeight=*/ 0.381083,
/*for i=2, -14.5<vz<=-14, vzWeight=*/ 0.420849,
/*for i=3, -14<vz<=-13.5, vzWeight=*/ 0.457166,
/*for i=4, -13.5<vz<=-13, vzWeight=*/ 0.489082,
/*for i=5, -13<vz<=-12.5, vzWeight=*/ 0.526434,
/*for i=6, -12.5<vz<=-12, vzWeight=*/ 0.568316,
/*for i=7, -12<vz<=-11.5, vzWeight=*/ 0.631824,
/*for i=8, -11.5<vz<=-11, vzWeight=*/ 0.627796,
/*for i=9, -11<vz<=-10.5, vzWeight=*/ 0.690459,
/*for i=10, -10.5<vz<=-10, vzWeight=*/0.714526,
/*for i=11, -10<vz<=-9.5, vzWeight=*/ 0.750061,
/*for i=12, -9.5<vz<=-9, vzWeight=*/  0.799916,
/*for i=13, -9<vz<=-8.5, vzWeight=*/  0.816735,
/*for i=14, -8.5<vz<=-8, vzWeight=*/  0.85323 ,
/*for i=15, -8<vz<=-7.5, vzWeight=*/  0.877455,
/*for i=16, -7.5<vz<=-7, vzWeight=*/  0.895425,
/*for i=17, -7<vz<=-6.5, vzWeight=*/  0.92577 ,
/*for i=18, -6.5<vz<=-6, vzWeight=*/  0.960376,
/*for i=19, -6<vz<=-5.5, vzWeight=*/  1.00063 ,
/*for i=20, -5.5<vz<=-5, vzWeight=*/  1.00567 ,
/*for i=21, -5<vz<=-4.5, vzWeight=*/  1.01408 ,
/*for i=22, -4.5<vz<=-4, vzWeight=*/  1.03837 ,
/*for i=23, -4<vz<=-3.5, vzWeight=*/  1.05649 ,
/*for i=24, -3.5<vz<=-3, vzWeight=*/  1.07153 ,
/*for i=25, -3<vz<=-2.5, vzWeight=*/  1.07986 ,
/*for i=26, -2.5<vz<=-2, vzWeight=*/  1.09942 ,
/*for i=27, -2<vz<=-1.5, vzWeight=*/  1.1102  ,
/*for i=28, -1.5<vz<=-1, vzWeight=*/  1.12039 ,
/*for i=29, -1<vz<=-0.5, vzWeight=*/  1.116   ,
/*for i=30, -0.5<vz<=0, vzWeight=*/   1.1358  ,
/*for i=31, 0<vz<=0.5, vzWeight=*/    1.13279 ,
/*for i=32, 0.5<vz<=1, vzWeight=*/    1.1387  ,
/*for i=33, 1<vz<=1.5, vzWeight=*/    1.13261 ,
/*for i=34, 1.5<vz<=2, vzWeight=*/    1.14604 ,
/*for i=35, 2<vz<=2.5, vzWeight=*/    1.12861 ,
/*for i=36, 2.5<vz<=3, vzWeight=*/    1.14573 ,
/*for i=37, 3<vz<=3.5, vzWeight=*/    1.12039 ,
/*for i=38, 3.5<vz<=4, vzWeight=*/    1.11589 ,
/*for i=39, 4<vz<=4.5, vzWeight=*/    1.11258 ,
/*for i=40, 4.5<vz<=5, vzWeight=*/    1.10108 ,
/*for i=41, 5<vz<=5.5, vzWeight=*/    1.09245 ,
/*for i=42, 5.5<vz<=6, vzWeight=*/    1.08902 ,
/*for i=43, 6<vz<=6.5, vzWeight=*/    1.06721 ,
/*for i=44, 6.5<vz<=7, vzWeight=*/    1.0261  ,
/*for i=45, 7<vz<=7.5, vzWeight=*/    1.01827 ,
/*for i=46, 7.5<vz<=8, vzWeight=*/    0.989024,
/*for i=47, 8<vz<=8.5, vzWeight=*/    0.962674,
/*for i=48, 8.5<vz<=9, vzWeight=*/    0.936773,
/*for i=49, 9<vz<=9.5, vzWeight=*/    0.907251,
/*for i=50, 9.5<vz<=10, vzWeight=*/   0.891311,
/*for i=51, 10<vz<=10.5, vzWeight=*/  0.864917,
/*for i=52, 10.5<vz<=11, vzWeight=*/  0.807905,
/*for i=53, 11<vz<=11.5, vzWeight=*/  0.787968,
/*for i=54, 11.5<vz<=12, vzWeight=*/  0.725249,
/*for i=55, 12<vz<=12.5, vzWeight=*/  0.709917,
/*for i=56, 12.5<vz<=13, vzWeight=*/  0.661411,
/*for i=57, 13<vz<=13.5, vzWeight=*/  0.628832,
/*for i=58, 13.5<vz<=14, vzWeight=*/  0.580118,
/*for i=59, 14<vz<=14.5, vzWeight=*/  0.550339,
/*for i=60, 14.5<vz<=15, vzWeight=*/  0.479087
};
const int nbins_vzWeights = sizeof(vzWeights)/sizeof(float);//should be 60
const float minbinValue_vzWeights = -15., maxbinValue_vzWeights = 15.;
const float binsize_vzWeights = (maxbinValue_vzWeights-minbinValue_vzWeights)/nbins_vzWeights;//should be 0.5
//const float binsize_vzWeights = 0.5;

const float JEC_ptbins[] = {
  17., 22., 27.,    //15-30
  33., 39., 47.,    //30-50
  55., 64., 74.,    //50-80
  84., 97., 114.,   //80-120
  133., 153.,      //120-170
  174., 196.,      //170-220
  220., 245., 272., //220-300
  300., 350., 400., //300-500
  550., 790., 1000. //500-1000
};
const int nbins_JEC_ptbins = sizeof(JEC_ptbins)/sizeof(float)-1;

const float etabins[] = {
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
const int nbins_eta = sizeof(etabins)/sizeof(float)-1;

// tree names+directories
const std::string treeNames[]={ 
  "GARBAGE ENTRY" , //akPFJetAna of arb. radius, form strin
  //  "GARBAGE ENTRY" , //akCaloJet of same radius
  "hiEvtAnalyzer/HiTree" ,
  "skimanalysis/HltTree" ,
  "hltanalysis/HltTree" 
//  "hltobject/"+HLTBitStrings[0]+"_v" , 
//  "hltobject/"+HLTBitStrings[1]+"_v" , 
//  "hltobject/"+HLTBitStrings[2]+"_v" , 
//  "hltobject/"+HLTBitStrings[3]+"_v"   
}; 
const int N_trees = sizeof(treeNames)/sizeof(std::string);

// variable names for QA Plots
const std::string var[] = {   
  "jtpt" ,  "rawpt",  //jet weights only
  "jteta", "jtphi", //require normalization
  "trkMax", "trkSum", "trkHardSum", 
  "chMax",  "chSum",  "chHardSum", 
  "phMax",  "phSum",  "phHardSum", 
  "neMax",  "neSum", 
  "eMax",   "eSum", 
  "muMax",  "muSum", 
  "Aj",     "xj",      "dphi", //dijets
  "leadJetPt", "subleadJetPt"
};
const int N_vars = sizeof(var)/sizeof(std::string);


// do not use yet, for when we have MC with hlt emulation!
////L1
//const std::string L1BitStrings[]={//this array is a good idea
//  "L1_SingleJet28_BptxAND",
//  "L1_SingleJet40_BptxAND",
//  "L1_SingleJet48_BptxAND",   
//  "L1_SingleJet52_BptxAND"
//};
//const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);
//
//const std::string L1PresclBranches[]={//this array is a good idea
//  L1BitStrings[0]+"_Prescl",
//  L1BitStrings[1]+"_Prescl",
//  L1BitStrings[2]+"_Prescl",
//  L1BitStrings[3]+"_Prescl"
//};
//
//// HLT                           
//const std::string HLTJetType="Calo";//"PF";
//const std::string HLTBitStrings[]={	
//  "HLT_AK4"+HLTJetType+"Jet40_Eta5p1",
//  "HLT_AK4"+HLTJetType+"Jet60_Eta5p1",			
//  "HLT_AK4"+HLTJetType+"Jet80_Eta5p1",			
//  "HLT_AK4"+HLTJetType+"Jet100_Eta5p1"    		
//};						
////const int N_HLTBits=sizeof(HLTBitStrings)/sizeof(std::string);
//const int N_HLTBits=N_L1Bits;
//
//const std::string HLTPresclBranches[]={	
//  HLTBitStrings[0]+"_v1_Prescl",
//  HLTBitStrings[1]+"_v1_Prescl",
//  HLTBitStrings[2]+"_v1_Prescl",
//  HLTBitStrings[3]+"_v1_Prescl"
//};			
//			
//const std::string HLTBranches[]={	
//  HLTBitStrings[0]+"_v1",
//  HLTBitStrings[1]+"_v1",
//  HLTBitStrings[2]+"_v1",
//  HLTBitStrings[3]+"_v1"
//};						
//
//// tree names+directories
//const std::string treeNames[]={ 
//  "GARBAGE ENTRY" , //use jet ana of choice later
//  "hiEvtAnalyzer/HiTree" ,
//  "skimanalysis/HltTree" ,
//  "hltanalysis/HltTree" ,
//  "hltobject/"+HLTBitStrings[0]+"_v" , 
//  "hltobject/"+HLTBitStrings[1]+"_v" , 
//  "hltobject/"+HLTBitStrings[2]+"_v" , 
//  "hltobject/"+HLTBitStrings[3]+"_v"   
//}; 
//const int N_trees = sizeof(treeNames)/sizeof(std::string);

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


double trigComb(bool *trg, int *pscl, double pt){
  float weight=0.;
  if( trg[3] && pt>=100.            ) weight = pscl[3];
  if( trg[2] && pt>=80.  && pt<100. ) weight = pscl[2];
  if( trg[1] && pt>=60.  && pt<80.  ) weight = pscl[1];
  if( trg[0] && pt>=40.  && pt<60.  ) weight = pscl[0];
  return weight;
}


float deltaphi(float phi1, float phi2){
  float pi=TMath::Pi(); 
  float dphi=TMath::Abs(phi1-phi2);
  if(dphi > pi)dphi -= 2.*pi;
  return TMath::Abs(dphi);
}


//random old binning i found in raghavs code, here for safe keeping
//const float ptbins[] = {
//  3,   4,   5,   7,   9,  12,  15,   18,  21,  24,  28,
//  28,  32,  37,  43,  49,  56,  64,   74,  84,  97, 114,
//  133, 153, 174, 196, 220, 245, 272,  300, 330, 362, 395, 
//  430, 468, 507, 548, 592, 638, 686, 1000 
//}; 

//const int pthatbins[]={ 15, 30, 50, 80, 120, 170, 220, 280, 9999 }; //2015 bins?

