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
const std::string defInFilelist = "filelists/5p02TeV_Py8_CUETP8M1_QCDjet15_20Files_debug_forests.txt";
const int defRadius=4;
const std::string defJetType="PF";
const std::string defOutputName = "readForests_ppMC_defOut.root";
const bool defDebugMode = true;
int readForests_ppMC(int startfile = defStartFile , int endfile = defEndFile ,
                     std::string inFilelist = defInFilelist , std::string outfile = defOutputName ,
                     int radius = defRadius , std::string jetType=defJetType , bool debugMode = defDebugMode );
const int readFilesArgCount=7+minArgs;

//// helper functions
double trigComb(bool *trg, int *pscl, double pt);
void divideBinWidth(TH1 *h);
float deltaphi(float phi1, float phi2);

//// CONSTANTS
// ---------------------------------------------------------------------------------------------------------------

// misc 
const char *etaWidth = (char*)"20_eta_20";

// binning arrays
//const double ptbins[] = {
//  3,   4,   5,   7,   9,  12,  15,   18,  21,  24,  28,
//  28,  32,  37,  43,  49,  56,  64,   74,  84,  97, 114,
//  133, 153, 174, 196, 220, 245, 272,  300, 330, 362, 395, 
//  430, 468, 507, 548, 592, 638, 686, 1000 
//}; //random old binning i found in raghavs code, here for safe keeping

const int ptbins[] = { 15, 30, 50, 80, 120, 170, 220, 300, 500 };
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one

//const int pthatbins[]={ 15, 30, 50, 80, 120, 170, 220, 280, 9999 }; //2015 bins?
const int pthatbins[]={ 15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 9999 };
const int nbins_pthat=sizeof(pthatbins)/sizeof(int)-1;

const double pthatWeights[]={//from jetWeights/jetWeights_Py8_CUETP8M1_QCDjetAllPtBins.txt
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

const double vzWeights[]={//from doVzWeights on jtpt15_ldJetPt50 HPtJetTrigData/Py8QCDMC samples in printQAPlots.C 8.31.16
/*or i=1, -15<vz<=-14.5, vzWeight=*/   0.381095,
/*  or i=2, -14.5<vz<=-14, vzWeight=*/ 0.420839,
/*  or i=3, -14<vz<=-13.5, vzWeight=*/ 0.457158,
/*  or i=4, -13.5<vz<=-13, vzWeight=*/ 0.489084,
/*  or i=5, -13<vz<=-12.5, vzWeight=*/ 0.526421,
/*  or i=6, -12.5<vz<=-12, vzWeight=*/ 0.568266,
/*  or i=7, -12<vz<=-11.5, vzWeight=*/ 0.631775,
/*  or i=8, -11.5<vz<=-11, vzWeight=*/ 0.627823,
/*  or i=9, -11<vz<=-10.5, vzWeight=*/ 0.690401,
/*  or i=10, -10.5<vz<=-10, vzWeight=*/0.714522,
/*  or i=11, -10<vz<=-9.5, vzWeight=*/ 0.750037,
/*  or i=12, -9.5<vz<=-9, vzWeight=*/  0.799898,
/*  or i=13, -9<vz<=-8.5, vzWeight=*/  0.816712,
/*  or i=14, -8.5<vz<=-8, vzWeight=*/  0.8532  ,
/*  or i=15, -8<vz<=-7.5, vzWeight=*/  0.877429,
/*  or i=16, -7.5<vz<=-7, vzWeight=*/  0.895462,
/*  or i=17, -7<vz<=-6.5, vzWeight=*/  0.925748,
/*  or i=18, -6.5<vz<=-6, vzWeight=*/  0.960379,
/*  or i=19, -6<vz<=-5.5, vzWeight=*/  1.00065 ,
/*  or i=20, -5.5<vz<=-5, vzWeight=*/  1.00565 ,
/*  or i=21, -5<vz<=-4.5, vzWeight=*/  1.01403 ,
/*  or i=22, -4.5<vz<=-4, vzWeight=*/  1.03837 ,
/*  or i=23, -4<vz<=-3.5, vzWeight=*/  1.05647 ,
/*  or i=24, -3.5<vz<=-3, vzWeight=*/  1.07155 ,
/*  or i=25, -3<vz<=-2.5, vzWeight=*/  1.07985 ,
/*  or i=26, -2.5<vz<=-2, vzWeight=*/  1.09945 ,
/*  or i=27, -2<vz<=-1.5, vzWeight=*/  1.1102  ,
/*  or i=28, -1.5<vz<=-1, vzWeight=*/  1.12035 ,
/*  or i=29, -1<vz<=-0.5, vzWeight=*/  1.11602 ,
/*  or i=30, -0.5<vz<=0, vzWeight=*/   1.1358  ,
/*  or i=31, 0<vz<=0.5, vzWeight=*/    1.1328  ,
/*  or i=32, 0.5<vz<=1, vzWeight=*/    1.13871 ,
/*  or i=33, 1<vz<=1.5, vzWeight=*/    1.13266 ,
/*  or i=34, 1.5<vz<=2, vzWeight=*/    1.14604 ,
/*  or i=35, 2<vz<=2.5, vzWeight=*/    1.12859 ,
/*  or i=36, 2.5<vz<=3, vzWeight=*/    1.14574 ,
/*  or i=37, 3<vz<=3.5, vzWeight=*/    1.12041 ,
/*  or i=38, 3.5<vz<=4, vzWeight=*/    1.11591 ,
/*  or i=39, 4<vz<=4.5, vzWeight=*/    1.11259 ,
/*  or i=40, 4.5<vz<=5, vzWeight=*/    1.1011  ,
/*  or i=41, 5<vz<=5.5, vzWeight=*/    1.09243 ,
/*  or i=42, 5.5<vz<=6, vzWeight=*/    1.08899 ,
/*  or i=43, 6<vz<=6.5, vzWeight=*/    1.06719 ,
/*  or i=44, 6.5<vz<=7, vzWeight=*/    1.02613 ,
/*  or i=45, 7<vz<=7.5, vzWeight=*/    1.0183  ,
/*  or i=46, 7.5<vz<=8, vzWeight=*/    0.989017,
/*  or i=47, 8<vz<=8.5, vzWeight=*/    0.962692,
/*  or i=48, 8.5<vz<=9, vzWeight=*/    0.936776,
/*  or i=49, 9<vz<=9.5, vzWeight=*/    0.90725 ,
/*  or i=50, 9.5<vz<=10, vzWeight=*/   0.89135 ,
/*  or i=51, 10<vz<=10.5, vzWeight=*/  0.864868,
/*  or i=52, 10.5<vz<=11, vzWeight=*/  0.807878,
/*  or i=53, 11<vz<=11.5, vzWeight=*/  0.788   ,
/*  or i=54, 11.5<vz<=12, vzWeight=*/  0.725265,
/*  or i=55, 12<vz<=12.5, vzWeight=*/  0.709917,
/*  or i=56, 12.5<vz<=13, vzWeight=*/  0.661436,
/*  or i=57, 13<vz<=13.5, vzWeight=*/  0.628792,
/*  or i=58, 13.5<vz<=14, vzWeight=*/  0.580058,
/*  or i=59, 14<vz<=14.5, vzWeight=*/  0.550358,
/*  or i=60, 14.5<vz<=15, vzWeight=*/  0.479061//,
///*  or i=61, 15.<vz<inf, vzWeight=*/  0.
};
const int nbins_vzWeights = sizeof(vzWeights)/sizeof(double);//should be 60
const double minbinValue_vzWeights = -15.;
const double maxbinValue_vzWeights = 15.;
const double binsize_vzWeights = (maxbinValue_vzWeights-minbinValue_vzWeights)/nbins_vzWeights;//should be 0.5
//const double binsize_vzWeights = 0.5;

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

// tree names+directories
const std::string treeNames[]={ 
  "GARBAGE ENTRY" , //akPfJetAna of arb. radius, form strin
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
  "leadJetPt"
};
const int N_vars = sizeof(var)/sizeof(std::string);


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
  double weight=0;
  if(trg[3] && pt>=100 )          weight = pscl[3];
  if(trg[2] && pt>=80  && pt<100) weight = pscl[2];
  if(trg[1] && pt>=60  && pt<80 ) weight = pscl[1];
  if(trg[0] && pt>=40  && pt<60 ) weight = pscl[0];
  return weight;
}


float deltaphi(float phi1, float phi2){
  float pi=TMath::Pi(); 
  float dphi=TMath::Abs(phi1-phi2);
  if(dphi > pi)dphi -= 2*pi;
  return TMath::Abs(dphi);
}
