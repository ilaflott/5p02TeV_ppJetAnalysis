// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

////////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
// ROOTSYS
#include <TSystem.h>
#include <TProfile.h>
#include <TStyle.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
// utilities
#include <TMath.h>
#include <TStopwatch.h>
// histos
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
// custom
#include "L2L3ResidualWFits.h"


//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;//, fastDebugMode = true;
const std::string defMCOutputName="readForests_ppMC_defOut";//.root";

//int readForests_ppMC_MCJEC( std::string inFilelist="filelists/test_readForests_ppMC_local.txt",
int readForests_ppMC_MCJEC( std::string inFilelist=defMCInFilelist,
			    int startfile=0, int endfile=20, 
			    int radius=defRadius, std::string jetType=defJetType, 
			    bool debugMode=defDebugMode,
			    std::string outfile=(defMCOutputName+"_MCJEC.root")      );

const int readForestsArgCount=7+minArgs;

const float jtPtCut=15.;
const float jtEtaCutLo=0.0, jtEtaCutHi=4.7;//make htEtaCutHi 100. for inf.
const float jetQAPtCut=15.;//50.;
const float ldJetPtCut=20., subldJetPtCut=10., ptAveCut=15., dPhiCut=2./3.*TMath::Pi();//dijet cuts

//// HELPER FUNCTIONS
// -------------------------------------------------------------------------------------------------------------

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

float trigComb(bool *trgDec, int *treePrescl, double triggerPt){
  float weight_eS=0.;
  if(trgDec[0] && triggerPt>=40.  && triggerPt<60. ) weight_eS=treePrescl[0];
  if(trgDec[1] && triggerPt>=60.  && triggerPt<80. ) weight_eS=treePrescl[1];
  if(trgDec[2] && triggerPt>=80.  && triggerPt<100.) weight_eS=treePrescl[2];
  if(trgDec[3] && triggerPt>=100. )                  weight_eS=treePrescl[3];
  return weight_eS;
}

float deltaphi(float phi1, float phi2){
  float pi=TMath::Pi();
  float dphi=TMath::Abs(phi1-phi2);
  if(dphi > pi)dphi -=2.*pi;
  return TMath::Abs(dphi);
}



//// STRINGS
// ------------------------------------------------------------------------------------------------------------

// eta width
//const char *etaWidth=(char*)"20_eta_20";
const std::string etaWidth="20_eta_20";

// variable names for QA Plots
const std::string var[]={
  "jtpt"  , "rawpt" ,  //jets, 0-4=5
  "jteta" , "jtphi" , "jty",
  "trkMax", "trkSum", "trkHardSum", //jet constituents, 5-22=18
  "chMax" , "chSum" , "chHardSum",
  "phMax" , "phSum" , "phHardSum",
  "neMax" , "neSum" ,
  "eMax"  , "eSum"  ,
  "muMax" , "muSum" ,
  "neN", "chN", "sumN",
  "Aj" , "xj" , "dphi", //dijet variables, 23-27=5
  "leadJetPt", "subleadJetPt"
};
const int N_vars=sizeof(var)/sizeof(std::string);

//L1
const std::string L1BitStrings[]={//this array is a good idea
  "L1_SingleJet28_BptxAND",
  "L1_SingleJet40_BptxAND",
  "L1_SingleJet48_BptxAND",
  "L1_SingleJet52_BptxAND"//,
  //  "L1_SingleJet52_BptxAND",
  //  "L1_SingleJet68_BptxAND"
};
//const int N_L1Bits=N_HLTBits;
const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);

// HLT, both Calo and PF available
const std::string Calo_HLTBitStrings[]={
  "HLT_AK4CaloJet40_Eta5p1",
  "HLT_AK4CaloJet60_Eta5p1",
  "HLT_AK4CaloJet80_Eta5p1",
  "HLT_AK4CaloJet100_Eta5p1" //,
  //  "HLT_AK4CaloJet110_Eta5p1" ,
  //  "HLT_AK4CaloJet120_Eta5p1"
};
const int N_HLTBits=sizeof(Calo_HLTBitStrings)/sizeof(std::string);

const std::string PF_HLTBitStrings[]={
  "HLT_AK4PFJet40_Eta5p1",
  "HLT_AK4PFJet60_Eta5p1",
  "HLT_AK4PFJet80_Eta5p1",
  "HLT_AK4PFJet100_Eta5p1" //,
  //  "HLT_AK4PFJet110_Eta5p1" ,
  //  "HLT_AK4PFJet120_Eta5p1"
};

// data tree name array
const std::string dataTreeNames[]={
  "GARBAGE ENTRY",
  "hiEvtAnalyzer/HiTree", "skimanalysis/HltTree", "hltanalysis/HltTree",
  "hltobject/"+Calo_HLTBitStrings[0]+"_v" ,
  "hltobject/"+Calo_HLTBitStrings[1]+"_v" ,
  "hltobject/"+Calo_HLTBitStrings[2]+"_v" ,
  "hltobject/"+Calo_HLTBitStrings[3]+"_v" //, 
  //  "hltobject/"+Calo_HLTBitStrings[4]+"_v"  
  //  "hltobject/"+Calo_HLTBitStrings[5]+"_v"  
};
const int N_dataTrees=sizeof(dataTreeNames)/sizeof(std::string);

const std::string trgCombType="Calo";

// MC tree name array (NOT FINISHED, WILL ERROR IN CURRENT FORM)
const std::string MCTreeNames[]={
  "GARBAGE ENTRY",
  "hiEvtAnalyzer/HiTree", "skimanalysis/HltTree"//, "hltanalysis/HltTree",
  //  "hltobject/"+Calo_HLTBitStrings[0]+"_v" ,
  //  "hltobject/"+Calo_HLTBitStrings[1]+"_v" ,
  //  "hltobject/"+Calo_HLTBitStrings[2]+"_v" ,
  //  "hltobject/"+Calo_HLTBitStrings[3]+"_v" //, 
  //  "hltobject/"+Calo_HLTBitStrings[4]+"_v"  
  //  "hltobject/"+Calo_HLTBitStrings[5]+"_v"  
};
const int N_MCTrees=sizeof(MCTreeNames)/sizeof(std::string);
  


//// BINNING
// ------------------------------------------------------------------------------------------------------


const float ptbins[]={ 
  15., 30., 50., 80., 
  120., 170., 220., 300., 
  500. 
}; // original genpt JER binning i inherited from raghav
const int nbins_pt=sizeof(ptbins)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one

const float ptbins2[]={
     3., 4., 5., 7., 9., 12., 
     15., 18., 21., 24., 28.,
     32., 37., 43., 49., 56.,
     64., 74., 84., 97., 114.,
     133., 153., 174., 196.,
     220., 245., 272., 300., 
     330., 362., 395., 430.,
     468., 507., 548., 592.,
     638., 686., 1000.//, 1500
}; //raghavs suggested genpt binning for JER
const int nbins_pt2=sizeof(ptbins2)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one

const float JEC_ptbins[]={
  17., 22., 27.,    //15-30
  33., 39., 47.,    //30-50
  55., 64., 74.,    //50-80
  84., 97., 114.,   //80-120
  133., 153.,      //120-170
  174., 196.,      //170-220
  220., 245., 272., //220-300
  300., 350., 400., //300-500
  550., 790., 1000. //500-inf
};
const int nbins_JEC_ptbins=sizeof(JEC_ptbins)/sizeof(float)-1;

const float etabins[]={
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
const int nbins_eta=sizeof(etabins)/sizeof(float)-1;

const float rapbins[]={ 
  0.0, 0.5,
  1.0, 1.5,
  2.0, 2.5,
  3.0, 3.2,
  4.7	     
};	     
const int nbins_rap=sizeof(rapbins)/sizeof(float)-1;

//// WEIGHTS FOR MC
// ------------------------------------------------------------------------------------------------------

// pthat binning and weights
const float pthatbins[]={ 15., 30., 50., 80., 120., 170., 220., 280., 370., 460., 540., 9999. };
//const int pthatbins[]={ 15, 30, 50, 80, 120, 170, 220, 280, 9999 }; //2015 bins?
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

// vz binning and weights
const float vzWeights[]={//           HLTAK4CaloJets, HLTAK4PFJets
  /*for i=1, -15<vz<=-14.5, vzWeight=*/ 0.381083,  /* 0.35749,  */
  /*for i=2, -14.5<vz<=-14, vzWeight=*/ 0.420849,  /* 0.400545, */
  /*for i=3, -14<vz<=-13.5, vzWeight=*/ 0.457166,  /* 0.433741, */
  /*for i=4, -13.5<vz<=-13, vzWeight=*/ 0.489082,  /* 0.469437, */
  /*for i=5, -13<vz<=-12.5, vzWeight=*/ 0.526434,  /* 0.512839, */
  /*for i=6, -12.5<vz<=-12, vzWeight=*/ 0.568316,  /* 0.556416, */
  /*for i=7, -12<vz<=-11.5, vzWeight=*/ 0.631824,  /* 0.612334, */
  /*for i=8, -11.5<vz<=-11, vzWeight=*/ 0.627796,  /* 0.622496, */
  /*for i=9, -11<vz<=-10.5, vzWeight=*/ 0.690459,  /* 0.672646, */
  /*for i=10, -10.5<vz<=-10, vzWeight=*/0.714526,  /* 0.706209, */
  /*for i=11, -10<vz<=-9.5, vzWeight=*/ 0.750061,  /* 0.740972, */
  /*for i=12, -9.5<vz<=-9, vzWeight=*/  0.799916,  /* 0.792563, */
  /*for i=13, -9<vz<=-8.5, vzWeight=*/  0.816735,  /* 0.815354, */
  /*for i=14, -8.5<vz<=-8, vzWeight=*/  0.85323 ,  /* 0.853215, */
  /*for i=15, -8<vz<=-7.5, vzWeight=*/  0.877455,  /* 0.875179, */
  /*for i=16, -7.5<vz<=-7, vzWeight=*/  0.895425,  /* 0.891705, */
  /*for i=17, -7<vz<=-6.5, vzWeight=*/  0.92577 ,  /* 0.921846, */
  /*for i=18, -6.5<vz<=-6, vzWeight=*/  0.960376,  /* 0.958741, */
  /*for i=19, -6<vz<=-5.5, vzWeight=*/  1.00063 ,  /* 1.00373,  */
  /*for i=20, -5.5<vz<=-5, vzWeight=*/  1.00567 ,  /* 1.00595,  */
  /*for i=21, -5<vz<=-4.5, vzWeight=*/  1.01408 ,  /* 1.01643,  */
  /*for i=22, -4.5<vz<=-4, vzWeight=*/  1.03837 ,  /* 1.04327,  */
  /*for i=23, -4<vz<=-3.5, vzWeight=*/  1.05649 ,  /* 1.06326,  */
  /*for i=24, -3.5<vz<=-3, vzWeight=*/  1.07153 ,  /* 1.07991,  */
  /*for i=25, -3<vz<=-2.5, vzWeight=*/  1.07986 ,  /* 1.08791,  */
  /*for i=26, -2.5<vz<=-2, vzWeight=*/  1.09942 ,  /* 1.10771,  */
  /*for i=27, -2<vz<=-1.5, vzWeight=*/  1.1102  ,  /* 1.11622,  */
  /*for i=28, -1.5<vz<=-1, vzWeight=*/  1.12039 ,  /* 1.12078,  */
  /*for i=29, -1<vz<=-0.5, vzWeight=*/  1.116   ,  /* 1.12084,  */
  /*for i=30, -0.5<vz<=0, vzWeight=*/   1.1358  ,  /* 1.14016,  */
  /*for i=31, 0<vz<=0.5, vzWeight=*/    1.13279 ,  /* 1.14058,  */
  /*for i=32, 0.5<vz<=1, vzWeight=*/    1.1387  ,  /* 1.14684,  */
  /*for i=33, 1<vz<=1.5, vzWeight=*/    1.13261 ,  /* 1.14148,  */
  /*for i=34, 1.5<vz<=2, vzWeight=*/    1.14604 ,  /* 1.15419,  */
  /*for i=35, 2<vz<=2.5, vzWeight=*/    1.12861 ,  /* 1.1392,   */
  /*for i=36, 2.5<vz<=3, vzWeight=*/    1.14573 ,  /* 1.15569,  */
  /*for i=37, 3<vz<=3.5, vzWeight=*/    1.12039 ,  /* 1.13269,  */
  /*for i=38, 3.5<vz<=4, vzWeight=*/    1.11589 ,  /* 1.12478,  */
  /*for i=39, 4<vz<=4.5, vzWeight=*/    1.11258 ,  /* 1.11415,  */
  /*for i=40, 4.5<vz<=5, vzWeight=*/    1.10108 ,  /* 1.10534,  */
  /*for i=41, 5<vz<=5.5, vzWeight=*/    1.09245 ,  /* 1.09492,  */
  /*for i=42, 5.5<vz<=6, vzWeight=*/    1.08902 ,  /* 1.08401,  */
  /*for i=43, 6<vz<=6.5, vzWeight=*/    1.06721 ,  /* 1.06167,  */
  /*for i=44, 6.5<vz<=7, vzWeight=*/    1.0261  ,  /* 1.02492,  */
  /*for i=45, 7<vz<=7.5, vzWeight=*/    1.01827 ,  /* 1.01539,  */
  /*for i=46, 7.5<vz<=8, vzWeight=*/    0.989024,  /* 0.985653, */
  /*for i=47, 8<vz<=8.5, vzWeight=*/    0.962674,  /* 0.952256, */
  /*for i=48, 8.5<vz<=9, vzWeight=*/    0.936773,  /* 0.930472, */
  /*for i=49, 9<vz<=9.5, vzWeight=*/    0.907251,  /* 0.894218, */
  /*for i=50, 9.5<vz<=10, vzWeight=*/   0.891311,  /* 0.875605, */
  /*for i=51, 10<vz<=10.5, vzWeight=*/  0.864917,  /* 0.847596, */
  /*for i=52, 10.5<vz<=11, vzWeight=*/  0.807905,  /* 0.78767,  */
  /*for i=53, 11<vz<=11.5, vzWeight=*/  0.787968,  /* 0.768839, */
  /*for i=54, 11.5<vz<=12, vzWeight=*/  0.725249,  /* 0.700327, */
  /*for i=55, 12<vz<=12.5, vzWeight=*/  0.709917,  /* 0.687202, */
  /*for i=56, 12.5<vz<=13, vzWeight=*/  0.661411,  /* 0.629434, */
  /*for i=57, 13<vz<=13.5, vzWeight=*/  0.628832,  /* 0.597895, */
  /*for i=58, 13.5<vz<=14, vzWeight=*/  0.580118,  /* 0.555087, */
  /*for i=59, 14<vz<=14.5, vzWeight=*/  0.550339,  /* 0.523324, */
  /*for i=60, 14.5<vz<=15, vzWeight=*/  0.479087   /* 0.448861  */
};
const int nbins_vzWeights=sizeof(vzWeights)/sizeof(float);//should be 60
const float minbinValue_vzWeights=-15., maxbinValue_vzWeights=15.;
const float binsize_vzWeights=(maxbinValue_vzWeights-minbinValue_vzWeights)/nbins_vzWeights;//should be 0.5

//random old binning i found in raghavs code, here for safe keeping
//const float ptbins[]={
//  3,   4,   5,   7,   9,  12,  15,   18,  21,  24,  28,
//  28,  32,  37,  43,  49,  56,  64,   74,  84,  97, 114,
//  133, 153, 174, 196, 220, 245, 272,  300, 330, 362, 395,
//  430, 468, 507, 548, 592, 638, 686, 1000
//};

//float trigComb_sanityCheck(bool *trgDec, int *treePrescl, double triggerPt){
//  std::cout<<"trigComb_sanityCheck Called..."<<std::endl<<std::endl;
//
//  float weight_eS=0.;
//  std::cout<<"before if statements..."<<std::endl;
//  for(int k=0;k<4;k++){if(trgDec[k]) {
//      std::cout<<"trgDec["<<k<<"]="<<trgDec[k]<<std::endl;
//      std::cout<<"treePrescl["<<k<<"]="<<treePrescl[k]<<std::endl;
//      std::cout<<"triggerPt="<<triggerPt<<std::endl;
//      std::cout<<"weight_eS=="<<weight_eS<<std::endl<<std::endl;
//    }}
//
//  if(trgDec[0] && triggerPt>=40.  && triggerPt<60. ) weight_eS=treePrescl[0];
//  if(trgDec[1] && triggerPt>=60.  && triggerPt<80. ) weight_eS=treePrescl[1];
//  if(trgDec[2] && triggerPt>=80.  && triggerPt<100.) weight_eS=treePrescl[2];
//  if(trgDec[3] && triggerPt>=100. )                  weight_eS=treePrescl[3];
//
//  std::cout<<"after if statements..."<<std::endl;
//  for(int k=0;k<4;k++){if(trgDec[k]){
//      std::cout<<"trgDec["<<k<<"]="<<trgDec[k]<<std::endl;
//      std::cout<<"treePrescl["<<k<<"]="<<treePrescl[k]<<std::endl;
//      std::cout<<"triggerPt="<<triggerPt<<std::endl;
//      std::cout<<"weight_eS=="<<weight_eS<<std::endl<<std::endl;
//    }}
//
//  std::cout<<"trigComb_sanityCheck exiting..."<<std::endl<<std::endl;
//  return weight_eS;
//}

