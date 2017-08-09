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
//// custom
//#include "L2L3ResidualWFits.h"
//#include "MCTruthResidual.h"

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
  //jets, 0-3=4 vars
  "jtpt"  , "rawpt" ,
  "jteta" , "jtphi" ,
  //jet constituents, 4-30=27 vars
  "trkN", "trkSum", "trkMax", //trk
  "trkHardN", "trkHardSum",
  "phN", "phSum", "phMax",    //ph
  "phHardN", "phHardSum",
  "chN", "chSum", "chMax",    //ch
  "chHardN", "chHardSum",
  "neN" ,  "neSum" , "neMax" , //ne
  "eN"  ,  "eSum"  , "eMax"  ,    //e
  "muN" ,  "muSum" , "muMax" , //mu
  "neuMult", "chMult", "numConst",
  //"hcalSum", "ecalSum",       //hcal and ecal
  // dijets, 31-35=5 vars
  "Aj" , "xj" , "dphi",    //JERS might not like me
  "leadJetPt", "subleadJetPt"
};
const int N_vars=sizeof(var)/sizeof(std::string);
//"hcalSum", "ecalSum",       //hcal and ecal
// dijets, 31-35=5 vars
//"Aj" , "xj" , "dphi",
//"leadJetPt", "subleadJetPt"
const int jtInd=0, jtConInd=4, dijtInd=31;

//L1
const std::string L1BitStrings[]={//this array is a good idea
  "L1_SingleJet28_BptxAND",
  "L1_SingleJet40_BptxAND",
  "L1_SingleJet48_BptxAND",

  "L1_SingleJet52_BptxAND"
};
//const int N_L1Bits=N_HLTBits;
const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);


// HLT, both Calo and PF available
const std::string Calo_HLTBitStrings[]={
  "HLT_AK4CaloJet40_Eta5p1",
  "HLT_AK4CaloJet60_Eta5p1",
  "HLT_AK4CaloJet80_Eta5p1",

  "HLT_AK4CaloJet100_Eta5p1" 
};
const int N_HLTBits=sizeof(Calo_HLTBitStrings)/sizeof(std::string);

const std::string PF_HLTBitStrings[]={
  "HLT_AK4PFJet40_Eta5p1",
  "HLT_AK4PFJet60_Eta5p1",
  "HLT_AK4PFJet80_Eta5p1",

  "HLT_AK4PFJet100_Eta5p1" 
};


// data tree name array
const std::string dataTreeNames[]={
  "GARBAGE ENTRY",
  "hiEvtAnalyzer/HiTree", "skimanalysis/HltTree", "hltanalysis/HltTree",
  "hltobject/"+Calo_HLTBitStrings[0]+"_v" ,
  "hltobject/"+Calo_HLTBitStrings[1]+"_v" ,
  "hltobject/"+Calo_HLTBitStrings[2]+"_v" ,

  "hltobject/"+Calo_HLTBitStrings[3]+"_v" 
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
  //  "hltobject/"+Calo_HLTBitStrings[3]+"_v" 
};
const int N_MCTrees=sizeof(MCTreeNames)/sizeof(std::string);




// for 5 tev pp jets 2k17
const float absetabins[]={
  0.0, 0.5,
  1.0, 1.5,
  2.0, 2.5,
  3.0, 3.2,
  4.7
};
const int nbins_abseta=sizeof(absetabins)/sizeof(float)-1;


const float etabins[]={
  -5.0 ,
  -4.7 , -4.4 , -4.1 , -3.8 , -3.5 ,
  -3.2 , -3.0 , -2.7 , -2.4 , -2.1 ,
  -1.8 , -1.5 , -1.2 ,
  -0.9 , -0.6 , -0.3 ,
  0.0 ,
  0.3 ,  0.6 ,  0.9 ,
  1.2 ,  1.5 ,  1.8 ,
  2.1 ,  2.4 ,  2.7 ,  3.0 ,  3.2 ,
  3.5 ,  3.8 ,  4.1 ,  4.4 ,  4.7 ,
  5.0
};
const int nbins_eta=sizeof(etabins)/sizeof(float)-1;

//raghavs suggested genpt binning for JER
const float ptbins[]={
  //1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 
  //28,   32,   37, 
  43, 
  49, 
  56, 
  64, 
  74, 
  84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
  507, 548, 592, 638, 686, 737, 790, 846, 905, 967
  //,
  //1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
  //4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000
};
const int nbins_pt=sizeof(ptbins)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one


//// for jet RECO meeting 2k17
//const float absetabins[]={
//  0.00, 0.25, 0.50,  0.75, 
//  1.00, 1.25, 1.50, 1.75,
//  2.10, 3.00, 4.00, 5.00,
//};
//const int nbins_abseta=sizeof(absetabins)/sizeof(float)-1;
//
//const float etabins[]={
//  -5.00, 
//  -4.00, 
//  -3.00, -2.10,
//  -1.75, -1.50, -1.25, -1.00,
//  -0.75, -0.50, -0.25,
//  0.00, 
//  0.25, 0.50,  0.75, 
//  1.00, 1.25, 1.50, 1.75,
//  2.10, 3.00, 
//  4.00, 
//  5.00
//};
//const int nbins_eta=sizeof(etabins)/sizeof(float)-1;
//
//// binning i'm using for jet RECO JER studies 2k17
//const float ptbins[]={
//  30., 35., 40., 45., 50., 55., 65., 70.
//}; 
//const int nbins_pt=sizeof(ptbins)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one

//// WEIGHTS FOR MC
// ------------------------------------------------------------------------------------------------------

//pthat bins and weights for semi-private MC in millibarns
const float pthatbins[]={ 15., 30., 50., 80., 120., 170., 220., 280., 370., 460., 540., 9999. };
const int nbins_pthat=sizeof(pthatbins)/sizeof(int)-1;
const float pthatWeights[]={
  /*(pthat>=15)&&(pthat<30)    , n[0]=932778   , xsDiff=0.49235 mb    , weight= */ 5.27832e-07,
  /*(pthat>=30)&&(pthat<50)    , n[1]=903567   , xsDiff=0.030482 mb   , weight= */ 3.37352e-08,
  /*(pthat>=50)&&(pthat<80)    , n[2]=983531   , xsDiff=0.0035721 mb  , weight= */ 3.63191e-09,
  /*(pthat>=80)&&(pthat<120)   , n[3]=1820782  , xsDiff=0.00042494 mb , weight= */ 2.33383e-10,
  /*(pthat>=120)&&(pthat<170)  , n[4]=1080554  , xsDiff=5.873e-05 mb  , weight= */ 5.43517e-11,
  /*(pthat>=170)&&(pthat<220)  , n[5]=836152   , xsDiff=9.199e-06 mb  , weight= */ 1.10016e-11,
  /*(pthat>=220)&&(pthat<280)  , n[6]=954396   , xsDiff=2.2564e-06 mb , weight= */ 2.36422e-12,
  /*(pthat>=280)&&(pthat<370)  , n[7]=1083994  , xsDiff=6.336e-07 mb  , weight= */ 5.84505e-13,
  /*(pthat>=370)&&(pthat<460)  , n[8]=948240   , xsDiff=1.0884e-07 mb , weight= */ 1.14781e-13,
  /*(pthat>=460)&&(pthat<540)  , n[9]=1558268  , xsDiff=2.215e-08 mb  , weight= */ 1.42145e-14,
  /*(pthat>=540)&&(pthat<9999) , n[10]=2597338 , xsDiff=1.001e-08 mb  , weight= */ 3.85395e-15
};


////pthat bins and weights for official MC in millibarns
//const float pthatbins[]={ 15., 30., 50., 80., 120., 170., 220., 280., 370., 9999. }; 
//const int nbins_pthat=sizeof(pthatbins)/sizeof(int)-1;
//const float pthatWeights[]={
//  /*(pthat>=15)&&(pthat<30)     */   2.67718E-06,
//  /*(pthat>=30)&&(pthat<50)     */   1.59143E-07,
//  /*(pthat>=50)&&(pthat<80)     */   1.69053E-08,
//  /*(pthat>=80)&&(pthat<120)    */   1.95826E-09,
//  /*(pthat>=120)&&(pthat<170)   */   2.67247E-10,
//  /*(pthat>=170)&&(pthat<220)   */   4.29455E-11,
//  /*(pthat>=220)&&(pthat<280)   */   9.59090E-12,
//  /*(pthat>=280)&&(pthat<370)   */   2.32869E-12,
//  /*(pthat>=370)&&(pthat<9999)  */   4.38436E-13 
//};





const bool doVzWeights=true;
// non-event triggered vz weights
const float vzWeights[]={           
0.0579229 ,         // vzLow=-24
0.0622996 ,
0.061942  ,
0.0902018 ,
0.0744494 ,
0.0785317 ,                // vzLow=-21.5
0.088784  ,
0.106222  ,
0.113585  ,
0.142753  ,
0.147792  ,               // vzLow=-19
0.155734  ,
0.192015  ,
0.216288  ,
0.245479  ,
0.283993  ,               // vzLow=-16.5
0.31437	  ,
0.352618  ,
0.390549  ,
0.43129	  ,
0.468473  ,               // vzLow=-14
0.501128  ,
0.539377  ,
0.582421  ,
0.647345  ,
0.643217  ,               // vzLow=-11.5
0.707504  ,
0.73216	  ,
0.768601  ,
0.819635  ,
0.836843  ,               // vzLow=-9
0.874235  ,
0.89897	  ,
0.917435  ,
0.948561  ,
0.983956  ,               // vzLow=-6.5
1.02518	  ,
1.03036	  ,
1.039	  ,
1.06393	  ,
1.0824 	  ,            // vzLow=-4
1.09788	  ,
1.1064	  ,
1.12646	  ,
1.13744	  ,
1.14794   ,              // vzLow=-1.5
1.14349	  ,
1.16373	  ,
1.16068	  ,
1.16672	  ,
1.16045   ,              // vzLow=1
1.17421	  ,
1.15623	  ,
1.1739	  ,
1.14798	  ,
1.14334   ,             // vzLow=3.5
1.13985	  ,
1.12823	  ,
1.11934	  ,
1.11583	  ,
1.09352   ,              // vzLow=6
1.05137	  ,
1.04333	  ,
1.01349	  ,
0.986353  ,
0.959904  ,              // vzLow=8.5
0.929557  ,
0.913313  ,
0.88634	  ,
0.827843  ,
0.807459  ,              // vzLow=11
0.743207  ,
0.727486  ,
0.677889  ,
0.644392  ,
0.594507  ,               // vzLow=13.5
0.56391	  ,
0.491046  ,
0.465995  ,
0.429629  ,
0.390355  ,               // vzLow=16
0.337874  ,
0.31006	  ,
0.274059  ,
0.238072  ,
0.217408  ,               // vzLow=18.5
0.189755  ,
0.169915  ,
0.164071  ,
0.137212  ,
0.136876  ,               // vzLow=21
0.0945742 ,
0.101419  ,
0.0768486 ,
0.075515  ,
0.0863076                 // vzLow=23.5
};
const int nbins_vzWeights=sizeof(vzWeights)/sizeof(float);//should be 96
const float minbinValue_vzWeights=-24., maxbinValue_vzWeights=24.;
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


//const float rapbins[]={
//  0.0, 0.5,
//  1.0, 1.5,
//  2.0, 2.5,
//  3.0, 3.2,
//  4.7
//};
//const int nbins_rap=sizeof(rapbins)/sizeof(float)-1;











//const float JEC_ptbins[]={
//  17., 22., 27.,    //15-30
//  33., 39., 47.,    //30-50
//  55., 64., 74.,    //50-80
//  84., 97., 114.,   //80-120
//  133., 153.,      //120-170
//  174., 196.,      //170-220
//  220., 245., 272., //220-300
//  300., 350., 400., //300-500
//  550., 790., 1000. //500-inf
//};
//const int nbins_JEC_ptbins=sizeof(JEC_ptbins)/sizeof(float)-1;
//
//const float ptbins2[]={
//  //3., 4., 
//  5., 7., 9., 12.,
//  15., 18., 21., 24., 28.,
//  32., 37., 43., 49., 
//  56.,
//  64., 74., 84., 97., 114.,
//  133., 153., 174., 196.,
//  220., 245., 272., 300.,
//  330., 362., 395., 430.,
//  468., 507., 548., 592.,
//  638., 686., 1000.//, 1500
//}; //raghavs suggested genpt binning for JER
//const int nbins_pt2=sizeof(ptbins2)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one








//const float etabins[]={
//  -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013,
//  -3.839, -3.664, -3.489, -3.314, -3.139,
//  -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -2.043,
//  -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044,
//  -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174,
//  -0.087, +0.000, +0.087,
//  +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, +0.879, +0.957,
//  +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, +1.653, +1.740, +1.830, +1.930,
//  +2.043, +2.172, +2.322, +2.500, +2.650, +2.853, +2.964,
//  +3.139, +3.314, +3.489, +3.664, +3.839,
//  +4.013, +4.191, +4.363, +4.538, +4.716, +4.889, +5.191
//};
//const int nbins_eta=sizeof(etabins)/sizeof(float)-1;


// original genpt JER binning i inherited from raghav
//const float ptbins[]={
//  15., 30., 50., 80.,
//  120., 170., 220., 300.,
//  500.
//}; // original genpt JER binning i inherited from raghav
//const int nbins_pt=sizeof(ptbins)/sizeof(float)-1;//above values define edges of bins, not centers, so subtract one























//semiofficial MC/{no trig, HLTAK4CaloJets,HLTAK4PFJets}, official MC/HLTak4CaloJets
  //1.0, //vz=-24cm
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //0.381083, /*  0.35749,   //for i=1 ,   1.8995   ,   	 vz=-15cm */
  //0.420849, /*  0.400545,  //for i=2 ,   1.89802  ,   	 */
  //0.457166, /*  0.433741,  //for i=3 ,   1.75632  ,   	 */
  //0.489082, /*  0.469437,  //for i=4 ,   1.72277  ,   	 */
  //0.526434, /*  0.512839,  //for i=5 ,   1.70519  ,   	 */
  //0.568316, /*  0.556416,  //for i=6 ,   1.68102  ,   	 */
  //0.631824, /*  0.612334,  //for i=7 ,   1.60677  ,   	 */
  //0.627796, /*  0.622496,  //for i=8 ,   1.46807  ,   	 */
  //0.690459, /*  0.672646,  //for i=9 ,   1.46982  ,   	 */
  //0.714526, /*  0.706209,  //for i=10,   1.42762  ,   	 */
  //0.750061, /*  0.740972,  //for i=11,   1.38758  ,   	 */
  //0.799916, /*  0.792563,  //for i=12,   1.37528  ,   	 */
  //0.816735, /*  0.815354,  //for i=13,   1.33446  ,   	 */
  //0.85323 , /*  0.853215,  //for i=14,   1.29506  ,   	 */
  //0.877455, /*  0.875179,  //for i=15,   1.27491  ,   	 */
  //0.895425, /*  0.891705,  //for i=16,   1.18277  ,   	 */
  //0.92577 , /*  0.921846,  //for i=17,   1.17951  ,   	 */
  //0.960376, /*  0.958741,  //for i=18,   1.18097  ,   	 */
  //1.00063 , /*  1.00373,   //for i=19,   1.1312   ,   	 */
  //1.00567 , /*  1.00595,   //for i=20,   1.08347  ,   	 */
  //1.01408 , /*  1.01643,   //for i=21,   1.05749  ,   	 */
  //1.03837 , /*  1.04327,   //for i=22,   1.05493  ,   	 */
  //1.05649 , /*  1.06326,   //for i=23,   1.00842  ,   	 */
  //1.07153 , /*  1.07991,   //for i=24,   1.01667  ,   	 */
  //1.07986 , /*  1.08791,   //for i=25,   0.968532 ,   	 */
  //1.09942 , /*  1.10771,   //for i=26,   0.951276 ,   	 */
  //1.1102  , /*  1.11622,   //for i=27,   0.943379 ,   	 */
  //1.12039 , /*  1.12078,   //for i=28,   0.941838 ,   	 */
  //1.116   , /*  1.12084,   //for i=29,   0.919547 ,   	 */
  //1.1358  , /*  1.14016,   //for i=30,   0.921949 ,   	 */
  //1.13279 , /*  1.14058,   //for i=31,   0.907791 ,   	 */
  //1.1387  , /*  1.14684,   //for i=32,   0.898317 ,   	 */
  //1.13261 , /*  1.14148,   //for i=33,   0.900794 ,   	 */
  //1.14604 , /*  1.15419,   //for i=34,   0.897549 ,   	 */
  //1.12861 , /*  1.1392,    //for i=35,   0.911964 ,   	 */
  //1.14573 , /*  1.15569,   //for i=36,   0.898306 ,   	 */
  //1.12039 , /*  1.13269,   //for i=37,   0.902429 ,   	 */
  //1.11589 , /*  1.12478,   //for i=38,   0.899196 ,   	 */
  //1.11258 , /*  1.11415,   //for i=39,   0.899973 ,   	 */
  //1.10108 , /*  1.10534,   //for i=40,   0.930686 ,   	 */
  //1.09245 , /*  1.09492,   //for i=41,   0.937519 ,   	 */
  //1.08902 , /*  1.08401,   //for i=42,   0.931296 ,   	 */
  //1.06721 , /*  1.06167,   //for i=43,   0.962955 ,   	 */
  //1.0261  , /*  1.02492,   //for i=44,   0.978917 ,   	 */
  //1.01827 , /*  1.01539,   //for i=45,   0.989318 ,   	 */
  //0.989024, /*  0.985653,  //for i=46,   0.99599  ,   	 */
  //0.962674, /*  0.952256,  //for i=47,   1.00077  ,   	 */
  //0.936773, /*  0.930472,  //for i=48,   1.02948  ,   	 */
  //0.907251, /*  0.894218,  //for i=49,   0.983398 ,   	 */
  //0.891311, /*  0.875605,  //for i=50,   1.06389  ,   	 */
  //0.864917, /*  0.847596,  //for i=51,   1.08784  ,   	 */
  //0.807905, /*  0.78767,   //for i=52,   1.11999  ,   	 */
  //0.787968, /*  0.768839,  //for i=53,   1.09857  ,   	 */
  //0.725249, /*  0.700327,  //for i=54,   1.14199  ,   	 */
  //0.709917, /*  0.687202,  //for i=55,   1.12368  ,   	 */
  //0.661411, /*  0.629434,  //for i=56,   1.2143   ,   	 */
  //0.628832, /*  0.597895,  //for i=57,   1.28173  ,   	 */
  //0.580118, /*  0.555087,  //for i=58,   1.16107  ,   	 */
  //0.550339, /*  0.523324,  //for i=59,   1.07179  ,   	 */
  //0.479087,  /*  0.448861   //for i=60,   1.28413      	 vz=15cm */
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0,
  //1.0//, vz=25cm
