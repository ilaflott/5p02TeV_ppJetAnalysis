// C++, C, etc.  
#include <cstdlib>  
#include <cstdio>   
#include <cmath>   


// gen variable names for Th2 profiles//MCEff Ploits
const std::string genVars[] = {   
  "pt","eta","phi" //,"drjt"
};
const int N_genVars = sizeof(genVars)/sizeof(std::string);

const std::string genVars_ptrat[] = {   
  "eta","phi" ,"drjt"
};
const int N_genVars_ptrat = sizeof(genVars_ptrat)/sizeof(std::string);


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

// binning arrays
//const int ptbins[] = { 15, 30, 50, 80, 120, 170, 220, 300, 500 };
const int ptbins[]={
     3, 4, 5, 7, 9, 12, 
     15, 18, 21, 24, 28,
     32, 37, 43, 49, 56,
     64, 74, 84, 97, 114,
     133, 153, 174, 196,
     220, 245, 272, 300, 
     330, 362, 395, 430,
     468, 507, 548, 592,
     638, 686, 1000//, 1500
}; //raghavs suggested genpt binning for JER
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

const double xmin=ptbins[0];
const double xmax=ptbins[nbins_pt];

//// tree names+directories
//const std::string treeNames_ppMC[]={ 
//  "GARBAGE ENTRY" , //akPfJetAna of arb. radius, form strin
//  "GARBAGE ENTRY" , //akCaloJet of same radius
//  "hiEvtAnalyzer/HiTree" ,
//  "skimanalysis/HltTree" ,
//  "hltanalysis/HltTree" //,
//  //  "hltobject/"+HLTBitStrings[0]+"_v" , 
//  //  "hltobject/"+HLTBitStrings[1]+"_v" , 
//  //  "hltobject/"+HLTBitStrings[2]+"_v" , 
//  //  "hltobject/"+HLTBitStrings[3]+"_v"   
//}; 
//const int Ntrees_ppMC = sizeof(treeNames)/sizeof(std::string);

