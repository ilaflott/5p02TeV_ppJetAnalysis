#include <cstdlib>
#include <string>
#include <TMath.h>

//#define NOBJECT_MAX 16384

//highPtJetTrig 5.02TeV 2015 prompt reco

//const double intgrtdLumi=25.8*pow(10.,6.);// 25.8 pb^-1 to \microbarns^-1 //old lumi
const double intgrtdLumi=27.4*pow(10.,3.);// pb^-1 to \microbarns^-1 // new lumi
const double intgrtdLumi_unc=2.3/100.;
const int jetTrigQABinning=10;






//// momentum bins/boundaries
//const double pthat[]={ 15, 30, 50, 80, 120, 170, 
//		       220, 280, 370, 460, 540, 2000 };
////const double pthat[]={ 15, 30, 50, 80, 120, 170, 
////		       220, 280, 370, 9999 };
//const int Nbins_pthat=sizeof(pthat)/sizeof(double);

// giving a TH1F an array for binning is fussy; demands root-type numbers or else it gives 
// Error in <TAxis::TAxis::Set>: bins must be in increasing order
const Double_t ptbins_debug[]={  
//  5.,
//  6.,
//  8.,
//  10.,
//  12.,
//  15.,
//  18.,
//  21.,
//  24.,
//  28.,
//  32.,
//  37.,  //40., //this junk bin added in for cosmetics
//  43., 
//  49., //junk bins from here an above
  56.,   64.,   74.,   84., 97., 114., 
  133., 153., 174., 196., 220., 
  245., 272., 300., 
  330., 362., 395., 430., 
  468., 507., 
  548., 592., 638., 
  686., 
  737., 
  790., 
  846., 
  905., 
  967., //, // end of SMP Official bin limits 
  1032.//1000. //,  // to catch overflow or other entries?
  //1050. // junk bin, cosmetics only
}; 
const int nbins_pt_debug = sizeof(ptbins_debug)/sizeof(Double_t)-1;


// eta bins/boundarie
const Double_t absetabins[]={
  0.0, 
  0.5,  
  1.0, 
  1.5,  
  2.0, 
  2.5,  
  3.0//, 3.2,  
  //  4.7
};
std::string absetabins_str[]={
  "0.0", 
  "0.5",  
  "1.0", 
  "1.5",  
  "2.0", 
  "2.5",  
  "3.0"//, 3.2,  
  //  4.7
};
const int nbins_abseta = sizeof(absetabins)/sizeof(Double_t)-1;






const int pthatbins[]={ 15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 9999 };
//const int pthatbins[]={ 15, 30, 50, 80, 120, 170, 220, 280, 370, 9999};
const int nbins_pthat=sizeof(pthatbins)/sizeof(int)-1;

const double pthatWeights[]={//from ppMC_Py8_CUETP8M1_QCDjetAllPtBins
  /*(pthat>=15)&&(pthat<30),   */ 5.27832e-07,
  /*(pthat>=30)&&(pthat<50),   */ 3.37352e-08,
  /*(pthat>=50)&&(pthat<80),   */ 3.63191e-09,
  /*(pthat>=80)&&(pthat<120),  */ 2.33383e-10,
  /*(pthat>=120)&&(pthat<170), */ 5.43517e-11,
  /*(pthat>=170)&&(pthat<220), */ 1.10016e-11,
  /*(pthat>=220)&&(pthat<280), */ 2.36422e-12,
  /*(pthat>=280)&&(pthat<370), */ 5.84505e-13,
  /*(pthat>=370)&&(pthat<460), */ 1.14781e-13,
  /*(pthat>=460)&&(pthat<540), */ 1.42145e-14,
  /*(pthat>=540)&&(pthat<9999),*/ 3.85395e-15
};

//const double pthatWeights[]={//from ppMC_Py8_CUETP8M1_QCDjetAllPtBins_Official
//    /*(pthat>=15)&&(pthat<30)     */   2.67718E-06,
//    /*(pthat>=30)&&(pthat<50)     */   1.59143E-07,
//    /*(pthat>=50)&&(pthat<80)     */   1.69053E-08,
//    /*(pthat>=80)&&(pthat<120)    */   1.95826E-09,
//    /*(pthat>=120)&&(pthat<170)   */   2.67247E-10,
//    /*(pthat>=170)&&(pthat<220)   */   4.29455E-11,
//    /*(pthat>=220)&&(pthat<280)   */   9.59090E-12,
//    /*(pthat>=280)&&(pthat<370)   */   2.32869E-12,
//    /*(pthat>=370)&&(pthat<9999)  */   4.38436E-13
//};

// binning arrays

const int ptbins[]={
  //  43,
  //  49,
  56,
  64,
  74,
  84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
  507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
  1032,
  //1000   //SMP "official", last bin junk
}; 
// 
//  5,
//  15, 18, 21, 24, 28,
//  32, 37, 43, 49, 56,
//  64, 74, 84, 97, 114,
//  133, 153, 174, 196,
//  220, 245, 272, 300, 
//  330, 362, 395, 430,
//  468, 507, 548, 592,
//  638, 686, 1000//, 1500
//}; //my tweaked version of raghav's original

//  7, 9, 12, 
//  15, 18, 21, 24, 28,
//  32, 37, 43, 49, 56,
//  64, 74, 84, 97, 114,
//  133, 153, 174, 196,
//  220, 245, 272, 300, 
//  330, 362, 395, 430,
//  468, 507, 548, 592,
//  638, 686, 1000//, 1500
//}; //raghavs original suggested genpt binning for JER

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
//my tweaked eta binning
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
//};//original eta binning
const int nbins_eta = sizeof(etabins)/sizeof(double)-1;

//const double xmin=ptbins[0];
//const double xmax=ptbins[nbins_pt];

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
