// C++, C, etc.  
#include <cstdlib>  
#include <cstdio>   
#include <cmath>   


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

