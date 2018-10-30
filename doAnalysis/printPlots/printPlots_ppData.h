// C++, C, etc.  
#include <cstdlib>  
#include <cstdio>   
#include <cmath>   


//highPtJetTrig 5.02TeV 2015 prompt reco

//const double intgrtdLumi=25.8*pow(10.,6.);// 25.8 pb^-1 to \microbarns^-1 //old lumi
const double intgrtdLumi=27.4*pow(10.,3.);// pb^-1 to \microbarns^-1 // new lumi
const double intgrtdLumi_unc=2.3/100.;

// string arrays
//L1
const std::string L1BitStrings[]={
  "L1_SingleJet28_BptxAND",
  "L1_SingleJet40_BptxAND",
  "L1_SingleJet48_BptxAND",
  "L1_SingleJet52_BptxAND"
};
const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);

//HLT
const std::string HLTBitStrings[N_L1Bits]={
  "HLT_AK4CaloJet40_Eta5p1",
  "HLT_AK4CaloJet60_Eta5p1",
  "HLT_AK4CaloJet80_Eta5p1",
  "HLT_AK4CaloJet100_Eta5p1"    
};
const int N_HLTBits=sizeof(HLTBitStrings)/sizeof(std::string);

//// tree names+directories
//const std::string ppDataTrees[]={ 
//  "GARBAGE ENTRY" , //use jet ana of choice later
//  "hiEvtAnalyzer/HiTree" ,
//  "skimanalysis/HltTree" ,
//  "hltanalysis/HltTree" ,
//  "hltobject/"+HLTBitStrings[0]+"_v" , 
//  "hltobject/"+HLTBitStrings[1]+"_v" , 
//  "hltobject/"+HLTBitStrings[2]+"_v" , 
//  "hltobject/"+HLTBitStrings[3]+"_v"   
//}; 
//const int N_trees_ppData = sizeof(treeNames)/sizeof(std::string);


// trigs
const std::string HLTName[] = {
  "HLTComb","HLT40","HLT60","HLT80","HLT100"
};

// trigs
const std::string HLTCaloName_Leg[] = {
  "HLTak4CaloJetsCombo","HLTak4CaloJets40","HLTak4CaloJets60","HLTak4CaloJets80","HLTak4CaloJets100"
};
// trigs
const std::string HLTPFName_Leg[] = {
  "HLTPFJets Combo","HLTak4PFJets40","HLTak4PFJets60","HLTak4PFJets80","HLTak4PFJets100"
};
const int N_trigs=sizeof(HLTName)/sizeof(std::string);

//// combos
const std::string etaWidth = "20_eta_20";

const int jetTrigQABinning=10;





//// other useful things
//const std::string HLThNames[]={"HLT40_","HLT60_","HLT80_","HLT100_","HLTComb_","TrgCombTest_"};
//const std::string HLThTitle[]={"HLT40 ","HLT60 ","HLT80 ","HLT100 ","HLTComb ","Kurt Meth. HLTComb "};
//const int N_HLThNames=sizeof(HLThNames)/sizeof(std::string);

