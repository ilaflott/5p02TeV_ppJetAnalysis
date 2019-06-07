// custom
#include "readForests_common.h"

//#include "L2L3ResidualWFits.h"
//#include "L2L3ResidualWFits_2015.h"
//#include "MCTruthResidual.h"
//#include "readForests_JEC_twiki.h"
//#include "readForests_JEC_v7_ak4PF_75X_dataRun2_v13.h"
#include "readForests_JEC_v12_ak4PF_75X_data_Run2_v13.h"
#include "readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0.h"
#include "readForests_JEC_Fall15_25nsV2_DATA_Uncertainties.h"

//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//// readForests_ppData
//const std::string defDataInFilelist="filelists/5p02TeV_HighPtJet80_forests.txt";
//const std::string defDataOutputName="readForests_ppData_defOut_Jet80_jetPlots.root";
const std::string defDataInFilelist="filelists/5p02TeV_HighPtLowerJets_forests.txt";
const std::string defDataOutputName="readForests_ppData_defOut_LowerJets_jetPlots.root";
//const std::string defDataInFilelist="filelists/5p02TeV_MinBias6_758p3_forests.txt";
//const std::string defDataOutputName="readForests_ppData_defOut_MinBias_jetPlots.root";

//const std::string defDataOutputName="readForests_ppData_defOut_LowerJets_jetPlots_trgPtL2L3Res.root";
//const std::string defDataOutputName="readForests_ppData_defOut_LowerJets_jetPlots_trgPtNoRes.root";

//const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_forests_local.txt";
const std::string defMCOutputName="readForests_ppMC_defOut_jetPlots.root";
const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;//, fastDebugMode = true;
const float defEtaCutLo=0.0, defEtaCutHi=2.0;//really absetacut
//const float defEtaCutLo=0.0, defEtaCutHi=5.191;//really absetacut
//const float defEtaCutLo=0.0, defEtaCutHi=2.4;//really absetacut
//const float defEtaCutLo=2.4, defEtaCutHi=2.7;//really absetacut
//const float defEtaCutLo=2.7, defEtaCutHi=3.0;//really absetacut
//const float defEtaCutLo=3.2, defEtaCutHi=4.7;//really absetacut

int readForests_ppData_jetPlots( std::string inFilelist=defDataInFilelist, 
				 int startfile=0, int endfile=9,
				 int radius=defRadius, std::string jetType=defJetType, 
				 bool debugMode=defDebugMode,
				 std::string outfile=defDataOutputName, 
				 float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

int readForests_ppMC_jetPlots( std::string inFilelist=defMCInFilelist,
			       int startfile=0, int endfile=40, 
			       int radius=defRadius, std::string jetType=defJetType, 
			       bool debugMode=defDebugMode,			
			       std::string outfile=(defMCOutputName),
			       float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

const int readForestsArgCount=9+minArgs;

// extended eta range for jetID Eff, or more QA in diff region... etc.



//const float jtPtCut=30; // 49 or 56 or 64 or 74...
const float jtPtCut=43.; // 49 or 56 or 64 or 74...
//const float jtPtCut=5.; 
//const float jtPtCut=18; 
//  5.;  6.;  8.;   10.;   12.;   15.;   18.;   21.;   24.;   28.;  32.;   37.;  43.;   49.;   56.; //other low pT cuts
const float jtPtCut_Hi=1410.; 
const float jetQAPtCut=jtPtCut;
const float genJetPtCut=12.;//49.; 
const float genJetPtCut_Hi=jtPtCut_Hi; 
const float gendrCut=0.1;
const float ldJetPtCut=76., subldJetPtCut=56., ptAveCut=subldJetPtCut, dPhiCut=2./3.*TMath::Pi();//dijet cuts







// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

