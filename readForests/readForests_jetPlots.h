// custom
#include "readForests_common.h"
#include "L2L3ResidualWFits.h"


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
//const std::string defDataOutputName="readForests_ppData_defOut_LowerJets_jetPlots_trgPtL2L3Res.root";
//const std::string defDataOutputName="readForests_ppData_defOut_LowerJets_jetPlots_trgPtNoRes.root";

//const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_forests_local.txt";
const int defStartFile=5;
const int defEndFile=6;
const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;//, fastDebugMode = true;
const std::string defMCOutputName="readForests_ppMC_defOut_jetPlots.root";
const float defEtaCutLo=0.0, defEtaCutHi=2.0;//really absetacut
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
			       int startfile=0, int endfile=10, 
			       int radius=defRadius, std::string jetType=defJetType, 
			       bool debugMode=defDebugMode,			
			       std::string outfile=(defMCOutputName),
			       float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

const int readForestsArgCount=9+minArgs;

// extended eta range for jetID Eff, or more QA in diff region... etc.


//const float jtPtCut=30; // 49 or 56 or 64 or 74...
const float jtPtCut=56.; // 49 or 56 or 64 or 74...
const float jtPtCut_Hi=1410.; 

const float jetQAPtCut=jtPtCut;

const float genJetPtCut=49.; 
const float genJetPtCut_Hi=jtPtCut_Hi; 


const float ldJetPtCut=76., subldJetPtCut=56., ptAveCut=subldJetPtCut, dPhiCut=2./3.*TMath::Pi();//dijet cuts







// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

