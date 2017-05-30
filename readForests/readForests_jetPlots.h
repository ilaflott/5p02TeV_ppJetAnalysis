// custom
#include "readForests_common.h"
#include "L2L3ResidualWFits.h"


//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//// readForests_ppData
//const std::string defDataInFilelist="filelists/test_readForests_ppData_Jet80_local.txt";
const std::string defDataInFilelist="filelists/test_readForests_ppData_LowerJets_local.txt";
//const std::string defDataInFilelist="filelists/5p02TeV_HighPtJet80_forests.txt";
//const std::string defDataInFilelist="filelists/5p02TeV_HighPtLowerJets_forests.txt";
const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
//const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_forests_local.txt";
const int defStartFile=0;
const int defEndFile=1;
const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;//, fastDebugMode = true;
const std::string defDataOutputName="readForests_ppData_defOut.root";
const std::string defMCOutputName="readForests_ppMC_defOut";//.root";
const float defEtaCutLo=0.0, defEtaCutHi=4.7;//really absetacut

int readForests_ppData_jetPlots( std::string inFilelist=defDataInFilelist, 
				 int startfile=0, int endfile=1,
				 int radius=defRadius, std::string jetType=defJetType, 
				 bool debugMode=defDebugMode,
				 std::string outfile=defDataOutputName, 
				 float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

int readForests_ppMC_jetPlots( std::string inFilelist=defMCInFilelist,
			       int startfile=0, int endfile=10, 
			       int radius=defRadius, std::string jetType=defJetType, 
			       bool debugMode=defDebugMode,			
			       std::string outfile=(defMCOutputName+"_jetPlots.root"),
			       float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

const int readForestsArgCount=9+minArgs;

// extended eta range for jetID Eff, or more QA in diff region... etc.
const float jtPtCut=56.;
const float jetQAPtCut=jtPtCut;
const float genJetPtCut=30.;

const float ldJetPtCut=80., subldJetPtCut=56., ptAveCut=68., dPhiCut=2./3.*TMath::Pi();//dijet cuts







// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe
