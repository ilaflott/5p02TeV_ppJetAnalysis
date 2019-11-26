// custom
#include "readForests_common.h"
#include "L2L3ResidualWFits.h"
#include "readForests_JEC_v12_ak4PF_75X_data_Run2_v13.h" //these are the JECs 5 TeV pp uses
#include "readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0.h" //another set of slightly older 5 TeV JECs for comparison
#include "readForests_JEC_Fall15_25nsV2_DATA_Uncertainties.h" //uncertainties from 13 TeV, Run Era 2015D (Fall)
#include "johns_log_run262163_evt_lumisToRead.h"


//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//// readForests_ppData
//const std::string defDataJ80InFilelist="filelists/test_readForests_ppData_Jet80_local.txt";
const std::string defDataJ80InFilelist_forjohncomp="filelists/defDataJet80InFilelist_forjohncomp.txt";
//const std::string defDataLOJInFilelist="filelists/test_readForests_ppData_LowerJets_local.txt";
//const std::string defDataInFilelist="filelists/5p02TeV_HighPtLowerJets_forests.txt";
const std::string defDataInFilelist="filelists/5p02TeV_HighPtJet80_forests.txt";

//const std::string defDataInFilelist="filelists/5p02TeV_HighPtJet80_forests.txt";
//const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_forests_local.txt";
const int defStartFile=0;
const int defEndFile=99999;
const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;//, fastDebugMode = true;
const std::string defDataOutputName="readForests_ppData_findEvt_defOut.root";
const std::string defMCOutputName="readForests_ppMC_defOut";//.root";
const float defEtaCutLo=0.0, defEtaCutHi=5.1;//really absetacut

int readForests_ppData_findEvt( std::string inFilelist=defDataInFilelist, 
				int startfile=39, int endfile=39,
				int radius=defRadius, std::string jetType=defJetType, 
				bool debugMode=defDebugMode,
				std::string outfile=defDataOutputName, 
				float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

int readForests_ppData_findEvt_v2( std::string inFilelist=defDataInFilelist, 
				   int startfile=309, int endfile=309,
				   int radius=defRadius, std::string jetType=defJetType, 
				   bool debugMode=defDebugMode,
				   std::string outfile=defDataOutputName, 
				   float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );


int readForests_ppData_findEvt_v3( std::string inFilelist=defDataJ80InFilelist_forjohncomp, 
				   int startfile=0, int endfile=5,
				   int radius=defRadius, std::string jetType=defJetType, 
				   bool debugMode=defDebugMode,
				   std::string outfile=defDataOutputName, 
				   float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

//int readForests_ppMC_findEvt( std::string inFilelist=defMCInFilelist,
//			       int startfile=0, int endfile=4, 
//			       int radius=defRadius, std::string jetType=defJetType, 
//			       bool debugMode=defDebugMode,			
//			       std::string outfile=(defMCOutputName+"_findEvt.root"),
//			       float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

const int readForestsArgCount=9+minArgs;

// extended eta range for jetID Eff, or more QA in diff region... etc.
const float jtPtCut=10.; 
const float jtPtCut_Hi=1400.; 
const float jetQAPtCut=1400.;//need a "QA" pt cut to make sure jets at low pT for JEC Unc can migrate to bins up above 56 GeV + so constituent hists for jet ID QA look alright

//for MC jets only
const float genJetPtCut=43.;//49.; 
const float genJetPtCut_Hi=jtPtCut_Hi; 
const float gendrCut=0.1;

const float ldJetPtCut=74., subldJetPtCut=56., ptAveCut=64., dPhiCut=2./3.*TMath::Pi();//dijet cuts

// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

