// custom
#include "readForests_common.h"
#include "readForests_JEC_v12_ak4PF_75X_data_Run2_v13.h" //these are the JECs 5 TeV pp uses
#include "readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0.h" //another set of slightly older 5 TeV JECs for comparison
#include "readForests_JEC_Fall15_25nsV2_DATA_Uncertainties.h" //uncertainties from 13 TeV, Run Era 2015D (Fall)

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

//const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_forests_local.txt";
const std::string defMCOutputName="readForests_ppMC_defOut_jetPlots.root";
const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;
const float defEtaCutLo=0.0, defEtaCutHi=2.0;//really absetacut, don't use negative numbers

int readForests_ppData_jetPlots( std::string inFilelist=defDataInFilelist, 
				 int startfile=0, int endfile=1,
				 int radius=defRadius, std::string jetType=defJetType, 
				 bool debugMode=defDebugMode,
				 std::string outfile=defDataOutputName, 
				 float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

int readForests_ppMC_jetPlots( std::string inFilelist=defMCInFilelist,
			       int startfile=0, int endfile=4, 
			       int radius=defRadius, std::string jetType=defJetType, 
			       bool debugMode=defDebugMode,			
			       std::string outfile=(defMCOutputName),
			       float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );

const int readForestsArgCount=9+minArgs;

////CUTS: DO NOT CHANGE FOR NOW (6/25/19)
///  5.;  6.;  8.;   10.;   12.;   15.;   18.;   21.;   24.;   28.;  32.;   37.;  43.;   49.;   56.; //other potential low pT cuts, from SMP binning
const float jtPtCut=43.; 
const float jtPtCut_Hi=1410.; 
const float jetQAPtCut=56.;//need a "QA" pt cut to make sure jets at low pT for JEC Unc can migrate to bins up above 56 GeV + so constituent hists for jet ID QA look alright

//for MC jets only
const float genJetPtCut=43.;//49.; 
const float genJetPtCut_Hi=jtPtCut_Hi; 
const float gendrCut=0.1;

//for dijet histograms only
const float ldJetPtCut=76., subldJetPtCut=56., ptAveCut=subldJetPtCut, dPhiCut=2./3.*TMath::Pi();//dijet cuts

// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

