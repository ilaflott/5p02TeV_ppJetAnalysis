// custom
#include "readForests_common.h"

//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//// readForests_ppData
//const std::string defDataInFilelist="filelists/5p02TeV_HighPtJet80_forests.txt";
//const std::string defDataOutputName="readForests_ppData_defOut_Jet80_jetTrigEff.root";
const std::string defDataInFilelist="filelists/5p02TeV_HighPtLowerJets_forests.txt";
const std::string defDataOutputName="readForests_ppData_defOut_LowerJets_jetTrigEff.root";
//const std::string defDataInFilelist="filelists/5p02TeV_MinBias6_758p3_forests.txt";
//const std::string defDataOutputName="readForests_ppData_defOut_MinBias_jetTrigEff.root";

const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;
const float defEtaCutLo=0.0, defEtaCutHi=2.0;//really absetacut, don't use negative numbers

int readForests_ppData_jetTrigEff( std::string inFilelist=defDataInFilelist, 
				   int startfile=20, int endfile=21,
				   int radius=defRadius, std::string jetType=defJetType, 
				   bool debugMode=defDebugMode,
				   std::string outfile=defDataOutputName, 
				   float jtEtaCutLo=defEtaCutLo, float jtEtaCutHi=defEtaCutHi      );


const int readForestsArgCount=9+minArgs;

////CUTS: DO NOT CHANGE FOR NOW (6/25/19)
///  5.;  6.;  8.;   10.;   12.;   15.;   18.;   21.;   24.;   28.;  32.;   37.;  43.;   49.;   56.; //other potential low pT cuts, from SMP binning
const float jtPtCut=43.; 
const float jtPtCut_Hi=1410.; 
const float jetQAPtCut=56.;//need a "QA" pt cut to make sure jets at low pT for JEC Unc can migrate to bins up above 56 GeV + so constituent hists for jet ID QA look alright

//for dijet histograms only
const float ldJetPtCut=76., subldJetPtCut=56., ptAveCut=subldJetPtCut, dPhiCut=2./3.*TMath::Pi();//dijet cuts

// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

