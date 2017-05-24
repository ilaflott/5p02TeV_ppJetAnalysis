////////// (initializa/declara)tions //////////
// custom
#include "readForests_common.h"
#include "L2L3ResidualWFits.h"
//#include "MCTruthResidual.h"


//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1+7;//jtetalo/hi optional, default is 0to4.7

//const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
//const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_forests_local.txt";
const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_Official_forests_acrossBins.txt";
const std::string defMCOutputName="readForests_ppMC_JERS.root";

int readForests_ppMC_JERS( std::string inFilelist=defMCInFilelist, int startfile=0, int endfile=20,
			   int radius=4, std::string jetType="PF",
			   bool debugMode=true, 
			   std::string outfile=defMCOutputName,
			   float jtEtaCutLo=0., float jtEtaCutHi=4.7      );

//const int readForestsArgCount=7+minArgs;
const float jtPtCut=56.;
const float genJetPtCut=5.;

const int genptBin1Low=56, genptBin1High=64;	 
const int genptBin2Low=196, genptBin2High=220;


// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

