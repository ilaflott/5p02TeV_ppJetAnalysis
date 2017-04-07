////////// (initializa/declara)tions //////////
// custom
#include "readForests_common.h"
#include "L2L3ResidualWFits.h"
//#include "MCTruthResidual.h"


//// FUNCTIONS
// ---------------------------------------------------------------------------------------------------------------
// define "main" functions, their default inputs, number of input arguments in this section

//// main
const int minArgs=1;

//const std::string defMCInFilelist="filelists/5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_forests.txt";
//const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_forests_local.txt";
//const defEndile=10;
const std::string defMCInFilelist="filelists/test_readForests_ppMC_Py8_CUETP8M1_Official_forests_acrossBins.txt";
const int defEndfile=55;
const int defRadius=4;
const std::string defJetType="PF";
const bool defDebugMode=true;//, fastDebugMode = true;
const std::string defMCOutputName="readForests_ppMC_";//.root";

int readForests_ppMC_JERS( std::string inFilelist=defMCInFilelist,
			   int startfile=0, int endfile=defEndfile,
			   int radius=defRadius, std::string jetType=defJetType,
			   bool debugMode=defDebugMode,
			   std::string outfile=(defMCOutputName+"JERS.root"),
			   float jtEtaCutLo=0., float jtEtaCutHi=4.7      );

const int readForestsArgCount=7+minArgs;

const float jtPtCut=15.;
const float jtGenPtCut=5.;
//const float jtEtaCutLo=0.0, jtEtaCutHi=4.7;














// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ readForests.C $(root-config --cflags --libs) -Werror -Wall -O2 -o readForests.exe

