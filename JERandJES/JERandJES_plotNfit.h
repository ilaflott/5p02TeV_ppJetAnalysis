// C/C++
#include <cstdlib>
#include <cstring>
#include <iostream>

// file IO/settings
#include <TROOT.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TFile.h>

// fit + boundaries + plots
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

// draw/print
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLine.h>
#include <TLatex.h>

// methods, functions, numbers, etc.
//#include "readFiles_ppMC_JERandJES_binsNconsts.h"
//#include "readFiles_ppMC_JERandJES_fitsNfuncs.h"

// run settings
const bool printDebug=true;
const bool doBjets = false; //true;
const int iSave=1;
const int wJetID=1;  //int id = wJetID;
const int digi = 3;
const double fracRMS = 1.00;


// files/IO/strings
const std::string CMSSW_BASE=
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/";

const std::string inFile_MC_dir=
  "readFiles/ppMC/saved_outputCondor/readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_2016-07-09/";
const std::string inFile_MC_name="QCDjetAllPtBins_ak4PF-allFiles.root";
const std::string inFile_MC=CMSSW_BASE+inFile_MC_dir+inFile_MC_name;

const std::string inFile_Data_dir=
  "readFiles/ppData/saved_outputCondor/readFiles_ppData_5p02TeV_HighPtJetTrig_2016-06-10_allFiles/";
const std::string inFile_Data_name="HighPtJetTrig_ak4PF-allFiles.root";
const std::string inFile_Data=CMSSW_BASE+inFile_Data_dir+inFile_Data_name;



