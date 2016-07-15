// this header for code(s):
// norm_ppMCResponseMatrix.C, ppMC_JER_plotUncAndResiduals.C 
// based on pp_JER_plot.C and pp_normResponseMatrix.C writen by
// Raghav Kunnawalkam Elayavalli June 5th 2014
// tweaks+rearrangements: Ian Laflotte July 15 2016

// C++, C, etc.
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <utility>

// ROOTSYS
#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TStyle.h>

// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>

// utilities
#include <TEventList.h>
#include <TCut.h>
#include <TMath.h>
#include <TError.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSpline.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>

// histos
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
//#include <TH1D.h>
//#include <TH2D.h>

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

const double ptbins_ana[] = {
  50, 60, 70, 80, 90,
  100, 110, 130, 150, 170, 190,
  210, 240, 270, 300};
const int ptbins[] = {
  50, 60, 70, 80, 90,
  100, 110, 130, 150, 170, 190,
  210, 240, 270, 300};
const int nbins_ana=(sizeof(ptbins_ana)/sizeof(double))-1;
//const int nbins_ana = 14;

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}
