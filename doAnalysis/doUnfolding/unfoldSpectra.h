// unfoldSpectra.h written by Ian W. Laflotte
// June 28 2016
// heavily based on unfold.C and utilities.h by Raghav Kunnawalkam Elayavalli
// provides functions and binning common to doMCClosureTests and unfoldDataSpectra
// compile doMCClosureTests or unfoldDataSpectra with...
// g++ -L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold -I"${ROOUNFOLD_PATH}/src" <target>.C $(root-config --cflags --libs) -Werror -Wall -O2 -o <target>.exe
// OR, if rooUnfoldCompile function is set in .bashrc
// rooUnfoldCompile <target>.C

// C++, C, etc. 
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <vector>

// ROOTSYS
#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>

// I/O
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
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
#include <TMatrix.h>
#include <TMatrixT.h>
#include <TMatrixD.h>
#include <TSVDUnfold.h>

// plotting
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

// histos
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>

// RooUnfold
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"

// custom
#include "unfoldSpectra_TH1_funcs.h"
#include "unfoldSpectra_TH2_funcs.h"
#include "unfoldSpectra_TMatrix_funcs.h"
#include "unfoldSpectra_simpbins.h"
#include "unfoldSpectra_anabins.h"



// -----------------------------------------------------------------------------------------------------------------------
const std::string CMSSW_BASE=
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/readForests/outputCondor/";
const std::string SCRATCH_BASE=
  "/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/readForests/";
//const std::string unfoldSpectra_outdir="output/";
const std::string unfoldDataSpectra_outdir="output/unfoldDataSpectra/";
const std::string unfoldMCSpectra_outdir="output/unfoldMCSpectra/";





//Bayes setting that don't change too much
const bool doBayes=true; 
const int kIter = 4; //recommended is 4, default is 4
//,kIterRange=4, kIterDraw = 3, kIterCenter=21;

//SVD setting that don't change too much
const bool doSVD=true; //!(doBayes); 
const int nKregMax  = 9 , kRegRange=(nKregMax-1)/2 ;//max num of diff kregs to do

//other options, useful
const bool useSimplePtBinning=true;//bin by ten everywhere instead of custom binning
const bool doOverUnderflows=false;
const bool clearOverUnderflows=true;
const bool doMCIntegralScaling=true;//unfoldMCSpectra only

//other options, useful...?
const bool zeroBins=false;
const bool normalizedMCMatrix=false;
const bool fillRespHists=false;


//useful strings, numbers
const double integratedLuminosity=27.4*pow(10.,9.);//+/-2.4%
const std::string MCdesc= "QCD PY8 Tune CUETP8M1"; 
const std::string Datadesc1= "pp promptReco, #sqrt{s}=5.02 TeV"; 
const std::string Datadesc2= "L_{int}=27.4 pb^{-1} +/- 2.4%"; 



// color/styles
const Int_t color[13]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,kRed+2,kMagenta+1,kRed+1,kGreen+3};
// fill/fill-border color/style
const Int_t fcolor[5]={kRed-10,kBlue-10,kGreen-8,kOrange-4, kRed-6};
const Int_t lcolor[5]={kRed+1,kBlue-3,kGreen+3,kOrange+3, kRed+2};
// marker color/style
const Int_t fmstyle[6] = {20,21,22,23,29,3};
const Int_t emstyle[6] = {24,25,26,27,30,28};


void drawText(const char *text, float xp, float yp, int size){
  std::cout<<std::endl<<"in drawText"<<std::endl<<std::endl;

  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();

  std::cout<<std::endl<<"drawText done"<<std::endl<<std::endl;
  return;
}
