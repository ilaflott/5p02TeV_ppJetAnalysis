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
  "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/readForests/outputCondor/";
const std::string SCRATCH_BASE=
  "/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/";
//const std::string unfoldSpectra_outdir="output/";
const std::string unfoldDataSpectra_outdir="output/unfoldDataSpectra/";
const std::string unfoldMCSpectra_outdir="output/unfoldMCSpectra/";


//Bayes settings
const int kIter = 4; // recommended is 4, default is 4


//SVD settings 
//    max diff # of kreg to do  /  "width" of kreg from center; i.e. kregs looped over will be kRegCenter +/- kRegRange
const int nKregMax  = 9 , kRegRange=(nKregMax-1)/2 ;      

//other options, useful or interesting to change
//const int kRegDraw  = 4 ; // array entries w/ arguments 0-8. 4 -> middle hist on 3x3 SVDplot
const bool doToyErrs=false;
const bool clearOverUnderflows=false;
const bool doMCIntegralScaling=false;//unfoldMCSpectra only
const bool fillRespHists=false;
const bool doOverUnderflows=false;    // if normalizedMCMatrix=true; this should be false... why?
const bool zeroBins=false;
//const bool normalizedMCMatrix=false;   // do not change for now; 6.24.17

//other options, questionable usefulness, thing before changing...


//useful strings, numbers
//const double integratedLuminosity=27.4*pow(10.,9.);//+/-2.4%
const float integratedLuminosity=27.4*pow(10.,9.);//+/-2.4%
const std::string CMSPRELIM= "CMS PRELIMINARY"; 
const std::string MCdesc= "Py8 Tune CUETP8M1 QCD"; 
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


float computeEffLumi(TFile* finData){
  
  float effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  
  std::cout<<std::endl<<"dataset integrated Luminosity (microbarns) ="<<integratedLuminosity<<std::endl;
  TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_vzCut");
  TH1F *h_NEvents_read    = (TH1F*)finData->Get("NEvents_read");
  
  LumiEff_vz = (float) h_NEvents_vzCut->GetEntries()/h_NEvents_read->GetEntries();
  effIntgrtdLumi_vz=integratedLuminosity*LumiEff_vz;
  
  std::cout<<std::endl<<"lumi efficiency, vz cuts="<<LumiEff_vz<<std::endl;
  std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;
  //LumiEff_vz=0.999;//temp
  return (effIntgrtdLumi_vz);
}



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








void matStylePrint(TH2F* mat, std::string hTitle, TCanvas* canv, std::string outPdfFile, bool useSimpBins){
  
  
  // general for drawRespMatrix ---------------
  //canv->cd();
  
  if(!useSimpBins){
    mat->GetYaxis()->SetMoreLogLabels(true);
    mat->GetYaxis()->SetNoExponent(true);	
    mat->GetXaxis()->SetMoreLogLabels(true);
    mat->GetXaxis()->SetNoExponent(true);      }
  
  mat->GetZaxis()->SetLabelSize(0.025);      
  mat->GetYaxis()->SetLabelSize(0.02);
  mat->GetYaxis()->SetTitleSize(0.023);
  mat->GetYaxis()->SetTitle("gen p_{t}");      
  mat->GetXaxis()->SetLabelSize(0.02);
  mat->GetXaxis()->SetTitleSize(0.025);
  mat->GetXaxis()->SetTitle("reco p_{t}");
  
  // input resp matrix w/ full range ---------------
  
  mat->SetTitle(hTitle.c_str());


  if( hTitle.find("% Errors") != std::string::npos )
    mat->SetAxisRange(0.1,1000.,"Z");
  else if (hTitle.find("Errors") != std::string::npos )
    mat->SetAxisRange(0.000000000000000001,.0001,"Z");
  else if (hTitle.find("Column")  != std::string::npos)
    mat->SetAxisRange(0.000001,1.,"Z");
  else if (hTitle.find("Row")  != std::string::npos)
    mat->SetAxisRange(0.000001,1.,"Z");
  else
    mat->SetAxisRange(0.000000000000000001,.001,"Z");
  
//  if(useSimpBins){
//    if( hTitle.find("%errs") != std::string::npos )
//      mat->SetAxisRange(0.1,1000.,"Z");
//    else if (hTitle.find("errors") != std::string::npos )
//      mat->SetAxisRange(0.000000000000000001,.0001,"Z");
//    else if (hTitle.find("Column")  != std::string::npos)
//      mat->SetAxisRange(0.000001,1.,"Z");
//    else if (hTitle.find("Row")  != std::string::npos)
//      mat->SetAxisRange(0.000001,1.,"Z");
//    else
//      mat->SetAxisRange(0.000000000000000001,.001,"Z");}
//  else{
//    if( hTitle.find("%errs") != std::string::npos )
//      mat->SetAxisRange(0.1,1000.,"Z");
//    else if (hTitle.find("errors") != std::string::npos )
//      mat->SetAxisRange(0.000000000000000001,.0001,"Z");
//    else if (hTitle.find("Column")  != std::string::npos)
//      mat->SetAxisRange(0.000001,1.,"Z");
//    else if (hTitle.find("Row")  != std::string::npos)
//      mat->SetAxisRange(0.000001,1.,"Z");
//    else
//      mat->SetAxisRange(0.000000000000000001,.001,"Z");}
  
  canv->cd();
  
  mat->Draw("COLZ");           
  
  canv->Print(outPdfFile.c_str());
  
  
  return;
}

double* setBinning ( bool useSimpBins, std::string type){
  double* binarray=NULL;
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in setBinning"<<std::endl<<std::endl;
  
  if(useSimpBins){           
    if(type=="reco"){
      binarray   = (double*)simpbins_pt_reco ; }
    else if(type=="gen"){
      binarray   = (double*)simpbins_pt_gen ;}
    else{
      std::cout<<"type not found; bins not set!"<<std::endl;
      std::cout<<"SHUT"<<std::endl;
      std::cout<<"DOWN"<<std::endl;
      std::cout<<"EVERYTHING"<<std::endl;     
      assert(false);    }
  }  
  else{   
    if(type=="reco"){
      binarray   = (double*)anabins_pt_reco ; }
    else if(type=="gen"){
      binarray   = (double*)anabins_pt_gen ;}
    else{
      std::cout<<"type not found; bins not set!"<<std::endl;
      std::cout<<"SHUT"<<std::endl;
      std::cout<<"DOWN"<<std::endl;
      std::cout<<"EVERYTHING"<<std::endl;     
      assert(false);}        
  }
  
  if(funcDebug)std::cout<<std::endl<<"setBinning done, exiting setBinning."<<std::endl<<std::endl;

  return binarray;
  
}




int setNBins ( bool useSimpBins, std::string type){
  int nbins=-1;
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in set/nBins"<<std::endl<<std::endl;
  
  if(useSimpBins){           
    if(type=="reco"){
      nbins   = n_simpbins_pt_reco ; }
    else if(type=="gen"){
      nbins   = n_simpbins_pt_gen ; }
    else{
      std::cout<<"type not found; bins not set!"<<std::endl;
      std::cout<<"SHUT"<<std::endl;
      std::cout<<"DOWN"<<std::endl;
      std::cout<<"EVERYTHING"<<std::endl;     
      assert(false);    }
  }  
  else{   
    if(type=="reco"){
      nbins   = n_anabins_pt_reco ; }
    else if(type=="gen"){
      nbins   = n_anabins_pt_gen ; }
    else{
      std::cout<<"type not found; bins not set!"<<std::endl;
      std::cout<<"SHUT"<<std::endl;
      std::cout<<"DOWN"<<std::endl;
      std::cout<<"EVERYTHING"<<std::endl;     
      assert(false);}        
  }
  
  if(funcDebug)std::cout<<std::endl<<"setBinning done, exiting setBinning."<<std::endl<<std::endl;

  return nbins;
  
}


void setupRatioHist(TH1* h, bool useSimpBins, double* boundaries, int nbins){
  
  h->SetAxisRange(0.2, 1.8, "Y");
  h->GetXaxis()->SetMoreLogLabels(true);
  h->GetXaxis()->SetNoExponent(true);
  //if(!useSimpBins)h->GetXaxis()->SetMoreLogLabels(true);
  //if(!useSimpBins)h->GetXaxis()->SetNoExponent(true);
  h->SetAxisRange( boundaries[0] ,  
		   boundaries[nbins] + 
		   ( boundaries[nbins]-boundaries[nbins-1] ),
		   "X");           
  
  return;
}



void setupSpectraHist(TH1* h, bool useSimpBins, double* boundaries, int nbins){
  
  h->GetXaxis()->SetTitle("jet p_{T} (GeV)");
  //  if(!useSimpBins)h->GetXaxis()->SetMoreLogLabels(true);
  //  if(!useSimpBins)h->GetXaxis()->SetNoExponent(true);
  h->GetXaxis()->SetMoreLogLabels(true);
  h->GetXaxis()->SetNoExponent(true);
  h->SetAxisRange( boundaries[0] ,  
		   boundaries[nbins] + 
		   ( boundaries[nbins]-boundaries[nbins-1] ),
		   "X");           
//h->SetAxisRange( boundaries[0] - 
//		   (boundaries[1] - boundaries[0]), 
//		   boundaries[nbins] + 
//		   (boundaries[nbins]-boundaries[nbins-1]),"X");           
  
  return;
}


