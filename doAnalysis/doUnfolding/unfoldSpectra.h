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
#include <sys/stat.h>
#include <unistd.h>

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

// histos + fits
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
//#include <TH3.h>
//#include <TH1F.h>
//#include <TH2F.h>
//#include <TH1D.h>
//#include <TH2D.h>

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
#include "unfoldSpectra_ptbins.h"
#include "unfoldSpectra_simpbins.h"
#include "unfoldSpectra_anabins.h"
#include "unfoldSpectra_thySpectra.h"


// -----------------------------------------------------------------------------------------------------------------------
const std::string CMSSW_BASE=
  "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/readForests/outputCondor/";
const std::string SCRATCH_BASE=
  "/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/";
const std::string unfoldDataSpectra_outdir="output/unfoldDataSpectra/";
const std::string unfoldMCSpectra_outdir="output/unfoldMCSpectra/";

// -----------------------------------------------------------------------------------------------------------------------
//const int kIter = 4; // default is 4, Bayes
//const int DEF_kIter = 4; // default is 4, Bayes

// max diff # of kreg to do  /  "width" of kreg from center; i.e. kregs looped over will be kRegCenter +/- kRegRange
const int nKiterMax=9, kIterRange=(nKiterMax-1)/2; //Bayes
const int nKregMax  = 9 , kRegRange=(nKregMax-1)/2 ; //SVD
//const int KREGMAXN  = 9 , KREGRAMGE=(KREGMAXN-1)/2 ; //SVD TO DO
const int verbosity=1;



const double NLOMCscaling=1e+00; // i.e. for when i fuck up the scaling/normalization
const double MCscaling=1e+00; // i.e. for when i fuck up the scaling/normalization
const double etaBinWidth=1.;   // e.g. |y| < 2.0
//const double etaBinWidth=4.;   // e.g. |y| < 2.0
//const double etaBinWidth=2.; // e.g. |y| < 1.0
//const double etaBinWidth=1.; // e.g. |y| < 0.5

// for future way of dealing w/ mc response matrix th2
//const float TH2ProjRecPtLoCut=49.;
//const float TH2ProjRecPtHiCut=49.;
//const float TH2ProjGenPtLoCut=937.;
//const float TH2ProjGenPtHiCut=937.;

//const double integratedLuminosity=27.4;//+/-2.4% //int lumi is in pb-1
const double integratedLuminosity=27.4*pow(10.,3.);//+/-2.4% //int lumi is in nb^-1
//const double integratedLuminosity=27.4*pow(10.,6.);//+/-2.4% //int lumi is in microb^-1
//const double integratedLuminosity=27.4*pow(10.,9.);//+/-2.4% //int lumi is in mb^-1

// generally, fill the response hists
// generally, dont use over/under flows (try this again sometime)
// generally, cut off evts from resp matrix by rebinning. But if you do, make sure clear OverUnderflows = true 
// generally, MCIntegral scaling doesn't seem to help
const bool doToyErrs          =false; // error options: kCovToy, kCovariance, kErrors, kNoError
const bool doMCIntegralScaling=false;
const bool fillRespHists      =true;
//const bool useTH2ProjRespHist= fillRespHists&&false;
const bool trimRespMat=false&&fillRespHists;//only reason to do this is to adjust fake/misses correction. only will be done if fillRespHists is true.

// settings that don't really get used
const bool doOverUnderflows   =false;//leave false almost always.
const bool clearOverUnderflows=true; //leave me true almost always. 
const bool zeroBins=false; //leave false almost always.
//const bool normalizedMCMatrix=false;   

// -----------------------------------------------------------------------------------------------------------------------
//useful strings, numbers
const std::string CMSPRELIM= "CMS PRELIMINARY"; 
const std::string MCdesc= "PY8 CUETP8M1, LO QCD MC"; 
const std::string NLOMCdesc= "NLO QCD MC"; 
const std::string Datadesc1= "pp promptReco, #sqrt{s}=5.02 TeV"; 
const std::string Datadesc2= "L_{int}=27.4 pb^{-1} +/- 2.4%"; 

double computeEffLumi(TFile* finData){  
  bool funcDebug=false;
  double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency  

  TH1 *h_NEvents_vzCut   = (TH1*)finData->Get("NEvents_vzCut");
  TH1 *h_NEvents_read    = (TH1*)finData->Get("NEvents_read");
  
  //LumiEff_vz = (double) h_NEvents_vzCut->GetEntries()/h_NEvents_read->GetEntries();//old
  //LumiEff_vz = (double) h_NEvents_vzCut->GetBinContent(2)/h_NEvents_read->GetBinContent(2);//w/ prescale weights  
  LumiEff_vz = (double) h_NEvents_vzCut->GetBinContent(1)/h_NEvents_read->GetBinContent(1); //w/o prescale weights
  if(funcDebug)std::cout<<"eff if using entries = " << (double) h_NEvents_vzCut->GetEntries()/h_NEvents_read->GetEntries() << std::endl;
  if(funcDebug)std::cout<<"eff if using GetBinContent(1) = " << (double) h_NEvents_vzCut->GetBinContent(1)/h_NEvents_read->GetBinContent(1) << std::endl;
  if(funcDebug)std::cout<<"eff if using GetBinContent(2) = " << (double) h_NEvents_vzCut->GetBinContent(2)/h_NEvents_read->GetBinContent(2) << std::endl;
  
  effIntgrtdLumi_vz=integratedLuminosity*LumiEff_vz;
  std::cout<<std::endl<<"dataset integrated Luminosity (microbarns) ="<<integratedLuminosity<<std::endl;
  std::cout<<std::endl<<"lumi efficiency, vz cuts="<<LumiEff_vz<<std::endl;
  std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;
  return (effIntgrtdLumi_vz);
}

void drawText(const char *text, double xp, double yp, int size){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in drawText"<<std::endl<<std::endl;
  
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);//42
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
  
  if(funcDebug)std::cout<<std::endl<<"drawText done"<<std::endl<<std::endl;
  return;
}


void drawTLatex(float x=-1., float y=-1., float y_space=-1.,  TCanvas* canv=NULL, std::vector<std::string> input_strs={0}){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in drawTLatex"<<std::endl;
  if( (x<0.) || (y<0.) || (y_space<0.))return;
  if((!canv))return;
  if(input_strs.size()==0) return;
  
  canv->cd();
  
  for(unsigned int i=0;i<input_strs.size(); i++){
    TLatex *t=new TLatex(x, y, ((std::string)input_strs[i]).c_str());
    t->SetTextFont(63);
    t->SetTextColor(kBlack);
    t->SetTextSize(26);
    //t->SetLineWidth(1);
    t->SetNDC();
    t->Draw();
    y-=y_space;
  }
  
  return;  
}


void matStylePrint(TH2D * mat, std::string hTitle, TCanvas* canv, std::string outPdfFile, bool useSimpBins){
  //bool funcDebug=false;
  //canv->cd();  

  // general for drawRespMatrix ---------------
  if(!useSimpBins){
    mat->GetYaxis()->SetMoreLogLabels(true);
    mat->GetYaxis()->SetNoExponent(true);	
    mat->GetXaxis()->SetMoreLogLabels(true);
    mat->GetXaxis()->SetNoExponent(true);      }
  
  bool isTMatrixD=( (hTitle.find("Pearson")!=std::string::npos )    ||
		    (hTitle.find("Covariance")!=std::string::npos ) ||
		    (hTitle.find("Unfolding")!=std::string::npos ) 		    );
  
  mat->GetZaxis()->SetLabelSize(0.025);      
  
  mat->GetYaxis()->SetLabelSize(0.02);
  mat->GetYaxis()->SetTitleSize(0.023);
  std::string yaxisTitle="GEN Jet p_{T} (GeV)";  
  if(isTMatrixD)yaxisTitle="GEN Jet p_{T} Bin #";  
  mat->GetYaxis()->SetTitle(yaxisTitle.c_str());      
  
  mat->GetXaxis()->SetLabelSize(0.02);
  mat->GetXaxis()->SetTitleSize(0.025);
  std::string xaxisTitle="RECO Jet p_{T} (GeV)";  
  if(isTMatrixD)xaxisTitle="RECO Jet p_{T} Bin #";  
  mat->GetXaxis()->SetTitle(xaxisTitle.c_str());      
  
  // input resp matrix w/ full range ---------------  
  mat->SetTitle(hTitle.c_str());  
  if( hTitle.find("% Errors") != std::string::npos )
    mat->SetAxisRange(1.e-01,1.e+03,"Z");    //mat->SetAxisRange(0.,200.,"Z");
  else if (hTitle.find("Errors") != std::string::npos )
    mat->SetAxisRange(1.e-12,1.e+02,"Z");
  else if (hTitle.find("Column")  != std::string::npos)
    mat->SetAxisRange(1.e-06,1.e+0,"Z");
  else if (hTitle.find("Row")  != std::string::npos)
    mat->SetAxisRange(1.e-06,1.e+0,"Z");
  else if (hTitle.find("Scaled Data Covariance")!= std::string::npos)
    mat->SetAxisRange(1.e-10,1.e-01,"Z");
  else if (hTitle.find("Data Covariance")!= std::string::npos)
    mat->SetAxisRange(1.e+00,1.e+08,"Z");
  else if (hTitle.find("Covariance")  != std::string::npos)//for set data covariance matrix input
    mat->SetAxisRange(1.e-13,1.e-01,"Z");
  else if (hTitle.find("Pearson")  != std::string::npos)
    mat->SetAxisRange(-1.,1.,"Z");
  else if(hTitle.find("Unfolding")  != std::string::npos)
    mat->SetAxisRange(1.e-08,1.e+03,"Z");
  else
    mat->SetAxisRange(1.e-10,1.e+03,"Z");
//  else if (hTitle.find("Covariance")  != std::string::npos)//works if there's no set covariance matrix
//    mat->SetAxisRange(1.e-15,1.e+02,"Z");

  //mat->SetAxisRange(56.,1000.,"Y");
  //mat->SetAxisRange(56.,1000.,"X");
  
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
  //      mat->SetAxisRange(0.000000000000000001,.001,"Z");
  //}
  
  canv->cd();  
  mat->Draw("COLZ");             
  canv->Print(outPdfFile.c_str());
  return;
}

double* setBinning ( bool useSimpBins, std::string type){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in setBinning"<<std::endl<<std::endl;  
  
  if( type!="gen" && type!="reco" ){    
    std::cout<<"type not found; bins not set!"<<std::endl;
    std::cout<<"SHUT"<<std::endl;
    std::cout<<"DOWN"<<std::endl;
    std::cout<<"EVERYTHING"<<std::endl;     
    assert(false);    }
  
  double* binarray=NULL;  
  if(useSimpBins){           
    if(type=="reco")
      binarray   = (double*)simpbins_pt_reco ; 
    else if(type=="gen")
      binarray   = (double*)simpbins_pt_gen  ;  }  
  else{   
    if(type=="reco")
      binarray   = (double*)anabins_pt_reco ; 
    else if(type=="gen")
      binarray   = (double*)anabins_pt_gen ;  }
  
  if(funcDebug)std::cout<<std::endl<<"setBinning done, exiting setBinning."<<std::endl<<std::endl;
  return binarray;  
}


//returns binning array, puts # of bins in nbins
double* setBinning_etabin(int etabinint=0, std::string type="SMP", int* nbins=NULL    ){
  bool funcDebug=true;
  if(funcDebug)std::cout<<std::endl<<"in setBinning_etabin"<<std::endl<<std::endl;  
  
  double* binarray=NULL;   
  if(type=="default_SMP"){//default, unchanged, SMP pt bin edges out to last bin w/ content in it in data
    binarray=(double*)(default_SMP_ptbins[etabinint].data());
    *(nbins)=(int)(default_SMP_ptbins[etabinint].size()-1);
  }
  else if(type=="default2_SMP"){//default, unchanged, SMP pt bin edges out to last bin w/ content in it in data
    binarray=(double*)(default2_SMP_ptbins[etabinint].data());
    *(nbins)=(int)(default2_SMP_ptbins[etabinint].size()-1);
  }
  else if(type=="merged_SMP"){//SMP pt bin edges merged for statistics, for joao's NLO in |y| bins
    binarray=(double*)(merged_SMP_ptbins[etabinint].data());
    *(nbins)=(int)(merged_SMP_ptbins[etabinint].size()-1);
  }
  else if(type=="NLO_SMP"){//SMP pt bin edges for raghav's NLO spectra (EXCEPTION: 686 - 1000 GeV bin for first two |eta| bins, should be 686-1032 GeV)
    binarray=(double*)(NLO_SMP_ptbins[etabinint].data());
    *(nbins)=(int)(NLO_SMP_ptbins[etabinint].size()-1);    
  }
  else if(type=="john"){//for comparison with john OLD
    binarray=(double*)(john_ptbins[etabinint].data());
    *(nbins)=(int)(john_ptbins[etabinint].size()-1);
  }
  else if(type=="john2"){//for comparison with john OLD
    binarray=(double*)(john_ptbins2[etabinint].data());
    *(nbins)=(int)(john_ptbins2[etabinint].size()-1);
  }
  else if(type=="chris"){//for comparison with chris OLD
    binarray=(double*)(chris_ptbins[etabinint].data());
    *(nbins)=(int)(chris_ptbins[etabinint].size()-1);
  }
  else if(type=="chris2"){//for comparison with chris OLD
    if(etabinint!=0)assert(false);//00eta20 only
    binarray=(double*)(chris_ptbins2[etabinint].data());
    *(nbins)=(int)(chris_ptbins2[etabinint].size()-1);
  }
  else
    assert(false);
  
  if(funcDebug)std::cout<<std::endl<<"setBinning_etabin done, exiting."<<std::endl<<std::endl;
  return binarray;  

}


//void setBinning_etabin(
//		       double* binarr1,//gen
//		       double* binarr2,//reco
//		       double* binarr3,//genmat
//		       double* binarr4,//recomat
//		       int arrlength, TH1D* nlohist   ){
//  bool funcDebug=false;
//  if(funcDebug)std::cout<<"setting binning for unfolding single-etabin hist"<<std::endl;
//  int arrindex=0;  
//  int nbins_hist=nlohist->GetNbinsX();
//  if(funcDebug){
//    std::cout<<"arrlength="<<arrlength<<std::endl;
//    std::cout<<"nbins_hist="<<nbins_hist<<std::endl;  }
//  for(int i=1; i<=nbins_hist; i++){
//    if(funcDebug)std::cout<<"__________"<<std::endl;
//    if(funcDebug)std::cout<<"i="<<i<<std::endl;
//    if(funcDebug)std::cout<<"arrindex="<<arrindex<<std::endl;
//    double binedge=-1.;
//    binedge=nlohist->GetBinLowEdge(i);        
//    if(funcDebug)std::cout<<"binedge="<<binedge<<std::endl;
//    if(binedge<56.)continue;
//    binarr1[arrindex]=binedge;    binarr2[arrindex]=binedge;
//    binarr3[arrindex]=binedge;    binarr4[arrindex]=binedge;    
//    if(funcDebug)std::cout<<"binarr1["<<arrindex<<"]="<<binarr1[arrindex]<<std::endl;
//    arrindex++;
//    if(i==nbins_hist){//upper edge of last bin special case or else we run off the end of the array and/or hist. bad!
//      if(funcDebug)std::cout<<"i="<<i<<std::endl;
//      if(funcDebug)std::cout<<"(bool)(i==nbins_hist)="<<(i==nbins_hist)<<std::endl;
//      if(funcDebug)std::cout<<"hitting last bin!"<<std::endl;
//      binedge=nlohist->GetBinLowEdge(i)+nlohist->GetBinWidth(i);       
//      if(funcDebug)std::cout<<"arrindex="<<arrindex<<std::endl;
//      if(funcDebug)std::cout<<"binedge="<<binedge<<std::endl;
//      binarr1[arrindex]=binedge;    binarr2[arrindex]=binedge;
//      binarr3[arrindex]=binedge;    binarr4[arrindex]=binedge;                 
//      if(funcDebug)std::cout<<"binarr1["<<arrindex<<"]="<<binarr1[arrindex]<<std::endl;    }    
//  }
//  return;
//}

void setBinning_etabin_forPY8(
			      double* binarr1,//gen
			      double* binarr2,//reco
			      double* binarr3,//genmat
			      double* binarr4,//recomat
			      int arrlength, TH1D* nlohist   , int etabinint){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"setting binning for unfolding single-etabin hist"<<std::endl;
  int arrindex=0;  
  int nbins_hist=nlohist->GetNbinsX();
  //int nbins_hist=arrlength;//nlohist->GetNbinsX();
  if(funcDebug){
    std::cout<<"arrlength="<<arrlength<<std::endl;
    std::cout<<"nbins_hist="<<nbins_hist<<std::endl;  }
  for(int i=1; i<=nbins_hist; i++){
    
    if(funcDebug)std::cout<<"__________"<<std::endl;
    if(funcDebug)std::cout<<"i="<<i<<std::endl;
    if(funcDebug)std::cout<<"arrindex="<<arrindex<<std::endl;
    double binedge=-1.;
    binedge=nlohist->GetBinLowEdge(i);        
    if(funcDebug)std::cout<<"binedge="<<binedge<<std::endl;
    if(binedge<56.)continue;
    binarr1[arrindex]=binedge;    binarr2[arrindex]=binedge;
    binarr3[arrindex]=binedge;    binarr4[arrindex]=binedge;    
    if(funcDebug)std::cout<<"binarr1["<<arrindex<<"]="<<binarr1[arrindex]<<std::endl;
    arrindex++;
    if(i==nbins_hist){//upper edge of last bin special case or else we run off the end of the array and/or hist. bad!
      if(etabinint==0 || etabinint==1) {
	std::cout<<"not using last bin for etabin==0 or 1"<<std::endl; continue;
      }
      else {
	
	
	if(funcDebug)std::cout<<"i="<<i<<std::endl;
	if(funcDebug)std::cout<<"(bool)(i==nbins_hist)="<<(i==nbins_hist)<<std::endl;
	if(funcDebug)std::cout<<"hitting last bin!"<<std::endl;
	binedge=nlohist->GetBinLowEdge(i)+nlohist->GetBinWidth(i);       
	if(funcDebug)std::cout<<"arrindex="<<arrindex<<std::endl;
	if(funcDebug)std::cout<<"binedge="<<binedge<<std::endl;
	binarr1[arrindex]=binedge;    binarr2[arrindex]=binedge;
	binarr3[arrindex]=binedge;    binarr4[arrindex]=binedge;                 
	if(funcDebug)std::cout<<"binarr1["<<arrindex<<"]="<<binarr1[arrindex]<<std::endl;    }
    }    
  }
  return;
}


int setNBins ( bool useSimpBins, std::string type){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in set/nBins"<<std::endl<<std::endl;
  
  if( type!="gen" && type!="reco" ){    
    std::cout<<"type not found; bins not set!"<<std::endl;
    std::cout<<"SHUT"<<std::endl;
    std::cout<<"DOWN"<<std::endl;
    std::cout<<"EVERYTHING"<<std::endl;     
    assert(false);    }

  int nbins=-1;  
  if(useSimpBins){           
    if(type=="reco")
      nbins   = n_simpbins_pt_reco ; 
    else if(type=="gen")
      nbins   = n_simpbins_pt_gen ;   }  
  else{   
    if(type=="reco")
      nbins   = n_anabins_pt_reco ; 
    else if(type=="gen")
      nbins   = n_anabins_pt_gen ;   }
  
  if(funcDebug)std::cout<<std::endl<<"setBinning done, exiting setBinning."<<std::endl<<std::endl;
  return nbins;  
}

void setupRatioHist(TH1* h, bool useSimpBins, double* boundaries=NULL, int nbins=1){  
  //h->SetAxisRange(0.2, 1.8, "Y");
  //h->SetAxisRange(0.3, 1.7, "Y");
  //h->SetAxisRange(0.4, 1.6, "Y");
  h->SetAxisRange(0.6, 1.4, "Y");
  h->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  h->GetXaxis()->SetMoreLogLabels(true);
  h->GetXaxis()->SetNoExponent(true);
  h->GetYaxis()->CenterTitle(true);
  //if(!useSimpBins)h->GetXaxis()->SetMoreLogLabels(true);
  //if(!useSimpBins)h->GetXaxis()->SetNoExponent(true);
  //h->SetAxisRange( boundaries[0] ,  
  //		   boundaries[nbins] + 
  //		   ( boundaries[nbins]-boundaries[nbins-1] ),
  //		   "X");             
  return;
}
void setupSystPercErrsHist(TH1* h, bool useSimpBins, double* boundaries=NULL, int nbins=1){  
  //h->SetAxisRange(0.2, 1.8, "Y");
  //h->SetAxisRange(0.3, 1.7, "Y");
  //h->SetAxisRange(0.4, 1.6, "Y");
  h->SetAxisRange(-50., 50., "Y");
  h->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  h->GetXaxis()->SetMoreLogLabels(true);
  h->GetXaxis()->SetNoExponent(true);
  //if(!useSimpBins)h->GetXaxis()->SetMoreLogLabels(true);
  //if(!useSimpBins)h->GetXaxis()->SetNoExponent(true);
  //h->SetAxisRange( boundaries[0] ,  
  //		   boundaries[nbins] + 
  //		   ( boundaries[nbins]-boundaries[nbins-1] ),
  //		   "X");             
  return;
}

void setupSpectraHist(TH1* h, bool useSimpBins=false, double* boundaries=NULL, int nbins=1){  
  h->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#eta dp_{T}} [nb/GeV]");
  //h->GetYaxis()->SetTitle("N_{Jets}/L_{int}");
  //h->GetYaxis()->SetTitle("A.U.");
  h->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  h->GetXaxis()->SetMoreLogLabels(true);
  h->GetXaxis()->SetNoExponent(true);  
  //  if(!useSimpBins)h->GetXaxis()->SetMoreLogLabels(true);
  //  if(!useSimpBins)h->GetXaxis()->SetNoExponent(true);
  //  h->SetAxisRange( boundaries[0] ,  
  //		   boundaries[nbins] + 
  //		   ( boundaries[nbins]-boundaries[nbins-1] ),
  //		   "X");           
  //h->SetAxisRange( boundaries[0] - 
  //		   (boundaries[1] - boundaries[0]), 
  //		   boundaries[nbins] + 
  //		   (boundaries[nbins]-boundaries[nbins-1]),"X");             
  return;
}

inline bool fileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

//std::string getUnfDetails(TH1D* hunf){
void getUnfDetails(TH1D* hunf, std::string* unfdetails_str){
  (*unfdetails_str)="";
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in getUnfDetails"<<std::endl;
  std::string hunf_title=(std::string)hunf->GetTitle();
  if(funcDebug)std::cout<<"hunf_title="<<hunf_title<<std::endl;
  bool SVDUnf=false, BayesUnf=false;
  size_t details_pos;
  if(hunf_title.find("kReg")!=std::string::npos){
    SVDUnf=true; BayesUnf=false;
    details_pos=hunf_title.find("kReg");
  }
  else if (hunf_title.find("kIter")!=std::string::npos){
    SVDUnf=false; BayesUnf=true;
    details_pos=hunf_title.find("kIter");
  }
  else {
    if(funcDebug)std::cout<<"not bayesian nor SVD unfolding! exit!"<<std::endl;
    return ;//unfdetails_str;    
  }
  if(!SVDUnf && !BayesUnf){
    if(funcDebug)std::cout<<"not bayesian nor SVD unfolding! exit!"<<std::endl;
    return ;//unfdetails_str;    
  }
  
  //unsigned int title_length=hunf_title.length();
  //(*unfdetails_str)=hunf_title;
  std::string dummy_title=hunf_title;
  if(funcDebug)std::cout<<"before replacement, dummy_title="<<dummy_title<<std::endl;
  dummy_title.replace(0, details_pos, "");
  if(funcDebug)std::cout<<"after replacement, dummy_title="<<dummy_title<<std::endl;
  (*unfdetails_str)=dummy_title;
  return;// unfdetails_str;
}

void makeNPCorrName(std::string* NPs, std::string* NPs_name){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in makeNPCorrName"<<std::endl;
  std::string dummy=*NPs;
  if(funcDebug)std::cout<<"NPs_name="<<dummy<<std::endl;
  dummy=dummy.replace(0, (std::string(" #otimes ")).length(),"");
  if(funcDebug)std::cout<<"NPs_name="<<dummy<<std::endl;
  if(dummy.find("+")!=std::string::npos)    dummy=dummy.replace(dummy.find("+"), 1, "_");  
  else                                      dummy=dummy.replace(dummy.find(" "), 1, "_");
  if(funcDebug)std::cout<<"NPs_name="<<dummy<<std::endl;
  *NPs_name=dummy;
  return;
}

//void makeNPs(std::string* NPTH1_name, std::string* NPs){
//  
//}
void checkNRenameFiles (const std::string outFileName, std::string *outRespMatPdfFile, std::string *outPdfFile, std::string *outRootFile, std::string *out3x3PdfFile=NULL){
  bool funcDebug=false;

  int outputInd=1;  
  while( (bool)fileExists(*(outRootFile)) ) {
    if(funcDebug)
      std::cout<<"fileExists! adding ind="<<outputInd<<std::endl;
    (*outRespMatPdfFile)=outFileName+"_"+std::to_string(outputInd)+"_respMat.pdf";
    (*outPdfFile       )=outFileName+"_"+std::to_string(outputInd)+".pdf";
    (*outRootFile      )=outFileName+"_"+std::to_string(outputInd)+".root";      
    if((bool)out3x3PdfFile)
      (*out3x3PdfFile       )=outFileName+"_"+std::to_string(outputInd)+"_3x3.pdf";    
    outputInd++;  
    if(funcDebug)      std::cout<<"outPdffile="<<(*outPdfFile)<<std::endl;    
  }   
  return;
}

void drawRespMatrixFile(TH2D* hmat, TH2D* hmat_rebin, TH2D* hmat_errors,
			//TH2D* covmat_TH2,	//TH2D* Pearson,	//TH2D* unfmat_TH2,
			double* boundaries_pt_reco_mat , int nbins_pt_reco_mat,
			double* boundaries_pt_gen_mat , int nbins_pt_gen_mat,
			std::string outRespMatPdfFile   , bool  useSimpBins, TFile* fpp_MC=NULL){
  bool funcDebug=true;
  std::cout<<std::endl<<"drawing input response matrices..."<<std::endl;    
  if(!hmat)assert(false);
  std::string open_outRespMatPdfFile=outRespMatPdfFile+"[";      
  std::string close_outRespMatPdfFile=outRespMatPdfFile+"]";    
  
  //for TH2s, always make the canvas wider than it's height. 
  //The reason is the z-axis label in COLZ mode adds extra width, so for a symmetric canvas, the #'s on the axis are cut off and unreadable
  TCanvas* tempCanvForPdfPrint=new TCanvas("tempCanv_respMat","",1400,1000);    
  tempCanvForPdfPrint->cd();    
  if(useSimpBins){	tempCanvForPdfPrint->SetLogx(0);
    tempCanvForPdfPrint->SetLogy(0);       
    tempCanvForPdfPrint->SetLogz(1);          }
  else {	tempCanvForPdfPrint->SetLogx(1);
    tempCanvForPdfPrint->SetLogy(1);       
    tempCanvForPdfPrint->SetLogz(1);         }       
  
  // open file    
  if(funcDebug)std::cout<<"opening response matrix pdf file"<<std::endl;
  tempCanvForPdfPrint->Print(open_outRespMatPdfFile.c_str()); 
  

  // orig matrix ---------------          
  //bool isForNLOMC=(((std::string)hmat->GetName()).find("ynew")!=std::string::npos);
  bool isForNLOMC=(((std::string)hmat->GetName()).find("_th2")!=std::string::npos);
  double xLow,xHi;//=boundaries_pt_reco_mat[0];
  double yLow, yHi;//=boundaries_pt_gen_mat[0];
  if( (bool)fpp_MC && !isForNLOMC){
    xLow=(double)std::atof( ( (TH1*) (fpp_MC->Get("hJetPtCut") )
			      )->GetTitle() );
    xHi=(double)std::atof( ( (TH1*) (fpp_MC->Get("hJetPtCut_Hi") )
			     )->GetTitle() );
    yLow=(double)std::atof( ( (TH1*) (fpp_MC->Get("hGenJetPtCut") )
			      )->GetTitle() );
    yHi=(double)std::atof( ( (TH1*) (fpp_MC->Get("hGenJetPtCut_Hi") )
			     )->GetTitle() );	
    hmat->SetAxisRange(xLow,xHi,"X");
    hmat->SetAxisRange(yLow,yHi,"Y");        
  }
  else{
    xLow=boundaries_pt_reco_mat[0];
    xHi=boundaries_pt_reco_mat[nbins_pt_reco_mat];
    yLow=boundaries_pt_gen_mat[0];
    yHi=boundaries_pt_gen_mat[nbins_pt_gen_mat]; }
  
  
  
  //assert(false);
  if(funcDebug)std::cout<<"matStylePrint(hmat)"<<std::endl;
  matStylePrint(hmat, "MC Response Matrix, Original Histogram", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);      
  
  // orig matrix w/ used pt range ---------------    
  tempCanvForPdfPrint->cd();    
  hmat->SetTitle("MC Response Matrix, p_{T} Range Used for Unfolding");
  if(!isForNLOMC)hmat->SetAxisRange(boundaries_pt_reco_mat[0],boundaries_pt_reco_mat[nbins_pt_reco_mat],"X");
  if(!isForNLOMC)hmat->SetAxisRange(boundaries_pt_gen_mat[0],boundaries_pt_gen_mat[nbins_pt_gen_mat],"Y");      
  hmat->Draw("COLZ");               
  tempCanvForPdfPrint->Print(outRespMatPdfFile.c_str());      
  
  
  // matrix rebinned ---------------    
  if(hmat_rebin){
    if(funcDebug)std::cout<<"matStylePrint(hmat_rebin)"<<std::endl;
    matStylePrint(hmat_rebin, "MC Response Matrix, Rebinned", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  }
  
  if(hmat_errors){
    // error matrix in binning of interest ---------------    
    if(funcDebug)std::cout<<"matStylePrint(hmat_errors)"<<std::endl;
    matStylePrint(hmat_errors, "MC Response Matrix, Bin Errors", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  }
  
  TH2D* hmat_percenterrs=NULL;
  if( (hmat_rebin) && (hmat_errors)){
    // percent error matrix in binning of interest ---------------    
    if(funcDebug)std::cout<<"makeRespMatrixPercentErrs(hmat_errors, hmat_rebin)"<<std::endl;
    hmat_percenterrs= makeRespMatrixPercentErrs( (TH2D*) hmat_errors, (TH2D*) hmat_rebin,
						 (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						 (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
    if(funcDebug) hmat_percenterrs->Print("base");
    if(funcDebug)std::cout<<"matStylePrint(hmat_percenterrs)"<<std::endl;
    matStylePrint(hmat_percenterrs, "MC Response Matrix, Bin % Errors", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  }
  
  TH2D* hmat_rebin_colnormd=NULL, *hmat_rebin_rownormd=NULL;
  if(hmat_rebin){
    // col normd matrix in binning of interest  ---------------    
    // POTENTIAL ISSUE: using the resp matrix post rebinning/divide bin width/clearing overflows... should i be using the original matrix?
    if(funcDebug)std::cout<<"normalizeCol_RespMatrix(hmat_rebin)"<<std::endl;
    hmat_rebin_colnormd = normalizeCol_RespMatrix( (TH2D*)  hmat_rebin,
						   (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						   (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
    if(funcDebug)  hmat_rebin_colnormd->Print("base");  
    
    
    if(funcDebug)std::cout<<"matStylePrint(hmat_rebin_colnormd)"<<std::endl;
    matStylePrint(hmat_rebin_colnormd, "MC Response Matrix, Columns Sum to 1", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
    
    // row normd matrix in binning of interest  ---------------    
    if(funcDebug)std::cout<<"normalizeRow_RespMatrix(hmat_rebin)"<<std::endl;
    hmat_rebin_rownormd = normalizeRow_RespMatrix( (TH2D*)  hmat_rebin,
						   (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						   (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
    if(funcDebug)  hmat_rebin_rownormd->Print("base");
    if(funcDebug)std::cout<<"matStylePrint(hmat_rebin_rownormd)"<<std::endl;
    matStylePrint(hmat_rebin_rownormd, "MC Response Matrix, Rows Sum to 1", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  }
  
  // close file     
  if(funcDebug)std::cout<<"closing outRespMatPdfFile"<<std::endl;
  tempCanvForPdfPrint->Print(close_outRespMatPdfFile.c_str());   

  if(funcDebug)std::cout<<"done. exiting drawRespMatrixFile"<<std::endl;
  return;
}


TH1* makeThyHist_00eta20(std::string filename, bool applyNPCorrFactor=true){
  bool funcDebug=false;
  
  TFile* thyFile=TFile::Open(filename.c_str());
  if(!thyFile){
    std::cout<<"err no file"<<std::endl;
    return NULL;  }
  
  TH1* h_00eta05=(TH1*)thyFile->Get("h0100100");   //  TH1* h_00eta05=(TH1*)thyFile->Get("h1100100");
  TH1* h_05eta10=(TH1*)thyFile->Get("h0100200");   //  TH1* h_05eta10=(TH1*)thyFile->Get("h1100200");
  TH1* h_10eta15=(TH1*)thyFile->Get("h0100300");   //  TH1* h_10eta15=(TH1*)thyFile->Get("h1100300");
  TH1* h_15eta20=(TH1*)thyFile->Get("h0100400");   //  TH1* h_15eta20=(TH1*)thyFile->Get("h1100400");
  
  //highest eta bin has smallest pt range (this is how raghav rigged the files, no choice)
  //all eta bin hists start at the same pt; 40-something. they only differ in where they end.
  int nbinsx_forSum = h_15eta20->GetNbinsX();
  double finalBinLowEdge_forSum=h_15eta20->GetBinLowEdge(nbinsx_forSum);
  if(funcDebug)std::cout<< "nbinsx_forSum="         << nbinsx_forSum          << std::endl;
  if(funcDebug)std::cout<< "finalBinLowEdge_forSum="<< finalBinLowEdge_forSum << std::endl;
  
  int nbinsx_max =h_00eta05->GetNbinsX();
  double finalBinLowEdge_max=h_00eta05->GetBinLowEdge(nbinsx_max);
  if(funcDebug)std::cout<< "nbinsx_max="         << nbinsx_max          << std::endl;
  if(funcDebug)std::cout<< "finalBinLowEdge_max="<< finalBinLowEdge_max << std::endl;
  
  //this is the hist that will be returned, clone for binning
  TH1* thyHist=(TH1*)h_00eta05->Clone("sumThyHistClone_forBins");
  thyHist->Reset("ICES");//resets integ, contents, errs, stats
  thyHist->Reset("M");//resets max/min
  //thyHist->SumW2();
  if(funcDebug)thyHist->Print("base");
  
  TFile* NPCorrFile=NULL;
  TF1* fNPCorr_00eta05=NULL;
  TF1* fNPCorr_05eta10=NULL;
  TF1* fNPCorr_10eta15=NULL;
  TF1* fNPCorr_15eta20=NULL;

  if(applyNPCorrFactor){
    NPCorrFile=(TFile*)TFile::Open(NPCorr_filename.c_str());
    if(NPCorr_filename.find("nnlo")!=std::string::npos){
      fNPCorr_00eta05=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin0");
      fNPCorr_05eta10=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin1");
      fNPCorr_10eta15=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin2");
      fNPCorr_15eta20=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin3");    }
    else{
      fNPCorr_00eta05=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin0");
      fNPCorr_05eta10=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin1");
      fNPCorr_10eta15=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin2");
      fNPCorr_15eta20=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin3");    }
  }
  
  
  for(int i=1;i<=nbinsx_max;i++){
    if(funcDebug)std::cout<<"i="<<i<<std::endl;
    
    if(i<=nbinsx_forSum){
      
      double NPCF_00eta05=1.;
      double NPCF_05eta10=1.;
      double NPCF_10eta15=1.;
      double NPCF_15eta20=1.;
      if(applyNPCorrFactor&&((bool)NPCorrFile)){
	double bincenter=(double)thyHist->GetXaxis()->GetBinCenter(i);
	NPCF_00eta05=fNPCorr_00eta05->Eval(bincenter) ;
	NPCF_05eta10=fNPCorr_05eta10->Eval(bincenter) ;
	NPCF_10eta15=fNPCorr_10eta15->Eval(bincenter) ;
	NPCF_15eta20=fNPCorr_15eta20->Eval(bincenter) ;  
	if(funcDebug){
	  std::cout << "NPCF_00eta05 = " << NPCF_00eta05 << std::endl;
	  std::cout << "NPCF_05eta10 = " << NPCF_05eta10 << std::endl;
	  std::cout << "NPCF_10eta15 = " << NPCF_10eta15 << std::endl;
	  std::cout << "NPCF_15eta20 = " << NPCF_15eta20 << std::endl;	}
	//if(NPCF_00eta05>1.)NPCF_00eta05=1.;
	//if(NPCF_05eta10>1.)NPCF_05eta10=1.;
	//if(NPCF_10eta15>1.)NPCF_10eta15=1.;
	//if(NPCF_15eta20>1.)NPCF_15eta20=1.;	
      }
      
      double sumcontent=0.;     
      sumcontent+=(h_00eta05->GetBinContent(i)*NPCF_00eta05 );
      sumcontent+=(h_05eta10->GetBinContent(i)*NPCF_05eta10 );
      sumcontent+=(h_10eta15->GetBinContent(i)*NPCF_10eta15 );
      sumcontent+=(h_15eta20->GetBinContent(i)*NPCF_15eta20 );
      sumcontent/=4.;//etabin width
      sumcontent/=1000.;//picobarns to nanobarns
      if(funcDebug)std::cout<<"sumcontent="<<sumcontent<<std::endl;
      thyHist->SetBinContent(i,sumcontent);
      
      double sumerr=0.;
      sumerr+=h_00eta05->GetBinError(i)*h_00eta05->GetBinError(i)*NPCF_00eta05*NPCF_00eta05;      
      sumerr+=h_05eta10->GetBinError(i)*h_05eta10->GetBinError(i)*NPCF_05eta10*NPCF_05eta10;
      sumerr+=h_10eta15->GetBinError(i)*h_10eta15->GetBinError(i)*NPCF_10eta15*NPCF_10eta15;
      sumerr+=h_15eta20->GetBinError(i)*h_15eta20->GetBinError(i)*NPCF_15eta20*NPCF_15eta20;
      sumerr=sqrt(sumerr);//errs
      sumerr/=4.;//etabin width
      sumerr/=1000.;//picobarns to nanobarns      
      if(funcDebug)std::cout<<"sumerr="<<sumerr<<std::endl;
      thyHist->SetBinError(i,sumerr);
    }
    if(funcDebug)std::cout<<"thyHist BinConent="<<thyHist->GetBinContent(i)<<std::endl;
    if(funcDebug)std::cout<<"thyHist BinError=" <<thyHist->GetBinError(i)<<std::endl;          
  }
  
  if(funcDebug)std::cout<<"end of makeThyHist_00eta20"<<std::endl;
  if(funcDebug)thyHist->Print("base");
  return thyHist;
}



TH1* makeThyHist_00eta20_v2(std::string filename, bool applyNPCorrFactor=true, std::string* NPCorrString=NULL){
  bool funcDebug=false;
  
  TFile* thyFile=TFile::Open(filename.c_str());
  if(!thyFile){
    std::cout<<"err no file"<<std::endl;
    return NULL;  }
  
  TH1* h_00eta05=(TH1*)thyFile->Get("h0100100");   //  TH1* h_00eta05=(TH1*)thyFile->Get("h1100100");
  TH1* h_05eta10=(TH1*)thyFile->Get("h0100200");   //  TH1* h_05eta10=(TH1*)thyFile->Get("h1100200");
  TH1* h_10eta15=(TH1*)thyFile->Get("h0100300");   //  TH1* h_10eta15=(TH1*)thyFile->Get("h1100300");
  TH1* h_15eta20=(TH1*)thyFile->Get("h0100400");   //  TH1* h_15eta20=(TH1*)thyFile->Get("h1100400");
  int nbinsx_00eta05=h_00eta05->GetNbinsX();
  int nbinsx_05eta10=h_05eta10->GetNbinsX();
  int nbinsx_10eta15=h_10eta15->GetNbinsX();
  int nbinsx_15eta20=h_15eta20->GetNbinsX();
  
  int nbinsx =h_00eta05->GetNbinsX();
  double finalBinHighEdge=h_00eta05->GetBinLowEdge(nbinsx) + h_00eta05->GetBinWidth(nbinsx);
  double firstBinLowEdge=h_00eta05->GetBinLowEdge(1);
  if(funcDebug)std::cout<< "nbinsx="         << nbinsx         << std::endl;
  if(funcDebug)std::cout<< "finalBinHighEdge="         << finalBinHighEdge         << std::endl;
  if(funcDebug)std::cout<< "firstBinLowEdge="          << firstBinLowEdge         << std::endl;
  
  //this is the hist that will be returned, clone the most central |y| hist for binning (less-central |y| hists' pt binning are subsets of this one) 
  TH1* thyHist=(TH1*)h_00eta05->Clone("sumThyHistClone_forBins");
  thyHist->Reset("ICES");//resets integ, contents, errs, stats
  thyHist->Reset("M");//resets max/min
  thyHist->Sumw2(true);
  if(funcDebug)thyHist->Print("base");
  
  TFile* NPCorrFile=NULL;
  TF1* fNPCorr_00eta05=NULL;
  TF1* fNPCorr_05eta10=NULL;
  TF1* fNPCorr_10eta15=NULL;
  TF1* fNPCorr_15eta20=NULL;
  
  if(applyNPCorrFactor){
    NPCorrFile=(TFile*)TFile::Open(NPCorr_filename.c_str());
    //    if(NPCorr_filename.find("nnlo")!=std::string::npos){
    if(filename.find("nnlo")!=std::string::npos){
      fNPCorr_00eta05=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin0");
      fNPCorr_05eta10=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin1");
      fNPCorr_10eta15=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin2");
      fNPCorr_15eta20=(TF1*)NPCorrFile->Get("fNPC_POWPY8_R4_etabin3");    
      //if(NPCorrString)(*NPCorrString)=" #otimes POW+PY8 NPs";
      if(NPCorrString)(*NPCorrString)=" #otimes POW+PY8";
    }
    else{
      fNPCorr_00eta05=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin0");
      fNPCorr_05eta10=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin1");
      fNPCorr_10eta15=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin2");
      fNPCorr_15eta20=(TF1*)NPCorrFile->Get("fNPC_HerwigEE4C_R4_etabin3");    
      //      if(NPCorrString)(*NPCorrString)=" #otimes HERWIG EE4C NPs";
      if(NPCorrString)(*NPCorrString)=" #otimes HERWIG EE4C";
    }
    if(!((bool)NPCorrFile)){
      std::cout<<"warning; want to open NP correction file + apply corrections but NP correction file not open nor present. Exit."<<std::endl;
      assert(false);
    }
  }
  
  
  for(int i=1;i<=nbinsx;i++){
    
    
    double NPCF_00eta05=1.;
    double NPCF_05eta10=1.;
    double NPCF_10eta15=1.;
    double NPCF_15eta20=1.;    
    
    if(applyNPCorrFactor){
      double bincenter=(double)thyHist->GetXaxis()->GetBinCenter(i);
      if(i<=nbinsx_00eta05)NPCF_00eta05=fNPCorr_00eta05->Eval(bincenter) ;
      if(i<=nbinsx_05eta10)NPCF_05eta10=fNPCorr_05eta10->Eval(bincenter) ;
      if(i<=nbinsx_10eta15)NPCF_10eta15=fNPCorr_10eta15->Eval(bincenter) ;
      if(i<=nbinsx_15eta20)NPCF_15eta20=fNPCorr_15eta20->Eval(bincenter) ;  
      
      if(funcDebug){
	std::cout << "NPCF_00eta05 = " << NPCF_00eta05 << std::endl;
	std::cout << "NPCF_05eta10 = " << NPCF_05eta10 << std::endl;
	std::cout << "NPCF_10eta15 = " << NPCF_10eta15 << std::endl;
	std::cout << "NPCF_15eta20 = " << NPCF_15eta20 << std::endl;	}      
      
      //if(NPCF_00eta05>1.)NPCF_00eta05=1.;
      //if(NPCF_05eta10>1.)NPCF_05eta10=1.;
      //if(NPCF_10eta15>1.)NPCF_10eta15=1.;
      //if(NPCF_15eta20>1.)NPCF_15eta20=1.;
      
      if(funcDebug){
	std::cout << "NPCF_00eta05 = " << NPCF_00eta05 << std::endl;
	std::cout << "NPCF_05eta10 = " << NPCF_05eta10 << std::endl;
	std::cout << "NPCF_10eta15 = " << NPCF_10eta15 << std::endl;
	std::cout << "NPCF_15eta20 = " << NPCF_15eta20 << std::endl;	}      
    }
    
    
    if(funcDebug)std::cout<<"i="<<i<<std::endl;
    double sumcontent=0.;     
    
    if(i<=nbinsx_00eta05){
      sumcontent+=(h_00eta05->GetBinContent(i)*NPCF_00eta05 );}
    if(i<=nbinsx_05eta10){
      sumcontent+=(h_05eta10->GetBinContent(i)*NPCF_05eta10 );}
    if(i<=nbinsx_10eta15){
      sumcontent+=(h_10eta15->GetBinContent(i)*NPCF_10eta15 );}
    if(i<=nbinsx_15eta20){
      sumcontent+=(h_15eta20->GetBinContent(i)*NPCF_15eta20 );}
    
    sumcontent/=4.;//etabin width
    sumcontent/=1000.;//picobarns to nanobarns
    if(funcDebug)std::cout<<"sumcontent="<<sumcontent<<std::endl;
    thyHist->SetBinContent(i,sumcontent);
    
    double sumerr=0.;
    sumerr+=h_00eta05->GetBinError(i)*h_00eta05->GetBinError(i)*NPCF_00eta05*NPCF_00eta05;      
    sumerr+=h_05eta10->GetBinError(i)*h_05eta10->GetBinError(i)*NPCF_05eta10*NPCF_05eta10;
    sumerr+=h_10eta15->GetBinError(i)*h_10eta15->GetBinError(i)*NPCF_10eta15*NPCF_10eta15;
    sumerr+=h_15eta20->GetBinError(i)*h_15eta20->GetBinError(i)*NPCF_15eta20*NPCF_15eta20;
    sumerr=sqrt(sumerr);//errs
    sumerr/=4.;//etabin width
    sumerr/=1000.;//picobarns to nanobarns      
    if(funcDebug)std::cout<<"sumerr="<<sumerr<<std::endl;
    thyHist->SetBinError(i,sumerr);
    
    
    if(funcDebug)std::cout<<"thyHist BinConent="<<thyHist->GetBinContent(i)<<std::endl;
    if(funcDebug)std::cout<<"thyHist BinError=" <<thyHist->GetBinError(i)<<std::endl;          
  }

  if(funcDebug)std::cout<<"end of makeThyHist_00eta20"<<std::endl;
  if(funcDebug)thyHist->Print("base");
  return thyHist;
}

void applyNPCorr_etabin(std::string nlofile_str="", TH1D* hnlo=NULL, std::string* NPsstr=NULL, int etabinint=0, std::string sysdirec=""){
  bool funcDebug=true;
  
  
  //open NP file
  TFile* NPCorrFile=NULL;
  NPCorrFile=(TFile*)TFile::Open(NPCorr_filename.c_str());
  if(NPCorrFile->IsZombie())assert(false);
  
  //set NP TF1 name + NP desc str if not null
  std::string NPCorrName="";  
  if(NPCorr_filename.find("fromJohn")!=std::string::npos){
    //NPCorrName="f"+std::to_string(etabinint)+sysdirec;
    NPCorrName="fNPC_AVG_R4_etabin"+std::to_string(etabinint)+sysdirec;
    //if(NPsstr)(*NPsstr)=" #otimes PY8";
    if(NPsstr)(*NPsstr)=" #otimes AVG(PY8,HWG)";
  }
  else if(nlofile_str.find("nnlo")!=std::string::npos){    //if nnpdf nnlo, then use these nps
    //    NPCorrName="fNPC_POWPY8_R4_etabin";
    //if(NPsstr)(*NPsstr)=" #otimes POW+PY8";  
    NPCorrName="fNPC_HerwigEE5C_R4_etabin"+std::to_string(etabinint);
    if(NPsstr)(*NPsstr)=" #otimes HERWIG EE5C";  
  }
  else{    
    NPCorrName="fNPC_HerwigEE4C_R4_etabin"+std::to_string(etabinint);
    if(NPsstr)(*NPsstr)=" #otimes HERWIG EE4C";  }  
  //NPCorrName+=

  //apply the NP to the vals+errs in hnlo
  
  //grab the TF1
  TF1* fNPCorr=(TF1*)NPCorrFile->Get(NPCorrName.c_str());
  if(funcDebug)std::cout<<"NPCorr_filename is: " << NPCorr_filename<<std::endl;
  if(funcDebug)std::cout<<"contents of the file are..."<<std::endl;
  if(funcDebug)NPCorrFile->ls();
  if(funcDebug)std::cout<<"NP Fit Obj Name is: " << NPCorrName<<std::endl;
  if(funcDebug)std::cout<<"hnlo Name is: " << hnlo->GetName()<<std::endl;
  if(funcDebug)std::cout<<"fNPCorr->GetName()="<< fNPCorr->GetName() <<std::endl;
  if(!fNPCorr)assert(false);
  
  int nbins=hnlo->GetNbinsX();
  for(int j=1;j<=nbins;j++){
    double hnlo_NPcorrval=(hnlo->GetBinContent(j))*(fNPCorr->Eval(hnlo->GetBinCenter(j)));
    double hnlo_NPcorrval_err=(hnlo->GetBinError(j))*(fNPCorr->Eval(hnlo->GetBinCenter(j)));
    hnlo->SetBinContent(j, hnlo_NPcorrval);
    hnlo->SetBinError(j, hnlo_NPcorrval_err);  }
  
  if(funcDebug)std::cout<<"done with applyNPCorr_etabin"<<std::endl;
  return;
}

TH2D* absVal_TH2Content(TH2D* hist){
  //bool funcDebug=false;  
  int  nbinsx=(int)hist->GetNbinsX();
  int  nbinsy=(int)hist->GetNbinsY();  
  std::string newname=(std::string)hist->GetName();
  newname+="_absvalclone";  
  TH2D* clonehist=(TH2D*)hist->Clone(newname.c_str());
  for(int i = 1; i<=nbinsx; i++){
    for(int j = 1;j<=nbinsy; j++){
      double val=hist->GetBinContent(i,j);
      if(val<0.)val*=-1;
      clonehist->SetBinContent(i,j,val);    
    }  }
  return (TH2D*)clonehist;
}

void divBinWidth_DiAndSVals(double * binning, int binninglength, TH1D* hdisval){
  //bool funcDebug=false;
  int numbins=hdisval->GetNbinsX();
  for(int i=1; i<=numbins; i++){
    double val= hdisval->GetBinContent(i);
    double width=binning[i]-binning[i-1];
    val/=width;
    hdisval->SetBinContent(i, val);
  }
  return;
}


void init_kRegParamArray(int* kReg, int nbins, int kRegInput){
  bool funcDebug=true;
  
  int kReg_lo=kRegInput-kRegRange+1,kReg_hi=kRegInput+kRegRange;
  if(funcDebug) std::cout<<std::endl<<"length of kReg parameter array = "<<nKregMax<<" elements"<<std::endl;
  if(funcDebug) std::cout<<"kRegInput = "<< kRegInput<<std::endl;  
  if(funcDebug) std::cout<<(kReg_lo)<<" <= kReg[i] <= "<<(kReg_hi)<<std::endl;
  if(funcDebug) std::cout<<" nbins/2 = "<<nbins<<"/2 = "<< nbins/2<<std::endl;
  if(kReg_hi>nbins){
    if(funcDebug) std::cout<<"kRegInput was = "<< kRegInput<<std::endl;  
    kRegInput=nbins-1-kRegRange;
    kReg_lo=kRegInput-kRegRange+1;
    kReg_hi=kRegInput+kRegRange;
    std::cout<<"kRegInput too big. Setting to largest possible value."<<std::endl;
    if(funcDebug) std::cout<<"kRegInput is now = "<< kRegInput<<std::endl;      
    if(funcDebug) std::cout<<(kReg_lo)<<" <= kReg[i] <= "<<(kReg_hi)<<std::endl;
    if(funcDebug) std::cout<<" nbins/2 = "<<nbins<<"/2 = "<< nbins/2<<std::endl;      }
  bool addOne_avdOne=false; //in case the kRegParam=1 (which is a useless case to check, the normalization is consistently wrong in this casw
  bool addOne_avdkr0=false;
  for(int i=0; i<nKregMax; ++i)   {
    int kRegParam=kReg_lo+i-1;
    if( i==0 )  
      kReg[0]=(nbins/2);
    
    else if( i>=1){
      if(i==1 && kRegParam==1)
	addOne_avdOne=true;
      if(kRegParam==kReg[0])
	addOne_avdkr0=true;
      
      
      kReg[i]=kRegParam;

      if(addOne_avdOne)
	kReg[i]+=1;
      
      if(addOne_avdkr0)//avoid nbins/2
	kReg[i]+=1;
      
      
    }
    if(funcDebug) std::cout<<"kReg["<<i<<"]="<<kReg[i]<<std::endl;
  }
  //assert(false);
  return;
}



void draw_di_sv_canv(TCanvas* di_sv_canv, TH1D* hSVal, TH1D* hdi, int kRegInput){
  
  bool funcDebug=true;
  if(funcDebug)std::cout<<std::endl<<"drawing singular values on di_sv_canv canvas.."<<std::endl<<std::endl;
  di_sv_canv->cd(1);
  di_sv_canv->cd(1)->SetLogy();  
  
  //hSVal->SetAxisRange(1.,(double)(hSVal->GetNbinsX()-1),"X");
  //hSVal->SetAxisRange(1.,(double)(hSVal->GetNbinsX()),"X");
  hSVal->SetTitle("Singular Values (AC^{-1})");        
  hSVal->SetXTitle("index i");        
  hSVal->SetYTitle("s_{i}");        
  hSVal->DrawClone("HIST E");
  
  //double tau=hSVal->GetBinContent(kRegInput+1)*hSVal->GetBinContent(kRegInput+1);
  double tau=hSVal->GetBinContent(kRegInput)*hSVal->GetBinContent(kRegInput);
  std::cout<<"hSVal->GetBinContent("<<kRegInput<<")="<<hSVal->GetBinContent(kRegInput)<<std::endl;
  std::cout<<"tau                     ="<<tau<<std::endl;
  
  if(funcDebug){
    printf("(orig)tau=%f\n",tau);
    tau*=100.;
    printf("tau=%f\n",tau);
    tau=(int)tau;
    printf("tau=%f\n",tau);
    tau/=100;
    printf("tau=%f\n",tau);
    std::cout<<std::endl;
  }
  
  float x=0.550173, y=0.8459761;
  drawText( "5.02 TeV ak4PFJets",                                 x, y, 19);y-=0.03;
  drawText( "2015 Prompt Reco"  ,                                 x, y, 19);y-=0.03;
  drawText( MCdesc.c_str()      ,                                 x, y, 19);y-=0.03;
  drawText( ("Current kReg="+std::to_string(kRegInput)).c_str() , x, y, 19);y-=0.03;	
  drawText( ("#tau = "+std::to_string( tau ) ).c_str() , x, y, 19);	      
  
  // di vector values
  di_sv_canv->cd(2);
  di_sv_canv->cd(2)->SetLogy(1);    
  
  //hdi->SetAxisRange(1.,(double)(hdi->GetNbinsX()-1),"X");
  //hdi->SetAxisRange(1.,(double)(hdi->GetNbinsX()),"X");
  hdi->SetTitle("Divector Values (#||{d_{i}}) ");
  hdi->SetXTitle("index i");
  hdi->SetYTitle("#||{d_{i}}");
  hdi->DrawClone("HIST E"); 
  
  ((TPad*)(di_sv_canv->cd(2)))->Update();      //note, if you do this to a pad that is log scaled, this returns the log of the max+min, not the max+min
  double ymax= gPad->GetUymax();
  ymax=pow(10,ymax);
  ((TPad*)(di_sv_canv->cd(2)))->Update();            
  double ymin= gPad->GetUymin();
  ymin=pow(10,ymin);
  
  di_sv_canv->cd(2);      
  hdi->DrawClone("HIST E"); 
  
  double xcoord= ( ((double)kRegInput) - 0.5);	
  if(funcDebug){std::cout<<"ymax="<<ymax<<std::endl;
    std::cout<<"ymin="<<ymin<<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"hdi->GetNbinsX()="<<hdi->GetNbinsX()<<std::endl;
  std::cout<<"kRegInput="<<kRegInput<<std::endl;
  //std::cout<<"((double)kRegInput) = "<< (((double)kRegInput)) <<std::endl;
  std::cout<<"xcoord="<<xcoord<<std::endl;
  }
  
  TLine* theLineAtOne_hdi=new TLine(0., 1., (double)(hdi->GetNbinsX()), 1.);
  theLineAtOne_hdi->SetLineWidth(1);
  theLineAtOne_hdi->SetLineStyle(2);
  theLineAtOne_hdi->SetLineColor(36);
  theLineAtOne_hdi->Draw();
  
  TLine* kRegLine_hdi=new TLine(xcoord, ymin, xcoord, ymax);      
  kRegLine_hdi->SetLineWidth(1);
  kRegLine_hdi->SetLineStyle(2);
  kRegLine_hdi->SetLineColor(36);
  kRegLine_hdi->Draw();
  
  float x1=0.550173, y1=0.8459761;
  drawText( "5.02 TeV ak4PFJets",                                 x1, y1, 19);y1-=0.03;
  drawText( "2015 Prompt Reco"  ,                                 x1, y1, 19);y1-=0.03;
  drawText( MCdesc.c_str()      ,                                 x1, y1, 19);y1-=0.03;
  drawText( ("Current kReg="+std::to_string(kRegInput)).c_str() , x1, y1, 19);y1-=0.03;	
  drawText( ("#tau="+std::to_string(tau) ).c_str() , x1, y1, 19);	
  
  float hdi_signif_mean=0;
  int signif_count=0;
  float hdi_insignif_mean=0;
  int insignif_count=0;
  for(int k=1; k<= (hdi->GetNbinsX());k++){
    if(funcDebug){
      std::cout<<"k="<<k<<std::endl;
      std::cout<<"x range is; "<<hdi->GetBinLowEdge(k)<<" to "<< hdi->GetBinLowEdge(k)+1<<std::endl;
      std::cout<<"hdi->GetBinContent("<<k<<")="<<hdi->GetBinContent(k)<<std::endl;
      std::cout<<"hSVal->GetBinContent("<<k<<")="<<hSVal->GetBinContent(k)<<std::endl;    }
    int lowedge = (int)hdi->GetBinLowEdge(k);
    if(lowedge<1){
      if(funcDebug)std::cout<<"k="<<k<<", GARBAGE!!!! SKIP!!!"<<std::endl;
      continue;
    }
    else if(lowedge>=1 && lowedge <= kRegInput){
      if(funcDebug)std::cout<<"k="<<k<<", SIGNIFICANT!!!!"<<std::endl;
      signif_count++;
      hdi_signif_mean+=hdi->GetBinContent(k);
    }
    else{
      if(funcDebug)std::cout<<"k="<<k<<", INSIGNIFICANT!!!!"<<std::endl;
      insignif_count++;
      hdi_insignif_mean+=hdi->GetBinContent(k);
    }
  }
  hdi_signif_mean/=  (float)signif_count;
  hdi_insignif_mean/=  (float)insignif_count;
	
  float hdi_signif_stddev=0;
  float hdi_insignif_stddev=0;
  for(int k=1; k<= (hdi->GetNbinsX());k++){
    int lowedge = (int)hdi->GetBinLowEdge(k);
    if(lowedge<1){
      std::cout<<"k="<<k<<", GARBAGE!!!! SKIP!!!"<<std::endl;
      continue;	}
    else if(lowedge>=1 && lowedge <= kRegInput){
      std::cout<<"k="<<k<<", SIGNIFICANT!!!!"<<std::endl;
      hdi_signif_stddev+=(hdi->GetBinContent(k)-hdi_signif_mean)*(hdi->GetBinContent(k)-hdi_signif_mean);	}
    else{
      std::cout<<"k="<<k<<", INSIGNIFICANT!!!!"<<std::endl;
      hdi_insignif_stddev+=(hdi->GetBinContent(k)-hdi_insignif_mean)*(hdi->GetBinContent(k)-hdi_insignif_mean);	}
  }
  
  hdi_signif_stddev/=(float)(signif_count -1 );
  hdi_signif_stddev=sqrt(hdi_signif_stddev);
  hdi_insignif_stddev/=(float)(insignif_count -1 );
  hdi_insignif_stddev=sqrt(hdi_insignif_stddev);
  
  std::cout<<"fo divectors for significant principal components of AC^-1's SVD decomposition"<<std::endl;
  std::cout<<"mean = "<<hdi_signif_mean<<std::endl;
  std::cout<<"std dev = "<<hdi_signif_stddev<<std::endl;
  
  std::cout<<"fo divectors for INsignificant principal components of AC^-1's SVD decomposition"<<std::endl;
  std::cout<<"mean = "<<hdi_insignif_mean<<std::endl;
  std::cout<<"std dev = "<<hdi_insignif_stddev<<std::endl;
  

  
  return;
}




//void makeMissesHist(TH2D* hmat_rebin=NULL, TH1D* hgen_rebin=NULL, TH1D* misseshist=NULL){
TH1D* makeMissesHist(TH2D* hmat_rebin=NULL, TH1D* hgen_rebin=NULL){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in makeMissesHist!"<<std::endl;
  if( !((bool)hmat_rebin) || !((bool)hgen_rebin) )assert(false);
  TH1D* misseshist=(TH1D*)hgen_rebin->Clone("my_misses_hist");
  misseshist->Reset("MICES");
  int nbinsy=hgen_rebin->GetNbinsX();//should also be the # of bins of resp matrix
                                     //#xbins in hgen == #y bins in hmat_rebin
  TAxis* hmat_rebin_yax=(TAxis*)hmat_rebin->GetYaxis();
  int nbinsy_hmat=hmat_rebin_yax->GetNbins();
  assert(nbinsy_hmat=nbinsy);
  TAxis* hmat_rebin_xax=(TAxis*)hmat_rebin->GetXaxis();
  int nbinsx=hmat_rebin_xax->GetNbins();
  
  for(int ybin=1; ybin<=nbinsy;ybin++){
    double val   =hgen_rebin->GetBinContent(ybin);
    double valerr=hgen_rebin->GetBinError(ybin);
    valerr*=valerr;
    double histval=val;
    double rowsum=0.;
    for(int xbin=1; xbin<=nbinsx;xbin++){
      rowsum+=hmat_rebin->GetBinContent(xbin,ybin);
      val-=hmat_rebin->GetBinContent(xbin,ybin);
      valerr-=(hmat_rebin->GetBinError(xbin,ybin)*hmat_rebin->GetBinError(xbin,ybin));    }
    
    if(funcDebug){
      std::cout<<"---------------------------------------------------"<<std::endl;
      std::cout<<"ybin="<<ybin<<", "<<(int)misseshist->GetBinLowEdge(ybin)<<" GeV < pT < "<<
	(int)misseshist->GetBinLowEdge(ybin) + (int)misseshist->GetBinWidth(ybin)<<std::endl;}
    
    if(val<0.)      { 
      if(funcDebug)std::cout<<"WARNING, val<0.! val="<<val<<", setting equal to 0."<<std::endl;
      val=0.;        }
    if(valerr>0.) 
      valerr=sqrt(valerr);
    else             {  
      if(funcDebug)std::cout<<"WARNING, valerr<=0.! valerr="<<valerr<<", setting equal to sqrt(val)="<<sqrt(val)<<std::endl;
      valerr=sqrt(val);    }
    
    if(funcDebug){
      std::cout<<"misseshist content = "<< val<< " +/- "<<valerr<<std::endl;
      std::cout<<"hgen content      = "<<histval<<std::endl;
      std::cout<<"hmat content      = "<<rowsum<<std::endl;
      std::cout<<"hgen hmat diff    = "<<histval-rowsum<<std::endl;
      
    }
    
    misseshist->SetBinContent(ybin, val);
    misseshist->SetBinError(ybin, valerr);
  }
  misseshist->SetMarkerStyle(kOpenStar);
  misseshist->SetMarkerColor(kBlack);
  misseshist->SetMarkerSize(1.02);
  misseshist->SetLineColor(kBlack);
  //return;
  return (TH1D*)misseshist;
}

//void makeFakesHist(TH2D* hmat_rebin=NULL, TH1D* hrec_sameside_rebin=NULL, TH1D* fakeshist=NULL){

TH1D* makeFakesHist(TH2D* hmat_rebin=NULL, TH1D* hrec_sameside_rebin=NULL){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in makeFakesHist!"<<std::endl;
  if( !((bool)hmat_rebin) || !((bool)hrec_sameside_rebin) )assert(false);
  TH1D* fakeshist=(TH1D*)hrec_sameside_rebin->Clone("my_fake_hist");
  fakeshist->Reset("MICES");
  int nbinsx=hrec_sameside_rebin->GetNbinsX();//should also be the # of bins of resp matrix
  TAxis* hmat_rebin_xax=(TAxis*)hmat_rebin->GetXaxis();
  int nbinsx_hmat=hmat_rebin_xax->GetNbins();
  assert(nbinsx_hmat=nbinsx);
  TAxis* hmat_rebin_yax=(TAxis*)hmat_rebin->GetYaxis();
  int nbinsy=hmat_rebin_yax->GetNbins();
  
  for(int xbin=1; xbin<=nbinsx;xbin++){
    double val   =hrec_sameside_rebin->GetBinContent(xbin);
    double valerr=hrec_sameside_rebin->GetBinError(xbin);
    valerr*=valerr;
    double histval=val;
    double colsum=0.;
    for(int ybin=1; ybin<=nbinsy;ybin++){
      colsum+=hmat_rebin->GetBinContent(xbin,ybin);
      val-=hmat_rebin->GetBinContent(xbin,ybin);
      valerr-=(hmat_rebin->GetBinError(xbin,ybin)*hmat_rebin->GetBinError(xbin,ybin));    }
    
    if(funcDebug){
      std::cout<<"---------------------------------------------------"<<std::endl;
      std::cout<<"xbin="<<xbin<<", "<<(int)fakeshist->GetBinLowEdge(xbin)<<" GeV < pT < "<<
	(int)fakeshist->GetBinLowEdge(xbin) + (int)fakeshist->GetBinWidth(xbin)<<std::endl;}
    
    if(val<0.)      { 
      if(funcDebug)std::cout<<"WARNING, val<0.! val="<<val<<", setting equal to 0."<<std::endl;
      val=0.;        }
    if(valerr>0.) 
      valerr=sqrt(valerr);
    else             {  
      if(funcDebug)std::cout<<"WARNING, valerr<=0.! valerr="<<valerr<<", setting equal to sqrt(val)="<<sqrt(val)<<std::endl;
      valerr=sqrt(val);    }
    
    if(funcDebug){
      std::cout<<"fakeshist content = "<< val<< " +/- "<<valerr<<std::endl;      
      std::cout<<"hrec content      = "<<histval<<std::endl;
      std::cout<<"hmat content      = "<<colsum<<std::endl;
      std::cout<<"hrec hmat diff    = "<<histval-colsum<<std::endl;
    }
    
    fakeshist->SetBinContent(xbin, val);
    fakeshist->SetBinError(xbin, valerr);
  }
  fakeshist->SetMarkerStyle(kOpenStar);
  fakeshist->SetMarkerColor(kGreen);
  fakeshist->SetMarkerSize(1.02);
  fakeshist->SetLineColor(kGreen);
  //return;
  return (TH1D*)fakeshist;
}



void applyFakeCorrection(TH1D* hrec_rebin_fakecorr         =NULL, 			 TH1D* hrec_sameside_rebin_fakecorr=NULL, 
			 TH1D* hfold_fakecorr              =NULL, 			 TH1D* hfold_truth_fakecorr        =NULL,
			 TH1D* hfak=NULL){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in applyFakeCorrection"<<std::endl;
  if(!(hrec_rebin_fakecorr) ||
     !(hrec_sameside_rebin_fakecorr) ||
     !(hfak) ){
    std::cout<<"not enough hists given for fake correction. exiting applyFakeCorrection"<<std::endl;
    return;  }
  
  //normalization diff between data meas and MC meas
  Double_t fac=(hrec_rebin_fakecorr->Integral()/hrec_sameside_rebin_fakecorr->Integral());
  int nbins_fakcorr=hrec_rebin_fakecorr->GetNbinsX();
  for(int i=1;i<=nbins_fakcorr; i++){
    //get # of fakes in given bin
    Double_t hfak_i     = hfak->GetBinContent(i); Double_t hfak_i_err = hfak->GetBinError(i);
    
    //take away fakes in misses spectra + set error
    Double_t hrec_i     = hrec_rebin_fakecorr->GetBinContent(i); Double_t hrec_i_err = hrec_rebin_fakecorr->GetBinError(i);
    hrec_rebin_fakecorr->SetBinContent( i, hrec_i-fac*hfak_i );
    hrec_rebin_fakecorr->SetBinError(   i, sqrt(hrec_i_err*hrec_i_err + fac*hfak_i_err*fac*hfak_i_err));
    
    //fold(unf) == detector level spectra w/o fakes. put the fakes back in
    if(hfold_fakecorr){
      Double_t hfold_i     = hfold_fakecorr->GetBinContent(i); Double_t hfold_i_err = hfold_fakecorr->GetBinError(i);
      hfold_fakecorr->SetBinContent( i, hfold_i+fac*hfak_i );
      hfold_fakecorr->SetBinError(   i, sqrt(hfold_i_err*hfold_i_err + fac*hfak_i_err*fac*hfak_i_err));
    }
    
    //take away fakes in misses spectra + set error
    Double_t hrecss_i     = hrec_sameside_rebin_fakecorr->GetBinContent(i); Double_t hrecss_i_err = hrec_sameside_rebin_fakecorr->GetBinError(i);
    hrec_sameside_rebin_fakecorr->SetBinContent( i, hrecss_i-hfak_i );
    hrec_sameside_rebin_fakecorr->SetBinError(   i, sqrt(hrecss_i_err*hrecss_i_err + hfak_i_err*hfak_i_err));
    
    //fold(truth) == detector level spectra w/o fakes. put the fakes back in
    if(hfold_truth_fakecorr){
      Double_t hfold_truth_i     = hfold_truth_fakecorr->GetBinContent(i); Double_t hfold_truth_i_err = hfold_truth_fakecorr->GetBinError(i);
      hfold_truth_fakecorr->SetBinContent( i, hfold_truth_i+hfak_i );
      hfold_truth_fakecorr->SetBinError(   i, sqrt(hfold_truth_i_err*hfold_truth_i_err + hfak_i_err*hfak_i_err));
    }
  }

  return;
}




//miss correction does NOT work in an analagous way to fake correction....
void applyMissCorrection(TH1D* hunf_misscorr         =NULL, TH1D* hgen_rebin=NULL,
			 TH1D* hmiss=NULL){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in applyMissCorrection"<<std::endl;
  if(!(hunf_misscorr) ||
     !(hgen_rebin) ||
     !(hmiss) ){
    std::cout<<"not enough hists given for Miss correction. exiting applyMissCorrection"<<std::endl;
    return;  }
  


  //normalization diff between data meas and MC meas
  //Double_t fac=(hunf_misscorr->Integral()/hgen_rebin->Integral());//this might not make sense since the correct normalization factor would be with the fakes in the spectra already..
  int nbins_fakcorr=hunf_misscorr->GetNbinsX();
  for(int i=1;i<=nbins_fakcorr; i++){
    //get # of fakes in given bin
    Double_t hmiss_i     = hmiss->GetBinContent(i);     
    //Double_t hmiss_i_err = hmiss->GetBinError(i);//forget error calculation for now
    Double_t eff_i=(hgen_rebin->GetBinContent(i)-hmiss_i)/hgen_rebin->GetBinContent(i);
    //Double_t eff_i_err=(hgen_rebin->GetBinError(i)*hgen_rebin->GetBinError(i) - hmiss_i_err*hmiss_i_err);
    
    //add back into unfolded data the misses spectra + set error
    Double_t hunf_i     = hunf_misscorr->GetBinContent(i); Double_t hunf_i_err = hunf_misscorr->GetBinError(i);
    hunf_misscorr->SetBinContent( i, hunf_i*eff_i );
    hunf_misscorr->SetBinError(   i, hunf_i_err*eff_i);

  }
  //
  ////normalization diff between data meas and MC meas
  //Double_t fac=(hunf_misscorr->Integral()/hgen_rebin->Integral());//this might not make sense since the correct normalization factor would be with the fakes in the spectra already..
  //int nbins_fakcorr=hunf_misscorr->GetNbinsX();
  //for(int i=1;i<=nbins_fakcorr; i++){
  //  //get # of fakes in given bin
  //  Double_t hmiss_i     = hmiss->GetBinContent(i);     Double_t hmiss_i_err = hmiss->GetBinError(i);
  //  
  //  //add back into unfolded data the misses spectra + set error
  //  Double_t hunf_i     = hunf_misscorr->GetBinContent(i); Double_t hunf_i_err = hunf_misscorr->GetBinError(i);
  //  hunf_misscorr->SetBinContent( i, hunf_i+fac*hmiss_i );
  //  hunf_misscorr->SetBinError(   i, sqrt(hunf_i_err*hunf_i_err + fac*hmiss_i_err*fac*hmiss_i_err));
  //
  //}

  return;
}








//this funciton makes the ratios itself explicitly.
void multiratioplot( TCanvas* c=NULL, TLegend * legend=NULL, 
		     std::vector<TH1D*> hvect={0}, 
		     TH1D* hdenom=NULL, std::string ratiotitle=""
		     ) {
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in multiratio plot v1"<<std::endl;
  if(!c)    return;
  if(!hdenom)    return;
  if(hvect.size()==0) return;
  if(ratiotitle.length()==0)
    ratiotitle="ratio with"+(std::string)hdenom->GetName();
  

  c->cd();
  
  //setup pads + draw on canvases
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetLogx(1);    pad1->SetLogy(1);
  pad1->SetGridx(0);   pad1->SetGridy(0);       
  pad1->SetBottomMargin(0); // Upper and lower plot are joined  
  pad1->Draw();             // Draw the upper pad: pad1
  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLogx(1);    pad2->SetLogy(0);  
  pad2->SetGridx(0); pad2->SetGridy(0); 
  pad2->SetTopMargin(0);  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  
  //draw spectra on pad1
  pad1->cd();                 
  for(unsigned int i=0;i<hvect.size();i++){
    if(i==0)hvect[i]->DrawClone("P E");
    else hvect[i]->DrawClone("P E SAME");    }  
  hdenom->DrawClone("P E SAME");  
  
  if( ((bool)legend))legend->Draw();
  
  //draw ratio(s) on pad 2  
  pad2->cd();       
  std::vector<TH1D*> hvect_ratio;
  for(unsigned int i=0;i<hvect.size(); i++){
    hvect_ratio.push_back((TH1D*)hvect[i]->Clone());
    hvect_ratio[i]->Divide(hdenom);
    if(i!=0)
      hvect_ratio[i]->DrawClone("P E SAME");  
    else{ //if(i==0){
      hvect_ratio[i]->Draw("P E");      
      TLine* lineatone=new TLine ( hdenom->GetBinLowEdge(1),                                                             1.,
				   hdenom->GetBinLowEdge(hdenom->GetNbinsX())+hdenom->GetBinWidth(hdenom->GetNbinsX()) , 1.);
      lineatone->SetLineStyle(2);      lineatone->SetLineWidth(1);      lineatone->SetLineColor(36);
      lineatone->Draw();  
      TLine* lineA=(TLine*)(lineatone->Clone());
      lineA->SetY1(1.25);  lineA->SetY2(1.25);      lineA->Draw();  
      TLine* lineB=(TLine*)(lineatone->Clone());
      lineB->SetY1(0.75);  lineB->SetY2(0.75);            lineB->Draw();  
      hvect_ratio[i]->Draw("P E SAME");    
    }  }
  
  // Ratio plot settings
  hvect_ratio[0]->SetTitle(""); // Remove the ratio title
  hvect_ratio[0]->SetAxisRange(0.6,1.4,"Y");
  
  // Y axis ratio plot settings
  hvect_ratio[0]->GetYaxis()->SetTitle(ratiotitle.c_str());
  hvect_ratio[0]->GetYaxis()->SetTitleSize(20);
  hvect_ratio[0]->GetYaxis()->SetTitleFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetYaxis()->SetTitleOffset(1.55);
  hvect_ratio[0]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetYaxis()->SetLabelSize(15);
  
  // X axis ratio plot settings
  hvect_ratio[0]->GetXaxis()->SetTitleSize(20);
  hvect_ratio[0]->GetXaxis()->SetTitleFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetXaxis()->SetTitleOffset(4.);
  hvect_ratio[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetXaxis()->SetLabelSize(15);
  
  return;
}


////this function takes in a vector of ratio histograms already made. This may be desired because perhaps the ratios are between different histograms, rather than having a common denominator
void multiratioplot( TCanvas* c=NULL, TLegend * legend=NULL, 
		     std::vector<TH1D*> hvect={0}, 
		     std::vector<TH1D*> hvect_ratio={0}, 
		     std::string ratiotitle=""
		     ) {
  bool funcDebug=true;
  if(funcDebug)std::cout<<"in multiratio plot v2"<<std::endl;
  if(!c)    return;
  if(hvect.size()==0) return;
  if(hvect_ratio.size()==0) return;
  if(ratiotitle.length()==0)
    ratiotitle="ratio";// with"+(std::string)hdenom->GetName();
  
  c->cd();
  
  //setup pads + draw on canvases
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetLogx(1);    pad1->SetLogy(1);
  pad1->SetGridx(0);   pad1->SetGridy(0);       
  pad1->SetBottomMargin(0); // Upper and lower plot are joined  
  pad1->Draw();             // Draw the upper pad: pad1
  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLogx(1);    pad2->SetLogy(0);  
  pad2->SetGridx(0); pad2->SetGridy(0); 
  pad2->SetTopMargin(0);  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  
  //draw spectra on pad1
  pad1->cd();                 
  for(unsigned int i=0;i<hvect.size();i++){
    if(i==0)hvect[i]->DrawClone("P E");
    else hvect[i]->DrawClone("P E SAME");    }  
  if( ((bool)legend))legend->Draw();
  
  //draw ratio(s) on pad 2  
  pad2->cd();       
  for(unsigned int i=0;i<hvect_ratio.size(); i++){
    if(i!=0)
      hvect_ratio[i]->DrawClone("P E SAME");  
    else{ //if(i==0){
      hvect_ratio[i]->Draw("P E");      
      TLine* lineatone=new TLine ( hvect_ratio[i]->GetBinLowEdge(1), 1.,
				   hvect_ratio[i]->GetBinLowEdge(hvect_ratio[i]->GetNbinsX())+hvect_ratio[i]->GetBinWidth(hvect_ratio[i]->GetNbinsX()) , 1.);
      lineatone->SetLineStyle(2);      lineatone->SetLineWidth(1);      lineatone->SetLineColor(36);
      lineatone->Draw();  
      TLine* lineA=(TLine*)(lineatone->Clone());
      lineA->SetY1(1.25);  lineA->SetY2(1.25);      lineA->Draw();  
      TLine* lineB=(TLine*)(lineatone->Clone());
      lineB->SetY1(0.75);  lineB->SetY2(0.75);            lineB->Draw();  
      hvect_ratio[i]->Draw("P E SAME");    
    }  }
  
  // Ratio plot settings
  hvect_ratio[0]->SetTitle(""); // Remove the ratio title
  hvect_ratio[0]->SetAxisRange(0.6,1.4,"Y");
  
  // Y axis ratio plot settings
  hvect_ratio[0]->GetYaxis()->SetTitle(ratiotitle.c_str());
  hvect_ratio[0]->GetYaxis()->SetTitleSize(20);
  hvect_ratio[0]->GetYaxis()->SetTitleFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetYaxis()->SetTitleOffset(1.55);
  hvect_ratio[0]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetYaxis()->SetLabelSize(15);
  
  // X axis ratio plot settings
  hvect_ratio[0]->GetXaxis()->SetTitleSize(20);
  hvect_ratio[0]->GetXaxis()->SetTitleFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetXaxis()->SetTitleOffset(4.);
  hvect_ratio[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hvect_ratio[0]->GetXaxis()->SetLabelSize(15);
  
  return;
}





//void makeCombinedPlots(std::string outRootFile="", TCanvas* canvForPrint=NULL, std::string outPdfFile="", bool applyNPCorrs=false){
void makeCombinedPlots(std::string outRootFile="", TCanvas* canvForPrint=NULL, std::string outPdfFile="", bool compareToNLOThy=false, bool applyNPCorrs=false){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"makeCombinedPlots"<<std::endl;
  if(outRootFile.length()==0){std::cout<<"ERROR no outRootFile string. exit."<<std::endl;
    return; }
  if(!canvForPrint){std::cout<<"ERROR canvForPrint pointer not found. exit."<<std::endl;
    return; }
  if(outPdfFile.length()==0){std::cout<<"ERROR no outPdfFile string. exit."<<std::endl;
    return; }
  std::vector<std::string> desclines={ "#sqrt{s} = 5.02 TeV",
				       "ak4PF Jets",
				       "#||{y} < 2.0"};
  
  TFile* fout= new TFile(outRootFile.c_str(), "UPDATE");
  if(!fout){std::cout<<"ERROR output file pointer not found. exit."<<std::endl;
    return; }
  
  canvForPrint->Clear();


  // combined plots, I/O spectra ----------------------------------------------
  TH1D* hunf		     = (TH1D*)fout->Get("Data_unf");
  hunf->SetMarkerStyle(kOpenCircle);
  hunf->SetMarkerColor(kRed);
  hunf->SetLineColor(kRed);
  hunf->SetMarkerSize(1.02);     
  
  std::cout<<"hunf->GetTitle()="<<hunf->GetTitle()<<std::endl;
  std::string descString="";
  getUnfDetails(hunf, &descString);
  std::string methodString="";
  //bool SVDUnf=false, BayesUnf=false;
  
  if(descString.find("kReg")!=std::string::npos){
    methodString=", SVD ";
  }
  else if(descString.find("kIter")!=std::string::npos){
    methodString=", Bayesian ";
  }
  
  
  TH1D* hrec_rebin	     = (TH1D*)fout->Get("Data_meas");
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetLineColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02); 

  TH1D* hrec_sameside_rebin  = (TH1D*)fout->Get("MC_meas");
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin->SetLineColor(kBlue-3);     
  hrec_sameside_rebin->SetMarkerSize(1.02);   

  TH1D* hgen_rebin	     = (TH1D*)fout->Get("MC_truth");
  hgen_rebin->SetMarkerStyle(kOpenStar);
  hgen_rebin->SetMarkerColor(kMagenta);
  hgen_rebin->SetLineColor(kMagenta);
  hgen_rebin->SetMarkerSize(1.02);     
  
  TLegend* legend_in1=new TLegend(0.7,0.68,0.92,0.9);
  legend_in1->AddEntry(hrec_rebin	    , "Data Meas." , "lp");
  legend_in1->AddEntry(hunf		    , "Data Unf." , "lp");
  legend_in1->AddEntry(hrec_sameside_rebin   , "PY8 Meas." , "lp");
  legend_in1->AddEntry(hgen_rebin	    , "PY8 Truth" , "lp");
  legend_in1->SetFillStyle(0);
  legend_in1->SetBorderSize(0);
  
  TCanvas* canv_spectra_genratio=NULL;
  //  hrec_rebin->SetTitle("Bayesian, PY8 and Data Spectra");    
  hrec_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hrec_rebin,  hrec_sameside_rebin, hunf}),
		  hgen_rebin, "Ratio w/ PY8 Truth");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_genratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  
  TCanvas* canv_spectra_recratio=NULL;  
  //  hgen_rebin->SetTitle("Bayesian, PY8 and Data Spectra");
  hgen_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hgen_rebin,  hrec_sameside_rebin, hunf}),
		  hrec_rebin, "Ratio w/ Data Meas.");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_recratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  
  
  TCanvas* canv_spectra_unfratio=NULL;
  //  hgen_rebin->SetTitle("Bayesian, PY8 and Data Spectra");
  hgen_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hgen_rebin,  hrec_sameside_rebin, hrec_rebin}),
		  hunf, "Ratio w/ Data Unf.");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_unfratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();

  
  // combined plots, meas. spectra, folded spectra + fakes, ratios ----------------------------------------------
  TH1D* hfold_fakecorr       = (TH1D*)fout->Get("Data_foldfakcorr");//check
  hfold_fakecorr->SetMarkerStyle(kOpenCircle);
  hfold_fakecorr->SetMarkerColor(kGreen-5);
  hfold_fakecorr->SetLineColor(  kGreen-5);
  hfold_fakecorr->SetMarkerSize(1.02);
  
  TH1D* hfold_truth_fakecorr = (TH1D*)fout->Get("MC_truth_foldfakcorr");//check
  hfold_truth_fakecorr->SetMarkerStyle(kOpenStar);
  hfold_truth_fakecorr->SetMarkerColor(kGreen+3);
  hfold_truth_fakecorr->SetLineColor(  kGreen+3);
  hfold_truth_fakecorr->SetMarkerSize(1.02);     

  TH1D* h_foldratio_datafold = (TH1D*)fout->Get("ppData_Meas_Ratio_DataFoldpFakes");//check
  h_foldratio_datafold->SetMarkerStyle(kOpenCircle);
  h_foldratio_datafold->SetMarkerColor(kGreen-5);
  h_foldratio_datafold->SetLineColor(  kGreen-5);
  h_foldratio_datafold->SetMarkerSize(1.02);

  TH1D* h_foldratio_mcfold   = (TH1D*)fout->Get("ppData_Meas_Ratio_TruthFoldpFakes");//check
  h_foldratio_mcfold->SetMarkerStyle(kOpenStar);
  h_foldratio_mcfold->SetMarkerColor(kGreen+3);
  h_foldratio_mcfold->SetLineColor(  kGreen+3);
  h_foldratio_mcfold->SetMarkerSize(1.02);     

  TLegend* legendfold2 = new TLegend( 0.7,0.7,0.9,0.9);
  legendfold2->AddEntry(hfold_truth_fakecorr,"Fold(PY8 Truth) + Fakes",  "lp");
  legendfold2->AddEntry(hrec_sameside_rebin, "PY8 Meas."  ,  "lp");
  legendfold2->AddEntry(hfold_fakecorr ,     "Fold(Data Unf.) + Fakes" , "lp");
  legendfold2->AddEntry(hrec_rebin, "Data Meas." , "lp");    
  legendfold2->SetFillStyle(0);
  legendfold2->SetBorderSize(0);    
  
  TCanvas* canv_spectra_foldratio=NULL;
  //setupSpectraHist(hfold_truth_fakecorr);
  //hfold_truth_fakecorr->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  setupSpectraHist(hrec_sameside_rebin);
  hrec_sameside_rebin->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendfold2,
  		  ((std::vector<TH1D*>){hrec_sameside_rebin, hfold_truth_fakecorr,  hfold_fakecorr,  hrec_rebin}),
  		  ((std::vector<TH1D*>){h_foldratio_datafold, h_foldratio_mcfold}),
  		  "Ratio w/ Meas. ");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());       
  canv_spectra_foldratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  

  
  
  // combined plots, meas spectra - fakes, folded spectra, ratios ----------------------------------------------
  TH1D* hfold       = (TH1D*)fout->Get("Data_fold");//check
  hfold->SetMarkerStyle(kOpenCircle);
  hfold->SetMarkerColor(kGreen-5);
  hfold->SetLineColor(  kGreen-5);
  hfold->SetMarkerSize(1.02);

  TH1D* hrec_rebin_fakecorr = (TH1D*)fout->Get("Data_measfakcorr");//check
  hrec_rebin_fakecorr->SetMarkerStyle(kOpenCircle);
  hrec_rebin_fakecorr->SetMarkerColor(kBlue);     
  hrec_rebin_fakecorr->SetLineColor(kBlue);     
  hrec_rebin_fakecorr->SetMarkerSize(1.02);     
  
  TH1D* hfold_truth = (TH1D*)fout->Get("MC_truth_fold");//check
  hfold_truth->SetMarkerStyle(kOpenStar);
  hfold_truth->SetMarkerColor(kGreen+3);
  hfold_truth->SetLineColor(  kGreen+3);
  hfold_truth->SetMarkerSize(1.02);     
  
  TH1D* hrec_sameside_rebin_fakecorr = (TH1D*)fout->Get("MC_measfakcorr");//check
  hrec_sameside_rebin_fakecorr->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin_fakecorr->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin_fakecorr->SetLineColor(kBlue-3);     
  hrec_sameside_rebin_fakecorr->SetMarkerSize(1.02);       
  
  TH1D* h_recratio_oppfold = (TH1D*)fout->Get("ppData_MeasmFakes_Ratio_DataFold");//check
  h_recratio_oppfold->SetMarkerStyle(kOpenCircle);
  h_recratio_oppfold->SetMarkerColor(kGreen-5);
  h_recratio_oppfold->SetLineColor(  kGreen-5);
  h_recratio_oppfold->SetMarkerSize(1.02);
  
  TH1D* h_recratio_truthfold   = (TH1D*)fout->Get("ppData_MeasmFakes_Ratio_TruthFold");//check
  h_recratio_truthfold->SetMarkerStyle(kOpenStar);
  h_recratio_truthfold->SetMarkerColor(kGreen+3);
  h_recratio_truthfold->SetLineColor(  kGreen+3);
  h_recratio_truthfold->SetMarkerSize(1.02);     
  
  TLegend* legendfold = new TLegend( 0.7,0.7,0.9,0.9);
  legendfold->AddEntry(hfold_truth,"Fold(PY8 Truth)",  "lp");
  legendfold->AddEntry(hrec_sameside_rebin_fakecorr, "PY8 Meas. - Fakes"  ,  "lp");
  legendfold->AddEntry(hfold ,     "Fold(Data Unf.)" , "lp");
  legendfold->AddEntry(hrec_rebin_fakecorr, "Data Meas. - Fakes" , "lp");    
  legendfold->SetFillStyle(0);
  legendfold->SetBorderSize(0);    
  
  TCanvas* canv_spectra_foldratio2=NULL;
  setupSpectraHist(hrec_sameside_rebin_fakecorr);
  hrec_sameside_rebin_fakecorr->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendfold,
  		  ((std::vector<TH1D*>){hrec_sameside_rebin_fakecorr, hfold_truth,  hfold,  hrec_rebin_fakecorr}),
  		  ((std::vector<TH1D*>){h_recratio_oppfold, h_recratio_truthfold}),
  		  "Ratio w/ (Meas. - Fakes) ");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());       
  canv_spectra_foldratio2=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  

  

  TCanvas* canv_thy_spectra_1=NULL;      TCanvas* canv_thy_spectra_2=NULL;    
  if(compareToNLOThy){
    // thy spectra, ratios #1 etc. ----------------------------------------------  
    TH1D* CT10nlo	 = (TH1D*)fout->Get("NLO_CT10_NLO_R04_jtpt")	 ;
    CT10nlo->SetMarkerSize(0);
    CT10nlo->SetLineColor(kBlack);  
    TH1D* CT14nlo	 = (TH1D*)fout->Get("NLO_CT14_NLO_R04_jtpt")	 ;
    CT14nlo->SetMarkerSize(0);
    CT14nlo->SetLineColor(kGreen);
    TH1D* NNPDFnnlo= (TH1D*)fout->Get("NLO_NNPDF_NLO_R04_jtpt")    ;
    NNPDFnnlo->SetMarkerSize(0);
    NNPDFnnlo->SetLineColor(kCyan-6);  
    
    TLegend* legendThy1 =new TLegend( 0.6,0.7,1.0,0.9 );    
    if(applyNPCorrs){
      legendThy1->AddEntry(CT10nlo  ,"CT10 NLO #otimes HERWIG EE4C","l");
      legendThy1->AddEntry(CT14nlo  ,"CT14 NLO #otimes HERWIG EE4C","l");
      legendThy1->AddEntry(NNPDFnnlo,"NNPDF NNLO #otimes POW+PY8","l");   }
    else{
      legendThy1->AddEntry(CT10nlo  ,"CT10 NLO","l");
      legendThy1->AddEntry(CT14nlo  ,"CT14 NLO","l");
      legendThy1->AddEntry(NNPDFnnlo,"NNPDF NNLO","l");   }  
    legendThy1->AddEntry(hunf,"Data Unf.","lp");
    legendThy1->AddEntry(hgen_rebin,"PY8 Truth", "lp");
    legendThy1->SetFillStyle(0);
    legendThy1->SetBorderSize(0);    
    

    if(!applyNPCorrs)CT10nlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());    
    else CT10nlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());    
    multiratioplot( canvForPrint, legendThy1,
		    ((std::vector<TH1D*>){CT10nlo, CT14nlo, NNPDFnnlo, hgen_rebin}),
		    hunf, "Ratio w/ Data Unf.");    
    drawTLatex( 0.48, 0.88, 0.035, canvForPrint, desclines);
    canvForPrint->Print(outPdfFile.c_str()); 
    canv_thy_spectra_1=(TCanvas*)canvForPrint->DrawClone(); 
    canvForPrint->Clear();
    
    // thy spectra, ratios #2 etc. ----------------------------------------------  
    TH1D* HERAPDF	 = (TH1D*)fout->Get("NLO_HERAPDF105_NLO_R04_jtpt")	 ;
    HERAPDF->SetMarkerSize(0);
    HERAPDF->SetLineColor(kViolet-5);  
    
    TH1D* MMHTnlo	 = (TH1D*)fout->Get("NLO_MMHT2014_NLO_R04_jtpt")	 ;
    MMHTnlo->SetMarkerSize(0);
    MMHTnlo->SetLineColor(kOrange+7);    
    
    TLegend* legendThy2 =new TLegend( 0.6,0.7,1.0,0.9 );    
    if(applyNPCorrs){
      legendThy2->AddEntry(HERAPDF  ,"HEARPDF NLO #otimes HERWIG EE4C","l");
      legendThy2->AddEntry(MMHTnlo  ,"MMHT NLO #otimes HERWIG EE4C","l");}
    else{
      legendThy2->AddEntry(CT10nlo  ,"CT10 NLO","l");
      legendThy2->AddEntry(CT14nlo  ,"CT14 NLO","l");}
    legendThy2->AddEntry(hunf,"Data Unf.","lp");
    legendThy2->AddEntry(hgen_rebin,"PY8 Truth", "lp");
    legendThy2->SetFillStyle(0);
    legendThy2->SetBorderSize(0);    
    

    if(!applyNPCorrs)MMHTnlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());    
    else MMHTnlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());    
    multiratioplot( canvForPrint, legendThy2,
		    ((std::vector<TH1D*>){MMHTnlo, HERAPDF,  hgen_rebin}),
		    hunf, "Ratio w/ Data Unf.");    
    drawTLatex( 0.48, 0.88, 0.035, canvForPrint, desclines);
    canvForPrint->Print(outPdfFile.c_str()); 
    canv_thy_spectra_2=(TCanvas*)canvForPrint->DrawClone(); 
    canvForPrint->Clear();
  }

  
  fout->cd();  
  canv_spectra_genratio->SetTitle("PY8 Truth Ratios Combined Canvas");	      canv_spectra_genratio->Write("combcanv_spectra_genratio");		  
  canv_spectra_recratio->SetTitle("Data Meas Ratios Combined Canvas");	      canv_spectra_recratio->Write("combcanv_spectra_recratio");		  
  canv_spectra_unfratio->SetTitle("Data Unf. Ratios Combined Canvas");	      canv_spectra_unfratio->Write("combcanv_spectra_unfratio");		  
  canv_spectra_foldratio->SetTitle("Folding Test Add Fakes Ratios Combined Canvas");        canv_spectra_foldratio->Write("combcanv_spectra_pfakes_foldratio");    
  canv_spectra_foldratio2->SetTitle("Folding Test Minus Fakes Ratios Combiend Canvas");	    canv_spectra_foldratio2->Write("combcanv_spectra_mfakes_foldratio");

  if(compareToNLOThy){
    canv_thy_spectra_1->SetTitle("NLO Thy Spectra Combined Canvas");                          canv_thy_spectra_1->Write("combcanv_thy_spectra_1");     
    canv_thy_spectra_2->SetTitle("NLO Thy Spectra Combined Canvas");                          canv_thy_spectra_2->Write("combcanv_thy_spectra_2");          
  }
  
  return;
}

void makeCombinedPlots_NLO_etabin(
				  std::string outRootFile="", 
				  TCanvas* canvForPrint=NULL,
				  std::string outPdfFile="", 
				  bool applyNPCorrs=false, 
				  std::string NLOMCtitle="", 
				  int etabinint=0){
  
  bool funcDebug=true;
  if(funcDebug)std::cout<<"makeCombinedPlots"<<std::endl;
  if(outRootFile.length()==0){std::cout<<"ERROR no outRootFile string. exit."<<std::endl;
    return; }
  if(!canvForPrint){std::cout<<"ERROR canvForPrint pointer not found. exit."<<std::endl;
    return; }
  if(outPdfFile.length()==0){std::cout<<"ERROR no outPdfFile string. exit."<<std::endl;
    return; }
  
  
  std::string etabindesc_str="#||{y}";
  if(     etabinint==0)etabindesc_str="0.0 < "+etabindesc_str+" < 0.5";
  else if(etabinint==1)etabindesc_str="0.5 < "+etabindesc_str+" < 1.0";
  else if(etabinint==2)etabindesc_str="1.0 < "+etabindesc_str+" < 1.5";
  else if(etabinint==3)etabindesc_str="1.5 < "+etabindesc_str+" < 2.0";
  else if(etabinint==4)etabindesc_str="2.0 < "+etabindesc_str+" < 2.5";
  else if(etabinint==5)etabindesc_str="2.5 < "+etabindesc_str+" < 3.0";
  
  
  std::vector<std::string> desclines={ "#sqrt{s} = 5.02 TeV",
				       "ak4PF Jets",
				       etabindesc_str};
  
  TFile* fout= new TFile(outRootFile.c_str(), "UPDATE");
  if(!fout){std::cout<<"ERROR output file pointer not found. exit."<<std::endl;
    return; }
  
  canvForPrint->Clear();
  
  
  // combined plots, I/O spectra ----------------------------------------------
  TH1D* hunf		     = (TH1D*)fout->Get("Data_unf");
  hunf->SetMarkerStyle(kOpenCircle);
  hunf->SetMarkerColor(kRed);
  hunf->SetLineColor(kRed);
  hunf->SetMarkerSize(1.02);     
  
  std::cout<<"hunf->GetTitle()="<<hunf->GetTitle()<<std::endl;
  std::string descString="";
  getUnfDetails(hunf, &descString);
  std::string methodString="";
  //bool SVDUnf=false, BayesUnf=false;
  
  if(descString.find("kReg")!=std::string::npos){
    methodString=", SVD ";
  }
  else if(descString.find("kIter")!=std::string::npos){
    methodString=", Bayesian ";
  }
  
  
  TH1D* hrec_rebin	     = (TH1D*)fout->Get("Data_meas");
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetLineColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02); 
  
  TH1D* hrec_sameside_rebin  = (TH1D*)fout->Get("MC_meas");
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin->SetLineColor(kBlue-3);     
  hrec_sameside_rebin->SetMarkerSize(1.02);   
  
  TH1D* hgen_rebin	     = (TH1D*)fout->Get("MC_truth");
  hgen_rebin->SetMarkerStyle(kOpenStar);
  hgen_rebin->SetMarkerColor(kMagenta);
  hgen_rebin->SetLineColor(kMagenta);
  hgen_rebin->SetMarkerSize(1.02);     
  
  TLegend* legend_in1=new TLegend(0.7,0.68,0.92,0.9);
  legend_in1->AddEntry(hrec_rebin	    , "Data Meas." , "lp");
  legend_in1->AddEntry(hunf		    , "Data Unf." , "lp");
  legend_in1->AddEntry(hrec_sameside_rebin   , ("Smeared "+NLOMCtitle).c_str() , "lp");
  legend_in1->AddEntry(hgen_rebin	    , (NLOMCtitle+" Toy Truth").c_str() , "lp");
  legend_in1->SetFillStyle(0);
  legend_in1->SetBorderSize(0);
  
  TCanvas* canv_spectra_genratio=NULL;
  //  hrec_rebin->SetTitle("Bayesian, PY8 and Data Spectra");    
  hrec_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hrec_rebin,  hrec_sameside_rebin, hunf}),
		  hgen_rebin, "Ratio w/ "+NLOMCtitle+" Toy Truth");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_genratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  
  TCanvas* canv_spectra_recratio=NULL;  
  //  hgen_rebin->SetTitle("Bayesian, PY8 and Data Spectra");
  hgen_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hgen_rebin,  hrec_sameside_rebin, hunf}),
		  hrec_rebin, "Ratio w/ Data Meas.");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_recratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  
  
  TCanvas* canv_spectra_unfratio=NULL;
  //  hgen_rebin->SetTitle("Bayesian, PY8 and Data Spectra");
  hgen_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hgen_rebin,  hrec_sameside_rebin, hrec_rebin}),
		  hunf, "Ratio w/ Data Unf.");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_unfratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  
  
  // combined plots, meas. spectra, folded spectra + fakes, ratios ----------------------------------------------
  TH1D* hfold_fakecorr       = (TH1D*)fout->Get("Data_foldfakcorr");//check
  hfold_fakecorr->SetMarkerStyle(kOpenCircle);
  hfold_fakecorr->SetMarkerColor(kGreen-5);
  hfold_fakecorr->SetLineColor(  kGreen-5);
  hfold_fakecorr->SetMarkerSize(1.02);
  
  TH1D* hfold_truth_fakecorr = (TH1D*)fout->Get("MC_truth_foldfakcorr");//check
  hfold_truth_fakecorr->SetMarkerStyle(kOpenStar);
  hfold_truth_fakecorr->SetMarkerColor(kGreen+3);
  hfold_truth_fakecorr->SetLineColor(  kGreen+3);
  hfold_truth_fakecorr->SetMarkerSize(1.02);     

  TH1D* h_foldratio_datafold = (TH1D*)fout->Get("ppData_Meas_Ratio_DataFoldpFakes");//check
  h_foldratio_datafold->SetMarkerStyle(kOpenCircle);
  h_foldratio_datafold->SetMarkerColor(kGreen-5);
  h_foldratio_datafold->SetLineColor(  kGreen-5);
  h_foldratio_datafold->SetMarkerSize(1.02);

  TH1D* h_foldratio_mcfold   = (TH1D*)fout->Get("ppData_Meas_Ratio_TruthFoldpFakes");//check
  h_foldratio_mcfold->SetMarkerStyle(kOpenStar);
  h_foldratio_mcfold->SetMarkerColor(kGreen+3);
  h_foldratio_mcfold->SetLineColor(  kGreen+3);
  h_foldratio_mcfold->SetMarkerSize(1.02);     

  TLegend* legendfold2 = new TLegend( 0.7,0.7,0.9,0.9);
  legendfold2->AddEntry(hfold_truth_fakecorr,"Fold(NLO Toy Truth) + Fakes",  "lp");
  legendfold2->AddEntry(hrec_sameside_rebin, "Smeared NLO"  ,  "lp");
  legendfold2->AddEntry(hfold_fakecorr ,     "Fold(Data Unf.) + Fakes" , "lp");
  legendfold2->AddEntry(hrec_rebin, "Data Meas." , "lp");    
  legendfold2->SetFillStyle(0);
  legendfold2->SetBorderSize(0);    
  
  TCanvas* canv_spectra_foldratio=NULL;
  //setupSpectraHist(hfold_truth_fakecorr);
  //hfold_truth_fakecorr->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  setupSpectraHist(hrec_sameside_rebin);
  hrec_sameside_rebin->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendfold2,
  		  ((std::vector<TH1D*>){hrec_sameside_rebin, hfold_truth_fakecorr,  hfold_fakecorr,  hrec_rebin}),
  		  ((std::vector<TH1D*>){h_foldratio_datafold, h_foldratio_mcfold}),
  		  "Ratio w/ Meas. ");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());       
  canv_spectra_foldratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  

  
  
  // combined plots, meas spectra - fakes, folded spectra, ratios ----------------------------------------------
  TH1D* hfold       = (TH1D*)fout->Get("Data_fold");//check
  hfold->SetMarkerStyle(kOpenCircle);
  hfold->SetMarkerColor(kGreen-5);
  hfold->SetLineColor(  kGreen-5);
  hfold->SetMarkerSize(1.02);
  
  TH1D* hrec_rebin_fakecorr = (TH1D*)fout->Get("Data_measfakcorr");//check
  hrec_rebin_fakecorr->SetMarkerStyle(kOpenCircle);
  hrec_rebin_fakecorr->SetMarkerColor(kBlue);     
  hrec_rebin_fakecorr->SetLineColor(kBlue);     
  hrec_rebin_fakecorr->SetMarkerSize(1.02);     
  
  TH1D* hfold_truth = (TH1D*)fout->Get("MC_truth_fold");//check
  hfold_truth->SetMarkerStyle(kOpenStar);
  hfold_truth->SetMarkerColor(kGreen+3);
  hfold_truth->SetLineColor(  kGreen+3);
  hfold_truth->SetMarkerSize(1.02);     
  
  TH1D* hrec_sameside_rebin_fakecorr = (TH1D*)fout->Get("MC_measfakcorr");//check
  hrec_sameside_rebin_fakecorr->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin_fakecorr->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin_fakecorr->SetLineColor(kBlue-3);     
  hrec_sameside_rebin_fakecorr->SetMarkerSize(1.02);       
  
  TH1D* h_recratio_oppfold = (TH1D*)fout->Get("ppData_MeasmFakes_Ratio_DataFold");//check
  h_recratio_oppfold->SetMarkerStyle(kOpenCircle);
  h_recratio_oppfold->SetMarkerColor(kGreen-5);
  h_recratio_oppfold->SetLineColor(  kGreen-5);
  h_recratio_oppfold->SetMarkerSize(1.02);
  
  TH1D* h_recratio_truthfold   = (TH1D*)fout->Get("ppData_MeasmFakes_Ratio_TruthFold");//check
  h_recratio_truthfold->SetMarkerStyle(kOpenStar);
  h_recratio_truthfold->SetMarkerColor(kGreen+3);
  h_recratio_truthfold->SetLineColor(  kGreen+3);
  h_recratio_truthfold->SetMarkerSize(1.02);     
  
  TLegend* legendfold = new TLegend( 0.7,0.7,0.9,0.9);
  legendfold->AddEntry(hfold_truth,"Fold(Toy NLO Truth)",  "lp");
  legendfold->AddEntry(hrec_sameside_rebin_fakecorr, "Smeared NLO - Fakes"  ,  "lp");
  legendfold->AddEntry(hfold ,     "Fold(Data Unf.)" , "lp");
  legendfold->AddEntry(hrec_rebin_fakecorr, "Data Meas. - Fakes" , "lp");    
  legendfold->SetFillStyle(0);
  legendfold->SetBorderSize(0);    
  
  TCanvas* canv_spectra_foldratio2=NULL;
  setupSpectraHist(hrec_sameside_rebin_fakecorr);
  hrec_sameside_rebin_fakecorr->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendfold,
  		  ((std::vector<TH1D*>){hrec_sameside_rebin_fakecorr, hfold_truth,  hfold,  hrec_rebin_fakecorr}),
  		  ((std::vector<TH1D*>){h_recratio_oppfold, h_recratio_truthfold}),
  		  "Ratio w/ (Meas. - Fakes) ");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());       
  canv_spectra_foldratio2=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  

  

  TCanvas* canv_thy_spectra_1=NULL;      TCanvas* canv_thy_spectra_2=NULL;    

  // thy spectra, ratios #1 etc. ----------------------------------------------  
  TH1D* CT10nlo	 = (TH1D*)fout->Get("NLO_CT10_NLO_R04_jtpt")	 ;
  CT10nlo->SetMarkerSize(0);
  CT10nlo->SetLineColor(kBlack);  
  TH1D* CT14nlo	 = (TH1D*)fout->Get("NLO_CT14_NLO_R04_jtpt")	 ;
  CT14nlo->SetMarkerSize(0);
  CT14nlo->SetLineColor(kGreen);
  TH1D* NNPDFnnlo= (TH1D*)fout->Get("NLO_NNPDF_NLO_R04_jtpt")    ;
  NNPDFnnlo->SetMarkerSize(0);
  NNPDFnnlo->SetLineColor(kCyan-6);  
  
  TLegend* legendThy1 =new TLegend( 0.6,0.7,1.0,0.9 );    
  if(applyNPCorrs){
    legendThy1->AddEntry(CT10nlo  ,"CT10 NLO #otimes HERWIG EE4C","l");
    legendThy1->AddEntry(CT14nlo  ,"CT14 NLO #otimes HERWIG EE4C","l");
    legendThy1->AddEntry(NNPDFnnlo,"NNPDF NNLO #otimes POW+PY8","l");   }
  else{
    legendThy1->AddEntry(CT10nlo  ,"CT10 NLO","l");
    legendThy1->AddEntry(CT14nlo  ,"CT14 NLO","l");
    legendThy1->AddEntry(NNPDFnnlo,"NNPDF NNLO","l");   }  
  legendThy1->AddEntry(hunf,"Data Unf.","lp");
  legendThy1->AddEntry(hgen_rebin,(NLOMCtitle+" Toy Truth").c_str(), "lp");
  legendThy1->SetFillStyle(0);
  legendThy1->SetBorderSize(0);    
    
  
  if(!applyNPCorrs)CT10nlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());    
  else CT10nlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendThy1,
		  ((std::vector<TH1D*>){CT10nlo, CT14nlo, NNPDFnnlo, hgen_rebin}),
		  hunf, "Ratio w/ Data Unf.");    
  drawTLatex( 0.48, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str()); 
  canv_thy_spectra_1=(TCanvas*)canvForPrint->DrawClone(); 
  canvForPrint->Clear();
  
  // thy spectra, ratios #2 etc. ----------------------------------------------  
  TH1D* HERAPDF	 = (TH1D*)fout->Get("NLO_HERAPDF105_NLO_R04_jtpt")	 ;
  HERAPDF->SetMarkerSize(0);
  HERAPDF->SetLineColor(kViolet-5);  
  
  TH1D* MMHTnlo	 = (TH1D*)fout->Get("NLO_MMHT2014_NLO_R04_jtpt")	 ;
  MMHTnlo->SetMarkerSize(0);
  MMHTnlo->SetLineColor(kOrange+7);    
  
  TLegend* legendThy2 =new TLegend( 0.6,0.7,1.0,0.9 );    
  if(applyNPCorrs){
    legendThy2->AddEntry(HERAPDF  ,"HERAPDF NLO #otimes HERWIG EE4C","l");
    legendThy2->AddEntry(MMHTnlo  ,"MMHT NLO #otimes HERWIG EE4C","l");}
  else{
    legendThy2->AddEntry(CT10nlo  ,"CT10 NLO","l");
    legendThy2->AddEntry(CT14nlo  ,"CT14 NLO","l");}
  legendThy2->AddEntry(hunf,"Data Unf.","lp");
  legendThy2->AddEntry(hgen_rebin,(NLOMCtitle+" Toy Truth").c_str(), "lp");
  legendThy2->SetFillStyle(0);
  legendThy2->SetBorderSize(0);    
  
  
  if(!applyNPCorrs)MMHTnlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());    
  else MMHTnlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendThy2,
		    ((std::vector<TH1D*>){MMHTnlo, HERAPDF,  hgen_rebin}),
		  hunf, "Ratio w/ Data Unf.");    
  drawTLatex( 0.48, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str()); 
  canv_thy_spectra_2=(TCanvas*)canvForPrint->DrawClone(); 
  canvForPrint->Clear();
  
  
  
  fout->cd();  
  canv_spectra_genratio->SetTitle((NLOMCtitle+" Toy Truth Ratios Combined Canvas").c_str());	      canv_spectra_genratio->Write("combcanv_spectra_genratio");		  
  canv_spectra_recratio->SetTitle("Data Meas Ratios Combined Canvas");	      canv_spectra_recratio->Write("combcanv_spectra_recratio");		  
  canv_spectra_unfratio->SetTitle("Data Unf. Ratios Combined Canvas");	      canv_spectra_unfratio->Write("combcanv_spectra_unfratio");		  
  canv_spectra_foldratio->SetTitle("Folding Test Add Fakes Ratios Combined Canvas");        canv_spectra_foldratio->Write("combcanv_spectra_pfakes_foldratio");    
  canv_spectra_foldratio2->SetTitle("Folding Test Minus Fakes Ratios Combiend Canvas");	    canv_spectra_foldratio2->Write("combcanv_spectra_mfakes_foldratio");
  
  canv_thy_spectra_1->SetTitle("NLO Thy Spectra Combined Canvas");                          canv_thy_spectra_1->Write("combcanv_thy_spectra_1");     
  canv_thy_spectra_2->SetTitle("NLO Thy Spectra Combined Canvas");                          canv_thy_spectra_2->Write("combcanv_thy_spectra_2");          
    
  
  return;
}



void makeCombinedPlots_closure(std::string outRootFile="", TCanvas* canvForPrint=NULL, std::string outPdfFile="", bool applyNPCorrs=false){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"makeCombinedPlots_closure"<<std::endl;
  if(outRootFile.length()==0){std::cout<<"ERROR no outRootFile string. exit."<<std::endl;
    return; }
  if(!canvForPrint){std::cout<<"ERROR canvForPrint pointer not found. exit."<<std::endl;
    return; }
  if(outPdfFile.length()==0){std::cout<<"ERROR no outPdfFile string. exit."<<std::endl;
    return; }
  std::vector<std::string> desclines={ "#sqrt{s} = 5.02 TeV",
				       "ak4PF Jets",
				       "#||{y} < 2.0"};
  
  TFile* fout= new TFile(outRootFile.c_str(), "UPDATE");
  if(!fout){std::cout<<"ERROR output file pointer not found. exit."<<std::endl;
    return; }
  
  canvForPrint->Clear();


  // combined plots, I/O spectra ----------------------------------------------
  TH1D* hunf		     = (TH1D*)fout->Get("Data_unf");
  hunf->SetMarkerStyle(kOpenCircle);
  hunf->SetMarkerColor(kRed);
  hunf->SetLineColor(kRed);
  hunf->SetMarkerSize(1.02);     
  
  std::cout<<"hunf->GetTitle()="<<hunf->GetTitle()<<std::endl;
  std::string descString="";
  getUnfDetails(hunf, &descString);
  std::string methodString="";
  //bool SVDUnf=false, BayesUnf=false;
  
  if(descString.find("kReg")!=std::string::npos){
    methodString=", SVD ";
  }
  else if(descString.find("kIter")!=std::string::npos){
    methodString=", Bayesian ";
  }
  
  
  TH1D* hrec_rebin	     = (TH1D*)fout->Get("Data_meas");
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetLineColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02); 

  TH1D* hrec_sameside_rebin  = (TH1D*)fout->Get("MC_meas");
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin->SetLineColor(kBlue-3);     
  hrec_sameside_rebin->SetMarkerSize(1.02);   

  TH1D* hgen_rebin	     = (TH1D*)fout->Get("MC_truth");
  hgen_rebin->SetMarkerStyle(kOpenStar);
  hgen_rebin->SetMarkerColor(kMagenta);
  hgen_rebin->SetLineColor(kMagenta);
  hgen_rebin->SetMarkerSize(1.02);     
  
  TLegend* legend_in1=new TLegend(0.7,0.68,0.92,0.9);
  legend_in1->AddEntry(hrec_rebin	    , "Test PY8 Meas." , "lp");
  legend_in1->AddEntry(hunf		    , "Test PY8 Unf." , "lp");
  legend_in1->AddEntry(hrec_sameside_rebin   , "PY8 Meas." , "lp");
  legend_in1->AddEntry(hgen_rebin	    , "PY8 Truth" , "lp");
  legend_in1->SetFillStyle(0);
  legend_in1->SetBorderSize(0);
  
  TCanvas* canv_spectra_genratio=NULL;
  //  hrec_rebin->SetTitle("Bayesian, PY8 and Data Spectra");    
  hrec_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hrec_rebin,  hrec_sameside_rebin, hunf}),
		  hgen_rebin, "Ratio w/ PY8 Truth");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_genratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  
  TCanvas* canv_spectra_recratio=NULL;  
  //  hgen_rebin->SetTitle("Bayesian, PY8 and Data Spectra");
  hgen_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hgen_rebin,  hrec_sameside_rebin, hunf}),
		  hrec_rebin, "Ratio w/ Test PY8 Meas.");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_recratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  
  
  TCanvas* canv_spectra_unfratio=NULL;
  //  hgen_rebin->SetTitle("Bayesian, PY8 and Data Spectra");
  hgen_rebin->SetTitle(("Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legend_in1,
		  ((std::vector<TH1D*>){hgen_rebin,  hrec_sameside_rebin, hrec_rebin}),
		  hunf, "Ratio w/ Test PY8 Unf.");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());      
  canv_spectra_unfratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();

  
  // combined plots, meas. spectra, folded spectra + fakes, ratios ----------------------------------------------
  TH1D* hfold_fakecorr       = (TH1D*)fout->Get("Data_foldfakcorr");//check
  hfold_fakecorr->SetMarkerStyle(kOpenCircle);
  hfold_fakecorr->SetMarkerColor(kGreen-5);
  hfold_fakecorr->SetLineColor(  kGreen-5);
  hfold_fakecorr->SetMarkerSize(1.02);
  
  TH1D* hfold_truth_fakecorr = (TH1D*)fout->Get("MC_truth_foldfakcorr");//check
  hfold_truth_fakecorr->SetMarkerStyle(kOpenStar);
  hfold_truth_fakecorr->SetMarkerColor(kGreen+3);
  hfold_truth_fakecorr->SetLineColor(  kGreen+3);
  hfold_truth_fakecorr->SetMarkerSize(1.02);     

  TH1D* h_foldratio_datafold = (TH1D*)fout->Get("ppData_Meas_Ratio_DataFoldpFakes");//check
  h_foldratio_datafold->SetMarkerStyle(kOpenCircle);
  h_foldratio_datafold->SetMarkerColor(kGreen-5);
  h_foldratio_datafold->SetLineColor(  kGreen-5);
  h_foldratio_datafold->SetMarkerSize(1.02);

  TH1D* h_foldratio_mcfold   = (TH1D*)fout->Get("ppData_Meas_Ratio_TruthFoldpFakes");//check
  h_foldratio_mcfold->SetMarkerStyle(kOpenStar);
  h_foldratio_mcfold->SetMarkerColor(kGreen+3);
  h_foldratio_mcfold->SetLineColor(  kGreen+3);
  h_foldratio_mcfold->SetMarkerSize(1.02);     

  TLegend* legendfold2 = new TLegend( 0.7,0.7,0.9,0.9);
  legendfold2->AddEntry(hfold_truth_fakecorr,"Fold(PY8 Truth) + Fakes",  "lp");
  legendfold2->AddEntry(hrec_sameside_rebin, "PY8 Meas."  ,  "lp");
  legendfold2->AddEntry(hfold_fakecorr ,     "Fold(Test PY8 Unf.) + Fakes" , "lp");
  legendfold2->AddEntry(hrec_rebin, "Test PY8 Meas." , "lp");    
  legendfold2->SetFillStyle(0);
  legendfold2->SetBorderSize(0);    
  
  TCanvas* canv_spectra_foldratio=NULL;
  //setupSpectraHist(hfold_truth_fakecorr);
  //hfold_truth_fakecorr->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  setupSpectraHist(hrec_sameside_rebin);
  hrec_sameside_rebin->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendfold2,
  		  ((std::vector<TH1D*>){hrec_sameside_rebin, hfold_truth_fakecorr,  hfold_fakecorr,  hrec_rebin}),
  		  ((std::vector<TH1D*>){h_foldratio_datafold, h_foldratio_mcfold}),
  		  "Ratio w/ Meas.");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());       
  canv_spectra_foldratio=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  

  
  
  // combined plots, meas spectra - fakes, folded spectra, ratios ----------------------------------------------
  TH1D* hfold       = (TH1D*)fout->Get("Data_fold");//check
  hfold->SetMarkerStyle(kOpenCircle);
  hfold->SetMarkerColor(kGreen-5);
  hfold->SetLineColor(  kGreen-5);
  hfold->SetMarkerSize(1.02);

  TH1D* hrec_rebin_fakecorr = (TH1D*)fout->Get("Data_measfakcorr");//check
  hrec_rebin_fakecorr->SetMarkerStyle(kOpenCircle);
  hrec_rebin_fakecorr->SetMarkerColor(kBlue);     
  hrec_rebin_fakecorr->SetLineColor(kBlue);     
  hrec_rebin_fakecorr->SetMarkerSize(1.02);     
  
  TH1D* hfold_truth = (TH1D*)fout->Get("MC_truth_fold");//check
  hfold_truth->SetMarkerStyle(kOpenStar);
  hfold_truth->SetMarkerColor(kGreen+3);
  hfold_truth->SetLineColor(  kGreen+3);
  hfold_truth->SetMarkerSize(1.02);     
  
  TH1D* hrec_sameside_rebin_fakecorr = (TH1D*)fout->Get("MC_measfakcorr");//check
  hrec_sameside_rebin_fakecorr->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin_fakecorr->SetMarkerColor(kBlue-3);     
  hrec_sameside_rebin_fakecorr->SetLineColor(kBlue-3);     
  hrec_sameside_rebin_fakecorr->SetMarkerSize(1.02);       
  
  TH1D* h_recratio_oppfold = (TH1D*)fout->Get("ppData_MeasmFakes_Ratio_DataFold");//check
  h_recratio_oppfold->SetMarkerStyle(kOpenCircle);
  h_recratio_oppfold->SetMarkerColor(kGreen-5);
  h_recratio_oppfold->SetLineColor(  kGreen-5);
  h_recratio_oppfold->SetMarkerSize(1.02);
  
  TH1D* h_recratio_truthfold   = (TH1D*)fout->Get("ppData_MeasmFakes_Ratio_TruthFold");//check
  h_recratio_truthfold->SetMarkerStyle(kOpenStar);
  h_recratio_truthfold->SetMarkerColor(kGreen+3);
  h_recratio_truthfold->SetLineColor(  kGreen+3);
  h_recratio_truthfold->SetMarkerSize(1.02);     
  
  TLegend* legendfold = new TLegend( 0.7,0.7,0.9,0.9);
  legendfold->AddEntry(hfold_truth,"Fold(PY8 Truth)",  "lp");
  legendfold->AddEntry(hrec_sameside_rebin_fakecorr, "PY8 Meas. - Fakes"  ,  "lp");
  legendfold->AddEntry(hfold ,     "Fold(Test PY8 Unf.)" , "lp");
  legendfold->AddEntry(hrec_rebin_fakecorr, "Test PY8 Meas. - Fakes" , "lp");    
  legendfold->SetFillStyle(0);
  legendfold->SetBorderSize(0);    
  
  TCanvas* canv_spectra_foldratio2=NULL;
  setupSpectraHist(hrec_sameside_rebin_fakecorr);
  hrec_sameside_rebin_fakecorr->SetTitle(("Folded Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendfold,
  		  ((std::vector<TH1D*>){hrec_sameside_rebin_fakecorr, hfold_truth,  hfold,  hrec_rebin_fakecorr}),
  		  ((std::vector<TH1D*>){h_recratio_oppfold, h_recratio_truthfold}),
  		  "Ratio w/ (Meas. - Fakes) ");
  drawTLatex( 0.58, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str());       
  canv_spectra_foldratio2=(TCanvas*)canvForPrint->DrawClone();
  canvForPrint->Clear();
  

  

  
  
  // thy spectra, ratios etc. ----------------------------------------------  
  TH1D* CT10nlo	 = (TH1D*)fout->Get("NLO_CT10_NLO_R04_jtpt")	 ;
  CT10nlo->SetMarkerSize(0);
  CT10nlo->SetLineColor(kBlack);  
  TH1D* CT14nlo	 = (TH1D*)fout->Get("NLO_CT14_NLO_R04_jtpt")	 ;
  CT14nlo->SetMarkerSize(0);
  CT14nlo->SetLineColor(kGreen);
  TH1D* NNPDFnnlo= (TH1D*)fout->Get("NLO_NNPDF_NLO_R04_jtpt")    ;
  NNPDFnnlo->SetMarkerSize(0);
  NNPDFnnlo->SetLineColor(kCyan-6);  
  
  TLegend* legendThy1 =new TLegend( 0.6,0.7,1.0,0.9 );    
  if(applyNPCorrs){
    legendThy1->AddEntry(CT10nlo  ,"CT10 NLO #otimes HERWIG EE4C","l");
    legendThy1->AddEntry(CT14nlo  ,"CT14 NLO #otimes HERWIG EE4C","l");
    legendThy1->AddEntry(NNPDFnnlo,"NNPDF NNLO #otimes POW+PY8","l");   }
  else{
    legendThy1->AddEntry(CT10nlo  ,"CT10 NLO","l");
    legendThy1->AddEntry(CT14nlo  ,"CT14 NLO","l");
    legendThy1->AddEntry(NNPDFnnlo,"NNPDF NNLO","l");   }  
  legendThy1->AddEntry(hunf,"Test PY8 Unf.","lp");
  legendThy1->AddEntry(hgen_rebin,"PY8 Truth", "lp");
  legendThy1->SetFillStyle(0);
  legendThy1->SetBorderSize(0);    

  hunf->Scale(2.);
  hgen_rebin->Scale(2.);

  TCanvas* canv_thy_spectra_1=NULL;  
  if(!applyNPCorrs)CT10nlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());    
  else CT10nlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendThy1,
		  ((std::vector<TH1D*>){CT10nlo, CT14nlo, NNPDFnnlo, hgen_rebin}),
		  hunf, "Ratio w/ Test PY8 Unf.");    
  drawTLatex( 0.48, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str()); 
  canv_thy_spectra_1=(TCanvas*)canvForPrint->DrawClone(); 
  canvForPrint->Clear();


  // thy spectra, ratios #2 etc. ----------------------------------------------  
  TH1D* HERAPDF	 = (TH1D*)fout->Get("NLO_HERAPDF105_NLO_R04_jtpt")	 ;
  HERAPDF->SetMarkerSize(0);
  HERAPDF->SetLineColor(kViolet-5);  
  
  TH1D* MMHTnlo	 = (TH1D*)fout->Get("NLO_MMHT2014_NLO_R04_jtpt")	 ;
  MMHTnlo->SetMarkerSize(0);
  MMHTnlo->SetLineColor(kOrange+7);    
  
  TLegend* legendThy2 =new TLegend( 0.6,0.7,1.0,0.9 );    
  if(applyNPCorrs){
    legendThy2->AddEntry(HERAPDF  ,"HEARPDF NLO #otimes HERWIG EE4C","l");
    legendThy2->AddEntry(MMHTnlo  ,"MMHT NLO #otimes HERWIG EE4C","l");}
  else{
    legendThy2->AddEntry(CT10nlo  ,"CT10 NLO","l");
    legendThy2->AddEntry(CT14nlo  ,"CT14 NLO","l");}
  legendThy2->AddEntry(hunf,"Test PY8 Unf.","lp");
  legendThy2->AddEntry(hgen_rebin,"PY8 Truth", "lp");
  legendThy2->SetFillStyle(0);
  legendThy2->SetBorderSize(0);    
  
  TCanvas* canv_thy_spectra_2=NULL;  
  if(!applyNPCorrs)MMHTnlo->SetTitle(("NLO Jet Spectra"+methodString+descString).c_str());    
  else MMHTnlo->SetTitle(("NLO+NPs Jet Spectra"+methodString+descString).c_str());    
  multiratioplot( canvForPrint, legendThy2,
		  ((std::vector<TH1D*>){MMHTnlo, HERAPDF,  hgen_rebin}),
		  hunf, "Ratio w/ Data Unf.");    
  drawTLatex( 0.48, 0.88, 0.035, canvForPrint, desclines);
  canvForPrint->Print(outPdfFile.c_str()); 
  canv_thy_spectra_2=(TCanvas*)canvForPrint->DrawClone(); 
  canvForPrint->Clear();
  

  fout->cd();  
  canv_spectra_genratio->SetTitle("PY8 Truth Ratios Combined Canvas");	      canv_spectra_genratio->Write("combcanv_spectra_genratio");		  
  canv_spectra_recratio->SetTitle("Test PY8 Meas. Ratios Combined Canvas");	      canv_spectra_recratio->Write("combcanv_spectra_recratio");		  
  canv_spectra_unfratio->SetTitle("Test PY8 Unf. Ratios Combined Canvas");	      canv_spectra_unfratio->Write("combcanv_spectra_unfratio");		  
  canv_spectra_foldratio->SetTitle("Folding Test Add Fakes Ratios Combined Canvas");        canv_spectra_foldratio->Write("combcanv_spectra_pfakes_foldratio");    
  canv_spectra_foldratio2->SetTitle("Folding Test Minus Fakes Ratios Combiend Canvas");	    canv_spectra_foldratio2->Write("combcanv_spectra_mfakes_foldratio");	  
  canv_thy_spectra_1->SetTitle("NLO Thy Spectra Combined Canvas");                          canv_thy_spectra_1->Write("combcanv_thy_spectra_1");                 
  canv_thy_spectra_2->SetTitle("NLO Thy Spectra Combined Canvas");                          canv_thy_spectra_2->Write("combcanv_thy_spectra_2");                 
  
  return;
}


//i want the ratio of bias/no bias to look like a pedestal at 1 w/ a mild linear rise, totaling 1%
void addLinearBias(TH1* hrec){
  bool funcDebug=true;
  int nbinsx=hrec->GetNbinsX();
  float pTmin=hrec->GetBinCenter(1);
  float pTmax=hrec->GetBinCenter(nbinsx);
  float maxperc=0.03;
  float slope= -1.*maxperc/(pTmax-pTmin);
  float yint=-1.*slope*pTmax;
  
  TF1* linbias= new TF1 ("linbias", "[0]+x*[1]");
  linbias->SetParameter(0, yint);
  linbias->SetParameter(1, slope);
  
  for(int i=1; i<=nbinsx; i++){
    float x=hrec->GetBinCenter(i);
    if(funcDebug)std::cout<<"bin center for i="<<i<<" is "<<x<<std::endl;
    if(funcDebug)std::cout<<"bin width for i="<<i<<" is "<<hrec->GetBinWidth(i)<<std::endl;
    float val=hrec->GetBinContent(i);
    float valerr=hrec->GetBinError(i);
    if(funcDebug)std::cout<<" old val is "<< val<<std::endl;
    if(funcDebug)std::cout<<" old err is "<< valerr<<std::endl;
    float fact=linbias->Eval(x)+1.;
    if(funcDebug)std::cout<<" fact is "<< fact<<std::endl;    
    val=(val)*(fact) ;
    valerr=(valerr)*(fact) ;    
    if(funcDebug)std::cout<<" new val is "<< val<<std::endl;
    if(funcDebug)std::cout<<" new err is "<< valerr<<std::endl;    
    hrec->SetBinContent(i, val);
    hrec->SetBinError(i, valerr);
    if(funcDebug)std::cout<<"sanity check; bin content is = "<< hrec->GetBinContent(i, val)   <<std::endl;
    if(funcDebug)std::cout<<"sanity check; bin content is = "<< hrec->GetBinContent(i, valerr)<<std::endl;
    if(funcDebug)std::cout<<std::endl;
  }
  
  return;
}


void JERsystStyle(TH1D* JERsystHist=NULL){
  if(!JERsystHist)return;
  JERsystHist->SetFillStyle(1001);//1001-->Solid, 3ijk-->some kind of hatch pattern, 0--> hollow
  JERsystHist->SetFillColorAlpha(kRed, 0.30);
  JERsystHist->SetMarkerColorAlpha(kRed, 0.);//transparent
  JERsystHist->SetLineColorAlpha(kRed,0.);    
  return    ;
}


void makeJERsystSpectra(TH1D* JERsysterr=NULL, TFile* NLOfile=NULL, TH1D* histForError=NULL, bool doJERsyst=false){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"making JER systematics spectra error"<<std::endl;
  if(!NLOfile) {
    std::cout<<"NLOfile not passed. exit."<<std::endl;
    return;}
  if(!JERsysterr){
    std::cout<<"no hist for JER sys error. exit."<<std::endl;
    return;}
  if(! (NLOfile->TFile::IsOpen()) ){
    std::cout<<"TFile not open. exit."<<std::endl;
    return;}
  if(!histForError){
    std::cout<<"no spectra to make errors from. exit."<<std::endl;
  }
  //  TF1* JERsigmaFit=(TF1*)NLOfile->Get("SigmaFit_h");
  std::string fitname="SigmaFit_f";
  if(doJERsyst)fitname+="";//put fit string stuff here
  TF1* JERsigmaFit=(TF1*)NLOfile->Get(fitname.c_str());
  if(!JERsigmaFit){
    std::cout<<"TF1 for JER not found. exit"<<std::endl;
    return;
  }
  
  int nbinsx=JERsysterr->GetNbinsX();    
  JERsystStyle(JERsysterr);
  
  for(int i=1; i<=nbinsx; i++){
    double val=histForError->GetBinContent(i);
    double center=histForError->GetBinCenter(i);
    double JERerr=val*JERsigmaFit->Eval(center);//*1.079;    //JER scaling factor ?!? see makeToySpectra in smearTheory.h
    //double JERerr=val*JERsigmaFit->Eval(center)*5.;    //debug
    JERsysterr->SetBinContent(i,val);
    JERsysterr->SetBinError(i, JERerr);
    if(funcDebug)std::cout<<"i= "<<i<<"  bin content is = "<<val<<std::endl;
    if(funcDebug)std::cout<<"             JERerr = "<<JERerr<<std::endl;
    if(funcDebug)std::cout<<" hist stat error is = "<<histForError->GetBinError(i)<<std::endl;
    
  }
  
  
  return;
}


void makeJERsystRatio(TH1D* JERsysterr=NULL, TFile* NLOfile=NULL){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"making JER systematic error ratio"<<std::endl;
  if(!NLOfile) {
    std::cout<<"NLOfile not passed. exit."<<std::endl;
    return;}
  if(!JERsysterr){
    std::cout<<"no hist for JER sys error. exit."<<std::endl;
    return;}
  if(! (NLOfile->TFile::IsOpen()) ){
    std::cout<<"TFile not open. exit."<<std::endl;
    return;}
  
  std::string fitname="SigmaFit_f";//put jer sysfit up here
  TF1* JERsigmaFit=(TF1*)NLOfile->Get(fitname.c_str());
  if(!JERsigmaFit){
    std::cout<<"TF1 for JER not found. exit"<<std::endl;
    return;
  }
  
  int nbinsx=JERsysterr->GetNbinsX();    
  JERsystStyle(JERsysterr);
  
  for(int i=1; i<=nbinsx; i++){
    double val=0.;//histForError->GetBinContent(i);
    double center=JERsysterr->GetBinCenter(i);
    //double JERerr=val*JERsigmaFit->Eval(center);//*1.079;    //JER scaling factor ?!? see makeToySpectra in smearTheory.h
    double JERerr=JERsigmaFit->Eval(center);//*1.079;    //JER scaling factor ?!? see makeToySpectra in smearTheory.h
    JERerr*=100.;//
    //double JERerr=val*JERsigmaFit->Eval(center)*5.;    //debug
    JERsysterr->SetBinContent(i,val);
    JERsysterr->SetBinError(i, JERerr);
    if(funcDebug)std::cout<<"i= "<<i<<"  bin content is = "<<val<<std::endl;
    if(funcDebug)std::cout<<"             JERerr = "<<JERerr<<std::endl;
    //if(funcDebug)std::cout<<" hist stat error is = "<<histForError->GetBinError(i)<<std::endl;
    
  }
  
  
  return;
}

void StatErrStyle(TH1D* StatErrHist=NULL){
  if(!StatErrHist)return;
  StatErrHist->SetFillStyle(1001);//1001->Solid, 3ijk-->some kind of hatch pattern, 0--> hollow
  StatErrHist->SetFillColorAlpha(  kBlack, 0.);
  StatErrHist->SetMarkerStyle(kFullSquare);//transparent
  
  StatErrHist->SetMarkerColorAlpha(kBlack, 1.);//transparent
  StatErrHist->SetLineColorAlpha(  kBlack,.8);    
  return;
}




void makeStatErrRatio(TH1D* hratio=NULL, TH1D* hratio_staterr=NULL){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"making stat error ratio"<<std::endl;
  if(!hratio) {
    std::cout<<"ratio hist for stat err values not passed. exit."<<std::endl;
    return;}
  if(!hratio_staterr){
    std::cout<<"ratio hist to DEPICT stat err values not passed. exiit."<<std::endl;
    return;}
  
  int nbinsx=hratio_staterr->GetNbinsX();    
  StatErrStyle(hratio_staterr);
  
  for(int i=1; i<=nbinsx; i++){
    double val=0.;
    double err=hratio->GetBinError(i);
    hratio_staterr->SetBinContent(i,val);
    hratio_staterr->SetBinError(i, err*100.);
    if(funcDebug)std::cout<<"i= "<<i<<"  bin content is = "<<val<<std::endl;
    if(funcDebug)std::cout<<"                err = "<<err<<std::endl;
  
    
  }
  
  
  return;
}

float getnonzeromin(TH1* th1){
  float nonzeromin=1.e+40;
  int nbinsx=th1->GetNbinsX();
  for(int i=1;i<=nbinsx;i++){
    float content=th1->GetBinContent(i);
    if(!(content>0.)&&!(content<0.))//check for zero
      continue;
    else if(content < nonzeromin)//check that new content is less than current minium
      nonzeromin=content;
    else continue;//if the content is greater than current minimum
  }
  return nonzeromin;

}

//void divideThyByUnf(TH1* hthyratio, TH1* hunf){
//    
//  bool funcDebug=true;
//  
//  int nthybins=(int)hthyratio->GetNbinsX();
//  int nunfbin=(int)hunf->GetNBinsX();
//  for(int i=1;i<=nthybins; i++){
//    double thybin_lowedge=hthyratio->GetBinLowEdge(i);
//    
//    
//    for(int j=1; j<=nunfbins; j++){
//      double hunfbin_lowedge=hunf->GetBinLowEdge(j);
//      bool edgeMatch=safeDoubleEq(thybin_lowedge,hunfbin_lowedge);
//      
//      if(edgeMatch){
//	double thybinval=hthyratio->GetBinContent(i);
//	double thybinerr=hthyratio->GetBinError(i);
//
//	double hunfbinval=hunf->GetBinContent(i);
//	double hunfbinerr=hunf->GetBinError(i);
//
//	if(!(hunfbinval>0.))
//	  continue;
//	
//	
//	
//      }
//      else 
//	continue;
//      
//    }
//    
//    
//  }
//  
//  
//  return;
//}






//TH1* cloneTH1_oneBinWidth(TH1* hist){
//  int nbinsx=hist->GetNBinsx();
//  
//  TH1* hOneBinWidth = new TH1( (hist->GetName()+"_nubins").c_str(),
//				 (hist->GetTitle()+" one bin width only").c_str(),
//				 (int)nbinsx,0,(double)nbinsx) ;
//  
//  for(int i=1;i<=nbinsx;i++){
//    double content=hist->GetBinContent(i);
//    double err=hist->GetBinError(i);
//    hOneBinWidth->SetBinContent(i,content);
//    hOneBinWidth->SetBinError(i,err);    
//  }
//  
//  
//  return (TH1*)hOneBinWidth;
//}
//
//
//TH1* cloneTH2_oneBinWidth(TH2D* hist){
//  int nbinsx=hist->GetNBinsx();
//  int nbinsy=hist->GetNBinsy();
//  
//  TH2D* hOneBinWidth = new TH2D( (hist->GetName()+"_nubins").c_str(),
//				 (hist->GetTitle()+" one bin width only").c_str(), 
//				 (int)nbinsx,0,(double)nbinsx, 
//				 (int)nbinsy,0,(double)nbinsy 			 ) ;
//  
//  for(int i=1;i<=nbinsx;i++){
//    for(int j=1;j<=nbinsy;j++){
//      double content=hist->GetBinContent(i,j);
//      double err=hist->GetBinError(i,j);
//      hOneBinWidth->SetBinContent(i,j,content);
//      hOneBinWidth->SetBinError(i,j,err);    
//    }
//  }
//  
//  
//  return (TH2D*)hOneBinWidth;
//}


//void checkNRenameFiles (const std::string outFileName, std::string *outRespMatPdfFile, std::string *outPdfFile, std::string *outRootFile){
//  bool funcDebug=false;
//  int outputInd=1;
//  
//  while( (bool)fileExists(*(outRootFile)) ) {
//    if(funcDebug)std::cout<<"fileExists! adding ind="<<outputInd<<std::endl;
//    (*outRespMatPdfFile)=outFileName+"_"+std::to_string(outputInd)+"_respMat.pdf";
//    (*outPdfFile       )=outFileName+"_"+std::to_string(outputInd)+".pdf";
//    (*outRootFile      )=outFileName+"_"+std::to_string(outputInd)+".root";      
//    if(funcDebug)std::cout<<"outPdffile="<<(*outPdfFile)<<std::endl;    
//    outputInd++;
//  }
//  
//  
//  return;
//}



  ////---------------      
  //tempCanvForPdfPrint->cd();
  //tempCanvForPdfPrint->SetLogx(0);
  //tempCanvForPdfPrint->SetLogy(0);
  //tempCanvForPdfPrint->SetLogz(1);  
  ////always use simp bins for covmat to avoid log scaling the x/y axes
  //matStylePrint( (TH2*)covmat_TH2, "Bayesian Unfolding Covariance Matrix", tempCanvForPdfPrint, outRespMatPdfFile, true);
  //
  ////---------------      
  //tempCanvForPdfPrint->cd();
  //tempCanvForPdfPrint->SetLogx(0);
  //tempCanvForPdfPrint->SetLogy(0);
  //tempCanvForPdfPrint->SetLogz(0);
  //
  ////always use simp bins for covmat to avoid log scaling the x/y axes
  //matStylePrint( (TH2*)Pearson, "Bayesian Unfolding Pearson Matrix", tempCanvForPdfPrint, outRespMatPdfFile, true);
  //
  //
  ////---------------      
  //tempCanvForPdfPrint->cd();
  //tempCanvForPdfPrint->SetLogx(0);
  //tempCanvForPdfPrint->SetLogy(0);
  //tempCanvForPdfPrint->SetLogz(1);
  //
  ////always use simp bins for covmat to avoid log scaling the x/y axes
  //matStylePrint( (TH2*)unfmat_TH2, "Bayesian Unfolding Matrix", tempCanvForPdfPrint, outRespMatPdfFile, true);
  
// color/styles
//const Int_t color[13]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,kRed+2,kMagenta+1,kRed+1,kGreen+3};
//// fill/fill-border color/style
//const Int_t fcolor[5]={kRed-10,kBlue-10,kGreen-8,kOrange-4, kRed-6};
//const Int_t lcolor[5]={kRed+1,kBlue-3,kGreen+3,kOrange+3, kRed+2};
//// marker color/style
//const Int_t fmstyle[6] = {20,21,22,23,29,3};
//const Int_t emstyle[6] = {24,25,26,27,30,28};





//void testwrite(TFile* fout){
//  bool funcDebug=true;
//  if(funcDebug)std::cout<<"test write hist"<<std::endl;
//  TH1D* testHistWrite=(TH1D*)fout->Get("Data_meas");
//  fout->cd();
//  testHistWrite->Write("testHistWrite");
//  return;
//}


