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
const int kIter = 4; // default is 4, Bayes

// max diff # of kreg to do  /  "width" of kreg from center; i.e. kregs looped over will be kRegCenter +/- kRegRange
const int nKregMax  = 9 , kRegRange=(nKregMax-1)/2 ; //SVD

const double NLOMCscaling=1e+00; // i.e. for when i fuck up the scaling/normalization
const double MCscaling=1e+00; // i.e. for when i fuck up the scaling/normalization
const double etaBinWidth=4.;   // e.g. |y| < 2.0
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
const bool useTH2ProjRespHist= fillRespHists&&false;

// settings that don't really get used
const bool doOverUnderflows   =false;//leave false almost always
const bool clearOverUnderflows=true; 
const bool zeroBins=false; //leave false almost always
//const bool normalizedMCMatrix=false;   

// -----------------------------------------------------------------------------------------------------------------------
//useful strings, numbers
const std::string CMSPRELIM= "CMS PRELIMINARY"; 
const std::string MCdesc= "Py8 Tune CUETP8M1 QCD"; 
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
  else if (hTitle.find("Covariance")  != std::string::npos)
    mat->SetAxisRange(1.e-15,1.e+02,"Z");
  else if (hTitle.find("Pearson")  != std::string::npos)
    mat->SetAxisRange(-1.,1.,"Z");
  else if(hTitle.find("Unfolding")  != std::string::npos)
    mat->SetAxisRange(1.e-08,1.e+03,"Z");
  else
    mat->SetAxisRange(1.e-10,1.e+03,"Z");
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
  h->SetAxisRange(0.5, 1.6, "Y");
  h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
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

void setupSpectraHist(TH1* h, bool useSimpBins, double* boundaries=NULL, int nbins=1){  
  h->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{d#eta dp_{T}} [nb/GeV]");
  //h->GetYaxis()->SetTitle("N_{Jets}/L_{int}");
  //h->GetYaxis()->SetTitle("A.U.");
  h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
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
  bool isForNLOMC=(((std::string)hmat->GetName()).find("ynew")!=std::string::npos);
  double xLow,xHi;//=boundaries_pt_reco_mat[0];
  double yLow, yHi;//=boundaries_pt_gen_mat[0];
  if( (bool)fpp_MC && !isForNLOMC){
    xLow=(double)std::atof( ( (TH1*) (fpp_MC->Get("hJetPtCut_unf_lo") )
			      )->GetTitle() );
    xHi=(double)std::atof( ( (TH1*) (fpp_MC->Get("hJetPtCut_unf_hi") )
			     )->GetTitle() );
    yLow=(double)std::atof( ( (TH1*) (fpp_MC->Get("hGenJetPtCut_unf_lo") )
			      )->GetTitle() );
    yHi=(double)std::atof( ( (TH1*) (fpp_MC->Get("hGenJetPtCut_unf_hi") )
			     )->GetTitle() );	
    hmat->SetAxisRange(xLow,xHi,"X");
    hmat->SetAxisRange(yLow,yHi,"Y");        
  }
//  else{
//    xLow=boundaries_pt_reco_mat[0];
//    xHi=boundaries_pt_reco_mat[nbins_pt_reco_mat];
//    yLow=boundaries_pt_gen_mat[0];
//    yHi=boundaries_pt_gen_mat[nbins_pt_gen_mat]; }
  
  
  
  
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
  if(funcDebug)std::cout<<"matStylePrint(hmat_rebin)"<<std::endl;
  matStylePrint(hmat_rebin, "MC Response Matrix, Rebinned", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  
  // error matrix in binning of interest ---------------    
  if(funcDebug)std::cout<<"matStylePrint(hmat_errors)"<<std::endl;
  matStylePrint(hmat_errors, "MC Response Matrix, Bin Errors", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  
  // percent error matrix in binning of interest ---------------    
  if(funcDebug)std::cout<<"makeRespMatrixPercentErrs(hmat_errors, hmat_rebin)"<<std::endl;
  TH2D* hmat_percenterrs= makeRespMatrixPercentErrs( (TH2D*) hmat_errors, (TH2D*) hmat_rebin,
						     (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						     (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
  if(funcDebug) hmat_percenterrs->Print("base");
  if(funcDebug)std::cout<<"matStylePrint(hmat_percenterrs)"<<std::endl;
  matStylePrint(hmat_percenterrs, "MC Response Matrix, Bin % Errors", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  
  
  // col normd matrix in binning of interest  ---------------    
  // POTENTIAL ISSUE: using the resp matrix post rebinning/divide bin width/clearing overflows... should i be using the original matrix?
  if(funcDebug)std::cout<<"normalizeCol_RespMatrix(hmat_rebin)"<<std::endl;
  TH2D* hmat_rebin_colnormd = normalizeCol_RespMatrix( (TH2D*)  hmat_rebin,
						       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  if(funcDebug)  hmat_rebin_colnormd->Print("base");  
  
  if(funcDebug)std::cout<<"matStylePrint(hmat_rebin_colnormd)"<<std::endl;
  matStylePrint(hmat_rebin_colnormd, "MC Response Matrix, Columns Sum to 1", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  
  // row normd matrix in binning of interest  ---------------    
  if(funcDebug)std::cout<<"normalizeRow_RespMatrix(hmat_rebin)"<<std::endl;
  TH2D*  hmat_rebin_rownormd = normalizeRow_RespMatrix( (TH2D*)  hmat_rebin,
							(double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
							(double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  if(funcDebug)  hmat_rebin_rownormd->Print("base");
  if(funcDebug)std::cout<<"matStylePrint(hmat_rebin_rownormd)"<<std::endl;
  matStylePrint(hmat_rebin_rownormd, "MC Response Matrix, Rows Sum to 1", tempCanvForPdfPrint, outRespMatPdfFile, useSimpBins);
  
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



TH1* makeThyHist_00eta20_v2(std::string filename, bool applyNPCorrFactor=true){
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



