///////// (initializa/declara)tions //////////
// C++, C, etc.
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
// ROOTSYS
#include <TSystem.h>
#include <TProfile.h>
#include <TROOT.h>
// I/O
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TEventList.h>
#include <TCut.h>
// utilities
#include <TMath.h>
#include <TRandom3.h>
#include <TStopwatch.h>
// plotting
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVirtualFitter.h>
#include <Minuit2/Minuit2Minimizer.h>
// histos
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>


// I/O ------------------------
// env ------------------------

//CMSSW_BASE
const std::string CMSSW_BASE= 
  "/home/obaron/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string SCRATCH_BASE= 
  "/cms/heavyion/ilaflott/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/";

const std::string input_dir=SCRATCH_BASE;//+"4.07.16_outputCondor/readForests_HIJetRECO_4-07-17/";
//const std::string readForests_condor_dir=CMSSW_BASE+"src/readForests/outputCondor";

//const std::string ppData_inCondorDir=input_dir+"10.11.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_10-11-17_jetPlots_0.0eta2.0_rec56_HLTCalo/";
//const std::string ppData_inCondorDir="02.04.18_outputCondor/ppData_HighPtJetTrig_ak4PFJets_02-04-18_jetPlots_0.0eta4.7/";
const std::string ppData_inCondorDir="04.13.18_outputCondor/ppData_HighPtJetTrig_ak4PFJets_04-13-18_jetPlots_0.0eta0.5/";
const std::string input_ppData_Filename=input_dir+ppData_inCondorDir+"HighPtJetTrig_ak4PF-allFiles.root";
//const std::string input_ppData_Filename="HighPtJetTrig_ak4PF-allFiles.root";

const std::string ppMC_inCondorDir  ="04.13.18_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_04-13-18_jetPlots_0.0eta0.5/";
const std::string input_ppMC_Filename=input_dir+ppMC_inCondorDir+"Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
//const std::string input_ppMC_Filename="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";






const std::string ppData_inHist=
  //  "hVz"; //weights=1 always, no trigger requirement as of 4/11/18 sample
  //"hVz_BSshift"; //weights=1 always, no trigger requirement as of 4/11/18 sample
  //  "hVz_baryshift_pos"; //weights=1 always, no trigger requirement as of 4/11/18 sample
"hVz_baryshift_neg"; //weights=1 always, no trigger requirement as of 4/11/18 sample
  //"hWeightedVz"; //trigger weighted; events w/ weight=0 included as of 4/11/18 sample
  //"hTriggerVz_noWeights"; //triggered events (weight>0) but weight=1 always
  //"hTriggerVz_withWeights"; //triggered events (weight>0) with weights

const std::string outfilename="Py8_vzEvtWeights_"+ppData_inHist+"_04_11_18.root";

const std::string ppMC_inHist=
  //"hVz"; //vz w/ weights=1 always
  //"hvzWeightedVz";      //vz weights applied only
  "hpthatWeightedVz";   //pthat weights applied only
  //"hWeightedVz";        // vz and pthat weights applied. equiv to pthatweights only if vz weights not actually used

const std::string polyFit="pol8";
const std::string fitOptions="FRVM";
// F--> use minuit fitter for polN
// R--> use range specified in function range
// V--> verbose
// W--> Set weights to one for non-empty bins, ignore error bars
// WW--> Set weights to one for all bins (including empty), ignore error bars
// L--> use log likelihood method
// LL--> use improved log likelihood fit in case of low stat and when bin contents aren't integers. DO NOT USE FOR BIN CONTENT GREATER THAN 100
// E--> perform better errors with Minos technique
// M--> improve fit results (whatever this means)






const bool debugMode=true;

const bool doRootFileOut=false;

int main (int argc, char *argv[]){

  TFile* fout=NULL;

  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat();
  gStyle->SetStatW(0.095);//this sets the fit box W/H too, not just the stat box
  gStyle->SetStatH(0.095);
  
  gROOT->ForceStyle();
  
  if(debugMode) std::cout<<"   input_dir = "<<input_dir<<std::endl; 
  std::cout<<" data condor dir: "<<ppData_inCondorDir<<std::endl;
  std::cout<<" MC condor dir: "<<ppMC_inCondorDir<<std::endl;
  std::cout<<std::endl;
  
  
  std::cout<<"now opening data file"<<std::endl;
  if(debugMode)  std::cout<<"input_ppData_Filename= "<<input_ppData_Filename<<std::endl;
  TFile* finData = new TFile(input_ppData_Filename.c_str());
  
  TH1F *h_NEvents         = (TH1F*)finData->Get("NEvents_read");
  TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_vzCut");
  
  const Int_t NEvts=h_NEvents->GetEntries();
  const Int_t NEvtsPassvzCut=h_NEvents_vzCut->GetEntries();
  std::cout<<std::endl;
  std::cout<<"NEvents in dataset="<<NEvts<<std::endl;
  std::cout<<"NEvents passing vz/skim cut="<<NEvtsPassvzCut<<std::endl;
  std::cout<<std::endl;
  
  const double LumiEff_vz = (double) NEvtsPassvzCut/NEvts;//h_NEvents_vzCut->GetEntries()/h_NEvents->GetEntries(); 
  std::cout<<std::endl;
  std::cout<<"Lumi Efficiency post vz+skim cuts=NEvents_vzCut/NEvents ";
  std::cout<<"=LumiEff_vz ="<<LumiEff_vz<<std::endl;
  std::cout<<std::endl;

  const double intgrtdLumi=27.8*pow(10.,6.);// 25.8 pb^-1 to \microbarns^-1
  const double effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz; 
  std::cout<<std::endl;
  std::cout<<"effective Intgrted Lumi ="<<effIntgrtdLumi_vz<<std::endl;
  std::cout<<std::endl;




  std::cout<<"now opening MC file"<<std::endl;
  if(debugMode)  std::cout<<"input_ppMC_Filename= "<<input_ppMC_Filename<<std::endl;
  TFile* finMC = new TFile(input_ppMC_Filename.c_str());  
  
  std::cout<<" now opening input histograms from input files "<<std::endl;
  
  
  TH1F* theDataEvtInputHist= (TH1F*)finData->Get( ppData_inHist.c_str() );  
  std::cout<<"NBins Orig Data Hist     = "<< theDataEvtInputHist->GetNbinsX() << std::endl;
  std::cout<<"Bin Width Orig Data Hist = "<< theDataEvtInputHist->GetBinWidth(1) << std::endl;
  std::cout<<"xLow Orig Data Hist      = "<< theDataEvtInputHist->GetBinLowEdge(1) << std::endl;
  std::cout<<"xHigh Orig Data Hist     = "<< theDataEvtInputHist->GetBinLowEdge(theDataEvtInputHist->GetNbinsX()) + theDataEvtInputHist->GetBinWidth(1) << std::endl;
  std::cout<<std::endl;
  
  TH1F* theMCEvtInputHist= (TH1F*)finMC->Get( ppMC_inHist.c_str());
  std::cout<<"NBins Orig MC Hist     = "<< theMCEvtInputHist->GetNbinsX() << std::endl;
  std::cout<<"Bin Width Orig MC Hist = "<< theMCEvtInputHist->GetBinWidth(1) << std::endl;
  std::cout<<"xLow Orig MC Hist      = "<< theMCEvtInputHist->GetBinLowEdge(1) << std::endl;
  std::cout<<"xHigh Orig MC Hist     = "<< theMCEvtInputHist->GetBinLowEdge(theMCEvtInputHist->GetNbinsX()) + theMCEvtInputHist->GetBinWidth(1) << std::endl;
  std::cout<<std::endl;
  
  const int dataRebinInt=10*2;//works hVz/hWeightedVz (ppData sample from feb 4 2018)  
  //const int dataRebinInt=2;//works hTriggerVz_noWeights (ppData sample from feb 4 2018)  
  
  TH1F* theDataEvtQAHist = (TH1F*)theDataEvtInputHist->Rebin(dataRebinInt,"theDataEvtQAHist");  
  theDataEvtQAHist->Scale(1./theDataEvtQAHist->Integral());
  //  theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(1) );
  //  theDataEvtQAHist->Scale( 1/effIntgrtdLumi_vz );
  std::cout<<"NBins rebinned Data Hist     = "<< theDataEvtQAHist->GetNbinsX() << std::endl;
  std::cout<<"Bin Width rebinned Data Hist = "<< theDataEvtQAHist->GetBinWidth(1) << std::endl;
  std::cout<<"xLow rebinned Data Hist      = "<< theDataEvtQAHist->GetBinLowEdge(1) << std::endl;
  std::cout<<"xHigh rebinned Data Hist     = "<< theDataEvtQAHist->GetBinLowEdge(theDataEvtQAHist->GetNbinsX()) + theDataEvtQAHist->GetBinWidth(1) << std::endl;
  std::cout<<std::endl;
  
  
  
  TH1F* theRatio = (TH1F*)theDataEvtQAHist->Clone("theRatioHist");
  //  theRatio->Scale( 1/theRatio->GetBinWidth(1) );
  //  theRatio->Scale( 1/effIntgrtdLumi_vz );
  //  std::cout<<"NBins rebinned Data Hist     = "<< theRatio->GetNbinsX() << std::endl;
  //  std::cout<<"Bin Width rebinned Data Hist = "<< theRatio->GetBinWidth(1) << std::endl;
  //  std::cout<<"xLow rebinned Data Hist      = "<< theRatio->GetBinLowEdge(1) << std::endl;
  //  std::cout<<"xHigh rebinned Data Hist     = "<< theRatio->GetBinLowEdge(theRatio->GetNbinsX()) + theRatio->GetBinWidth(1) << std::endl;
  //  std::cout<<std::endl;
  
  //int MCRebinInt=1*2;//works hpthatWeightedVz (ppMC sample from feb 4 2018)
  const int MCRebinInt=10*2;//works hVz (ppMC sample from feb 4 2018)  
  TH1F* theMCEvtQAHist= (TH1F*)theMCEvtInputHist->Rebin(MCRebinInt,"theMCEvtQAHist");  
  //theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(1) );
  theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
  std::cout<<"NBins rebinned MC Hist     = "<< theMCEvtQAHist->GetNbinsX() << std::endl;
  std::cout<<"Bin Width rebinned MC Hist = "<< theMCEvtQAHist->GetBinWidth(1) << std::endl;
  std::cout<<"xLow rebinned MC Hist      = "<< theMCEvtQAHist->GetBinLowEdge(1) << std::endl;
  std::cout<<"xHigh rebinned MC Hist     = "<< theMCEvtQAHist->GetBinLowEdge(theMCEvtQAHist->GetNbinsX()) + theMCEvtQAHist->GetBinWidth(1) << std::endl;
  std::cout<<std::endl;
  
  assert( (theDataEvtQAHist->GetNbinsX()) == (theMCEvtQAHist->GetNbinsX()) );
  
  
  float xLow=theMCEvtQAHist->GetBinLowEdge(1);
  float xHigh = theMCEvtQAHist->GetBinLowEdge(theMCEvtQAHist->GetNbinsX()) + theMCEvtQAHist->GetBinWidth(1);
  float fit_xLow=-24.,fit_xHigh=24.;  
  
  int nbins=theMCEvtQAHist->GetNbinsX();
  

  //theDataEvtQAHist->SetTitle (    h_Title.c_str() );
  //theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  //theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  //
  //theDataEvtQAHist->SetMarkerStyle(kDot);
  //theDataEvtQAHist->SetMarkerSize(1.1);
  //theDataEvtQAHist->SetMarkerColor( kBlack);
  //theDataEvtQAHist->SetLineColor( theRatioLineColor );
  //theDataEvtQAHist->SetAxisRange(0.,1.5,"Y");



  
  TCanvas* TCan = NULL;//general canvas
  TCanvas* TCanMC = NULL;
  TCanvas* TCanData = NULL;
  
  TCanvas* TCanWeightRat = NULL;
  TCanvas* TCanWeightfn = NULL;
  TCanvas* TCanWeightbin = NULL;
  TCanvas* TCanWeightpoly = NULL;
  
  TCanvas* TCanDatMCRat = NULL;
  TCanvas* TCanBinPoly = NULL;
  TCanvas* TCanGausPoly = NULL;
  
  //std::cout<<"Input Hist Stats:"<<std::endl;
  //std::cout<<"Data hist nbins = "<<theDataEvtQAHist->GetNbinsX()<<"  Data hist bin width = "<<theDataEvtQAHist->GetBinWidth(1)<<std::endl;
  //std::cout<<"MC hist nbins = "<<theMCEvtInputHist->GetNbinsX()<<"  MC hist bin width = "<<theMCEvtInputHist->GetBinWidth(1)<<std::endl;

  

  TCan = new TCanvas("TCan","c",900,900);   
  TCanMC = new TCanvas("TCanMC","cMC",900,900);   
  TCanData = new TCanvas("TCanData","cData",900,900);   
  TCanWeightRat = new TCanvas("TCanWeightRat","Weight: Bin/Fn",900,900);   
  TCanWeightfn = new TCanvas("TCanWeightfn","Weight from Function",900,900);
  TCanWeightbin = new TCanvas("TCanWeightbin","Weight from Bin",900,900);
  TCanDatMCRat = new TCanvas("TCanDatMCRat","Data / MC ratio",900,900);
  TCanWeightpoly = new TCanvas("TCanWeightpoly","Weight from Polynomial",900,900);
  TCanBinPoly = new TCanvas("TCanBinPoly","Bin / Poly ratio",900,900);
  TCanGausPoly = new TCanvas("TCanGausPoly","Gaus / Poly ratio",900,900);
  
  
  //TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("theDataHistClone"); I'm going to replace this with the rebin which should make its own clone
  //double norm = theRatio->GetMaximumStored();

  //  TH1F* binWeight = new TH1F("binWeight","Bin-based Weight",100,-24,24);
  //  TH1F* fnWeight = new TH1F("fnWeight","Fn-based Weight",100,-24,24);
  //  TH1F* tPolyweight = new TH1F("tPolyweight","Polynomial-based Weight",100,-24,24);  
  
  //  TH1F* binWeight = new TH1F("binWeight","Bin-based Weight" ,           96, -24, 24);
  //  TH1F* fnWeight = new TH1F("fnWeight","Fn-based Weight"    ,           96, -24, 24);
  //  TH1F* tPolyweight = new TH1F("tPolyweight","Polynomial-based Weight", 96, -24, 24); 
  
  TH1F* binWeight = new TH1F("binWeight","Bin-based Weight" ,           nbins, xLow, xHigh);
  TH1F* fnWeight = new TH1F("fnWeight","Fn-based Weight"    ,           nbins, xLow, xHigh);
  TH1F* tPolyweight = new TH1F("tPolyweight","Polynomial-based Weight", nbins, xLow, xHigh); 
 
  //make fit functions
  TF1* fgaussData = NULL;
  TF1* fgaussMC = NULL;
  fgaussData = new TF1("fgaussData","gaus", fit_xLow, fit_xHigh);
  fgaussMC = new TF1("fgaussMC","gaus", fit_xLow, fit_xHigh);

  TF1* fpolyRat = NULL;
  fpolyRat = new TF1("fpolyRat",polyFit.c_str(), fit_xLow, fit_xHigh);
  
  int fitstatusdata = -99; 
  std::cout<<"fitting data w/ gaussian"<<std::endl;
  fitstatusdata = theDataEvtQAHist->Fit(fgaussData, fitOptions.c_str());  
  std::cout<< "Data Fit Status: "<< fitstatusdata<< ", Fit Error: "<< fgaussData->GetParError(1)<< std::endl;
  
  int fitstatusMC = -99; 
  std::cout<<"fitting MC w/ gaussian"<<std::endl;
  fitstatusMC = theMCEvtQAHist->Fit(fgaussMC, fitOptions.c_str());
  std::cout<< "MC Fit Status: "<< fitstatusMC<< ", Fit Error: "<< fgaussMC->GetParError(1)<< std::endl;
 
 
  //this gives us the bin-by-bin weights
  theRatio->Divide(theMCEvtQAHist);
  
  TH1F* tRatPoly = (TH1F*)theRatio->Clone("tRatPoly");
  std::cout<<"fitting ratio w/ polynomial"<<std::endl;
  int fitstatusPoly = tRatPoly->Fit(fpolyRat, fitOptions.c_str());
  std::cout<< "Polynomial Fit Status: "<< fitstatusPoly<< ", Fit Error: "<< fpolyRat->GetParError(1)<< std::endl;
  
  std::cout<<"now grabbing vzWeights for "<<nbins<<" bins for ( "<<xLow<<"< vz <"<<xHigh<<" )"<<std::endl;
  
  
  
  for (int i=1;i<=nbins;++i){//binsX loop
    
    
    //    //when binWeight hist had 96 bins from -24 to +24 
    //    if(hist_xLow > -24.5 && hist_xLow < 24.)
    //      {	
    //	Float_t binWeight_xLow = binWeight->TH1::GetBinLowEdge(i-2);
    //	std::cout<<"binWeight Low Bin Edge = "<<binWeight_xLow<<std::endl;
    //	binWeight->SetBinContent(i-2,vzWeight);      	
    //	std::cout<<"binWeight content = "<<binWeight->GetBinContent(i-2)<<std::endl;
    //      }
    //Double_t test = 0.4;

    Float_t hist_xLow = theRatio->TH1::GetBinLowEdge(i);
    Float_t binWeight_xLow = binWeight->TH1::GetBinLowEdge(i);

    std::cout<<std::endl;
    std::cout<<"i = "<< i << std::endl;
    std::cout<<"MC input hist LowEdge = "<<theMCEvtQAHist->GetBinLowEdge(i)<<std::endl;
    std::cout<<"theRatio Low Bin Edge = "<<hist_xLow<<std::endl;
    std::cout<<"binWeight Low Bin Edge = "<<binWeight_xLow<<std::endl;            
    std::cout<<std::endl;

    float binCenter = theMCEvtQAHist->GetBinLowEdge(i)+(theMCEvtQAHist->GetBinWidth(i)/2.);
    std::cout<<"binCenter="<< binCenter << std::endl;
    
    Float_t vzWeight = theRatio->TH1::GetBinContent(i);    //TH1 bin counting starts at i=1?! why?!    
    
    if( (binWeight_xLow >= fit_xLow) &&
	(binWeight_xLow < fit_xHigh) ){
      
      //function fit ratio weights 
      //No no no - this needs to be the x value, not the bin content! //so do I want center, low edge, or high edge? Or something else?
      //Double_t gaussMC =     fgaussMC->Eval(theMCEvtQAHist->GetBinLowEdge(i)+(theMCEvtQAHist->GetBinWidth(i)/2.));
      //Double_t gaussData = fgaussData->Eval(theMCEvtQAHist->GetBinLowEdge(i)+(theMCEvtQAHist->GetBinWidth(i)/2.));
      //Double_t polyFit =     fpolyRat->Eval(theMCEvtQAHist->GetBinLowEdge(i)+(theMCEvtQAHist->GetBinWidth(i)/2.));
      
      Double_t gaussMC =     fgaussMC->Eval(binCenter);
      Double_t gaussData = fgaussData->Eval(binCenter);
      Double_t polyFit =     fpolyRat->Eval(binCenter);
      
      Double_t gaussFit = (gaussData/gaussMC);     
      
      binWeight->SetBinContent(i,vzWeight);      	
      fnWeight->SetBinContent(i,gaussFit);
      tPolyweight->SetBinContent(i,polyFit);
      
      std::cout<<"theRatio bin content                = "<<theRatio->GetBinContent(i)<<std::endl;
      std::cout<<"binWeights weight                   = "<<binWeight->GetBinContent(i)<<std::endl;
      std::cout<<"fnWeight (gaussFits) weight         = "<<fnWeight->GetBinContent(i)<<std::endl;
      std::cout<<"tPolyWeight (polynomial fit) weight = "<<tPolyweight->GetBinContent(i)<<std::endl;      
      
      std::cout<<std::endl;
      
    }
    else 
      std::cout<<"outside of vz fit range ("<<fit_xLow<<" to "<< fit_xHigh << ")" << std::endl;

  }

  if(doRootFileOut)
    fout=new TFile(outfilename.c_str(),"RECREATE");
  
  
  
  //orig histos
  TCan->cd();
  TCan->SetLogy(1);
  theDataEvtInputHist = (TH1F*)theDataEvtInputHist->Rebin(dataRebinInt,"theDataEvtInputHist4Draw");  
  theDataEvtInputHist->SetLineColor(kBlack);
  theDataEvtInputHist->Scale( 1./effIntgrtdLumi_vz );
  theDataEvtInputHist->Scale( 1./theDataEvtInputHist->GetBinWidth(1) );  
  theDataEvtInputHist->Draw("E");  
  
  theMCEvtInputHist= (TH1F*)theMCEvtInputHist->Rebin(MCRebinInt,"theMCEvtInputHist4Draw");
  theMCEvtInputHist->SetLineColor(kRed);  
  theMCEvtInputHist->Scale( 1./theMCEvtInputHist->GetBinWidth(1) );  
  theMCEvtInputHist->Scale( theDataEvtInputHist->Integral()/theMCEvtInputHist->Integral() );
  theMCEvtInputHist->Draw("HIST E SAME");

  if(fout){
    fout->cd();
    TCan->SetTitle("MC v. Data evt vz Canvas");
    TCan->Write("MCDataVz_canv");
  }
  else
    TCan->Print("MCDataVz.png","png");
  
  //Compare fit to histogram
  TCanMC->cd();
  TCanMC->SetLogy(1);
  theMCEvtQAHist->SetTitle("Gaussian Fit of MC");
  theMCEvtQAHist->Draw();
  fgaussMC->Draw("same");
  if(fout){
    fout->cd();
    TCanMC->SetTitle("MC evt vz Gauss Fit Canvas");
    TCanMC->Write("gaussfitMC_canv");
  }
  else  
    TCanMC->Print("gaussfitMC.png","png");
  
  TCanData->cd();
  TCanData->SetLogy(1);
  theDataEvtQAHist->SetTitle("Gaussian Fit of Data");
  theDataEvtQAHist->Draw();
  fgaussData->Draw("same");
  if(fout){
    fout->cd();
    TCanData->SetTitle("Data evt vz Gauss Fit Canvas");
    TCanData->Write("gaussfitData_canv");
  }
  else  
    TCanData->Print("gaussfitData.png","png");
  
  
  
  
  TCanWeightfn->cd();
  fnWeight->Draw();
  if(fout){
    fout->cd();
    TCanWeightfn->SetTitle("Data/MC Gauss Fit Weights Canvas");
    TCanWeightfn->Write("weightFn_canv");
  }
  else  
    TCanWeightfn->Print("weightFn.png","png");
  
  TCanWeightbin->cd();
  binWeight->Draw();
  if(fout){
    fout->cd();
    TCanWeightbin->SetTitle("Data/MC Bin-by-Bin Weights Canvas");
    TCanWeightbin->Write("weightBin_canv");
  }
  else  
    TCanWeightbin->Print("weightBin.png","png"); 
  
  TCanWeightpoly->cd();
  tPolyweight->Draw();
  if(fout){
    fout->cd();
    TCanWeightpoly->SetTitle("Data/MC Polynomial Fit Weights Canvas");
    TCanWeightpoly->Write("weightPoly_canv");
  }
  else  
    TCanWeightpoly->Print("weightPoly.png","png");
  
  
  TH1F* BinOverGauss = (TH1F*)binWeight->Clone("tweightRatio");
  BinOverGauss->Divide(fnWeight);
  BinOverGauss->SetTitle("Bin Weights Over Gaussian Weights");
  BinOverGauss->SetAxisRange(.4,2.8,"Y");
  TCanWeightRat->cd();
  BinOverGauss->Draw();
  if(fout){
    fout->cd();
    TCanWeightRat->SetTitle("Bin-by-Bin v. Gaussian Fit Weights Ratio Canvas");
    TCanWeightRat->Write("BbB_v_GaussFit_wgtsRatio_canv");
  }
  else  
    TCanWeightRat->Print("BbB_v_GaussFit_wgtsRatio.png","png");
  
  TH1F* BinOverPoly = (TH1F*)binWeight->Clone("tBinPolyRatio");
  BinOverPoly->Divide(tPolyweight);
  BinOverPoly->SetTitle("Bin Weights Over Polynomial Weights");
  BinOverPoly->SetAxisRange(.4,2.8,"Y");
  TCanBinPoly->cd();
  BinOverPoly->Draw();
  if(fout){
    fout->cd();
    TCanBinPoly->SetTitle("Bin-by-Bin v. Polynomial Fit Weights Ratio Canvas");
    TCanBinPoly->Write("BbB_v_Polynomial_wgtsRatio_canv");
  }
  else      
    TCanBinPoly->Print("BbB_v_Polynomial_wgtsRatio.png","png");
  
  
  TH1F* GaussOverPoly = (TH1F*)fnWeight->Clone("tGausPolyRatio");
  GaussOverPoly->Divide(tPolyweight);
  GaussOverPoly->SetTitle("Gaussian Weights Over Polynomial Weights");
  GaussOverPoly->SetAxisRange(.4,2.8,"Y");
  TCanGausPoly->cd();
  GaussOverPoly->Draw();
  if(fout){
    fout->cd();
    TCanGausPoly->SetTitle("Gaussian Fit v. Polynomial Fit Weights Ratio Canvas");
    TCanGausPoly->Write("Gauss_v_Polynomial_wgtsRatio_canv");
  }
  else      
    TCanGausPoly->Print("Gauss_v_Polynomial_wgtsRatio.png","png");
  
  
  TCanDatMCRat->cd();
  tRatPoly->SetTitle("Polynomial Fit of Data/MC");
  tRatPoly->Draw();
  fpolyRat->Draw("same");
  if(fout){
    fout->cd();
    TCanDatMCRat->SetTitle("Poly Fit Data/MC evt vz Canvas");
    TCanDatMCRat->Write("PolyFit_DataMC_evtvz_wgts_canv");
  }
  else      
    TCanDatMCRat->Print("PolyFit_DataMC_evtvz_wgts_canv.png","png");
  
  if(fout)
    fout->Close();
  std::cout<<"program end"<<std::endl; 
  
  return 0 ;
}
