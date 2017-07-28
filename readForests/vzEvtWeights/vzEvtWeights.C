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
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string SCRATCH_BASE= 
  "/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/readForests/";

const std::string input_dir=SCRATCH_BASE;//+"4.07.16_outputCondor/readForests_HIJetRECO_4-07-17/";
//const std::string readForests_condor_dir=CMSSW_BASE+"src/readForests/outputCondor";

const std::string ppData_inCondorDir=input_dir+"7.25.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_07-25-17_jetPlots_0.0eta0.5/";
const std::string input_ppData_Filename=ppData_inCondorDir+"HighPtJetTrig_ak4PF-allFiles.root";
//const std::string input_ppData_Filename="HighPtJetTrig_ak4PF-allFiles.root";

const std::string ppMC_inCondorDir  =input_dir+"7.25.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_07-25-17_jetPlots_0.0eta0.5/";
const std::string input_ppMC_Filename=ppMC_inCondorDir+"Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
//const std::string input_ppMC_Filename="Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";

int main (int argc, char *argv[]){
  
  
  std::cout<<" now opening Data File "<<std::endl<<input_ppData_Filename<<std::endl<<std::endl;
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

  TH1F* theDataEvtQAHist= (TH1F*)finData->Get( "hWeightedVz" );
  theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(1) );
  theDataEvtQAHist->Scale( 1/effIntgrtdLumi_vz );
  
  //theDataEvtQAHist->SetTitle (    h_Title.c_str() );
  //theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  //theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  //
  //theDataEvtQAHist->SetMarkerStyle(kDot);
  //theDataEvtQAHist->SetMarkerSize(1.1);
  //theDataEvtQAHist->SetMarkerColor( kBlack);
  //theDataEvtQAHist->SetLineColor( theRatioLineColor );
  //theDataEvtQAHist->SetAxisRange(0.,1.5,"Y");
  
  
  std::cout<<" now opening MC File "<<std::endl<<input_ppMC_Filename<<std::endl<<std::endl;
  TFile* finMC = new TFile(input_ppMC_Filename.c_str());  
  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( "hpthatWeightedVz" );
  theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
  theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
  
  TH1F *theRatio=theDataEvtQAHist;
  theRatio->Divide(theMCEvtQAHist);
  //theRatio->Draw();  
  //theRatio->SetLineColor( altRatioLineColor1 );
  //theEvtQALeg->AddEntry(theRatio,"MC not vz-weighted","lp");
  
  Float_t theVzBinWidth=theRatio->TH1::GetBinWidth(1);
  Float_t xLow=-24., xHigh=24.;
  Float_t NvzWeightBins_F=(xHigh-xLow)/(theVzBinWidth);
  Int_t NvzWeightBins=(Int_t)NvzWeightBins_F;//should be 60
  
  std::cout<<"float-division says NBins="<<NvzWeightBins_F<<" (exactly)"<<std::endl;
  std::cout<<"int-typecast says NBins="<<NvzWeightBins<<std::endl;
  
  std::cout<<"now grabbing vzWeights for "<<NvzWeightBins<<" bins for ( "<<xLow<<"< vz <"<<xHigh<<" )"<<std::endl;
  for (int i=0;i<NvzWeightBins;++i){//binsX loop

    Float_t vzWeight = theRatio->TH1::GetBinContent(i+1);    //TH1 bin counting starts at i=1?! why?!
    if(i%5==0)
      
      std::cout<<"i=="<<i<<", vzWeight="<<vzWeight <<" , vzLow="<<xLow<<std::endl;
    else
      std::cout<<"i=="<<i<<", vzWeight="<<vzWeight<< std::endl;

    //Float_t leftSideOfBin=xLow+(i)*theVzBinWidth;
    //Float_t rightSideOfBin=xLow+(i+1)*theVzBinWidth;
    //std::cout<<"i=="<<i<<", "<<leftSideOfBin<<"<vz<="<<rightSideOfBin<<", vzWeight="<<vzWeight<< std::endl;
    xLow=xLow+theRatio->TH1::GetBinWidth(1);
    
    
  }
  
  return 0 ;
}

































































