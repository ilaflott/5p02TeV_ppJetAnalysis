#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>

const std::string thePDFFileName="printPlots_defOut.pdf";

const std::string radius="4";
const std::string etaWidth ="20_eta_20";

const std::string HLTName[]={"HLT40_","HLT60_","HLT80_","HLT100_","HLTComb_","TrgCombTest_"};
const int N_HLTNames=sizeof(HLTName)/sizeof(std::string);

const std::string var[] = {
  "jtpt" ,  "rawpt",  "jteta", "jtphi",
  "trkMax", "trkSum", "trkHardSum",
  "chMax",  "chSum",  "chHardSum",
  "phMax",  "phSum",  "phHardSum",
  "neMax",  "neSum",
  "eMax",   "eSum",
  "muMax",  "muSum",
  "Aj",     "xj"
};
const int N_vars=sizeof(var)/sizeof(std::string);

const float  var_xLow[N_vars]  = {0, };
const float var_xHigh[N_vars]  = {};
//const float var_yHigh[N_vars]  = {};
//const float  var_yLow[N_vars]  = {};


int main(int argc, char *argv[]){

  const bool debugMode=true;
  int rStatus = -1;

  if( argc==1 ) {//no input arguments, error
    std::cout<<"give me root file output from readFiles!"<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl<<std::endl;
    return rStatus;
  }
  rStatus=1;

  //input file
  std::string filename=argv[1];
  std::cout<< "filename is " << filename <<std::endl;
  TFile *fin = new TFile(filename.c_str());
  
  //get output file ready
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
  std::string open_thePDFFileName=thePDFFileName+"[";
  std::string close_thePDFFileName=thePDFFileName+"]";
  temp_canv->Print( open_thePDFFileName.c_str() );      //open the .pdf file we want to print to
  temp_canv->cd();


  //----------------------
  //first, raw QA plots, 2*N_vars plots
  for(int i=0;i<2;i++){
    for(int j=0;j<N_vars;j++){      
      //e.g. hJetQA_{0,1}wJetID_{jtpt,rawpt,jteta...}
      std::string theHistName="hJetQA_"+std::to_string(i)+"wJetID_"+var[j];
      if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
      TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
      theJetQAHist->Draw();
      temp_canv->Print( thePDFFileName.c_str() );
    }
  }
  //----------------------


  //----------------------
  //then the vz hist, one plot, e.g. hVz
  {
    std::string theHistName="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //----------------------


  //----------------------
  //then trig and trigCombo plots, 12 plots
  for(int i=0;i<2;i++){
    //hpp_{HLT40,HLT60,HLT80,HLT100,HLTComb}_{no,}JetID_R4_20_eta_20
    for(int j=0;j<N_HLTNames;j++){
      std::string theHistName="hpp_"+HLTName[j];
      if(i==0)theHistName+="no";
      theHistName+="JetID_R"+radius+"_"+etaWidth;
      if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
      TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
      theJetQAHist->Draw();
      temp_canv->Print( thePDFFileName.c_str() );
    }
  }
  //----------------------

  fin->Close();
  temp_canv->Print( close_thePDFFileName.c_str() );//close the output file
  return 0;
}

