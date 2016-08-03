#include "readFiles_ppData.h"

const std::string CMSSW_BASE =
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string inputDir =
  "src/readFiles/ppData/saved_outputCondor/";
const std::string inputFilename =
  "readFiles_ppData_5p02TeV_HighPtJetTrig_2016-06-10_allFiles/HighPtJetTrig_ak4PF-allFiles.root";
const std::string fullFilename = CMSSW_BASE+inputDir+inputFilename;

const std::string thePDFFileName="readFiles_ppData_printPlots_HighPtJetTrigs_JetQA-vz-TrigComb_6.11.16.pdf";

const std::string radius="4";
const std::string HLTName[]={"HLT40_","HLT60_","HLT80_","HLT100_","HLTComb_","TrgCombTest_"};
const int N_HLTNames=sizeof(HLTName)/sizeof(std::string);

//const float  var_xLow[N_vars]  = {};
//const float var_xHigh[N_vars]  = {};
//const float var_yHigh[N_vars]  = {};
//const float  var_yLow[N_vars]  = {};

int main(int argc, char *argv[]){

  const bool debugMode=true;
  int rStatus = -1;
  if( argc!=1 ) {//no input arguments, error
    std::cout<<"settings hard coded, just do ./readFiles_ppMC_printPlots.exe"<<std::endl;
    return rStatus;
  }
  rStatus=1;

  //input file
  std::string filename=fullFilename;
  std::cout<< "filename is " << filename <<std::endl;
  TFile *fin = new TFile(filename.c_str());
  
  //get output file ready
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
  std::string open_thePDFFileName=thePDFFileName+"[";
  std::string close_thePDFFileName=thePDFFileName+"]";
  temp_canv->Print( open_thePDFFileName.c_str() );      //open the .pdf file we want to print to
  temp_canv->cd();


  // evt plots
  //vz
  {
    std::string theHistName="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //----------------------


  // jetQA plots ----------------------
  //{0,1}wJetID
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


  // HLT/L1 and Comb plots ----------------------
  //hpp_{HLT40,HLT60,HLT80,HLT100,HLTComb}_{no,}JetID_R4_20_eta_20
  for(int i=0;i<2;i++){
    for(int j=0;j<N_HLTNames;j++){
      std::string theHistName="hpp_"+HLTName[j];
      if(i==0)theHistName+="no";
      theHistName+="JetID_R"+radius+"_"+(std::string)etaWidth;
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

