#include "readFiles_ppMC.h"

const std::string CMSSW_BASE = 
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string inputDir = 
  "src/readFiles/ppMC/saved_outputCondor/";
const std::string inputFilename = 
  "readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_2016-07-09/QCDjetAllPtBins_ak4PF-allFiles.root";
const std::string fullFilename = CMSSW_BASE+inputDir+inputFilename;

const std::string thePDFFileName = "readFiles_printPlots_ppMC_PY6_TuneZ2-QCDDijet_JetQA-vz-pthat_7.12.16.pdf";
const std::string radius="4";
const std::string pfRad_etaWidth="R"+std::to_string(radius)"_"+(std::string)etaWidth;

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


  // evt plots----------------------
  //vz
  {
    std::string theHistName="hVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //weighted vz
  {
    std::string theHistName="hWeightedVz";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //pthat
  {
    std::string theHistName="hpthat";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //weighted pthat
  {
    std::string theHistName="hWeightedpthat";
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //----------------------



  // JetQA plots ----------------------
  //{1,0}wJetID
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



  // unfolding plots----------------------
  //unfolding inputs
  {  
    std::string theHistName="hpp_gen_"+pfRad_etaWidth;
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {  
    std::string theHistName="hpp_reco_"+pfRad_etaWidth;
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {  
    std::string theHistName="hpp_matrix_"+pfRad_etaWidth;
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH2F* theJetQAHist= (TH2F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  //MCClosure inputs
  {  
    std::string theHistName="hpp_mcclosure_gen_"+pfRad_etaWidth;
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {  
    std::string theHistName="hpp_mcclosure_data_"+pfRad_etaWidth;
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {  
    std::string theHistName="hpp_mcclosure_data_train_"+pfRad_etaWidth;
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }
  {  
    std::string theHistName="hpp_mcclosure_matrix_"+pfRad_etaWidth;
    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
    TH2F* theJetQAHist= (TH2F*)fin->Get( theHistName.c_str() );
    theJetQAHist->Draw();
    temp_canv->Print( thePDFFileName.c_str() );
  }



//  // Jet Corrections----------------------
//  //JEC
//  {  
//    std::string theHistName="";
//    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
//    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
//    theJetQAHist->Draw();
//    temp_canv->Print( thePDFFileName.c_str() );
//  }
//  //JER
//  {  
//    std::string theHistName="";
//    if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
//    TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
//    theJetQAHist->Draw();
//    temp_canv->Print( thePDFFileName.c_str() );
//  }



  fin->Close();
  temp_canv->Print( close_thePDFFileName.c_str() );//close the output file
  return 0;
}

