#include "DijetResponse.h"

const std::string CMSSW_BASE =
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";
const std::string inputDir =
  "src/DijetResponse/saved_outputCondor/";

const std::string inputDataFilename = 
  "deriveResponse_5p02TeV_HighPtJetTrig_2016-06-11__Jet80PFBug_allResponses/HighPtJetTrig_ak4PF-sumResponse_RelMPFandMPFAbsPho.root";
const std::string inputDataFullFilename = CMSSW_BASE+inputDir+inputDataFilename;

const std::string inputMCFilename   = 
  "deriveResponse_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_merged_2016-06-22__allResponses/QCDJetAllPtBins_ak4PF-sumResponse_allResponses.root";
const std::string inputMCFullFilename = CMSSW_BASE+inputDir+inputMCFilename;

const bool printDataPlots=true;
const bool printMCPlots=true;

const std::string theDataPDFFileName="DijetResponse_printPlots_HighPtJetTrigs_allResponses_7.12.16.pdf";
const std::string theMCPDFFileName="DijetResponse_printPlots_PY6_TuneZ2-QCDDijet_allResponses_7.12.16.pdf";

int main(int argc, char *argv[]){

  const bool debugMode=true;
  int rStatus = -1;
  if( argc!=1 ) {
    std::cout<<"settings hard coded, just do ./DijetResponse_printPlots.exe"<<std::endl;
    return rStatus;
  }
  rStatus=1;

  // -----------------------------------------------------------------------
  // Data ------------------------------------------------------------------
  // -----------------------------------------------------------------------
  if(printDataPlots){
    //input Data file
    std::string filename=inputDataFullFilename;
    std::cout<< "filename is " << filename <<std::endl;
    TFile *fin = new TFile(filename.c_str());
  
    //create canvas for plots + open output PDF file
    TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
    std::string thePDFFileName=theDataPDFFileName;
    std::string open_thePDFFileName=theDataPDFFileName+"[";
    std::string close_thePDFFileName=theDataPDFFileName+"]";
    temp_canv->Print( open_thePDFFileName.c_str() );      //open the .pdf file we wan t to print to, no blank page at start
    // draw each plot on the temp canvas, then print to the pdf file

    // relative response ( individual pt bins, and all together)
    //----------------------
    {
      temp_canv->cd();    
      for(int i=0;i<5;++i){
	std::string theHistName="hRelResponse_pt"+std::to_string(i);
	if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
	TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
	theJetQAHist->Draw();
	temp_canv->Print( thePDFFileName.c_str() );
      }
    }
    {
      std::string theCanvName="hRelResponse_all";
      if(debugMode)std::cout<<"theCanvName="<<theCanvName<<std::endl;
      TCanvas* theAllRespCanv= (TCanvas*)fin->Get( theCanvName.c_str() );
      theAllRespCanv->cd();
      theAllRespCanv->Print( thePDFFileName.c_str() );
    }

    // MPFResponse ( individual pt bins, and all together)
    //----------------------
    {
      temp_canv->cd();    
      for(int i=0;i<5;++i){
	std::string theHistName="hMPFResponse_pt"+std::to_string(i);
	if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
	TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
	theJetQAHist->Draw();
	temp_canv->Print( thePDFFileName.c_str() );
      }
    }
    {
      std::string theCanvName="hMPFResponse_all";
      if(debugMode)std::cout<<"theCanvName="<<theCanvName<<std::endl;
      TCanvas* theAllRespCanv= (TCanvas*)fin->Get( theCanvName.c_str() );
      theAllRespCanv->cd();
      theAllRespCanv->Print( thePDFFileName.c_str() );
    }

    // MPFAbsResponse ( individual pt bins, and all together)
    //----------------------
    {
      temp_canv->cd();    
      for(int i=0;i<5;++i){
	std::string theHistName="hMPFAbsPhoResponse_pt"+std::to_string(i);
	if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
	TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
	theJetQAHist->Draw();
	temp_canv->Print( thePDFFileName.c_str() );
      }
    }
    {
      std::string theCanvName="hMPFAbsPhoResponse_all";
      if(debugMode)std::cout<<"theCanvName="<<theCanvName<<std::endl;
      TCanvas* theAllRespCanv= (TCanvas*)fin->Get( theCanvName.c_str() );
      theAllRespCanv->cd();
      theAllRespCanv->Print( thePDFFileName.c_str() );
    }

    fin->Close();
    temp_canv->Print( close_thePDFFileName.c_str() );//close the output file
  }//end print Data plots condition




  // ---------------------------------------------------------------------
  // MC ------------------------------------------------------------------
  // ---------------------------------------------------------------------
  if(printMCPlots){
    //input Data file
    std::string filename=inputMCFullFilename;
    std::cout<< "filename is " << filename <<std::endl;
    TFile *fin = new TFile(filename.c_str());
  
    //get output file ready
    TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);
    std::string thePDFFileName=theMCPDFFileName;
    std::string open_thePDFFileName=theMCPDFFileName+"[";
    std::string close_thePDFFileName=theMCPDFFileName+"]";
    temp_canv->Print( open_thePDFFileName.c_str() );      //open the .pdf file we want to print to
    //temp_canv->cd();
    
    
    // relative response ( individual pt bins, and all together)
    //----------------------
    {
      temp_canv->cd();
      for(int i=0;i<5;++i){
	std::string theHistName="hRelResponse_pt"+std::to_string(i);
	if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
	TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
	theJetQAHist->Draw();
	temp_canv->Print( thePDFFileName.c_str() );
      }
    }
    {
      std::string theCanvName="hRelResponse_all";
      if(debugMode)std::cout<<"theCanvName="<<theCanvName<<std::endl;
      TCanvas* theAllRespCanv= (TCanvas*)fin->Get( theCanvName.c_str() );
      theAllRespCanv->cd();
      theAllRespCanv->Print( thePDFFileName.c_str() );
    }
    //----------------------

    // MPFResponse ( individual pt bins, and all together)
    //----------------------
    {
      temp_canv->cd();
      for(int i=0;i<5;++i){
	std::string theHistName="hMPFResponse_pt"+std::to_string(i);
	if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
	TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
	theJetQAHist->Draw();
	temp_canv->Print( thePDFFileName.c_str() );
      }
    }
    {
      std::string theCanvName="hMPFResponse_all";
      if(debugMode)std::cout<<"theCanvName="<<theCanvName<<std::endl;
      TCanvas* theAllRespCanv= (TCanvas*)fin->Get( theCanvName.c_str() );
      theAllRespCanv->cd();
      theAllRespCanv->Print( thePDFFileName.c_str() );
    }
    //----------------------

    // MPFAbsResponse ( individual pt bins, and all together)
    //----------------------
    {
      temp_canv->cd();
      for(int i=0;i<5;++i){
	std::string theHistName="hMPFAbsPhoResponse_pt"+std::to_string(i);
	if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;
	TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );
	theJetQAHist->Draw();
	temp_canv->Print( thePDFFileName.c_str() );
      }
    }
    {
      std::string theCanvName="hMPFAbsPhoResponse_all";
      if(debugMode)std::cout<<"theCanvName="<<theCanvName<<std::endl;
      TCanvas* theAllRespCanv= (TCanvas*)fin->Get( theCanvName.c_str() );
      theAllRespCanv->cd();
      theAllRespCanv->Print( thePDFFileName.c_str() );
    }
    //----------------------    
    
    fin->Close();
    temp_canv->Print( close_thePDFFileName.c_str() );//close the output file
  }// end do MC plots condition
  // MC ------------------------------------------------------------------
  
  return 0;
}

