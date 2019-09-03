#include "printPlots.h"
#include "TVirtualFitter.h"
#include <Minuit2/Minuit2Minimizer.h>



// code/job switches ------------------------
const bool debugMode=true;

// the macro ------------------------
int printPlots_jetTrigEff ( const std::string input_ppData_condorDir ,   const std::string input_ppData_filename, 
			    const std::string output_PDFname_base, const std::string targTrig="", const std::string refTrig=""
			    ){
  
  //globalTrigHistStyle();
  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  TH1::SetDefaultSumw2(true);

  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  std::string fullJetType, radius;
  makeFullJetTypeString(&fullJetType, &radius, input_ppData_condorDir);
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl<<std::endl;;

  // OPEN INPUT SECTION
  //TFile *finDataJet80=NULL;    
  //TFile *finDataLowJets=NULL;    
  TFile *finData=NULL;    
  if(input_ppData_filename.find("MinBias")!=std::string::npos)
    std::cout<<"running jetTrigEff for MinBias PD"<<std::endl;
  else if(input_ppData_filename.find("LowerJets")!=std::string::npos)
    std::cout<<"running jetTrigEff for LowerJets PD"<<std::endl;
  else if(input_ppData_filename.find("Jet80")!=std::string::npos)
    std::cout<<"running jetTrigEff for Jet80 PD"<<std::endl;
  else if(input_ppData_filename.find("HighPtJetTrig")!=std::string::npos)
    std::cout<<"running jetTrigEff for HighPtJetTrig"<<std::endl;

  if(input_ppData_filename.find("wMB")!=std::string::npos)
    std::cout<<"MinBias sample is in the input file"<<std::endl;  
  std::cout<<std::endl;
  
  const std::string ppData_fullFilename=inputDir+input_ppData_condorDir+input_ppData_filename;
  
  //open the input files
  std::cout<<" now opening input data filename: "<<std::endl<<input_ppData_filename<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_ppData_condorDir<<std::endl;
  //if(debugMode)std::cout<<" ppData_fullFilename = "<<ppData_fullFilename	   <<std::endl;
  finData = new TFile(ppData_fullFilename.c_str());      
  if(!finData){
    std::cout << " DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppData_condorDir =" << input_ppData_condorDir << std::endl;
    std::cout << "input_ppData_filename  =" << input_ppData_filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);
  }
  std::cout<<std::endl;
  
  ////// GET OUTPUT FILE(S) READY
  std::cout<<" now opening output files in directory: "<<outputDir<<std::endl;
  //root file
  std::string theROOTFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetTrigEff.root";
  std::cout<<"theROOTFileName="<<theROOTFileName<<std::endl;
  TFile* fout = new TFile(theROOTFileName.c_str(), "RECREATE");
  //pdf file
  std::string thePDFFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetTrigEff.pdf";
  std::cout<<"thePDFFileName="<<thePDFFileName<<std::endl;
  std::string open_thePDFFileName=thePDFFileName+"[";    std::string close_thePDFFileName=thePDFFileName+"]";
  std::cout<<std::endl;

  // just for opening the pdf
  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 1000);
  temp_canv->Print( open_thePDFFileName.c_str() );    
  std::cout<<std::endl;
  
  //i want to look at eff curves for all possible combos of etabins
  //e.g. i=0,j=0 --> just 00eta05
  //e.g. i=0,j=2 --> 00eta15, etc.
  
  
  //----------------------------------------------------------------------//
//  //REFERENCE TRIGGER EFFICIENCY 0.0<|Y|<0.5
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    0, 0, targTrig, refTrig,
//  		    "reference");
//  //REFERENCE TRIGGER EFFICIENCY 0.5<|Y|<1.0
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    1, 1, targTrig, refTrig,
//  		    "reference");
//  //REFERENCE TRIGGER EFFICIENCY 1.0<|Y|<1.5
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    2, 2, targTrig, refTrig,
//  		    "reference");
//  //REFERENCE TRIGGER EFFICIENCY 1.5<|Y|<2.0
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    3, 3, targTrig, refTrig,
//  		    "reference");
  //REFERENCE TRIGGER EFFICIENCY 0.0<|Y|<1.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
  		    0, 1, targTrig, refTrig,
  		    "reference");
  //REFERENCE TRIGGER EFFICIENCY 1.0<|Y|<2.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
		    2, 3, targTrig, refTrig,
		    "reference");
  //REFERENCE TRIGGER EFFICIENCY |Y|<2.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
		    0, 3, targTrig, refTrig,
		    "reference");
  //----------------------------------------------------------------------//
//  //TRIGGER EMULATION EFFICIENCY 0.0<|Y|<0.5
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    0, 0, targTrig, refTrig,
//  		    "emulation");
//  //TRIGGER EMULATION EFFICIENCY 0.5<|Y|<1.0
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    1, 1, targTrig, refTrig,
//  		    "emulation");
//  //TRIGGER EMULATION EFFICIENCY 1.0<|Y|<1.5
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    2, 2, targTrig, refTrig,
//  		    "emulation");
//  //TRIGGER EMULATION EFFICIENCY 1.5<|Y|<2.0
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    3, 3, targTrig, refTrig,
//  		    "emulation");
  //TRIGGER EMULATION EFFICIENCY 0.0<|Y|<1.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
  		    0, 1, targTrig, refTrig,
  		    "emulation");
  //TRIGGER EMULATION EFFICIENCY 1.0<|Y|<2.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
  		    2, 3, targTrig, refTrig,
  		    "emulation");
  //TRIGGER EMULATION EFFICIENCY |Y|<2.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
		    0, 3, targTrig, refTrig,
		    "emulation");
  //----------------------------------------------------------------------//
//  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 0.0<|Y|<0.5
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    0, 0, targTrig, refTrig,
//		    "emulation","wL1Thresh");
//  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 0.5<|Y|<1.0
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    1, 1, targTrig, refTrig,
//		    "emulation","wL1Thresh");
//  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 1.0<|Y|<1.5
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    2, 2, targTrig, refTrig,
//  		    "emulation","wL1Thresh");
//  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 1.5<|Y|<2.0
//  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
//  		    3, 3, targTrig, refTrig,
//  		    "emulation","wL1Thresh");
  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 0.0<|Y|<1.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
  		    0, 1, targTrig, refTrig,
		    "emulation","wL1Thresh");
  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 1.0<|Y|<2.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
  		    2, 3, targTrig, refTrig,
  		    "emulation","wL1Thresh");
  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD |Y|<2.0
  drawAndFitTrigEff(finData, thePDFFileName, fullJetType, radius, fout,
		    0, 3, targTrig, refTrig,
		    "emulation","wL1Thresh");
  //----------------------------------------------------------------------//  
  
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  temp_canv->Print( close_thePDFFileName.c_str() );  
  temp_canv->Close();    
  
  
  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;

  finData->Close();  
  fout->Close();
  
  return 0;
}// end printplots




// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=6 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetTrigEff.exe <ppDataCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  
  if(argc==6) rStatus=printPlots_jetTrigEff( (const std::string) argv[1], (const std::string) argv[2], (const std::string)argv[3], (const std::string)argv[4], (const std::string)argv[5]);
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
  
}
