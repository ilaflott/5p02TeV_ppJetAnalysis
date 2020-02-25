#include "printPlots.h"
#include "TVirtualFitter.h"
#include <Minuit2/Minuit2Minimizer.h>



// code/job switches ------------------------
const bool debugMode=true;

// the macro ------------------------
int printPlots_jetTrigEff ( const std::string input_ppDataTrgEff_condorDir ,   const std::string input_ppData_condorDir ,   
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
  makeFullJetTypeString(&fullJetType, &radius, input_ppDataTrgEff_condorDir);
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl<<std::endl;;

  // OPEN INPUT SECTION
  TFile *finDataJet80=NULL;    
  TFile *finDataLowJets=NULL;    
  //TFile *finDataMinBias=NULL;    
  TFile *finData=NULL;    
  //  if(input_ppData_filename.find("MinBias")!=std::string::npos)
  //    std::cout<<"running jetTrigEff for MinBias PD"<<std::endl;
  //  else if(input_ppData_filename.find("LowerJets")!=std::string::npos)
  //    std::cout<<"running jetTrigEff for LowerJets PD"<<std::endl;
  //  else if(input_ppData_filename.find("Jet80")!=std::string::npos)
  //    std::cout<<"running jetTrigEff for Jet80 PD"<<std::endl;
  //  std::cout<<std::endl;
  
  const std::string input_ppDataJet80_filename="HighPtJet80_ak4PF-allFiles.root";
  const std::string ppDataJet80_fullFilename=inputDir+input_ppDataTrgEff_condorDir+input_ppDataJet80_filename;
  const std::string input_ppDataLowJets_filename="HighPtLowerJets_ak4PF-allFiles.root";
  const std::string ppDataLowJets_fullFilename=inputDir+input_ppDataTrgEff_condorDir+input_ppDataLowJets_filename;
  //const std::string input_ppDataMinBias_filename="HighPtMinBias_ak4PF-allFiles.root";
  //const std::string ppDataMinBias_fullFilename=inputDir+input_ppDataTrgEff_condorDir+input_ppDataMinBias_filename;
  const std::string input_ppData_filename="HighPtJetTrig_noMB_ak4PF-allFiles.root";
  const std::string ppData_fullFilename=inputDir+input_ppData_condorDir+input_ppData_filename;
  
  //open the input files
  std::cout<<" in directory: "<<inputDir+input_ppDataTrgEff_condorDir<<std::endl  <<std::endl;
  std::cout<<" now opening input data Jet80 PD filename: "<<std::endl<<input_ppDataJet80_filename<<std::endl;
  //finDataJet80 = new TFile(ppDataJet80_fullFilename.c_str());      
  finDataJet80 = TFile::Open(ppDataJet80_fullFilename.c_str(),"READ");
  if(finDataJet80->IsZombie()){
    std::cout << " JET80 DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppDataJet80_condorDir =" << input_ppDataTrgEff_condorDir << std::endl;
    std::cout << "input_ppDataJet80_filename  =" << input_ppDataJet80_filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);
  }
  std::cout<<std::endl;
  std::cout<<" now opening input data LowerJets PD filename: "<<std::endl<<input_ppDataLowJets_filename<<std::endl;
  //finDataLowJets = new TFile(ppDataLowJets_fullFilename.c_str());     
  finDataLowJets = TFile::Open(ppDataLowJets_fullFilename.c_str(),"READ"); 
  if(finDataLowJets->IsZombie()){
    std::cout << " LOWERJETS DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppDataLowJets_condorDir =" << input_ppDataTrgEff_condorDir << std::endl;
    std::cout << "input_ppDataLowJets_filename  =" << input_ppDataLowJets_filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);
  }
  std::cout<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_ppData_condorDir<<std::endl  <<std::endl;
  std::cout<<" now opening input data Trigger Combined PD filename: "<<std::endl<<input_ppData_filename<<std::endl;
  finData = TFile::Open(ppData_fullFilename.c_str(),"READ");
  //finData = new TFile(ppData_fullFilename.c_str());      
  if(finData->IsZombie()){
    std::cout << " COMBINED DATA file not found" << std::endl;
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

  
  int trigpt_rebinfact=1;
  if(     targTrig=="HLT80")
    trigpt_rebinfact=1;
  else if(targTrig=="HLT60")
    trigpt_rebinfact=1;






  //----------------------------------------------------------------------//
  //----------------- TURN ONS USING LEAD TRIGGER JET PT -----------------//
  //----------------------------------------------------------------------//
  //REFERENCE TRIGGER EFFICIENCY 0.0<|Y|<5.1
  drawAndFitTrigEff_trigPt(finDataJet80, finDataLowJets, 
			   fullJetType, radius, trigpt_rebinfact,
			   targTrig, refTrig, "reference", "", 
			   thePDFFileName, fout);  
  //REFERENCE TRIGGER EFFICIENCY 0.0<|Y|<2.0
  drawAndFitTrigEff_trigPt(finDataJet80, finDataLowJets, 
			   fullJetType, radius, trigpt_rebinfact,
			   targTrig, refTrig, "reference", "", 
			   thePDFFileName, fout, "_00eta20");  
  //REFERENCE TRIGGER EFFICIENCY 2.0<|Y|<5.1
  drawAndFitTrigEff_trigPt(finDataJet80, finDataLowJets, 
			   fullJetType, radius, trigpt_rebinfact,
			   targTrig, refTrig, "reference", "", 
			   thePDFFileName, fout, "_20eta51");  


  //////TRIGGER EMULATION W L1 THRESHOLD EFFICIENCY 0.0<|Y|<5.1//too much self-correlation to be useful
  //drawAndFitTrigEff_trigPt(finDataJet80, finDataLowJets, 
  //			   fullJetType, radius, trigpt_rebinfact,
  //			   targTrig, refTrig, "emulation", "wL1Thresh", 
  //			   thePDFFileName, fout);  
  //----------------------------------------------------------------------//  
  
  
  
  //----------------------------------------------------------------------//
  //----------------- TURN ONS USING INCL JET XSEC -----------------------//
  //----------------------------------------------------------------------//
  //i want to look at eff curves for all possible combos of etabins
  //e.g. i=0,j=0 --> just 00eta05
  //e.g. i=0,j=2 --> 00eta15, etc.
  //----------------------------------------------------------------------//
  //REFERENCE TRIGGER EFFICIENCY  0.0<|Y|<1.0
  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
		    fullJetType, radius, 0, 1, 
		    targTrig, refTrig, "reference","",
		    thePDFFileName, fout);
  //REFERENCE TRIGGER EFFICIENCY  1.0<|Y|<2.0
  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
		    fullJetType, radius, 2, 3, 
		    targTrig, refTrig, "reference","",
		    thePDFFileName, fout);
  //REFERENCE TRIGGER EFFICIENCY  |Y|<2.0
  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
		    fullJetType, radius, 0, 3, 
		    targTrig, refTrig, "reference","",
		    thePDFFileName, fout);
  //----------------------------------------------------------------------//
  //  //TRIGGER EMULATION EFFICIENCY  0.0<|Y|<1.0
  //  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
  //		    fullJetType, radius, 0, 1, 
  //		    targTrig, refTrig, "emulation","",
  //		    thePDFFileName, fout);
  //  //TRIGGER EMULATION EFFICIENCY  1.0<|Y|<2.0
  //  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
  //		    fullJetType, radius, 2, 3, 
  //		    targTrig, refTrig, "emulation","",
  //		    thePDFFileName, fout);
  //  //TRIGGER EMULATION EFFICIENCY  |Y|<2.0
  //  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
  //		    fullJetType, radius, 0, 3, 
  //		    targTrig, refTrig, "emulation","",
  //		    thePDFFileName, fout);
  //  //----------------------------------------------------------------------//
  //  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 0.0<|Y|<1.0
  //  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
  //		    fullJetType, radius, 0, 1, 
  //		    targTrig, refTrig, "emulation","wL1Thresh",
  //		    thePDFFileName, fout);
  //  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD 1.0<|Y|<2.0
  //  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
  //		    fullJetType, radius, 2, 3, 
  //		    targTrig, refTrig, "emulation","wL1Thresh",
  //		    thePDFFileName, fout);
  //  //TRIGGER EMULATION EFFICIENCY W L1 THRESHOLD |Y|<2.0
  //  drawAndFitTrigEff(finDataJet80, finDataLowJets, 
  //		    fullJetType, radius, 0, 3, 
  //		    targTrig, refTrig, "emulation","wL1Thresh",
  //		    thePDFFileName, fout);
  //----------------------------------------------------------------------//  

  

  //globalTrigHistStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  //TH1::SetDefaultSumw2(true);

  //-------------------------------------------------------------------------//
  //----------------- LEAD TRIGGER JET OBJECT SPECTRA -----------------------//
  //-------------------------------------------------------------------------//
  //COMBO w/ EXCLUSIVE TRIGGER SPECTRA
  printTrigPtHist(finDataLowJets, finDataJet80, finData, 
    		  fullJetType, radius, "excl", 
    		  thePDFFileName, (TFile*) fout);
  printTrigPtHist(finDataLowJets, finDataJet80, finData, 
    		  fullJetType, radius, "excl", 
    		  thePDFFileName, (TFile*) fout,"_00eta20");
  printTrigPtHist(finDataLowJets, finDataJet80, finData, 
    		  fullJetType, radius, "excl", 
    		  thePDFFileName, (TFile*) fout,"_20eta51");
  //-------------------------------------------------------------------------------//
  //COMBO w/ INCLUSIVE TRIGGER SPECTRA
  printTrigPtHist(finDataLowJets, finDataJet80, finData, 
  		  fullJetType, radius, "incl", 
  		  thePDFFileName, (TFile*) fout);  
  //printTrigPtHist(finDataLowJets, finDataJet80, finData,  //these don't make sense for... reasons
  //fullJetType, radius, "incl", 
  //thePDFFileName, (TFile*) fout,"_00eta20");
  //printTrigPtHist(finDataLowJets, finDataJet80, finData, 
  //fullJetType, radius, "incl", 
  //thePDFFileName, (TFile*) fout,"_20eta51");
  //----------------------------------------------------------------------//  


  //-----------------------------------------------------------------------//
  //----------------- INLCUSIVE TRIGGER JET SPECTRA -----------------------//
  //-----------------------------------------------------------------------//
  //COMBO w/ EXCLUSIVE TRIGGER SPECTRA
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  for(int etabin=0;etabin<4;etabin++)
    printJetTrigHist_wRatio(finDataLowJets, finDataJet80, finData, 
  			    fullJetType, radius, "excl", etabin,			    
  			    thePDFFileName, (TFile*) fout);
  //-----------------------------------------------------------------------//
  //  //COMBO w/ INCLUSIVE TRIGGER SPECTRA
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  for(int etabin=0;etabin<4;etabin++)
    printJetTrigHist_wRatio(finDataLowJets, finDataJet80, finData, 
			    fullJetType, radius, "incl", etabin,			    
    			    thePDFFileName, (TFile*) fout);
  //-----------------------------------------------------------------------//

  //-----------------------------------------------------------------------------//
  //----------------- INLCUSIVE TRIGGER JET SPECTRA V RUN -----------------------//
  //-----------------------------------------------------------------------------//
  //COMBO w/ EXCLUSIVE TRIGGER SPECTRA
  //for(int runind=0;runind<Nruns;runind++){
  //  printJetTrigHist_wRatio_vRun(finData, false,  runind,
  //				   thePDFFileName, fullJetType, "excl", radius, usedHLTPF,  (TFile*) fout);
  //}
  //-----------------------------------------------------------------------------//
  //COMBO w/ INCLUSIVE TRIGGER SPECTRA
  //for(int runind=0;runind<Nruns;runind++){
  //  printJetTrigHist_wRatio_vRun(finData,  false,  runind,
  //				   thePDFFileName, fullJetType, "incl", radius, usedHLTPF,  (TFile*) fout);	
  //}
  //-----------------------------------------------------------------------------//
  
  //----------------------------------------------------------------------------------------//
  //----------------- LEAD TRIGGER JET OBJECT ETA DIST AND ETA ASYMM -----------------------//
  //----------------------------------------------------------------------------------------//
  //COMBO w/ EXCLUSIVE TRIGGER ETA DIST
  //    printTrigEtaHist(finData, usedHLT100,
  //		     thePDFFileName, fullJetType, "excl", radius, usedHLTPF , (TFile*) fout);    
  //    printTrigEtaAsymmHist(finData, usedHLT100,
  //			  thePDFFileName, fullJetType, "excl", radius, usedHLTPF , (TFile*) fout);            
  //----------------------------------------------------------------------------------------//
  //COMBO w/ INCLUSIVE TRIGGER ETA DIST
  //    printTrigEtaHist(finData, usedHLT100, 
  //		     thePDFFileName, fullJetType, "incl", radius, usedHLTPF , (TFile*) fout);    
  //    printTrigEtaAsymmHist(finData, usedHLT100,
  //			  thePDFFileName, fullJetType, "incl", radius, usedHLTPF , (TFile*) fout);    
  //----------------------------------------------------------------------------------------//
  
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  temp_canv->Print( close_thePDFFileName.c_str() );  
  temp_canv->Close();    
  
  
  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;

  finDataJet80->Close();  
  finDataLowJets->Close();  
  finData->Close();  
  fout->Close();
  
  return 0;
}// end printplots




// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=6 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetTrigEff.exe [actually read the code i'm lazy right now sorry future ian]"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  
  if(argc==6) rStatus=printPlots_jetTrigEff( (const std::string) argv[1], (const std::string) argv[2], (const std::string)argv[3], (const std::string)argv[4], (const std::string)argv[5]);
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
  
}
