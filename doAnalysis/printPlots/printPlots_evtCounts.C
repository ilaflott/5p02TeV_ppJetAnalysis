//const int Netabins=4;
#include "printPlots.h"


// code/job switches ------------------------

//options
const bool debugMode=false;


// the macro ------------------------
int printPlots_evtCounts(const std::string input_ppData_condorDir , const std::string input_ppMC_condorDir , 
			const std::string output_PDFname_base){
  

  //  globalHistStyle();
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=input_ppData_condorDir.find("_ak")+3;  
  const std::string radiusInt= input_ppData_condorDir.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=input_ppData_condorDir.find("Jets");;
  
  //make final jetType strings
  const std::string jetType=input_ppData_condorDir.substr( jetTypePos,(jetsPos-jetTypePos) );
  const std::string fullJetType="ak"+radiusInt+jetType;
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;



  //put together input file strings
  const std::string input_ppData_Filename_noMB="HighPtJetTrig_noMB_" +fullJetType+ "-allFiles.root";
  std::string input_ppData_Filename_Jet80="HighPtJet80_" +fullJetType+ "-allFiles.root";
  std::string input_ppData_Filename_LowJets="HighPtLowerJets_" +fullJetType+ "-allFiles.root";
  
  //const std::string input_ppData_Filename="HighPtJet80_" +fullJetType+ "-allFiles.root";
  //const std::string input_ppData_Filename="HighPtLowerJets_" +fullJetType+ "-allFiles.root";
  const std::string input_ppMC_Filename  ="Py8_CUETP8M1_QCDjetAllPtBins_" +fullJetType+ "-allFiles.root";
  
  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;  
  const std::string ppData_fullFilename_Jet80  =inputDir+input_ppData_condorDir+input_ppData_Filename_Jet80;
  const std::string ppData_fullFilename_LowJets=inputDir+input_ppData_condorDir+input_ppData_Filename_LowJets;
  const std::string ppMC_fullFilename  =inputDir+input_ppMC_condorDir+input_ppMC_Filename;
  
  // OPEN INPUT SECTION    
  TFile *finData_Jet80=NULL;  
  std::cout<<" now opening ppData: "<<std::endl<<input_ppData_Filename_Jet80<<std::endl;
  std::cout<<" in directory: "<<input_ppData_condorDir<<std::endl<<std::endl;
  finData_Jet80 = new TFile(ppData_fullFilename_Jet80.c_str(), "READ");      
  if(finData_Jet80->IsZombie()){
    std::cout << " DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppData_condorDir =" << input_ppData_condorDir << std::endl;
    std::cout << "input_ppData_filename  =" << input_ppData_Filename_Jet80 << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);      }  


  TFile *finData_LowJets=NULL;  
  std::cout<<" now opening ppData: "<<std::endl<<input_ppData_Filename_LowJets<<std::endl;
  std::cout<<" in directory: "<<input_ppData_condorDir<<std::endl<<std::endl;
  finData_LowJets = new TFile(ppData_fullFilename_LowJets.c_str(), "READ");      
  if(finData_LowJets->IsZombie()){
    std::cout << " DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppData_condorDir =" << input_ppData_condorDir << std::endl;
    std::cout << "input_ppData_filename  =" << input_ppData_Filename_LowJets << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);      }  
  
  TFile *finMC=NULL;
  std::cout<<" now opening ppMC: "<<std::endl<<input_ppMC_Filename<<std::endl;
  std::cout<<" in directory: "<<input_ppMC_condorDir<<std::endl<<std::endl;  
  finMC = new TFile(ppMC_fullFilename.c_str(), "READ");      
  if(finMC->IsZombie()){
    std::cout << " MC file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppMC_condorDir =" << input_ppMC_condorDir << std::endl;
    std::cout << "input_ppMC_filename  =" << input_ppMC_Filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);     }
  //assert(false);
  
  
    

  
  // GET OUTPUT ROOT FILE READY
  std::string theROOTFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_evtCounts.root";
  TFile* fout=new TFile(theROOTFileName.c_str(),"RECREATE");
  
  TH1D* h_NEvents_Jet80_read        =(TH1D*)finData_Jet80->Get("NEvents_read"   );
  TH1D* h_NEvents_Jet80_skimCut     =(TH1D*)finData_Jet80->Get("NEvents_skimCut");
  TH1D* h_NEvents_Jet80_vzCut       =(TH1D*)finData_Jet80->Get("NEvents_vzCut"  );
  TH1D* h_NEvents_Jet80_PFMETfracCut=(TH1D*)finData_Jet80->Get("NEvents_PFMETfracCut"  );
  TH1D* h_NEvents_Jet80_is80        =(TH1D*)finData_Jet80->Get("NEvents_is80"   );
  h_NEvents_Jet80_read   ->SetName("NEvents_Jet80_read"   );
  h_NEvents_Jet80_skimCut->SetName("NEvents_Jet80_skimCut");
  h_NEvents_Jet80_vzCut  ->SetName("NEvents_Jet80_vzCut"  );
  h_NEvents_Jet80_PFMETfracCut  ->SetName("NEvents_Jet80_PFMETfracCut"  );  
  h_NEvents_Jet80_is80   ->SetName("NEvents_Jet80_is80"   );

  std::cout<<std::endl;
  std::cout<<" ----- Jet80 PD Event Count Report ----- "<<std::endl;
  std::cout<<"# Events ............................................ in Dataset: "<< h_NEvents_Jet80_read->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<< h_NEvents_Jet80_skimCut->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after ..................................|vz|< 24 cm cut: "<< h_NEvents_Jet80_vzCut->GetBinContent(1)<<std::endl;
  std::cout<<"               Evt. Selection Efficiency 1 =        NEvents_vzCut/NEvents = "<< h_NEvents_Jet80_vzCut->GetBinContent(1)/h_NEvents_Jet80_read->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after ..................................PFMETfrac < 0.3: "<< h_NEvents_Jet80_PFMETfracCut->GetBinContent(1)<<std::endl;
  std::cout<<"               Evt. Selection Efficiency 2 = NEvents_PFMETfracCut/NEvents = "<< h_NEvents_Jet80_PFMETfracCut->GetBinContent(1)/h_NEvents_Jet80_read->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after ............... passing trigger jet req for HLT80: "<< h_NEvents_Jet80_is80 ->GetBinContent(1)<<std::endl;


  TH1D* h_NEvents_LowJets_read   =(TH1D*)finData_LowJets->Get("NEvents_read"   );
  TH1D* h_NEvents_LowJets_skimCut=(TH1D*)finData_LowJets->Get("NEvents_skimCut");
  TH1D* h_NEvents_LowJets_vzCut  =(TH1D*)finData_LowJets->Get("NEvents_vzCut"  );
  TH1D* h_NEvents_LowJets_PFMETfracCut=(TH1D*)finData_LowJets->Get("NEvents_PFMETfracCut"  );
  TH1D* h_NEvents_LowJets_is40   =(TH1D*)finData_LowJets->Get("NEvents_is40"   );
  TH1D* h_NEvents_LowJets_is60   =(TH1D*)finData_LowJets->Get("NEvents_is60"   );
  h_NEvents_LowJets_read   ->SetName("NEvents_LowJets_read"   );
  h_NEvents_LowJets_skimCut->SetName("NEvents_LowJets_skimCut");
  h_NEvents_LowJets_vzCut  ->SetName("NEvents_LowJets_vzCut"  );
  h_NEvents_LowJets_PFMETfracCut  ->SetName("NEvents_LowJets_PFMETfracCut"  );
  h_NEvents_LowJets_is40   ->SetName("NEvents_LowJets_is40"   );
  h_NEvents_LowJets_is60   ->SetName("NEvents_LowJets_is60"   );
  std::cout<<std::endl;
  std::cout<<" ----- LowerJets PD Event Count Report ----- "<<std::endl;
  std::cout<<"# Events ............................................ in Dataset: "<< h_NEvents_LowJets_read->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<< h_NEvents_LowJets_skimCut->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after ..................................|vz|< 24 cm cut: "<< h_NEvents_LowJets_vzCut->GetBinContent(1)<<std::endl;
  std::cout<<"               Evt. Selection Efficiency 1 =        NEvents_vzCut/NEvents = "<< h_NEvents_LowJets_vzCut->GetBinContent(1)/h_NEvents_LowJets_read->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after ..................................PFMETfrac < 0.3: "<< h_NEvents_LowJets_PFMETfracCut->GetBinContent(1)<<std::endl;
  std::cout<<"               Evt. Selection Efficiency 2 = NEvents_PFMETfracCut/NEvents = "<< h_NEvents_LowJets_PFMETfracCut->GetBinContent(1)/h_NEvents_LowJets_read->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after ...... passing trigger jet req for HLT60 or HLT40: "<< h_NEvents_LowJets_is40 ->GetBinContent(1)+h_NEvents_LowJets_is60 ->GetBinContent(1)<<std::endl;
  std::cout<<"# Events ......................... passing trigger jet req for HLT40 only: "<< h_NEvents_LowJets_is40 ->GetBinContent(1)<<std::endl;
  std::cout<<"# Events ......................... passing trigger jet req for HLT60 only: "<< h_NEvents_LowJets_is60 ->GetBinContent(1)<<std::endl;
  
  TH1D* h_NEvents_MC_read   =(TH1D*)finMC->Get("NEvents_read"   );
  TH1D* h_NEvents_MC_skimCut=(TH1D*)finMC->Get("NEvents_skimCut");
  TH1D* h_NEvents_MC_vzCut  =(TH1D*)finMC->Get("NEvents_vzCut"  );
  h_NEvents_MC_read   ->SetName("NEvents_MC_read"   );
  h_NEvents_MC_skimCut->SetName("NEvents_MC_skimCut");
  h_NEvents_MC_vzCut  ->SetName("NEvents_MC_vzCut"  );

  std::cout<<std::endl;
  std::cout<<" ----- PYTHIA8 PD Event Count Report ----- "<<std::endl;
  std::cout<<"# Events ............................................ in Dataset: "<< h_NEvents_MC_read->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after HBHENoise, BeamScraping, and PV  quality  filters: "<< h_NEvents_MC_skimCut->GetBinContent(1)<<std::endl;
  std::cout<<"# Events after ..................................|vz|< 24 cm cut: "<< h_NEvents_MC_vzCut->GetBinContent(1)<<std::endl;

  fout->cd();

  h_NEvents_Jet80_read   ->Write();
  h_NEvents_Jet80_skimCut->Write();  
  h_NEvents_Jet80_vzCut  ->Write();  
  h_NEvents_Jet80_PFMETfracCut  ->Write();  
  h_NEvents_Jet80_is80   ->Write();

  h_NEvents_LowJets_read   ->Write();
  h_NEvents_LowJets_skimCut->Write();
  h_NEvents_LowJets_vzCut  ->Write();
  h_NEvents_LowJets_PFMETfracCut  ->Write();  
  h_NEvents_LowJets_is40   ->Write();
  h_NEvents_LowJets_is60   ->Write();

  h_NEvents_MC_read   ->Write();
  h_NEvents_MC_skimCut->Write();
  h_NEvents_MC_vzCut  ->Write();



//  // --------------- JET COUNTS SECTION
//  
//  
//  TH1D* h_NJets_Jet80	=(TH1D*)finData_Jet80->Get("NJets"   );
//  TH1D* h_NJets_Jet80_00y20_jtptLo56_0wjetID=NULL;
//  TH1D* h_NJets_Jet80_00y20_jtptLo56_1wjetID;
//  
//  for(int i=0; i<4;i++){
//  }
//
//
//
//
//
//
//
//  fout ->cd();
//  
//  h_NJets_Jet80_jtptCut	  ->Write();
//  h_NJets_Jet80_jtptCut_Hi->Write();
//  h_NJets_Jet80_jtetaCut1 ->Write();
//  h_NJets_Jet80_jtetaCut2 ->Write();
//  h_NJets_Jet80_jetQAPtCut->Write();
//  h_NJets_Jet80_JetIDCut  ->Write();
//
//
//  h_NJets_LowJets_jtptCut   ->Write();
//  h_NJets_LowJets_jtptCut_Hi->Write();
//  h_NJets_LowJets_jtetaCut1 ->Write();
//  h_NJets_LowJets_jtetaCut2 ->Write();
//  h_NJets_LowJets_jetQAPtCut->Write();
//  h_NJets_LowJets_JetIDCut  ->Write();
//





  
  
  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData_Jet80->Close();  
  finData_LowJets->Close();  

  if(debugMode)std::cout<<std::endl<<"closing output ROOT"<<std::endl;
  fout->Close();
  return 0;
  
}// end printplots

				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=4 ) {//no input arguments, error
    std::cout<<"do ./printPlots_evtCount.exe <ppDataCondorDir> <ppMCCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  if(argc==4) rStatus=printPlots_evtCounts( (const std::string) argv[1], (const std::string) argv[2],
					(const std::string) argv[3]);
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
   
}



///////////////////////
//// PROBLEM SECTION, do NOT delete comments or other commented code
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//
//if(debugMode)std::cout<<std::endl<<"N_trigs="<<N_trigs<<" and j="<<j<<std::endl;
//if(debugMode)std::cout<<"j<(N_trigs-1)="<<(bool)(j<(N_trigs-1))<<std::endl;
//if(debugMode)std::cout<<"j<N_trigs="<<(bool)(j<(N_trigs))<<std::endl;
//if(debugMode)std::cout<<"j<(N_trigs+1)="<<(bool)(j<(N_trigs+1))<<std::endl;	
//
////assert(j!=N_trigs);// NEVER FIRES
////assert(j<N_trigs); // NEVER FIRES
//
////open the hists + do scaling 
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//
//if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
//  std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
//  continue;}  
//else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
//
////if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
////  std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
////  continue;}	
////else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
//  
////open the hists + do scaling 
//if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
//// PROBLEM SECTION
///////////////////////








//const bool drawJetTrigQAPlots=true, drawJetRapBinsPlot=false;
//const bool drawDataMCOverlaysInput=false; 		 //this should always be false, but needs to be cleaned up later.
//const bool drawDataMCOverlays   = drawDataMCOverlaysInput;
//const bool drawDataMCRatios     = true;   //!drawDataMCOverlays;
//const bool drawTrigDataOverlays = false;  //drawDataMCOverlays && drawJetTrigQAPlots;
//const bool drawTrigDataRatios   = false; //drawDataMCRatios   && drawJetTrigQAPlots;
//const int theTrigOverlayLineColor[]  ={  1,  1,  1,  1,  1,  1 };
//const int theTrigCombMarkerColor=1, altTrigCombMarkerColor=12;
//const int theTrigOverlayMarkerColor[]={  2,  3,  6,  7,  1,  4 };
//const int theRapOverlayMarkerColor[] ={  2,  3,  6,  7,  1,  4 , 8, 9};
//const int theTrigOverlayMarker[]     ={ 20, 20, 20, 20, 20, 32 };
//const std::string crossSectionAxTitle="#sigma (#mub)";
//const std::string AUAxTitle="A.U.";
//const std::string ratioTitle="MC/Data";
