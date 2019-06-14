#include "printPlots.h"


// code/job switches ------------------------

//options
const bool debugMode=false, doEventCounts=false;//, doJetIDPlots=true;


//draw switches
const bool drawEvtQAPlots=true;
const bool drawJetQAPlots=true;
const bool drawTupelJetQAPlots=false;
const bool drawJetConstituentPlots=drawJetQAPlots&&false, drawDijetPlots=drawJetQAPlots&&false;
const bool drawJetRapBinsPlot=false;//, drawGENJetRapBinsPlot=true;  

const bool doJetIDPlots=true;
const bool useMB=false;
const int Netabins=4;
const bool print_incjetana_tupelequiv=false;

// the macro ------------------------
int printPlots_jetPlots(const std::string input_ppData_condorDir , const std::string input_ppMC_condorDir , 
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
  const std::string input_ppData_Filename_MB="HighPtJetTrig_" +fullJetType+ "-allFiles.root";
  const std::string input_ppData_Filename_noMB="HighPtJetTrig_noMB_" +fullJetType+ "-allFiles.root";
  std::string input_ppData_Filename;
  if(useMB)input_ppData_Filename=input_ppData_Filename_MB;
  else input_ppData_Filename=input_ppData_Filename_noMB;
  
  //const std::string input_ppData_Filename="HighPtJet80_" +fullJetType+ "-allFiles.root";
  //const std::string input_ppData_Filename="HighPtLowerJets_" +fullJetType+ "-allFiles.root";
  const std::string input_ppMC_Filename  ="Py8_CUETP8M1_QCDjetAllPtBins_" +fullJetType+ "-allFiles.root";

  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;  
  const std::string ppData_fullFilename=inputDir+input_ppData_condorDir+input_ppData_Filename;
  const std::string ppMC_fullFilename  =inputDir+input_ppMC_condorDir+input_ppMC_Filename;
  
  // OPEN INPUT SECTION    
  TFile *finData=NULL;  
  std::cout<<" now opening ppData: "<<std::endl<<input_ppData_Filename<<std::endl;
  std::cout<<" in directory: "<<input_ppData_condorDir<<std::endl<<std::endl;
  finData = new TFile(ppData_fullFilename.c_str(), "READ");      
  if(!finData){
    std::cout << " DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppData_condorDir =" << input_ppData_condorDir << std::endl;
    std::cout << "input_ppData_filename  =" << input_ppData_Filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);      }  
  
  TFile *finMC=NULL;
  std::cout<<" now opening ppMC: "<<std::endl<<input_ppMC_Filename<<std::endl;
  std::cout<<" in directory: "<<input_ppMC_condorDir<<std::endl<<std::endl;  
  finMC = new TFile(ppMC_fullFilename.c_str(), "READ");      
  //finMC = TFile::Open(ppMC_fullFilename.c_str(), "READ");      
  //if(finMC->TFile::IsOpen()){
  //if(!((bool)finMC)){
  if(finMC->IsZombie()){
    std::cout << " MC file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppMC_condorDir =" << input_ppMC_condorDir << std::endl;
    std::cout << "input_ppMC_filename  =" << input_ppMC_Filename << std::endl;
    //std::cout << "exiting." << std::endl;
    //assert(false);      
  }
  //assert(false);
  
  
  

  //TH1F* testJetID_fin=(TH1F*)finData->Get("hJetQA_1wJetID_jtpt_etabin0");
  //TH1F* testNOJetID_fin=(TH1F*)finData->Get("hJetQA_0wJetID_jtpt_etabin0");
  //if( !testJetID_fin && testNOJetID_fin)
  //  doJetIDPlots=false;
  //else if(testJetID_fin && !testNOJetID_fin)
  //  doJetIDPlots=true;
  //
  //std::cout<<"doJetIDPlots="<<doJetIDPlots<<std::endl;
  //assert(false);
  
  
  
  
  // GET OUTPUT PDF FILE READY
  std::string theROOTFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetPlots.root";
  TFile* fout=new TFile(theROOTFileName.c_str(),"RECREATE");
  
  std::string thePDFFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetPlots.pdf";
  std::string open_thePDFFileName=thePDFFileName+"[";    std::string close_thePDFFileName=thePDFFileName+"]";  
  
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1000, 800);
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  
  //temp_canvOpen->UseCurrentStyle();
  temp_canvOpen->Close();  
  
  
  // evtcounts/effective integrated luminosity ----------------------  
  long double theLumi;
  if(doEventCounts){
    //old_printDataEventCountReport((TFile*) finData);//this + next 3 lines compatible w/ readForests samples pre 11-27-18, deprecated
    //old_printMCEventCountReport( (TFile*)finMC);
    //old_printDataJetCountReport( (TFile*)finData);
    //old_printMCJetCountReport( (TFile*)finMC);    
    
    theLumi=computeEffLumi( (TFile*) finData);  	    
    //printDataEventCountReport((TFile*) finData);
    //printDataJetCountReport( (TFile*)finData);
    //
    //if((bool)finMC)printMCEventCountReport( (TFile*)finMC);    
    //if((bool)finMC)printMCJetCountReport( (TFile*)finMC);    
    //assert(false);

    //assert(false);    
    

  }
  else {

    printDataEventCountReport((TFile*) finData);
    printDataJetCountReport( (TFile*)finData);
    
    if(!(finMC->IsZombie())){
      std::cout<<"calling printMCEventCountReport"<<std::endl;
      printMCEventCountReport( (TFile*)finMC);    
    }
    if(!(finMC->IsZombie())){
      std::cout<<"calling printMCJetCountReport"<<std::endl;
      printMCJetCountReport( (TFile*)finMC);    
    }
    //std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
    theLumi=intgrtdLumi;
  }
  
  
  
  // evt plots ----------------------
  if(drawEvtQAPlots && !(finMC->IsZombie()) ){
    
    std::cout<<std::endl<<" drawing evt QA Plots now! "<<std::endl<<std::endl;
    
    printMCEvtQAHist( (TFile*) finMC   , "hpthat" ,
		      (std::string) thePDFFileName , (TFile*) fout);
    
    printMCEvtQAHist( (TFile*) finMC   , "hWeightedpthat" ,
		      (std::string) thePDFFileName , (TFile*) fout);
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVz" , 
		       (TFile*) finMC   , "hpthatWeightedVz" ,
		       (int) 10, (std::string) thePDFFileName  , (long double) theLumi , (TFile*) fout ) ;

    printEvtVtxQAHist( (TFile*) finData , "hWeightedVz" , 
		       (TFile*) finMC   , "hWeightedVz" ,
		       (int) 10, (std::string) thePDFFileName  , (long double) theLumi , (TFile*) fout ) ;
    //assert(false);
    
//    printEvtVtxQAHist( (TFile*) finData , "hWeightedVy" , 
//		       (TFile*) finMC   , "hWeightedVy" ,
//		       (int) 20, (std::string) thePDFFileName  , (long double) theLumi , (TFile*) fout  ) ;
//    
//    printEvtVtxQAHist( (TFile*) finData , "hWeightedVx" , 
//		       (TFile*) finMC   , "hWeightedVx" ,
//		       (int) 20, (std::string) thePDFFileName  , (long double) theLumi , (TFile*) fout ) ;
    
//    printEvtVtxQAHist( (TFile*) finData , "hWeightedVr" , 
//		       (TFile*) finMC   , "hWeightedVr" ,
//		       (int) 20, (std::string) thePDFFileName  , (long double) theLumi , (TFile*) fout ) ;
    
    //printEvtNrefQAHist( (TFile*) finData , "hWNref" , 
    //			(TFile*) finMC   , "hWNref" ,
    //			(std::string) thePDFFileName  , (std::string) fullJetType, 
    //			(long double) theLumi  );
    //
    //printEvtNrefQAHist( (TFile*) finData , "hWjetsPEvt" , 
    //			(TFile*) finMC   , "hWjetsPEvt" ,
    //			(std::string) thePDFFileName  , (std::string) fullJetType, 
    //			(long double) theLumi  );
    //
    //printEvtNrefQAHist( (TFile*) finData , "hLeadJetPt" , 
    //			(TFile*) finMC   , "hLeadJetPt" ,
    //			(std::string) thePDFFileName  , (std::string) fullJetType, 
    //			(long double) theLumi  );

    
    
  }  
  else std::cout<<std::endl<<"skipping evt QA plots..."<<std::endl<<std::endl; 
  
  
  
  // jet plots ----------------------  
  if(drawJetQAPlots){
    std::cout<<" drawing jet QA Plots..."<<std::endl;    
    
    std::string jetIDInt;
    if(doJetIDPlots)jetIDInt="1";
    else jetIDInt="0";  
    
    //TH1s
    //int etabin=0;
    for(int j=0;j<N_vars;j++){ 
      for(int etabin=0; etabin<Netabins; etabin++){
	std::cout<<std::endl;
	if(debugMode)std::cout<<std::endl<<" var ="<<var[j]<<", j="<<j<<std::endl;
	
	//constituent/dijet skips
	bool skipConstitPlot= (!drawJetConstituentPlots && j>=jetConstits_varStart && j<dijet_varStart);
	bool skipDijetPlot= (!drawDijetPlots && j>=dijet_varStart);
	
	//skip plot
	bool skipPlot=skipConstitPlot||skipDijetPlot;
	if(skipPlot) {
	  if(debugMode)std::cout<<"skipping jet plot for "<<var[j]<<std::endl;
	continue;}
	
	//printPlot      
	//std::string inHistName="hJetQA_"+jetIDInt+"wJetID_"+var[j];          
	std::string inHistName="hJetQA_"+jetIDInt+"wJetID_"+var[j]+"_etabin"+std::to_string(etabin);          
	printJetQAHist( (TFile*) finData , (TFile*) finMC   ,  (int) j, (bool)doJetIDPlots, (int)etabin, 
			(std::string) inHistName , (std::string) thePDFFileName  , (std::string) fullJetType, 
			(long double) theLumi  , (TFile*) fout );      
      }
    }
    
  }
  


  if(drawTupelJetQAPlots){
    std::cout<<" drawing tupel jet QA Plots..."<<std::endl;    
    
    std::string jetIDInt;
    if(doJetIDPlots)jetIDInt="1";
    else jetIDInt="0";  
    
    //TH1s
    //int etabin=0;
    for(int j=0;j<N_tupelvars;j++){ 
      for(int etabin=0; etabin<Netabins; etabin++){
	std::cout<<std::endl;
	if(debugMode)std::cout<<std::endl<<" tupelvar ="<< tupelvar[j]<<", j="<<j<<std::endl;
	
	//constituent/dijet skips
	//bool skipConstitPlot= (!drawJetConstituentPlots && j>=jetConstits_varStart && j<dijet_varStart);
	//bool skipDijetPlot= (!drawDijetPlots && j>=dijet_varStart);
	
	//skip plot
	//bool skipPlot=skipConstitPlot||skipDijetPlot;
	//if(skipPlot) {
	//  if(debugMode)std::cout<<"skipping jet plot for "<<var[j]<<std::endl;
	//continue;}
	
	//printPlot      
	
	//"hTupelJetQA_"+std::to_string(i)+"jetID_etabin"+std::to_string(j)+"_"+tupelvar[k]).c_str()
	std::string inHistName="hTupelJetQA_"+jetIDInt+"jetID_etabin"+std::to_string(etabin)+"_"+tupelvar[j];          
	if(print_incjetana_tupelequiv){
	  if(tupelvar_incjetanaequiv[j]=="NULL")continue;	  
	  inHistName="hJetQA_"+jetIDInt+"wJetID_"+tupelvar_incjetanaequiv[j]+"_etabin"+std::to_string(etabin);          
	}
	//	printTupelJetQAHist( (TFile*) finData , (TFile*) NULL   ,  (int) j, (bool)doJetIDPlots, (int)etabin, 
	printTupelJetQAHist( (TFile*) finData ,  (int) j, (bool)doJetIDPlots, (int)etabin, 
			     (std::string) inHistName , (std::string) thePDFFileName  , (std::string) fullJetType, 
			     (long double) theLumi  , (TFile*) fout );      
      }
    }
    
  }
  

  
  
  
  
  // jet plots ----------------------  
  if(drawJetRapBinsPlot){
    
    std::cout<<" drawing dual-diff xsec plot(s)..."<<std::endl;
    
    //print dual diff xsec w/ data and MC (RECO)
    //printJetSpectraRapHists( (TFile*) finData , (TFile*) finMC   ,  (bool)doJetIDPlots, 
    //			     (std::string) thePDFFileName  , (std::string) fullJetType, 
    //			     (long double) theLumi  );   
    
    //print dual diff xsec w/ data only (RECO)
    printJetSpectraRapHists( (TFile*) finData , NULL   ,  (bool)doJetIDPlots, 
    			     (std::string) thePDFFileName  , (std::string) fullJetType, 
			     (long double) theLumi  );   
    
    ////print dual diff xsec w/ MC only (RECO v GEN)
    //printMCJetSpectraRapHists( (TFile*) finMC   ,  (bool)doJetIDPlots, 
    //			       (std::string) thePDFFileName  , (std::string) fullJetType);
    //assert(false);
  }
  




  
  
  
  
  
  
  


  
  
  // for DEBUG ONLY
  if(debugMode)std::cout<<std::endl<<"closing output PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1200, 600);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
  
  
  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData->Close();  
  if(debugMode)std::cout<<std::endl<<"closing output ROOT"<<std::endl;
  fout->Close();
  return 0;
  
}// end printplots

				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=4 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetPlots.exe <ppDataCondorDir> <ppMCCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  //if(argc==1)rStatus=printPlots_jetPlots();
  if(argc==4) rStatus=printPlots_jetPlots( (const std::string) argv[1], (const std::string) argv[2],
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
