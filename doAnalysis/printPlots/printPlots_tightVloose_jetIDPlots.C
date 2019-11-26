#include "printPlots.h"


// code/job switches ------------------------
//options
const bool debugMode=true, doEventCounts=true;
//draw switches
const bool drawJetIDPlots=true;
const bool drawJetConstituentPlots=true, drawDijetPlots=false;

const bool drawJetQAPlots=false;
const bool drawTupelJetQAPlots=true&&!drawJetQAPlots;
//bool useMB=false;
const int Netabins=4;
const bool print_incjetana_tupelequiv=false;

//int printPlots_tightVloose_jetIDPlots(const std::string input_condorDir ,
int printPlots_tightVloose_jetIDPlots(const std::string input_condorDir_loose_jetID , 
				 const std::string input_condorDir_tight_jetID , 
				 const std::string output_PDFname_base ){
  
  
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=input_condorDir_loose_jetID.find("_ak")+3;  
  const std::string radiusInt= input_condorDir_loose_jetID.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=input_condorDir_loose_jetID.find("Jets");;
  
  //make final jetType strings
  const std::string jetType=input_condorDir_loose_jetID.substr( jetTypePos,(jetsPos-jetTypePos) );
  const std::string fullJetType="ak"+radiusInt+jetType;
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  //bool useMB=false;
  std::string input_ppData_Filename="HighPtJetTrig_noMB_" +fullJetType+ "-allFiles.root";
  std::cout<<"input_ppData_Filename="<<input_ppData_Filename<<std::endl;
  //assert(false);
  
  
  bool isData=true;
  //bool isMC=false;
  
  std::string input_Filename;//, input_Filename_nojetID;
  if(isData) input_Filename=input_ppData_Filename;
  //else if (isMC)input_Filename=input_ppMC_Filename;
  else assert(false);
  //else  input_Filename=input_ppMC_Filename;
  const std::string fullFilename_loose_jetID=inputDir+input_condorDir_loose_jetID+input_Filename;
  const std::string fullFilename_tight_jetID=inputDir+input_condorDir_tight_jetID+input_Filename;
  
  // OPEN INPUT SECTION
  TFile *fin_loose_jetID=NULL;
  TFile *fin_tight_jetID=NULL;
  std::cout<<std::endl<<"printing jetID Plots, now opening input files!!"<<std::endl<<std::endl;
  
  //open the input files
  std::cout<<" now opening HIN Jet ID file "<<input_Filename<<std::endl;//fullFilename_jetID<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_condorDir_loose_jetID<<std::endl<<std::endl;
  fin_loose_jetID = new TFile(fullFilename_loose_jetID.c_str());      
  
  std::cout<<" now opening SMP Jet ID file "<<input_Filename<<std::endl;//fullFilename_jetID<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_condorDir_tight_jetID<<std::endl<<std::endl;
  fin_tight_jetID = new TFile(fullFilename_tight_jetID.c_str());      
  
  //  std::cout<<" now opening non-Jet ID file "<<input_Filename<<std::endl;//fullFilename_jetID<<std::endl;
  //  std::cout<<" in directory: "<<inputDir+input_condorDir_nojetID<<std::endl<<std::endl;
  //  fin_nojetID = new TFile(fullFilename_nojetID.c_str());      
  
  //  if(!fin_jetID || !fin_nojetID)assert(false);
  
  
  // GET OUTPUT PDF FILE READY
  std::string theROOTFileName=outputDir+fullJetType;
  theROOTFileName+="_"+output_PDFname_base+"_tightVloose_jetIDPlots.root";
  TFile* fout=new TFile(theROOTFileName.c_str(), "RECREATE");
  
  std::string thePDFFileName=outputDir+fullJetType;
  thePDFFileName+="_"+output_PDFname_base+"_tightVloose_jetIDPlots.pdf";
  
  std::string open_thePDFFileName, close_thePDFFileName;
  open_thePDFFileName=thePDFFileName+"[";  close_thePDFFileName=thePDFFileName+"]";
  
  // just for opening the pdf
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1000, 800);  
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  
  temp_canvOpen->Close();
  
  
  
  // evtcounts/effective integrated luminosity ----------------------
  //long double theLumi=1.;
  //if(doEventCounts){
  //  if(isData)getEventCounts( (TFile*)fin_jetID, isData );
  //  if(isData)theLumi=computeEffLumi( (TFile*) fin_jetID);   }
  //else {
  //  std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
  //  if(isData)theLumi=intgrtdLumi;
  //  else theLumi=1.;
  //}
  //theLumi=1.;
  
  


  std::cout<<" drawing tightVloose jetID Plots..."<<std::endl;

  
  for( int etabin=0; etabin<Netabins ; etabin++ ){
    std::cout<<std::endl;
    print_tightVlooseJetIDHist_pt( (TFile*) fin_loose_jetID , (TFile*) fin_tight_jetID , 
				   (int)etabin,  (std::string)"jtpt", 
				   (bool) isData,
				   (std::string) thePDFFileName  , 
				   (std::string) fullJetType, (TFile*) fout);  }//etabin loop
  
  
  
  
  // for DEBUG ONLY
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1000, 800);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
  
  
  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  fin_loose_jetID->Close();  
  fin_tight_jetID->Close();  
  fout->Close();
  return 0;
}// end printplots

				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=4 ) {//no input arguments, error
    std::cout<<"do ./printPlots_tightVloose_jetIDPlots.exe <inputCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  if(argc==4) rStatus=printPlots_tightVloose_jetIDPlots( (const std::string) argv[1], (const std::string) argv[2] , (const std::string) argv[3] );
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
   
}
