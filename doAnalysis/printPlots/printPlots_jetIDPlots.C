#include "printPlots.h"


// code/job switches ------------------------
//options
const bool debugMode=true, doEventCounts=true;
//draw switches
const bool drawJetIDPlots=true;
const bool drawJetConstituentPlots=true, drawDijetPlots=false;


//int printPlots_jetIDPlots(const std::string input_condorDir ,
int printPlots_jetIDPlots(const std::string input_condorDir_jetID ,  const std::string input_condorDir_nojetID ,
			  const std::string output_PDFname_base ){
  
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=input_condorDir_jetID.find("_ak")+3;  
  const std::string radiusInt= input_condorDir_jetID.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=input_condorDir_jetID.find("Jets");;
  
  //make final jetType strings
  const std::string jetType=input_condorDir_jetID.substr( jetTypePos,(jetsPos-jetTypePos) );
  const std::string fullJetType="ak"+radiusInt+jetType;
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  //put together input file strings
  const std::string input_ppData_Filename="HighPtJetTrig_" +fullJetType+ "-allFiles.root";
  const std::string input_ppMC_Filename  ="Py8_CUETP8M1_QCDjetAllPtBins_" +fullJetType+ "-allFiles.root";
  
  
  std::size_t isDataString=input_condorDir_jetID.find("ppData");
  std::size_t isMCString=input_condorDir_jetID.find("ppMC");
  bool isData, isMC;
  //if(isDataString==std::string::npos ) isData=false;
  //  else if (isDataString != std::string::npos &&
  //	   (isMCString>isDataString) ) isData=false;
  if(isDataString==std::string::npos ) isData=false;
  else isData=true;
  if(isMCString==std::string::npos ) isMC=false;
  else isMC=true;
  if( (isData && isMC ) || (!isData && !isMC)) assert(false);
  
  
  std::string input_Filename;//, input_Filename_nojetID;
  if(isData) input_Filename=input_ppData_Filename;
  else  input_Filename=input_ppMC_Filename;
  const std::string fullFilename_jetID=inputDir+input_condorDir_jetID+input_Filename;
  const std::string fullFilename_nojetID=inputDir+input_condorDir_nojetID+input_Filename;
  
  
  
  
  // OPEN INPUT SECTION
  TFile *fin_jetID=NULL, *fin_nojetID=NULL;  
  std::cout<<std::endl<<"printing jetID Plots, now opening input files!!"<<std::endl<<std::endl;
  
  //open the input files
  std::cout<<" now opening Jet ID file "<<input_Filename<<std::endl;//fullFilename_jetID<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_condorDir_jetID<<std::endl<<std::endl;
  fin_jetID = new TFile(fullFilename_jetID.c_str());      

  std::cout<<" now opening non-Jet ID file "<<input_Filename<<std::endl;//fullFilename_jetID<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_condorDir_nojetID<<std::endl<<std::endl;
  fin_nojetID = new TFile(fullFilename_nojetID.c_str());      
  
  if(!fin_jetID || !fin_nojetID)assert(false);
  
  
  // GET OUTPUT PDF FILE READY
  std::string theROOTFileName=outputDir+fullJetType;
  theROOTFileName+="_"+output_PDFname_base+"_jetIDPlots.root";
  TFile* fout=new TFile(theROOTFileName.c_str(), "RECREATE");

  std::string thePDFFileName=outputDir+fullJetType;
  thePDFFileName+="_"+output_PDFname_base+"_jetIDPlots.pdf";
  
  std::string open_thePDFFileName, close_thePDFFileName;
  open_thePDFFileName=thePDFFileName+"[";  close_thePDFFileName=thePDFFileName+"]";
  
  // just for opening the pdf
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1000, 800);  
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  
  temp_canvOpen->Close();
  
  
  
  // evtcounts/effective integrated luminosity ----------------------
  long double theLumi=1.;
  if(doEventCounts){
    if(isData)getEventCounts( (TFile*)fin_jetID, isData );
    if(isData)theLumi=computeEffLumi( (TFile*) fin_jetID);   }
  else {
    std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
    if(isData)theLumi=intgrtdLumi;
    else theLumi=1.;
  }
  
  


  // jet ID plots ----------------------
  if(drawJetIDPlots){
    std::cout<<" drawing jet ID Plots..."<<std::endl;

    //std::string jetIDInt;
    //if(doJetIDPlots)jetIDInt="1";
    //else jetIDInt="0";
    //for(int jetIDInt=0;jetIDInt<2;jetIDInt++){
    
    for(int j=0;j<N_vars;j++){

      std::cout<<std::endl;
      if(debugMode)std::cout<<std::endl<<" var ="<<var[j]<<", j="<<j<<std::endl;
      
      //constituent/dijet skips
      bool skipConstitPlot= (!drawJetConstituentPlots && j>=jetConstits_varStart && j<dijet_varStart);
      bool skipDijetPlot= (!drawDijetPlots && j>=dijet_varStart);
      
      //skip plot
      bool skipPlot=skipConstitPlot||skipDijetPlot;
      if(skipPlot) std::cout<<"skipping jet plot for "<<var[j]<<std::endl;
      if(skipPlot) continue;
      
      //printPlot
      //std::string inHistName="hJetQA_"+jetIDInt+"wJetID_"+var[j];
      //      printJetIDHist( (TFile*) fin , (int) j, (bool) isData,
      //		      (std::string) thePDFFileName  , 
      //		      (std::string) fullJetType, (long double) theLumi  );
      printJetIDHist( (TFile*) fin_jetID , (TFile*) fin_nojetID , 
		      (int) j, (bool) isData,
		      (std::string) thePDFFileName  , 
		      (std::string) fullJetType, (long double) theLumi,  
		      (TFile*) fout);
	}

  }  



  // for DEBUG ONLY
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1000, 800);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
  

  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  fin_jetID->Close();  
  fin_nojetID->Close();  
  fout->Close();
  return 0;
}// end printplots
				     
				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=4 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetIDPlots.exe <inputCondorDir_jetID> <inputCondorDir_noJetID> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  if(argc==4) rStatus=printPlots_jetIDPlots( (const std::string) argv[1], (const std::string) argv[2] , (const std::string) argv[3] );
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
   
}
