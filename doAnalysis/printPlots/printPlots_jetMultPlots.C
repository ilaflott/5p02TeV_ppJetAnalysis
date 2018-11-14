#include "printPlots.h"


// code/job switches ------------------------

//options
const bool debugMode=true, doEventCounts=true, doJetIDPlots=true;

//draw switches
const bool drawEvtQAPlots=false;
const bool drawJetQAPlots=false;
const bool drawJetConstituentPlots=false, drawDijetPlots=false;

const bool drawJetMultQAPlots=true;
//const bool drawJetTrackCorrPlots=true;

const bool determineMultClasses=true;  
const int maxNumMultClasses=10+2;
const float multClassFracSize=(1./((float)(maxNumMultClasses-2))); 





//// hist painting ------------------------
//
//// line colors
//const int theDataOverlayLineColor=1, theMCOverlayLineColor=1;//, altOverlayLineColor=3; 
//const int theRatioLineColor=1;//,altRatioLineColor1=8, altRatioLineColor2=7;
//
//// marker colors
//const int theDataOverlayMarkerColor=2, theMCOverlayMarkerColor=4;//,theRatioMarkerColor=9;












// the macro ------------------------
int printPlots_jetMultPlots(const std::string input_ppData_condorDir , const std::string input_ppMC_condorDir , 
			    const std::string output_PDFname_base){
  
  globalHistStyle();
  
  //int printPlots_jetMultPlots(const std::string input_ppData_condorDir , const std::string output_PDFname_base){
  //const std::string input_ppMC_condorDir="";
  
  ////format of the filename is always HighPt{filelist}_ak{3,4}{PF,Calo}-allFiles.root
  //std::size_t jobTypePos=input_ppMC_condorDir.find("__")+2;  
  //std::string jobType;
  //if(jobTypePos!=std::string::npos) jobType=input_ppMC_condorDir.substr( jobTypePos);
  //else                              jobType="noJobType";
  //if(debugMode)std::cout<<"jobType string is = "<<jobType<<std::endl;  
  
  
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
  const std::string input_ppData_Filename="HighPtJetTrig_" +fullJetType+ "-allFiles.root";
  const std::string input_ppMC_Filename  ="Py8_CUETP8M1_QCDjetAllPtBins_" +fullJetType+ "-allFiles.root";
  
  const std::string ppData_fullFilename=inputDir+input_ppData_condorDir+input_ppData_Filename;
  const std::string ppMC_fullFilename  =inputDir+input_ppMC_condorDir+input_ppMC_Filename;
  
  
  
  
  
  
  
  // OPEN INPUT SECTION
  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;  
  
  
  
  TFile *finData=NULL;  
  std::cout<<" now opening ppData: "<<std::endl<<input_ppData_Filename<<std::endl;
  std::cout<<" in directory: "<<input_ppData_condorDir<<std::endl;
  finData = new TFile(ppData_fullFilename.c_str(), "READ");      
  std::cout<<std::endl;
  if(!finData){
    std::cout << " DATA file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppData_condorDir =" << input_ppData_condorDir << std::endl;
    std::cout << "input_ppData_filename  =" << input_ppData_Filename << std::endl;
    std::cout << "exiting." << std::endl;
    assert(false);    
  }  
  
  TFile *finMC=NULL;
  std::cout<<" now opening ppMC: "<<std::endl<<input_ppMC_Filename<<std::endl;
  std::cout<<" in directory: "<<input_ppMC_condorDir<<std::endl<<std::endl;  
  finMC = new TFile(ppMC_fullFilename.c_str(), "READ");      
  std::cout<<std::endl;
  if(!finMC){
    std::cout << " MC file not found" << std::endl;
    std::cout << "inputDir               =" << inputDir << std::endl;
    std::cout << "input_ppMC_condorDir =" << input_ppMC_condorDir << std::endl;
    std::cout << "input_ppMC_filename  =" << input_ppMC_Filename << std::endl;
    //std::cout << "exiting." << std::endl;
    //assert(false);    
  }
  
  
  
  
  

  
  
  // GET OUTPUT PDF FILE READY
  std::string thePDFFileName=outputDir+fullJetType+"_"+output_PDFname_base+"_jetMultPlots.pdf";
  std::string open_thePDFFileName=thePDFFileName+"[";    std::string close_thePDFFileName=thePDFFileName+"]";  
  
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1000, 800);
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  
  //temp_canvOpen->UseCurrentStyle();
  temp_canvOpen->Close();  
  
  
  // evtcounts/effective integrated luminosity ----------------------  
  long double theLumi;
  if(doEventCounts){
    std::cout<<std::endl;
    std::cout<<"doing evt counts"<<std::endl;
    std::cout<<std::endl;
    getEventCountsTrkMult( (TFile*)finData, true );
    //getEventCounts( (TFile*)finMC,   false);
    theLumi=computeEffLumiTrkMult( (TFile*) finData);   }
  else {
    std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
    theLumi=intgrtdLumi;}
  




  if(determineMultClasses){
    std::cout<<std::endl;
    std::cout<<"determining multiplicity classes w/ track counts"<<std::endl;   
    std::cout<<std::endl;
    
    //hJetMultQA_1wJetID_Ntrks
    //hJetMultQA_1wJetID_Ntrks_passCuts
    //hJetMultQA_1wJetID_Nvtxtrks
    //hJetMultQA_1wJetID_Nvtxtrks_vtx0
    
    //std::string inHistName="hJetMultQA_"+ std::to_string( ((int)doJetIDPlots)) + "wJetID_Ntrks";
    //std::string inHistName="hJetMultQA_"+ std::to_string( ((int)doJetIDPlots)) + "wJetID_Ntrks_passCuts";
    std::string inHistName="hJetMultQA_"+ std::to_string( ((int)doJetIDPlots)) + "wJetID_Nvtxtrks";
    //std::string inHistName="hJetMultQA_"+ std::to_string( ((int)doJetIDPlots)) + "wJetID_Nvtxtrks_vtx0";
    
    std::cout<<"using hist " << inHistName << "for determining multplicity classes"<<std::endl;
    TH1F* nTrkHist=(TH1F*)(finData->Get(inHistName.c_str()));
    
    
    int numbins=nTrkHist->GetNbinsX();
    float binWidth=nTrkHist->GetBinWidth(1);
    float xmin=nTrkHist->GetBinLowEdge(1);
    float xmax=nTrkHist->GetBinLowEdge(numbins) + binWidth;
    
    std::cout<< "# of bins = " << numbins  <<std::endl;
    std::cout<< "bin width = " << binWidth << std::endl;
    std::cout<< "xmin = " << xmin << std::endl;
    std::cout<< "xmax = " << xmax << std::endl;
    std::cout<<std::endl;
    float totalIntegral=nTrkHist->Integral();
    std::cout<<"Total Integral = "<< totalIntegral<<std::endl;
    
    
    float multClassIntegral=-1.;
    int nTrks_lo=-1, nTrks_hi = -1;


    ////most peripheral first
    //int multClass=1; //mult class 1 --> most peripheral collision class
    //int multClass_firstBin=1, multClass_secondBin=2;
    //int perc=0, percHi=100, percLo=100;
    
    //most central first // TO DO
    int multClass=maxNumMultClasses; 
    int multClass_firstBin=numbins-1, multClass_secondBin=numbins;
    float perc_F=0., percHi_F=0., percLo_F=0.;
    int perc=0, percHi=0, percLo=0;
    
    
    int multClass_arr[maxNumMultClasses]={0};
    int firstBin_arr[maxNumMultClasses]={0}, secondBin_arr[maxNumMultClasses]={0};
    int nTrksLo_arr[maxNumMultClasses]={0} , nTrksHi_arr[maxNumMultClasses]={0};
    float frac_arr[maxNumMultClasses]={0.};
    int perc_arr[maxNumMultClasses]={0},percHi_arr[maxNumMultClasses]={0},percLo_arr[maxNumMultClasses]={0};
    //    int perc_arr[maxNumMultClasses]={0},percHi_arr[maxNumMultClasses]={100},percLo_arr[maxNumMultClasses]={100};
    
    do{
      
      multClassIntegral=nTrkHist->Integral(multClass_firstBin, multClass_secondBin);
      
      
      if((multClassIntegral/totalIntegral) < multClassFracSize ){
	
	//if(multClass_secondBin%25==0){// most peripheral first
	//  std::cout<<"(multClassIntegral/totalIntegral) < 10%, increasing integral range"<<std::endl;
	//}
	//multClass_secondBin++; 

	if(multClass_firstBin%10==0){ //most central first
	  std::cout<<"(multClassIntegral/totalIntegral) < 10%, increasing integral range"<<std::endl;
	}
	multClass_firstBin--; // most central first // TO DO
      }
      else{
	multClass_firstBin++;//most central first
	//multClass_secondBin--;//most peripheral first
	multClassIntegral=nTrkHist->Integral(multClass_firstBin, multClass_secondBin);
	std::cout<<std::endl;
	std::cout<<"mult class # " << multClass << " done"<<std::endl;
	nTrks_lo=nTrkHist->GetBinLowEdge(multClass_firstBin);
	nTrks_hi=nTrkHist->GetBinLowEdge(multClass_secondBin);// + binWidth;
	std::cout<<"fraction of events w/ # of tracks between bins "<< multClass_firstBin << " and " << multClass_secondBin << " is..."<<std::endl;	
	std::cout<<"fraction of events w/ # of tracks between "<< nTrks_lo << " and " << nTrks_hi << " is..."<<std::endl;	
	std::cout<<"(multClassIntegral/totalIntegral) ="<<(multClassIntegral/totalIntegral)<<std::endl;	
	std::cout<<std::endl;
	
	multClass_arr[multClass-1]=multClass;
	firstBin_arr[multClass-1]=multClass_firstBin;
	secondBin_arr[multClass-1]=multClass_secondBin;
	nTrksLo_arr[multClass-1]=nTrks_lo;
	nTrksHi_arr[multClass-1]=nTrks_hi;
	frac_arr[multClass-1]=(multClassIntegral/totalIntegral);
	
	////most peripheral first
	//perc=((int)((multClassIntegral/totalIntegral)*100.));
	//percHi=percLo;
	//percLo=percLo-perc;	
	//perc_arr[multClass-1]  =perc;
	//percHi_arr[multClass-1]=percHi;
	//percLo_arr[multClass-1]=percLo;	
	//multClass++;
	//multClass_firstBin=multClass_secondBin+1;
	//multClass_secondBin=multClass_firstBin+1;


	//most central first  // TO DO
	perc_F=( ( (multClassIntegral/totalIntegral)*100.) );
	percLo_F=percHi_F;
	percHi_F=percHi_F+perc_F;	

	perc  =((int)(perc_F));
	percLo=((int)(percLo_F));
	percHi=((int)(percHi_F));
	
	perc_arr[multClass-1]  =perc;
	percHi_arr[multClass-1]=percHi;
	percLo_arr[multClass-1]=percLo;	
	
	std::cout<<"multClass = "<<multClass <<std::endl;
	std::cout<<"array ind = "<<(multClass-1) <<std::endl;
	std::cout<<"perc = "   << perc   <<std::endl;
	std::cout<<"percLo = " << percLo <<std::endl;
	std::cout<<"percHi = " << percHi <<std::endl;
	std::cout<<"perc_arr[multClass-1] = "   << perc_arr[multClass-1] <<std::endl;
	std::cout<<"percLo_arr[multClass-1] = " << percLo_arr[multClass-1] <<std::endl;
	std::cout<<"percHi_arr[multClass-1] = " << percHi_arr[multClass-1] <<std::endl;


	multClass--;
	multClass_secondBin=multClass_firstBin-1;
	multClass_firstBin=multClass_secondBin-1;
	

      }
      
      
      
      
      
    }
    while(multClass >= 1 && multClass_firstBin > 1 ) ; //fill peripheral first, then central  // TO DO
    //while(multClass <= maxNumMultClasses && multClass_secondBin < numbins) ; //fill peripheral first, then central

    
    
    std::cout<<std::endl;    
    std::cout<<"mult class # " << multClass << " done"<<std::endl;
    nTrks_lo=nTrkHist->GetBinLowEdge(multClass_firstBin);
    nTrks_hi=nTrkHist->GetBinLowEdge(multClass_secondBin) + binWidth;
    std::cout<<"fraction of events w/ # of tracks between bins "<< multClass_firstBin << " and " << multClass_secondBin << " is..."<<std::endl;	
    //std::cout<<"fraction of events w/ # of tracks between "<< nTrks_lo << " and " << nTrks_hi << " is..."<<std::endl;	
    std::cout<<"(multClassIntegral/totalIntegral) ="<<(multClassIntegral/totalIntegral)<<std::endl;	
    multClass_arr[multClass-1]=multClass;
    firstBin_arr[multClass-1]=multClass_firstBin;
    secondBin_arr[multClass-1]=multClass_secondBin;
    nTrksLo_arr[multClass-1]=nTrks_lo;
    nTrksHi_arr[multClass-1]=nTrks_hi;
    
//    //most peripheral first
//    frac_arr[multClass-1]=(multClassIntegral/totalIntegral);
//    perc=((int)((multClassIntegral/totalIntegral)*100.));
//    percHi=percLo;
//    percLo=percLo-perc;
//    perc_arr[multClass-1]  =perc;
//    percHi_arr[multClass-1]=percHi;
//    percLo_arr[multClass-1]=0;//percLo;


    //most central first
    frac_arr[multClass-1]=(multClassIntegral/totalIntegral);    
    perc_F=( ( (multClassIntegral/totalIntegral)*100.) );
    percLo_F=percHi_F;
    percHi_F=percHi_F+perc_F;	
    perc=((int)(perc_F));
    percLo=((int)(percHi_F));
    percHi=((int)(percHi_F+perc_F));	
    
    perc_arr[multClass-1]  =perc;
    percHi_arr[multClass-1]=percHi;//100;
    percLo_arr[multClass-1]=percLo;

    std::cout<<std::endl;    
    std::cout<<"done determining multiplicity classes. printing out summary table.."<<std::endl;
    std::cout<<std::endl;    
    
    
    //TPERC - total percent of events between NTRKL and NTRKH
    //LPERC/HPERC - low/high percentile value of mult class
    //NTRKL/NTRKH - low/high NTRK limits for specific mult class
    std::cout <<"FORMAT     ::MULTI::ARAYI::TPERC::LPERC::HPERC::NTRKL::NTRKH" << std::endl;    
    //for(int i =0; i<maxNumMultClasses; ++i){
    for(int i =(maxNumMultClasses-1); i>=0; --i){
      

      if(multClass_arr[i]==maxNumMultClasses)
	std::cout<<"CENTRAL    :: ";
      else if(multClass_arr[i]==1)
	std::cout<<"PERIPHERAL :: ";
      else
	std::cout<<"           :: ";

      std::cout<< multClass_arr[i] ;
      if(multClass_arr[i] < 10) std::cout<< "   :: " ;
      else if(multClass_arr[i] < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;
      
      std::cout<< frac_arr[i] ;
      //if(frac_arr[i] < 10.) std::cout<< "   :: " ;
      //else if(frac_arr[i] < 100.) std::cout<< "  :: " ;
      //else std::cout<< " :: " ;
      std::cout<< " :: " ;
      
      std::cout<< i ;
      if(i < 10) std::cout<< "   :: " ;
      else if(i < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;

      std::cout<< perc_arr[i] ;
      if(perc_arr[i] < 10) std::cout<< "   :: " ;
      else if(perc_arr[i] < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;

      std::cout<< percLo_arr[i] ;
      if(percLo_arr[i] < 10) std::cout<< "   :: " ;
      else if(percLo_arr[i] < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;

      std::cout<< percHi_arr[i] ;
      if(percHi_arr[i] < 10) std::cout<< "   :: " ;
      else if(percHi_arr[i] < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;

      std::cout<< firstBin_arr[i]  ;
      if(firstBin_arr[i] < 10) std::cout<< "   :: " ;
      else if(firstBin_arr[i] < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;

      std::cout<< secondBin_arr[i] ;
      if(secondBin_arr[i] < 10) std::cout<< "   :: " ;
      else if(secondBin_arr[i] < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;
      
      std::cout<< nTrksLo_arr[i]   ;
      if(nTrksLo_arr[i] < 10) std::cout<< "   :: " ;
      else if(nTrksLo_arr[i] < 100) std::cout<< "  :: " ;
      else std::cout<< " :: " ;
      
      std::cout<< nTrksHi_arr[i]   ;
      std::cout<<std::endl;    
      
    }
    std::cout << "done getting results using the histogram " << inHistName << std::endl;

  }
  
  
  // evt plots ----------------------
  if(drawEvtQAPlots){
    
    std::cout<<std::endl<<" drawing evt QA Plots now! "<<std::endl<<std::endl;
    
    printMCEvtQAHist( (TFile*) finMC   , "hpthat" ,
		      (std::string) thePDFFileName );
    
    printMCEvtQAHist( (TFile*) finMC   , "hWeightedpthat" ,
		      (std::string) thePDFFileName );
    

    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVz" , 
		       (TFile*) finMC   , "hWeightedVz" ,
		       (int) 10, (std::string) thePDFFileName  , (long double) theLumi  ) ;
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVy" , 
		       (TFile*) finMC   , "hWeightedVy" ,
		       (int) 5, (std::string) thePDFFileName  , (long double) theLumi  ) ;
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVx" , 
		       (TFile*) finMC   , "hWeightedVx" ,
		       (int) 5, (std::string) thePDFFileName  , (long double) theLumi  ) ;
    
    printEvtVtxQAHist( (TFile*) finData , "hWeightedVr" , 
		       (TFile*) finMC   , "hWeightedVr" ,
		       (int) 5, (std::string) thePDFFileName  , (long double) theLumi  ) ;
    



    
    printEvtNrefQAHist( (TFile*) finData , "hWNref" , 
			(TFile*) finMC   , "hWNref" ,
			(std::string) thePDFFileName  , (std::string) fullJetType, 
			(long double) theLumi  );
    
    printEvtNrefQAHist( (TFile*) finData , "hWjetsPEvt" , 
			(TFile*) finMC   , "hWjetsPEvt" ,
			(std::string) thePDFFileName  , (std::string) fullJetType, 
			(long double) theLumi  );

    printEvtNrefQAHist( (TFile*) finData , "hLeadJetPt" , 
			(TFile*) finMC   , "hLeadJetPt" ,
			(std::string) thePDFFileName  , (std::string) fullJetType, 
			(long double) theLumi  );

    
    
  }  
  else std::cout<<std::endl<<"skipping evt QA plots..."<<std::endl<<std::endl; 
  
    // evt plots ----------------------
  if(drawJetMultQAPlots){
    
    std::cout<<std::endl<<" drawing jetMult QA Plots now! "<<std::endl<<std::endl;
    std::string jetMultHist_str="hJetMultQA_1wJetID_";

    //VTX TH1s
    std::cout<<"begin printJetMultQATH1_vtx"<<std::endl;
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_x"          ).c_str() , "x_{vtx}    (cm)" ,(int) 20,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_y"          ).c_str() , "y_{vtx}    (cm)" ,(int) 20,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_z"          ).c_str() , "z_{vtx}    (cm)" ,(int) 20,(std::string) thePDFFileName, (long double) theLumi);
    
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_r"          ).c_str() , "r_{vtx}    (cm)" ,(int) 20,(std::string) thePDFFileName, (long double) theLumi); //cylindrical rad

    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_rad"          ).c_str() , "R_{vtx} (cm)" ,      (int) 1,(std::string) thePDFFileName, (long double) theLumi);//spherical rad
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_phi"          ).c_str() , "#phi_{vtx} (rad)" ,      (int) 1,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_theta"          ).c_str() , "#theta_{vtx} (rad)",     (int) 1,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_tantheta2"          ).c_str() , "tan #theta/2" ,(int) 1,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_lntantheta2"          ).c_str() , "ln tan #theta/2" ,    (int) 1,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_eta"          ).c_str() , "#eta" ,    (int) 1,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_yOverx"          ).c_str() , "y/x" ,    (int) 1,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_rOverz"          ).c_str() , "r/z" ,    (int) 1,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_zOverrad"          ).c_str() , "z/R" ,    (int) 1,(std::string) thePDFFileName, (long double) theLumi);
    
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Nvtx"           ).c_str() , "N_{evt}^{vtx} "       ,(int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Ntrks"          ).c_str() , "N_{evt}^{trks}"      ,(int) 3 ,(std::string) thePDFFileName, (long double) theLumi);
    //printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Ntrks_vtxs"       ).c_str() , "N_{evt}^{trks} assoc. w/ vtx"      ,(int) 2 ,(std::string) thePDFFileName, (long double) theLumi);        
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Ntrks_passCuts"          ).c_str() , "N_{evt}^{trks} passing cuts"      ,(int) 3 ,(std::string) thePDFFileName, (long double) theLumi);    
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Nvtxtrks"       ).c_str() , "N_{vtx}^{trks}"      ,(int) 2 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Nvtxtrks_vtx0"  ).c_str() , "N_{vtx0}^{trks}"     ,(int) 2 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Nvtxtrks_vtxgt0").c_str() , "N_{vtxGT0}^{trks}"   ,(int) 2 ,(std::string) thePDFFileName, (long double) theLumi);  

    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_maxPt"           ).c_str() , "max p_{T} vtx "       ,(int) 1 ,(std::string) thePDFFileName, (long double) theLumi);     
    printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"vtx_maxMult"           ).c_str() , "max mult vtx "       ,(int) 1 ,(std::string) thePDFFileName, (long double) theLumi);    

    //printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Ntrks_vtx0"       ).c_str() , "N_{vtx0}^{trks}"      ,(int) 2 ,(std::string) thePDFFileName, (long double) theLumi);
    //printJetMultQATH1_vtx( (TFile*) finData, (TFile*) finMC,  (jetMultHist_str+"Ntrks_passCuts_vtx0"       ).c_str() , "N_{vtx0}^{trks} passing cuts"      ,(int) 2 ,(std::string) thePDFFileName, (long double) theLumi);
    
    
    

    
    //TRK TH1s
    std::cout<<"begin printJetMultQATH1_trk"<<std::endl;
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkPt" ).c_str() , "p_{T}^{trk} (GeV)", (int) 4 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkPtError" ).c_str() , "#sigma(p_{T}^{trk}) (GeV)", (int) 4 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkPtPercError" ).c_str() , "% #sigma(p_{T}^{trk})/p_{T}^{trk}", (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);

    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkEta").c_str() , "#eta^{trk}" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkPhi").c_str() , "#phi^{trk} (rad)" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);

    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkDz").c_str() , "D_{z}^{trk} (cm)" , (int) 5 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkDzError").c_str() , "#sigma(D_{z}^{trk}) (cm)" , (int) 5 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkDzSig").c_str() , "(D_{z}^{trk})/#sigma(D_{z}^{trk})" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);

    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkDxy").c_str() , "D_{xy}^{trk} (cm)" , (int) 5 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkDxyError").c_str() , "#sigma(D_{xy}^{trk}) (cm)" , (int) 5 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkDxySig").c_str() , "(D_{xy}^{trk})/#sigma(D_{xy}^{trk})" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);

    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkNVtx").c_str() , "N_{vtx} assoc. w/ each trk" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkChi2").c_str() , "#chi^{2} of trk" , (int) 5 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkNdof").c_str() , "N_{d.o.f.}^{trk}" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkNHit").c_str() , "N_{hits}^{trk}" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkNLayer").c_str() , "N_{layers}^{trk}" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);

    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkCharge").c_str() , "Charge of Track" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkHighPur").c_str() , "High Purity Track Flag" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkTight").c_str() , "Tight Track Flag" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkLoose").c_str() , "Loose Track Flag" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    printJetMultQATH1_trk( (TFile*) finData, (TFile*) finMC, (jetMultHist_str+"trkFake").c_str() , "Fake Track Flag" , (int) 1 ,(std::string) thePDFFileName, (long double) theLumi);
    
    

    

    ////DATA VTX TH2s
    std::cout<<"begin printJetMultQATH2"<<std::endl;
    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"vtx_rVz"          ).c_str(), "r_{vtx} (cm)", (int) 5,  "z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"vtx_xVy"          ).c_str(), "x_{vtx} (cm)", (int) 5,  "y_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );

    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"vtx_maxPt_v_maxMult"          ).c_str(), "max p_{T} vtx", (int) 1,  "max mult vtx", (int) 1, (std::string) thePDFFileName );

//    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"vtx_maxPt_v_maxMult"          ).c_str(), "max p_{T} vtx", (int) 1,  "max mult vtx", (int) 1, (std::string) thePDFFileName );
//    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"vtx_maxPt_v_maxMult"          ).c_str(), "max p_{T} vtx", (int) 1,  "max mult vtx", (int) 1, (std::string) thePDFFileName );
//    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"vtx_maxPt_v_maxMult"          ).c_str(), "max p_{T} vtx", (int) 1,  "max mult vtx", (int) 1, (std::string) thePDFFileName );

    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"NvtxtrksVz"       ).c_str(), "N_{vtx}^{trks}"  , (int) 1,  "z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"NvtxtrksVz_vtx0"  ).c_str(), "N_{vtx0}^{trks}"  , (int) 1,  "z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"NvtxtrksVz_vtxgt0").c_str(), "N_{vtxGT0}^{trks}"  , (int) 1,  "z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"NvtxtrksVr"       ).c_str(), "N_{vtx}^{trks}"  , (int) 1,  "r_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"NvtxtrksVr_vtx0"  ).c_str(), "N_{vtx0}^{trks}"  , (int) 1,  "r_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finData , (jetMultHist_str+"NvtxtrksVr_vtxgt0").c_str(), "N_{vtxGT0}^{trks}"  , (int) 1,  "r_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );

    ////DATA TRK TH2s
    printJetMultQATH2( (TFile*) finData   , (jetMultHist_str+"trkPtvEta" ).c_str() , "p_{T}^{trk}" , (int) 2, "#eta^{trk}" , (int) 1,  (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finData   , (jetMultHist_str+"trkEtaVPhi").c_str() , "#eta^{trk}" , (int) 1, "#phi^{trk}" , (int) 1,   (std::string) thePDFFileName );


    ////MC VTX TH2s
    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"vtx_rVz"          ).c_str(), "MC r_{vtx} (cm)"      , (int) 5,  "MC z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"vtx_xVy"          ).c_str(), "MC x_{vtx} (cm)"      , (int) 5,  "MC y_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );

    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"vtx_maxPt_v_maxMult").c_str(), "MC max p_{T} vtx", (int) 1,  "MC max mult vtx", (int) 1, (std::string) thePDFFileName );

    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"NvtxtrksVz"       ).c_str(), "MC N_{vtx}^{trks}"    , (int) 1,  "MC z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"NvtxtrksVz_vtx0"  ).c_str(), "MC N_{vtx0}^{trks}"   , (int) 1,  "MC z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"NvtxtrksVz_vtxgt0").c_str(), "MC N_{vtxGT0}^{trks}" , (int) 1,  "MC z_{vtx} (cm)", (int) 20, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"NvtxtrksVr"       ).c_str(), "MC N_{vtx}^{trks}"    , (int) 1,  "MC r_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"NvtxtrksVr_vtx0"  ).c_str(), "MC N_{vtx0}^{trks}"   , (int) 1,  "MC r_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finMC , (jetMultHist_str+"NvtxtrksVr_vtxgt0").c_str(), "MC N_{vtxGT0}^{trks}" , (int) 1,  "MC r_{vtx} (cm)", (int) 5, (std::string) thePDFFileName );
    
    ////MC TRK TH2s
    printJetMultQATH2( (TFile*) finMC   , (jetMultHist_str+"trkPtvEta" ).c_str() , "MC p_{T}^{trk}" , (int) 2, "MC #eta^{trk}" , (int) 1,  (std::string) thePDFFileName );
    printJetMultQATH2( (TFile*) finMC   , (jetMultHist_str+"trkEtaVPhi").c_str() , "MC #eta^{trk}"  , (int) 1, "MC #phi^{trk}" , (int) 1,   (std::string) thePDFFileName );
        
    ////MIX TH2s
    //printJetMultQATH2_mix( (TFile*) finData   , "Nvtxs_v_Nref" , (std::string) thePDFFileName );
    //printJetMultQATH2_mix( (TFile*) finData   , "Ntrks_v_Nref" , (std::string) thePDFFileName );
    //printJetMultQATH2_mix( (TFile*) finData   , "Ntrks_v_Nvtxs" , (std::string) thePDFFileName );
    
  }
  
  // jet plots ----------------------  
  if(drawJetQAPlots){
    std::cout<<" drawing jet QA Plots..."<<std::endl;    
    
    std::string jetIDInt;
    if(doJetIDPlots)jetIDInt="1";
    else jetIDInt="0";  

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
      std::string inHistName="hJetQA_"+jetIDInt+"wJetID_"+var[j];          
      printJetQAHist( (TFile*) finData , (TFile*) finMC   ,  (int) j, (bool)doJetIDPlots, 
		      (std::string) inHistName , (std::string) thePDFFileName  , (std::string) fullJetType, 
		      (long double) theLumi  );      
    }
    
  }
  
  
  
  // for DEBUG ONLY
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1200, 600);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
  
  
  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData->Close();  
  
  return 0;
  
}// end printplots

				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=4 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetMultPlots.exe <ppDataCondorDir> <ppMCCondorDir> <outputNameBase>"<<std::endl;
    //std::cout<<"do ./printPlots_jetMultPlots.exe <ppDataCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  //if(argc==1)rStatus=printPlots_jetMultPlots();
  if(argc==4) rStatus=printPlots_jetMultPlots( (const std::string) argv[1], (const std::string) argv[2], 
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
