#include "printQAPlots.h"


// env ------------------------

//CMSSW_BASE
const std::string CMSSW_BASE= 
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";

const std::string readForests_ppData_dir=CMSSW_BASE+
  "src/readFiles/readForests_ppData/saved_outputCondor/";

const std::string readForests_ppMC_dir=CMSSW_BASE+
  "src/readFiles/readForests_ppMC/saved_outputCondor/";

//output
const std::string printQAPlots_dir=CMSSW_BASE+"src/printQAPlots/output/";

// code/job switches ------------------------
const bool doEventCounts=true, drawEvtQAPlots=true;
const bool drawJetQAPlots=true, drawJetTrigQAPlots=false;
const bool drawMCEffPlots=false;//, drawJECandJERPlots=true;//MC Specific for now

const bool debugMode=true;//, dataDebugMode=true, MCDebugMode=true; //debug



// the macro ------------------------
int printQAPlots(const std::string input_ppData_condorDir, const std::string input_ppData_Filename,
		 const std::string input_ppMC_condorDir, const std::string input_ppMC_Filename, 
		 const std::string output_PDFname_base_full, const bool drawDataMCOverlaysInput=true ){
 
  //plot style
  const bool drawDataMCOverlays=drawDataMCOverlaysInput, drawDataMCRatios =!(drawDataMCOverlays);
  const bool drawTrigDataOverlays=drawDataMCOverlays,    drawTrigDataRatios=!(drawTrigDataOverlays);

  //for efficiency later if/when I get a chance 
  const bool openDataFile=true;
  const bool openMCFile  =true;


  //figure out what radius/jetcollection we are looking at using the ppData filename
  //format of the filename is always HighPt{filelist}_ak{3,4}{PF,Calo}-allFiles.root
  std::size_t radPos=input_ppData_Filename.find("_ak")+3;  
  const std::string radiusIntOnly=input_ppData_Filename.substr( radPos,1 );

  const std::string radius="R"+
    input_ppData_Filename.substr( radPos,1 )+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;

  std::size_t jetTypePos=radPos+1;
  std::size_t hyphenPos=input_ppData_Filename.find("-allFiles");;
  const std::string jetType=input_ppData_Filename.substr( jetTypePos,(hyphenPos-jetTypePos) );
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;

  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;

  //input file(s)
  TFile *finMC=NULL, *finData=NULL;
  TFile* finData_LowerJets=NULL;
  assert( openMCFile||openDataFile );

  

  const std::string ppData_fullFilename=readForests_ppData_dir+
    input_ppData_condorDir+input_ppData_Filename;
  const std::string ppData_fullFilename_LowerJets=readForests_ppData_dir+
    input_ppData_condorDir+"HighPtLowerJets_ak"+radiusIntOnly+jetType+"-allFiles.root";
  const std::string ppMC_fullFilename  =readForests_ppMC_dir+ 
    input_ppMC_condorDir+input_ppMC_Filename;;  
  
  if(openDataFile){
    if(debugMode)std::cout<<" now opening Data File "<<std::endl<<ppData_fullFilename<<std::endl<<std::endl;
    else std::cout<<" now opening Data File "<<std::endl<<input_ppData_Filename<<std::endl<<std::endl;
    finData = new TFile(ppData_fullFilename.c_str());      
    finData_LowerJets = new TFile(ppData_fullFilename_LowerJets.c_str());        }
  
  if(openMCFile){  
    if(debugMode) std::cout<<" now opening MC File "<<std::endl<<ppMC_fullFilename<<std::endl<<std::endl;
    else std::cout<<" now opening MC File "<<std::endl<<input_ppMC_Filename<<std::endl<<std::endl;
    finMC = new TFile(ppMC_fullFilename.c_str());}
  assert(finData||finMC);
  
  // Lumi Efficiency Calculation ------------------------
  long double effIntgrtdLumi_vz=1.   , LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  long double effIntgrtdLumi_kmat=1. , LumiEff_kmat=1.;// eff-> effective, Eff->efficiency
  long double effIntgrtdLumi_JetID=1., LumiEff_JetID=1.;// eff-> effective, Eff->efficiency

  long double alt_effIntgrtdLumi_vz=1.   ;
  long double alt_effIntgrtdLumi_kmat=1. ;
  long double alt_effIntgrtdLumi_JetID=1.;

  // evt counts----------------------
  if(!doEventCounts)std::cout<<"skipping evt QA plots..."<<std::endl<<std::endl;
  else{ std::cout<<" grabbing raw event counts "<<std::endl;

    const bool doDataCounts=true, doMCCounts=true;
    if(doDataCounts){ std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;

      TH1F *h_NEvents           = (TH1F*)finData->Get("NEvents");
      TH1F *h_NEvents_LowerJets = (TH1F*)finData_LowerJets->Get("NEvents");

      TH1F *h_NEvents_read      = (TH1F*)finData->Get("NEvents_read");

      std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;


      TH1F *h_NEvents_skimCut = (TH1F*)finData->Get("NEvents_skimCut");

      TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_vzCut");
      TH1F *h_NEvents_vzCut_LowerJets   = (TH1F*)finData_LowerJets->Get("NEvents_vzCut");
      
      std::cout<<"Total Num of Events read passing skimCuts                              = " <<
	h_NEvents_skimCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read passing vzCuts and skimCuts                   = " <<
	h_NEvents_vzCut->GetEntries()<<std::endl;
      std::cout<<std::endl;


      //TH1F *h_NEvents_trgDup   = (TH1F*)finData->Get("NEvents_trgDup");   
      //TH1F *h_NEvents_trgDupCut= (TH1F*)finData->Get("NEvents_trgDupCut");
      //
      //std::cout<<"Total Num of duplicate Events read post skim+vz cuts                 = " <<
      //	h_NEvents_trgDup->GetEntries()<<std::endl;
      //std::cout<<"Total Num of Events read passing skim+vz+duplicate (aka good events) = " <<
      //	h_NEvents_trgDupCut->GetEntries()<<std::endl;
      //std::cout<<std::endl;


      TH1F *h_NEvents_withJets           = (TH1F*)finData->Get("NEvents_withJets");

      TH1F *h_NEvents_withJets_kmatCut   = (TH1F*)finData->Get("NEvents_withJets_kmatCut");
      TH1F *h_NEvents_withJets_kmatCut_LowerJets   = (TH1F*)finData_LowerJets->Get("NEvents_withJets_kmatCut");

      TH1F *h_NEvents_withJets_JetIDCut  = (TH1F*)finData->Get("NEvents_withJets_JetIDCut");
      TH1F *h_NEvents_withJets_JetIDCut_LowerJets  = (TH1F*)finData_LowerJets->Get("NEvents_withJets_JetIDCut");
      
      std::cout<<"Total Num of good Events w/ jets                             = " <<
	h_NEvents_withJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
	h_NEvents_withJets_kmatCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
	h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;      
      
      LumiEff_vz = (h_NEvents_vzCut->GetEntries()+h_NEvents_vzCut_LowerJets->GetEntries())/
	(h_NEvents->GetEntries()+h_NEvents_LowerJets->GetEntries());
      effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;
      LumiEff_kmat = (h_NEvents_withJets_kmatCut->GetEntries()+h_NEvents_withJets_kmatCut_LowerJets->GetEntries())/
	(h_NEvents->GetEntries()+h_NEvents_LowerJets->GetEntries());
      effIntgrtdLumi_kmat=intgrtdLumi*LumiEff_kmat;
      LumiEff_JetID = (h_NEvents_withJets_JetIDCut->GetEntries()+h_NEvents_withJets_JetIDCut_LowerJets->GetEntries())/
	(h_NEvents->GetEntries()+h_NEvents_LowerJets->GetEntries());
      effIntgrtdLumi_JetID=intgrtdLumi*LumiEff_JetID;
      
      alt_effIntgrtdLumi_vz    =0.99*intgrtdLumi*LumiEff_vz;
      alt_effIntgrtdLumi_kmat  =0.99*intgrtdLumi*LumiEff_kmat;
      alt_effIntgrtdLumi_JetID =0.99*intgrtdLumi*LumiEff_JetID;

      std::cout<<std::endl<<"dataset integrated Luminosity (millibarns) ="<<intgrtdLumi<<std::endl;
      std::cout<<std::endl<<"lumi efficiency, vz cuts="<<LumiEff_vz<<std::endl;
      std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;
      std::cout<<std::endl<<"lumi efficiency, kmatCuts="<<LumiEff_kmat<<std::endl;
      std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_kmat<<std::endl;
      std::cout<<std::endl<<"lumi efficiency, JetID="<<LumiEff_JetID<<std::endl;
      std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_JetID<<std::endl;


      TH1F *h_NJets          = (TH1F*)finData->Get("NJets");
      TH1F *h_NJets_kmatCut  = (TH1F*)finData->Get("NJets_kmatCut");
      TH1F *h_NJets_JetIDCut = (TH1F*)finData->Get("NJets_JetIDCut");
      
      std::cout<<std::endl<<"/// Data Jet-Count Summary Summary ///"<<std::endl<<std::endl;
      std::cout<<"Total Num of Jets read from good events                          = " <<
	h_NJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut              = " <<
	h_NJets_kmatCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
	h_NJets_JetIDCut->GetEntries()<<std::endl;
      std::cout<<std::endl;    }

    
    //MC
    if(doMCCounts){ std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;
      TH1F *h_NEvents         = (TH1F*)finMC->Get("NEvents");
      TH1F *h_NEvents_read    = (TH1F*)finMC->Get("NEvents_read");

      std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;


      TH1F *h_NEvents_skimCut = (TH1F*)finMC->Get("NEvents_skimCut");
      TH1F *h_NEvents_vzCut   = (TH1F*)finMC->Get("NEvents_vzCut");
      
      std::cout<<"Total Num of Events read passing skimCuts                              = " <<
	h_NEvents_skimCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
	h_NEvents_vzCut->GetEntries()<<std::endl;
      std::cout<<std::endl;


      TH1F *h_NEvents_withJets           = (TH1F*)finMC->Get("NEvents_withJets");
      TH1F *h_NEvents_withJets_kmatCut   = (TH1F*)finMC->Get("NEvents_withJets_kmatCut");
      TH1F *h_NEvents_withJets_JetIDCut = (TH1F*)finMC->Get("NEvents_withJets_JetIDCut");
      
      std::cout<<"Total Num of good Events w/ jets                             = " <<
	h_NEvents_withJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
	h_NEvents_withJets_kmatCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
	h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;      


      TH1F *h_NJets          = (TH1F*)finMC->Get("NJets");
      TH1F *h_NJets_kmatCut  = (TH1F*)finMC->Get("NJets_kmatCut");
      TH1F *h_NJets_JetIDCut = (TH1F*)finMC->Get("NJets_JetIDCut");
      
      std::cout<<std::endl<<"/// Data Jet-Count Summary Summary ///"<<std::endl<<std::endl;
      std::cout<<"Total Num of Jets read from good events                          = " <<
	h_NJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut              = " <<
	h_NJets_kmatCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
	h_NJets_JetIDCut->GetEntries()<<std::endl;

      std::cout<<std::endl;    }
  }// end evt counts
  
  


  std::string thePDFFileName=printQAPlots_dir;
  if(drawDataMCOverlaysInput)thePDFFileName+=output_PDFname_base_full+"DataMCOverlaid.pdf";
  else thePDFFileName+=output_PDFname_base_full+"DataMCRatios.pdf";
  
  std::string open_thePDFFileName, close_thePDFFileName;
  open_thePDFFileName=thePDFFileName+"[";  close_thePDFFileName=thePDFFileName+"]";
  
  // just for opening the pdf
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1200, 600);
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  temp_canvOpen->Close();



  // evt plots ----------------------
  if(!drawEvtQAPlots)std::cout<<std::endl<<"skipping evt QA plots..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;
    std::cout<<" drawing evt QA Plots now! "<<std::endl<<std::endl;

    if(debugMode)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
    TCanvas *temp_canvEvt_wLog = new TCanvas("tempEvt", "temp Evt Canv", 1200, 600);
    temp_canvEvt_wLog->SetLogy(1);
    temp_canvEvt_wLog->cd(); 
  
    if(debugMode)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
    TCanvas *temp_canvEvt = new TCanvas("tempEvt_wLog", "temp Evt Canv withLog", 1200, 600);
    temp_canvEvt->cd();

    if(drawDataMCOverlays){    
      if(debugMode)std::cout<<"drawing Data/MC overlays for Evt QA..."<<std::endl;
      
      {//MC and data total weighted vz
	std::string inHistName="hWeightedVz";
	long double theLumi, alt_theLumi;
	if(doEventCounts)theLumi=effIntgrtdLumi_vz;
	else theLumi=intgrtdLumi;
	if(doEventCounts)alt_theLumi=alt_effIntgrtdLumi_vz;
	
	std::cout<<" drawing Data+Evt Vz hist(s) "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title    ="Evt QA, vz";
	std::string h_XAx_Title="vz (cm)"                          ;
	std::string h_YAx_Title="millibarns/bin"                   ;
	
	TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1/theLumi );

	TH1F* theDataEvtQAHist_LowerJets= (TH1F*)finData_LowerJets->Get( inHistName.c_str() );
	theDataEvtQAHist_LowerJets->Scale( 1/theDataEvtQAHist_LowerJets->GetBinWidth(0) );
	theDataEvtQAHist_LowerJets->Scale( 1/alt_theLumi );
	
	theDataEvtQAHist->TH1::Add(theDataEvtQAHist_LowerJets);
	theDataEvtQAHist->SetMarkerStyle(kDot);
	theDataEvtQAHist->SetMarkerSize(1.3);
	theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
	theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );
	
	theDataEvtQAHist->SetTitle (    h_Title.c_str() );
	theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	theDataEvtQAHist->Draw();
	
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
	theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
	
	theMCEvtQAHist->SetMarkerStyle(kDot);
	theMCEvtQAHist->SetMarkerSize(1.3);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	
	theMCEvtQAHist->Draw("same");
	
	TLegend* theEvtQALeg=new TLegend(0.72,0.84,0.98,0.98, NULL,"brNDC");//new TLegend(0.75,0.75,0.90,0.85, NULL,"brNDC");
	theEvtQALeg->AddEntry(theDataEvtQAHist, "2015 promptReco 5.02 TeV ppData" ,"lp");
	theEvtQALeg->AddEntry(theMCEvtQAHist,   "Py8 CUETP8M1 QCD, ppMC"       ,"lp");
	

 


	//{//pthat weighted only MC, can use this for Vz weighting
	//  TH1F* theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	//  theMCEvtQAHist2->Scale( 1/theMCEvtQAHist2->GetBinWidth(0) );
	//  theMCEvtQAHist2->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist2->Integral() );
	//  theMCEvtQAHist2->SetMarkerStyle(kDot);
	//  theMCEvtQAHist2->SetMarkerSize(1.3);
	//  theMCEvtQAHist2->SetMarkerColor( kBlack );
	//  theMCEvtQAHist2->SetLineColor( altOverlayLineColor );
	//  theMCEvtQAHist2->Draw("SAME");            
	//  theEvtQALeg->AddEntry(theMCEvtQAHist2,"MC w/o vz-weights","lp");	
	//}
	
	theEvtQALeg->Draw();
	temp_canvEvt->Print( thePDFFileName.c_str() ); 
      }// endtotal weighted vz
      
      
      
      
      {//MC pthat, unnweighted
	temp_canvEvt_wLog->cd();
	std::string inHistName="hpthat";
	std::cout<<" drawing MC pthat hist "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title    ="Evt QA, pthat, MC ONLY" ;
	std::string h_XAx_Title="pthat (GeV)"       ;
	std::string h_YAx_Title="NEvents/bin"       ;
	
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
	
	theMCEvtQAHist->SetTitle (    h_Title.c_str() );
	theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	theMCEvtQAHist->SetMarkerStyle(kDot);
	theMCEvtQAHist->SetMarkerSize(1.3);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	
	theMCEvtQAHist->Draw();
	
	//TLegend* theEvtQALeg=new TLegend(0.10,0.70,0.30,0.85, NULL,"brNDC");
	TLegend* theEvtQALeg=new TLegend(0.72,0.84,0.98,0.98, NULL,"brNDC");//new TLegend(0.75,0.75,0.90,0.85, NULL,"brNDC");
	theEvtQALeg->AddEntry(theMCEvtQAHist,"Py8 QCD Dijet MC","lp");
	
	theEvtQALeg->Draw();
	temp_canvEvt_wLog->Print( thePDFFileName.c_str() );	
      }
      
      
      
      {//MC pthat, weighted
	temp_canvEvt_wLog->cd();	
	std::string inHistName="hWeightedpthat";     
	std::cout<<" drawing MC weighted pthat hist "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title    ="Evt QA, weighted pthat, MC ONLY" ;
	std::string h_XAx_Title="pthat (GeV)"                ;
	std::string h_YAx_Title="millibarns/bin"             ;
	
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
	
	theMCEvtQAHist->SetTitle (    h_Title.c_str() );
	theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	theMCEvtQAHist->SetMarkerStyle(kDot);
	theMCEvtQAHist->SetMarkerSize(1.3);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	
	theMCEvtQAHist->Draw();
	
	TLegend* theEvtQALeg=new TLegend(0.72,0.84,0.98,0.98, NULL,"brNDC");//new TLegend(0.75,0.75,0.90,0.85, NULL,"brNDC");
	//TLegend* theEvtQALeg=new TLegend(0.10,0.70,0.30,0.85, NULL,"brNDC");
	theEvtQALeg->AddEntry(theMCEvtQAHist,"Py8 QCD Dijet MC","lp");
	
	theEvtQALeg->Draw();
	temp_canvEvt_wLog->Print( thePDFFileName.c_str() );	
      }
    }
    
    else if(drawDataMCRatios)
      {	std::string inHistName="hWeightedVz";
	
	std::cout<<" drawing DataMC Ratio hist "<<inHistName<<std::endl<<std::endl;
	
	long double theLumi, alt_theLumi;
	if(doEventCounts)theLumi=effIntgrtdLumi_vz;
	else theLumi=intgrtdLumi;
	if(doEventCounts)alt_theLumi=alt_effIntgrtdLumi_vz;

	std::string h_Title   ="Evt QA, vz, Data/MC";
	std::string h_XAx_Title="vz (cm)", h_YAx_Title="Data/MC";
	TLegend* theEvtQALeg=new TLegend(0.72,0.84,0.98,0.98, NULL,"brNDC");//new TLegend(0.10,0.70,0.30,0.85, NULL,"brNDC");
	
	//TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1/theLumi );

	TH1F* theDataEvtQAHist_LowerJets= (TH1F*)finData_LowerJets->Get( inHistName.c_str() );
	theDataEvtQAHist_LowerJets->Scale( 1/theDataEvtQAHist_LowerJets->GetBinWidth(0) );
	theDataEvtQAHist_LowerJets->Scale( 1/alt_theLumi );
	
	theDataEvtQAHist->TH1::Add(theDataEvtQAHist_LowerJets);
	theDataEvtQAHist->SetTitle (    h_Title.c_str() );
	theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	theDataEvtQAHist->SetMarkerStyle(kDot);
	theDataEvtQAHist->SetMarkerSize(1.3);
	theDataEvtQAHist->SetMarkerColor( kBlack);
	theDataEvtQAHist->SetLineColor( theRatioLineColor );
	theDataEvtQAHist->SetAxisRange(0.,1.5,"Y");

	//for ratio w/ 1/vzweights
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
	theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
	
	TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
	theRatio->Divide(theMCEvtQAHist);
	theRatio->Draw();
	
	theEvtQALeg->Draw();	
	
	
	////for ratio w/ pre-vz-weights MC, use this for Vz weighting
	//TH1F* theDataEvtQAHist2=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio2");
	//
	//TH1F* theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	//theMCEvtQAHist2->Scale( 1/theMCEvtQAHist2->GetBinWidth(0) );
	//theMCEvtQAHist2->Scale( theDataEvtQAHist2->Integral()/theMCEvtQAHist2->Integral() );
	//
	//TH1F *theRatio2=theDataEvtQAHist2;
	//theRatio2->SetLineColor( altRatioLineColor1 );	
	//theRatio2->Divide(theMCEvtQAHist2);	
	//theEvtQALeg->AddEntry(theRatio2,"MC not vz-weighted","lp"); 		
	//
	//theRatio2->Draw("SAME");
	
	
	////for ratio w/ vzweights
	//TH1F* theDataEvtQAHist3=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio3");
	//
	//TH1F* theMCEvtQAHist3= (TH1F*)finMC->Get( "hWeightedAnnaVz" );
	//theMCEvtQAHist3->Scale( 1/theMCEvtQAHist3->GetBinWidth(0) );
	//theMCEvtQAHist3->Scale( theDataEvtQAHist3->Integral()/theMCEvtQAHist3->Integral() );
   	//
	//TH1F *theRatio3=theDataEvtQAHist3;
	//theRatio3->SetLineColor( altRatioLineColor2 );	
	//theRatio3->Divide(theMCEvtQAHist3);	
	//theEvtQALeg->AddEntry(theRatio3,"Annas VzWeights","lp"); 	
	//theRatio3->Draw("SAME");
	
	temp_canvEvt->Print( thePDFFileName.c_str() );
      }
    temp_canvEvt->Close();  
  }//end total weighted vz
  
  
  
  
  
  
  // jet plots----------------------
  if(!drawJetQAPlots)std::cout<<std::endl<<"skipping jet QA plots..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing Jet plots..."<<std::endl;
    TCanvas *temp_canvJet = new TCanvas("tempJet", "tempJet", 1200, 600);
    temp_canvJet->cd();
    temp_canvJet->SetLogy(1);
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing no-log Jet plots..."<<std::endl;
    TCanvas *temp_canvJet_noLog = new TCanvas("tempJet_noLog", "tempJet_noLog", 1200, 600);
    temp_canvJet_noLog->cd();
    temp_canvJet_noLog->SetLogy(0);
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing ratio Jet plots..."<<std::endl;
    TCanvas *temp_canvJetRat = new TCanvas("tempJetRat", "tempJetRat", 1200, 600);
    temp_canvJetRat->cd();
    temp_canvJetRat->SetLogy(0);
    
    
    std::cout<<" drawing jet QA Plots..."<<std::endl;
    for(int j=0;j<N_vars;j++){ std::cout<<std::endl;

      std::cout<<std::endl<<" var ="<<var[j]<<std::endl;
      for(int i=0;i<2;i++){
	
	long double theLumi, alt_theLumi;
	if(i==0)theLumi=effIntgrtdLumi_kmat;
	else if(i==1)theLumi=effIntgrtdLumi_JetID;

	if(i==0)alt_theLumi=alt_effIntgrtdLumi_kmat;
	else if(i==1)alt_theLumi=alt_effIntgrtdLumi_JetID;
	


	std::string inHistName="hJetQA_"+std::to_string(i)+"wJetID_"+var[j];

	//open the hists
	if(debugMode)std::cout<<" opening input Data "<<inHistName<<std::endl<<std::endl;

	TH1F* theDataJetQAHist= (TH1F*)
	  ( (TH1*)finData->Get(inHistName.c_str()) )
	  ->TH1::Rebin(var_xAx_reBin[j]);	
	theDataJetQAHist->Print("base"); 
	std::cout<<std::endl;
	
	theDataJetQAHist->Scale( 1/theDataJetQAHist->GetBinWidth(1) );
	theDataJetQAHist->Scale( 1/theLumi);
	theDataJetQAHist->Print("base"); 
	std::cout<<std::endl;

	TH1F* theDataJetQAHist_LowerJets= (TH1F*)
	  ( (TH1*)finData_LowerJets->Get(inHistName.c_str()) )
	  ->TH1::Rebin(var_xAx_reBin[j]);	
	theDataJetQAHist_LowerJets->Print("base"); 
	std::cout<<std::endl;
	
	theDataJetQAHist_LowerJets->Scale( 1/theDataJetQAHist_LowerJets->GetBinWidth(1) );
	theDataJetQAHist_LowerJets->Scale( 1/alt_theLumi);
	theDataJetQAHist_LowerJets->Print("base"); 
	std::cout<<std::endl;

	theDataJetQAHist->TH1::Add(theDataJetQAHist_LowerJets);
	theDataJetQAHist->SetMarkerStyle(kDot);
	theDataJetQAHist->SetMarkerSize(1.3);
	theDataJetQAHist->SetMarkerColor( theDataOverlayMarkerColor );
	theDataJetQAHist->SetLineColor( theDataOverlayLineColor );
	

	if(debugMode)std::cout<<" opening input MC "<<inHistName<<std::endl<<std::endl;

	TH1F* theMCJetQAHist= (TH1F*)
	  ( (TH1*)finMC->Get(inHistName.c_str()) )
	  ->TH1::Rebin(var_xAx_reBin[j]);	
	theMCJetQAHist->Print("base");
	std::cout<<std::endl;

	theMCJetQAHist->Scale( 1/theMCJetQAHist->GetBinWidth(1) );	  
	if( var[j]!="jtpt"&&var[j]!="rawpt"&&
	    var[j]!="leadJetPt"&&var[j]!="subleadJetPt" )  
	  theMCJetQAHist->Scale( theDataJetQAHist->Integral()/theMCJetQAHist->Integral() );	  
	theMCJetQAHist->Print("base");
	std::cout<<std::endl;
 
	theMCJetQAHist->SetMarkerStyle(kDot);
	theMCJetQAHist->SetMarkerSize(1.3);
	theMCJetQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCJetQAHist->SetLineColor( theMCOverlayLineColor );
	
	//assert(false);
	
	

	if(drawDataMCOverlays){  
	  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title="millibarns/bin"; 
	  if(debugMode)std::cout<<" x-axis title= "<<h_XAx_Title<<std::endl<<std::endl;

	  std::string h_Title="Jet QA, "+var[j];
	  if(i==1)h_Title+=", w/ JetIDCut";      
	  //else h_Title+=", no JetIDCut";
	  
	  theDataJetQAHist->SetTitle (    h_Title.c_str() );
	  theDataJetQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theDataJetQAHist->SetYTitle( h_YAx_Title.c_str() );

	  theMCJetQAHist->SetTitle (    h_Title.c_str() );
	  theMCJetQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theMCJetQAHist->SetYTitle( h_YAx_Title.c_str() );


	  TLegend* theJetQALeg=new TLegend(0.72,0.84,0.98,0.98, NULL,"brNDC");
	  theJetQALeg->AddEntry(theDataJetQAHist,"2015 promptReco 5.02 TeV ppData","lp");
	  theJetQALeg->AddEntry(theMCJetQAHist,  "Py8 CUETP8M1 QCD, ppMC"      ,"lp");
	 
	  theDataJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 

	  // variable specific draw orders, canv, etc.
	  if(var[j]=="jtpt"||var[j]=="rawpt") { 
	    temp_canvJet->cd();
	    theMCJetQAHist->Draw(); theDataJetQAHist->Draw("same");  
	    theJetQALeg->Draw(); 
	    temp_canvJet->Print( thePDFFileName.c_str() ); }		  
	  else if(var[j]=="xj"||var[j]=="Aj"||
		  var[j]=="jtphi"||var[j]=="jteta") { 
	    temp_canvJet_noLog->cd();
	    theDataJetQAHist->Draw(); theMCJetQAHist->Draw("SAME"); 
	    theJetQALeg->Draw();
	    temp_canvJet_noLog->Print( thePDFFileName.c_str() ); }		  
	  //else if(var[j]=="dphi") { temp_canvJet->cd();
	  //  theDataJetQAHist->Draw(); theMCJetQAHist->Draw("SAME");  
	  //  theJetQALeg->Draw();
	  //  temp_canvJet->Print( thePDFFileName.c_str() ); }		  
	  else { temp_canvJet->cd();
	    theDataJetQAHist->Draw(); theMCJetQAHist->Draw("SAME");  
	    theJetQALeg->Draw();
	    temp_canvJet->Print( thePDFFileName.c_str() ); }		  	  
	  std::cout<<std::endl; }
	
	
	if(drawDataMCRatios){ temp_canvJetRat->cd();

	  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title="Data/MC"    ;
	  std::string h_Title   ="Jet QA, Data/MC, "+var[j];
	  if(i==1)h_Title+=", w/ JetIDCut";      
	  //else h_Title+=", no JetIDCut";
	  

	  TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");

	  theRatio->SetTitle (    h_Title.c_str() );
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle( h_YAx_Title.c_str() );
	  theRatio->SetMarkerColor( kBlack );
	  theRatio->SetLineColor( theRatioLineColor );

	  theRatio->SetAxisRange(0,2.,"Y");
	  theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");

	  TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
	  TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
	  TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
	           lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	       lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);

	  theRatio->Divide(theMCJetQAHist);
	  theRatio->Draw();

	  lineAtOne->Draw("same");
	  lineAtOneHalf->Draw("same");
	  lineAtOneEtOneHalf->Draw("same");

	  
	  TLegend* theJetQALeg=new TLegend(0.72,0.84,0.98,0.98, NULL,"brNDC");//new TLegend(0.81,0.84,0.96,0.94, NULL,"brNDC");
	  theJetQALeg->AddEntry(theRatio,"ppData/Py8 QCD ppMC","lp");
	  theJetQALeg->Draw();
	  
	  temp_canvJetRat->Print( thePDFFileName.c_str() );	
	  std::cout<<std::endl;	    	
	} 
      } }
    //if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl<<std::endl;
    //temp_canvJet->Print( close_thePDFFileName.c_str() );//close the output file
    
    temp_canvJet_noLog->Close();
    temp_canvJetRat->Close();
    temp_canvJet->Close();
  }// end Jet QA Plots
  
  




  // jet trig plots----------------------
  if(!drawJetTrigQAPlots)std::cout<<std::endl<<"skipping jet trig QA plots..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;
    std::cout<<" drawing jet trig QA Plots..."<<std::endl;

    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
    TCanvas *temp_canvJetTrig = new TCanvas("tempJetTrig", "tempJetTrig", 1200, 600);
    temp_canvJetTrig->cd();
    temp_canvJetTrig->SetLogy(1);    

    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
    TCanvas *temp_canvJetRatio = new TCanvas("tempJetRatio", "tempJetRatio", 1200, 600);
    temp_canvJetRatio->cd();
    temp_canvJetRatio->SetLogy(0);

    
    for(int i=0; i<2; i++){ if(debugMode)std::cout<<std::endl<<"i= "<<i<<std::endl;
      
      long double theLumi;
      if(i==0)theLumi=effIntgrtdLumi_kmat;
      else theLumi=effIntgrtdLumi_JetID;
      
      TLegend* JetTrigLegend=new TLegend(0.84,0.72,0.98,0.98, NULL,"brNDC");
      
      if(drawTrigDataOverlays){   temp_canvJetTrig->cd();
	if(debugMode)std::cout<<"drawing dataTrig QA Overlay Plots"<<std::endl;

	for(int j=0; j<N_trigs; j++){ 
	  
	  /////////////////////
	  // PROBLEM SECTION, do NOT delete comments or other commented code
	  if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
	  
	  if(debugMode)std::cout<<std::endl<<"N_trigs="<<N_trigs<<" and j="<<j<<std::endl;
	  if(debugMode)std::cout<<"j<(N_trigs-1)="<<(bool)(j<(N_trigs-1))<<std::endl;
	  if(debugMode)std::cout<<"j<N_trigs="<<(bool)(j<(N_trigs))<<std::endl;
	  if(debugMode)std::cout<<"j<(N_trigs+1)="<<(bool)(j<(N_trigs+1))<<std::endl;	
	  
	  //assert(j!=N_trigs);// NEVER FIRES
	  //assert(j<N_trigs); // NEVER FIRES
	  
	  //open the hists + do scaling 
	  if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
	  
	  if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
	    std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
	    continue;}  
	  else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
	  
	  //if(j==N_trigs){  // uncommenting this if/else statement fixes the pathological behavior
	  //  std::cout<<"warning! j==N_trigs, about to seg fault..1."<<std::endl;
	  //  continue;}	
	  //else std::cout<<" HLTName ="<<HLTName[j]<<std::endl;
	    
	  //open the hists + do scaling 
	  if(j==N_trigs)std::cout<<"warning! j=="<< N_trigs<<"=N_trigs, about to seg fault..."<<std::endl;// NEVER FIRES
	  // PROBLEM SECTION
	  /////////////////////
	  
	  std::string inHistName="hpp_"+HLTName[j];
	  if(i==1)inHistName+="_wJetID";
	  inHistName+="_"+radius+etaWidth;
	  
	  if(debugMode)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
	  
	  //open the hists + do scaling 
	  TH1F* theJetTrigQAHist= (TH1F*)
	    ( (TH1*)finData->Get(inHistName.c_str()) )
	    ->TH1::Rebin(jetTrigQABinning);
	  
	  //theJetTrigQAHist->Print("base");
	  theJetTrigQAHist->Scale( 1/theJetTrigQAHist->GetBinWidth(1) );
	  theJetTrigQAHist->Scale( 1/theLumi);
	  
	  theJetTrigQAHist->SetMarkerStyle(kDot);
	  theJetTrigQAHist->SetMarkerSize(3.0);
	  theJetTrigQAHist->SetMarkerColor( theTrigOverlayLineColor[j]+j+1 );//theTrigOverlayMarkerColor );
	  theJetTrigQAHist->SetLineColor(   theTrigOverlayLineColor[j] );
	  
	  JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTName[j]).c_str(),"lp");
	  
	  if(j==0){
	    std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title="millibarns/bin"    ; 
	    std::string h_Title   ="Jet Trig QA, DATA ONLY";
	    if(i==1)h_Title+=", w/ JetIDCut";      
	    //else h_Title+=", no JetIDCut";
	    theJetTrigQAHist->SetAxisRange(0.,800.,"X");
	    theJetTrigQAHist->SetTitle (    h_Title.c_str() );
	    theJetTrigQAHist->SetXTitle( h_XAx_Title.c_str() );
	    theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
	    theJetTrigQAHist->Draw(); 	    }
	  else theJetTrigQAHist->Draw("SAME");   	    
	}//end loop over trigs
	
	JetTrigLegend->Draw();
	temp_canvJetTrig->Print(thePDFFileName.c_str());
      }//end DRAW TRIG OVERLAY
      

      
      if(drawTrigDataRatios){		temp_canvJetRatio->cd();
	if(debugMode)std::cout<<"drawing dataTrig QA Overlay Plots"<<std::endl;
	
	std::string comboPlotName="hpp_"+HLTName[4];
	if(i==1)comboPlotName+="_wJetID";
	comboPlotName+="_"+radius+etaWidth;
	
	const int jetTrigQARatioBinning=jetTrigQABinning/2;
	
	TH1F* theComboHist=(TH1F*)
	  ( (TH1*)finData->Get(comboPlotName.c_str()) )
	  ->TH1::Rebin(jetTrigQARatioBinning);
	
	for(int j=0; j<N_trigs; j++){ 
	  if(j==N_trigs-1)continue;//need to divide each curve by HLT Comb
	  
	  std::string inHistName="hpp_"+HLTName[j];
	  if(i==1)inHistName+="_wJetID";
	  inHistName+="_"+radius+etaWidth;
	  
	  TH1F* comboHistClone=(TH1F*)theComboHist->Clone("comboHistClone4Div");
	  
	  TH1F* theRatio=(TH1F*)
	    ( (TH1*)finData->Get(inHistName.c_str()) )
	    ->TH1::Rebin(jetTrigQARatioBinning);
	  
	  theRatio->SetMarkerStyle(kDot);
	  theRatio->SetMarkerSize(1.3);
	  theRatio->SetMarkerColor( theTrigOverlayMarkerColor );
	  theRatio->SetLineColor( theTrigOverlayLineColor[j] );
	  
	  theRatio->Divide(comboHistClone);
	  
	  JetTrigLegend->AddEntry(theRatio,(HLTName[j]+"/Trig Comb.").c_str(),"lp");
	  
	  if(j==0){
	    std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title="JetTrig/JetTrigCombo"    ; 
	    std::string h_Title   ="Jet Trig QA, DATA ONLY, Ratios";
	    if(i==1)h_Title+=", w/ JetIDCut";      
	    
	    //else h_Title+=", no JetIDCut";
	    theRatio->SetAxisRange(0.,1000.,"X");
	    theRatio->SetTitle (    h_Title.c_str() );
	    theRatio->SetXTitle( h_XAx_Title.c_str() );
	    theRatio->SetYTitle( h_YAx_Title.c_str() );
	    theRatio->Draw(); }
	  else theRatio->Draw("SAME");
	}//end trig name loop	  

	JetTrigLegend->Draw();
	temp_canvJetRatio->Print(thePDFFileName.c_str());
      }//end DRAW TRIG Ratios

    }//end jet ID loop, i=0,1
    temp_canvJetTrig->Close();        
    temp_canvJetRatio->Close();
  }//end else statement for het trig plots
  
  







  // MCEffQAPlots plots----------------------
  if(!drawMCEffPlots||!drawDataMCOverlays)std::cout<<std::endl<<"skipping MC Eff. QA plots..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;

    std::cout<<" drawing MC Eff. QA Plots..."<<std::endl;

    std::cout<<std::endl<<"creating temporary canvas for printing MCEff plots..."<<std::endl;
    TCanvas *temp_canvMCEff = new TCanvas("tempMCEff", "tempMCEff", 1200, 600);
    temp_canvMCEff->cd();

    if(drawDataMCOverlays){
      std::cout<<"drawDataMCOverlays must be true; it's value is "<<drawDataMCOverlays<<std::endl;

      for(int j=0; j<(N_genVars); j++){     
	std::cout<<std::endl<<"j= "<<j<<std::endl;
	for(int i=0; i<2; i++){      	    
	  
	  std::string inHistName="hpp_mceff_"+genVars[j];
	  //std::string inHistName="hpp_mceff_gen"+genVars[j];//for some slightly older samples, 9.7.2016
	  
	  if(i==1)inHistName+="_wJetID";
	  inHistName+="_"+radius+etaWidth;
	  
	  std::cout<<"opening TH2F "<<inHistName << std::endl;
	  TH2F* the2DMCEffQAHist= (TH2F*)finMC->Get( inHistName.c_str() );
	  the2DMCEffQAHist->Print("base");
	  
	  std::string h_Title   ="MC Eff. QA, TH2 profile";
	  if(i==1)h_Title+=", w/ JetIDCut";      
	  
	  //std::string h_ZAx_Title="Entries"; 
	  std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
	  std::string h_YAx_Title="avg "+stringRatio+"/bin";
	  
	  std::string h_XAx_Title="";	    
	  if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
	  if(j==1)h_XAx_Title="gen Jet Eta"       ;
	  if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
	  
	  std::cout<<"taking the profile of "<<inHistName << std::endl;
	  TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX("MCEffProfile",1,-1, "o");
	  
	  theMCEffQAHist->SetTitle (    h_Title.c_str() );
	  theMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
	  
	  theMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
	  if(j==0)theMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
	  if(j==1)theMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
	  if(j==2)theMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
	  
	  theMCEffQAHist->Draw(); 
	 
	  temp_canvMCEff->Print(thePDFFileName.c_str());   
	} // end JetID loop
      }// gen genvar loop

      for(int j=0; j<(N_genVars_ptrat); j++){     
	std::cout<<std::endl<<"j= "<<j<<std::endl;      
	for(int i=0; i<2; i++){      	    
	  //open the hists + do scaling 
	  std::string inHistName="hpp_mceff_ptrat_"+genVars_ptrat[j];
	  if(i==1)inHistName+="_wJetID";
	  inHistName+="_"+radius+etaWidth;
	  
	    TH2F* the2DMCEffQAHist= (TH2F*)finMC->Get( inHistName.c_str() );
	    
	    std::string h_Title   ="MC Eff. QA, TH2 profile";
	    if(i==1)h_Title+=", w/ JetIDCut";      
	    
	    std::string h_YAx_Title="avg (recpt)/(genpt)/bin";	    
	    std::string h_XAx_Title="";	    
	    if(j==0)h_XAx_Title="gen Jet Eta (GeV)"  ;
	    if(j==1)h_XAx_Title="gen Jet Phi (rad)"  ;
	    if(j==2)h_XAx_Title="gen Jet dRJet"      ; 
	    
	    TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX("MCEffProfile",1,-1, "o");
	    
	    theMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
	    if(j==0)theMCEffQAHist->SetAxisRange( -3.0, 3.0, "X");
	    if(j==1)theMCEffQAHist->SetAxisRange( -4.0, 4.0, "X");
	    if(j==2)theMCEffQAHist->SetAxisRange(   0., 0.5, "X");
	    
	    theMCEffQAHist->SetTitle (    h_Title.c_str() );
	    theMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
	    theMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
	    
	    theMCEffQAHist->Draw(); 
	    
	    temp_canvMCEff->Print(thePDFFileName.c_str());   
	}// end JetID loop
      }//end genvar ptrat loops

    }//end draw data/MCoverlays    
    temp_canvMCEff->Close();        
  }//end  MCEff Plots else
  
  
  
  
  
  
    
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1200, 600);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
    

  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData->Close();
  finData_LowerJets->Close();
  
  return 0;
}// end printplots
				     
				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=7 ) {//no input arguments, error
    std::cout<<"do ./printQAPlots.exe <ppDataCondorDir> <ppDataFilename> <ppMCCondorDir> <ppMCFilename> <outputNameBase> <OverlayFlag>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  //if(argc==1)rStatus=printQAPlots();
  if(argc==7) rStatus=printQAPlots( (const std::string) argv[1], (const std::string) argv[2],
			     (const std::string) argv[3], (const std::string) argv[4],
			     (const std::string) argv[5], (bool) std::atoi(argv[6])    );
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
}



//// I/O ------------------------
////input
//const std::string input_ppData_condorDir=readForests_ppData_dir+
////"ppData_HighPtJetTrig_ak4CaloJets_HLT.ak4CaloJets_09-17-16__newJetID/";
////"ppData_HighPtJetTrig_ak4PFJets_HLT.ak4CaloJets_09-16-16__newJetID/";
////"ppData_HighPtJetTrig_ak3PFJets_HLT.ak4CaloJets_09-18-16__newJetID/";
////"ppData_HighPtJetTrig_ak4CaloJets_HLT.ak4PFJets_09-18-16__newJetID/";
//  "ppData_HighPtJetTrig_ak4PFJets_HLT.ak4PFJets_09-18-16__newJetID/";
////"ppData_HighPtJetTrig_ak3PFJets_HLT.ak4PFJets_09-18-16__newJetID/";
//
//const std::string input_ppData_Filename=input_ppData_condorDir+
////"HighPtJetTrig_ak4Calo-allFiles.root";
//  "HighPtJetTrig_ak4PF-allFiles.root";
////"HighPtJetTrig_ak3PF-allFiles.root";
//
//const std::string input_ppMC_condorDir=readForests_ppMC_dir+
////"ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4CaloJets_09-17-16__newJetID/";
//  "ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_09-16-16__newJetID/";  
////"ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_09-18-16__newJetID/";
//
//const std::string input_ppMC_Filename=input_ppMC_condorDir+
////"Py8_CUETP8M1_QCDjetAllPtBins_ak4Calo-allFiles.root";
//  "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
////"Py8_CUETP8M1_QCDjetAllPtBins_ak3PF-allFiles.root";

//const std::string output_PDFname_base_full=printQAPlots_dir+
////"printQAplots_9.18.16_newJetID_ak4CaloJets_HLT.ak4PFJets_QAPlots_";
//  "printQAplots_9.18.16_newJetID_ak4PFJets_HLT.ak4PFJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak3PFJets_HLT.ak4PFJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak4CaloJets_HLT.ak4CaloJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak4PFJets_HLT.ak4CaloJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak3PFJets_HLT.ak4CaloJets_QAPlots_";
  
 

