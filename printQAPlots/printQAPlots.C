#include "printQAPlots.h"


// I/O ------------------------

//input
const std::string input_ppData_Filename=
  //"readForests_ppData_5p02TeV_HighPtJetTrig_08-29-16__wEvtJetCounts_AjxjBugFix/HighPtJetTrig_ak4PF-allFiles.root";
  //"readForests_ppData_5p02TeV_HighPtJetTrig_08-30-16__wdphi_recpt2rawpt/HighPtJetTrig_ak4PF-allFiles.root";
  //"readForests_ppData_5p02TeV_HighPtJetTrig_08-31-16__jtpt60/HighPtJetTrig_ak4PF-allFiles.root";
  "readForests_ppData_5p02TeV_HighPtJetTrig_08-31-16__jtpt15_ldJetPt60/HighPtJetTrig_ak4PF-allFiles.root";
//const std::string ppData_fullFilename=CMSSW_BASE+input_ppData_dir+condor_dir+input_ppData_Filename;
const std::string ppData_fullFilename=CMSSW_BASE+input_ppData_dir+"test_ppData.root";//debug

const std::string input_ppMC_Filename=
  //"readForests_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-29-16__noJECorJER_wmceff_wEvtJetCounts_AjxjBugFix/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
  //"readForests_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-29-16__divByVzWeightTry/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
  "readForests_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_08-31-16__jtpt15_ldJetPt60/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
//const std::string ppMC_fullFilename=CMSSW_BASE+input_ppMC_dir+condor_dir+input_ppMC_Filename;
const std::string ppMC_fullFilename=CMSSW_BASE+input_ppMC_dir+"test_ppMC.root";//debug


//output
const std::string outputDir=CMSSW_BASE+"src/printQAPlots/output/";
const std::string output_ppDataMCOverlaid_PDFname=outputDir+
  //"spectraOnly_sanityCheckForVZ_printQAPlots_ppDataMC_QA-jets-evts-trig_8.31.16.pdf";
  "spectraOnly_TEST_ldJtPtPlot_jteta.et.al_jtpt50_printQAPlots_ppDataMC_QA-jets-evts-trig_8.31.16.pdf";
const std::string output_ppDataMCRatio_PDFname=outputDir+
  //"ratiosOnly_sanityCheckForVZ_printQAPlots_ppDataMC_QA-jets-evts-trig_8.31.16.pdf";
  "ratiosOnly_TEST_ldJtPtPlot_jteta.et.al._jtpt50_printQAPlots_ppDataMC_QA-jets-evts-trig_8.31.16.pdf";


// switches ------------------------
const bool doEventCounts=true, doVzWeights=true;

//plot style
const bool drawDataMCOverlays= false;
const bool drawDataMCRatios =!(drawDataMCOverlays);//temporary? unsure
//const bool drawTrigOverlays=true, drawTrigRatios=false;//data only for now

//plot content
const bool drawEvtQAPlots=true, drawJetQAPlots=true;

////data-specific (for now)
//const bool drawTrigPlots=drawTrigOverlays||drawTrigRatios;
////mc-specific (fornow)
//const bool drawMCEffPlots=true, drawJECandJERPlots=true;//MC Specific for now

//for efficiency later if/when I get a chance 
const bool openDataFile=//doEventCounts||drawTrigPlots||
  drawDataMCOverlays||drawDataMCRatios;
const bool openMCFile  =//doEventCounts||//drawMCEffPlots||drawJECandJERPlots||
  drawDataMCOverlays||drawDataMCRatios;



// hist painting ------------------------
const int theDataOverlayMarkerColor=1, theMCOverlayMarkerColor=1,theRatioMarkerColor=1;//black
const int theDataOverlayLineColor=2, theMCOverlayLineColor=4,theRatioLineColor=9;//red, blue, olive
const int altOverlayLineColor=3; 
const int altRatioLineColor1=8, altRatioLineColor2=7;//green, violet or something



// Lumi Efficiency Calculation ------------------------
long double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
long double effIntgrtdLumi_kmat=1., LumiEff_kmat=1.;// eff-> effective, Eff->efficiency
long double effIntgrtdLumi_JetID=1., LumiEff_JetID=1.;// eff-> effective, Eff->efficiency


// the macro ------------------------
const bool debugMode=true;//, dataDebugMode=true, MCDebugMode=true; //debug

int printQAPlots(){
 
  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;

  //input file(s)
  TFile *finMC=NULL, *finData=NULL;
  assert( openMCFile||openDataFile );

  if(openDataFile&&debugMode)  std::cout<<" now opening Data File "<<std::endl<<ppData_fullFilename<<std::endl<<std::endl;
  if(openDataFile) finData = new TFile(ppData_fullFilename.c_str());

  if(openMCFile&&debugMode)  std::cout<<" now opening MC File "<<std::endl<<ppMC_fullFilename<<std::endl<<std::endl;
  if(openMCFile) finMC = new TFile(ppMC_fullFilename.c_str());
  

  //get output file ready
  std::string thePDFFileName, open_thePDFFileName, close_thePDFFileName;

  if(drawDataMCOverlays)thePDFFileName=output_ppDataMCOverlaid_PDFname;  
  else if(drawDataMCRatios)thePDFFileName=output_ppDataMCRatio_PDFname; 
  else thePDFFileName="GENERIC_printQAPlots_OUTPUTNAME.pdf";

  open_thePDFFileName=thePDFFileName+"[";  close_thePDFFileName=thePDFFileName+"]";

  TCanvas *temp_canv = new TCanvas("temp", "temp", 1200, 600);// just for opening the pdf
                                                              // canvases get weird if i abuse them
                                                              // so i close+make new ones as i go along
  temp_canv->Print( open_thePDFFileName.c_str() );        
  temp_canv->Close();
  // evt counts----------------------

  if(!doEventCounts)std::cout<<"skipping evt QA plots..."<<std::endl<<std::endl;
  else{ std::cout<<" grabbing raw event counts "<<std::endl;

    const bool doDataCounts=true, doMCCounts=false;
    if(doDataCounts){ std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;

      TH1F *h_NEvents         = (TH1F*)finData->Get("NEvents");
      TH1F *h_NEvents_read    = (TH1F*)finData->Get("NEvents_read");

      std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;


      TH1F *h_NEvents_skimCut = (TH1F*)finData->Get("NEvents_skimCut");
      TH1F *h_NEvents_vzCut   = (TH1F*)finData->Get("NEvents_vzCut");
      
      std::cout<<"Total Num of Events read passing skimCuts                              = " <<
	h_NEvents_skimCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
	h_NEvents_vzCut->GetEntries()<<std::endl;
      std::cout<<std::endl;


      TH1F *h_NEvents_withJets           = (TH1F*)finData->Get("NEvents_withJets");
      TH1F *h_NEvents_withJets_kmatCut   = (TH1F*)finData->Get("NEvents_withJets_kmatCut");
      TH1F *h_NEvents_withJets_JetIDCut = (TH1F*)finData->Get("NEvents_withJets_JetIDCut");
      
      std::cout<<"Total Num of good Events w/ jets                             = " <<
	h_NEvents_withJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
	h_NEvents_withJets_kmatCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
	h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;      
      
      LumiEff_vz = h_NEvents_vzCut->GetEntries()/h_NEvents->GetEntries();
      effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;
      LumiEff_kmat = h_NEvents_withJets_kmatCut->GetEntries()/h_NEvents->GetEntries();
      effIntgrtdLumi_kmat=intgrtdLumi*LumiEff_kmat;
      LumiEff_JetID = h_NEvents_withJets_JetIDCut->GetEntries()/h_NEvents->GetEntries();
      effIntgrtdLumi_JetID=intgrtdLumi*LumiEff_JetID;

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
  




  // evt plots ----------------------
  if(debugMode)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas("tempEvt", "tempEvt", 1200, 600);
  temp_canvEvt->cd();



  if(!drawEvtQAPlots)std::cout<<std::endl<<"skipping evt QA plots..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;
    std::cout<<" drawing evt QA Plots now! "<<std::endl<<std::endl;

    {//MC and data total weighted vz
      std::string inHistName="hWeightedVz";
      long double theLumi;
      if(doEventCounts)theLumi=effIntgrtdLumi_vz;
      else theLumi=intgrtdLumi;
      
      if(drawDataMCOverlays){  if(debugMode)std::cout<<" drawing DataMC Overlay hist "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title   ="Data+MC, Evt QA, vz";
	std::string h_XAx_Title="vz (cm)"                          ;
	std::string h_YAx_Title="millibarns/bin"                   ;
	

	TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1/theLumi );

	theDataEvtQAHist->SetMarkerStyle(kDot);
	theDataEvtQAHist->SetMarkerSize(1.1);
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
	theMCEvtQAHist->SetMarkerSize(1.1);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );

	theMCEvtQAHist->Draw("SAME");


	TLegend* theEvtQALeg=new TLegend(0.10,0.70,0.30,0.85, NULL,"brNDC");
	theEvtQALeg->AddEntry(theDataEvtQAHist,"HighPt Jet Trig Data","lp");
	theEvtQALeg->AddEntry(theMCEvtQAHist,"Py8 QCD Dijet MC","lp");
	
	{//pthat weighted only MC, can use this for Vz weighting
	  TH1F* theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	  theMCEvtQAHist2->Scale( 1/theMCEvtQAHist2->GetBinWidth(0) );
	  theMCEvtQAHist2->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist2->Integral() );
	  theMCEvtQAHist2->SetMarkerStyle(kDot);
	  theMCEvtQAHist2->SetMarkerSize(1.1);
	  theMCEvtQAHist2->SetMarkerColor( kBlack );
	  theMCEvtQAHist2->SetLineColor( altOverlayLineColor );
	  theMCEvtQAHist2->Draw("SAME");            
	  theEvtQALeg->AddEntry(theMCEvtQAHist2,"MC w/o vz-weights","lp");	
	}
	theEvtQALeg->Draw();
	temp_canvEvt->Print( thePDFFileName.c_str() );	

      }
      
      
      if(drawDataMCRatios){  if(debugMode)std::cout<<" drawing DataMC Ratio hist "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title   ="Data/MC Ratio, Evt QA, vz";
	std::string h_XAx_Title="vz (cm)", h_YAx_Title="Data/MC";
	TLegend* theEvtQALeg=new TLegend(0.10,0.70,0.30,0.85, NULL,"brNDC");

	//TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1/theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1/theLumi );

	theDataEvtQAHist->SetTitle (    h_Title.c_str() );
	theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );

	theDataEvtQAHist->SetMarkerStyle(kDot);
	theDataEvtQAHist->SetMarkerSize(1.1);
	theDataEvtQAHist->SetMarkerColor( kBlack);
	theDataEvtQAHist->SetLineColor( theRatioLineColor );
	theDataEvtQAHist->SetAxisRange(0.,1.5,"Y");

	//for ratio w/ 1/vzweights
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1/theMCEvtQAHist->GetBinWidth(0) );
	theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
	
	TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
	theRatio->Divide(theMCEvtQAHist);
	theEvtQALeg->AddEntry(theRatio,"MC wght~1/vz-wght","lp");	
	theRatio->Draw();
	
	
	//for ratio w/ pre-vz-weights MC, use this for Vz weighting
	TH1F* theDataEvtQAHist2=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio2");

	TH1F* theMCEvtQAHist2= (TH1F*)finMC->Get( "hpthatWeightedVz" );
	theMCEvtQAHist2->Scale( 1/theMCEvtQAHist2->GetBinWidth(0) );
	theMCEvtQAHist2->Scale( theDataEvtQAHist2->Integral()/theMCEvtQAHist2->Integral() );
	
	TH1F *theRatio2=theDataEvtQAHist2;
	theRatio2->SetLineColor( altRatioLineColor1 );	
	theRatio2->Divide(theMCEvtQAHist2);	
	theEvtQALeg->AddEntry(theRatio2,"MC not vz-weighted","lp"); 		
	
	if(doVzWeights){
	  Int_t NvzWeightBins=60;
	  Float_t theVzBinWidth=theRatio2->TH1::GetBinWidth(0);
	  Float_t xLow=-15., xHigh=15.;

	  std::cout<<"now grabbing vzWeights for "<< NvzWeightBins<<"bins for "<<xLow<<"< vz <"<<xHigh << std::endl;
	  for (int i=1;i<NvzWeightBins+1;++i){//binsX loop
	    Float_t leftSideOfBin=-15.+(i-1)*theVzBinWidth;
	    Float_t rightSideOfBin=-15.+(i)*theVzBinWidth;
	    Float_t vzWeight = theRatio2->TH1::GetBinContent(i);
	    
	    if(i!=NvzWeightBins)
	      std::cout<<"  or i="<<i<<", "<<leftSideOfBin<<<<"<vz<="<<rightSideOfBin<<", vzWeight="<<vzWeight<< std::endl; 
	    else
	      std::cout<<"  or i="<<i<<", "<<leftSideOfBin<<<<"<vz<="<<rightSideOfBin<<", vzWeight="<<vzWeight<< std::endl; 
	  }}//end vz weights
	
	
	theRatio2->Draw("SAME");
	
	
	//for ratio w/ vzweights
	TH1F* theDataEvtQAHist3=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio3");
	
	TH1F* theMCEvtQAHist3= (TH1F*)finMC->Get( "hWeightedProdVz" );
	theMCEvtQAHist3->Scale( 1/theMCEvtQAHist3->GetBinWidth(0) );
	theMCEvtQAHist3->Scale( theDataEvtQAHist3->Integral()/theMCEvtQAHist3->Integral() );
   	
	TH1F *theRatio3=theDataEvtQAHist3;
	theRatio3->SetLineColor( altRatioLineColor2 );	
	theRatio3->Divide(theMCEvtQAHist3);	
	theEvtQALeg->AddEntry(theRatio3,"MC wgth~1*vz-wght","lp"); 	
	theRatio3->Draw("SAME");
	
	theEvtQALeg->Draw();	

	////if(debugMode)std::cout<<"trying out DrawLineNDC.."<<std::endl;
	//TLine* lineAtOne=new TLine();
	//lineAtOne->TLine::DrawLineNDC( 0., 1.*(1./1.5), 1., 1.*(1./1.5) );
	//lineAtOne->Draw();

	temp_canvEvt->Print( thePDFFileName.c_str() );	
      }    
    }//end total weighted vz
  }  //end evt Ratio, Overlaid plots
  temp_canvEvt->Close();
  



  // jet plots----------------------
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
  


  if(!drawJetQAPlots)std::cout<<std::endl<<"skipping jet QA plots..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;

    std::cout<<" drawing jet QA Plots..."<<std::endl;
    for(int j=0;j<N_vars;j++){ std::cout<<std::endl;

      std::cout<<std::endl<<" var ="<<var[j]<<std::endl;
      for(int i=0;i<2;i++){

	//if(var[j]=="dphi") continue; //for redoing a few plots 8.31.16

	long double theLumi;
	if(i==0)theLumi=effIntgrtdLumi_kmat;
	else theLumi=effIntgrtdLumi_JetID;

	std::string inHistName="hJetQA_"+std::to_string(i)+"wJetID_"+var[j];
	if(debugMode)std::cout<<" opening input Data+MC histos "<<inHistName<<std::endl<<std::endl;


	//open the hists + do scaling + set paint style
	TH1F* theDataJetQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	//theDataJetQAHist->Print("base");
	theDataJetQAHist->Scale( 1/theDataJetQAHist->GetBinWidth(0) );
	theDataJetQAHist->Scale( 1/theLumi);

	theDataJetQAHist->SetMarkerStyle(kDot);
	theDataJetQAHist->SetMarkerSize(1.1);
	theDataJetQAHist->SetMarkerColor( theDataOverlayMarkerColor );
	theDataJetQAHist->SetLineColor( theDataOverlayLineColor );


	TH1F* theMCJetQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	//theMCJetQAHist->Print("base");
	theMCJetQAHist->Scale( 1/theMCJetQAHist->GetBinWidth(0) );	  
	if(var[j]!="jtpt"&&var[j]!="rawpt")  
	  theMCJetQAHist->Scale( theDataJetQAHist->Integral()/theMCJetQAHist->Integral() );	  

	theMCJetQAHist->SetMarkerStyle(kDot);
	theMCJetQAHist->SetMarkerSize(1.1);
	theMCJetQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCJetQAHist->SetLineColor( theMCOverlayLineColor );

	

	if(drawDataMCOverlays){	
	  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title="millibarns/bin"    ; 
	  std::string h_Title   ="Data+MC, Jet QA, "+var[j];
	  if(i==1)h_Title+=", w/ JetIDCut";      
	  else h_Title+=", no JetIDCut";
	  
	  theDataJetQAHist->SetTitle (    h_Title.c_str() );
	  theDataJetQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theDataJetQAHist->SetYTitle( h_YAx_Title.c_str() );
	  theMCJetQAHist->SetTitle (    h_Title.c_str() );
	  theMCJetQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theMCJetQAHist->SetYTitle( h_YAx_Title.c_str() );


	  TLegend* theJetQALeg=new TLegend(0.81,0.84,0.96,0.99, NULL,"brNDC");
	  theJetQALeg->AddEntry(theDataJetQAHist,"HighPt Jet Trig Data","lp");
	  theJetQALeg->AddEntry(theMCJetQAHist,"Py8 QCD Dijet MC","lp");
	  
	  // variable specific ranges, draw orders, etc.
	  if(var[j]=="xj"||var[j]=="Aj") { temp_canvJet_noLog->cd();
	    theDataJetQAHist->SetAxisRange(0.,1.05,"X");
	    theDataJetQAHist->Draw();  
	    theMCJetQAHist->Draw("SAME"); 
	    theJetQALeg->Draw();
	    temp_canvJet_noLog->Print( thePDFFileName.c_str() ); }		  

	  else if(var[j]=="jtphi") { temp_canvJet_noLog->cd();
	    theDataJetQAHist->SetAxisRange(-4.0,4.0,"X"); 
	    theDataJetQAHist->Draw();  
	    theMCJetQAHist->Draw("SAME");  
	    theJetQALeg->Draw();
	    temp_canvJet_noLog->Print( thePDFFileName.c_str() ); }		  

	  else if(var[j]=="dphi") { temp_canvJet->cd();
	    theDataJetQAHist->SetAxisRange(0.,4.0,"X"); 
	    theDataJetQAHist->Draw();  
	    theMCJetQAHist->Draw("SAME");  
	    theJetQALeg->Draw();
	    temp_canvJet_noLog->Print( thePDFFileName.c_str() ); }		  

	  else if(var[j]=="jteta") { temp_canvJet_noLog->cd();
	    theDataJetQAHist->SetAxisRange(-2.5,2.5,"X"); 
	    theDataJetQAHist->Draw();  
	    theMCJetQAHist->Draw("SAME");  
	    theJetQALeg->Draw();
	    temp_canvJet_noLog->Print( thePDFFileName.c_str() ); }		  

	  else { temp_canvJet->cd();
	    theMCJetQAHist->Draw();  
	    theDataJetQAHist->Draw("SAME");  
	    theJetQALeg->Draw();
	    temp_canvJet->Print( thePDFFileName.c_str() ); }		  
	  
	  std::cout<<std::endl;	
	}
	
	
	if(drawDataMCRatios){ temp_canvJetRat->cd();
	  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title="Data/MC"    ;
	  std::string h_Title   ="Data/MC Ratio, Jet QA, "+var[j];
	  if(i==1)h_Title+=", w/ JetIDCut";      
	  else h_Title+=", no JetIDCut";
	  
	  TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
	  theRatio->SetTitle (    h_Title.c_str() );
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle( h_YAx_Title.c_str() );
	  theRatio->SetMarkerColor( kBlack );
	  theRatio->SetLineColor( theRatioLineColor );
	  theRatio->SetAxisRange(0,2.,"Y");
	  
	  
	  // variable specific ranges
	  if(var[j]=="jteta") theRatio->SetAxisRange(-2.5,2.5,"X"); 
	  else if(var[j]=="jtphi") theRatio->SetAxisRange(-4.0,4.0,"X"); 
	  else if(var[j]=="dphi") { theRatio->SetAxisRange(0.,4.0,"X"); }
	  else if(var[j]=="xj"||var[j]=="Aj") theRatio->SetAxisRange(0.,1.05,"X");

	  theRatio->Divide(theMCJetQAHist);

	  theRatio->Draw();
	  
	  //TLine* lineAtOne=new TLine();
	  //lineAtOne->TLine::DrawLineNDC( );
	  //lineAtOne->Draw();

	  //TLegend* theJetQALeg=new TLegend(0.81,0.84,0.96,0.99, NULL,"brNDC");
	  //theJetQALeg->AddEntry(theRatio,"HighPtJets Data/Py8 QCD MC","lp");
	  //theJetQALeg->Draw();
	  
	  temp_canvJetRat->Print( thePDFFileName.c_str() );	
	  std::cout<<std::endl;	    	
	}
      }    
    }
  }// end Jet QA Plots
  
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl<<std::endl;
  temp_canvJet->Print( close_thePDFFileName.c_str() );//close the output file

  temp_canvJet_noLog->Close();
  temp_canvJetRat->Close();
  temp_canvJet->Close();



  /*
  // trig plots ----------------------
  if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas("tempJetTrig", "tempJetTrig", 1200, 600);
  temp_canvJetTrig->cd();
  temp_canvJetTrig->SetLogy(1);
    // jet plots----------------------
  if(!drawJetTrigPlots)std::cout<<std::endl<<"skipping jetTrig QA plots..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;
    
    if()

  }


  temp_canv->SetLogy(1);
  // Data Only, HLT/L1/TrigComb Plots
  for(int j=0;j<N_HLTNames;j++){
    for(int i=0;i<2;i++){

      std::string theHistName="hpp_"+HLTName[j];
      if(i==0)theHistName+="no";
      theHistName+="JetID_R"+radius+"_"+(std::string)etaWidth;
      if(debugMode)std::cout<<"theHistName="<<theHistName<<std::endl;

      TH1F* theJetQAHist= (TH1F*)fin->Get( theHistName.c_str() );

      std::string theHistTitle=HLThTitle[j];
      if(i==0)theHistTitle+="w/o JetID R";
      else theHistTitle+="w/ JetID R";      
      theHistTitle+=radius+" "+(std::string)etaWidth;
      theJetQAHist->SetTitle( theHistTitle.c_str() );
      theJetQAHist->SetYTitle("Cross-section (millibarn)");

      theJetQAHist->SetAxisRange(0,1300.,"X");
      theJetQAHist->Scale(1/theJetQAHist->GetBinWidth(0));
      theJetQAHist->Scale(1/integratedLuminosity);

      theJetQAHist->Draw();
      temp_canv->Print( thePDFFileName.c_str() );
    }
  }
  //----------------------
  */



  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData->Close();

  return 0;
}// end printplots




// steering ----------------------

int main(int argc, char *argv[]){

  int rStatus = -1;
  if( argc!=1 ) {//no input arguments, error
    std::cout<<"debugging, for now just do ./printQAPlots.exe"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  rStatus=printQAPlots();
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
}
