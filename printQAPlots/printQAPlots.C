#include "printQAPlots.h"

// env and I/O ------------------------
const std::string CMSSW_BASE="/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/";//CMSSW_BASE
const std::string SCRATCH_BASE="/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis_archivedCondorOutput/readForests/10.18.16_outputCondor/";

//const std::string readForests_dir=CMSSW_BASE+"src/readForests/saved_outputCondor/";//input dir
const std::string readForests_dir=SCRATCH_BASE;

const std::string outputDir=CMSSW_BASE+"src/printQAPlots/output/";//output dir

// code/job switches ------------------------
//evt switches
const bool doEventCounts=true, drawEvtQAPlots=true;
//jet switches
const bool doJetIDPlots=true;
const bool drawJetQAPlots=true, drawJetTrigQAPlots=true, drawJetRapBinsPlot=true;
const bool drawMCEffPlots=false;//, drawJECandJERPlots=true;//MC Specific for now
//minor jet switches
const bool drawDijetPlots=false;
const bool doBasicJetQAOnly=true;

const bool debugMode=true;

// hist painting ------------------------
const int theDataOverlayMarkerColor=2, theMCOverlayMarkerColor=4,theRatioMarkerColor=9;//red, blue, olive
const int theDataOverlayLineColor=1, theMCOverlayLineColor=1,theRatioLineColor=1;//black

const int altOverlayLineColor=3; 
const int altRatioLineColor1=8, altRatioLineColor2=7;//green, violet or something

const int theTrigCombMarkerColor=1;
const int altTrigCombMarkerColor=12;
const int theTrigOverlayMarker[]={20,20,20,20,20,32};
const int theTrigOverlayMarkerColor[]={2,3,6,7,1,4};
const int theTrigOverlayLineColor[]={1,1,1,1,1,1};

const int theRapOverlayMarkerColor[]={2,3,6,7,1,4};

const std::string crossSectionYAxis="#sigma (#mub)";

// the macro ------------------------
int printQAPlots(const std::string input_ppData_condorDir , const std::string input_ppMC_condorDir , 
		 const std::string output_PDFname_base ,    const bool drawDataMCOverlaysInput=true ){
 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

  // PLOT STYLE
  const bool drawDataMCOverlays   = drawDataMCOverlaysInput;
  const bool drawDataMCRatios     = !drawDataMCOverlays;
  const bool drawMCEffProfiles       = drawMCEffPlots && drawDataMCOverlays;

  const bool drawTrigDataOverlays = drawDataMCOverlays && drawJetTrigQAPlots;
  const bool drawTrigDataRatios   = drawDataMCRatios   && drawJetTrigQAPlots;


  // MAKE STRINGS
  //format of the filename is always HighPt{filelist}_ak{3,4}{PF,Calo}-allFiles.root
  std::size_t jobTypePos=input_ppMC_condorDir.find("__")+2;  
  std::string jobType;
  if(jobTypePos!=std::string::npos) jobType=input_ppMC_condorDir.substr( jobTypePos);
  else                         jobType="noJobType";
  if(debugMode)std::cout<<"jobType string is = "<<jobType<<std::endl;

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
  const std::string input_ppData_Filename="/HighPtJetTrig_" +fullJetType+ "-allFiles.root";
  const std::string input_ppMC_Filename  ="/Py8_CUETP8M1_QCDjetAllPtBins_" +fullJetType+ "-allFiles.root";
  const std::string ppData_fullFilename=readForests_dir+input_ppData_condorDir+input_ppData_Filename;
  const std::string ppMC_fullFilename  =readForests_dir+input_ppMC_condorDir+input_ppMC_Filename;;  

  if(debugMode)std::cout<<"ppData Filename is = "<<input_ppData_Filename<<std::endl;
  if(debugMode)std::cout<<"ppMC Filename is = "<<input_ppMC_Filename<<std::endl;



  // OPEN INPUT SECTION
  TFile *finMC=NULL, *finData=NULL;  
  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;

  //open the input files
  if(debugMode)std::cout<<" now opening Data File "<<std::endl<<ppData_fullFilename<<std::endl<<std::endl;
  else std::cout<<" now opening Data File "<<std::endl<<input_ppData_Filename<<std::endl<<std::endl;
  finData = new TFile(ppData_fullFilename.c_str());      
  
  if(debugMode) std::cout<<" now opening MC File "<<std::endl<<ppMC_fullFilename<<std::endl<<std::endl;
  else std::cout<<" now opening MC File "<<std::endl<<input_ppMC_Filename<<std::endl<<std::endl;
  finMC = new TFile(ppMC_fullFilename.c_str());

  assert( finData && finMC );
  



  // Lumi Eff + evt + jet counts ----------------------
  long double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  long double effIntgrtdLumi_kmat=1., LumiEff_kmat=1.;
  long double effIntgrtdLumi_JetID=1., LumiEff_JetID=1.;
  
  bool didJetIDCut;
  
  if(!doEventCounts)std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
  else{ std::cout<<" grabbing raw event counts "<<std::endl;
    
    for(int file=0;file<2;file++){ 
      
      TFile * fin=NULL;
      if(file==0)fin=finData;
      else fin=finMC;
      
      if(!fin){	std::cout<<"input file not found, cannot look at event counts"<<std::endl;
	std::cout<<"cannot compute integrated luminosity for file=="<<file<<std::endl;  continue;}
      else{
	
	TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
	TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");      
	TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
	TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");
	//TH1F *h_NEvents_trgDup   = (TH1F*)fin->Get("NEvents_trgDup");   
	//TH1F *h_NEvents_trgDupCut= (TH1F*)fin->Get("NEvents_trgDupCut");      

	TH1F *h_NEvents_withJets          = (TH1F*) fin->Get("NEvents_withJets");
	TH1F *h_NEvents_withJets_kmatCut  = (TH1F*) fin->Get("NEvents_withJets_kmatCut");	
	TH1F *h_NEvents_withJets_JetIDCut = (TH1F*) fin->Get("NEvents_withJets_JetIDCut");

	didJetIDCut=false;
	if(h_NEvents_withJets_JetIDCut)didJetIDCut=true&&doJetIDPlots;
	if (doJetIDPlots) assert(didJetIDCut);

	if(file==0)std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;      
	else std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;      
	std::cout<<"Total Num of Events in file(s) opened       = " <<
	  h_NEvents->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events read from those file(s) = " <<
	  h_NEvents_read->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events read passing skimCuts   = " <<
	  h_NEvents_skimCut->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<
	  h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
	//std::cout<<"Total Num of duplicate Events read post skim+vz cuts                 = " <<
	//	h_NEvents_trgDup->GetEntries()<<std::endl;
	//std::cout<<"Total Num of Events read passing skim+vz+duplicate (aka good events) = " <<
	//	h_NEvents_trgDupCut->GetEntries()<<std::endl;
	//std::cout<<std::endl;
	
      std::cout<<"Total Num of good Events w/ jets                             = " <<
	h_NEvents_withJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of good Events, w/ jets, post kmatCuts             = " <<
	h_NEvents_withJets_kmatCut->GetEntries()<<std::endl;
      if(didJetIDCut)std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " <<
	h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;      
      
      TH1F *h_NJets          = (TH1F*)fin->Get("NJets");
      TH1F *h_NJets_kmatCut  = (TH1F*)fin->Get("NJets_kmatCut");
      TH1F *h_NJets_JetIDCut = (TH1F*)fin->Get("NJets_JetIDCut");
      
      std::cout<<std::endl<<"/// Data Jet-Count Summary Summary ///"<<std::endl<<std::endl;
      std::cout<<"Total Num of Jets read from good events                          = " <<
	h_NJets->GetEntries()<<std::endl;
      std::cout<<"Total Num of Jets read from good events post kinCut              = " <<
	h_NJets_kmatCut->GetEntries()<<std::endl;
      if(didJetIDCut)std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " <<
	h_NJets_JetIDCut->GetEntries()<<std::endl;
      std::cout<<std::endl;    

      if(file==0){//data file evt counts

	LumiEff_vz = h_NEvents_vzCut->GetEntries()/h_NEvents->GetEntries();
	LumiEff_kmat = h_NEvents_withJets_kmatCut->GetEntries()/h_NEvents->GetEntries();
	if(didJetIDCut) LumiEff_JetID = h_NEvents_withJets_JetIDCut->GetEntries()/h_NEvents->GetEntries();
	
	effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;
	effIntgrtdLumi_kmat=intgrtdLumi*LumiEff_kmat;
	if(didJetIDCut)effIntgrtdLumi_JetID=intgrtdLumi*LumiEff_JetID;
	
	std::cout<<std::endl<<"dataset integrated Luminosity (microbarns) ="<<intgrtdLumi<<std::endl;
	std::cout<<std::endl<<"lumi efficiency, vz cuts="<<LumiEff_vz<<std::endl;
	std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;
	std::cout<<std::endl<<"lumi efficiency, kmatCuts="<<LumiEff_kmat<<std::endl;
	std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_kmat<<std::endl;
	if(didJetIDCut) std::cout<<std::endl<<"lumi efficiency, JetID="<<LumiEff_JetID<<std::endl;
	if(didJetIDCut) std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_JetID<<std::endl;            }
      }
    }// end data/MC file loop
  }// end evt counts
  

  

  // OUTPUT PDF FILE
  std::string thePDFFileName=outputDir+fullJetType+"_"+jobType;
  if(drawDataMCOverlaysInput)thePDFFileName+="_"+output_PDFname_base+"DataMCOverlaid.pdf";
  else thePDFFileName+="_"+output_PDFname_base+"DataMCRatios.pdf";//one of these two will ALWAYS fire
  
  std::string open_thePDFFileName, close_thePDFFileName;
  open_thePDFFileName=thePDFFileName+"[";  close_thePDFFileName=thePDFFileName+"]";
  
  // just for opening the pdf
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1200, 600);
  temp_canvOpen->Print( open_thePDFFileName.c_str() );  temp_canvOpen->Close();


  
  // evt plots ----------------------
  if(!drawEvtQAPlots)std::cout<<std::endl<<"skipping evt QA plots..."<<std::endl<<std::endl;
  else if(doJetIDPlots)std::cout<<std::endl<<"doing JetID plots, skipping evt QA plots..."<<std::endl<<std::endl;
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
	long double theLumi;
	if(doEventCounts)theLumi=effIntgrtdLumi_vz;
	else theLumi=intgrtdLumi;
	
	std::cout<<" drawing Data+Evt Vz hist(s) "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title    ="EvtQA";
	std::string h_XAx_Title="v_{z}^{evt} (cm)"                          ;
	std::string h_YAx_Title=crossSectionYAxis;//"#sigma (#mub/bin)";
	
	TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1./theLumi );
	
	theDataEvtQAHist->SetMarkerStyle(kFullCircle);
	//theDataEvtQAHist->SetMarkerSize(1.8);
	theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
	theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );
	
	
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(0) );
	theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
	
	theMCEvtQAHist->SetMarkerStyle(kFullCircle);
	//	theMCEvtQAHist->SetMarkerSize(1.8);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	
	theMCEvtQAHist->SetTitle (    h_Title.c_str() );
	theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	//theMCEvtQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
	
	theMCEvtQAHist->Draw();
	theDataEvtQAHist->Draw("same");	

	float legx1=0.78, legx2=legx1+0.07;
	float legy1=0.76, legy2=legy1+0.11;
	TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
	theEvtQALeg->AddEntry(theDataEvtQAHist, "data" ,"lp");
	theEvtQALeg->AddEntry(theMCEvtQAHist,   "MC"       ,"lp");

	theEvtQALeg->Draw();
	
	float pp1x=0.14,pp1y=0.85;
	TLatex *pp1=new TLatex(pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets");
	pp1->SetTextFont(63);
	pp1->SetTextColor(kBlack);
	pp1->SetTextSize(19);
	pp1->SetLineWidth(1);
	pp1->SetNDC();
	pp1->Draw();
	TLatex *pp2=new TLatex(pp1x,(pp1y-0.05),"L_{int} = 25.8 pb^{-1}, #sqrt{s} = 5.02 TeV");
	pp2->SetTextFont(63);
	pp2->SetTextColor(kBlack);
	pp2->SetTextSize(19);
	pp2->SetLineWidth(1);
	pp2->SetNDC();
	pp2->Draw();
	//TLatex *pp3=new TLatex(pp1x,(pp1y-0.10),(fullJetType+"Jets").c_str() );
	//pp3->SetTextFont(63);
	//pp3->SetTextColor(kBlack);
	//pp3->SetTextSize(19);
	//pp3->SetLineWidth(1);
	//pp3->SetNDC();
	//pp3->Draw();
	TLatex *pp4=new TLatex(pp1x,(pp1y-0.10),"Py8 Tune CUETP8M1" );
	pp4->SetTextFont(63);
	pp4->SetTextColor(kBlack);
	pp4->SetTextSize(19);
	pp4->SetLineWidth(1);
	pp4->SetNDC();
	pp4->Draw();

	temp_canvEvt->Print( thePDFFileName.c_str() ); 

      }// endtotal weighted vz

      
      {//MC pthat, unnweighted
	temp_canvEvt_wLog->cd();
	std::string inHistName="hpthat";
	std::cout<<" drawing MC pthat hist "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title    ="MC EvtQA";
	std::string h_XAx_Title="MC evt #hat{p}_{t} (GeV)"       ;
	std::string h_YAx_Title="A.U."       ;
	
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(0) );
	
	theMCEvtQAHist->SetTitle (    h_Title.c_str() );
	theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	//theMCEvtQAHist->SetMarkerStyle(kFullCircle);
	theMCEvtQAHist->SetMarkerStyle(kDot);
	//theMCEvtQAHist->SetMarkerSize(1.8);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	
	theMCEvtQAHist->Draw();
	
	
	float t1x=.70,t1y=.84;
	TLatex *t1=new TLatex(t1x,t1y,"QCD MC, unweighted");
	t1->SetTextFont(63);
	t1->SetTextColor(kBlack);
	t1->SetTextSize(22);
	t1->SetLineWidth(1);
	t1->SetNDC();
	t1->Draw();	

	TLatex *t2=new TLatex(t1x,t1y-0.05,"Py8 Tune CUETP8M1");
	t2->SetTextFont(63);
	t2->SetTextColor(kBlack);
	t2->SetTextSize(22);
	t2->SetLineWidth(1);
	t2->SetNDC();
	t2->Draw();	

	//theEvtQALeg->Draw();
	temp_canvEvt_wLog->Print( thePDFFileName.c_str() );	
      }
      
      {//MC pthat, weighted
	temp_canvEvt_wLog->cd();	
	std::string inHistName="hWeightedpthat";     
	std::cout<<" drawing MC weighted pthat hist "<<inHistName<<std::endl<<std::endl;
	
	std::string h_Title    ="MC EvtQA" ;
	std::string h_XAx_Title="MC evt #hat{p}_{t} (GeV)"                ;
	std::string h_YAx_Title=crossSectionYAxis;//"#sigma (#mub/bin)"             ;
	
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(0) );
	
	theMCEvtQAHist->SetTitle (    h_Title.c_str() );
	theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	//theMCEvtQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
	
	//theMCEvtQAHist->SetMarkerStyle(kFullCircle);
	theMCEvtQAHist->SetMarkerStyle(kDot);
	//theMCEvtQAHist->SetMarkerSize(1.8);
	theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	
	theMCEvtQAHist->Draw();
	
	
	float t1x=.70,t1y=.84;
	TLatex *t1=new TLatex(t1x,t1y,"QCD MC, weighted");
	t1->SetTextFont(63);
	t1->SetTextColor(kBlack);
	t1->SetTextSize(22);
	t1->SetLineWidth(1);
	t1->SetNDC();
	t1->Draw();	

	TLatex *t2=new TLatex(t1x,t1y-0.05,"Py8 Tune CUETP8M1");
	t2->SetTextFont(63);
	t2->SetTextColor(kBlack);
	t2->SetTextSize(22);
	t2->SetLineWidth(1);
	t2->SetNDC();
	t2->Draw();	

	temp_canvEvt_wLog->Print( thePDFFileName.c_str() );	
      }
    }
    
    else if(drawDataMCRatios)
      {	std::string inHistName="hWeightedVz";
	
	std::cout<<" drawing DataMC Ratio hist "<<inHistName<<std::endl<<std::endl;
	
	long double theLumi;
	if(doEventCounts)theLumi=effIntgrtdLumi_vz;
	else theLumi=intgrtdLumi;

	std::string h_Title   ="EvtQA, Data/MC";
	std::string h_XAx_Title="v_{z}^{evt} (cm)", h_YAx_Title="Data/MC";
	
	TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
	theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(0) );
	theDataEvtQAHist->Scale( 1./theLumi );
	
	theDataEvtQAHist->SetTitle (    h_Title.c_str() );
	theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
	theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	theDataEvtQAHist->SetMarkerStyle(kFullCircle);
	theDataEvtQAHist->SetMarkerSize(1.8);
	theDataEvtQAHist->SetMarkerColor( kBlack);
	theDataEvtQAHist->SetLineColor( theRatioLineColor );
	theDataEvtQAHist->SetAxisRange(0.,1.5,"Y");

	//for ratio w/ 1./vzweights
	TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
	theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(0) );
	theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );
	
	TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
	theRatio->Divide(theMCEvtQAHist);
	theRatio->Draw();
	
	//theEvtQALeg->Draw();	
	
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

      if( doBasicJetQAOnly && j>4 )continue;
      if( !drawDijetPlots  && j>20 )continue;

      std::cout<<std::endl<<" var ="<<var[j]<<std::endl;
      
      for(int i=0;i<2;i++){

	if( i==1 && (!doJetIDPlots||!didJetIDCut) ) {
	  if(j==0)std::cout<<"drawing non-JetID Plots ONLY for JetQA!"<<std::endl;
	  continue; }

	else if( i==0 && (doJetIDPlots&&didJetIDCut) ) {
	  if(j==0)std::cout<<"drawing JetID Plots ONLY for JetQA!"<<std::endl;
	  continue; }
	
	long double theLumi;
	if(i==0)theLumi=effIntgrtdLumi_kmat;
	else theLumi=effIntgrtdLumi_JetID;

	std::string inHistName="hJetQA_"+std::to_string(i)+"wJetID_"+var[j];

	//open the hists
	if(debugMode)std::cout<<" opening input Data "<<inHistName<<std::endl<<std::endl;

	TH1F* theDataJetQAHist= (TH1F*)
	  ( (TH1*)finData->Get(inHistName.c_str()) )
	  ->TH1::Rebin(var_xAx_reBin[j]);	
	theDataJetQAHist->Print("base"); 
	std::cout<<std::endl;
	
	theDataJetQAHist->Scale( 1./theDataJetQAHist->GetBinWidth(1) );
	theDataJetQAHist->Scale( 1./theLumi);
	theDataJetQAHist->Scale( 1./4.);//eta bin width scaling
	theDataJetQAHist->Print("base"); 
	std::cout<<std::endl;

	theDataJetQAHist->SetMarkerStyle(kFullCircle);
	theDataJetQAHist->SetMarkerSize(1.2);
	theDataJetQAHist->SetMarkerColor( theDataOverlayMarkerColor );
	theDataJetQAHist->SetLineColor( theDataOverlayLineColor );
	

	if(debugMode)std::cout<<" opening input MC "<<inHistName<<std::endl<<std::endl;

	TH1F* theMCJetQAHist= (TH1F*)
	  ( (TH1*)finMC->Get(inHistName.c_str()) )
	  ->TH1::Rebin(var_xAx_reBin[j]);	
	theMCJetQAHist->Scale( 1./theMCJetQAHist->GetBinWidth(1) );	  
	theMCJetQAHist->Scale(1000.);
	theMCJetQAHist->Scale(1./4.);//eta bin width scaling
	theMCJetQAHist->Print("base");
	std::cout<<std::endl;

	if( var[j]!="jtpt"&&var[j]!="rawpt"&&
	    var[j]!="leadJetPt"&&var[j]!="subleadJetPt" )  
	  theMCJetQAHist->Scale( theDataJetQAHist->Integral()/theMCJetQAHist->Integral() );	  

	std::cout<<std::endl;
 
	theMCJetQAHist->SetMarkerStyle(kFullCircle);
	theMCJetQAHist->SetMarkerSize(1.2);
	theMCJetQAHist->SetMarkerColor( theMCOverlayMarkerColor );
	theMCJetQAHist->SetLineColor( theMCOverlayLineColor );
	
	//assert(false);
	
	

	if(drawDataMCOverlays){  
	  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title=crossSectionYAxis;//"#sigma (#mub/bin)";
	  if(debugMode)std::cout<<" x-axis title= "<<h_XAx_Title<<std::endl<<std::endl;

	  std::string h_Title="JetQA";//, "+var[j];
	  if(i==1)h_Title+=", w/ JetIDCut";      
	  //else h_Title+=", no JetIDCut";
	  
	  theDataJetQAHist->SetTitle (    h_Title.c_str() );
	  theDataJetQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theDataJetQAHist->SetYTitle( h_YAx_Title.c_str() );
	  //theDataJetQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );

	  theMCJetQAHist->SetTitle (    h_Title.c_str() );
	  theMCJetQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theMCJetQAHist->SetYTitle( h_YAx_Title.c_str() );
	  //theMCJetQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );

	  float legx1=0.77, legx2=legx1+0.07;
	  float legy1=0.76, legy2=legy1+0.11;
	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	  theJetQALeg->AddEntry(theDataJetQAHist,"data","lp");
	  theJetQALeg->AddEntry(theMCJetQAHist,  "MC"  ,"lp");
	  
	  theDataJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 

	  // variable specific draw orders, canv, etc.
	  if(var[j]=="jtpt"||var[j]=="rawpt") { temp_canvJet->cd();
	    theMCJetQAHist->Draw(); 
	    theDataJetQAHist->Draw("same");  
	    theJetQALeg->Draw(); 

	    float t1Loc1=0.73, t1Loc2=0.68; 	    
	    TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	    t1->SetTextFont(63);
	    t1->SetTextColor(kBlack);
	    t1->SetTextSize(20);
	    t1->SetLineWidth(1);
	    t1->SetNDC();
	    t1->Draw();
	    
	    TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}>15 GeV, |#eta|<2.0");
	    t2->SetTextFont(63);
	    t2->SetTextColor(kBlack);
	    t2->SetTextSize(20);
	    t2->SetLineWidth(1);
	    t2->SetNDC();
	    t2->Draw();	
	    
	    temp_canvJet->Print( thePDFFileName.c_str() ); }

	  else if(var[j]=="jteta"||var[j]=="jtphi"){ temp_canvJet->cd();
	    theDataJetQAHist->Draw(); 
	    theMCJetQAHist->Draw("SAME"); 
	    theJetQALeg->Draw();

	    float t1Loc1=0.71, t1Loc2=0.74; 	    
	    TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	    t1->SetTextFont(63);
	    t1->SetTextColor(kBlack);
	    t1->SetTextSize(20);
	    t1->SetLineWidth(1);
	    t1->SetNDC();
	    t1->Draw();
	    
	    TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}>50 GeV   |#eta|<2.0");
	    t2->SetTextFont(63);
	    t2->SetTextColor(kBlack);
	    t2->SetTextSize(20);
	    t2->SetLineWidth(1);
	    t2->SetNDC();
	    t2->Draw();	
	    //separate dijet from jetqa, put dijet cuts here
	    temp_canvJet->Print( thePDFFileName.c_str() ); }

	  else if(var[j]=="xj"||var[j]=="Aj"||var[j]=="dphi"){ temp_canvJet_noLog->cd();
	    

	    theDataJetQAHist->Draw(); 
	    theMCJetQAHist->Draw("SAME"); 
	    theJetQALeg->Draw();

	    float t1Loc1=0.18, t1Loc2=0.72; 	    

	    if(var[j]=="Aj") t1Loc1=.7;

	    TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	    t1->SetTextFont(63);
	    t1->SetTextColor(kBlack);
	    t1->SetTextSize(20);
	    t1->SetLineWidth(1);
	    t1->SetNDC();
	    t1->Draw();
	    
	    TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}^{ld} > 60 GeV");
	    t2->SetTextFont(63);
	    t2->SetTextColor(kBlack);
	    t2->SetTextSize(20);
	    t2->SetLineWidth(1);
	    t2->SetNDC();
	    t2->Draw();	
	   
	    TLatex *t3=new TLatex((t1Loc1),(t1Loc2-.10),"p_{t}^{subld} > 30GeV");
	    t3->SetTextFont(63);
	    t3->SetTextColor(kBlack);
	    t3->SetTextSize(20);
	    t3->SetLineWidth(1);
	    t3->SetNDC();
	    t3->Draw();	

	    TLatex *t4=new TLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");	    
	    t4->SetTextFont(63);
	    t4->SetTextColor(kBlack);
	    t4->SetTextSize(20);
	    t4->SetLineWidth(1);
	    t4->SetNDC();
	    t4->Draw();	
	    //separate dijet from jetqa, put dijet cuts here
	    temp_canvJet_noLog->Print( thePDFFileName.c_str() ); }

	  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt"){ temp_canvJet->cd();

	    theDataJetQAHist->Draw(); 
	    theMCJetQAHist->Draw("SAME"); 
	    theJetQALeg->Draw();

	    float t1Loc1=0.18, t1Loc2=0.72; 	    

	    if(var[j]=="Aj"||
	       var[j]=="leadJetPt"||var[j]=="subleadJetPt") t1Loc1=.7;

	    TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	    t1->SetTextFont(63);
	    t1->SetTextColor(kBlack);
	    t1->SetTextSize(20);
	    t1->SetLineWidth(1);
	    t1->SetNDC();
	    t1->Draw();
	    
	    TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}^{ld} > 60 GeV");
	    t2->SetTextFont(63);
	    t2->SetTextColor(kBlack);
	    t2->SetTextSize(20);
	    t2->SetLineWidth(1);
	    t2->SetNDC();
	    t2->Draw();	
	   
	    TLatex *t3=new TLatex((t1Loc1),(t1Loc2-.10),"p_{t}^{subld} > 30GeV");
	    t3->SetTextFont(63);
	    t3->SetTextColor(kBlack);
	    t3->SetTextSize(20);
	    t3->SetLineWidth(1);
	    t3->SetNDC();
	    t3->Draw();	

	    TLatex *t4=new TLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");	    
	    t4->SetTextFont(63);
	    t4->SetTextColor(kBlack);
	    t4->SetTextSize(20);
	    t4->SetLineWidth(1);
	    t4->SetNDC();
	    t4->Draw();	
	    //separate dijet from jetqa, put dijet cuts here
	    temp_canvJet->Print( thePDFFileName.c_str() ); }
	  
	  //else if(var[j]=="dphi") { temp_canvJet->cd();
	  //  theDataJetQAHist->Draw(); 
	  //  theMCJetQAHist->Draw("SAME");  
	  //  theJetQALeg->Draw();
	  //  
	  //  float t1Loc1=0.73, t1Loc2=0.80; 	    
	  //  TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	  //  t1->SetTextFont(63);
	  //  t1->SetTextColor(kBlack);
	  //  t1->SetTextSize(20);
	  //  t1->SetLineWidth(1);
	  //  t1->SetNDC();
	  //  t1->Draw();
	  //  TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_t^1>60 GeV, p_t^2> 30GeV, #delta#phi<#frac{2}{3}#pi");
	  //  t2->SetTextFont(63);
	  //  t2->SetTextColor(kBlack);
	  //  t2->SetTextSize(20);
	  //  t2->SetLineWidth(1);
	  //  t2->SetNDC();
	  //  t2->Draw();	
	  //  temp_canvJet->Print( thePDFFileName.c_str() ); }		  

	  else { temp_canvJet->cd();
	    theDataJetQAHist->Draw(); 
	    theMCJetQAHist->Draw("SAME");  
	    theJetQALeg->Draw();

	    float t1Loc1=0.73, t1Loc2=0.80; 	    
	    TLatex *t1=new TLatex(t1Loc1,t1Loc2, (fullJetType+"Jets").c_str());
	    t1->SetTextFont(63);
	    t1->SetTextColor(kBlack);
	    t1->SetTextSize(20);
	    t1->SetLineWidth(1);
	    t1->SetNDC();
	    t1->Draw();
	    
	    TLatex *t2=new TLatex(t1Loc1, (t1Loc2-.05),"p_{t}>50 GeV, |#eta|<2.0");
	    t2->SetTextFont(63);
	    t2->SetTextColor(kBlack);
	    t2->SetTextSize(20);
	    t2->SetLineWidth(1);
	    t2->SetNDC();
	    t2->Draw();	
	    temp_canvJet->Print( thePDFFileName.c_str() ); 	  }		  	  

	  std::cout<<std::endl; 

	}//end draw overlays

	else if(drawDataMCRatios){ 
	  temp_canvJetRat->cd();	  
	  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title="Data/MC"    ;
	  std::string h_Title   ="JetQA, Data/MC";//, "+var[j];
	  if(i==1)h_Title+=", w/ JetIDCut";      
	  //else h_Title+=", no JetIDCut";
	  
	  TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
	  
	  theRatio->SetTitle (    h_Title.c_str() );
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle( h_YAx_Title.c_str() );
	  theRatio->SetMarkerColor( kRed );
	  theRatio->SetLineColor( theRatioLineColor );
	  
	  theRatio->SetAxisRange(0,2.,"Y");
	  theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");

	  if(var[j]=="dphi") 
	    theRatio->SetAxisRange(2.0,3.2,"X");
	  
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
	  //theJetQALeg->Draw();
	  
	  temp_canvJetRat->Print( thePDFFileName.c_str() );	
	  std::cout<<std::endl;	    	
	}//end draw ratios

      } }//end loop over jetQA Plots     
    temp_canvJet_noLog->Close();//clean up neatly
    temp_canvJetRat->Close();
    temp_canvJet->Close();
  }
  
  // jet trig plots----------------------
  if(!drawJetRapBinsPlot || drawDataMCRatios )
    std::cout<<std::endl<<"skipping jet rap bin plot..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;
    std::cout<<" drawing jet rap bin Plots..."<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing Jet rapBins plots..."<<std::endl;
    TCanvas *temp_canvJetRapBins = new TCanvas("tempJetRapBins", "tempJetRapBins", 1200, 600);
    temp_canvJetRapBins->cd();
    temp_canvJetRapBins->SetLogy(1);    

    for(int i=0; i<2; i++){       if(debugMode)std::cout<<std::endl<<"i= "<<i<<std::endl;

      if( i==1 && (!doJetIDPlots||!didJetIDCut) ) {
	std::cout<<"drawing non-JetID Plots ONLY for JetTrigQA!"<<std::endl;
	continue; }
      else if( i==0 && (doJetIDPlots&&didJetIDCut) ) {
	std::cout<<"drawing JetID Plots ONLY for JetTrigQA!"<<std::endl;
	continue; }
      
      long double theLumi;
      if(i==0)theLumi=effIntgrtdLumi_kmat;
      else theLumi=effIntgrtdLumi_JetID;

      TLegend* theJetRapHistLeg=new TLegend(0.75,0.60,0.88,0.86, NULL,"brNDC");
      //TLegend* theMCJetRapHistLeg=new TLegend(0.61,0.60,0.74,0.86, NULL,"brNDC");
      
      if(drawDataMCOverlays){ float power=std::pow(10., 0.);
	for(int rapbin=0; rapbin<(nbins_rap); ++rapbin){
	  if(rapbin==0)std::cout<<"drawing jet pt in rapBin plot...."<<std::endl;
	  
	  std::string inHistName="hJetSpectraRap";
          if(i==1)inHistName+="_wJetID";
          inHistName+="_bin_"+std::to_string(rapbin);
          if(debugMode)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;

	  //std::string inMCHistName="hJetSpectraRap";
          //if(i==1)inMCHistName+="_wJetID";
          //inMCHistName+="_bin"+std::to_string(rapbin);
          //if(debugMode)std::cout<<"inMCHistName="<<inMCHistName<<std::endl<<std::endl;

          //open the hists + do scaling
          TH1F* theJetSpectraRapHist= (TH1F*)
            ( (TH1*)finData->Get(inHistName.c_str()) )
	    ->TH1::Rebin(jetSpectraRapBinning);
          theJetSpectraRapHist->Scale( 1./theJetSpectraRapHist->GetBinWidth(1) );
          theJetSpectraRapHist->Scale( 1./theLumi);
          theJetSpectraRapHist->Scale( power );
	  
          theJetSpectraRapHist->SetMarkerStyle( 20 ) ;//theRapOverlayMarker[j]);
          theJetSpectraRapHist->SetMarkerSize(1.1);
          theJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
          theJetSpectraRapHist->SetLineColor( kBlack );//theRapOverlayLineColor[j]   );
	  //theJetSpectraRapHist->SetAxisRange(0.,800.,"X");


          ////open the hists + do scaling
          //TH1F* theMCJetSpectraRapHist= (TH1F*)
          //  ( (TH1*)finMC->Get(inMCHistName.c_str()) )
	  //  ->TH1::Rebin(jetSpectraRapBinning);
          //theMCJetSpectraRapHist->Scale( 1./theMCJetSpectraRapHist->GetBinWidth(1) );
          //theMCJetSpectraRapHist->Scale( power );
          //theMCJetSpectraRapHist->Scale( theJetSpectraRapHist->Integral()/theMCJetSpectraRapHist->Integral() );
          //theMCJetSpectraRapHist->SetMarkerStyle( 27 ) ;//theRapOverlayMarker[j]);
          //theMCJetSpectraRapHist->SetMarkerSize(1.1);
          //theMCJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
          //theMCJetSpectraRapHist->SetLineColor( -1 );//theRapOverlayLineColor[j]   );
	  ////theJetSpectraRapHist->SetAxisRange(0.,800.,"X");	 


	  std::string rapHistLegDesc;
	  std::stringstream stream1, stream2;
	  stream1.precision(1); stream1 << std::fixed << rapbins[rapbin];
	  stream2.precision(1); stream2 << std::fixed << rapbins[rapbin+1];
	  rapHistLegDesc=stream1.str()+"<|y|<"+stream2.str()+" x10^{"+std::to_string(rapbin)+"}";
	  
          theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
          //theMCJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
	  
          if(rapbin==0){
	    std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title=crossSectionYAxis;//"#sigma (#mub/bin)"    ;
	    std::string h_Title   ="Jet Pt in |y| bins";//, DATA ONLY";
            if(i==1)h_Title+=", w/ JetIDCut";
            
            theJetSpectraRapHist->SetAxisRange(std::pow(10.,-8.),std::pow(10.,5.),"Y");
            theJetSpectraRapHist->SetTitle (    h_Title.c_str() );
            theJetSpectraRapHist->SetXTitle( h_XAx_Title.c_str() );
            theJetSpectraRapHist->SetYTitle( h_YAx_Title.c_str() );
            theJetSpectraRapHist->Draw();     	  
	    //            theMCJetSpectraRapHist->Draw("SAME");     	  
	    theJetRapHistLeg->Draw(); }
	    //theMCJetRapHistLeg->Draw("SAME"); 	  }
	  else { 
	    theJetSpectraRapHist->Draw("SAME");}
	    //theMCJetSpectraRapHist->Draw("SAME");   }
	  power*=10.;
	}//end rapbin hist loop
      }//end draw mcdata overlays

      float t1Loc1=0.50, t1Loc2=0.86; 	    
      TLatex *t1=new TLatex(t1Loc1,t1Loc2,"pp 2015 promptReco, #sqrt{s} = 5.02 TeV");
      t1->SetTextFont(63);
      t1->SetTextColor(kBlack);
      t1->SetTextSize(19);
      t1->SetLineWidth(1);
      t1->SetNDC();
      t1->Draw();
      
      TLatex *pp1=new TLatex(t1Loc1,(t1Loc2-0.05),"Jet80+LowerJets, L_{int} = 25.8 pb^{-1}");
      pp1->SetTextFont(63);
      pp1->SetTextColor(kBlack);
      pp1->SetTextSize(19);
      pp1->SetLineWidth(1);
      pp1->SetNDC();
      pp1->Draw();
      
      TLatex *pp2=new TLatex(t1Loc1,(t1Loc2-0.10),(fullJetType+"Jets, p_{t}>15 GeV").c_str());
      pp2->SetTextFont(63);
      pp2->SetTextColor(kBlack);
      pp2->SetTextSize(19);
      pp2->SetLineWidth(1);
      pp2->SetNDC();
      pp2->Draw();
      
      //TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}>15 GeV, |#eta|<2.0");
      //t2->SetTextFont(63);
      //t2->SetTextColor(kBlack);
      //t2->SetTextSize(19);
      //t2->SetLineWidth(1);
      //t2->SetNDC();
      //t2->Draw();	

      temp_canvJetRapBins->Print( thePDFFileName.c_str() );
    }//end JetID loop
    
  }//end jet pt rap bins plot
  
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
      

      if( i==1 && (!doJetIDPlots||!didJetIDCut) ) {
	std::cout<<"drawing non-JetID Plots ONLY for JetTrigQA!"<<std::endl;
	continue; }

      else if( i==0 && (doJetIDPlots&&didJetIDCut) ) {
	std::cout<<"drawing JetID Plots ONLY for JetTrigQA!"<<std::endl;
	continue; }
      

      long double theLumi;
      if(i==0)theLumi=effIntgrtdLumi_kmat;
      else theLumi=effIntgrtdLumi_JetID;
      
      TLegend* JetTrigLegend=new TLegend(0.75,0.60,0.88,0.86, NULL,"brNDC");
      TLegend* JetTrigRatioLegend=new TLegend(0.12,0.59,0.30,0.87, NULL,"brNDC");
      
      if(drawTrigDataOverlays){   temp_canvJetTrig->cd();
	if(debugMode)std::cout<<"drawing dataTrig QA Overlay Plots"<<std::endl;

	for(int j=0; j<N_trigs; j++){ 
	  
	  std::string inHistName="hpp_"+HLTName[j];
	  if(i==1)inHistName+="_wJetID";
	  inHistName+="_"+radius+etaWidth;
	  
	  if(debugMode)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
	  
	  //open the hists + do scaling 
	  TH1F* theJetTrigQAHist= (TH1F*)
	    ( (TH1*)finData->Get(inHistName.c_str()) )
	    ->TH1::Rebin(jetTrigQABinning);
	  
	  //theJetTrigQAHist->Print("base");
	  theJetTrigQAHist->Scale( 1./theJetTrigQAHist->GetBinWidth(1) );
	  theJetTrigQAHist->Scale( 1./theLumi);
	  theJetTrigQAHist->Scale( 1./4.);//eta bin width
	  
	  theJetTrigQAHist->SetMarkerStyle(theTrigOverlayMarker[j]);
	  theJetTrigQAHist->SetMarkerSize(1.2);
	  theJetTrigQAHist->SetMarkerColor( theTrigOverlayMarkerColor[j] );
	  theJetTrigQAHist->SetLineColor(   theTrigOverlayLineColor[j]   );
	  
	  JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTName_Leg[j]).c_str(),"lp");
	  
	  if(j==0){
	    std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title=crossSectionYAxis;//"#sigma (#mub/bin)"    ; 
	    std::string h_Title   ="JetTrigQA";//, DATA ONLY";
	    if(i==1)h_Title+=", w/ JetIDCut";      
	    //else h_Title+=", no JetIDCut";
	    theJetTrigQAHist->SetAxisRange(0.,800.,"X");
	    theJetTrigQAHist->SetTitle (    h_Title.c_str() );
	    theJetTrigQAHist->SetXTitle( h_XAx_Title.c_str() );
	    theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
	    //theJetTrigQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
	    theJetTrigQAHist->Draw(); 	    }
	  else theJetTrigQAHist->Draw("SAME");   	    

	}//end loop over trigs
	JetTrigLegend->Draw();
	temp_canvJetTrig->Print(thePDFFileName.c_str());
	


	
      }//end DRAW TRIG OVERLAY      

      else if(drawTrigDataRatios){		temp_canvJetRatio->cd();
	if(debugMode)std::cout<<"drawing dataTrig QA Overlay Plots"<<std::endl;
	
	std::string comboPlotName="hpp_"+HLTName[N_trigs-1];
	if(i==1)comboPlotName+="_wJetID";
	comboPlotName+="_"+radius+etaWidth;
	
	const int jetTrigQARatioBinning=jetTrigQABinning/2;
	
	TH1F* theComboHist=(TH1F*)
	  ( (TH1*)finData->Get(comboPlotName.c_str()) )
	  ->TH1::Rebin(jetTrigQARatioBinning);
	
	for(int j=0; j<(N_trigs-1); j++){ 
	  
	  std::string inHistName="hpp_"+HLTName[j];
	  if(i==1)inHistName+="_wJetID";
	  inHistName+="_"+radius+etaWidth;
	  
	  TH1F* comboHistClone=(TH1F*)theComboHist->Clone("comboHistClone4Div");
	  
	  TH1F* theRatio=(TH1F*)
	    ( (TH1*)finData->Get(inHistName.c_str()) )
	    ->TH1::Rebin(jetTrigQARatioBinning);
	  
	  theRatio->SetMarkerStyle(theTrigOverlayMarker[j]);
	  theRatio->SetMarkerSize(1.2);
	  theRatio->SetMarkerColor( theTrigOverlayMarkerColor[j] );
	  theRatio->SetLineColor( theTrigOverlayLineColor[j] );
	  
	  theRatio->Divide(comboHistClone);
	  
	  JetTrigRatioLegend->AddEntry(theRatio,(HLTName_Leg[j]+"/TrigComb.").c_str(),"lp");
	  
	  if(j==0){
	    std::string h_Title   ="JetTrigQA Ratios";
	    if(i==1)    h_Title+=", w/ JetIDCut";      
	    std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title="Trig/TrigCombo"    ; 

	    theRatio->SetAxisRange(0.,500.,"X");
	    theRatio->SetAxisRange(0.,2.,"Y");
	    theRatio->SetTitle (    h_Title.c_str() );
	    theRatio->SetXTitle( h_XAx_Title.c_str() );
	    theRatio->SetYTitle( h_YAx_Title.c_str() );
	    theRatio->Draw(); }

	  else theRatio->Draw("SAME");
	}
	
	JetTrigRatioLegend->Draw();
	temp_canvJetRatio->Print(thePDFFileName.c_str());
      }//end DRAW TRIG Ratios
      
    }//end jet ID loop, i=0,1    
    temp_canvJetTrig->Close();        
    temp_canvJetRatio->Close();    
  }//end else statement for het trig plots
  
  
  

    // MCEffQAPlots plots----------------------
  if(!drawMCEffProfiles)std::cout<<std::endl<<"skipping MC Eff. QA plots..."<<std::endl<<std::endl;
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
	  if(i==1)inHistName+="_wJetID";
	  inHistName+="_"+radius+etaWidth;
	  
	  std::cout<<"opening TH2F "<<inHistName << std::endl;
	  TH2F* the2DMCEffQAHist= (TH2F*)finMC->Get( inHistName.c_str() );
	  if(!the2DMCEffQAHist) {std::cout<<"no MCEff plot, continuing..."<<std::endl; continue;}
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
  
  
  
// for DEBUG ONLY
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1200, 600);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
    

  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  finMC->Close();
  finData->Close();  
  //assert(false);
  
  return 0;
}// end printplots
				     
				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=5 ) {//no input arguments, error
    std::cout<<"do ./printQAPlots.exe <ppDataCondorDir> <ppMCCondorDir> <outputNameBase> <OverlayFlag>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  //if(argc==1)rStatus=printQAPlots();
  if(argc==5) rStatus=printQAPlots( (const std::string) argv[1], (const std::string) argv[2],
				    (const std::string) argv[3], (bool) std::atoi(argv[4])    );
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
}



//// I/O ------------------------
////input
//const std::string input_ppData_condorDir=readForests_dir+
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
//const std::string input_ppMC_condorDir=readForests_dir+
////"ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4CaloJets_09-17-16__newJetID/";
//  "ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_09-16-16__newJetID/";  
////"ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak3PFJets_09-18-16__newJetID/";
//
//const std::string input_ppMC_Filename=input_ppMC_condorDir+
////"Py8_CUETP8M1_QCDjetAllPtBins_ak4Calo-allFiles.root";
//  "Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root";
////"Py8_CUETP8M1_QCDjetAllPtBins_ak3PF-allFiles.root";

//const std::string output_PDFname_base=outputDir+
////"printQAplots_9.18.16_newJetID_ak4CaloJets_HLT.ak4PFJets_QAPlots_";
//  "printQAplots_9.18.16_newJetID_ak4PFJets_HLT.ak4PFJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak3PFJets_HLT.ak4PFJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak4CaloJets_HLT.ak4CaloJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak4PFJets_HLT.ak4CaloJets_QAPlots_";
////"printQAplots_9.18.16_newJetID_ak3PFJets_HLT.ak4CaloJets_QAPlots_";
  
 

//// trigger closure/turn on plot
//temp_canvJetRatio->cd();
//for(int trigTH2=0;trigTH2<3;trigTH2++){
//
//  std::string inHistName;//="hTrgQA_";//+HLTName[j];
//  //if(i==1)inHistName+="_wJetID";
//  //inHistName+="_"+radius+etaWidth;	   
//
//  if(trigTH2==0)inHistName="hTrgQA_recpt_rawpt";
//  else if(trigTH2==1)inHistName="hTrgQA_recpt_trgpt";
//  else if(trigTH2==2)inHistName="hTrgQA_rawpt_trgpt";
//  else break;
//
//  
//  std::string h_Title   ="JetTrigQA";//, DATA ONLY";
//  if(i==1)h_Title+=", w/ JetIDCut";      
//
//  std::string h_XAx_Title, h_YAx_Title;
//  if(trigTH2==0){      h_XAx_Title=" raw JetPt (GeV)";  h_YAx_Title="reco JetPt (GeV)"; }
//  else if(trigTH2==1){ h_XAx_Title="trig JetPt (GeV)";  h_YAx_Title="reco JetPt (GeV)"; }
//  else if(trigTH2==2){ h_XAx_Title=" raw JetPt (GeV)";  h_YAx_Title="trig JetPt (GeV)"; }
//  else break;
//
//  if(debugMode)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
//
//  //open the hists + do scaling
//  TH2F* theJetTrigQAHist=(TH2F*)((TH2*)finData->Get(inHistName.c_str()));//->TH1::Rebin(jetTrigQABinning);	  
//
//  float jetTrigTh2_ptLow=10,jetTrigTh2_ptHigh=300;
//  theJetTrigQAHist->SetAxisRange(jetTrigTh2_ptLow,jetTrigTh2_ptHigh,"X");
//  theJetTrigQAHist->SetAxisRange(jetTrigTh2_ptLow,jetTrigTh2_ptHigh,"Y");
//
//  theJetTrigQAHist->SetTitle (    h_Title.c_str() );
//  theJetTrigQAHist->SetXTitle( h_XAx_Title.c_str() );
//  theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
//
//  theJetTrigQAHist->Print("base");
//  //theJetTrigQAHist->Scale( 1./theJetTrigQAHist->GetBinWidth(1) );
//  //theJetTrigQAHist->Scale( 1./theLumi);
//
//  theJetTrigQAHist->Draw("COLZ");
//  temp_canvJetRatio->Print(thePDFFileName.c_str());
//}//end TH2Loop for trigpt/rawpt/recpt plots




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


