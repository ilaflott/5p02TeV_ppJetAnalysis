#include "printPlots_combine.h"


// code/job switches ------------------------

//options
const bool debugMode=true, doEventCounts=true, doJetIDPlots=false;

//draw switches
const bool drawEvtQAPlots=true;
const bool drawJetQAPlots=true;
const bool drawJetConstituentPlots=true, drawDijetPlots=true;
const bool drawJetTrigQAPlots=true, drawJetRapBinsPlot=false;

const bool drawDataMCOverlaysInput=false; 		 //this should always be false, but needs to be cleaned up later.
const bool drawDataMCOverlays   = drawDataMCOverlaysInput;
const bool drawDataMCRatios     = !drawDataMCOverlays;
const bool drawTrigDataOverlays = drawDataMCOverlays && drawJetTrigQAPlots;
const bool drawTrigDataRatios   = drawDataMCRatios   && drawJetTrigQAPlots;
  
 


// hist painting ------------------------

// line colors
const int theDataOverlayLineColor=1, theMCOverlayLineColor=1, altOverlayLineColor=3; 
const int theRatioLineColor=1,altRatioLineColor1=8, altRatioLineColor2=7;
const int theTrigOverlayLineColor[]  ={  1,  1,  1,  1,  1,  1 };

// marker colors
const int theDataOverlayMarkerColor=2, theMCOverlayMarkerColor=4,theRatioMarkerColor=9;
const int theTrigCombMarkerColor=1, altTrigCombMarkerColor=12;
const int theTrigOverlayMarkerColor[]={  2,  3,  6,  7,  1,  4 };
const int theRapOverlayMarkerColor[] ={  2,  3,  6,  7,  1,  4 , 8, 9};

// markers
const int theTrigOverlayMarker[]     ={ 20, 20, 20, 20, 20, 32 };



// misc ------------------------
const std::string crossSectionYAxis="#sigma (#mub)";


// the macro ------------------------
int printPlots_jetPlots(const std::string input_ppData_condorDir , const std::string input_ppMC_condorDir , 
		 const std::string output_PDFname_base){
  
  //// global style settings
  //gStyle->SetTitleSize(0.035,"t");
  //gStyle->SetTitleOffset(1.0,"t");
  //
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  //
  gStyle->SetPadTickY(1);   
  gStyle->SetTickLength(0.015,"y");
  gStyle->SetLabelSize(0.02,"y");
  gStyle->SetTitleOffset(0.89,"y");
  gStyle->SetTitleSize(0.03,"y");
  //
  gStyle->SetPadTickX(1);
  //gStyle->SetTickLength(0.018,"x");
  gStyle->SetLabelSize(0.02,"x");
  gStyle->SetTitleOffset(0.94,"x");
  gStyle->SetTitleSize(.03,"x");
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();





  // MAKE STRINGS
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
  TFile *finMC=NULL, *finData=NULL;  
  std::cout<<std::endl<<"printing QA Plots, now opening input files!!"<<std::endl<<std::endl;
  
  //open the input files
  std::cout<<" now opening ppData: "<<std::endl<<input_ppData_Filename<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_ppData_condorDir<<std::endl<<std::endl;
  finData = new TFile(ppData_fullFilename.c_str());      
  
  std::cout<<" now opening ppMC: "<<std::endl<<input_ppMC_Filename<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_ppMC_condorDir<<std::endl<<std::endl;
  finMC = new TFile(ppMC_fullFilename.c_str());
  

  











  // GET OUTPUT PDF FILE READY
  std::string thePDFFileName=outputDir+fullJetType;;
  thePDFFileName+="_"+output_PDFname_base+"_combine.pdf";
  
  std::string open_thePDFFileName, close_thePDFFileName;
  open_thePDFFileName=thePDFFileName+"[";  close_thePDFFileName=thePDFFileName+"]";
  
  // just for opening the pdf
  TCanvas *temp_canvOpen = new TCanvas("temp", "temp", 1000, 800);

  temp_canvOpen->Print( open_thePDFFileName.c_str() );  
  //temp_canvOpen->Print( thePDFFileName.c_str() );  
  //temp_canvOpen->UseCurrentStyle();
  temp_canvOpen->Close();
  
  

  // Do Lumi Eff + evt + jet counts ----------------------
  long double effIntgrtdLumi_vz=1., LumiEff_vz=1.;// eff-> effective, Eff->efficiency
  //long double effIntgrtdLumi_kmat=1., LumiEff_kmat=1.;
  //long double effIntgrtdLumi_JetID=1., LumiEff_JetID=1.;
  //bool didJetIDCut=false;
  
  if(!doEventCounts)std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
  else{ 
    
    std::cout<<" grabbing raw event counts "<<std::endl;    
    for(int file=0;file<2;file++){ 
      
      TFile * fin=NULL;
      if(file==0)fin=finData;
      else fin=finMC;
      
      if(!fin){	std::cout<<"input file (int file="<<file<<") not found, cannot look at event counts"<<std::endl; continue; }
      
      
      TH1F *h_NEvents         = (TH1F*)fin->Get("NEvents");
      TH1F *h_NEvents_read    = (TH1F*)fin->Get("NEvents_read");      
      TH1F *h_NEvents_skimCut = (TH1F*)fin->Get("NEvents_skimCut");
      TH1F *h_NEvents_vzCut   = (TH1F*)fin->Get("NEvents_vzCut");

      TH1F *h_NEvents_skipped = NULL;
      TH1F *h_NEvents_trigd=NULL;
      TH1F *h_NEvents_jet100=NULL,*h_NEvents_jet80=NULL,*h_NEvents_jet60=NULL,*h_NEvents_jet40=NULL;
      if(file==0)
	{ h_NEvents_skipped=(TH1F*)fin->Get("NEvents_skipped");      
	  h_NEvents_trigd=(TH1F*)fin->Get("NEvents_trigd");      
	  h_NEvents_jet40=(TH1F*)fin->Get("NEvents_jet40Trigd");      
	  h_NEvents_jet60=(TH1F*)fin->Get("NEvents_jet60Trigd");      
	  h_NEvents_jet80=(TH1F*)fin->Get("NEvents_jet80Trigd");      
	  h_NEvents_jet100=(TH1F*)fin->Get("NEvents_jet100Trigd");      
	}
      
      if(file==0)std::cout<<std::endl<<"/// Data Evt-Count Summary ///"<<std::endl<<std::endl;      
      else std::cout<<std::endl<<"/// MC Evt-Count Summary ///"<<std::endl<<std::endl;      
      
      std::cout<<"Total Num of Events in file(s) opened       = " <<
	h_NEvents->GetEntries()<<std::endl;
      if(file==0)
	std::cout<<"Total Num of Events skipped from those file(s) = " <<
	  h_NEvents_skipped->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read from those file(s) = " <<
	h_NEvents_read->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read passing skimCuts   = " <<
	h_NEvents_skimCut->GetEntries()<<std::endl;
      std::cout<<"Total Num of Events read passing vzCuts and skimCuts    = " <<
	h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;

      if(file==0){
	std::cout<<"Total Num of Events trigd in those file(s) = " <<
	  h_NEvents_trigd->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events jet40 in those file(s) = " <<
	  h_NEvents_jet40->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events jet60 in those file(s) = " <<
	  h_NEvents_jet60->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events jet80 in those file(s) = " <<
	  h_NEvents_jet80->GetEntries()<<std::endl;
	std::cout<<"Total Num of Events jet100 in those file(s) = " <<
	  h_NEvents_jet100->GetEntries()<<std::endl;}

      
      /*
	TH1F *h_NEvents_withJets_kmatCut  = (TH1F*) fin->Get("NEvents_withJets_kmatCut");	
	TH1F *h_NEvents_withJets_JetIDCut = (TH1F*) fin->Get("NEvents_withJets_JetIDCut");
	
	didJetIDCut= ((bool)h_NEvents_withJets_JetIDCut);// && doJetIDPlots;
	if(doJetIDPlots && !didJetIDCut){
	std::cout<<std::endl<<"cannot make hists w/ jetID cut, jetID hists not there! exiting..."<<std::endl<<std::endl;
	assert(didJetIDCut);}
	
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
      */
      
      if(file==0){ //if data, do luminosity/integrated luminosity
	
	std::cout<<std::endl<<"dataset integrated Luminosity (microbarns) ="<<intgrtdLumi<<std::endl;
	
	LumiEff_vz = h_NEvents_vzCut->GetEntries()/h_NEvents_read->GetEntries();
	//LumiEff_kmat = h_NEvents_withJets_kmatCut->GetEntries()/h_NEvents->GetEntries();
	//if(didJetIDCut) LumiEff_JetID = h_NEvents_withJets_JetIDCut->GetEntries()/h_NEvents->GetEntries();
	
	effIntgrtdLumi_vz=intgrtdLumi*LumiEff_vz;
	//effIntgrtdLumi_kmat=intgrtdLumi*LumiEff_kmat;
	//if(didJetIDCut)effIntgrtdLumi_JetID=intgrtdLumi*LumiEff_JetID;
	
	std::cout<<std::endl<<"lumi efficiency, vz cuts="<<LumiEff_vz<<std::endl;
	std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_vz<<std::endl;
	
	//std::cout<<std::endl<<"lumi efficiency, kmatCuts="<<LumiEff_kmat<<std::endl;
	//std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_kmat<<std::endl;
	
	//if(didJetIDCut) std::cout<<std::endl<<"lumi efficiency, JetID="<<LumiEff_JetID<<std::endl;
	//if(didJetIDCut) std::cout<<"effective ingrtdLumi for dataset+cuts="<<effIntgrtdLumi_JetID<<std::endl;            
      }
      
    }// end data/MC file loop
  }// end evt counts
  
  








  // effective integrated luminosity ----------------------
  long double theLumi;
  if(doEventCounts)theLumi=effIntgrtdLumi_vz;
  else theLumi=intgrtdLumi;
  


  // evt plots ----------------------
  if(!drawEvtQAPlots)std::cout<<std::endl<<"skipping evt QA plots..."<<std::endl<<std::endl;
  //else if(doJetIDPlots)std::cout<<std::endl<<"doing JetID plots, skipping evt QA plots..."<<std::endl<<std::endl;

  else { 
    std::cout<<std::endl<<" drawing evt QA Plots now! "<<std::endl<<std::endl;
    
    if(debugMode)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
    TCanvas *temp_canvEvt_wLog = new TCanvas("tempEvt_wLog", "temp Evt Canv", 1000, 800);//This is used for the Jet QA Evt plots with log plots.
    temp_canvEvt_wLog->cd(); 
    //temp_canvEvt_wLog->UseCurrentStyle();
    temp_canvEvt_wLog->SetLogy(1);
    
    
    if(debugMode)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
    TCanvas *temp_canvEvt = new TCanvas("tempEvt", "temp Evt Canv withLog", 1000, 800); //This is for the Jet QA Evt plots without logs
    temp_canvEvt->cd();

    //temp_canvEvt->UseCurrentStyle();

    //if(drawDataMCOverlays){    
    //if(debugMode)std::cout<<"drawing Data/MC overlays for Evt QA..."<<std::endl; //makes Evt QA plots
    
    {// data, MC unweighted
      temp_canvEvt->cd();
      //	temp_canvEvt->Divide(1,2);
      
      TPad *evtpad1 = new TPad("vz1evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
      TPad *evtpad2 = new TPad("vz1evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
      
      //temp_canvEvt->cd(1);
      evtpad1->SetBottomMargin(0); 
      evtpad1->Draw();
      evtpad1->cd();
      // open hist
      
      std::string inHistName  ="hVz";
      std::cout<<" drawing "<<inHistName<<std::endl<<std::endl;
      TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inHistName.c_str() );
      
      theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
      theDataEvtQAHist->Scale( 1./theLumi );
      theDataEvtQAHist->SetMarkerStyle(kCircle);
      theDataEvtQAHist->SetMarkerSize(0.99);
      theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
      theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );
      

      std::string inMCHistName  ="hVz";
      TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
      
      theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
      theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
      theMCEvtQAHist->SetMarkerStyle(kMultiply);
      theMCEvtQAHist->SetMarkerSize(0.99);
      theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
      
      // title+axes
      std::string h_Title     ="EvtQA";
      std::string h_XAx_Title ="v_{z}^{evt} (cm)";
      std::string h_YAx_Title =crossSectionYAxis;
      theDataEvtQAHist->SetTitle (    h_Title.c_str() );
      theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
      theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theDataEvtQAHist->SetAxisRange(0.0,2.0, "Y");  
      
      theDataEvtQAHist->Draw();
      theMCEvtQAHist->Draw("same");	
          
      
      // legend
      float legx1=0.78, legx2=legx1+0.11;
      float legy1=0.74, legy2=legy1+0.09;
      TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
      
      theEvtQALeg->AddEntry(theDataEvtQAHist, "data" ,"lp");
      theEvtQALeg->AddEntry(theMCEvtQAHist,   "MC"       ,"lp");
      theEvtQALeg->Draw();
      
      // text
      float pp1x=0.14,pp1y=0.82;
      TLatex *pp1=new TLatex(pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets");
      pp1->SetTextFont(63);
      pp1->SetTextColor(kBlack);
      pp1->SetTextSize(15);
      pp1->SetLineWidth(1);
      pp1->SetNDC();
      pp1->Draw();
      TLatex *pp2=new TLatex(pp1x,(pp1y-0.04),"L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV");
      pp2->SetTextFont(63);
      pp2->SetTextColor(kBlack);
      pp2->SetTextSize(15);
      pp2->SetLineWidth(1);
      pp2->SetNDC();
      pp2->Draw();
      //TLatex *pp3=new TLatex(pp1x,(pp1y-0.10),(fullJetType+"Jets").c_str() );
      //pp3->SetTextFont(63);
      //pp3->SetTextColor(kBlack);
      //pp3->SetTextSize(17);
      //pp3->SetLineWidth(1);
      //pp3->SetNDC();
      //pp3->Draw();
      TLatex *pp4=new TLatex(pp1x,(pp1y-0.08),"Py8 Tune CUETP8M1" );
      pp4->SetTextFont(63);
      pp4->SetTextColor(kBlack);
      pp4->SetTextSize(15);
      pp4->SetLineWidth(1);
      pp4->SetNDC();
      pp4->Draw();


      ////tlatex for cut
      //TLatex *pp5=new TLatex(pp1x,(pp1y-0.11),EvtPtCut.c_str());
      //pp5->SetTextFont(63);
      //pp5->SetTextColor(kBlack);
      //pp5->SetTextSize(15);
      //pp5->SetLineWidth(1);
      //pp5->SetNDC();
      //pp5->Draw();
      
      
      
      
      //ratios go here
      temp_canvEvt->cd(); //change back to main canvas
      evtpad2->Draw();
      evtpad2->cd();
      //hist for ratio; use clone
      
      TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
      evtpad2->SetTopMargin(0);
      evtpad2->SetBottomMargin(0.3);
      
      theRatio->GetYaxis()->SetTitleSize(15);
      theRatio->GetYaxis()->SetTitleFont(43);
      theRatio->GetYaxis()->SetTitleOffset(2);
      theRatio->GetYaxis()->SetLabelFont(43); 
      theRatio->GetYaxis()->SetLabelSize(13);
      theRatio->GetXaxis()->SetTitleSize(20);
      theRatio->GetXaxis()->SetTitleFont(43);
      theRatio->GetXaxis()->SetTitleOffset(4.);
      theRatio->GetXaxis()->SetLabelFont(43); 
      theRatio->GetXaxis()->SetLabelSize(15);
      
      theRatio->SetTitle("");
      theRatio->SetXTitle( h_XAx_Title.c_str() );
      theRatio->SetYTitle("Data/MC");
      theRatio->Divide(theMCEvtQAHist);
      theRatio->Draw("h");
      
	
      // print to PDF file
      temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.

    }// end vz




    {// MC pthat weighted, data trigger weighted
      temp_canvEvt->cd();
      //	temp_canvEvt->Divide(1,2);
      
      TPad *evtpad1 = new TPad("vz2evtpad1", "Overlay Pad",0.0,0.3,1.0,1.0);
      TPad *evtpad2 = new TPad("vz2evtpad2", "Ratio Pad",0.0,0.05,1.0,0.3);

      
      //temp_canvEvt->cd(1);
      evtpad1->SetBottomMargin(0); 
      evtpad1->Draw();
      evtpad1->cd();
      // open hist
      
      std::string inDataHistName  ="hWeightedVz";
      std::cout<<" drawing "<<inDataHistName<<std::endl<<std::endl;
      TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
      
      theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
      theDataEvtQAHist->Scale( 1./theLumi );
      theDataEvtQAHist->SetMarkerStyle(kCircle);
      theDataEvtQAHist->SetMarkerSize(0.99);
      theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
      theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );
      
      //get jet pt cut
      //std::string jetPTcut ="hJetPtCut";
      //TH1F* hJetPTCut= (TH1F*)finData->Get( jetPTcut.c_str() );
      //hJetPTCut->GetMean();
      //std::string EvtPt = ""; //std::to_string((int)hJetPTCut->GetMean()).c_str();
      //std::string EvtPtCut = "p_{T} > " + EvtPt + " GeV";
      //std::cout<<"EvtPtCut = "<<EvtPtCut.c_str()<<std::endl;
      
      std::string inMCHistName  ="hpthatWeightedVz";
      TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
      
      theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
      theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
      theMCEvtQAHist->SetMarkerStyle(kMultiply);
      theMCEvtQAHist->SetMarkerSize(0.99);
      theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
      
      // title+axes
      std::string h_Title     ="EvtQA";
      std::string h_XAx_Title ="v_{z}^{evt} (cm)";
      std::string h_YAx_Title =crossSectionYAxis;
      theDataEvtQAHist->SetTitle (    h_Title.c_str() );
      theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
      theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theDataEvtQAHist->SetAxisRange(0.0,2.0, "Y");  
      
      theDataEvtQAHist->Draw();	
      theMCEvtQAHist->Draw("same");
      
      
      // legend
      float legx1=0.78, legx2=legx1+0.11;
      float legy1=0.74, legy2=legy1+0.09;
      TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
      
      theEvtQALeg->AddEntry(theDataEvtQAHist, "data" ,"lp");
      theEvtQALeg->AddEntry(theMCEvtQAHist,   "MC"       ,"lp");
      theEvtQALeg->Draw();
      
      // text
      float pp1x=0.14,pp1y=0.82;
      TLatex *pp1=new TLatex(pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets");
      pp1->SetTextFont(63);
      pp1->SetTextColor(kBlack);
      pp1->SetTextSize(15);
      pp1->SetLineWidth(1);
      pp1->SetNDC();
      pp1->Draw();
      TLatex *pp2=new TLatex(pp1x,(pp1y-0.04),"L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV");
      pp2->SetTextFont(63);
      pp2->SetTextColor(kBlack);
      pp2->SetTextSize(15);
      pp2->SetLineWidth(1);
      pp2->SetNDC();
      pp2->Draw();
      //TLatex *pp3=new TLatex(pp1x,(pp1y-0.10),(fullJetType+"Jets").c_str() );
      //pp3->SetTextFont(63);
      //pp3->SetTextColor(kBlack);
      //pp3->SetTextSize(17);
      //pp3->SetLineWidth(1);
      //pp3->SetNDC();
      //pp3->Draw();
      TLatex *pp4=new TLatex(pp1x,(pp1y-0.08),"Py8 Tune CUETP8M1" );
      pp4->SetTextFont(63);
      pp4->SetTextColor(kBlack);
      pp4->SetTextSize(15);
      pp4->SetLineWidth(1);
      pp4->SetNDC();
      pp4->Draw();
    
    
      ////tlatex for cut
      //TLatex *pp5=new TLatex(pp1x,(pp1y-0.11),EvtPtCut.c_str());
      //pp5->SetTextFont(63);
      //pp5->SetTextColor(kBlack);
      //pp5->SetTextSize(15);
      //pp5->SetLineWidth(1);
      //pp5->SetNDC();
      //pp5->Draw();
      
      
      
      
      //ratios go here
      temp_canvEvt->cd(); //change back to main canvas
      evtpad2->Draw();
      evtpad2->cd();
      //hist for ratio; use clone
      
      TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
      evtpad2->SetTopMargin(0);
      evtpad2->SetBottomMargin(0.3);
      
      theRatio->GetYaxis()->SetTitleSize(15);
      theRatio->GetYaxis()->SetTitleFont(43);
      theRatio->GetYaxis()->SetTitleOffset(2);
      theRatio->GetYaxis()->SetLabelFont(43); 
      theRatio->GetYaxis()->SetLabelSize(13);
      theRatio->GetXaxis()->SetTitleSize(20);
      theRatio->GetXaxis()->SetTitleFont(43);
      theRatio->GetXaxis()->SetTitleOffset(4.);
      theRatio->GetXaxis()->SetLabelFont(43); 
      theRatio->GetXaxis()->SetLabelSize(15);
      
      theRatio->SetTitle("");
      theRatio->SetXTitle( h_XAx_Title.c_str() );
      theRatio->SetYTitle("Data/MC");
      theRatio->Divide(theMCEvtQAHist);
      theRatio->Draw("h");
      
    	
      // print to PDF file
      temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
    
    }// end trig/pthat weighted vz


    {// MC weighted, data weighted
      temp_canvEvt->cd();
      //	temp_canvEvt->Divide(1,2);
      
      TPad *evtpad1 = new TPad("vz3evtpad1", "Overlay Pad",0.0,0.3,1.0,1.0);
      TPad *evtpad2 = new TPad("vz3evtpad2", "Ratio Pad",0.0,0.05,1.0,0.3);

      
      //temp_canvEvt->cd(1);
      evtpad1->SetBottomMargin(0); 
      evtpad1->Draw();
      evtpad1->cd();
      // open hist
      
      std::string inDataHistName  ="hWeightedVz";
      std::cout<<" drawing "<<inDataHistName<<std::endl<<std::endl;
      TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
      
      theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
      theDataEvtQAHist->Scale( 1./theLumi );
      theDataEvtQAHist->SetMarkerStyle(kCircle);
      theDataEvtQAHist->SetMarkerSize(0.99);
      theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
      theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );
      
      //get jet pt cut
      //std::string jetPTcut ="hJetPtCut";
      //TH1F* hJetPTCut= (TH1F*)finData->Get( jetPTcut.c_str() );
      //hJetPTCut->GetMean();
      //std::string EvtPt = ""; //std::to_string((int)hJetPTCut->GetMean()).c_str();
      //std::string EvtPtCut = "p_{T} > " + EvtPt + " GeV";
      //std::cout<<"EvtPtCut = "<<EvtPtCut.c_str()<<std::endl;
      
      std::string inMCHistName  ="hWeightedVz";
      TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
      
      theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
      theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
      theMCEvtQAHist->SetMarkerStyle(kMultiply);
      theMCEvtQAHist->SetMarkerSize(0.99);
      theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
      
      // title+axes
      std::string h_Title     ="EvtQA";
      std::string h_XAx_Title ="v_{z}^{evt} (cm)";
      std::string h_YAx_Title =crossSectionYAxis;
      theDataEvtQAHist->SetTitle (    h_Title.c_str() );
      theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
      theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theDataEvtQAHist->SetAxisRange(0.0,2.0, "Y");  
      
      theDataEvtQAHist->Draw();	
      theMCEvtQAHist->Draw("same");
      
      
      // legend
      float legx1=0.78, legx2=legx1+0.11;
      float legy1=0.74, legy2=legy1+0.09;
      TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
      
      theEvtQALeg->AddEntry(theDataEvtQAHist, "data" ,"lp");
      theEvtQALeg->AddEntry(theMCEvtQAHist,   "MC"       ,"lp");
      theEvtQALeg->Draw();
      
      // text
      float pp1x=0.14,pp1y=0.82;
      TLatex *pp1=new TLatex(pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets");
      pp1->SetTextFont(63);
      pp1->SetTextColor(kBlack);
      pp1->SetTextSize(15);
      pp1->SetLineWidth(1);
      pp1->SetNDC();
      pp1->Draw();
      TLatex *pp2=new TLatex(pp1x,(pp1y-0.04),"L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV");
      pp2->SetTextFont(63);
      pp2->SetTextColor(kBlack);
      pp2->SetTextSize(15);
      pp2->SetLineWidth(1);
      pp2->SetNDC();
      pp2->Draw();
      //TLatex *pp3=new TLatex(pp1x,(pp1y-0.10),(fullJetType+"Jets").c_str() );
      //pp3->SetTextFont(63);
      //pp3->SetTextColor(kBlack);
      //pp3->SetTextSize(17);
      //pp3->SetLineWidth(1);
      //pp3->SetNDC();
      //pp3->Draw();
      TLatex *pp4=new TLatex(pp1x,(pp1y-0.08),"Py8 Tune CUETP8M1" );
      pp4->SetTextFont(63);
      pp4->SetTextColor(kBlack);
      pp4->SetTextSize(15);
      pp4->SetLineWidth(1);
      pp4->SetNDC();
      pp4->Draw();
    
    
      ////tlatex for cut
      //TLatex *pp5=new TLatex(pp1x,(pp1y-0.11),EvtPtCut.c_str());
      //pp5->SetTextFont(63);
      //pp5->SetTextColor(kBlack);
      //pp5->SetTextSize(15);
      //pp5->SetLineWidth(1);
      //pp5->SetNDC();
      //pp5->Draw();
      
      
      
      
      //ratios go here
      temp_canvEvt->cd(); //change back to main canvas
      evtpad2->Draw();
      evtpad2->cd();
      //hist for ratio; use clone
      
      TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
      evtpad2->SetTopMargin(0);
      evtpad2->SetBottomMargin(0.3);
      
      theRatio->GetYaxis()->SetTitleSize(15);
      theRatio->GetYaxis()->SetTitleFont(43);
      theRatio->GetYaxis()->SetTitleOffset(2);
      theRatio->GetYaxis()->SetLabelFont(43); 
      theRatio->GetYaxis()->SetLabelSize(13);
      theRatio->GetXaxis()->SetTitleSize(20);
      theRatio->GetXaxis()->SetTitleFont(43);
      theRatio->GetXaxis()->SetTitleOffset(4.);
      theRatio->GetXaxis()->SetLabelFont(43); 
      theRatio->GetXaxis()->SetLabelSize(15);
      
      theRatio->SetTitle("");
      theRatio->SetXTitle( h_XAx_Title.c_str() );
      theRatio->SetYTitle("Data/MC");
      theRatio->Divide(theMCEvtQAHist);
      theRatio->Draw("h");
      
    	
      // print to PDF file
      temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
    
    }// endtotal weighted vz



    
    {// NUM JETS PER EVENT
      temp_canvEvt->cd();
      //temp_canvEvt_wLog->cd();
      //	temp_canvEvt->Divide(1,2);
      
      TPad *evtpad1 = new TPad("nrefevtpad1", "Overlay Pad",0.0,0.3,1.0,1.0);
      TPad *evtpad2 = new TPad("nrefevtpad2", "Ratio Pad",0.0,0.05,1.0,0.3);
    
      
      //temp_canvEvt->cd(1);
      evtpad1->SetBottomMargin(0); 
      evtpad1->SetLogy(1);
      evtpad1->Draw();
      evtpad1->cd();
      // open hist
      
      std::string inDataHistName  ="hNref";
      std::cout<<" drawing "<<inDataHistName<<std::endl<<std::endl;
      TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
      
      //theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
      theDataEvtQAHist->Scale( 1./theLumi );
      theDataEvtQAHist->SetMarkerStyle(kCircle);
      theDataEvtQAHist->SetMarkerSize(0.99);
      theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
      theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );
      
      //get jet pt cut
      //std::string jetPTcut ="hJetPtCut";
      //TH1F* hJetPTCut= (TH1F*)finData->Get( jetPTcut.c_str() );
      //hJetPTCut->GetMean();
      //std::string EvtPt = ""; //std::to_string((int)hJetPTCut->GetMean()).c_str();
      //std::string EvtPtCut = "p_{T} > " + EvtPt + " GeV";
      //std::cout<<"EvtPtCut = "<<EvtPtCut.c_str()<<std::endl;
      
      std::string inMCHistName  ="hNref";
      TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
      
      //theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
      theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
      theMCEvtQAHist->SetMarkerStyle(kMultiply);
      theMCEvtQAHist->SetMarkerSize(0.99);
      theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
      
      // title+axes
      std::string h_Title     ="EvtQA";
      std::string h_XAx_Title ="N_{jets}/{evt}";
      std::string h_YAx_Title =crossSectionYAxis;
      theDataEvtQAHist->SetTitle (    h_Title.c_str() );
      theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
      theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theDataEvtQAHist->SetAxisRange(0.0,2.0, "Y");  
      
      theDataEvtQAHist->Draw();	
      theMCEvtQAHist->Draw("same");
      
      
      // legend
      float legx1=0.78, legx2=legx1+0.11;
      float legy1=0.74, legy2=legy1+0.09;
      TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
      
      theEvtQALeg->AddEntry(theDataEvtQAHist, "data" ,"lp");
      theEvtQALeg->AddEntry(theMCEvtQAHist,   "MC"       ,"lp");
      theEvtQALeg->Draw();
      
      //// text
      //float pp1x=0.14,pp1y=0.82;
      //TLatex *pp1=new TLatex(pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets");
      //pp1->SetTextFont(63);
      //pp1->SetTextColor(kBlack);
      //pp1->SetTextSize(15);
      //pp1->SetLineWidth(1);
      //pp1->SetNDC();
      //pp1->Draw();
      //TLatex *pp2=new TLatex(pp1x,(pp1y-0.04),"L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV");
      //pp2->SetTextFont(63);
      //pp2->SetTextColor(kBlack);
      //pp2->SetTextSize(15);
      //pp2->SetLineWidth(1);
      //pp2->SetNDC();
      //pp2->Draw();
      //TLatex *pp3=new TLatex(pp1x,(pp1y-0.10),(fullJetType+"Jets").c_str() );
      //pp3->SetTextFont(63);
      //pp3->SetTextColor(kBlack);
      //pp3->SetTextSize(17);
      //pp3->SetLineWidth(1);
      //pp3->SetNDC();
      //pp3->Draw();
      //TLatex *pp4=new TLatex(pp1x,(pp1y-0.08),"Py8 Tune CUETP8M1" );
      //pp4->SetTextFont(63);
      //pp4->SetTextColor(kBlack);
      //pp4->SetTextSize(15);
      //pp4->SetLineWidth(1);
      //pp4->SetNDC();
      //pp4->Draw();
      //TLatex *pp5=new TLatex(pp1x,(pp1y-0.11),EvtPtCut.c_str());
      //pp5->SetTextFont(63);
      //pp5->SetTextColor(kBlack);
      //pp5->SetTextSize(15);
      //pp5->SetLineWidth(1);
      //pp5->SetNDC();
      //pp5->Draw();
      
      
      
      
      //ratios go here
      temp_canvEvt->cd(); //change back to main canvas
      evtpad2->Draw();
      evtpad2->cd();
      //hist for ratio; use clone
      
      TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
      evtpad2->SetTopMargin(0);
      evtpad2->SetBottomMargin(0.3);
      
      theRatio->GetYaxis()->SetTitleSize(15);
      theRatio->GetYaxis()->SetTitleFont(43);
      theRatio->GetYaxis()->SetTitleOffset(2);
      theRatio->GetYaxis()->SetLabelFont(43); 
      theRatio->GetYaxis()->SetLabelSize(13);
      theRatio->GetXaxis()->SetTitleSize(20);
      theRatio->GetXaxis()->SetTitleFont(43);
      theRatio->GetXaxis()->SetTitleOffset(4.);
      theRatio->GetXaxis()->SetLabelFont(43); 
      theRatio->GetXaxis()->SetLabelSize(15);
      
      theRatio->SetTitle("");
      theRatio->SetXTitle( h_XAx_Title.c_str() );
      theRatio->SetYTitle("Data/MC");
      theRatio->Divide(theMCEvtQAHist);
      theRatio->SetAxisRange(0.0, 2.0, "Y");
      theRatio->Draw("h");
      
    	
      // print to PDF file
      temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
    
    }// end numjets/event





    {// NUM JETS PER EVENT
      temp_canvEvt->cd();
      //temp_canvEvt_wLog->cd();
      //	temp_canvEvt->Divide(1,2);
      
      TPad *evtpad1 = new TPad("nrefevtpad1", "Overlay Pad",0.0,0.3,1.0,1.0);
      TPad *evtpad2 = new TPad("nrefevtpad2", "Ratio Pad",0.0,0.05,1.0,0.3);
    
      
      //temp_canvEvt->cd(1);
      evtpad1->SetBottomMargin(0); 
      evtpad1->SetLogy(1);
      evtpad1->Draw();
      evtpad1->cd();
      // open hist
      
      std::string inDataHistName  ="hWNref";
      std::cout<<" drawing "<<inDataHistName<<std::endl<<std::endl;
      TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
      
      //theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
      theDataEvtQAHist->Scale( 1./theLumi );
      theDataEvtQAHist->SetMarkerStyle(kCircle);
      theDataEvtQAHist->SetMarkerSize(0.99);
      theDataEvtQAHist->SetMarkerColor( theDataOverlayMarkerColor );
      theDataEvtQAHist->SetLineColor( theDataOverlayLineColor );
      
      //get jet pt cut
      //std::string jetPTcut ="hJetPtCut";
      //TH1F* hJetPTCut= (TH1F*)finData->Get( jetPTcut.c_str() );
      //hJetPTCut->GetMean();
      //std::string EvtPt = ""; //std::to_string((int)hJetPTCut->GetMean()).c_str();
      //std::string EvtPtCut = "p_{T} > " + EvtPt + " GeV";
      //std::cout<<"EvtPtCut = "<<EvtPtCut.c_str()<<std::endl;
      
      std::string inMCHistName  ="hWNref";
      TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
      
      //theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
      theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
      theMCEvtQAHist->SetMarkerStyle(kMultiply);
      theMCEvtQAHist->SetMarkerSize(0.99);
      theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
      
      // title+axes
      std::string h_Title     ="EvtQA";
      std::string h_XAx_Title ="N_{jets}/{evt}";
      std::string h_YAx_Title =crossSectionYAxis;
      theDataEvtQAHist->SetTitle (    h_Title.c_str() );
      theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
      theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theDataEvtQAHist->SetAxisRange(0.0,2.0, "Y");  
      
      theDataEvtQAHist->Draw();	
      theMCEvtQAHist->Draw("same");
      
      
      // legend
      float legx1=0.78, legx2=legx1+0.11;
      float legy1=0.74, legy2=legy1+0.09;
      TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
      
      theEvtQALeg->AddEntry(theDataEvtQAHist, "data" ,"lp");
      theEvtQALeg->AddEntry(theMCEvtQAHist,   "MC"       ,"lp");
      theEvtQALeg->Draw();
      
      //// text
      //float pp1x=0.14,pp1y=0.82;
      //TLatex *pp1=new TLatex(pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets");
      //pp1->SetTextFont(63);
      //pp1->SetTextColor(kBlack);
      //pp1->SetTextSize(15);
      //pp1->SetLineWidth(1);
      //pp1->SetNDC();
      //pp1->Draw();
      //TLatex *pp2=new TLatex(pp1x,(pp1y-0.04),"L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV");
      //pp2->SetTextFont(63);
      //pp2->SetTextColor(kBlack);
      //pp2->SetTextSize(15);
      //pp2->SetLineWidth(1);
      //pp2->SetNDC();
      //pp2->Draw();
      //TLatex *pp3=new TLatex(pp1x,(pp1y-0.10),(fullJetType+"Jets").c_str() );
      //pp3->SetTextFont(63);
      //pp3->SetTextColor(kBlack);
      //pp3->SetTextSize(17);
      //pp3->SetLineWidth(1);
      //pp3->SetNDC();
      //pp3->Draw();
      //TLatex *pp4=new TLatex(pp1x,(pp1y-0.08),"Py8 Tune CUETP8M1" );
      //pp4->SetTextFont(63);
      //pp4->SetTextColor(kBlack);
      //pp4->SetTextSize(15);
      //pp4->SetLineWidth(1);
      //pp4->SetNDC();
      //pp4->Draw();
    
    
      ////tlatex for cut
      //TLatex *pp5=new TLatex(pp1x,(pp1y-0.11),EvtPtCut.c_str());
      //pp5->SetTextFont(63);
      //pp5->SetTextColor(kBlack);
      //pp5->SetTextSize(15);
      //pp5->SetLineWidth(1);
      //pp5->SetNDC();
      //pp5->Draw();
      
      
      
      
      //ratios go here
      temp_canvEvt->cd(); //change back to main canvas
      evtpad2->Draw();
      evtpad2->cd();
      //hist for ratio; use clone
      
      TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
      evtpad2->SetTopMargin(0);
      evtpad2->SetBottomMargin(0.3);
      
      theRatio->GetYaxis()->SetTitleSize(15);
      theRatio->GetYaxis()->SetTitleFont(43);
      theRatio->GetYaxis()->SetTitleOffset(2);
      theRatio->GetYaxis()->SetLabelFont(43); 
      theRatio->GetYaxis()->SetLabelSize(13);
      theRatio->GetXaxis()->SetTitleSize(20);
      theRatio->GetXaxis()->SetTitleFont(43);
      theRatio->GetXaxis()->SetTitleOffset(4.);
      theRatio->GetXaxis()->SetLabelFont(43); 
      theRatio->GetXaxis()->SetLabelSize(15);
      
      theRatio->SetTitle("");
      theRatio->SetXTitle( h_XAx_Title.c_str() );
      theRatio->SetYTitle("Data/MC");
      theRatio->Divide(theMCEvtQAHist);
      theRatio->SetAxisRange(0.0, 2.0, "Y");
      theRatio->Draw("h");
      
    	
      // print to PDF file
      temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
    
    }// end numjets/event
    
    
    {//MC pthat, unnweighted
      temp_canvEvt_wLog->cd();	//I'm probably gonna need to divide here, too. (nope)
      
      // open hist
      std::string inHistName="hpthat";
      std::cout<<" drawing MC pthat hist "<<inHistName<<std::endl<<std::endl;
      TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
      
      theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
      theMCEvtQAHist->SetMarkerStyle(kDot);
      theMCEvtQAHist->SetMarkerSize(0.95);
      theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
	
      // title+axes
      std::string h_Title    ="MC EvtQA";
      std::string h_XAx_Title="MC evt #hat{p}_{t} (GeV)"       ;
      std::string h_YAx_Title="A.U."       ;
      
      theMCEvtQAHist->SetTitle (    h_Title.c_str() );
      theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
      theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
      
      theMCEvtQAHist->Draw();
      
      // test
      float t1x=.65,t1y=.81;
      TLatex *t1=new TLatex(t1x,t1y,"QCD MC, unweighted");
      t1->SetTextFont(63);
      t1->SetTextColor(kBlack);
      t1->SetTextSize(15);
      t1->SetLineWidth(1);
      t1->SetNDC();
      t1->Draw();	
      
      TLatex *t2=new TLatex(t1x,t1y-0.04,"Py8 Tune CUETP8M1");
      t2->SetTextFont(63);
      t2->SetTextColor(kBlack);
      t2->SetTextSize(15);
      t2->SetLineWidth(1);
      t2->SetNDC();
      t2->Draw();	
      
      // print to PDF file
      temp_canvEvt_wLog->Print( thePDFFileName.c_str() );	
    }
    
    {//MC pthat, weighted
      temp_canvEvt_wLog->cd();	//Don't divide here either
      
      // open hist
      std::string inHistName="hWeightedpthat";     
      std::cout<<" drawing MC weighted pthat hist "<<inHistName<<std::endl<<std::endl;
      TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inHistName.c_str() );
      
      theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
      theMCEvtQAHist->SetMarkerStyle(kDot);
      theMCEvtQAHist->SetMarkerSize(0.99);
      theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCEvtQAHist->SetLineColor( theMCOverlayLineColor );
      
      // title+axzes
      std::string h_Title    ="MC EvtQA" ;
      std::string h_XAx_Title="MC evt #hat{p}_{t} (GeV)"                ;
      std::string h_YAx_Title=crossSectionYAxis;//"#sigma (#mub/bin)"             ;
      
      theMCEvtQAHist->SetTitle (    h_Title.c_str() );
      theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
      theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theMCEvtQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
	
      theMCEvtQAHist->Draw();
      
      // text
      float t1x=.65,t1y=.81;
      TLatex *t1=new TLatex(t1x,t1y,"QCD MC, weighted");
      t1->SetTextFont(63);
      t1->SetTextColor(kBlack);
      t1->SetTextSize(15);
      t1->SetLineWidth(1);
      t1->SetNDC();
      t1->Draw();	
      
      TLatex *t2=new TLatex(t1x,t1y-0.04,"Py8 Tune CUETP8M1");
      t2->SetTextFont(63);
      t2->SetTextColor(kBlack);
      t2->SetTextSize(15);
      t2->SetLineWidth(1);
      t2->SetNDC();
      t2->Draw();	
      
      // print to pdf
      temp_canvEvt_wLog->Print( thePDFFileName.c_str() );	

    }
    
    temp_canvEvt->Close();  
  }//end total weighted vz
  
  




  std::string jetIDInt;
  if(doJetIDPlots)jetIDInt="1";
  else jetIDInt="0";
  
  // jet plots----------------------// creates jet QA plots
  if(!drawJetQAPlots)std::cout<<std::endl<<"skipping jet QA plots..."<<std::endl<<std::endl;
  else { 
    
    std::cout<<std::endl<<std::endl;    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing Jet plots..."<<std::endl;
    TCanvas *temp_canvJet = new TCanvas("tempJet", "tempJet", 700, 700);
    temp_canvJet->cd();
    
    TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad",0.0,0.3,1.0,1.0);
    jetpad1->SetLogy(1);
    jetpad1->SetGridx(1);
    jetpad1->SetGridy(1);
    jetpad1->SetBottomMargin(0);
    
    TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad",0.0,0.05,1.0,0.3);
    jetpad2->SetTopMargin(0);
    jetpad2->SetBottomMargin(0.3);
    jetpad2->SetLogy(0);
    

    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing no-log Jet plots..."<<std::endl;
    TCanvas *temp_canvJet_noLog = new TCanvas("tempJet_noLog", "tempJet_noLog", 600, 600);
    temp_canvJet_noLog->cd();
    temp_canvJet_noLog->SetLogy(0);

    TPad *jetpad_noLog1 = new TPad("jetpad_noLog1", "Overlay Pad",0.0,0.3,1.0,1.0);
    jetpad_noLog1->SetLogy(0);
    jetpad_noLog1->SetBottomMargin(0);

    TPad *jetpad_noLog2 = new TPad("jetpad_noLog2", "Ratio Pad",0.0,0.05,1.0,0.3);
    jetpad_noLog2->SetTopMargin(0);
    jetpad_noLog2->SetBottomMargin(0.3);

    // PUT JET CUTS SECTION IN HERE

    TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
    std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
    std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;

    TH1F* jtptQACut_h= (TH1F*)finData->Get( "hJetQAPtCut" );
    std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
    std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;

    TH1F* jtetaLoCut_h= (TH1F*)finData->Get( "hJetEtaCutLo" );
    std::stringstream etaLo; etaLo.precision(1);
    float jtetaLoCut_F=jtetaLoCut_h->GetMean();
    etaLo << std::fixed << jtetaLoCut_F;
    std::string jtetaLoCut_str = etaLo.str();
    std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;

    TH1F* jtetaHiCut_h= (TH1F*)finData->Get( "hJetEtaCutHi" );
    std::stringstream etaHi; etaHi.precision(1);
    float jtetaHiCut_F=jtetaHiCut_h->GetMean();
    etaHi << std::fixed << jtetaHiCut_F;
    std::string jtetaHiCut_str = etaHi.str();
    std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;

    std::string jetCutString="p_{t}>"+jtptCut_str+" GeV, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
    std::cout<<"jetCutString="<<jetCutString<<std::endl;

    float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
    std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;

    // DIJET CUT HISTS
    //TH1F* dijetCut_h= (TH1F*)finData->Get( "" );
    //std::string dijetCut_str = std::to_string((int)dijetCut_h->GetMean()).c_str();
    //std::cout<<"dijetCut_str = "<<dijetCut_str<<std::endl;
    //TH1F* dijetCut_h= (TH1F*)finData->Get( "" );
    //std::string dijetCut_str = std::to_string((int)dijetCut_h->GetMean()).c_str();
    //std::cout<<"dijetCut_str = "<<dijetCut_str<<std::endl;
    //TH1F* dijetCut_h= (TH1F*)finData->Get( "" );
    //std::string dijetCut_str = std::to_string((int)dijetCut_h->GetMean()).c_str();
    //std::cout<<"dijetCut_str = "<<dijetCut_str<<std::endl;



    //hJetEtaCutHi;12.400000
    //hJetEtaCutLo;10.000000
    //hJetQAPtCut;115.000000
    //hLdJetPtCut;155.000000
    //hSubldJetPtCut;120.000000
    //hPtAveCut;125.000000
	    
        
    std::cout<<" drawing jet QA Plots..."<<std::endl;
    for(int j=0;j<N_vars;j++){ std::cout<<std::endl;



      std::cout<<std::endl<<" var ="<<var[j]<<", j="<<j<<std::endl;
      
      if(!drawJetConstituentPlots && 
	 j>=jetConstits_varStart && j<dijet_varStart){
	if(debugMode) std::cout<<"skipping jet consituent plots..."<<std::endl;
	continue;}
      else if(!drawDijetPlots && 
	 j>=dijet_varStart){
	if(debugMode) std::cout<<"skipping dijet plots..."<<std::endl;
	continue;}
      
      std::string inHistName="hJetQA_"+jetIDInt+"wJetID_"+var[j];      
      
      if(debugMode)std::cout<<" opening input Data "<<inHistName<<std::endl<<std::endl;      

      //if(j==27)continue;//temporary; issue with muMax plot?
      //assert(false);

      TH1F* theDataJetQAHist= (TH1F*)
	( (TH1*)finData->Get(inHistName.c_str()) )
	->TH1::Rebin(var_xAx_reBin[j]);	
      theDataJetQAHist->Print("base"); 
      std::cout<<std::endl;
      
      theDataJetQAHist->Scale( 1./theDataJetQAHist->GetBinWidth(1) );
      theDataJetQAHist->Scale( 1./theLumi);
      theDataJetQAHist->Scale( 1./etaBinWidth);//eta bin width scaling

      theDataJetQAHist->SetMarkerStyle(kCircle);
      theDataJetQAHist->SetMarkerSize(0.90);
      theDataJetQAHist->SetMarkerColor( theDataOverlayMarkerColor );
      theDataJetQAHist->SetLineColor( theDataOverlayLineColor );
      
      
      if(debugMode)std::cout<<" opening input MC "<<inHistName<<std::endl<<std::endl;      
      TH1F* theMCJetQAHist= (TH1F*)
	( (TH1*)finMC->Get(inHistName.c_str()) )
	->TH1::Rebin(var_xAx_reBin[j]);	
      theMCJetQAHist->Print("base");
      std::cout<<std::endl;
      
      theMCJetQAHist->Scale( 1./theMCJetQAHist->GetBinWidth(1) );	  
      theMCJetQAHist->Scale(1000.);
      theMCJetQAHist->Scale(1./etaBinWidth);//eta bin width scaling
      
      if( var[j] != "jtpt"      &&
      	  var[j] != "rawpt"     &&
      	  var[j] != "leadJetPt" &&
      	  var[j] != "subleadJetPt" ) theMCJetQAHist->Scale( theDataJetQAHist->Integral()/theMCJetQAHist->Integral() );	        
      
      theMCJetQAHist->SetMarkerStyle(kMultiply);
      theMCJetQAHist->SetMarkerSize(0.90);
      theMCJetQAHist->SetMarkerColor( theMCOverlayMarkerColor );
      theMCJetQAHist->SetLineColor( theMCOverlayLineColor );
      
      
      std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title=crossSectionYAxis;//"#sigma (#mub/bin)";
      if(debugMode)std::cout<<" x-axis title= "<<h_XAx_Title<<std::endl<<std::endl;
      
      std::string h_Title="JetQA";//, "+var[j];
      if(doJetIDPlots) h_Title+=", w/ JetIDCut";      
      
      theDataJetQAHist->SetTitle (    h_Title.c_str() );
      theDataJetQAHist->SetXTitle( h_XAx_Title.c_str() );
      theDataJetQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theDataJetQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
      
      theMCJetQAHist->SetTitle (    h_Title.c_str() );
      theMCJetQAHist->SetXTitle( h_XAx_Title.c_str() );
      theMCJetQAHist->SetYTitle( h_YAx_Title.c_str() );
      //theMCJetQAHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
      
      theDataJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
      theMCJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
      


      // variable specific draw orders, canv, etc.
      if(var[j]=="jtpt"||var[j]=="rawpt") {
	
	
	temp_canvJet->cd();			//here is the draw command for pt 

	jetpad1->SetBottomMargin(0);
	jetpad1->SetLogy(1);
	jetpad1->Draw();
	jetpad1->cd();

	theMCJetQAHist->Draw(); 
	theDataJetQAHist->Draw("same");  
	
	//float t1Loc1=0.73, t1Loc2=0.68; 	    
	float t1Loc1=0.54, t1Loc2=0.82;
	TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	t1->SetTextFont(63);
	t1->SetTextColor(kBlack);
	t1->SetTextSize(14);
	t1->SetLineWidth(1);
	t1->SetNDC();
	t1->Draw();
	

	TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),(jetCutString).c_str());
	t2->SetTextFont(63);
	t2->SetTextColor(kBlack);
	t2->SetTextSize(14);
	t2->SetLineWidth(1);
	t2->SetNDC();
	t2->Draw();	
	
	float legx1=0.77, legx2=legx1+0.11;
	float legy1=0.76, legy2=legy1+0.09;
	TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	theJetQALeg->SetFillStyle(0);
	theJetQALeg->AddEntry(theDataJetQAHist,"data","lp");
	theJetQALeg->AddEntry(theMCJetQAHist,  "MC"  ,"lp");
	theJetQALeg->Draw(); 
	
	
	temp_canvJet->cd();
	jetpad2->Draw();
	jetpad2->cd();
	jetpad2->SetTopMargin(0);
	TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
	
	//theRatio->SetMarkerStyle( kFullCross );
	theRatio->SetMarkerSize(0.90 );
	theRatio->SetMarkerColor( kRed );
	theRatio->SetLineColor( theRatioLineColor );
	
	theRatio->SetAxisRange(0.,2.,"Y");
	theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	
	theRatio->GetYaxis()->SetTitleSize(15);
	theRatio->GetYaxis()->SetTitleFont(43);
	theRatio->GetYaxis()->SetTitleOffset(2);
	theRatio->GetYaxis()->SetLabelFont(43); 
	theRatio->GetYaxis()->SetLabelSize(13);
	theRatio->GetXaxis()->SetTitleSize(20);
	theRatio->GetXaxis()->SetTitleFont(43);
	theRatio->GetXaxis()->SetTitleOffset(4.);
	theRatio->GetXaxis()->SetLabelFont(43); 
	theRatio->GetXaxis()->SetLabelSize(15);
	
	theRatio->SetTitle("");
	theRatio->SetXTitle( h_XAx_Title.c_str() );
	theRatio->SetYTitle("Data/MC");
	theRatio->Divide(theMCJetQAHist);
	theRatio->Draw("HIST E");
	
	TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
	lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
	lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
	lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	
	lineAtOne->Draw("same");
	lineAtOneHalf->Draw("same");
	lineAtOneEtOneHalf->Draw("same");
	
	
	temp_canvJet->Print( thePDFFileName.c_str() ); }
	
	else if(var[j]=="jteta" ||
		var[j]=="jtphi"	||
		( j>=jetConstits_varStart && j<dijet_varStart )){ 

		
	  temp_canvJet->cd();			//here is the draw command for pt 
	  jetpad1->SetBottomMargin(0);
	  jetpad1->SetLogy(1);
	  jetpad1->Draw();
	  jetpad1->cd();

	  theDataJetQAHist->Draw(); 
	  theMCJetQAHist->Draw("SAME"); 

	  float t1Loc1=0.54, t1Loc2=0.82;
	  //float t1Loc1=0.50, t1Loc2=0.72; 	    
	  TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	  t1->SetTextFont(63);
	  t1->SetTextColor(kBlack);
	  t1->SetTextSize(14);
	  t1->SetLineWidth(1);
	  t1->SetNDC();
	  t1->Draw();
	  
	  TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),(jetCutString).c_str());
	  t2->SetTextFont(63);
	  t2->SetTextColor(kBlack);
	  t2->SetTextSize(14);
	  t2->SetLineWidth(1);
	  t2->SetNDC();
	  t2->Draw();	

	  float legx1=0.77, legx2=legx1+0.11;
	  float legy1=0.76, legy2=legy1+0.09;
	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	  theJetQALeg->SetFillStyle(0);
	  theJetQALeg->AddEntry(theDataJetQAHist,"data","lp");
	  theJetQALeg->AddEntry(theMCJetQAHist,  "MC"  ,"lp");
	  theJetQALeg->Draw();

	  //Ratio Plots
	  //temp_canvJet->cd();	  
	  jetpad2->cd();
	  jetpad2->Draw();
	  TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
	  
	  theRatio->SetMarkerColor( kRed );
	  theRatio->SetLineColor( theRatioLineColor );
	  
	  theRatio->SetAxisRange(0.,2.,"Y");
	  theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  
	  theRatio->GetYaxis()->SetTitleSize(15);
	  theRatio->GetYaxis()->SetTitleFont(43);
	  theRatio->GetYaxis()->SetTitleOffset(2);
	  theRatio->GetYaxis()->SetLabelFont(43); 
	  theRatio->GetYaxis()->SetLabelSize(13);
	  theRatio->GetXaxis()->SetTitleSize(20);
	  theRatio->GetXaxis()->SetTitleFont(43);
	  theRatio->GetXaxis()->SetTitleOffset(4.);
	  theRatio->GetXaxis()->SetLabelFont(43); 
	  theRatio->GetXaxis()->SetLabelSize(15);
	  
	  theRatio->SetTitle("");
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle("Data/MC");
	  theRatio->Divide(theMCJetQAHist);
	  theRatio->Draw("h");
	  
	  TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
	  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	  TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
	  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	  TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
	  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	  
	  lineAtOne->Draw("same");
	  lineAtOneHalf->Draw("same");
	  lineAtOneEtOneHalf->Draw("same");
	  
	  
	  temp_canvJet->Print( thePDFFileName.c_str() ); }
      
	else if(var[j]=="Aj"||
		var[j]=="xj"||
		var[j]=="dphi"){ 
	  //continue;
	  temp_canvJet_noLog->cd();	  
	  jetpad_noLog1->SetBottomMargin(0);
	  jetpad_noLog1->Draw();
	  jetpad_noLog1->cd();
	  theDataJetQAHist->Draw(); 
	  theMCJetQAHist->Draw("SAME"); 
	  
	  float t1Loc1=0.15, t1Loc2=0.7; 	    	  
	  if(var[j]=="Aj") t1Loc1=.7;

	  
	  TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	  t1->SetTextFont(63);
	  t1->SetTextColor(kBlack);
	  t1->SetTextSize(14);
	  t1->SetLineWidth(0.95);
	  t1->SetNDC();
	  t1->Draw();
	  
	  TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.04),"p_{t}^{ld} > 60 GeV");
	  t2->SetTextFont(63);
	  t2->SetTextColor(kBlack);
	  t2->SetTextSize(14);
	  t2->SetLineWidth(0.95);
	  t2->SetNDC();
	  t2->Draw();	
	  
	  TLatex *t3=new TLatex((t1Loc1),(t1Loc2-.08),"p_{t}^{subld} > 30GeV");
	  t3->SetTextFont(63);
	  t3->SetTextColor(kBlack);
	  t3->SetTextSize(14);
	  t3->SetLineWidth(0.95);
	  t3->SetNDC();
	  t3->Draw();	
	  
	  TLatex *t4=new TLatex(t1Loc1, (t1Loc2-0.14), "#Delta#phi > #frac{2}{3}#pi");	    
	  t4->SetTextFont(63);
	  t4->SetTextColor(kBlack);
	  t4->SetTextSize(14);
	  t4->SetLineWidth(0.95);
	  t4->SetNDC();
	  t4->Draw();	

	  float legx1=0.15, legx2=legx1+0.11;
	  if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
	  float legy1=0.74, legy2=legy1+0.09;

	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	  theJetQALeg->SetFillStyle(0);
	  theJetQALeg->AddEntry(theDataJetQAHist,"data","lp");
	  theJetQALeg->AddEntry(theMCJetQAHist,  "MC"  ,"lp");
	  theJetQALeg->Draw(); 

	  //Ratio Plots
	  temp_canvJet_noLog->cd();	  
	  jetpad_noLog2->Draw();
	  jetpad_noLog2->cd();
	  TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
	  //I need to redefine the axis titles here
	  std::string h_XAx_ratTitle=var_xAx_Titles[j], h_YAx_ratTitle="Data/MC"    ;
	  std::string h_ratTitle   ="JetQA, Data/MC";//, "+var[j];
	  if(doJetIDPlots)h_ratTitle+=", w/ JetIDCut";      
	  
	  theRatio->SetMarkerColor( kRed );
	  theRatio->SetLineColor( theRatioLineColor );
	  
	  theRatio->SetAxisRange(0.,2.,"Y");
	  theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  
	  theRatio->GetYaxis()->SetTitleSize(15);
	  theRatio->GetYaxis()->SetTitleFont(43);
	  theRatio->GetYaxis()->SetTitleOffset(2);
	  theRatio->GetYaxis()->SetLabelFont(43); 
	  theRatio->GetYaxis()->SetLabelSize(13);
	  theRatio->GetXaxis()->SetTitleSize(20);
	  theRatio->GetXaxis()->SetTitleFont(43);
	  theRatio->GetXaxis()->SetTitleOffset(4.);
	  theRatio->GetXaxis()->SetLabelFont(43); 
	  theRatio->GetXaxis()->SetLabelSize(15);
	  
	  theRatio->SetTitle("");
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle("Data/MC");
	  theRatio->Divide(theMCJetQAHist);
	  theRatio->Draw();
	  
	  TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
	  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	  TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
	  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	  TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
	  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	  
	  lineAtOne->Draw("same");
	  lineAtOneHalf->Draw("same");
	  lineAtOneEtOneHalf->Draw("same");
	  
	  temp_canvJet_noLog->Print( thePDFFileName.c_str() ); }
      
	else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt"){ 
	  //continue;
	  temp_canvJet->cd(0);
	  jetpad1->Draw();
	  jetpad1->cd();
	  theDataJetQAHist->Draw(); 
	  theMCJetQAHist->Draw("SAME"); 

	  float t1Loc1=0.7, t1Loc2=0.72; 	    	  
	  TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	  t1->SetTextFont(63);
	  t1->SetTextColor(kBlack);
	  t1->SetTextSize(14);
	  t1->SetLineWidth(1);
	  t1->SetNDC();
	  t1->Draw();
	  
	  TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}^{ld} > 60 GeV");
	  t2->SetTextFont(63);
	  t2->SetTextColor(kBlack);
	  t2->SetTextSize(14);
	  t2->SetLineWidth(1);
	  t2->SetNDC();
	  t2->Draw();	
	  
	  TLatex *t3=new TLatex((t1Loc1),(t1Loc2-.10),"p_{t}^{subld} > 30GeV");
	  t3->SetTextFont(63);
	  t3->SetTextColor(kBlack);
	  t3->SetTextSize(14);
	  t3->SetLineWidth(1);
	  t3->SetNDC();
	  t3->Draw();	
	  
	  TLatex *t4=new TLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");	    
	  t4->SetTextFont(63);
	  t4->SetTextColor(kBlack);
	  t4->SetTextSize(14);
	  t4->SetLineWidth(1);
	  t4->SetNDC();
	  t4->Draw();	

	  float legx1=0.70, legx2=legx1+0.11;
	  float legy1=0.76, legy2=legy1+0.09;
	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	theJetQALeg->SetFillStyle(0);
	  theJetQALeg->AddEntry(theDataJetQAHist,"data","lp");
	  theJetQALeg->AddEntry(theMCJetQAHist,  "MC"  ,"lp");
	  theJetQALeg->Draw(); 
	  
	  //Ratio Plots
	  temp_canvJet->cd();	  
	  jetpad2->Draw();
	  jetpad2->cd();
	  TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
	  
	  theRatio->SetMarkerColor( kRed );
	  theRatio->SetLineColor( theRatioLineColor );
	  
	  theRatio->SetAxisRange(0.,2.,"Y");
	  theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  
	  theRatio->GetYaxis()->SetTitleSize(15);
	  theRatio->GetYaxis()->SetTitleFont(43);
	  theRatio->GetYaxis()->SetTitleOffset(2);
	  theRatio->GetYaxis()->SetLabelFont(43); 
	  theRatio->GetYaxis()->SetLabelSize(13);
	  theRatio->GetXaxis()->SetTitleSize(20);
	  theRatio->GetXaxis()->SetTitleFont(43);
	  theRatio->GetXaxis()->SetTitleOffset(4.);
	  theRatio->GetXaxis()->SetLabelFont(43); 
	  theRatio->GetXaxis()->SetLabelSize(15);
	  
	  theRatio->SetTitle ("");
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle( "Data/MC" );
	  
	  theRatio->Divide(theMCJetQAHist);
	  theRatio->Draw("h");
	  
	  TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
	  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	  TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
	  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	  TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
	  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	  
	  lineAtOne->Draw("same");
	  lineAtOneHalf->Draw("same");
	  lineAtOneEtOneHalf->Draw("same");
	  
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
      //  t1->SetTextSize(10);
      //  t1->SetLineWidth(1);
      //  t1->SetNDC();
      //  t1->Draw();
      //  TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_t^1>60 GeV, p_t^2> 30GeV, #delta#phi<#frac{2}{3}#pi");
      //  t2->SetTextFont(63);
      //  t2->SetTextColor(kBlack);
      //  t2->SetTextSize(10);
      //  t2->SetLineWidth(1);
      //  t2->SetNDC();
      //  t2->Draw();	
      //  temp_canvJet->Print( thePDFFileName.c_str() ); }		  
      
	//else { temp_canvJet->cd(1); //Else what?
	//
	//  theDataJetQAHist->Draw(); 
	//  theMCJetQAHist->Draw("SAME");  
	//  
	//  float t1Loc1=0.55, t1Loc2=0.83; 	    
	//  TLatex *t1=new TLatex(t1Loc1,t1Loc2, (fullJetType+"Jets").c_str());
	//  t1->SetTextFont(63);
	//  t1->SetTextColor(kBlack);
	//  t1->SetTextSize(10);
	//  t1->SetLineWidth(1);
	//  t1->SetNDC();
	//  t1->Draw();
	//  
	//  TLatex *t2=new TLatex(t1Loc1, (t1Loc2-.05),"p_{t}>50 GeV, |#eta|<2.0");
	//  t2->SetTextFont(63);
	//  t2->SetTextColor(kBlack);
	//  t2->SetTextSize(10);
	//  t2->SetLineWidth(1);
	//  t2->SetNDC();
	//  t2->Draw();	
	//  
	//  float legx1=0.77, legx2=legx1+0.11;
	//  float legy1=0.77, legy2=legy1+0.09;
	//  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	//  theJetQALeg->SetFillStyle(0);
	//  theJetQALeg->AddEntry(theDataJetQAHist,"data","lp");
	//  theJetQALeg->AddEntry(theMCJetQAHist,  "MC"  ,"lp");
	//  theJetQALeg->Draw(); 
	//
	//  //temp_canvJet->Print( thePDFFileName.c_str() ); }
	//}

      std::cout<<std::endl; 	

      //      }//end draw overlays
      
      //else if(drawDataMCRatios){ temp_canvJetRat->cd();	  //redirect this part to .cd(2)
      //      else if(drawDataMCRatios){
      /*temp_canvJet->cd(2);	  
	TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone("dataJetHistClone4Ratio");
	//I need to redefine the axis titles here
	std::string h_XAx_ratTitle=var_xAx_Titles[j], h_YAx_ratTitle="Data/MC"    ;
	std::string h_ratTitle   ="JetQA, Data/MC";//, "+var[j];
	if(doJetIDPlots)h_ratTitle+=", w/ JetIDCut";      
	
	theRatio->SetTitle (    h_ratTitle.c_str() );
	theRatio->SetXTitle( h_XAx_ratTitle.c_str() );
	theRatio->SetYTitle( h_YAx_ratTitle.c_str() );
	theRatio->SetMarkerColor( kRed );
	theRatio->SetLineColor( theRatioLineColor );
	
	theRatio->SetAxisRange(0,2.,"Y");
	theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	
	
	theRatio->Divide(theMCJetQAHist);
	theRatio->Draw();
	
	TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j],1.0); 
	lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j],0.5); 
	lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j],1.5); 
	lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	
	lineAtOne->Draw("same");
	lineAtOneHalf->Draw("same");
	lineAtOneEtOneHalf->Draw("same");
	
	temp_canvJet->Print( thePDFFileName.c_str() );	*/
      //      }//end draw ratios 
      
      std::cout<<std::endl;	    	
      
    }//end loop over jetQA Plots     

    temp_canvJet_noLog->Close();//clean up neatly
    //    temp_canvJetRat->Close();
    temp_canvJet->Close();
  }
  

  // jets rap-bins plots----------------------//this is set to "false" right now
  if( !drawJetRapBinsPlot || drawDataMCRatios )
    std::cout<<std::endl<<"skipping jet rap bin plot..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl;
    std::cout<<" drawing jet rap bin Plots..."<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing Jet rapBins plots..."<<std::endl;
    TCanvas *temp_canvJetRapBins = new TCanvas("tempJetRapBins", "tempJetRapBins", 900, 600);
    temp_canvJetRapBins->SetLogy(1);    

    TLegend* theJetRapHistLeg=new TLegend(0.72,0.60,0.85,0.86, NULL,"brNDC");
    //TLegend* theJetRapHistLeg=new TLegend();

    
    if(drawDataMCOverlays){     temp_canvJetRapBins->cd();
      
      float power=std::pow(10., 6.);
      //float power=std::pow(10., 0.);
      std::cout<<"drawing dual-differential cross section plot"<<std::endl;	
      std::cout<<"looping over "<< nbins_abseta <<" rapidity bins"<<std::endl;
      
      for(int rapbin=0; rapbin<8; ++rapbin){	
	
	// assert(rapbin!=0);
	// assert(rapbin!=-1);
	
	int bincounter = nbins_abseta-rapbin;
	std::cout<<"bincounter = "<<bincounter<<" "<<std::endl;

	//open the hists + rebin
	std::string inHistName="hJetSpectraRap";	
	if(doJetIDPlots) inHistName+="_wJetID"; 
	inHistName+="_bin"+std::to_string(rapbin);
	if(debugMode)std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl<<std::endl;	

	if(rapbin==6){
	  std::cout<<"rapbin==6, continuing on..."<<std::endl;
	  continue;//still not sure why i need this...
	}
	
	TH1F* theJetSpectraRapHist= (TH1F*)
	  ( (TH1*)finData->Get(inHistName.c_str()) )
	  ->TH1::Rebin(jetSpectraRapBinning);

	theJetSpectraRapHist->Print( "base" );

	theJetSpectraRapHist->Scale( 1./theJetSpectraRapHist->GetBinWidth(1) );
	theJetSpectraRapHist->Scale( 1./theLumi);
	theJetSpectraRapHist->Scale( power );
	
	theJetSpectraRapHist->SetMarkerStyle( kCircle ) ;//theRapOverlayMarker[j]);
	theJetSpectraRapHist->SetMarkerSize(0.99);
	theJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
	theJetSpectraRapHist->SetLineColor( kBlack );//theRapOverlayLineColor[j]   );
	theJetSpectraRapHist->SetAxisRange(0.,800.,"X");
	
	
	//TH1F* theMCJetSpectraRapHist= (TH1F*)
	//  ( (TH1*)finMC->Get(inMCHistName.c_str()) )
	//  ->TH1::Rebin(jetSpectraRapBinning);
	//
	//theMCJetSpectraRapHist->Scale( 1./theMCJetSpectraRapHist->GetBinWidth(1) );
	//theMCJetSpectraRapHist->Scale( power );
	//theMCJetSpectraRapHist->Scale( theJetSpectraRapHist->Integral()/theMCJetSpectraRapHist->Integral() );
	//
	//theMCJetSpectraRapHist->SetMarkerStyle( 27 ) ;//theRapOverlayMarker[j]);
	//theMCJetSpectraRapHist->SetMarkerSize(0.99);
	//theMCJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
	//theMCJetSpectraRapHist->SetLineColor( -1 );//theRapOverlayLineColor[j]   );
	////theJetSpectraRapHist->SetAxisRange(0.,800.,"X");	 
	
	std::stringstream stream1, stream2;
	stream1.precision(1); stream1 << std::fixed << absetabins[rapbin];
	stream2.precision(1); stream2 << std::fixed << absetabins[rapbin+1];
	
	std::string rapHistLegDesc=stream1.str()+"<|y|<"+stream2.str()+" x10^{"+std::to_string(bincounter-1)+"}";	
	theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
	
	if(rapbin==0){

	  std::string h_Title   ="Jet Pt in |y| bins";
	  if(doJetIDPlots) h_Title+=", w/ JetIDCut";
	  theJetSpectraRapHist->SetTitle (    h_Title.c_str() );

	  std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title=crossSectionYAxis;
	  theJetSpectraRapHist->SetXTitle( h_XAx_Title.c_str() );
	  theJetSpectraRapHist->SetYTitle( h_YAx_Title.c_str() );
	  theJetSpectraRapHist->SetAxisRange(std::pow(10.,-9.3),std::pow(10.,5.),"Y");
	  theJetSpectraRapHist->SetAxisRange(10.,800.,"X");

	  theJetSpectraRapHist->Draw();     	  
	  theJetRapHistLeg->Draw(); }
	else { 
	  theJetSpectraRapHist->Draw("SAME");}
	power/=10.;	
      }//end rapbin hist loop
    
    float t1Loc1=0.42, t1Loc2=0.83; 	    
    TLatex *t1=new TLatex(t1Loc1,t1Loc2,"pp 2015 promptReco, #sqrt{s} = 5.02 TeV");
    t1->SetTextFont(63);
    t1->SetTextColor(kBlack);
    t1->SetTextSize(10);
    t1->SetLineWidth(1);
    t1->SetNDC();
    t1->Draw();
    
    TLatex *pp1=new TLatex(t1Loc1,(t1Loc2-0.04),"Jet80+LowerJets, L_{int} = 27.4 pb^{-1}");
    pp1->SetTextFont(63);
    pp1->SetTextColor(kBlack);
    pp1->SetTextSize(10);
    pp1->SetLineWidth(1);
    pp1->SetNDC();
    pp1->Draw();
    
    TLatex *pp2=new TLatex(t1Loc1,(t1Loc2-0.08),(fullJetType+"Jets, p_{t}>15 GeV").c_str());
    pp2->SetTextFont(63);
    pp2->SetTextColor(kBlack);
    pp2->SetTextSize(10);
    pp2->SetLineWidth(1);
    pp2->SetNDC();
    pp2->Draw();
    
    //TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.05),"p_{t}>15 GeV, |#eta|<2.0");
    //t2->SetTextFont(63);
    //t2->SetTextColor(kBlack);
    //t2->SetTextSize(10);
    //t2->SetLineWidth(1);
    //t2->SetNDC();
    //t2->Draw();	
    
    temp_canvJetRapBins->Print( thePDFFileName.c_str() );

    }//end draw mcdata overlays    
  }//end jet pt rap bins plot
  
  
  
  
  // jet trig plots----------------------
  if(!drawJetTrigQAPlots)
    std::cout<<std::endl<<"skipping jet trig QA plots..."<<std::endl<<std::endl;
  else { 
    std::cout<<std::endl<<std::endl;
    std::cout<<" drawing jet trig QA Plots..."<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
    TCanvas *temp_canvJetTrig = new TCanvas("tempJetTrig", "tempJetTrig", 1200, 600);
    temp_canvJetTrig->cd();
    temp_canvJetTrig->SetLogy(1);    
    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
    TCanvas *temp_canvJetRatio = new TCanvas("tempJetRatio", "tempJetRatio", 1200, 600);
    temp_canvJetRatio->cd();
    temp_canvJetRatio->SetLogy(0);
    
    
    if(drawTrigDataOverlays) { temp_canvJetTrig->cd();
      
      TLegend* JetTrigLegend=new TLegend(0.75,0.60,0.88,0.86, NULL,"brNDC");
      
      if(debugMode)	std::cout<<"drawing dataTrig QA Overlay Plots, looping over trigs"<<std::endl;
      for(int j=0; j<N_trigs; j++){ 
	
	//open the hists + do scaling 
	std::string inHistName="hpp_"+HLTName[j];
	if(doJetIDPlots)inHistName+="_wJetID";
	inHistName+="_"+radius+etaWidth;	
	if(debugMode)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
	TH1F* theJetTrigQAHist= (TH1F*)
	  ( (TH1*)finData->Get(inHistName.c_str()) )
	  ->TH1::Rebin(jetTrigQABinning);
		
	theJetTrigQAHist->Scale( 1./theJetTrigQAHist->GetBinWidth(1) );
	theJetTrigQAHist->Scale( 1./theLumi);
	theJetTrigQAHist->Scale( 1./4.);//eta bin width

	theJetTrigQAHist->SetMarkerStyle(theTrigOverlayMarker[j]);
	theJetTrigQAHist->SetMarkerSize(0.99);
	theJetTrigQAHist->SetMarkerColor( theTrigOverlayMarkerColor[j] );
	theJetTrigQAHist->SetLineColor(   theTrigOverlayLineColor[j]   );
	
	JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTName_Leg[j]).c_str(),"lp");
	
	if(j==0){
	  std::string h_Title   ="JetTrigQA";//, DATA ONLY";
	  if(doJetIDPlots) h_Title+=", w/ JetIDCut";      
	  //else h_Title+=", no JetIDCut";
	  theJetTrigQAHist->SetTitle (    h_Title.c_str() );

	  std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title=crossSectionYAxis;
	  theJetTrigQAHist->SetAxisRange(0.,1000.,"X");
	  theJetTrigQAHist->SetXTitle( h_XAx_Title.c_str() );
	  theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );

	  theJetTrigQAHist->Draw(); 	    	}
	else 
	  theJetTrigQAHist->Draw("SAME");   	    
	
      }//end loop over trigs
      JetTrigLegend->Draw();
      temp_canvJetTrig->Print(thePDFFileName.c_str());
    }//end DRAW TRIG OVERLAY      
    
    else if(drawTrigDataRatios){ temp_canvJetRatio->cd();

      TLegend* JetTrigRatioLegend=new TLegend(0.12,0.59,0.30,0.87, NULL,"brNDC");    
      
      if(debugMode)std::cout<<"drawing dataTrig QA Overlay Plots"<<std::endl;
      
      std::string comboPlotName="hpp_"+HLTName[N_trigs-1];
      if(doJetIDPlots)comboPlotName+="_wJetID";
      comboPlotName+="_"+radius+etaWidth;
      TH1F* theComboHist=(TH1F*)
	( (TH1*)finData->Get(comboPlotName.c_str()) )
	  ->TH1::Rebin(jetTrigQABinning/2);
      
      for(int j=0; j<(N_trigs-1); j++){ 
	
	std::string inHistName="hpp_"+HLTName[j];
	if(doJetIDPlots)inHistName+="_wJetID";
	inHistName+="_"+radius+etaWidth;	
	TH1F* theRatio=(TH1F*)
	  ( (TH1*)finData->Get(inHistName.c_str()) )
	  ->TH1::Rebin(jetTrigQABinning/2);
	
	theRatio->SetMarkerStyle(theTrigOverlayMarker[j]);
	theRatio->SetMarkerSize(0.99);
	theRatio->SetMarkerColor( theTrigOverlayMarkerColor[j] );
	theRatio->SetLineColor( theTrigOverlayLineColor[j] );	

	JetTrigRatioLegend->AddEntry(theRatio,(HLTName_Leg[j]+"/TrigComb.").c_str(),"lp");

	TH1F* comboHistClone=(TH1F*)theComboHist->Clone("comboHistClone4Div");
       	theRatio->Divide(comboHistClone);
	
	if(j==0){
	  std::string h_Title   ="JetTrigQA Ratios";
	  if(doJetIDPlots)    h_Title+=", w/ JetIDCut";      
	  theRatio->SetTitle (    h_Title.c_str() );

	  std::string h_XAx_Title="reco JetPt (GeV)", h_YAx_Title="Trig/TrigCombo"    ; 
	  theRatio->SetAxisRange(0.,1000.,"X");
	  theRatio->SetAxisRange(0.,2.,"Y");
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle( h_YAx_Title.c_str() );
	  theRatio->Draw(); }	
	else 
	  theRatio->Draw("SAME");
      }      

      JetTrigRatioLegend->Draw();
      temp_canvJetRatio->Print(thePDFFileName.c_str());
    }//end DRAW TRIG Ratios       

    temp_canvJetTrig->Close();        
    temp_canvJetRatio->Close();    
  }//end else statement for het trig plots

  
  

  
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


