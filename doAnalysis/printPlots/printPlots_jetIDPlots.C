#include "printPlots.h"


// code/job switches ------------------------

//options
const bool debugMode=true, doEventCounts=true;

//draw switches
const bool drawJetIDPlots=true;
const bool drawJetConstituentPlots=true, drawDijetPlots=false;


int printPlots_jetIDPlots(const std::string input_condorDir ,
			  const std::string output_PDFname_base ){
  
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=input_condorDir.find("_ak")+3;  
  const std::string radiusInt= input_condorDir.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=input_condorDir.find("Jets");;
  
  //make final jetType strings
  const std::string jetType=input_condorDir.substr( jetTypePos,(jetsPos-jetTypePos) );
  const std::string fullJetType="ak"+radiusInt+jetType;
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  //put together input file strings
  const std::string input_ppData_Filename="HighPtJetTrig_" +fullJetType+ "-allFiles.root";
  const std::string input_ppMC_Filename  ="Py8_CUETP8M1_QCDjetAllPtBins_" +fullJetType+ "-allFiles.root";
  
  
  std::size_t isDataString=input_condorDir.find("ppData");
  //std::size_t isMCString=input_condorDir.find("ppMC");
  bool isData;
  //if(isDataString==std::string::npos ) isData=false;
  //  else if (isDataString != std::string::npos &&
  //	   (isMCString>isDataString) ) isData=false;
  if(isDataString==std::string::npos ) isData=false;
  else isData=true;
  
  
  std::string input_Filename;
  if(isData) input_Filename=input_ppData_Filename;
  else  input_Filename=input_ppMC_Filename;
  const std::string fullFilename=inputDir+input_condorDir+input_Filename;
  
  
  
  
  // OPEN INPUT SECTION
  TFile *fin=NULL;  
  std::cout<<std::endl<<"printing jetID Plots, now opening input files!!"<<std::endl<<std::endl;
  
  //open the input files
  std::cout<<" now opening "<<std::endl<<input_Filename<<std::endl;
  std::cout<<" in directory: "<<inputDir+input_condorDir<<std::endl<<std::endl;
  fin = new TFile(fullFilename.c_str());      
  
  
  // GET OUTPUT PDF FILE READY
  std::string thePDFFileName=outputDir+fullJetType;;
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
    getEventCounts( (TFile*)fin, isData );
    if(isData)theLumi=computeEffLumi( (TFile*) fin);   }
  else {
    std::cout<<"skipping evt/jet QA counts + plots..."<<std::endl<<std::endl;
    theLumi=intgrtdLumi;
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
      printJetIDHist( (TFile*) fin , (int) j, (bool) isData,
		      (std::string) thePDFFileName  , 
		      (std::string) fullJetType, (long double) theLumi  );
    }

  }  








  /*
  // jet ID plots ----------------------
  if(!drawJetIDPlots)std::cout<<std::endl<<"skipping jet QA plots..."<<std::endl<<std::endl;
  else { 
    
    std::cout<<std::endl<<std::endl;    
    if(debugMode)std::cout<<std::endl<<"creating temporary canvas for printing Jet plots..."<<std::endl;
    TCanvas *temp_canvJet = new TCanvas("tempJet", "tempJet", 1000, 800);
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
    TCanvas *temp_canvJet_noLog = new TCanvas("tempJet_noLog", "tempJet_noLog", 1000, 800);
    temp_canvJet_noLog->cd();
    temp_canvJet_noLog->SetLogy(0);
    
    TPad *jetpad_noLog1 = new TPad("jetpad_noLog1", "Overlay Pad",0.0,0.3,1.0,1.0);
    jetpad_noLog1->SetLogy(0);
    jetpad_noLog1->SetBottomMargin(0);
    
    TPad *jetpad_noLog2 = new TPad("jetpad_noLog2", "Ratio Pad",0.0,0.05,1.0,0.3);
    jetpad_noLog2->SetTopMargin(0);
    jetpad_noLog2->SetBottomMargin(0.3);
    
    // PUT JET CUTS SECTION IN HERE
    
    TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
    std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
    std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
    
    TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
    std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
    std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;

    TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
    std::stringstream etaLo; etaLo.precision(1);
    float jtetaLoCut_F=jtetaLoCut_h->GetMean();
    etaLo << std::fixed << jtetaLoCut_F;
    std::string jtetaLoCut_str = etaLo.str();
    std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;

    TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
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
    //TH1F* dijetCut_h= (TH1F*)fin->Get( "" );
    //std::string dijetCut_str = std::to_string((int)dijetCut_h->GetMean()).c_str();
    //std::cout<<"dijetCut_str = "<<dijetCut_str<<std::endl;
    //TH1F* dijetCut_h= (TH1F*)fin->Get( "" );
    //std::string dijetCut_str = std::to_string((int)dijetCut_h->GetMean()).c_str();
    //std::cout<<"dijetCut_str = "<<dijetCut_str<<std::endl;
    //TH1F* dijetCut_h= (TH1F*)fin->Get( "" );
    //std::string dijetCut_str = std::to_string((int)dijetCut_h->GetMean()).c_str();
    //std::cout<<"dijetCut_str = "<<dijetCut_str<<std::endl;

    //hJetEtaCutHi;12.400000
    //hJetEtaCutLo;10.000000
    //hJetQAPtCut;115.000000
    //hLdJetPtCut;155.000000
    //hSubldJetPtCut;120.000000
    //hPtAveCut;125.000000
	    
        
    std::cout<<" drawing jet QA Plots..."<<std::endl;
    for(int j=0;j<N_vars;j++){ 
      
      std::cout<<std::endl;
      
      std::cout<<std::endl<<" var ="<<var[j]<<", j="<<j<<std::endl;
      
      if(!drawJetConstituentPlots && 
	 j>=jetConstits_varStart && j<dijet_varStart){
	if(debugMode) std::cout<<"skipping jet consituent plots..."<<std::endl;
	continue;}
      else if(!drawDijetPlots && 
	      j>=dijet_varStart){
	if(debugMode) std::cout<<"skipping dijet plots..."<<std::endl;
	continue;}
      
      
      // I O
      //std::string inHistName=
      
      std::string nonJetIDHistName="hJetQA_0wJetID_"+var[j];      
      if(debugMode)std::cout<<" opening input non jet ID hist "<<nonJetIDHistName<<std::endl<<std::endl;      
      TH1F* nonJetIDHist=(TH1F*) ( (TH1*)fin->Get(nonJetIDHistName.c_str()) );
      
      std::string theJetIDHistName="hJetQA_1wJetID_"+var[j];      
      if(debugMode)std::cout<<" opening input jetID hist "<<   theJetIDHistName<<std::endl<<std::endl;      
      TH1F* theJetIDHist= (TH1F*) ( (TH1*)fin->Get(theJetIDHistName.c_str()) );
      
      // REBIN AND SCALING
      
      if(!((bool)theJetIDHist) || !((bool)nonJetIDHist) )
	{  std::cout<<"jet ID hist not found. skipping hist."<<std::endl<<std::endl;
	  continue;}
      
      //nonJetIDHist->Scale( 1./theLumi);
      nonJetIDHist->Scale( 1./etaBinWidth);//eta bin width scaling
      
      //theJetIDHist->Scale(1000.);
      theJetIDHist->Scale(1./etaBinWidth);//eta bin width scaling
      
      if(var[j]=="jtpt" || var[j]=="rawpt" ||
	 var[j]=="leadJetPt" || var[j]=="subleadJetPt") {
	nonJetIDHist=(TH1F*)nonJetIDHist->TH1::Rebin(nbins_pt_debug, "nonJetIDHist_rebin_jtPt", ptbins_debug);
	divideBinWidth(nonJetIDHist);      
	
	theJetIDHist=(TH1F*)theJetIDHist->TH1::Rebin(nbins_pt_debug, "wJetIDHist_rebin_jtPt", ptbins_debug);
	divideBinWidth(theJetIDHist);      
      }
      else {
	nonJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
	nonJetIDHist->Scale( 1./nonJetIDHist->GetBinWidth(1) );          
	
	theJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
	theJetIDHist->Scale( 1./theJetIDHist->GetBinWidth(1) );          
	
	
      }
      
      if(debugMode)nonJetIDHist->Print("base");       std::cout<<std::endl;
      if(debugMode)theJetIDHist->Print("base");       std::cout<<std::endl;
      
      // STYLE
      
      nonJetIDHist->SetMarkerStyle(kCircle);
      nonJetIDHist->SetMarkerSize(0.90);
      nonJetIDHist->SetMarkerColor( theDataOverlayMarkerColor );
      nonJetIDHist->SetLineColor( theDataOverlayLineColor );
      
      
      theJetIDHist->SetMarkerStyle(kMultiply);
      theJetIDHist->SetMarkerSize(0.90);
      theJetIDHist->SetMarkerColor( theMCOverlayMarkerColor );
      theJetIDHist->SetLineColor( theMCOverlayLineColor );
      
      
      std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title=AUAxTitle;//crossSectionAxTitle;
      if(debugMode)std::cout<<" x-axis title= "<<h_XAx_Title<<std::endl<<std::endl;
      
      std::string h_Title="jetID Eff Plot for "+var[j];
      
      nonJetIDHist->SetTitle (    h_Title.c_str() );
      nonJetIDHist->SetXTitle( h_XAx_Title.c_str() );
      nonJetIDHist->SetYTitle( h_YAx_Title.c_str() );
      //nonJetIDHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
      
      theJetIDHist->SetTitle (    h_Title.c_str() );
      theJetIDHist->SetXTitle( h_XAx_Title.c_str() );
      theJetIDHist->SetYTitle( h_YAx_Title.c_str() );
      //theJetIDHist->SetYTitle( "#frac{d\sigma}{dp_t} (\{mu}b/bin)" );
      
      
      // variable specific draw orders, canv, etc.
      if(var[j]=="jtpt"||var[j]=="rawpt") {
	
	
	temp_canvJet->cd();			//here is the draw command for pt 
	
	jetpad1->SetBottomMargin(0);
	jetpad1->SetLogy(1);
	jetpad1->SetLogx(1);
	jetpad1->Draw();
	jetpad1->cd();

	theJetIDHist->Draw("E"); 
	nonJetIDHist->Draw("E same");  
	
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


	TLatex *t3=NULL;
	if(isData)t3=new TLatex((t1Loc1),(t1Loc2-.10),sqrts2k15ppString.c_str());
	else t3=new TLatex((t1Loc1),(t1Loc2-.10),Py8TuneString.c_str());
	t3->SetTextFont(63);
	t3->SetTextColor(kBlack);
	t3->SetTextSize(14);
	t3->SetLineWidth(1);
	t3->SetNDC();
	t3->Draw();	
	
	float legx1=0.77, legx2=legx1+0.11;
	float legy1=0.76, legy2=legy1+0.09;
	TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	theJetQALeg->SetFillStyle(0);
	theJetQALeg->AddEntry(nonJetIDHist, "no JetID "  ,"lp");
	theJetQALeg->AddEntry(theJetIDHist, "with jetID"  ,"lp");
	theJetQALeg->Draw(); 
	
	
	temp_canvJet->cd();
	jetpad2->Draw();
	jetpad2->SetLogx(1);
	jetpad2->cd();
	jetpad2->SetTopMargin(0);
	TH1F* theRatio=(TH1F*)theJetIDHist->Clone("MCJetHistClone4Ratio");
	
	theRatio->SetMarkerStyle( kCircle );
	theRatio->SetMarkerSize(0.99 );
	theRatio->SetMarkerColor( kRed );
	theRatio->SetLineColor( theRatioLineColor );
	
	theRatio->SetAxisRange(0.,2.,"Y");
	//theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	
	theRatio->GetYaxis()->SetTitleSize(15);
	theRatio->GetYaxis()->SetTitleFont(43);
	theRatio->GetYaxis()->SetTitleOffset(1.1);
	theRatio->GetYaxis()->SetLabelFont(43); 
	theRatio->GetYaxis()->SetLabelSize(13);
	theRatio->GetXaxis()->SetTitleSize(20);
	theRatio->GetXaxis()->SetTitleFont(43);
	theRatio->GetXaxis()->SetTitleOffset(4.);
	theRatio->GetXaxis()->SetLabelFont(43); 
	theRatio->GetXaxis()->SetLabelSize(15);
	
	theRatio->SetTitle("");
	theRatio->SetXTitle( h_XAx_Title.c_str() );
	theRatio->SetYTitle(jetIDratioTitle.c_str() );
	theRatio->Divide(nonJetIDHist);
	theRatio->Draw("HIST E");
	
	TLine* lineAtOne          = NULL;
	TLine* lineAtOneHalf      = NULL;
	TLine* lineAtOneEtOneHalf = NULL;

	lineAtOne          = new TLine(ptbins_debug[0],1.0,jetQAxmax[j],1.0); 
	lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,jetQAxmax[j],0.5); 
	lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,jetQAxmax[j],1.5); 
	
	lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	
	lineAtOne->Draw("same");
	lineAtOneHalf->Draw("same");
	lineAtOneEtOneHalf->Draw("same");
	
	
	temp_canvJet->Print( thePDFFileName.c_str() ); 

	jetpad1->SetLogx(0);
	jetpad2->SetLogx(0);    }
	
	else if(var[j]=="jteta" ||
		var[j]=="jtphi"	||
		( j>=jetConstits_varStart && j<dijet_varStart )){ 

		
	  temp_canvJet->cd();			//here is the draw command for pt 
	  jetpad1->SetBottomMargin(0);
	  jetpad1->SetLogy(1);
	  jetpad1->Draw();
	  jetpad1->cd();

	  
	  nonJetIDHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  theJetIDHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	

	  nonJetIDHist->Draw("E"); 
	  theJetIDHist->Draw("E SAME"); 

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

	  TLatex *t3=NULL;
	  if(isData)t3=new TLatex((t1Loc1),(t1Loc2-.10),sqrts2k15ppString.c_str());
	  else t3=new TLatex((t1Loc1),(t1Loc2-.10),Py8TuneString.c_str());
	  t3->SetTextFont(63);
	  t3->SetTextColor(kBlack);
	  t3->SetTextSize(14);
	  t3->SetLineWidth(1);
	  t3->SetNDC();
	  t3->Draw();	

	  float legx1=0.77, legx2=legx1+0.11;
	  float legy1=0.76, legy2=legy1+0.09;
	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	  theJetQALeg->SetFillStyle(0);
	  theJetQALeg->AddEntry(nonJetIDHist,"no JetID ","lp");
	  theJetQALeg->AddEntry(theJetIDHist,"with jetID"  ,"lp");
	  theJetQALeg->Draw();

	  //Ratio Plots
	  //temp_canvJet->cd();	  
	  jetpad2->cd();
	  jetpad2->Draw();
	  TH1F* theRatio=(TH1F*)theJetIDHist->Clone("mcJetHistClone4Ratio");
	  
	  theRatio->SetMarkerColor( kRed );
	  theRatio->SetLineColor( theRatioLineColor );
	  
	  theRatio->SetAxisRange(0.,2.,"Y");
	  
	  theRatio->GetYaxis()->SetTitleSize(15);
	  theRatio->GetYaxis()->SetTitleFont(43);
	  theRatio->GetYaxis()->SetTitleOffset(1.1);
	  theRatio->GetYaxis()->SetLabelFont(43); 
	  theRatio->GetYaxis()->SetLabelSize(13);
	  theRatio->GetXaxis()->SetTitleSize(20);
	  theRatio->GetXaxis()->SetTitleFont(43);
	  theRatio->GetXaxis()->SetTitleOffset(4.);
	  theRatio->GetXaxis()->SetLabelFont(43); 
	  theRatio->GetXaxis()->SetLabelSize(15);
	  
	  theRatio->SetTitle("");
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle(jetIDratioTitle.c_str());
	  theRatio->Divide(nonJetIDHist);
	  theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  theRatio->Draw("E");
	  
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
	  jetpad_noLog1->SetGridx(1);
	  jetpad_noLog1->SetGridy(1);
	  jetpad_noLog1->Draw();
	  jetpad_noLog1->cd();


	  nonJetIDHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  theJetIDHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	

	  nonJetIDHist->Draw("E"); 
	  theJetIDHist->Draw("E SAME"); 
	  
	  float t1Loc1=0.15, t1Loc2=0.7; 	    	  
	  if(var[j]=="Aj") t1Loc1=.7;

	  
	  TLatex *t1=new TLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
	  t1->SetTextFont(63);
	  t1->SetTextColor(kBlack);
	  t1->SetTextSize(14);
	  t1->SetLineWidth(1);
	  t1->SetNDC();
	  t1->Draw();
	  
	  TLatex *t2=new TLatex((t1Loc1),(t1Loc2-.04),"p_{t}^{ld} > 60 GeV");
	  t2->SetTextFont(63);
	  t2->SetTextColor(kBlack);
	  t2->SetTextSize(14);
	  t2->SetLineWidth(1);
	  t2->SetNDC();
	  t2->Draw();	
	  
	  TLatex *t3=new TLatex((t1Loc1),(t1Loc2-.08),"p_{t}^{subld} > 30GeV");
	  t3->SetTextFont(63);
	  t3->SetTextColor(kBlack);
	  t3->SetTextSize(14);
	  t3->SetLineWidth(1);
	  t3->SetNDC();
	  t3->Draw();	
	  
	  TLatex *t4=new TLatex(t1Loc1, (t1Loc2-0.13), "#Delta#phi > #frac{2}{3}#pi");	    
	  t4->SetTextFont(63);
	  t4->SetTextColor(kBlack);
	  t4->SetTextSize(14);
	  t4->SetLineWidth(1);
	  t4->SetNDC();
	  t4->Draw();	

	  TLatex *t5=NULL;
	  if(isData)t5=new TLatex((t1Loc1),(t1Loc2-.18),sqrts2k15ppString.c_str());
	  else t5=new TLatex((t1Loc1),(t1Loc2-.18),Py8TuneString.c_str());
	  t5->SetTextFont(63);
	  t5->SetTextColor(kBlack);
	  t5->SetTextSize(14);
	  t5->SetLineWidth(1);
	  t5->SetNDC();
	  t5->Draw();	

	  float legx1=0.15, legx2=legx1+0.11;
	  if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
	  float legy1=0.74, legy2=legy1+0.09;

	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	  theJetQALeg->SetFillStyle(0);
	  theJetQALeg->AddEntry(nonJetIDHist,"no JetID ","lp");
	  theJetQALeg->AddEntry(theJetIDHist,"with jetID"  ,"lp");
	  theJetQALeg->Draw(); 

	  //Ratio Plots
	  temp_canvJet_noLog->cd();	  
	  jetpad_noLog2->Draw();
	  jetpad_noLog2->cd();
	  TH1F* theRatio=(TH1F*)theJetIDHist->Clone("MCJetHistClone4Ratio");
	  //I need to redefine the axis titles here
	  std::string h_XAx_ratTitle=var_xAx_Titles[j], h_YAx_ratTitle="Data/MC"    ;
	  	  
	  theRatio->SetMarkerColor( kRed );
	  theRatio->SetLineColor( theRatioLineColor );
	  
	  theRatio->SetAxisRange(0.,2.,"Y");
	  
	  theRatio->GetYaxis()->SetTitleSize(15);
	  theRatio->GetYaxis()->SetTitleFont(43);
	  theRatio->GetYaxis()->SetTitleOffset(1.1);
	  theRatio->GetYaxis()->SetLabelFont(43); 
	  theRatio->GetYaxis()->SetLabelSize(13);
	  theRatio->GetXaxis()->SetTitleSize(20);
	  theRatio->GetXaxis()->SetTitleFont(43);
	  theRatio->GetXaxis()->SetTitleOffset(4.);
	  theRatio->GetXaxis()->SetLabelFont(43); 
	  theRatio->GetXaxis()->SetLabelSize(15);
	  
	  theRatio->SetTitle("");
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle(jetIDratioTitle.c_str());
	  theRatio->Divide(nonJetIDHist);
	  theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  theRatio->Draw("E");
	  
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
	  jetpad1->SetLogy(1);
	  jetpad1->SetLogx(1);
	  jetpad1->Draw();
	  jetpad1->cd();

	  nonJetIDHist->Draw("E"); 
	  theJetIDHist->Draw("E SAME"); 

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
	  
	  TLatex *t4=new TLatex(t1Loc1, (t1Loc2-0.16), "#Delta#phi > #frac{2}{3}#pi");	    
	  t4->SetTextFont(63);
	  t4->SetTextColor(kBlack);
	  t4->SetTextSize(14);
	  t4->SetLineWidth(1);
	  t4->SetNDC();
	  t4->Draw();	

	  TLatex *t5=NULL;
	  if(isData)t5=new TLatex((t1Loc1),(t1Loc2-.23),sqrts2k15ppString.c_str());
	  else t5=new TLatex((t1Loc1),(t1Loc2-.23),Py8TuneString.c_str());
	  t5->SetTextFont(63);
	  t5->SetTextColor(kBlack);
	  t5->SetTextSize(14);
	  t5->SetLineWidth(1);
	  t5->SetNDC();
	  t5->Draw();	

	  float legx1=0.70, legx2=legx1+0.11;
	  float legy1=0.76, legy2=legy1+0.09;
	  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
	  theJetQALeg->SetFillStyle(0);
	  theJetQALeg->AddEntry(nonJetIDHist,"no JetID ","lp");
	  theJetQALeg->AddEntry(theJetIDHist,"with jetID"  ,"lp");
	  theJetQALeg->Draw(); 
	  
	  //Ratio Plots
	  temp_canvJet->cd();	  
	  jetpad2->SetLogx(1);
	  jetpad2->Draw();
	  jetpad2->cd();
	  TH1F* theRatio=(TH1F*)theJetIDHist->Clone("MCJetHistClone4Ratio");
	  
	  theRatio->SetMarkerColor( kRed );
	  theRatio->SetLineColor( theRatioLineColor );
	  
	  theRatio->SetAxisRange(0.,2.,"Y");
	  //theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
	  
	  theRatio->GetYaxis()->SetTitleSize(15);
	  theRatio->GetYaxis()->SetTitleFont(43);
	  theRatio->GetYaxis()->SetTitleOffset(1.1);
	  theRatio->GetYaxis()->SetLabelFont(43); 
	  theRatio->GetYaxis()->SetLabelSize(13);
	  theRatio->GetXaxis()->SetTitleSize(20);
	  theRatio->GetXaxis()->SetTitleFont(43);
	  theRatio->GetXaxis()->SetTitleOffset(4.);
	  theRatio->GetXaxis()->SetLabelFont(43); 
	  theRatio->GetXaxis()->SetLabelSize(15);
	  
	  theRatio->SetTitle ("");
	  theRatio->SetXTitle( h_XAx_Title.c_str() );
	  theRatio->SetYTitle( jetIDratioTitle.c_str() );
	  
	  theRatio->Divide(nonJetIDHist);
	  theRatio->Draw("E");
	  
	  TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,jetQAxmax[j],1.0); 
	  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	  TLine* lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,jetQAxmax[j],0.5); 
	  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	  TLine* lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,jetQAxmax[j],1.5); 
	  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	  
	  lineAtOne->Draw("same");
	  lineAtOneHalf->Draw("same");
	  lineAtOneEtOneHalf->Draw("same");
	  
	  temp_canvJet->Print( thePDFFileName.c_str() ); 
	  jetpad1->SetLogx(0);
	  jetpad2->SetLogx(0);
	}
      
      
      std::cout<<std::endl;	    	
      

    }//end loop over jetQA Plots     

    temp_canvJet_noLog->Close();//clean up neatly
    //    temp_canvJetRat->Close();
    temp_canvJet->Close();
  }
  */  


  
  // for DEBUG ONLY
  if(debugMode)std::cout<<std::endl<<"closing the PDF file"<<std::endl;
  TCanvas *temp_canvClose = new TCanvas("tempClose", "tempClose", 1000, 800);
  temp_canvClose->Print( close_thePDFFileName.c_str() );  
  temp_canvClose->Close();    
  

  if(debugMode)std::cout<<std::endl<<"closing input files"<<std::endl;
  fin->Close();  
  
  return 0;
}// end printplots
				     
				     


// steering ----------------------

int main(int argc, char *argv[]){
  
  int rStatus = -1;
  if( argc!=3 ) {//no input arguments, error
    std::cout<<"do ./printPlots_jetIDPlots.exe <inputCondorDir> <outputNameBase>"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  
  if(argc==3) rStatus=printPlots_jetIDPlots( (const std::string) argv[1], (const std::string) argv[2] );
  
  std::cout<<"done, rStatus="<<rStatus<<std::endl;
  return rStatus;
   
}
