




void printEvtVtxQAHist( TFile* finData , std::string inDataHistName, 
			TFile* finMC   , std::string inMCHistName  ,
			int rebin, std::string thePDFFileName, long double theLumi , TFile* fout=NULL) {
  
  if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("vtxCanv"+inDataHistName).c_str(), "temp Evt Canv withLog", CANVX,CANVY);
  temp_canvEvt->cd();
  
  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  
  evtpad1->SetBottomMargin(0); 
  
  evtpad1->Draw();
  evtpad1->cd();
  
  // open hist
  std::cout<<" getting hist(s): "<<inDataHistName<<std::endl;
  
  TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
  if(!theDataEvtQAHist){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  theDataEvtQAHist=(TH1F*)theDataEvtQAHist->TH1::Rebin( (int)rebin, (inDataHistName+"_rebin").c_str());
  theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
  theDataEvtQAHist->Scale( 1./theLumi );
  dataHistStyle(theDataEvtQAHist);
  
  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
  if(!theMCEvtQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  theMCEvtQAHist=(TH1F*)theMCEvtQAHist->TH1::Rebin( (int)rebin, (inMCHistName+"_rebin").c_str());
  theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
  theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
  MCHistStyle(theMCEvtQAHist);
  
  


  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  std::string h_XAx_Title;// ="v_{z}^{evt} (cm)";
  float vxyr_ymin=1.4e-02, vxyr_ymax=1.e+06;
  
  if( inDataHistName.find("Vz") != std::string::npos ){
    h_XAx_Title ="v_{z}^{evt} (cm)";}
  else if( inDataHistName.find("Vx") != std::string::npos ){
    evtpad1->SetLogy(1);
    h_XAx_Title ="v_{x}^{evt} (cm)";
    theDataEvtQAHist->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
    theMCEvtQAHist  ->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
  }
  else if( inDataHistName.find("Vy") != std::string::npos ){
    evtpad1->SetLogy(1);
    h_XAx_Title ="v_{y}^{evt} (cm)";
    theDataEvtQAHist->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
    theMCEvtQAHist  ->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
  }
  else if( inDataHistName.find("Vr") != std::string::npos ){
    evtpad1->SetLogy(1);
    h_XAx_Title ="v_{r}^{evt} (cm)";
    theDataEvtQAHist->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
    theMCEvtQAHist  ->SetAxisRange(vxyr_ymin,vxyr_ymax,"Y");
  }
  std::string h_YAx_Title ="A.U.";
  



  
  Double_t dataMax = theDataEvtQAHist->GetMaximum(), MCMax=theMCEvtQAHist->GetMaximum();
  std::cout<<"dataMax="<<dataMax<<std::endl;
  std::cout<<"MCMax="<<    MCMax<<std::endl;
  if(dataMax>MCMax){
    std::cout<<"dataMax > MC Max"<<std::endl;
    theDataEvtQAHist->SetTitle (    h_Title.c_str() );
    theDataEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
    theDataEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
    
    theDataEvtQAHist->Draw("E"); 
    theMCEvtQAHist->Draw("HIST E SAME"); 
    theDataEvtQAHist->Draw("E SAME"); //because i want data points written over MC's red line hist
  }
  else {
    theMCEvtQAHist->SetTitle (    h_Title.c_str() );
    theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
    theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
    //theMCEvtQAHist->SetAxisRange(.001,10.*MCMax,"Y");
    theMCEvtQAHist->Draw("HIST E"); 
    theDataEvtQAHist->Draw("E SAME"); 
  }
  
  
  // legend
  float legx1=0.65, legx2=legx1+0.08;
  float legy1=0.74, legy2=legy1+0.09;

  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->SetBorderSize(0.);

  theEvtQALeg->AddEntry(theDataEvtQAHist, " data" ,"lp");
  theEvtQALeg->AddEntry(theMCEvtQAHist,   " MC"       ,"lp");
  theEvtQALeg->Draw();
  
  // text
  float pp1x=0.14,pp1y=0.82;
  TLatex *pp1=makeTLatex(pp1x,pp1y,promptRecoPDString.c_str());  pp1->Draw();
  TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),PDStatsString.c_str());  pp2->Draw();  
  TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();
  


  //ratios go here
  temp_canvEvt->cd(); //change back to main canvas
  
  evtpad2->Draw();
  evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  TH1F *theRatio=(TH1F*)theMCEvtQAHist->Clone("MCClone4Ratio");
  TH1F *theDenom=(TH1F*)theDataEvtQAHist->Clone("DataClone4Ratio");
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
  theRatio->SetYTitle( ratioTitle.c_str() );

  theRatio->Divide(theDenom);
  theRatio->SetAxisRange(0.0,2.0, "Y");  
  
  theRatio->Draw("HIST E");
  
  float min = theRatio->GetBinLowEdge( 1 );
  float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
    + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  
  TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  lineAtOne->Draw("same");
  lineAtOneHalf->Draw("same");
  lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvEvt->DrawClone();
    //outcanv->SetTitle(("Evt Vtx QA Plots Canvas "+inDataHistName).c_str());
    //outcanv->Write(("EvtVtxQAPlots_"+inDataHistName+"_canv").c_str());    
    outcanv->SetTitle(("Evt Vtx QA Plots Canvas "+inMCHistName).c_str());
    outcanv->Write(("EvtVtxQAPlots_"+inMCHistName+"_canv").c_str());    
    
  }
  
  temp_canvEvt->Close();
  return;
}

































void printEvtNrefQAHist( TFile* finData , std::string inDataHistName, 
			 TFile* finMC   , std::string inMCHistName  ,
			 std::string thePDFFileName, std::string fullJetType, 
			 long double theLumi ) {
  
  if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas("tempEvt", "temp Evt Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  
  evtpad1->SetBottomMargin(0); 
  evtpad1->SetLogy(1);
  evtpad1->Draw();
  evtpad1->cd();
  
  if(inDataHistName=="hLeadJetPt"){
    evtpad1->SetGridx(1);
    evtpad1->SetGridy(1);
    evtpad1->SetLogx(1);

    evtpad2->SetGridx(1);
    evtpad2->SetGridy(0);
    evtpad2->SetLogx(1);
  }


  

  // PUT JET CUTS SECTION IN HERE
  TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
  std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
  
  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{gen}>"+genjtptCut_str+" GeV, "+"p_{T}>"+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;

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
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  if(inMCHistName=="hWNref"||inMCHistName=="hNref" ){
    jetEtaCutString="|y|<5.1";
  }
  if(inMCHistName=="hWjetsPEvt"||inMCHistName=="hjetPEvt" ){
    jetEtaCutString="|y|<4.7";
  }

  
  
  // open hist
  std::cout<<" drawing "<<inDataHistName<<std::endl<<std::endl;
  TH1F* theDataEvtQAHist= (TH1F*)finData->Get( inDataHistName.c_str() );
  if(!theDataEvtQAHist){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  //theDataEvtQAHist->Scale( 1./theDataEvtQAHist->GetBinWidth(1) );
  dataHistStyle(theDataEvtQAHist);


  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
  if(!theMCEvtQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  //theMCEvtQAHist->Scale( 1./theMCEvtQAHist->GetBinWidth(1) );
  //theMCEvtQAHist->SetMarkerStyle(kMultiply);
  //theMCEvtQAHist->SetMarkerSize(0.90);
  //theMCEvtQAHist->SetMarkerColor( theMCOverlayMarkerColor );
  theMCEvtQAHist->SetLineColor( theMCLineColor );  
  
  //hist-type specific scaling
  if(inDataHistName!=inMCHistName)std::cout<<"WARNING: inDataHistName!=inMCHistName!=\"hLeadJetPt\"!!!!!!"<<std::endl;
  
  if(inMCHistName!="hLeadJetPt"){
    //if(inMCHistName=="hWjetsPEvt"||inMCHistName=="hjetsPEvt") 
    theDataEvtQAHist->SetBinContent(1, 0.);
    theMCEvtQAHist->SetBinContent(1, 0.);


    theMCEvtQAHist->Scale( theDataEvtQAHist->Integral(1, (int)theMCEvtQAHist->GetNbinsX() )/
			   theMCEvtQAHist->Integral(   1, (int)theMCEvtQAHist->GetNbinsX() ) 
			   );	
      //else theMCEvtQAHist->Scale( theDataEvtQAHist->Integral()/theMCEvtQAHist->Integral() );	
  }
  else{// if(inMCHistName=="hLeadJetPt"){
    theDataEvtQAHist=(TH1F*)theDataEvtQAHist->
      TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_evtLeadJetPt", ptbins_debug);
    divideBinWidth(theDataEvtQAHist);
    //theDataEvtQAHist->Scale( 1./theLumi );
    
    theMCEvtQAHist=(TH1F*)theMCEvtQAHist->
      TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_evtLeadJetPt", ptbins_debug);
    divideBinWidth(theMCEvtQAHist);
    //theMCEvtQAHist->Scale(1000.);
    theMCEvtQAHist->Scale( theDataEvtQAHist->Integral("width")/
			   theMCEvtQAHist->Integral("width") );	
  }
  

  

  // title+axes
  std::string h_Title=PDStatsString_2;
  std::string h_XAx_Title;
  std::string h_YAx_Title ="A.U.";
  
  if(inMCHistName=="hLeadJetPt") h_XAx_Title="Event Maximum Jet p_{T}^{RECO} (GeV)";
  else h_XAx_Title="N_{jets}/Evt";
  
  theDataEvtQAHist->SetTitle( h_Title.c_str() );
  theMCEvtQAHist->SetTitle( h_Title.c_str() );
  theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  theDataEvtQAHist->Draw("E"); 
  theMCEvtQAHist->Draw("HIST E SAME"); 
  theDataEvtQAHist->Draw("E SAME"); 

  // legend
  float legx1=0.74, legx2=legx1+0.11;
  float legy1=0.74, legy2=legy1+0.09;

  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");  
  theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->SetBorderSize(0.);

  theEvtQALeg->AddEntry(theDataEvtQAHist, "Data" ,"lp");
  theEvtQALeg->AddEntry(theMCEvtQAHist,   "P8 MC"       ,"lp");
  theEvtQALeg->Draw();
  
  // text
  float pp1x=0.48,pp1y=0.84;
  TLatex *pp1=makeTLatex(  pp1x,pp1y,"pp 2015 promptReco, Jet80+LowerJets" );  pp1->Draw();
  TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),"L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV");  pp2->Draw();
  TLatex *pp3=makeTLatex(pp1x,(pp1y-0.08),"Py8 Tune CUETP8M1" );  pp3->Draw();
  TLatex *pp4=makeTLatex(pp1x,(pp1y-0.12), fullJetType+"Jets" );  pp4->Draw();
  TLatex *pp5=makeTLatex(pp1x,(pp1y-0.16), jetCutString );  pp5->Draw();
  TLatex *pp6=makeTLatex(pp1x,(pp1y-0.20), jetEtaCutString );  pp6->Draw();
  
  //ratios go here
  temp_canvEvt->cd(); //change back to main canvas
  evtpad2->Draw();
  evtpad2->cd();
  
  //hist for ratio; use clone

  //TH1F *theRatio=(TH1F*)theDataEvtQAHist->Clone("DataVzClone4Ratio");
  TH1F *theRatio=(TH1F*)theMCEvtQAHist->Clone("MCClone4Ratio");
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
  theRatio->SetYTitle( ratioTitle.c_str() );
  theRatio->Divide(theDataEvtQAHist);
  theRatio->SetAxisRange(0.0,2.0, "Y");  
  
  theRatio->Draw("HIST E");
  
  float min = theRatio->GetBinLowEdge( 1 );
  float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
    + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  
  TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  lineAtOne->Draw("same");
  lineAtOneHalf->Draw("same");
  lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  return;
}









void printMCEvtQAHist( TFile* finMC   , std::string inMCHistName , std::string thePDFFileName  , TFile* fout=NULL){
  
  if(!finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing Evt plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("tempEvt_"+inMCHistName).c_str(), "temp Evt Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->SetLogy(1);
  temp_canvEvt->SetLogx(1);
  temp_canvEvt->cd();  


  std::cout<<" drawing "<<inMCHistName<<std::endl<<std::endl;  
  TH1F* theMCEvtQAHist= (TH1F*)finMC->Get( inMCHistName.c_str() );
  if(!theMCEvtQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  MCHistStyle(theMCEvtQAHist);
  
  std::string h_Title=PDStatsString_2;
  std::string h_XAx_Title="MC evt #hat{p}_{t} (GeV)"       ;
  std::string h_YAx_Title="A.U."       ;

  theMCEvtQAHist->SetTitle( h_Title.c_str() );
  theMCEvtQAHist->SetXTitle( h_XAx_Title.c_str() );
  theMCEvtQAHist->SetYTitle( h_YAx_Title.c_str() );
  theMCEvtQAHist->SetAxisRange(15.,1000.,"X");
  
  
  if (inMCHistName=="hWeightedpthat") 
    theMCEvtQAHist->SetAxisRange(1.e-7,5.e+5,"Y");
  else//pthat, unweighted
    theMCEvtQAHist->SetAxisRange(1.e+1,5.e+5,"Y");
  
  theMCEvtQAHist->Draw("HIST E");
  
  // test
  float t1x=.65,t1y=.81;
  TLatex *t1;
  if (inMCHistName=="hWeightedpthat") {
    t1=makeTLatex(t1x,t1y,"QCD MC, weighted");
  }
  else{
    t1=makeTLatex(t1x,t1y,"QCD MC, unweighted");
  }
  t1->Draw();	
  
  TLatex *t2=makeTLatex(t1x,t1y-0.04,"Py8 Tune CUETP8M1");
  t2->Draw();	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); 
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvEvt->DrawClone();
    outcanv->SetTitle(("MC Evt QA Plots Canvas "+inMCHistName).c_str());
    outcanv->Write(("MCEvtQAPlots_"+inMCHistName+"_canv").c_str());    
    //temp_canvEvt->Write();
  }
  
  temp_canvEvt->Close();  
  return;
}




void printJetQAHist( TFile* finData , TFile* finMC, int j, bool doJetIDPlots, int etabin,
		     std::string inHistName , std::string thePDFFileName , std::string fullJetType, 
		     long double theLumi  ,TFile* fout=NULL) {

  
  std::string jetIDInt;
  if(doJetIDPlots)jetIDInt="1";
  else jetIDInt="0";
  
  
  if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  
  
  // STRINGS FOR HISTS
  TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
  std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
  
  std::string genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
  std::cout<<"genJetCutString="<<genJetCutString<<std::endl;
  
  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)finData->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)finData->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
//  TH1F* jtetaHiCut_h= (TH1F*)finData->Get( "hJetEtaCutHi" );
//  std::stringstream etaHi; etaHi.precision(1);
//  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
//  etaHi << std::fixed << jtetaHiCut_F;
//  std::string jtetaHiCut_str = etaHi.str();
//  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString;
  if(etabin==0) jetEtaCutString = "0.0 < #||{y} < 0.5";
  if(etabin==1) jetEtaCutString = "0.5 < #||{y} < 1.0";
  if(etabin==2) jetEtaCutString = "1.0 < #||{y} < 1.5";
  if(etabin==3) jetEtaCutString = "1.5 < #||{y} < 2.0";
  if(etabin==4) jetEtaCutString = "2.0 < #||{y} < 2.5";
  if(etabin==5) jetEtaCutString = "2.5 < #||{y} < 3.0";
//  if(jtetaLoCut_str=="0.0") jetEtaCutString="|y| < "+jtetaHiCut_str;
//  else jetEtaCutString=jtetaLoCut_str+" < |y| < "+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  //float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  float etaBinWidth=1.;//2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  
  
  
  
  // OPEN FILES
  std::string inDataHist=inHistName;
  std::cout<<" opening input Data "<<inDataHist<<std::endl<<std::endl;      
  
  TH1F* theDataJetQAHist=(TH1F*) ( (TH1*)finData->Get(inDataHist.c_str()) );
  theDataJetQAHist=(TH1F*) (theDataJetQAHist->Clone((inDataHist+"_"+std::to_string(j)+"clone").c_str()) );
  if(!theDataJetQAHist){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  
  theDataJetQAHist->Scale( 1./theLumi);
  theDataJetQAHist->Scale( 1./etaBinWidth);//eta bin width scaling
  
  
  std::string inMCHist=inHistName;
  float MCscale4Units=1./etaBinWidth;
  if(var[j]=="jtpt_forRes" || var[j]=="jtpt_L2Res" || var[j]=="jtpt_L3Res"){
    inMCHist="hJetQA_"+std::to_string((int)doJetIDPlots)+"wJetID_jtpt";
    //MCscale4Units=1.;
  }
  
  std::cout<<" opening input MC "<<inMCHist<<std::endl<<std::endl;      
  TH1F* theMCJetQAHist= (TH1F*) ( (TH1*)finMC->Get(inMCHist.c_str()) );  
  theMCJetQAHist=(TH1F*) (theMCJetQAHist->Clone((inMCHist+"_"+std::to_string(j)+"clone").c_str()) );
  if(!theMCJetQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  theMCJetQAHist->Scale(MCscale4Units);
  
  
  
  // MORE SCALING
  if(var[j]=="jtpt" || var[j]=="rawpt" ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res" ||
     var[j]=="leadJetPt" || var[j]=="subleadJetPt") {
    
    theDataJetQAHist=(TH1F*)theDataJetQAHist->TH1::Rebin(nbins_pt_debug, ("ppData_analysisBins_"+var[j]).c_str(), ptbins_debug);
    divideBinWidth(theDataJetQAHist);      
    
    theMCJetQAHist=(TH1F*)theMCJetQAHist->TH1::Rebin(nbins_pt_debug, ("ppMC_analysisBins_"+var[j]).c_str(), ptbins_debug);
    divideBinWidth(theMCJetQAHist); 
      
  }
  
  else {
    
    theDataJetQAHist->TH1::Rebin(var_xAx_reBin[j]);
    theDataJetQAHist->Scale( 1./theDataJetQAHist->GetBinWidth(1) );          
    theDataJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
    
    theMCJetQAHist->TH1::Rebin(var_xAx_reBin[j]);
    theMCJetQAHist->Scale( 1./theMCJetQAHist->GetBinWidth(1) );          
    theMCJetQAHist->Scale( theDataJetQAHist->Integral()/theMCJetQAHist->Integral() );	        
    theMCJetQAHist->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
    
  }
  
  if(funcDebug)theDataJetQAHist->Print("base");       std::cout<<std::endl;  
  if(funcDebug)theMCJetQAHist->Print("base");       std::cout<<std::endl;
  


  
  // STYLE  
  dataHistStyle(theDataJetQAHist);
  MCHistStyle(theMCJetQAHist);
  
  
  
  // TITLES
  std::string h_XAx_Title= var_xAx_Titles[j];
  //std::string h_YAx_Title= dcrossSectionAxTitle;
  //std::string h_YAx_Title= AUAxTitle;//for APS DNP
  std::string h_YAx_Title= ddcrossSectionAxTitle;//for APS DNP
  std::string h_Title = PDStatsString_2; //for APS DNP
  
  theDataJetQAHist->SetTitle (    h_Title.c_str() );
  theDataJetQAHist->SetXTitle( h_XAx_Title.c_str() );
  theDataJetQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  theMCJetQAHist->SetTitle (    h_Title.c_str() );
  theMCJetQAHist->SetXTitle( h_XAx_Title.c_str() );
  theMCJetQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  
  
  TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone(  ("MCJetHistClone4Ratio_"+var[j]).c_str());  
  
  
  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvJet->cd();
  
    
  TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad1->SetGridy(0);
  jetpad1->SetGridx(0);
  jetpad1->SetLogy(1);


  TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  jetpad2->SetGridx(1);
  jetpad2->SetLogy(0);
  
  
  if(var[j]=="jtpt"      || var[j]=="rawpt"          ||
     var[j]=="jtpt_forRes"|| var[j]=="jtpt_L2Res"|| var[j]=="jtpt_L3Res"|| 
     var[j]=="leadJetPt" || var[j]=="subleadJetPt"     ){
    
    jetpad1->SetLogx(1);  
    jetpad2->SetLogx(1);  
  }
  else {
    jetpad1->SetLogx(0);
    jetpad2->SetLogx(0);
  }
  
  jetpad1->SetBottomMargin(0);
  
  jetpad2->SetTopMargin(0);
  jetpad2->SetBottomMargin(0.3);
  
  //if(var[j]=="jtpt"      || var[j]=="rawpt"          ||
  //   var[j]=="leadJetPt" || var[j]=="subleadJetPt"     )
  //  jetpad2->SetLogx(1);  
  //else jetpad2->SetLogx(0);
  
  
  
  
  jetpad1->Draw();
  jetpad2->Draw();  
  

  
  // variable specific draw orders, canv, etc.
  if(var[j]=="jtpt"||var[j]=="rawpt" ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      )         {
    
    jetpad1->cd();
    if(var[j]=="jtpt"||var[j]=="rawpt"){
      //theMCJetQAHist->Draw("HIST E"); 
      //theDataJetQAHist->Draw("E SAME"); 
      theDataJetQAHist->Draw("E"); 
      theMCJetQAHist->Draw("HIST E SAME ]["); 
 
    } 
    
    else if(var[j]=="jtpt_forRes"||var[j]=="jtpt_L2Res"||var[j]=="jtpt_L3Res"){
      theDataJetQAHist->Draw("E"); 
      theMCJetQAHist->Draw("HIST E SAME"); 
    } 
    
    float t1Loc1=0.50, t1Loc2=0.84;
    TLatex* t1= makeTLatex(  t1Loc1  ,  t1Loc2      ,  fullJetType+"Jets" );    t1->Draw();
    TLatex *t2= makeTLatex( (t1Loc1) , (t1Loc2-.05) , (jetCutString)      );    t2->Draw();	
    TLatex *t3= makeTLatex( (t1Loc1) , (t1Loc2-.10) , (genJetCutString)      ); t3->Draw();	
    TLatex *t4= makeTLatex( (t1Loc1) , (t1Loc2-.15) , (jetEtaCutString)   );    t4->Draw();	    
    if(doJetIDPlots){
      TLatex *t5= makeTLatex( (t1Loc1), (t1Loc2-.20), "Jet ID Applied" );      t5->Draw();	}
    
    float legx1=0.70, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->SetBorderSize(0.);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw(); 
    
    
    jetpad2->cd();
    
    
    TH1F* theDenom=(TH1F*)theDataJetQAHist->Clone(("DataHistClone4Ratio_"+var[j]).c_str());
    
    theRatio->SetLineColor( theMCLineColor );
    
    theRatio->SetAxisRange(0.5,1.5,"Y");
    
    theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);
    
    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(18);
    
    theRatio->SetTitle("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle(ratioTitle.c_str() );
    
    theRatio->Divide(theDenom);
    
    theRatio->GetXaxis()->SetMoreLogLabels(true);
    theRatio->GetXaxis()->SetNoExponent(true);
    
    theRatio->Draw("HIST E");
    
    TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,ptbins_debug[nbins_pt_debug],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,ptbins_debug[nbins_pt_debug],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,ptbins_debug[nbins_pt_debug],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
    
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");



    
    
  }
  
  else if(var[j]=="jteta" ||
	  var[j]=="jtphi"	||
	  ( j>=jetConstits_varStart && j<dijet_varStart )){ 
    
    jetpad1->cd();
    //if(var[j]=="jtphi"){
      //theMCJetQAHist->SetAxisRange(  .02,.8,"Y");
      //theDataJetQAHist->SetAxisRange(.02,.8,"Y");
    //};
    
    theMCJetQAHist->Draw("HIST E"); 
    theDataJetQAHist->Draw("E same");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    //float t1Loc1=0.50, t1Loc2=0.72; 	    
    
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(genJetCutString).c_str());    t2->Draw();
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();		
    if(doJetIDPlots){
      TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.20),"Jet ID Applied");      t5->Draw();	
    }
    
    float legx1=0.70, legx2=legx1+0.09;
    float legy1=0.76, legy2=legy1+0.09;

    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->SetBorderSize(0.);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw();
    
    //Ratio Plots
    
    jetpad2->cd();

    //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone("mcJetHistClone4Ratio");
    //    theRatio->SetMarkerStyle( kMultiply );
    //    theRatio->SetMarkerSize(0.99 );
    //    theRatio->SetMarkerColor( theMCColor );
    theRatio->SetLineColor( theMCLineColor );
    //theRatio->SetLineWidth( 2 );

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
    theRatio->GetXaxis()->SetLabelSize(18);
    
    theRatio->SetTitle("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle(ratioTitle.c_str());
    theRatio->Divide(theDataJetQAHist);
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
    
    
  }
  
  else if(var[j]=="Aj"||
	  var[j]=="xj"||
	  var[j]=="dphi"){ 
    
    jetpad1->cd();
    theMCJetQAHist->Draw("HIST E"); 
    theDataJetQAHist->Draw("E SAME"); 
    
    float t1Loc1=0.15, t1Loc2=0.7; 	    	  
    if(var[j]=="Aj") t1Loc1=.7;
    
    
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),"p_{T}^{ld} > 60 GeV");    t2->Draw();	
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),"p_{T}^{subld} > 30GeV");    t3->Draw();	
    TLatex *t4=makeTLatex(t1Loc1, (t1Loc2-0.14), "#Delta#phi > #frac{2}{3}#pi");      t4->Draw();	
    if(doJetIDPlots){
      TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.18),"Jet ID Applied");      t5->Draw();	}
    
    float legx1=0.15, legx2=legx1+0.11;
    if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
    float legy1=0.74, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->SetBorderSize(0.);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();

    //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone("MCJetHistClone4Ratio");
    //I need to redefine the axis titles here
    std::string h_XAx_ratTitle=var_xAx_Titles[j], h_YAx_ratTitle="Data/MC"    ;
    std::string h_ratTitle   ="JetQA, Data/MC";//, "+var[j];
    if(doJetIDPlots)h_ratTitle+=", w/ JetIDCut";      
    
    //theRatio->SetMarkerStyle( kMultiply );
    //theRatio->SetMarkerSize(0.99 );
    //theRatio->SetMarkerColor( theMCColor );
    theRatio->SetLineColor( theMCLineColor );
    //theRatio->SetLineWidth( 2 );
    
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
    theRatio->SetYTitle(ratioTitle.c_str());
    theRatio->Divide(theDataJetQAHist);
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
    
  }
      
  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt"){ 
    jetpad1->cd();
    theMCJetQAHist->Draw("HIST E"); 
    theDataJetQAHist->Draw("E SAME"); 
    
    float t1Loc1=0.7, t1Loc2=0.72; 	    	  
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),"p_{T}^{ld} > 60 GeV");    t2->Draw();	
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),"p_{T}^{subld} > 30GeV");    t3->Draw();	
    TLatex *t4=makeTLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");	        t4->Draw();	
    if(doJetIDPlots){
      TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.22),"Jet ID Applied");      t5->Draw();	}
    
    float legx1=0.70, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->SetBorderSize(0.);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();
    
    //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone("MCJetHistClone4Ratio");
    //theRatio->SetMarkerStyle( kMultiply );
    //theRatio->SetMarkerSize(0.99 );
    //theRatio->SetMarkerColor( theMCColor );
    theRatio->SetLineColor( theMCLineColor );
    //theRatio->SetLineWidth( 2 );
    
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
    theRatio->SetYTitle( ratioTitle.c_str() );
    
    theRatio->Divide(theDataJetQAHist);
    theRatio->Draw("HIST E");
    
    TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,jetQAxmax[j],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,jetQAxmax[j],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,jetQAxmax[j],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
    
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");
  }
  
  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() );
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    outcanv->SetTitle((inHistName+" Canvas").c_str());
    outcanv->Write((inHistName+"_canv").c_str());    
    
    theDataJetQAHist->Write(("DATATH1_"+inHistName).c_str());
    theMCJetQAHist  ->Write(("MCTH1_"+inHistName).c_str());    
    theRatio        ->Write(("RATIOTH1_"+inHistName).c_str());
  }
  temp_canvJet->Close();  
  return;
}



void printJetSpectraRapHists( TFile* finData , TFile* finMC, bool doJetIDPlots,
			      std::string thePDFFileName , std::string fullJetType, 
			      long double theLumi  ) {
  bool funcDebug=true;
  if(funcDebug)std::cout<<"drawing dual-differential cross section plot"<<std::endl;	
  
  bool drawMC=(bool)finMC;
  if(funcDebug && drawMC)std::cout<<"drawing dual-differential cross section plot w/ MC"<<std::endl;	
  
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing Jet rapBins plots..."<<std::endl;
  TCanvas *temp_canvJetRapBins = new TCanvas("tempJetRapBins", "tempJetRapBins", 1600, 1200);
  temp_canvJetRapBins->SetLogy(1);    
  temp_canvJetRapBins->SetLogx(1);    
  temp_canvJetRapBins->cd();

  int latexsize=26;
  float t1Loc1=0.30, t1Loc2=0.85, deltaLoc=0.027;//for text later
  
  float TLegx2=0.88 ,TLegx1=TLegx2-.3;
  float TLegy2=0.88, TLegy1=TLegy2-.25;

  TLegend* theJetRapHistLeg=new TLegend(TLegx1, TLegy1, TLegx2, TLegy2, NULL,"brNDC");
  theJetRapHistLeg->SetFillStyle(0);
  theJetRapHistLeg->SetBorderSize(0);
  
  

  std::cout<<"looping over "<<nbins_abseta<<" rapidity bins"<<std::endl;      
  

  float totalIntXsec=0.0;
  const int bincounter_perm=nbins_abseta-1;  
  int bincounter=bincounter_perm;
  float power=std::pow(10., ((float)bincounter)); 
  
  
  float min=-1., max=-1.;
  //find min/max for making the range of the y axis appropriate
  if(funcDebug)std::cout<<"finding min/max of hists"<<std::endl;
  for(int rapbin=0; rapbin<nbins_abseta; ++rapbin){	

    if(rapbin==6)
      continue;  
    
    //only care about most central and most peripheral |y| bins
    std::string rapbin_str=absetabins_str[rapbin]+" < |y| <"+absetabins_str[rapbin+1];
    if(funcDebug)std::cout<<"rapbin="<<rapbin<<std::endl;
    if(funcDebug)std::cout<<"rapbin_str="<<rapbin_str<<std::endl;
    
    //open the hists + rebin
    std::string inHistName="hJetSpectraRap";	
    if(doJetIDPlots) inHistName+="_wJetID"; 
    inHistName+="_bin"+std::to_string(rapbin);
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;	    
    
    TH1F* theJetSpectraRapHist= (TH1F*) finData->Get( inHistName.c_str() );
    theJetSpectraRapHist=(TH1F*)theJetSpectraRapHist->
      TH1::Rebin(nbins_pt_debug, (inHistName+"_rebin").c_str(), ptbins_debug);
    //theJetSpectraRapHist->Scale(1./(1.e+6));//to nanobarns
    
    //rescaling
    divideBinWidth(theJetSpectraRapHist);    
    if(rapbin==7) theJetSpectraRapHist->Scale(1/3.); //3.2 - 4.7
    theJetSpectraRapHist->Scale( 1./theLumi);
    theJetSpectraRapHist->Scale( power );
    
    
    
    power/=10.;
    bincounter=bincounter-1;

    
    if(rapbin==0 || rapbin==(nbins_abseta-1))      {
      if(funcDebug)std::cout<<"this is one of the bins i'm looking for!"<<std::endl;}
    else{
      if(funcDebug)std::cout<<"not the most central nor most peripheral bin; continue"<<std::endl;
      continue;
    }

    if(rapbin==0){
      //theJetSpectraRapHist->Reset("M");//reset min/max
      max=theJetSpectraRapHist->GetMaximum()*1000.;
      min=1.e-7;
    }
    //else if(rapbin==(nbins_abseta-1))min=theJetSpectraRapHist->GetMaximum();///(1.e+8);//*(1.e-8);//because root thinks the minimum is 0 because empty bins
    
    if(funcDebug)std::cout<<"max="<<max<<std::endl;
    if(funcDebug)std::cout<<"min="<<min<<std::endl;
    
    
  }
  if(min<0. || max<0.)assert(false);
  //assert(false);
  bincounter=bincounter_perm;  
  power=std::pow(10., ((float)bincounter)); 
  for(int rapbin=0; rapbin<nbins_abseta; ++rapbin){	
    
    
    if(rapbin==6){  std::cout<<"rapbin==6, continuing on..."<<std::endl;
      continue;  }
    

    //open the hists + rebin
    std::string inHistName="hJetSpectraRap";	
    if(doJetIDPlots) inHistName+="_wJetID"; 
    inHistName+="_bin"+std::to_string(rapbin);
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;	
    
    
    TH1F* theJetSpectraRapHist= (TH1F*) finData->Get( inHistName.c_str() );
    theJetSpectraRapHist=(TH1F*)theJetSpectraRapHist->
      TH1::Rebin(nbins_pt_debug, (inHistName+"_rebin").c_str(), ptbins_debug);
    //theJetSpectraRapHist->Scale(1./(1.e+6));//to nanobarns
    //rescaling
    divideBinWidth(theJetSpectraRapHist);    
    if(rapbin==7) theJetSpectraRapHist->Scale(1/3.); //3.2 - 4.7
    theJetSpectraRapHist->Scale( 1./theLumi);
    std::cout<<"bincounter for power scaling=="<<bincounter<<std::endl;
    theJetSpectraRapHist->Scale( power );
    
    //integrated xsec calculation
    float integral=theJetSpectraRapHist->Integral(); //correct this, integral doesnt account for bin width
    std::cout<<"integral="<<integral/power<<std::endl;
    totalIntXsec+=(integral/power);
    std::cout<<"totalIntXsec="<<totalIntXsec<<std::endl;    
    
    theJetSpectraRapHist->SetMarkerStyle( kFullDotLarge ) ;
    theJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
    theJetSpectraRapHist->SetLineColor( kBlack );
    
    
    // MC Hist
    TH1F* theMCJetSpectraRapHist=NULL;
    if(drawMC){
      std::string inMCHistName="hJetSpectraRap";	
      if(doJetIDPlots) inMCHistName+="_wJetID"; 
      inMCHistName+="_bin"+std::to_string(rapbin);
      if(funcDebug)std::cout<<"inMCHistName="<<inMCHistName<<std::endl;	          
      theMCJetSpectraRapHist= (TH1F*)finMC->Get(inMCHistName.c_str()) ;
      theMCJetSpectraRapHist=(TH1F*)theMCJetSpectraRapHist->
	TH1::Rebin(nbins_pt_debug, (inMCHistName+"_rebin").c_str(), ptbins_debug);
      divideBinWidth(theMCJetSpectraRapHist);
      if(rapbin==7) theMCJetSpectraRapHist->Scale(1./3.); //3.2 - 4.7
      theMCJetSpectraRapHist->Scale( 1000. ); 
      theMCJetSpectraRapHist->Scale( power ); //i think the MC pthat weights are in units of nanobarns? multiply by 1000 
      //theMCJetSpectraRapHist->Scale( theJetSpectraRapHist->Integral()/theMCJetSpectraRapHist->Integral() );      
      theMCJetSpectraRapHist->SetLineColor( theRapOverlayMarkerColor[rapbin] );    }
    
    std::string rapHistLegDesc="("+absetabins_str[rapbin]+" < #||{y} < "+absetabins_str[rapbin+1]+") x 10^{  "+std::to_string(bincounter)+"}";	
    theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
    if(drawMC)theJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,"PY8 QCD MC","l");

    
    power/=10.;
    bincounter=bincounter-1;
    
    if(rapbin==0){
      
      std::string h_Title   ="Detector Level Jet p_{  T}";
      //if(doJetIDPlots) h_Title+=", w/ JetIDCut";
      theJetSpectraRapHist->SetTitle (    h_Title.c_str() );
      
      std::string h_XAx_Title="Jet p^{  RECO}_{T} (GeV)", h_YAx_Title="d^{2}#sigma/d#||{y}dp_{T} [nb/GeV]";//crossSectionAxTitle;
      theJetSpectraRapHist->SetXTitle( h_XAx_Title.c_str() );
      theJetSpectraRapHist->SetYTitle( h_YAx_Title.c_str() );
      theJetSpectraRapHist->SetAxisRange(min, max,"Y");
      //theJetSpectraRapHist->SetAxisRange(std::pow(10.,-10.),std::pow(10.,5.5),"Y");
      //theJetSpectraRapHist->SetAxisRange(std::pow(10.,-9.),std::pow(10.,6.5),"Y");
      //theJetSpectraRapHist->SetAxisRange(60.,840.,"X");
      
      theJetSpectraRapHist->GetXaxis()->SetMoreLogLabels(true);
      theJetSpectraRapHist->GetXaxis()->SetNoExponent(true);
      theJetSpectraRapHist->Draw("P E");     	  
      if(drawMC)theMCJetSpectraRapHist->Draw("HIST E SAME");     	  
    }
    else { 
      theJetSpectraRapHist->Draw("P E SAME");
      if(drawMC)theMCJetSpectraRapHist->Draw("HIST E SAME");
    }
    
  }//end rapbin hist loop

  theJetRapHistLeg->Draw(); 

  

  // STRINGS FOR HISTS

  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;

  std::cout<<"FINAL TOTAL INT JET XSEC="<<totalIntXsec<<std::endl;
  TLatex *pp1=makeTLatex(t1Loc1,t1Loc2,"2015 pp, #sqrt{s} = 5.02 TeV",latexsize);  
  TLatex *pp2=makeTLatex(t1Loc1,(t1Loc2-deltaLoc),"Jet80+LowerJets, L_{   int} = 27.4 pb^{ -1}",latexsize);  
  TLatex *t1=makeTLatex(t1Loc1,(t1Loc2-2*deltaLoc),(fullJetType+"Jets").c_str(),latexsize);   
  TLatex *t2=makeTLatex(t1Loc1,(t1Loc2-3*deltaLoc),(jetCutString).c_str(),latexsize);    
  pp1->Draw();  
  pp2->Draw();
  t1->Draw();
  t2->Draw();
  
  if(drawMC)    { 
    TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
    std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
    std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
    
    std::string genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
    std::cout<<"genJetCutString="<<genJetCutString<<std::endl;    
    
    TLatex *t3=makeTLatex(t1Loc1,(t1Loc2-4*deltaLoc),(genJetCutString).c_str());    
    t3->Draw(); 
  }
  if(doJetIDPlots)    { 
    TLatex *t4=makeTLatex(t1Loc1,(t1Loc2-4.5*deltaLoc),"Jet ID Applied",latexsize);      
    t4->Draw();	  
  }
  
  temp_canvJetRapBins->Print( thePDFFileName.c_str() );
  temp_canvJetRapBins->Close();  
  return;
}



void printMCJetSpectraRapHists(  TFile* finMC, bool doJetIDPlots,
				 std::string thePDFFileName , std::string fullJetType  ) {
  bool funcDebug=true;
  
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing Jet rapBins plots..."<<std::endl;
  TCanvas *temp_canvJetRapBins = new TCanvas("tempJetRapBins", "tempJetRapBins", 2000, 1500);
  temp_canvJetRapBins->SetLogy(1);    
  temp_canvJetRapBins->SetLogx(1);    
  
  TLegend* theJetRapHistLeg=new TLegend(0.72,0.60,0.85,0.86, NULL,"brNDC");
  
  temp_canvJetRapBins->cd();
  
  std::cout<<"drawing dual-differential cross section plot"<<std::endl;	
  std::cout<<"looping over "<<nbins_abseta<<" rapidity bins"<<std::endl;      
  
  float power=std::pow(10., 8.);
  int bincounter=8;
  float totalIntXsec=0.0;
  
  for(int rapbin=0; rapbin<8; ++rapbin){	
    
    
    if(rapbin!=6){
      power/=10;
      bincounter=bincounter-1;
      if(rapbin>5)continue;
    }
    
    else{  std::cout<<"rapbin==6, continuing on..."<<std::endl;
      continue;  }
    
    //open the hists + rebin
    std::string inHistName="hJetSpectraRap";	
    if(doJetIDPlots) inHistName+="_wJetID"; 
    inHistName+="_bin"+std::to_string(rapbin);
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;	
    
        
    TH1F* theJetSpectraRapHist= (TH1F*) finMC->Get( inHistName.c_str() );
    theJetSpectraRapHist=(TH1F*)theJetSpectraRapHist->
      TH1::Rebin(nbins_pt_debug, (inHistName+"_rebin").c_str(), ptbins_debug);
    divideBinWidth(theJetSpectraRapHist);
    if(rapbin==7) theJetSpectraRapHist->Scale(1/3.); //3.2 - 4.7
    theJetSpectraRapHist->Scale( power*1000 );
    
    float integral=theJetSpectraRapHist->Integral(); //correct this, integral doesnt account for bin width
    std::cout<<"integral="<<integral/power<<std::endl;
    totalIntXsec+=(integral/power);
    std::cout<<"totalIntXsec="<<totalIntXsec<<std::endl;
    
    
    theJetSpectraRapHist->SetMarkerStyle( kFullDotLarge ) ;
    theJetSpectraRapHist->SetMarkerColor( theRapOverlayMarkerColor[rapbin] );
    theJetSpectraRapHist->SetLineColor( kBlack );
    
    
    // MC Hist//hJetSpectraRap_1wJetID_bin5_genpt;1
    std::string inMCHistName="hJetSpectraRap";	
    if(doJetIDPlots) inMCHistName+="_1wJetID"; 
    else inMCHistName+="_0wJetID"; 
    inMCHistName+="_bin"+std::to_string(rapbin)+"_genpt";
    if(funcDebug)std::cout<<"inMCHistName="<<inMCHistName<<std::endl;	
    
    
    TH1F* theMCJetSpectraRapHist=NULL;
    theMCJetSpectraRapHist= (TH1F*)finMC->Get(inMCHistName.c_str()) ;
    theMCJetSpectraRapHist=(TH1F*)theMCJetSpectraRapHist->
      TH1::Rebin(nbins_pt_debug, (inMCHistName+"_rebin").c_str(), ptbins_debug);
    divideBinWidth(theMCJetSpectraRapHist);
    if(rapbin==7) theMCJetSpectraRapHist->Scale(1/3.); //3.2 - 4.7
    theMCJetSpectraRapHist->Scale( power * 1000 ); //i think the MC pthat weights are in units of nanobarns? multiply by 1000 
    //theMCJetSpectraRapHist->Scale( theJetSpectraRapHist->Integral()/theMCJetSpectraRapHist->Integral() );
    
    theMCJetSpectraRapHist->SetLineColor( theRapOverlayMarkerColor[rapbin] );
    
    
    std::stringstream stream1, stream2;
    stream1.precision(1); stream1 << std::fixed << absetabins[rapbin];
    stream2.precision(1); stream2 << std::fixed << absetabins[rapbin+1];
    
    if(rapbin==0){
      std::string rapHistLegDesc="(|y| < "+stream2.str()+") x 10^{"+std::to_string(bincounter)+"}";	
      theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
      theJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,"PY8 QCD MC GEN LEVEL","l");
    }
    else{
      std::string rapHistLegDesc="("+stream1.str()+" < |y| <"+stream2.str()+") x10^{"+std::to_string(bincounter)+"}";	
      theJetRapHistLeg->AddEntry(theJetSpectraRapHist,(rapHistLegDesc).c_str(),"lp");
      theJetRapHistLeg->AddEntry(theMCJetSpectraRapHist,"PY8 QCD MC GEN LEVEL","l");
    }
    
    if(rapbin==0){
      
      std::string h_Title   ="Jet Pt in |y| bins";
      if(doJetIDPlots) h_Title+=", w/ JetIDCut";
      theJetSpectraRapHist->SetTitle (    h_Title.c_str() );
      
      std::string h_XAx_Title="Jet p^{RECO}_{T} (GeV)", h_YAx_Title="d^{2}#sigma / d|y| dp_{T} (#mu b)";//crossSectionAxTitle;
      theJetSpectraRapHist->SetXTitle( h_XAx_Title.c_str() );
      theJetSpectraRapHist->SetYTitle( h_YAx_Title.c_str() );
      theJetSpectraRapHist->SetAxisRange(std::pow(10.,-9.),std::pow(10.,6.5),"Y");
      //theJetSpectraRapHist->SetAxisRange(60.,840.,"X");
      
      theJetSpectraRapHist->GetXaxis()->SetMoreLogLabels(true);
      theJetSpectraRapHist->GetXaxis()->SetNoExponent(true);
      theJetSpectraRapHist->Draw("E");     	  
      theMCJetSpectraRapHist->Draw("HIST E SAME");     	  
      theJetRapHistLeg->Draw(); 
    }
    else { 
      theJetSpectraRapHist->Draw("E SAME");
      theMCJetSpectraRapHist->Draw("HIST E SAME");
      }
    
  }//end rapbin hist loop
  
  std::cout<<"FINAL TOTAL INT JET XSEC="<<totalIntXsec<<std::endl;

  // STRINGS FOR HISTS
  TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
  std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
  
  std::string genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
  std::cout<<"genJetCutString="<<genJetCutString<<std::endl;
  
  TH1F* jtptCut_h= (TH1F*)finMC->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;

  float t1Loc1=0.62, t1Loc2=0.83, deltaLoc=.03; 	    
  TLatex *pp1=makeTLatex(t1Loc1,(t1Loc2),"PYTHIA8 ONLY, #sqrt{s} = 5.02 TeV"); pp1->Draw();
  TLatex *t1=makeTLatex(t1Loc1,(t1Loc2-deltaLoc),(fullJetType+"Jets").c_str());    t1->Draw();
  TLatex *t2=makeTLatex(t1Loc1,(t1Loc2-2*deltaLoc),(genJetCutString).c_str());    t2->Draw();
  TLatex *t3=makeTLatex(t1Loc1,(t1Loc2-3*deltaLoc),(jetCutString).c_str());    t3->Draw();
  if(doJetIDPlots){
    TLatex *t4=makeTLatex(t1Loc1,(t1Loc2-4*deltaLoc),"Jet ID Applied");      t4->Draw();	
  }
  
  temp_canvJetRapBins->Print( thePDFFileName.c_str() );
  temp_canvJetRapBins->Close();  
  return;
}












//void printJetQAHist( TFile* finData , TFile* finMC, int j, bool doJetIDPlots, int etabin,
void printTupelJetQAHist( TFile* finData , int j, bool doJetIDPlots, int etabin,
		     std::string inHistName , std::string thePDFFileName , std::string fullJetType, 
		     long double theLumi  ,TFile* fout=NULL) {

  
  std::string jetIDInt;
  if(doJetIDPlots)jetIDInt="1";
  else jetIDInt="0";
  
  
  //  if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
  if(!finData ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  
  
  // STRINGS FOR HISTS
  //  TH1F* genjtptCut_h= (TH1F*)finMC->Get( "hGenJetPtCut" );
  //  std::string genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  //  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;
  //  
  //  std::string genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
  //  std::cout<<"genJetCutString="<<genJetCutString<<std::endl;
  
  TH1F* jtptCut_h= (TH1F*)finData->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
//  TH1F* jtptQACut_h= (TH1F*)finData->Get( "hJetQAPtCut" );
//  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
//  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
//  TH1F* jtetaLoCut_h= (TH1F*)finData->Get( "hJetEtaCutLo" );
//  std::stringstream etaLo; etaLo.precision(1);
//  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
//  etaLo << std::fixed << jtetaLoCut_F;
//  std::string jtetaLoCut_str = etaLo.str();
//  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
  
  //  TH1F* jtetaHiCut_h= (TH1F*)finData->Get( "hJetEtaCutHi" );
  //  std::stringstream etaHi; etaHi.precision(1);
  //  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  //  etaHi << std::fixed << jtetaHiCut_F;
  //  std::string jtetaHiCut_str = etaHi.str();
//  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString;
  if(etabin==0) jetEtaCutString = "0.0 < #||{y} < 0.5";
  if(etabin==1) jetEtaCutString = "0.5 < #||{y} < 1.0";
  if(etabin==2) jetEtaCutString = "1.0 < #||{y} < 1.5";
  if(etabin==3) jetEtaCutString = "1.5 < #||{y} < 2.0";
  if(etabin==4) jetEtaCutString = "2.0 < #||{y} < 2.5";
  if(etabin==5) jetEtaCutString = "2.5 < #||{y} < 3.0";
//  if(jtetaLoCut_str=="0.0") jetEtaCutString="|y| < "+jtetaHiCut_str;
//  else jetEtaCutString=jtetaLoCut_str+" < |y| < "+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  //float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  float etaBinWidth=1.;//2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  
  
  
  
  // OPEN FILES
  std::string inDataHist=inHistName;
  std::cout<<" opening input Data "<<inDataHist<<std::endl<<std::endl;      
  
  TH1F* theDataJetQAHist=(TH1F*) ( (TH1*)finData->Get(inDataHist.c_str()) );
  theDataJetQAHist=(TH1F*) (theDataJetQAHist->Clone((inDataHist+"_"+std::to_string(j)+"clone").c_str()) );
  if(!theDataJetQAHist){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  
  theDataJetQAHist->Scale( 1./theLumi);
  theDataJetQAHist->Scale( 1./etaBinWidth);//eta bin width scaling
  
  
//  std::string inMCHist=inHistName;
//  float MCscale4Units=1./etaBinWidth;
//  if(tupelvar[j]=="jtpt_forRes" || tupelvar[j]=="jtpt_L2Res" || tupelvar[j]=="jtpt_L3Res"){
//    inMCHist="hJetQA_"+std::to_string((int)doJetIDPlots)+"wJetID_jtpt";
//    //MCscale4Units=1.;
//  }
  
//  std::cout<<" opening input MC "<<inMCHist<<std::endl<<std::endl;      
//  TH1F* theMCJetQAHist= (TH1F*) ( (TH1*)finMC->Get(inMCHist.c_str()) );  
//  theMCJetQAHist=(TH1F*) (theMCJetQAHist->Clone((inMCHist+"_"+std::to_string(j)+"clone").c_str()) );
//  if(!theMCJetQAHist){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
//    return;}
//  theMCJetQAHist->Scale(MCscale4Units);
  
  
  
  // MORE SCALING
  if(tupelvar[j]=="Pt" || tupelvar[j]=="RawPt" ||
     tupelvar[j]=="E" ||     tupelvar[j]=="RawE"){
    
    theDataJetQAHist=(TH1F*)theDataJetQAHist->TH1::Rebin(nbins_pt_debug, ("ppData_analysisBins_"+tupelvar[j]).c_str(), ptbins_debug);
    divideBinWidth(theDataJetQAHist);      
    
    //theMCJetQAHist=(TH1F*)theMCJetQAHist->TH1::Rebin(nbins_pt_debug, ("ppMC_analysisBins_"+tupelvar[j]).c_str(), ptbins_debug);
    //divideBinWidth(theMCJetQAHist); 
      
  }
  
  else {
    
    theDataJetQAHist->TH1::Rebin(tupelvar_xAx_reBin[j]);
    theDataJetQAHist->Scale( 1./theDataJetQAHist->GetBinWidth(1) );          
    theDataJetQAHist->SetAxisRange(tupeljetQAxmin[j],tupeljetQAxmax[j],"X"); 
    
    //theMCJetQAHist->TH1::Rebin(tupelvar_xAx_reBin[j]);
    //theMCJetQAHist->Scale( 1./theMCJetQAHist->GetBinWidth(1) );          
    //theMCJetQAHist->Scale( theDataJetQAHist->Integral()/theMCJetQAHist->Integral() );	        
    //theMCJetQAHist->SetAxisRange(tupeljetQAxmin[j],tupeljetQAxmax[j],"X"); 
    
  }
  
  if(funcDebug)theDataJetQAHist->Print("base");       std::cout<<std::endl;  
  //if(funcDebug)theMCJetQAHist->Print("base");       std::cout<<std::endl;
  


  
  // STYLE  
  dataHistStyle(theDataJetQAHist);
  //MCHistStyle(theMCJetQAHist);
  
  
  
  // TITLES
  std::string h_XAx_Title= tupelvar_xAx_Titles[j];
  //std::string h_YAx_Title= dcrossSectionAxTitle;
  //std::string h_YAx_Title= AUAxTitle;//for APS DNP
  
  std::string h_YAx_Title;
  if(tupelvar[j]=="Pt" || tupelvar[j]=="RawPt")h_YAx_Title= h_YAx_Title=ddcrossSectionAxTitle;//for APS DNP
  else h_YAx_Title="A.U.";
  
  std::string h_Title = PDStatsString_2; //for APS DNP
  
  theDataJetQAHist->SetTitle (    h_Title.c_str() );
  theDataJetQAHist->SetXTitle( h_XAx_Title.c_str() );
  theDataJetQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  //theMCJetQAHist->SetTitle (    h_Title.c_str() );
  //theMCJetQAHist->SetXTitle( h_XAx_Title.c_str() );
  //theMCJetQAHist->SetYTitle( h_YAx_Title.c_str() );
  
  
  
  //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone(  ("MCJetHistClone4Ratio_"+tupelvar[j]).c_str());  
  TH1F* theRatio=(TH1F*)theDataJetQAHist->Clone(  ("DataJetHistClone4Ratio_"+tupelvar[j]).c_str());  //place holder
  
  
  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvJet->cd();
  
    
  TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad1->SetGridy(0);
  jetpad1->SetGridx(0);
  jetpad1->SetLogy(1);


  TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  jetpad2->SetGridx(1);
  jetpad2->SetLogy(0);
  
  
  if(tupelvar[j]=="Pt" || tupelvar[j]=="RawPt" ||
     tupelvar[j]=="E" ||     tupelvar[j]=="RawE" || tupelvar[j]=="HfHadE" || tupelvar[j] == "HfEmE" ){
    
    
    jetpad1->SetLogx(1);  
    jetpad2->SetLogx(1);  
  }
  else {
    jetpad1->SetLogx(0);
    jetpad2->SetLogx(0);
  }
  
  jetpad1->SetBottomMargin(0);
  
  jetpad2->SetTopMargin(0);
  jetpad2->SetBottomMargin(0.3);
  
  
  jetpad1->Draw();
  jetpad2->Draw();  
  
  
  
  // variable specific draw orders, canv, etc.
  if(     tupelvar[j]=="Pt" || tupelvar[j]=="RawPt" ||
	  tupelvar[j]=="E" ||     tupelvar[j]=="RawE" || tupelvar[j]=="HfHadE" || tupelvar[j] == "HfEmE" ){
    
    jetpad1->cd();
    //if(tupelvar[j]=="jtpt"||tupelvar[j]=="rawpt"){
    //theMCJetQAHist->Draw("HIST E"); 
    //theDataJetQAHist->Draw("E SAME"); 
    theDataJetQAHist->Draw("HIST ][ E"); 
    //theMCJetQAHist->Draw("HIST E SAME ]["); 
    
    
    
    //  else if(tupelvar[j]=="jtpt_forRes"||tupelvar[j]=="jtpt_L2Res"||tupelvar[j]=="jtpt_L3Res"){
    //    theDataJetQAHist->Draw("E"); 
    //    theMCJetQAHist->Draw("HIST E SAME"); 
    //  } 
    
    float t1Loc1=0.50, t1Loc2=0.84;
    TLatex* t1= makeTLatex(  t1Loc1  ,  t1Loc2      ,  fullJetType+"Jets" );    t1->Draw();
    TLatex *t2= makeTLatex( (t1Loc1) , (t1Loc2-.05) , (jetCutString)      );    t2->Draw();	
    //    TLatex *t3= makeTLatex( (t1Loc1) , (t1Loc2-.10) , (genJetCutString)      ); t3->Draw();	
    //TLatex *t4= makeTLatex( (t1Loc1) , (t1Loc2-.15) , (jetEtaCutString)   );    t4->Draw();	    
    //if(doJetIDPlots){
    //  TLatex *t5= makeTLatex( (t1Loc1), (t1Loc2-.20), "Jet ID Applied" );      t5->Draw();	}
    TLatex *t4= makeTLatex( (t1Loc1) , (t1Loc2-.10) , (jetEtaCutString)   );    t4->Draw();	    
    if(doJetIDPlots){
      TLatex *t5= makeTLatex( (t1Loc1), (t1Loc2-.15), "Jet ID Applied" );      t5->Draw();	}
    
    float legx1=0.70, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->SetBorderSize(0.);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    //theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw(); 
    
    
    jetpad2->cd();
    
    
    TH1F* theDenom=(TH1F*)theDataJetQAHist->Clone(("DataHistClone4Ratio_"+tupelvar[j]).c_str());
    
    theRatio->SetLineColor( theMCLineColor );
    
    theRatio->SetAxisRange(0.5,1.5,"Y");
    
    theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);
    
    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(18);
    
    theRatio->SetTitle("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle(ratioTitle.c_str() );
    
    theRatio->Divide(theDenom);
    
    theRatio->GetXaxis()->SetMoreLogLabels(true);
    theRatio->GetXaxis()->SetNoExponent(true);
    
    theRatio->Draw("HIST ][ E");
    
    TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,ptbins_debug[nbins_pt_debug],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,ptbins_debug[nbins_pt_debug],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,ptbins_debug[nbins_pt_debug],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
    
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");



    
    
  }
  
  //  else if(tupelvar[j]=="Eta" ||
  //	  tupelvar[j]=="Phi"){
  else {
    
    jetpad1->cd();
    //if(tupelvar[j]=="jtphi"){
    //theMCJetQAHist->SetAxisRange(  .02,.8,"Y");
    //theDataJetQAHist->SetAxisRange(.02,.8,"Y");
    //};
    
    //theMCJetQAHist->Draw("HIST E"); 
    //theDataJetQAHist->Draw("E same");  
    theDataJetQAHist->Draw("HIST E ][ ");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    //float t1Loc1=0.50, t1Loc2=0.72; 	    
    
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
    //TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(genJetCutString).c_str());    t2->Draw();
    //TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();
    //TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();		
    //if(doJetIDPlots){
    //  TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.20),"Jet ID Applied");      t5->Draw();	
    //}
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.05),(jetCutString).c_str());    t3->Draw();
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.10),(jetEtaCutString).c_str());    t4->Draw();		
    if(doJetIDPlots){
      TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.15),"Jet ID Applied");      t5->Draw();	
    }
    
    float legx1=0.70, legx2=legx1+0.09;
    float legy1=0.76, legy2=legy1+0.09;

    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->SetBorderSize(0.);
    theJetQALeg->AddEntry(theDataJetQAHist,"Data","lp");
    //theJetQALeg->AddEntry(theMCJetQAHist,  "P8 MC"  ,"lp");
    theJetQALeg->Draw();
    
    //Ratio Plots
    
    jetpad2->cd();

    //TH1F* theRatio=(TH1F*)theMCJetQAHist->Clone("mcJetHistClone4Ratio");
    //    theRatio->SetMarkerStyle( kMultiply );
    //    theRatio->SetMarkerSize(0.99 );
    //    theRatio->SetMarkerColor( theMCColor );
    theRatio->SetLineColor( theMCLineColor );
    //theRatio->SetLineWidth( 2 );

    theRatio->SetAxisRange(0.,2.,"Y");
    theRatio->SetAxisRange(tupeljetQAxmin[j],tupeljetQAxmax[j],"X");	
    
    theRatio->GetYaxis()->SetTitleSize(15);
    theRatio->GetYaxis()->SetTitleFont(43);
    theRatio->GetYaxis()->SetTitleOffset(2);
    theRatio->GetYaxis()->SetLabelFont(43); 
    theRatio->GetYaxis()->SetLabelSize(13);

    theRatio->GetXaxis()->SetTitleSize(20);
    theRatio->GetXaxis()->SetTitleFont(43);
    theRatio->GetXaxis()->SetTitleOffset(4.);
    theRatio->GetXaxis()->SetLabelFont(43); 
    theRatio->GetXaxis()->SetLabelSize(18);
    
    theRatio->SetTitle("");
    theRatio->SetXTitle( h_XAx_Title.c_str() );
    theRatio->SetYTitle(ratioTitle.c_str());
    theRatio->Divide(theDataJetQAHist);
    theRatio->Draw("HIST E");
	  
    TLine* lineAtOne          = new TLine(tupeljetQAxmin[j],1.0,tupeljetQAxmax[j],1.0); 
    lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
    TLine* lineAtOneHalf      = new TLine(tupeljetQAxmin[j],0.5,tupeljetQAxmax[j],0.5); 
    lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
    TLine* lineAtOneEtOneHalf = new TLine(tupeljetQAxmin[j],1.5,tupeljetQAxmax[j],1.5); 
    lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
    
    lineAtOne->Draw("same");
    lineAtOneHalf->Draw("same");
    lineAtOneEtOneHalf->Draw("same");
    
    
  }
  
  
  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() );
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    outcanv->SetTitle((inHistName+" Canvas").c_str());
    outcanv->Write((inHistName+"_canv").c_str());    
    
    theDataJetQAHist->Write(("DATATH1_"+inHistName).c_str());
    //theMCJetQAHist  ->Write(("MCTH1_"+inHistName).c_str());    
    theRatio        ->Write(("RATIOTH1_"+inHistName).c_str());
  }
  temp_canvJet->Close();  
  return;
}
