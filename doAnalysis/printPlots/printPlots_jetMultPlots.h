





void printJetMultQATH1_vtx( TFile* finData ,  TFile* finMC, std::string inHistName, 
			    //			    TFile* finMC   , std::string inHistName  ,
			    std::string h_XAx_Title , int rebin, std::string thePDFFileName , long double theLumi) {
  
  if(!finData || !finMC ){    std::cout<<"input file(s) not found, cannot run printJetMultQATH1_vtx"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing jetMultQA plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("jetMultCanv"+inHistName).c_str(), "temp jetMult Canv", 1000, 800); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  //if( ((std::size_t)inHistName.find("_N"))!=std::string::npos || 
  //    ((std::size_t)inHistName.find(""))!=std::string::npos 
  //    )
  //  temp_canvEvt->SetLogy(1);
  //if(findsubstr(inHistName, "_N") ||
  //findsubstr(inHistName, "vtx_r") )
  
  if( findsubstr(inHistName, "vtx_z") )
    temp_canvEvt->SetLogy(0);
  else 
    temp_canvEvt->SetLogy(1);  
  //else if( findsubstr(inHistName, "Nvtxtrks_vtx0") )
  //  temp_canvEvt->SetLogy(0);
  //else if( findsubstr(inHistName, "Nvtxtrks_vtxgt0") )
  //  temp_canvEvt->SetLogy(0);
  //else if( inHistName=="Nvtxtrks" )
  //  temp_canvEvt->SetLogy(0);
  
  
  //  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  //  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  //  
  //  evtpad1->SetBottomMargin(0); 
  //  evtpad1->Draw();
  //  evtpad1->cd();
  
  
  // open hist
  std::cout<<" getting hist(s): "<<inHistName<<std::endl;
  
  TH1F* theDataJetMultTH1= (TH1F*)finData->Get( inHistName.c_str() );
  if(!theDataJetMultTH1) 
    std::cout<<"input data hist not found!"<<std::endl;
  
  TH1F* theMCJetMultTH1= (TH1F*)finMC->Get( inHistName.c_str() );
  if(!theMCJetMultTH1) 
    std::cout<<"input MC hist not found!"<<std::endl;
  
  if(!theDataJetMultTH1 || !theMCJetMultTH1) { 
    std::cout<<"skipping hist!"<<std::endl<<std::endl;
    return;   }
  
  theDataJetMultTH1=(TH1F*)theDataJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theDataJetMultTH1->Scale( 1./theDataJetMultTH1->GetBinWidth(1) );
  //theDataJetMultTH1->Scale( 1./theLumi );
  dataHistStyle(theDataJetMultTH1);
  
  
  theMCJetMultTH1=(TH1F*)theMCJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theMCJetMultTH1->Scale( 1./theMCJetMultTH1->GetBinWidth(1) );
  //theMCJetMultTH1->Scale( 1000. );
  theMCJetMultTH1->Scale( theDataJetMultTH1->Integral()/theMCJetMultTH1->Integral() );	
  MCHistStyle(theMCJetMultTH1);
  
  
  
  
  
  


  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  std::string h_YAx_Title ="A.U.";  

  //std::string h_XAx_Title="";
  //if( inHistName.find("Vz") != std::string::npos ){
  //  h_XAx_Title ="v_{z}^{evt} (cm)";}
  //else if( inHistName.find("Vx") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{x}^{evt} (cm)";}
  //else if( inHistName.find("Vy") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{y}^{evt} (cm)";}
  //else if( inHistName.find("Vr") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{r}^{evt} (cm)";}
  
  //Double_t dataMax = theDataJetMultTH1->GetMaximum();//, MCMax=theMCJetMultTH1->GetMaximum();
  //std::cout<<"dataMax="<<dataMax<<std::endl;
  //std::cout<<"MCMax="<<    MCMax<<std::endl;
  //if(dataMax>MCMax){
  //  std::cout<<"dataMax > MC Max"<<std::endl;
  //  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  
  //  theDataJetMultTH1->Draw("E"); 
  //  theMCJetMultTH1->Draw("HIST E SAME"); 
  //  theDataJetMultTH1->Draw("E SAME"); //because i want data points written over MC's red line hist
  //}
  //else {
  //  theMCJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theMCJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theMCJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  theMCJetMultTH1->SetAxisRange(.001,10.*MCMax,"Y");
  //  theMCJetMultTH1->Draw("HIST E"); 
  //  theDataJetMultTH1->Draw("E SAME"); 
  //}
  
  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  
  theDataJetMultTH1->Draw("E"); 
  theMCJetMultTH1->Draw("HIST E SAME");

  // legend
  float legx1=0.80, legx2=legx1+0.10;
  float legy1=0.80, legy2=legy1+0.10;
  
  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  //theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->AddEntry(theDataJetMultTH1, " data" ,"lp");
  theEvtQALeg->AddEntry(theMCJetMultTH1,   " MC"       ,"lp");
  theEvtQALeg->Draw();
  
  //// text
  //float pp1x=0.55,pp1y=0.85;
  //TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str());  pp1->Draw();
  //TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str());  pp2->Draw();  
  //TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();  
  
  
  //ratios go here
  //temp_canvEvt->cd(); //change back to main canvas
  
  //evtpad2->Draw();
  //evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  //TH1F *theRatio=(TH1F*)theMCJetMultTH1->Clone("MCClone4Ratio");
  //TH1F *theDenom=(TH1F*)theDataJetMultTH1->Clone("DataClone4Ratio");
  //evtpad2->SetTopMargin(0);
  //evtpad2->SetBottomMargin(0.3);
  //
  //theRatio->GetYaxis()->SetTitleSize(15);
  //theRatio->GetYaxis()->SetTitleFont(43);
  //theRatio->GetYaxis()->SetTitleOffset(2);
  //theRatio->GetYaxis()->SetLabelFont(43); 
  //theRatio->GetYaxis()->SetLabelSize(13);
  //theRatio->GetXaxis()->SetTitleSize(20);
  //theRatio->GetXaxis()->SetTitleFont(43);
  //theRatio->GetXaxis()->SetTitleOffset(4.);
  //theRatio->GetXaxis()->SetLabelFont(43); 
  //theRatio->GetXaxis()->SetLabelSize(15);
  //
  //theRatio->SetTitle("");
  //theRatio->SetXTitle( h_XAx_Title.c_str() );
  //theRatio->SetYTitle( ratioTitle.c_str() );
  //
  //theRatio->Divide(theDenom);
  //theRatio->SetAxisRange(0.0,2.0, "Y");  
  //
  //theRatio->Draw("HIST E");
  //
  //float min = theRatio->GetBinLowEdge( 1 );
  //float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
  //  + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  //
  //TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  //lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  //TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  //TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  //
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  std::cout<<std::endl;
  return;
}







//void printJetMultQATH1_trk( TFile* finData , std::string inDataHistName, 
void printJetMultQATH1_trk( TFile* finData , TFile* finMC, std::string inHistName, 
			    std::string h_XAx_Title , int rebin, std::string thePDFFileName , long double theLumi) {
  
  if(!finData || !finMC ){    std::cout<<"input file(s) not found, cannot run printJetMultQATH1_vtx"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing jetMultQA plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("jetMultCanv"+inHistName).c_str(), "temp jetMult Canv", 1000, 800); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  //if( ((std::size_t)inHistName.find("_N"))!=std::string::npos || 
  //    ((std::size_t)inHistName.find(""))!=std::string::npos 
  //    )
  //  temp_canvEvt->SetLogy(1);
  //if(findsubstr(inHistName, "_N") ||
  //findsubstr(inHistName, "vtx_r") )
  temp_canvEvt->SetLogy(1);
  if(findsubstr(inHistName, "trkPtPercError"))
    temp_canvEvt->SetLogx(0);
  else if(findsubstr(inHistName, "trkPtError"))
    temp_canvEvt->SetLogx(0);
  else if(findsubstr(inHistName, "trkPt"))
    temp_canvEvt->SetLogx(0);
  //  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  //  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  //  
  //  evtpad1->SetBottomMargin(0); 
  //  evtpad1->Draw();
  //  evtpad1->cd();
  
  
  // open hist
  std::cout<<" getting hist(s): "<<inHistName<<std::endl;
  
  TH1F* theDataJetMultTH1= (TH1F*)finData->Get( inHistName.c_str() );
  if(!theDataJetMultTH1){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  theDataJetMultTH1=(TH1F*)theDataJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theDataJetMultTH1->Scale( 1./theDataJetMultTH1->GetBinWidth(1) );
  //theDataJetMultTH1->Scale( 1./theLumi );
  dataHistStyle(theDataJetMultTH1);
  
  
  TH1F* theMCJetMultTH1= (TH1F*)finMC->Get( inHistName.c_str() );
  if(!theMCJetMultTH1){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
    return;}
  theMCJetMultTH1=(TH1F*)theMCJetMultTH1->TH1::Rebin( (int)rebin, (inHistName+"_rebin").c_str());
  theMCJetMultTH1->Scale( 1./theMCJetMultTH1->GetBinWidth(1) );
  //theMCJetMultTH1->Scale(1000.);
  theMCJetMultTH1->Scale( theDataJetMultTH1->Integral()/theMCJetMultTH1->Integral() );	
  MCHistStyle(theMCJetMultTH1);
  
  
  
  
  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  std::string h_YAx_Title ="A.U.";  
  //std::string h_XAx_Title="";
  //if( inHistName.find("Vz") != std::string::npos ){
  //  h_XAx_Title ="v_{z}^{evt} (cm)";}
  //else if( inHistName.find("Vx") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{x}^{evt} (cm)";}
  //else if( inHistName.find("Vy") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{y}^{evt} (cm)";}
  //else if( inHistName.find("Vr") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{r}^{evt} (cm)";}
  
  //Double_t dataMax = theDataJetMultTH1->GetMaximum();//, MCMax=theMCJetMultTH1->GetMaximum();
  //std::cout<<"dataMax="<<dataMax<<std::endl;
  //std::cout<<"MCMax="<<    MCMax<<std::endl;
  //if(dataMax>MCMax){
  //  std::cout<<"dataMax > MC Max"<<std::endl;
  //  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  
  //  theDataJetMultTH1->Draw("E"); 
  //  theMCJetMultTH1->Draw("HIST E SAME"); 
  //  theDataJetMultTH1->Draw("E SAME"); //because i want data points written over MC's red line hist
  //}
  //else {
  //  theMCJetMultTH1->SetTitle (    h_Title.c_str() );
  //  theMCJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  //  theMCJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  //  theMCJetMultTH1->SetAxisRange(.001,10.*MCMax,"Y");
  //  theMCJetMultTH1->Draw("HIST E"); 
  //  theDataJetMultTH1->Draw("E SAME"); 
  //}
  
  if(findsubstr(inHistName, "trkPtPercError"))
    theDataJetMultTH1->SetAxisRange(0, 20., "X");
  
  theDataJetMultTH1->SetTitle (    h_Title.c_str() );
  theDataJetMultTH1->SetXTitle( h_XAx_Title.c_str() );
  theDataJetMultTH1->SetYTitle( h_YAx_Title.c_str() );
  
  theDataJetMultTH1->Draw("E"); 
  theMCJetMultTH1->Draw("HIST E SAME"); 
  
  // legend
  float legx1=0.80, legx2=legx1+0.10;
  float legy1=0.80, legy2=legy1+0.10;
  
  TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  theEvtQALeg->SetFillStyle(0);
  theEvtQALeg->AddEntry(theDataJetMultTH1, " data" ,"lp");
  theEvtQALeg->AddEntry(theMCJetMultTH1,   " MC"       ,"lp");
  theEvtQALeg->Draw();
  
  // text
  float pp1x=0.55,pp1y=0.85;
  TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str());  pp1->Draw();
  TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str());  pp2->Draw();  
  //TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();
  


  //ratios go here
  //temp_canvEvt->cd(); //change back to main canvas
  
  //evtpad2->Draw();
  //evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  //TH1F *theRatio=(TH1F*)theMCJetMultTH1->Clone("MCClone4Ratio");
  //TH1F *theDenom=(TH1F*)theDataJetMultTH1->Clone("DataClone4Ratio");
  //evtpad2->SetTopMargin(0);
  //evtpad2->SetBottomMargin(0.3);
  //
  //theRatio->GetYaxis()->SetTitleSize(15);
  //theRatio->GetYaxis()->SetTitleFont(43);
  //theRatio->GetYaxis()->SetTitleOffset(2);
  //theRatio->GetYaxis()->SetLabelFont(43); 
  //theRatio->GetYaxis()->SetLabelSize(13);
  //theRatio->GetXaxis()->SetTitleSize(20);
  //theRatio->GetXaxis()->SetTitleFont(43);
  //theRatio->GetXaxis()->SetTitleOffset(4.);
  //theRatio->GetXaxis()->SetLabelFont(43); 
  //theRatio->GetXaxis()->SetLabelSize(15);
  //
  //theRatio->SetTitle("");
  //theRatio->SetXTitle( h_XAx_Title.c_str() );
  //theRatio->SetYTitle( ratioTitle.c_str() );
  //
  //theRatio->Divide(theDenom);
  //theRatio->SetAxisRange(0.0,2.0, "Y");  
  //
  //theRatio->Draw("HIST E");
  //
  //float min = theRatio->GetBinLowEdge( 1 );
  //float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
  //  + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  //
  //TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  //lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  //TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  //TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  //
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  std::cout<<std::endl;
  return;
}












void printJetMultQATH2( TFile* finData , std::string inDataHistName, 
			//TFile* finMC   , std::string inMCHistName  ,
			std::string h_XAx_Title, int xrebin, 
			std::string h_YAx_Title , int yrebin,
			std::string thePDFFileName ) {
  
  //if(!finData || !finMC ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
  if(!finData ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(funcDebug)std::cout<<"creating temporary canvas for printing jetMultQA plots..."<<std::endl;
  TCanvas *temp_canvEvt = new TCanvas(("jetMultCanv"+inDataHistName).c_str(), "temp jetMult Canv", 1200, 800); //This is for the Jet QA Evt plots without logs
  temp_canvEvt->cd();
  
  //if( ((std::size_t)inDataHistName.find("_N"))!=std::string::npos || 
  //    ((std::size_t)inDataHistName.find(""))!=std::string::npos 
  //    )
  //  temp_canvEvt->SetLogy(1);
  //if(findsubstr(inDataHistName, "_N") ||
  //findsubstr(inDataHistName, "vtx_r") )
  //  temp_canvEvt->SetLogy(1);
  
  if(findsubstr(inDataHistName,"trkEtaVPhi"))
    temp_canvEvt->SetLogz(0);
  else
    temp_canvEvt->SetLogz(1);
  
  //  TPad *evtpad1 = new TPad("evtpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  //  TPad *evtpad2 = new TPad("evtpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  //  
  //  evtpad1->SetBottomMargin(0); 
  //  evtpad1->Draw();
  //  evtpad1->cd();
  
  
  // open hist
  std::cout<<" getting hist(s): "<<inDataHistName<<std::endl;
  
  TH2F* theDataJetMultTH2= (TH2F*)finData->Get( inDataHistName.c_str() );
  if(!theDataJetMultTH2){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
    return;}
  theDataJetMultTH2=(TH2F*)theDataJetMultTH2->TH2::Rebin2D( (int)xrebin, (int)yrebin, (inDataHistName+"_rebin").c_str());
  
  float scaleFactor = theDataJetMultTH2->GetXaxis()->GetBinWidth(1);
  scaleFactor *= theDataJetMultTH2->GetYaxis()->GetBinWidth(1);
  theDataJetMultTH2->Scale( 1./scaleFactor );
  //theDataJetMultTH2->Scale( 1./theLumi );
  //dataHistStyle(theDataJetMultTH2);
  
  
  //TH2F* theMCJetMultTH2= (TH2F*)finMC->Get( inMCHistName.c_str() );
  //if(!theMCJetMultTH2){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
  //  return;}
  //theMCJetMultTH2=(TH2F*)theMCJetMultTH2->TH1::Rebin( (int)rebin, (inMCHistName+"_rebin").c_str());
  //theMCJetMultTH2->Scale( 1./theMCJetMultTH2->GetBinWidth(1) );
  //theMCJetMultTH2->Scale( theDataJetMultTH2->Integral()/theMCJetMultTH2->Integral() );	
  //MCHistStyle(theMCJetMultTH2);
  
  


  // title+axes
  
  std::string h_Title     =PDStatsString_2;
  //std::string h_YAx_Title ="A.U.";  
  //std::string h_XAx_Title="";
  //if( inDataHistName.find("Vz") != std::string::npos ){
  //  h_XAx_Title ="v_{z}^{evt} (cm)";}
  //else if( inDataHistName.find("Vx") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{x}^{evt} (cm)";}
  //else if( inDataHistName.find("Vy") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{y}^{evt} (cm)";}
  //else if( inDataHistName.find("Vr") != std::string::npos ){
  //  evtpad1->SetLogy(1);
  //  h_XAx_Title ="v_{r}^{evt} (cm)";}
  
  //Double_t dataMax = theDataJetMultTH2->GetMaximum();//, MCMax=theMCJetMultTH2->GetMaximum();
  //std::cout<<"dataMax="<<dataMax<<std::endl;
  //std::cout<<"MCMax="<<    MCMax<<std::endl;
  //if(dataMax>MCMax){
  //  std::cout<<"dataMax > MC Max"<<std::endl;
  //  theDataJetMultTH2->SetTitle (    h_Title.c_str() );
  //  theDataJetMultTH2->SetXTitle( h_XAx_Title.c_str() );
  //  theDataJetMultTH2->SetYTitle( h_YAx_Title.c_str() );
  //  
  //  theDataJetMultTH2->Draw("E"); 
  //  theMCJetMultTH2->Draw("HIST E SAME"); 
  //  theDataJetMultTH2->Draw("E SAME"); //because i want data points written over MC's red line hist
  //}
  //else {
  //  theMCJetMultTH2->SetTitle (    h_Title.c_str() );
  //  theMCJetMultTH2->SetXTitle( h_XAx_Title.c_str() );
  //  theMCJetMultTH2->SetYTitle( h_YAx_Title.c_str() );
  //  theMCJetMultTH2->SetAxisRange(.001,10.*MCMax,"Y");
  //  theMCJetMultTH2->Draw("HIST E"); 
  //  theDataJetMultTH2->Draw("E SAME"); 
  //}
  
  theDataJetMultTH2->SetTitle (    h_Title.c_str() );
  theDataJetMultTH2->SetXTitle( h_XAx_Title.c_str() );
  theDataJetMultTH2->SetYTitle( h_YAx_Title.c_str() );
  
  theDataJetMultTH2->Draw("COLZ"); 
  
  //// legend
  //float legx1=0.80, legx2=legx1+0.10;
  //float legy1=0.80, legy2=legy1+0.10;
  //
  //TLegend* theEvtQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");
  //theEvtQALeg->SetFillStyle(0);
  //theEvtQALeg->AddEntry(theDataJetMultTH2, " data" ,"lp");
  //theEvtQALeg->AddEntry(theMCJetMultTH2,   " MC"       ,"lp");
  //theEvtQALeg->Draw();
  
  // text
  //float pp1x=0.55,pp1y=0.85;
  //TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str());  //pp1->Draw();
  //TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str()); // pp2->Draw();  
  //TLatex *pp4=makeTLatex(pp1x,(pp1y-0.08),Py8TuneString.c_str() );  pp4->Draw();
  
  
  
  //ratios go here
  //temp_canvEvt->cd(); //change back to main canvas
  
  //evtpad2->Draw();
  //evtpad2->cd();
  
  
  
  //hist for ratio; use clone
  
  //TH1F *theRatio=(TH1F*)theMCJetMultTH1->Clone("MCClone4Ratio");
  //TH1F *theDenom=(TH1F*)theDataJetMultTH1->Clone("DataClone4Ratio");
  //evtpad2->SetTopMargin(0);
  //evtpad2->SetBottomMargin(0.3);
  //
  //theRatio->GetYaxis()->SetTitleSize(15);
  //theRatio->GetYaxis()->SetTitleFont(43);
  //theRatio->GetYaxis()->SetTitleOffset(2);
  //theRatio->GetYaxis()->SetLabelFont(43); 
  //theRatio->GetYaxis()->SetLabelSize(13);
  //theRatio->GetXaxis()->SetTitleSize(20);
  //theRatio->GetXaxis()->SetTitleFont(43);
  //theRatio->GetXaxis()->SetTitleOffset(4.);
  //theRatio->GetXaxis()->SetLabelFont(43); 
  //theRatio->GetXaxis()->SetLabelSize(15);
  //
  //theRatio->SetTitle("");
  //theRatio->SetXTitle( h_XAx_Title.c_str() );
  //theRatio->SetYTitle( ratioTitle.c_str() );
  //
  //theRatio->Divide(theDenom);
  //theRatio->SetAxisRange(0.0,2.0, "Y");  
  //
  //theRatio->Draw("HIST E");
  //
  //float min = theRatio->GetBinLowEdge( 1 );
  //float max = theRatio->GetBinLowEdge( (int)theRatio->GetNbinsX() ) 
  //  + theRatio->GetBinWidth( (int)theRatio->GetNbinsX() ) ;
  //
  //TLine* lineAtOne          = new TLine(min,1.0,max,1.0); 
  //lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  //TLine* lineAtOneHalf      = new TLine(min,0.5,max,0.5); 
  //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  //TLine* lineAtOneEtOneHalf = new TLine(min,1.5,max,1.5); 
  //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  //
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");	
  
  // print to PDF file
  temp_canvEvt->Print( thePDFFileName.c_str() ); //ratio plots moved above.
  temp_canvEvt->Close();
  std::cout<<std::endl;
  return;
}
