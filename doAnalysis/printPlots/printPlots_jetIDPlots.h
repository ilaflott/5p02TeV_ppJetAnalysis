





void printJetIDHist( TFile* fin_jetID , TFile* fin_nojetID, 
		     int j, bool isData,
		     std::string thePDFFileName , std::string fullJetType, 
		     long double theLumi,  
		     TFile* fout=NULL) {
  
  
  if(!fin_jetID || !fin_nojetID ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
    
  // PUT JET CUTS SECTION IN HERE
  TH1F* genjtptCut_h=NULL; 
  std::string genjtptCut_str;
  if(!isData){
    genjtptCut_h=(TH1F*)fin_jetID->Get( "hGenJePttCut" );
    genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
    std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;}
  
  TH1F* jtptCut_h= (TH1F*)fin_jetID->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString;
  if(!isData)jetCutString="p_{T}^{gen}>"+genjtptCut_str+" GeV, "+"p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  else jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin_jetID->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin_jetID->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str(); 
 std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
 
 TH1F* jtetaHiCut_h= (TH1F*)fin_jetID->Get( "hJetEtaCutHi" );
  std::stringstream etaHi; etaHi.precision(1);
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  etaHi << std::fixed << jtetaHiCut_F;
  std::string jtetaHiCut_str = etaHi.str();
  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString;
  if(isData)sampleDescString=sqrts2k15ppString;
  else      sampleDescString=Py8TuneString;
  
  std::string inHistName="hJetQA_0wJetID_"+var[j];
  if(funcDebug)std::cout<<" opening input hist "<<inHistName<<std::endl<<std::endl;      
  
  TH1F* theNonJetIDHist=(TH1F*) ( (TH1*)fin_nojetID->Get(inHistName.c_str()) );
  if(!theNonJetIDHist){ std::cout<<"input hist not found! skipping hist"<<std::endl;
    return;}
  
  std::string inJetIDHistName="hJetQA_1wJetID_"+var[j];
  if(funcDebug)std::cout<<" opening input jetID hist "<<inJetIDHistName<<std::endl<<std::endl;      
  TH1F* theJetIDHist= (TH1F*) ( (TH1*)fin_jetID->Get(inJetIDHistName.c_str()) );  
  if(!theJetIDHist){ std::cout<<"input jetID hist not found! skipping hist"<<std::endl;
    return;}
  
  // MORE SCALING
  if(var[j]=="jtpt" || var[j]=="rawpt" ||
     var[j]=="leadJetPt" || var[j]=="subleadJetPt"      ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
    
    theNonJetIDHist=(TH1F*)theNonJetIDHist->TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_jtPt", ptbins_debug);
    divideBinWidth(theNonJetIDHist);      
    if(isData)theNonJetIDHist->Scale(1./theLumi);
    theNonJetIDHist->Scale(1./etaBinWidth);
    
    theJetIDHist=(TH1F*)theJetIDHist->TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_jtPt", ptbins_debug);
    divideBinWidth(theJetIDHist);      
    if(isData)theJetIDHist->Scale(1./theLumi);
    theJetIDHist->Scale(1./etaBinWidth);
  }
  
  else {
    
    theNonJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
    theNonJetIDHist->Scale( 1./theNonJetIDHist->GetBinWidth(1) );          
    if(isData)theNonJetIDHist->Scale(1./theLumi);
    theNonJetIDHist->Scale(1./etaBinWidth);
    
    theJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
    theJetIDHist->Scale( 1./theJetIDHist->GetBinWidth(1) );         
    if(isData)theJetIDHist->Scale(1./theLumi);
    theJetIDHist->Scale(1./etaBinWidth);
    //    theJetIDHist->Scale( theNonJetIDHist->Integral()/theJetIDHist->Integral() );	        
    
  }
  
  
  
  
  // TITLES
  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title=AUAxTitle;
  std::string h_Title="HIN Jet ID Efficiency";//, "+var[j];
  
  theNonJetIDHist->SetTitle (    h_Title.c_str() );
  theNonJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theNonJetIDHist->SetYTitle( h_YAx_Title.c_str() );
  
  theJetIDHist->SetTitle (    h_Title.c_str() );
  theJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theJetIDHist->SetYTitle( h_YAx_Title.c_str() );  
  
  
  // STYLE  
  TH1F* theRatio=(TH1F*)theJetIDHist->Clone("MCJetHistClone4Ratio");
  
  dataHistStyle(theNonJetIDHist);
  MCHistStyle(theJetIDHist);
  //ratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  jetIDratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  
  if(funcDebug)theNonJetIDHist->Print("base");       std::cout<<std::endl;  
  if(funcDebug)theJetIDHist->Print("base");       std::cout<<std::endl;
  
  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvJet->cd();
  
  TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad1->SetLogy(1); 
  jetpad1->SetGridx(1);
  jetpad1->SetGridy(1);
  jetpad1->SetBottomMargin(0);
  
    
  TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  jetpad2->SetTopMargin(0);
  jetpad2->SetBottomMargin(0.3);
  //jetpad2->SetLogy(0);
  jetpad2->SetLogy(1);
  
  
  if(var[j]=="jtpt"      || var[j]=="rawpt"          ||
     var[j]=="leadJetPt" || var[j]=="subleadJetPt"   ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
    
    jetpad1->SetLogx(1);  
    jetpad2->SetLogx(1);  
    
    theJetIDHist->GetXaxis()->SetNoExponent(true);
    theJetIDHist->GetXaxis()->SetMoreLogLabels(true);
    
    theNonJetIDHist->GetXaxis()->SetNoExponent(true);
    theNonJetIDHist->GetXaxis()->SetMoreLogLabels(true);

    theRatio->GetXaxis()->SetNoExponent(true);
    theRatio->GetXaxis()->SetMoreLogLabels(true);
  }
 
  jetpad1->Draw();
  jetpad2->Draw();    


  // variable specific draw orders, canv, etc.
  if(var[j]=="jtpt"||var[j]=="rawpt" ||
     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
    
    jetpad1->cd();    
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E SAME");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	    
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	    
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
    
    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    theJetQALeg->Draw(); 
    
    jetpad2->cd();
    
    theRatio->Divide(theNonJetIDHist); //jetid/nojetid
    fracSubtracted(theRatio); //1-jetid/nojetid


    theRatio->Draw("HIST E");
    
    

    
  }
  
  else if(var[j]=="jteta" ||
	  var[j]=="jtphi"	||
	  ( j>=jetConstits_varStart && j<dijet_varStart )){ 
    
    
    jetpad1->cd();
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E same");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
    
    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    
    theJetQALeg->Draw();
    
    jetpad2->cd();
    
    
    theRatio->Divide(theNonJetIDHist);
    fracSubtracted(theRatio); //1-jetid/nojetid
    theRatio->Draw("HIST E");
    
    
  }
  
  else if(var[j]=="Aj"||
	  var[j]=="xj"||
	  var[j]=="dphi"){ 
    
    jetpad1->cd();
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E SAME"); 
    
    float t1Loc1=0.15, t1Loc2=0.7; 	    	  
    if(var[j]=="Aj") t1Loc1=.7;    
    
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
    TLatex *t2=makeTLatex(t1Loc1,t1Loc2-.04,(fullJetType+"Jets").c_str());    t2->Draw();    
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),"p_{T}^{ld} > 60 GeV");    t3->Draw();		  
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.12),"p_{T}^{subld} > 30GeV");    t4->Draw();	    
    TLatex *t5=makeTLatex(t1Loc1, (t1Loc2-0.16), "#Delta#phi > #frac{2}{3}#pi");	        t5->Draw();	

    float legx1=0.15, legx2=legx1+0.11;
    if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
    float legy1=0.74, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();

    
    theRatio->Divide(theNonJetIDHist);
    fracSubtracted(theRatio); //1-jetid/nojetid
    theRatio->Draw("E");
    
    
  }
      
  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt"){ 
    
    jetpad1->cd();    
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E SAME"); 
    
    float t1Loc1=0.7, t1Loc2=0.72; 	    	  

    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();
    TLatex *t2=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t2->Draw();
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.05),"p_{T}^{ld} > 60 GeV");    t3->Draw();	
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.10),"p_{T}^{subld} > 30GeV");    t4->Draw();	
    TLatex *t5=makeTLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");    t5->Draw();	

    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    theJetQALeg->Draw(); 
    
    //Ratio Plots
    
    jetpad2->cd();
    
    theRatio->Divide(theNonJetIDHist);
    fracSubtracted(theRatio); //1-jetid/nojetid
    theRatio->Draw("HIST E");
    
  }
  
  
  //TLines for ratio plots
  float lowLim, highLim;
  lowLim=theRatio->GetBinLowEdge(1);
  highLim=theRatio->GetBinLowEdge(theRatio->GetNbinsX())+theRatio->GetBinWidth(theRatio->GetNbinsX());
  
  TLine* lineAtOne          = new TLine(lowLim,1.0,highLim,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  TLine* lineAtOneHalf      = new TLine(lowLim,0.5,highLim,0.5); 
  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  TLine* lineAtOneEtOneHalf = new TLine(lowLim,1.5,highLim,1.5); 
  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");

  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() ); 
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    
    std::string outcanvtitle="JetID Eff. "+var[j]+" Canvas";
    if(isData)outcanvtitle="Data "+outcanvtitle;
    else outcanvtitle="MC "+outcanvtitle;    
    outcanv->SetTitle(outcanvtitle.c_str());
    
    std::string outcanvname="JetIDEff_"+var[j];
    if(isData)outcanvname="Data_"+outcanvname;
    else outcanvname="MC_"+outcanvname;
    outcanv->Write(outcanvname.c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}







void printTupelJetIDHist( TFile* fin_jetID , int etabin, bool print_incjetana_tupelequiv,
			  int j, bool isData,
			  std::string thePDFFileName , std::string fullJetType, 
			  long double theLumi,  
			  TFile* fout=NULL) {
  
  
  if(!fin_jetID ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;

  if(!isData)assert(false);

  //// PUT JET CUTS SECTION IN HERE
  //TH1F* genjtptCut_h=NULL; 
  //std::string genjtptCut_str;
  //if(!isData){
  //  genjtptCut_h=(TH1F*)fin_jetID->Get( "hGenJetPtCut" );
  //  genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  //  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;}
  
  TH1F* jtptCut_h= (TH1F*)fin_jetID->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString;
  //if(!isData)jetCutString="p_{T}^{gen}>"+genjtptCut_str+" GeV, "+"p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  //else jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtptQACut_h= (TH1F*)fin_jetID->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  
//  TH1F* jtetaLoCut_h= (TH1F*)fin_jetID->Get( "hJetEtaCutLo" );
//  std::stringstream etaLo; etaLo.precision(1);
//  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
//  etaLo << std::fixed << jtetaLoCut_F;
//  std::string jtetaLoCut_str = etaLo.str(); 
// std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
 
 // TH1F* jtetaHiCut_h= (TH1F*)fin_jetID->Get( "hJetEtaCutHi" );
 //  std::stringstream etaHi; etaHi.precision(1);
 //  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
 //  etaHi << std::fixed << jtetaHiCut_F;
 //  std::string jtetaHiCut_str = etaHi.str();
 //  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
 //  
 //  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
 //  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  std::string jetEtaCutString;
  if(etabin==0) jetEtaCutString = "0.0 < #||{y} < 0.5";
  if(etabin==1) jetEtaCutString = "0.5 < #||{y} < 1.0";
  if(etabin==2) jetEtaCutString = "1.0 < #||{y} < 1.5";
  if(etabin==3) jetEtaCutString = "1.5 < #||{y} < 2.0";
  if(etabin==4) jetEtaCutString = "2.0 < #||{y} < 2.5";
  if(etabin==5) jetEtaCutString = "2.5 < #||{y} < 3.0";  
  //float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  float etaBinWidth=1.;
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString;
  if(isData)sampleDescString=sqrts2k15ppString;
  else      sampleDescString=Py8TuneString;
  
  
  std::string inHistName="hTupelJetQA_0jetID_etabin"+std::to_string(etabin)+"_"+tupelvar[j];          
  if(print_incjetana_tupelequiv){
    if(tupelvar_incjetanaequiv[j]=="NULL")return;	  
    inHistName="hJetQA_0wJetID_"+tupelvar_incjetanaequiv[j]+"_etabin"+std::to_string(etabin);          
  }
  
  if(funcDebug)std::cout<<" opening input hist "<<inHistName<<std::endl<<std::endl;      
  
  TH1F* theNonJetIDHist=(TH1F*) ( (TH1*)fin_jetID->Get(inHistName.c_str()) );
  if(!theNonJetIDHist){ std::cout<<"input hist not found! skipping hist"<<std::endl;
    return;}
  
  //std::string inJetIDHistName="hJetQA_1wJetID_"+tupelvar[j];
  std::string inJetIDHistName="hTupelJetQA_1jetID_etabin"+std::to_string(etabin)+"_"+tupelvar[j];          
  if(print_incjetana_tupelequiv){
    if(tupelvar_incjetanaequiv[j]=="NULL")return;	  
    inJetIDHistName="hJetQA_1wJetID_"+tupelvar_incjetanaequiv[j]+"_etabin"+std::to_string(etabin);          
  }
  if(funcDebug)std::cout<<" opening input jetID hist "<<inJetIDHistName<<std::endl<<std::endl;      
  TH1F* theJetIDHist= (TH1F*) ( (TH1*)fin_jetID->Get(inJetIDHistName.c_str()) );  
  if(!theJetIDHist){ std::cout<<"input jetID hist not found! skipping hist"<<std::endl;
    return;}
  
  // MORE SCALING
  if(tupelvar[j]=="Pt" || tupelvar[j]=="RawPt" ||
     tupelvar[j]=="E" ||     tupelvar[j]=="RawE"){
    
    theNonJetIDHist=(TH1F*)theNonJetIDHist->TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_jtPt", ptbins_debug);
    divideBinWidth(theNonJetIDHist);      
    if(isData)theNonJetIDHist->Scale(1./theLumi);
    theNonJetIDHist->Scale(1./etaBinWidth);
    
    theJetIDHist=(TH1F*)theJetIDHist->TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_jtPt", ptbins_debug);
    divideBinWidth(theJetIDHist);      
    if(isData)theJetIDHist->Scale(1./theLumi);
    theJetIDHist->Scale(1./etaBinWidth);
  }
  
  else {
    
    theNonJetIDHist->TH1::Rebin(tupelvar_xAx_reBin[j]);
    theNonJetIDHist->Scale( 1./theNonJetIDHist->GetBinWidth(1) );          
    if(isData)theNonJetIDHist->Scale(1./theLumi);
    theNonJetIDHist->Scale(1./etaBinWidth);
    
    theJetIDHist->TH1::Rebin(tupelvar_xAx_reBin[j]);
    theJetIDHist->Scale( 1./theJetIDHist->GetBinWidth(1) );         
    if(isData)theJetIDHist->Scale(1./theLumi);
    theJetIDHist->Scale(1./etaBinWidth);
    //    theJetIDHist->Scale( theNonJetIDHist->Integral()/theJetIDHist->Integral() );	        
    
  }
  
  
  
  
  // TITLES
  std::string h_XAx_Title=tupelvar_xAx_Titles[j];//, h_YAx_Title=AUAxTitle;
  std::string h_YAx_Title;
  if(tupelvar[j]=="Pt" || tupelvar[j]=="RawPt")h_YAx_Title= h_YAx_Title=ddcrossSectionAxTitle;//for APS DNP
  else h_YAx_Title="A.U.";
  //  std::string h_Title="HIN Jet ID Efficiency";//, "+var[j];
  std::string h_Title="SMP Jet ID Efficiency";//, "+var[j];
  if(print_incjetana_tupelequiv)
    h_Title="HIN Jet ID Efficiency";
  
  theNonJetIDHist->SetTitle (    h_Title.c_str() );
  theNonJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theNonJetIDHist->SetYTitle( h_YAx_Title.c_str() );
  
  theJetIDHist->SetTitle (    h_Title.c_str() );
  theJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theJetIDHist->SetYTitle( h_YAx_Title.c_str() );  
  
  
  // STYLE  
  TH1F* theRatio=(TH1F*)theJetIDHist->Clone("MCJetHistClone4Ratio");
  
  dataHistStyle(theNonJetIDHist);
  MCHistStyle(theJetIDHist);
  //ratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  jetIDratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  
  if(funcDebug)theNonJetIDHist->Print("base");       std::cout<<std::endl;  
  if(funcDebug)theJetIDHist->Print("base");       std::cout<<std::endl;
  
  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
  temp_canvJet->cd();
  
  TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad1->SetLogy(1); 
  jetpad1->SetGridx(1);
  jetpad1->SetGridy(1);
  jetpad1->SetBottomMargin(0);
  
    
  TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  jetpad2->SetTopMargin(0);
  jetpad2->SetBottomMargin(0.3);
  //jetpad2->SetLogy(0);
  jetpad2->SetLogy(1);
  
  
  if(tupelvar[j]=="Pt" || tupelvar[j]=="RawPt" ||
     tupelvar[j]=="E" ||     tupelvar[j]=="RawE" || tupelvar[j]=="HfHadE" || tupelvar[j] == "HfEmE" ){

    jetpad1->SetLogx(1);  
    jetpad2->SetLogx(1);  
    
    theJetIDHist->GetXaxis()->SetNoExponent(true);
    theJetIDHist->GetXaxis()->SetMoreLogLabels(true);
    
    theNonJetIDHist->GetXaxis()->SetNoExponent(true);
    theNonJetIDHist->GetXaxis()->SetMoreLogLabels(true);
    
    theRatio->GetXaxis()->SetNoExponent(true);
    theRatio->GetXaxis()->SetMoreLogLabels(true);
  }
  
  jetpad1->Draw();
  jetpad2->Draw();    
  
  
  // variable specific draw orders, canv, etc.
  if(     tupelvar[j]=="Pt" || tupelvar[j]=="RawPt" ||
	  tupelvar[j]=="E" ||     tupelvar[j]=="RawE" || tupelvar[j]=="HfHadE" || tupelvar[j] == "HfEmE" ){
    
    
    jetpad1->cd();    
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E SAME");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	    
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	    
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
    
    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    theJetQALeg->Draw(); 
    
    jetpad2->cd();
    
    theRatio->Divide(theNonJetIDHist); //jetid/nojetid
    fracSubtracted(theRatio); //1-jetid/nojetid


    theRatio->Draw("HIST E");
    
    

    
  }
  
//  else if(var[j]=="jteta" ||
//	  var[j]=="jtphi"	||
//	  ( j>=jetConstits_varStart && j<dijet_varStart )){ 
  else { 
    
    
    jetpad1->cd();
    theJetIDHist->Draw("HIST E"); 
    theNonJetIDHist->Draw("E same");  
    
    float t1Loc1=0.54, t1Loc2=0.82;
    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();
    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	
    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	
    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
    
    float legx1=0.73, legx2=legx1+0.11;
    float legy1=0.76, legy2=legy1+0.09;
    
    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
    theJetQALeg->SetFillStyle(0);
    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
    
    theJetQALeg->Draw();
    
    jetpad2->cd();
    
    
    theRatio->Divide(theNonJetIDHist);
    fracSubtracted(theRatio); //1-jetid/nojetid
    theRatio->Draw("HIST E");
    
    
  }

  
  
  
  //TLines for ratio plots
  float lowLim, highLim;
  lowLim=theRatio->GetBinLowEdge(1);
  highLim=theRatio->GetBinLowEdge(theRatio->GetNbinsX())+theRatio->GetBinWidth(theRatio->GetNbinsX());
  
  TLine* lineAtOne          = new TLine(lowLim,1.0,highLim,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  TLine* lineAtOneHalf      = new TLine(lowLim,0.5,highLim,0.5); 
  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  TLine* lineAtOneEtOneHalf = new TLine(lowLim,1.5,highLim,1.5); 
  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  //lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");

  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() ); 
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    
    std::string outcanvtitle="JetID Eff. "+tupelvar[j]+" Canvas";
    if(isData)outcanvtitle="Data "+outcanvtitle;
    else outcanvtitle="MC "+outcanvtitle;    
    outcanv->SetTitle(outcanvtitle.c_str());
    
    std::string outcanvname="JetIDEff_"+tupelvar[j];
    if(isData)outcanvname="Data_"+outcanvname;
    else outcanvname="MC_"+outcanvname;
    outcanv->Write(outcanvname.c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}



