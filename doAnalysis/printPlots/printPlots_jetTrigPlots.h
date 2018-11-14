

const int rebinfactor=1;//for simpbins aka not-analysis-bins


void printJetTrigHist_wRatio( TFile* fin , bool usedHLT100, bool analysisRebin,
			      std::string thePDFFileName , std::string fullJetType , 
			      std::string trigType, std::string radius , bool usedHLTPF, bool didJetID , TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=true;
  //bool analysisRebin=true;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  
  
  
  
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
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
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigSpec"+trigType).c_str(), ("blahTrigPt"+trigType).c_str(), 1200, 1000);
  temp_canvJetTrig->cd();

  

  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING JET SPECTRA BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  TPad *jetpad_excsp = new TPad("jetpad_excsp", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad_excsp->SetLogx(1);
  jetpad_excsp->SetLogy(1);
  jetpad_excsp->SetGridx(1);
  jetpad_excsp->SetGridy(1);
  
  jetpad_excsp->SetBottomMargin(0);
  jetpad_excsp->Draw();
  
  TPad *jetpad_excrat = new TPad("jetpad_excrat", "Ratio Pad", 0.0, 0.05, 1.0, 0.3);
  jetpad_excrat->SetLogx(1);
  jetpad_excrat->SetLogy(0);
  jetpad_excrat->SetGridx(1);
  jetpad_excrat->SetGridy(0);
  
  jetpad_excrat->SetTopMargin(0);
  jetpad_excrat->SetBottomMargin(0.3);
  
  jetpad_excrat->Draw();
  
  jetpad_excsp->cd();
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.66,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  for(int j=0; j<(N_trigs); j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    //open the hists + do scaling
    std::string inHistName;
    if(trigType=="excl")inHistName="hpp_exc";
    else if(trigType=="incl")inHistName="hpp_";
    
    if(j==0)inHistName="hpp_";
    inHistName+=HLTName[j];
    if(didJetID)inHistName+="_wJetID";
    inHistName+="_"+radius+etaWidth;
    
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base"); std::cout<<std::endl;
    theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
    theJetTrigQAHist->GetXaxis()->SetNoExponent(true);
    
    if(analysisRebin) {
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_spec_anabins").c_str(), ptbins_debug    );
      divideBinWidth(theJetTrigQAHist);
    }
    else {
      //theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(10, (inHistName+"_spec_rebin10").c_str() );
      //theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(1, (inHistName+"_spec_rebin1").c_str() );
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_spec_rebin1").c_str() );
      theJetTrigQAHist->Scale(1./((float)rebinfactor));
    }
    theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");
	  
    
    if(funcDebug)std::cout<<"setting legend entry titles"<<std::endl;
    if(usedHLTPF){
      if(funcDebug)std::cout<<"using HLTPF Trig titles"<<std::endl;
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    }
    else {
      if(funcDebug)std::cout<<"using HLTCalo Trig titles"<<std::endl;
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    }
    
    if( j==0 ){
      
      std::string h_Title;
      if(trigType=="excl")h_Title="Excl. Jet Spectra by Trigger";
      else if(trigType=="incl")h_Title="Incl. Jet Spectra by Trigger";      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      
      std::string h_YAx_Title= AUAxTitle;
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );      
      
      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw("E");
      
      
      float t1Loc1=0.54, t1Loc2=0.82;
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      t1->Draw();
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),(jetCutString).c_str());      t2->Draw();
      TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),(jetEtaCutString).c_str());      t3->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.12;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );      pp1->Draw();
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );    pp2->Draw();      
    }
    else{
      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw("E SAME");}
    
  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();
  
  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------ "<<std::endl;
  std::cout<<" --------------- DRAWING JET RATIOS BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------ "<<std::endl<<std::endl<<std::endl;
  
  
  jetpad_excrat->cd();
  
  
  TH1F* theDenominator=NULL;
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    
    if( j==0 ){//combo hlt, the denominator
      
      //open the hists + do scaling
      inHistName="hpp_"+HLTName[j];
      if(didJetID)inHistName+="_wJetID";      
      inHistName+="_"+radius+etaWidth;
      
      if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
      
      theDenominator = (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
      if(funcDebug)theDenominator->Print("base");
      
      if(analysisRebin) 
	theDenominator=(TH1F*)theDenominator->TH1::Rebin(nbins_pt_debug,(inHistName+"_rat_anabins").c_str(), ptbins_debug    );
      else {
	theDenominator=(TH1F*)theDenominator->TH1::Rebin(rebinfactor, (inHistName+"_rat_rebin10").c_str() );
	theDenominator->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");
      }
      
      trigRatioHistStyle(theDenominator,j);      
      
      
    }
    
    else if(j>0){//individial hlt paths, numerators
      
      if(trigType=="excl")      inHistName="hpp_exc"+HLTName[j];//+"_"+radius+etaWidth;
      else if(trigType=="incl") inHistName="hpp_"+HLTName[j];//+"_"+radius+etaWidth;
      if(didJetID)inHistName+="_wJetID";
      inHistName+="_"+radius+etaWidth;
      
      if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
      
      TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
      if(funcDebug)theJetTrigQAHist->Print("base");
      
      if(analysisRebin) 
	theJetTrigQAHist= (TH1F*) theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_rat_anabins").c_str(), ptbins_debug    );
      else {
	theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_rat_rebin10").c_str() );
	theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");}
      
      trigRatioHistStyle(theJetTrigQAHist, j);
      
      theJetTrigQAHist->Divide(theDenominator);
      
      theJetTrigQAHist->SetAxisRange(0.0,2.0,"Y");

      std::string h_XAx_Title="jet p_{T} (GeV)     ";
      theJetTrigQAHist->SetXTitle( h_XAx_Title.c_str() );
      theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
      theJetTrigQAHist->GetXaxis()->SetNoExponent(true);
      
      std::string h_YAx_Title="Combo/Indiv.";
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
      
      theJetTrigQAHist->SetTitle("");	
      
      if(j==1)theJetTrigQAHist->Draw("E");
      else theJetTrigQAHist->Draw("E SAME");


      if( j==1 )  {
	
	float min = theJetTrigQAHist->GetBinLowEdge( 1 );
	float max = theJetTrigQAHist->GetBinLowEdge( (int)theJetTrigQAHist->GetNbinsX() ) 
	  + theJetTrigQAHist->GetBinWidth( (int)theJetTrigQAHist->GetNbinsX() ) ;
	
	TLine* lineAtOne          = new TLine(min,1.0,max,1.0);
	lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	TLine* lineAtOneHalf      = new TLine(min,0.2,max,0.2);
	lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
	TLine* lineAtOneEtOneHalf = new TLine(min,1.8,max,1.8);
	lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
	
	lineAtOne->Draw("same");
	lineAtOneHalf->Draw("same");
	lineAtOneEtOneHalf->Draw("same");
	
      }



      


    }
    
  }//end loop over trigs for ratio
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());
  
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();
    
    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Jet p_{T}";    
    if(analysisRebin)outcanvtitle+=" Ana. Bins ";
    else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_jtpt_";
    if(analysisRebin) outcanvname+="anabins_";
    else              outcanvname+="simpbins_";
    outcanvname +=trigType;
  
    outcanv->Write((outcanvname+"_canv").c_str());    

  }
  
  //jetpad_excsp->Close();
  //jetpad_excrat->Close();
  temp_canvJetTrig->Close();
  
  
  
  
  return;
}






void printPFvCaloTrigHist_wRatio( TFile* fin , TFile* fin2, bool usedHLT100, bool usedHLTPF, bool analysisRebin, bool comboOnly,
				  std::string thePDFFileName , std::string fullJetType , 
				  std::string trigType, std::string radius , bool didJetID){
  
  if(!fin || !fin2){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  //bool analysisRebin=true;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  
  
  
  
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
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
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas(("tempPFvCaloTrigSpec"+trigType).c_str(), ("blahTrigPt"+trigType).c_str(), 1200, 1000);
  temp_canvJetTrig->cd();
  
  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING PFvCalo SPECTRA BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  TPad *jetpad_excsp = new TPad("PFvCalopad_excsp", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  jetpad_excsp->SetLogx(1);
  jetpad_excsp->SetLogy(1);
  jetpad_excsp->SetGridx(1);
  jetpad_excsp->SetGridy(1);
  
  jetpad_excsp->SetBottomMargin(0);
  jetpad_excsp->Draw();
  
  TPad *jetpad_excrat = new TPad("PFvCalopad_excrat", "Ratio Pad", 0.0, 0.05, 1.0, 0.3);
  jetpad_excrat->SetLogx(1);
  jetpad_excrat->SetLogy(0);
  jetpad_excrat->SetGridx(1);
  jetpad_excrat->SetGridy(0);
  
  jetpad_excrat->SetTopMargin(0);
  jetpad_excrat->SetBottomMargin(0.3);
  
  jetpad_excrat->Draw();
  
  jetpad_excsp->cd();
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.66,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);

  for(int j=0; j<(N_trigs); j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    //open the hists + do scaling
    std::string inHistName;
    if(trigType=="excl")inHistName="hpp_exc";
    else if(trigType=="incl")inHistName="hpp_";
    
    if(j==0)inHistName="hpp_";
    else {if(comboOnly)continue; }//temporary
    inHistName+=HLTName[j];
    if(didJetID)inHistName+="_wJetID";
    inHistName+="_"+radius+etaWidth;
    
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base"); std::cout<<std::endl;

    TH1F* theJetTrigQAHist2= (TH1F*) ( (TH1*)fin2->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist2->Print("base"); std::cout<<std::endl;
    
    if(analysisRebin) {
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_spec_anabins").c_str(), ptbins_debug    );
      theJetTrigQAHist2=(TH1F*)theJetTrigQAHist2->TH1::Rebin(nbins_pt_debug,(inHistName+"2_spec_anabins").c_str(), ptbins_debug    );}
    else {
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(10, (inHistName+"_spec_rebin10").c_str() );
      theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");      
      
      theJetTrigQAHist2=(TH1F*)theJetTrigQAHist2->TH1::Rebin(10, (inHistName+"2_spec_rebin10").c_str() );
      theJetTrigQAHist2->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          }
    divideBinWidth(theJetTrigQAHist);
    divideBinWidth(theJetTrigQAHist2);
    
    trigSpectraHistStyle(theJetTrigQAHist, j);    
    trigSpectra2HistStyle(theJetTrigQAHist2, j);    
    
    if(usedHLTPF){
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
      JetTrigLegend->AddEntry(theJetTrigQAHist2,(HLTCaloName_Leg[j]).c_str(),"lp");}
    
    else {
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
      JetTrigLegend->AddEntry(theJetTrigQAHist2,(HLTPFName_Leg[j]).c_str(),"lp");}
    
    if( j==0 ){
      
      std::string h_Title;
      if(trigType=="excl")h_Title="Excl. Jet Spectra by Trigger";
      else if(trigType=="incl")h_Title="Incl. Jet Spectra by Trigger";      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      
      std::string h_YAx_Title= AUAxTitle;
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
      
      theJetTrigQAHist->Draw("E");      
      theJetTrigQAHist2->Draw("E SAME");      
      
      float t1Loc1=0.54, t1Loc2=0.82;
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      t1->Draw();      
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),(jetCutString).c_str());      t2->Draw();      
      TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),(jetEtaCutString).c_str());      t3->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.12;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );      pp1->Draw();      
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );      pp2->Draw();  
    }
    else{
      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw("E SAME");
      trigSpectra2HistStyle(theJetTrigQAHist2, j);
      theJetTrigQAHist2->Draw("E SAME");
    }
    
  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();
  
  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------ "<<std::endl;
  std::cout<<" --------------- DRAWING JET RATIOS BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------ "<<std::endl<<std::endl<<std::endl;
  
  
  jetpad_excrat->cd();
  
  
  TH1F* theDenominator=NULL;
  TH1F* theNumerator=NULL;
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    
    
    
    //open the hists + do scaling
    if(trigType=="excl")      inHistName="hpp_exc"+HLTName[j];//+"_"+radius+etaWidth;
    else if(trigType=="incl") inHistName="hpp_"+HLTName[j];//+"_"+radius+etaWidth;
    
    if(j==0)inHistName="hpp_"+HLTName[j];    
    else{ if(comboOnly)continue;} 
    
    if(didJetID)inHistName+="_wJetID";      
    inHistName+="_"+radius+etaWidth;
    
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
    
    
    
    theDenominator = (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theDenominator->Print("base");
    
    if(analysisRebin) 
      theDenominator=(TH1F*)theDenominator->TH1::Rebin(nbins_pt_debug,(inHistName+"_rat_anabins").c_str(), ptbins_debug    );
    else {
      theDenominator=(TH1F*)theDenominator->TH1::Rebin(10, (inHistName+"_rat_rebin10").c_str() );    
      theDenominator->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");      }
    
    trigRatioHistStyle(theDenominator,j);      
    
    
    

    theNumerator = (TH1F*) ( (TH1*)fin2->Get(inHistName.c_str()) );
    if(funcDebug)theNumerator->Print("base");
    
    if(analysisRebin) 
      theNumerator=(TH1F*)theNumerator->TH1::Rebin(nbins_pt_debug,(inHistName+"2_rat_anabins").c_str(), ptbins_debug    );
    else {
      theNumerator=(TH1F*)theNumerator->TH1::Rebin(10, (inHistName+"2_rat_rebin10").c_str() );    
      theNumerator->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          
    }

    trigRatioHistStyle2(theNumerator,j);      
    
    
    theNumerator->Divide(theDenominator);
    
    
    
    theNumerator->SetAxisRange(0.0,2.0,"Y");
    
    std::string h_XAx_Title="jet p_{T} (GeV)     ";
    theNumerator->SetXTitle( h_XAx_Title.c_str() );
    theNumerator->GetXaxis()->SetMoreLogLabels(true);
    theNumerator->GetXaxis()->SetNoExponent(true);
    
    std::string h_YAx_Title;
    if(usedHLTPF)h_YAx_Title="Calo/PF";
    else h_YAx_Title="PF/Calo";
    
    
    theNumerator->SetTitle("");	
    theNumerator->SetYTitle( (h_YAx_Title).c_str() );    
    theNumerator->GetYaxis()->SetTitleOffset( .40 );        
    theNumerator->GetYaxis()->SetTitleSize( 0.08 );        
    if(j==0){
      theNumerator->Draw("E");}
    else theNumerator->Draw("E SAME");    



    if( j==0 )  {     
      float min = theNumerator->GetBinLowEdge( 1 );
      float max = theNumerator->GetBinLowEdge( (int)theNumerator->GetNbinsX() ) 
	+ theNumerator->GetBinWidth( (int)theNumerator->GetNbinsX() ) ;
      
      TLine* lineAtOne          = new TLine(min,1.0,max,1.0);
      lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
      TLine* lineAtOneHalf      = new TLine(min,0.2,max,0.2);
      lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
      TLine* lineAtOneEtOneHalf = new TLine(min,1.8,max,1.8);
      lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
      
      lineAtOne->Draw("same");
      lineAtOneHalf->Draw("same");
      lineAtOneEtOneHalf->Draw("same");
      
    }
    
    
  }
  
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());

  
  temp_canvJetTrig->Close();
  
  
  
  
  return;
}





void printTrigPtHist( TFile* fin , bool usedHLT100, bool analysisRebin,
		      std::string thePDFFileName , std::string fullJetType ,
		      std::string trigType, std::string radius , bool usedHLTPF, TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
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
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;

  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigPt"+trigType).c_str(),("blahTrigPt"+trigType).c_str() , 1200, 1000);
  temp_canvJetTrig->SetLogx(1);
  temp_canvJetTrig->SetLogy(1);
  temp_canvJetTrig->SetGridx(1);
  temp_canvJetTrig->SetGridy(1);
  temp_canvJetTrig->cd();
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING TRIGGER JET PT SPECTRA FOR TRIGGERS --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  
  //TH1F* jettrigcomb_clone=NULL;
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.70,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  if(funcDebug)     std::cout<<"drawing evt trigPt Spectra, looping over trigs"<<std::endl<<std::endl;
  
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    //if(HLTName[j]=="HLT80")continue;//debug
    //if(HLTName[j]=="HLT60")continue;//debug
    
    std::string inHistName;
    if(trigType=="excl") inHistName="is"+HLTName[j]+"_trgPt";
    else if (trigType=="incl") inHistName="Inc"+HLTName[j]+"_trgPt";    
    
    if(HLTName[j]=="HLTComb")inHistName=HLTName[j]+"_trgPt";
    
    std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    theJetTrigQAHist->GetXaxis()->SetNoExponent(true);
    theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
    theJetTrigQAHist->Print("base");
    std::cout<<std::endl<<std::endl;
    
    if(analysisRebin){
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_spec_anabins").c_str(), ptbins_debug    );
      divideBinWidth(theJetTrigQAHist);
    }
    else{
      theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_rebin").c_str());
      theJetTrigQAHist->Scale(1./((float)rebinfactor));
      //theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          
      theJetTrigQAHist->SetAxisRange(40,1400,"X");          
    }
    
    
    trigPtHistStyle(theJetTrigQAHist,j);
    
    if(usedHLTPF) JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    else JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    
    //if(j==0){
    if(HLTName[j]=="HLTComb"){
      
      std::string h_Title;
      if(trigType=="excl")h_Title="is[40,60,80] Event Trigger Jet Pt";
      else if(trigType=="incl")h_Title="Trigger Jet Pt for HLT[40,60,80]";
      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      theJetTrigQAHist->SetYTitle( (AUAxTitle).c_str() );
      theJetTrigQAHist->SetXTitle( "trigger jet p_{T} (GeV)   ");
      
      //theJetTrigQAHist->SetAxisRange(0.1,100000000, "Y");
      
      theJetTrigQAHist->Draw("E");
      //jettrigcomb_clone=(TH1F*)theJetTrigQAHist->Clone("HLTComb_cloneHist");
      
      
      float t1Loc1=0.54, t1Loc2=0.82;
      
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());
      t1->SetTextFont(63);
      t1->SetTextColor(kBlack);
      t1->SetTextSize(14);
      t1->SetLineWidth(1);
      t1->SetNDC();
      //t1->Draw();

      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0< #||{y} <5.1"));
      t2->SetTextFont(63);
      t2->SetTextColor(kBlack);
      t2->SetTextSize(14);
      t2->SetLineWidth(1);
      t2->SetNDC();
      t2->Draw();

      float pp1x=t1Loc1,pp1y=t1Loc2-0.08;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );
      pp1->SetTextFont(63);
      pp1->SetTextColor(kBlack);
      pp1->SetTextSize(15);
      pp1->SetLineWidth(1);
      pp1->SetNDC();
      pp1->Draw();
      
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );
      pp2->SetTextFont(63);
      pp2->SetTextColor(kBlack);
      pp2->SetTextSize(15);
      pp2->SetLineWidth(1);
      pp2->SetNDC();
      pp2->Draw();  
    }//if HLTComb    
    else {
      theJetTrigQAHist->Draw("E SAME");
    }//if not HLTComb
    
    
    
  }//end loop over trigs for spectra
  
  //jettrigcomb_clone->Draw("E SAME");
  JetTrigLegend->Draw();
  




  temp_canvJetTrig->Print(thePDFFileName.c_str());
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();

    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Trig. Jet p_{T}";    
    if(analysisRebin)outcanvtitle+=" Ana. Bins ";
    else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_trgjtpt_";
    if(analysisRebin) outcanvname+="anabins_";
    else              outcanvname+="simpbins_";
    outcanvname +=trigType;

    outcanv->Write((outcanvname+"_canv").c_str());    

    //std::string outcanvtitle="HLTak4PFJet",outcanvname="HLTak4PFJet";//= "HLT"fullJetType+
    //if(usedHLT100)outcanvname+="40to100_";
    //else outcanvname+="40to80_";
    //if(usedHLT100)outcanvtitle+="40 to 100 ";
    //else outcanvtitle+="40 to 80 ";
    //outcanvname+=fullJetType+"Jets_";
    //outcanvtitle+=fullJetType+"Jets ";    
    //outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    //outcanv->Write((outcanvname+"canv").c_str());    
  }  
  temp_canvJetTrig->Close();

  return;
}







void printTrigEtaHist( TFile* fin , bool usedHLT100,
		       std::string thePDFFileName , std::string fullJetType ,
		       std::string trigType, std::string radius , bool usedHLTPF, TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
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
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  
  
  
  
  
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigEta"+trigType).c_str(),("blahTrigEta"+trigType).c_str() , 1200, 1000);
  temp_canvJetTrig->SetLogx(0);
  temp_canvJetTrig->SetLogy(1);
  temp_canvJetTrig->SetGridx(1);
  temp_canvJetTrig->SetGridy(1);
  temp_canvJetTrig->cd();

  

  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING TRIGGER JET ETA SPECTRA FOR TRIGGERS --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  
  
  
  TLegend* JetTrigLegend=new TLegend(0.50,0.35,0.65,0.46, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  if(funcDebug)     std::cout<<"drawing evt trigEta Spectra, looping over trigs"<<std::endl<<std::endl;
  
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    if(trigType=="excl") inHistName="is"+HLTName[j]+"_trgEta";
    else if (trigType=="incl") inHistName="Inc"+HLTName[j]+"_trgEta";
    if(j==0)inHistName=HLTName[j]+"_trgEta";
    
    if(funcDebug)std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base");
    std::cout<<std::endl<<std::endl;
    
    
    theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(2, (inHistName+"_rebin").c_str());
    theJetTrigQAHist->Scale(1./2.);    
    


    trigEtaHistStyle(theJetTrigQAHist,j);
    
    theJetTrigQAHist->SetAxisRange(800,40000000, "Y");

    if(usedHLTPF) JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    else JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    
    if(j==0){
      
      std::string h_Title;//   ="Exclusive Event Trigger Jet Eta";//, Data ONLY";
      if(trigType=="excl")h_Title="is[40,60,80] Event Trigger Jet Eta";
      else if(trigType=="incl")h_Title="Trigger Jet Eta for HLT[40,60,80]";
      
      
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );
      theJetTrigQAHist->SetYTitle( (AUAxTitle).c_str() );
      theJetTrigQAHist->SetXTitle( "trigger jet #||{y}   ");

      theJetTrigQAHist->Draw("E");



      float t1Loc1=0.35, t1Loc2=0.50;
      //TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      //t1->Draw();      
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0< #||{y} <5.1"));      t2->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.08;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );      pp1->Draw();
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );      pp2->Draw();  
    }
    
    else
      theJetTrigQAHist->Draw("E SAME");
    
    
    
  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();
    
    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Trig. Jet #eta ";    
    //if(analysisRebin)outcanvtitle+=" Ana. Bins ";
    //else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_trgjteta_";
    //if(analysisRebin) outcanvname+="anabins_";
    //else              outcanvname+="simpbins_";
    outcanvname +=trigType;
    
    outcanv->Write((outcanvname+"_canv").c_str());    

  }

  temp_canvJetTrig->Close();

  return;
}





TH1F* makeAsymmEtaHist(TH1F* theEtaHist,std::string trigType) {
  //  bool funcDebug=false;
  std::string hTitle=theEtaHist->GetName();
  int NBinsNewHist=theEtaHist->GetNbinsX()/2;
  float min = 0. ; // = theEtaHist->GetBinLowEdge( 1 );   
  float max = theEtaHist->GetBinLowEdge( (int)theEtaHist->GetNbinsX() )
    + theEtaHist->GetBinWidth( (int)theEtaHist->GetNbinsX() ) ;
  
  if(theEtaHist->GetNbinsX()%2 != 0) 
    std::cout<<std::endl<<"WARNING: AsymmEtaHist will have odd # of bins! Will not look good!"<<std::endl<<std::endl;
  
  TH1F* newAsymmEtaHist= new TH1F( (hTitle+"_asymmHist_"+trigType).c_str(), "",   NBinsNewHist, min,max);
  //newAsymmEtaHist->Print("base");
  
  for(int i = 1; i<=NBinsNewHist; i++){
    

    float negSideContent, negSideErr;
    float posiSideContent, posiSideErr;
    
    negSideContent=theEtaHist->GetBinContent(NBinsNewHist-i+1);
    negSideErr=theEtaHist      ->GetBinError(NBinsNewHist-i+1);
    
    posiSideContent=theEtaHist->GetBinContent( NBinsNewHist+i );
    posiSideErr=theEtaHist      ->GetBinError( NBinsNewHist+i );
    
    if(i==(NBinsNewHist*2-i+1) )continue;
    
    float asymmBinContent,asymmBinErr; 
    asymmBinContent=posiSideContent/negSideContent;
    asymmBinErr=
      (posiSideContent/negSideContent)* std::sqrt( 
						  (negSideErr/negSideContent)*(negSideErr/negSideContent) +
						  (posiSideErr/posiSideContent)*(posiSideErr/posiSideContent)    
						   );   
    
        
    if( asymmBinContent != asymmBinContent ||
	asymmBinErr != asymmBinErr
	){
      continue;
    }
    else{

      newAsymmEtaHist->SetBinContent( i, asymmBinContent );
      newAsymmEtaHist->SetBinError(   i, asymmBinErr   );
    }
    
    
    
    
  }
  
  
  return newAsymmEtaHist;
  
}




void printTrigEtaAsymmHist( TFile* fin , bool usedHLT100,
			    std::string thePDFFileName , std::string fullJetType ,
			    std::string trigType, std::string radius , bool usedHLTPF, TFile* fout=NULL){
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return;}
  
  // PUT JET CUTS SECTION IN HERE
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
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
  
  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
  
  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  std::string sampleDescString=sqrts2k15ppString;
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  
  
  
  
  
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigAsymmEta"+trigType).c_str(),("blahTrigAsymmEta"+trigType).c_str() , 1200, 1000);
  temp_canvJetTrig->SetLogx(0);
  temp_canvJetTrig->SetLogy(0);
  temp_canvJetTrig->SetGridx(1);
  temp_canvJetTrig->SetGridy(1);
  temp_canvJetTrig->cd();

  
  
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING TRIGGER JET ASYMM ETA HISTS FOR TRIGGERS --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;
  
  
  
  
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.70,0.85,0.86, NULL,"brNDC");
  JetTrigLegend->SetFillStyle(0);
  
  if(funcDebug)     std::cout<<"drawing evt trigEta AsymmSpectra, looping over trigs"<<std::endl<<std::endl;
  
  for(int j=0; j<N_trigs; j++){
    
    //open the hists + do scaling
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    
    std::string inHistName;
    if(trigType=="excl") inHistName="is"+HLTName[j]+"_trgEta";
    else if (trigType=="incl") inHistName="Inc"+HLTName[j]+"_trgEta";    
    if(j==0)inHistName=HLTName[j]+"_trgEta";
    if(funcDebug)std::cout<<std::endl<<"inHistName="<<inHistName<<std::endl;
    
    TH1F* theJetTrigQAHist= (TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base");
    std::cout<<std::endl<<std::endl;
    
    theJetTrigQAHist=(TH1F*)theJetTrigQAHist->TH1::Rebin(5, (inHistName+"_rebin").c_str());
    theJetTrigQAHist->Scale(1./5.);                
    
    TH1F* theJetTrigAsymmQAHist= makeAsymmEtaHist(theJetTrigQAHist,trigType) ; //(TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigAsymmQAHist->Print("base");
    trigEtaHistStyle(theJetTrigAsymmQAHist,j);    
    
    if(usedHLTPF) JetTrigLegend->AddEntry(theJetTrigAsymmQAHist,(HLTPFName_Leg[j]).c_str(),"lp");
    else JetTrigLegend->AddEntry(theJetTrigAsymmQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");
    
    theJetTrigAsymmQAHist->SetAxisRange(0.0,2.0,"Y");
    theJetTrigAsymmQAHist->SetAxisRange(0.0,5.0,"X");
    if(j==0){
      
      std::string h_Title;//   ="Exclusive Event Trigger Jet Eta";//, Data ONLY";
      if(trigType=="excl")h_Title="is[40,60,80] Event Trigger Jet Eta Asymm";
      else if(trigType=="incl")h_Title="Trigger Jet Eta Asymm for HLT[40,60,80]";
      
      
      theJetTrigAsymmQAHist->SetTitle (    h_Title.c_str() );
      theJetTrigAsymmQAHist->SetYTitle( "+y bin / -y bin" );
      theJetTrigAsymmQAHist->SetXTitle( "trigger jet #||{y}  ");
      
      //theJetTrigAsymmQAHist->SetAxisRange(0.1,100000000, "Y");
      theJetTrigAsymmQAHist->Draw("E");
      
      float t1Loc1=0.54, t1Loc2=0.82;
      //TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      //t1->Draw();
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0< #||{y} <5.1"));      t2->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.08;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );       pp1->Draw();
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );   pp2->Draw();  
    }    
    else
      theJetTrigAsymmQAHist->Draw("E SAME");    

  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();  

  temp_canvJetTrig->Print(thePDFFileName.c_str());  
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();
    
    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Trig. Jet #eta Asymm. ";    
//    if(analysisRebin)outcanvtitle+=" Ana. Bins ";
//    else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;

    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_trgjteta_asymm_";
//    if(analysisRebin) outcanvname+="anabins_";
//    else              outcanvname+="simpbins_";
    outcanvname +=trigType;
    
    outcanv->Write((outcanvname+"_canv").c_str());    


//    std::string outcanvtitle="HLTak4PFJet",outcanvname="HLTak4PFJet";//= "HLT"fullJetType+
//    if(usedHLT100)outcanvname+="40to100_";
//    else outcanvname+="40to80_";
//    if(usedHLT100)outcanvtitle+="40 to 100 ";
//    else outcanvtitle+="40 to 80 ";
//    outcanvname+=fullJetType+"Jets_";
//    outcanvtitle+=fullJetType+"Jets ";    
//    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
//    outcanv->Write((outcanvname+"canv").c_str());    
  }
  temp_canvJetTrig->Close();
  return;

}




