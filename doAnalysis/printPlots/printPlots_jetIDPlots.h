
//const double JETQACANVX=1300.;
//const double JETQACANVY=JETQACANVX;
//void printJetIDHist( TFile* fin , 
//		     int j, bool isData,
//		     std::string thePDFFileName , std::string fullJetType, 
//		     long double theLumi,  
//		     TFile* fout=NULL) {

void printJetIDHist( TFile* fin, int j, bool isData, int NetabinstoDraw,
		     std::string thePDFFileName , std::string fullJetType, 
		     long double theLumi  ,TFile* fout=NULL) {
  bool funcDebug=false;
  
  //std::string jetIDInt;
  //if(doJetIDPlots)jetIDInt="1";
  //else jetIDInt="0";
  
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  bool isConstitHist= (bool)( j>=jetConstits_varStart && j<dijet_varStart );  
  
  // STRINGS FOR HISTS
  TH1F* genjtptCut_h= NULL;
  std::string genjtptCut_str = "";
  std::string genJetCutString="";
  if(!isData){
    genjtptCut_h=(TH1F*)fin->Get( "hGenJetPtCut" );
    genjtptCut_str=std::to_string( (int) genjtptCut_h->GetMean() );
    std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;  
    genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
    std::cout<<"genJetCutString="<<genJetCutString<<std::endl;
  }
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  std::string jtptCut_str=jtptQACut_str;
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;  
  
  float etaBinWidth=1.;
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  std::string jetIDEfftitle="Jet ID Eff.";

  bool isybin=((var[j]=="jtpt" && j==0)
	       ||var[j]=="jtE"||var[j]=="jtm");
  
  TH1F* theJetIDHist[Netabins]={};
  TH1F* theNoJetIDHist[Netabins]={};
  TH1F* theRatio[Netabins]={};
  std::string inHistName="hJetQA_JetIDEff_"+var[j];
  if(isybin)
    inHistName+="_"+absetabins_str[0]+"y"+absetabins_str[NetabinstoDraw];
  else
    inHistName+="_"+absetabins_str[0]+"eta"+absetabins_str[NetabinstoDraw];
  
  for(int etabin=0; etabin<NetabinstoDraw;etabin++){

    //hJetQA_0wJetID_muSum_etabin4
    std::string inNoJetIDHist="hJetQA_0wJetID_"+var[j];
    std::string inJetIDHist  ="hJetQA_1wJetID_"+var[j];
    if(isybin){
      inNoJetIDHist+="_ybin"+std::to_string(etabin);
      inJetIDHist  +="_ybin"+std::to_string(etabin);
    }
    else {
      inNoJetIDHist+="_etabin"+std::to_string(etabin);
      inJetIDHist  +="_etabin"+std::to_string(etabin);
    }

    std::cout<<" opening input Data "<<inNoJetIDHist<<std::endl<<std::endl;          
    theNoJetIDHist[etabin] = (TH1F*) ( (TH1*)fin->Get(inNoJetIDHist.c_str()) );
    theNoJetIDHist[etabin] = (TH1F*) (theNoJetIDHist[etabin]->Clone((inNoJetIDHist+"_"+std::to_string(j)+"clone").c_str()) );
    if(!theNoJetIDHist[etabin]){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
      return;}
    
    theNoJetIDHist[etabin]->Scale( 1./theLumi);
    theNoJetIDHist[etabin]->Scale( 1./etaBinWidth);//eta bin width scaling
    theNoJetIDHist[etabin]->ClearUnderflowAndOverflow();  
    
    //std::string inJetIDHist=inHistName+std::to_string(etabin);
    //if(var[j]=="jtpt_forRes" || var[j]=="jtpt_L2Res" || var[j]=="jtpt_L3Res"){
    //  inJetIDHist="hJetQA_"+std::to_string((int)doJetIDPlots)+"wJetID_jtpt_etabin"+std::to_string(etabin);
    //}
    
    std::cout<<" opening input MC "<<inJetIDHist<<std::endl<<std::endl;      
    theJetIDHist[etabin]= (TH1F*) ( (TH1*)fin->Get(inJetIDHist.c_str()) );  
    theJetIDHist[etabin]= (TH1F*) (theJetIDHist[etabin]->Clone((inJetIDHist+"_"+std::to_string(j)+"clone").c_str()) );
    if(!theJetIDHist[etabin]){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
      return;}
    //theJetIDHist[etabin]->Scale(MCscale4Units)
    theJetIDHist[etabin]->Scale( 1./theLumi);
    theJetIDHist[etabin]->Scale(1./etaBinWidth);  
    theJetIDHist[etabin]->ClearUnderflowAndOverflow();      
    
    // MORE SCALING
    if(var[j]=="jtpt" || var[j]=="rawpt" ||
       var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res" ||
       var[j]=="leadJetPt" || var[j]=="subleadJetPt") {
      
      theNoJetIDHist[etabin]=(TH1F*)theNoJetIDHist[etabin]->TH1::Rebin(nbins_pt_debug, (
									    ((std::string)theNoJetIDHist[etabin]->GetName()+"_rebin")).c_str(), ptbins_debug);
      divideBinWidth(theNoJetIDHist[etabin]);      
      
      theJetIDHist[etabin]=(TH1F*)theJetIDHist[etabin]->TH1::Rebin(nbins_pt_debug, (
									((std::string)theJetIDHist[etabin]->GetName()+"_rebin")).c_str(), ptbins_debug);
      divideBinWidth(theJetIDHist[etabin]); 
      
    }//if jet pT hist
    
    else {
      
      if(j>=jetConstits_varStart && j<dijet_varStart){
	//i'm sorry, C++ gods, this is because the hists' bins in readForests_ppMC/Data_jetPlots are defined differently
	if(funcDebug)std::cout<<"BEFORE REBIN theNoJetIDHist nbins="<<theJetIDHist[etabin]->GetNbinsX()<<std::endl;      
	if(funcDebug)std::cout<<"BEFORE REBIN theNoJetIDHist xmin="<<theJetIDHist[etabin]->GetBinLowEdge(1)<<std::endl;      
	if(funcDebug)std::cout<<"BEFORE REBIN theNoJetIDHist xmax="<<
		       theJetIDHist[etabin]->GetBinLowEdge(theJetIDHist[etabin]->GetNbinsX()) +
		       theJetIDHist[etabin]->GetBinWidth(theJetIDHist[etabin]->GetNbinsX()) 
			      <<std::endl;      
	if(funcDebug)std::cout<<std::endl;
	theJetIDHist[etabin]->TH1::Rebin(var_xAx_reBin[j]);
	theJetIDHist[etabin]->Scale( 1./theJetIDHist[etabin]->GetBinWidth(1) );          
	theJetIDHist[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
	//theJetIDHist[etabin]->GetXaxis()->SetRange(jetQAxmin[j],jetQAxmax[j]); 
	std::cout<<"theNoJetIDHist nbins="<<theJetIDHist[etabin]->GetNbinsX()<<std::endl;      
	std::cout<<"theNoJetIDHist xmin="<<theJetIDHist[etabin]->GetBinLowEdge(1)<<std::endl;      
	std::cout<<"theNoJetIDHist xmax="<<
	  theJetIDHist[etabin]->GetBinLowEdge(theJetIDHist[etabin]->GetNbinsX()) +
	  theJetIDHist[etabin]->GetBinWidth(theJetIDHist[etabin]->GetNbinsX()) 
		 <<std::endl;
	std::cout<<std::endl;      
	theJetIDHist[etabin]->ClearUnderflowAndOverflow();     
	
	const int mcjetqa_nbins=theJetIDHist[etabin]->GetNbinsX();
	double mcjetqa_bins[mcjetqa_nbins+1]={0.};
	for(int i=1; i<=mcjetqa_nbins; i++)
	  mcjetqa_bins[i-1]=theJetIDHist[etabin]->GetBinLowEdge(i);
	mcjetqa_bins[mcjetqa_nbins]=theJetIDHist[etabin]->GetBinLowEdge(mcjetqa_nbins) + theJetIDHist[etabin]->GetBinWidth(mcjetqa_nbins);
	if(funcDebug)      std::cout<<"mcjetqa_nbins="<<mcjetqa_nbins<<std::endl;
	if(funcDebug)      std::cout<<"mcjetqa_bins[0]="<<mcjetqa_bins[0]<<std::endl;
	if(funcDebug)      std::cout<<"mcjetqa_bins["<< mcjetqa_nbins <<"]="<<mcjetqa_bins[mcjetqa_nbins]<<std::endl;
	if(funcDebug)      std::cout<<std::endl;
	
	if(funcDebug)std::cout<<"BEFORE REBIN theJetIDHist nbins="<<theNoJetIDHist[etabin]->GetNbinsX()<<std::endl;      
	if(funcDebug)std::cout<<"BEFORE REBIN theJetIDHist xmin="<<theNoJetIDHist[etabin]->GetBinLowEdge(1)<<std::endl;      
	if(funcDebug)std::cout<<"BEFORE REBIN theJetIDHist xmax="<<
		       theNoJetIDHist[etabin]->GetBinLowEdge(theNoJetIDHist[etabin]->GetNbinsX()) +
		       theNoJetIDHist[etabin]->GetBinWidth(theNoJetIDHist[etabin]->GetNbinsX()) 
			      <<std::endl;      
	if(funcDebug)std::cout<<std::endl;
	//theNoJetIDHist[etabin]->TH1::Rebin(var_xAx_reBin[j]);
	theNoJetIDHist[etabin]=(TH1F*)theNoJetIDHist[etabin]->TH1::Rebin(mcjetqa_nbins, 
							     ( ((std::string)theNoJetIDHist[etabin]->GetName())+"_mcjetqarebin").c_str(), 
							     (const Double_t*)mcjetqa_bins);
	theNoJetIDHist[etabin]->ClearUnderflowAndOverflow();  
	theNoJetIDHist[etabin]->Scale( 1./theNoJetIDHist[etabin]->GetBinWidth(1) );          //its a constant binning so this is fine
	
	theNoJetIDHist[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");   
	//theNoJetIDHist[etabin]->GetXaxis()->SetRange(jetQAxmin[j],jetQAxmax[j]);
	
	if(funcDebug)      std::cout<<"theJetIDHist nbins="<<theNoJetIDHist[etabin]->GetNbinsX()<<std::endl;      
	if(funcDebug)      std::cout<<"theJetIDHist xmin="<<theNoJetIDHist[etabin]->GetBinLowEdge(1)<<std::endl;      
	if(funcDebug)      std::cout<<"theJetIDHist xmax="<<
			     theNoJetIDHist[etabin]->GetBinLowEdge(theNoJetIDHist[etabin]->GetNbinsX()) +
			     theNoJetIDHist[etabin]->GetBinWidth(theNoJetIDHist[etabin]->GetNbinsX()) 
				    <<std::endl;      
	theJetIDHist[etabin]->Scale( theNoJetIDHist[etabin]->Integral()/theJetIDHist[etabin]->Integral() );	        
	
	//assert(false);
	
      }//end if constituent hist
      else{//if not jet spectra, or constituent, e.g. jeteta, jetphi, jety, jetm, jetE
	theNoJetIDHist[etabin]->TH1::Rebin(var_xAx_reBin[j]);
	theNoJetIDHist[etabin]->Scale( 1./theNoJetIDHist[etabin]->GetBinWidth(1) );          
	theNoJetIDHist[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
	//theNoJetIDHist[etabin]->GetXaxis()->SetRange(jetQAxmin[j],jetQAxmax[j]);
	
	theJetIDHist[etabin]->TH1::Rebin(var_xAx_reBin[j]);
	theJetIDHist[etabin]->Scale( 1./theJetIDHist[etabin]->GetBinWidth(1) );          
	theJetIDHist[etabin]->Scale( theNoJetIDHist[etabin]->Integral()/theJetIDHist[etabin]->Integral() );	        
	theJetIDHist[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
	//theJetIDHist[etabin]->GetXaxis()->SetRange(jetQAxmin[j],jetQAxmax[j]); 
      }//end if dijet hist
      
    }//end if not a jet pt hist
    
    
    // STYLE  
    //dataHistStyle(theNoJetIDHist[etabin]);
    //MCHistStyle(theJetIDHist[etabin]);
    nojetIDHistStyle(theNoJetIDHist[etabin]);
    jetIDHistStyle(theJetIDHist[etabin]);
    
    
    


    
    // TITLES
    std::string h_Title = PDStatsString_2; //for APS DNP
    std::string h_XAx_Title= var_xAx_Titles[j];
    //std::string h_YAx_Title= dcrossSectionAxTitle;
    //std::string h_YAx_Title= AUAxTitle;//for APS DNP
    std::string h_YAx_Title= ddcrossSectionAxTitle;//for APS DNP
    if((var[j]=="jtE"||var[j]=="jtm"||var[j]=="jty"||var[j]=="jteta" || var[j]=="jtphi" || isConstitHist) )h_YAx_Title=AUAxTitle;
    
    theNoJetIDHist[etabin]->SetTitle (    h_Title.c_str() );
    theNoJetIDHist[etabin]->SetXTitle( h_XAx_Title.c_str() );
    theNoJetIDHist[etabin]->SetYTitle( h_YAx_Title.c_str() );
    
    theJetIDHist[etabin]->SetTitle (    h_Title.c_str() );
    theJetIDHist[etabin]->SetXTitle( h_XAx_Title.c_str() );
    theJetIDHist[etabin]->SetYTitle( h_YAx_Title.c_str() );
    
    
    // GET RATIO PLOT READY FOR BOTTOM PANEL
    theRatio[etabin]=(TH1F*)theJetIDHist[etabin]->Clone(  ("MCJetHistClone4Ratio_"+var[j]+"_etabin"+std::to_string(etabin)).c_str());  
    
    //general
    theRatio[etabin]->SetTitle("");
    theRatio[etabin]->SetLineColor( theMCLineColor );
    theRatio[etabin]->SetFillColorAlpha(kRed, 0.);
    theRatio[etabin]->SetFillStyle(0);

    //yaxis, range, title, labels
    theRatio[etabin]->SetAxisRange(0.25,1.75,"Y");  
    theRatio[etabin]->GetYaxis()->CenterTitle(true);
    theRatio[etabin]->GetYaxis()->SetTitleSize(17);
    theRatio[etabin]->GetYaxis()->SetTitleFont(43);
    theRatio[etabin]->GetYaxis()->SetTitleOffset(2.5);
    theRatio[etabin]->SetYTitle(jetIDEfftitle.c_str() );    
    theRatio[etabin]->GetYaxis()->SetLabelFont(43); 
    theRatio[etabin]->GetYaxis()->SetLabelSize(13);  
    
    //xaxis, range, title, labels
    //  if((var[j]=="jteta" || var[j]=="jtphi" || isConstitHist) )
    //    theRatio[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
    theRatio[etabin]->GetXaxis()->SetTitleSize(19);
    theRatio[etabin]->GetXaxis()->SetTitleFont(43);
    theRatio[etabin]->GetXaxis()->SetTitleOffset(7.0);
    theRatio[etabin]->SetXTitle( h_XAx_Title.c_str() );  
    theRatio[etabin]->GetXaxis()->SetLabelFont(43); 
    theRatio[etabin]->GetXaxis()->SetLabelSize(15);
    
  }//end etabin loop for grabbing hists and processing the #'s for the histograms
  
  
  //this etabin loop for drawing stuff only
  
  TPad* jetpad1[Netabins]={};
  TPad* jetpad2[Netabins]={};
  
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", (int)JETQACANVX, (int)JETQACANVY); //This is for the Jet QA Evt plots without logs
  //temp_canvJet->SetCanvasSize(JETQACANVX,JETQACANVY);
  //temp_canvJet->SetWindowSize(JETQACANVX+(JETQACANVX-temp_canvJet->GetWw()),JETQACANVY+(JETQACANVY-temp_canvJet->GetWh()));
  temp_canvJet->cd();
  temp_canvJet->Divide(2,2);//lotta wasted space
  //temp_canvJet->TPad::Divide(2,2,0.005,0.005);//no change
  //temp_canvJet->TPad::Divide(2,2,0.002,0.002);//no change
  //temp_canvJet->TPad::Divide(2,2,-.1,-.1);//ROOT optimizes for me, no change

  //temp_canvJet->Divide(2,2);
  
  
  for(int etabin=0; etabin<NetabinstoDraw;etabin++){    

    
    // TITLES
    std::string h_Title = PDStatsString_2; //for APS DNP
    std::string h_XAx_Title= var_xAx_Titles[j];
    //std::string h_YAx_Title= dcrossSectionAxTitle;
    //std::string h_YAx_Title= AUAxTitle;//for APS DNP
    std::string h_YAx_Title= ddcrossSectionAxTitle;//for APS DNP
    if((var[j]=="jtm"||var[j]=="jtE"||var[j]=="jty"||var[j]=="jteta" || var[j]=="jtphi" || isConstitHist) )h_YAx_Title=AUAxTitle;
    
    std::string jetEtaCutString=absetabins_str[etabin]+" <#||{";
    if(isybin) jetEtaCutString+="y";
    else       jetEtaCutString+="#eta";
    jetEtaCutString+="} < "+absetabins_str[etabin+1];
    
    std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;    
    std::cout<<"JETQACANVX="<<JETQACANVX<<std::endl;
    std::cout<<"JETQACANVY="<<JETQACANVY<<std::endl;
    
    bool isJetSpectra=      (var[j]=="jtpt"      || var[j]=="rawpt"          ||
			     var[j]=="jtpt_forRes"|| var[j]=="jtpt_L2Res"|| var[j]=="jtpt_L3Res"|| 
			     var[j]=="leadJetPt" || var[j]=="subleadJetPt"     );
    
    if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
    temp_canvJet->cd(etabin+1);
    
    
    jetpad1[etabin] = new TPad(("jetpad1_"+std::to_string(etabin)).c_str(), "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
    //std::cout<<"DEFAULT MARGINS ARE...."<<std::endl;//always 0.1
    //std::cout<<"LEFT  : "<<jetpad1[etabin]->GetLeftMargin()<<std::endl;
    //std::cout<<"RIGHT : "<<jetpad1[etabin]->GetRightMargin()<<std::endl;
    //std::cout<<"TOP   : "<<jetpad1[etabin]->GetTopMargin()<<std::endl;
    //std::cout<<"BOTTOM: "<<jetpad1[etabin]->GetBottomMargin()<<std::endl;
    jetpad1[etabin]->SetGridx(0);
    jetpad1[etabin]->SetGridy(0);  
    jetpad1[etabin]->SetLeftMargin(0.09);
    jetpad1[etabin]->SetRightMargin(0.04);
    jetpad1[etabin]->SetTopMargin(0.07);//good idea?
    jetpad1[etabin]->SetBottomMargin(0);//do not change
    if(isJetSpectra){
      jetpad1[etabin]->SetLogx(1);  
      jetpad1[etabin]->SetLogy(1);  
    }
    else if(var[j]=="jtE"){
      jetpad1[etabin]->SetLogx(1);  
      jetpad1[etabin]->SetLogy(1);  
    } 
    else {
      jetpad1[etabin]->SetLogx(0);  
      jetpad1[etabin]->SetLogy(1);  
    }
    jetpad1[etabin]->Draw(); 
 
    jetpad2[etabin] = new TPad(("jetpad2_"+std::to_string(etabin)).c_str(), "Ratio Pad"  , 0.0, 0.0, 1.0, 0.3);
    jetpad2[etabin]->SetGridx(1);
    jetpad2[etabin]->SetGridy(0);  
    jetpad2[etabin]->SetLeftMargin(0.09);
    jetpad2[etabin]->SetRightMargin(0.04);
    jetpad2[etabin]->SetTopMargin(0);//do not change
    jetpad2[etabin]->SetBottomMargin(0.3);//maybe change
    if(isJetSpectra){
      jetpad2[etabin]->SetLogx(1);  
      jetpad2[etabin]->SetLogy(0);  
    }
    else if(var[j]=="jtE"){
      jetpad2[etabin]->SetLogx(1);  
      jetpad2[etabin]->SetLogy(0);  
    }
    else{
      jetpad2[etabin]->SetLogx(0);  
      jetpad2[etabin]->SetLogy(0);  
    }
    jetpad2[etabin]->Draw();
  

  
    // variable specific draw orders, canv, etc.
    if(var[j]=="jtpt"||var[j]=="rawpt" || 
       var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      )         {
      
      jetpad1[etabin]->cd();
      
      float histmin1=getnonzeromin(theNoJetIDHist[etabin]);
      float histmin2=getnonzeromin(  theJetIDHist[etabin]);
      float histmax1=theNoJetIDHist[etabin]->GetMaximum();
      float histmax2=  theJetIDHist[etabin]->GetMaximum();
      if(histmin1<histmin2)
	theNoJetIDHist[etabin]->SetMinimum(histmin1*0.7);
      else 
	theNoJetIDHist[etabin]->SetMinimum(histmin2*0.7);
      if(histmax1>histmax2)
	theNoJetIDHist[etabin]->SetMaximum(histmax1*3.);
      else 
	theNoJetIDHist[etabin]->SetMaximum(histmax1*3.);

      theNoJetIDHist[etabin]->Draw("HIST E"); 
      theJetIDHist[etabin]->Draw("HIST E SAME "); 
      
      
      float t1Loc1=0.50, t1Loc2=0.84;
      TLatex* t1= makeTLatex(  t1Loc1  , t1Loc2 ,  fullJetType+"Jets" );    t1->Draw(); t1Loc2-=0.05;
      TLatex *t2= makeTLatex( (t1Loc1) , t1Loc2 , (jetCutString)      );    t2->Draw(); t1Loc2-=0.05;
      if(!isData){
      TLatex *t3= makeTLatex( (t1Loc1) , t1Loc2 , (genJetCutString)   );    t3->Draw(); t1Loc2-=0.05; }
      TLatex *t4= makeTLatex( (t1Loc1) , t1Loc2 , (jetEtaCutString)   );    t4->Draw();	    
      //if(doJetIDPlots){
      //	TLatex *t5= makeTLatex( (t1Loc1), (t1Loc2-.20), "Jet ID Applied" );      t5->Draw();	}
      
      float legx1=0.70, legx2=0.89;
      float legy1=0.77, legy2=0.89;
      TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
      theJetQALeg->SetFillStyle(0);
      theJetQALeg->SetBorderSize(0.);
      std::string legdesc="";
      if(isData)legdesc="Data ";
      else      legdesc="PY8 MC ";
      theJetQALeg->AddEntry(theNoJetIDHist[etabin],(legdesc + "No Jet ID").c_str(),"lp");
      theJetQALeg->AddEntry(theJetIDHist[etabin],  (legdesc + "With Jet ID").c_str()  ,"lpf");
      theJetQALeg->Draw(); 
      
      
      jetpad2[etabin]->cd();
      
      
      //TH1F* theDenom=(TH1F*)theNoJetIDHist[etabin]->Clone(("DataHistClone4Ratio_"+var[j]).c_str());
      
      theRatio[etabin]->GetXaxis()->SetMoreLogLabels(true);
      theRatio[etabin]->GetXaxis()->SetNoExponent(true);        
      theRatio[etabin]->Divide(theNoJetIDHist[etabin]);    
      theRatio[etabin]->Draw("L P E");
      
      TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,ptbins_debug[nbins_pt_debug],1.0); 
      lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
      //TLine* lineAtOneHalf      = new TLine(ptbins_debug[0],0.5,ptbins_debug[nbins_pt_debug],0.5); 
      //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
      //TLine* lineAtOneEtOneHalf = new TLine(ptbins_debug[0],1.5,ptbins_debug[nbins_pt_debug],1.5); 
      //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
      
      lineAtOne->Draw("same");
      //lineAtOneHalf->Draw("same");
      //lineAtOneEtOneHalf->Draw("same");
      
      
      
      
    
    }//end if jet pt spectra
    
    else if(var[j]=="jteta" || var[j]=="jtm"||var[j]=="jtE"||
	    var[j]=="jtphi"	|| var[j]=="jty"||
	    isConstitHist){ 

      if(var[j]=="jtE"){
	theRatio[etabin]->GetXaxis()->SetMoreLogLabels(true);
	theRatio[etabin]->GetXaxis()->SetNoExponent(true);        
      }

      
      float histmin1=getnonzeromin(theNoJetIDHist[etabin]);
      float histmin2=getnonzeromin(  theJetIDHist[etabin]);
      float histmax1=theNoJetIDHist[etabin]->GetMaximum();
      float histmax2=  theJetIDHist[etabin]->GetMaximum();
      if(histmin1<histmin2)
	theJetIDHist[etabin]->SetMinimum(histmin1*0.7);
      else 
	theJetIDHist[etabin]->SetMinimum(histmin2*0.7);
      if(histmax1>histmax2)
	theJetIDHist[etabin]->SetMaximum(histmax1*3.);
      else 
	theJetIDHist[etabin]->SetMaximum(histmax1*3.);

      jetpad1[etabin]->cd();
      
      theJetIDHist[etabin]->Draw("HIST E"); 
      theNoJetIDHist[etabin]->Draw("E SAME");  
      
      float t1Loc1=0.54, t1Loc2=0.82;
      //float t1Loc1=0.50, t1Loc2=0.72; 	    
      TLatex* t1= makeTLatex(  t1Loc1  , t1Loc2 ,  fullJetType+"Jets" );    t1->Draw(); t1Loc2-=0.05;
      TLatex *t2= makeTLatex( (t1Loc1) , t1Loc2 , (jetCutString)      );    t2->Draw(); t1Loc2-=0.05;
      if(!isData){
      TLatex *t3= makeTLatex( (t1Loc1) , t1Loc2 , (genJetCutString)   );    t3->Draw(); t1Loc2-=0.05; }
      TLatex *t4= makeTLatex( (t1Loc1) , t1Loc2 , (jetEtaCutString)   );    t4->Draw();	    
      //if(doJetIDPlots){
      //	TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.20),"Jet ID Applied");      t5->Draw();	
      //}
      
      float legx1=0.70, legx2=0.89;
      float legy1=0.77, legy2=0.89;

      TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
      theJetQALeg->SetFillStyle(0);
      theJetQALeg->SetBorderSize(0.);
      std::string legdesc="";
      if(isData)legdesc="Data ";
      else      legdesc="PY8 MC ";
      theJetQALeg->AddEntry(theNoJetIDHist[etabin],(legdesc + "No Jet ID").c_str(),"lp");
      theJetQALeg->AddEntry(theJetIDHist[etabin],  (legdesc + "With Jet ID").c_str()  ,"lpf");
      //theJetQALeg->AddEntry(theNoJetIDHist[etabin],"Data","lp");
      //theJetQALeg->AddEntry(theJetIDHist[etabin],  "P8 MC"  ,"lp");
      theJetQALeg->Draw();
    
      //Ratio Plots    
      jetpad2[etabin]->cd();
      theRatio[etabin]->Divide(theNoJetIDHist[etabin]);        
      theRatio[etabin]->Draw("L P E");
      
      //lines on ratio plot for y axis
      TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j]+theRatio[etabin]->GetBinWidth(1),1.0); 
      lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
      //TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j]+theRatio[etabin]->GetBinWidth(1),0.5); 
      //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
      //TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j]+theRatio[etabin]->GetBinWidth(1),1.5); 
      //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
      
      lineAtOne->Draw("same");
      //lineAtOneHalf->Draw("same");
      //lineAtOneEtOneHalf->Draw("same");
      
      
    }
    
    else if(var[j]=="Aj"||
	    var[j]=="xj"||
	    var[j]=="dphi"){ 
      
      jetpad1[etabin]->cd();
      theJetIDHist[etabin]->Draw("HIST E"); 
      theNoJetIDHist[etabin]->Draw("E SAME"); 
      
      float t1Loc1=0.15, t1Loc2=0.7; 	    	  
      if(var[j]=="Aj") t1Loc1=.7;
      
      
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),"p_{T}^{ld} > 60 GeV");    t2->Draw();	
      TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),"p_{T}^{subld} > 30GeV");    t3->Draw();	
      TLatex *t4=makeTLatex(t1Loc1, (t1Loc2-0.14), "#Delta#phi > #frac{2}{3}#pi");      t4->Draw();	
      //if(doJetIDPlots){
      //	TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.18),"Jet ID Applied");      t5->Draw();	}
      
      float legx1=0.15, legx2=legx1+0.11;
      if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
      float legy1=0.74, legy2=legy1+0.09;
      
      TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
      theJetQALeg->SetFillStyle(0);
      theJetQALeg->SetBorderSize(0.);
      std::string legdesc="";
      if(isData)legdesc="Data ";
      else      legdesc="PY8 MC ";
      theJetQALeg->AddEntry(theNoJetIDHist[etabin],(legdesc + "No Jet ID").c_str(),"lp");
      theJetQALeg->AddEntry(theJetIDHist[etabin],  (legdesc + "With Jet ID").c_str()  ,"lp");
      //theJetQALeg->AddEntry(theNoJetIDHist[etabin],"Data","lp");
      //theJetQALeg->AddEntry(theJetIDHist[etabin],  "P8 MC"  ,"lp");
      theJetQALeg->Draw(); 
      
      //Ratio Plots
      
      jetpad2[etabin]->cd();
      
      //TH1F* theRatio[etabin]=(TH1F*)theJetIDHist[etabin]->Clone("MCJetHistClone4Ratio");
      //I need to redefine the axis titles here
      std::string h_XAx_ratTitle=var_xAx_Titles[j], h_YAx_ratTitle="Data/MC"    ;
      std::string h_ratTitle   ="JetQA, Data/MC";//, "+var[j];
      //if(doJetIDPlots)h_ratTitle+=", w/ JetIDCut";      
      
      //theRatio[etabin]->SetMarkerStyle( kMultiply );
      //theRatio[etabin]->SetMarkerSize(0.99 );
      //theRatio[etabin]->SetMarkerColor( theMCColor );
      theRatio[etabin]->SetLineColor( theMCLineColor );
      //theRatio[etabin]->SetLineWidth( 2 );
      
      theRatio[etabin]->SetAxisRange(0.,2.,"Y");
      theRatio[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
      
      theRatio[etabin]->GetYaxis()->SetTitleSize(15);
      theRatio[etabin]->GetYaxis()->SetTitleFont(43);
      theRatio[etabin]->GetYaxis()->SetTitleOffset(2);
      theRatio[etabin]->GetYaxis()->SetLabelFont(43); 
      theRatio[etabin]->GetYaxis()->SetLabelSize(13);
      theRatio[etabin]->GetXaxis()->SetTitleSize(20);
      theRatio[etabin]->GetXaxis()->SetTitleFont(43);
      theRatio[etabin]->GetXaxis()->SetTitleOffset(4.);
      theRatio[etabin]->GetXaxis()->SetLabelFont(43); 
      theRatio[etabin]->GetXaxis()->SetLabelSize(15);
      
      theRatio[etabin]->SetTitle("");
      theRatio[etabin]->SetXTitle( h_XAx_Title.c_str() );
      theRatio[etabin]->SetYTitle(jetIDEfftitle.c_str());
      theRatio[etabin]->Divide(theNoJetIDHist[etabin]);
      theRatio[etabin]->Draw("E");
      
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
      jetpad1[etabin]->cd();
      theJetIDHist[etabin]->Draw("HIST E"); 
      theNoJetIDHist[etabin]->Draw("E SAME"); 
      
      float t1Loc1=0.7, t1Loc2=0.72; 	    	  
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t1->Draw();
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),"p_{T}^{ld} > 60 GeV");    t2->Draw();	
      TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),"p_{T}^{subld} > 30GeV");    t3->Draw();	
      TLatex *t4=makeTLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");	        t4->Draw();	
      //if(doJetIDPlots){
      //	TLatex *t5=makeTLatex((t1Loc1),(t1Loc2-.22),"Jet ID Applied");      t5->Draw();	}
      
      float legx1=0.70, legx2=legx1+0.11;
      float legy1=0.76, legy2=legy1+0.09;
      
      TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
      theJetQALeg->SetFillStyle(0);
      theJetQALeg->SetBorderSize(0.);
      std::string legdesc="";
      if(isData)legdesc="Data ";
      else      legdesc="PY8 MC ";
      theJetQALeg->AddEntry(theNoJetIDHist[etabin],(legdesc + "No Jet ID").c_str(),"lp");
      theJetQALeg->AddEntry(theJetIDHist[etabin],  (legdesc + "With Jet ID").c_str()  ,"lp");
      //theJetQALeg->AddEntry(theNoJetIDHist[etabin],"Data","lp");
      //theJetQALeg->AddEntry(theJetIDHist[etabin],  "P8 MC"  ,"lp");
      theJetQALeg->Draw(); 
      
      //Ratio Plots
      
      jetpad2[etabin]->cd();
      
      //TH1F* theRatio[etabin]=(TH1F*)theJetIDHist[etabin]->Clone("MCJetHistClone4Ratio");
      //theRatio[etabin]->SetMarkerStyle( kMultiply );
      //theRatio[etabin]->SetMarkerSize(0.99 );
      //theRatio[etabin]->SetMarkerColor( theMCColor );
      theRatio[etabin]->SetLineColor( theMCLineColor );
      //theRatio[etabin]->SetLineWidth( 2 );
      
      theRatio[etabin]->SetAxisRange(0.,2.,"Y");
      theRatio[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
      
      theRatio[etabin]->GetYaxis()->SetTitleSize(15);
      theRatio[etabin]->GetYaxis()->SetTitleFont(43);
      theRatio[etabin]->GetYaxis()->SetTitleOffset(2);
      theRatio[etabin]->GetYaxis()->SetLabelFont(43); 
      theRatio[etabin]->GetYaxis()->SetLabelSize(13);
      theRatio[etabin]->GetXaxis()->SetTitleSize(20);
      theRatio[etabin]->GetXaxis()->SetTitleFont(43);
      theRatio[etabin]->GetXaxis()->SetTitleOffset(4.);
      theRatio[etabin]->GetXaxis()->SetLabelFont(43); 
      theRatio[etabin]->GetXaxis()->SetLabelSize(15);
      
      theRatio[etabin]->SetTitle ("");
      theRatio[etabin]->SetXTitle( h_XAx_Title.c_str() );
      theRatio[etabin]->SetYTitle( jetIDEfftitle.c_str() );
      
      theRatio[etabin]->Divide(theNoJetIDHist[etabin]);
      theRatio[etabin]->Draw("HIST E");
      
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
    
  }//end etabin loop

  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() );
  
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    outcanv->SetTitle((inHistName+" Canvas").c_str());
    outcanv->Write((inHistName+"_canv").c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}




void printJetIDEtaHist( TFile* fin, int j, bool isData, int NetabinstoDraw,
			std::string thePDFFileName , std::string fullJetType, 
			long double theLumi  ,TFile* fout=NULL) {
  bool funcDebug=false;
  
  if(!fin){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  bool isConstitHist= (bool)( j>=jetConstits_varStart && j<dijet_varStart );  
  
  // STRINGS FOR HISTS
  TH1F* genjtptCut_h= NULL;
  std::string genjtptCut_str ="";
  std::string genJetCutString="";
  if(!isData){
    genjtptCut_h=(TH1F*)fin->Get( "hGenJetPtCut" );
    genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
    std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;  
    genJetCutString="p_{T}^{GEN} > "+genjtptCut_str+" GeV";
    std::cout<<"genJetCutString="<<genJetCutString<<std::endl;
  }
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
  std::string jtptCut_str=jtptQACut_str;
  std::string jetCutString="p_{T}^{RECO} > "+jtptCut_str+" GeV";
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  std::stringstream etaLo; etaLo.precision(1);
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  etaLo << std::fixed << jtetaLoCut_F;
  std::string jtetaLoCut_str = etaLo.str();
  std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;  
  
  float etaBinWidth=1.;
  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
  
  bool isybin=(var[j]=="jty");
  std::string inHistName="hJetQA_JetIDEff_"+var[j];
  if(isybin)
    inHistName+="_"+absetabins_str[0]+"y"+absetabins_str[NetabinstoDraw];
  else
    inHistName+="_"+absetabins_str[0]+"eta"+absetabins_str[NetabinstoDraw];
  
  TH1F* theNoJetIDHist[Netabins]={};
  TH1F* theJetIDHist[Netabins]={};
  TH1F* theRatio=NULL;
  TH1F* theNoJetIDEtaHist=NULL;
  TH1F* theJetIDEtaHist=NULL;
  
  for(int etabin=0; etabin<NetabinstoDraw;etabin++){
    
    std::string inNoJetIDHist="hJetQA_0wJetID_"+var[j];
    std::string inJetIDHist  ="hJetQA_1wJetID_"+var[j];
    if(isybin){
      inNoJetIDHist+="_ybin"+std::to_string(etabin);
      inJetIDHist  +="_ybin"+std::to_string(etabin);
    }
    else {
      inNoJetIDHist+="_etabin"+std::to_string(etabin);
      inJetIDHist  +="_etabin"+std::to_string(etabin);
    }
    std::cout<<" opening input Data "<<inNoJetIDHist<<std::endl<<std::endl;      
    
    theNoJetIDHist[etabin] = (TH1F*) ( (TH1*)fin->Get(inNoJetIDHist.c_str()) );
    theNoJetIDHist[etabin] = (TH1F*) (theNoJetIDHist[etabin]->Clone((inNoJetIDHist+"_"+std::to_string(j)+"clone").c_str()) );
    if(!theNoJetIDHist[etabin]){ std::cout<<"input data hist not found! skipping hist"<<std::endl;
      return;}
    
    theNoJetIDHist[etabin]->Scale( 1./theLumi);
    theNoJetIDHist[etabin]->Scale( 1./etaBinWidth);//eta bin width scaling
    
    std::cout<<" opening input MC "<<inJetIDHist<<std::endl<<std::endl;      
    theJetIDHist[etabin]= (TH1F*) ( (TH1*)fin->Get(inJetIDHist.c_str()) );  
    theJetIDHist[etabin]= (TH1F*) (theJetIDHist[etabin]->Clone((inJetIDHist+"_"+std::to_string(j)+"clone").c_str()) );
    if(!theJetIDHist[etabin]){ std::cout<<"input MC hist not found! skipping hist"<<std::endl;
      return;}
    theJetIDHist[etabin]->Scale(1./etaBinWidth);
    theJetIDHist[etabin]->Scale( 1./theLumi);  
    theJetIDHist[etabin]->ClearUnderflowAndOverflow();      theNoJetIDHist[etabin]->ClearUnderflowAndOverflow();  
    
    //theNoJetIDHist[etabin]->TH1::Rebin(var_xAx_reBin[j]);
    //theNoJetIDHist[etabin]->Scale( 1./theNoJetIDHist[etabin]->GetBinWidth(1) );          
    theNoJetIDHist[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
    //theNoJetIDHist[etabin]->GetXaxis()->SetRange(jetQAxmin[j],jetQAxmax[j]);
    
    //theJetIDHist[etabin]->TH1::Rebin(var_xAx_reBin[j]);
    //theJetIDHist[etabin]->Scale( 1./theJetIDHist[etabin]->GetBinWidth(1) );          
    theJetIDHist[etabin]->Scale( theNoJetIDHist[etabin]->Integral()/theJetIDHist[etabin]->Integral() );	        
    theJetIDHist[etabin]->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X"); 
    //theJetIDHist[etabin]->GetXaxis()->SetRange(jetQAxmin[j],jetQAxmax[j]); 

    // STYLE  
    //dataHistStyle(theNoJetIDHist[etabin]);
    //MCHistStyle(theJetIDHist[etabin]);
    nojetIDHistStyle(theNoJetIDHist[etabin]);
    jetIDHistStyle(theJetIDHist[etabin]);
    
    
    
    // TITLES
    std::string h_Title = PDStatsString_2; //for APS DNP
    std::string h_XAx_Title= var_xAx_Titles[j];
    //std::string h_YAx_Title= dcrossSectionAxTitle;
    //std::string h_YAx_Title= AUAxTitle;//for APS DNP
    std::string h_YAx_Title= ddcrossSectionAxTitle;//for APS DNP
    if((var[j]=="jty"||var[j]=="jteta" || var[j]=="jtphi" || isConstitHist) )h_YAx_Title=AUAxTitle;
    
    theNoJetIDHist[etabin]->SetTitle (    h_Title.c_str() );
    theNoJetIDHist[etabin]->SetXTitle( h_XAx_Title.c_str() );
    theNoJetIDHist[etabin]->SetYTitle( h_YAx_Title.c_str() );
    
    theJetIDHist[etabin]->SetTitle (    h_Title.c_str() );
    theJetIDHist[etabin]->SetXTitle( h_XAx_Title.c_str() );
    theJetIDHist[etabin]->SetYTitle( h_YAx_Title.c_str() );
    
    
  }

  std::string nojetidhistclonename="NoJetID";
  if(inHistName.find("etabin")!=std::string::npos)
    nojetidhistclonename+="EtaHist";
  else 
    nojetidhistclonename+="YHist";
  theNoJetIDEtaHist=(TH1F*)theNoJetIDHist[0]->Clone(nojetidhistclonename.c_str());
  theNoJetIDEtaHist->Reset("MICES");
  std::string jetidhistclonename="JetID";
  if(inHistName.find("etabin")!=std::string::npos)
    jetidhistclonename+="EtaHist";
  else 
    jetidhistclonename+="YHist";
  theJetIDEtaHist=(TH1F*)theJetIDHist[0]->Clone(jetidhistclonename.c_str());
  theJetIDEtaHist->Reset("MICES");
  for(int etabin=0; etabin<NetabinstoDraw;etabin++){
    theNoJetIDEtaHist->Add(theNoJetIDHist[etabin]);
    theJetIDEtaHist->Add(theJetIDHist[etabin]);
  }
  
  // TITLES
  std::string h_Title = PDStatsString_2; //for APS DNP
  std::string h_XAx_Title= var_xAx_Titles[j];
  //std::string h_YAx_Title= dcrossSectionAxTitle;
  //std::string h_YAx_Title= AUAxTitle;//for APS DNP
  std::string h_YAx_Title= ddcrossSectionAxTitle;//for APS DNP
  if((var[j]=="jtE"||var[j]=="jtm"||var[j]=="jteta" || var[j]=="jtphi" || isConstitHist) )h_YAx_Title=AUAxTitle;
  
  std::string jetEtaCutString="0.0 < #||{";
  if(inHistName.find("etabin")!=std::string::npos)
    jetEtaCutString+="#eta} < "+absetabins_str[NetabinstoDraw];
  else
    jetEtaCutString+="y} < "+absetabins_str[NetabinstoDraw];
  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;    

  // GET RATIO PLOT READY FOR BOTTOM PANEL
  std::string ratiojethistclonename="RatioJetQA";
  if(inHistName.find("etabin")!=std::string::npos)
    ratiojethistclonename+="EtaHist";
  else 
    ratiojethistclonename+="YHist";
  theRatio=(TH1F*)theJetIDEtaHist->Clone(ratiojethistclonename.c_str());  
  
  //general
  theRatio->SetTitle("");
  theRatio->SetLineColor( theMCLineColor );
  theRatio->SetFillColorAlpha(kRed, 0.);
  theRatio->SetFillStyle(0);

  //yaxis, range, title, labels
  theRatio->SetAxisRange(0.25,1.75,"Y");  
  theRatio->GetYaxis()->CenterTitle(true);
  theRatio->GetYaxis()->SetTitleSize(17);
  theRatio->GetYaxis()->SetTitleFont(43);
  theRatio->GetYaxis()->SetTitleOffset(1.4);
  theRatio->SetYTitle("JetID Eff." );    
  theRatio->GetYaxis()->SetLabelFont(43); 
  theRatio->GetYaxis()->SetLabelSize(13);  
  
  //xaxis, range, title, labels
  //  if((var[j]=="jteta" || var[j]=="jtphi" || isConstitHist) )
  //    theRatio->SetAxisRange(jetQAxmin[j],jetQAxmax[j],"X");	
  theRatio->GetXaxis()->SetTitleSize(20);
  theRatio->GetXaxis()->SetTitleFont(43);
  theRatio->GetXaxis()->SetTitleOffset(4.);
  theRatio->SetXTitle( h_XAx_Title.c_str() );  
  theRatio->GetXaxis()->SetLabelFont(43); 
  theRatio->GetXaxis()->SetLabelSize(18);
    
  
  
  
  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", (int)JETQACANVX, (int)JETQACANVY); //This is for the Jet QA Evt plots without logs
  TPad* jetpad1=NULL;
  TPad* jetpad2=NULL;

  std::cout<<"JETQACANVX="<<JETQACANVX<<std::endl;
  std::cout<<"JETQACANVY="<<JETQACANVY<<std::endl;
  
  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
  temp_canvJet->cd();
  
  
  jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  //std::cout<<"DEFAULT MARGINS ARE...."<<std::endl;//always 0.1
  //std::cout<<"LEFT  : "<<jetpad1->GetLeftMargin()<<std::endl;
  //std::cout<<"RIGHT : "<<jetpad1->GetRightMargin()<<std::endl;
  //std::cout<<"TOP   : "<<jetpad1->GetTopMargin()<<std::endl;
  //std::cout<<"BOTTOM: "<<jetpad1->GetBottomMargin()<<std::endl;
  jetpad1->SetGridx(0);
  jetpad1->SetGridy(0);  
  //jetpad1->SetLeftMargin(0.09);
  //jetpad1->SetRightMargin(0.04);
  //jetpad1->SetTopMargin(0.07);//good idea?
  jetpad1->SetBottomMargin(0);//do not change
  jetpad1->SetLogx(0);  
  jetpad1->SetLogy(1);  
  jetpad1->Draw(); 
  
  jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
  jetpad2->SetGridx(1);
  jetpad2->SetGridy(0);  
  //jetpad2->SetLeftMargin(0.09);
  //jetpad2->SetRightMargin(0.04);
  jetpad2->SetTopMargin(0);//do not change
  jetpad2->SetBottomMargin(0.3);//maybe change
  
  jetpad2->SetLogx(0);  
  jetpad2->SetLogy(0);  
  jetpad2->Draw();
  
  jetpad1->cd();
  
  theJetIDEtaHist->Draw("HIST E"); 
  theNoJetIDEtaHist->Draw("E SAME");  
  
  float t1Loc1=0.54, t1Loc2=0.82;
  //float t1Loc1=0.50, t1Loc2=0.72; 	    
  TLatex* t1= makeTLatex(  t1Loc1  , t1Loc2 ,  fullJetType+"Jets" );    t1->Draw(); t1Loc2-=0.05;
  TLatex *t2= makeTLatex( (t1Loc1) , t1Loc2 , (jetCutString)      );    t2->Draw(); t1Loc2-=0.05;
  if(!isData){
    TLatex *t3= makeTLatex( (t1Loc1) , t1Loc2 , (genJetCutString)   );    t3->Draw(); t1Loc2-=0.05; }
  TLatex *t4= makeTLatex( (t1Loc1) , t1Loc2 , (jetEtaCutString)   );    t4->Draw();	    

  
  float legx1=0.8, legx2=0.89;
  float legy1=0.8, legy2=0.89;
  
  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
  theJetQALeg->SetFillStyle(0);
  theJetQALeg->SetBorderSize(0.);
  std::string legdesc="";
  if(isData)legdesc="Data ";
  else      legdesc="PY8 MC ";
  theJetQALeg->AddEntry(theNoJetIDEtaHist,(legdesc + "No Jet ID").c_str(),"lp");
  theJetQALeg->AddEntry(theJetIDEtaHist  ,  (legdesc + "With Jet ID").c_str()  ,"lp");
  theJetQALeg->Draw();
  
  //Ratio Plots    
  jetpad2->cd();
  theRatio->Divide(theNoJetIDEtaHist);        
  theRatio->Draw("HIST E");
  
  //lines on ratio plot for y axis
  TLine* lineAtOne          = new TLine(jetQAxmin[j],1.0,jetQAxmax[j]+theRatio->GetBinWidth(1),1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
  //TLine* lineAtOneHalf      = new TLine(jetQAxmin[j],0.5,jetQAxmax[j]+theRatio[etabin]->GetBinWidth(1),0.5); 
  //lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
  //TLine* lineAtOneEtOneHalf = new TLine(jetQAxmin[j],1.5,jetQAxmax[j]+theRatio[etabin]->GetBinWidth(1),1.5); 
  //lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
  
  lineAtOne->Draw("same");
  //lineAtOneHalf->Draw("same");
  //lineAtOneEtOneHalf->Draw("same");
  
  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() );
  
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    outcanv->SetTitle((inHistName+" Canvas").c_str());
    outcanv->Write((inHistName+"_canv").c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}

////void printJetIDHist( TFile* fin_jetID , TFile* fin_nojetID, 
//void printJetIDHist( TFile* fin , 
//		     int j, bool isData,
//		     std::string thePDFFileName , std::string fullJetType, 
//		     long double theLumi,  
//		     TFile* fout=NULL) {
//  
//  
//  if(!fin ){    std::cout<<"input file not found, cannot look at jetID selection"<<std::endl; 
//    return; }
//  bool funcDebug=false;
//    
//  // PUT JET CUTS SECTION IN HERE
//  TH1F* genjtptCut_h=NULL; 
//  std::string genjtptCut_str;
//  if(!isData){
//    genjtptCut_h=(TH1F*)fin->Get( "hGenJePttCut" );
//    genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
//    std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;}
//  
//  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
//  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
//  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
//  
//  std::string jetCutString;
//  if(!isData)jetCutString="p_{T}^{gen}>"+genjtptCut_str+" GeV, "+"p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
//  else jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
//  std::cout<<"jetCutString="<<jetCutString<<std::endl;
//  
//  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
//  std::string jtptQACut_str = std::to_string( (int) jtptQACut_h->GetMean() );
//  std::cout<<"jtptQACut_str = "<<jtptQACut_str<<std::endl;
//  
//  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
//  std::stringstream etaLo; etaLo.precision(1);
//  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
//  etaLo << std::fixed << jtetaLoCut_F;
//  std::string jtetaLoCut_str = etaLo.str(); 
// std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;
// 
// TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
//  std::stringstream etaHi; etaHi.precision(1);
//  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
//  etaHi << std::fixed << jtetaHiCut_F;
//  std::string jtetaHiCut_str = etaHi.str();
//  std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;
//  
//  std::string jetEtaCutString=jtetaLoCut_str+"<|y|<"+jtetaHiCut_str;
//  std::cout<<"jetEtaCutString="<<jetEtaCutString<<std::endl;
//  
//  float etaBinWidth=2.*(jtetaHiCut_F-jtetaLoCut_F);//factor of two, because abseta
//  std::cout<<"etaBinWidth="<<etaBinWidth<<std::endl;
//  
//  std::string sampleDescString;
//  if(isData)sampleDescString=sqrts2k15ppString;
//  else      sampleDescString=Py8TuneString;
//  
//  std::string inHistName="hJetQA_0wJetID_"+var[j];
//  if(funcDebug)std::cout<<" opening input hist "<<inHistName<<std::endl<<std::endl;      
//  
//  TH1F* theNonJetIDHist=(TH1F*) ( (TH1*)fin->Get(inHistName.c_str()) );
//  if(!theNonJetIDHist){ std::cout<<"input hist not found! skipping hist"<<std::endl;
//    return;}
//  
//  std::string inJetIDHistName="hJetQA_1wJetID_"+var[j];
//  if(funcDebug)std::cout<<" opening input jetID hist "<<inJetIDHistName<<std::endl<<std::endl;      
//  TH1F* theJetIDHist= (TH1F*) ( (TH1*)fin->Get(inJetIDHistName.c_str()) );  
//  if(!theJetIDHist){ std::cout<<"input jetID hist not found! skipping hist"<<std::endl;
//    return;}
//  
//  // MORE SCALING
//  if(var[j]=="jtpt" || var[j]=="rawpt" ||
//     var[j]=="leadJetPt" || var[j]=="subleadJetPt"      ||
//     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
//    
//    theNonJetIDHist=(TH1F*)theNonJetIDHist->TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_jtPt", ptbins_debug);
//    divideBinWidth(theNonJetIDHist);      
//    if(isData)theNonJetIDHist->Scale(1./theLumi);
//    theNonJetIDHist->Scale(1./etaBinWidth);
//    
//    theJetIDHist=(TH1F*)theJetIDHist->TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_jtPt", ptbins_debug);
//    divideBinWidth(theJetIDHist);      
//    if(isData)theJetIDHist->Scale(1./theLumi);
//    theJetIDHist->Scale(1./etaBinWidth);
//  }
//  
//  else {
//    
//    theNonJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
//    theNonJetIDHist->Scale( 1./theNonJetIDHist->GetBinWidth(1) );          
//    if(isData)theNonJetIDHist->Scale(1./theLumi);
//    theNonJetIDHist->Scale(1./etaBinWidth);
//    
//    theJetIDHist->TH1::Rebin(var_xAx_reBin[j]);
//    theJetIDHist->Scale( 1./theJetIDHist->GetBinWidth(1) );         
//    if(isData)theJetIDHist->Scale(1./theLumi);
//    theJetIDHist->Scale(1./etaBinWidth);
//    //    theJetIDHist->Scale( theNonJetIDHist->Integral()/theJetIDHist->Integral() );	        
//    
//  }
//  
//  
//  
//  
//  // TITLES
//  std::string h_XAx_Title=var_xAx_Titles[j], h_YAx_Title=AUAxTitle;
//  std::string h_Title="HIN Jet ID Efficiency";//, "+var[j];
//  
//  theNonJetIDHist->SetTitle (    h_Title.c_str() );
//  theNonJetIDHist->SetXTitle( h_XAx_Title.c_str() );
//  theNonJetIDHist->SetYTitle( h_YAx_Title.c_str() );
//  
//  theJetIDHist->SetTitle (    h_Title.c_str() );
//  theJetIDHist->SetXTitle( h_XAx_Title.c_str() );
//  theJetIDHist->SetYTitle( h_YAx_Title.c_str() );  
//  
//  
//  // STYLE  
//  TH1F* theRatio=(TH1F*)theJetIDHist->Clone("MCJetHistClone4Ratio");
//  
//  dataHistStyle(theNonJetIDHist);
//  MCHistStyle(theJetIDHist);
//  //ratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
//  jetIDratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
//  
//  if(funcDebug)theNonJetIDHist->Print("base");       std::cout<<std::endl;  
//  if(funcDebug)theJetIDHist->Print("base");       std::cout<<std::endl;
//  
//  if(funcDebug) std::cout<<"creating temporary canvas for printing Jet plots..."<<std::endl;
//  TCanvas *temp_canvJet = new TCanvas("tempEvt", "temp Jet Canv withLog", CANVX, CANVY); //This is for the Jet QA Evt plots without logs
//  temp_canvJet->cd();
//  
//  TPad *jetpad1 = new TPad("jetpad1", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
//  jetpad1->SetLogy(1); 
//  jetpad1->SetGridx(1);
//  jetpad1->SetGridy(1);
//  jetpad1->SetBottomMargin(0);
//  
//    
//  TPad *jetpad2 = new TPad("jetpad2", "Ratio Pad"  , 0.0, 0.05, 1.0, 0.3);
//  jetpad2->SetTopMargin(0);
//  jetpad2->SetBottomMargin(0.3);
//  //jetpad2->SetLogy(0);
//  jetpad2->SetLogy(1);
//  
//  
//  if(var[j]=="jtpt"      || var[j]=="rawpt"          ||
//     var[j]=="leadJetPt" || var[j]=="subleadJetPt"   ||
//     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
//    
//    jetpad1->SetLogx(1);  
//    jetpad2->SetLogx(1);  
//    
//    theJetIDHist->GetXaxis()->SetNoExponent(true);
//    theJetIDHist->GetXaxis()->SetMoreLogLabels(true);
//    
//    theNonJetIDHist->GetXaxis()->SetNoExponent(true);
//    theNonJetIDHist->GetXaxis()->SetMoreLogLabels(true);
//
//    theRatio->GetXaxis()->SetNoExponent(true);
//    theRatio->GetXaxis()->SetMoreLogLabels(true);
//  }
// 
//  jetpad1->Draw();
//  jetpad2->Draw();    
//
//
//  // variable specific draw orders, canv, etc.
//  if(var[j]=="jtpt"||var[j]=="rawpt" ||
//     var[j]=="jtpt_forRes" ||     var[j]=="jtpt_L2Res" ||     var[j]=="jtpt_L3Res"      ) {
//    
//    jetpad1->cd();    
//    theJetIDHist->Draw("HIST E"); 
//    theNonJetIDHist->Draw("E SAME");  
//    
//    float t1Loc1=0.54, t1Loc2=0.82;
//    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
//    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	    
//    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	    
//    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
//    
//    float legx1=0.73, legx2=legx1+0.11;
//    float legy1=0.76, legy2=legy1+0.09;
//    
//    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
//    theJetQALeg->SetFillStyle(0);
//    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
//    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
//    theJetQALeg->Draw(); 
//    
//    jetpad2->cd();
//    
//    theRatio->Divide(theNonJetIDHist); //jetid/nojetid
//    fracSubtracted(theRatio); //1-jetid/nojetid
//
//
//    theRatio->Draw("HIST E");
//    
//    
//
//    
//  }
//  
//  else if(var[j]=="jteta" ||
//	  var[j]=="jtphi"	||
//	  ( j>=jetConstits_varStart && j<dijet_varStart )){ 
//    
//    
//    jetpad1->cd();
//    theJetIDHist->Draw("HIST E"); 
//    theNonJetIDHist->Draw("E same");  
//    
//    float t1Loc1=0.54, t1Loc2=0.82;
//    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();
//    TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	
//    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	
//    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
//    
//    float legx1=0.73, legx2=legx1+0.11;
//    float legy1=0.76, legy2=legy1+0.09;
//    
//    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
//    theJetQALeg->SetFillStyle(0);
//    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
//    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
//    
//    theJetQALeg->Draw();
//    
//    jetpad2->cd();
//    
//    
//    theRatio->Divide(theNonJetIDHist);
//    fracSubtracted(theRatio); //1-jetid/nojetid
//    theRatio->Draw("HIST E");
//    
//    
//  }
//  
//  else if(var[j]=="Aj"||
//	  var[j]=="xj"||
//	  var[j]=="dphi"){ 
//    
//    jetpad1->cd();
//    theJetIDHist->Draw("HIST E"); 
//    theNonJetIDHist->Draw("E SAME"); 
//    
//    float t1Loc1=0.15, t1Loc2=0.7; 	    	  
//    if(var[j]=="Aj") t1Loc1=.7;    
//    
//    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
//    TLatex *t2=makeTLatex(t1Loc1,t1Loc2-.04,(fullJetType+"Jets").c_str());    t2->Draw();    
//    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),"p_{T}^{ld} > 60 GeV");    t3->Draw();		  
//    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.12),"p_{T}^{subld} > 30GeV");    t4->Draw();	    
//    TLatex *t5=makeTLatex(t1Loc1, (t1Loc2-0.16), "#Delta#phi > #frac{2}{3}#pi");	        t5->Draw();	
//
//    float legx1=0.15, legx2=legx1+0.11;
//    if(var[j]=="Aj") {legx1=.7; legx2=legx1+0.11;}
//    float legy1=0.74, legy2=legy1+0.09;
//    
//    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
//    theJetQALeg->SetFillStyle(0);
//    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
//    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
//    theJetQALeg->Draw(); 
//    
//    //Ratio Plots
//    
//    jetpad2->cd();
//
//    
//    theRatio->Divide(theNonJetIDHist);
//    fracSubtracted(theRatio); //1-jetid/nojetid
//    theRatio->Draw("E");
//    
//    
//  }
//      
//  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt"){ 
//    
//    jetpad1->cd();    
//    theJetIDHist->Draw("HIST E"); 
//    theNonJetIDHist->Draw("E SAME"); 
//    
//    float t1Loc1=0.7, t1Loc2=0.72; 	    	  
//
//    TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();
//    TLatex *t2=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());    t2->Draw();
//    TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.05),"p_{T}^{ld} > 60 GeV");    t3->Draw();	
//    TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.10),"p_{T}^{subld} > 30GeV");    t4->Draw();	
//    TLatex *t5=makeTLatex(t1Loc1, (t1Loc2-0.17), "#Delta#phi > #frac{2}{3}#pi");    t5->Draw();	
//
//    float legx1=0.73, legx2=legx1+0.11;
//    float legy1=0.76, legy2=legy1+0.09;
//    
//    TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
//    theJetQALeg->SetFillStyle(0);
//    theJetQALeg->AddEntry(theNonJetIDHist, "no JetID",   "lp");
//    theJetQALeg->AddEntry(theJetIDHist,    "with JetID", "lp");
//    theJetQALeg->Draw(); 
//    
//    //Ratio Plots
//    
//    jetpad2->cd();
//    
//    theRatio->Divide(theNonJetIDHist);
//    fracSubtracted(theRatio); //1-jetid/nojetid
//    theRatio->Draw("HIST E");
//    
//  }
//  
//  
//  //TLines for ratio plots
//  float lowLim, highLim;
//  lowLim=theRatio->GetBinLowEdge(1);
//  highLim=theRatio->GetBinLowEdge(theRatio->GetNbinsX())+theRatio->GetBinWidth(theRatio->GetNbinsX());
//  
//  TLine* lineAtOne          = new TLine(lowLim,1.0,highLim,1.0); 
//  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
//  TLine* lineAtOneHalf      = new TLine(lowLim,0.5,highLim,0.5); 
//  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
//  TLine* lineAtOneEtOneHalf = new TLine(lowLim,1.5,highLim,1.5); 
//  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);
//  
//  //lineAtOne->Draw("same");
//  //lineAtOneHalf->Draw("same");
//  //lineAtOneEtOneHalf->Draw("same");
//
//  
//  // print to PDF file
//  temp_canvJet->Print( thePDFFileName.c_str() ); 
//  if(fout){
//    fout->cd();
//    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
//    
//    std::string outcanvtitle="JetID Eff. "+var[j]+" Canvas";
//    if(isData)outcanvtitle="Data "+outcanvtitle;
//    else outcanvtitle="MC "+outcanvtitle;    
//    outcanv->SetTitle(outcanvtitle.c_str());
//    
//    std::string outcanvname="JetIDEff_"+var[j];
//    if(isData)outcanvname="Data_"+outcanvname;
//    else outcanvname="MC_"+outcanvname;
//    //outcanvname+="_etabin"+std::to_string(etabin);
//    outcanv->Write(outcanvname.c_str());    
//  }
//  
//  temp_canvJet->Close();  
//  return;
//
//}







void printTupelJetIDHist( TFile* fin_jetID , int etabin, bool print_incjetana_tupelequiv,
			  int j, bool isData,
			  std::string thePDFFileName , std::string fullJetType, 
			  long double theLumi,  
			  TFile* fout=NULL) {
  
  
  if(!fin_jetID ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;

  //if(!isData)assert(false);

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
  jetCutString="p_{T}>"+jtptQACut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  
  
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
  else if(etabin==1) jetEtaCutString = "0.5 < #||{y} < 1.0";
  else if(etabin==2) jetEtaCutString = "1.0 < #||{y} < 1.5";
  else if(etabin==3) jetEtaCutString = "1.5 < #||{y} < 2.0";
  else if(etabin==4) jetEtaCutString = "2.0 < #||{y} < 2.5";
  else if(etabin==5) jetEtaCutString = "2.5 < #||{y} < 3.0";  
  else if(etabin==6) jetEtaCutString = "3.2 < #||{y} < 4.7";  
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
  std::string h_XAx_Title;
  if(!print_incjetana_tupelequiv)h_XAx_Title=tupelvar_xAx_Titles[j];//, h_YAx_Title=AUAxTitle;
  else h_XAx_Title=tupelvar_incjetanaequiv_xAx_Titles[j];

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
    outcanvname+="_etabin"+std::to_string(etabin);
    outcanv->Write(outcanvname.c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}














void printHINvSMPJetIDHist_pt( TFile* fin_HINjetID=NULL , TFile* fin_SMPjetID=NULL , 
			       int etabin=0, std::string jetptvar="jtpt",
			       bool isData=true,
			       std::string thePDFFileName="" , std::string fullJetType="", 
			       TFile* fout=NULL) {
  
  
  if(fin_HINjetID->IsZombie() ){    std::cout<<"HIN jetID input file not found, cannot look at event counts"<<std::endl; 
    return; }
  if(fin_SMPjetID->IsZombie() ){    std::cout<<"SMP jetID input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;

  //if(!isData)assert(false);

  //// PUT JET CUTS SECTION IN HERE
  //TH1F* genjtptCut_h=NULL; 
  //std::string genjtptCut_str;
  //if(!isData){
  //  genjtptCut_h=(TH1F*)fin_jetID->Get( "hGenJetPtCut" );
  //  genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  //  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;}
  
  TH1F* jtptCut_h= (TH1F*)fin_HINjetID->Get( "hJetQAPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  
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
  
  //hJetQA_1wJetID_jtpt_etabin0
  std::string HINinHistName="hJetQA_1wJetID_";
  if(jetptvar=="jtpt") HINinHistName+="jtpt_";
  else if(jetptvar=="rawpt") HINinHistName+="rawpt_";
  HINinHistName+="etabin"+std::to_string(etabin)+"";
  if(funcDebug)std::cout<<" opening input HIN Jet ID hist "<<HINinHistName<<std::endl<<std::endl;      
  TH1F* theHINJetIDHist=(TH1F*) ( (TH1*)fin_HINjetID->Get(HINinHistName.c_str()) );
  theHINJetIDHist=(TH1F*)theHINJetIDHist->TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_jtPt", ptbins_debug);
  divideBinWidth(theHINJetIDHist);      
  theHINJetIDHist->Scale(1./etaBinWidth);
  
  //hTupelJetQA_1jetID_etabin0_Pt
  std::string SMPinHistName="hTupelJetQA_1jetID_etabin"+std::to_string(etabin)+"_";
  if(jetptvar=="jtpt")SMPinHistName+="Pt";
  else if(jetptvar=="rawpt")SMPinHistName+="RawPt";
  if(funcDebug)std::cout<<" opening input SMP Jet ID hist "<<SMPinHistName<<std::endl<<std::endl;      
  TH1F* theSMPJetIDHist= (TH1F*) ( (TH1*)fin_SMPjetID->Get(SMPinHistName.c_str()) );  
  theSMPJetIDHist=(TH1F*)theSMPJetIDHist->TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_jtPt", ptbins_debug);
  divideBinWidth(theSMPJetIDHist);      
  theSMPJetIDHist->Scale(1./etaBinWidth);
  
  
 
  
  // TITLES
  std::string h_XAx_Title="Jet p_{T}";
  if(jetptvar=="jtpt")h_XAx_Title="RECO "+h_XAx_Title;
  else if(jetptvar=="rawpt")h_XAx_Title="RAW "+h_XAx_Title;
  else if(jetptvar=="genpt")h_XAx_Title="GEN "+h_XAx_Title;
  
  std::string  h_YAx_Title=ddcrossSectionAxTitle;//for APS DNP
  std::string h_Title="SMP v. HIN Jet ID";
  
  theHINJetIDHist->SetTitle (    h_Title.c_str() );
  theHINJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theHINJetIDHist->SetYTitle( h_YAx_Title.c_str() );
  
  theSMPJetIDHist->SetTitle (    h_Title.c_str() );
  theSMPJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theSMPJetIDHist->SetYTitle( h_YAx_Title.c_str() );  
  
  
  // STYLE  
  TH1F* theRatio=(TH1F*)theSMPJetIDHist->Clone("MCJetHistClone4Ratio");
  
  dataHistStyle(theHINJetIDHist);
  MCHistStyle(theSMPJetIDHist);
  //ratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  jetIDratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  
  if(funcDebug)theHINJetIDHist->Print("base");       std::cout<<std::endl;  
  if(funcDebug)theSMPJetIDHist->Print("base");       std::cout<<std::endl;
  
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
  jetpad2->SetGridx(1);
  //jetpad2->SetGridy(1);
  //jetpad2->SetLogy(0);
  jetpad2->SetLogy(0);  
  
  jetpad1->SetLogx(1);  
  jetpad2->SetLogx(1);  
  
  theSMPJetIDHist->GetXaxis()->SetNoExponent(true);
  theSMPJetIDHist->GetXaxis()->SetMoreLogLabels(true);
  
  theHINJetIDHist->GetXaxis()->SetNoExponent(true);
  theHINJetIDHist->GetXaxis()->SetMoreLogLabels(true);
  
  theRatio->GetYaxis()->SetTitle("SMP/HIN Jet ID Ratio");
  theRatio->SetMaximum(1.05);
  theRatio->SetMinimum(0.85);
  theRatio->GetXaxis()->SetNoExponent(true);
  theRatio->GetXaxis()->SetMoreLogLabels(true);
  
  jetpad1->Draw();
  jetpad2->Draw();    
  
  
  
  jetpad1->cd();    
  theSMPJetIDHist->Draw("HIST E"); 
  theHINJetIDHist->Draw("E SAME");  
  
  float t1Loc1=0.54, t1Loc2=0.82;
  TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
  TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	    
  TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	    
  TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
  
  float legx1=0.73, legx2=legx1+0.11;
  float legy1=0.76, legy2=legy1+0.09;
  
  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
  theJetQALeg->SetFillStyle(0);
  theJetQALeg->AddEntry(theHINJetIDHist, "HIN JetID",   "lp");
  theJetQALeg->AddEntry(theSMPJetIDHist,    "SMP JetID", "lp");
  theJetQALeg->Draw(); 
  
  jetpad2->cd();
  
  theRatio->Divide(theHINJetIDHist); //jetid/nojetid
  //fracSubtracted(theRatio); //1-jetid/nojetid
  


  float lowLim, highLim;
  lowLim=theRatio->GetBinLowEdge(1);
  highLim=theRatio->GetBinLowEdge(theRatio->GetNbinsX())+theRatio->GetBinWidth(theRatio->GetNbinsX());
  
  TLine* lineAtOne          = new TLine(lowLim,1.0,highLim,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
//  TLine* lineAtOneHalf      = new TLine(lowLim,0.5,highLim,0.5); 
//  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
//  TLine* lineAtOneEtOneHalf = new TLine(lowLim,1.5,highLim,1.5); 
//  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);


  
  theRatio->Draw("HIST");
  lineAtOne->Draw();
  
  
  

  
  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() ); 
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    
    std::string outcanvtitle="HIN v SMP JetID";
    if (jetptvar=="jtpt")outcanvtitle+=", RECO Jet p_{T}";
    else if (jetptvar=="rawpt")outcanvtitle+=", RAW Jet p_{T}";
    if(isData)outcanvtitle="Data "+outcanvtitle;
    else outcanvtitle="MC "+outcanvtitle;    
    outcanv->SetTitle(outcanvtitle.c_str());
    
    std::string outcanvname="HINvSMPJetID_";
    if (jetptvar=="jtpt")outcanvname+="jtpt";
    else if (jetptvar=="rawpt")outcanvname+="rawpt";
    if(isData)outcanvname="Data_"+outcanvname;
    else outcanvname="MC_"+outcanvname;
    outcanvname+="_etabin"+std::to_string(etabin);
    outcanv->Write(outcanvname.c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}









void print_tightVlooseJetIDHist_pt( TFile* fin_loose_jetID=NULL , TFile* fin_tight_jetID=NULL , 
				    int etabin=0, std::string jetptvar="jtpt",
				    bool isData=true,
				    std::string thePDFFileName="" , std::string fullJetType="", 
				    TFile* fout=NULL) {
  
  
  if(fin_loose_jetID->IsZombie() ){    std::cout<<"HIN jetID input file not found, cannot look at event counts"<<std::endl; 
    return; }
  if(fin_tight_jetID->IsZombie() ){    std::cout<<"SMP jetID input file not found, cannot look at event counts"<<std::endl; 
    return; }
  bool funcDebug=false;

  //if(!isData)assert(false);

  //// PUT JET CUTS SECTION IN HERE
  //TH1F* genjtptCut_h=NULL; 
  //std::string genjtptCut_str;
  //if(!isData){
  //  genjtptCut_h=(TH1F*)fin_jetID->Get( "hGenJetPtCut" );
  //  genjtptCut_str = std::to_string( (int) genjtptCut_h->GetMean() );
  //  std::cout<<"genjtptCut_str = "<<genjtptCut_str<<std::endl;}
  
  TH1F* jtptCut_h= (TH1F*)fin_loose_jetID->Get( "hJetQAPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;
  
  std::string jetCutString="p_{T}>"+jtptCut_str+" GeV";//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  std::cout<<"jetCutString="<<jetCutString<<std::endl;
  
  
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
  
  //hJetQA_1wJetID_jtpt_etabin0
  std::string loose_inHistName="hJetQA_1wJetID_";
  if(jetptvar=="jtpt") loose_inHistName+="jtpt_";
  else if(jetptvar=="rawpt") loose_inHistName+="rawpt_";
  loose_inHistName+="etabin"+std::to_string(etabin)+"";
  if(funcDebug)std::cout<<" opening input HIN Jet ID hist "<<loose_inHistName<<std::endl<<std::endl;      
  TH1F* theLooseJetIDHist=(TH1F*) ( (TH1*)fin_loose_jetID->Get(loose_inHistName.c_str()) );
  theLooseJetIDHist=(TH1F*)theLooseJetIDHist->TH1::Rebin(nbins_pt_debug, "ppData_analysisBins_jtPt", ptbins_debug);
  divideBinWidth(theLooseJetIDHist);      
  theLooseJetIDHist->Scale(1./etaBinWidth);
  
  //hTupelJetQA_1jetID_etabin0_Pt
  std::string tight_inHistName=loose_inHistName;
  //std::string tight_inHistName="hTupelJetQA_1jetID_etabin"+std::to_string(etabin)+"_";
  //if(jetptvar=="jtpt")tight_inHistName+="Pt";
  //else if(jetptvar=="rawpt")tight_inHistName+="RawPt";
  //if(funcDebug)std::cout<<" opening input SMP Jet ID hist "<<tight_inHistName<<std::endl<<std::endl;      
  TH1F* theTightJetIDHist= (TH1F*) ( (TH1*)fin_tight_jetID->Get(tight_inHistName.c_str()) );  
  theTightJetIDHist=(TH1F*)theTightJetIDHist->TH1::Rebin(nbins_pt_debug, "ppMC_analysisBins_jtPt", ptbins_debug);
  divideBinWidth(theTightJetIDHist);      
  theTightJetIDHist->Scale(1./etaBinWidth);

 
  
  // TITLES
  std::string h_XAx_Title="Jet p_{T}";
  if(jetptvar=="jtpt")h_XAx_Title="RECO "+h_XAx_Title;
  else if(jetptvar=="rawpt")h_XAx_Title="RAW "+h_XAx_Title;
  else if(jetptvar=="genpt")h_XAx_Title="GEN "+h_XAx_Title;
  
  std::string  h_YAx_Title=ddcrossSectionAxTitle;//for APS DNP
  std::string h_Title="loose v. togjt Jet ID";
  
  theLooseJetIDHist->SetTitle (    h_Title.c_str() );
  theLooseJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theLooseJetIDHist->SetYTitle( h_YAx_Title.c_str() );
  
  theTightJetIDHist->SetTitle (    h_Title.c_str() );
  theTightJetIDHist->SetXTitle( h_XAx_Title.c_str() );
  theTightJetIDHist->SetYTitle( h_YAx_Title.c_str() );  
  
  
  // STYLE  
  TH1F* theRatio=(TH1F*)theTightJetIDHist->Clone("MCJetHistClone4Ratio");
  
  dataHistStyle(theLooseJetIDHist);
  MCHistStyle(theTightJetIDHist);
  //ratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  jetIDratioHistStyle(theRatio, h_XAx_Title, jetIDratioTitle);
  
  if(funcDebug)theLooseJetIDHist->Print("base");       std::cout<<std::endl;  
  if(funcDebug)theTightJetIDHist->Print("base");       std::cout<<std::endl;
  
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
  jetpad2->SetGridx(1);
  //jetpad2->SetGridy(1);
  //jetpad2->SetLogy(0);
  jetpad2->SetLogy(0);  
  
  jetpad1->SetLogx(1);  
  jetpad2->SetLogx(1);  
  
  theTightJetIDHist->GetXaxis()->SetNoExponent(true);
  theTightJetIDHist->GetXaxis()->SetMoreLogLabels(true);
  
  theLooseJetIDHist->GetXaxis()->SetNoExponent(true);
  theLooseJetIDHist->GetXaxis()->SetMoreLogLabels(true);
  
  theRatio->GetYaxis()->SetTitle("tight/loose Jet ID Ratio");
  theRatio->SetMaximum(1.05);
  theRatio->SetMinimum(0.85);
  theRatio->GetXaxis()->SetNoExponent(true);
  theRatio->GetXaxis()->SetMoreLogLabels(true);
  
  jetpad1->Draw();
  jetpad2->Draw();    
  
  
  
  jetpad1->cd();    
  theTightJetIDHist->Draw("HIST E"); 
  theLooseJetIDHist->Draw("E SAME");  
  
  float t1Loc1=0.54, t1Loc2=0.82;
  TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(sampleDescString).c_str());    t1->Draw();    
  TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.05),(fullJetType+"Jets").c_str());    t2->Draw();	    
  TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.10),(jetCutString).c_str());    t3->Draw();	    
  TLatex *t4=makeTLatex((t1Loc1),(t1Loc2-.15),(jetEtaCutString).c_str());    t4->Draw();	
  
  float legx1=0.73, legx2=legx1+0.11;
  float legy1=0.76, legy2=legy1+0.09;
  
  TLegend* theJetQALeg=new TLegend(legx1,legy1,legx2,legy2, NULL,"brNDC");	  
  theJetQALeg->SetFillStyle(0);
  theJetQALeg->AddEntry(theLooseJetIDHist, "loose JetID",   "lp");
  theJetQALeg->AddEntry(theTightJetIDHist,    "tight JetID", "lp");
  theJetQALeg->Draw(); 
  
  jetpad2->cd();
  
  theRatio->Divide(theLooseJetIDHist); //jetid/nojetid
  //fracSubtracted(theRatio); //1-jetid/nojetid
  


  float lowLim, highLim;
  lowLim=theRatio->GetBinLowEdge(1);
  highLim=theRatio->GetBinLowEdge(theRatio->GetNbinsX())+theRatio->GetBinWidth(theRatio->GetNbinsX());
  
  TLine* lineAtOne          = new TLine(lowLim,1.0,highLim,1.0); 
  lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
//  TLine* lineAtOneHalf      = new TLine(lowLim,0.5,highLim,0.5); 
//  lineAtOneHalf->SetLineColor(12);      lineAtOneHalf->SetLineStyle(2);
//  TLine* lineAtOneEtOneHalf = new TLine(lowLim,1.5,highLim,1.5); 
//  lineAtOneEtOneHalf->SetLineColor(12); lineAtOneEtOneHalf->SetLineStyle(2);


  
  theRatio->Draw("HIST");
  lineAtOne->Draw();
  
  
  

  
  
  // print to PDF file
  temp_canvJet->Print( thePDFFileName.c_str() ); 
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJet->DrawClone();
    
    std::string outcanvtitle="loose v tight JetID";
    if (jetptvar=="jtpt")outcanvtitle+=", RECO Jet p_{T}";
    else if (jetptvar=="rawpt")outcanvtitle+=", RAW Jet p_{T}";
    if(isData)outcanvtitle="Data "+outcanvtitle;
    else outcanvtitle="MC "+outcanvtitle;    
    outcanv->SetTitle(outcanvtitle.c_str());
    
    std::string outcanvname="looseVtightJetID_";
    if (jetptvar=="jtpt")outcanvname+="jtpt";
    else if (jetptvar=="rawpt")outcanvname+="rawpt";
    if(isData)outcanvname="Data_"+outcanvname;
    else outcanvname="MC_"+outcanvname;
    outcanvname+="_etabin"+std::to_string(etabin);
    outcanv->Write(outcanvname.c_str());    
  }
  
  temp_canvJet->Close();  
  return;
}

