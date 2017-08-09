{
  // TFile * fData = TFile::Open("PP_5p02TeV_Data_ak3PF_20_eta_20.root");
  //TFile * fMC;
  TFile * fData;
//  std::string etabin = "eta32to47";
  std::string etabin = "eta50to50";

  
  //fMC = TFile::Open(Form("ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_12-15-16_jetPlots_jtIDv2_%s/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root", etabin.c_str()));
//  fMC = TFile::Open("/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/readForests/1.10.17_outputCondor/ppMC_Py8_CUETP8M1_QCDjetAllPtBins_ak4PFJets_01-10-17_jetPlots_jtIDv3/Py8_CUETP8M1_QCDjetAllPtBins_ak4PF-allFiles.root");
  fData = TFile::Open("/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/readForests/1.16.17_outputCondor/ppData_HighPtJetTrig_ak4PFJets_01-16-17_jetPlots_jtIDv3L_eta20/HighPtJetTrig_ak4PF-allFiles.root");
  
  std::string etabinlabel = "|#eta| < 2.0";			//check eta
  /*
  std::string etabinlabel = "";
  if(etabin == "eta27to30")
    etabinlabel = "2.7 < |#eta| < 3.0";
  if(etabin == "eta32to47")
    etabinlabel = "3.2 < |#eta| < 4.7";
    */
	
  const int N = 20;
  
  std::string var[N] = {"jtpt" ,"rawpt", "jteta", "jtphi", "jty", "trkMax", "trkSum", "chMax", "chSum","phMax", "phSum", "neMax", "neSum", "eMax", "eSum", "muMax", "muSum", "neN", "chN", "sumN"};
  std::string lable[N] = {"Reco Jet p_{T} [GeV/c]" ,"Raw Jet p_{T} [GeV/c]", "Jet #eta", "Jet #phi", "Jet rapidity (y)", "Track Max/Raw p_{T}", "Track Sum/Raw p_{T}", "Charged Max/Raw p_{T}", "Charged Sum/Raw p_{T}","Photon Max/Raw p_{T}", "Photon Sum/Raw p_{T}", "Neutral Max/Raw p_{T}", "Neutral Sum/Raw p_{T}", "Electron Max/Raw p_{T}", "Electron Sum/Raw p_{T}", "Muon Max/Raw p_{T}", "Muon Sum/Raw p_{T}", "Number of Neutrals", "Number of Charged", "Total N in jet"};
  TH1F * hJetQA[2][2][N];
  TH1F * hRatio[N];

  TCanvas * cJetQA[N];
  TLegend* leg0 = new TLegend(0.53,0.63,0.83,0.85);
  TLegend* leg1 = new TLegend(0.53,0.63,0.83,0.85);

  // variable loop
  for(int i = 0; i<N; ++i){ 
    // jet ID loop, 0- no jetID, 1- yes JetID
    for(int k = 0; k<2; ++k){
      
      cout<<"Running for Variable: "<<Form("hJetQA_%dwJetID_%s",k,var[i].c_str())<<endl;
      hJetQA[k][0][i] = (TH1F*)fData->Get(Form("hJetQA_%dwJetID_%s",k,var[i].c_str()));
	  hJetQA[k][0][i]->SetYTitle("Arbitrary Units");
	  hJetQA[k][0][i]->GetYaxis()->SetTitleOffset(1.2);
      hJetQA[k][0][i]->Print("base");
      if(i<=16)
	hJetQA[k][0][i]->Rebin(5);
      if(k == 0){
	hJetQA[k][0][i]->SetFillColor(4);
	hJetQA[k][0][i]->SetLineColor(4);
	hJetQA[k][0][i]->SetFillStyle(3354);
      }else {
	hJetQA[k][0][i]->SetLineColor(2);
	hJetQA[k][0][i]->SetMarkerColor(2);
	hJetQA[k][0][i]->SetMarkerStyle(33);
      }
      hJetQA[k][0][i]->GetXaxis()->CenterTitle();
      hJetQA[k][0][i]->GetYaxis()->CenterTitle();

      // if(i>=4){
      // 	hJetQA[k][0][i]->Scale(1./hJetQA[k][0][i]->Integral());
      // 	// hJetQA[k][0][i]->Scale(1./hJetQA[k][0][i]->Integral());
      // }
      // else {
      // 	// hJetQA[k][0][i]->Scale(1./(12.71 * 1e9));
      // 	// hJetQA[k][1][i]->Scale(1./(1e3));
      // }

      // if(i==2 || i==3){
      // 	hJetQA[k][0][i]->Scale((float)10/1.1);	
      // }
    }// k loop
    
    cJetQA[i] = new TCanvas(Form("cJetQA_%d",i),"",800,800);
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->SetGridy();         // horizontal grid
    if(i!=2 || i!=3) pad1->SetLogy();          // Vertical Logy
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    hJetQA[0][0][i]->SetStats(0);          // No statistics on upper plot
    hJetQA[0][0][i]->Draw("h");               // Draw h1
    hJetQA[1][0][i]->Draw("psame");               // Draw h1
    	
    if(i == 0){
	leg0->AddEntry("","CMS Data","");
	leg0->AddEntry("","PP #sqrt{s_{NN}}=5.02 TeV, 2015","");
	leg0->AddEntry("","ak4PF Jets","");											//check radius
	leg0->AddEntry("",Form("%s",etabinlabel.c_str()),"");
	leg0->AddEntry("","Modified JetID from 5.02 TeV","");
	leg0->AddEntry(hJetQA[1][0][0],"w/ Jet ID","lpf");
	leg0->AddEntry(hJetQA[0][0][0],"no Jet ID","lpf");
	leg0->SetTextSize(0.025);  
    }
 
    leg0->Draw();
  
    // lower plot will be in pad
    cJetQA[i]->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->SetGridy(); // horizontal grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    
    // Define the ratio plot
    hRatio[i] = (TH1F*)hJetQA[1][0][i]->Clone(Form("hRatio_%s",var[i].c_str()));
    hRatio[i]->Divide(hJetQA[0][0][i]);
    hRatio[i]->SetTitle("");
    hRatio[i]->SetYTitle("Ratio: w/JetID / w/o JetID");
    hRatio[i]->SetXTitle(Form("%s",lable[i].c_str()));
    hRatio[i]->SetMinimum(0);  // Define Y ..
    hRatio[i]->SetMaximum(2); // .. range
    hRatio[i]->SetStats(0);      // No statistics on lower plot
    hRatio[i]->SetLineColor(kBlack);
    hRatio[i]->SetMarkerStyle(21);
    hRatio[i]->Draw("h");
    hRatio[i]->GetYaxis()->SetNdivisions(505);
    hRatio[i]->GetYaxis()->SetTitleSize(15);
    hRatio[i]->GetYaxis()->SetTitleFont(43);
    hRatio[i]->GetYaxis()->SetTitleOffset(1.55); //was 1.55
    hRatio[i]->GetYaxis()->SetLabelFont(43); 
    hRatio[i]->GetYaxis()->SetLabelSize(15);
    hRatio[i]->GetXaxis()->SetTitleSize(20);
    hRatio[i]->GetXaxis()->SetTitleFont(43);
    hRatio[i]->GetXaxis()->SetTitleOffset(4.);
    hRatio[i]->GetXaxis()->SetLabelFont(43); 
    hRatio[i]->GetXaxis()->SetLabelSize(15);

	
//    cJetQA[i]->SaveAs(Form("%s/JetQAPlots_pp_5TeV_%s.pdf", etabin.c_str(),var[i].c_str()),"RECREATE");    
    cJetQA[i]->SaveAs(Form("JetIDEff_ppData_AllPtBins_ak4PFJets_eta20_01-17-17_jetPlots_jtIDv3_5TeV_%s.pdf", var[i].c_str()),"RECREATE");    //check radius and eta
  }// var loop

}// jet ID loop
