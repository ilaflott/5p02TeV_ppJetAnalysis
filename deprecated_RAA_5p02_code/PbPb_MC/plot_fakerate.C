{

  const int nRad = 3;
  int Radii[nRad] = {2,3,4};
  int MarkerStyle = 24;
  int MarkerColor[nRad] = {1, 2, 4};
  
  const int nbins_cent = 6;
  int boundaries_cent[nbins_cent+1] = {0, 5, 10, 30, 50, 70, 90};
  //Int_t nbins_jtpt[nbins_cent]       = {60, 60, 60, 60, 40, 24};
  // Int_t nbins_jtpt[nbins_cent]       = {30, 30, 30, 30, 20, 12};
  // Double_t binLow_jtpt[nbins_cent]   = {50.0, 50.0, 50.0, 50.0, 50.0, 50.0};
  // Double_t binHigh_jtpt[nbins_cent]  = {350.0, 350.0, 350.0, 350.0, 250.0, 170.0};

  // double boundaries_truncate_cent0[nbins_jtpt[0]+1] = {50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0};
  // double boundaries_truncate_cent1[nbins_jtpt[1]+1] = {50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0};
  // double boundaries_truncate_cent2[nbins_jtpt[2]+1] = {50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0};
  // double boundaries_truncate_cent3[nbins_jtpt[3]+1] = {50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0};
  // double boundaries_truncate_cent4[nbins_jtpt[4]+1] = {50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0};
  // double boundaries_truncate_cent5[nbins_jtpt[5]+1] = {50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0};

  double boundaries[] = {50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0};
  int nbins = 6;
  
  TFile * fin[nRad];
  TH1F * hFake[nRad][nbins_cent];
  TH1F * hReco[nRad][nbins_cent];

  TH1F * hFakeRate[nRad][nbins_cent];
  TH1F * hMisIDCorrFactor[nRad][nbins_cent];

  TCanvas * cFakeRate = new TCanvas("cFakeRate","",1200, 1000);
  cFakeRate->Divide(3,2);
  TCanvas * cMisIDCorrFactor = new TCanvas("cMisIDCorrFactor","",1200, 1000);
  cMisIDCorrFactor->Divide(3,2);
  
  for(int rad = 0; rad<nRad; rad++){
    fin[rad] = TFile::Open(Form("PbPb_MC_new_akPu%d_20_eta_20.root", Radii[rad]));
    
    for(int i = 0; i<nbins_cent; ++i){

      hFake[rad][i] = (TH1F*)fin[rad]->Get(Form("hpbpb_DijetBKGSelection_LeadPFJet_lt_A_cent%d",i));
      hReco[rad][i] = (TH1F*)fin[rad]->Get(Form("hpbpb_reco_R%d_20_eta_20_cent%d",Radii[rad], i));

      hFake[rad][i]->Rebin(10);
      hReco[rad][i]->Rebin(10);
      // if(i == 0) hFake[rad][i] = (TH1F*)hFake[rad][i]->Rebin(nbins_jtpt[i],Form("hpbpb_fake_R%d_cent%d",Radii[rad],i),boundaries_truncate_cent0);
      // if(i == 1) hFake[rad][i] = (TH1F*)hFake[rad][i]->Rebin(nbins_jtpt[i],Form("hpbpb_fake_R%d_cent%d",Radii[rad],i),boundaries_truncate_cent1);
      // if(i == 2) hFake[rad][i] = (TH1F*)hFake[rad][i]->Rebin(nbins_jtpt[i],Form("hpbpb_fake_R%d_cent%d",Radii[rad],i),boundaries_truncate_cent2);
      // if(i == 3) hFake[rad][i] = (TH1F*)hFake[rad][i]->Rebin(nbins_jtpt[i],Form("hpbpb_fake_R%d_cent%d",Radii[rad],i),boundaries_truncate_cent3);
      // if(i == 4) hFake[rad][i] = (TH1F*)hFake[rad][i]->Rebin(nbins_jtpt[i],Form("hpbpb_fake_R%d_cent%d",Radii[rad],i),boundaries_truncate_cent4);
      // if(i == 5) hFake[rad][i] = (TH1F*)hFake[rad][i]->Rebin(nbins_jtpt[i],Form("hpbpb_fake_R%d_cent%d",Radii[rad],i),boundaries_truncate_cent5);

      hFake[rad][i]->Print("base");
      hReco[rad][i]->Print("base");
      
      hFakeRate[rad][i] = (TH1F*)hFake[rad][i]->Clone(Form("FakeRate_R%d_cent%d",Radii[rad],i));
      hFakeRate[rad][i]->Divide(hReco[rad][i]);

      hMisIDCorrFactor[rad][i] = (TH1F*)hFakeRate[rad][i]->Clone(Form("hMisIDCorrFactor_R%d_cent%d",Radii[rad], i));
      
      for(int bin = 1; bin<=hMisIDCorrFactor[rad][i]->GetNbinsX(); ++bin){
	double val = hMisIDCorrFactor[rad][i]->GetBinContent(bin);
	val = (double)1./(1.0 - val);
	hMisIDCorrFactor[rad][i]->SetBinContent(bin, val);
      }
      
      cFakeRate->cd(nbins_cent-i);
      cFakeRate->cd(nbins_cent-i)->SetLogy();
      hFakeRate[rad][i] = (TH1F*)hFakeRate[rad][i]->Rebin(nbins, Form("FakeRate_R%d_cent%d",Radii[rad],i),boundaries);
      hFakeRate[rad][i]->SetTitle(Form("%d-%d% %",boundaries_cent[i], boundaries_cent[i+1]));
      hFakeRate[rad][i]->SetXTitle("Jet p_{T} (GeV/c)");
      hFakeRate[rad][i]->SetYTitle("PY+HD, BKG Selection RealJet Mis ID Fake Rate");
      hFakeRate[rad][i]->SetMarkerStyle(MarkerStyle);
      hFakeRate[rad][i]->SetMarkerColor(MarkerColor[rad]);      
      //hFakeRate[rad][i]->SetAxisRange(50, 110, "X");
      hFakeRate[rad][i]->SetAxisRange(1e-2, 1, "Y");
      if(rad == 0) hFakeRate[rad][i]->Draw();
      else hFakeRate[rad][i]->Draw("same");

      cMisIDCorrFactor->cd(nbins_cent-i);
      hMisIDCorrFactor[rad][i] = (TH1F*)hMisIDCorrFactor[rad][i]->Rebin(nbins, Form("MisIDCorrFactor_R%d_cent%d",Radii[rad],i),boundaries);
      hMisIDCorrFactor[rad][i]->SetTitle(Form("%d-%d% %",boundaries_cent[i], boundaries_cent[i+1]));
      hMisIDCorrFactor[rad][i]->SetXTitle("Jet p_{T} (GeV/c)");
      hMisIDCorrFactor[rad][i]->SetYTitle("PY+HD, MB RealJet Mis ID CorrFactor");
      hMisIDCorrFactor[rad][i]->SetMarkerStyle(MarkerStyle);
      hMisIDCorrFactor[rad][i]->SetMarkerColor(MarkerColor[rad]);      
      //hMisIDCorrFactor[rad][i]->SetAxisRange(50, 110, "X");
      hMisIDCorrFactor[rad][i]->SetAxisRange(1, 1.15, "Y");
      if(rad == 0) hMisIDCorrFactor[rad][i]->Draw();
      else hMisIDCorrFactor[rad][i]->Draw("same");

      
    }

  }

  cFakeRate->SaveAs("FakeRate_PYTHAHYDJET.pdf","RECREATE");
  cMisIDCorrFactor->SaveAs("MisIDCorrFactor_PYTHAHYDJET.pdf","RECREATE");

  TFile fout("CombinatorialFakeCut_RealJets_CorrFactor.root","RECREATE");
  fout.cd();
  
  for(int rad = 0; rad<nRad; rad++){
    
    for(int i = 0; i<nbins_cent; ++i){

      hFakeRate[rad][i]->Write();
      hMisIDCorrFactor[rad][i]->Write();
      
    }

  }
  fout.Close();
}
