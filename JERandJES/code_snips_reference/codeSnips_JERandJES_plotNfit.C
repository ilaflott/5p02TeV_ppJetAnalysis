  /*
  // PbPb plot
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  TLegend *leg0 = getLegend(0.1950635,0.4625641,0.8242194,0.7497436);
  //TLegend *leg0 = getLegend(0.4418803,0.3887016,0.8692251,0.6763714);  
  if(wJetID) leg0->SetHeader("w JetID");
  else leg0->SetHeader("w/o JetID");
  leg0->SetTextSize(0.10);
  leg0->AddEntry( hMean[0][0], Form("pp %s",algname.c_str()) , "p" );
  //leg0->AddEntry( hMean[0][0]   , Form("PbPb ak%s",(algo+srad[nj]+jetType).c_str()), "p" );
  //leg0->AddEntry(hMean_resd[0][0],"Resd. Correc.","p");
  
  TLine *l0 = new TLine(xmin,1.00,xmax,1.0);   l0->SetLineStyle(2);  l0->SetLineWidth(2);
  TLine *l1 = new TLine(xmin,0.99,xmax,0.99);  l1->SetLineStyle(2);
  TLine *l2 = new TLine(xmin,1.01,xmax,1.01);  l2->SetLineStyle(2);

  TCanvas *c11[Nrad];     int ipad=0;  int maxr=2;  int maxc=7;
  for(int nj=0; nj<Nrad; nj++){
    c11[nj] = new TCanvas( ( "c11_"+std::to_string(nj) ).c_str() , algname.c_str(), 1695, 530);
    makeMultiPanelCanvas(c11[nj], maxc, maxr, 0.0,0.0, 0.22,0.22, 0.02,0);

    for(int ic=ncen; ic>=0; --ic){    //for(int ic=ncen; ic>=0; ic--){
      c11[nj]->cd(++ipad);
      gPad->SetLogx();

      hSigma[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hSigma[nj][ic]->Draw("p");
      hSigma[nj][ic]->SetMarkerColor(2);
      hSigma[nj][ic]->SetLineColor(2);
      hSigma[nj][ic]->Draw("psame");

      // if(ipad==2 && ic==ncen-1){
      // 	drawText2("PYTHIA+HYDJET",0.19,0.85,21);
      // 	leg0->Draw();
      // }
      // if(ipad==3 && ic==4){
      // 	drawText2("#sqrt{s_{NN}} = 5.02 TeV,|#eta| < 2",0.10,0.85,19);
      // }
      // if(ipad==1 && ic==5){
      // 	if(wJetID)drawText2(Form("%s, w JetID",(srad[]+jetType).c_str()),0.28,0.60,21);
      // 	else drawText2(Form("%s, w/o JetID",(srad[]+jetType).c_str()),0.28,0.60,21);
      // }

      drawText2(ccent[ic],0.28,0.70,21);
      c11[nj]->cd(ipad+maxc);
      gPad->SetLogx();
    
      // if(ipad!=1)
      hMean[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hMean[nj][ic]->GetYaxis()->SetRangeUser(0.8,1.05);
      hMean[nj][ic]->Draw("p");
      
      hMean[nj][ncen]->SetMarkerColor(2);
      hMean[nj][ncen]->SetLineColor(2);
      hMean[nj][ncen]->Draw("psame");

      l1->Draw();
      l0->Draw();
      l2->Draw();
    }
    if(iSave) c11[nj]->SaveAs(Form("Plots_JEC_%d_PbPb_%s.pdf",wJetID,algname.c_str()));
  }

  
  // PbPb plot
  TCanvas *c99[Nrad][ncen];
  for(int nj=0;nj<Nrad;nj++){
    for(int ic=ncen-1;ic>=0; --ic){
      c99[nj][ic] = new TCanvas(Form("c99%d_%d",nj,ic),Form("%s %s ",srad[nj].c_str(),ccent[ic]),100,102,1399,942);
      c99[nj][ic]->Divide(4,4,0,0);

      int ipad=0;
      for(int ip=0;ip<nbins_pt;ip++){      
	c99[nj][ic]->cd(++ipad);
	
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	
	// hrsp[nj][ic][ip]->Scale(1./hrsp[nj][ic][ip]->Integral());
	// if(ic==ncen-3)hrsp[nj][ncen-1][ip]->Scale(1./hrsp[nj][ncen-1][ip]->Integral());	

	double ymax = TMath::Max(hrsp[nj][ic][ip]->GetMaximum(),hrsp[nj][ncen-1][ip]->GetMaximum());
	//double ymax = hrsp[nj][ic][ip]->GetMaximum();

	hrsp[nj][ic][ip]->SetMaximum(1.45*ymax);
	hrsp[nj][ic][ip]->SetMinimum(0);
	hrsp[nj][ic][ip]->SetTitle(0);

	hrsp[nj][ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
	hrsp[nj][ic][ip]->GetXaxis()->SetTitleFont(42);
	hrsp[nj][ic][ip]->GetXaxis()->SetLabelFont(42);
	hrsp[nj][ic][ip]->GetXaxis()->SetLabelSize(0.08);
	hrsp[nj][ic][ip]->GetXaxis()->SetTitleSize(0.07);

	hrsp[nj][ic][ip]->GetYaxis()->SetTitle("");
	hrsp[nj][ic][ip]->GetYaxis()->SetTitleFont(42);
	hrsp[nj][ic][ip]->GetYaxis()->SetLabelFont(42);
	hrsp[nj][ic][ip]->GetYaxis()->SetLabelSize(0.08);
      
	hrsp[nj][ic][ip]->SetMarkerStyle(20);
	hrsp[nj][ic][ip]->SetMarkerColor(1);
	hrsp[nj][ic][ip]->SetLineColor(1);
	hrsp[nj][ic][ip]->SetMarkerSize(1.3);
	hrsp[nj][ic][ip]->Draw("p");  

	hrsp[nj][ncen-1][ip]->SetMarkerStyle(20);
	hrsp[nj][ncen-1][ip]->SetMarkerColor(2);
	hrsp[nj][ncen-1][ip]->SetLineColor(2);
	hrsp[nj][ncen-1][ip]->SetMarkerSize(1.3);
	// hrsp[nj][ncen-1][ip]->Draw("psame");  
      
	//std::ostringstream strs; 
	//strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
	//std::string spt = strs.str();
	std::string spt=std::to_string(ptbins[ip])+"< p_{T} (GeV/c) <"+std::to_string(ptbins[ip+1]);
	if(ipad==1){
	  //drawText2(Form("%s",(srad[nj]+jetType).c_str()),0.28,0.90,21);      
	  drawText2(ccent[ic],0.75,0.87,21);	  
	  drawText2(spt.c_str(),0.22,0.80,21);		
	  if(wJetID)drawText2("w JetID",0.22,0.65,21);
	  else drawText2("w/o JetID",0.22,0.66,21);	  
	} 
	else drawText2(spt.c_str(),0.17,0.80,21);
      }
      if(iSave)	c99[nj][ic]->SaveAs(Form("Plots_IndDist_PbPb_%s_%s.pdf",ccent[ic],algname.c_str()));
    }
  }


  // output file
  fout->cd( );
  const std::string dirname=algname+"JetAnalyzer";
  fout->mkdir( dirname.c_str() );
  fout->cd( dirname.c_str() );

  // PbPb write output
  for(int nj=0; nj<Nrad; nj++){
    for(int ic=0; ic<ncen; ic++){
      hMean[nj][ic]->SetName(  Form("hMean_PbPb_R%s_PuPF_%s", srad[nj].c_str(), cdir[ic]));
      hMean[nj][ic]->SetTitle( Form("hMean_PbPb_R%s_PuPF_%s", srad[nj].c_str(), cdir[ic]));
      hMean[nj][ic]->Write();
      hSigma[nj][ic]->SetName(  Form("hSigma_PbPb_R%s_PuPF_%s", srad[nj].c_str(), cdir[ic]));
      hSigma[nj][ic]->SetTitle( Form("hSigma_PbPb_R%s_PuPF_%s", srad[nj].c_str(), cdir[ic]));
      hSigma[nj][ic]->Write();
    }
  }
  */
