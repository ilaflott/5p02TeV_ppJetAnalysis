void setupJetTrigSpectraRatioCanvas(TCanvas* canv, 
				    TPad* specpad, 
				    TPad* ratpad){

  canv->cd();

  specpad->SetLogx(1);
  specpad->SetLogy(1);
  specpad->SetGridx(1);
  specpad->SetGridy(1);
  
  specpad->SetBottomMargin(0);
  specpad->Draw();
  
  ratpad->SetLogx(1);
  ratpad->SetLogy(0);
  ratpad->SetGridx(1);
  ratpad->SetGridy(0);
  
  ratpad->SetTopMargin(0);
  ratpad->SetBottomMargin(0.3);
  
  ratpad->Draw();
  return;
}



