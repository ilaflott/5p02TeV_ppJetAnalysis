#include "unfoldSpectra.h"

// general settings
const int radius=4;
const bool doToyErrs=false;
const bool debugMode=true;

// Bayesian settings -------------------------------------------------- 
const bool doBayes=true; 
const int kIter = 4;

// SVD settings --------------------------------------------------
const bool doSVD    = true; // testing now
const bool drawPDFs = true; //untested
const int  nKregMax = 9, kRegDraw = 3, param=21;


// ----------------------------------------------------------------------------------------------------------------------
int doMCClosureTests( ){

  // performance tracking
  TStopwatch timer; 
  timer.Start();
  
  // set default error calculation for TH1/2's from here on out + gStyle
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  //// ppData input histos --------------------------------------------------
  //std::cout<<std::endl<<std::endl<<"opening INPUT histos from DATA file"; 
  //if(debugMode)std::cout<<": "<<inFile_Data_name; 
  //std::cout<<std::endl<<std::endl;
  //TFile *fpp_Data = TFile::Open(inFile_Data.c_str());

  // ppMC input histos --------------------------------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"; 
  if(debugMode)std::cout<<": "<<inFile_MC_name; 
  std::cout<<std::endl<<std::endl;
  TFile *fpp_MC = TFile::Open(inFile_MC.c_str());

  //mat
  //TH2F *hmat, *hmat_anabin;
  TH2F *hmat = (TH2F*)fpp_MC->Get(Form("hpp_mcclosure_matrix_R%d_20_eta_20",radius));
  hmat->Print("base");
  TH2F* hmat_anabin = (TH2F*)hmat->Clone(Form("hpp_anaBin_Trans_matrix_HLT_R%d_20_eta_20",radius));
  hmat_anabin->Rebin2D(10, 10);
  hmat_anabin->Print("base");    
  
  //gen
  //TH1F *hgen, *hgen_anabin;
  TH1F *hgen= (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_gen_R%d_20_eta_20",radius));
  hgen->Print("base");
  TH1F *hgen_anabin = (TH1F*)hgen->Clone(Form("hpp_anaBin_mcclosure_gen_JetComb_R%d_20_eta_20",radius));
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins, Form("hpp_anaBin_mcclosure_gen_JetComb_R%d_20_eta_20",radius), ptbins);
  //hgen_anabin = (TH1F*)hgen_anabin->Rebin(10);
  divideBinWidth(hgen_anabin);
  hgen_anabin->Print("base");

  //rec
  //TH1F *hrec, *hrec_anabin;
  TH1F *hrec = (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_data_R%d_20_eta_20",radius));
  hrec->Print("base");
  TH1F *hrec_anabin = (TH1F*)hrec->Clone(Form("hpp_anaBin_mcclosure_rec_JetComb_R%d_20_eta_20",radius));
  hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins, Form("hpp_anaBin_mcclosure_rec_JetComb_R%d_20_eta_20",radius), ptbins);
  //hrec_anabin = (TH1F*)hrec_anabin->Rebin(10);
  divideBinWidth(hrec_anabin);
  hrec_anabin->Print("base");
  //if(dotrigcor) doTrigCorr(hrec_anabin, hMC_turnon);

  //rec_check
  //TH1F *hrec_check=NULL,*hrec_check_anabin=NULL;
  TH1F *hrec_check = (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_data_train_R%d_20_eta_20",radius));
  hrec_check->Print("base");
  //if(dotrigcor) doTrigCorr(hrec_check, hMC_turnon);
  TH1F *hrec_check_anabin = (TH1F*)hrec_check->Clone(Form("hpp_mcclosure_JetComb_data_train_R%d_20_eta_20",radius));
  hrec_check_anabin = (TH1F*)hrec_check_anabin->Rebin(nbins, Form("hpp_mcclosure_JetComb_data_train_R%d_20_eta_20",radius), ptbins);
  //hrec_check_anabin = (TH1F*)hrec_check_anabin->Rebin(10);
  divideBinWidth(hrec_check_anabin);


  // output histos --------------------------------------------------
  std::cout<<std::endl<<std::endl<<"creating NEW OUTPUT histos..."<<std::endl;

  //gen
  if(debugMode)std::cout<<"opening  gen/genresponse histograms..."<<std::endl;
  //TH1F *hgen_resp, *hgen_resp_anabin;
  TH1F *hgen_resp = new TH1F( Form("hpp_gen_responseR%d_20_eta_20",radius),"", 
			hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  hgen_resp->Sumw2();
  hgen_resp->Print("base");
  TH1F *hgen_resp_anabin = new TH1F(Form("hpp_gen_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hgen_resp_anabin->Sumw2();
  hgen_resp_anabin->Print("base");
  
  //rec
  if(debugMode)std::cout<<"opening  rec/recresponse histograms..."<<std::endl;
  //TH1F *hrec_resp, *hrec_resp_anabin; 
  TH1F *hrec_resp = new TH1F(Form("hpp_rec_response_R%d_20_eta_20",radius),"", 
		       hrec->GetNbinsX(), hrec->GetXaxis()->GetXmin(), hrec->GetXaxis()->GetXmax());
  hrec_resp->Sumw2();
  hrec_resp->Print("base");
  TH1F *hrec_resp_anabin = new TH1F(Form("hpp_rec_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hrec_resp_anabin->Sumw2();
  hrec_resp_anabin->Print("base");
  
  // open output file for later
  const std::string outFile="PP_5p02_doMCClosureTests.root";
  TFile* fout = new TFile(outFile.c_str(),"RECREATE");
  fout->cd();

  // unfolding settings --------------------------------------------------
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
  if(doToyErrs) errorTreatment = RooUnfold::kCovToy;
  std::cout<<std::endl<<"errorTreatment: "<<errorTreatment<<std::endl;
  
  // Bayesian unfolding -------------------------------------------------- 
  if(doBayes){
    std::cout<<std::endl<<" beginning Bayesian unfolding..."<<std::endl;
    std::cout<<"calling RooUnfoldResponse "<<std::endl;
    RooUnfoldResponse roo_resp( hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_R%d",radius));

    std::cout<<"calling RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );

    TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);
    hunf->SetName("PP_Bayesian_Unfolded_Spectra");

    TH1F *hratio = (TH1F*)hunf->Clone(Form("MCClosure_test_oppside_Bayes_R%d_20_eta_20", radius));
    //hratio = (TH1F*)hunf->Clone(Form("MCClosure_test_oppside_Bayes_R%d_20_eta_20", radius));
    hratio->Print("base");

    std::cout<<"dividing hunf and hgen_anabin..."<<std::endl;
    //hgen_anabin->Print("base");
    hratio->Divide(hgen_anabin);
    hratio->SetMarkerStyle(24);
    hratio->SetMarkerColor(kRed);

    std::cout<<"calling RooUnfoldBayes for MC Closure..."<<std::endl;
    RooUnfoldBayes unf_bayes_check(&roo_resp, hrec_check_anabin, kIter);
    
    TH1F *hunf_check = (TH1F*)unf_bayes_check.Hreco(errorTreatment);
    hunf_check->Print("base");

    TH1F *hratio_check = (TH1F*)hunf_check->Clone(Form("MCClosure_test_sameside_Bayes_R%d_20_eta_20", radius));
    hratio_check->Print("base");

    std::cout<<"dividing hunf_check and hgen_anabin..."<<std::endl;
    //hgen_anabin->Print("base");
    hratio_check->Divide(hgen_anabin);
    hratio_check->SetMarkerStyle(24);
    hratio_check->SetMarkerColor(kBlack);
    
    std::cout<<"writing bayesian unfolding output..."<<std::endl;
    //fout->cd();

    hmat->Write();
    hunf->Write();
    hratio->Write();
    hunf_check->Write();
    hratio_check->Write();

    std::cout<<"Bayesian MC Closure test done!"<<std::endl;
  } // end bayesian specifics
  else  std::cout<<std::endl<<" skipping Bayesian MC Closure test..."<<std::endl;


  // SVD unfolding --------------------------------------------------
  if(doSVD){
    std::cout<<std::endl<<" beginning SVD unfolding..."<<std::endl;

    // not sure what this is for, some kinda regularization parameter?
    int kReg[nKregMax];  //{param-4, param-3, param-2, param-1, param, param+1, param+2, param+3, param+4};
    for(int i=(-4); (i+4)<nKregMax; ++i) kReg[i+4]=param+i;
    
    TH1F *hunf_svd[nKregMax];
    TH1F *hratio_svd[nKregMax];
    TH1F *hunf_svd_check[nKregMax];
    TH1F *hratio_svd_check[nKregMax];

    TH1 *hrec_folded_ratio[nKregMax]; 
    TH1 *hrec_unfolded_ratio[nKregMax];
    TH2 *hPearsonSVDPriorMeas[nKregMax]; 
    TH1 *hFoldedSVDPriorMeas[nKregMax];
    TH2 *hPearsonSVDPriorTruth[nKregMax];
    TH1 *hFoldedSVDPriorTruth[nKregMax];

    std::cout<<"calling RooUnfoldResponse"<<std::endl;
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_anabin_R%d",radius));

    // declare the canvases+legends for later
    std::cout<<"opening TCanvases..."<<std::endl;
    TCanvas *cRatio = new TCanvas("cRatio","",1200,1000);    cRatio->Divide(3,3);      
    TCanvas *cSpectra = new TCanvas("cSpectra","",1200,1000);    cSpectra->Divide(3,3);
    TCanvas *cSpectraCheck = new TCanvas("cSpectraCheck","",1200, 1000);    
    TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","",1200,1000);    cPearsonMatrixIter->Divide(3,3);
    TCanvas *c11 = new TCanvas("c11","c11: Singular Values", 1200, 1000);
    TCanvas *cdi = new TCanvas("cdi","cdi: di vectors", 1200, 1000);
    TLegend *leg[nKregMax],*leg1[nKregMax];      

    std::cout<<"entering SVD Unfolding Loop..."<<std::endl;
    for(int kr = 0; kr < nKregMax; ++kr){
      std::cout<<"int kr= "<<kr<<std::endl;
  
      std::cout<<"calling RooUnfoldSvd"<<std::endl;
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
  
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      hunf_svd[kr]->SetName(Form("PP_SVD_Unfolded_Spectra_kReg%d",kReg[kr]));
      //hunf_svd[kr]->SetName(Form("h_Unfolded_mcclosure_spectra_SVD_kReg%d_R%d_20_eta_20", kReg[kr],radius));
      hunf_svd[kr]->Print("base");
  
      hratio_svd[kr] = (TH1F*)hunf_svd[kr]->Clone(Form("MCClosure_test_oppside_SVD_kReg%d_R%d_20_eta_20", kReg[kr], radius));
      hratio_svd[kr]->Divide(hgen_anabin);
      hratio_svd[kr]->SetMarkerStyle(33);
      hratio_svd[kr]->SetMarkerColor(kRed);
  
      // get covariance matrix and calculate pearson coefficients
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);
      //std::cout<<"&covmat="<<&covmat<<std::endl;
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat);
      //std::cout<<"pearson="<<pearson<<std::endl;
      //std::cout<<"&pearson="<<&pearson<<std::endl;

      hPearsonSVDPriorMeas[kr] = new TH2D(*pearson);
      hPearsonSVDPriorMeas[kr]->SetTitle(" ");
      hPearsonSVDPriorMeas[kr]->SetName(Form("hPearsonSVDPriorMeas_kReg%d", kReg[kr]));
      hPearsonSVDPriorMeas[kr]->SetMinimum(-1.);
      hPearsonSVDPriorMeas[kr]->SetMaximum(1.);
      hPearsonSVDPriorMeas[kr]->GetZaxis()->SetLabelSize(0.06);
      
      //Apply to Truth
      std::cout<<"applying roo_resp to truth dist hunf_svd[kr="<<kr<<"]"<<std::endl;
      hFoldedSVDPriorMeas[kr] = roo_resp.ApplyToTruth(hunf_svd[kr]);
      hFoldedSVDPriorMeas[kr]->SetName(Form("hFoldedSVDPriorMeas_kReg%d", kReg[kr]));  
  
      cPearsonMatrixIter->cd(kr+1);
  
      hPearsonSVDPriorMeas[kr]->SetTitle(Form("kReg = %d",kReg[kr]));
      hPearsonSVDPriorMeas[kr]->SetMinimum(-1.);
      hPearsonSVDPriorMeas[kr]->SetMaximum(1.);
      hPearsonSVDPriorMeas[kr]->SetAxisRange(0, 35, "X");
      hPearsonSVDPriorMeas[kr]->SetAxisRange(0, 35, "Y");
      hPearsonSVDPriorMeas[kr]->Draw("colz");
  
      cSpectra->cd(kr+1);
      cSpectra->cd(kr+1)->SetLogy();
  
      hrec_anabin->SetTitle(Form("kReg = %d",kReg[kr]));
      hrec_anabin->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_anabin->SetMarkerStyle(24);
      hrec_anabin->SetMarkerColor(kBlack);
      hrec_anabin->SetAxisRange(45, 1000, "X");
      hrec_anabin->Draw();
  
      hunf_svd[kr]->SetMarkerStyle(33);
      hunf_svd[kr]->SetMarkerColor(kBlue);
      hunf_svd[kr]->Draw("same");
  
      hFoldedSVDPriorMeas[kr]->SetLineStyle(33);
      hFoldedSVDPriorMeas[kr]->SetLineColor(kRed);
      hFoldedSVDPriorMeas[kr]->Draw("same");
  
      leg[kr] = new TLegend(0.1, 0.1, 0.30, 0.30, NULL,"BRNDC");
      leg[kr]->AddEntry(hrec_anabin,"Measured","pl");
      leg[kr]->AddEntry(hunf_svd[kr],Form("Unfolded, kr = %d", kReg[kr]),"pl");
      leg[kr]->AddEntry(hFoldedSVDPriorMeas[kr],Form("Folded kr = %d", kReg[kr]),"pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg[kr]->SetTextSize(0.04); 
      leg[kr]->Draw();
      	
      cRatio->cd(kr+1);
  
      hrec_folded_ratio[kr] = (TH1F*)hFoldedSVDPriorMeas[kr]->Clone(Form("ratio_folded_with_measured_kReg%d",kReg[kr]));
      hrec_folded_ratio[kr]->SetTitle(Form("kReg = %d",kReg[kr]));
      hrec_folded_ratio[kr]->Divide(hrec_anabin);
      hrec_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_folded_ratio[kr]->SetMarkerColor(kRed);
      hrec_folded_ratio[kr]->SetAxisRange(45, 1000, "X");
      hrec_folded_ratio[kr]->SetAxisRange(0.5, 1.5, "Y");
      hrec_folded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_folded_ratio[kr]->Draw();
  
      hrec_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone(Form("ratio_unfolded_with_measured_kReg%d",kr));
      hrec_unfolded_ratio[kr]->Divide(hrec_anabin);
      hrec_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_unfolded_ratio[kr]->SetMarkerColor(kBlue);
      hrec_unfolded_ratio[kr]->SetTitle(Form("kReg = %d",kReg[kr]));
      hrec_unfolded_ratio[kr]->Print("base");
      hrec_unfolded_ratio[kr]->Draw("same");
  
      leg1[kr] = new TLegend(0.1, 0.1, 0.30, 0.3, NULL,"BRNDC");
      leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],Form("Unfolded/Measured, kr = %d", kReg[kr]),"pl");
      leg1[kr]->AddEntry(hrec_folded_ratio[kr],Form("Folded/Measured kr = %d", kReg[kr]),"pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg1[kr]->SetTextSize(0.04); 
      leg1[kr]->Draw();
  	
      //Get and print singular values and d_i vector, print also on screen
      //Note that these do not depend on the regularization. 
      //The opposite: they tell you which regularization to use!
      if(kr == 0){
  
  	TSVDUnfold *svdUnfold = unf_svd.Impl();
  	TH1D *hdi = (TH1D*)svdUnfold->GetD();
  	TH1 *hSVal = (TH1*)svdUnfold->GetSV();
  	//TH1D  *hSValClone = (TH1D*)hSVal->Clone("hSValClone");
  
  	std::cout << "print singular values: " << std::endl;
  	for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)std::cout<<"bin: "<<bin<<", SV: "<<hSVal->GetBinContent(bin)<< std::endl;
  	
  	std::cout << "print di vector: " <<  std::endl;
  	for(int bin=1; bin<=hdi->GetNbinsX(); bin++)std::cout<<"i: "<<bin<<", di: "<<hdi->GetBinContent(bin)<<std::endl;
  	
  	c11->cd();
  	gPad->SetLogy();
  	
  	hSVal->SetTitle(" ");
  	hSVal->SetAxisRange(0,35,"X");
  	hSVal->SetXTitle("singular values");
  	hSVal->DrawCopy();
  	
  	drawText("PP 5 TeV", 0.608173, 0.8459761, 22);
  
  	cdi->cd();
  	gPad->SetLogy();
  	
  	hdi->SetTitle(" ");
  	hdi->SetAxisRange(0,35,"X");
  	hdi->SetXTitle("|d_{i}^{kreg}|");
  	hdi->DrawCopy();
  	
  	drawText("PP 5 TeV", 0.608173, 0.8459761, 22);
      } // end kr==0 specific
  
      hrec_unfolded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Print("base");

      // the closure tests for svd unfolding
      std::cout<<"calling RooUnfoldSVD for MCClosure tests"<<std::endl;
      RooUnfoldSvd unf_svd_check(&roo_resp, hrec_check_anabin, kReg[kr]);
      
      hunf_svd_check[kr] = (TH1F*)unf_svd_check.Hreco(errorTreatment);
      hunf_svd_check[kr]->SetName(Form("h_Unfolded_mcclosure_sameside_SVD_kReg%d_check_R%d_20_eta_20", kReg[kr],radius));
      hunf_svd_check[kr]->Print("base");
      
      hratio_svd_check[kr] = (TH1F*)hunf_svd_check[kr]->Clone(Form("MCClosure_test_sameside_SVD_kReg%d_R%d_20_eta_20", 
								   kReg[kr], radius));
      hratio_svd_check[kr]->Print("base");
      hratio_svd_check[kr]->Divide(hgen_anabin);
      hratio_svd_check[kr]->SetMarkerStyle(33);
      hratio_svd_check[kr]->SetMarkerColor(kBlack);
      
      // get covariance matrix and calculate pearson coefficients
      TMatrixD covmat_truth = unf_svd_check.Ereco(errorTreatment);
      TMatrixD *pearson_truth = CalculatePearsonCoefficients(&covmat_truth);
      
      std::cout << "print pearson coefficients" << std::endl;      

      hPearsonSVDPriorTruth[kr] = new TH2D(*pearson_truth);
      hPearsonSVDPriorTruth[kr]->SetName(Form("hPearsonSVDPriorTruth_kReg%d",kReg[kr]));
      hPearsonSVDPriorTruth[kr]->SetMinimum(-1.);
      hPearsonSVDPriorTruth[kr]->SetMaximum(1.);
      hPearsonSVDPriorTruth[kr]->GetZaxis()->SetLabelSize(0.06);
      
      std::cout<<"applying roo_resp to truth dist hunf_svd_check[kr="<<kr<<"]"<<std::endl;
      hFoldedSVDPriorTruth[kr] = roo_resp.ApplyToTruth(hunf_svd_check[kr]);
      hFoldedSVDPriorTruth[kr]->SetName(Form("hFoldedSVDPriorTruth_kReg%d", kReg[kr]));
  
      std::cout<<std::endl<<std::endl<<std::endl<<"MC: kReg = "<<kReg[kr]<<std::endl<<std::endl<<std::endl;

      TSVDUnfold *svdUnfold = unf_svd_check.Impl();      
      TH1 *hSVal = svdUnfold->GetSV();
      TH1D *hdi = svdUnfold->GetD();
      //TH1D  *hSValClone = (TH1D*)hSVal->Clone("hSValClone");
      
      std::cout << "print singular values " << std::endl;
      for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
	std::cout << "bin: " << bin << "  SV: " << hSVal->GetBinContent(bin) << std::endl;
      
      std::cout << "print di vector " <<  std::endl;
      for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
	std::cout << "i: " << bin << "  di: " << hdi->GetBinContent(bin) << std::endl;

    }// kReg loop

    // loop over histos pointers in arrays to write
    //fout->cd();
    for(int kr = 0; kr<nKregMax; ++kr){
      // for svd unfolding
      hunf_svd[kr]->Write();
      hratio_svd[kr]->Write();
      hPearsonSVDPriorMeas[kr]->Write();
      hFoldedSVDPriorMeas[kr]->Write();
      // for MC CLosure tests
      hunf_svd_check[kr]->Write();
      hratio_svd_check[kr]->Write();
      hPearsonSVDPriorTruth[kr]->Write();
      hFoldedSVDPriorTruth[kr]->Write();
    }

    // draw PDF files if desired
    if(drawPDFs){

      std::string type="MCClosure";//should only be Data or MC
      cPearsonMatrixIter->SaveAs(Form("Pearson_Matrix_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");
      cSpectra->SaveAs(Form("Spectra_meas_unf_fol_%s_R%d_plots.pdf",type.c_str(), radius),"RECREATE");
      cRatio->SaveAs(Form("Ratio_meas_unf_fol_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");          
      cdi->SaveAs(Form("di_Vectors_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
      c11->SaveAs(Form("Singular_Values_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
      
      cSpectraCheck->cd();
      
      hrec_folded_ratio[kRegDraw]->SetTitle(" ");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(0.4, 1.2, "Y");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(45, 1000, "X");
      
      TLegend * leg2 = new TLegend(0.1, 0.1, 0.30, 0.3, NULL,"BRNDC");
      leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Unfolded/Measured","pl");
      leg2->AddEntry(hrec_folded_ratio[kRegDraw],"Folded/Measured","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2->SetTextSize(0.04); 
      
      hrec_folded_ratio[kRegDraw]->Draw();
      hrec_unfolded_ratio[kRegDraw]->Draw("same");
      drawText("PP 5 TeV", 0.608173, 0.8459761, 22);
      drawText(Form("kReg = %d",kReg[kRegDraw]), 0.60, 0.75, 22);
      leg2->Draw();    
      
      cSpectraCheck->SaveAs(Form("Spectra_Ratio_withkRegDraw_%s_R%d.pdf", type.c_str(), radius),"RECREATE");	

      // MCClosure Plot
      std::string c1_filename="PP_5p02TeV_"+type+"_R"+std::to_string(radius)+"_SVDOnly.pdf";
      TCanvas *c1 = new TCanvas("c1","Spectra",1300,1000);
      c1->cd();
      
      TH1F *hDum = new TH1F("hDum","",250, 30, 900);
      hDum->SetYTitle("Unfolding Closure");//y-axis
      hDum->GetYaxis()->SetNdivisions(610); 
      hDum->GetYaxis()->SetLabelFont(43);
      hDum->GetYaxis()->SetTitleFont(43);
      hDum->GetYaxis()->SetLabelSize(20);
      hDum->GetYaxis()->SetTitleSize(22);
      hDum->GetYaxis()->SetTitleOffset(2.6);
      
      hDum->SetXTitle("Gen Jet p_{T} (GeV/c)");//x-axis
      hDum->GetXaxis()->SetLabelFont(43);
      hDum->GetXaxis()->SetTitleFont(43);
      hDum->GetXaxis()->SetLabelSize(20);
      hDum->GetXaxis()->SetTitleSize(22);
      hDum->GetXaxis()->SetTitleOffset(2.3);
      hDum->GetXaxis()->SetNoExponent();
      hDum->GetXaxis()->SetMoreLogLabels();
      hDum->SetAxisRange(0.8,1.2,"Y");
      
      hDum->Draw("hist");
      
      TLegend *leg0 = new TLegend(0.25, 0.15, 0.60, 0.35, NULL,"BRNDC");
      leg0->SetTextSize(0.07);

      if(doBayes){
	c1_filename="PP_5p02TeV_"+type+"_R"+std::to_string(radius)+"_BayesAndSVD.pdf";

	TH1F* hratio=(TH1F*)fout->Get(Form("MCClosure_test_oppside_Bayes_R%d_20_eta_20", radius));
	hratio->Print("base");
	hratio->Draw("psame");

	TH1F* hratio_check=(TH1F*)fout->Get(Form("MCClosure_test_sameside_Bayes_R%d_20_eta_20", radius));
	hratio_check->Print("base");
	hratio_check->Draw("psame");

	leg0->SetTextSize(0.04);
	leg0->AddEntry(hratio,"OppSide Bayesian","p");
	leg0->AddEntry(hratio_check,"SameSide Bayesian","p");
      }

      hratio_svd[kRegDraw]->Draw("psame");
      hratio_svd_check[kRegDraw]->Draw("psame");
      
      drawText("PP, 5 TeV", 0.648173, 0.8459761, 22);
      
      TLine *line = new TLine(30,1,1000,1);
      line->SetLineStyle(2); 
      line->SetLineWidth(2);
      line->Draw();
      
      drawText(Form("kReg = %d", kReg[kRegDraw]), 0.64, 0.7, 21);
      drawText("PYTHIA",0.22,0.88,21);
      
      leg0->AddEntry(hratio_svd[kRegDraw],"OppSide SVD","p");
      leg0->AddEntry(hratio_svd_check[kRegDraw],"SameSide SVD","p");
      leg0->Draw();
      
      c1->SaveAs( c1_filename.c_str(), "RECREATE");
    }// end drawPDFs
  }// end SVD specific
  else  std::cout<<std::endl<<" skipping SVD unfolding..."<<std::endl;

  std::cout<<"writing histos to output file"<<std::endl;
  //fout->cd();
  hgen->Write();
  hgen_resp->Write();
  hrec->Write();    
  hrec_resp->Write();
  hrec_check->Write();
  //fout->Write();

  std::cout<<"cleaning up"<<std::endl;
  fout->Close();
  fpp_MC->Close();  

  return 0;
} // end doMCClosureTests

//-----------------------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]){
  int rStatus=-1;//-1 until input right
  if(argc!=1){
    std::cout<<" bad usage, do ";
    std::cout<<"./doMCClosureTests"<<std::endl<<std::endl;
    return rStatus;
  }
  rStatus=1;//1 until a function returns 0 (normal exit by my convention)
  rStatus=doMCClosureTests();
  std::cout<<"done! return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}











