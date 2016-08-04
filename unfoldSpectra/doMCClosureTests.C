#include "unfoldSpectra.h"

// general settings -------------------------------------------------- 
const std::string baseName="doMCClosureTests_defOut"; // output name
const std::string outRootFile=baseName+".root";

const int radius=4;
const bool doToyErrs=false;
const bool defDebugMode=false;
RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

// Bayesian settings -------------------------------------------------- 
const bool doBayes=true; 
const int kIter = 4;

// SVD settings --------------------------------------------------
const bool doSVD    = true; // testing now
const bool drawPDFs = true; //untested
const int  nKregMax = 9, kRegRange=4, kRegDraw = 3, kRegCenter=21;


int doMCClosureTests( const bool debugMode=defDebugMode){

  // script performance tracking if needed
  TStopwatch timer; 
  timer.Start();
  
  // set error handing + def. stat. info on plots
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  // other settings that are easily forgott
  if(debugMode)std::cout<<std::endl<<"debugMode=TRUE"<<std::endl;
  if(debugMode&&doToyErrs)std::cout<<"using toy errors"<<std::endl;
  if(debugMode&&doToyErrs) errorTreatment = RooUnfold::kCovToy;
  const bool debugPearson=false&&debugMode;

  // ppMC input histos --------------------------------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"; 
  if(debugMode)std::cout<<": "<<inFile_MC_name; 
  std::cout<<std::endl<<std::endl;
  TFile *fpp_MC = TFile::Open(inFile_MC.c_str());

  //mat
  TH2F *hmat,*hmat_anabin;
  hmat = (TH2F*)fpp_MC->Get(Form("hpp_mcclosure_matrix_R%d_20_eta_20",radius));
  hmat->Print("base");
  hmat_anabin = (TH2F*)hmat->Clone(Form("hpp_anaBin_Trans_matrix_HLT_R%d_20_eta_20",radius));
  hmat_anabin->Rebin2D(10, 10);
  hmat_anabin->Print("base");    
  
  //gen
  TH1F*hgen,*hgen_anabin;
  hgen= (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_gen_R%d_20_eta_20",radius));
  hgen->Print("base");
  hgen_anabin = (TH1F*)hgen->Clone(Form("hpp_anaBin_mcclosure_gen_JetComb_R%d_20_eta_20",radius));
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins, Form("hpp_anaBin_mcclosure_gen_JetComb_R%d_20_eta_20",radius), ptbins);
  //hgen_anabin = (TH1F*)hgen_anabin->Rebin(10);
  divideBinWidth(hgen_anabin);
  hgen_anabin->Print("base");

  //rec
  TH1F*hrec,*hrec_anabin;
  hrec = (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_data_R%d_20_eta_20",radius));
  hrec->Print("base");
  hrec_anabin = (TH1F*)hrec->Clone(Form("hpp_anaBin_mcclosure_rec_JetComb_R%d_20_eta_20",radius));
  hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins, Form("hpp_anaBin_mcclosure_rec_JetComb_R%d_20_eta_20",radius), ptbins);
  //hrec_anabin = (TH1F*)hrec_anabin->Rebin(10);
  divideBinWidth(hrec_anabin);
  hrec_anabin->Print("base");
  //if(dotrigcor) doTrigCorr(hrec_anabin, hMC_turnon);

  //rec_check
  TH1F *hrec_check, *hrec_check_anabin;
  hrec_check = (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_data_train_R%d_20_eta_20",radius));
  hrec_check->Print("base");
  //if(dotrigcor) doTrigCorr(hrec_check, hMC_turnon);
  hrec_check_anabin = (TH1F*)hrec_check->Clone(Form("hpp_mcclosure_JetComb_data_train_R%d_20_eta_20",radius));
  hrec_check_anabin = (TH1F*)hrec_check_anabin->Rebin(nbins, Form("hpp_mcclosure_JetComb_data_train_R%d_20_eta_20",radius), ptbins);
  //hrec_check_anabin = (TH1F*)hrec_check_anabin->Rebin(10);
  divideBinWidth(hrec_check_anabin);
  hrec_check_anabin->Print("base");

  // new histos --------------------------------------------------
  std::cout<<std::endl<<std::endl<<"creating NEW histos..."<<std::endl<<std::endl;

  if(debugMode)std::cout<<"opening  gen/genresponse histograms..."<<std::endl;
  TH1F *hgen_resp, *hgen_resp_anabin;
  hgen_resp = new TH1F( Form("hpp_gen_responseR%d_20_eta_20",radius),"", 
			hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  hgen_resp->Sumw2();
  hgen_resp->Print("base");
  hgen_resp_anabin = new TH1F(Form("hpp_gen_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hgen_resp_anabin->Sumw2();
  hgen_resp_anabin->Print("base");
  

  if(debugMode)std::cout<<"opening  rec/recresponse histograms..."<<std::endl;
  TH1F *hrec_resp, *hrec_resp_anabin;
  hrec_resp = new TH1F(Form("hpp_rec_response_R%d_20_eta_20",radius),"", 
		       hrec->GetNbinsX(), hrec->GetXaxis()->GetXmin(), hrec->GetXaxis()->GetXmax());
  hrec_resp->Sumw2();
  hrec_resp->Print("base");
  hrec_resp_anabin = new TH1F(Form("hpp_rec_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hrec_resp_anabin->Sumw2();
  hrec_resp_anabin->Print("base");
  
  // open output file for writing output to a root file as we go along
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE"); fout->cd();

  // start with input histos for convenient access later
  hgen->Write();  hgen_resp->Write();
  hrec->Write();  hrec_resp->Write();  hrec_check->Write();  

  // Bayesian unfolding -------------------------------------------------- 
  if(!doBayes) std::cout<<std::endl<<"   skipping Bayesian MC Closure test..."<<std::endl;
  else { std::cout<<std::endl<<std::endl<<"   beginning Bayesian unfolding..."<<std::endl;

    if(debugMode)std::cout<<std::endl<<"errorTreatment: "<<errorTreatment<<std::endl;

    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldResponse..."<<std::endl;
    RooUnfoldResponse roo_resp( hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_R%d",radius));

    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );

    if(debugMode)std::cout<<std::endl<<"calling Hreco..."<<std::endl<<std::endl;
    TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);
    hunf->SetName("PP_Bayesian_Unfolded_Spectra");
    TH1F *hratio = (TH1F*)hunf->Clone(Form("MCClosure_test_oppside_Bayes_R%d_20_eta_20", radius));
    hratio->Print("base");

    if(debugMode)std::cout<<std::endl<<"dividing hunf and hgen_anabin..."<<std::endl;
    //hgen_anabin->Print("base");
    hratio->Divide(hgen_anabin);

    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldBayes for MC Closure..."<<std::endl;
    RooUnfoldBayes unf_bayes_check(&roo_resp, hrec_check_anabin, kIter);
    
    if(debugMode)std::cout<<std::endl<<"calling Hreco..."<<std::endl<<std::endl;
    TH1F *hunf_check = (TH1F*)unf_bayes_check.Hreco(errorTreatment);
    hunf_check->Print("base");
    TH1F *hratio_check = (TH1F*)hunf_check->Clone(Form("MCClosure_test_sameside_Bayes_R%d_20_eta_20", radius));
    hratio_check->Print("base");

    if(debugMode)std::cout<<std::endl<<"dividing hunf_check and hgen_anabin..."<<std::endl;
    hratio_check->Divide(hgen_anabin);
    hratio_check->SetMarkerStyle(24);
    hratio_check->SetMarkerColor(kBlack);
    
    std::cout<<std::endl<<"writing bayesian unfolding output to file..."<<std::endl;
    fout->cd();
    hmat->Write();
    hunf->Write();
    hunf_check->Write();
    hratio->SetMarkerStyle(24);    hratio->SetMarkerColor(kRed);    hratio->Write();
    hratio_check->Write();
    if(drawPDFs){
      std::cout<<std::endl<<"drawing PDFs for Bayesian Unfolding..."<<std::endl;
      std::string outPdfFile=baseName+"_BayesianUnfoldingPlots.pdf";
      std::string open_outPdfFile=outPdfFile+"[";
      std::string close_outPdfFile=outPdfFile+"]";
      TCanvas* tempCanvForPdfPrint=new TCanvas("","",800,800);
      tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); tempCanvForPdfPrint->cd();
      hmat->Draw(); tempCanvForPdfPrint->Print(outPdfFile.c_str());
      hunf->Draw(); tempCanvForPdfPrint->Print(outPdfFile.c_str());
      hunf_check->Draw(); tempCanvForPdfPrint->Print(outPdfFile.c_str());
      hunf_check->Draw(); tempCanvForPdfPrint->Print(outPdfFile.c_str());
      hratio->Draw(); tempCanvForPdfPrint->Print(outPdfFile.c_str());
      hratio_check->Draw(); tempCanvForPdfPrint->Print(outPdfFile.c_str());
      tempCanvForPdfPrint->Print(close_outPdfFile.c_str());
    }
    std::cout<<"Bayesian MC Closure test done and all output written!"<<std::endl;
  } // end bayesian specifics


  // SVD unfolding --------------------------------------------------
  if(!doSVD) std::cout<<std::endl<<"   skipping SVD unfolding..."<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<"   beginning SVD unfolding..."<<std::endl;
    
    std::cout<<std::endl<<"initializing kReg parameter array w/ "<<nKregMax<<" elements"<<std::endl;
    int kReg[nKregMax];  
    std::cout<<kRegCenter-kRegRange<<" <= kReg[i] <= "<<kRegCenter+kRegRange<<std::endl;
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) kReg[i+kRegRange]=kRegCenter+i;
    if(debugMode)for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) 
      std::cout<<"kReg["<<i+kRegRange<<"]="<<kReg[i+kRegRange]<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating histo-arrays..."<<std::endl;
    TH2 *hPearsonSVDPriorMeas[nKregMax]; TH2 *hPearsonSVDPriorTruth[nKregMax];
    TH1 *hFoldedSVDPriorMeas[nKregMax]; TH1 *hFoldedSVDPriorTruth[nKregMax];
    TH1 *hrec_folded_ratio[nKregMax]; TH1 *hrec_unfolded_ratio[nKregMax];
    TH1F *hunf_svd[nKregMax];    TH1F *hunf_svd_check[nKregMax];
    TH1F *hratio_svd[nKregMax];  TH1F *hratio_svd_check[nKregMax];

    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldResponse..."<<std::endl;
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_anabin_R%d",radius));

    if(debugMode)std::cout<<std::endl<<"opening TCanvases..."<<std::endl;
    TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","", 1200, 1000);  cPearsonMatrixIter->Divide(3,3);
    TCanvas *cRatio             = new TCanvas("cRatio","",             1200, 1000);  cRatio->Divide(3,3);      
    TCanvas *cSpectra           = new TCanvas("cSpectra","",           1200, 1000);  cSpectra->Divide(3,3);
    TCanvas *cSpectraCheck      = new TCanvas("cSpectraCheck","",      1200, 1000);    

    TCanvas *c11 = new TCanvas("c11","c11: Singular Values", 1200, 1000);
    TCanvas *cdi = new TCanvas("cdi","cdi: di vectors",      1200, 1000);
    TLegend *leg[nKregMax],*leg1[nKregMax];      

    std::cout<<std::endl<<"entering SVD Unfolding Loop..."<<std::endl;
    for(int kr = 0; kr < nKregMax; ++kr){

      std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl<<std::endl;
  
      if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd"<<std::endl;
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
  
      if(debugMode)std::cout<<std::endl<<"calling Hreco..."<<std::endl<<std::endl;
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      hunf_svd[kr]->SetName(Form("PP_SVD_Unfolded_mcclosure_Spectra_kReg%d_20_eta_20",kReg[kr]));
      //hunf_svd[kr]->SetName(Form("h_Unfolded_mcclosure_spectra_SVD_kReg%d_R%d_20_eta_20", kReg[kr],radius));
      hunf_svd[kr]->Print("base");
  
      hratio_svd[kr] = (TH1F*)hunf_svd[kr]->Clone(Form("MCClosure_test_oppside_SVD_kReg%d_R%d_20_eta_20", kReg[kr], radius));
      hratio_svd[kr]->Divide(hgen_anabin);
      hratio_svd[kr]->SetMarkerStyle(33);
      hratio_svd[kr]->SetMarkerColor(kRed);
  
      // get covariance matrix and calculate pearson coefficients
      if(debugMode)std::cout<<std::endl<<"calling Ereco"<<std::endl;
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);

      if(debugMode)std::cout<<std::endl<<"calling CalculatePearsonCoefficients..."<<std::endl<<std::endl;
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat,debugPearson);

      if(debugMode)std::cout<<std::endl<<"writing pearson matrix as TH2D..."<<std::endl;
      hPearsonSVDPriorMeas[kr] = new TH2D(*pearson);
      //      hPearsonSVDPriorMeas[kr]->SetTitle(" ");
      hPearsonSVDPriorMeas[kr]->SetName(Form("hPearsonSVDPriorMeas_kReg%d", kReg[kr]));

      // draw Pearson Coefficient Matrix
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter canvas..."<<std::endl;
      cPearsonMatrixIter->cd(kr+1);  
      hPearsonSVDPriorMeas[kr]->SetTitle( ( "kReg="+std::to_string(kReg[kr]) ).c_str() );
      hPearsonSVDPriorMeas[kr]->SetMinimum(-1.);
      hPearsonSVDPriorMeas[kr]->SetMaximum(1.);
      hPearsonSVDPriorMeas[kr]->SetAxisRange(0, 35, "X");
      hPearsonSVDPriorMeas[kr]->SetAxisRange(0, 35, "Y");
      hPearsonSVDPriorMeas[kr]->GetZaxis()->SetLabelSize(0.06);
      hPearsonSVDPriorMeas[kr]->Draw("colz");

      // Apply to Truth
      if(debugMode)std::cout<<std::endl<<"applying roo_resp to truth dist hunf_svd[kr="<<kr<<"]"<<std::endl;
      hFoldedSVDPriorMeas[kr] = roo_resp.ApplyToTruth(hunf_svd[kr]);
      hFoldedSVDPriorMeas[kr]->SetName(Form("hFoldedSVDPriorMeas_kReg%d", kReg[kr]));  

      // draw on cSpectra canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl;
      cSpectra->cd(kr+1);
      cSpectra->cd(kr+1)->SetLogy();


      std::cout <<"CHECK: kr="<<kr<<"  and kReg[kr]="<<kReg[kr]<<std::endl<<std::endl;
      TH1F* hrec_anabin_clone=(TH1F*)hrec_anabin->Clone("untitledClone");
      hrec_anabin_clone->SetTitle( ("kReg = "+std::to_string(kReg[kr])).c_str() );//Form("kReg = %d",kReg[kr]) );
      hrec_anabin_clone->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_anabin_clone->SetMarkerStyle(24);
      hrec_anabin_clone->SetMarkerColor(kBlack);
      hrec_anabin_clone->SetAxisRange(45, 1000, "X");
      hrec_anabin_clone->Print("base");
      hrec_anabin_clone->Draw();
  
      hunf_svd[kr]->SetMarkerStyle(33);
      hunf_svd[kr]->SetMarkerColor(kBlue);
      hunf_svd[kr]->Draw("same");
  
      hFoldedSVDPriorMeas[kr]->SetLineStyle(33);
      hFoldedSVDPriorMeas[kr]->SetLineColor(kRed);
      hFoldedSVDPriorMeas[kr]->Draw("same");
  
      leg[kr] = new TLegend(0.1, 0.1, 0.30, 0.30, NULL,"BRNDC");
      leg[kr]->AddEntry(hrec_anabin_clone,"Measured","pl");
      leg[kr]->AddEntry(hunf_svd[kr],Form("Unfolded, kr = %d", kReg[kr]),"pl");
      leg[kr]->AddEntry(hFoldedSVDPriorMeas[kr],Form("Folded kr = %d", kReg[kr]),"pl");
      leg[kr]->SetTextSize(0.04); 
      leg[kr]->Draw();
      	
      // draw on cRatio canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio canvas..."<<std::endl;
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
      leg1[kr]->Draw("same");//leg1[kr]->Draw();
  	
      //Get and print singular values and d_i vector, print also on screen
      //Note that these do not depend on the regularization. 
      //The opposite: they tell you which regularization to use!
      if(kr == 0){
	std::cout<<std::endl<<"kr==0, getting singular values and di vectors!"<<std::endl<<std::endl;
	
  	TSVDUnfold *svdUnfold = unf_svd.Impl();

  	TH1 *hSVal = (TH1*)svdUnfold->GetSV();
  	if(debugMode)std::cout << "  printing singular values: " << std::endl;
  	if(debugMode)for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
		       std::cout<<"bin: "<<bin<<",  SV: "<<hSVal->GetBinContent(bin)<< std::endl;
  	TH1D *hdi = (TH1D*)svdUnfold->GetD();
  	if(debugMode)std::cout << std::cout<<"  printing di vector: " <<  std::endl;
  	if(debugMode)for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
		       std::cout<<"i: "<<bin<<",  di: "<<hdi->GetBinContent(bin)<<std::endl;

	if(debugMode)std::cout<<std::endl<<"drawing stuff on c11 canvas..."<<std::endl;
  	c11->cd();
	gPad->SetLogy();	//c11->SetLogy();  	
	hSVal->SetTitle(" ");
  	hSVal->SetAxisRange(0,35,"X");
  	hSVal->SetXTitle("singular values");
  	hSVal->DrawCopy();  	
  	drawText("PP 5 TeV", 0.608173, 0.8459761, 22);
	
	if(debugMode)std::cout<<std::endl<<"drawing stuff on cdi canvas..."<<std::endl;
  	cdi->cd();
	gPad->SetLogy();	//cdi->SetLogy();
  	hdi->SetTitle(" ");
  	hdi->SetAxisRange(0,35,"X");
  	hdi->SetXTitle("|d_{i}^{kreg}|");
  	hdi->DrawCopy();
  	
  	drawText("PP 5 TeV", 0.608173, 0.8459761, 22);
      } // end kr==0 specific
  
      hrec_unfolded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Print("base");

      // the closure tests for svd unfolding
      if(debugMode)std::cout<<"calling RooUnfoldSVD for MCClosure tests"<<std::endl;
      RooUnfoldSvd unf_svd_check(&roo_resp, hrec_check_anabin, kReg[kr]);
      
      if(debugMode)std::cout<<std::endl<<"calling Hreco..."<<std::endl<<std::endl;
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
      if(debugMode)std::cout<<std::endl<<"calling Ereco..."<<std::endl;
      TMatrixD covmat_truth = unf_svd_check.Ereco(errorTreatment);

      if(debugMode)std::cout<<std::endl<<"calling CalculatePearsonCoefficients..."<<std::endl;
      std::cout<<std::endl;
      TMatrixD *pearson_truth = CalculatePearsonCoefficients(&covmat_truth,debugPearson);      

      if(debugMode)std::cout<<"writing pearson matrix as TH2D..."<<std::endl;
      hPearsonSVDPriorTruth[kr] = new TH2D(*pearson_truth);
      hPearsonSVDPriorTruth[kr]->SetName(Form("hPearsonSVDPriorTruth_kReg%d",kReg[kr]));
      hPearsonSVDPriorTruth[kr]->SetMinimum(-1.);
      hPearsonSVDPriorTruth[kr]->SetMaximum(1.);
      hPearsonSVDPriorTruth[kr]->GetZaxis()->SetLabelSize(0.06);
      
      if(debugMode)std::cout<<std::endl<<"applying roo_resp to truth dist hunf_svd_check[kr="<<kr<<"]"<<std::endl;
      hFoldedSVDPriorTruth[kr] = roo_resp.ApplyToTruth(hunf_svd_check[kr]);
      hFoldedSVDPriorTruth[kr]->SetName(Form("hFoldedSVDPriorTruth_kReg%d", kReg[kr]));
  
      TSVDUnfold *svdUnfold = unf_svd_check.Impl();      
      TH1D *hdi = svdUnfold->GetD();
      TH1 *hSVal = svdUnfold->GetSV();
      //TH1D  *hSValClone = (TH1D*)hSVal->Clone("hSValClone");
      
      if(debugMode)std::cout << "  printing singular values: " << std::endl;
      if(debugMode)for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)//did start at 1, now starts at 0?
	std::cout << "bin: " << bin << "  SV: " << hSVal->GetBinContent(bin) << std::endl;
      if(debugMode)std::cout << std::endl<<"  printing di vector: " <<  std::endl;
      if(debugMode)for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
	std::cout << "i: " << bin << "  di: " << hdi->GetBinContent(bin) << std::endl;
    }// end kReg loop

    // loop over histos pointers in arrays to write
    if(debugMode)std::cout<<std::endl<<"writing histo arrays to file..."<<std::endl;
    fout->cd();
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
    if(drawPDFs){ std::cout<<std::endl<<"drawing PDFs for SVD Unfolding"<<std::endl;
      //std::string type="MCClosure";//should only be Data or MC

      std::string outPdfFile=baseName+"_SVDUnfoldingPlots";
      if(doBayes)outPdfFile+="_wBayes";
      outPdfFile+=".pdf";

      std::string open_outPdfFile=outPdfFile+"[";
      std::string close_outPdfFile=outPdfFile+"]";
      cPearsonMatrixIter->Print(open_outPdfFile.c_str()); 
      cPearsonMatrixIter->cd();   cPearsonMatrixIter->Print(outPdfFile.c_str()); 
      cSpectra->cd();   cSpectra->Print(outPdfFile.c_str());
      cRatio->cd();   cRatio->Print(outPdfFile.c_str()); 
      cdi->cd();   cdi->Print(outPdfFile.c_str()); 
      c11->cd();   c11->Print(outPdfFile.c_str()); 

      //cPearsonMatrixIter->SaveAs(Form("Pearson_Matrix_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");
      //cSpectra->SaveAs(Form("Spectra_meas_unf_fol_%s_R%d_plots.pdf",type.c_str(), radius),"RECREATE");
      //cRatio->SaveAs(Form("Ratio_meas_unf_fol_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");          
      //cdi->SaveAs(Form("di_Vectors_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
      //c11->SaveAs(Form("Singular_Values_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
      
      cSpectraCheck->cd();
      
      hrec_folded_ratio[kRegDraw]->SetTitle(" ");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(0.4, 1.2, "Y");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(45, 1000, "X");
      
      hrec_folded_ratio[kRegDraw]->Draw();
      hrec_unfolded_ratio[kRegDraw]->Draw("same");
      drawText("PP 5 TeV", 0.608173, 0.8459761, 22);
      drawText(Form("kReg = %d",kReg[kRegDraw]), 0.60, 0.75, 22);

      TLegend * leg2 = new TLegend(0.1, 0.1, 0.30, 0.3, NULL,"BRNDC");
      leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Unfolded/Measured","pl");
      leg2->AddEntry(hrec_folded_ratio[kRegDraw],"Folded/Measured","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2->SetTextSize(0.04); 
      leg2->Draw();    
      
      //cSpectraCheck->SaveAs(Form("Spectra_Ratio_withkRegDraw_%s_R%d.pdf", type.c_str(), radius),"RECREATE");	
      cSpectraCheck->Print(outPdfFile.c_str());
      
      // MCClosure Plot
      //std::string c1_filename="PP_5p02TeV_"+type+"_R"+std::to_string(radius)+"_SVDOnly.pdf";
      TCanvas *c1 = new TCanvas("c1","Spectra",1300,1000);  c1->cd();
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

      if(doBayes){ // if we did bayes as well, put it on the same MCClosure plot as the SVD
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
      
      drawText("PP, 5 TeV", 0.648173, 0.8459761,22);
      TLine *line = new TLine(30,1,1000,1);
      line->SetLineStyle(2); 
      line->SetLineWidth(2);
      line->Draw();
      
      drawText(Form("kReg = %d", kReg[kRegDraw]), 0.64, 0.7, 21);
      drawText("PYTHIA",0.22,0.88,21);
      
      leg0->AddEntry(hratio_svd[kRegDraw],"OppSide SVD","p");
      leg0->AddEntry(hratio_svd_check[kRegDraw],"SameSide SVD","p");
      leg0->Draw();
      
      c1->Print(outPdfFile.c_str());
      c1->Print(close_outPdfFile.c_str());
      //c1->SaveAs( c1_filename.c_str(), "RECREATE");
    }// end drawPDFs
  }// end SVD specific

  std::cout<<"writing histos to output file"<<std::endl;
  fout->Write();
  fout->Close();
  fpp_MC->Close();  

  return 0;
} // end doMCClosureTests




//-----------------------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]){
  int rStatus = -1;
  if( argc!=1 && argc!=2 ) {//no input arguments, error
    std::cout<<"do ./doMCClosureTests.exe <debugMode,def0> to run"<<std::endl;
    return rStatus;
  }

  rStatus=1;//1 until a function returns 0 (normal exit by my convention)
  if(argc==1)   rStatus=doMCClosureTests();
  if(argc==2)   rStatus=doMCClosureTests( (const bool) atoi(argv[1]) ); 
  std::cout<<"done! return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}











  //// ppData input histos --------------------------------------------------
  //std::cout<<std::endl<<std::endl<<"opening INPUT histos from DATA file"; 
  //if(debugMode)std::cout<<": "<<inFile_Data_name; 
  //std::cout<<std::endl<<std::endl;
  //TFile *fpp_Data = TFile::Open(inFile_Data.c_str());

      //if(debugMode)std::cout<<std::endl<<""<<std::endl;
