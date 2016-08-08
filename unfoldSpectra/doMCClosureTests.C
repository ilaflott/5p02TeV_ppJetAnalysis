#include "unfoldSpectra.h"


//  global settings -------------------------------------------------- 
const int radius=4;
const bool doToyErrs=false;

// Bayesian settings 
const bool doBayes=true; 
const int kIter = 4;

// SVD settings 
const bool doSVD    = true; // testing now
const int  nKregMax = 9, kRegRange=4, kRegDraw = 3, kRegCenter=21;

// I/O
const std::string baseName="doMCClosureTests_defOut"; // output name
const std::string outRootFile=baseName+".root";
const std::string outBayesPdfFile=baseName+"_BayesianUnfoldingPlots.pdf";
const std::string outSVDPdfFile=baseName+"_SVDUnfoldingPlots"; // see drawPDFs part for rest of string

const bool drawPDFs = true; 
const bool drawPDFs_BayesInputHistos= true&&drawPDFs;
const bool drawPDFs_SVDInputHistos  = true&&drawPDFs;
const bool defDebugMode=false;
bool debugPearson=false;   // for easy debugging of pearson matrix stuff

//  the code --------------------------------------------------
int doMCClosureTests( const bool debugMode=defDebugMode){

  // other settings that are easily forgot
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl; 
  debugPearson=(bool)(debugPearson && debugMode);
  if(debugMode)std::cout<<"debugPearson="<<debugPearson<<std::endl;

  // script performance tracking if needed
  TStopwatch timer; 
  timer.Start();
  
  // set error handing, stat info, other settings  // fix me
  RooUnfold::ErrorTreatment errorTreatment;
  if(!doToyErrs) errorTreatment = RooUnfold::kCovariance;
  else errorTreatment = RooUnfold::kCovToy; 
  if(debugMode)std::cout<<"doToyErrs="<<doToyErrs<<std::endl; 

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);


  //  ppMC input histos --------------------------------------------------
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

  
  //  new histos --------------------------------------------------
  std::cout<<std::endl<<std::endl<<"creating NEW histos..."<<std::endl<<std::endl;

  TH1F *hgen_resp, *hgen_resp_anabin;
  hgen_resp = new TH1F( Form("hpp_gen_response_R%d_20_eta_20",radius),"", 
			hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  hgen_resp->Sumw2();
  hgen_resp->Print("base");
  hgen_resp_anabin = new TH1F(Form("hpp_gen_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hgen_resp_anabin->Sumw2();
  hgen_resp_anabin->Print("base");
  
  TH1F *hrec_resp, *hrec_resp_anabin;
  hrec_resp = new TH1F(Form("hpp_rec_response_R%d_20_eta_20",radius),"", 
		       hrec->GetNbinsX(), hrec->GetXaxis()->GetXmin(), hrec->GetXaxis()->GetXmax());
  hrec_resp->Sumw2();
  hrec_resp->Print("base");
  hrec_resp_anabin = new TH1F(Form("hpp_rec_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hrec_resp_anabin->Sumw2();
  hrec_resp_anabin->Print("base");
  
  // open output file for writing output to a root file as we go along
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE"); 

  //  Bayesian unfolding -------------------------------------------------- 
  if(!doBayes) std::cout<<std::endl<<"   skipping Bayesian MC Closure test..."<<std::endl<<std::endl;
  else { std::cout<<std::endl<<std::endl<<std::endl<<"   beginning Bayesian MC Closure Test..."<<std::endl;

    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldResponse..."<<std::endl;
    RooUnfoldResponse roo_resp( hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_R%d",radius));

    if(debugMode)std::cout<<"calling RooUnfoldBayes and Hreco..."<<std::endl<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );

    TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);
    hunf->SetName("PP_Bayesian_Unfolded_Spectra");
    std::cout<<std::endl; hunf->Print("base");

    TH1F *hratio = (TH1F*)hunf->Clone(Form("MCClosure_test_oppside_Bayes_R%d_20_eta_20", radius));
    //hratio->SetName("PP_Bayesian_Unfolded_Ratio");
    hratio->Print("base");
    hratio->Divide(hgen_anabin);

    std::cout<<std::endl<<"doing unfolding checks..."<<std::endl;
    if(debugMode)std::cout<<"calling RooUnfoldBayes and Hreco for MC Closure..."<<std::endl; 
    
    std::cout<<std::endl;    
    RooUnfoldBayes unf_bayes_check(&roo_resp, hrec_check_anabin, kIter);
    TH1F *hunf_check = (TH1F*)unf_bayes_check.Hreco(errorTreatment);
    hunf_check->SetName("PP_Bayesian_Unfolded_Spectra_Check");
    std::cout<<std::endl;        hunf_check->Print("base");

    TH1F *hratio_check = (TH1F*)hunf_check->Clone(Form("MCClosure_test_sameside_Bayes_R%d_20_eta_20", radius));
    //hratio_check->SetName("PP_Bayesian_Unfolded_Ratio_Check");
    hratio_check->Print("base");
    hratio_check->Divide(hgen_anabin);
    
    std::cout<<std::endl<<"writing bayesian unfolding output to file..."<<std::endl;
    fout->cd();
    hunf->Write();    hunf_check->Write();
    hgen_resp->Write();  hgen_resp_anabin->Write();
    hrec_resp->Write();	 hrec_resp_anabin->Write();

    hratio->SetMarkerStyle(24);    hratio_check->SetMarkerStyle(24); 
    hratio->SetMarkerColor(kRed);  hratio_check->SetMarkerColor(kBlack);    
    hratio->Write();               hratio_check->Write();

    std::cout<<"writing input histos to output file for easy access later..."<<std::endl;
    hmat->Write(); hmat_anabin->Write(); 
    hgen->Write(); hgen_anabin->Write(); 
    hrec->Write(); hrec_anabin->Write(); 
    hrec_check->Write(); hrec_check_anabin->Write();

    //  drawPDFS -------------------------------------------------- 
    if(!drawPDFs)std::cout<<std::endl<<"NOT drawing PDFs for Bayesian Unfolding!"<<std::endl<<std::endl;
    else{ std::cout<<std::endl<<"drawing output PDFs for Bayesian Unfolding..."<<std::endl;
      std::string outPdfFile=outBayesPdfFile;
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";

      // temp canvas for printing
      TCanvas* tempCanvForPdfPrint=new TCanvas("","",800,800);    
      tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); 
      tempCanvForPdfPrint->cd();

      // draw and print unfold/ratio plots
      hunf->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
      hratio->Draw();       tempCanvForPdfPrint->Print(outPdfFile.c_str());

      hunf_check->Draw();   tempCanvForPdfPrint->Print(outPdfFile.c_str());
      hratio_check->Draw(); tempCanvForPdfPrint->Print(outPdfFile.c_str());

      // draw and print responses
      //hgen_resp->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
      //hrec_resp->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
      //hgen_resp_anabin->DrawCopy();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
      //hrec_resp_anabin->DrawCopy();         tempCanvForPdfPrint->Print(outPdfFile.c_str());

      if(drawPDFs_BayesInputHistos){
	std::cout<<std::endl<<"drawing input histos to Bayesian Unfolding..."<<std::endl;
	hmat->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
	//hmat_anabin->Draw();   tempCanvForPdfPrint->Print(outPdfFile.c_str());
	hgen->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
	//hgen_anabin->Draw();   tempCanvForPdfPrint->Print(outPdfFile.c_str());
	hrec->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
	//hrec_anabin->Draw();   tempCanvForPdfPrint->Print(outPdfFile.c_str());
	hrec_check->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());
	//hrec_check_anabin->Draw();   tempCanvForPdfPrint->Print(outPdfFile.c_str());   
      }
      tempCanvForPdfPrint->Print(close_outPdfFile.c_str());      
    }// end draw pdfs
  }  // end doBayes
  std::cout<<std::endl<<"Bayesian MC Closure test done and all output written!"<<std::endl;
  
  
  //  SVD unfolding --------------------------------------------------
  if(!doSVD) std::cout<<std::endl<<"   skipping SVD MCClosure Test..."<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<std::endl<<"   beginning SVD MCClosure Test..."<<std::endl;
    
    int kReg[nKregMax];  
    std::cout<<std::endl<<"initializing kReg parameter array w/ "<<nKregMax<<" elements"<<std::endl;
    std::cout<<kRegCenter-kRegRange<<" <= kReg[i] <= "<<kRegCenter+kRegRange<<std::endl;
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) kReg[i+kRegRange]=kRegCenter+i;
    if(debugMode)for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) 
      std::cout<<"kReg["<<i+kRegRange<<"]="<<kReg[i+kRegRange]<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating histo-arrays..."<<std::endl;
    TH1 *hrec_folded_ratio[nKregMax]; TH1 *hrec_unfolded_ratio[nKregMax];

    TH2 *hPearsonSVDPriorMeas[nKregMax]; TH2 *hPearsonSVDPriorTruth[nKregMax];
    TH1 *hFoldedSVDPriorMeas[nKregMax]; TH1 *hFoldedSVDPriorTruth[nKregMax];
    TH1F *hunf_svd[nKregMax];    TH1F *hunf_svd_check[nKregMax];
    TH1F *hratio_svd[nKregMax];  TH1F *hratio_svd_check[nKregMax];

    if(debugMode)std::cout<<"calling RooUnfoldResponse..."<<std::endl;
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_anabin_R%d",radius));

    if(debugMode)std::cout<<"opening TCanvases..."<<std::endl;
    TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","", 1200, 1000);  cPearsonMatrixIter->Divide(3,3);
    TCanvas *cRatio             = new TCanvas("cRatio","",             1200, 1000);  cRatio->Divide(3,3);      
    TCanvas *cSpectra           = new TCanvas("cSpectra","",           1200, 1000);  cSpectra->Divide(3,3);
    TCanvas *cRatioCheck      = new TCanvas("cRatioCheck","",      1200, 1000);    
    TCanvas *c11 = new TCanvas("c11"," Singular Values and divectors", 1200, 1000); c11->Divide(2);
    //TCanvas *cdi = new TCanvas("cdi","cdi: di vectors",      1200, 1000);

    TLegend *leg[nKregMax],*leg1[nKregMax];      

    std::cout<<"entering SVD Unfolding Loop..."<<std::endl;
    for(int kr = 0; kr < nKregMax; ++kr){
      std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl;
      
      // unfold the MC dist, then create a ratio with the gen-level MC dist
      if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd and Hreco..."<<std::endl<<std::endl;
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      
      std::cout<<std::endl;   
      hunf_svd[kr]->SetName(Form("PP_SVD_Unfolded_mcclosure_Spectra_kReg%d_20_eta_20",kReg[kr]));
      hunf_svd[kr]->SetLineStyle(33);
      hunf_svd[kr]->SetLineColor(kBlue);
      hunf_svd[kr]->Print("base");
  
      hratio_svd[kr] = (TH1F*)hunf_svd[kr]->Clone(Form("MCClosure_test_oppside_SVD_kReg%d_R%d_20_eta_20", kReg[kr], radius));
      hratio_svd[kr]->Divide(hgen_anabin);
      hratio_svd[kr]->SetMarkerStyle(33);
      hratio_svd[kr]->SetMarkerColor(kRed);
      hratio_svd[kr]->Print("base");

      // get covariance matrix and calculate pearson coefficients
      if(debugMode)std::cout<<std::endl<<"calling Ereco..."<<std::endl;
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);

      if(debugPearson)std::cout<<std::endl;
      if(debugMode)std::cout<<"calling CalculatePearsonCoefficients..."<<std::endl;
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat, debugPearson );
      if(debugPearson)std::cout<<std::endl;
      
      // drawing/visuals
      if(debugMode)std::cout<<"creating TH2D for pearson matrix..."<<std::endl;
      hPearsonSVDPriorMeas[kr] = new TH2D(*pearson);
      hPearsonSVDPriorMeas[kr]->SetName(Form("hPearsonSVDPriorMeas_kReg%d", kReg[kr]));
      hPearsonSVDPriorMeas[kr]->SetTitle( ( "pearson matrix, kReg="+std::to_string(kReg[kr]) ).c_str() ); 
      
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter canvas..."<<std::endl;
      cPearsonMatrixIter->cd(kr+1);  
      
      float minNmax=0.3; 
      hPearsonSVDPriorMeas[kr]->SetMinimum(-1*minNmax);  hPearsonSVDPriorMeas[kr]->SetMaximum(minNmax);
      hPearsonSVDPriorMeas[kr]->SetAxisRange(0, 35, "X");
      hPearsonSVDPriorMeas[kr]->SetAxisRange(0, 35, "Y");
      hPearsonSVDPriorMeas[kr]->GetZaxis()->SetLabelSize(0.035);
      hPearsonSVDPriorMeas[kr]->Draw("colz");

      if(debugMode)std::cout<<"applying roo_resp to \"truth\" histo hunf_svd[kr="<<kr<<"]..."<<std::endl;
      hFoldedSVDPriorMeas[kr] = roo_resp.ApplyToTruth(hunf_svd[kr]);
      hFoldedSVDPriorMeas[kr]->SetName(Form("hFoldedSVDPriorMeas_kReg%d", kReg[kr]));  
      hFoldedSVDPriorMeas[kr]->SetLineStyle(33);
      hFoldedSVDPriorMeas[kr]->SetLineColor(kRed);

      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
      cSpectra->cd(kr+1);
      cSpectra->cd(kr+1)->SetLogy();

      // what's happening here exactly?
      //std::cout <<"CHECK: kr="<<kr<<"  and kReg[kr]="<<kReg[kr]<<std::endl<<std::endl;
      TH1F* hrec_anabin_clone=(TH1F*)hrec_anabin->Clone("hrec_anabin_clone");
      hrec_anabin_clone->SetTitle( ("spectra, kReg = "+std::to_string(kReg[kr])).c_str() );//Form("kReg = %d",kReg[kr]) );
      hrec_anabin_clone->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_anabin_clone->SetMarkerStyle(24);
      hrec_anabin_clone->SetMarkerColor(kBlack);
      hrec_anabin_clone->SetAxisRange(45, 1000, "X");
      hrec_anabin_clone->Print("base");
      hrec_anabin_clone->Draw();

      hunf_svd[kr]->Draw("same");
      hFoldedSVDPriorMeas[kr]->Draw("same");
  
      leg[kr] = new TLegend(0.1, 0.1, 0.40, 0.30, NULL,"NBNDC");//x1,y1,x2,y2,header,option
      leg[kr]->AddEntry(hrec_anabin_clone,"Measured","pl");
      leg[kr]->AddEntry(hunf_svd[kr],"Unfolded","pl");
      leg[kr]->AddEntry(hFoldedSVDPriorMeas[kr],"Folded","pl");
      leg[kr]->SetTextSize(0.04); 
      leg[kr]->Draw();
      	
      // draw on cRatio canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio canvas..."<<std::endl;
      cRatio->cd(kr+1);
  
      hrec_folded_ratio[kr] = (TH1F*)hFoldedSVDPriorMeas[kr]->Clone(Form("ratio_folded_with_measured_kReg%d",kReg[kr]));
      hrec_folded_ratio[kr]->SetTitle(Form("ratio, kReg = %d",kReg[kr]));
      hrec_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_folded_ratio[kr]->SetMarkerColor(kRed);
      hrec_folded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_folded_ratio[kr]->SetAxisRange(45, 1000, "X");
      hrec_folded_ratio[kr]->SetAxisRange(0.5, 1.5, "Y");
      hrec_folded_ratio[kr]->Divide(hrec_anabin);
      hrec_folded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Draw();
  
      hrec_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone(Form("ratio_unfolded_with_measured_kReg%d",kReg[kr]));
      hrec_unfolded_ratio[kr]->SetTitle(Form("div. by measured, kReg = %d",kReg[kr]));
      hrec_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_unfolded_ratio[kr]->SetMarkerColor(kBlue);
      hrec_unfolded_ratio[kr]->Divide(hrec_anabin);
      hrec_unfolded_ratio[kr]->Print("base");
      hrec_unfolded_ratio[kr]->Draw("same");

      leg1[kr] = new TLegend(0.1, 0.1, 0.40, 0.3, NULL,"NBNDC");
      leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],"Unfold/Meas","pl");
      leg1[kr]->AddEntry(hrec_folded_ratio[kr],"Fold/Meas","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg1[kr]->SetTextSize(0.04); 
      leg1[kr]->Draw("same");//leg1[kr]->Draw();
  	

      //  singular values and d_i vector ---------------------------
      //Note that these do not depend on the regularization. 
      //The opposite: they tell you which regularization to use! (ian note: how?)
      if(kr == 0){   	TSVDUnfold *svdUnfold = unf_svd.Impl();
	if(debugMode) std::cout<<std::endl<<"    !!! kr==0 !!!"<<std::endl<<std::endl; 

	// singular values
  	if(debugMode)std::cout << "  getting singular values... " << std::endl;
  	TH1 *hSVal = (TH1*)svdUnfold->GetSV();
  	if(debugMode)for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
		       std::cout<<"bin: "<<bin<<",  SV: "<<hSVal->GetBinContent(bin)<< std::endl;

	if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
  	c11->cd(1);	gPad->SetLogy();	//c11->SetLogy();  	

	hSVal->SetTitle(" singular values ");
  	hSVal->SetXTitle(" singular values ");
  	hSVal->SetAxisRange(0,35,"X");
  	hSVal->DrawCopy();  	
	drawText( "5.02 TeV ppMC, QCD Py6 Tune Z2",0.358173, 0.8459761, 19);
	drawText( ("kReg="+std::to_string(kReg[kr])).c_str(),0.408173, 0.8059761, 19);

	// di vector values
  	if(debugMode)std::cout << "  getting di vector... " <<  std::endl;
  	TH1D *hdi = (TH1D*)svdUnfold->GetD();
  	if(debugMode)for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
		       std::cout<<"i: "<<bin<<",  di: "<<hdi->GetBinContent(bin)<<std::endl;

	if(debugMode)std::cout<<std::endl<<"drawing di vectors on c11 canvas.."<<std::endl;
	c11->cd(2);	gPad->SetLogy();	//cdi->SetLogy();

  	hdi->SetTitle(" di vectors ");
  	hdi->SetXTitle(" |d_{i}^{kreg}| ");
  	hdi->SetAxisRange(0,35,"X");
  	hdi->DrawCopy();
	drawText( "5.02 TeV ppMC, QCD Py6 Tune Z2",0.358173, 0.8459761, 19);
	drawText( ("kReg="+std::to_string(kReg[kr])).c_str(),0.408173, 0.8059761, 19);
	if(debugMode)std::cout<<std::endl<<"done with kr==0 specifics"<<std::endl<<std::endl;
      }
        
      // the closure tests for svd unfolding
      if(debugMode)std::cout<<"calling RooUnfoldSVD for MCClosure tests"<<std::endl;
      RooUnfoldSvd unf_svd_check(&roo_resp, hrec_check_anabin, kReg[kr]);
      
      if(debugMode)std::cout<<std::endl<<"calling Hreco..."<<std::endl;
      hunf_svd_check[kr] = (TH1F*)unf_svd_check.Hreco(errorTreatment);

      std::string hunf_svd_check_name="h_Unfolded_mcclosure_sameside_SVD";
      hunf_svd_check_name+="_kReg"+std::to_string(kReg[kr]);
      hunf_svd_check_name+="_R"+std::to_string(radius)+"_20_eta_20";
      hunf_svd_check[kr]->SetName(hunf_svd_check_name.c_str());
      std::cout<<std::endl;      hunf_svd_check[kr]->Print("base");

      std::string hratio_svd_check_name="MCClosure_test_sameside_SVD";
      hratio_svd_check_name+="_kReg"+std::to_string(kReg[kr]);
      hratio_svd_check_name+="_R"+std::to_string(radius)+"_20_eta_20";
      hratio_svd_check[kr] = (TH1F*)hunf_svd_check[kr]->Clone( hratio_svd_check_name.c_str());
      hratio_svd_check[kr]->SetMarkerStyle(33);
      hratio_svd_check[kr]->SetMarkerColor(kBlack);
      hratio_svd_check[kr]->Divide(hgen_anabin);
      hratio_svd_check[kr]->Print("base");
      
      // get covariance matrix and calculate pearson coefficients
      if(debugMode)std::cout<<std::endl<<"calling Ereco..."<<std::endl;
      TMatrixD covmat_truth = unf_svd_check.Ereco(errorTreatment);
      
      if(debugPearson)std::cout<<std::endl<<"pearsonDebug is ON! ";
      if(debugMode)std::cout<<"calling CalculatePearsonCoefficients..."<<std::endl;
      TMatrixD *pearson_truth = CalculatePearsonCoefficients(&covmat_truth,debugPearson);      

      if(debugMode)std::cout<<"writing pearson matrix as TH2D..."<<std::endl;
      hPearsonSVDPriorTruth[kr] = new TH2D(*pearson_truth);
      hPearsonSVDPriorTruth[kr]->SetName(Form("hPearsonSVDPriorTruth_kReg%d",kReg[kr]));
      hPearsonSVDPriorTruth[kr]->SetMinimum(-1.);
      hPearsonSVDPriorTruth[kr]->SetMaximum(1.);
      hPearsonSVDPriorTruth[kr]->GetZaxis()->SetLabelSize(0.06);
      
      if(debugMode)std::cout<<std::endl<<"applying roo_resp to \"truth\" histo, hunf_svd_check["<<kr<<"]"<<std::endl;
      hFoldedSVDPriorTruth[kr] = roo_resp.ApplyToTruth(hunf_svd_check[kr]);
      hFoldedSVDPriorTruth[kr]->SetName(Form("hFoldedSVDPriorTruth_kReg%d", kReg[kr]));
  
      TSVDUnfold *svdUnfold = unf_svd_check.Impl();      
      
      if(debugMode)std::cout << std::endl<<std::endl<<"  printing singular values: " << std::endl<<std::endl;
      TH1 *hSVal = svdUnfold->GetSV();
      if(debugMode)for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)//did start at 1, now starts at 0?
	std::cout << "bin: " << bin << "  SV: " << hSVal->GetBinContent(bin) << std::endl;

      if(debugMode)std::cout << std::endl<<std::endl<<"  printing di vector: " <<  std::endl<<std::endl;
      TH1D *hdi = svdUnfold->GetD();
      if(debugMode)for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
	std::cout << "i: " << bin << "  di: " << hdi->GetBinContent(bin) << std::endl;
    }// end kReg loop

    std::cout<<std::endl<<"writing SVD histo arrays to file..."<<std::endl;    fout->cd();
    for(int kr = 0; kr<nKregMax; ++kr){
      hunf_svd[kr]->Write();  hunf_svd_check[kr]->Write();
      hratio_svd[kr]->Write();  hratio_svd_check[kr]->Write();
      hPearsonSVDPriorMeas[kr]->Write();  hPearsonSVDPriorTruth[kr]->Write();
      hFoldedSVDPriorMeas[kr]->Write();  hFoldedSVDPriorTruth[kr]->Write();    
    }

    if(!drawPDFs){ if(debugMode)std::cout<<std::endl<<"NOT drawing PDFs for SVD Unfolding"<<std::endl<<std::endl; }
    else{ std::cout<<std::endl<<"drawing PDFs for SVD Unfolding..."<<std::endl;
      
      // form filename string
      std::string outPdfFile=outSVDPdfFile;      
      if(doBayes)outPdfFile+="_wBayes"; 
      outPdfFile+=".pdf";
      
      // open pdf file and draw the canvases we have so far
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
      cPearsonMatrixIter->Print(open_outPdfFile.c_str()); 


      //  already-drawn canvases --------------------------------------------------      
      cSpectra->cd()           ;   cSpectra->Print(outPdfFile.c_str());
      cRatio->cd()             ;   cRatio->Print(outPdfFile.c_str()); 
      cPearsonMatrixIter->cd() ;   cPearsonMatrixIter->Print(outPdfFile.c_str()); 
      c11->cd()                ;   c11->Print(outPdfFile.c_str()); 
      

      //  cRatio Check --------------------------------------------------
      cRatioCheck->cd();
      
      hrec_folded_ratio[kRegDraw]->SetTitle(" ");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(0.4, 1.2, "Y");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(45, 1000, "X");      
      hrec_folded_ratio[kRegDraw]->Draw();
      
      TLegend * leg2 = new TLegend(0.1, 0.1, 0.40, 0.3, NULL,"NBNDC");
      leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Unfold/Meas","pl");
      leg2->AddEntry(hrec_folded_ratio[kRegDraw],"Fold/Meas","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2->SetTextSize(0.04); 
      leg2->Draw();    
      
      hrec_unfolded_ratio[kRegDraw]->Draw("same");
      drawText( "Ratio plot check",	0.508173, 0.8659761, 22);
      drawText( "5.02 TeV ppMC, QCD Py6 Tune Z2",	0.508173, 0.8359761, 22);
      drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.508173, 0.8059761, 22);
      
      cRatioCheck->Print(outPdfFile.c_str());      
      

      //  MCClosure Plot --------------------------------------------------
      if(debugMode)std::cout<<std::endl<<"drawing MCClosure plot..."<<std::endl;      
      TCanvas *c1 = new TCanvas("c1","Spectra",1300,1000);  c1->cd();

      TH1F *hDum = new TH1F("MC_Unfold_Closure_genpt","MCClosure Plot",250, 30, 900);      
      hDum->SetYTitle("Unfolding Closure");//y-axis
      hDum->GetYaxis()->SetNdivisions(610); 
      hDum->GetYaxis()->SetLabelFont(43);
      hDum->GetYaxis()->SetLabelSize(20);
      hDum->GetYaxis()->SetTitleFont(43);
      hDum->GetYaxis()->SetTitleSize(22);
      hDum->GetYaxis()->SetTitleOffset(2.0);
      
      hDum->SetXTitle("Gen Jet p_{T} (GeV/c)");//x-axis
      hDum->GetXaxis()->SetLabelFont(43);
      hDum->GetXaxis()->SetLabelSize(20);
      hDum->GetXaxis()->SetTitleFont(43);
      hDum->GetXaxis()->SetTitleSize(22);
      hDum->GetXaxis()->SetTitleOffset(2.0);
      hDum->GetXaxis()->SetNoExponent();
      hDum->GetXaxis()->SetMoreLogLabels();
      hDum->SetAxisRange(0.8,1.2,"Y");
      
      hDum->Draw("hist");
      
      TLegend *leg0 = new TLegend(0.25, 0.15, 0.60, 0.35, NULL,"NBNDC");
      leg0->SetTextSize(0.07);
      
      if(doBayes){  
	if(debugMode)std::cout<<"Bayesian MCClosure was done too! adding to plot..."<<std::endl<<std::endl;

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
      
      TLine *line = new TLine(30,1,1000,1);
      line->SetLineStyle(2); 
      line->SetLineWidth(2);
      line->Draw();
      
      drawText( "MCClosure Tests", 0.608173, 0.8659761, 22);
      drawText( "5.02 TeV ppMC, QCD Py6 Tune Z2",	0.608173, 0.8359761, 21);
      drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.608173, 0.8059761, 21);
      
      leg0->AddEntry(hratio_svd[kRegDraw],"OppSide SVD","p");
      leg0->AddEntry(hratio_svd_check[kRegDraw],"SameSide SVD","p");
      leg0->Draw();
      
      c1->Print(outPdfFile.c_str());
      c1->Print(close_outPdfFile.c_str());
      std::cout<<std::endl<<"done drawing SVD PDFs!"<<std::endl<<std::endl;
    }// end drawPDFs
  }// end SVD specific
  
  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();

  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_MC->Close();  
  
  return 0;
} // end doMCClosureTests




//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){
  int rStatus = -1;
  if( argc!=1 && argc!=2 ) {//no input arguments, error
    std::cout<<"do ./doMCClosureTests.exe <debugMode,def0> to run"<<std::endl;
    return rStatus;
  }
  
  rStatus=1;//1 until a function returns 0 (normal exit by my convention)
  if(argc==1)   rStatus=doMCClosureTests();
  if(argc==2)   rStatus=doMCClosureTests( (const bool) atoi(argv[1]) ); 
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}











  //// ppData input histos --------------------------------------------------
  //std::cout<<std::endl<<std::endl<<"opening INPUT histos from DATA file"; 
  //if(debugMode)std::cout<<": "<<inFile_Data_name; 
  //std::cout<<std::endl<<std::endl;
  //TFile *fpp_Data = TFile::Open(inFile_Data.c_str());

      //if(debugMode)std::cout<<std::endl<<""<<std::endl;
      //cPearsonMatrixIter->SaveAs(Form("Pearson_Matrix_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");
      //cSpectra->SaveAs(Form("Spectra_meas_unf_fol_%s_R%d_plots.pdf",type.c_str(), radius),"RECREATE");
      //cRatio->SaveAs(Form("Ratio_meas_unf_fol_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");          
      //cdi->SaveAs(Form("di_Vectors_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
      //c11->SaveAs(Form("Singular_Values_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
      //cRatioCheck->SaveAs(Form("Spectra_Ratio_withkRegDraw_%s_R%d.pdf", type.c_str(), radius),"RECREATE");	      
