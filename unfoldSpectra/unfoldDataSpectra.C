#include "unfoldSpectra.h"

const int radius=4;
const bool doToyErrs=false;
const bool debugMode=true;


// ----------------------------------------------------------------------------------------------------------------------
int unfoldDataSpectra( ){

  // performance tracking
  TStopwatch timer; 
  timer.Start();
  
  // set default error calculation for TH1/2's from here on out + gStyle
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);


  // ppData input histos -------------------------

  std::cout<<std::endl<<std::endl<<"opening INPUT histos from DATA file"; 
  if(debugMode)std::cout<<": "<<inFile_Data_name; 
  std::cout<<std::endl<<std::endl;
  TFile *fpp_Data = TFile::Open(inFile_Data.c_str());

  //rec
  TH1F *hrec, *hrec_anabin;
  hrec = (TH1F*)fpp_Data->Get(Form("hpp_TrgCombTest_JetID_R%d_20_eta_20",radius));
  hrec->Print("base");
  
  //cloned from input histo, not "real" input histo
  hrec_anabin = (TH1F*)hrec->Clone(Form("hpp_anaBin_HLTComb_R%d_20_eta_20",radius));
  hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins, Form("hpp_anaBin_HLTComb_R%d_20_eta_20",radius), ptbins);
  divideBinWidth(hrec_anabin);
  hrec_anabin->Print("base");    
  

  // ppMC input histos -------------------------

  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"; 
  if(debugMode)std::cout<<": "<<inFile_MC_name; 
  std::cout<<std::endl<<std::endl;
  TFile *fpp_MC = TFile::Open(inFile_MC.c_str());

  //mat?
  TH2F *hmat, *hmat_anabin;
  hmat = (TH2F*)fpp_MC->Get(Form("hpp_matrix_R%d_20_eta_20",radius));
  hmat->Print("base");
  
  //cloned from input histo, not "real" input histo
  hmat_anabin = (TH2F*)hmat->Clone(Form("hpp_anaBin_Trans_matrix_HLT_R%d_20_eta_20",radius));
  hmat_anabin->Rebin2D(10, 10);
  hmat_anabin->Print("base");    
  
  //gen
  TH1F *hgen, *hgen_anabin;
  hgen = (TH1F*)fpp_MC->Get(Form("hpp_gen_R%d_20_eta_20",radius));
  hgen->Print("base");
  
  //cloned from input histo, not "real" input histo
  hgen_anabin = (TH1F*)hgen->Clone(Form("hpp_anaBin_JetComb_gen_R%d_20_eta_20",radius));
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins, Form("hpp_anaBin_JetComb_gen_R%d_20_eta_20",radius), ptbins);
  divideBinWidth(hgen_anabin);
  hgen_anabin->Print("base");


  // output histos -------------------------
  std::cout<<std::endl<<std::endl<<"creating new OUTPUT histos..."<<std::endl;
  //gen
  if(debugMode)std::cout<<"opening output gen/genresponse histograms..."<<std::endl;
  TH1F *hgen_resp, *hgen_resp_anabin;
  hgen_resp = new TH1F( Form("hpp_gen_responseR%d_20_eta_20",radius),"", 
			hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  hgen_resp->Sumw2();
  hgen_resp->Print("base");
  hgen_resp_anabin = new TH1F(Form("hpp_gen_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hgen_resp_anabin->Sumw2();
  hgen_resp_anabin->Print("base");
  
  //rec
  if(debugMode)std::cout<<"opening output rec/recresponse histograms..."<<std::endl;
  TH1F *hrec_resp, *hrec_resp_anabin; 
  hrec_resp = new TH1F(Form("hpp_rec_response_R%d_20_eta_20",radius),"", 
		       hrec->GetNbinsX(), hrec->GetXaxis()->GetXmin(), hrec->GetXaxis()->GetXmax());
  hrec_resp->Sumw2();
  hrec_resp->Print("base");
  hrec_resp_anabin = new TH1F(Form("hpp_rec_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hrec_resp_anabin->Sumw2();
  hrec_resp_anabin->Print("base");
  
  // open output file for later
  const std::string outFile="PP_5p02_Spectra_Unfolded.root";
  TFile* fout = new TFile(outFile.c_str(),"RECREATE");


  // unfolding settings ------------------------- 
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
  if(doToyErrs) errorTreatment = RooUnfold::kCovToy;
  std::cout<<std::endl<<"errorTreatment: "<<errorTreatment<<std::endl;
  
  // Bayesian unfolding ------------------------- 
  const bool doBayes=true; 
  const int kIter = 4;
  if(doBayes){
    std::cout<<std::endl<<" beginning Bayesian unfolding..."<<std::endl;
    std::cout<<"calling RooUnfoldResponse "<<std::endl;
    RooUnfoldResponse roo_resp( hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_R%d",radius));

    std::cout<<"calling RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );

    TH1F *hunf; 
    hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);
    hunf->SetName("PP_Bayesian_Unfolded_Spectra");

    TH1F *hratio;
    hratio = (TH1F*)hunf->Clone(Form("MCClosure_test_oppside_Bayes_R%d_20_eta_20", radius));
    hratio->Print("base");
    hgen_anabin->Print("base");

    std::cout<<"dividing hunf and hgen_anabin..."<<std::endl;
    hratio->Divide(hgen_anabin);
    hratio->SetMarkerStyle(24);
    hratio->SetMarkerColor(kRed);

    std::cout<<"writing bayesian unfolding output..."<<std::endl;
    fout->cd();
    hmat->Write();
    hunf->Write();
    hratio->Write();

    std::cout<<"Bayesian Unfolding done!"<<std::endl;
  } // end bayesian specifics
  else  std::cout<<std::endl<<" skipping Bayesian unfolding..."<<std::endl;


  // SVD unfolding -------------------------
  const bool doSVD   =true; // testing now
  const bool drawPDFs=true; //untested
  const int  nKregMax = 9, kRegDraw = 3, param=21;
  if(doSVD){

    std::cout<<std::endl<<" beginning SVD unfolding..."<<std::endl;

    // not sure what this is for, some kinda regularization parameter?
    int kReg[nKregMax];  //{param-4, param-3, param-2, param-1, param, param+1, param+2, param+3, param+4};
    for(int i=(-4); (i+4)<nKregMax; ++i) kReg[i+4]=param+i;
    
    TH1F *hunf_svd[nKregMax];
    TH1F *hratio_svd[nKregMax];
    TH1 *hrec_folded_ratio[nKregMax], *hrec_unfolded_ratio[nKregMax];
    
    // pearson coefficients
    TH2 *hPearsonSVDPriorMeas[nKregMax]; 
    TH1 *hFoldedSVDPriorMeas[nKregMax];
    
    std::cout<<"calling RooUnfoldResponse"<<std::endl;
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_anabin_R%d",radius));

    // declare the canvases
    std::cout<<"opening TCanvases..."<<std::endl;
    TCanvas *cPearsonMatrixIter; 
    cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","",1200,1000);
    cPearsonMatrixIter->Divide(3,3);
    
    TCanvas * cSpectra;
    cSpectra = new TCanvas("cSpectra","",1200,1000);
    cSpectra->Divide(3,3);
    
    TCanvas * cRatio;
    cRatio = new TCanvas("cRatio","",1200,1000);
    cRatio->Divide(3,3);  
    
    TCanvas *c11 = new TCanvas("c11","c11: Singular Values", 1200, 1000);
    TCanvas *cdi = new TCanvas("cdi","cdi: di vectors", 1200, 1000);
    TCanvas *cSpectraCheck = new TCanvas("cSpectraCheck","",1200, 1000);
    
    TLegend *leg[nKregMax],*leg1[nKregMax];      

    std::cout<<"entering SVD Unfolding Loop..."<<std::endl;
    for(int kr = 0; kr < nKregMax; ++kr){
      std::cout<<"int kr= "<<kr<<std::endl;
  
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
  
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      hunf_svd[kr]->SetName(Form("PP_SVD_Unfolded_Spectra_kReg%d",kReg[kr]));
      //hunf_svd[kr]->SetName(Form("h_Unfolded_mcclosure_spectra_SVD_kReg%d_R%d_20_eta_20", kReg[kr],radius));
      hunf_svd[kr]->Print("base");
  
      hratio_svd[kr] = (TH1F*)hunf_svd[kr]->Clone(Form("MCClosure_test_oppside_SVD_kReg%d_R%d_20_eta_20", kReg[kr], radius));
      hratio_svd[kr]->Divide(hgen_anabin);
      hratio_svd[kr]->SetMarkerStyle(33);
      hratio_svd[kr]->SetMarkerColor(kRed);
  
      //Get covariance matrix and calculate corresponding Pearson coefficients    
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);
      //      std::cout<<"&covmat="<<&covmat<<std::endl;
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
    }// kReg loop
    
    // loop over histos pointers in arrays to write
    fout->cd();
    for(int kr = 0; kr<nKregMax; ++kr){
      hunf_svd[kr]->Write();
      hPearsonSVDPriorMeas[kr]->Write();
      hFoldedSVDPriorMeas[kr]->Write();
      hratio_svd[kr]->Write();
    }

    if(drawPDFs){

      std::string type="unfoldData";//should only be Data or MC
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
    }// end drawPDFs
  }// end SVD specific
  else  std::cout<<std::endl<<" skipping SVD unfolding..."<<std::endl;

  std::cout<<"writing histos to output file"<<std::endl;
  fout->cd();
  hgen->Write();
  hgen_resp->Write();
  hrec->Write();    
  hrec_resp->Write();
  fout->Write();

  std::cout<<"cleaning up"<<std::endl;
  fout->Close();
  fpp_MC->Close();  
  fpp_Data->Close();

  return 0;
} // end unfoldDataSpectra


//-----------------------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]){
  int rStatus = -1;
  if( argc!=1 ) {
    std::cout<<"settings hard coded, just do ./unfoldDataSpectra.exe"<<std::endl;
    return rStatus;
  }
  rStatus=1;
  rStatus=unfoldDataSpectra();
  std::cout<<"done! return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}
