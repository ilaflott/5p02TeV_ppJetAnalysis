// Raghav Kunnawalkam Elayavalli
// June 24 2015
// Rutgers
// Macro to perform Bayes and SVD Unfolding direct from the official RooUnfold classes 

#if !defined(__CINT__) || defined(__MAKECINT__)
// C++, C, etc.
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
// ROOTSYS
#include <TStyle.h>
#include <TSystem.h>
#include <TProfile.h>
// I/O
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TEventList.h>
// utilities
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TMath.h>
// plotting
#include <TCanvas.h>
// histos
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
// roounfold
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/RooUnfold-1.1.1/src/RooUnfold.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#endif
// custom header
#include "utilities.h"

// CalculatePearsonCoefficients
//-----------------------------------------------------------------------------------------------------------------------
// explanation goes here, may be self explanatory
TMatrixT<double> *CalculatePearsonCoefficients(TMatrixT<double> *covmat) {

  TMatrixT<double> *pearsonCoefs = (TMatrixT<double>*)covmat->Clone("pearsonCoefs");
  //pearsonCoefs->Clear();

  Int_t nrows = covmat->GetNrows();
  Int_t ncols = covmat->GetNcols();
  Double_t pearson = 0.;

  for(int row = 0; row<nrows; row++) {
    for(int col = 0; col<ncols; col++) {
      pearson = covmat(row,col)/TMath::Sqrt(covmat(row,row)*covmat(col,col));
      pearsonCoefs(row,col) = pearson;
      //cout << "(" << row << "," << col << ") = " << pearson << endl;
    }
  }
  return pearsonCoefs;
} // end CalculatePearsonCoefficients


// PrintMatrix
//-----------------------------------------------------------------------------------------------------------------------
// Print the matrix as a table of elements.
// Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
// By default, format ="%11.4g", is used to print one element.
// One can specify an alternative format, e.g. format ="%6.2f  " 
// <ian note: white space in format assignment above on purpose?>
PrintMatrix( const TMatrixD& m, 
	     const char* format, const char* name, Int_t cols_per_sheet ){

  if (!m.IsValid()) {
    m.Error("PrintMatrix","%s is invalid",name);
    return;
  }

  const Int_t ncols  = m.GetNcols();
  const Int_t nrows  = m.GetNrows();
  const Int_t collwb = m.GetColLwb();
  const Int_t rowlwb = m.GetRowLwb();

  if (!(format && format[0])) format= "%11.4g ";
  char topbar[1000];
  snprintf(topbar,1000,format,123.456789);
  Int_t nch = strlen(topbar)+1;`
  if (nch > 18) nch = 18;
  char ftopbar[20];
  for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
  Int_t nk = 1 + Int_t(log10(ncols));
  snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
  Int_t nch2 = strlen(ftopbar);
  for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
  ftopbar[nch] = '|';
  ftopbar[nch+1] = 0;

  printf("\n%dx%d %s is as follows",nrows,ncols,name);

  if (cols_per_sheet <= 0) {
    cols_per_sheet = 5;
    if (nch <= 8) cols_per_sheet =10;
  }
  nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
  for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
  topbar[nk] = 0;
  for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
    printf("\n\n     |");
    for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
      printf(ftopbar,j+collwb-1);
    printf("\n%s\n",topbar);
    if (m.GetNoElements() <= 0) continue;
    for (Int_t i=1 ; i<=nrows ; i++) {
      printf("%4d |",i+rowlwb-1);
      for( Int_t j=sheet_counter ; j<sheet_counter+cols_per_sheet && j<=ncols ; j++) printf(format,m(i+rowlwb-1,j+collwb-1));
      printf("\n");
    }
  }
  printf("\n");
} // end PrintMatrix


//-----------------------------------------------------------------------------------------------------------------------
//static const int nKregMax = 9;
const int nKregMax = 9;
const int kIter = 4;
const int kRegDraw = 3;

void unfold(int radius = 4,
	    bool doToy = true, bool doBayes = true, bool doSVD = true,
	    int param = 21, std::string type = "Data"){

  TStopwatch timer; 
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);


// is this necessary?
#ifdef __CINT__
  gSystem->Load("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/RooUnfold-1.1.1/libRooUnfold.so");
#endif


  int kReg[nKregMax];  //{param-4, param-3, param-2, param-1, param, param+1, param+2, param+3, param+4};
  for(int i=(-4); (i+4)<=nKregMax; ++i)kReg[(i+4)]=param+i;

  TFile* fout;  
  if( type=="Data" )    fout = new TFile("PP_5p02_Spectra_Unfolded.root","RECREATE");
  else if( type=="MC" ) fout = new TFile("PP_5p02_ClosureTest.root","RECREATE");
  //fout = new TFile(outfileName,"RECREATE");
  fout->cd();
        
  // are these pointers defined correctly?
  TH1F *hgen, *hrec, *hrec_check;
  TH2F *hmat;

  TH1F *hgen_anabin, *hrec_anabin, * hrec_check_anabin;
  TH2F *hmat_anabin;
  TH1F *hgen_resp, *hrec_resp, *hgen_resp_anabin, *hrec_resp_anabin;
  
  TH1F *hunf, *hunf_check;
  TH1F *hunf_svd[nKregMax], * hunf_svd_check[nKregMax];

  TH1F *hratio, *hratio_check;
  TH1F *hratio_svd[nKregMax], *hratio_svd_check[nKregMax];  

  TH1 * hFoldedSVDPriorMeas[nKregMax];
  TH2 * hPearsonSVDPriorMeas[nKregMax];
  TH1 * hFoldedSVDPriorTruth[nKregMax];
  TH2 * hPearsonSVDPriorTruth[nKregMax];

  TH1 * hrec_folded_ratio[nKregMax];
  TH1 * hrec_unfolded_ratio[nKregMax];
  
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
  if(doToy) errorTreatment = RooUnfold::kCovToy;
  cout << "errorTreatment: " << errorTreatment << endl;
  
  TFile *fpp_MC, *fpp_Data; 
  // input files, i think this is hadd'd readFiles output specifically
  fpp_MC = TFile::Open("PP_5p02TeV_MC_ak4PF_20_eta_20.root");
  // for MC, we need to add the kinematic correction factors
  fpp_Data = TFile::Open("PP_5p02TeV_Data_ak4PF_20_eta_20.root");
  
  if(type == "Data"){
    // hKineticEff = (TH1F*)fPP_KineticEff->Get("hKineticEfficiency");    
    hgen = (TH1F*)fpp_MC->Get(Form("hpp_gen_R%d_20_eta_20",radius));
    hgen->Print("base");
    hrec = (TH1F*)fpp_Data->Get(Form("hpp_TrgCombTest_JetID_R%d_20_eta_20",radius));
    hrec->Print("base");

    hmat = (TH2F*)fpp_MC->Get(Form("hpp_matrix_R%d_20_eta_20",radius));
    hmat->Print("base");
    hmat_anabin = (TH2F*)hmat->Clone(Form("hpp_anaBin_Trans_matrix_HLT_R%d_20_eta_20",radius));
    hmat_anabin->Rebin2D(10, 10);
    hmat_anabin->Print("base");    

    hgen_anabin = (TH1F*)hgen->Clone(Form("hpp_anaBin_JetComb_gen_R%d_20_eta_20",radius));
    hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins, Form("hpp_anaBin_JetComb_gen_R%d_20_eta_20",radius), ptbins);
    divideBinWidth(hgen_anabin);
    hgen_anabin->Print("base");
    hrec_anabin = (TH1F*)hrec->Clone(Form("hpp_anaBin_HLTComb_R%d_20_eta_20",radius));
    hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins, Form("hpp_anaBin_HLTComb_R%d_20_eta_20",radius), ptbins);
    divideBinWidth(hrec_anabin);
    hrec_anabin->Print("base");    
  }
  
  if(type == "MC"){
    hgen = (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_gen_R%d_20_eta_20",radius));
    hgen->Print("base");
    hrec = (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_data_R%d_20_eta_20",radius));
    hrec->Print("base");
    hrec_check = (TH1F*)fpp_MC->Get(Form("hpp_mcclosure_data_train_R%d_20_eta_20",radius));

    hmat = (TH2F*)fpp_MC->Get(Form("hpp_mcclosure_matrix_R%d_20_eta_20",radius));
    hmat->Print("base");
    hmat_anabin = (TH2F*)hmat->Clone(Form("hpp_anaBin_Trans_matrix_HLT_R%d_20_eta_20",radius));
    hmat_anabin->Rebin2D(10, 10);
    hmat_anabin->Print("base");
    
    hgen_anabin = (TH1F*)hgen->Clone(Form("hpp_anaBin_mcclosure_gen_JetComb_R%d_20_eta_20",radius));
    hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins, Form("hpp_anaBin_mcclosure_gen_JetComb_R%d_20_eta_20",radius), ptbins);
    //hgen_anabin = (TH1F*)hgen_anabin->Rebin(10);
    divideBinWidth(hgen_anabin);
    hgen_anabin->Print("base");
    hrec_anabin = (TH1F*)hrec->Clone(Form("hpp_anaBin_mcclosure_rec_JetComb_R%d_20_eta_20",radius));
    hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins, Form("hpp_anaBin_mcclosure_rec_JetComb_R%d_20_eta_20",radius), ptbins);
    //hrec_anabin = (TH1F*)hrec_anabin->Rebin(10);
    divideBinWidth(hrec_anabin);
    hrec_anabin->Print("base");
    hrec_check->Print("base");
    //if(dotrigcor) doTrigCorr(hrec_check, hMC_turnon);
    hrec_check_anabin = (TH1F*)hrec_check->Clone(Form("hpp_mcclosure_JetComb_data_train_R%d_20_eta_20",radius));
    hrec_check_anabin = (TH1F*)hrec_check_anabin->Rebin(nbins, Form("hpp_mcclosure_JetComb_data_train_R%d_20_eta_20",radius), ptbins);
    //hrec_check_anabin = (TH1F*)hrec_check_anabin->Rebin(10);
    divideBinWidth(hrec_check_anabin);
    //if(dotrigcor) doTrigCorr(hrec_anabin, hMC_turnon);
  }
  
  hrec_resp = new TH1F(Form("hpp_rec_response_R%d_20_eta_20",radius),"", hrec->GetNbinsX(), hrec->GetXaxis()->GetXmin(), hrec->GetXaxis()->GetXmax());
  hrec_resp->Sumw2();

  hrec_resp_anabin = new TH1F(Form("hpp_rec_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hrec_resp_anabin->Sumw2();

  hgen_resp = new TH1F(Form("hpp_gen_responseR%d_20_eta_20",radius),"", hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  hgen_resp->Sumw2();

  hgen_resp_anabin = new TH1F(Form("hpp_gen_response_anabin_R%d_20_eta_20",radius),"", nbins, ptbins);
  hgen_resp_anabin->Sumw2();
  
  // declare the canvas
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
 
  if(doBayes){
    // BAYESIAN unfolding!
      
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_R%d",radius));
    RooUnfoldBayes unf_bayes(&roo_resp, hrec_anabin, kIter);
    hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);
    hunf->SetName("PP_Bayesian_Unfolded_Spectra");
    hratio = (TH1F*)hunf->Clone(Form("MCClosure_test_oppside_Bayes_R%d_20_eta_20", radius));
    hratio->Print("base");
    hgen_anabin->Print("base");
    hratio->Divide(hgen_anabin);
    hratio->SetMarkerStyle(24);
    hratio->SetMarkerColor(kRed);

    if(type == "MC"){
      RooUnfoldBayes unf_bayes_check(&roo_resp, hrec_check_anabin, kIter);
      hunf_check = (TH1F*)unf_bayes_check.Hreco(errorTreatment);
      hunf_check->Print("base");
      hratio_check = (TH1F*)hunf_check->Clone(Form("MCClosure_test_sameside_Bayes_R%d_20_eta_20", radius));
      hratio_check->Print("base");
      hgen_anabin->Print("base");
      hratio_check->Divide(hgen_anabin);
      hratio_check->SetMarkerStyle(24);
      hratio_check->SetMarkerColor(kBlack);

    }     

  }


  TLegend* leg[9];
  TLegend* leg1[9];  
    
  if(doSVD){
    // // SVD unfolding!
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, Form("Response_matrix_anabin_R%d",radius));
      
    for(int kr = 0; kr < 9; ++kr){
      cout<<kr<<endl;
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
      TMatrixD *pearson = (TMatrixD*)CalculatePearsonCoefficients(&covmat);
      hPearsonSVDPriorMeas[kr] = new TH2D(*pearson);
      hPearsonSVDPriorMeas[kr]->SetTitle(" ");
      hPearsonSVDPriorMeas[kr]->SetName(Form("hPearsonSVDPriorMeas_kReg%d", kReg[kr]));
      hPearsonSVDPriorMeas[kr]->SetMinimum(-1.);
      hPearsonSVDPriorMeas[kr]->SetMaximum(1.);
      hPearsonSVDPriorMeas[kr]->GetZaxis()->SetLabelSize(0.06);
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
      //Note that these do not depend on the regularization. The opposite: they tell you which regularization to use!

      if(kr == 0){
	
	TSVDUnfold *svdUnfold = unf_svd.Impl();
	TH1 *hSVal = svdUnfold->GetSV();
	TH1D *hdi = svdUnfold->GetD();
	TH1D  *hSValClone = (TH1D*)hSVal->Clone("hSValClone");
	cout << "print singular values " << endl;
	for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
	  cout << "bin: " << bin << "  SV: " << hSVal->GetBinContent(bin) << endl;
	cout << "print di vector " <<  endl;
	for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
	  cout << "i: " << bin << "  di: " << hdi->GetBinContent(bin) << endl;

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
      }
	
      if(type == "MC"){
	  
	RooUnfoldSvd unf_svd_check(&roo_resp, hrec_check_anabin, kReg[kr]);
	hunf_svd_check[kr] = (TH1F*)unf_svd_check.Hreco(errorTreatment);
	hunf_svd_check[kr]->SetName(Form("h_Unfolded_mcclosure_sameside_SVD_kReg%d_check_R%d_20_eta_20", kReg[kr],radius));
	hunf_svd_check[kr]->Print("base");
    
	hratio_svd_check[kr] = (TH1F*)hunf_svd_check[kr]->Clone(Form("MCClosure_test_sameside_SVD_kReg%d_R%d_20_eta_20",  kReg[kr], radius));
	hratio_svd_check[kr]->Divide(hgen_anabin);
	hratio_svd_check[kr]->SetMarkerStyle(33);
	hratio_svd_check[kr]->SetMarkerColor(kBlack);

	//Get covariance matrix and calculate corresponding Pearson coefficients    
	TMatrixD covmat_truth = unf_svd_check.Ereco(errorTreatment);
	TMatrixD *pearson_truth = (TMatrixD*)CalculatePearsonCoefficients(&covmat_truth);

	cout << "print pearson coefficients" << endl;      
	hPearsonSVDPriorTruth[kr] = new TH2D(*pearson_truth);
	
	hPearsonSVDPriorTruth[kr]->SetName(Form("hPearsonSVDPriorTruth_kReg%d",kReg[kr]));
	hPearsonSVDPriorTruth[kr]->SetMinimum(-1.);
	hPearsonSVDPriorTruth[kr]->SetMaximum(1.);
	hPearsonSVDPriorTruth[kr]->GetZaxis()->SetLabelSize(0.06);

	hFoldedSVDPriorTruth[kr] = roo_resp.ApplyToTruth(hunf_svd_check[kr]);
	hFoldedSVDPriorTruth[kr]->SetName(Form("hFoldedSVDPriorTruth_kReg%d", kReg[kr]));

	cout<<endl<<endl<<endl<<"MC: kReg = "<<kReg[kr]<<endl<<endl<<endl;
	TSVDUnfold *svdUnfold = unf_svd_check.Impl();
	TH1 *hSVal = svdUnfold->GetSV();
	TH1D *hdi = svdUnfold->GetD();
	TH1D  *hSValClone = (TH1D*)hSVal->Clone("hSValClone");
	cout << "print singular values " << endl;
	for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
	  cout << "bin: " << bin << "  SV: " << hSVal->GetBinContent(bin) << endl;

	cout << "print di vector " <<  endl;
	for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
	  cout << "i: " << bin << "  di: " << hdi->GetBinContent(bin) << endl;
	  
      }// MC conditional

      hrec_unfolded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Print("base");
      
    }// kReg loop
    
  }// SVD conditional

  cPearsonMatrixIter->SaveAs(Form("Pearson_Matrix_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");
  cSpectra->SaveAs(Form("Spectra_meas_unf_fol_%s_R%d_plots.pdf",type.c_str(), radius),"RECREATE");
  cRatio->SaveAs(Form("Ratio_meas_unf_fol_%s_R%d_plots.pdf", type.c_str(), radius),"RECREATE");          
  cdi->SaveAs(Form("di_Vectors_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
  c11->SaveAs(Form("Singular_Values_%s_R%d.pdf",type.c_str(), radius),"RECREATE");
  cSpectraCheck->cd();
  cout<<"Passed here !"<<endl<<endl;
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
  cSpectraCheck->SaveAs(Form("Spectra_Ratio_withkRegDraw_R%d.pdf", radius),"RECREATE");	

  fout->cd();

  hgen->Write();
  hrec->Write();    
  if(type == "MC")  hrec_check->Write();
  hgen_resp->Write();
  hrec_resp->Write();
    
  if(doBayes){ 
    hmat->Write();
    hunf->Write();
    if(type == "MC") hunf_check->Write();
    hratio->Write();
    if(type == "MC") hratio_check->Write();    
  }
    
  if(doSVD){
    for(int kr = 0; kr<7; ++kr){
      hunf_svd[kr]->Write();
      hPearsonSVDPriorMeas[kr]->Write();
      hFoldedSVDPriorMeas[kr]->Write();
      hratio_svd[kr]->Write();
      if(type == "MC"){
	hunf_svd_check[kr]->Write();
	hratio_svd_check[kr]->Write();
	hPearsonSVDPriorTruth[kr]->Write();
	hFoldedSVDPriorTruth[kr]->Write();
      }
    }
  }

  
  if(type == "MC"){
    //hgen[0]->Draw();
    TLegend *leg0 = new TLegend(0.25, 0.15, 0.60, 0.35, NULL,"BRNDC");
    // leg0->SetHeader("Bayesian Unfolding");
    leg0->SetTextSize(0.07);
    leg0->AddEntry(hratio,"OppSide Bayesian","p");
    leg0->AddEntry(hratio_check,"SameSide Bayesian","p");
    leg0->AddEntry(hratio_svd[kRegDraw],"OppSide SVD","p");
    leg0->AddEntry(hratio_svd_check[kRegDraw],"SameSide SVD","p");
    leg0->SetTextSize(0.04);
    
    TLine *line = new TLine(30,1,1000,1);
    line->SetLineStyle(2); 
    line->SetLineWidth(2);
    
    TH1F *hDum = new TH1F("hDum","",250, 30, 900);
    hDum->SetXTitle("Gen Jet p_{T} (GeV/c)");
    hDum->SetYTitle("Unfolding Closure");
    //TH1F *hDum = GetDummyHist(50, 300.0, 0, 2,"hDum","Gen Jet p_{T} (GeV/c)","PbPb Unfolding Closure", false);
    hDum->GetYaxis()->SetNdivisions(610); 
    hDum->GetYaxis()->SetLabelFont(43);
    hDum->GetYaxis()->SetTitleFont(43);
    hDum->GetYaxis()->SetLabelSize(20);
    hDum->GetYaxis()->SetTitleSize(22);
    hDum->GetYaxis()->SetTitleOffset(2.6);
    hDum->GetXaxis()->SetLabelFont(43);
    hDum->GetXaxis()->SetTitleFont(43);
    hDum->GetXaxis()->SetLabelSize(20);
    hDum->GetXaxis()->SetTitleSize(22);
    hDum->GetXaxis()->SetTitleOffset(2.3);
    hDum->GetXaxis()->SetNoExponent();
    hDum->GetXaxis()->SetMoreLogLabels();
    
    int ipad=0;
    TCanvas *c1 = new TCanvas("c1","Spectra",1300,1000);
    c1->cd();

    hDum->SetAxisRange(0.8,1.2,"Y");
    hDum->Draw("hist");
    hratio->Draw("psame");
    hratio_check->Draw("psame");
    hratio_svd[kRegDraw]->Draw("psame");
    hratio_svd_check[kRegDraw]->Draw("psame");
    drawText("PP, 5 TeV", 0.648173, 0.8459761, 22);
    line->Draw();
    drawText(Form("kReg = %d", kReg[kRegDraw]), 0.64, 0.7, 21);
    drawText("PYTHIA",0.22,0.88,21);
    leg0->Draw();

    c1->SaveAs(Form("PP_MCClosure_R%d.pdf",radius),"RECREATE");
    fout->Write();
    fout->Close();
  }
} // end unfold
