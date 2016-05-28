// Raghav Kunnawalkam Elayavalli
// June 5th 2014
// CERN
// for questions or comments: raghav.k.e at CERN dot CH

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include "../../Headers/plot.h"


using namespace std;

void PP_MC_normResponseMatrix(int radius = 3){

  gStyle->SetOptStat(0);
  
  bool printDebug = true;
  
  TFile * fin = TFile::Open(Form("pp_MC_new_ak%d_20_eta_20.root",radius));
  
  TH2F * mPP_Matrix;
  TH2F * mPP_Response;
  TH2F * mPP_ResponseNorm;
  TH2F * mPP_ResponseNorm_recoTrunc;

  // TH2F * mPbPb_Matrix_recoTrunc[nbins_cent];
  // TH2F * mPbPb_Response_recoTrunc[nbins_cent];
  // TH2F * mPbPb_ResponseNorm_recoTrunc[nbins_cent];

  TFile fout(Form("PP_kineticEfficiency_correctionFactors_R%d.root",radius),"RECREATE");
  fout.cd();
  
  TH1F * hGen_Projection;
  TH1F * hGen_Projection_recoTrunc;
  TH1F * hKineticEfficiency;

  TCanvas * cGenProj = new TCanvas("cGenProj","",800,600);
  TCanvas * cKineticEfficiency = new TCanvas("cKineticEfficiency","",800,600);  

  TCanvas * cMatrix = new TCanvas("cMatrix","",800,600);
  TCanvas * cMatrix_recoTrunc = new TCanvas("cMatrix_recoTrunc","",800,600);
  
  mPP_Matrix = (TH2F*)fin->Get(Form("hpp_matrix_HLT_R%d_20_eta_20",radius));
  mPP_Matrix->Rebin2D(10, 10);

  for (int y=1;y<=mPP_Matrix->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_Matrix->GetNbinsX();x++) {
      if (mPP_Matrix->GetBinContent(x,y)<=1*mPP_Matrix->GetBinError(x,y)) {
	mPP_Matrix->SetBinContent(x,y,0);
	mPP_Matrix->SetBinError(x,y,0);
      }
      sum+=mPP_Matrix->GetBinContent(x,y);
    }
    
    for (int x=1;x<=mPP_Matrix->GetNbinsX();x++) {	   
      double ratio = 1;
      // if (hGenSpectraCorrPP->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorrPP->GetBinContent(x);
      mPP_Matrix->SetBinContent(x,y,mPP_Matrix->GetBinContent(x,y)*ratio);
      mPP_Matrix->SetBinError(x,y,mPP_Matrix->GetBinError(x,y)*ratio);
    }
  }
  // mPbPb_Matrix[i]->Smooth(0);
  
  // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  

  if(printDebug) cout<<"getting the response matrix"<<endl;

  mPP_Response = (TH2F*)mPP_Matrix->Clone("mPP_Response");
  TH1F *hProjPP = (TH1F*)mPP_Response->ProjectionY()->Clone("hProjPP");
  
  
  for (int y=1;y<=mPP_Response->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_Response->GetNbinsX();x++) {
      if (mPP_Response->GetBinContent(x,y)<=1*mPP_Response->GetBinError(x,y)) {
	// in the above if statement, kurt has 1*error and my old has 0*error
	mPP_Response->SetBinContent(x,y,0);
	mPP_Response->SetBinError(x,y,0);
      }
      sum+=mPP_Response->GetBinContent(x,y);
    }
    
    for (int x=1;x<=mPP_Response->GetNbinsX();x++) {  	
      if (sum==0) continue;
      double ratio = 1;
      //if(dPbPb_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
      // else ratio = dPbPb_TrgComb[i]->GetBinContent(y)/sum
      ratio = 1./sum;
      if (hProjPP->GetBinContent(y)==0) ratio = 1e-100/sum;
      else ratio = hProjPP->GetBinContent(y)/sum;
      mPP_Response->SetBinContent(x,y,mPP_Response->GetBinContent(x,y)*ratio);
      mPP_Response->SetBinError(x,y,mPP_Response->GetBinError(x,y)*ratio);
    }
  }
  if(printDebug) cout<<"getting the normalized response matrix"<<endl;
  mPP_ResponseNorm = (TH2F*)mPP_Matrix->Clone("mPP_ResponseNorm");
  for (int x=1;x<=mPP_ResponseNorm->GetNbinsX();x++) {
    double sum=0;
    for (int y=1;y<=mPP_ResponseNorm->GetNbinsY();y++) {
      if (mPP_ResponseNorm->GetBinContent(x,y)<=1*mPP_ResponseNorm->GetBinError(x,y)) {
	mPP_ResponseNorm->SetBinContent(x,y,0);
	mPP_ResponseNorm->SetBinError(x,y,0);
      }
      sum+=mPP_ResponseNorm->GetBinContent(x,y);
    }
    
    for (int y=1;y<=mPP_ResponseNorm->GetNbinsY();y++) {  	
      if (sum==0) continue;
      double ratio = 1./sum;
      mPP_ResponseNorm->SetBinContent(x,y,mPP_ResponseNorm->GetBinContent(x,y)*ratio);
      mPP_ResponseNorm->SetBinError(x,y,mPP_ResponseNorm->GetBinError(x,y)*ratio);
    }
    
    
  }

  
  hGen_Projection = (TH1F*)mPP_ResponseNorm->ProjectionX();
  hGen_Projection->SetName("hGen_Projection");
  // hGen_Projection->Rebin(10);

  cout<<"Finished getting histograms for full untruncated response matrix"<<endl;

  mPP_ResponseNorm_recoTrunc = (TH2F*)mPP_ResponseNorm->Clone("mPP_ResponseNorm_recoTrunc");

  for(int x = 0; x < mPP_ResponseNorm->GetNbinsX(); ++x){

    for(int y = 0; y < mPP_ResponseNorm->GetNbinsY(); ++y){

      if(y < mPP_ResponseNorm->GetYaxis()->FindBin(50) || y > mPP_ResponseNorm->GetYaxis()->FindBin(350)){
	mPP_ResponseNorm_recoTrunc->SetBinContent(x, y, 0.0);
	mPP_ResponseNorm_recoTrunc->SetBinError(x, y, 0.0);
      }
	
    }
      
  }
    
  hGen_Projection_recoTrunc = (TH1F*)mPP_ResponseNorm_recoTrunc->ProjectionX();
  hGen_Projection_recoTrunc->SetName("hGen_Projection_recoTrunc");
  // hGen_Projection_recoTrunc->Rebin(10);

  cout<<"finished getting histograms for truncted response matrix"<<endl;
    
  hKineticEfficiency = (TH1F*)hGen_Projection_recoTrunc->Clone("hKineticEfficiency");
  hKineticEfficiency->Divide(hGen_Projection);

  for(int bin = 1; bin<=hKineticEfficiency->GetNbinsX(); ++bin){
    double val = hKineticEfficiency->GetBinContent(bin);
    val = 1./val;
    hKineticEfficiency->SetBinContent(bin, val);
  }

  cout<<"Plotting kinetic efficiency"<<endl;
    
  cKineticEfficiency->cd();

  hKineticEfficiency->SetTitle(" ");
  hKineticEfficiency->SetXTitle("Gen p_{T} (GeV/c)");
  hKineticEfficiency->SetYTitle("Kinetic Efficiency");

  hKineticEfficiency->SetMarkerStyle(20);
  hKineticEfficiency->SetMarkerColor(kBlack);
  hKineticEfficiency->SetAxisRange(60, 299, "X");
  hKineticEfficiency->SetAxisRange(0.98, 1.05, "Y");
  hKineticEfficiency->Draw();

  drawText("PYTHIA, Kinetic Efficiency", 0.25,0.7,14);
  drawText(Form("ak%dPF Jets",radius), 0.25, 0.65, 14);
    
  cout<<"Finished plotting kinetic efficiency"<<endl;

  cGenProj->cd();
  hGen_Projection->SetTitle(" ");
  hGen_Projection->SetXTitle("Gen p_{T} (GeV/c)");
  hGen_Projection->SetYTitle("GenProjection of Normalized Response Matrix");
  hGen_Projection->Scale(1./10);
  hGen_Projection->SetMarkerStyle(33);
  hGen_Projection->SetMarkerColor(kBlack);
  hGen_Projection->SetAxisRange(20, 350, "X");
  hGen_Projection->SetAxisRange(0.6, 1.4, "Y");
  hGen_Projection->Draw();

  hGen_Projection_recoTrunc->Scale(1./10);
  hGen_Projection_recoTrunc->SetMarkerStyle(25);
  hGen_Projection_recoTrunc->SetMarkerColor(kRed);
  hGen_Projection_recoTrunc->Draw("same");

  drawText("PYTHIA, Projection onto Gen axis", 0.25,0.7,14);
  drawText(Form("ak%dPF Jets",radius), 0.25, 0.65, 14);
    
  TLegend * leg = myLegend(0.3,0.6,0.7,0.7);
  leg->AddEntry(hGen_Projection,"No Truncation in Reco pT","pl");
  leg->AddEntry(hGen_Projection_recoTrunc,"reco pT truncated 40 to 350 GeV","pl");
  leg->SetTextSize(0.04);
  leg->Draw();


  cMatrix->cd();
  // mPP_ResponseNorm->Rebin2D(10, 10);
  // mPP_ResponseNorm->Scale(1./TMath::Sqrt(10));
  makeHistTitle(mPP_ResponseNorm," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");
  mPP_ResponseNorm->SetAxisRange(0, 500, "X");
  mPP_ResponseNorm->SetAxisRange(0, 500, "Y");
  mPP_ResponseNorm->SetAxisRange(0.001, 1, "Z");
  mPP_ResponseNorm->Draw("colz");  
  drawText("PYTHIA", 0.15,0.8,14);
  drawText(Form("ak%dPF Jets",radius), 0.15, 0.75, 14);

  
  cMatrix_recoTrunc->cd();
  // mPP_ResponseNorm_recoTrunc->Rebin2D(10, 10);
  // mPP_ResponseNorm_recoTrunc->Scale(1./TMath::Sqrt(10));
  makeHistTitle(mPP_ResponseNorm_recoTrunc," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");
  mPP_ResponseNorm_recoTrunc->SetAxisRange(0, 500, "X");
  mPP_ResponseNorm_recoTrunc->SetAxisRange(0, 500, "Y");
  mPP_ResponseNorm_recoTrunc->SetAxisRange(0.001, 1, "Z");
  mPP_ResponseNorm_recoTrunc->Draw("colz");  
  drawText("PYTHIA, reco p_{T} truncated", 0.15,0.8,14);
  drawText(Form("ak%dPF Jets",radius), 0.15, 0.75, 14);
  
  cKineticEfficiency->SaveAs(Form("KineticEfficiency_R%d_PPMC.pdf",radius),"RECREATE");
  cGenProj->SaveAs(Form("GenProjection_R%d_PPMC.pdf",radius),"RECREATE");
  cMatrix->SaveAs(Form("ResponseMatrix_R%d_PPMC.pdf",radius),"RECREATE");
  cMatrix_recoTrunc->SaveAs(Form("ResponseMatrix_recoTrunc_R%d_PPMC.pdf",radius),"RECREATE");
  
  
  fout.Write();
  fout.Close();
  
}
