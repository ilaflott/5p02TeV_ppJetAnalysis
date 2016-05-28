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

void PbPb_MC_normResponseMatrix(int radius = 3){

  gStyle->SetOptStat(0);
  
  const int nbins_cent = 6;
  const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
  bool printDebug = true;
  
  TFile * fin = TFile::Open(Form("PbPb_MC_new_akPu%d_20_eta_20.root",radius));
  
  TH2F * mPbPb_Matrix[nbins_cent];
  TH2F * mPbPb_Response[nbins_cent];
  TH2F * mPbPb_ResponseNorm[nbins_cent];

  // TH2F * mPbPb_Matrix_recoTrunc[nbins_cent];
  // TH2F * mPbPb_Response_recoTrunc[nbins_cent];
  TH2F * mPbPb_ResponseNorm_recoTrunc[nbins_cent];

  TFile fout(Form("PbPb_kineticEfficiency_correctionFactors_R%d.root",radius),"RECREATE");
  fout.cd();
  
  TH1F * hGen_Projection[nbins_cent];
  TH1F * hGen_Projection_recoTrunc[nbins_cent];
  TH1F * hKineticEfficiency[nbins_cent];

  TCanvas * cGenProj = new TCanvas("cGenProj","",1200,1000);
  makeMultiPanelCanvas(cGenProj,3,2,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cKineticEfficiency = new TCanvas("cKineticEfficiency","",1200,1000);
  makeMultiPanelCanvas(cKineticEfficiency,3,2,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMatrix = new TCanvas("cMatrix","",1200,1000);
  makeMultiPanelCanvas(cMatrix,3,2,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMatrix_recoTrunc = new TCanvas("cMatrix_recoTrunc","",1200,1000);
  makeMultiPanelCanvas(cMatrix_recoTrunc,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0;i<nbins_cent;++i){
    if(printDebug) cout<<"centrality bin iteration = "<<i<<endl;

    mPbPb_Matrix[i] = (TH2F*)fin->Get(Form("hpbpb_matrix_HLT_R%d_20_eta_20_cent%d",radius,i));
    mPbPb_Matrix[i]->Print("base");
    mPbPb_Matrix[i]->Rebin2D(10, 10);

    //mPbPb_Matrix_recoTrunc[i] = (TH2F*)fin->Get(Form("hpbpb_matrix_recoTrunc_HLT_R%d_20_eta_20_cent%d",radius,i));
    //mPbPb_Matrix_recoTrunc[i]->Print("base");
    
    // TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
    // f->SetParameters(1e10,-8.8,40);
    // TH1F *hGenSpectraCorr = (TH1F*)mPbPb_Matrix[i]->ProjectionX()->Clone(Form("hGenSpectraCorr_cent%d",i));
    // hGenSpectraCorr->Fit("f"," ");
    // hGenSpectraCorr->Fit("f","","");
    // hGenSpectraCorr->Fit("f","LL");
    // TH1F *fHist = functionHist(f,hGenSpectraCorr,Form("fHist_cent%d",i));// function that you get from the fitting 
    // hGenSpectraCorr->Divide(fHist);

    for (int y=1;y<=mPbPb_Matrix[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_Matrix[i]->GetNbinsX();x++) {
	if (mPbPb_Matrix[i]->GetBinContent(x,y)<=1*mPbPb_Matrix[i]->GetBinError(x,y)) {
	  //in the above line mine had 0*getbinerror while Kurt's had 1*. 
	  mPbPb_Matrix[i]->SetBinContent(x,y,0);
	  mPbPb_Matrix[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_Matrix[i]->GetBinContent(x,y);
      }
      
      for (int x=1;x<=mPbPb_Matrix[i]->GetNbinsX();x++) {	   
	double ratio = 1;
	// if (hGenSpectraCorr->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorr->GetBinContent(x);
	mPbPb_Matrix[i]->SetBinContent(x,y,mPbPb_Matrix[i]->GetBinContent(x,y)*ratio);
	mPbPb_Matrix[i]->SetBinError(x,y,mPbPb_Matrix[i]->GetBinError(x,y)*ratio);
      }
    }

    cout<<"passed segment 1"<<endl;
    
    //mPbPb_Matrix[i]->Smooth(0);
    // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  
    
    mPbPb_Response[i] = (TH2F*)mPbPb_Matrix[i]->Clone(Form("mPbPb_Response_cent%d",i));
    TH1F *hProj = (TH1F*)mPbPb_Response[i]->ProjectionY()->Clone(Form("hProj_cent%d",i));

    for (int y=1;y<=mPbPb_Response[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_Response[i]->GetNbinsX();x++) {
	if (mPbPb_Response[i]->GetBinContent(x,y)<=1*mPbPb_Response[i]->GetBinError(x,y)) {
	  // in the above if loop, kurt has 1*error and my old had 0*error
	  mPbPb_Response[i]->SetBinContent(x,y,0);
	  mPbPb_Response[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_Response[i]->GetBinContent(x,y);
      }
      
      for (int x=1;x<=mPbPb_Response[i]->GetNbinsX();x++) {  	
	if (sum==0) continue;
	double ratio = 1;
	//if(dPbPb_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
	// else ratio = dPbPb_TrgComb[i]->GetBinContent(y)/sum
	ratio = 1./sum;
	if (hProj->GetBinContent(y)==0) ratio = 1e-100/sum;
	else ratio = hProj->GetBinContent(y)/sum;
	mPbPb_Response[i]->SetBinContent(x,y,mPbPb_Response[i]->GetBinContent(x,y)*ratio);
	mPbPb_Response[i]->SetBinError(x,y,mPbPb_Response[i]->GetBinError(x,y)*ratio);
      }
    }
    cout<<"passed segment 2"<<endl;
    
    mPbPb_ResponseNorm[i] = (TH2F*)mPbPb_Matrix[i]->Clone(Form("mPbPb_ResponseNorm_cent%d",i));
    for (int x=1;x<=mPbPb_ResponseNorm[i]->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=mPbPb_ResponseNorm[i]->GetNbinsY();y++) {
	if (mPbPb_ResponseNorm[i]->GetBinContent(x,y)<=1*mPbPb_ResponseNorm[i]->GetBinError(x,y)) {
	  mPbPb_ResponseNorm[i]->SetBinContent(x,y,0);
	  mPbPb_ResponseNorm[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_ResponseNorm[i]->GetBinContent(x,y);
      }
      
      for (int y=1;y<=mPbPb_ResponseNorm[i]->GetNbinsY();y++) {  	
	if (sum==0) continue;
	double ratio = 1./sum;
	mPbPb_ResponseNorm[i]->SetBinContent(x,y,mPbPb_ResponseNorm[i]->GetBinContent(x,y)*ratio);
	mPbPb_ResponseNorm[i]->SetBinError(x,y,mPbPb_ResponseNorm[i]->GetBinError(x,y)*ratio);
      }
      
    }    
    cout<<"passed segment 3"<<endl;
    mPbPb_ResponseNorm[i]->Print("base");

    hGen_Projection[i] = (TH1F*)mPbPb_ResponseNorm[i]->ProjectionX();
    hGen_Projection[i]->SetName(Form("hGen_Projection_cent%d",i));
    // hGen_Projection[i]->Rebin(10);

    cout<<"Finished getting histograms for full untruncated response matrix"<<endl;

    mPbPb_ResponseNorm_recoTrunc[i] = (TH2F*)mPbPb_ResponseNorm[i]->Clone(Form("mPbPb_ResponseNorm_recoTrunc_cent%d",i));

    for(int x = 0; x < mPbPb_ResponseNorm_recoTrunc[i]->GetNbinsX(); ++x){

      for(int y = 0; y < mPbPb_ResponseNorm_recoTrunc[i]->GetNbinsY(); ++y){

	if(y < mPbPb_ResponseNorm[i]->GetYaxis()->FindBin(50) || y > mPbPb_ResponseNorm[i]->GetYaxis()->FindBin(350)){
	  mPbPb_ResponseNorm_recoTrunc[i]->SetBinContent(x, y, 0.0);
	  mPbPb_ResponseNorm_recoTrunc[i]->SetBinError(x, y, 0.0);
	}
	
      }
      
    }
    
    
    // for (int y=1;y<=mPbPb_Matrix_recoTrunc[i]->GetNbinsY();y++) {
    //   double sum=0;
    //   for (int x=1;x<=mPbPb_Matrix_recoTrunc[i]->GetNbinsX();x++) {
    // 	if (mPbPb_Matrix_recoTrunc[i]->GetBinContent(x,y)<=1*mPbPb_Matrix_recoTrunc[i]->GetBinError(x,y)) {
    // 	  //in the above line mine had 0*getbinerror while Kurt's had 1*. 
    // 	  mPbPb_Matrix_recoTrunc[i]->SetBinContent(x,y,0);
    // 	  mPbPb_Matrix_recoTrunc[i]->SetBinError(x,y,0);
    // 	}
    // 	sum+=mPbPb_Matrix_recoTrunc[i]->GetBinContent(x,y);
    //   }
      
    //   for (int x=1;x<=mPbPb_Matrix_recoTrunc[i]->GetNbinsX();x++) {	   
    // 	double ratio = 1;
    // 	// if (hGenSpectraCorr->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorr->GetBinContent(x);
    // 	mPbPb_Matrix_recoTrunc[i]->SetBinContent(x,y,mPbPb_Matrix_recoTrunc[i]->GetBinContent(x,y)*ratio);
    // 	mPbPb_Matrix_recoTrunc[i]->SetBinError(x,y,mPbPb_Matrix_recoTrunc[i]->GetBinError(x,y)*ratio);
    //   }
    // }
    // //mPbPb_Matrix_recoTrunc[i]->Smooth(0);
    // // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  
    // cout<<"passed segment 1"<<endl;

    // mPbPb_Response_recoTrunc[i] = (TH2F*)mPbPb_Matrix_recoTrunc[i]->Clone(Form("mPbPb_Response_recoTrunc_cent%d",i));
    // TH1F *hProj_recoTrunc = (TH1F*)mPbPb_Response[i]->ProjectionY()->Clone(Form("hProj_recoTrunc_cent%d",i));

    // for (int y=1;y<=mPbPb_Response_recoTrunc[i]->GetNbinsY();y++) {
    //   double sum=0;
    //   for (int x=1;x<=mPbPb_Response_recoTrunc[i]->GetNbinsX();x++) {
    // 	if (mPbPb_Response_recoTrunc[i]->GetBinContent(x,y)<=1*mPbPb_Response_recoTrunc[i]->GetBinError(x,y)) {
    // 	  // in the above if loop, kurt has 1*error and my old had 0*error
    // 	  mPbPb_Response_recoTrunc[i]->SetBinContent(x,y,0);
    // 	  mPbPb_Response_recoTrunc[i]->SetBinError(x,y,0);
    // 	}
    // 	sum+=mPbPb_Response_recoTrunc[i]->GetBinContent(x,y);
    //   }
      
    //   for (int x=1;x<=mPbPb_Response_recoTrunc[i]->GetNbinsX();x++) {  	
    // 	if (sum==0) continue;
    // 	double ratio = 1;
    // 	//if(dPbPb_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
    // 	// else ratio = dPbPb_TrgComb[i]->GetBinContent(y)/sum
    // 	ratio = 1./sum;
    // 	if (hProj_recoTrunc->GetBinContent(y)==0) ratio = 1e-100/sum;
    // 	else ratio = hProj_recoTrunc->GetBinContent(y)/sum;
    // 	mPbPb_Response_recoTrunc[i]->SetBinContent(x,y,mPbPb_Response_recoTrunc[i]->GetBinContent(x,y)*ratio);
    // 	mPbPb_Response_recoTrunc[i]->SetBinError(x,y,mPbPb_Response_recoTrunc[i]->GetBinError(x,y)*ratio);
    //   }
    // }
    
    // cout<<"passed segment 2"<<endl;

    // mPbPb_ResponseNorm_recoTrunc[i] = (TH2F*)mPbPb_Matrix_recoTrunc[i]->Clone(Form("mPbPb_ResponseNorm_recoTrunc_cent%d",i));
    // for (int x=1;x<=mPbPb_ResponseNorm_recoTrunc[i]->GetNbinsX();x++) {
    //   double sum=0;
    //   for (int y=1;y<=mPbPb_ResponseNorm_recoTrunc[i]->GetNbinsY();y++) {
    // 	if (mPbPb_ResponseNorm_recoTrunc[i]->GetBinContent(x,y)<=1*mPbPb_ResponseNorm_recoTrunc[i]->GetBinError(x,y)) {
    // 	  mPbPb_ResponseNorm_recoTrunc[i]->SetBinContent(x,y,0);
    // 	  mPbPb_ResponseNorm_recoTrunc[i]->SetBinError(x,y,0);
    // 	}
    // 	sum+=mPbPb_ResponseNorm_recoTrunc[i]->GetBinContent(x,y);
    //   }
      
    //   for (int y=1;y<=mPbPb_ResponseNorm_recoTrunc[i]->GetNbinsY();y++) {  	
    // 	if (sum==0) continue;
    // 	double ratio = 1./sum;
    // 	mPbPb_ResponseNorm_recoTrunc[i]->SetBinContent(x,y,mPbPb_ResponseNorm_recoTrunc[i]->GetBinContent(x,y)*ratio);
    // 	mPbPb_ResponseNorm_recoTrunc[i]->SetBinError(x,y,mPbPb_ResponseNorm_recoTrunc[i]->GetBinError(x,y)*ratio);
    //   }
      
    // }    
    // cout<<"passed segment 3"<<endl;

    // mPbPb_ResponseNorm_recoTrunc[i]->Print("base");
 
    hGen_Projection_recoTrunc[i] = (TH1F*)mPbPb_ResponseNorm_recoTrunc[i]->ProjectionX();
    hGen_Projection_recoTrunc[i]->SetName(Form("hGen_Projection_recoTrunc_cent%d",i));
    // hGen_Projection_recoTrunc[i]->Rebin(10);

    cout<<"finished getting histograms for truncted response matrix"<<endl;
    
    hKineticEfficiency[i] = (TH1F*)hGen_Projection_recoTrunc[i]->Clone(Form("hKineticEfficiency_cent%d",i));
    hKineticEfficiency[i]->Divide(hGen_Projection[i]);

    for(int bin = 1; bin<=hKineticEfficiency[i]->GetNbinsX(); ++bin){
      double val = hKineticEfficiency[i]->GetBinContent(bin);
      val = (double)1./val;
      hKineticEfficiency[i]->SetBinContent(bin, val);
    }

    cout<<"Plotting kinetic efficiency"<<endl;
    
    cKineticEfficiency->cd(nbins_cent - i);

    hKineticEfficiency[i]->SetTitle(" ");
    hKineticEfficiency[i]->SetXTitle("Gen p_{T} (GeV/c)");
    hKineticEfficiency[i]->SetYTitle("Kinetic Efficiency");

    hKineticEfficiency[i]->SetMarkerStyle(20);
    hKineticEfficiency[i]->SetMarkerColor(kBlack);
    hKineticEfficiency[i]->SetAxisRange(60, 300, "X");
    hKineticEfficiency[i]->SetAxisRange(0.5, 1.5, "Y");
    hKineticEfficiency[i]->Draw();

    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);    
    if(i == nbins_cent-1) drawText("PYTHIA+HYDJET, Kinetic Efficiency", 0.25,0.7,14);
    if(i == nbins_cent-1) drawText(Form("akPu%dPF Jets",radius), 0.25, 0.65, 14);
    
    cout<<"Finished plotting kinetic efficiency"<<endl;

    cGenProj->cd(nbins_cent - i);
    hGen_Projection[i]->SetTitle(" ");
    hGen_Projection[i]->SetXTitle("Gen p_{T} (GeV/c)");
    hGen_Projection[i]->SetYTitle("GenProjection of Normalized Response Matrix");
    hGen_Projection[i]->Scale(1./10);
    hGen_Projection[i]->SetMarkerStyle(33);
    hGen_Projection[i]->SetMarkerColor(kBlack);
    hGen_Projection[i]->SetAxisRange(20, 350, "X");
    hGen_Projection[i]->SetAxisRange(0.6, 1.4, "Y");
    hGen_Projection[i]->Draw();

    hGen_Projection_recoTrunc[i]->Scale(1./10);
    hGen_Projection_recoTrunc[i]->SetMarkerStyle(25);
    hGen_Projection_recoTrunc[i]->SetMarkerColor(kRed);
    hGen_Projection_recoTrunc[i]->Draw("same");

    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);    
    if(i == nbins_cent-1) drawText("PYTHIA+HYDJET, Projection onto Gen axis", 0.25,0.7,14);
    if(i == nbins_cent-1) drawText(Form("akPu%dPF Jets",radius), 0.25, 0.65, 14);

    cMatrix->cd(nbins_cent-i);
    // mPbPb_ResponseNorm[i]->Rebin2D(10,10);
    // mPbPb_ResponseNorm[i]->Scale(1./TMath::Sqrt(10));
    makeHistTitle(mPbPb_ResponseNorm[i]," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");
    mPbPb_ResponseNorm[i]->SetAxisRange(0, 500, "X");
    mPbPb_ResponseNorm[i]->SetAxisRange(0, 500, "Y");
    mPbPb_ResponseNorm[i]->SetAxisRange(0, 1, "Z");
    mPbPb_ResponseNorm[i]->Draw("colz");
    
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.2,0.85,20);    
    if(i == nbins_cent-1) drawText("PYTHIA+HYDJET", 0.2,0.8,14);
    if(i == nbins_cent-1) drawText(Form("akPu%dPF Jets",radius), 0.2, 0.75, 14);

    cMatrix_recoTrunc->cd(nbins_cent-i);
    // mPbPb_ResponseNorm_recoTrunc[i]->Rebin2D(10,10);
    // mPbPb_ResponseNorm_recoTrunc[i]->Scale(1./TMath::Sqrt(10));
    makeHistTitle(mPbPb_ResponseNorm_recoTrunc[i]," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");
    mPbPb_ResponseNorm_recoTrunc[i]->SetAxisRange(0, 500, "X");
    mPbPb_ResponseNorm_recoTrunc[i]->SetAxisRange(0, 500, "Y");
    mPbPb_ResponseNorm_recoTrunc[i]->SetAxisRange(0, 1, "Z");
    mPbPb_ResponseNorm_recoTrunc[i]->Draw("colz");
    
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.2,0.85,20);    
    if(i == nbins_cent-1) drawText("PYTHIA+HYDJET, reco p_{T} truncated", 0.2,0.8,14);
    if(i == nbins_cent-1) drawText(Form("akPu%dPF Jets",radius), 0.2, 0.75, 14);

    
  }
  cGenProj->cd(1);
  TLegend * leg = myLegend(0.3,0.6,0.7,0.7);
  leg->AddEntry(hGen_Projection[0],"No Truncation in Reco pT","pl");
  leg->AddEntry(hGen_Projection_recoTrunc[0],"reco pT truncated 50 to 350 GeV","pl");
  leg->SetTextSize(0.04);
  leg->Draw();
  
  cKineticEfficiency->SaveAs(Form("KineticEfficiency_R%d_PbPbMC.pdf",radius),"RECREATE");
  cGenProj->SaveAs(Form("GenProjection_R%d_PbPbMC.pdf",radius),"RECREATE");
  cMatrix->SaveAs(Form("ResponseMatrix_R%d_PbPbMC.pdf",radius),"RECREATE");
  cMatrix_recoTrunc->SaveAs(Form("ResponseMatrix_recoTrunc_R%d_PbPbMC.pdf",radius),"RECREATE");

  fout.Write();
  fout.Close();
  
}
