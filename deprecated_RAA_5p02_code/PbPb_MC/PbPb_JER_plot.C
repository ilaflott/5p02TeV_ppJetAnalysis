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

const int ptbins_ana[] = {50, 80, 120, 170, 220, 300};
const int nbins_ana = sizeof(ptbins_ana)/sizeof(int) -1;
const int ptbins[] = {50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
const int nbins_cent = 6;
const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 

using namespace std;

void PbPb_JER_plot(int radius = 4){

  bool printDebug = true;
  
  TFile * fin = TFile::Open(Form("PbPb_5p02TeV_MC_akPu%dPF_20_eta_20.root",radius));

  TH1F * hJER[nbins_ana][nbins_cent];
  
  TF1 * F1[nbins_ana][nbins_cent];
  float sigma[nbins_ana][nbins_cent];
  TCanvas * cJER[nbins_cent];

  for(int i = 0; i<nbins_cent; ++i){
    
    cJER[i] = new TCanvas(Form("cJER_cent%d",i),"",1200,1000);
    cJER[i]->Divide(5,3);

    for(int bin = 1; bin < nbins_ana; ++bin){

      cJER[i]->cd(bin);
      cJER[i]->cd(bin)->SetLogy();

      F1[bin][i] = new TF1(Form("F1_bin%d_cent%d",bin,i),"gaus", 0.4, 1.6);

      hJER[bin][i] = (TH1F*)fin->Get(Form("hJER_%d_pt_%d_cent%d",ptbins[bin], ptbins[bin+1], i));  
      hJER[bin][i]->SetTitle(Form("%d < Gen p_{T} < %d; |#eta|<2", ptbins[bin], ptbins[bin+1]));
      hJER[bin][i]->SetXTitle("p_{T}^{detector jet}/p_{T}^{particle jet}");
      hJER[bin][i]->Print("base");
      hJER[bin][i]->Fit(F1[bin][i],"R");
      hJER[bin][i]->Draw();
      drawText(Form("Mean = %f", F1[bin][i]->GetParameter(1)),0.15,0.75,14);
      drawText(Form("Sigma = %f", F1[bin][i]->GetParameter(2)),0.15,0.7,14);
      sigma[bin][i] = F1[bin][i]->GetParameter(2);
    }

    cJER[i]->cd(14);
    drawText(Form("PYTHIA+HYDJET, akPu%dPF Jets, |#eta|<2",radius),0.3,0.5,18);
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);    
    
  
    cJER[i]->SaveAs(Form("PbPb_JER_ptrange_R%d_cent%d.pdf",radius, i),"RECREATE");
  }
  //define the residual correction factors
  TF1 * fResidual = new TF1("fResidual","1 - [0]/pow(x,[1])");

  if(radius == 2){
    fResidual->SetParameter(0, -0.05756);
    fResidual->SetParameter(1,  0.42750);
  }
  if(radius == 3){
    fResidual->SetParameter(0, -0.66229);
    fResidual->SetParameter(1,  1.02119);
  }
  if(radius == 4){
    fResidual->SetParameter(0, -1.28178);
    fResidual->SetParameter(1,  1.17348);
  }

  TFile fout(Form("PbPb_JER_Uncert_R%d.root",radius),"RECREATE");
  fout.cd();

  
  TH1F * hPbPb_JER_Uncert[nbins_cent];
  TCanvas * cPbPb_JER;
  cPbPb_JER = new TCanvas("cPbPb_JER","",800,600);
  cPbPb_JER->Divide(3,2);
  
  for(int i = 0; i<nbins_cent; ++i){
    hPbPb_JER_Uncert[i] = new TH1F(Form("hPbPb_JER_Uncert_cent%d",i),"",nbins_ana, ptbins_ana);

    for(int bin = 1; bin <= hPbPb_JER_Uncert[i]->GetNbinsX(); ++bin){
      
      float pt = hPbPb_JER_Uncert[i]->GetBinCenter(bin);
      float res_val = (float)1./fResidual->Eval(pt) + sigma[bin-1][i];
      
      hPbPb_JER_Uncert[i]->SetBinContent(bin, res_val);
      
    }
    hPbPb_JER_Uncert[i]->Write();

    cPbPb_JER->cd(nbins_cent-i);

    hPbPb_JER_Uncert[i]->SetLineColor(kBlue);
    hPbPb_JER_Uncert[i]->SetAxisRange(0.5, 1.5, "Y");
    hPbPb_JER_Uncert[i]->SetAxisRange(60, 300, "X");
    hPbPb_JER_Uncert[i]->SetTitle(" ");
    hPbPb_JER_Uncert[i]->SetXTitle("Jet p_{T} (GeV/c)");
    hPbPb_JER_Uncert[i]->SetYTitle("JER Uncertainty (%)");
    hPbPb_JER_Uncert[i]->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);

    if(i == 0) drawText(Form("akPu%dPF Jets PY+HD", radius), 0.2, 0.2, 16);
    
  }

  cPbPb_JER->SaveAs(Form("PbPb_JER_Uncert_R%d.pdf",radius),"RECREATE");
  
  fout.Close();

}
