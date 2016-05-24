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

const double ptbins_ana[] = {50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
const int ptbins[] = {50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
const int nbins_ana = 14;

using namespace std;

void pp_JER_plot(int radius = 3){

  bool printDebug = true;
  
  TFile * fin = TFile::Open(Form("pp_MC_new_ak%d_20_eta_20.root",radius));

  TH1F * hJER[nbins_ana];
  
  TF1 * F1[nbins_ana];
  float sigma[nbins_ana];
  TCanvas * cJER = new TCanvas("cJER","",1200,1000);
  cJER->Divide(5,3);

  for(int bin = 1; bin < nbins_ana; ++bin){

    cJER->cd(bin);
    cJER->cd(bin)->SetLogy();

    F1[bin] = new TF1(Form("F1_bin%d",bin),"gaus", 0.4, 1.6);

    hJER[bin] = (TH1F*)fin->Get(Form("hJER_%d_pt_%d",ptbins[bin], ptbins[bin+1]));  
    hJER[bin]->SetTitle(Form("%d < Gen p_{T} < %d; |#eta|<2", ptbins[bin], ptbins[bin+1]));
    hJER[bin]->SetXTitle("p_{T}^{detector jet}/p_{T}^{particle jet}");
    hJER[bin]->Print("base");
    hJER[bin]->Fit(F1[bin],"R");
    hJER[bin]->Draw();
    drawText(Form("Mean = %f", F1[bin]->GetParameter(1)),0.15,0.75,14);
    drawText(Form("Sigma = %f", F1[bin]->GetParameter(2)),0.15,0.7,14);
    sigma[bin] = F1[bin]->GetParameter(2);
  }

  cJER->cd(14);
  drawText(Form("PYTHIA Z2, ak%dPF Jets",radius),0.3,0.5,18);
  
  cJER->SaveAs(Form("pp_JER_ptrange_R%d.pdf",radius),"RECREATE");

  // define the residual correction factors
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
  
  TH1F * hPP_JER_Uncert = new TH1F("hPP_JER_Uncert","",nbins_ana, ptbins_ana);

  for(int bin = 1; bin <= hPP_JER_Uncert->GetNbinsX(); ++bin){

    float pt = hPP_JER_Uncert->GetBinCenter(bin);
    float res_val = (float)1./fResidual->Eval(pt) + sigma[bin-1];

    hPP_JER_Uncert->SetBinContent(bin, res_val);
    
  }

  TCanvas * cPP_JER = new TCanvas("cPP_JER","",800,600);
  hPP_JER_Uncert->SetLineColor(kBlue);
  hPP_JER_Uncert->SetAxisRange(0.5, 1.5, "Y");
  hPP_JER_Uncert->SetAxisRange(60, 300, "X");
  hPP_JER_Uncert->SetTitle(" ");
  hPP_JER_Uncert->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_JER_Uncert->SetYTitle("JER Uncertainty (%)");
  hPP_JER_Uncert->Draw();
  drawText(Form("ak%dPF Jets, PYTHIA", radius), 0.2, 0.2, 16);
  cPP_JER->SaveAs(Form("PP_JER_Uncert_R%d.pdf",radius),"RECREATE");

  TFile fout(Form("PP_JER_Uncert_R%d.root",radius),"RECREATE");
  fout.cd();
  hPP_JER_Uncert->Write();
  fout.Close();

}
