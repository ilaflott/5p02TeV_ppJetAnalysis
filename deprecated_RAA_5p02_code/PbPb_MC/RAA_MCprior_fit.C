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
#include "../../Headers/utilities.h"

using namespace std;

void RAA_MCprior_fit(int radius = 2){
    
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  TFile * fin = TFile::Open(Form("PbPb_MC_histograms_FromForest_akPu%d_20_eta_20.root",radius));

  TFile * fout = new TFile(Form("response_fit_weight_R%d.root",radius),"RECREATE");
  fout->cd();
  TH2F * hResponse[nbins_cent];
  TH1F * hPrior[nbins_cent];
  TF1 * fPowerLaw[nbins_cent];
  TH1F * hHistPrior_fit[nbins_cent];
  float fitStart = 40.0;
  float fitEnd = 500.0;

  for(int i = 0; i<nbins_cent; ++i){

    fPowerLaw[i] = new TF1(Form("fPowerLaw_cent%d",i),"[0]*pow(x+[2],[1])");
    fPowerLaw[i]->SetParameters(1e2,-5,0);
    
    hResponse[i] = (TH2F*)fin->Get(Form("hpbpb_Trans_matrix_HLT_R%d_20_eta_20_cent%d",radius,i));
    hPrior[i] = (TH1F*)hResponse[i]->ProjectionY();
    hPrior[i]->SetName(Form("hPrior_cent%d",i));
    hPrior[i]->Fit(Form("fPowerLaw_cent%d",i),"","", fitStart, fitEnd);
    hPrior[i]->Fit(Form("fPowerLaw_cent%d",i),"","", fitStart, fitEnd);
    hPrior[i]->Fit(Form("fPowerLaw_cent%d",i),"","", fitStart, fitEnd);
    hPrior[i]->Fit(Form("fPowerLaw_cent%d",i),"","", fitStart, fitEnd);
    hPrior[i]->Fit(Form("fPowerLaw_cent%d",i),"","", fitStart, fitEnd);

    hHistPrior_fit[i] = new TH1F(Form("hHistPrior_fit_cent%d",i),"",1000, 0, 1000);
    hHistPrior_fit[i] = functionHist(fPowerLaw[i],hHistPrior_fit[i],Form("hHistPrior_fit_cent%d",i));
    hHistPrior_fit[i]->Divide(hPrior[i]);
    
  }

  fout->Write();
  fout->Close();

}
