
//Plotting macro to simply sum all the stuff in the deriveDijetResponse macro
//K. Jung, December 2015


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
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
#include <TRandom.h>
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


using namespace std;

void sumDijetResponse(std::string filename="relDijetResponse_data_calo.root", int isMC= 0){

  bool doDraw = false;

  const int nbins_pt = 5;
  double xbins_pt[nbins_pt+1] = {40,60,80,110,200,1000};
  double xbins_eta[] = {-5.191, -4.716, -4.363, -4.013,
			-3.664, -3.314, -2.964,
			-2.650, -2.322, -2.043, -1.830,
			-1.653, -1.479, -1.305,
			-1.131, -0.957, -0.783, -0.609,
			-0.435, -0.261, -0.087,
			0.000,
			0.087,  0.261,  0.435, 0.609,
			0.783,  0.957,  1.131, 
			1.305,  1.479,  1.653,  1.830,
			2.043,  2.322,  2.650, 
			2.964,  3.314,   3.664, 4.013,
			4.363,  4.716,  5.191};
  const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;
  
	
  TH1F *hRelResponse[nbins_pt];
  TH1F *hMPFResponse[nbins_pt];

  TH1F *hAbsPhoResponse[nbins_pt];
  TH1F *hMPFAbsPhoResponse[nbins_pt];

  TH1F *avgAHisto[nbins_pt][nbins_eta];
  TH1F *avgBHisto[nbins_pt][nbins_eta];
  TH1F *hAvgAbsPhoResponse[nbins_pt][nbins_eta];

  std::string type = "";
  if(isMC) type = "PYTHIA_CUETP8M1";
  else type = "Data";

  TFile *fin = new TFile(filename.c_str());

  for(int i=0; i<nbins_pt; i++){
    cout<<i<<endl;
    hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hAbsPhoResponse[i] = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFAbsPhoResponse[i] = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    for(int j=0; j<nbins_eta; j++){
      cout<<j<<endl;
      avgAHisto[i][j] = (TH1F*)fin->Get(Form("avgAHisto_pt%d_eta%d",i,j))->Clone(Form("avgAHisto_pt%d_eta%d",i,j));
      avgBHisto[i][j] = (TH1F*)fin->Get(Form("avgBHisto_pt%d_eta%d",i,j))->Clone(Form("avgBHisto_pt%d_eta%d",i,j));
      //hAvgAbsPhoResponse[i][j] = (TH1F*)fin->Get(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j))->Clone(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j));
    }
  }

  for(int i=0; i<nbins_pt; i++){
    for(int j=0; j<nbins_eta; j++){
      cout<<i<<" "<<j<<endl;
    
      hRelResponse[i]->SetBinContent(j+1,(1+avgAHisto[i][j]->GetMean())/(1-avgAHisto[i][j]->GetMean()));
      hRelResponse[i]->SetBinError(j+1,hRelResponse[i]->GetBinContent(j+1)*(1./sqrt(avgAHisto[i][j]->GetEntries())));

      hMPFResponse[i]->SetBinContent(j+1,(1+avgBHisto[i][j]->GetMean())/(1-avgBHisto[i][j]->GetMean()));
      hMPFResponse[i]->SetBinError(j+1,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(avgBHisto[i][j]->GetEntries())));

      //hMPFAbsPhoResponse[i]->SetBinContent(j+1,hAvgAbsPhoResponse[i][j]->GetMean());
      //hMPFAbsPhoResponse[i]->SetBinError(j+1,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(hAvgAbsPhoResponse[i][j]->GetEntries())));

      int totEntriesA=0, totEntriesB=0, totEntriesAbs=0;
      for(int j=0; j<nbins_eta; j++){ 
	totEntriesA+=avgAHisto[i][j]->GetEntries(); 
	totEntriesB+=avgBHisto[i][j]->GetEntries();
	//totEntriesAbs+=hAvgAbsPhoResponse[i][j]->GetEntries();
      }
      hRelResponse[i]->SetEntries(totEntriesA);
      hMPFResponse[i]->SetEntries(totEntriesB);
      //hMPFAbsPhoResponse[i]->SetEntries(totEntriesAbs);
    }
  }
  
  int color[5] = {1,2,4,8,20};
  TFile *fout = new TFile(Form("simplePlots_%disMC.root",isMC),"recreate");
  fout->cd();
  for(int i=0; i<nbins_pt; i++){
    cout<<i<<endl;
    hRelResponse[i]->SetMarkerColor(color[i]);
    hRelResponse[i]->SetLineColor(color[i]);
    hRelResponse[i]->SetTitle(Form("Rel, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hRelResponse[i]->Write();

    hMPFResponse[i]->SetMarkerColor(color[i]);
    hMPFResponse[i]->SetLineColor(color[i]);
    hMPFResponse[i]->SetMarkerStyle(21);
    hMPFResponse[i]->SetLineStyle(2);
    hRelResponse[i]->SetTitle(Form("MPF, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hMPFResponse[i]->Write();

    // hMPFAbsPhoResponse[i]->SetMarkerColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetLineColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetTitle(Form("MPF Abs, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hMPFAbsPhoResponse[i]->Write();
  }

  cout<<"going to draw histograms"<<endl;

  if(doDraw){

    TCanvas *c0 = new TCanvas("c0","",600,600);
    TLegend *l0 = new TLegend(0.4,0.8,0.7,0.9);
    hRelResponse[0]->Draw();
    l0->AddEntry(hRelResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
    for(int i=1; i<nbins_pt; i++){
      cout<<i<<endl;
      hRelResponse[i]->Draw("same");
      if(isMC) l0->AddEntry("",Form("%s",type.c_str()),"");
      l0->AddEntry(hRelResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    }
    l0->Draw("Same");
    cout<<"saving c0"<<endl;
    //c0->SaveAs("test.png");
    //sprintf("Relresponse_%s.pdf",type.c_str());
    //cout<<type.c_str()<<endl;
    c0->SaveAs("test.pdf");
    cout<<"after savng"<<endl;
    //Printf("c0 saved");
    TCanvas *c1 = new TCanvas("c1","",600,600);
    c1->cd();
    TLegend *l1 = new TLegend(0.4,0.8,0.7,0.9);
    hMPFResponse[0]->Draw();
    l1->AddEntry(hMPFResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
    for(int i=1; i<nbins_pt; i++){
      cout<<i<<endl;
      hMPFResponse[i]->Draw("same");
      if(isMC) l1->AddEntry("",Form("%s",type.c_str()),"");
      l1->AddEntry(hMPFResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    }
    l1->Draw("Same");
    cout<<"saving c1"<<endl;
    c1->SaveAs(Form("MPFresponse_%s.pdf",type.c_str()),"RECREATE");

    // TCanvas *c2 = new TCanvas("c2","",600,600);
    // c2->cd();
    // hMPFAbsPhoResponse[0]->Draw();
    // TLegend *l2 = new TLegend(0.4,0.8,0.7,0.9);
    // l2->AddEntry(hMPFAbsPhoResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
    // for(int i=1; i<nbins_pt; i++){
    //   cout<<i<<endl;
    //   hMPFAbsPhoResponse[i]->Draw("same");
    //   if(isMC) l1->AddEntry("",Form("%s",type.c_str()),"");
    //   l2->AddEntry(hMPFAbsPhoResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    // }
    // l2->Draw("Same");
    // cout<<"saving c2"<<endl;
    // c2->SaveAs(Form("MPFABSPhoresponse_%s.pdf",type.c_str()),"RECREATE");
		
  }

  fout->Close();


}
