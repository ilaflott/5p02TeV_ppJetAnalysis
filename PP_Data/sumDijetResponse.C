// Orginally Written by K. Jung
// December 2015
// Destroyed By Ian Laflotte
// May 1st 2016
// Plotting macro to simply sum all the stuff in the deriveDijetResponse macro
// usage: hadd deriveDijetResponse output files, feed hadd'd root file into this macro

//C includes			    
#include <cstdlib>		    
#include <stdio.h>		    
#include <iostream>		    
#include <fstream>		    
#include <sstream>		    
#include <cmath>		    
				    
//load/read/write		    
#include <TSystem.h>		    
#include <TProfile.h>		    
#include <TFile.h>		    
#include <TTree.h>		    
#include <TNtuple.h>		    
				    
//utilities			    
#include <TChain.h>		    
#include <TCut.h>		    
#include <TMath.h>		    
#include <TRandom.h>
#include <TRandom3.h>		    
#include <TStopwatch.h>		    
#include <TEventList.h>		    
				    
//plotting			    
#include <TStyle.h>		    
#include <TCanvas.h>		    
#include <TGraphErrors.h>	    
#include <TGraphAsymmErrors.h>	    
#include <TLegend.h>		    
#include <TLatex.h>		    
#include <TLine.h>		    

//histos				    
#include <TH1.h>		    
#include <TF1.h>		    
#include <TH2.h>		    
#include <TH3.h>		    
#include <TH1F.h>		    
#include <TH1F.h>		    
#include <TH2F.h>                   

using namespace std;

// takes in a file, and an isMC integer flag, why integer?
void sumDijetResponse(std::string filename="relDijetResponse_data_calo.root", int isMC= 0){

  // decide if histograms are to be drawn
  bool doDraw = false;

  //pt and eta binning arrays

  //  const int nbins_pt = 5;
  //  double xbins_pt[nbins_pt+1] = {40,60,80,110,200,1000};
  double xbins_pt[] = {40,60,80,110,200,1000};
  const int nbins_pt=sizeof(xbins_pt)/sizeof(double)-1;

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
  
  TH1F *hRelResponse[nbins_pt];// Rel
  TH1F *hMPFResponse[nbins_pt];// MPF

  TH1F *hAbsPhoResponse[nbins_pt]; // AbsPho
  TH1F *hMPFAbsPhoResponse[nbins_pt]; // MPFAbsPho
  TH1F *hAvgAbsPhoResponse[nbins_pt][nbins_eta]; //AvgAbsPho

  TH1F *avgAHisto[nbins_pt][nbins_eta]; // avg A
  TH1F *avgBHisto[nbins_pt][nbins_eta]; // avg B

  //indicate the type of thing we're reading
  std::string type = "";
  if(!isMC) type = "Data"; 
  else type = "PYTHIA_CUETP8M1";

  //open the file
  TFile *fin = new TFile(filename.c_str());
  
  //array initializing
  for(int i=0; i<nbins_pt; i++){ //pt loop initializes response histo arrays
    cout<<i<<endl;

    hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hAbsPhoResponse[i] = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFAbsPhoResponse[i] = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);

    for(int j=0; j<nbins_eta; j++){//eta loop gets+clones avgAB histos from the input file
      cout<<j<<endl;

      avgAHisto[i][j] = (TH1F*)fin->Get( Form("avgAHisto_pt%d_eta%d",i,j) )->Clone(Form("avgAHisto_pt%d_eta%d",i,j));
      avgBHisto[i][j] = (TH1F*)fin->Get( Form("avgBHisto_pt%d_eta%d",i,j) )->Clone(Form("avgBHisto_pt%d_eta%d",i,j));
      //hAvgAbsPhoResponse[i][j] = (TH1F*)fin->Get(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j))->Clone(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j));

    }//end eta loop
  }//end pt loop

  //
  for(int i=0; i<nbins_pt; i++){//pt loop
    for(int j=0; j<nbins_eta; j++){//eta loop
      cout<<i<<" "<<j<<endl;
    
      //Ahisto stuff
      hRelResponse[i]->SetBinContent( j+1, ( 1+avgAHisto[i][j]->GetMean() )/( 1-avgAHisto[i][j]->GetMean() ) );
      hRelResponse[i]->SetBinError( j+1 , (hRelResponse[i]->GetBinContent(j+1)) * (1./sqrt(avgAHisto[i][j]->GetEntries()) ) );

      //Bhisto stuff
      hMPFResponse[i]->SetBinContent( j+1 ,(1+avgBHisto[i][j]->GetMean())/(1-avgBHisto[i][j]->GetMean()));
      hMPFResponse[i]->SetBinError( j+1 ,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(avgBHisto[i][j]->GetEntries())));

      //PhoResponse stuff
      //hMPFAbsPhoResponse[i]->SetBinContent( j+1 ,hAvgAbsPhoResponse[i][j]->GetMean());
      //hMPFAbsPhoResponse[i]->SetBinError( j+1 ,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(hAvgAbsPhoResponse[i][j]->GetEntries())));

      //count entries in avg AB histos
      int totEntriesA=0, totEntriesB=0, totEntriesAbs=0;

      for(int j=0; j<nbins_eta; j++){//eta loop

	totEntriesA+=avgAHisto[i][j]->GetEntries(); 
	totEntriesB+=avgBHisto[i][j]->GetEntries();
	//totEntriesAbs+=hAvgAbsPhoResponse[i][j]->GetEntries();

      }//eta loop

      //set entries in Rel/MPF response histo arrays after filling because... why?
      hRelResponse[i]->SetEntries(totEntriesA);
      hMPFResponse[i]->SetEntries(totEntriesB);
      //hMPFAbsPhoResponse[i]->SetEntries(totEntriesAbs);

    }//end eta loop
  }//end pt loop
  


  //open an output file for output plots if desired
  TFile *fout = new TFile(Form("simplePlots_%disMC.root",isMC),"recreate");
  fout->cd();

  //color matrix for later
  int color[5] = {1,2,4,8,20};

  for(int i=0; i<nbins_pt; i++){//pt loop
    cout<<i<<endl;

    hRelResponse[i]->SetMarkerColor(color[i]);
    hRelResponse[i]->SetLineColor(color[i]);
    hRelResponse[i]->SetTitle(Form("Rel, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hRelResponse[i]->Write();

    hMPFResponse[i]->SetMarkerColor(color[i]);
    hMPFResponse[i]->SetLineColor(color[i]);
    hMPFResponse[i]->SetMarkerStyle(21);
    hMPFResponse[i]->SetLineStyle(2);
    //    hRelResponse[i]->SetTitle(Form("MPF, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hMPFResponse[i]->SetTitle(Form("MPF, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hMPFResponse[i]->Write();

    // hMPFAbsPhoResponse[i]->SetMarkerColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetLineColor(color[i]);
    // hMPFAbsPhoResponse[i]->SetTitle(Form("MPF Abs, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    // hMPFAbsPhoResponse[i]->Write();
  }//end pt loop


  if(!doDraw) cout<<"not going to draw histograms"<<endl;
  else{
    cout<<"going to draw histograms"<<endl;

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
    //
    // TLegend *l2 = new TLegend(0.4,0.8,0.7,0.9);
    // l2->AddEntry(hMPFAbsPhoResponse[0],Form("%g<p_{T}<%g",xbins_pt[0],xbins_pt[1]));
    //
    // for(int i=1; i<nbins_pt; i++){
    //   cout<<i<<endl;
    //   hMPFAbsPhoResponse[i]->Draw("same");
    //   if(isMC) l1->AddEntry("",Form("%s",type.c_str()),"");
    //   l2->AddEntry(hMPFAbsPhoResponse[i],Form("%g<p_{T}<%g",xbins_pt[i],xbins_pt[i+1]));
    // }
    //
    // l2->Draw("Same");
    // cout<<"saving c2"<<endl;
    // c2->SaveAs(Form("MPFABSPhoresponse_%s.pdf",type.c_str()),"RECREATE");
		
  }

  //close the output file
  fout->Close();


}
