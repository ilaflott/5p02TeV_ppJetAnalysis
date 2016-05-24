#include <iostream>
#include <string>

#include <TH1F.h>
#include <TH1D.h>
#include <TFile.h>

//#include <TStyle.h>
//#include <TLine.h>
//#include <TGraph.h>
//#include <TGraphAsymmErrors.h>
//#include <TCanvas.h>
//#include <TAxis.h>
//#include <TMath.h>
//#include <TLegend.h>
//#include <TLatex.h>
//#include <TMathText.h>
#include "DijetResponse.h"

void analyze_derive_DijetResponseOutput(){

  const bool debugMode=false;
  
  std::cout<<"nbins_pt = " << nbins_pt << std::endl;
  std::cout<<"nbins_eta = " << nbins_eta << std::endl;

  //open the file
  std::string filename="deriveResponseOut_1file.root";
  std::string filename="deriveResponseOut_files0.4_HighPtJet80.root";
  std::string filename="deriveResponseOut_files5.9_HighPtJet80.root";
  std::string filename="deriveResponseOut_hadd_files0.9_HighPtJet80.root";

  std::cout<< "filename is " << filename <<std::endl;
  TFile *fin = new TFile(filename.c_str());
  
  //named in root file like
  //TH1F avgAHisto   _pti_etaj;
  //TH1D h_avgA      _i_j;
  //TH1D h_nEntriesA _i_j;  
  TH1F *avgAHist[nbins_pt][nbins_eta];
  TH1D *h_avgA[nbins_pt][nbins_eta];
  TH1D *h_nEntriesA[nbins_pt][nbins_eta];
  
  //loop over histos, grab #entries, mean, etc.
  int AHistCount_tot=0;
  int zeroAHistCount_tot=0;
  int zeroAHistCount_avgtot=0;
  int zeroAHistCount_Ntot=0;

  std::cout<< "looping over histos..."<<std::endl;
  for(int i=0; i<nbins_pt ; ++i){//pt loop
    for(int j=0; j<nbins_eta; ++j){//eta loop
      
    

      std::string avgAHistTitle   ="avgAHisto_pt"+ std::to_string(i)+ "_eta"+ std::to_string(j);
      std::string h_avgATitle     ="h_avgA_"+      std::to_string(i)+ "_"+    std::to_string(j);
      std::string h_nEntriesTitle ="h_nEntriesA_"+ std::to_string(i)+ "_"+    std::to_string(j);
 
      //avgAHist[i][j]    =(TH1F*)fin->Get( avgAHistTitle.c_str() )->Clone( avgAHistTitle.c_str() );
      //h_avgA[i][j]      =(TH1D*)fin->Get( h_avgATitle.c_str() )->Clone( h_avgATitle.c_str() );
      //h_nEntriesA[i][j] =(TH1D*)fin->Get( h_nEntriesTitle.c_str() )->Clone( h_nEntriesTitle.c_str() );
      avgAHist[i][j]    =(TH1F*)fin->Get( avgAHistTitle.c_str() );
      h_avgA[i][j]      =(TH1D*)fin->Get( h_avgATitle.c_str() );
      h_nEntriesA[i][j] =(TH1D*)fin->Get( h_nEntriesTitle.c_str() );
      
      AHistCount_tot++;

      int N_avgAHisto  =avgAHist[i][j]    ->GetEntries();
      int N_h_avgA     =h_avgA[i][j]      ->GetEntries();
      int N_h_nEntriesA=h_nEntriesA[i][j] ->GetEntries();
      double avg_avgAHisto  =avgAHist[i][j]    ->GetMean();
      double avg_h_avgA     =h_avgA[i][j]      ->GetMean();
      double avg_h_nEntriesA=h_nEntriesA[i][j] ->GetMean();

      if( ((N_avgAHisto==0)||(N_h_avgA==0)||(N_h_nEntriesA==0)) && 
	  ((avg_avgAHisto==0.)||(avg_h_avgA==0.)||(avg_h_nEntriesA==0.)) ){
	std::cout<<"looking at histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "no entries and a mean of zero!" <<std::endl;
	zeroAHistCount_tot++;
	zeroAHistCount_Ntot++;
	zeroAHistCount_avgtot++;
	if(debugMode)std::cout<< "avgAHistTitle   is " << avgAHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgATitle     is " << h_avgATitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Entries in avgAHisto   " << N_avgAHisto  <<std::endl;
	if(debugMode)std::cout<<"Entries in h_avgA      " << N_h_avgA     <<std::endl;
	if(debugMode)std::cout<<"Entries in h_nEntriesA " << N_h_nEntriesA<<std::endl;
	if(debugMode)std::cout<<"Mean of avgAHisto   " << avg_avgAHisto  <<std::endl;
	if(debugMode)std::cout<<"Mean of h_avgA      " << avg_h_avgA     <<std::endl;
	if(debugMode)std::cout<<"Mean of h_nEntriesA " << avg_h_nEntriesA<<std::endl;
	std::cout<<std::endl;
	continue;
      }

      if( (N_avgAHisto) ==0 || (N_h_avgA) ==0 || (N_h_nEntriesA) ==0 ){
	std::cout<<"looking at histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "no entries!" <<std::endl;
	zeroAHistCount_Ntot++;
	if(debugMode)std::cout<< "avgAHistTitle   is " << avgAHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgATitle     is " << h_avgATitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Entries in avgAHisto   " << N_avgAHisto  <<std::endl;
	if(debugMode)std::cout<<"Entries in h_avgA      " << N_h_avgA     <<std::endl;
	if(debugMode)std::cout<<"Entries in h_nEntriesA " << N_h_nEntriesA<<std::endl;
	std::cout<<std::endl;
      }

      if( (avg_avgAHisto) ==0 || (avg_h_avgA) ==0 || (avg_h_nEntriesA) ==0 ){
	zeroAHistCount_avgtot++;
	std::cout<<"looking at histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "zero mean!" <<std::endl;
	if(debugMode)std::cout<< "avgAHistTitle   is " << avgAHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgATitle     is " << h_avgATitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Mean of avgAHisto   " << avg_avgAHisto  <<std::endl;
	if(debugMode)std::cout<<"Mean of h_avgA      " << avg_h_avgA     <<std::endl;
	if(debugMode)std::cout<<"Mean of h_nEntriesA " << avg_h_nEntriesA<<std::endl;
	std::cout<<std::endl;
      }
    }
  }

  std::cout<<"number of hists   " << AHistCount_tot  <<std::endl;
  std::cout<<"number of hists w/ zero entries is " << zeroAHistCount_Ntot  <<std::endl;
  std::cout<<"number of hists w/ zero mean is "<< zeroAHistCount_avgtot <<std::endl;
  std::cout<<"number of hists w/ both above traits is is "<< zeroAHistCount_tot <<std::endl;

  fin->Close();
  return;
}
