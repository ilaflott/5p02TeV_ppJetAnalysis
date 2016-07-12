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

void analyze_deriveResponseOutput(){

  const bool debugMode=true;
  //const bool debugMode=false;
  
  std::cout<<"nbins_pt = " << nbins_pt << std::endl;
  std::cout<<"nbins_eta = " << nbins_eta << std::endl;

  // histos in root file named like
  // avg{A,B,AbsPho}Histo   _pti_etaj;
  // h_avg{A,B,AbsPho}      _i_j;
  // h_nEntries{A,B,AbsPho} _i_j;  
  std::string filename="deriveResponse_defOut.root";
  std::cout<< "filename is " << filename <<std::endl;
  TFile *fin = new TFile(filename.c_str());
  
  // avgAhistos loop
  std::cout<< "looping over avg A histos..."<<std::endl;
  TH1F *avgAHist[nbins_pt][nbins_eta];
  TH1D *h_avgA[nbins_pt][nbins_eta];
  TH1D *h_nEntriesA[nbins_pt][nbins_eta];
  int AHistCount_tot=0;
  int zeroAHistCount_tot=0;
  int zeroAHistCount_avgtot=0;
  int zeroAHistCount_Ntot=0;
  for(int i=0; i<nbins_pt ; ++i){//pt loop
    for(int j=0; j<nbins_eta; ++j){//eta loop
      std::string avgAHistTitle   ="avgAHisto_pt"+ std::to_string(i)+ "_eta"+ std::to_string(j);
      std::string h_avgATitle     ="h_avgA_"+      std::to_string(i)+ "_"+    std::to_string(j);
      std::string h_nEntriesTitle ="h_nEntriesA_"+ std::to_string(i)+ "_"+    std::to_string(j);
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
	std::cout<<"looking at avgA histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
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
	std::cout<<"looking at avgA histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
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
	std::cout<<"looking at avg A histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
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
  } // end avgAhistos loop

  // avgBhistos loop
  std::cout<< "looping over avg B histos..."<<std::endl;
  TH1F *avgBHist[nbins_pt][nbins_eta];
  TH1D *h_avgB[nbins_pt][nbins_eta];
  TH1D *h_nEntriesB[nbins_pt][nbins_eta];
  int BHistCount_tot=0;
  int zeroBHistCount_tot=0;
  int zeroBHistCount_avgtot=0;
  int zeroBHistCount_Ntot=0;
  for(int i=0; i<nbins_pt ; ++i){//pt loop
    for(int j=0; j<nbins_eta; ++j){//eta loop
      std::string avgBHistTitle   ="avgBHisto_pt"+ std::to_string(i)+ "_eta"+ std::to_string(j);
      std::string h_avgBTitle     ="h_avgB_"+      std::to_string(i)+ "_"+    std::to_string(j);
      std::string h_nEntriesTitle ="h_nEntriesB_"+ std::to_string(i)+ "_"+    std::to_string(j);
      avgBHist[i][j]    =(TH1F*)fin->Get( avgBHistTitle.c_str() );
      h_avgB[i][j]      =(TH1D*)fin->Get( h_avgBTitle.c_str() );
      h_nEntriesB[i][j] =(TH1D*)fin->Get( h_nEntriesTitle.c_str() );      
      BHistCount_tot++;
      int N_avgBHisto  =avgBHist[i][j]    ->GetEntries();
      int N_h_avgB     =h_avgB[i][j]      ->GetEntries();
      int N_h_nEntriesB=h_nEntriesB[i][j] ->GetEntries();
      double avg_avgBHisto  =avgBHist[i][j]    ->GetMean();
      double avg_h_avgB     =h_avgB[i][j]      ->GetMean();
      double avg_h_nEntriesB=h_nEntriesB[i][j] ->GetMean();
      if( ((N_avgBHisto==0)||(N_h_avgB==0)||(N_h_nEntriesB==0)) && 
	  ((avg_avgBHisto==0.)||(avg_h_avgB==0.)||(avg_h_nEntriesB==0.)) ){
	std::cout<<"looking at avgB histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "no entries and a mean of zero!" <<std::endl;
	zeroBHistCount_tot++;
	zeroBHistCount_Ntot++;
	zeroBHistCount_avgtot++;
	if(debugMode)std::cout<< "avgBHistTitle   is " << avgBHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgBTitle     is " << h_avgBTitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Entries in avgBHisto   " << N_avgBHisto  <<std::endl;
	if(debugMode)std::cout<<"Entries in h_avgB      " << N_h_avgB     <<std::endl;
	if(debugMode)std::cout<<"Entries in h_nEntriesB " << N_h_nEntriesB<<std::endl;
	if(debugMode)std::cout<<"Mean of avgBHisto   " << avg_avgBHisto  <<std::endl;
	if(debugMode)std::cout<<"Mean of h_avgB      " << avg_h_avgB     <<std::endl;
	if(debugMode)std::cout<<"Mean of h_nEntriesB " << avg_h_nEntriesB<<std::endl;
	std::cout<<std::endl;
	continue;
      }
      if( (N_avgBHisto) ==0 || (N_h_avgB) ==0 || (N_h_nEntriesB) ==0 ){
	std::cout<<"looking at avgB histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "no entries!" <<std::endl;
	zeroBHistCount_Ntot++;
	if(debugMode)std::cout<< "avgBHistTitle   is " << avgBHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgBTitle     is " << h_avgBTitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Entries in avgBHisto   " << N_avgBHisto  <<std::endl;
	if(debugMode)std::cout<<"Entries in h_avgB      " << N_h_avgB     <<std::endl;
	if(debugMode)std::cout<<"Entries in h_nEntriesB " << N_h_nEntriesB<<std::endl;
	std::cout<<std::endl;
      }
      if( (avg_avgBHisto) ==0 || (avg_h_avgB) ==0 || (avg_h_nEntriesB) ==0 ){
	zeroBHistCount_avgtot++;
	std::cout<<"looking at avgB histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "zero mean!" <<std::endl;
	if(debugMode)std::cout<< "avgBHistTitle   is " << avgBHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgBTitle     is " << h_avgBTitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Mean of avgBHisto   " << avg_avgBHisto  <<std::endl;
	if(debugMode)std::cout<<"Mean of h_avgB      " << avg_h_avgB     <<std::endl;
	if(debugMode)std::cout<<"Mean of h_nEntriesB " << avg_h_nEntriesB<<std::endl;
	std::cout<<std::endl;
      }
    }
  }// end avgBhistos loop
  
  // avgAbsPhohistos loop
  std::cout<< "looping over avg AbsPho histos..."<<std::endl;
  TH1F *avgAbsPhoHist[nbins_pt][nbins_eta];
  TH1D *h_avgAbsPho[nbins_pt][nbins_eta];
  TH1D *h_nEntriesAbsPho[nbins_pt][nbins_eta];
  int AbsPhoHistCount_tot=0;
  int zeroAbsPhoHistCount_tot=0;
  int zeroAbsPhoHistCount_avgtot=0;
  int zeroAbsPhoHistCount_Ntot=0;
  for(int i=0; i<nbins_pt ; ++i){//pt loop
    for(int j=0; j<nbins_eta; ++j){//eta loop
      std::string avgAbsPhoHistTitle   ="avgAbsPhoHisto_pt"+ std::to_string(i)+ "_eta"+ std::to_string(j);
      std::string h_avgAbsPhoTitle     ="h_avgAbsPho_"+      std::to_string(i)+ "_"+    std::to_string(j);
      std::string h_nEntriesTitle ="h_nEntriesAbsPho_"+ std::to_string(i)+ "_"+    std::to_string(j);
      avgAbsPhoHist[i][j]    =(TH1F*)fin->Get( avgAbsPhoHistTitle.c_str() );
      h_avgAbsPho[i][j]      =(TH1D*)fin->Get( h_avgAbsPhoTitle.c_str() );
      h_nEntriesAbsPho[i][j] =(TH1D*)fin->Get( h_nEntriesTitle.c_str() );      
      AbsPhoHistCount_tot++;
      int N_avgAbsPhoHisto  =avgAbsPhoHist[i][j]    ->GetEntries();
      int N_h_avgAbsPho     =h_avgAbsPho[i][j]      ->GetEntries();
      int N_h_nEntriesAbsPho=h_nEntriesAbsPho[i][j] ->GetEntries();
      double avg_avgAbsPhoHisto  =avgAbsPhoHist[i][j]    ->GetMean();
      double avg_h_avgAbsPho     =h_avgAbsPho[i][j]      ->GetMean();
      double avg_h_nEntriesAbsPho=h_nEntriesAbsPho[i][j] ->GetMean();
      if( ((N_avgAbsPhoHisto==0)||(N_h_avgAbsPho==0)||(N_h_nEntriesAbsPho==0)) && 
	  ((avg_avgAbsPhoHisto==0.)||(avg_h_avgAbsPho==0.)||(avg_h_nEntriesAbsPho==0.)) ){
	std::cout<<"looking at avgAbsPho histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "no entries and a mean of zero!" <<std::endl;
	zeroAbsPhoHistCount_tot++;
	zeroAbsPhoHistCount_Ntot++;
	zeroAbsPhoHistCount_avgtot++;
	if(debugMode)std::cout<< "avgAbsPhoHistTitle   is " << avgAbsPhoHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgAbsPhoTitle     is " << h_avgAbsPhoTitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Entries in avgAbsPhoHisto   " << N_avgAbsPhoHisto  <<std::endl;
	if(debugMode)std::cout<<"Entries in h_avgAbsPho      " << N_h_avgAbsPho     <<std::endl;
	if(debugMode)std::cout<<"Entries in h_nEntriesAbsPho " << N_h_nEntriesAbsPho<<std::endl;
	if(debugMode)std::cout<<"Mean of avgAbsPhoHisto   " << avg_avgAbsPhoHisto  <<std::endl;
	if(debugMode)std::cout<<"Mean of h_avgAbsPho      " << avg_h_avgAbsPho     <<std::endl;
	if(debugMode)std::cout<<"Mean of h_nEntriesAbsPho " << avg_h_nEntriesAbsPho<<std::endl;
	std::cout<<std::endl;
	continue;
      }
      if( (N_avgAbsPhoHisto) ==0 || (N_h_avgAbsPho) ==0 || (N_h_nEntriesAbsPho) ==0 ){
	std::cout<<"looking at avgAbsPho histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "no entries!" <<std::endl;
	zeroAbsPhoHistCount_Ntot++;
	if(debugMode)std::cout<< "avgAbsPhoHistTitle   is " << avgAbsPhoHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgAbsPhoTitle     is " << h_avgAbsPhoTitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Entries in avgAbsPhoHisto   " << N_avgAbsPhoHisto  <<std::endl;
	if(debugMode)std::cout<<"Entries in h_avgAbsPho      " << N_h_avgAbsPho     <<std::endl;
	if(debugMode)std::cout<<"Entries in h_nEntriesAbsPho " << N_h_nEntriesAbsPho<<std::endl;
	std::cout<<std::endl;
      }
      if( (avg_avgAbsPhoHisto) ==0 || (avg_h_avgAbsPho) ==0 || (avg_h_nEntriesAbsPho) ==0 ){
	zeroAbsPhoHistCount_avgtot++;
	std::cout<<"looking at avgAbsPho histo in bins [ptBin][etaBin] = "<<"["<<i<<"]["<<j<<"]"<<std::endl;
	std::cout<< "zero mean!" <<std::endl;
	if(debugMode)std::cout<< "avgAbsPhoHistTitle   is " << avgAbsPhoHistTitle   <<std::endl;
	if(debugMode)std::cout<< "h_avgAbsPhoTitle     is " << h_avgAbsPhoTitle     <<std::endl;
	if(debugMode)std::cout<< "h_nEntriesTitle is " << h_nEntriesTitle <<std::endl;
	if(debugMode)std::cout<<"Mean of avgAbsPhoHisto   " << avg_avgAbsPhoHisto  <<std::endl;
	if(debugMode)std::cout<<"Mean of h_avgAbsPho      " << avg_h_avgAbsPho     <<std::endl;
	if(debugMode)std::cout<<"Mean of h_nEntriesAbsPho " << avg_h_nEntriesAbsPho<<std::endl;
	std::cout<<std::endl;
      }
    }
  }


  std::cout<<"number of avg A histos   " << AHistCount_tot  <<std::endl;
  std::cout<<"number of avg A histos w/ zero entries is " << zeroAHistCount_Ntot  <<std::endl;
  std::cout<<"number of avg A histos w/ zero mean is "<< zeroAHistCount_avgtot <<std::endl;
  std::cout<<"number of avg A histos w/ both above traits is is "<< zeroAHistCount_tot <<std::endl;

  std::cout<<"number of avg B histos   " << BHistCount_tot  <<std::endl;
  std::cout<<"number of avg B histos w/ zero entries is " << zeroBHistCount_Ntot  <<std::endl;
  std::cout<<"number of avg B histos w/ zero mean is "<< zeroBHistCount_avgtot <<std::endl;
  std::cout<<"number of avg B histos w/ both above traits is is "<< zeroBHistCount_tot <<std::endl;

  std::cout<<"number of avg AbsPho histos   " << AbsPhoHistCount_tot  <<std::endl;
  std::cout<<"number of avg AbsPho histos w/ zero entries is " << zeroAbsPhoHistCount_Ntot  <<std::endl;
  std::cout<<"number of avg AbsPho histos w/ zero mean is "<< zeroAbsPhoHistCount_avgtot <<std::endl;
  std::cout<<"number of avg AbsPho histos w/ both above traits is is "<< zeroAbsPhoHistCount_tot <<std::endl;

  fin->Close();
  return;
}
