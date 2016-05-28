/*****************
 * Macro for deriving the dijet relative and absolute response from data (or MC)
 * Uses JME-13-004 as a reference - "A" and "B" formulae are defined there
 * Kurt Jung, December 2015
 *******************/

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "TLorentzVector.h"

using namespace std;

int findEtaBin(float eta, int nbins_eta, double* xbins_eta){

  for(int i=0; i<nbins_eta; i++){
    if(eta<xbins_eta[i+1]) return i;
  }
  return nbins_eta-1;
}

int findPtBin(float pt, int nbins_pt, double* xbins_pt){

  for(int i=0; i<nbins_pt; i++){
    if(pt<xbins_pt[i+1]) return i;
  }
  return nbins_pt-1;
}

TLorentzVector findMissEt(int nPFpart, int* pfId, float* pfPt, float* pfEta, float* pfPhi, int nref, float* pt_F, float* rawpt_F, float* eta_F, float* phi_F, float* m_F, float* eSum, float *phoSum){
	
  TLorentzVector missEt(0,0,0,0);
  for(int i=0; i<nPFpart; i++){
    TLorentzVector pfp;
    double pfMass=0;
			
    if(pfId[i]==1 || pfId[i]==4){ pfMass = 0.1395702; }
    pfp.SetPtEtaPhiM(pfPt[i],pfEta[i],pfPhi[i],pfMass); 
    pfp.SetPz(0.); //2d projection in x-y plane 
    missEt += pfp;
  }
  //include jet residuals
  for(int i=0; i<nref; i++){
    if(pt_F[i]>15){
      TLorentzVector jtTmp(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
      if((eSum[i]+phoSum[i])/jtTmp.E() < 0.9){
	TLorentzVector jt;
	jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]); 
	jt.SetPz(0.); //project to x-y plane
	missEt += jt;
      }
    }
  }
  //rotate to get *missing* ET from residual ET
  missEt*=-1;
  missEt.SetE(0.);  //really only care about the eT (not full E)
  //cout << "missEt:  Mag: " << missEt.Mag() << endl;
  //missEt.Print();
  return missEt;
}

void deriveDijetResponse(int startfile=0,
			 int endfile=1,
			 int isMC=1,
			 int radius = 4,
			 std::string out="test.root")
{
  
  const double minJetPt = 100;
  
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
  
  for(int i=0; i<nbins_pt; i++){
    hRelResponse[i] = new TH1F(Form("hRelResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFResponse[i] = new TH1F(Form("hMPFResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hAbsPhoResponse[i] = new TH1F(Form("hAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    hMPFAbsPhoResponse[i] = new TH1F(Form("hMPFAbsPhoResponse_pt%d",i),"",nbins_eta,xbins_eta);
    for(int j=0; j<nbins_eta; j++){
      avgAHisto[i][j] = new TH1F(Form("avgAHisto_pt%d_eta%d",i,j),"",50,-1,1);
      avgBHisto[i][j] = new TH1F(Form("avgBHisto_pt%d_eta%d",i,j),"",50,-1,1);
      hAvgAbsPhoResponse[i][j] = new TH1F(Form("hAvgAbsPhoResponse_pt%d_eta%d",i,j),"",50,-1,1);
    }
  }

  double avgA[nbins_pt][nbins_eta];
  TH1D *h_avgA[nbins_pt][nbins_eta];
  int nEntries[nbins_pt][nbins_eta];
  TH1D *h_nEntries[nbins_pt][nbins_eta];
  double avgB[nbins_pt][nbins_eta];
  TH1D *h_avgB[nbins_pt][nbins_eta];
  int nEntriesB[nbins_pt][nbins_eta];
  TH1D *h_nEntriesB[nbins_pt][nbins_eta];
  double avgAbsPhoResponse[nbins_pt][nbins_eta];
  TH1D *h_avgAbsPhoResponse[nbins_pt][nbins_eta];
  int nEntriesAbs[nbins_pt][nbins_eta];
  TH1D *h_nEntriesAbs[nbins_pt][nbins_eta];
  
  for(int i=0; i<nbins_pt; i++){
    for(int j=0; j<nbins_eta; j++){
      avgA[i][j]=0;
      nEntries[i][j]=0;
      avgB[i][j]=0;
      nEntriesB[i][j]=0;
      avgAbsPhoResponse[i][j]=0;
      nEntriesAbs[i][j]=0;
      h_avgA[i][j]=new TH1D(Form("h_avgA_%d_%d",i,j),"",2, 0, 2);;
      h_nEntries[i][j]=new TH1D(Form("h_nEntries_%d_%d",i,j),"",2, 0, 2);;
      h_avgB[i][j]=new TH1D(Form("h_avgB_%d_%d",i,j),"",2, 0, 2);;
      h_nEntriesB[i][j]=new TH1D(Form("h_nEntriesB_%d_%d",i,j),"",2, 0, 2);;
      h_avgAbsPhoResponse[i][j]=new TH1D(Form("h_avgAbsPhoResponse_%d_%d",i,j),"",2, 0, 2);;
      h_nEntriesAbs[i][j]=new TH1D(Form("h_nEntriesAbs_%d_%d",i,j),"",2, 0, 2);;      
    }
  }
  
  std::string infile_Forest;
  if(isMC) infile_Forest = "pp_mc_forests_new.txt";
  else infile_Forest = "pp_data_forests.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;

  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  int totalEntries = 0;
  for(int ifile = startfile; ifile<endfile; ++ifile){
    
    //signal(SIGINT,break);		
    instr_Forest>>filename_Forest;
    cout << "opening file " << filename_Forest << endl;
    TFile *fin = TFile::Open(filename_Forest.c_str());
    TTree *jtTree = (TTree*)fin->Get(Form("ak%dPFJetAnalyzer/t",radius));
    TTree *pfTree = (TTree*)fin->Get("pfcandAnalyzer/pfTree");
    TTree *phoTree = (TTree*)fin->Get("ggHiNtuplizer/EventTree");

    float pt_F[1000], eta_F[1000], phi_F[1000], rawpt_F[1000], m_F[1000], eSum[1000], phoSum[1000];
    float pfPt[10000], pfEta[10000], pfPhi[10000];
    int nref, nPFpart, pfId[10000], nPho;
    vector<float> *phoEt=0, *phoEta=0, *phoPhi=0; 

    jtTree->SetBranchAddress("jtpt",pt_F);
    jtTree->SetBranchAddress("jteta",eta_F);
    jtTree->SetBranchAddress("jtphi",phi_F);
    jtTree->SetBranchAddress("jtm",m_F);
    jtTree->SetBranchAddress("rawpt",rawpt_F);
    jtTree->SetBranchAddress("nref",&nref);
    jtTree->SetBranchAddress("eSum",&eSum);
    jtTree->SetBranchAddress("photonSum",&phoSum);
    pfTree->SetBranchAddress("nPFpart",&nPFpart);
    pfTree->SetBranchAddress("pfPt",pfPt);
    pfTree->SetBranchAddress("pfEta",pfEta);
    pfTree->SetBranchAddress("pfPhi",pfPhi);
    pfTree->SetBranchAddress("pfId",pfId);
    phoTree->SetBranchAddress("nPho",&nPho);
    phoTree->SetBranchAddress("phoEt",&phoEt);
    phoTree->SetBranchAddress("phoEta",&phoEta);
    phoTree->SetBranchAddress("phoPhi",&phoPhi);

    cout << "file entries: "<< jtTree->GetEntries() << endl;

    for(int ientry=0; ientry<jtTree->GetEntries(); ientry++){
    //for(int ientry=0; ientry<2000000; ientry++){
      jtTree->GetEntry(ientry);
      pfTree->GetEntry(ientry);
      phoTree->GetEntry(ientry);

      //signal(SIGINT,break);
      totalEntries++;
      if(nPFpart > 10000 || nref > 1000) cout << " warning! nPF: "<< nPFpart << " and njet: "<< nref << endl;
      if(totalEntries && totalEntries%100000==0) cout << "entry: "<< ientry << endl;

      //first look for photon balance (AN2013-179 / JME-13-004)
      for(int i=0; i<nPho; i++){
	if(phoEt->at(i)>30 && abs(phoEta->at(i))<1.4442){
	  //cout << "photon found, Et: " << phoEt->at(i) << endl;
	  for(int ijet=0; ijet<nref-1; ijet++){
	    if(pt_F[ijet]<minJetPt) continue;
	    if(abs(phoPhi->at(i)-phi_F[ijet])>2.95 && pt_F[ijet+1]/phoEt->at(i)<0.2){
	      //cout << "jet photon pair found" << endl;

	      TLorentzVector missEt = findMissEt(nPFpart, pfId, pfPt, pfEta, pfPhi, nref, pt_F, rawpt_F, eta_F, phi_F, m_F, eSum, phoSum);

	      int etaBin = findEtaBin(phoEta->at(i),nbins_eta,xbins_eta);
	      int ptBin = findPtBin(phoEt->at(i),nbins_pt,xbins_pt);
	      TLorentzVector phoVec(phoEt->at(i),phoEta->at(i),phoPhi->at(i),0.);
	      double num = missEt.Dot(phoVec);
	      avgAbsPhoResponse[ptBin][etaBin] += (1+(num/phoVec.Dot(phoVec)));
	      h_avgAbsPhoResponse[ptBin][etaBin]->Fill((1+(num/phoVec.Dot(phoVec))));
	      hAvgAbsPhoResponse[ptBin][etaBin]->Fill(1+(num/phoVec.Dot(phoVec)));
	      nEntriesAbs[ptBin][etaBin]++;
	      h_nEntriesAbs[ptBin][etaBin]->Fill(1);
	    }
	  }
	}
      }

      //Doing standard Rrel method
      if(nref>1){
	int rJet, pJet;
	if(abs(eta_F[0])<1.3 && abs(eta_F[1])<1.3){
	  pJet = rand()%2;
	  rJet = (pJet+1)%2;
	}
	else if(abs(eta_F[0])<1.3 && abs(eta_F[1])<5.){
	  pJet = 1; rJet = 0;
	}
	else if(abs(eta_F[1])<1.3 && abs(eta_F[0])<5.){
	  pJet = 0; rJet = 1;
	}
	else continue;
	double avgPt = 0.5*(pt_F[rJet]+pt_F[pJet]);
	if(nref>2){
	  if(pt_F[2]/avgPt > 0.2){
	    continue;
	  }
	}
	if(pt_F[rJet]<minJetPt) continue;
	double dphi = abs(phi_F[rJet]-phi_F[pJet]);
	if(dphi>(2*3.14159)) dphi-=(2*3.14159);
	if(dphi < 2.7) continue;
	int etaBin = findEtaBin(eta_F[pJet],nbins_eta,xbins_eta);
	int ptBin = findPtBin(avgPt,nbins_pt,xbins_pt);
	//if(etaBin>xbins_eta[nbins_eta] || ptBin>xbins_pt[nbins_pt]) cout << "OH NO! Bin Mismatch!" << endl;
	avgA[ptBin][etaBin] += 0.5*(pt_F[pJet]-pt_F[rJet])/avgPt;
	h_avgA[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);
	avgAHisto[ptBin][etaBin]->Fill(0.5*(pt_F[pJet]-pt_F[rJet])/avgPt);
	nEntries[ptBin][etaBin]++;
	h_nEntries[ptBin][etaBin]->Fill(1);

	//Starting MPF Method
	TLorentzVector missEt = findMissEt(nPFpart, pfId, pfPt, pfEta, pfPhi, nref, pt_F, rawpt_F, eta_F, phi_F, m_F, eSum, phoSum);

	for(int i=0; i<nref; i++){
	  int etaBin = findEtaBin(eta_F[pJet],nbins_eta,xbins_eta);
	  int ptBin = findPtBin(avgPt,nbins_pt,xbins_pt);
	  TLorentzVector jetVec(pt_F[rJet],eta_F[rJet],phi_F[rJet],m_F[rJet]);
	  avgB[ptBin][etaBin] += (missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag()));
	  h_avgB[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	  avgBHisto[ptBin][etaBin]->Fill((missEt.Dot(jetVec) / ( 2*avgPt * jetVec.Mag())));
	  nEntriesB[ptBin][etaBin]++;
	  h_nEntriesB[ptBin][etaBin]->Fill(1);
	}
	
      }
    } //close entries
    fin->Close();
  } //close file loop

  for(int i=0; i<nbins_pt; i++){
    for(int j=0; j<nbins_eta; j++){
      avgA[i][j] = avgA[i][j]/(double)nEntries[i][j];
      avgB[i][j] = avgB[i][j]/(double)nEntriesB[i][j];
      if(nEntries[i][j]) hRelResponse[i]->SetBinContent(j+1,(1+avgA[i][j])/(1-avgA[i][j]));
      //if(nEntries[i][j]) hRelResponse[i]->SetBinError(j+1,hRelResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntries[i][j])));
      if(nEntries[i][j]) hRelResponse[i]->SetBinError(j+1,avgAHisto[i][j]->GetRMS()*(1./sqrt(nEntries[i][j])));
      else hRelResponse[i]->SetBinContent(j+1,0);
      
      if(nEntriesB[i][j]) hMPFResponse[i]->SetBinContent(j+1,(1+avgB[i][j])/(1-avgB[i][j]));
      //if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,hMPFResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesB[i][j])));
      if(nEntriesB[i][j]) hMPFResponse[i]->SetBinError(j+1,avgBHisto[i][j]->GetRMS()*(1./sqrt(nEntriesB[i][j])));
      else hMPFResponse[i]->SetBinContent(j+1,0);
   
      if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinContent(j+1,avgAbsPhoResponse[i][j]/nEntriesAbs[i][j]);
      //if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hMPFAbsPhoResponse[i]->GetBinContent(j+1)*(1./sqrt(nEntriesAbs[i][j])));
      if(nEntriesAbs[i][j]) hMPFAbsPhoResponse[i]->SetBinError(j+1,hAvgAbsPhoResponse[i][j]->GetRMS()*(1./sqrt(nEntriesAbs[i][j])));
      else hMPFAbsPhoResponse[i]->SetBinContent(j+1,0);

      int totEntriesA=0, totEntriesB=0, totEntriesAbs=0;
      for(int j=0; j<nbins_eta; j++){ 
  	totEntriesA+=nEntries[i][j]; 
  	totEntriesB+=nEntriesB[i][j];
  	totEntriesAbs+=nEntriesAbs[i][j];
      }
      hRelResponse[i]->SetEntries(totEntriesA);
      hMPFResponse[i]->SetEntries(totEntriesB);
      hMPFAbsPhoResponse[i]->SetEntries(totEntriesAbs);
    }
  }
	
  int color[5] = {1,2,4,8,20};
  TFile *fout;
  fout = new TFile(Form("%s",out.c_str()),"recreate");

  fout->cd();
  for(int i=0; i<nbins_pt; i++){
    hRelResponse[i]->SetMarkerColor(color[i]);
    hRelResponse[i]->SetLineColor(color[i]);
    hRelResponse[i]->SetTitle(Form("Rel, %g<p^{avg}_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hRelResponse[i]->Write();
    //hRelResponse[i]->Scale(1./hRelResponse[i]->Integral());

    hMPFResponse[i]->SetMarkerColor(color[i]);
    hMPFResponse[i]->SetLineColor(color[i]);
    hMPFResponse[i]->SetMarkerStyle(21);
    hMPFResponse[i]->SetLineStyle(2);
    hMPFResponse[i]->SetTitle(Form("MPF dijet, %g<p^{avg}_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hMPFResponse[i]->Write();
    //hMPFResponse[i]->Scale(1./hMPFResponse[i]->Integral());

    hMPFAbsPhoResponse[i]->SetMarkerColor(color[i]);
    hMPFAbsPhoResponse[i]->SetLineColor(color[i]);
    hMPFAbsPhoResponse[i]->SetTitle(Form("MPF #gamma+Jet, %g<p_{T}<%g GeV",xbins_pt[i],xbins_pt[i+1]));
    hMPFAbsPhoResponse[i]->Write();
    for(int j=0; j<nbins_eta; j++){
      avgAHisto[i][j]->Write();
      avgBHisto[i][j]->Write();
      hAvgAbsPhoResponse[i][j]->Write();
      h_avgA[i][j]->Write();
      h_nEntries[i][j]->Write();
      h_avgB[i][j]->Write();
      h_nEntriesB[i][j]->Write();
      h_avgAbsPhoResponse[i][j]->Write();
      h_nEntriesAbs[i][j]->Write();
    }
  }

  fout->Close();

}
