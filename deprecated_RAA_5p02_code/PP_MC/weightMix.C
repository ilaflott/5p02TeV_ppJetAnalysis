#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TChain.h"

#include "TCut.h"
#include "TNtuple.h"

#include "THStack.h"

using namespace std;

static bool subtract = 0;
static int QID = 3;

static TString weightString;

static bool normLead = 0;

static int mixColor = 2;
static int dataColor = 1;
static int ppColor = 4;

static int centralBin = 8;
static int leadCut = 120;
static int subleadCut = 30;

static double sideMin = 0.1;
static double sideMax = TMath::Pi()/3 + 0.1;

static double sideCorrect = 1;

static const char* LUM = "#int L dt=150";

static bool plotSubtraction = 0;

static bool reweightCentrality = 1;

static const double pi = TMath::Pi();
static const int Nfiles = 11;

void weightMix(){

  TH1::SetDefaultSumw2();

  int Npt = 11;
  //int Npt = 10;

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest2015#Dijet_Cross_Sections_for_reweigh
  double pthatBins2015[] = {15       , 30       , 50       , 80       , 120      , 170      , 220      , 280      , 370      , 460      , 540      , 9999};
  double xs2015[]        = {5.269E-01, 3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05, 1.223E-05, 3.031E-06, 7.746E-07, 1.410E-07, 3.216E-08, 1.001E-08, 0.0};
  //double pthatBins2015[] = {15       , 30       , 50       , 80       , 120      , 170      , 220      , 280      , 460      , 540      , 9999};
  //double xs2015[]        = {5.269E-01, 3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05, 1.223E-05, 3.031E-06, 7.746E-07, 3.216E-08, 1.001E-08, 0.0};
  
  double *xs;
  double *pthats;
  int n[11];

  Npt = 11;
  xs = xs2015;
  pthats = pthatBins2015;

  for(int i = 0; i < Npt-1; ++i){
    xs[i] -= xs[i+1];
  }

  TChain* nt;
  // TChain* evt;
  TFile* outf[Nfiles];

  nt = new TChain("ak4PFJetAnalyzer/t");
  //evt = new TChain("hiEvtAnalyzer/HiTree");

  std::string infile_Forest;
  infile_Forest = "pp_mc_forests_new.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  std::string filename_array[Nfiles];

  for(int ifile = 0; ifile<Nfiles; ++ifile){
    instr_Forest>>filename_Forest;
    nt->AddFile(filename_Forest.c_str());
    filename_array[ifile] = filename_Forest.c_str();
  }

  cout<<"loaded all the files"<<endl;
  
  for(int i = 0; i < Npt; ++i){
    TCut pthatCut(Form("pthat >= %f && pthat < %f",pthats[i],pthats[i+1]));
    n[i] = nt->GetEntries(pthatCut);
    cout<<"no of events in pthat = "<<pthats[i]<<" = "<<n[i]<<endl;
  }

  
  for(int i = 0; i<Nfiles; ++i){

    cout<<"Now looping over individual file : "<<filename_array[i].c_str()<<endl;

    TFile * fin = new TFile(filename_array[i].c_str());
    TTree * newTree = (TTree*)fin->Get("ak4PFJetAnalyzer/t");
    //TTree * evtTree = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");

    float vz, vz_n;
    double pthatweight = 0;
    int evnt, evnt_n, lumi, lumi_n;
    float varpthat, varvz;
  
    newTree->SetBranchAddress("pthat",&varpthat);
    //evtTree->SetBranchAddress("evt", &evnt);
    //evtTree->SetBranchAddress("lumi", &lumi);
    //evtTree->SetBranchAddress("vz", &vz);
    //newTree->AddFriend(evtTree);

    outf[i] = new TFile(Form("weights_pp_%d.root",i),"recreate");
    
    TTree *ntw = new TTree("weights","");
    ntw->Branch("pthatweight",&pthatweight,"pthatweight/D");
    //ntw->Branch("evt",&evnt_n,"evt/I");
    //ntw->Branch("lumi",&lumi_n,"lumi/I");
    //ntw->Branch("vz",&vz_n,"vz/F");
    
    for(int ie = 0; ie < newTree->GetEntries(); ++ie){

      if(ie%50000 == 0) cout<<ie<<"/"<<newTree->GetEntries()<<endl;
    
      newTree->GetEntry(ie);
    
      for(int i = 0; i < Npt; ++i){
	if(n[i] > 0 && varpthat >= pthats[i]) pthatweight = (float)(xs[i]/n[i]);
      }

      // vz_n = vz;
      // lumi_n = lumi;
      // evnt_n = evnt;
      
      ntw->Fill();

    }

    outf[i]->cd();
    ntw->Write();
    outf[i]->Write();
    outf[i]->Close();

  }

}

