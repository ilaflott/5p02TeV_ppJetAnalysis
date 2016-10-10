// Raghav Kunnawalkam Elayavalli
// Nov 24th 2015
// Rutgers @CERN
// for questions or comments: raghav.k.e at CERN dot CH

// read all the MC files for pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 

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

// static const int nbins_pt = 39;
// static const double boundaries_pt[nbins_pt+1] = {
//   3, 4, 5, 7, 9, 12, 
//   15, 18, 21, 24, 28,
//   32, 37, 43, 49, 56,
//   64, 74, 84, 97, 114,
//   133, 153, 174, 196,
//   220, 245, 272, 300, 
//   330, 362, 395, 430,
//   468, 507, 548, 592,
//   638, 686, 1000 
// };

const int ptbins[] = {15, 30, 50, 80, 120, 170, 220, 300, 500};
const int nbins_pt = sizeof(ptbins)/sizeof(int) -1;

const double ptbins_jec[] = {17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 350, 400, 550, 790, 1000};
const int nbins_pt_jec = sizeof(ptbins_jec)/sizeof(double) -1;

const double etabins[] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, 
            -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
            -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
            -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
            -0.435, -0.348, -0.261, -0.174, -0.087, 
            +0.000, 
            +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
            +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
            +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, 
            +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, 
            +4.363, +4.538, +4.716, +4.889, +5.191
};
const int nbins_eta = sizeof(etabins)/sizeof(double) -1;

static const int trigValue = 5;
static const char trigName [trigValue][256] = {"HLT40","HLT60","HLT80","HLT100","Combined"};
static const char * etaWidth = (char*)"20_eta_20";

// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

float deltaphi(float phi1, float phi2)
{
  float pi=TMath::Pi();
 
  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}



using namespace std;

void RAA_read_mc_pp(int startfile = 0,
		    int endfile = 8,
		    int radius = 4,
		    std::string algo="",
		    std::string jetType= "PF",
		    std::string kFoname="test_output.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  double pthatBins2015[] = {15, 30, 50, 80, 120, 170, 220, 280, 2000};

  infile_Forest = "pp_mc_forests_new.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 4; //6

  TChain * jetpp[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str());
  dir[3] = "ppTrack";
  //dir[3] = "hiEvtAnalyzer";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "trackTree"
    //"HiTree"
    // , "jetObjTree"
  };

  for(int t = 0;t<N;t++){
    jetpp[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  TFile * fweight = new TFile(Form("weights_pp_%d.root",startfile), "r");
  TTree * weight_Tree = (TTree*)fweight->Get("weights");

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    jetpp[0]->Add(filename_Forest.c_str());
    jetpp[1]->Add(filename_Forest.c_str());
    jetpp[2]->Add(filename_Forest.c_str());
    jetpp[3]->Add(filename_Forest.c_str());
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[3]->GetEntries() << endl;

  }
  
  jetpp[2]->AddFriend(jetpp[0]);
  jetpp[2]->AddFriend(jetpp[1]);
  jetpp[2]->AddFriend(jetpp[3]);
  
  // Forest files 
  int nref_F;
  float pt_F[1000];
  float refpt_F[1000];
  float refeta_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float chMax_F[1000];
  int chN_F[1000];
  float trkMax_F[1000];
  float chSum_F[1000];
  float phSum_F[1000];
  float neSum_F[1000];
  float trkSum_F[1000];
  float phMax_F[1000];
  float neMax_F[1000];

  float chHardMax_F[1000];
  float trkHardMax_F[1000];
  float chHardSum_F[1000];
  float phHardSum_F[1000];
  float trkHardSum_F[1000];
  float phHardMax_F[1000];

  float eMax_F[1000];
  float muMax_F[1000];
  float eSum_F[1000];
  float muSum_F[1000];
  float jtpu_F[1000];
  int subid_F[1000];
  float refdrjt_F[1000];
  int refparton_F[1000];

  float pthat_F;
  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet100_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  int jet100_p_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int pcollisionEventSelection_F;

  int nVtx;
  float vz_F[100];
  float vx_F[100];
  float vy_F[100];
  int nTrk;
  
  jetpp[3]->SetBranchAddress("nTrk",&nTrk);
  jetpp[3]->SetBranchAddress("nVtx",&nVtx);
  jetpp[3]->SetBranchAddress("zVtx",vz_F);
  jetpp[3]->SetBranchAddress("yVtx",vy_F);
  jetpp[3]->SetBranchAddress("xVtx",vx_F);
  
  // jetpp[4]->SetBranchAddress("evt",&evt_F);
  // jetpp[4]->SetBranchAddress("run",&run_F);
  // jetpp[4]->SetBranchAddress("lumi",&lumi_F);
  // jetpp[4]->SetBranchAddress("vz",&vz_F);
  //jetpp[1]->SetBranchAddress("pPAcollisionEventSelectionPA",&pcollisionEventSelection_F);
  //jetpp[0]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  //jetpp[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  //jetpp[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpp[2]->SetBranchAddress("nref",&nref_F);
  jetpp[2]->SetBranchAddress("subid",&subid_F);
  jetpp[2]->SetBranchAddress("refdrjt",refdrjt_F);
  jetpp[2]->SetBranchAddress("refparton_flavor",refparton_F);
  jetpp[2]->SetBranchAddress("refpt", &refpt_F);
  jetpp[2]->SetBranchAddress("refeta", &refeta_F);
  jetpp[2]->SetBranchAddress("jtpt",&pt_F);
  jetpp[2]->SetBranchAddress("jteta",&eta_F);
  jetpp[2]->SetBranchAddress("jtphi",&phi_F);
  jetpp[2]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[2]->SetBranchAddress("jtpu",&jtpu_F);
  jetpp[2]->SetBranchAddress("pthat", & pthat_F);
  jetpp[2]->SetBranchAddress("chargedMax",&chMax_F);
  jetpp[2]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[2]->SetBranchAddress("chargedN",&chN_F);
  //jetpp[2]->SetBranchAddress("chargedHardMax",&chMax_F);
  jetpp[2]->SetBranchAddress("chargedHardSum",&chSum_F);
  jetpp[2]->SetBranchAddress("trackMax",&trkMax_F);
  jetpp[2]->SetBranchAddress("trackSum",&trkSum_F);
  //jetpp[2]->SetBranchAddress("trackHardMax",&trkMax_F);
  jetpp[2]->SetBranchAddress("trackHardSum",&trkSum_F);
  jetpp[2]->SetBranchAddress("photonMax",&phMax_F);
  jetpp[2]->SetBranchAddress("photonSum",&phSum_F);
  //jetpp[2]->SetBranchAddress("photonHardMax",&phMax_F);
  jetpp[2]->SetBranchAddress("photonHardSum",&phSum_F);
  jetpp[2]->SetBranchAddress("neutralMax",&neMax_F);
  jetpp[2]->SetBranchAddress("neutralSum",&neSum_F);

  jetpp[2]->SetBranchAddress("eSum",&eSum_F);
  jetpp[2]->SetBranchAddress("eMax",&eMax_F);
  jetpp[2]->SetBranchAddress("muSum",&muSum_F);
  jetpp[2]->SetBranchAddress("muMax",&muMax_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&jet40_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1_Prescl",&jet40_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&jet60_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1_Prescl",&jet60_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&jet80_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1_Prescl",&jet80_p_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1",&jet100_F);
  jetpp[0]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1_Prescl",&jet100_p_F);

  double pthatweight;
  // // int evt_w;
  // // int lumi_w;
  //float vz_w;

  weight_Tree->SetBranchAddress("pthatweight",&pthatweight);
  // // weight_Tree->SetBranchAddress("evt",&evt_w);
  // // weight_Tree->SetBranchAddress("lumi",&lumi_w);
  //weight_Tree->SetBranchAddress("vz",&vz_w);

  jetpp[2]->AddFriend(weight_Tree);
  
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();
  
  // TH1F *hpp_Jet100_gen,*hpp_Jet100_reco;
  // TH1F *hpp_Jet80_gen,*hpp_Jet80_reco;
  // TH1F *hpp_Jet60_gen,*hpp_Jet60_reco;
  // TH1F *hpp_Jet40_gen,*hpp_Jet40_reco;
  // TH1F *hpp_JetComb_gen,*hpp_JetComb_reco;

  // TH1F *hpp_anaBin_Jet100_gen,*hpp_anaBin_Jet100_reco;
  // TH1F *hpp_anaBin_Jet80_gen,*hpp_anaBin_Jet80_reco;
  // TH1F *hpp_anaBin_Jet60_gen,*hpp_anaBin_Jet60_reco;
  // TH1F *hpp_anaBin_Jet40_gen,*hpp_anaBin_Jet40_reco;
  // TH1F *hpp_anaBin_JetComb_gen,*hpp_anaBin_JetComb_reco;
  
  TH1F *hpp_gen,*hpp_reco;
  TH2F *hpp_matrix;
  TH2F *hpp_mcclosure_matrix;
  TH1F *hpp_mcclosure_data;
  TH1F *hpp_mcclosure_gen;

  // TH2F *hpp_matrix_HLT;
  // TH2F *hpp_Trans_matrix_HLT;
  // TH2F *hpp_anaBin_matrix_HLT;
  // TH2F *hpp_anaBin_Trans_matrix_HLT;
  // TH2F *hpp_mcclosure_matrix_HLT;
  // TH2F *hpp_mcclosure_Trans_matrix_HLT;
  // TH2F *hpp_anaBin_mcclosure_matrix_HLT;
  // TH2F *hpp_anaBin_mcclosure_Trans_matrix_HLT;

  // TH2F *hpp_matrix_HLT_Jet80;
  // TH2F *hpp_Trans_matrix_HLT_Jet80;
  // TH2F *hpp_anaBin_matrix_HLT_Jet80;
  // TH2F *hpp_anaBin_Trans_matrix_HLT_Jet80;
  // TH2F *hpp_mcclosure_matrix_HLT_Jet80;
  // TH2F *hpp_mcclosure_Trans_matrix_HLT_Jet80;
  // TH2F *hpp_anaBin_mcclosure_matrix_HLT_Jet80;
  // TH2F *hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet80;

  // TH2F *hpp_matrix_HLT_Jet60;
  // TH2F *hpp_Trans_matrix_HLT_Jet60;
  // TH2F *hpp_anaBin_matrix_HLT_Jet60;
  // TH2F *hpp_anaBin_Trans_matrix_HLT_Jet60;
  // TH2F *hpp_mcclosure_matrix_HLT_Jet60;
  // TH2F *hpp_mcclosure_Trans_matrix_HLT_Jet60;
  // TH2F *hpp_anaBin_mcclosure_matrix_HLT_Jet60;
  // TH2F *hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet60;

  // TH2F *hpp_matrix_HLT_Jet40;
  // TH2F *hpp_Trans_matrix_HLT_Jet40;
  // TH2F *hpp_anaBin_matrix_HLT_Jet40;
  // TH2F *hpp_anaBin_Trans_matrix_HLT_Jet40;
  // TH2F *hpp_mcclosure_matrix_HLT_Jet40;
  // TH2F *hpp_mcclosure_Trans_matrix_HLT_Jet40;
  // TH2F *hpp_anaBin_mcclosure_matrix_HLT_Jet40;
  // TH2F *hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet40;

  // TH2F * hpp_Trans_matrix_HLT_GenSmear;
  // //TH2F * hpp_Trans_matrix_HLT_gen2pSmear;
  // TH2F * hpp_Trans_matrix_HLT_RecoSmear;
  // TH2F * hpp_Trans_matrix_HLT_BothSmear;

  // TH2F * hpp_matrix_HLT_GenSmear;
  // //TH2F * hpp_matrix_HLT_gen2pSmear;
  // TH2F * hpp_matrix_HLT_RecoSmear;
  // TH2F * hpp_matrix_HLT_BothSmear;
  
  //TH2F *hpp_response;
  // TH1F *hpp_mcclosure_JetComb_data;
  TH1F *hpp_mcclosure_data_train;
  // TH1F *hpp_mcclosure_JetComb_data_train;
  // TH1F *hpp_mcclosure_Jet80_data_train;
  // TH1F *hpp_mcclosure_Jet60_data_train;
  // TH1F *hpp_mcclosure_Jet40_data_train;
  // TH1F *hpp_mcclosure_Jet80_data;
  // TH1F *hpp_mcclosure_Jet60_data;
  // TH1F *hpp_mcclosure_Jet40_data;
  // TH1F *hpp_mcclosure_JetComb_gen;
  // TH1F *hpp_mcclosure_Jet80_gen;
  // TH1F *hpp_mcclosure_Jet60_gen;
  // TH1F *hpp_mcclosure_Jet40_gen;

  TH1F * hpthat = new TH1F("hpthat","",1000,0,1000);

  // macros for quark vs gluon uncertainty from the MC.
  // TH1F * hpp_quark_Jet80_gen;
  // TH1F * hpp_quark_Jet60_gen;
  // TH1F * hpp_quark_Jet40_gen;
  // TH1F * hpp_quark_JetComb_gen;
  // TH1F * hpp_quark_Jet80_reco;
  // TH1F * hpp_quark_Jet60_reco;
  // TH1F * hpp_quark_Jet40_reco;
  // TH1F * hpp_quark_JetComb_reco;

  // TH1F * hpp_gluon_Jet80_gen;
  // TH1F * hpp_gluon_Jet60_gen;
  // TH1F * hpp_gluon_Jet40_gen;
  // TH1F * hpp_gluon_JetComb_gen;
  // TH1F * hpp_gluon_Jet80_reco;
  // TH1F * hpp_gluon_Jet60_reco;
  // TH1F * hpp_gluon_Jet40_reco;
  // TH1F * hpp_gluon_JetComb_reco;

  // TH1F * hpp_mcclosure_quark_Jet80_gen;
  // TH1F * hpp_mcclosure_quark_Jet60_gen;
  // TH1F * hpp_mcclosure_quark_Jet40_gen;
  // TH1F * hpp_mcclosure_quark_JetComb_gen;
  // TH1F * hpp_mcclosure_quark_Jet80_reco;
  // TH1F * hpp_mcclosure_quark_Jet60_reco;
  // TH1F * hpp_mcclosure_quark_Jet40_reco;
  // TH1F * hpp_mcclosure_quark_JetComb_reco;

  // TH1F * hpp_mcclosure_gluon_Jet80_gen;
  // TH1F * hpp_mcclosure_gluon_Jet60_gen;
  // TH1F * hpp_mcclosure_gluon_Jet40_gen;
  // TH1F * hpp_mcclosure_gluon_JetComb_gen;
  // TH1F * hpp_mcclosure_gluon_Jet80_reco;
  // TH1F * hpp_mcclosure_gluon_Jet60_reco;
  // TH1F * hpp_mcclosure_gluon_Jet40_reco;
  // TH1F * hpp_mcclosure_gluon_JetComb_reco;
  
  // TH2F * hpp_Trans_matrix_quark_HLT;
  // TH2F * hpp_Trans_matrix_gluon_HLT;

  // TH2F * hpp_mcclosure_Trans_matrix_quark_HLT;
  // TH2F * hpp_mcclosure_Trans_matrix_gluon_HLT;  

  // hpp_quark_Jet80_gen = new TH1F(Form("hpp_quark_Jet80_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_quark_Jet60_gen = new TH1F(Form("hpp_quark_Jet60_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_quark_Jet40_gen = new TH1F(Form("hpp_quark_Jet40_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_quark_JetComb_gen = new TH1F(Form("hpp_quark_JetComb_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_quark_Jet80_reco = new TH1F(Form("hpp_quark_Jet80_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_quark_Jet60_reco = new TH1F(Form("hpp_quark_Jet60_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_quark_Jet40_reco = new TH1F(Form("hpp_quark_Jet40_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_quark_JetComb_reco = new TH1F(Form("hpp_quark_JetComb_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_Jet80_gen = new TH1F(Form("hpp_gluon_Jet80_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_Jet60_gen = new TH1F(Form("hpp_gluon_Jet60_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_Jet40_gen = new TH1F(Form("hpp_gluon_Jet40_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_JetComb_gen = new TH1F(Form("hpp_gluon_JetComb_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_Jet80_reco = new TH1F(Form("hpp_gluon_Jet80_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_Jet60_reco = new TH1F(Form("hpp_gluon_Jet60_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_Jet40_reco = new TH1F(Form("hpp_gluon_Jet40_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_gluon_JetComb_reco = new TH1F(Form("hpp_gluon_JetComb_reco_R%d", radius),"",1000, 0, 1000);

  // hpp_mcclosure_quark_Jet80_gen = new TH1F(Form("hpp_mcclosure_quark_Jet80_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_quark_Jet60_gen = new TH1F(Form("hpp_mcclosure_quark_Jet60_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_quark_Jet40_gen = new TH1F(Form("hpp_mcclosure_quark_Jet40_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_quark_JetComb_gen = new TH1F(Form("hpp_mcclosure_quark_JetComb_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_quark_Jet80_reco = new TH1F(Form("hpp_mcclosure_quark_Jet80_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_quark_Jet60_reco = new TH1F(Form("hpp_mcclosure_quark_Jet60_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_quark_Jet40_reco = new TH1F(Form("hpp_mcclosure_quark_Jet40_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_quark_JetComb_reco = new TH1F(Form("hpp_mcclosure_quark_JetComb_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_Jet80_gen = new TH1F(Form("hpp_mcclosure_gluon_Jet80_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_Jet60_gen = new TH1F(Form("hpp_mcclosure_gluon_Jet60_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_Jet40_gen = new TH1F(Form("hpp_mcclosure_gluon_Jet40_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_JetComb_gen = new TH1F(Form("hpp_mcclosure_gluon_JetComb_gen_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_Jet80_reco = new TH1F(Form("hpp_mcclosure_gluon_Jet80_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_Jet60_reco = new TH1F(Form("hpp_mcclosure_gluon_Jet60_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_Jet40_reco = new TH1F(Form("hpp_mcclosure_gluon_Jet40_reco_R%d", radius),"",1000, 0, 1000);
  // hpp_mcclosure_gluon_JetComb_reco = new TH1F(Form("hpp_mcclosure_gluon_JetComb_reco_R%d", radius),"",1000, 0, 1000);

  // hpp_Trans_matrix_quark_HLT = new TH2F(Form("hpp_Trans_matrix_quark_HLT_R%d",radius),"",1000, 0, 1000, 1000, 0, 1000);
  // hpp_Trans_matrix_gluon_HLT = new TH2F(Form("hpp_Trans_matrix_gluon_HLT_R%d",radius),"",1000, 0, 1000, 1000, 0, 1000);
  // hpp_mcclosure_Trans_matrix_quark_HLT = new TH2F(Form("hpp_mcclosure_Trans_matrix_quark_HLT_R%d",radius),"",1000, 0, 1000, 1000, 0, 1000);
  // hpp_mcclosure_Trans_matrix_gluon_HLT = new TH2F(Form("hpp_mcclosure_Trans_matrix_gluon_HLT_R%d",radius),"",1000, 0, 1000, 1000, 0, 1000);
  

  // hpp_Trans_matrix_HLT = new TH2F(Form("hpp_Trans_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000, 0, 1000, 1000, 0, 1000);
  // hpp_Trans_matrix_HLT_GenSmear = new TH2F(Form("hpp_Trans_matrix_HLT_GenSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);
  // hpp_Trans_matrix_HLT_RecoSmear = new TH2F(Form("hpp_Trans_matrix_HLT_RecoSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);
  // hpp_Trans_matrix_HLT_BothSmear = new TH2F(Form("hpp_Trans_matrix_HLT_BothSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);
  // //hpp_Trans_matrix_HLT_gen2pSmear = new TH2F(Form("hpp_Trans_matrix_HLT_gen2pSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);

  // hpp_matrix_HLT_GenSmear = new TH2F(Form("hpp_matrix_HLT_GenSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);
  // hpp_matrix_HLT_RecoSmear = new TH2F(Form("hpp_matrix_HLT_RecoSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);
  // hpp_matrix_HLT_BothSmear = new TH2F(Form("hpp_matrix_HLT_BothSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);
  // //hpp_matrix_HLT_gen2pSmear = new TH2F(Form("hpp_matrix_HLT_gen2pSmear_R%d_%s", radius, etaWidth),"",1000, 0, 1000, 1000, 0, 1000);
  
  
  hpp_gen = new TH1F(Form("hpp_gen_R%d_%s",radius,etaWidth),Form("Gen refpt R%d %s ",radius,etaWidth),1000,0,1000);
  //cout<<"A"<<endl;
  hpp_reco = new TH1F(Form("hpp_reco_R%d_%s",radius,etaWidth),Form("Reco jtpt R%d %s ",radius,etaWidth),1000,0,1000);
  //cout<<"B"<<endl;
  hpp_matrix = new TH2F(Form("hpp_matrix_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt R%d %s ; reco pT; gen pT",radius,etaWidth),1000,0,1000,1000,0,1000);

  // hpp_matrix_HLT = new TH2F(Form("hpp_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000,0,1000,1000,0,1000);
  // hpp_anaBin_matrix_HLT = new TH2F(Form("hpp_anaBin_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_Trans_matrix_HLT = new TH2F(Form("hpp_anaBin_Trans_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

  // hpp_matrix_HLT_Jet80 = new TH2F(Form("hpp_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000,0,1000,1000,0,1000);
  // hpp_Trans_matrix_HLT_Jet80 = new TH2F(Form("hpp_Trans_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000,0,1000,1000,0,1000);
  // hpp_anaBin_matrix_HLT_Jet80 = new TH2F(Form("hpp_anaBin_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_Trans_matrix_HLT_Jet80 = new TH2F(Form("hpp_anaBin_Trans_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

  // hpp_matrix_HLT_Jet60 = new TH2F(Form("hpp_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000,0,1000,1000,0,1000);
  // hpp_Trans_matrix_HLT_Jet60 = new TH2F(Form("hpp_Trans_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000,0,1000,1000,0,1000);
  // hpp_anaBin_matrix_HLT_Jet60 = new TH2F(Form("hpp_anaBin_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_Trans_matrix_HLT_Jet60 = new TH2F(Form("hpp_anaBin_Trans_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

  // hpp_matrix_HLT_Jet40 = new TH2F(Form("hpp_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000,0,1000,1000,0,1000);
  // hpp_Trans_matrix_HLT_Jet40 = new TH2F(Form("hpp_Trans_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),1000,0,1000,1000,0,1000);
  // hpp_anaBin_matrix_HLT_Jet40 = new TH2F(Form("hpp_anaBin_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_Trans_matrix_HLT_Jet40 = new TH2F(Form("hpp_anaBin_Trans_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);


  // hpp_JetComb_gen = new TH1F(Form("hpp_JetComb_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from HLT trigger combined R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_JetComb_reco = new TH1F(Form("hpp_JetComb_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from HLT trigger combined R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_Jet80_gen = new TH1F(Form("hpp_Jet80_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet80 trigger R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_Jet80_reco = new TH1F(Form("hpp_Jet80_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet80 trigger R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_Jet60_gen = new TH1F(Form("hpp_Jet60_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_Jet60_reco = new TH1F(Form("hpp_Jet60_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_Jet40_gen = new TH1F(Form("hpp_Jet40_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_Jet40_reco = new TH1F(Form("hpp_Jet40_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),1000,0,1000);
  // hpp_anaBin_JetComb_gen = new TH1F(Form("hpp_anaBin_JetComb_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from HLT trigger combined R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_JetComb_reco = new TH1F(Form("hpp_anaBin_JetComb_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from HLT trigger combined R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_Jet80_gen = new TH1F(Form("hpp_anaBin_Jet80_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_Jet80_reco = new TH1F(Form("hpp_anaBin_Jet80_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_Jet60_gen = new TH1F(Form("hpp_anaBin_Jet60_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_Jet60_reco = new TH1F(Form("hpp_anaBin_Jet60_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_Jet40_gen = new TH1F(Form("hpp_anaBin_Jet40_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  // hpp_anaBin_Jet40_reco = new TH1F(Form("hpp_anaBin_Jet40_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);  

 
  // hpp_anaBin_mcclosure_Trans_matrix_HLT = new TH2F(Form("hpp_anaBin_mcclosure_Trans_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_mcclosure_matrix_HLT = new TH2F(Form("hpp_anaBin_mcclosure_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_mcclosure_matrix_HLT = new TH2F(Form("hpp_mcclosure_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);
  // hpp_mcclosure_Trans_matrix_HLT = new TH2F(Form("hpp_mcclosure_Trans_matrix_HLT_R%d_%s",radius,etaWidth),Form("Trans Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);

  // hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet40 = new TH2F(Form("hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_mcclosure_matrix_HLT_Jet40 = new TH2F(Form("hpp_anaBin_mcclosure_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_mcclosure_matrix_HLT_Jet40 = new TH2F(Form("hpp_mcclosure_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);
  // hpp_mcclosure_Trans_matrix_HLT_Jet40 = new TH2F(Form("hpp_mcclosure_Trans_matrix_HLT_Jet40_R%d_%s",radius,etaWidth),Form("Trans Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);

  // hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet60 = new TH2F(Form("hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_mcclosure_matrix_HLT_Jet60 = new TH2F(Form("hpp_anaBin_mcclosure_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_mcclosure_matrix_HLT_Jet60 = new TH2F(Form("hpp_mcclosure_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);
  // hpp_mcclosure_Trans_matrix_HLT_Jet60 = new TH2F(Form("hpp_mcclosure_Trans_matrix_HLT_Jet60_R%d_%s",radius,etaWidth),Form("Trans Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);

  // hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet80 = new TH2F(Form("hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_anaBin_mcclosure_matrix_HLT_Jet80 = new TH2F(Form("hpp_anaBin_mcclosure_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  // hpp_mcclosure_matrix_HLT_Jet80 = new TH2F(Form("hpp_mcclosure_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);
  // hpp_mcclosure_Trans_matrix_HLT_Jet80 = new TH2F(Form("hpp_mcclosure_Trans_matrix_HLT_Jet80_R%d_%s",radius,etaWidth),Form("Trans Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);


  hpp_mcclosure_matrix = new TH2F(Form("hpp_mcclosure_matrix_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt R%d %s;reco pT; gen pT ",radius,etaWidth),1000, 0, 1000,1000, 0, 1000);
  hpp_mcclosure_data = new TH1F(Form("hpp_mcclosure_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test R%d %s ",radius,etaWidth),1000, 0, 1000);
  hpp_mcclosure_gen = new TH1F(Form("hpp_mcclosure_gen_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test R%d %s ",radius,etaWidth),1000, 0, 1000);

  // hpp_mcclosure_JetComb_data = new TH1F(Form("hpp_mcclosure_JetComb_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger combined  R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet80_data = new TH1F(Form("hpp_mcclosure_Jet80_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 80  R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet60_data = new TH1F(Form("hpp_mcclosure_Jet60_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 60  R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet40_data = new TH1F(Form("hpp_mcclosure_Jet40_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 40  R%d %s ",radius,etaWidth),1000, 0, 1000);
  hpp_mcclosure_data_train = new TH1F(Form("hpp_mcclosure_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_JetComb_data_train = new TH1F(Form("hpp_mcclosure_JetComb_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger combined  R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet80_data_train = new TH1F(Form("hpp_mcclosure_Jet80_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger 80  R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet60_data_train = new TH1F(Form("hpp_mcclosure_Jet60_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger 60  R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet40_data_train = new TH1F(Form("hpp_mcclosure_Jet40_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger 40  R%d %s ",radius,etaWidth),1000, 0, 1000);  
  // hpp_mcclosure_JetComb_gen = new TH1F(Form("hpp_mcclosure_gen_JetComb_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger combined R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet80_gen = new TH1F(Form("hpp_mcclosure_gen_Jet80_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet60_gen = new TH1F(Form("hpp_mcclosure_gen_Jet60_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 60 R%d %s ",radius,etaWidth),1000, 0, 1000);
  // hpp_mcclosure_Jet40_gen = new TH1F(Form("hpp_mcclosure_gen_Jet40_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 40 R%d %s ",radius,etaWidth),1000, 0, 1000);


  TH3F * hJEC= new TH3F("hJEC",";raw p_{T};#eta;JEC",500, 0, 500, 200, -5, +5, 300, 0, 5);
  TH1F * hJER[nbins_pt];
  TH1F * hJER_eta_30pt50[nbins_eta];
  TH1F * hJER_eta_150pt200[nbins_eta];
  TH1F * hJEC_check[nbins_pt_jec][nbins_eta];

  for(int x = 0; x<nbins_pt_jec; ++x){
    for(int y = 0; y<nbins_eta; ++y){
      hJEC_check[x][y] = new TH1F(Form("hJEC_check_ptbin%d_etabin%d",x,y),Form("rawpt/genpt %2.0f < genpt < %2.0f, %2.4f < geneta < %2.4f",ptbins_jec[x], ptbins_jec[x+1], etabins[y], etabins[y+1]),100, 0, 3);
    }
  }
  
  TH1F * hVz = new TH1F("hVz","",200, -20, 20);
  
  for(int bin = 0; bin<nbins_pt; ++bin){
    hJER[bin] = new TH1F(Form("hJER_%d_pt_%d", ptbins[bin], ptbins[bin+1]),"",100, 0, 2);
  }

  for(int bin = 0; bin<nbins_eta; ++bin){
    hJER_eta_30pt50[bin] = new TH1F(Form("hJER_etabin%d_30_pt_50", bin),Form("rawpt/genpt 30 < genpt < 50, %2.4f < geneta < %2.4f", etabins[bin], etabins[bin+1]),100, 0, 2);
    hJER_eta_150pt200[bin] = new TH1F(Form("hJER_etabin%d_150_pt_200", bin),Form("rawpt/genpt 150 < genpt < 200, %2.4f < geneta < %2.4f", etabins[bin], etabins[bin+1]),100, 0, 2);
  }

  // Add the Jet ID plots:
  // list of variables:
  std::string var[21] = {"jtpt" ,"rawpt", "jteta", "jtphi", "trkMax", "trkSum", "trkHardSum", "chMax", "chSum", "chHardSum","phMax", "phSum", "phHardSum", "neMax", "neSum", "eMax", "eSum", "muMax", "muSum" ,"Aj","xj"};
  TH1F * hJetQA[2][21];

  for(int k = 0; k<2; ++k){
    for(int j = 0; j<21; ++j){
      if(j==2) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),100, -5, +5);
      else if(j==3) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),100, -4, +4);
      else if(j<=1)hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),500, 0, 500);
      else if(j>=4)hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),Form(";%s;",var[j].c_str()),200, 0, 2);
    }
  }
  // double RMS = 0.22;
  
  // TH1F * hpp_DijetSignalSelection_LeadPFJet_gt_A = new TH1F("hpp_DijetSignalSelection_LeadPFJet_gt_A","",1000, 0, 1000);
  // TH2F * hdphiptcent_DijetSignalSelection_LeadPFJet_gt_A = new TH2F("hdphiptcent_DijetSignalSelection_LeadPFJet_gt_A","",100, 0, 200,160,0,3.2);  
  // TH1F * hpp_DijetBKGSelection_LeadPFJet_lt_A = new TH1F("hpp_DijetBKGSelection_LeadPFJet_lt_A","",1000, 0, 1000);;
  // TH2F * hdphiptcent_DijetBKGSelection_LeadPFJet_lt_A = new TH2F("hdphiptcent_DijetBKGSelection_LeadPFJet_lt_A","",100, 0, 200, 160,0,3.2); 
  
  // now start the event loop for each file.
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpp[0]->GetEntries();
  //Long64_t nGoodEvt = 0;
  if(printDebug) nentries = 10;
  TRandom rnd; 
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  // my weights 
  //fVzPP->SetParameters(1.05602e00,5.74688e-03,-3.37288e-03,-1.44764e-05,8.59060e-07);
  // anna's weights
  fVzPP->SetParameters(+0.941, -0.0173, +3.23e-3, +3.61e-6, -1.04e-5);
  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    jetpp[0]->GetEntry(nEvt);
    jetpp[1]->GetEntry(nEvt);
    jetpp[2]->GetEntry(nEvt);
    // jetpp[4]->GetEntry(nEvt);
    jetpp[3]->GetEntry(nEvt);
    weight_Tree->GetEntry(nEvt);
    
    //if(printDebug) cout<<"forest values = "<<hiBin_F<<", "<<evt_F<<", "<<run_F<<", "<<lumi_F<<", "<<vz_F<<endl;
    
    double weight_eS = pthatweight;
    if(printDebug) cout<<"weight_eS ="<<weight_eS<<endl;
    
    if(fabs(vz_F[0])>15) continue;

    //if(pcollisionEventSelection_F==0) continue; 
    // aply proxy for collision eventselection : position.Rho <= 2 && tracksSize >= 2
    double rho = TMath::Sqrt((vx_F[0]*vx_F[0]) + (vy_F[0]*vy_F[0]));
    if(rho > 2 || nTrk < 2) continue;
    
    double vz_weight = fVzPP->Eval(vz_F[0]);
    weight_eS = pthatweight * vz_weight;

    hpthat->Fill(pthat_F, weight_eS);

    hVz->Fill(vz_F[0], weight_eS);
    //cout<<"after filling vz"<<endl;
    
    // std::vector <float> pt;
    // std::vector <float> phi;
    // std::vector <float> eta;
    
    for(int jet = 0; jet<nref_F; ++jet){

      if(fabs(eta_F[jet]) > 2) continue;
      if(subid_F[jet] != 0) continue;
      if(pt_F[jet] > 3 * pthat_F) continue;
      if(refdrjt_F[jet] > 0.3) continue; 

      float genpt = refpt_F[jet];
      float geneta = refeta_F[jet];
      float recpt = pt_F[jet];
      float rawpt = rawpt_F[jet];
	
      hJEC->Fill(rawpt, eta_F[jet], (float)(recpt/rawpt));

      if(fabs(geneta)>=5.191) continue;
      int etabin = -1;
      for(int bin = 0; bin<nbins_eta; ++bin){
	if(geneta > etabins[bin]) etabin = bin;
      }
      if(etabin == -1) continue; 
      if(genpt >= 30 && genpt<50) hJER_eta_30pt50[etabin]->Fill((float)recpt/genpt, weight_eS);
      if(genpt >= 150 && genpt<200) hJER_eta_150pt200[etabin]->Fill((float)recpt/genpt, weight_eS);

      int binx = -1;
      for(int bin = 0; bin<nbins_pt_jec; ++bin){
	if(genpt > ptbins_jec[bin]) binx = bin;
      }
      if(binx == -1) continue;
      hJEC_check[binx][etabin]->Fill((float)rawpt/genpt);
      
      int ptbin = -1;
      for(int bin = 0; bin<nbins_pt; ++bin){
	if(genpt > ptbins[bin]) ptbin = bin;
      }
      if(ptbin == -1) continue;
      hJER[ptbin]->Fill((float)recpt/genpt, weight_eS);

      // for unfolding later
      hpp_gen->Fill(genpt, weight_eS);
      hpp_reco->Fill(recpt, weight_eS);
      hpp_matrix->Fill(recpt, genpt, weight_eS);

      if(nEvt%2 == 0){
	hpp_mcclosure_data->Fill(recpt, weight_eS);
      }else {
	hpp_mcclosure_matrix->Fill(recpt, genpt, weight_eS);
	hpp_mcclosure_gen->Fill(genpt, weight_eS);
	hpp_mcclosure_data_train->Fill(recpt, weight_eS);
      }

      if(recpt>30){
	hJetQA[0][0]->Fill(recpt, weight_eS);
	hJetQA[0][1]->Fill(rawpt_F[jet], weight_eS);
	hJetQA[0][2]->Fill(eta_F[jet], weight_eS);
	hJetQA[0][3]->Fill(phi_F[jet], weight_eS);
	hJetQA[0][4]->Fill(trkMax_F[jet]/recpt, weight_eS);
	hJetQA[0][5]->Fill(trkSum_F[jet]/recpt, weight_eS);
	hJetQA[0][6]->Fill(trkHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][7]->Fill(chMax_F[jet]/recpt, weight_eS);
	hJetQA[0][8]->Fill(chSum_F[jet]/recpt, weight_eS);
	hJetQA[0][9]->Fill(chHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][10]->Fill(phMax_F[jet]/recpt, weight_eS);
	hJetQA[0][11]->Fill(phSum_F[jet]/recpt, weight_eS);
	hJetQA[0][12]->Fill(phHardSum_F[jet]/recpt, weight_eS);
	hJetQA[0][13]->Fill(neMax_F[jet]/recpt, weight_eS);
	hJetQA[0][14]->Fill(neSum_F[jet]/recpt, weight_eS);
	hJetQA[0][15]->Fill(eMax_F[jet]/recpt, weight_eS);
	hJetQA[0][16]->Fill(eSum_F[jet]/recpt, weight_eS);
	hJetQA[0][17]->Fill(muMax_F[jet]/recpt, weight_eS);
	hJetQA[0][18]->Fill(muSum_F[jet]/recpt, weight_eS);

	// apply JetID
	// charged hadron fraction > 0
	// charged hadron multiplicity > 0
	// charged EM fraction < 0.99 
	
	if(chSum_F[jet]/recpt>0 && neSum_F[jet]/recpt<0.99 && chN_F[jet]>0 && phSum_F[jet]/recpt<0.99){
	  hJetQA[1][0]->Fill(recpt, weight_eS);
	  hJetQA[1][1]->Fill(rawpt_F[jet], weight_eS);
	  hJetQA[1][2]->Fill(eta_F[jet], weight_eS);
	  hJetQA[1][3]->Fill(phi_F[jet], weight_eS);
	  hJetQA[1][4]->Fill(trkMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][5]->Fill(trkSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][6]->Fill(trkHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][7]->Fill(chMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][8]->Fill(chSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][9]->Fill(chHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][10]->Fill(phMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][11]->Fill(phSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][12]->Fill(phHardSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][13]->Fill(neMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][14]->Fill(neSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][15]->Fill(eMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][16]->Fill(eSum_F[jet]/recpt, weight_eS);
	  hJetQA[1][17]->Fill(muMax_F[jet]/recpt, weight_eS);
	  hJetQA[1][18]->Fill(muSum_F[jet]/recpt, weight_eS);
	}
      }
      /*
      if(jet40_F == 1 && jet60_F==0 && jet80_F == 0){
 
	hpp_Jet40_gen->Fill(genpt, weight_eS * jet40_p_F);
	hpp_Jet40_reco->Fill(recpt, weight_eS * jet40_p_F);

	hpp_matrix_HLT_Jet40->Fill(genpt, recpt, weight_eS * jet40_p_F);
	hpp_Trans_matrix_HLT_Jet40->Fill(recpt, genpt, weight_eS * jet40_p_F);
	hpp_anaBin_Trans_matrix_HLT_Jet40->Fill(recpt, genpt, weight_eS * jet40_p_F);
	hpp_anaBin_matrix_HLT_Jet40->Fill(genpt, recpt, weight_eS * jet40_p_F);

	//hpp_Trans_matrix_HLT_gen2pSmear->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS * jet55_p_F);
	hpp_Trans_matrix_HLT_GenSmear->Fill( recpt , genpt + rnd.Gaus(0,1), weight_eS * jet40_p_F);
	hpp_Trans_matrix_HLT_RecoSmear->Fill( recpt + rnd.Gaus(0,1), genpt, weight_eS * jet40_p_F);
	hpp_Trans_matrix_HLT_BothSmear->Fill( recpt + rnd.Gaus(0,1), genpt + rnd.Gaus(0,1), weight_eS * jet40_p_F);

	//hpp_matrix_HLT_gen2pSmear->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS * jet55_p_F);
	hpp_matrix_HLT_GenSmear->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS * jet40_p_F);
	hpp_matrix_HLT_RecoSmear->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS * jet40_p_F);
	hpp_matrix_HLT_BothSmear->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS * jet40_p_F);

	
	hpp_anaBin_Jet40_gen->Fill(genpt, weight_eS * jet40_p_F);
	hpp_anaBin_Jet40_reco->Fill(recpt, weight_eS * jet40_p_F);
	
	if(nEvt%2==0) {
	  hpp_anaBin_mcclosure_matrix_HLT_Jet40->Fill(genpt, recpt, weight_eS * jet40_p_F);
	  hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet40->Fill(recpt, genpt, weight_eS * jet40_p_F);
	  hpp_mcclosure_Trans_matrix_HLT_Jet40->Fill(recpt, genpt, weight_eS * jet40_p_F);
	  hpp_mcclosure_matrix_HLT_Jet40->Fill(genpt, recpt, weight_eS * jet40_p_F);
	  
	  hpp_mcclosure_Jet40_gen->Fill(genpt, weight_eS * jet40_p_F);
	  hpp_mcclosure_Jet40_data_train->Fill(recpt, weight_eS * jet40_p_F);
	}
	if(nEvt%2==1) {
	  hpp_mcclosure_Jet40_data->Fill(recpt, weight_eS * jet40_p_F);
	}


	if(refparton_F[jet] < 9){
	  hpp_quark_Jet40_gen->Fill(genpt, weight_eS * jet40_p_F);
	  hpp_quark_Jet40_reco->Fill(recpt, weight_eS * jet40_p_F);
	  hpp_Trans_matrix_quark_HLT->Fill(recpt, genpt, weight_eS * jet40_p_F);
	  if(nEvt%2==0){
	    hpp_mcclosure_Trans_matrix_quark_HLT->Fill(recpt, genpt, weight_eS * jet40_p_F);
	  }else if(nEvt%2==1){
	    hpp_mcclosure_quark_Jet40_gen->Fill(genpt, weight_eS * jet40_p_F);
	    hpp_mcclosure_quark_Jet40_reco->Fill(recpt, weight_eS * jet40_p_F);
	  }
	}

	if(refparton_F[jet] == 21){
	  hpp_gluon_Jet40_gen->Fill(genpt, weight_eS * jet40_p_F);
	  hpp_gluon_Jet40_reco->Fill(recpt, weight_eS * jet40_p_F);
	  hpp_Trans_matrix_gluon_HLT->Fill(recpt, genpt, weight_eS * jet40_p_F);
	  if(nEvt%2==0){
	    hpp_mcclosure_Trans_matrix_gluon_HLT->Fill(recpt, genpt, weight_eS * jet40_p_F);
	  }else if(nEvt%2==1){
	    hpp_mcclosure_gluon_Jet40_gen->Fill(genpt, weight_eS * jet40_p_F);
	    hpp_mcclosure_gluon_Jet40_reco->Fill(recpt, weight_eS * jet40_p_F);
	  }
	}
	
      }
      
      if(jet60_F == 1 && jet80_F == 0){

	hpp_Jet60_gen->Fill(genpt, weight_eS);
	hpp_Jet60_reco->Fill(recpt, weight_eS);
	hpp_anaBin_Jet60_gen->Fill(genpt, weight_eS);
	hpp_anaBin_Jet60_reco->Fill(recpt, weight_eS);

	//hpp_Trans_matrix_HLT_gen2pSmear->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpp_Trans_matrix_HLT_GenSmear->Fill( recpt , genpt + rnd.Gaus(0,1), weight_eS);
	hpp_Trans_matrix_HLT_RecoSmear->Fill( recpt + rnd.Gaus(0,1), genpt, weight_eS);
	hpp_Trans_matrix_HLT_BothSmear->Fill( recpt + rnd.Gaus(0,1), genpt + rnd.Gaus(0,1), weight_eS);

	//hpp_matrix_HLT_gen2pSmear->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS);
	hpp_matrix_HLT_GenSmear->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS);
	hpp_matrix_HLT_RecoSmear->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS);
	hpp_matrix_HLT_BothSmear->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS);

	
	hpp_matrix_HLT_Jet60->Fill(genpt, recpt, weight_eS);
	hpp_anaBin_matrix_HLT_Jet60->Fill(genpt, recpt, weight_eS);
	hpp_Trans_matrix_HLT_Jet60->Fill(recpt, genpt, weight_eS);
	hpp_anaBin_Trans_matrix_HLT_Jet60->Fill(recpt, genpt, weight_eS);

	if(nEvt%2==0) {
	  hpp_anaBin_mcclosure_matrix_HLT_Jet60->Fill(genpt, recpt, weight_eS);
	  hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet60->Fill(recpt, genpt, weight_eS);
	  hpp_mcclosure_Trans_matrix_HLT_Jet60->Fill(recpt, genpt, weight_eS);
	  hpp_mcclosure_matrix_HLT_Jet60->Fill(genpt, recpt, weight_eS);
	  hpp_mcclosure_Jet60_gen->Fill(genpt, weight_eS);
	  hpp_mcclosure_Jet60_data_train->Fill(recpt, weight_eS);
	}
	if(nEvt%2==1) {
	  hpp_mcclosure_Jet60_data->Fill(recpt, weight_eS);
	}


	if(refparton_F[jet] < 9){
	  hpp_quark_Jet60_gen->Fill(genpt, weight_eS);
	  hpp_quark_Jet60_reco->Fill(recpt, weight_eS);
	  hpp_Trans_matrix_quark_HLT->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpp_mcclosure_Trans_matrix_quark_HLT->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpp_mcclosure_quark_Jet60_gen->Fill(genpt, weight_eS);
	    hpp_mcclosure_quark_Jet60_reco->Fill(recpt, weight_eS);
	  }
	}

	if(refparton_F[jet] == 21){
	  hpp_gluon_Jet60_gen->Fill(genpt, weight_eS);
	  hpp_gluon_Jet60_reco->Fill(recpt, weight_eS);
	  hpp_Trans_matrix_gluon_HLT->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpp_mcclosure_Trans_matrix_gluon_HLT->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpp_mcclosure_gluon_Jet60_gen->Fill(genpt, weight_eS);
	    hpp_mcclosure_gluon_Jet60_reco->Fill(recpt, weight_eS);
	  }
	}
	
      }

      if(jet80_F == 1){

	hpp_Jet80_gen->Fill(genpt, weight_eS);
	hpp_Jet80_reco->Fill(recpt, weight_eS);
	hpp_anaBin_Jet80_gen->Fill(genpt, weight_eS);
	hpp_anaBin_Jet80_reco->Fill(recpt, weight_eS);
	hpp_matrix_HLT_Jet80->Fill(genpt, recpt, weight_eS);
	hpp_anaBin_matrix_HLT_Jet80->Fill(genpt, recpt, weight_eS);
	hpp_Trans_matrix_HLT_Jet80->Fill(recpt, genpt, weight_eS);
	hpp_anaBin_Trans_matrix_HLT_Jet80->Fill(recpt, genpt, weight_eS);
	
	//hpp_Trans_matrix_HLT_gen2pSmear->Fill(recpt, genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpp_Trans_matrix_HLT_GenSmear->Fill( recpt , genpt + rnd.Gaus(0,1), weight_eS);
	hpp_Trans_matrix_HLT_RecoSmear->Fill( recpt + rnd.Gaus(0,1), genpt, weight_eS);
	hpp_Trans_matrix_HLT_BothSmear->Fill( recpt + rnd.Gaus(0,1), genpt + rnd.Gaus(0,1), weight_eS);

	//hpp_matrix_HLT_gen2pSmear->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS);
	hpp_matrix_HLT_GenSmear->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS);
	hpp_matrix_HLT_RecoSmear->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS);
	hpp_matrix_HLT_BothSmear->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS);

	if(nEvt%2==0) {
	  hpp_anaBin_mcclosure_matrix_HLT_Jet80->Fill(genpt, recpt, weight_eS);
	  hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet80->Fill(recpt, genpt, weight_eS);
	  hpp_mcclosure_Trans_matrix_HLT_Jet80->Fill(recpt, genpt, weight_eS);
	  hpp_mcclosure_matrix_HLT_Jet80->Fill(genpt, recpt, weight_eS);
	  hpp_mcclosure_Jet80_gen->Fill(genpt, weight_eS);
	  hpp_mcclosure_Jet80_data_train->Fill(recpt, weight_eS);
	}
	if(nEvt%2==1) {
	  hpp_mcclosure_Jet80_data->Fill(recpt, weight_eS);
	}

	if(refparton_F[jet] < 9){
	  hpp_quark_Jet80_gen->Fill(genpt, weight_eS);
	  hpp_quark_Jet80_reco->Fill(recpt, weight_eS);
	  hpp_Trans_matrix_quark_HLT->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpp_mcclosure_Trans_matrix_quark_HLT->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpp_mcclosure_quark_Jet80_gen->Fill(genpt, weight_eS);
	    hpp_mcclosure_quark_Jet80_reco->Fill(recpt, weight_eS);
	  }
	}

	if(refparton_F[jet] == 21){
	  hpp_gluon_Jet80_gen->Fill(genpt, weight_eS);
	  hpp_gluon_Jet80_reco->Fill(recpt, weight_eS);
	  hpp_Trans_matrix_gluon_HLT->Fill(recpt, genpt, weight_eS);
	  if(nEvt%2==0){
	    hpp_mcclosure_Trans_matrix_gluon_HLT->Fill(recpt, genpt, weight_eS);
	  }else if(nEvt%2==1){
	    hpp_mcclosure_gluon_Jet80_gen->Fill(genpt, weight_eS);
	    hpp_mcclosure_gluon_Jet80_reco->Fill(recpt, weight_eS);
	  }
	}
	
      }

      */}
    // jet loop

    hJetQA[0][19]->Fill((float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]), weight_eS);
    hJetQA[0][20]->Fill((float)(pt_F[1]/pt_F[0]), weight_eS);
	    

    //if(pt.size() <=1) continue;

    // double delphi_lead_sublead = deltaphi(phi[0], phi[1]);
    // bool isSignal = true;

    // if((pt[0] > 3 * RMS) && (delphi_lead_sublead < 2 * (float)TMath::Pi()/3))
    //   isSignal = false;
    // if(pt[0] > 3 * RMS && delphi_lead_sublead > 2 * (float)TMath::Pi()/3 && pt[1] < 1.6 * RMS ) isSignal = false;
    // if(pt[0] < 3 * RMS ) isSignal = false;

    // for(unsigned j = 1; j<pt.size(); ++j){
    //   float delphi = deltaphi (phi[0], phi[j]);
    //   if(isSignal) hdphiptcent_DijetSignalSelection_LeadPFJet_gt_A->Fill(pt[0], delphi, weight_eS);
    //   if(!isSignal) hdphiptcent_DijetBKGSelection_LeadPFJet_lt_A->Fill(pt[0], delphi, weight_eS);      
    // }

    // for(unsigned j = 0; j<pt.size(); ++j){
    //   if(isSignal) hpp_DijetSignalSelection_LeadPFJet_gt_A->Fill(pt[j], weight_eS);
    //   if(!isSignal) hpp_DijetBKGSelection_LeadPFJet_lt_A->Fill(pt[j], weight_eS);
    // }

    //pt.clear();
    //eta.clear();
    //phi.clear();
    
    if(printDebug)cout<<endl;

  }// event loop

  // hpp_quark_JetComb_gen->Add(hpp_quark_Jet80_gen);
  // hpp_quark_JetComb_gen->Add(hpp_quark_Jet60_gen);
  // hpp_quark_JetComb_gen->Add(hpp_quark_Jet40_gen);

  // hpp_gluon_JetComb_gen->Add(hpp_gluon_Jet80_gen);
  // hpp_gluon_JetComb_gen->Add(hpp_gluon_Jet60_gen);
  // hpp_gluon_JetComb_gen->Add(hpp_gluon_Jet40_gen);

  // hpp_quark_JetComb_reco->Add(hpp_quark_Jet80_reco);
  // hpp_quark_JetComb_reco->Add(hpp_quark_Jet60_reco);
  // hpp_quark_JetComb_reco->Add(hpp_quark_Jet40_reco);

  // hpp_gluon_JetComb_reco->Add(hpp_gluon_Jet80_reco);
  // hpp_gluon_JetComb_reco->Add(hpp_gluon_Jet60_reco);
  // hpp_gluon_JetComb_reco->Add(hpp_gluon_Jet40_reco);

  // hpp_mcclosure_quark_JetComb_gen->Add(hpp_mcclosure_quark_Jet80_gen);
  // hpp_mcclosure_quark_JetComb_gen->Add(hpp_mcclosure_quark_Jet60_gen);
  // hpp_mcclosure_quark_JetComb_gen->Add(hpp_mcclosure_quark_Jet40_gen);

  // hpp_mcclosure_gluon_JetComb_gen->Add(hpp_mcclosure_gluon_Jet80_gen);
  // hpp_mcclosure_gluon_JetComb_gen->Add(hpp_mcclosure_gluon_Jet60_gen);
  // hpp_mcclosure_gluon_JetComb_gen->Add(hpp_mcclosure_gluon_Jet40_gen);

  // hpp_mcclosure_quark_JetComb_reco->Add(hpp_mcclosure_quark_Jet80_reco);
  // hpp_mcclosure_quark_JetComb_reco->Add(hpp_mcclosure_quark_Jet60_reco);
  // hpp_mcclosure_quark_JetComb_reco->Add(hpp_mcclosure_quark_Jet40_reco);

  // hpp_mcclosure_gluon_JetComb_reco->Add(hpp_mcclosure_gluon_Jet80_reco);
  // hpp_mcclosure_gluon_JetComb_reco->Add(hpp_mcclosure_gluon_Jet60_reco);
  // hpp_mcclosure_gluon_JetComb_reco->Add(hpp_mcclosure_gluon_Jet40_reco);
    
  
  // hpp_matrix_HLT->Add(hpp_matrix_HLT_Jet80);
  // hpp_matrix_HLT->Add(hpp_matrix_HLT_Jet60);
  // hpp_matrix_HLT->Add(hpp_matrix_HLT_Jet40);
  
  // hpp_anaBin_matrix_HLT->Add(hpp_anaBin_matrix_HLT_Jet80);
  // hpp_anaBin_matrix_HLT->Add(hpp_anaBin_matrix_HLT_Jet60);
  // hpp_anaBin_matrix_HLT->Add(hpp_anaBin_matrix_HLT_Jet40);
  
  // hpp_Trans_matrix_HLT->Add(hpp_Trans_matrix_HLT_Jet80);
  // hpp_Trans_matrix_HLT->Add(hpp_Trans_matrix_HLT_Jet60);
  // hpp_Trans_matrix_HLT->Add(hpp_Trans_matrix_HLT_Jet40);
  
  // hpp_anaBin_Trans_matrix_HLT->Add(hpp_anaBin_Trans_matrix_HLT_Jet80);
  // hpp_anaBin_Trans_matrix_HLT->Add(hpp_anaBin_Trans_matrix_HLT_Jet60);
  // hpp_anaBin_Trans_matrix_HLT->Add(hpp_anaBin_Trans_matrix_HLT_Jet40);
  
  // hpp_mcclosure_matrix_HLT->Add(hpp_mcclosure_matrix_HLT_Jet80);
  // hpp_mcclosure_matrix_HLT->Add(hpp_mcclosure_matrix_HLT_Jet60);
  // hpp_mcclosure_matrix_HLT->Add(hpp_mcclosure_matrix_HLT_Jet40);
  
  // hpp_anaBin_mcclosure_matrix_HLT->Add(hpp_anaBin_mcclosure_matrix_HLT_Jet80);
  // hpp_anaBin_mcclosure_matrix_HLT->Add(hpp_anaBin_mcclosure_matrix_HLT_Jet60);
  // hpp_anaBin_mcclosure_matrix_HLT->Add(hpp_anaBin_mcclosure_matrix_HLT_Jet40);
  
  // hpp_mcclosure_Trans_matrix_HLT->Add(hpp_mcclosure_Trans_matrix_HLT_Jet80);
  // hpp_mcclosure_Trans_matrix_HLT->Add(hpp_mcclosure_Trans_matrix_HLT_Jet60);
  // hpp_mcclosure_Trans_matrix_HLT->Add(hpp_mcclosure_Trans_matrix_HLT_Jet40);
  
  // hpp_anaBin_mcclosure_Trans_matrix_HLT->Add(hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet80);
  // hpp_anaBin_mcclosure_Trans_matrix_HLT->Add(hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet60);
  // hpp_anaBin_mcclosure_Trans_matrix_HLT->Add(hpp_anaBin_mcclosure_Trans_matrix_HLT_Jet40);
  
  // hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet80_data);
  // hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet60_data);
  // hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet40_data);

  // divideBinWidth(hpp_mcclosure_JetComb_data);
  
  // hpp_mcclosure_JetComb_data_train->Add(hpp_mcclosure_Jet80_data_train);
  // hpp_mcclosure_JetComb_data_train->Add(hpp_mcclosure_Jet60_data_train);
  // hpp_mcclosure_JetComb_data_train->Add(hpp_mcclosure_Jet40_data_train);

  // divideBinWidth(hpp_mcclosure_JetComb_data_train);

  // hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet80_gen);
  // hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet60_gen);
  // hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet40_gen);

  // divideBinWidth(hpp_mcclosure_JetComb_gen);
  
  // hpp_JetComb_reco->Add(hpp_Jet80_reco);
  // hpp_JetComb_reco->Add(hpp_Jet60_reco);
  // hpp_JetComb_reco->Add(hpp_Jet40_reco);

  // divideBinWidth(hpp_JetComb_reco);
  
  // hpp_JetComb_gen->Add(hpp_Jet80_gen);
  // hpp_JetComb_gen->Add(hpp_Jet60_gen);
  // hpp_JetComb_gen->Add(hpp_Jet40_gen);

  // divideBinWidth(hpp_JetComb_gen);

  // hpp_anaBin_JetComb_reco->Add(hpp_anaBin_Jet80_reco);
  // hpp_anaBin_JetComb_reco->Add(hpp_anaBin_Jet60_reco);
  // hpp_anaBin_JetComb_reco->Add(hpp_anaBin_Jet40_reco);

  // divideBinWidth(hpp_anaBin_JetComb_reco);
  
  // hpp_anaBin_JetComb_gen->Add(hpp_anaBin_Jet80_gen);
  // hpp_anaBin_JetComb_gen->Add(hpp_anaBin_Jet60_gen);
  // hpp_anaBin_JetComb_gen->Add(hpp_anaBin_Jet40_gen);

  // divideBinWidth(hpp_anaBin_JetComb_gen);

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
