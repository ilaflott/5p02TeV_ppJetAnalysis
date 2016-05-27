// original author: Raghav Kunnawalkam Elayavalli
// Nov 24th 2015
// Overhaul, Ian Laflotte
// Apr 29th 2016
// reads and writes jets from pp data forest files
// for producing quality assurance spectra, and for unfolding later
// compile with...
// g++ RAA_read_data_pp.C $(root-config --cflags --libs) -Werror -Wall -O2 -o RAA_read_data_pp.exe

///////////////////////////////////////////
////////// includes/preprocessor //////////
///////////////////////////////////////////

// C++, C, etc.			    
#include <cstdlib>		    
#include <cstdio>		    
#include <iostream>		    
#include <fstream>		    
#include <cmath>		    
// ROOTSYS
#include <TSystem.h>		    
#include <TProfile.h>		    
// I/O
#include <TChain.h>		    
#include <TFile.h>		    
#include <TTree.h>		    
#include <TNtuple.h>		    
#include <TEventList.h>		    
#include <TCut.h>		    
// utilities			    
#include <TMath.h>		    
#include <TRandom3.h>		    
#include <TStopwatch.h>		    
// plotting			    
#include <TStyle.h>		    
#include <TCanvas.h>		    
#include <TGraphErrors.h>	    
#include <TGraphAsymmErrors.h>	    
#include <TLegend.h>		    
#include <TLatex.h>		    
#include <TLine.h>		    
// histos				    
#include <TH1.h>		    
#include <TH2.h>		    
#include <TH3.h>		    
#include <TH1F.h>		    
#include <TH2F.h>                   
#include <TF1.h>		    
// custom header
#include "readFiles_ppData.h"

//////////////////////////////////////
////////// readFiles_ppData //////////
//////////////////////////////////////

// default inputs to the routine
/////
const int defStartFile=0;
const int defEndFile=1; //inclusive boundary
const std::string defInFilelist = "../filelists/5p02TeV_HighPtJet80_9Files_debug_forests.txt";
/////
const int defRadius=4;
const std::string defJetType="PF";
/////
const std::string defOutputName = "readFiles_ppData_defOut.root"; 
const bool defDebugMode = true; 
/////
//int iii=0;
//std::cout<<"here"<<iii<<std::endl;iii++;
// the macro

int readFiles_ppData(int startfile = defStartFile , int endfile = defEndFile , std::string inFilelist = defInFilelist , std::string outfile = defOutputName ,
		     int radius = defRadius , std::string jetType = defJetType , bool debugMode = defDebugMode ){ 

  // for monitoring performance
  TStopwatch timer;  timer.Start();

  // debug mode settings+warnings
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  const bool fastDebugMode = (debugMode)&&true; //if debugMode is off, fastDebugMode shouldn't be on
  if(fastDebugMode)std::cout<<"fastDebugMode is ON"<<std::endl;

  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading files "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius =  " << radius<<std::endl;
  std::cout<<"jetType =  " << jetType<<std::endl;
  if(debugMode)std::cout<<"N trees = " << N_trees<<std::endl;
  if(debugMode)std::cout<<"N vars  = " << N_vars<<std::endl;
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  // form jetTreeName, initialize tree name array
  std::string trees[N_trees];
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  trees[0]=jetTreeName;
  if(debugMode)std::cout<<"looking at jetTree "<<trees[0]<<std::endl;
  for(int i=1;i<N_trees;++i)trees[i]=treeNames[i];
  
  // chain the trees together
  TChain* jetpp[N_trees];
  for(int t = 0;t<N_trees;++t)  jetpp[t] = new TChain( trees[t].data() );
  
  // open filelist, loop to the starting file
  std::cout<<"filelist used is "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  for(int ifile = 0;ifile<startfile;++ifile) instr_Forest >> filename_Forest;
  
  // add input files to the ppjet tree, inclusive boundary
  for(int ifile = startfile; ifile<=endfile; ++ifile){//input file loop

    // grab filename + create full path to the filename, add file to the TCh
    instr_Forest>>filename_Forest; 
    //if(atMIT) filename_Forest=hadoopDir+filename_Forest;
    //else filename_Forest=xrootdDirect+filename_Forest;

    // add files to the trees
    std::cout << "adding file " <<filename_Forest <<" to TChain..."<<std::endl;
    for(int t = 0;t<N_trees;++t){
      jetpp[t]->Add(filename_Forest.c_str());
      if(debugMode)std::cout<<"read file# "<<ifile<<" and loaded "<<trees[t]<<std::endl ;
      if(debugMode&&!(fastDebugMode))std::cout<<"this tree has Entries: "<<jetpp[t]->GetEntries()<<std::endl ;
    }
  }//end input file loop
  
  // friend all trees to jet analyzer tree in jetpp[0]
  for(int t = 1;t<N_trees;++t)jetpp[0]->AddFriend( jetpp[t] );
  
  // declare variables/arrays for input files
  // ak${radius}${jetType}JetAnalyzer
  //jets
  int nref_F;  
  float pt_F[1000];     float eta_F[1000];  float phi_F[1000]; 
  float rawpt_F[1000];  float jtpu_F[1000];
  //charged particles
  int chN_F[1000];  
  float chSum_F[1000]; float chHardSum_F[1000];  float chMax_F[1000];  //float chHardMax_F[1000];   
  //tracks
  float trkSum_F[1000];  float trkHardSum_F[1000];  float trkMax_F[1000];  //float trkHardMax_F[1000];
  //photons
  float phSum_F[1000];  float phHardSum_F[1000];  float phMax_F[1000];  //float phHardMax_F[1000];   
  //e, mu, neutral particles
  float eSum_F[1000];  float eMax_F[1000];  
  float muSum_F[1000];  float muMax_F[1000];  
  float neSum_F[1000];  float neMax_F[1000];   

  // HiEvtAnalyzer
  ULong64_t evt_F;   UInt_t run_F;   UInt_t lumi_F;   float vz_F;

  // skimanalysis
  int pBeamScrapingFilter_F;   int pHBHENoiseFilter_F;   int pprimaryvertexFilter_F;

  // hltanalysis
  //L1, one prescale branch per L1 path we're curious about
  int jet40_l1seed_p_F;  int jet60_l1seed_p_F;  int jet80_l1seed_p_F;  int jet100_l1seed_p_F;  //prescale
  //HLT, two branches, one decision one prescale, per HLT path we're curious about
  int jet40_F;           int jet60_F;           int jet80_F;           int jet100_F;           //decision
  int jet40_p_F;         int jet60_p_F;         int jet80_p_F;         int jet100_p_F;         //prescale

  // four trees, four specific HLT trigger pt objects
  //currently (5/26/16): HLT_AK${radius}CaloJet${trgObjpt}_Eta${}_v
  std::vector<float> *trgObjpt_40  = 0;   
  std::vector<float> *trgObjpt_60  = 0;   
  std::vector<float> *trgObjpt_80  = 0;   
  std::vector<float> *trgObjpt_100 = 0;

  // set branch addresses for the input file treesp, map them to above variables
  //ak${radius}${jetType}JetAnalyzer
  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_F);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jtpu",&jtpu_F);
  //charged particles
  jetpp[0]->SetBranchAddress("chargedN",&chN_F);
  jetpp[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[0]->SetBranchAddress("chargedHardSum",&chHardSum_F);
  jetpp[0]->SetBranchAddress("chargedMax",&chMax_F);  //  jetpp[0]->SetBranchAddress("chargedHardMax",&chMax_F);
  //tracks
  jetpp[0]->SetBranchAddress("trackSum",&trkSum_F);
  jetpp[0]->SetBranchAddress("trackHardSum",&trkHardSum_F);
  jetpp[0]->SetBranchAddress("trackMax",&trkMax_F);  //  jetpp[0]->SetBranchAddress("trackHardMax",&trkMax_F);
  //photons
  jetpp[0]->SetBranchAddress("photonSum",&phSum_F);
  jetpp[0]->SetBranchAddress("photonHardSum",&phHardSum_F);
  jetpp[0]->SetBranchAddress("photonMax",&phMax_F);  //  jetpp[0]->SetBranchAddress("photonHardMax",&phMax_F);
  //leptons
  jetpp[0]->SetBranchAddress("eSum",&eSum_F);
  jetpp[0]->SetBranchAddress("eMax",&eMax_F);
  jetpp[0]->SetBranchAddress("muSum",&muSum_F);
  jetpp[0]->SetBranchAddress("muMax",&muMax_F);
  //neutral partcles assoc. w/ jet
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F);
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);

  //hiEvtAnalyzer
  jetpp[1]->SetBranchAddress("evt",&evt_F);
  jetpp[1]->SetBranchAddress("run",&run_F);
  jetpp[1]->SetBranchAddress("lumi",&lumi_F);
  jetpp[1]->SetBranchAddress("vz",&vz_F);

  //skimanalysis
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_F);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_F);

  //hltanalysis
  //L1 prescales
  jetpp[3]->SetBranchAddress("L1_SingleJet28_BptxAND_Prescl",&jet40_l1seed_p_F);
  jetpp[3]->SetBranchAddress("L1_SingleJet40_BptxAND_Prescl",&jet60_l1seed_p_F);
  jetpp[3]->SetBranchAddress("L1_SingleJet48_BptxAND_Prescl",&jet80_l1seed_p_F);
  jetpp[3]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&jet100_l1seed_p_F);
  //HLT decisions+prescales
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&jet40_F);
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1_Prescl",&jet40_p_F);
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&jet60_F);
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1_Prescl",&jet60_p_F);
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&jet80_F);
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1_Prescl",&jet80_p_F);
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1",&jet100_F);
  jetpp[3]->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1_Prescl",&jet100_p_F);
  //specific HLT path trees, one path one tree, one pt branch (see trees[] in header)
  // trgObjpt_40 is filled with jet pt?
  jetpp[4]->SetBranchAddress("pt",&trgObjpt_40);
  jetpp[5]->SetBranchAddress("pt",&trgObjpt_60);  
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_80);  
  jetpp[7]->SetBranchAddress("pt",&trgObjpt_100);

  // Declare the output File and the necessary histograms after that:
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();

  // book jet variable histograms
  TH1F *hJetQA[2][N_vars];
  for(int k = 0; k<2; ++k){
    for(int j = 0; j<N_vars; ++j){
      //jtpt and rawpt special binning for QA
      if(j<=1) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
				       Form(";%s;",var[j].c_str()),500,  0, 500 );
      //jteta special binning for QA
      else if(j==2) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
					    Form(";%s;",var[j].c_str()),100, -5, +5  );
      //jtphi special binning for QA
      else if(j==3) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
					    Form(";%s;",var[j].c_str()),100, -4, +4  );
      //same binning for all other QA, alter later if needed
      else if(j>=4) hJetQA[k][j] = new TH1F(Form("hJetQA_%dwJetID_%s",k,var[j].c_str()),
					    Form(";%s;",var[j].c_str()),200,  0, 2   );
    }
  }

  // for vertex position, can use for vz weighting
  TH1F *hVz;
  hVz = new TH1F("hVz","",200, -20, 20);

  // for triggerpt spectra before/after jetID
  // consider renaming h_ppTrgObj, 'hpp' confusing
  TH1F *hpp_TrgObj100[2];
  TH1F *hpp_TrgObj80[2];
  TH1F *hpp_TrgObj60[2];
  TH1F *hpp_TrgObj40[2];
  TH1F *hpp_TrgObjComb[2];
  TH1F *hpp_CombJetpT[2];
  
  hpp_TrgObj100[0]  = new TH1F(Form("hpp_HLT100_noJetID_R%d_%s"     , radius, etaWidth), 
			       Form("Spectra from  Jet 100 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj80[0]   = new TH1F(Form("hpp_HLT80_noJetID_R%d_%s"     , radius, etaWidth), 
			       Form("Spectra from  Jet 80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj60[0]   = new TH1F(Form("hpp_HLT60_noJetID_R%d_%s"               , radius, etaWidth), 
			       Form("Spectra from  Jet 60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj40[0]   = new TH1F(Form("hpp_HLT40_noJetID_R%d_%s"                        , radius, etaWidth), 
			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObjComb[0] = new TH1F(Form("hpp_HLTComb_noJetID_R%d_%s"    , radius, etaWidth), 
			       Form("Trig Combined Spectra R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_CombJetpT[0]  = new TH1F(Form("hpp_TrgCombTest_noJetID_R%d_%s"           , radius, etaWidth), 
			       Form("Trig Combined Spectra KurtMethod R%d %s " , radius, etaWidth), 2000, 0, 2000);

  hpp_TrgObj100[1]  = new TH1F(Form("hpp_HLT100_JetID_R%d_%s"       , radius, etaWidth), 
			       Form("Spectra from  Jet 100 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj80[1]   = new TH1F(Form("hpp_HLT80_JetID_R%d_%s"       , radius, etaWidth), 
			       Form("Spectra from  Jet 80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj60[1]   = new TH1F(Form("hpp_HLT60_JetID_R%d_%s"                 , radius, etaWidth), 
			       Form("Spectra from  Jet 60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObj40[1]   = new TH1F(Form("hpp_HLT40_JetID_R%d_%s"                          , radius, etaWidth), 
			       Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_TrgObjComb[1] = new TH1F(Form("hpp_HLTComb_JetID_R%d_%s"      , radius, etaWidth), 
			       Form("Trig Combined Spectra R%d %s " , radius, etaWidth), 2000, 0, 2000);
  hpp_CombJetpT[1]  = new TH1F(Form("hpp_TrgCombTest_JetID_R%d_%s"             , radius, etaWidth), 
			       Form("Trig Combined Spectra KurtMethod R%d %s " , radius, etaWidth), 2000, 0, 2000);


  // trigger counts for weighting later
  Long64_t NEvents = 0;
  Long64_t NEvents_100 = 0;
  Long64_t NEvents_80 = 0;
  Long64_t NEvents_60 = 0;
  Long64_t NEvents_40 = 0;    
  //Long64_t nGoodEvt = 0;
  
  Long64_t nentries = jetpp[0]->GetEntries();
  if(debugMode) nentries = 1000;  if(fastDebugMode) nentries = 10;
  std::cout<<"Running through "<<nentries<<" events"<<std::endl;  

  // now loop over the events
  for(int nEvt = 0; nEvt < nentries; ++nEvt) {//event loop   
    
    // grab entry from each tree 
    for(int i=0;i<N_trees; ++i)jetpp[i]->GetEntry(nEvt);
    if(debugMode&&nEvt%1000==0)std::cout<<"nEvt = "<<nEvt<<std::endl;
    
    // skim/HiEvtAnalysis criteria
    for(int i = 0; i<N_trees;++i) jetpp[i]->GetEntry(nEvt);
    if( pHBHENoiseFilter_F==0     || 
        pBeamScrapingFilter_F==0  || 
        pprimaryvertexFilter_F==0 ||
        fabs(vz_F)>15              ) 
      continue;

    // fill vz histo
    hVz->Fill(vz_F);
    
    bool is40  = false;
    bool is60  = false;
    bool is80  = false;
    bool is100 = false;

    // trigger decision boolean array w/ prescale + l1 decision
    bool trgDec[4]={//l1 jetpt bits
      (bool)jet40_F, 
      (bool)jet60_F, 
      (bool)jet80_F, 
      (bool)jet100_F
    };
    int treePrescl[4]={//if bit 0, prescale for that l1 bit is 0
      (jet40_p_F*jet40_l1seed_p_F), 
      (jet60_p_F*jet60_l1seed_p_F), 
      (jet80_p_F*jet80_l1seed_p_F), 
      (jet100_p_F*jet100_l1seed_p_F)
    };    
    
    // still don't know why this was here in original code
    // scared to get rid of it completely
    //int maxtrg= -1;      
    //for(int ii=4; ii>=0; ii--){
    //  if(trgDec[ii]==1){
    //	maxtrg=ii;
    //	break;
    //  }
    //}

    // prefer highPt l1 trigger for jet if it matches multiple l1jet pt criteria
    double triggerPt=0.;
    if(jet40_F){
      for(unsigned int itt=0; itt<trgObjpt_40->size(); ++itt){
	if(trgObjpt_40->at(itt) > triggerPt) triggerPt = trgObjpt_40->at(itt);
      }
    }
    if(jet60_F){
      for(unsigned int itt=0; itt<trgObjpt_60->size(); ++itt){
	if(trgObjpt_60->at(itt) > triggerPt) triggerPt = trgObjpt_60->at(itt);
      }
    }
    if(jet80_F){
      for(unsigned int itt=0; itt<trgObjpt_80->size(); ++itt){
	if(trgObjpt_80->at(itt) > triggerPt) triggerPt = trgObjpt_80->at(itt);
      }
    }
    if(jet100_F){
      for(unsigned int itt=0; itt<trgObjpt_100->size(); ++itt){
	if(trgObjpt_100->at(itt) > triggerPt) triggerPt = trgObjpt_100->at(itt);
      }
    }    

    // come up with event weight for later use, trigComb returns weight for that event, given it's HLT, prescale(s), and trigPt
    double weight_eS = trigComb(trgDec, treePrescl, triggerPt);
    
    // check trigger decisions for events + exclusivity between them
    if(trgDec[3] && triggerPt>=100                ) is100 = true;
    if(trgDec[2] && triggerPt>=80 && triggerPt<100) is80  = true;
    if(trgDec[1] && triggerPt>=60 && triggerPt<80 ) is60  = true;
    if(trgDec[0] && triggerPt>=40 && triggerPt<60 ) is40  = true;

    // count events and events per each trigger
    NEvents++;
    if(is100) NEvents_100++;
    if(is80)  NEvents_80++;
    if(is60)  NEvents_60++;
    if(is40)  NEvents_40++;
    
    // loop over jets in the event
    for(int jet = 0; jet<nref_F; ++jet){

      // pt is too low or not in barrel, toss jet
      if( fabs(eta_F[jet]) > 2 ||
          pt_F[jet] <= 15       ) continue;     
      
      float recpt = pt_F[jet];//reco pt for this jet,
      
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
      
      if(is40)  hpp_TrgObj40[0]->Fill(recpt, treePrescl[0]);//recpt here+below prescale weighted
      if(is60)  hpp_TrgObj60[0]->Fill(recpt, treePrescl[1]);
      if(is80)  hpp_TrgObj80[0]->Fill(recpt, treePrescl[2]);
      if(is100) hpp_TrgObj100[0]->Fill(recpt, treePrescl[3]);
      if(true)  hpp_CombJetpT[0]->Fill(recpt, weight_eS);//recpt here is weighted-weighted
      
      // dijet imbalance, dijet pt ratio computation
      if(jet == 0 ) {//feel like i should add &&nref>1 req
	  hJetQA[0][19]->Fill( (float) (pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]) , weight_eS ); //A_j
	  hJetQA[0][20]->Fill( (float) (pt_F[1]/pt_F[0])                   , weight_eS ); //x_j
      }

      // apply JetID
      //I should really do this with just a flag instead of a 2xN_vars array, but no time to fix it...
      if( chSum_F[jet]/recpt > 0    &&    // if  charged hadron pt frac > 0 and 
	  neSum_F[jet]/recpt < 0.99 &&    // neutral hadron pt frac < 99%   and 
	  phSum_F[jet]/recpt < 0.99 &&    // photon hadron pt frac < 99%	   and 
	  chN_F[jet] > 0             ) {  // there's a charged particle   
	
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
	
	if(is40)  hpp_TrgObj40[1]->Fill(recpt, treePrescl[0]);
	if(is60)  hpp_TrgObj60[1]->Fill(recpt, treePrescl[1]);
	if(is80)  hpp_TrgObj80[1]->Fill(recpt, treePrescl[2]);
	if(is100) hpp_TrgObj100[1]->Fill(recpt, treePrescl[3]);
	if(true)hpp_CombJetpT[1]->Fill(recpt, weight_eS);
	
	if(jet == 0){//feel like i should add &&nref>1 req
	  hJetQA[1][19]->Fill((float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]), weight_eS);
	  hJetQA[1][20]->Fill((float)(pt_F[1]/pt_F[0]), weight_eS);
	}  
      }//end jet id cut      
    }//end jet loop
  }//end event loop
  std::cout<<std::endl;

  // before JetID
  hpp_TrgObjComb[0]->Add(hpp_TrgObj100[0]);
  hpp_TrgObjComb[0]->Add(hpp_TrgObj80[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj60[0] );
  hpp_TrgObjComb[0]->Add(hpp_TrgObj40[0] );

  // after JetID
  hpp_TrgObjComb[1]->Add(hpp_TrgObj100[1]);
  hpp_TrgObjComb[1]->Add(hpp_TrgObj80[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj60[1] );
  hpp_TrgObjComb[1]->Add(hpp_TrgObj40[1] );
  
  // before JetID
  hpp_TrgObj40[0]->Print("base");
  hpp_TrgObj60[0]->Print("base");
  hpp_TrgObj80[0]->Print("base");
  hpp_TrgObj100[0]->Print("base");
  hpp_TrgObjComb[0]->Print("base");//diff between this
  hpp_CombJetpT[0]->Print("base");//and this?
  
  // after JetID
  hpp_TrgObj40[1]->Print("base");
  hpp_TrgObj60[1]->Print("base");
  hpp_TrgObj80[1]->Print("base");
  hpp_TrgObj100[1]->Print("base");
  hpp_TrgObjComb[1]->Print("base");
  hpp_CombJetpT[1]->Print("base");
  
  std::cout<<std::endl<<"Event Summary"<<std::endl;
  std::cout<<"Total Number of Events read   = "<<NEvents<<std::endl;
  std::cout<<"Total Number of HLT100 events = "<<NEvents_100<<std::endl;
  std::cout<<"Total Number of HLT80 events = "<<NEvents_80<<std::endl;
  std::cout<<"Total Number of HLT60 events = "<<NEvents_60<<std::endl;
  std::cout<<"Total Number of HLT40 events = "<<NEvents_40<<std::endl<<std::endl;
  
  std::cout<<"writing output file..."<<std::endl;
  fout->Write();

  if(debugMode)std::cout<<"clearing trgObjpt..."<<std::endl<<std::endl;
  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();
  
  timer.Stop();

  std::cout<<std::endl<<"readFiles_ppData finished."<<std::endl;
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl<<std::endl;
  
  return 0;
}// end readFiles_ppData

////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  int rStatus = -1;
  // error, not enough arguments
  if(argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readFiles_ppData.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readFiles_ppData.exe ";
    std::cout<<"<startFile> <endFile> <inputFileList> <outputFilename> ";
    std::cout<<"<jetRadius> <jetType> <debugMode>";
    std::cout<<std::endl<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readFiles_ppData();
  else{//read input argument vector
    int startfile= atoi(argv[1]); int endfile= atoi(argv[2]); std::string inputFileList=argv[3]; std::string outputFileName=argv[4];
    int jetRadius= atoi(argv[4]); std::string jetType=argv[6];     bool debug=(bool)atoi(argv[7]);      
    rStatus = readFiles_ppData( startfile, endfile, inputFileList, outputFileName,
				jetRadius, jetType, debug);
  }
  return rStatus;
}
