// custom header
#include "readForests_jetPlots.h"

// ppData switches
const bool fillDataEvtQAHists=true;

const bool fillDataJetQAHists=true;
const bool fillDataTupelJetQAHists=false&&!fillDataJetQAHists;

const bool fillDataJetTrigQAHists=true; //data-specific
const bool fillDataJetCovMatrix=true;
const bool fillDataJetSpectraRapHists=true; //other

const bool fillDataDijetHists=false;//hardly used options

//const bool fillDataJetIDHists=true;//, tightJetID=false;
//const int jetIDint=(int)fillDataJetIDHists;

const bool fillDataJetJECQAHists=false;
const bool useHLT100=false;
//const bool ResCorrTrigPtIfApp=false;
//const bool fillDataVtxTrkQAHists=true; //in the works
const bool useTupel=fillDataTupelJetQAHists;
const bool useIncJetAnalyzer=fillDataJetQAHists;


//const std::string trgCombType="Calo";
const std::string trgCombType="PF";


//// readForests_ppData_jetPlots
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData_jetPlots( std::string inFilelist , int startfile , int endfile , 
				 int radius , std::string jetType , bool debugMode ,
				 std::string outfile, float jtEtaCutLo, float jtEtaCutHi){
  
  assert(fillDataJetQAHists!=fillDataTupelJetQAHists);//just dont use them at the same time. just. dont.
  
  // for monitoring performance + debugging
  TStopwatch timer;  timer.Start();
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;

  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading filelist "<< inFilelist<<std::endl;
  std::cout<<"reading files #'s "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius = " << radius;
  std::cout<<", jetType = " << jetType;
  std::cout<<", jtEtaCutLo = "<<jtEtaCutLo<<std::endl;
  std::cout<<", jtEtaCutHi = "<<jtEtaCutHi<<std::endl;
  std::cout<<", jtPtCut = "<<jtPtCut<<std::endl;
  std::cout<<", jtPtCut_Hi = "<<jtPtCut_Hi<<std::endl;  
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  

  //check out what kind of filelist we're using
  bool filelistIsLowerJets=false; //to know if i should apply the duplicate check or not
  bool filelistIsJet80=false;  
  bool filelistIsMinBias=false;
  if( inFilelist.find("LowerJets")!=std::string::npos ) {
    filelistIsLowerJets=true; 
    std::cout<<"running LowerJets filelist!"<<std::endl; }
  else if( inFilelist.find("Jet80")    !=std::string::npos ) {
    filelistIsJet80=true;
    std::cout<<"running Jet80 filelist!"<<std::endl; }
  else if( inFilelist.find("MinBias")    !=std::string::npos ) {//minbias
    filelistIsMinBias=true;
    std::cout<<"running MinBias filelist!"<<std::endl; }
  else assert(false);
  assert( !(filelistIsJet80&&filelistIsLowerJets) );
  assert( !(filelistIsJet80&&filelistIsMinBias) );//minbias
  assert( !(filelistIsMinBias&&filelistIsLowerJets) );//minbias
  
  // jet tree name+directory formation
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;

  // initialize tree name array
  std::string trees[N_dataTrees];
  trees[0]=jetTreeName;
  for(int i=1;i<N_dataTrees;++i){
    if(i<6)trees[i]=dataTreeNames[i];
    else{
      if(trgCombType=="PF")
	trees[i]=dataTreeNames[i]+PF_HLTBitStrings[i-6]+"_v";
      else if(trgCombType=="Calo")
	trees[i]=dataTreeNames[i]+Calo_HLTBitStrings[i-6]+"_v";
      else assert(false);      
    }
  }
  if(debugMode)for(int i=0;i<N_dataTrees;++i)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
  //  assert(false);

  // declare TChains for each tree + friend them
  TChain* jetpp[N_dataTrees];
  for(int t = 0;t<N_dataTrees;++t)  { 
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0)jetpp[0]->AddFriend( jetpp[t] );      }
  
  // open filelist, add files, including startfile+endfile, to chains
  std::cout<<"opening filelist: "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  std::string lastFileAdded=""; bool filesAdded=false; 
  for(int ifile = 0; ifile<=endfile; ++ifile){
    
    // grab a filename, check startfile and end of filelist condition
    instr_Forest>>filename_Forest; 
    if(ifile<startfile){ lastFileAdded=filename_Forest; 
      continue; }
    if(filename_Forest==lastFileAdded){ std::cout<<"end of filelist!"<<std::endl; 
      break; }
    
    std::cout<<"adding file #"<<ifile; 
    if(debugMode)std::cout<<", "<<filename_Forest; 
    
    std::cout<<", to each TChain in array"<<std::endl;
    for(int t = 0;t<N_dataTrees;++t) {
      filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
      if(debugMode)
	jetpp[t]->TChain::Lookup(true);
    }
    lastFileAdded=filename_Forest;  
  }
  assert(filesAdded);//avoid segfault later
  

  // Declare the output file, declare hists after for easy write
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();
  
  
  // declare cut hists
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"      ,(std::to_string(jtPtCut)).c_str()   ,    100, 0.,100.); hJetPtCut->Fill(jtPtCut);           
  TH1F *hJetPtCut_Hi        =new TH1F("hJetPtCut_Hi"      ,(std::to_string(jtPtCut_Hi)).c_str()   ,    1000,1500.,2500.); hJetPtCut_Hi->Fill(jtPtCut_Hi);           
  TH1F *hJetEtaCutHi       =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6.  ); hJetEtaCutHi->Fill(jtEtaCutHi); 
  TH1F *hJetEtaCutLo       =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6.  ); hJetEtaCutLo->Fill(jtEtaCutLo);     
  TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);
  TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  
  
  TH1F *hTrigComb=new TH1F("TrigComb", (trgCombType).c_str(), 1,0,1); 
  
  //evt counts, 0-1 bin --> fill w/ unweighted entry, 1-2 bin --> fill w/ trigger prescale weighted entry
  TH1D *h_NEvents          = new TH1D("NEvents"          , "NEvents",         2,0.,2.); // total event count in dataset
  TH1D *h_NEvents_skipped  = new TH1D("NEvents_skipped"  , "NEvents skipped", 2,0.,2.); // # events skipped 
  TH1D *h_NEvents_read     = new TH1D("NEvents_read"     , "NEvents read",    2,0.,2.); // # events looped over (should = total - skipped)
  TH1D *h_NEvents_trigd_1  = new TH1D("NEvents_trigd_1"  , "NEvents excl trigd no quality",   2,0.,2.);  //i want this to represent # excl trigd evts regardless of event quality (for sanity)
  TH1D *h_NEvents_skimCut  = new TH1D("NEvents_skimCut"  , "NEvents read post skimCut",  2,0.,2.); // # events passing skim cut
  TH1D *h_NEvents_vzCut    = new TH1D("NEvents_vzCut"    , "NEvents read post vzCut AND skimCut", 2,0.,2.); // # events passing vz cut
  TH1D *h_NEvents_trigd_2  = new TH1D("NEvents_trigd_2", "NEvents excl trigd w quality",   2,0.,2.);  //i want this to represent # of trigd excl evts that pass event quality
  
  TH1D *h_NEvents_isMB    = new TH1D("NEvents_isMB"   , "NEvents isMB ",    2,0.,2.); // # events passing isMB, w event quality//minbias
  TH1D *h_NEvents_is40    = new TH1D("NEvents_is40"   , "NEvents is40 ",    2,0.,2.); // # events passing is40, w event quality
  TH1D *h_NEvents_is60    = new TH1D("NEvents_is60"   , "NEvents is60 ",    2,0.,2.); // similar
  TH1D *h_NEvents_is80    = new TH1D("NEvents_is80"   , "NEvents is80 ",    2,0.,2.); // similar
  TH1D *h_NEvents_is100   = new TH1D("NEvents_is100"  , "NEvents is100",   2,0.,2.); // similar  
  
  TH1D *h_NEvents_HLTMB    = new TH1D("NEvents_HLTMB"   , "NEvents HLTMB ",    2,0.,2.); // # events passing HLTMB, w event quality//minbias
  TH1D *h_NEvents_HLT40    = new TH1D("NEvents_HLT40"   , "NEvents HLT40 ",    2,0.,2.); // # events passing HLT40, w event quality
  TH1D *h_NEvents_HLT60    = new TH1D("NEvents_HLT60"   , "NEvents HLT60 ",    2,0.,2.); // similar
  TH1D *h_NEvents_HLT80    = new TH1D("NEvents_HLT80"   , "NEvents HLT80 ",    2,0.,2.); // similar
  TH1D *h_NEvents_HLT100   = new TH1D("NEvents_HLT100"  , "NEvents HLT100",   2,0.,2.); // similar  
  
  TH1D *h_NEvents_HLT40thresh    = new TH1D("NEvents_HLT40thresh"   , "NEvents HLT40 + Theshold",    2,0.,2.); // # events passing HLT40+threshold, w event quality
  TH1D *h_NEvents_HLT60thresh    = new TH1D("NEvents_HLT60thresh"   , "NEvents HLT60 + Theshold",    2,0.,2.); // similar
  TH1D *h_NEvents_HLT80thresh    = new TH1D("NEvents_HLT80thresh"   , "NEvents HLT80 + Theshold",    2,0.,2.); // similar
  TH1D *h_NEvents_HLT100thresh   = new TH1D("NEvents_HLT100thresh"  , "NEvents HLT100+ Theshold",   2,0.,2.); // similar  

  TH1D* h_NEvents_Jet80_PDexc   = new TH1D( "NEvents_Jet80_PDexc" , "Events Unique to Jet80 PD" , 2 , 0., 2.);
  TH1D* h_NEvents_LowJets_PDexc = new TH1D( "NEvents_LowJets_PDexc" , "Events Unique to LowerJets PD" , 2 , 0., 2.);
  TH1D* h_NEvents_MinBias_PDexc = new TH1D( "NEvents_MinBias_PDexc" , "Events Unique to MinBias PD" , 2 , 0., 2.);

  TH1D* h_NEvents_Jet80LowJets_NOTMinBias = new TH1D( "NEvents_Jet80LowJets_NOTMinBias" , "Events in Jet80+LowJets   PDs but not MinBias PD" , 2 , 0., 2.);
  TH1D* h_NEvents_MinBiasJet80_NOTLowJets = new TH1D( "NEvents_MinBiasJet80_NOTLowJets" , "Events in MinBias+Jet80   PDs but not LowJets PD" , 2 , 0., 2.);
  TH1D* h_NEvents_LowJetsMinBias_NOTJet80 = new TH1D( "NEvents_LowJetsMinBias_NOTJet80" , "Events in LowJets+MinBias PDs but not Jet80 PD" , 2 , 0., 2.);
  
  TH1D* h_NEvents_Jet80LowJetsMinBias = new TH1D( "NEvents_Jet80LowJetsMinBias" , "Events in Jet80+LowJets+MinBias PDs" , 2 , 0., 2.);

  
  //EvtQA, i.e. hists filled once per event
  TH1D *hVr=NULL, *hWVr=NULL, *hTrgVr_noW=NULL;
  TH1D *hVx=NULL, *hWVx=NULL, *hTrgVx_noW=NULL;
  TH1D *hVy=NULL, *hWVy=NULL, *hTrgVy_noW=NULL;
  
  
  TH1D *hVz=NULL, *hWVz=NULL, *hTrgVz_noW=NULL, *hTrgVz_wW=NULL;
  
  const float pp5TeV2k15_BS=0.2811;//in cm
  TH1D *hVz_BSshift=NULL, *hWVz_BSshift=NULL, *hTrgVz_noW_BSshift=NULL;
  
  const float pp5TeV2k15_baryshift=-0.5160;//in cm, that of PbPb run 262922, tracker alignment record TrackerAlignment_prompt
  TH1D *hVz_baryshift_pos=NULL, *hWVz_baryshift_pos=NULL, *hTrgVz_noW_baryshift_pos=NULL;
  TH1D *hVz_baryshift_neg=NULL, *hWVz_baryshift_neg=NULL, *hTrgVz_noW_baryshift_neg=NULL;
  
  TH1D* hNref=NULL, *hWNref=NULL;
  TH1D* hjetsPEvt=NULL,* hWjetsPEvt=NULL;
  TH1D* hLeadJetPt=NULL, *hLeadJetPt_wCuts=NULL;
  
  if(fillDataEvtQAHists){
    
    hVr = new TH1D("hVr","vr, no trig, no weights",  2000,0.,0.60); 
    hWVr = new TH1D("hWeightedVr","vr, trigd, with weights",  2000,0.,0.60); 
    hTrgVr_noW = new TH1D("hTriggerVr_noWeights","vr, trigd, no weights",2000,0.,0.60); 
    
    hVx = new TH1D("hVx","vx, no trig, no weights",  2000,-0.30,0.30); //bin size of .0025
    hWVx = new TH1D("hWeightedVx","vx, trigd, with weights",  2000,-0.30,0.30);    
    hTrgVx_noW = new TH1D("hTriggerVx_noWeights","vx, trigd, no weights",  2000,-0.30,0.30);    
    
    hVy = new TH1D("hVy","vy, no trig, no weights",  2000,-0.30,0.30); 
    hWVy = new TH1D("hWeightedVy","vy, trigd, with weights",  2000,-0.30,0.30);    
    hTrgVy_noW = new TH1D("hTriggerVy_noWeights","vy, trigd, no weights",  2000,-0.30,0.30);    
    
    hVz = new TH1D("hVz","vz, no trig, no weights",  1000,-25.,25.); //50/1000=5*10^-2 cm sized bins
    hWVz = new TH1D("hWeightedVz","vz, w/ weights",  1000,-25.,25.);    
    hTrgVz_noW = new TH1D("hTriggerVz_noWeights","vz, trigd, no weights",  1000,-25.,25.);  
    hTrgVz_wW = new TH1D("hTriggerVz_withWeights","vz, trigd, w/ weights",  1000,-25.,25.);  

           hVz_BSshift = new TH1D(                 "hVz_BSshift","vz-BS, no trig, no weights",  1000,-25.,25.); //50/1000=5*10^-2 cm sized bins
          hWVz_BSshift = new TH1D(         "hWeightedVz_BSshift","vz-BS, trigd, with weights",  1000,-25.,25.);    
    hTrgVz_noW_BSshift = new TH1D("hTriggerVz_noWeights_BSshift","vz-BS, trigd, no weights  ",  1000,-25.,25.);  

            hVz_baryshift_pos = new TH1D(                 "hVz_baryshift_pos","vz+barycenter, no trig, no weights",  1000,-25.,25.); 
	   hWVz_baryshift_pos = new TH1D(         "hWeightedVz_baryshift_pos","vz+barycenter, trigd, with weights",  1000,-25.,25.);    
     hTrgVz_noW_baryshift_pos = new TH1D("hTriggerVz_noWeights_baryshift_pos","vz+barycenter, trigd, no weights  ",  1000,-25.,25.);  
     
            hVz_baryshift_neg = new TH1D(                 "hVz_baryshift_neg","vz-barycenter, no trig, no weights",  1000,-25.,25.); 
	   hWVz_baryshift_neg = new TH1D(         "hWeightedVz_baryshift_neg","vz-barycenter, trigd, with weights",  1000,-25.,25.);    
     hTrgVz_noW_baryshift_neg = new TH1D("hTriggerVz_noWeights_baryshift_neg","vz-barycenter, trigd, no weights  ",  1000,-25.,25.);  
    
    hNref      = new TH1D("hNref","nref each evt",30,0,30);         
    hWNref     = new TH1D("hWNref","weighted nref each evt",30,0,30);
    //hTrgNref_noW     = new TH1D("hTrgNref_noW","trig'd evts w/o weights nref each evt",30,0,30);
    
    hjetsPEvt  = new TH1D("hjetsPEvt","numJets each evt",30,0,30);  
    hWjetsPEvt = new TH1D("hWjetsPEvt","weighted numJets each evt",30,0,30);
    //hTrgjetsPEvt_noW     = new TH1D("hTrgjetsPEvt_noW","trig'd evts w/o weights numJets each evt",30,0,30);
    
    hLeadJetPt      = new TH1D("hLeadJetPt","lead jet pt for each evt",2500,0,2500);
    hLeadJetPt_wCuts      = new TH1D("hLeadJetPt_wCuts","lead jet pt for each evt",2500,0,2500);
    
  }
  
  
  
  //jet counts

  TH1D *h_NJets          = new TH1D("NJets","NJets read", 2,0.,2.);
  TH1D *h_NJets_jtptCut  = new TH1D("NJets_jtptCut ","NJets read post jtptCut ", 2,0.,2.);
  TH1D *h_NJets_jtptCut_Hi  = new TH1D("NJets_jtptCut_Hi ","NJets read post jtptCut_Hi ", 2,0.,2.);
  TH1D *h_NJets_jtetaCut1  = new TH1D("NJets_jtetaCut1 ","NJets read post jtetaCut1 ", 2,0.,2.);
  TH1D *h_NJets_jtetaCut2  = new TH1D("NJets_jtetaCut2 ","NJets read post jtetaCut2 ", 2,0.,2.);
  TH1D *h_NJets_JetIDCut= new TH1D("NJets_JetIDCut","NJets read post JetIDCut", 2,0.,2.);



  // JEC QA
  TH1D *hjtpt[nbins_abseta]={};
  TH1D *hL2L3Res_rawpt[nbins_abseta]={};
  TH1D *hL2L3Res_rawpt_jtpt_ratio[nbins_abseta]={};
  TH1D *hjtpt_diff[nbins_abseta]={};
  TH2D *hjtpt_v_L2L3Res_rawpt[nbins_abseta]={};
  TH2D *hjtpt_v_diff_L2L3Res_rawpt[nbins_abseta]={};  
  if(fillDataJetJECQAHists){
    for(int k = 0; k<nbins_abseta; ++k){
      hjtpt[k]            = new TH1D( ("jtpt_etabin"+std::to_string(k)).c_str()               , (                                   "Corr. Raw Jet p_{T}, |y| bin "+std::to_string(k)+";HiForest Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 2000,0,2000);
      hL2L3Res_rawpt[k]            = new TH1D( ("L2L3Res_rawpt_etabin"+std::to_string(k)).c_str()               , (                                   "Corr. Raw Jet p_{T}, |y| bin "+std::to_string(k)+";Corr. Raw Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 2000,0,2000);
      hL2L3Res_rawpt_jtpt_ratio[k] = new TH1D( ("L2L3Res_rawpt_jtpt_ratio_etabin"+std::to_string(k)).c_str()    , ("Corr. Raw Jet p_{T} Spectra/HiForest Jet p_{T} Spectra, |y| bin "+std::to_string(k)+";Jet p_{T} [GeV];Ratio of Spectra").c_str(), 2500,0,2500);      
      hjtpt_diff[k]                = new TH1D( ("jtpt_diff_etabin"+std::to_string(k)).c_str()               , ("#Delta Jet p_{T}=HiForest Jet p_{t} - Corr. Raw Jet p_{T}, |y| bin "+std::to_string(k)+";#Delta Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 800,-200,200);
      
      hjtpt_v_L2L3Res_rawpt[k] = new TH2D(("jtpt_v_L2L3Res_rawpt_etabin"+std::to_string(k)).c_str()    , ("HiForest Jet p_{T} v. Corr. Raw Jet p_{T}, |y| bin "+std::to_string(k)+";HiForest Jet p_{T} [GeV];Corr. Raw Jet p_{T} [GeV]").c_str(), 2000,0,2000, 2000,0,2000);
      
      hjtpt_v_diff_L2L3Res_rawpt[k] = new TH2D(("jtpt_v_diff_L2L3Res_rawpt_etabin"+std::to_string(k)).c_str()    , ("HiForest Jet p_{T} v. #Delta Jet p_{T}, |y| bin "+std::to_string(k)+";HiForest Jet p_{T} [GeV];#Delta Jet p_{T} [GeV]").c_str(), 2000,0,2000, 800,-200,200);      
    }
  }


  TH1D *hJetQA_jtptEntries[nbins_abseta][2]={};
  if(fillDataJetQAHists || fillDataTupelJetQAHists){
    for(int k = 0; k<nbins_abseta; ++k){      
      for(int i=0; i<2; i++){
	hJetQA_jtptEntries[k][i] = new TH1D( Form("hJetQA_%dwJetID_jtptEntries_etabin%d", i,k) , "N_{Jets}, no Pscl. Weights;N_{Jets};Jet p_{T} [GeV]" , 2500,0,2500);                 
      }
    }
  }
  
  //jet QA, i.e. hists filled once per jet
  TH1D *hJetQA[nbins_abseta][N_vars][2]={};    

//  TH1D *hJetQA_ChHadFrac[nbins_abseta][2]=NULL;
//  TH1D *hJetQA_NeutralHadAndHfFrac[nbins_abseta][2]=NULL;
//  TH1D *hJetQA_ChEmFrac[nbins_abseta][2]=NULL;
//  TH1D *hJetQA_NeutralEmFrac[nbins_abseta][2]=NULL;
//  TH1D *hJetQA_ChMult[nbins_abseta][2]=NULL;
//  TH1D *hJetQA_ConstCnt[nbins_abseta][2]=NULL;
//  TH1D *hJetQA_NeuMult[nbins_abseta][2]=NULL;
  
  if(fillDataJetQAHists){        
    for(int j = 0; j<N_vars; ++j){		      
      for(int k = 0; k<nbins_abseta; ++k){      
	for(int i=0; i<2; i++){
	  //jets
	  if(var[j]=="jtpt"||var[j]=="rawpt")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;", var[j].c_str()) , 2500,0,2500);       
	  else if(var[j]=="jteta") 
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;",var[j].c_str()) , 100,-5,+5);
	  else if(var[j]=="jtphi") 
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;",var[j].c_str()) , 100,-4,+4);
	  //jetconst. counts
	  else if(var[j]=="trkN"|| var[j]=="phN"|| var[j]=="chN"|| var[j]=="neN"|| var[j]=="eN"|| var[j]=="muN")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 100,0,100);         
	  else if(var[j]=="trkHardN"|| var[j]=="phHardN"|| var[j]=="chHardN")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 100,0,100);         
	  else if(var[j]=="neuMult"|| var[j]=="chMult"|| var[j]=="numConst")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 100,0,100);         
	  //dijets
	  else if (fillDataDijetHists){
	    if(var[j]=="dphi") 
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;",var[j].c_str()) , 50,0,4);
	    else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt")
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;", var[j].c_str()), 2500,0,2500);
	    else //xj and Aj binnings
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 200,0,2);         
	  }//dijet hists
	  //jetpt fractions binnings
	  else 
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 200,0,2);         
	}//jet ID loop
      }//nvars loop 
    }//nbins_abseta loop

    //trying to line up HIN Jet ID vars w/ SMP Jet ID vars in tupel
    //for(int k = 0; k<nbins_abseta; ++k){      
    //  for(int i=0; i<2; i++){	
    //	std::string absetastr=absetabins_str[j]+ "< #||{y} <"+absetabins_str[j+1];
    //	hJetQA_ChHadFrac[k][i]=new TH1D(           Form("hJetQA_%dwJetID_ChHadFrac_etabin%d", i, k), ("%s;ChHadFrac;N_{Jets} x Pscl. Weight",absetastr.c_str()), 200,0,2);         
    //	hJetQA_NeutralHadAndHfFrac[k][i]=new TH1D( Form("hJetQA_%dwJetID_NeutralHadAndHfFrac_etabin%d", i, k), ("%s;NeutralHadAndHfFrac;N_{Jets} x Pscl. Weight",absetastr.c_str()), 200,0,2);         
    //	hJetQA_ChEmFrac[k][i]=new TH1D(            Form("hJetQA_%dwJetID_ChEmFrac_etabin%d", i, k), ("%s;ChEmFrac;N_{Jets} x Pscl. Weight",absetastr.c_str()), 200,0,2);         
    //	hJetQA_NeutralEmFrac[k][i]=new TH1D(       Form("hJetQA_%dwJetID_NeutralEmFrac_etabin%d", i, k), ("%s;NeutralEmFrac;N_{Jets} x Pscl. Weight",absetastr.c_str()), 200,0,2);         
    //	hJetQA_ChMult[k][i]=new TH1D(              Form("hJetQA_%dwJetID_ChMult_etabin%d", i, k), ("%s;ChMult;N_{Jets} x Pscl. Weight",absetastr.c_str()), 100,0,100);         
    //	hJetQA_ConstCnt[k][i]=new TH1D(            Form("hJetQA_%dwJetID_ConstCnt_etabin%d", i, k), ("%s;ConstCnt;N_{Jets} x Pscl. Weight",absetastr.c_str()), 100,0,100);         
    //	hJetQA_NeuMult[k][i]=new TH1D(             Form("hJetQA_%dwJetID_NeuMult_etabin%d", i, k), ("%s;NeuMult;N_{Jets} x Pscl. Weight",absetastr.c_str()), 100,0,100);         
    //  }
    //}
    
  }//fill datajetqa hists
  
  TH1D *hJetSpectraRap[nbins_abseta][2]={};
  if(fillDataJetSpectraRapHists)
    for(int j = 0; j<nbins_abseta; ++j){
      for(int i=0; i<2; i++){
	std::string h_Title="hJetSpectraRap_";
	//if(fillDataJetIDHists)h_Title+="wJetID_";
	if(i==1)h_Title+="wJetID_";
	h_Title+="bin"+std::to_string(j);
	std::string h_Desc;
	//      std::stringstream stream1, stream2;	
	//      stream1.precision(1); stream1 << std::fixed << absetabins[j];
	//      stream2.precision(1); stream2 << std::fixed << absetabins[j+1];
	h_Desc="JetPt Spectra for "+absetabins_str[j]+ "< #||{y} <"+absetabins_str[j+1];
	hJetSpectraRap[j][i]=new TH1D(h_Title.c_str(),h_Desc.c_str(), 2000,0,2000);  
      } 
    }
  
  //jet QA, i.e. hists filled once per jet
  TH1D *hTupelJetQA[nbins_abseta][N_tupelvars][2]={};    
  if(fillDataTupelJetQAHists)
    for(int k = 0; k<N_tupelvars; ++k){		      
      for(int j = 0; j<nbins_abseta; ++j){	
	for (int i=0; i<2; i++)	  
	  hTupelJetQA[j][k][i]= new TH1D( ("hTupelJetQA_"+std::to_string(i)+"jetID_etabin"+std::to_string(j)+"_"+tupelvar[k]).c_str(),
					  (absetabins_str[j]+ "< #||{y} <"+absetabins_str[j+1]+";"+tupelvar[k]+";N_{Jets} x Pscl. Weight").c_str(), 
					  tupel_nbins[k], tupel_xlo[k], tupel_xhi[k] );	
      }
    }
  //assert(false);

    
    
  // trigger (also jet) level

  
  // plots for the trigger jet pt spectra (THIS IS NOT THE INCL JET SPECTRA SIGNAL)
  TH1D  *hpp_IncHLT40trgPt =NULL ,   *hpp_IncHLT40trgEta =NULL;
  TH1D  *hpp_IncHLT60trgPt =NULL ,   *hpp_IncHLT60trgEta =NULL;
  TH1D  *hpp_IncHLT80trgPt =NULL ,   *hpp_IncHLT80trgEta =NULL;
  TH1D  *hpp_IncHLT100trgPt =NULL ,  *hpp_IncHLT100trgEta =NULL;
  
  TH1D  *hpp_HLT40trgPt =NULL ,   *hpp_HLT40trgEta =NULL;
  TH1D  *hpp_HLT60trgPt =NULL ,   *hpp_HLT60trgEta =NULL;
  TH1D  *hpp_HLT80trgPt =NULL ,   *hpp_HLT80trgEta =NULL;
  TH1D *hpp_HLT100trgPt =NULL ,  *hpp_HLT100trgEta =NULL;
  
  TH1D *hpp_HLTCombtrgPt=NULL, *hpp_HLTCombtrgEta=NULL;
  //TH1D *hpp_HLTCombtrgPt_wgt1=NULL;
  //TH1D *hpp_HLTCombtrgPt_wgt0=NULL;
  //TH1D *hpp_HLTCombtrgPt_weight_eS=NULL;
  if(fillDataJetTrigQAHists){
    
    //no exclusion trigger pt
    hpp_IncHLT40trgPt   = new TH1D( "IncHLT40_trgPt"  , "inc trgPt for HLT40   ", 2500, 0, 2500 ); 
    hpp_IncHLT60trgPt   = new TH1D( "IncHLT60_trgPt"  , "inc trgPt for HLT60   ", 2500, 0, 2500 ); 
    hpp_IncHLT80trgPt   = new TH1D( "IncHLT80_trgPt"  , "inc trgPt for HLT80   ", 2500, 0, 2500 ); 
    hpp_IncHLT100trgPt  = new TH1D( "IncHLT100_trgPt" , "inc trgPt for HLT100  ", 2500, 0, 2500 );       
    
    //with exclusion trigger pt
    hpp_HLT40trgPt   = new TH1D( "isHLT40_trgPt"  , "exc trgPt for HLT40   ", 2500, 0, 2500 ); 
    hpp_HLT60trgPt   = new TH1D( "isHLT60_trgPt"  , "exc trgPt for HLT60   ", 2500, 0, 2500 ); 
    hpp_HLT80trgPt   = new TH1D( "isHLT80_trgPt"  , "exc trgPt for HLT80   ", 2500, 0, 2500 ); 
    hpp_HLT100trgPt  = new TH1D( "isHLT100_trgPt" , "exc trgPt for HLT100  ", 2500, 0, 2500 );       

    hpp_HLTCombtrgPt = new TH1D( "HLTComb_trgPt"  , "trgPt for HLTComb ", 2500, 0, 2500 );       

  //  hpp_HLTCombtrgPt_wgt1 = new TH1D( "HLTComb_trgPt_wgt1"  , "trgPt for HLTComb weight=1.", 2500, 0, 2500 );       
  //  hpp_HLTCombtrgPt_wgt0 = new TH1D( "HLTComb_trgPt_wgt0"  , "trgPt for HLTComb weight=0.", 2500, 0, 2500 );       
  //  hpp_HLTCombtrgPt_weight_eS = new TH1D( "HLTComb_trgPt_weight_eS"  , "trgPt for HLTComb weight=weight_eS=0.", 2500, 0, 2500 );       
    
    //no exclusion trigger eta
    hpp_IncHLT40trgEta   = new TH1D( "IncHLT40_trgEta"  , "inc trgEta for HLT40   ", 100, -5., 5. ); 
    hpp_IncHLT60trgEta   = new TH1D( "IncHLT60_trgEta"  , "inc trgEta for HLT60   ", 100, -5., 5. ); 
    hpp_IncHLT80trgEta   = new TH1D( "IncHLT80_trgEta"  , "inc trgEta for HLT80   ", 100, -5., 5. ); 
    hpp_IncHLT100trgEta  = new TH1D( "IncHLT100_trgEta" , "inc trgEta for HLT100  ", 100, -5., 5. );       

    //with exclusion trigger eta
    hpp_HLT40trgEta   = new TH1D( "isHLT40_trgEta"  , "exc trgEta for HLT40   ", 100, -5., 5. ); 
    hpp_HLT60trgEta   = new TH1D( "isHLT60_trgEta"  , "exc trgEta for HLT60   ", 100, -5., 5. ); 
    hpp_HLT80trgEta   = new TH1D( "isHLT80_trgEta"  , "exc trgEta for HLT80   ", 100, -5., 5. ); 
    hpp_HLT100trgEta  = new TH1D( "isHLT100_trgEta" , "exc trgEta for HLT100  ", 100, -5., 5. );       

    hpp_HLTCombtrgEta = new TH1D( "HLTComb_trgEta"  , "trgEta for HLTComb ", 100, -5., 5. );       
  }



  TH1D  *hpp_incHLT40[nbins_abseta][2]={}, *hpp_incHLT60[nbins_abseta][2]={}, *hpp_incHLT80[nbins_abseta][2]={}, *hpp_incHLT100[nbins_abseta][2]={};// plots for incl spectra for each individual trigger
  TH1D *hpp_incHLTMB[nbins_abseta][2]={};
  
  TH1D  *hpp_excHLT40[nbins_abseta][2]={}, *hpp_excHLT60[nbins_abseta][2]={}, *hpp_excHLT80[nbins_abseta][2]={}, *hpp_excHLT100[nbins_abseta][2]={};// plots for trigger-exclusive inclusive jet spectra
  TH1D* hpp_excHLTMB[nbins_abseta][2]={};
  TH1D* hpp_excHLTMB_missedHLT40_etabin0=NULL;
  TH1D *hpp_HLTComb[nbins_abseta][2]={};//combination of triggers w/ exclusion    
  
  
  
  if(fillDataJetTrigQAHists)    {
    hpp_excHLTMB_missedHLT40_etabin0=new TH1D("hppt_excHLTMB_R4_etabin0_MissedHLT40", "", 2500,0,2500);
    
    std::string hJetTrigQATitleArray[]={ "hpp_HLTMB","hpp_HLT40", "hpp_HLT60",  "hpp_HLT80", "hpp_HLT100", "hpp_HLTComb"};
    std::string hJetTrigQADescArray[]={  "HLTMBJet","HLT40Jet","HLT60Jet",  "HLT80Jet", "HLT100Jet", "HLTCombo"};
    const int nJetTrigQATitles=sizeof(hJetTrigQATitleArray)/sizeof(std::string);
    
    std::string hJetExcTrigQATitleArray[]={ "hpp_excHLTMB","hpp_excHLT40", "hpp_excHLT60",  "hpp_excHLT80", "hpp_excHLT100"};
    std::string hJetExcTrigQADescArray[]={ "ExcHLTMBJet","ExcHLT40Jet","ExcHLT60Jet",  "ExcHLT80Jet", "ExcHLT100Jet"};
    const int nJetExcTrigQATitles=sizeof(hJetExcTrigQATitleArray)/sizeof(std::string);
    
    //this loop inits hists for the incl. jet *spectra* under each trigger path decision
    for(int k=0; k<nJetTrigQATitles; k++){	
      if(debugMode) std::cout<<"k=="<<k<<std::endl;	         
      for(int j=0;j<nbins_abseta;j++){    
	if(debugMode) std::cout<<"j=="<<j<<std::endl;	  	  
	for(int i=0; i<2; i++){
	  
	  //incl trigg plots
	  std::string hDesc=hJetTrigQADescArray[k];
	  if(i==1)hDesc+=" w/ JetID";
	  hDesc+= ", "+absetabins_str[j]+ "< #||{y} <"+absetabins_str[j+1];
	  if(debugMode)std::cout<<"inclusive HLT plot hDesc="<<hDesc<<std::endl;
	  
	  std::string hTitle=hJetTrigQATitleArray[k];
	  if(i==1)hTitle+="_wJetID";
	  hTitle+="_R"+std::to_string(radius)+"_etabin"+std::to_string(j);//etaWidth;	
	  if(debugMode)std::cout<<"inclusive HLT plot hTitle="<<hTitle<<std::endl;
	  
	  if     (hJetTrigQATitleArray[k]=="hpp_HLTMB"  )  hpp_incHLTMB[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetTrigQATitleArray[k]=="hpp_HLT40"  )  hpp_incHLT40[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetTrigQATitleArray[k]=="hpp_HLT60"  )  hpp_incHLT60[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetTrigQATitleArray[k]=="hpp_HLT80"  )  hpp_incHLT80[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetTrigQATitleArray[k]=="hpp_HLT100" )  hpp_incHLT100[j][i]    = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);	
	  else if(hJetTrigQATitleArray[k]=="hpp_HLTComb")  hpp_HLTComb[j][i]      = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	
	  //excl trigg plots	
	  if(k<nJetExcTrigQATitles){
	    std::string hExcTitle=hJetExcTrigQATitleArray[k];
	    //if(fillDataJetIDHists)hExcTitle+="_wJetID";
	    if(i==1)hExcTitle+="_wJetID";
	    hExcTitle+="_R"+std::to_string(radius)+"_etabin"+std::to_string(j);//etaWidth;
	    if(debugMode)std::cout<<" exclusive HLT plothExcTitle="<<hExcTitle<<std::endl;	  
	    
	    std::string hExcDesc=hJetExcTrigQADescArray[k];
	    //if(fillDataJetIDHists)hExcDesc+=" w/ JetID";
	    if(i==1)hExcDesc+=" w/ JetID";
	    hExcDesc+=", "+absetabins_str[j] +"< #||{y} <"+absetabins_str[j+1];	    
	    if(debugMode)std::cout<<" exclusive HLT plothExcDesc="<<hExcDesc<<std::endl<<std::endl;	    
	    
	    if     (hJetExcTrigQATitleArray[k]=="hpp_excHLTMB"  )  {
	      std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLTMB[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if     (hJetExcTrigQATitleArray[k]=="hpp_excHLT40"  )  {
	      std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT40[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT60"  )  {
	      std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT60[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT80"  )  {
	      std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT80[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT100" )  {
	      std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT100[j][i]    = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}	  
	  }
	}//jet id loop
      }//jettrig qa titles 
    }//absetabins
  }//end fill datajettrigqa plots
  
  
  //TH2D* hpp_covmat=NULL;
  TH2D* hpp_covmat_eta_arr[nbins_abseta]={};
  TH1D* hpp_covmat_eta_arr_helpers[nbins_abseta]={};
  bool etabin_bool_helper[nbins_abseta]={0};
  if(fillDataJetCovMatrix){
    for(int i=0;i<nbins_abseta;i++){
      std::string covmat_eta_title="Data Jet Covariance Matrix, "+absetabins_str[i]+" < #||{y} < "+absetabins_str[i+1];    
      if(debugMode)std::cout<<"covmat_eta_title =" << covmat_eta_title<< std::endl;      
      std::string covmat_eta_name="hpp_covmat_";
      //if(fillDataJetIDHists)covmat_eta_name+="wJetID_";
      covmat_eta_name+="wJetID_";
      covmat_eta_name+="etabin_"+std::to_string(i);      
      if(debugMode)std::cout<<"covmat_eta_name  =" << covmat_eta_name << std::endl;      
      hpp_covmat_eta_arr[i]=new TH2D(covmat_eta_name.c_str(), covmat_eta_title.c_str(), nbins_pt, ptbins, nbins_pt, ptbins);     
      hpp_covmat_eta_arr_helpers[i]=new TH1D(("covmat_helper_etabin"+std::to_string(i)).c_str(), "helper only", nbins_pt, ptbins);      
    }//nbinsabseta
  }//fill datajetcovmatrix
  
  // EVENT LOOP + PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  float rawpt_F[1000];
  float eta_F[1000];   float phi_F[1000];  
  
  //tracks
  int trkN_I[1000], trkHardN_I[1000]; 
  float trkSum_F[1000],  trkMax_F[1000], trkHardSum_F[1000];

  //charged hadrons
  int chN_I[1000], chHardN_I[1000]; 
  float chSum_F[1000],  chMax_F[1000], chHardSum_F[1000];
  
  //photons  
  int phN_I[1000], phHardN_I[1000]; 
  float phSum_F[1000],  phMax_F[1000], phHardSum_F[1000];

  //leptons
  int eN_I[1000], muN_I[1000];
  float eSum_F[1000];   float eMax_F[1000];        
  float muSum_F[1000];  float muMax_F[1000];    

  //neutrals
  int neN_I[1000];
  float neSum_F[1000];  float neMax_F[1000];        

  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  
  //tracks
  jetpp[0]->SetBranchAddress("trackN",&trkN_I);  
  jetpp[0]->SetBranchAddress("trackSum",&trkSum_F);  
  jetpp[0]->SetBranchAddress("trackMax",&trkMax_F);  
  jetpp[0]->SetBranchAddress("trackHardN",&trkHardN_I);  
  jetpp[0]->SetBranchAddress("trackHardSum",&trkHardSum_F);

  //charged hadrons from PF (no HF-Calo ctontributions)
  jetpp[0]->SetBranchAddress("chargedN",&chN_I);  
  jetpp[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[0]->SetBranchAddress("chargedMax",&chMax_F);  
  jetpp[0]->SetBranchAddress("chargedHardN",&chHardN_I);  
  jetpp[0]->SetBranchAddress("chargedHardSum",&chHardSum_F);

  //photons from PF (no HF-Calo contributions)
  jetpp[0]->SetBranchAddress("photonN",&phN_I);  
  jetpp[0]->SetBranchAddress("photonSum",&phSum_F);  
  jetpp[0]->SetBranchAddress("photonMax",&phMax_F);  
  jetpp[0]->SetBranchAddress("photonHardN",&phHardN_I);  
  jetpp[0]->SetBranchAddress("photonHardSum",&phHardSum_F);

  //electrons from PF
  jetpp[0]->SetBranchAddress("eN",&eN_I);  
  jetpp[0]->SetBranchAddress("eSum",&eSum_F);  
  jetpp[0]->SetBranchAddress("eMax",&eMax_F);

  //muons from PF
  jetpp[0]->SetBranchAddress("muN",&muN_I);
  jetpp[0]->SetBranchAddress("muSum",&muSum_F);
  jetpp[0]->SetBranchAddress("muMax",&muMax_F);

  //neutral hadrons from PF (no HF-Calo contributions)
  jetpp[0]->SetBranchAddress("neutralN",&neN_I);  
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F); 
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);


  // HiEvtAnalyzer
  float vx_F,vy_F,vz_F;
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I;   
  jetpp[1]->SetBranchAddress("vz",&vz_F);
  jetpp[1]->SetBranchAddress("vy",&vy_F);
  jetpp[1]->SetBranchAddress("vx",&vx_F);
  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);


  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilter_I, pprimaryvertexFilter_I,puvertexFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);
  jetpp[2]->SetBranchAddress("pVertexFilterCutGtight",&puvertexFilter_I);


  // ppTrack/trackTree
  //vtxs
  int nVtx_I;
  int nTrkVtx_I[1000];
  float xVtx_F[1000];
  float yVtx_F[1000];
  float zVtx_F[1000];

  jetpp[3]->SetBranchAddress("nVtx",&nVtx_I);
  jetpp[3]->SetBranchAddress("nTrkVtx",&nTrkVtx_I);
  jetpp[3]->SetBranchAddress("xVtx",&xVtx_F);
  jetpp[3]->SetBranchAddress("yVtx",&yVtx_F);
  jetpp[3]->SetBranchAddress("zVtx",&zVtx_F);  
  
  //trks
  int     nTrk_I[1000];
  float  trkPt_F[10000];
  float trkEta_F[10000];
  float trkPhi_F[10000];

  jetpp[3]->SetBranchAddress("nTrk"   ,   &nTrk_I  );
  jetpp[3]->SetBranchAddress("trkPt"   ,   &trkPt_F  );
  jetpp[3]->SetBranchAddress("trkEta"   ,   &trkEta_F  );
  jetpp[3]->SetBranchAddress("trkPhi"   ,   &trkPhi_F  );
  //jetpp[3]->SetBranchAddress(""   ,   &  );  

  // hltanalysis
  // L1 Trigger Branches
  std::string L1Branches[N_L1Bits], L1PresclBranches[N_L1Bits];
  for(int i=0; i<N_L1Bits; ++i) {
    L1Branches[i] = L1BitStrings[i];
    L1PresclBranches[i] = L1BitStrings[i] + "_Prescl";  
    if(debugMode){std::cout<<std::endl;
      std::cout<<"L1Branches      [i="<<i<<"]= "<<L1Branches[i]<<std::endl;
      std::cout<<"L1PresclBranches      [i="<<i<<"]= "<<L1PresclBranches[i]<<std::endl;    }    }
  
  int jet40_l1s_I, jet60_l1s_I,  jet80_l1s_I,  jet100_l1s_I;   //L1 
  int jet40_l1s_ps_I, jet60_l1s_ps_I,  jet80_l1s_ps_I,  jet100_l1s_ps_I;     
  jetpp[4]->SetBranchAddress( L1Branches[0].c_str()  , &jet40_l1s_I);
  jetpp[4]->SetBranchAddress( L1Branches[1].c_str()  , &jet60_l1s_I);
  jetpp[4]->SetBranchAddress( L1Branches[2].c_str()  , &jet80_l1s_I);
  jetpp[4]->SetBranchAddress( L1Branches[3].c_str()  , &jet100_l1s_I);
  
  jetpp[4]->SetBranchAddress( L1PresclBranches[0].c_str()  , &jet40_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[1].c_str()  , &jet60_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[2].c_str()  , &jet80_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[3].c_str()  , &jet100_l1s_ps_I);

  int mb_l1s_I;
  int mb_l1s_ps_I;
  jetpp[4]->SetBranchAddress( "L1_MinimumBiasHF1_OR", &mb_l1s_I);
  jetpp[4]->SetBranchAddress( "L1_MinimumBiasHF1_OR_Prescl", &mb_l1s_ps_I);
  
  
  
  // make sure # of L1/HLT bits is correct
  assert(N_HLTBits==N_L1Bits);  
  
  std::string PF_HLTBranches[N_HLTBits],PF_HLTPresclBranches[N_HLTBits];
  std::string Calo_HLTBranches[N_HLTBits], Calo_HLTPresclBranches[N_HLTBits];
  double HLTthresh[N_HLTBits]={-1.};
  for(int i=0; i<N_HLTBits; ++i) {
    PF_HLTBranches[i] = PF_HLTBitStrings[i] + "_v1"; 
    PF_HLTPresclBranches[i] = PF_HLTBitStrings[i] + "_v1_Prescl";   
    Calo_HLTBranches[i] = Calo_HLTBitStrings[i] + "_v1"; 
    Calo_HLTPresclBranches[i] = Calo_HLTBitStrings[i] + "_v1_Prescl";     
    if(trgCombType=="PF") HLTthresh[i]=HLTPFthresh[i];	
    else if (trgCombType=="Calo")HLTthresh[i]=HLTCalothresh[i];
    else assert(false);
    if(debugMode){
      std::cout<<std::endl;
      std::cout<<"trgCombType            =="<<trgCombType<<std::endl;
      std::cout<<"HLTthresh             [i="<<i<<"]= "<<HLTthresh[i]<<std::endl;
      std::cout<<std::endl;
      std::cout<<"PF_HLTBranches        [i="<<i<<"]= "<<PF_HLTBranches[i]<<std::endl;
      std::cout<<"PF_HLTPresclBranches  [i="<<i<<"]= "<<PF_HLTPresclBranches[i]<<std::endl;
      std::cout<<"HLTPFthresh           [i="<<i<<"]= "<<HLTPFthresh[i]<<std::endl;
      std::cout<<std::endl;
      std::cout<<"Calo_HLTBranches      [i="<<i<<"]= "<<Calo_HLTBranches[i]<<std::endl;
      std::cout<<"Calo_HLTPresclBranches[i="<<i<<"]= "<<Calo_HLTPresclBranches[i]<<std::endl;
      std::cout<<"HLTCalothresh         [i="<<i<<"]= "<<HLTCalothresh[i]<<std::endl;//<<std::endl; 
      std::cout<<std::endl; }      
  }
  
  int Jet40_I, Jet60_I, Jet80_I, Jet100_I;   //primary HLT bits, either Calo or PF
  int Jet40_p_I, Jet60_p_I, Jet80_p_I, Jet100_p_I;    //prescales

  int Jet40_2_I, Jet60_2_I, Jet80_2_I, Jet100_2_I; //the other 40/60 bits from whichever one wasn't picked (for duplicate skip)
  if(trgCombType=="PF"){
    jetpp[4]->SetBranchAddress( PF_HLTBranches[0].c_str() , &Jet40_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[1].c_str() , &Jet60_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[2].c_str() , &Jet80_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[3].c_str() , &Jet100_I);    
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[0].c_str() , &Jet40_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[1].c_str() , &Jet60_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[2].c_str() , &Jet80_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[3].c_str() , &Jet100_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[0].c_str() , &Jet40_2_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[1].c_str() , &Jet60_2_I);  
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[2].c_str() , &Jet80_2_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[3].c_str() , &Jet100_2_I);  }
  else if(trgCombType=="Calo"){
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[0].c_str() , &Jet40_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[1].c_str() , &Jet60_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[2].c_str() , &Jet80_I);
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[3].c_str() , &Jet100_I);    
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[0].c_str() , &Jet40_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[1].c_str() , &Jet60_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[2].c_str() , &Jet80_p_I);
    jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[3].c_str() , &Jet100_p_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[0].c_str() , &Jet40_2_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[1].c_str() , &Jet60_2_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[2].c_str() , &Jet80_2_I);
    jetpp[4]->SetBranchAddress( PF_HLTBranches[3].c_str() , &Jet100_2_I);  }
  else assert(false);

  int MB_HF1ORp5_I, MB_HF1ORp5_p_I;
  jetpp[4]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1"        , &MB_HF1ORp5_I);
  jetpp[4]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1_Prescl" , &MB_HF1ORp5_p_I);
  

  //tupel/EventTree
  std::vector<float> *jetPt, *jetEta, *jetPhi, *jetE, *jetID, *jetRawPt, *jetRawE;
  std::vector<float> *jetHfHadE, *jetHfEmE, *jetChHadFrac, *jetNeutralHadAndHfFrac, *jetChEmFrac, *jetNeutralEmFrac, *jetChMult, *jetConstCnt;  
  jetpp[5]->SetBranchAddress( "JetAk04Pt",                  &jetPt                  );
  jetpp[5]->SetBranchAddress( "JetAk04Eta",                 &jetEta                 );
  jetpp[5]->SetBranchAddress( "JetAk04Phi",                 &jetPhi                 );
  jetpp[5]->SetBranchAddress( "JetAk04E",                   &jetE                   );
  jetpp[5]->SetBranchAddress( "JetAk04Id",                  &jetID                  );
  jetpp[5]->SetBranchAddress( "JetAk04RawPt",               &jetRawPt               );
  jetpp[5]->SetBranchAddress( "JetAk04RawE",                &jetRawE                );
  jetpp[5]->SetBranchAddress( "JetAk04HfHadE",              &jetHfHadE              );
  jetpp[5]->SetBranchAddress( "JetAk04HfEmE",               &jetHfEmE               );
  jetpp[5]->SetBranchAddress( "JetAk04ChHadFrac",           &jetChHadFrac           );
  jetpp[5]->SetBranchAddress( "JetAk04NeutralHadAndHfFrac", &jetNeutralHadAndHfFrac );
  jetpp[5]->SetBranchAddress( "JetAk04ChEmFrac",            &jetChEmFrac            );
  jetpp[5]->SetBranchAddress( "JetAk04NeutralEmFrac",       &jetNeutralEmFrac       );
  jetpp[5]->SetBranchAddress( "JetAk04ChMult",              &jetChMult              );
  jetpp[5]->SetBranchAddress( "JetAk04ConstCnt",            &jetConstCnt            );
  





  //ONE HLT path ONE tree ONE trig obj pt branch
  std::vector<double>  *trgObjpt_40, *trgObjpt_60, *trgObjpt_80, *trgObjpt_100;
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_40);
  jetpp[7]->SetBranchAddress("pt",&trgObjpt_60);  
  jetpp[8]->SetBranchAddress("pt",&trgObjpt_80);  
  jetpp[9]->SetBranchAddress("pt",&trgObjpt_100);
  
  std::vector<double>  *trgObjeta_40, *trgObjeta_60, *trgObjeta_80, *trgObjeta_100;
  jetpp[6]->SetBranchAddress("eta",&trgObjeta_40);
  jetpp[7]->SetBranchAddress("eta",&trgObjeta_60);  
  jetpp[8]->SetBranchAddress("eta",&trgObjeta_80);  
  jetpp[9]->SetBranchAddress("eta",&trgObjeta_100);
  

  
  
  //evtcounts check
  UInt_t NEvents_jetAnalyzr=jetpp[0]->GetEntries();   // preskim event count from files
  UInt_t NEvents_skimAnalyzr=jetpp[3]->GetEntries();   // preskim event count from files
  bool jetEventsExist= (bool)NEvents_jetAnalyzr, skimEventsExist=(bool)NEvents_skimAnalyzr;
  bool EventCountsEqual= (NEvents_jetAnalyzr==NEvents_skimAnalyzr);
  if( !EventCountsEqual || 
      !jetEventsExist || !skimEventsExist ){//if this fires, we're gonna segfault
    std::cout<<"jetAnalyzr # evts  = " << NEvents_jetAnalyzr  << std::endl;
    std::cout<<"skimAnalyzr # evts = " << NEvents_skimAnalyzr << std::endl; }
  assert( jetEventsExist && skimEventsExist && EventCountsEqual );
  
  UInt_t NEvents_allFiles=NEvents_jetAnalyzr;   // preskim event count from files  
  UInt_t NEvents_toRead= NEvents_allFiles;  
  //if(debugMode)NEvents_toRead=100000;
  std::cout<<"reading "<<NEvents_toRead<<" events"<<std::endl;  
  

  UInt_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;  //trigger event counts
  UInt_t NEvents_MB=0;//minbias


  //JetEnergy Residuals
//  int etacutForResid=3;
//  if(radius==3)etacutForResid=4;

  //  //readForests_JEC_twiki
  //  TF1* L2RelJEC_arr[n_JEC_etabins]={};
  //  for(int i=0; i<n_JEC_etabins; i++){
  //    L2RelJEC_arr[i]=new TF1(("L2RecJEC_etabin"+std::to_string(i)).c_str(), L2RelativeJEC_str.c_str(),L2RelJEC_ptlo_arr[i],L2RelJEC_pthi_arr[i]);
  //    L2RelJEC_arr[i]->SetParameter(0,L2RelJEC_par0_arr[i]);
  //    L2RelJEC_arr[i]->SetParameter(1,L2RelJEC_par1_arr[i]);
  //    L2RelJEC_arr[i]->SetParameter(2,L2RelJEC_par2_arr[i]);
  //    L2RelJEC_arr[i]->SetParameter(3,L2RelJEC_par3_arr[i]);
  //    L2RelJEC_arr[i]->SetParameter(4,L2RelJEC_par4_arr[i]);
  //    L2RelJEC_arr[i]->SetParameter(5,L2RelJEC_par5_arr[i]);
  //  }
  //  //readForests_JEC_v1
  //  TF1* L3AbsJEC=new TF1 ("L3AbsJEC",L3AbsoluteJEC_str.c_str(),L3AbsJEC_ptlo, L3AbsJEC_pthi);
  //  L3AbsJEC->SetParameter(0,L3Abs_par0);
  //  L3AbsJEC->SetParameter(1,L3Abs_par1);
  //  L3AbsJEC->SetParameter(2,L3Abs_par2);
  //  L3AbsJEC->SetParameter(3,L3Abs_par3);
  //  L3AbsJEC->SetParameter(4,L3Abs_par4);
  //  L3AbsJEC->SetParameter(5,L3Abs_par5);
  //  L3AbsJEC->SetParameter(6,L3Abs_par6);
  //  //assert(false);

//  //  //readForests_JEC_v7_ak4PF_75X_dataRun2_v13
//  TF1* L2RelL3AbsJEC_arr[n_L2RelL3AbsJEC_etabins]={};
//  for(int i=0; i<n_L2RelL3AbsJEC_etabins; i++){
//    L2RelL3AbsJEC_arr[i]=new TF1(("L2RelL3AbsJEC_etabin"+std::to_string(i)).c_str(), L2RelL3AbsJEC_str.c_str(),L2RelL3AbsJEC_ptlo[i],L2RelL3AbsJEC_pthi[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(0, L2RelL3Abs_par0_arr[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(1, L2RelL3Abs_par1_arr[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(2, L2RelL3Abs_par2_arr[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(3, L2RelL3Abs_par3_arr[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(4, L2RelL3Abs_par4_arr[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(5, L2RelL3Abs_par5_arr[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(6, L2RelL3Abs_par6_arr[i]);
//    L2RelL3AbsJEC_arr[i]->SetParameter(7, L2RelL3Abs_par7_arr[i]);
//  }
//  
//  TF1* L2L3ResJEC_arr[n_L2L3ResJEC_etabins]={};
//  for(int i=0; i<n_L2L3ResJEC_etabins; i++){
//    L2L3ResJEC_arr[i]=new TF1(("L2L3ResJEC_etabin"+std::to_string(i)).c_str(), L2L3ResJEC_str.c_str(),L2L3ResJEC_ptlo,L2L3ResJEC_pthi);
//    L2L3ResJEC_arr[i]->SetParameter(0, L2L3Res_par0_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(1, L2L3Res_par1_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(2, L2L3Res_par2_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(3, L2L3Res_par3_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(4, L2L3Res_par4_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(5, L2L3Res_par5_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(6, L2L3Res_par6_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(7, L2L3Res_par7_arr[i]);
//    L2L3ResJEC_arr[i]->SetParameter(8, L2L3Res_par8_arr[i]);
//  }
  //assert(false);

  //  //readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0
  TF1* L2RelL3AbsJEC_arr[n_L2RelL3AbsJEC_etabins]={};
  TF1* L2L3ResJEC_arr[n_L2L3ResJEC_etabins]={};
  if(fillDataJetJECQAHists){
    for(int i=0; i<n_L2RelL3AbsJEC_etabins; i++){
      L2RelL3AbsJEC_arr[i]=new TF1(("L2RelL3AbsJEC_etabin"+std::to_string(i)).c_str(), L2RelL3AbsJEC_str.c_str(),L2RelL3AbsJEC_ptlo,L2RelL3AbsJEC_pthi);
      L2RelL3AbsJEC_arr[i]->SetParameter(0, L2RelL3Abs_par0_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(1, L2RelL3Abs_par1_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(2, L2RelL3Abs_par2_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(3, L2RelL3Abs_par3_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(4, L2RelL3Abs_par4_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(5, L2RelL3Abs_par5_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(6, L2RelL3Abs_par6_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(7, L2RelL3Abs_par7_arr[i]);
      L2RelL3AbsJEC_arr[i]->SetParameter(8, L2RelL3Abs_par8_arr[i]);
    }
    

    for(int i=0; i<n_L2L3ResJEC_etabins; i++){
      L2L3ResJEC_arr[i]=new TF1(("L2L3ResJEC_etabin"+std::to_string(i)).c_str(), L2L3ResJEC_str.c_str(),L2L3ResJEC_ptlo,L2L3ResJEC_pthi);
      L2L3ResJEC_arr[i]->SetParameter(0, L2L3Res_par0_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(1, L2L3Res_par1_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(2, L2L3Res_par2_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(3, L2L3Res_par3_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(4, L2L3Res_par4_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(5, L2L3Res_par5_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(6, L2L3Res_par6_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(7, L2L3Res_par7_arr[i]);
      L2L3ResJEC_arr[i]->SetParameter(8, L2L3Res_par8_arr[i]);
    }
  }


  
  //begin event loop
  int NJets_largeJECdiff=0;
  //int assertcount=0;
  //float currentIntegral=0.;
  //int problemCount=0;
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {
    
    //grab an entry
    if( debugMode && (nEvt%10000==0) ) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if (nEvt%100000==0) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    jetpp[0]->GetEntry(nEvt);
    
    bool trgDec[N_HLTBits]={ (bool) Jet40_I , (bool)  Jet60_I ,
			     (bool) Jet80_I , (bool) Jet100_I };    
    int trgPscl[N_HLTBits]={ Jet40_p_I * jet40_l1s_ps_I ,  Jet60_p_I *  jet60_l1s_ps_I ,   
                             Jet80_p_I * jet80_l1s_ps_I , Jet100_p_I * jet100_l1s_ps_I };    

    bool trgDec2[N_HLTBits]={ (bool) Jet40_2_I , (bool)  Jet60_2_I ,
			      (bool) Jet80_2_I , (bool) Jet100_2_I };    
    
    //    bool MBtrgDec= (
    //		    (bool) ( MB_HF1ORp5_I * mb_l1s_I )
    //		    ) && 
    //      !( ((bool) Jet40_p_I) || ((bool)Jet60_p_I) || ((bool)Jet80_p_I) );
    //bool MBtrgDec =false;// ( (bool) ( MB_HF1ORp5_I * mb_l1s_I ) );
    //int  MBtrgPscl=0.;// MB_HF1ORp5_p_I * mb_l1s_ps_I;//set to -1 if it's not the minBias filelist
    //if(filelistIsMinBias){
    bool MBtrgDec= ( (bool) ( MB_HF1ORp5_I));// * mb_l1s_I ) );
    int  MBtrgPscl= MB_HF1ORp5_p_I * mb_l1s_ps_I;//set to -1 if it's not the minBias filelist
    //}

    bool isInJet80PD= (trgDec[2] || trgDec2[2]);
    bool isInLowJetsPD=(trgDec[0] || trgDec[1] || trgDec2[0] || trgDec2[1]);
    bool isInMinBiasPD=(MBtrgDec);
    if(isInJet80PD && !(isInLowJetsPD || isInMinBiasPD) ) h_NEvents_Jet80_PDexc->Fill(0.,1.);
    if(isInLowJetsPD && !(isInMinBiasPD || isInJet80PD) ) h_NEvents_LowJets_PDexc->Fill(0.,1.);
    if(isInMinBiasPD && !(isInJet80PD || isInLowJetsPD ) ) h_NEvents_MinBias_PDexc->Fill(0.,1.);
    
    if((isInJet80PD && isInLowJetsPD) && !(isInMinBiasPD) ) h_NEvents_Jet80LowJets_NOTMinBias->Fill(0.,1.);
    if((isInLowJetsPD && isInMinBiasPD ) && !(isInJet80PD) ) h_NEvents_LowJetsMinBias_NOTJet80->Fill(0.,1.);
    if((isInMinBiasPD && isInJet80PD) && !(isInLowJetsPD) ) h_NEvents_MinBiasJet80_NOTLowJets->Fill(0.,1.);
    
    if(isInJet80PD && isInLowJetsPD && isInMinBiasPD ) h_NEvents_Jet80LowJetsMinBias->Fill(0.,1.);
    
    // grab largest triggerPt among HLT paths' trigger object collectiosn
    double maxTrgPt=0.,maxTrgEta=0.; 
    unsigned int trgObj40_size=trgObjpt_40->size(), trgObj60_size=trgObjpt_60->size();
    unsigned int trgObj80_size=trgObjpt_80->size(), trgObj100_size=trgObjpt_100->size();    
    //double trgpt_noresid, trigpt_L2Res, trgpt_L2L3Res;//TO DO; RESID CORRECTED TRG PT HISTOS
    if(trgDec[3]&&useHLT100)
      for(unsigned int itt=0; itt<trgObj100_size; ++itt){
	double trgpt=trgObjpt_100->at(itt);
	double trgeta=trgObjeta_100->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta = trgeta;
	  maxTrgPt = trgpt;
	}
      }
    
    if(trgDec[2])
      for(unsigned int itt=0; itt<trgObj80_size; ++itt){
	double trgpt=trgObjpt_80->at(itt);
	double trgeta=trgObjeta_80->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta = trgeta;
	  maxTrgPt  = trgpt;
	}	  	
      }
    
    if(trgDec[1]) 
      for(unsigned int itt=0; itt<trgObj60_size; ++itt){
	double trgpt=trgObjpt_60->at(itt);
	double trgeta=trgObjeta_60->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta =trgeta;
	  maxTrgPt  =trgpt;
	}
      }
    
    
    if(trgDec[0])  
      for(unsigned int itt=0; itt<trgObj40_size; ++itt)	{
	double trgpt=trgObjpt_40->at(itt);
	double trgeta=trgObjeta_40->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta = trgeta;
	  maxTrgPt  = trgpt;
	} 
      }
    
    double trgPt=maxTrgPt; 
    double trgEta=maxTrgEta;
    
    // check trigger decisions for events + exclusivity between them, count events, assign prescale weight
    double weight_eS=0.;
    //float weight_eS = trigComb(trgDec, trgPscl, trgPt); // trig comb function replicates the procedure below
    
    //TO DO put this is40/60/80/100 stuff in a sep void function; feed it the is-flags, the trg pt, the trgDec array, and the trgPscl array
    bool is40  = false, is60  = false, is80  = false, is100 = false;
    bool isMB=false;//minbias
    if(useHLT100){
      if(      trgDec[3] && !(trgPt<HLTthresh[3])               ) 
	{ is100 = true;  weight_eS=(double)trgPscl[3]; }
      else if( trgDec[2] && !(trgPt<HLTthresh[2])  && trgPt<HLTthresh[3] ) 
	{ is80  = true;  weight_eS=(double)trgPscl[2]; }
      else if( trgDec[1] && !(trgPt<HLTthresh[1])  && trgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=(double)trgPscl[1]; }
      else if( trgDec[0] && !(trgPt<HLTthresh[0])  && trgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=(double)trgPscl[0]; }                             
      else if( MBtrgDec ) { //i want absolutely only TRUE min bias events; ones that happen to fire a jet trigger to are not going to be what i need.	  
	if( !( ((bool) Jet40_I) || ((bool)Jet60_I) || ((bool)Jet80_I ) || ((bool) Jet40_2_I) || ((bool)Jet60_2_I) || ((bool)Jet80_2_I )  )) //PF and Calo triggers
	  { isMB=true; weight_eS=(double)MBtrgPscl; }
      }
    }
    else {
      if( trgDec[2] && !(trgPt<HLTthresh[2])  ) 
	{ is80  = true;  weight_eS=(double)trgPscl[2]; }
      else if( trgDec[1] && !(trgPt<HLTthresh[1])  && trgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=(double)trgPscl[1]; }
      else if( trgDec[0] && !(trgPt<HLTthresh[0])  && trgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=(double)trgPscl[0]; }                             
      else if( MBtrgDec ) { //i want absolutely only TRUE min bias events; ones that happen to fire a jet trigger to are not going to be what i need.	  
	if( !( ((bool) Jet40_I) || ((bool)Jet60_I) || ((bool)Jet80_I ) || ((bool) Jet40_2_I) || ((bool)Jet60_2_I) || ((bool)Jet80_2_I )  )) //PF and Calo triggers
	  { 
	    bool missedJet40Evt=false;
	    for(int jet = 0; jet<nref_I; ++jet)
	      if(pt_F[jet]>(HLTthresh[0]-15.))
		{missedJet40Evt=true;break; }	   
	    if(!missedJet40Evt)
	      { isMB=true; weight_eS=(double)MBtrgPscl; }
	    else {
	      for(int jet = 0; jet<nref_I; ++jet)
		if(fabs(eta_F[jet]<0.5))
		  hpp_excHLTMB_missedHLT40_etabin0->Fill(pt_F[jet],MBtrgPscl);
	    }
	  }
      }
    }
    
    

    
    //if(debugMode) std::cout<<"weight_eS="<<weight_eS<<std::endl;    
    bool firedJetTrigger=(is40 || is60 || is80 || is100 );
    bool firedMBTrigger=(isMB);
    bool firedTrigger=firedJetTrigger || firedMBTrigger;
    
    //bool firedTrigger=!((bool)(weight_eS<1.)) ;//including exclusivity req. 
    // if an event doesn't exclusively fire the trigger, don't skip it yet; there's histograms below that must be filled first
    

    //event counts, pre-duplicate skip, pre-evt cuts, pre-trigger reqs.
    h_NEvents->Fill(0.);   h_NEvents->Fill(1.,weight_eS);      
    
    //duplicate skipping between LowerJets and Jet80
    if(filelistIsMinBias)
      if( (bool)Jet40_I   || (bool)Jet60_I   || (bool)Jet80_I   || 
	  (bool)Jet40_2_I || (bool)Jet60_2_I || (bool)Jet80_2_I    ) {
	if(debugMode)std::cout<<"this event is in Jet80 AND/OR LowerJets dataset!"<<std::endl;
	if(debugMode)std::cout<<"Skipping event, will read it in Jet80 OR LowerJets instead!"<<std::endl;
	h_NEvents_skipped->Fill(0.);		h_NEvents_skipped->Fill(1.,weight_eS);  
	continue;
      }
    
    if(filelistIsJet80)//EXCLUSION FOR READING MIN BIAS DATA MAY NEED WORK TODO
      if( (bool)Jet40_I || (bool)Jet60_I || (bool)Jet40_2_I || (bool)Jet60_2_I ) {
	//if(debugMode)std::cout<<"this event is in Jet80 AND LowerJets dataset!"<<std::endl;
	//if(debugMode)std::cout<<"Skipping event, will read it in LowerJets instead!"<<std::endl;	
	h_NEvents_skipped->Fill(0.);		h_NEvents_skipped->Fill(1.,weight_eS);  
	continue;      }        
    

    
    h_NEvents_read->Fill(0.);    h_NEvents_read->Fill(1.,weight_eS);        
    
    //event counts, post-duplicate skip, pre-evt cuts, having fired a trigger
    if(firedTrigger)
      h_NEvents_trigd_1->Fill(0.);   h_NEvents_trigd_1->Fill(1.,weight_eS );    
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0    || 
	pBeamScrapingFilter_I==0 || 
	pprimaryvertexFilter_I==0  ) continue;    //	//	puvertexFilter_I==0  	) continue;    
    h_NEvents_skimCut->Fill(0.);  h_NEvents_skimCut->Fill(1.,weight_eS);
    
    if( fabs(vz_F)>24. )     continue;
    h_NEvents_vzCut->Fill(0.); h_NEvents_vzCut->Fill(1.,weight_eS);    
    
    if( firedTrigger )
      h_NEvents_trigd_2->Fill(0.);      h_NEvents_trigd_2->Fill(1.,weight_eS);        
    
    //event counts excl trigd, post dup skip + event quality
    if     ( is100 )  { NEvents_100++ ;  h_NEvents_is100->Fill(0.)  ;  h_NEvents_is100->Fill(1.,weight_eS ) ; }
    else if( is80  )  { NEvents_80++  ;  h_NEvents_is80 ->Fill(0.)  ;  h_NEvents_is80 ->Fill(1.,weight_eS ) ; }
    else if( is60  )  { NEvents_60++  ;  h_NEvents_is60 ->Fill(0.)  ;  h_NEvents_is60 ->Fill(1.,weight_eS ) ; }
    else if( is40  )  { NEvents_40++  ;  h_NEvents_is40 ->Fill(0.)  ;  h_NEvents_is40 ->Fill(1.,weight_eS ) ; }                        
    else if( isMB  )  { NEvents_MB++  ;  h_NEvents_isMB ->Fill(0.)  ;  h_NEvents_isMB ->Fill(1.,weight_eS ) ; }//minbias
    
    if(trgDec[3])   {h_NEvents_HLT100 -> Fill(0.); h_NEvents_HLT100 -> Fill(1.,(double)trgPscl[0]);}
    if(trgDec[2])   {h_NEvents_HLT80  -> Fill(0.); h_NEvents_HLT80  -> Fill(1.,(double)trgPscl[1]);}
    if(trgDec[1])   {h_NEvents_HLT60  -> Fill(0.); h_NEvents_HLT60  -> Fill(1.,(double)trgPscl[2]);}
    if(trgDec[0])   {h_NEvents_HLT40  -> Fill(0.); h_NEvents_HLT40  -> Fill(1.,(double)trgPscl[3]);}
    if(MBtrgDec  )  {h_NEvents_HLTMB  -> Fill(0.); h_NEvents_HLTMB  -> Fill(1.,(double)MBtrgPscl  );}//minbias
    
    if(trgDec[3]&&!(trgPt<HLTthresh[3])) {h_NEvents_HLT100thresh -> Fill(0.); h_NEvents_HLT100thresh -> Fill(1.,weight_eS);}
    if(trgDec[2]&&!(trgPt<HLTthresh[2])) {h_NEvents_HLT80thresh  -> Fill(0.); h_NEvents_HLT80thresh  -> Fill(1.,weight_eS);}
    if(trgDec[1]&&!(trgPt<HLTthresh[1])) {h_NEvents_HLT60thresh  -> Fill(0.); h_NEvents_HLT60thresh  -> Fill(1.,weight_eS);}
    if(trgDec[0]&&!(trgPt<HLTthresh[0])) {h_NEvents_HLT40thresh  -> Fill(0.); h_NEvents_HLT40thresh  -> Fill(1.,weight_eS);}
    
    
    if(fillDataJetTrigQAHists && firedJetTrigger){ //only want to fill these trigger jet plots if they pass all our quality criteria      

      if(trgDec[0]&&!(trgPt<HLTthresh[0]))   hpp_IncHLT40trgPt->Fill(  trgPt, (double)trgPscl[0] );
      if(trgDec[1]&&!(trgPt<HLTthresh[1]))   hpp_IncHLT60trgPt->Fill(  trgPt, (double)trgPscl[1] );
      if(trgDec[2]&&!(trgPt<HLTthresh[2]))   hpp_IncHLT80trgPt->Fill(  trgPt, (double)trgPscl[2] );
      if(trgDec[3]&&!(trgPt<HLTthresh[3]))   hpp_IncHLT100trgPt->Fill( trgPt, (double)trgPscl[3] );
      
      if(     is100 )  hpp_HLT100trgPt->Fill(   trgPt , weight_eS );   
      else if(is80  )   hpp_HLT80trgPt->Fill(   trgPt , weight_eS );
      else if(is60  )   hpp_HLT60trgPt->Fill(   trgPt , weight_eS );
      else if(is40  )   hpp_HLT40trgPt->Fill(   trgPt , weight_eS );      
      if(true)     {
      //if(firedJetTrigger)     {
	
	hpp_HLTCombtrgPt->Fill(  trgPt , weight_eS );   
	//if( hpp_HLTCombtrgPt->Integral(1,(int)HLTthresh[2]) > currentIntegral){
	//  problemCount++;
	//  std::cout<<"//--------------------//"<<std::endl;
	//  std::cout<<"WARNING!!!!"<<std::endl;
	//  std::cout<<"weight_eS="<<weight_eS<<std::endl;
	//  std::cout<<"problemCount="<<problemCount<<std::endl;
	//  std::cout<<"//--------------------//"<<std::endl;
	//}

//	if(filelistIsJet80){
//	  
//	  //if(hpp_HLTCombtrgPt->GetBinContent(1) > 0.){
//	  if(hpp_HLTCombtrgPt->Integral(1,(int)HLTthresh[2]) > currentIntegral &&
//	     trgPt>0.){
//	    //std::cout<<"hpp_HLTCombtrgPt->GetBinLowEdge(1)="<<hpp_HLTCombtrgPt->GetBinLowEdge(1)<<std::endl;
//	    //std::cout<<"hpp_HLTCombtrgPt->GetBinLowEdge("<<(int)HLTthresh[2]<<")="<<hpp_HLTCombtrgPt->GetBinLowEdge((int)HLTthresh[2])<<std::endl;
//	    std::cout<<std::endl;
//	    std::cout<<"//--------------------//"<<std::endl;
//	    std::cout<<"WARNING!!!!"<<std::endl;
//	    std::cout<<"weight_eS="<<weight_eS<<std::endl;
//	    std::cout<<"//--------------------//"<<std::endl;
//	    std::cout<<std::endl;
//	    std::cout<<"hpp_HLTCombtrgPt->Integral(1,(int)HLTthresh[2]) WAS    "<< currentIntegral<<std::endl;
//	    std::cout<<"hpp_HLTCombtrgPt->Integral(1,(int)HLTthresh[2]) IS NOW "<<  hpp_HLTCombtrgPt->Integral(1,(int)HLTthresh[2]) <<std::endl;
//	    std::cout<<std::endl;
//	    std::cout<<"trgPt="<<trgPt<<std::endl;
//	    std::cout<<"HLTthresh[2]="<<HLTthresh[2]<<std::endl;
//	    std::cout<<"is40/is60/is80/is100/isMB="<<is40<<"/"<<is60<<"/"<<is80<<"/"<<is100<<"/"<<isMB<<std::endl;
//	    std::cout<<"weight_eS="<<weight_eS<<std::endl;
//	    
//	    std::cout<<std::endl;
//	    std::cout<<"Jet100_I        ="<<Jet100_I<<std::endl;
//	    std::cout<<"Jet100_p_I      ="<<Jet100_p_I<<std::endl;
//	    std::cout<<"jet100_l1s_I    ="<<jet100_l1s_I<<std::endl;
//	    std::cout<<"jet100_l1s_ps_I ="<<jet100_l1s_ps_I<<std::endl;
//	    std::cout<<"trgObj100_size  ="<<trgObj100_size<<std::endl;
//	    for(unsigned int itt=0; itt<trgObj100_size; ++itt){
//	      std::cout<<"trgpt at "<<itt<<"="<<trgObjpt_100->at(itt)<<std::endl;
//	    }
//	    std::cout<<std::endl;
//	    std::cout<<"Jet80_I         ="<<Jet80_I<<std::endl;
//	    std::cout<<"Jet80_p_I       ="<<Jet80_p_I<<std::endl;
//	    std::cout<<"jet80_l1s_I     ="<<jet80_l1s_I<<std::endl;
//	    std::cout<<"jet80_l1s_ps_I  ="<<jet80_l1s_ps_I<<std::endl;
//	    std::cout<<"trgObj80_size   ="<<trgObj80_size<<std::endl;
//	    for(unsigned int itt=0; itt<trgObj80_size; ++itt){
//	      std::cout<<"trgpt at "<<itt<<"="<<trgObjpt_80->at(itt)<<std::endl;
//	    }
//	    std::cout<<std::endl;
//	    std::cout<<"Jet60_I         ="<<Jet60_I<<std::endl;
//	    std::cout<<"Jet60_p_I       ="<<Jet60_p_I<<std::endl;
//	    std::cout<<"jet60_l1s_I     ="<<jet60_l1s_I<<std::endl;
//	    std::cout<<"jet60_l1s_ps_I  ="<<jet60_l1s_ps_I<<std::endl;
//	    std::cout<<"trgObj60_size   ="<<trgObj60_size<<std::endl;
//	    for(unsigned int itt=0; itt<trgObj60_size; ++itt){
//	      std::cout<<"trgpt at "<<itt<<"="<<trgObjpt_60->at(itt)<<std::endl;
//	    }
//	    std::cout<<std::endl;
//	    std::cout<<"Jet40_I         ="<<Jet40_I<<std::endl;
//	    std::cout<<"Jet40_p_I       ="<<Jet40_p_I<<std::endl;
//	    std::cout<<"jet40_l1s_I     ="<<jet40_l1s_I<<std::endl;
//	    std::cout<<"jet40_l1s_ps_I  ="<<jet40_l1s_ps_I<<std::endl;
//	    std::cout<<"trgObj40_size   ="<<trgObj40_size<<std::endl;
//	    for(unsigned int itt=0; itt<trgObj40_size; ++itt){
//	      std::cout<<"trgpt at "<<itt<<"="<<trgObjpt_40->at(itt)<<std::endl;
//	    }
//	    std::cout<<std::endl;	    
//	    std::cout<<"MB_HF1ORp5_I    =" <<MB_HF1ORp5_I    <<std::endl;
//	    std::cout<<"MB_HF1ORp5_p_I  =" <<MB_HF1ORp5_p_I  <<std::endl;
//	    std::cout<<"mb_l1s_I        =" <<mb_l1s_I        <<std::endl;
//	    std::cout<<"mb_l1s_ps_I     =" <<mb_l1s_ps_I     <<std::endl;
//	    std::cout<<std::endl;
//
//	    currentIntegral=hpp_HLTCombtrgPt->Integral(1,(int)HLTthresh[2]) ;
//
//
//	    
//	    //assert(false);
//	  }
//
////	  if(trgPt<HLTthresh[2] && weight_eS>0.){
////	    std::cout<<std::endl;
////	    std::cout<<"WARNING!!!!"<<std::endl;
////	    std::cout<<"trgPt="<<trgPt<<std::endl;
////	    std::cout<<"HLTthresh[2]="<<HLTthresh[2]<<std::endl;
////	    std::cout<<"is40/is60/is80/is100/isMB="<<is40<<"/"<<is60<<"/"<<is80<<"/"<<is100<<"/"<<isMB<<std::endl;
////	    std::cout<<"weight_eS="<<weight_eS<<std::endl;
////	    hpp_HLTCombtrgPt_wgt1->Fill(  trgPt , 1. );   
////	    hpp_HLTCombtrgPt_wgt0->Fill(  trgPt , 0. );   
////	    hpp_HLTCombtrgPt_weight_eS->Fill(  trgPt , weight_eS );   
////	    //assertcount++;
////	    //if(assertcount==10)assert(false);
////	  }
//	}
      }
      
      // trgEta
      if(trgDec[0]&&!(trgPt<HLTthresh[0]))   hpp_IncHLT40trgEta->Fill(  trgEta, (double)trgPscl[0] );
      if(trgDec[1]&&!(trgPt<HLTthresh[1]))   hpp_IncHLT60trgEta->Fill(  trgEta, (double)trgPscl[1] );
      if(trgDec[2]&&!(trgPt<HLTthresh[2]))   hpp_IncHLT80trgEta->Fill(  trgEta, (double)trgPscl[2] );
      if(trgDec[3]&&!(trgPt<HLTthresh[3]))   hpp_IncHLT100trgEta->Fill( trgEta, (double)trgPscl[3] );
      
      if(     is100 )  hpp_HLT100trgEta->Fill(   trgEta , weight_eS );   
      else if(is80  )  hpp_HLT80trgEta->Fill(    trgEta , weight_eS );
      else if(is60  )  hpp_HLT60trgEta->Fill(    trgEta , weight_eS );
      else if(is40  )  hpp_HLT40trgEta->Fill(    trgEta , weight_eS );      
      if(true)         hpp_HLTCombtrgEta->Fill(  trgEta , weight_eS );                   
    }
    
    
    // fill evt vz histo
    if(fillDataEvtQAHists){
      
      hVx->Fill(vx_F);   
      hWVx->Fill(vx_F, weight_eS);  
      if(firedTrigger)hTrgVx_noW->Fill(vx_F);
      
      hVy->Fill(vy_F);   
      hWVy->Fill(vy_F, weight_eS);
      if(firedTrigger)hTrgVy_noW->Fill(vy_F);
      
      float vr_F=std::sqrt(vx_F*vx_F + vy_F*vy_F);
      hVr->Fill(vr_F);   
      hWVr->Fill(vr_F, weight_eS);  
      if(firedTrigger)hTrgVr_noW->Fill(vr_F);
      
      hVz->Fill(vz_F);                        //this will be # of events post cuts, no trigger, no weights, post skim/vz Cut 
      hWVz->Fill(vz_F, weight_eS);            //         ''       trigger events w/ weights, post skim/vz Cut
      if(firedTrigger)hTrgVz_noW->Fill(vz_F); //         ''       trigger events w/o weights post-skim/vz Cut
      if(firedTrigger)hTrgVz_wW->Fill(vz_F, weight_eS); // triggered events w/ weights post quality
      
      hVz_BSshift->Fill((vz_F-pp5TeV2k15_BS));                        
      hWVz_BSshift->Fill((vz_F-pp5TeV2k15_BS), weight_eS);            
      if(firedTrigger)hTrgVz_noW_BSshift->Fill((vz_F-pp5TeV2k15_BS)); 
      
      hVz_baryshift_pos->Fill((vz_F+pp5TeV2k15_baryshift));                        
      hWVz_baryshift_pos->Fill((vz_F+pp5TeV2k15_baryshift), weight_eS);            
      if(firedTrigger)hTrgVz_noW_baryshift_pos->Fill((vz_F+pp5TeV2k15_baryshift)); 
      
      hVz_baryshift_neg->Fill((vz_F-pp5TeV2k15_baryshift));                        
      hWVz_baryshift_neg->Fill((vz_F-pp5TeV2k15_baryshift), weight_eS);            
      if(firedTrigger)hTrgVz_noW_baryshift_neg->Fill((vz_F-pp5TeV2k15_baryshift)); 

      hNref->Fill(nref_I);
      hWNref->Fill(nref_I,weight_eS);
      //if(firedTrigger) hTrgNref_noW->Fill(nref_I);
    }
    
    //if the histograms are just be filled w/ zero weight entries... fugghetabouit! skip it! whatsamattayou??!
    if(!firedTrigger)continue;
    
    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;

    // for event counting + avoiding duplicate fills in dijet hists
    bool dijetHistsFilled=false;    
    float evt_leadJetPt=-1., evt_leadJetPt_wCuts=-1.;

    int jetsPerEvent;//=tupel_njets;    
    if(useTupel)jetsPerEvent=jetPt->size();
    else jetsPerEvent=nref_I;
    
    if(useTupel){
      int tupel_njets=jetPt->size();
      
      
      for(int jet=0; jet<(tupel_njets) ; jet++){
	
	float jtpt  = jetPt->at(jet); //this will be the jetpt pulled from the forest
	float receta = jetEta->at(jet);
	float absreceta = fabs(receta);
	float jetNeuMult=jetConstCnt->at(jet)-jetChMult->at(jet);

	// event+jet counting
	h_NJets->Fill(0.); h_NJets->Fill(1.,weight_eS);      
	
	if(jtpt>evt_leadJetPt){
	  evt_leadJetPt=jtpt;      }
	
	if( !(jtpt > jtPtCut) ){ jetsPerEvent--; continue;} //paranoia about float comparisons
	h_NJets_jtptCut->Fill(0.); h_NJets_jtptCut->Fill(1.,weight_eS);     
	
	if( !(jtpt < jtPtCut_Hi) ){ jetsPerEvent--; continue;}      
	h_NJets_jtptCut_Hi->Fill(0.); h_NJets_jtptCut_Hi->Fill(1.,weight_eS);      
	
	if( !(absreceta < 4.7) ) { jetsPerEvent--; continue;}      //this prett redundant; PDs dont seem to write jets past 4.7
	h_NJets_jtetaCut1->Fill(0.); h_NJets_jtetaCut1->Fill(1.,weight_eS);      
	
	if( absreceta < jtEtaCutLo ) { jetsPerEvent--;	continue;}
	else if( !(absreceta < jtEtaCutHi) ) { jetsPerEvent--;	continue;}      
	h_NJets_jtetaCut2->Fill(0.); h_NJets_jtetaCut2->Fill(1.,weight_eS);                        
	
	
	
	bool passesJetID=false;//this is current based on the 8 TeV criteria
	if (!(absreceta > 2.4)) 
	  passesJetID=(bool)jetID_00eta24_v2(
					     jetNeutralHadAndHfFrac->at(jet) ,
					     jetNeutralEmFrac->at(jet),
					     jetConstCnt->at(jet),
					     jetChHadFrac->at(jet),
					     jetChMult->at(jet),
					     jetChEmFrac->at(jet) );
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  passesJetID=(bool) jetID_24eta27_v2(
					      jetNeutralHadAndHfFrac->at(jet) ,	   
					      jetNeutralEmFrac->at(jet),		   
					      jetConstCnt->at(jet));               
	else if( !(absreceta>3.0) && absreceta>2.7 )
	  passesJetID=(bool) jetID_27eta30_v2(
					       jetNeutralEmFrac->at(jet),
					       jetNeuMult );	
	else  
	  passesJetID=(bool)jetID_32eta47_v2( 
					     jetNeutralEmFrac->at(jet),
					     jetNeuMult );
	
	if(passesJetID) { 
	  h_NJets_JetIDCut->Fill(0.);	h_NJets_JetIDCut->Fill(1.,weight_eS);	
	  // largest jet pt in each event
	  if(jtpt>evt_leadJetPt_wCuts){
	    evt_leadJetPt_wCuts=jtpt;      }
	}
	
	int theRapBin=-1;	
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin)
	  if( absetabins[rapbin]<=absreceta  && 		
	      absreceta<absetabins[rapbin+1]    	      ) {	    
	    theRapBin=rapbin;
	    break;
	  }
	

	
	



	
	for(int jtid=0;jtid<2;jtid++){	  
	  
	  if(jtid==1 && !passesJetID)continue;
	  
	  if( fillDataJetSpectraRapHists ) 
	    hJetSpectraRap[theRapBin][jtid]->Fill(jtpt,weight_eS);  	    	    
	  
	  if(jtid==1)
	    if(fillDataJetCovMatrix){
	      hpp_covmat_eta_arr_helpers[theRapBin]->Fill(jtpt);//,weight_eS);	    
	      etabin_bool_helper[theRapBin]=true;      }
	  
	  
	  if(fillDataJetTrigQAHists){       	
	    
	    hpp_HLTComb[theRapBin][jtid]->Fill(jtpt, weight_eS); //kurts method
	    
	    //no exclusion between triggers
	    if(MBtrgDec)  hpp_incHLTMB[theRapBin][jtid]->Fill(jtpt, (double)MBtrgPscl); //minbias
	    if(trgDec[0]&&!(trgPt<HLTthresh[0]))  hpp_incHLT40[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[0]);
	    if(trgDec[1]&&!(trgPt<HLTthresh[1]))  hpp_incHLT60[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[1]);
	    if(trgDec[2]&&!(trgPt<HLTthresh[2]))  hpp_incHLT80[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[2]);
	    if(trgDec[3]&&!(trgPt<HLTthresh[3]))  hpp_incHLT100[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[3]);  
	    
	    //with exclusion
	    if(isMB)   hpp_excHLTMB[theRapBin][jtid]->Fill( jtpt , weight_eS);	//minbias
	    else if(is40 )  hpp_excHLT40[theRapBin][jtid]->Fill( jtpt , weight_eS);
	    else if(is60 )  hpp_excHLT60[theRapBin][jtid]->Fill( jtpt , weight_eS);
	    else if(is80 )  hpp_excHLT80[theRapBin][jtid]->Fill( jtpt , weight_eS);
	    else if(is100)  hpp_excHLT100[theRapBin][jtid]->Fill( jtpt , weight_eS); 
	  }//end fill data jet trig qa hists
	  

	  if(fillDataTupelJetQAHists){
	    
	    
	    hJetQA_jtptEntries[theRapBin][jtid]->Fill(jtpt,1.);	  	
	    
	    
	    int varind=0;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetPt                  ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetEta                 ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetPhi                 ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetE                   ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetID                  ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetRawPt               ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetRawE                ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetHfHadE              ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetHfEmE               ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetChHadFrac           ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetNeutralHadAndHfFrac ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetChEmFrac            ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetNeutralEmFrac       ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetChMult              ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetConstCnt            ->at(jet),weight_eS);varind++;
	    hTupelJetQA[theRapBin][varind][jtid]->Fill(jetNeuMult                      ,weight_eS);varind++;
	  }//end fill data tupel jet qa hists
	}//end jtid loop
      }//end tupel jets loop
    }//end use tupel
    
    
    
    if(useIncJetAnalyzer){
      for(int jet = 0; jet<nref_I; ++jet){
	
	// event+jet counting
	h_NJets->Fill(0.); h_NJets->Fill(1.,weight_eS);      
	
      
	float jtpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
	float receta = eta_F[jet];
	float absreceta = fabs(receta);
	
	//old, Doga's L2L3Res corrections at Raghav's recommendation. May very well be redundant.
	//      if(doResidualCorr){
	//	
	//	recpt=L2JES->getCorrectedPt(recpt,receta);      
	//	jtpt_L2Res=recpt;
	//	//std::cout<<"jtpt_L2Res="<<jtpt_L2Res<<std::endl;
	//	
	//	recpt=L3JES->getCorrectedPt(recpt);        
	//	jtpt_L3Res=recpt;
	//	//std::cout<<"jtpt_L3Res="<<jtpt_L3Res<<std::endl;
	//
	//      }
	
	
	
	
	
	if(jtpt>evt_leadJetPt){
	  evt_leadJetPt=jtpt;      }
	
	if( !(jtpt > jtPtCut) ){ jetsPerEvent--; continue;} //paranoia about float comparisons
	h_NJets_jtptCut->Fill(0.); h_NJets_jtptCut->Fill(1.,weight_eS);     
	
	if( !(jtpt < jtPtCut_Hi) ){ jetsPerEvent--; continue;}      
	h_NJets_jtptCut_Hi->Fill(0.); h_NJets_jtptCut_Hi->Fill(1.,weight_eS);      
	
	if( !(absreceta < 4.7) ) { jetsPerEvent--; continue;}      //this prett redundant; PDs dont seem to write jets past 4.7
	h_NJets_jtetaCut1->Fill(0.); h_NJets_jtetaCut1->Fill(1.,weight_eS);      
	
	if( absreceta < jtEtaCutLo ) { jetsPerEvent--;	continue;}
	else if( !(absreceta < jtEtaCutHi) ) { jetsPerEvent--;	continue;}      
	h_NJets_jtetaCut2->Fill(0.); h_NJets_jtetaCut2->Fill(1.,weight_eS);                        
	
	
	
	// readForests_JEC_v7_ak4PF_75X_dataRun2_v13 
	// ... OR ... 
      // readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0
	float rawpt=rawpt_F[jet];      
	float L2RelL3Abs_rawpt=-1.;
	int L2RelL3AbsJECetabin=-1;
	float L2L3Res_rawpt=-1.;
	int L2L3ResJECetabin=-1;
	
	if(fillDataJetJECQAHists){
	  for(int uu=0; uu<1; uu++){//loop so that if the jet doesn't meet some criteria i can just quit while i'm ahead
	    
	    
	    //L2RelL3Abs
	    for(int i=0; i<n_L2RelL3AbsJEC_etabins;i++){
	      if(receta>L2RelL3AbsJEC_etabins[i] && receta < L2RelL3AbsJEC_etabins[i+1]) { 
		L2RelL3AbsJECetabin=i; break;}
	      else continue;
	    }
	    if(L2RelL3AbsJECetabin==-1)
	      continue;
	    
	    ////for 75X_dataRun2_v13 only
	    //if(rawpt < L2RelL3AbsJEC_ptlo[L2RelL3AbsJECetabin] || 
	    //   rawpt > L2RelL3AbsJEC_pthi[L2RelL3AbsJECetabin] 
	    //   ) continue;
	    
	    L2RelL3Abs_rawpt=rawpt*(L2RelL3AbsJEC_arr[L2RelL3AbsJECetabin]->Eval(rawpt));		
	    
	    //L2L3Res
	    for(int i=0; i<n_L2L3ResJEC_etabins;i++){
	      if(receta>L2L3ResJEC_etabins[i] && receta < L2L3ResJEC_etabins[i+1]) { 
		L2L3ResJECetabin=i; break;}
	      else continue;
	    }
	    if(L2L3ResJECetabin==-1)
	      continue;
	    
	    L2L3Res_rawpt=L2RelL3Abs_rawpt*(L2L3ResJEC_arr[L2L3ResJECetabin]->Eval(L2RelL3Abs_rawpt));
	    
	    
	    float DIFFERENCE=jtpt-L2L3Res_rawpt;
	    if(fabs(DIFFERENCE)>1.){
	      //	  std::cout<<"/////// -------------- ///////////" <<std::endl;
	      //	  std::cout<<std::endl;
	      //	  std::cout<<"receta="<<receta<<std::endl;
	      //	  //std::cout<<"... should be between L2RelL3AbsJEC_etabins["<<L2RelL3AbsJECetabin<<"]="<<L2RelL3AbsJEC_etabins[L2RelL3AbsJECetabin]<<" and L2RelL3AbsJEC_etabins["<<L2RelL3AbsJECetabin+1<<"]="<<L2RelL3AbsJEC_etabins[L2RelL3AbsJECetabin+1]<<std::endl;	
	      //	  //std::cout<<std::endl;
	      //	  //std::cout<<"receta="<<receta<<std::endl;
	      //	  //std::cout<<"... should be between L2L3ResJEC_etabins["<<L2L3ResJECetabin<<"]="<<L2L3ResJEC_etabins[L2L3ResJECetabin]<<" and L2L3ResJEC_etabins["<<L2L3ResJECetabin+1<<"]="<<L2L3ResJEC_etabins[L2L3ResJECetabin+1]<<std::endl;   
	      //	  //std::cout<<std::endl;
	      //	  std::cout<<"                         rawpt="<<rawpt<<std::endl;
	      //	  std::cout<<"from HiForest,           jtpt="<<jtpt<< " (~=~) L2RelL3Abs_rawpt ="<<L2RelL3Abs_rawpt<<std::endl;
	      //	  std::cout<<std::endl;
	      //	  std::cout<<"from HiForest,           jtpt="<<jtpt<< " (~=~) L2L3Res_rawpt    ="<<L2L3Res_rawpt<<std::endl;  
	      //	  std::cout<<"DIFFERENCE ="<<DIFFERENCE<<std::endl;
	      //	  std::cout<<std::endl;
	      NJets_largeJECdiff++;	      //	  //assert(false);
	    }	
	    
	  }
	  
	  
	  
	}//end fill DataJetJECQAHists
	
	
	
	
	
	//float rawpt  = rawpt_F[jet];
	float recphi = phi_F[jet];
	
	int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
	int neuMult = neN_I[jet] + phN_I[jet] ;
	int numConst  = chMult + neuMult;
      
	// 13 TeV JetID criterion, loose or tight
	bool passesJetID=false;
	//float jetIDpt=jtpt;//ala HIN jetID, jtpt is corrected w/ L2/L3 residuals
	
	if (!(absreceta > 2.4)) 
	  passesJetID=(bool)jetID_00eta24( jtpt, 
					   neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet],
					   numConst,  chMult);
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  passesJetID=(bool) jetID_24eta27( jtpt,
					    neSum_F[jet],  phSum_F[jet], 
					    numConst);
	else if( !(absreceta>3.0) && absreceta>2.7 )
	  passesJetID=(bool) jetID_27eta30( jtpt,
					    neSum_F[jet],  phSum_F[jet], 
					    neuMult);
	else  
	  passesJetID=(bool)jetID_32eta47( jtpt, 
					   phSum_F[jet]);
	
	//	  if(!passesJetID) {
	//	    jetsPerEvent--;
	//	    continue; }


	
	if(passesJetID){ 
	  h_NJets_JetIDCut->Fill(0.);	h_NJets_JetIDCut->Fill(1.,weight_eS);
	  // largest jet pt in each event
	  if(jtpt>evt_leadJetPt_wCuts){
	    evt_leadJetPt_wCuts=jtpt;      }
	}
	
	
	//fill jetspectraRapHists w/ passing jetID criterion
	int theRapBin=-1;	
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin)
	  if( absetabins[rapbin]<=absreceta  && 		
	      absreceta<absetabins[rapbin+1]    	      ) {	    
	    theRapBin=rapbin;
	    break;
	  }
	

	  if(fillDataJetJECQAHists){
	    hjtpt[theRapBin]->Fill(jtpt,weight_eS);
	    hL2L3Res_rawpt[theRapBin]->Fill(L2L3Res_rawpt,weight_eS);
	    hL2L3Res_rawpt_jtpt_ratio[theRapBin]->Fill(L2L3Res_rawpt,weight_eS);//to be divided later
	    hjtpt_diff[theRapBin]->Fill((jtpt-L2L3Res_rawpt),weight_eS);
	    hjtpt_v_L2L3Res_rawpt[theRapBin]->Fill(jtpt, L2L3Res_rawpt, weight_eS);
	    hjtpt_v_diff_L2L3Res_rawpt[theRapBin]->Fill(jtpt, (jtpt-L2L3Res_rawpt), weight_eS);
	  }



	  
	  
	  
	  for(int jtid=0; jtid<2; jtid++){
	    if(jtid==1 && !passesJetID)continue;
	    
	    if( fillDataJetSpectraRapHists ) 
	      hJetSpectraRap[theRapBin][jtid]->Fill(jtpt,weight_eS);  	    	    
	    
	    if(jtid==1)
	      if(fillDataJetCovMatrix){
		hpp_covmat_eta_arr_helpers[theRapBin]->Fill(jtpt);//,weight_eS);	    
		etabin_bool_helper[theRapBin]=true;      }
	    
	    
	    // trig plots
	    if(fillDataJetTrigQAHists){
	      
	      hpp_HLTComb[theRapBin][jtid]->Fill(jtpt, weight_eS); //kurts method
	      
	      //no exclusion between triggers
	      if(MBtrgDec)  hpp_incHLTMB[theRapBin][jtid]->Fill(jtpt, (double)MBtrgPscl); //minbias
	      if(trgDec[0]&&!(trgPt<HLTthresh[0]))  hpp_incHLT40[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[0]);
	      if(trgDec[1]&&!(trgPt<HLTthresh[1]))  hpp_incHLT60[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[1]);
	      if(trgDec[2]&&!(trgPt<HLTthresh[2]))  hpp_incHLT80[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[2]);
	      if(trgDec[3]&&!(trgPt<HLTthresh[3]))  hpp_incHLT100[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[3]);  
	      
	      //with exclusion
	      if(isMB)   hpp_excHLTMB[theRapBin][jtid]->Fill( jtpt , weight_eS);	//minbias
	      else if(is40 )  hpp_excHLT40[theRapBin][jtid]->Fill( jtpt , weight_eS);
	      else if(is60 )  hpp_excHLT60[theRapBin][jtid]->Fill( jtpt , weight_eS);
	      else if(is80 )  hpp_excHLT80[theRapBin][jtid]->Fill( jtpt , weight_eS);
	      else if(is100)  hpp_excHLT100[theRapBin][jtid]->Fill( jtpt , weight_eS); 	    //assert(false);
	    }//end fill data jettrig qa hists	
	    
	    
	    
	    
	    if(fillDataJetQAHists){
	      
	      hJetQA_jtptEntries[theRapBin][jtid]->Fill(jtpt,1.);	  	
	      
	      
	      int ind=0;
	      //jets	
	      hJetQA[theRapBin][ind][jtid]->Fill(jtpt, weight_eS); ind++;
	      hJetQA[theRapBin][ind][jtid]->Fill(rawpt_F[jet], weight_eS); ind++;

	      if(jtpt>jetQAPtCut){//second jet pt cut
		hJetQA[theRapBin][ind][jtid]->Fill(eta_F[jet], weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(phi_F[jet], weight_eS); ind++;
		
		//tracks
		hJetQA[theRapBin][ind][jtid]->Fill(trkN_I[jet], weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(trkSum_F[jet]/jtpt, weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(trkMax_F[jet]/jtpt, weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(trkHardN_I[jet], weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(trkHardSum_F[jet]/jtpt, weight_eS); ind++;
		
		//PF photons
		hJetQA[theRapBin][ind][jtid]->Fill(phN_I[jet], weight_eS); ind++;       //HIN PHSUM/JTPT ~= SMP NEUTRALEMFRAC
		hJetQA[theRapBin][ind][jtid]->Fill(phSum_F[jet]/jtpt, weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(phMax_F[jet]/jtpt, weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(phHardN_I[jet], weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(phHardSum_F[jet]/jtpt, weight_eS); ind++;
		
		//PF charged hadrons
		hJetQA[theRapBin][ind][jtid]->Fill(chN_I[jet], weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(chSum_F[jet]/jtpt, weight_eS); ind++;//HIN CHSUM/JTPT ~= SMP CHHADFRAC

		hJetQA[theRapBin][ind][jtid]->Fill(chMax_F[jet]/jtpt, weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(chHardN_I[jet], weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(chHardSum_F[jet]/jtpt, weight_eS); ind++;
		
		//PF neutral hadons
		hJetQA[theRapBin][ind][jtid]->Fill(neN_I[jet], weight_eS); ind++;   
		hJetQA[theRapBin][ind][jtid]->Fill(neSum_F[jet]/jtpt, weight_eS); ind++; //HIN NESUM/JTPT ~= SMP NEUTRALHADANDHFFRAC
		hJetQA[theRapBin][ind][jtid]->Fill(neMax_F[jet]/jtpt, weight_eS); ind++;
		
		//PF electrons
		hJetQA[theRapBin][ind][jtid]->Fill(eN_I[jet], weight_eS); ind++;   
		hJetQA[theRapBin][ind][jtid]->Fill(eSum_F[jet]/jtpt, weight_eS); ind++; //HIN ESUM/JTPT ~= SMP CHEMFRAC
		hJetQA[theRapBin][ind][jtid]->Fill(eMax_F[jet]/jtpt, weight_eS); ind++;
		
		//PF muons
		hJetQA[theRapBin][ind][jtid]->Fill(muN_I[jet], weight_eS); ind++;   
		hJetQA[theRapBin][ind][jtid]->Fill(muSum_F[jet]/jtpt, weight_eS); ind++;
		hJetQA[theRapBin][ind][jtid]->Fill(muMax_F[jet]/jtpt, weight_eS); ind++;
		
		//PF particle sums
		hJetQA[theRapBin][ind][jtid]->Fill(neuMult, weight_eS); ind++;	   //HIN NEUMULT ~= SMP NEUMULT 
		hJetQA[theRapBin][ind][jtid]->Fill(chMult, weight_eS); ind++; //HIN CHMULT ~= SMP CHMULT 
		hJetQA[theRapBin][ind][jtid]->Fill(numConst, weight_eS); ind++; //HIN NUMCONST ~= SMP CONSTCNT
		
//		float ChHadFrac=chSum_F[jet]/jtpt;
//		float NeutralHadAndHfFrac=neSum_F[jet]/jtpt;//should contain HF info but doesnt, in central rapidity, doesn't matter
//		float ChEmFrac=eSum_F[jet]/jtpt;
//		float NeutralEmFrac=phSum_F[jet]/jtpt;//should contain HF info but doesnt, in central rapidity, doesn't matter
//		int ChMult=muN_I[jet]+eN_I[jet]+chN_I[jet];
//		int NeuMult=phN_I[jet]+neN_I[jet];
//		int ConstCnt=ChMult+NeuMult;
//		
//		
//		hJetQA_ChHadFrac[theRapBin][jtid]->Fill( ChHadFrac, weight_eS);
//		hJetQA_NeutralHadAndHfFrac[theRapBin][jtid]->Fill( NeutralHadAndHfFrac, weight_eS);
//		hJetQA_ChEmFrac[theRapBin][jtid]->Fill( ChEmFrac, weight_eS);
//		hJetQA_NeutralEmFrac[theRapBin][jtid]->Fill( NeutralEmFrac, weight_eS);
//		hJetQA_ChMult[theRapBin][jtid]->Fill(ChMult, weight_eS);
//		hJetQA_ConstCntMult[theRapBin][jtid]->Fill(ConstCnt, weight_eS);
//		hJetQA_NeuMult[theRapBin][jtid]->Fill(NeuMult, weight_eS);
			      



	      }
	      
	      
	      //this dijet routine probably wont work since the rapidity bin changes. CAUTION!
	      //looking for the first two good jets that meet the criteria specified
	      if(fillDataDijetHists && passesJetID){
		
		if ( !firstGoodJetFound ){
		  if(jtpt>ldJetPtCut) { firstGoodJetFound=true;
		    firstGoodJetPt =jtpt; 
		    firstGoodJetPhi=recphi; 
		  } 	  }
		
		else if ( !secondGoodJetFound && 
			  firstGoodJetFound      ) {
		  if (jtpt>subldJetPtCut) { 
		    float checkdPhi=deltaphi(firstGoodJetPhi,recphi);
		    float ptAve=(firstGoodJetPt+jtpt)/2.;
		    if( checkdPhi>dPhiCut && ptAve>ptAveCut) {
		      if(jtpt>firstGoodJetPt){secondGoodJetFound=true;
			std::cout<<std::endl<<
			  "WARNING: picked wrong jet for lead jet! Swapping..."<<std::endl<<std::endl;
			secondGoodJetPt  = firstGoodJetPt;          
			secondGoodJetPhi = firstGoodJetPhi;
			firstGoodJetPt   = jtpt;          
			firstGoodJetPhi  = recphi; }
		      else { secondGoodJetFound=true;
			secondGoodJetPt  = jtpt;
			secondGoodJetPhi = recphi; }
		    }	    }	  }
		
		if (!dijetHistsFilled &&
		    firstGoodJetFound && 
		    secondGoodJetFound ) { dijetHistsFilled=true;
		  float A_j=(firstGoodJetPt-secondGoodJetPt)
		    /(firstGoodJetPt+secondGoodJetPt);
		  float x_j=(secondGoodJetPt/firstGoodJetPt); 	
		  float dphi=deltaphi(firstGoodJetPhi,secondGoodJetPhi);
		  int ind=dijtInd;
		  hJetQA[theRapBin][ind][jtid]->Fill( A_j , weight_eS ); ind++; 
		  hJetQA[theRapBin][ind][jtid]->Fill( x_j , weight_eS ); ind++;       
		  hJetQA[theRapBin][ind][jtid]->Fill( dphi , weight_eS ); ind++;       
		  hJetQA[theRapBin][ind][jtid]->Fill( firstGoodJetPt , weight_eS ); ind++;       
		  hJetQA[theRapBin][ind][jtid]->Fill( secondGoodJetPt , weight_eS ); ind++; }
	      }//end fillDataDijetHists	      
	    } //end fill datajetqahists	    
	  }//end jtid loop
      }//end incjetanalyzer jet loop
    }//end if useIncJetAnalyzer
    
    
    if(fillDataEvtQAHists){
      hjetsPEvt ->Fill(jetsPerEvent,1.0); 
      hWjetsPEvt->Fill(jetsPerEvent,weight_eS);     
      if(evt_leadJetPt>0.) 	  hLeadJetPt->Fill(evt_leadJetPt ,weight_eS);
      if(evt_leadJetPt_wCuts>0.)  hLeadJetPt_wCuts->Fill(evt_leadJetPt_wCuts ,weight_eS);
    }
    
    
    if(fillDataJetCovMatrix){      
      for(int i=0; i<nbins_abseta; i++){
	if(etabin_bool_helper[i]){//don't undergo expensive covmatrix calc if the etabin wasn't filled.
	  fillCovMatrix( (TH2D*)hpp_covmat_eta_arr[i], (TH1D*)hpp_covmat_eta_arr_helpers[i], nbins_pt , (double) weight_eS);
	  hpp_covmat_eta_arr_helpers[i]->Reset();//reset contents of TH1 without resetting binning, min/max, etc....
	  etabin_bool_helper[i]=false;	}      }      }
    
  }//end event loop

  std::cout<<std::endl;  
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"/// Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"# of Events (criteria/description) = count / pscl weighted count" <<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"# evts passing w cut applied = count / pscl weighted count" <<std::endl;
  std::cout<<"in file(s) opened            = " << NEvents_allFiles<<std::endl; 
  std::cout<<"looped over                  = " <<h_NEvents->GetBinContent(1)          << " / " << h_NEvents->GetBinContent(2)<<std::endl; 
  std::cout<<"skipped                      = " << h_NEvents_skipped->GetBinContent(1) << " / " << h_NEvents_skipped->GetBinContent(2)<<std::endl;
  std::cout<<"read = in files - skipped    = " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;  
  std::cout<<"skimCuts  = read - !skim     = " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;//same
  std::cout<<"vzCuts    = skimCuts - !vz   = " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl; //same
  std::cout<<"excl trig w/o quality = read - !trigd        = " << h_NEvents_trigd_1->GetBinContent(1) << " / " << h_NEvents_trigd_1->GetBinContent(2)<<std::endl;//non-duplicates
  std::cout<<"excl trig w/ quality = read - !qual - !trigd = " << h_NEvents_trigd_2->GetBinContent(1) << " / " << h_NEvents_trigd_2->GetBinContent(2)<<std::endl;//close to trigd_1
  std::cout << std::endl;
  
  std::cout<<std::endl;
  std::cout<<"Trigger Counts summary, post duplicate skip, post quality (unless otherwise noted)"<<std::endl;
  std::cout<<"Trigger Combination type = " << hTrigComb->GetTitle() << std::endl;
  std::cout<<std::endl;
  std::cout<<"Pass Trigger Inclusively; Evt Counts"<<std::endl;
  std::cout<<"# evts passing HLTMB      = " << h_NEvents_HLTMB->GetBinContent(1)  << " / " << h_NEvents_HLTMB->GetBinContent(2)  <<std::endl;//minbias
  std::cout<<"# evts passing HLT40      = " << h_NEvents_HLT40->GetBinContent(1)  << " / " << h_NEvents_HLT40->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT60      = " << h_NEvents_HLT60->GetBinContent(1)  << " / " << h_NEvents_HLT60->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT80      = " << h_NEvents_HLT80->GetBinContent(1)  << " / " << h_NEvents_HLT80->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT100     = " << h_NEvents_HLT100->GetBinContent(1) << " / " << h_NEvents_HLT100->GetBinContent(2) <<std::endl;
  std::cout<<std::endl;    
  std::cout<<"Pass Trigger + Threshold Inclusively; Evt Counts"<<std::endl;
  std::cout<<"Trigger Thresholds for 40/60/80/100 = " << HLTthresh[0] <<"/ "<< HLTthresh[1]<<"/ "<< HLTthresh[2]<<"/ "<< HLTthresh[3]<<std::endl;
  std::cout<<"# evts passing HLT40 +threshold     = " << h_NEvents_HLT40thresh->GetBinContent(1)  << " / " << h_NEvents_HLT40thresh->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT60 +threshold     = " << h_NEvents_HLT60thresh->GetBinContent(1)  << " / " << h_NEvents_HLT60thresh->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT80 +threshold     = " << h_NEvents_HLT80thresh->GetBinContent(1)  << " / " << h_NEvents_HLT80thresh->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing HLT100+threshold     = " << h_NEvents_HLT100thresh->GetBinContent(1) << " / " << h_NEvents_HLT100thresh->GetBinContent(2) <<std::endl;
  std::cout<<std::endl;    
  std::cout<<"Pass Trigger Exclusively; Evt Counts"<<std::endl;
  std::cout<<"excl trigd w/o quality   = " << h_NEvents_trigd_1->GetBinContent(1) << " / " << h_NEvents_trigd_1->GetBinContent(2)<<std::endl;//non-duplicates
  std::cout<<"excl trigd w/ quality    = " << h_NEvents_trigd_2->GetBinContent(1) << " / " << h_NEvents_trigd_2->GetBinContent(2)<<std::endl;//same
  std::cout<<"# evts passing isMB      = " << h_NEvents_isMB->GetBinContent(1)  << " / " << h_NEvents_isMB->GetBinContent(2)  <<std::endl;//minbias
  std::cout<<"# evts passing is40      = " << h_NEvents_is40->GetBinContent(1)  << " / " << h_NEvents_is40->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is60      = " << h_NEvents_is60->GetBinContent(1)  << " / " << h_NEvents_is60->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is80      = " << h_NEvents_is80->GetBinContent(1)  << " / " << h_NEvents_is80->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is100     = " << h_NEvents_is100->GetBinContent(1) << " / " << h_NEvents_is100->GetBinContent(2) <<std::endl;
  std::cout<<std::endl; 

  
  std::cout<<"//////////////////////////////"<<std::endl;
  std::cout<<"/// Jet-Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////////"<<std::endl;
    
  std::cout<<std::endl;    
  if(fillDataJetJECQAHists)std::cout<<"# of Jets w/ large corrected jet p_T diff = "<< NJets_largeJECdiff << std::endl;
  std::cout<<"# of Jets read                   = " << h_NJets->GetBinContent(1) << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Jets w/ jtpt > jtptCut      = " << h_NJets_jtptCut->GetBinContent(1) << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;   
  std::cout<<"# of Jets w/ jtpt > jtptCut_Hi   = " << h_NJets_jtptCut_Hi->GetBinContent(1) << " / " << h_NJets_jtptCut_Hi->GetBinContent(2) <<std::endl;   
  std::cout<<"# of Jets w/ |jteta|<4.7         = " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Jets w/ passing jtEtaCut(s) = " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Jets w/ jetIDCut            = " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl; 

  std::cout<<std::endl<<"writing output file... "<<std::endl;
  fout->Write(); 
  
  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();
  
  trgObjeta_40->clear();
  trgObjeta_60->clear();
  trgObjeta_80->clear();
  trgObjeta_100->clear();
  
  std::cout<<"readForests_ppData_jetPlots finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  std::cout<<"output file written is  "<<outfile<<std::endl<<std::endl;
  return 0;
}    // end readForests_ppData_jetPlots


////// main //////
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1 && argc!=10){
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppData_jetPlots.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> [<absEtaCutLo> <absEtaCutHi>]";
    std::cout<<std::endl<<std::endl;
    std::cout<<"where args in [] are optional."<<std::endl;
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    std::cout<<std::endl<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppData_jetPlots();
  else if (argc==8){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppData_jetPlots( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
				  outputFileName);
  }
  else if (argc==10){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
    rStatus = readForests_ppData_jetPlots( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
					   outputFileName, jtEtaCutLo, jtEtaCutHi);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}










//TH1D *hpp_HLT40InEff =NULL;
//TH1D *hpp_HLT60InEff =NULL;
//TH1D *hpp_HLT80InEff =NULL;
//TH1D *hpp_HLT100InEff=NULL;
//TH1D *hpp_HLTAllInEff=NULL;      

////ineff
//hpp_HLT40InEff  = new TH1D( "isHLT40_inEff"  , "isHLT40  inEff Check", 1000, 0,1000 ); 
//hpp_HLT60InEff  = new TH1D( "isHLT60_inEff"  , "isHLT60  inEff Check", 1000, 0,1000 ); 
//hpp_HLT80InEff  = new TH1D( "isHLT80_inEff"  , "isHLT80  inEff Check", 1000, 0,1000 ); 
//hpp_HLT100InEff = new TH1D( "isHLT100_inEff" , "isHLT100 inEff Check", 1000, 0,1000 );       
//
//hpp_HLTAllInEff = new TH1D( "isHLTAll_inEff" , "isHLTAll inEff Check", 1000, 0,1000 );                 

//bool L1trgDec[N_HLTBits]   ={ (bool)jet40_l1s_I, (bool)jet60_l1s_I,  (bool)jet80_l1s_I,  (bool)jet100_l1s_I    };

//if(      trgDec[3] && !L1trgDec[3])  hpp_HLT100InEff->Fill((float)trgPt,weight_eS);   
//else if( trgDec[2] && !L1trgDec[2])   hpp_HLT80InEff->Fill((float)trgPt,weight_eS);
//else if( trgDec[1] && !L1trgDec[1])   hpp_HLT60InEff->Fill((float)trgPt,weight_eS);
//else if( trgDec[0] && !L1trgDec[0])   hpp_HLT40InEff->Fill((float)trgPt,weight_eS);    

//if(fillDataJetTrigQAHists){
//  
//  hpp_HLTAllInEff->Add(hpp_HLT40InEff );
//  hpp_HLTAllInEff->Add(hpp_HLT60InEff );
//  hpp_HLTAllInEff->Add(hpp_HLT80InEff );
//  hpp_HLTAllInEff->Add(hpp_HLT100InEff);       
//  
//}
 
   //float jetIDCut_neSum, jetIDCut_phSum;
  //if(tightJetID){     jetIDCut_neSum=0.90;  jetIDCut_phSum=0.90;}
  //else{     jetIDCut_neSum=0.99;  jetIDCut_phSum=0.99;}
      
