// custom header
#include "readForests_jetPlots.h"

// ppData switches
const bool fillDataEvtQAHists=true;// leave me on always

const bool useIncJetAnalyzer=true;// leave me on almost always
const bool fillDataJetQAHists=true; // leave me on almost always


const bool fillDataJetTrigQAHists=true; // leave me on almost always
const bool fillDataJetCovMatrix=true; // leave me on almost always

const bool fillDataJetJECQAHists=false;//expensive computationally, use only if needed (i.e. if someone asks about JECs + wants QA)
const bool fillDataJetJECUncHists=true&&fillDataJetQAHists;// leave me on almost always

const bool fillDataDijetHists=false;// leave me off almost always
const bool useHLT100=false; // leave me off almost always
const bool useMBevts=false; // leave me off until Ian says to use me
const bool doRunExclStudy=false; // leave me off almost always

const bool useTupel=false&&!useIncJetAnalyzer;// leave me off almost always
const bool fillDataTupelJetQAHists=false&&useTupel;// leave me off almost always

//const std::string trgCombType="Calo"; // almost never used
const std::string trgCombType="PF"; // almost always used

const bool deepDebug=false; //put to true only if in trouble...

//// readForests_ppData_jetPlots
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData_jetPlots( std::string inFilelist , int startfile , int endfile , 
				 int radius , std::string jetType , bool debugMode ,
				 std::string outfile, float jtEtaCutLo, float jtEtaCutHi){
  
  assert(!fillDataJetQAHists || !fillDataTupelJetQAHists);//just dont use them at the same time. just. dont.
  
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
    if(i<6){
      trees[i]=dataTreeNames[i];
      if(debugMode){
	if(!useTupel&&i==5)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<" [NOT BEING USED THIS JOB]"<<std::endl;  
	else std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
    }
    else{
      if(trgCombType=="PF"){
	trees[i]=dataTreeNames[i]+PF_HLTBitStrings[i-6]+"_v";
	if(debugMode)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
      else if(trgCombType=="Calo"){
	trees[i]=dataTreeNames[i]+Calo_HLTBitStrings[i-6]+"_v";
	if(debugMode)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
      else assert(false);      
    }
  }
  //if(debugMode)for(int i=0;i<N_dataTrees;++i)
  //assert(false);

  // declare TChains for each tree + friend them
  TChain* jetpp[N_dataTrees];
  for(int t = 0;t<N_dataTrees;++t)  { 
    if(!useTupel && t==5)continue;
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0){
      jetpp[0]->AddFriend( jetpp[t] );      
    }    
  }
  
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
      if(!useTupel && t==5)continue;
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

  TH1D *h_NEvents_HLTMB_AND_L1jet52    = new TH1D("NEvents_HLTMB_AND_L1jet52"   , "NEvents HLTMB AND L1jet52",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet48    = new TH1D("NEvents_HLTMB_AND_L1jet48"   , "NEvents HLTMB AND L1jet48",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet40    = new TH1D("NEvents_HLTMB_AND_L1jet40"   , "NEvents HLTMB AND L1jet40",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet28    = new TH1D("NEvents_HLTMB_AND_L1jet28"   , "NEvents HLTMB AND L1jet28",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet24    = new TH1D("NEvents_HLTMB_AND_L1jet24"   , "NEvents HLTMB AND L1jet24",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet20    = new TH1D("NEvents_HLTMB_AND_L1jet20"   , "NEvents HLTMB AND L1jet20",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet16    = new TH1D("NEvents_HLTMB_AND_L1jet16"   , "NEvents HLTMB AND L1jet16",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet12    = new TH1D("NEvents_HLTMB_AND_L1jet12"   , "NEvents HLTMB AND L1jet12",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_HLTMB_AND_L1jet8     = new TH1D("NEvents_HLTMB_AND_L1jet8 "   , "NEvents HLTMB AND L1jet8 ",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  
  
  TH1D *h_NEvents_isMB_AND_L1jet52    = new TH1D("NEvents_isMB_AND_L1jet52"   , "NEvents isMB AND L1jet52",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet48    = new TH1D("NEvents_isMB_AND_L1jet48"   , "NEvents isMB AND L1jet48",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet40    = new TH1D("NEvents_isMB_AND_L1jet40"   , "NEvents isMB AND L1jet40",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet28    = new TH1D("NEvents_isMB_AND_L1jet28"   , "NEvents isMB AND L1jet28",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet24    = new TH1D("NEvents_isMB_AND_L1jet24"   , "NEvents isMB AND L1jet24",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet20    = new TH1D("NEvents_isMB_AND_L1jet20"   , "NEvents isMB AND L1jet20",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet16    = new TH1D("NEvents_isMB_AND_L1jet16"   , "NEvents isMB AND L1jet16",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet12    = new TH1D("NEvents_isMB_AND_L1jet12"   , "NEvents isMB AND L1jet12",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias
  TH1D *h_NEvents_isMB_AND_L1jet8     = new TH1D("NEvents_isMB_AND_L1jet8 "   , "NEvents isMB AND L1jet8 ",    2,0.,2.); // # events passing HLTMB, w event quality, w/ L1Jetseed //minbias





  
  //EvtQA, i.e. hists filled once per event
  TH1D* hRuns=NULL;
  TH1D* hLumiByRun[Nruns]={};

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
    
    hRuns=new TH1D("hRuns","Run #'s;Run #;N_{evts}",  1000, 261500, 262500);
    for(int i=0; i<Nruns;i++){      
      hLumiByRun[i]=new TH1D( ("hLumi_Run"+std::to_string(runarray[i])).c_str(), 
			      ("Events By Lumi, Run #"+std::to_string(runarray[i])+";Lumi Section;N_{evts}").c_str(), 
			      1500,0,1500);
    }
    
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
  TH1D *h_NJets_jetQAPtCut  = new TH1D("NJets_jetQAPtCut ","NJets read post jettQAPtCut ", 2,0.,2.);
  TH1D *h_NJets_JetIDCut= new TH1D("NJets_JetIDCut","NJets read post JetIDCut", 2,0.,2.);



  // JEC QA
  TH1D *hjtpt[nbins_abseta]={};//forest
  TH1D *hL2L3Res_rawpt[nbins_abseta]={};//rawpt w/ L2Rel+L3Abs+L2L3Res JEC I FOUND
  TH1D *hL2L3Res_2_rawpt[nbins_abseta]={};// CHRIS'S JEC
  
  TH1D *hjtpt_L2L3Res_rawpt_diff[nbins_abseta]={};
  TH1D *hjtpt_L2L3Res_2_rawpt_diff[nbins_abseta]={};
  TH1D *hL2L3Res_rawpt_L2L3Res_2_rawpt_diff[nbins_abseta]={};
  
  TH1D *hjtpt_L2L3Res_rawpt_ratio[nbins_abseta]={};
  TH1D *hjtpt_L2L3Res_2_rawpt_ratio[nbins_abseta]={};
  TH1D *hL2L3Res_rawpt_L2L3Res_2_rawpt_ratio[nbins_abseta]={};
  
  TH2D *hjtpt_v_diff_L2L3Res_rawpt[nbins_abseta]={};
  TH2D *hjtpt_v_diff_L2L3Res_2_rawpt[nbins_abseta]={};
  TH2D *hL2L3Res_rawpt_v_diff_L2L3Res_2_rawpt[nbins_abseta]={};


  if(fillDataJetJECQAHists){
    for(int k = 0; k<nbins_abseta; ++k){
      hjtpt[k]              = new TH1D( ("jtpt_etabin"+std::to_string(k)).c_str(), 
					("HiForest Jet p_{T}, |y| bin "+std::to_string(k)+";HiForest Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(),
					2500,0,2500);
      hL2L3Res_rawpt[k]     = new TH1D(  ("L2L3Res_rawpt_etabin"+std::to_string(k)).c_str(), 
					 ("74X_dataRun2_HLT_ppAt5TeV_v0 JEC Corr. Raw Jet p_{T}, |y| bin "+std::to_string(k)+";74X_dataRun2_HLT_ppAt5TeV_v0 JEC Corr. Raw Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 
					 2500,0,2500);
      hL2L3Res_2_rawpt[k]   = new TH1D(  ("L2L3Res_2_rawpt_etabin"+std::to_string(k)).c_str(), 
					 ("75X_data_Run2_v13 JEC Corr. Raw Jet p_{T}, |y| bin "+std::to_string(k)+";75X_data_Run2_v13 JEC Corr. Raw Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 
					 2500,0,2500);

      hjtpt_L2L3Res_rawpt_diff[k]              = new TH1D(  ("jtpt_L2L3Res_rawpt_diff_etabin"+std::to_string(k)).c_str(), 
							    ("(HiForest Jet p_{T}) - (74X_dataRun2_HLT_ppAt5TeV_v0 JEC Corr. Raw Jet p_{T} [GeV]), |y| bin "+std::to_string(k)+";#Delta Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 
							    800,-200,200);
      hjtpt_L2L3Res_2_rawpt_diff[k]          = new TH1D(  ("jtpt_L2L3Res_2_rawpt_diff_etabin"+std::to_string(k)).c_str(), 
							    ("(HiForest Jet p_{T}) - (75X_data_Run2_v13 JEC Corr. Raw Jet p_{T} [GeV]), |y| bin "+std::to_string(k)+";#Delta Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 
							    800,-200,200);      
      hL2L3Res_rawpt_L2L3Res_2_rawpt_diff[k] = new TH1D(  ("L2L3Res_rawpt_L2L3Res_2_rawpt_diff_etabin"+std::to_string(k)).c_str(), 
							    ("(74X_dataRun2_HLT_ppAt5TeV_v0 JEC Corr. Raw Jet p_{T} [GeV]) - (75X_data_Run2_v13 JEC Corr. Raw Jet p_{T} [GeV]), |y| bin "+std::to_string(k)+";#Delta Jet p_{T} [GeV];N_{jets} x Pscl. Weight").c_str(), 
							    800,-200,200);      
      
      hjtpt_L2L3Res_rawpt_ratio[k]              = new TH1D(  ("jtpt_L2L3Res_rawpt_ratio_etabin"+std::to_string(k)).c_str(), 
							     ("(HiForest Jet p_{T})/(74X_dataRun2_HLT_ppAt5TeV_v0 JEC Corr. Raw Jet p_{T}), |y| bin "+std::to_string(k)+";#Delta Jet p_{T} [GeV];Ratio").c_str(), 
							     2500,0,2500);
      hjtpt_L2L3Res_2_rawpt_ratio[k]          = new TH1D(  ("jtpt_L2L3Res_2_rawpt_ratio_etabin"+std::to_string(k)).c_str(), 
							     ("(HiForest Jet p_{T})/(75X_data_Run2_v13 JEC Corr. Raw Jet p_{T}), |y| bin "+std::to_string(k)+";#Delta Jet p_{T} [GeV];Ratio").c_str(), 
							     2500,0,2500);      
      hL2L3Res_rawpt_L2L3Res_2_rawpt_ratio[k] = new TH1D(  ("L2L3Res_rawpt_L2L3Res_2_rawpt_ratio_etabin"+std::to_string(k)).c_str(), 
							     ("(74X_dataRun2_HLT_ppAt5TeV_v0 JEC Raw Jet p_{T})/(75X_data_Run2_v13 JEC Corr. Raw Jet p_{T}), |y| bin "+std::to_string(k)+";#Delta Jet p_{T} [GeV];Ratio").c_str(), 
							     2500,0,2500);
            
      hjtpt_v_diff_L2L3Res_rawpt[k]              = new TH2D( ("jtpt_v_diff_L2L3Res_rawpt_etabin"+std::to_string(k)).c_str(), 
							     ("(HiForest Jet p_{T}) - (74X_dataRun2_HLT_ppAt5TeV_v0 JEC Raw Jet p_{T}), |y| bin "+std::to_string(k)+";HiForest Jet p_{T} [GeV];#Delta Jet p_{T} [GeV]").c_str(), 
							     2500,0,2500, 800,-200,200);      
      hjtpt_v_diff_L2L3Res_2_rawpt[k]          = new TH2D( ("jtpt_v_diff_L2L3Res_2_rawpt_etabin"+std::to_string(k)).c_str(), 
							     ("(HiForest Jet p_{T}) - (75X_data_Run2_v13 JEC Corr. Raw Jet p_{T}), |y| bin "+std::to_string(k)+";HiForest Jet p_{T} [GeV];#Delta Jet p_{T} [GeV]").c_str(), 
							     2500,0,2500, 800,-200,200);      
      hL2L3Res_rawpt_v_diff_L2L3Res_2_rawpt[k] = new TH2D( ("L2L3Res_rawpt_v_diff_L2L3Res_2_rawpt_etabin"+std::to_string(k)).c_str(), 
							     ("(74X_dataRun2_HLT_ppAt5TeV_v0 JEC Corr. Raw Jet p_{T}) - (75X_data_Run2_v13 JEC Corr. Raw Jet p_{T}), |y| bin "+std::to_string(k)+";74X_dataRun2_HLT_ppAt5TeV_v0 JEC Corr. Raw Jet p_{T} [GeV];#Delta Jet p_{T} [GeV]").c_str(), 
							     2500,0,2500, 800,-200,200);      
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
  
  //  //temp check JEC Unc Hists
  //  std::cout<<std::endl;
  //  std::cout<<"n_jecunc_etabins="<<n_jecunc_etabins<<std::endl;
  //  for(int i=0; i<n_jecunc_etabins; i++){
  //    std::cout<<"jecunc_etabins["<<i<<"]="<<jecunc_etabins[i]<<std::endl;
  //  }
  //  std::cout<<std::endl;
  //  std::cout<<"n_jecunc_ptbins="<<n_jecunc_ptbins<<std::endl;
  //  for(int i=0; i<n_jecunc_ptbins; i++){
  //    std::cout<<"jecunc_ptbins["<<i<<"]="<<jecunc_ptbins[i]<<std::endl;
  //  }
  //  std::cout<<std::endl;
  //  for(int i=0; i<n_jecunc_etabins;i++){
  //    for(int j=0; j<n_jecunc_ptbins;j++){
  //      if(j!=n_jecunc_ptbins-1)std::cout<<"jecunc["<<i<<"]["<<j<<"]="<<jecunc[i][j]<<", ";
  //      else                    std::cout<<"jecunc["<<i<<"]["<<j<<"]="<<jecunc[i][j]<<"////////////////////"<<std::endl<<std::endl;
  //    }
  //  }
  //  assert(false);
  //temp

  //speed up bin searching by starting at spot in array that reflects the cuts being used to run
  int jecunc_eta_start=-1;
  if(fillDataJetJECUncHists)  for(int i=0; i<n_jecunc_etabins; i++){
      if( (-1.*jtEtaCutHi) >= jecunc_etabins[i] && (-1.*jtEtaCutHi) < jecunc_etabins[i+1]){
	jecunc_eta_start=i;       std::cout<<"jecunc_etabins["<<jecunc_eta_start<<"]="<<jecunc_etabins[i]<<std::endl;
	break;}
      else continue;
    }
  int jecunc_eta_end=-1;
  if(fillDataJetJECUncHists)  for(int i=0; i<n_jecunc_etabins; i++){
      if( (jtEtaCutHi) >= jecunc_etabins[i] && (jtEtaCutHi) < jecunc_etabins[i+1]){ 
	jecunc_eta_end=i;  std::cout<<"jecunc_etabins["<<jecunc_eta_end<<"]="<<jecunc_etabins[i]<<std::endl;     
	break;}
      else continue;
    }
  int jecunc_pt_start=-1;
  if(fillDataJetJECUncHists)  for(int i=0; i<n_jecunc_ptbins; i++){
      if( jtPtCut >= jecunc_ptbins[i] && jtPtCut < jecunc_ptbins[i+1]){
	jecunc_pt_start=i; std::cout<<"jecunc_ptbins["<<jecunc_pt_start<<"]="<<jecunc_ptbins[i]<<std::endl;      
	break;}
      else continue;
    }
  int jecunc_pt_end=-1;
  if(fillDataJetJECUncHists)  for(int i=0; i<n_jecunc_ptbins; i++){
      if( (jtPtCut_Hi) >= jecunc_ptbins[i] && (jtPtCut_Hi) < jecunc_ptbins[i+1]){
	jecunc_pt_end=i;      std::cout<<"jecunc_ptbins["<<jecunc_pt_end<<"]="<<jecunc_ptbins[i]<<std::endl;
	break;}
      else continue;
  }
  //assert(false);

  
  
  
  //jet QA, i.e. hists filled once per jet
  TH1D *hJetQA[nbins_abseta][N_vars][2]={};    
  TH1D *hJetQA_jtpt_JEC_sysup[nbins_abseta]={};    
  TH1D *hJetQA_jtpt_JEC_sysdown[nbins_abseta]={};    
  if(fillDataJetQAHists){        
    for(int j = 0; j<N_vars; ++j){		      
      for(int k = 0; k<nbins_abseta; ++k){      
	for(int i=0; i<2; i++){
	  //jets
	  if(var[j]=="jtpt"||var[j]=="rawpt"){
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;", var[j].c_str()) , 2500,0,2500);       
	    if(i==1 && fillDataJetJECUncHists && var[j]=="jtpt"){
	      hJetQA_jtpt_JEC_sysup[k]   = new TH1D( Form("hJetQA_%dwJetID_%s_JEC_sysup_etabin%d"  , i,var[j].c_str(),k) , Form("JEC sysup Jet p_{T};%s;"  , var[j].c_str()) , 2500,0,2500);       
	      hJetQA_jtpt_JEC_sysdown[k] = new TH1D( Form("hJetQA_%dwJetID_%s_JEC_sysdown_etabin%d", i,var[j].c_str(),k) , Form("JEC sysdown Jet p_{T};%s;", var[j].c_str()) , 2500,0,2500);       
	    }
	  }
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
  }//fill datajetqa hists
  
  
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
  
  //// trigger object histograms
  
  // //these get filled 
  // // // post cuts and post skips 
  // // // if the given HLT path fires AND max trig pt across all fired HLT bits meets the given LOWER threshold 
  // // // with the max trig pt across all fired triggger objects
  // // // weighted by the given HLT objects prescale
  TH1D *hpp_IncHLTMBtrgPt=NULL;
  TH1D  *hpp_IncHLT40trgPt =NULL ,   *hpp_IncHLT40trgEta =NULL;
  TH1D  *hpp_IncHLT60trgPt =NULL ,   *hpp_IncHLT60trgEta =NULL;
  TH1D  *hpp_IncHLT80trgPt =NULL ,   *hpp_IncHLT80trgEta =NULL;
  TH1D  *hpp_IncHLT100trgPt =NULL ,  *hpp_IncHLT100trgEta =NULL;  
  
  // //these get filled 
  // // // post cuts and post skips 
  // // // POST EXCLUSION CRITERIA
  // // // with the max trig pt across all fired HLT bits
  // // // weighted by the given HLT objects prescale
  TH1D *hpp_HLTMBtrgPt=NULL;
  TH1D  *hpp_HLT40trgPt =NULL ,   *hpp_HLT40trgEta =NULL;
  TH1D  *hpp_HLT60trgPt =NULL ,   *hpp_HLT60trgEta =NULL;
  TH1D  *hpp_HLT80trgPt =NULL ,   *hpp_HLT80trgEta =NULL;
  TH1D *hpp_HLT100trgPt =NULL ,  *hpp_HLT100trgEta =NULL;  
  TH1D *hpp_HLTCombtrgPt=NULL, *hpp_HLTCombtrgEta=NULL;

  
  if(fillDataJetTrigQAHists){    
    
    //semi excluded, looks at lower threshold but only respective trigger path's decision + pscl
    hpp_IncHLTMBtrgPt   = new TH1D( "IncHLTMB_trgPt"  , "inc trgPt for HLTMB   ", 2500, 0, 2500 ); 
    hpp_IncHLT40trgPt   = new TH1D( "IncHLT40_trgPt"  , "inc trgPt for HLT40   ", 2500, 0, 2500 ); 
    hpp_IncHLT60trgPt   = new TH1D( "IncHLT60_trgPt"  , "inc trgPt for HLT60   ", 2500, 0, 2500 ); 
    hpp_IncHLT80trgPt   = new TH1D( "IncHLT80_trgPt"  , "inc trgPt for HLT80   ", 2500, 0, 2500 ); 
    hpp_IncHLT100trgPt  = new TH1D( "IncHLT100_trgPt" , "inc trgPt for HLT100  ", 2500, 0, 2500 );       
    
    //full exclusion between triggers
    hpp_HLTMBtrgPt   = new TH1D( "isHLTMB_trgPt"  , "exc trgPt for HLTMB   ", 2500, 0, 2500 ); 
    hpp_HLT40trgPt   = new TH1D( "isHLT40_trgPt"  , "exc trgPt for HLT40   ", 2500, 0, 2500 ); 
    hpp_HLT60trgPt   = new TH1D( "isHLT60_trgPt"  , "exc trgPt for HLT60   ", 2500, 0, 2500 ); 
    hpp_HLT80trgPt   = new TH1D( "isHLT80_trgPt"  , "exc trgPt for HLT80   ", 2500, 0, 2500 ); 
    hpp_HLT100trgPt  = new TH1D( "isHLT100_trgPt" , "exc trgPt for HLT100  ", 2500, 0, 2500 );       
    
    //combination after full exclusion criteria applied
    hpp_HLTCombtrgPt = new TH1D( "HLTComb_trgPt"  , "trgPt for HLTComb ", 2500, 0, 2500 );       


    //semi excluded, looks at lower threshold but only respective trigger path's decision + pscl
    hpp_HLT40trgEta   = new TH1D( "isHLT40_trgEta"  , "exc trgEta for HLT40   ", 100, -5., 5. ); 
    hpp_HLT60trgEta   = new TH1D( "isHLT60_trgEta"  , "exc trgEta for HLT60   ", 100, -5., 5. ); 
    hpp_HLT80trgEta   = new TH1D( "isHLT80_trgEta"  , "exc trgEta for HLT80   ", 100, -5., 5. ); 
    hpp_HLT100trgEta  = new TH1D( "isHLT100_trgEta" , "exc trgEta for HLT100  ", 100, -5., 5. );       

    //full exclusion between triggers
    hpp_IncHLT40trgEta   = new TH1D( "IncHLT40_trgEta"  , "inc trgEta for HLT40   ", 100, -5., 5. ); 
    hpp_IncHLT60trgEta   = new TH1D( "IncHLT60_trgEta"  , "inc trgEta for HLT60   ", 100, -5., 5. ); 
    hpp_IncHLT80trgEta   = new TH1D( "IncHLT80_trgEta"  , "inc trgEta for HLT80   ", 100, -5., 5. ); 
    hpp_IncHLT100trgEta  = new TH1D( "IncHLT100_trgEta" , "inc trgEta for HLT100  ", 100, -5., 5. );       

    //combined after full exclusion criteria applied
    hpp_HLTCombtrgEta = new TH1D( "HLTComb_trgEta"  , "trgEta for HLTComb ", 100, -5., 5. );       

  }

  //special trigpt hists
  
  // // these get fillied
  // // // post skips but pre-cuts
  // // // if the given HLT path fired 
  // // // with the given trig pt('s)
  // // // weighted by the given HLT objects prescale
  TH1D   *hpp_IncHLTMBpscl =NULL ;
  TH1D   *hpp_IncHLT40pscl =NULL ;
  TH1D   *hpp_IncHLT60pscl =NULL ;
  TH1D   *hpp_IncHLT80pscl =NULL ;
  TH1D  *hpp_IncHLT100pscl =NULL ;  

  TH1D   *hpp_IncHLT40trgPt_allobj =NULL ;
  TH1D   *hpp_IncHLT60trgPt_allobj =NULL ;
  TH1D   *hpp_IncHLT80trgPt_allobj =NULL ;
  TH1D  *hpp_IncHLT100trgPt_allobj =NULL ;  

  TH1D   *hpp_IncHLT40trgPt_leadobj =NULL ;
  TH1D   *hpp_IncHLT60trgPt_leadobj =NULL ;
  TH1D   *hpp_IncHLT80trgPt_leadobj =NULL ;
  TH1D  *hpp_IncHLT100trgPt_leadobj =NULL;
  TH1D   *hpp_IncHLT40MBtrgPt_leadobj =NULL ;
  
  TH2D   *hpp_IncHLT40trgPt_leadobj_leadjet =NULL ;
  TH2D   *hpp_IncHLT60trgPt_leadobj_leadjet =NULL ;
  TH2D   *hpp_IncHLT80trgPt_leadobj_leadjet =NULL ;
  TH2D  *hpp_IncHLT100trgPt_leadobj_leadjet =NULL ;  
  TH2D   *hpp_IncHLT40MBtrgPt_leadobj_leadjet =NULL ;
  
  TH2D   *hpp_IncHLT40trgPt_leadobj_globleadobj =NULL ;
  TH2D   *hpp_IncHLT60trgPt_leadobj_globleadobj =NULL ;
  TH2D   *hpp_IncHLT80trgPt_leadobj_globleadobj =NULL ;
  TH2D  *hpp_IncHLT100trgPt_leadobj_globleadobj =NULL ;  
  TH2D   *hpp_IncHLT40MBtrgPt_leadobj_globleadobj =NULL ;  
  
  if(fillDataJetTrigQAHists){

    hpp_IncHLTMBpscl   =new TH1D( "IncHLTMBpscl",    "HLTMB pscl;Prescale;N_{evts}" , 1000,0,1000);
    hpp_IncHLT40pscl   =new TH1D( "IncHLT40pscl",    "HLT40 pscl;Prescale;N_{evts}" , 1000,0,1000);
    hpp_IncHLT60pscl   =new TH1D( "IncHLT60pscl",    "HLT60 pscl;Prescale;N_{evts}" , 1000,0,1000);
    hpp_IncHLT80pscl   =new TH1D( "IncHLT80pscl",    "HLT80 pscl;Prescale;N_{evts}" , 1000,0,1000);
    hpp_IncHLT100pscl  =new TH1D( "IncHLT100pscl",  "HLT100 pscl;Prescale;N_{evts}", 1000,0,1000);    
    
    hpp_IncHLT40trgPt_allobj   =new TH1D( "IncHLT40trgPt_allobj",    "all trigger jet object for HLT40;Trigger Jet Obj p_{T};N_{jets}" , 2500,0,2500);
    hpp_IncHLT60trgPt_allobj   =new TH1D( "IncHLT60trgPt_allobj",    "all trigger jet object for HLT60;Trigger Jet Obj p_{T};N_{jets}" , 2500,0,2500);
    hpp_IncHLT80trgPt_allobj   =new TH1D( "IncHLT80trgPt_allobj",    "all trigger jet object for HLT80;Trigger Jet Obj p_{T};N_{jets}" , 2500,0,2500);
    hpp_IncHLT100trgPt_allobj  =new TH1D( "IncHLT100trgPt_allobj",  "all trigger jet object for HLT100;Trigger Jet Obj p_{T};N_{jets}", 2500,0,2500);    

    hpp_IncHLT40trgPt_leadobj   =new TH1D(  "IncHLT40trgPt_leadobj",  "lead trigger jet object for HLT40;Lead Trigger Jet Obj p_{T};N_{jets}" , 2500,0,2500);
    hpp_IncHLT60trgPt_leadobj   =new TH1D(  "IncHLT60trgPt_leadobj",  "lead trigger jet object for HLT60;Lead Trigger Jet Obj p_{T};N_{jets}" , 2500,0,2500);
    hpp_IncHLT80trgPt_leadobj   =new TH1D(  "IncHLT80trgPt_leadobj",  "lead trigger jet object for HLT80;Lead Trigger Jet Obj p_{T};N_{jets}" , 2500,0,2500);
    hpp_IncHLT100trgPt_leadobj  =new TH1D(  "IncHLT100trgPt_leadobj",  "lead trigger jet object for HLT100;Lead Trigger Jet Obj p_{T};N_{jets}", 2500,0,2500);
    hpp_IncHLT40MBtrgPt_leadobj =new TH1D(  "IncHLT40MBtrgPt_leadobj", "lead trigger jet obj for HLT40 && HLTMB" , 2500,0,2500);
    
    hpp_IncHLT40trgPt_leadobj_leadjet  =new TH2D(  "IncHLT40trgPt_leadobj_leadjet", "lead trigger jet obj v lead RECO jet obj, only for HLT40;Lead Trigger Obj p_{T};Lead RECO Obj p_{T}" , 
						   2500,0,2500, 2500,0,2500);
    hpp_IncHLT60trgPt_leadobj_leadjet  =new TH2D(  "IncHLT60trgPt_leadobj_leadjet", "lead trigger jet obj v lead RECO jet obj, only for HLT60;Lead Trigger Obj p_{T};Lead RECO Obj p_{T}" , 
						   2500,0,2500, 2500,0,2500);
    hpp_IncHLT80trgPt_leadobj_leadjet  =new TH2D(  "IncHLT80trgPt_leadobj_leadjet", "lead trigger jet obj v lead RECO jet obj, only for HLT80;Lead Trigger Obj p_{T};Lead RECO Obj p_{T}" , 
						   2500,0,2500, 2500,0,2500);
    hpp_IncHLT100trgPt_leadobj_leadjet =new TH2D(  "IncHLT100trgPt_leadobj_leadjet", "lead trigger jet obj v lead RECO jet obj, only for HLT100;Lead Trigger Obj p_{T};Lead RECO Obj p_{T}", 
						   2500,0,2500, 2500,0,2500);
    hpp_IncHLT40MBtrgPt_leadobj_leadjet  =new TH2D(  "IncHLT40MBtrgPt_leadobj_leadjet", "lead trigger jet obj v lead RECO jet obj, only for HLT40 && HLTMB" ,
						     2500,0,2500, 2500,0,2500);    
    
    
    hpp_IncHLT40trgPt_leadobj_globleadobj  =new TH2D(  "IncHLT40trgPt_leadobj_globleadobj", "lead trigger jet obj across HLT40 v lead trigger jet obj across all triggers, only for HLT40;Lead Trigger Obj p_{T} across HLT40;Lead Obj p_{T} across HLT40+60+80" , 
						       2500,0,2500, 2500,0,2500);
    hpp_IncHLT60trgPt_leadobj_globleadobj  =new TH2D(  "IncHLT60trgPt_leadobj_globleadobj", "lead trigger jet obj across HLT60 v lead trigger jet obj across all triggers, only for HLT60;Lead Trigger Obj p_{T} across HLT60;Lead Obj p_{T} across HLT40+60+80" , 
						       2500,0,2500, 2500,0,2500);
    hpp_IncHLT80trgPt_leadobj_globleadobj  =new TH2D(  "IncHLT80trgPt_leadobj_globleadobj", "lead trigger jet obj across HLT80 v lead trigger jet obj across all triggers, only for HLT80;Lead Trigger Obj p_{T} across HLT80;Lead Obj p_{T} across HLT40+60+80" , 
						       2500,0,2500, 2500,0,2500);
    hpp_IncHLT100trgPt_leadobj_globleadobj =new TH2D(  "IncHLT100trgPt_leadobj_globleadobj", "lead trigger jet obj across HLT100 v lead trigger jet obj across all triggers, only for HLT100;Lead Trigger Obj p_{T} across HLT100;Lead Obj p_{T} across HLT40+60+80", 
						       2500,0,2500, 2500,0,2500); 
    hpp_IncHLT40MBtrgPt_leadobj_globleadobj  =new TH2D(  "IncHLT40MBtrgPt_leadobj_globleadobj", "lead trigger jet obj across HLT40 v lead trigger jet obj across all triggers, only for HLT40+HLTMB;Lead Trigger Obj p_{T} across HLT40+HLTMB;Lead Obj p_{T} across HLT40+60+80" , 
						       2500,0,2500, 2500,0,2500);   


  }






  //inclusive jet spectra by-trigger histograms
  TH1D  *hpp_HLT40[nbins_abseta][2]={}, *hpp_HLT60[nbins_abseta][2]={}, *hpp_HLT80[nbins_abseta][2]={}, *hpp_HLT100[nbins_abseta][2]={}; //no exclusion, looks at each path on its own terms
  TH1D *hpp_HLTMB[nbins_abseta][2]={};
  
  TH1D  *hpp_incHLT40[nbins_abseta][2]={}, *hpp_incHLT60[nbins_abseta][2]={}, *hpp_incHLT80[nbins_abseta][2]={}, *hpp_incHLT100[nbins_abseta][2]={};// semi excluded, using max trig pt between all objects and looks at lower threshold
  TH1D *hpp_incHLTMB[nbins_abseta][2]={};
  
  TH1D  *hpp_excHLT40[nbins_abseta][2]={}, *hpp_excHLT60[nbins_abseta][2]={}, *hpp_excHLT80[nbins_abseta][2]={}, *hpp_excHLT100[nbins_abseta][2]={};// full exclusion, using max trig pt between all objects and looking at all thresholds
  TH1D* hpp_excHLTMB[nbins_abseta][2]={};  
  TH1D *hpp_HLTComb[nbins_abseta][2]={};//combination of triggers w/ exclusion    

  //special histograms

  //|y|<0.5 only, wJetID only
  TH1D* hpp_excHLTMB_NOT_l1jet=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet8=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet12=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet16=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet20=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet24=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet28=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet40=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet48=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet52=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet60=NULL;
  TH1D* hpp_excHLTMB_AND_l1jet68=NULL;

  TH1D* hpp_excHLTMB_NOT_HLTjet=NULL;
  TH1D* hpp_excHLTMB_AND_HLTjet40=NULL;
  TH1D* hpp_excHLTMB_AND_HLTjet60=NULL;
  TH1D* hpp_excHLTMB_AND_HLTjet80=NULL;
  //TH1D* hpp_excHLTMB_AND_HLTjet100=NULL;

  TH1D* hpp_excHLTMB_NOT_HFcoincidence=NULL;

  TH1D* hpp_excHLT40_NOT_HLTMB=NULL;
  TH1D* hpp_excHLT60_NOT_HLTMB=NULL;
  TH1D* hpp_excHLT80_NOT_HLTMB=NULL;
  TH1D* hpp_incHLT40_NOT_HLTMB=NULL;
  TH1D* hpp_incHLT60_NOT_HLTMB=NULL;
  TH1D* hpp_incHLT80_NOT_HLTMB=NULL;
  
  if(fillDataJetTrigQAHists)    {
    
    //no thresholds check, no exclusion checked
    std::string hJetHLTQATitleArray[]={ "hpp_HLTMB","hpp_HLT40", "hpp_HLT60",  "hpp_HLT80", "hpp_HLT100"};
    std::string hJetHLTQADescArray[]={  "HLTMBJet","HLT40Jet","HLT60Jet",  "HLT80Jet", "HLT100Jet"};
    const int nJetHLTQATitles=sizeof(hJetHLTQATitleArray)/sizeof(std::string);
    
    //meets lower threshold only, no exclusion checked. combination in this array simply by historical reasons
    std::string hJetIncTrigQATitleArray[]={ "hpp_incHLTMB","hpp_incHLT40", "hpp_incHLT60",  "hpp_incHLT80", "hpp_incHLT100", "hpp_HLTComb"};
    std::string hJetIncTrigQADescArray[]={  "HLTMBJet+Threshold","HLT40Jet+Threshold","HLT60Jet+Threshold",  "HLT80Jet+Threshold", "HLT100Jet+Threshold", "HLTCombo"};
    const int nJetIncTrigQATitles=sizeof(hJetIncTrigQATitleArray)/sizeof(std::string);
    
    //meets lower+upper thresholds, exclusion applied
    std::string hJetExcTrigQATitleArray[]={ "hpp_excHLTMB","hpp_excHLT40", "hpp_excHLT60",  "hpp_excHLT80", "hpp_excHLT100"};
    std::string hJetExcTrigQADescArray[]={ "ExcHLTMBJet","ExcHLT40Jet","ExcHLT60Jet",  "ExcHLT80Jet", "ExcHLT100Jet"};
    const int nJetExcTrigQATitles=sizeof(hJetExcTrigQATitleArray)/sizeof(std::string);
    
    //this loop inits hists for the incl. jet *spectra* under each trigger path decision

    for(int k=0; k<nJetIncTrigQATitles; k++){	      //if(debugMode) std::cout<<"k=="<<k<<std::endl;	         
      for(int j=0;j<nbins_abseta;j++){    	//if(debugMode) std::cout<<"j=="<<j<<std::endl;	  	  
	for(int i=0; i<2; i++){
	  
	  //incl trigg plots
	  std::string hDesc=hJetIncTrigQADescArray[k];
	  if(i==1)hDesc+=" w/ JetID";
	  hDesc+= ", "+absetabins_str[j]+ "< #||{y} <"+absetabins_str[j+1];
	  //if(debugMode)std::cout<<"inclusive HLT plot hDesc="<<hDesc<<std::endl;
	  
	  std::string hTitle=hJetIncTrigQATitleArray[k];
	  if(i==1)hTitle+="_wJetID";
	  hTitle+="_R"+std::to_string(radius)+"_etabin"+std::to_string(j);//etaWidth;	
	  //if(debugMode)std::cout<<"inclusive HLT plot hTitle="<<hTitle<<std::endl;
	  
	  if     (hJetIncTrigQATitleArray[k]=="hpp_incHLTMB"  )  hpp_incHLTMB[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetIncTrigQATitleArray[k]=="hpp_incHLT40"  )  hpp_incHLT40[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetIncTrigQATitleArray[k]=="hpp_incHLT60"  )  hpp_incHLT60[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetIncTrigQATitleArray[k]=="hpp_incHLT80"  )  hpp_incHLT80[j][i]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	  else if(hJetIncTrigQATitleArray[k]=="hpp_incHLT100" )  hpp_incHLT100[j][i]    = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);	
	  else if(hJetIncTrigQATitleArray[k]=="hpp_HLTComb")     hpp_HLTComb[j][i]      = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	

	  if(k<nJetHLTQATitles){
	    std::string hHLTTitle=hJetHLTQATitleArray[k];
	    //if(fillDataJetIDHists)hExcTitle+="_wJetID";
	    if(i==1)hHLTTitle+="_wJetID";
	    hHLTTitle+="_R"+std::to_string(radius)+"_etabin"+std::to_string(j);//etaWidth;
	    //if(debugMode)std::cout<<" exclusive HLT plothExcTitle="<<hExcTitle<<std::endl;	  
	    
	    std::string hHLTDesc=hJetHLTQADescArray[k];
	    //if(fillDataJetIDHists)hExcDesc+=" w/ JetID";
	    if(i==1)hHLTDesc+=" w/ JetID";
	    hHLTDesc+=", "+absetabins_str[j] +"< #||{y} <"+absetabins_str[j+1];	    
	    //if(debugMode)std::cout<<" exclusive HLT plothExcDesc="<<hExcDesc<<std::endl<<std::endl;	    
	    
	    if     (hJetHLTQATitleArray[k]=="hpp_HLTMB"  )  hpp_HLTMB[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT40"  )  hpp_HLT40[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT60"  )  hpp_HLT60[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT80"  )  hpp_HLT80[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT100" )  hpp_HLT100[j][i]    = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);	
	  }

	  

	  //excl trigg plots	
	  if(k<nJetExcTrigQATitles){

	    std::string hExcTitle=hJetExcTrigQATitleArray[k];
	    //if(fillDataJetIDHists)hExcTitle+="_wJetID";
	    if(i==1)hExcTitle+="_wJetID";
	    hExcTitle+="_R"+std::to_string(radius)+"_etabin"+std::to_string(j);//etaWidth;
	    //if(debugMode)std::cout<<" exclusive HLT plothExcTitle="<<hExcTitle<<std::endl;	  
	    
	    std::string hExcDesc=hJetExcTrigQADescArray[k];
	    //if(fillDataJetIDHists)hExcDesc+=" w/ JetID";
	    if(i==1)hExcDesc+=" w/ JetID";
	    hExcDesc+=", "+absetabins_str[j] +"< #||{y} <"+absetabins_str[j+1];	    
	    //if(debugMode)std::cout<<" exclusive HLT plothExcDesc="<<hExcDesc<<std::endl<<std::endl;	    
	    
	    if     (hJetExcTrigQATitleArray[k]=="hpp_excHLTMB"  )  {
	      //std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLTMB[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if     (hJetExcTrigQATitleArray[k]=="hpp_excHLT40"  )  {
	      //std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT40[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT60"  )  {
	      //std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT60[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT80"  )  {
	      //std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT80[j][i]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	    else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT100" )  {
	      //std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	      hpp_excHLT100[j][i]    = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}	  
	  }
	  
	}//jet id loop
      }//jettrig qa titles 
    }//absetabins

    
    hpp_excHLTMB_NOT_l1jet  =new TH1D(  "excHLTMB_NOT_l1jet"  ,  "isMB NOT l1jet  , #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet8 =new TH1D(  "excHLTMB_AND_l1jet8" ,  "isMB AND l1jet8 , #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet12=new TH1D(  "excHLTMB_AND_l1jet12",  "isMB AND l1jet12, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet16=new TH1D(  "excHLTMB_AND_l1jet16",  "isMB AND l1jet16, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet20=new TH1D(  "excHLTMB_AND_l1jet20",  "isMB AND l1jet20, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet24=new TH1D(  "excHLTMB_AND_l1jet24",  "isMB AND l1jet24, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet28=new TH1D(  "excHLTMB_AND_l1jet28",  "isMB AND l1jet28, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet40=new TH1D(  "excHLTMB_AND_l1jet40",  "isMB AND l1jet40, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet48=new TH1D(  "excHLTMB_AND_l1jet48",  "isMB AND l1jet48, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet52=new TH1D(  "excHLTMB_AND_l1jet52",  "isMB AND l1jet52, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet60=new TH1D(  "excHLTMB_AND_l1jet60",  "isMB AND l1jet60, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_l1jet68=new TH1D(  "excHLTMB_AND_l1jet68",  "isMB AND l1jet68, #||{y} < 0.5, w JetID" , 2500, 0, 2500);

    hpp_excHLTMB_NOT_HLTjet  =new TH1D(  "excHLTMB_NOT_HLTjet",    "isMB NOT   HLTjet, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_HLTjet40=new TH1D(  "excHLTMB_AND_HLTjet40",  "isMB AND HLTjet40, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_HLTjet60=new TH1D(  "excHLTMB_AND_HLTjet60",  "isMB AND HLTjet60, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLTMB_AND_HLTjet80=new TH1D(  "excHLTMB_AND_HLTjet80",  "isMB AND HLTjet80, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    
    hpp_excHLTMB_NOT_HFcoincidence=new TH1D(  "excHLTMB_NOT_HFcoincidence",  "isMB NOT HF coincidence, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    
    hpp_excHLT40_NOT_HLTMB=new TH1D(  "excHLT40_NOT_HLTMB",  "is40 NOT HLTMB, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLT60_NOT_HLTMB=new TH1D(  "excHLT60_NOT_HLTMB",  "is60 NOT HLTMB, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_excHLT80_NOT_HLTMB=new TH1D(  "excHLT80_NOT_HLTMB",  "is80 NOT HLTMB, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    
    hpp_incHLT40_NOT_HLTMB=new TH1D(  "incHLT40_NOT_HLTMB",  "HLT40 NOT HLTMB, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_incHLT60_NOT_HLTMB=new TH1D(  "incHLT60_NOT_HLTMB",  "HLT60 NOT HLTMB, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
    hpp_incHLT80_NOT_HLTMB=new TH1D(  "incHLT80_NOT_HLTMB",  "HLT60 NOT HLTMB, #||{y} < 0.5, w JetID" , 2500, 0, 2500);
  }//end fill datajettrigqa plots
  
  
  //TH2D* hpp_covmat=NULL;
  TH2D* hpp_covmat_eta_arr[nbins_abseta]={};
  TH1D* hpp_covmat_eta_arr_helpers[nbins_abseta]={};
  bool etabin_bool_helper[nbins_abseta]={0};
  
  TH2D* hpp_covmat_eta_arr_JEC_sysup[nbins_abseta]={};
  TH1D* hpp_covmat_eta_arr_helpers_JEC_sysup[nbins_abseta]={};
  TH2D* hpp_covmat_eta_arr_JEC_sysdown[nbins_abseta]={};
  TH1D* hpp_covmat_eta_arr_helpers_JEC_sysdown[nbins_abseta]={};
  bool JEC_etabin_bool_helper[nbins_abseta]={0};
  
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

      if(fillDataJetJECUncHists){
	std::string covmat_eta_name_JEC_sysup="hpp_covmat_";
	//if(fillDataJetIDHists)covmat_eta_name+="wJetID_";
	covmat_eta_name_JEC_sysup+="wJetID_JEC_sysup_";
	covmat_eta_name_JEC_sysup+="etabin_"+std::to_string(i);      
	if(debugMode)std::cout<<"covmat_eta_name_JEC_sysup  =" << covmat_eta_name_JEC_sysup << std::endl;      
	hpp_covmat_eta_arr_JEC_sysup[i]=new TH2D((covmat_eta_name_JEC_sysup).c_str(), ("JEC sysup "+covmat_eta_title).c_str(), nbins_pt, ptbins, nbins_pt, ptbins);     
	hpp_covmat_eta_arr_helpers_JEC_sysup[i]=new TH1D(("covmat_helper_JEC_sysup_etabin"+std::to_string(i)).c_str(), "helper only", nbins_pt, ptbins);      

	std::string covmat_eta_name_JEC_sysdown="hpp_covmat_";
	//if(fillDataJetIDHists)covmat_eta_name+="wJetID_";
	covmat_eta_name_JEC_sysdown+="wJetID_JEC_sysdown_";
	covmat_eta_name_JEC_sysdown+="etabin_"+std::to_string(i);      
	if(debugMode)std::cout<<"covmat_eta_name_JEC_sysdown  =" << covmat_eta_name_JEC_sysdown << std::endl;      
	hpp_covmat_eta_arr_JEC_sysdown[i]=new TH2D((covmat_eta_name_JEC_sysdown).c_str(), ("JEC sysdown "+covmat_eta_title).c_str(), nbins_pt, ptbins, nbins_pt, ptbins);     
	hpp_covmat_eta_arr_helpers_JEC_sysdown[i]=new TH1D(("covmat_helper_JEC_sysdown_etabin"+std::to_string(i)).c_str(), "helper only", nbins_pt, ptbins);      


      }
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
  
  //neutral hadrons from PF (no HF e/gamma contributions)
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
    //if(debugMode){std::cout<<std::endl;
    //  std::cout<<"L1Branches      [i="<<i<<"]= "<<L1Branches[i]<<std::endl;
    //  std::cout<<"L1PresclBranches      [i="<<i<<"]= "<<L1PresclBranches[i]<<std::endl;    }    
  }




  //L1 seed used for HLT seeding HLTMB
  int mb_l1s_I;
  int mb_l1s_ps_I;
  jetpp[4]->SetBranchAddress( "L1_MinimumBiasHF1_OR", &mb_l1s_I);
  jetpp[4]->SetBranchAddress( "L1_MinimumBiasHF1_OR_Prescl", &mb_l1s_ps_I);
  
  //technical bit to perhaps be used for cleaning up min bias signal
  int L1_HF_coincidence_I;
  jetpp[4]->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2",&L1_HF_coincidence_I);

  //other L1 Jets not used for seeding HLT40,60,80,100
  int l1jet24_l1s_I, l1jet20_l1s_I, l1jet16_l1s_I, l1jet12_l1s_I, l1jet8_l1s_I ;
  jetpp[4]->SetBranchAddress (  "L1_SingleJet8_BptxAND", & l1jet8_l1s_I  );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet12_BptxAND", & l1jet12_l1s_I );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet16_BptxAND", & l1jet16_l1s_I );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet20_BptxAND", & l1jet20_l1s_I );
  jetpp[4]->SetBranchAddress (  "L1_SingleJet24_BptxAND", & l1jet24_l1s_I );

  //L1 Jets used for HLT seeding HLT40,60,80,100
  int jet40_l1s_I, jet60_l1s_I,  jet80_l1s_I,  jet100_l1s_I;   //L1 
  int jet40_l1s_ps_I, jet60_l1s_ps_I,  jet80_l1s_ps_I,  jet100_l1s_ps_I;     
  jetpp[4]->SetBranchAddress( L1Branches[0].c_str()  , &jet40_l1s_I);   //l1jet28
  jetpp[4]->SetBranchAddress( L1Branches[1].c_str()  , &jet60_l1s_I);   //l1jet40
  jetpp[4]->SetBranchAddress( L1Branches[2].c_str()  , &jet80_l1s_I);   //l1jet48
  jetpp[4]->SetBranchAddress( L1Branches[3].c_str()  , &jet100_l1s_I);  //l1jet52
  
  jetpp[4]->SetBranchAddress( L1PresclBranches[0].c_str()  , &jet40_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[1].c_str()  , &jet60_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[2].c_str()  , &jet80_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[3].c_str()  , &jet100_l1s_ps_I);
  
  int l1jet60_l1s_I, l1jet68_l1s_I;
  jetpp[4]->SetBranchAddress( "L1_SingleJet60_BptxAND", &l1jet60_l1s_I);
  jetpp[4]->SetBranchAddress( "L1_SingleJet68_BptxAND", &l1jet68_l1s_I);




 
  
  
  // make sure # of L1/HLT bits is correct
  assert(N_HLTBits==N_L1Bits);  


  int MB_HF1ORp5_I, MB_HF1ORp5_p_I;
  jetpp[4]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1"        , &MB_HF1ORp5_I);
  jetpp[4]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1_Prescl" , &MB_HF1ORp5_p_I);
  
  
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
    //if(debugMode){
    //  std::cout<<std::endl;
    //  std::cout<<"trgCombType            =="<<trgCombType<<std::endl;
    //  std::cout<<"HLTthresh             [i="<<i<<"]= "<<HLTthresh[i]<<std::endl;
    //  std::cout<<std::endl;
    //  std::cout<<"PF_HLTBranches        [i="<<i<<"]= "<<PF_HLTBranches[i]<<std::endl;
    //  std::cout<<"PF_HLTPresclBranches  [i="<<i<<"]= "<<PF_HLTPresclBranches[i]<<std::endl;
    //  std::cout<<"HLTPFthresh           [i="<<i<<"]= "<<HLTPFthresh[i]<<std::endl;
    //  std::cout<<std::endl;
    //  std::cout<<"Calo_HLTBranches      [i="<<i<<"]= "<<Calo_HLTBranches[i]<<std::endl;
    //  std::cout<<"Calo_HLTPresclBranches[i="<<i<<"]= "<<Calo_HLTPresclBranches[i]<<std::endl;
    //  std::cout<<"HLTCalothresh         [i="<<i<<"]= "<<HLTCalothresh[i]<<std::endl;//<<std::endl; 
    //  std::cout<<std::endl; }      
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
  
  // other bits that i need to keep track of for full proper duplicate skip over jet80 events in the min bias dataset
  // pretty sure everyone of these bits are redundant w/ jet80 bits i.e. if they fire, so did the jet80 bits.
  // but i'm having issues... so i suspect no one and i suspect everyone.
  int CaloJet80_45_45_Eta2p1_I,  CaloJet80_Jet35_Eta0p7_I,  CaloJet80_Jet35_Eta1p1_I;
  int CaloJet100_Jet35_Eta0p7_I,  CaloJet100_Jet35_Eta1p1_I;
  int CaloJet110_I,  CaloJet120_I,  CaloJet150_I;
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet80_45_45_Eta2p1_v1", &CaloJet80_45_45_Eta2p1_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta0p7_v1", &CaloJet80_Jet35_Eta0p7_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta1p1_v1", &CaloJet80_Jet35_Eta1p1_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta0p7_v1", &CaloJet100_Jet35_Eta0p7_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta1p1_v1", &CaloJet100_Jet35_Eta1p1_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet110_Eta5p1_v1", &CaloJet110_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet120_Eta5p1_v1", &CaloJet120_I );
  jetpp[4]->SetBranchAddress("HLT_AK4CaloJet150_v1", &CaloJet150_I );
  
  int PFJet120_I, PFJet110_I;
  jetpp[4]->SetBranchAddress("HLT_AK4PFJet110_Eta5p1_v1", &PFJet110_I);
  jetpp[4]->SetBranchAddress("HLT_AK4PFJet120_Eta5p1_v1", &PFJet120_I);
  


  
  

  //tupel/EventTree
  std::vector<float> *jetPt, *jetEta, *jetPhi, *jetE, *jetID, *jetRawPt, *jetRawE;
  std::vector<float> *jetHfHadE, *jetHfEmE, *jetChHadFrac, *jetNeutralHadAndHfFrac, *jetChEmFrac, *jetNeutralEmFrac, *jetChMult, *jetConstCnt;  
  if(useTupel){
    if(deepDebug)std::cout<<"setting branch addresses for tupel/EventTree"<<std::endl;

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
  }
  

  


  
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
  //if(debugMode)NEvents_toRead=1666666;
  //if(debugMode)NEvents_toRead=(int)(NEvents_toRead/10);
  std::cout<<"reading "<<NEvents_toRead<<" out of "<< NEvents_allFiles<<" events"<<std::endl;  
  //std::cout<<"out of "<<NEvents_allFiles<<" events"<<std::endl;  
  

  UInt_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;  //trigger event counts
  UInt_t NEvents_MB=0;//minbias


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
      L2RelL3AbsJEC_arr[i]->SetParameter(8, L2RelL3Abs_par8_arr[i]);    }        
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
      L2L3ResJEC_arr[i]->SetParameter(8, L2L3Res_par8_arr[i]);    }
  }

  //  //readForests_JEC_v12_ak4PF_75X_data_Run2_v13
  TF1* L2RelL3AbsJEC_2_arr[n_L2RelL3AbsJEC_2_etabins]={};
  TF1* L2L3ResJEC_2_arr[n_L2L3ResJEC_2_etabins]={};
  if(fillDataJetJECQAHists){    
    for(int i=0; i<n_L2RelL3AbsJEC_2_etabins; i++){
      L2RelL3AbsJEC_2_arr[i]=new TF1(("L2RelL3AbsJEC_2_etabin"+std::to_string(i)).c_str(), L2RelL3AbsJEC_2_str.c_str(),L2RelL3AbsJEC_2_ptlo,L2RelL3AbsJEC_2_pthi);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(0, L2RelL3Abs_2_par0_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(1, L2RelL3Abs_2_par1_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(2, L2RelL3Abs_2_par2_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(3, L2RelL3Abs_2_par3_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(4, L2RelL3Abs_2_par4_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(5, L2RelL3Abs_2_par5_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(6, L2RelL3Abs_2_par6_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(7, L2RelL3Abs_2_par7_arr[i]);
      L2RelL3AbsJEC_2_arr[i]->SetParameter(8, L2RelL3Abs_2_par8_arr[i]);    }        
    for(int i=0; i<n_L2L3ResJEC_2_etabins; i++){
      L2L3ResJEC_2_arr[i]=new TF1(("L2L3ResJEC_2_etabin"+std::to_string(i)).c_str(), L2L3ResJEC_2_str.c_str(),L2L3ResJEC_2_ptlo,L2L3ResJEC_2_pthi);
      L2L3ResJEC_2_arr[i]->SetParameter(0, L2L3Res_2_par0_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(1, L2L3Res_2_par1_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(2, L2L3Res_2_par2_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(3, L2L3Res_2_par3_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(4, L2L3Res_2_par4_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(5, L2L3Res_2_par5_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(6, L2L3Res_2_par6_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(7, L2L3Res_2_par7_arr[i]);
      L2L3ResJEC_2_arr[i]->SetParameter(8, L2L3Res_2_par8_arr[i]);    }
  }



  
  //begin event loop
  int NJets_largeJECdiff=0;
  //int assertcount=0;
  //float currentIntegral=0.;
  //int problemCount=0;
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {
//    if(nEvt!=152)continue;
//    if(nEvt>152)assert(false);
    
    if (nEvt%100000==0 || deepDebug) {
      //if (true) {
      float nEvt_fracdone=((float)nEvt)/((float)NEvents_toRead);
      float nEvt_percdone=nEvt_fracdone*100.;
      int   nEvt_intpercdone=(int)nEvt_percdone;                
      std::cout<<"from trees, grabbing Evt # = "<<nEvt<<", ~"<<nEvt_intpercdone<<"% of job completed"<<std::endl;
    }
    
    //grab an entry    
    jetpp[0]->GetEntry(nEvt);
    
    if(doRunExclStudy){
      bool inGoodRunList=false;
      for(int i=0; i<NGoodRuns; i++){
	if(run_I==runarray_forExclRunStudy[i]){
	  inGoodRunList=true; break;}
      }
      if(!inGoodRunList) continue;	
    }
    
    if(fillDataEvtQAHists){
      


      hRuns->Fill(run_I, 1.);//fill this w/ no pscl wight for now, pre-skip/QA cuts
      
      for(int i=0;i<Nruns;i++){
	if(i!=(Nruns-1)){
	  //std::cout<<"run_I="<<run_I<<std::endl;
	  //std::cout<<"runarray["<<i<<"]="<<runarray[i]<<std::endl;
	  //std::cout<<"hLumiByRun["<<i<<"]->GetTitle()="<<hLumiByRun[i]->GetTitle()<<std::endl;
	  if(run_I==runarray[i]){
	    hLumiByRun[i]->Fill(lumi_I,1.);
	    break;
	  }
	  else continue;
	}
	else{//in case the run is not in the run list from gold JSON
	  //std::cout<<"run_I="<<run_I<<std::endl;
	  //std::cout<<"runarray["<<i<<"]="<<runarray[i]<<std::endl;
	  //std::cout<<"hLumiByRun["<<i<<"]->GetTitle()="<<hLumiByRun[i]->GetTitle()<<std::endl;
	  std::cout<<"WARNING! Evt not in gold JSON run list! details are:"<<std::endl;
	  std::cout<<"run_I="<<run_I<<std::endl;
	  std::cout<<"lumi_I="<<lumi_I<<std::endl;
	  std::cout<<"evt_I="<<evt_I<<std::endl;
	  hLumiByRun[i]->Fill(lumi_I,1.);
	}
      }
    }
    

    bool MBtrgDec= ( (bool) ( MB_HF1ORp5_I * mb_l1s_I ) ) && useMBevts;
    int  MBtrgPscl= MB_HF1ORp5_p_I * mb_l1s_ps_I;    
    
    bool trgDec[N_HLTBits]={ (bool) (Jet40_I * jet40_l1s_I) , (bool) ( Jet60_I * jet60_l1s_I ) ,
			     (bool) (Jet80_I * jet80_l1s_I) , (bool) (Jet100_I * jet100_l1s_I) };    
    int trgPscl[N_HLTBits]={ Jet40_p_I * jet40_l1s_ps_I ,  Jet60_p_I *  jet60_l1s_ps_I ,   
                             Jet80_p_I * jet80_l1s_ps_I , Jet100_p_I * jet100_l1s_ps_I };    

    //    bool trgDec2[N_HLTBits]={ (bool) (Jet40_2_I* jet40_l1s_I ), (bool) ( Jet60_2_I * jet60_l1s_I ) , //keep track of the other trigger decision bits
    //			      (bool) (Jet80_2_I* jet80_l1s_I ), (bool) (Jet100_2_I * jet100_l1s_I) };    
    
    
    //bool isInMinBiasPD=( (bool) mb_l1s_I );
    bool isInMinBiasPD=( (bool) MB_HF1ORp5_I );
    if(filelistIsMinBias)isInMinBiasPD=true;
    //bool isInMinBiasPD=( (bool) MB_HF1ORp5_I ||  (bool) mb_l1s_I);
    //bool isInMinBiasPD=( (bool) MB_HF1ORp5_I &&  (bool) mb_l1s_I);

    bool isInLowJetsPD=( (bool) Jet40_I   ||
			 (bool) Jet60_I   ||
			 (bool) Jet40_2_I ||
			 (bool) Jet60_2_I 
			 );
    if(filelistIsLowerJets)isInLowJetsPD=true;

    bool isInJet80PD  =( (bool) Jet80_I      || (bool) Jet100_I     || 
			 (bool) Jet80_2_I    || (bool) Jet100_2_I   || 
			 (bool) PFJet110_I   || (bool) PFJet120_I   ||
			 (bool) CaloJet110_I || (bool) CaloJet120_I || (bool) CaloJet150_I ||  
			 (bool) CaloJet100_Jet35_Eta0p7_I || (bool) CaloJet100_Jet35_Eta1p1_I  ||
			 (bool) CaloJet80_45_45_Eta2p1_I || (bool) CaloJet80_Jet35_Eta0p7_I || (bool) CaloJet80_Jet35_Eta1p1_I
			 );
    if(filelistIsJet80)isInJet80PD=true;
    
    if(isInJet80PD && !(isInLowJetsPD || isInMinBiasPD) )    h_NEvents_Jet80_PDexc->Fill(0.,1.);
    if(isInLowJetsPD && !(isInMinBiasPD || isInJet80PD) )    h_NEvents_LowJets_PDexc->Fill(0.,1.);
    if(isInMinBiasPD && !(isInJet80PD || isInLowJetsPD ) )   h_NEvents_MinBias_PDexc->Fill(0.,1.);    
    if((isInJet80PD && isInLowJetsPD) && !(isInMinBiasPD) )  h_NEvents_Jet80LowJets_NOTMinBias->Fill(0.,1.);
    if((isInLowJetsPD && isInMinBiasPD ) && !(isInJet80PD) ) h_NEvents_LowJetsMinBias_NOTJet80->Fill(0.,1.);
    if((isInMinBiasPD && isInJet80PD) && !(isInLowJetsPD) )  h_NEvents_MinBiasJet80_NOTLowJets->Fill(0.,1.);    
    if(isInJet80PD && isInLowJetsPD && isInMinBiasPD )       h_NEvents_Jet80LowJetsMinBias->Fill(0.,1.);

    

    // grab largest triggerPt among HLT paths' trigger object collectiosn
    double maxJetPt=-1.;
    for(int jet = 0; jet<nref_I; ++jet){
      if(fabs(eta_F[jet])<5.1){
	if(pt_F[jet]>maxJetPt)
	  maxJetPt=pt_F[jet];	      }
      else continue;
    }


    double maxTrgPt=-1.,maxTrgEta=0.; 
    double HLT40maxTrgPt=-1.,  HLT60maxTrgPt=-1.,  HLT80maxTrgPt=-1.,  HLT100maxTrgPt=-1.; 
    unsigned int trgObj40_size=trgObjpt_40->size(), trgObj60_size=trgObjpt_60->size();
    unsigned int trgObj80_size=trgObjpt_80->size(), trgObj100_size=trgObjpt_100->size();   
    if(trgDec[3]&&useHLT100){
      for(unsigned int itt=0; itt<trgObj100_size; ++itt){
	double trgpt=trgObjpt_100->at(itt);
	double trgeta=trgObjeta_100->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta = trgeta;
	  maxTrgPt = trgpt;
	}
	if(trgpt > HLT100maxTrgPt)
	  HLT100maxTrgPt=trgpt;
      }
    }
    

    if(trgDec[2]){
      for(unsigned int itt=0; itt<trgObj80_size; ++itt){
	
	double trgpt=trgObjpt_80->at(itt);
	double trgeta=trgObjeta_80->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta = trgeta;
	  maxTrgPt  = trgpt;
	}	  	
	if(trgpt > HLT80maxTrgPt)
	  HLT80maxTrgPt=trgpt;
      }
    }
    
    if(trgDec[1]) {
      for(unsigned int itt=0; itt<trgObj60_size; ++itt){
	double trgpt=trgObjpt_60->at(itt);
	//std::cout<<"trgpt HLT60 = "<< trgpt << std::endl;
	double trgeta=trgObjeta_60->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta =trgeta;
	  maxTrgPt  =trgpt;
	}
	if(trgpt > HLT60maxTrgPt)
	  HLT60maxTrgPt=trgpt;
      }    
    }
    
    if(trgDec[0])  {
      for(unsigned int itt=0; itt<trgObj40_size; ++itt)	{
	double trgpt=trgObjpt_40->at(itt);
	//std::cout<<"trgpt HLT40 = "<< trgpt << std::endl;
	double trgeta=trgObjeta_40->at(itt);
    	if(trgpt > maxTrgPt) { 
	  maxTrgEta = trgeta;
	  maxTrgPt  = trgpt;
	} 
	if(trgpt > HLT40maxTrgPt)
	  HLT40maxTrgPt=trgpt;
      }
    }
    //std::cout<<"hello5"<<std::endl;
    //assert(false);    

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
      else if( MBtrgDec ){ //&& trgPt<HLTthresh[0] ) { 
	//	// i want absolutely only TRUE min bias events; ones that happen to fire a jet trigger to are not going to be what i need.	 
	//	// if( !( ((bool) Jet40_I) || ((bool)Jet60_I) || ((bool)Jet80_I ) || ((bool) Jet40_2_I) || ((bool)Jet60_2_I) || ((bool)Jet80_2_I )  )) //PF and Calo triggers
	//	// if( !( ((bool) Jet40_I) || ((bool)Jet60_I) || ((bool)Jet80_I )   )) //just PF triggers
	//	// if(  !( ((bool) jet40_l1s_I) || ((bool) jet60_l1s_I)|| ((bool) jet80_l1s_I) || ((bool) jet100_l1s_I) )  )
	//	// if(  !(  ((bool) jet60_l1s_I)|| ((bool) jet80_l1s_I) || ((bool) jet100_l1s_I) )  )
	//	// if(  ! ( ((bool) jet80_l1s_I) || ((bool) jet100_l1s_I) )  )
	//	// if(  ! ( ((bool) jet100_l1s_I) )  )
	isMB=true; weight_eS=(double)MBtrgPscl;	    	
      }
    }
    
    

    
    //if(debugMode) std::cout<<"weight_eS="<<weight_eS<<std::endl;    
    bool firedJetTrigger=(is40 || is60 || is80 || is100 );
    bool firedMBTrigger=(isMB);
    bool firedTrigger=firedJetTrigger || firedMBTrigger;
    
    //event counts, pre-duplicate skip, pre-evt cuts, pre-trigger reqs.
    h_NEvents->Fill(0.);   h_NEvents->Fill(1.,weight_eS);      


    
    ////duplicate skipping between LowerJets and Jet80, old version extension
    if(filelistIsMinBias){
      if(isInJet80PD || isInLowJetsPD){
	h_NEvents_skipped->Fill(0.);		h_NEvents_skipped->Fill(1.,weight_eS);  
       	continue;
      }      
    }

    if(filelistIsJet80){//EXCLUSION FOR READING MIN BIAS DATA MAY NEED WORK TODO
      if(isInLowJetsPD){
	h_NEvents_skipped->Fill(0.);		h_NEvents_skipped->Fill(1.,weight_eS);  
	continue;      
      }
    }  
    
    if(filelistIsLowerJets){//do nothing in this duplicate skip version
    }
    
    
    
    ////// alternate duplicate skipping scheme; may be "more correct" TO DO
    //if(filelistIsMinBias)
    //  if(!isMB){
    //	//if(isInJet80PD || isInLowJetsPD){//some stuff i was string
    //	//      bool skipevt=false;
    //	//      if( (bool) Jet40_I   || 
    //	//	  (bool) Jet60_I   || 
    //	//	  (bool) Jet40_2_I || 
    //	//	  (bool) Jet60_2_I   ) skipevt=true;//in LowerJets
    //	//      else if(  (bool) Jet80_I   || 
    //	//		(bool) Jet100_I  ||
    //	//		(bool) Jet80_2_I || 
    //	//		(bool) Jet100_2_I  ) skipevt=true; //in Jet80
    //	//      else if( (bool) PFJet110_I || 
    //	//      	       (bool) PFJet120_I   ) skipevt=true; //in Jet80
    //	//      else if( (bool) CaloJet110_I || 
    //	//      	       (bool) CaloJet120_I || 
    //	//      	       (bool) CaloJet150_I    ) skipevt=true; //in Jet80
    //	//      else if( (bool) CaloJet100_Jet35_Eta0p7_I || 
    //	//      	       (bool) CaloJet100_Jet35_Eta1p1_I   ) skipevt=true; //in Jet80
    //	//      else if( (bool) CaloJet80_45_45_Eta2p1_I || 
    //	//      	       (bool) CaloJet80_Jet35_Eta0p7_I || 
    //	//      	       (bool) CaloJet80_Jet35_Eta1p1_I   ) skipevt=true; //in Jet80
    //	//      else skipevt=false;      
    //	//      if(skipevt){
    //	//	//if(debugMode)std::cout<<"this event is in Jet80 AND/OR LowerJets dataset!"<<std::endl;
    //	//	//if(debugMode)std::cout<<"Skipping event, will read it in Jet80 OR LowerJets instead!"<<std::endl;
    //	h_NEvents_skipped->Fill(0.);		h_NEvents_skipped->Fill(1.,weight_eS);  
    //   	continue;
    //  }      
    //
    //
    //
    //if(filelistIsJet80) //EXCLUSION FOR READING MIN BIAS DATA MAY NEED WORK TODO
    //  if( !(is80 || is100) ){
    //	//	//if(debugMode)std::cout<<"this event is in Jet80 AND LowerJets dataset!"<<std::endl;
    //	//	//if(debugMode)std::cout<<"Skipping event, will read it in LowerJets instead!"<<std::endl;	
    //	h_NEvents_skipped->Fill(0.);		h_NEvents_skipped->Fill(1.,weight_eS);  
    //	continue;      
    //  }
    //
    //
    //if(filelistIsLowerJets)
    //  if( !(is40 || is60) ){
    //	//	if(fillDataJetTrigQAHists){
    //	//	  //b.c. it makes the hist look scary to do this kind of duplicate skip. 
    //	//	  //so i'm gonna put this here temporarily to test my understanding + hopefully make sure it doesn't look scary...	  
    //	//	  if( trgDec[0] ) hpp_IncHLT40trgPt_2-> Fill(HLT40maxTrgPt, trgPscl[0]);      	
    //	//	  if( trgDec[0] && !(trgPt<HLTthresh[0]) )   hpp_IncHLT40trgPt->Fill(  trgPt, (double)trgPscl[0] );      
    //	//	}
    //	h_NEvents_skipped->Fill(0.);		h_NEvents_skipped->Fill(1.,weight_eS);  
    //	continue;      
    //  }
    //

  
    h_NEvents_read->Fill(0.);    h_NEvents_read->Fill(1.,weight_eS);        
    
    //event counts, post-duplicate skip, pre-evt cuts, having fired a trigger
    if(firedTrigger)
      h_NEvents_trigd_1->Fill(0.);   h_NEvents_trigd_1->Fill(1.,weight_eS );    
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I    ==0    || 
	pBeamScrapingFilter_I ==0 || 
	pprimaryvertexFilter_I==0  )  continue;  
    h_NEvents_skimCut->Fill(0.);  h_NEvents_skimCut->Fill(1.,weight_eS); 
    
    if( fabs(vz_F)>24. )     continue;
    h_NEvents_vzCut->Fill(0.); h_NEvents_vzCut->Fill(1.,weight_eS);    
    
    if( firedTrigger )
      h_NEvents_trigd_2->Fill(0.);      h_NEvents_trigd_2->Fill(1.,weight_eS);        
    
    //event counts excl trigd, post dup skip + event quality
    
    if(MBtrgDec  )  {h_NEvents_HLTMB  -> Fill(0.); h_NEvents_HLTMB  -> Fill(1.,(double)MBtrgPscl  );}//minbias
    if(trgDec[0])   {h_NEvents_HLT40  -> Fill(0.); h_NEvents_HLT40  -> Fill(1.,(double)trgPscl[0]);}
    if(trgDec[1])   {h_NEvents_HLT60  -> Fill(0.); h_NEvents_HLT60  -> Fill(1.,(double)trgPscl[1]);}
    if(trgDec[2])   {h_NEvents_HLT80  -> Fill(0.); h_NEvents_HLT80  -> Fill(1.,(double)trgPscl[2]);}
    if(trgDec[3])   {h_NEvents_HLT100 -> Fill(0.); h_NEvents_HLT100 -> Fill(1.,(double)trgPscl[3]);}

    if(trgDec[0]&&!(trgPt<HLTthresh[0])) {h_NEvents_HLT40thresh  -> Fill(0.); h_NEvents_HLT40thresh  -> Fill(1.,weight_eS);}
    if(trgDec[1]&&!(trgPt<HLTthresh[1])) {h_NEvents_HLT60thresh  -> Fill(0.); h_NEvents_HLT60thresh  -> Fill(1.,weight_eS);}
    if(trgDec[2]&&!(trgPt<HLTthresh[2])) {h_NEvents_HLT80thresh  -> Fill(0.); h_NEvents_HLT80thresh  -> Fill(1.,weight_eS);}
    if(trgDec[3]&&!(trgPt<HLTthresh[3])) {h_NEvents_HLT100thresh -> Fill(0.); h_NEvents_HLT100thresh -> Fill(1.,weight_eS);}

    if     ( isMB  )  { NEvents_MB++  ;  h_NEvents_isMB ->Fill(0.)  ;  h_NEvents_isMB ->Fill(1.,weight_eS ) ; }//minbias
    else if( is40  )  { NEvents_40++  ;  h_NEvents_is40 ->Fill(0.)  ;  h_NEvents_is40 ->Fill(1.,weight_eS ) ; }                        
    else if( is60  )  { NEvents_60++  ;  h_NEvents_is60 ->Fill(0.)  ;  h_NEvents_is60 ->Fill(1.,weight_eS ) ; }
    else if( is80  )  { NEvents_80++  ;  h_NEvents_is80 ->Fill(0.)  ;  h_NEvents_is80 ->Fill(1.,weight_eS ) ; }
    else if( is100 )  { NEvents_100++ ;  h_NEvents_is100->Fill(0.)  ;  h_NEvents_is100->Fill(1.,weight_eS ) ; }
    
    if( MBtrgDec && ((bool) l1jet8_l1s_I  )  )  {   h_NEvents_HLTMB_AND_L1jet8  -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet8  ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) l1jet12_l1s_I )  )  {   h_NEvents_HLTMB_AND_L1jet12 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet12 ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) l1jet16_l1s_I )  )  {   h_NEvents_HLTMB_AND_L1jet16 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet16 ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) l1jet20_l1s_I )  )  {   h_NEvents_HLTMB_AND_L1jet20 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet20 ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) l1jet24_l1s_I )  )  {   h_NEvents_HLTMB_AND_L1jet24 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet24 ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) jet40_l1s_I   )  )  {   h_NEvents_HLTMB_AND_L1jet28 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet28 ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) jet60_l1s_I   )  )  {   h_NEvents_HLTMB_AND_L1jet40 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet40 ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) jet80_l1s_I   )  )  {   h_NEvents_HLTMB_AND_L1jet48 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet48 ->Fill( 1., MBtrgPscl ) ; }//minbias
    if( MBtrgDec && ((bool) jet100_l1s_I  )  )  {   h_NEvents_HLTMB_AND_L1jet52 -> Fill( 0. )  ;  h_NEvents_HLTMB_AND_L1jet52 ->Fill( 1., MBtrgPscl ) ; }//minbias    

    if( isMB && ((bool) l1jet8_l1s_I  )  )  {   h_NEvents_isMB_AND_L1jet8  -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet8  ->Fill(1.,weight_eS ) ; }//minbias
    if( isMB && ((bool) l1jet12_l1s_I )  )  {   h_NEvents_isMB_AND_L1jet12 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet12 ->Fill(1.,weight_eS ) ; }//minbias
    if( isMB && ((bool) l1jet16_l1s_I )  )  {   h_NEvents_isMB_AND_L1jet16 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet16 ->Fill(1.,weight_eS ) ; }//minbias
    if( isMB && ((bool) l1jet20_l1s_I )  )  {   h_NEvents_isMB_AND_L1jet20 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet20 ->Fill(1.,weight_eS ) ; }//minbias
    if( isMB && ((bool) l1jet24_l1s_I )  )  {   h_NEvents_isMB_AND_L1jet24 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet24 ->Fill(1.,weight_eS ) ; }//minbias
    if( isMB && ((bool) jet40_l1s_I   )  )  {   h_NEvents_isMB_AND_L1jet28 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet28 ->Fill(1.,weight_eS ) ; }//minbias 
    if( isMB && ((bool) jet60_l1s_I   )  )  {   h_NEvents_isMB_AND_L1jet40 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet40 ->Fill(1.,weight_eS ) ; }//minbias
    if( isMB && ((bool) jet80_l1s_I   )  )  {   h_NEvents_isMB_AND_L1jet48 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet48 ->Fill(1.,weight_eS ) ; }//minbias
    if( isMB && ((bool) jet100_l1s_I  )  )  {   h_NEvents_isMB_AND_L1jet52 -> Fill(0.)  ;  h_NEvents_isMB_AND_L1jet52 ->Fill(1.,weight_eS ) ; }//minbias
    
    if(fillDataJetTrigQAHists){ //only want to fill these trigger jet plots if they pass all our quality criteria      

      //prescales
      if(MBtrgDec)hpp_IncHLTMBpscl  ->Fill( (double)MBtrgPscl, 1.);
      if(trgDec[0])hpp_IncHLT40pscl ->Fill( (double)trgPscl[0], 1.);
      if(trgDec[1])hpp_IncHLT60pscl ->Fill( (double)trgPscl[1], 1.);
      if(trgDec[2])hpp_IncHLT80pscl ->Fill( (double)trgPscl[2], 1.);
      if(trgDec[3])hpp_IncHLT100pscl->Fill( (double)trgPscl[3], 1.);
      //if(MBtrgDec && trgDec[0])hpp_IncHLT40MBpscl  ->Fill( (float)MBtrgPscl, 1.);


      //all trigger jet object pts
      if(trgDec[0])
	for(unsigned int itt=0; itt<trgObj40_size; ++itt)hpp_IncHLT40trgPt_allobj->Fill(trgObjpt_40->at(itt),1.);
      if(trgDec[1])
	for(unsigned int itt=0; itt<trgObj60_size; ++itt)hpp_IncHLT60trgPt_allobj->Fill(trgObjpt_60->at(itt),1.);
      if(trgDec[2])
	for(unsigned int itt=0; itt<trgObj80_size; ++itt)hpp_IncHLT80trgPt_allobj->Fill(trgObjpt_80->at(itt),1.);
      if(trgDec[3])
	for(unsigned int itt=0; itt<trgObj100_size; ++itt)hpp_IncHLT100trgPt_allobj->Fill(trgObjpt_100->at(itt),1.);

      //lead trigger jet object pt for each respective trigger
      if( trgDec[0] )  hpp_IncHLT40trgPt_leadobj-> Fill(HLT40maxTrgPt , 1.);      	
      if( trgDec[1] )  hpp_IncHLT60trgPt_leadobj-> Fill(HLT60maxTrgPt , 1.);
      if( trgDec[2] )  hpp_IncHLT80trgPt_leadobj-> Fill(HLT80maxTrgPt , 1.);
      if( trgDec[3] ) hpp_IncHLT100trgPt_leadobj->Fill(HLT100maxTrgPt, 1.);           
      if( MBtrgDec && trgDec[0] ) hpp_IncHLT40MBtrgPt_leadobj-> Fill(HLT40maxTrgPt , 1.);      	 

      //lead trigger jet object pt v lead reco jet object pt
      if( trgDec[0] ) hpp_IncHLT40trgPt_leadobj_leadjet->Fill( HLT40maxTrgPt , maxJetPt,  1.);      	
      if( trgDec[1] ) hpp_IncHLT60trgPt_leadobj_leadjet->Fill( HLT60maxTrgPt , maxJetPt,  1.);
      if( trgDec[2] ) hpp_IncHLT80trgPt_leadobj_leadjet->Fill( HLT80maxTrgPt , maxJetPt,  1.);
      if( trgDec[3] ) hpp_IncHLT100trgPt_leadobj_leadjet->Fill(HLT100maxTrgPt, maxJetPt,  1.);            
      if( trgDec[0] && MBtrgDec ) hpp_IncHLT40MBtrgPt_leadobj_leadjet->Fill( HLT40maxTrgPt , maxJetPt,  1.);      	

      if( trgDec[0] ) hpp_IncHLT40trgPt_leadobj_globleadobj->Fill( HLT40maxTrgPt , trgPt,  1.);      	
      if( trgDec[1] ) hpp_IncHLT60trgPt_leadobj_globleadobj->Fill( HLT60maxTrgPt , trgPt,  1.);
      if( trgDec[2] ) hpp_IncHLT80trgPt_leadobj_globleadobj->Fill( HLT80maxTrgPt , trgPt,  1.);
      if( trgDec[3] ) hpp_IncHLT100trgPt_leadobj_globleadobj->Fill(HLT100maxTrgPt, trgPt,  1.);            
      if( trgDec[0] && MBtrgDec ) hpp_IncHLT40MBtrgPt_leadobj_globleadobj->Fill( HLT40maxTrgPt , trgPt,  1.);      	

      //global lead trigger jet object (lead obj across all triggers) by trigger decision + lower threshold only
      if( MBtrgDec  &&  trgPt<HLTthresh[0]    )   hpp_IncHLTMBtrgPt->Fill(  trgPt, (double)MBtrgPscl  );
      if( trgDec[0] && !(trgPt<HLTthresh[0]) )   hpp_IncHLT40trgPt->Fill(  trgPt, (double)trgPscl[0] );
      if( trgDec[1] && !(trgPt<HLTthresh[1]) )   hpp_IncHLT60trgPt->Fill(  trgPt, (double)trgPscl[1] );
      if( trgDec[2] && !(trgPt<HLTthresh[2]) )   hpp_IncHLT80trgPt->Fill(  trgPt, (double)trgPscl[2] );
      if( trgDec[3] && !(trgPt<HLTthresh[3]) )   hpp_IncHLT100trgPt->Fill( trgPt, (double)trgPscl[3] );

      // trgEta [who cares really]
      if( trgDec[0] && !(trgPt<HLTthresh[0]) )   hpp_IncHLT40trgEta->Fill(  trgEta, (double)trgPscl[0] );
      if( trgDec[1] && !(trgPt<HLTthresh[1]) )   hpp_IncHLT60trgEta->Fill(  trgEta, (double)trgPscl[1] );
      if( trgDec[2] && !(trgPt<HLTthresh[2]) )   hpp_IncHLT80trgEta->Fill(  trgEta, (double)trgPscl[2] );
      if( trgDec[3] && !(trgPt<HLTthresh[3]) )   hpp_IncHLT100trgEta->Fill( trgEta, (double)trgPscl[3] );
      

      //global lead jet trigger object across exclusively firing triggers (i.e. meeting all exclusion criteria)
      if     ( isMB  )   hpp_HLTMBtrgPt->Fill(   trgPt , weight_eS );
      else if( is40  )   hpp_HLT40trgPt->Fill(   trgPt , weight_eS );      
      else if( is60  )   hpp_HLT60trgPt->Fill(   trgPt , weight_eS );
      else if( is80  )   hpp_HLT80trgPt->Fill(   trgPt , weight_eS );
      else if( is100 )  hpp_HLT100trgPt->Fill(   trgPt , weight_eS );   
      if(firedJetTrigger)	hpp_HLTCombtrgPt->Fill(  trgPt , weight_eS );   

      if     ( is40  )  hpp_HLT40trgEta->Fill(    trgEta , weight_eS );      
      else if( is60  )  hpp_HLT60trgEta->Fill(    trgEta , weight_eS );
      else if( is80  )  hpp_HLT80trgEta->Fill(    trgEta , weight_eS );
      else if( is100 )  hpp_HLT100trgEta->Fill(   trgEta , weight_eS );   
      if(firedJetTrigger)         hpp_HLTCombtrgEta->Fill(  trgEta , weight_eS );                         

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
    //if(useTupel)
    //    else jetsPerEvent=nref_I;
    
    
    if(useIncJetAnalyzer){
      
      jetsPerEvent=nref_I;
      for(int jet = 0; jet<nref_I; ++jet){
	
	float jtpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
	float rawpt=rawpt_F[jet];      
	float receta = eta_F[jet];
	float absreceta = fabs(receta);
	float recphi = phi_F[jet];
	
	int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
	int neuMult = neN_I[jet] + phN_I[jet] ;
	int numConst  = chMult + neuMult;
	
	
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
	
	
	// readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0
	float L2RelL3Abs_rawpt=-1.;//L2RelL3Abs JEC
	int L2RelL3AbsJECetabin=-1;
	float L2L3Res_rawpt=-1.;//L2L3Res JEC
	int L2L3ResJECetabin=-1;

	// readForests_JEC_v12_ak4PF_75X_data_Run2_v13
	float L2RelL3Abs_2_rawpt=-1.;//L2RelL3Abs JEC
	int L2RelL3Abs_2_JECetabin=-1;
	float L2L3Res_2_rawpt=-1.;//L2L3Res JEC
	int L2L3Res_2_JECetabin=-1;
	

	
	if(fillDataJetJECQAHists){

	  //readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0
	  for(int uu=0; uu<1; uu++){//loop so that if the jet doesn't meet some criteria i can just quit while i'm ahead
	    for(int i=0; i<n_L2RelL3AbsJEC_etabins;i++){	    //L2RelL3Abs
	      if(receta>=L2RelL3AbsJEC_etabins[i] && receta < L2RelL3AbsJEC_etabins[i+1]) { 
		L2RelL3AbsJECetabin=i; break;}
	      else continue;	    }
	    if(L2RelL3AbsJECetabin==-1)
	      continue;
	    
	    L2RelL3Abs_rawpt=rawpt*(L2RelL3AbsJEC_arr[L2RelL3AbsJECetabin]->Eval(rawpt));		
	    
	    for(int i=0; i<n_L2L3ResJEC_etabins;i++){	    //L2L3Res
	      if(receta>=L2L3ResJEC_etabins[i] && receta < L2L3ResJEC_etabins[i+1]) { 
		L2L3ResJECetabin=i; break;}
	      else continue;
	    }
	    if(L2L3ResJECetabin==-1)
	      continue;
	    
	    L2L3Res_rawpt=L2RelL3Abs_rawpt*(L2L3ResJEC_arr[L2L3ResJECetabin]->Eval(L2RelL3Abs_rawpt));	  
	    
	  }//end loop for easy continuing


	  ////readForests_JEC_v12_ak4PF_75X_data_Run2_v13
	  for(int uu=0; uu<1; uu++){//loop so that if the jet doesn't meet some criteria i can just quit while i'm ahead
	    for(int i=0; i<n_L2RelL3AbsJEC_2_etabins;i++){	    //L2RelL3Abs
	      if(receta>=L2RelL3AbsJEC_2_etabins[i] && receta < L2RelL3AbsJEC_2_etabins[i+1]) { 
		L2RelL3Abs_2_JECetabin=i; break;}
	      else continue;	    }
	    if(L2RelL3Abs_2_JECetabin==-1)
	      continue;
	    
	    L2RelL3Abs_2_rawpt=rawpt*(L2RelL3AbsJEC_2_arr[L2RelL3Abs_2_JECetabin]->Eval(rawpt));		
	    
	    for(int i=0; i<n_L2L3ResJEC_2_etabins;i++){	    //L2L3Res
	      if(receta>=L2L3ResJEC_2_etabins[i] && receta < L2L3ResJEC_2_etabins[i+1]) { 
		L2L3Res_2_JECetabin=i; break;}
	      else continue;
	    }
	    if(L2L3Res_2_JECetabin==-1)
	      continue;
	    
	    L2L3Res_2_rawpt=L2RelL3Abs_2_rawpt*(L2L3ResJEC_2_arr[L2L3Res_2_JECetabin]->Eval(L2RelL3Abs_2_rawpt));	  
	    
	  }//end loop for easy continuing
	}//end fill DataJetJECQAHists
	

	// 13 TeV JetID criterion, loose or tight
	bool passesJetID=false;
	//float jetIDpt=jtpt;//ala HIN jetID, jtpt is corrected w/ L2/L3 residuals
	
	if (!(absreceta > 2.4)) 
	  passesJetID=(bool)jetID_00eta24( jtpt, 
					   neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet], muSum_F[jet],
					   numConst,  chMult);
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  passesJetID=(bool) jetID_24eta27( jtpt,
					    neSum_F[jet],  phSum_F[jet], muSum_F[jet],
					    numConst);
	else if( !(absreceta>3.0) && absreceta>2.7 )
	  passesJetID=(bool) jetID_27eta30( jtpt,
					    neSum_F[jet],  phSum_F[jet], 
					    neuMult);
	else  
	  passesJetID=(bool)jetID_32eta47( jtpt, 
					   phSum_F[jet]);
	
	
	// get repidity bin
	int theRapBin=-1;	
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin)
	  if( absreceta>=absetabins[rapbin]  && 		
	      absreceta<absetabins[rapbin+1]    	      ) {	    
	    theRapBin=rapbin;
	    break;
	  }
	

	if(fillDataJetJECQAHists && passesJetID ){
	  hjtpt[theRapBin]->Fill(jtpt,weight_eS);//spectra
	  hL2L3Res_rawpt[theRapBin]->Fill(L2L3Res_rawpt,weight_eS);
	  hL2L3Res_2_rawpt[theRapBin]->Fill(L2L3Res_2_rawpt,weight_eS);
	  
	  hjtpt_L2L3Res_rawpt_diff[theRapBin]->Fill((jtpt-L2L3Res_rawpt),weight_eS);//diff
	  hjtpt_L2L3Res_2_rawpt_diff[theRapBin]->Fill((jtpt-L2L3Res_2_rawpt),weight_eS);
	  hL2L3Res_rawpt_L2L3Res_2_rawpt_diff[theRapBin]->Fill((L2L3Res_rawpt-L2L3Res_2_rawpt),weight_eS);
	  
	  hjtpt_L2L3Res_rawpt_ratio[theRapBin]->Fill((jtpt),weight_eS);//ratios (to be divided by corresponding spectra post processing)
	  hjtpt_L2L3Res_2_rawpt_ratio[theRapBin]->Fill((jtpt),weight_eS);
	  hL2L3Res_rawpt_L2L3Res_2_rawpt_ratio[theRapBin]->Fill((L2L3Res_rawpt),weight_eS);
	  
	  hjtpt_v_diff_L2L3Res_rawpt[theRapBin]->Fill(jtpt, (jtpt-L2L3Res_rawpt), weight_eS);//pt v. diff TH2
	  hjtpt_v_diff_L2L3Res_2_rawpt[theRapBin]->Fill(jtpt, (jtpt-L2L3Res_2_rawpt), weight_eS);
	  hL2L3Res_rawpt_v_diff_L2L3Res_2_rawpt[theRapBin]->Fill(L2L3Res_rawpt, (L2L3Res_rawpt-L2L3Res_2_rawpt), weight_eS);
	}





	double jtpt_sysup=-1., jtpt_sysdown=-1.;
	if(fillDataJetJECUncHists && passesJetID){
	  //std::cout<<std::endl;
	  //std::cout<<"doing JEC sysup/sysdown hists"<<std::endl;
	  int JECUnc_etabin=-1;
	  int JECUnc_ptbin=-1;
	  for(int uu=0; uu<1; uu++){//loop so that if the jet doesn't meet some criteria i can just quit while i'm ahead
	    
	    //find JEC eta bin
	    //for(int i=0; i<n_jecunc_etabins;i++){	  
	    for(int i=jecunc_eta_start; i<=jecunc_eta_end;i++){	  
	      if(receta>=jecunc_etabins[i] && receta< jecunc_etabins[i+1]){
		JECUnc_etabin=i;break;
	      }
	      else continue;}
	    //std::cout<<"receta="<<receta<<std::endl;
	    //std::cout<<"JECUnc_etabin="<<JECUnc_etabin<<std::endl;
	    if(JECUnc_etabin==-1)continue;
	    //std::cout<<"jecunc_etabins["<<JECUnc_etabin<<"]="<<jecunc_etabins[JECUnc_etabin]<<std::endl;
	    //std::cout<<"jecunc_etabins["<<JECUnc_etabin+1<<"]="<<jecunc_etabins[JECUnc_etabin+1]<<std::endl;
	    
	    
	    //find JEC pt bin
	    //	      for(int i=0; i<n_jecunc_ptbins;i++){	  
	    for(int i=jecunc_pt_start; i<=jecunc_pt_end;i++){	  
	      if(jtpt>=jecunc_ptbins[i] && jtpt< jecunc_ptbins[i+1]){
		JECUnc_ptbin=i;break;
	      }
	      else continue;}
	    //std::cout<<"jtpt="<<jtpt<<std::endl;
	    //std::cout<<"JECUnc_ptbin="<<JECUnc_ptbin<<std::endl;
	    if(JECUnc_ptbin==-1)continue;
	    //std::cout<<"jecunc_ptbins["<<JECUnc_ptbin<<"]="<<jecunc_ptbins[JECUnc_ptbin]<<std::endl;
	    //std::cout<<"jecunc_ptbins["<<JECUnc_ptbin+1<<"]="<<jecunc_ptbins[JECUnc_ptbin+1]<<std::endl;
	    
	    
	    //find JEC pt bin
	    double JECUnc=jecunc[JECUnc_etabin][JECUnc_ptbin];
	    //std::cout<<"JECUnc="<<JECUnc<<std::endl;	    
	    jtpt_sysup=jtpt+JECUnc*jtpt;	      
	    jtpt_sysdown=jtpt-JECUnc*jtpt;
	    //std::cout<<"jtpt_sysup="<<jtpt_sysup<<std::endl;
	    //std::cout<<"jtpt_sysdown="<<jtpt_sysdown<<std::endl;
	    hJetQA_jtpt_JEC_sysup[theRapBin]  ->Fill(jtpt_sysup,  weight_eS);
	    hJetQA_jtpt_JEC_sysdown[theRapBin]->Fill(jtpt_sysdown,weight_eS);
	    if(fillDataJetJECUncHists){
	      JEC_etabin_bool_helper[theRapBin]=true;      
	      hpp_covmat_eta_arr_helpers_JEC_sysup[theRapBin]->Fill(jtpt_sysup);//,weight_eS);	    
	      hpp_covmat_eta_arr_helpers_JEC_sysdown[theRapBin]->Fill(jtpt_sysdown);//,weight_eS);	    
	    }      
	    //assert(false);
	    
	  }//end JECUnc loop
	}//end JECUnc hists + JetID	
	
	
	if(!(jtpt>jetQAPtCut)){ jetsPerEvent--; continue;}//need JEC syst. stuff before this cut; so jets from below 56 GeV can migrate above	
	h_NJets_jetQAPtCut->Fill(0.);	  h_NJets_jetQAPtCut->Fill(1., weight_eS);
	
	
	//STUFF I WANT FILLED FOR 56 GEV AND ABOVE ONLY
	for(int jtid=0; jtid<2; jtid++){
	  
	  if(jtid==1 ){
	    
	    if(!passesJetID){ jetsPerEvent--; continue;	  }
	    else {  h_NJets_JetIDCut->Fill(0.);  h_NJets_JetIDCut->Fill(1., weight_eS); }
	    
	    if(jtpt>evt_leadJetPt_wCuts)
	      evt_leadJetPt_wCuts=jtpt;      
	    
	    if(fillDataJetCovMatrix){
	      etabin_bool_helper[theRapBin]=true;      
	      hpp_covmat_eta_arr_helpers[theRapBin]->Fill(jtpt);//,weight_eS);	    		      
	    }
	    
	  }//end jtid==1 specific. 	  
	  
	  
	  
	  
	  // trig plots
	  if(fillDataJetTrigQAHists){
	    
	    hpp_HLTComb[theRapBin][jtid]->Fill(jtpt, weight_eS); //kurts method
	    
	    //no exclusion between triggers
	    if(MBtrgDec ){
	      hpp_HLTMB[theRapBin][jtid]->Fill(jtpt, (double)MBtrgPscl); //minbias
	      if( trgPt < HLTthresh[0] ) {
		hpp_incHLTMB[theRapBin][jtid]->Fill(jtpt, (double)MBtrgPscl); //minbisa
	      }
	    }
	    
	    if(trgDec[0]){
	      hpp_HLT40[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[0]);  
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		if(!MBtrgDec)hpp_incHLT40_NOT_HLTMB->Fill(jtpt, (double)trgPscl[0]);
	      }		
	      if( (trgPt>HLTthresh[0]) && (trgPt<HLTthresh[1]) )  { 
		hpp_incHLT40[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[0]);  
	      }		
	    }
	    
	    if(trgDec[1] ){
	      hpp_HLT60[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[1]);  
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		if(!MBtrgDec)hpp_incHLT60_NOT_HLTMB->Fill(jtpt, (double)trgPscl[1]);
	      }		
	      if( (trgPt>HLTthresh[1]) && (trgPt<HLTthresh[2]))  { 
		hpp_incHLT60[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[1]);  
	      }
	    }

	    if(trgDec[2]){
	      hpp_HLT80[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[2]);  
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		if(!MBtrgDec)hpp_incHLT80_NOT_HLTMB->Fill(jtpt, (double)trgPscl[2]);
	      }				
	      if((trgPt>HLTthresh[2]))  {//this would have to change if useHLT100 was turned on
		hpp_incHLT80[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[2]);  
	      }
	    }
	    
	    if(trgDec[3]){
	      hpp_HLT100[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[3]); 
	      if(!(trgPt<HLTthresh[3]))  { //this would have to change if useHLT100 was turned on
		hpp_incHLT100[theRapBin][jtid]->Fill(jtpt, (double)trgPscl[3]); 
	      } 
	    }
	    
	    //with exclusion
	    if(isMB)   {
	      hpp_excHLTMB[theRapBin][jtid]->Fill( jtpt , weight_eS);	//minbias	
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		if(      (bool)l1jet68_l1s_I  ) hpp_excHLTMB_AND_l1jet68->Fill( jtpt, weight_eS);
		else if( (bool)l1jet60_l1s_I  ) hpp_excHLTMB_AND_l1jet60->Fill( jtpt, weight_eS);
		else if( (bool)jet100_l1s_I   ) hpp_excHLTMB_AND_l1jet52->Fill( jtpt, weight_eS);
		else if( (bool)jet80_l1s_I    ) hpp_excHLTMB_AND_l1jet48->Fill( jtpt, weight_eS);
		else if( (bool)jet60_l1s_I    ) hpp_excHLTMB_AND_l1jet40->Fill( jtpt, weight_eS);
		else if( (bool)jet40_l1s_I    ) hpp_excHLTMB_AND_l1jet28->Fill( jtpt, weight_eS);
		else if( (bool)l1jet24_l1s_I  ) hpp_excHLTMB_AND_l1jet24->Fill( jtpt, weight_eS);
		else if( (bool)l1jet20_l1s_I  ) hpp_excHLTMB_AND_l1jet20->Fill( jtpt, weight_eS);
		else if( (bool)l1jet16_l1s_I  ) hpp_excHLTMB_AND_l1jet16->Fill( jtpt, weight_eS);
		else if( (bool)l1jet12_l1s_I  ) hpp_excHLTMB_AND_l1jet12->Fill( jtpt, weight_eS);
		else if( (bool)l1jet8_l1s_I   ) hpp_excHLTMB_AND_l1jet8 ->Fill( jtpt, weight_eS);
		else                            hpp_excHLTMB_NOT_l1jet  ->Fill( jtpt, weight_eS);
		
		if (trgDec[2])      hpp_excHLTMB_AND_HLTjet80->Fill(jtpt,weight_eS);
		else if (trgDec[1]) hpp_excHLTMB_AND_HLTjet60->Fill(jtpt,weight_eS);
		else if (trgDec[0]) hpp_excHLTMB_AND_HLTjet40->Fill(jtpt,weight_eS);
		else                hpp_excHLTMB_NOT_HLTjet  ->Fill(jtpt,weight_eS);
		
		if(!((bool)L1_HF_coincidence_I)) hpp_excHLTMB_NOT_HFcoincidence->Fill(jtpt,weight_eS);
	      }		
	    }
	    else if(is40 )  {
	      hpp_excHLT40[theRapBin][jtid]->Fill( jtpt , weight_eS);
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		if(!MBtrgDec)hpp_excHLT40_NOT_HLTMB->Fill(jtpt, weight_eS);
	      }
	    }
	    else if(is60 ) { 
	      hpp_excHLT60[theRapBin][jtid]->Fill( jtpt , weight_eS);
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		if(!MBtrgDec)hpp_excHLT60_NOT_HLTMB->Fill(jtpt, weight_eS);
	      }
	    }
	    else if(is80 ) { 
	      hpp_excHLT80[theRapBin][jtid]->Fill( jtpt , weight_eS);
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		if(!MBtrgDec)hpp_excHLT80_NOT_HLTMB->Fill(jtpt, weight_eS);
	      }
	    }
	    else if(is100) { 
	      hpp_excHLT100[theRapBin][jtid]->Fill( jtpt , weight_eS); 	    //assert(false);
	    }
	    
	  }//end fill data jettrig qa hists	
	    
	    
	  
	  
	  if(fillDataJetQAHists){
	    
	    int ind=0;
	    //if(jtpt>jetQAPtCut){//second jet pt cut
	    
	    //jets	
	    hJetQA_jtptEntries[theRapBin][jtid]->Fill(jtpt,1.);	  	
	    hJetQA[theRapBin][ind][jtid]->Fill(jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(rawpt_F[jet], weight_eS); ind++;	    	    
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
	    //}
	    
	    
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
    
    
    if(useTupel){
      
      jetsPerEvent=jetPt->size();
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
	    if(isMB)   {
	      hpp_excHLTMB[theRapBin][jtid]->Fill( jtpt , weight_eS);	//minbias
	    }	    
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
	}//end normal etabin hlelper	  
	etabin_bool_helper[i]=false;	
      }//end normal etabin loop
      
      if(fillDataJetJECUncHists){
	for(int i=0; i<nbins_abseta; i++){
	  if(JEC_etabin_bool_helper[i]){//don't undergo expensive covmatrix calc if the etabin wasn't filled.
	    fillCovMatrix( (TH2D*)hpp_covmat_eta_arr_JEC_sysup[i], (TH1D*)hpp_covmat_eta_arr_helpers_JEC_sysup[i], nbins_pt , (double) weight_eS);
	    hpp_covmat_eta_arr_helpers_JEC_sysup[i]->Reset();//reset contents of TH1 without resetting binning, min/max, etc....	  	    
	    fillCovMatrix( (TH2D*)hpp_covmat_eta_arr_JEC_sysdown[i], (TH1D*)hpp_covmat_eta_arr_helpers_JEC_sysdown[i], nbins_pt , (double) weight_eS);
	    hpp_covmat_eta_arr_helpers_JEC_sysdown[i]->Reset();//reset contents of TH1 without resetting binning, min/max, etc....	  	    
	  }//end JECUnc specific etabin helper
	  JEC_etabin_bool_helper[i]=false;
	}//end JECUnc specific etabin loop	      
      }//end JECUnc specific covmatrix
      
    }//end covariance matrix


    
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
  //assert(false);
  
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








// OLD CODE IANS NOT READY TO LET GO OF YET BUT PROBABLY COULD: TOO BAD HES PARANOID

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
      











      //if(true)         hpp_HLTCombtrgEta->Fill(  trgEta , weight_eS );                   

      







      //if(true)    // {
      //if(firedJetTrigger)     {

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
      //}
      

	      




  
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









	//	if(isMB){
	//	  if( ( pt_F[jet]>119.      )        && 
	//	      ( pt_F[jet]<120.      )      &&
	//	      ( fabs(eta_F[jet])<0.5)       && 
	//	      ((bool)jet100_l1s_I)){
	//	    std::cout << "evt w/ crazy jet(s) found!"<<std::endl;
	//	    std::cout << "evt  =" << evt_I  << std::endl;
	//	    std::cout << "run  =" << run_I  << std::endl;
	//	    std::cout << "lumi =" << lumi_I << std::endl;
	//	    std::cout << "vz   =" << vz_F << std::endl;
	//	    std::cout << "vx   =" << vx_F << std::endl;
	//	    std::cout << "vy   =" << vy_F << std::endl;
	//	    std::cout << "nVtx =" << nVtx_I << std::endl;
	//	    std::cout << "nTrk =" << nTrk_I<< std::endl;		
	//	    std::cout << std::endl;
	//	    std::cout << "weight_eS = "<<weight_eS<<std::endl;
	//	    std::cout << "MBtrgDec  = "<<MBtrgDec<<std::endl;
	//	    std::cout << "MBtrgPscl = "<< MBtrgPscl<<std::endl;
	//	    std::cout << "MB_HF1ORp5_I   = "<< MB_HF1ORp5_I <<std::endl;
	//	    std::cout << "MB_HF1ORp5_p_I = "<< MB_HF1ORp5_p_I<<std::endl;
	//	    std::cout << "mb_l1s_I    = "<<mb_l1s_I <<std::endl;
	//	    std::cout << "mb_l1s_ps_I = "<< mb_l1s_ps_I<<std::endl;
	//	    std::cout << std::endl;
	//	    std::cout<<"here's the other jets in the event"<<std::endl;
	//	    std::cout << "nref_I="<<nref_I<<std::endl;
	//	    for(int j=0; j<nref_I; j++){
	//	      std::cout << std::endl;
	//	      std::cout << "// ------------------------------- //"<<std::endl;		  
	//	      std::cout << "jet # = "<< jet       <<std::endl;
	//	      std::cout << "rawpt = "<< rawpt_F[j]<<std::endl;
	//	      std::cout << "pt    = "<< pt_F[j]   <<std::endl;
	//	      std::cout << "eta   = "<< eta_F[j]  <<std::endl;
	//	      std::cout << "phi   = "<< phi_F[j]  <<std::endl;
	//	      std::cout <<  std::endl;
	//	      std::cout << "#constituents:"   <<std::endl;
	//	      std::cout << "trkN = "<< trkN_I[j] <<std::endl;
	//	      std::cout << "chN  = "<<  chN_I[j] <<std::endl;
	//	      std::cout << "phN  = "<<  phN_I[j] <<std::endl;
	//	      std::cout << "eN   = "<<   eN_I[j] <<std::endl;
	//	      std::cout << "muN  = "<<  muN_I[j] <<std::endl;
	//	      std::cout << "neN  = "<<  neN_I[j] <<std::endl;
	//	      std::cout <<  std::endl;
	//	      std::cout << "#constituent fractions:"<<std::endl;
	//	      std::cout << "trkSum = "<< trkSum_F[j]/pt_F[jet] <<std::endl;
	//	      std::cout << "chSum  = "<<  chSum_F[j]/pt_F[jet] <<std::endl;
	//	      std::cout << "phSum  = "<<  phSum_F[j]/pt_F[jet] <<std::endl;
	//	      std::cout << "eSum   = "<<   eSum_F[j]/pt_F[jet] <<std::endl;
	//	      std::cout << "muSum  = "<<  muSum_F[j]/pt_F[jet] <<std::endl;
	//	      std::cout << "neSum  = "<<  neSum_F[j]/pt_F[jet] <<std::endl;
	//	      
	//	    } 
	//	    assert(false);
	//	  }
	//	  else continue;
	//	}
