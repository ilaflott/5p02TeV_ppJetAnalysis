// custom header
#include "readForests_jetTrigEff.h"

// ppData switches
const bool fillDataEvtQAHists=true;// leave me on always
const bool fillDataJetQAHists=true; // leave me on almost always
const bool useIncJetAnalyzer=true;// leave me on almost always
const bool fillDataJetTrigQAHists=true; // leave me on almost always

const bool fillDataDijetHists=false;// leave me off almost always
const bool useHLT100=false; // leave me off almost always
const bool useMBevts=false; // leave me off until Ian says to use me

//const std::string trgCombType="Calo"; // almost never used
const std::string trgCombType="PF"; // almost always used

const bool deepDebug=false; //put to true only if in trouble...
//const bool deepDebug=true; //put to true only if in trouble...

//lumi in Jet80/LowerJets PDs according to brilcalc, accounting for event selection efficiency [inv picobarns]
//const long double effJet80Lumi  =(0.9974435 * 27.385867785) ;//w/ skim cuts [PVqual, HBHENoise, beamscrape] + vzcut
//const long double effLowJetsLumi=(0.9984391 * 27.108290680) ;
//const long double effJet80Lumi  =( * 27.385867785) ;//w/ skim cuts [PVqual, HBHENoise, beamscrape] + vzcut + PFMETfracCut
//const long double effLowJetsLumi=( * 27.108290680) ;
//const long double effJet80Lumi  =(0.996080 * 27.385867785) ;//w/ skim cuts [PVqual, HBHENoise, beamscrape, HBHEIsoNoise] + vzcut
//const long double effLowJetsLumi=(0.997251 * 27.108290680) ;
const long double effJet80Lumi  =(0.993674 * 27.385867785) ;//w/ skim cuts [PVqual, HBHENoise, beamscrape, HBHEIsoNoise] + vzcut + PFMETfracCut
const long double effLowJetsLumi=(0.994089 * 27.108290680) ;
const bool useTightJetID=false;

const bool usePFCandAnalyzer=true;
const bool useTupel=false;
//// readForests_ppData_jetTrigEff
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData_jetTrigEff( std::string inFilelist , int startfile , int endfile , 
				 int radius , std::string jetType , bool debugMode ,
				 std::string outfile, float jtEtaCutLo, float jtEtaCutHi){
  
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
    if(i<7){
      trees[i]=dataTreeNames[i];
      if(debugMode){
	if(!useTupel&&i==6)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<" [NOT BEING USED THIS JOB]"<<std::endl;  
	else if(!usePFCandAnalyzer&&i==5)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<" [NOT BEING USED THIS JOB]"<<std::endl;  
	else std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
    }
    else{
      if(trgCombType=="PF"){
	trees[i]=dataTreeNames[i]+PF_HLTBitStrings[i-7]+"_v";
	if(debugMode)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  
      }
      else if(trgCombType=="Calo"){
	trees[i]=dataTreeNames[i]+Calo_HLTBitStrings[i-7]+"_v";
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
    if(!usePFCandAnalyzer && t==5) continue;
    else if(!useTupel && t==6) continue;
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
      if(!usePFCandAnalyzer && t==5) continue;
      if(!useTupel && t==6)continue;
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
  TH1D *h_NJets_jetQAPtCut  = new TH1D("NJets_jetQAPtCut ","NJets read post jettQAPtCut ", 2,0.,2.);
  TH1D *h_NJets_JetIDCut= new TH1D("NJets_JetIDCut","NJets read post JetIDCut", 2,0.,2.);





  
  
  
  //jet QA, i.e. hists filled once per jet
  TH1D *hJetQA_jtptEntries[nbins_abseta][2]={};
  TH1D *hJetQA[nbins_abseta][N_vars][2]={};    
  if(fillDataJetQAHists){        
    for(int j = 0; j<N_vars; ++j){		      
      for(int k = 0; k<nbins_abseta; ++k){      
	for(int i=0; i<2; i++){

	  if(j==0){
	    hJetQA_jtptEntries[k][i] = new TH1D( Form("hJetQA_%dwJetID_jtptEntries_etabin%d", i,k) , "N_{Jets}, no Pscl. Weights;N_{Jets};Jet p_{T} [GeV]" , 2500,0,2500);    
	    //hJetQA_jtptEntries[k][i]->Print("BASE");
	  }
	  
	  //jets
	  if(var[j]=="jtpt"||var[j]=="rawpt"){
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;", var[j].c_str()) , 2500,0,2500);       
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
 	  else if(var[j]=="trkSum"||var[j]=="phSum"||var[j]=="chSum"||var[j]=="neSum"||var[j]=="eSum"||var[j]=="muSum"||
		  var[j]=="trkHardSum"||var[j]=="phHardSum"||var[j]=="chHardSum"||
		  var[j]=="trkMax"||var[j]=="phMax"||var[j]=="chMax"||var[j]=="neMax"||var[j]=="eMax"||var[j]=="muMax"		  )
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 200,0,2);         
	  //dijets
	  else if (fillDataDijetHists){
	    std::cout<<"booking dijet hists"<<std::endl;
	    if(var[j]=="dphi") 
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;",var[j].c_str()) , 50,0,4);
	    else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt")
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;", var[j].c_str()), 2500,0,2500);
	    else //xj and Aj binnings
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 200,0,2);         
	  }//dijet hists
	}//jet ID loop
      }//nvars loop 
    }//nbins_abseta loop        
  }//fill datajetqa hists
  
  
  
  //// trigger object histograms
  
  // // these get filled 
  // // // POST EXCLUSION CRITERIA
  // // // with the max trig pt across the final trigger path for exclusion decision
  TH1D    *hpp_is40pscl =NULL ;
  TH1D    *hpp_is60pscl =NULL ;
  TH1D    *hpp_is80pscl =NULL ;

  TH1D   *hpp_is40trgPt_allobj =NULL ;
  TH1D   *hpp_is60trgPt_allobj =NULL ;
  TH1D   *hpp_is80trgPt_allobj =NULL ;  

  TH1D  *hpp_is40trgPt =NULL ,   *hpp_is40trgEta =NULL;
  TH1D  *hpp_is60trgPt =NULL ,   *hpp_is60trgEta =NULL;
  TH1D  *hpp_is80trgPt =NULL ,   *hpp_is80trgEta =NULL;
  TH1D *hpp_HLTCombtrgPt=NULL ,   *hpp_HLTCombtrgEta=NULL;

  TH1D  * hpp_is40trgPt_00eta20 =NULL , *hpp_is40trgPt_20eta51 =NULL ;  
  TH1D  * hpp_is60trgPt_00eta20 =NULL , *hpp_is60trgPt_20eta51 =NULL ;  
  TH1D  * hpp_is80trgPt_00eta20 =NULL , *hpp_is80trgPt_20eta51 =NULL ;  
  TH1D *hpp_HLTCombtrgPt_00eta20 =NULL, *hpp_HLTCombtrgPt_20eta51 =NULL ;

  TH2D   *hpp_is40trgPt_leadobj_globleadobj =NULL ;
  TH2D   *hpp_is60trgPt_leadobj_globleadobj =NULL ;
  TH2D   *hpp_is80trgPt_leadobj_globleadobj =NULL ;

  TH2D   *hpp_is40trgPt_trgEta_leadobj =NULL ;
  TH2D   *hpp_is60trgPt_trgEta_leadobj =NULL ;
  TH2D   *hpp_is80trgPt_trgEta_leadobj =NULL ;

  TH1D   *hpp_is40trgPt_leadobj_finalSpectra =NULL ,  *hpp_is40trgEta_leadobj_finalSpectra =NULL ;
  TH1D   *hpp_is60trgPt_leadobj_finalSpectra =NULL ,  *hpp_is60trgEta_leadobj_finalSpectra =NULL ;
  TH1D   *hpp_is80trgPt_leadobj_finalSpectra =NULL ,  *hpp_is80trgEta_leadobj_finalSpectra =NULL ;
  
  // // these get fillied
  // // // post skips but pre-cuts
  // // // if the given HLT path fired 
  // // // with the given path's object w/ highest pt(s)/eta/prescale etc.
  // // // weighted by the given HLT objects prescale
  TH1D   *hpp_IncHLT40pscl =NULL ;
  TH1D   *hpp_IncHLT60pscl =NULL ;
  TH1D   *hpp_IncHLT80pscl =NULL ;
  
  TH1D   *hpp_IncHLT40trgPt_allobj =NULL ;
  TH1D   *hpp_IncHLT60trgPt_allobj =NULL ;
  TH1D   *hpp_IncHLT80trgPt_allobj =NULL ;

  TH1D   *hpp_IncHLT40trgPt_leadobj =NULL, *hpp_IncHLT40trgEta =NULL ;
  TH1D   *hpp_IncHLT60trgPt_leadobj =NULL, *hpp_IncHLT60trgEta =NULL ;
  TH1D   *hpp_IncHLT80trgPt_leadobj =NULL, *hpp_IncHLT80trgEta =NULL ;
  //
  
  TH1D   *hpp_IncHLT40trgPt_leadobj_00eta20 =NULL , *hpp_IncHLT40trgPt_leadobj_20eta51 =NULL ;
  TH1D   *hpp_IncHLT60trgPt_leadobj_00eta20 =NULL , *hpp_IncHLT60trgPt_leadobj_20eta51 =NULL ;
  TH1D   *hpp_IncHLT80trgPt_leadobj_00eta20 =NULL , *hpp_IncHLT80trgPt_leadobj_20eta51 =NULL ;
  //

  TH2D   *hpp_IncHLT40trgPt_leadobj_globleadobj =NULL ;
  TH2D   *hpp_IncHLT60trgPt_leadobj_globleadobj =NULL ;
  TH2D   *hpp_IncHLT80trgPt_leadobj_globleadobj =NULL ;

  TH2D   *hpp_IncHLT40trgPt_trgEta_leadobj =NULL ;
  TH2D   *hpp_IncHLT60trgPt_trgEta_leadobj =NULL ;
  TH2D   *hpp_IncHLT80trgPt_trgEta_leadobj =NULL ;
  
  TH1D   *hpp_IncHLT40trgPt_leadobj_finalSpectra =NULL ,  *hpp_IncHLT40trgEta_leadobj_finalSpectra =NULL ;
  TH1D   *hpp_IncHLT60trgPt_leadobj_finalSpectra =NULL ,  *hpp_IncHLT60trgEta_leadobj_finalSpectra =NULL ;
  TH1D   *hpp_IncHLT80trgPt_leadobj_finalSpectra =NULL ,  *hpp_IncHLT80trgEta_leadobj_finalSpectra =NULL ;


  //TRIGGER EFFICIENCY HISTS
  //for measuring trigger eff. of HLT80 --------
  TH1D* hpp_HLT80ref_IncHLT60trgPt_allobj =NULL;//ref method
  TH1D* hpp_HLT80ref_IncHLT60trgPt_leadobj=NULL;
  TH1D* hpp_HLT80ref_IncHLT60trgPt_leadobj_00eta20=NULL;
  TH1D* hpp_HLT80ref_IncHLT60trgPt_leadobj_20eta51=NULL;
  TH1D* hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj=NULL;//emulation method (basically ref + L1 req)
  TH1D* hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_00eta20=NULL;
  TH1D* hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_20eta51=NULL;

  //for measuring trigger eff. of HLT60 --------
  TH1D* hpp_HLT60ref_IncHLT40trgPt_allobj =NULL;//ref method
  TH1D* hpp_HLT60ref_IncHLT40trgPt_leadobj=NULL;
  TH1D* hpp_HLT60ref_IncHLT40trgPt_leadobj_00eta20=NULL;
  TH1D* hpp_HLT60ref_IncHLT40trgPt_leadobj_20eta51=NULL;
  TH1D* hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj=NULL;//emulation method (basically ref + L1 req)
  TH1D* hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_00eta20=NULL;
  TH1D* hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_20eta51=NULL;
  //for measuring trigger eff. of HLT40 --------
  //NULL? USE MB SOMEHOW?



  
  if(fillDataJetTrigQAHists){

    //FULL EXCLUSION CRITERIA APPLIED
    hpp_is40pscl   =new TH1D( "is40pscl",    "HLT40 pscl;Prescale;N_{evts}" , 200,0,200);
    hpp_is60pscl   =new TH1D( "is60pscl",    "HLT60 pscl;Prescale;N_{evts}" , 200,0,200);
    hpp_is80pscl   =new TH1D( "is80pscl",    "HLT80 pscl;Prescale;N_{evts}" , 200,0,200);

    hpp_is40trgPt_allobj   =new TH1D( "is40trgPt_allobj",    "all trig jet object for HLT40;Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_is60trgPt_allobj   =new TH1D( "is60trgPt_allobj",    "all trig jet object for HLT60;Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_is80trgPt_allobj   =new TH1D( "is80trgPt_allobj",    "all trig jet object for HLT80;Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    
    hpp_is40trgPt   = new TH1D( "is40_trgPt"  , "trgPt for is40   ", 2000,0,2000 ); 
    hpp_is60trgPt   = new TH1D( "is60_trgPt"  , "trgPt for is60   ", 2000,0,2000 ); 
    hpp_is80trgPt   = new TH1D( "is80_trgPt"  , "trgPt for is80   ", 2000,0,2000 ); 
    hpp_HLTCombtrgPt = new TH1D( "HLTComb_trgPt"  , "trgPt for HLTComb ", 2000,0,2000 );       
    
    hpp_is40trgEta   = new TH1D( "is40_trgEta"  , "exc trgEta for HLT40   ", 102, -5.1, 5.1 ); 
    hpp_is60trgEta   = new TH1D( "is60_trgEta"  , "exc trgEta for HLT60   ", 102, -5.1, 5.1 ); 
    hpp_is80trgEta   = new TH1D( "is80_trgEta"  , "exc trgEta for HLT80   ", 102, -5.1, 5.1 ); 
    hpp_HLTCombtrgEta = new TH1D( "HLTComb_trgEta"  , "trgEta for HLTComb ", 102, -5.1, 5.1 );           

    hpp_is40trgPt_00eta20  = new TH1D(  "is40_trgPt_00eta20"  , "trgPt for is40  #||{y}<2.0 ", 2000,0,2000 ); 
    hpp_is60trgPt_00eta20  = new TH1D(  "is60_trgPt_00eta20"  , "trgPt for is60  #||{y}<2.0 ", 2000,0,2000 ); 
    hpp_is80trgPt_00eta20  = new TH1D(  "is80_trgPt_00eta20"  , "trgPt for is80  #||{y}<2.0 ", 2000,0,2000 ); 
    hpp_HLTCombtrgPt_00eta20 = new TH1D( "HLTComb_trgPt_00eta20"  , "trgPt for HLTComb  #||{y}<2.0", 2000,0,2000 );       
    
    hpp_is40trgPt_20eta51  = new TH1D(  "is40_trgPt_20eta51"  , "trgPt for is40  #||{y}>2.0 ", 2000,0,2000 ); 
    hpp_is60trgPt_20eta51  = new TH1D(  "is60_trgPt_20eta51"  , "trgPt for is60  #||{y}>2.0 ", 2000,0,2000 ); 
    hpp_is80trgPt_20eta51  = new TH1D(  "is80_trgPt_20eta51"  , "trgPt for is80  #||{y}>2.0 ", 2000,0,2000 ); 
    hpp_HLTCombtrgPt_20eta51 = new TH1D( "HLTComb_trgPt_20eta51"  , "trgPt for HLTComb  #||{y}>2.0", 2000,0,2000 );       


    hpp_is40trgPt_leadobj_globleadobj  =new TH2D(  "is40trgPt_leadobj_globleadobj", "lead trig jet obj across HLT40 v lead trig jet obj across all triggers, only for HLT40;Lead Trig Obj p_{T} across HLT40;Lead Obj p_{T} across HLT40+60+80" , 
						       2000,0,2000, 2000,0,2000);
    hpp_is60trgPt_leadobj_globleadobj  =new TH2D(  "is60trgPt_leadobj_globleadobj", "lead trig jet obj across HLT60 v lead trig jet obj across all triggers, only for HLT60;Lead Trig Obj p_{T} across HLT60;Lead Obj p_{T} across HLT40+60+80" , 
						       2000,0,2000, 2000,0,2000);
    hpp_is80trgPt_leadobj_globleadobj  =new TH2D(  "is80trgPt_leadobj_globleadobj", "lead trig jet obj across HLT80 v lead trig jet obj across all triggers, only for HLT80;Lead Trig Obj p_{T} across HLT80;Lead Obj p_{T} across HLT40+60+80" , 
						       2000,0,2000, 2000,0,2000);


    hpp_is40trgPt_trgEta_leadobj  =new TH2D(  "is40trgPt_trgEta_leadobj", "lead trig jet pt HLT40 v lead trig jet y, only for HLT40;Lead Trig Obj p_{T} across HLT40;Lead Obj y across HLT40" , 
					      2000,0,2000, 102,-5.1,5.1);
    hpp_is60trgPt_trgEta_leadobj  =new TH2D(  "is60trgPt_trgEta_leadobj", "lead trig jet pt HLT60 v lead trig jet y, only for HLT60;Lead Trig Obj p_{T} across HLT60;Lead Obj y across HLT60" , 
					      2000,0,2000, 102,-5.1,5.1);
    hpp_is80trgPt_trgEta_leadobj  =new TH2D(  "is80trgPt_trgEta_leadobj", "lead trig jet pt HLT80 v lead trig jet y, only for HLT80;Lead Trig Obj p_{T} across HLT80;Lead Obj y across HLT80" , 
					      2000,0,2000, 102,-5.1,5.1);    
    
    hpp_is40trgPt_leadobj_finalSpectra    =new TH1D(  "is40trgPt_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_is60trgPt_leadobj_finalSpectra    =new TH1D(  "is60trgPt_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_is80trgPt_leadobj_finalSpectra    =new TH1D(  "is80trgPt_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT80;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);

    hpp_is40trgEta_leadobj_finalSpectra   =new TH1D(  "is40trgEta_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT40;Lead Trig Jet Obj y;N_{jets}" , 102,-5.1,5.1);
    hpp_is60trgEta_leadobj_finalSpectra   =new TH1D(  "is60trgEta_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT60;Lead Trig Jet Obj y;N_{jets}" , 102,-5.1,5.1);
    hpp_is80trgEta_leadobj_finalSpectra   =new TH1D(  "is80trgEta_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT80;Lead Trig Jet Obj y;N_{jets}" , 102,-5.1,5.1);


    



    //NO EXCLUSION CRITERIA APPLIED    
    hpp_IncHLT40pscl   =new TH1D( "IncHLT40pscl",    "HLT40 pscl;Prescale;N_{evts}" , 200,0,200);
    hpp_IncHLT60pscl   =new TH1D( "IncHLT60pscl",    "HLT60 pscl;Prescale;N_{evts}" , 200,0,200);
    hpp_IncHLT80pscl   =new TH1D( "IncHLT80pscl",    "HLT80 pscl;Prescale;N_{evts}" , 200,0,200);

    hpp_IncHLT40trgPt_allobj   =new TH1D( "IncHLT40trgPt_allobj",    "all trig jet object for HLT40;Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT60trgPt_allobj   =new TH1D( "IncHLT60trgPt_allobj",    "all trig jet object for HLT60;Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT80trgPt_allobj   =new TH1D( "IncHLT80trgPt_allobj",    "all trig jet object for HLT80;Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);

    hpp_IncHLT40trgEta   = new TH1D( "IncHLT40_trgEta"  , "inc trgEta for HLT40   ", 102, -5.1, 5.1 ); 
    hpp_IncHLT60trgEta   = new TH1D( "IncHLT60_trgEta"  , "inc trgEta for HLT60   ", 102, -5.1, 5.1 ); 
    hpp_IncHLT80trgEta   = new TH1D( "IncHLT80_trgEta"  , "inc trgEta for HLT80   ", 102, -5.1, 5.1 ); 
    
    hpp_IncHLT40trgPt_leadobj   =new TH1D(  "IncHLT40trgPt_leadobj",  "lead trig jet object for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT60trgPt_leadobj   =new TH1D(  "IncHLT60trgPt_leadobj",  "lead trig jet object for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT80trgPt_leadobj   =new TH1D(  "IncHLT80trgPt_leadobj",  "lead trig jet object for HLT80;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);

    hpp_IncHLT40trgPt_leadobj_00eta20   =new TH1D(  "IncHLT40trgPt_leadobj_00eta20",  "lead trig jet obj. for HLT40, #||{y}<2.0;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT60trgPt_leadobj_00eta20   =new TH1D(  "IncHLT60trgPt_leadobj_00eta20",  "lead trig jet obj. for HLT60, #||{y}<2.0;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT80trgPt_leadobj_00eta20   =new TH1D(  "IncHLT80trgPt_leadobj_00eta20",  "lead trig jet obj. for HLT80, #||{y}<2.0;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    
    hpp_IncHLT40trgPt_leadobj_20eta51   =new TH1D(  "IncHLT40trgPt_leadobj_20eta51",  "lead trig jet obj. for HLT40, 2.0<#||{y}<5.1;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT60trgPt_leadobj_20eta51   =new TH1D(  "IncHLT60trgPt_leadobj_20eta51",  "lead trig jet obj. for HLT60, 2.0<#||{y}<5.1;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT80trgPt_leadobj_20eta51   =new TH1D(  "IncHLT80trgPt_leadobj_20eta51",  "lead trig jet obj. for HLT80, 2.0<#||{y}<5.1;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    

    hpp_IncHLT40trgPt_leadobj_globleadobj  =new TH2D(  "IncHLT40trgPt_leadobj_globleadobj", "lead trig jet obj across HLT40 v lead trig jet obj across all triggers, only for HLT40;Lead Trig Obj p_{T} across HLT40;Lead Obj p_{T} across HLT40+60+80" , 
						       2000,0,2000, 2000,0,2000);
    hpp_IncHLT60trgPt_leadobj_globleadobj  =new TH2D(  "IncHLT60trgPt_leadobj_globleadobj", "lead trig jet obj across HLT60 v lead trig jet obj across all triggers, only for HLT60;Lead Trig Obj p_{T} across HLT60;Lead Obj p_{T} across HLT40+60+80" , 
						       2000,0,2000, 2000,0,2000);
    hpp_IncHLT80trgPt_leadobj_globleadobj  =new TH2D(  "IncHLT80trgPt_leadobj_globleadobj", "lead trig jet obj across HLT80 v lead trig jet obj across all triggers, only for HLT80;Lead Trig Obj p_{T} across HLT80;Lead Obj p_{T} across HLT40+60+80" , 
						       2000,0,2000, 2000,0,2000);

    hpp_IncHLT40trgPt_trgEta_leadobj  =new TH2D(  "IncHLT40trgPt_trgEta_leadobj", "lead trig jet pt HLT40 v lead trig jet y, only for HLT40;Lead Trig Obj p_{T} across HLT40;Lead Obj y across HLT40" , 
						  2000,0,2000, 102,-5.1,5.1);
    hpp_IncHLT60trgPt_trgEta_leadobj  =new TH2D(  "IncHLT60trgPt_trgEta_leadobj", "lead trig jet pt HLT60 v lead trig jet y, only for HLT60;Lead Trig Obj p_{T} across HLT60;Lead Obj y across HLT60" , 
						  2000,0,2000, 102,-5.1,5.1);
    hpp_IncHLT80trgPt_trgEta_leadobj  =new TH2D(  "IncHLT80trgPt_trgEta_leadobj", "lead trig jet pt HLT80 v lead trig jet y, only for HLT80;Lead Trig Obj p_{T} across HLT80;Lead Obj y across HLT80" , 
						  2000,0,2000, 102,-5.1,5.1);


    hpp_IncHLT40trgPt_leadobj_finalSpectra   =new TH1D(  "IncHLT40trgPt_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT60trgPt_leadobj_finalSpectra   =new TH1D(  "IncHLT60trgPt_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_IncHLT80trgPt_leadobj_finalSpectra   =new TH1D(  "IncHLT80trgPt_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT80;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);

    hpp_IncHLT40trgEta_leadobj_finalSpectra   =new TH1D(  "IncHLT40trgEta_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT40;Lead Trig Jet Obj y;N_{jets}" , 102,-5.1,5.1);
    hpp_IncHLT60trgEta_leadobj_finalSpectra   =new TH1D(  "IncHLT60trgEta_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT60;Lead Trig Jet Obj y;N_{jets}" , 102,-5.1,5.1);
    hpp_IncHLT80trgEta_leadobj_finalSpectra   =new TH1D(  "IncHLT80trgEta_leadobj_finalSpectra",  "in final spectra, lead trig jet obj. for HLT80;Lead Trig Jet Obj y;N_{jets}" , 102,-5.1,5.1);    





    
    
    //TRIGGER EFFICIENCY HISTS
    //for measuring trigger eff. of HLT80 --------
    hpp_HLT80ref_IncHLT60trgPt_allobj   =new TH1D( "HLT80ref_IncHLT60trgPt_allobj" ,   "HLT80 Ref., all trig jet object for HLT60;Trig Jet Obj p_{T};N_{jets}"      , 2000,0,2000);   
    hpp_HLT80ref_IncHLT60trgPt_leadobj  =new TH1D( "HLT80ref_IncHLT60trgPt_leadobj",   "HLT80 Ref., lead trig jet object for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT80ref_IncHLT60trgPt_leadobj_00eta20  =new TH1D( "HLT80ref_IncHLT60trgPt_leadobj_00eta20",   "HLT80 Ref., lead trig jet object for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT80ref_IncHLT60trgPt_leadobj_20eta51  =new TH1D( "HLT80ref_IncHLT60trgPt_leadobj_20eta51",   "HLT80 Ref., lead trig jet object for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj  =new TH1D( "HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj",   "HLT80 Ref., lead trig jet object for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_00eta20  =new TH1D( "HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_00eta20",   "HLT80 Ref., lead trig jet object for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_20eta51  =new TH1D( "HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_20eta51",   "HLT80 Ref., lead trig jet object for HLT60;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);    
    
    //for measuring trigger eff. of HLT60 --------
    hpp_HLT60ref_IncHLT40trgPt_allobj   =new TH1D( "HLT60ref_IncHLT40trgPt_allobj" ,   "HLT60 Ref., all trig jet object for HLT40;Trig Jet Obj p_{T};N_{jets}"      , 2000,0,2000);
    hpp_HLT60ref_IncHLT40trgPt_leadobj  =new TH1D( "HLT60ref_IncHLT40trgPt_leadobj",   "HLT60 Ref., lead trig jet object for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT60ref_IncHLT40trgPt_leadobj_00eta20  =new TH1D( "HLT60ref_IncHLT40trgPt_leadobj_00eta20",   "HLT60 Ref., lead trig jet object for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT60ref_IncHLT40trgPt_leadobj_20eta51  =new TH1D( "HLT60ref_IncHLT40trgPt_leadobj_20eta51",   "HLT60 Ref., lead trig jet object for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj  =new TH1D( "HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj",   "HLT60 Ref., lead trig jet object for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_00eta20  =new TH1D( "HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_00eta20",   "HLT60 Ref., lead trig jet object for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_20eta51  =new TH1D( "HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_20eta51",   "HLT60 Ref., lead trig jet object for HLT40;Lead Trig Jet Obj p_{T};N_{jets}" , 2000,0,2000);
    
    
  }
  
  
  
  
  
  
  //inclusive jet spectra by-trigger histograms
  TH1D  *hpp_HLT40[nbins_abseta][2]={}, *hpp_HLT60[nbins_abseta][2]={}, *hpp_HLT80[nbins_abseta][2]={}, *hpp_HLT100[nbins_abseta][2]={}; //no exclusion, looks at each path on its own terms
  TH1D *hpp_HLTMB[nbins_abseta][2]={};
  
  TH1D  *hpp_incHLT40[nbins_abseta][2]={}, *hpp_incHLT60[nbins_abseta][2]={}, *hpp_incHLT80[nbins_abseta][2]={}, *hpp_incHLT100[nbins_abseta][2]={};// semi excluded, using max trig pt between all objects and looks at lower threshold
  TH1D *hpp_incHLTMB[nbins_abseta][2]={};
  
  TH1D  *hpp_excHLT40[nbins_abseta][2]={}, *hpp_excHLT60[nbins_abseta][2]={}, *hpp_excHLT80[nbins_abseta][2]={}, *hpp_excHLT100[nbins_abseta][2]={};// full exclusion, using max trig pt between all objects and looking at all thresholds
  TH1D* hpp_excHLTMB[nbins_abseta][2]={};  
  TH1D *hpp_HLTComb[nbins_abseta][2]={};//combination of triggers w/ exclusion    



  //TRIGGER EFFICIENCY HISTS
  //for measuring trigger eff. of HLT80 --------  
  TH1D* hpp_HLT80ref_HLT60[nbins_abseta]={};
  TH1D* hpp_HLT80ref_HLT60_wL1Thresh[nbins_abseta]={};
  //TH1D* hpp_HLT80ref_HLT40[nbins_abseta]={};  
  //TH1D* hpp_HLT80ref_HLT40_wL1Thresh[nbins_abseta]={};  
  //for measuring trigger eff. of HLT60 --------
  TH1D* hpp_HLT60ref_HLT40[nbins_abseta]={};  
  TH1D* hpp_HLT60ref_HLT40_wL1Thresh[nbins_abseta]={};  
  //for measuring trigger eff. of HLT40 --------
  //NULL


  //special histograms

  //|y|<0.5 only, wJetID only
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
	    if(i==1)hHLTTitle+="_wJetID";
	    hHLTTitle+="_R"+std::to_string(radius)+"_etabin"+std::to_string(j);//etaWidth;
	    //if(debugMode)std::cout<<" exclusive HLT plothExcTitle="<<hExcTitle<<std::endl;	  
	    
	    std::string hHLTDesc=hJetHLTQADescArray[k];
	    if(i==1)hHLTDesc+=" w/ JetID";
	    hHLTDesc+=", "+absetabins_str[j] +"< #||{y} <"+absetabins_str[j+1];	    
	    //if(debugMode)std::cout<<" exclusive HLT plothExcDesc="<<hExcDesc<<std::endl<<std::endl;	    
	    
	    if     (hJetHLTQATitleArray[k]=="hpp_HLTMB"  )  hpp_HLTMB[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT40"  )  hpp_HLT40[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT60"  )  hpp_HLT60[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT80"  )  hpp_HLT80[j][i]     = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    else if(hJetHLTQATitleArray[k]=="hpp_HLT100" )  hpp_HLT100[j][i]    = new TH1D( hHLTTitle.c_str(), (hHLTDesc).c_str(), 2500,0,2500);
	    
	    if(i==1){
	      //HLT80 REF TRIGGERS
	      if(hJetHLTQATitleArray[k]=="hpp_HLT60"  )  hpp_HLT80ref_HLT60[j]     = new TH1D( ("hpp_HLT80ref_HLT60_wJetID_R"+
												std::to_string(radius)+"_etabin"+std::to_string(j)).c_str(),
											       ("HLT80 Ref., "+hHLTDesc).c_str(), 2500,0,2500);	
	      if(hJetHLTQATitleArray[k]=="hpp_HLT60"  )  hpp_HLT80ref_HLT60_wL1Thresh[j]     = new TH1D( ("hpp_HLT80ref_HLT60_wL1Thresh_wJetID_R"+
													  std::to_string(radius)+"_etabin"+std::to_string(j)).c_str(),
													 ("HLT80 Ref. && LeadL1JetEt>48, "+hHLTDesc).c_str(), 2500,0,2500);	
	      //if(hJetHLTQATitleArray[k]=="hpp_HLT40"  )  hpp_HLT80ref_HLT40[j]     = new TH1D( ("hpp_HLT80ref_HLT40_wJetID_R"+
	      //											std::to_string(radius)+"_etabin"+std::to_string(j)).c_str(), 
	      //										       ("HLT80 Ref., "+hHLTDesc).c_str(), 2500,0,2500);
	      //if(hJetHLTQATitleArray[k]=="hpp_HLT40"  )  hpp_HLT80ref_HLT40_wL1Thresh[j]     = new TH1D( ("hpp_HLT80ref_HLT40_wL1Thresh_wJetID_R"+
	      //												  std::to_string(radius)+"_etabin"+std::to_string(j)).c_str(), 
	      //												 ("HLT80 Ref. && LeadL1JetEt>48, "+hHLTDesc).c_str(), 2500,0,2500);
	      //HLT60 REF TRIGGERS
	      if(hJetHLTQATitleArray[k]=="hpp_HLT40"  )  hpp_HLT60ref_HLT40[j]     = new TH1D( ("hpp_HLT60ref_HLT40_wJetID_R"+
												std::to_string(radius)+"_etabin"+std::to_string(j)).c_str(), 
											       ("HLT60 Ref., "+hHLTDesc).c_str(), 2500,0,2500);
	      if(hJetHLTQATitleArray[k]=="hpp_HLT40"  )  hpp_HLT60ref_HLT40_wL1Thresh[j]     = new TH1D( ("hpp_HLT60ref_HLT40_wL1Thresh_wJetID_R"+
													  std::to_string(radius)+"_etabin"+std::to_string(j)).c_str(), 
													 ("HLT60 Ref. && LeadL1JetEt>40, "+hHLTDesc).c_str(), 2500,0,2500);
	    }
	    
	    
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
  
  
  
  
  TH1D* hpp_NL1ForJets=NULL, *hpp_NL1CenJets=NULL, *hpp_NL1Jets=NULL;
  TH1D* hpp_LeadL1ForJet_Et=NULL, *hpp_LeadL1CenJet_Et=NULL, *hpp_LeadL1Jet_Et=NULL;
  TH1D* hpp_L1ForJet_E=NULL, *hpp_L1ForJet_Et=NULL, *hpp_L1ForJet_Eta=NULL, *hpp_L1ForJet_Phi=NULL;
  TH1D* hpp_L1CenJet_E=NULL, *hpp_L1CenJet_Et=NULL, *hpp_L1CenJet_Eta=NULL, *hpp_L1CenJet_Phi=NULL;
  TH1D* hpp_L1Jet_E=NULL,    *hpp_L1Jet_Et=NULL, *hpp_L1Jet_Eta=NULL, *hpp_L1Jet_Phi=NULL;

  TH1D* hpp_SJ48_NL1ForJets=NULL, *hpp_SJ48_NL1CenJets=NULL, *hpp_SJ48_NL1Jets=NULL;
  TH1D* hpp_SJ48_LeadL1ForJet_Et=NULL, *hpp_SJ48_LeadL1CenJet_Et=NULL, *hpp_SJ48_LeadL1Jet_Et=NULL;
  TH1D* hpp_SJ48_L1ForJet_E=NULL, *hpp_SJ48_L1ForJet_Et=NULL, *hpp_SJ48_L1ForJet_Eta=NULL, *hpp_SJ48_L1ForJet_Phi=NULL;
  TH1D* hpp_SJ48_L1CenJet_E=NULL, *hpp_SJ48_L1CenJet_Et=NULL, *hpp_SJ48_L1CenJet_Eta=NULL, *hpp_SJ48_L1CenJet_Phi=NULL;
  TH1D* hpp_SJ48_L1Jet_E=NULL,    *hpp_SJ48_L1Jet_Et=NULL, *hpp_SJ48_L1Jet_Eta=NULL, *hpp_SJ48_L1Jet_Phi=NULL;
  
  TH1D* hpp_SJ40_NL1ForJets=NULL, *hpp_SJ40_NL1CenJets=NULL, *hpp_SJ40_NL1Jets=NULL;
  TH1D* hpp_SJ40_LeadL1ForJet_Et=NULL, *hpp_SJ40_LeadL1CenJet_Et=NULL, *hpp_SJ40_LeadL1Jet_Et=NULL;
  TH1D* hpp_SJ40_L1ForJet_E=NULL, *hpp_SJ40_L1ForJet_Et=NULL, *hpp_SJ40_L1ForJet_Eta=NULL, *hpp_SJ40_L1ForJet_Phi=NULL;
  TH1D* hpp_SJ40_L1CenJet_E=NULL, *hpp_SJ40_L1CenJet_Et=NULL, *hpp_SJ40_L1CenJet_Eta=NULL, *hpp_SJ40_L1CenJet_Phi=NULL;
  TH1D* hpp_SJ40_L1Jet_E=NULL,    *hpp_SJ40_L1Jet_Et=NULL, *hpp_SJ40_L1Jet_Eta=NULL, *hpp_SJ40_L1Jet_Phi=NULL;
  
  TH1D* hpp_SJ28_NL1ForJets=NULL, *hpp_SJ28_NL1CenJets=NULL, *hpp_SJ28_NL1Jets=NULL;
  TH1D* hpp_SJ28_LeadL1ForJet_Et=NULL, *hpp_SJ28_LeadL1CenJet_Et=NULL, *hpp_SJ28_LeadL1Jet_Et=NULL;
  TH1D* hpp_SJ28_L1ForJet_E=NULL, *hpp_SJ28_L1ForJet_Et=NULL, *hpp_SJ28_L1ForJet_Eta=NULL, *hpp_SJ28_L1ForJet_Phi=NULL;
  TH1D* hpp_SJ28_L1CenJet_E=NULL, *hpp_SJ28_L1CenJet_Et=NULL, *hpp_SJ28_L1CenJet_Eta=NULL, *hpp_SJ28_L1CenJet_Phi=NULL;
  TH1D* hpp_SJ28_L1Jet_E=NULL,    *hpp_SJ28_L1Jet_Et=NULL, *hpp_SJ28_L1Jet_Eta=NULL, *hpp_SJ28_L1Jet_Phi=NULL;
  
  if(fillDataJetTrigQAHists){
    
    //L1 Jet Counts
    hpp_NL1ForJets=new TH1D("NL1ForJets", "N L1 Forward Jets", 8,0,8);
    hpp_NL1CenJets=new TH1D("NL1CenJets", "N L1 Central Jets", 8,0,8);
    hpp_NL1Jets   =new TH1D("NL1Jets"   , "N L1 Jets"        , 8,0,8);
    //Lead L1 Jets
    hpp_LeadL1ForJet_Et=new TH1D("LeadL1ForJet_Et", "Lead L1 Forward Jet E_{T}", 1000,0.,1000.);
    hpp_LeadL1CenJet_Et=new TH1D("LeadL1CenJet_Et", "Lead L1 Central Jet E_{T}", 1000,0.,1000.);
    hpp_LeadL1Jet_Et   =new TH1D("LeadL1Jet_Et"   , "Lead L1 Jet E_{T}"        , 1000,0.,1000.);
    //all L1 Jets
    hpp_L1Jet_E  =new TH1D("L1Jet_E"   , "L1 Jet E"    , 1000, 0., 1000.);
    hpp_L1Jet_Et =new TH1D("L1Jet_Et"  , "L1 Jet E_{T}", 1000, 0., 1000.);
    hpp_L1Jet_Eta=new TH1D("L1Jet_Eta" , "L1 Jet #eta" , 102 , -5.1, 5.1);
    hpp_L1Jet_Phi=new TH1D("L1Jet_Phi" , "L1 Jet #phi" , 64  , -3.2,3.2 );
    //Forward L1 Jets
    hpp_L1ForJet_E  =new TH1D("L1ForJet_E"   , "L1 Forward Jet E"    , 1000, 0., 1000.);
    hpp_L1ForJet_Et =new TH1D("L1ForJet_Et"  , "L1 Forward Jet E_{T}", 1000, 0., 1000.);
    hpp_L1ForJet_Eta=new TH1D("L1ForJet_Eta" , "L1 Forward Jet #eta" , 102 , -5.1, 5.1);
    hpp_L1ForJet_Phi=new TH1D("L1ForJet_Phi" , "L1 Forward Jet #phi" , 64  , -3.2,3.2 );
    //Central Jets
    hpp_L1CenJet_E  =new TH1D("L1CenJet_E"   , "L1 Central Jet E"    , 1000, 0., 1000.);
    hpp_L1CenJet_Et =new TH1D("L1CenJet_Et"  , "L1 Central Jet E_{T}", 1000, 0., 1000.);
    hpp_L1CenJet_Eta=new TH1D("L1CenJet_Eta" , "L1 Central Jet #eta" , 102 , -5.1, 5.1);
    hpp_L1CenJet_Phi=new TH1D("L1CenJet_Phi" , "L1 Central Jet #phi" , 64  , -3.2,3.2 );
    
    //by L1 Trigger Path Firings; L1SingleJet48_BptxAND
    //L1 Jet Counts
    hpp_SJ48_NL1ForJets=new TH1D("SJ48_NL1ForJets", "L1_SingleJet48==1, N L1 Forward Jets", 8,0,8);
    hpp_SJ48_NL1CenJets=new TH1D("SJ48_NL1CenJets", "L1_SingleJet48==1, N L1 Central Jets", 8,0,8);
    hpp_SJ48_NL1Jets   =new TH1D("SJ48_NL1Jets"   , "L1_SingleJet48==1, N L1 Jets"        , 8,0,8);
    //Lead L1 Jets
    hpp_SJ48_LeadL1ForJet_Et=new TH1D("SJ48_LeadL1ForJet_Et", "L1_SingleJet48==1, Lead L1 Forward Jet E_{T}", 1000,0.,1000.);
    hpp_SJ48_LeadL1CenJet_Et=new TH1D("SJ48_LeadL1CenJet_Et", "L1_SingleJet48==1, Lead L1 Central Jet E_{T}", 1000,0.,1000.);
    hpp_SJ48_LeadL1Jet_Et   =new TH1D("SJ48_LeadL1Jet_Et"   , "L1_SingleJet48==1, Lead L1 Jet E_{T}"        , 1000,0.,1000.);
    //all L1 Jets
    hpp_SJ48_L1Jet_E  =new TH1D("SJ48_L1Jet_E"   , "L1_SingleJet48==1, L1 Jet E"    , 1000, 0., 1000.);
    hpp_SJ48_L1Jet_Et =new TH1D("SJ48_L1Jet_Et"  , "L1_SingleJet48==1, L1 Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ48_L1Jet_Eta=new TH1D("SJ48_L1Jet_Eta" , "L1_SingleJet48==1, L1 Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ48_L1Jet_Phi=new TH1D("SJ48_L1Jet_Phi" , "L1_SingleJet48==1, L1 Jet #phi" , 64  , -3.2,3.2 );
    //Forward L1 Jets
    hpp_SJ48_L1ForJet_E  =new TH1D("SJ48_L1ForJet_E"   , "L1_SingleJet48==1, L1 Forward Jet E"    , 1000, 0., 1000.);
    hpp_SJ48_L1ForJet_Et =new TH1D("SJ48_L1ForJet_Et"  , "L1_SingleJet48==1, L1 Forward Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ48_L1ForJet_Eta=new TH1D("SJ48_L1ForJet_Eta" , "L1_SingleJet48==1, L1 Forward Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ48_L1ForJet_Phi=new TH1D("SJ48_L1ForJet_Phi" , "L1_SingleJet48==1, L1 Forward Jet #phi" , 64  , -3.2,3.2 );
    //Central Jets
    hpp_SJ48_L1CenJet_E  =new TH1D("SJ48_L1CenJet_E"   , "L1_SingleJet48==1, L1 Central Jet E"    , 1000, 0., 1000.);
    hpp_SJ48_L1CenJet_Et =new TH1D("SJ48_L1CenJet_Et"  , "L1_SingleJet48==1, L1 Central Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ48_L1CenJet_Eta=new TH1D("SJ48_L1CenJet_Eta" , "L1_SingleJet48==1, L1 Central Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ48_L1CenJet_Phi=new TH1D("SJ48_L1CenJet_Phi" , "L1_SingleJet48==1, L1 Central Jet #phi" , 64  , -3.2,3.2 );


    //by L1 Trigger Path Firings; L1SingleJet40_BptxAND
    //L1 Jet Counts
    hpp_SJ40_NL1ForJets=new TH1D("SJ40_NL1ForJets", "L1_SingleJet40==1, N L1 Forward Jets", 8,0,8);
    hpp_SJ40_NL1CenJets=new TH1D("SJ40_NL1CenJets", "L1_SingleJet40==1, N L1 Central Jets", 8,0,8);
    hpp_SJ40_NL1Jets   =new TH1D("SJ40_NL1Jets"   , "L1_SingleJet40==1, N L1 Jets"        , 8,0,8);
    //Lead L1 Jets
    hpp_SJ40_LeadL1ForJet_Et=new TH1D("SJ40_LeadL1ForJet_Et", "L1_SingleJet40==1, Lead L1 Forward Jet E_{T}", 1000,0.,1000.);
    hpp_SJ40_LeadL1CenJet_Et=new TH1D("SJ40_LeadL1CenJet_Et", "L1_SingleJet40==1, Lead L1 Central Jet E_{T}", 1000,0.,1000.);
    hpp_SJ40_LeadL1Jet_Et   =new TH1D("SJ40_LeadL1Jet_Et"   , "L1_SingleJet40==1, Lead L1 Jet E_{T}"        , 1000,0.,1000.);
    //all L1 Jets
    hpp_SJ40_L1Jet_E  =new TH1D("SJ40_L1Jet_E"   , "L1_SingleJet40==1, L1 Jet E"    , 1000, 0., 1000.);
    hpp_SJ40_L1Jet_Et =new TH1D("SJ40_L1Jet_Et"  , "L1_SingleJet40==1, L1 Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ40_L1Jet_Eta=new TH1D("SJ40_L1Jet_Eta" , "L1_SingleJet40==1, L1 Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ40_L1Jet_Phi=new TH1D("SJ40_L1Jet_Phi" , "L1_SingleJet40==1, L1 Jet #phi" , 64  , -3.2,3.2 );
    //Forward L1 Jets
    hpp_SJ40_L1ForJet_E  =new TH1D("SJ40_L1ForJet_E"   , "L1_SingleJet40==1, L1 Forward Jet E"    , 1000, 0., 1000.);
    hpp_SJ40_L1ForJet_Et =new TH1D("SJ40_L1ForJet_Et"  , "L1_SingleJet40==1, L1 Forward Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ40_L1ForJet_Eta=new TH1D("SJ40_L1ForJet_Eta" , "L1_SingleJet40==1, L1 Forward Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ40_L1ForJet_Phi=new TH1D("SJ40_L1ForJet_Phi" , "L1_SingleJet40==1, L1 Forward Jet #phi" , 64  , -3.2,3.2 );
    //Central Jets
    hpp_SJ40_L1CenJet_E  =new TH1D("SJ40_L1CenJet_E"   , "L1_SingleJet40==1, L1 Central Jet E"    , 1000, 0., 1000.);
    hpp_SJ40_L1CenJet_Et =new TH1D("SJ40_L1CenJet_Et"  , "L1_SingleJet40==1, L1 Central Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ40_L1CenJet_Eta=new TH1D("SJ40_L1CenJet_Eta" , "L1_SingleJet40==1, L1 Central Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ40_L1CenJet_Phi=new TH1D("SJ40_L1CenJet_Phi" , "L1_SingleJet40==1, L1 Central Jet #phi" , 64  , -3.2,3.2 );


    //by L1 Trigger Path Firings; L1SingleJet28_BptxAND
    //L1 Jet Counts
    hpp_SJ28_NL1ForJets=new TH1D("SJ28_NL1ForJets", "L1_SingleJet28==1, N L1 Forward Jets", 8,0,8);
    hpp_SJ28_NL1CenJets=new TH1D("SJ28_NL1CenJets", "L1_SingleJet28==1, N L1 Central Jets", 8,0,8);
    hpp_SJ28_NL1Jets   =new TH1D("SJ28_NL1Jets"   , "L1_SingleJet28==1, N L1 Jets"        , 8,0,8);
    //Lead L1 Jets
    hpp_SJ28_LeadL1ForJet_Et=new TH1D("SJ28_LeadL1ForJet_Et", "L1_SingleJet28==1, Lead L1 Forward Jet E_{T}", 1000,0.,1000.);
    hpp_SJ28_LeadL1CenJet_Et=new TH1D("SJ28_LeadL1CenJet_Et", "L1_SingleJet28==1, Lead L1 Central Jet E_{T}", 1000,0.,1000.);
    hpp_SJ28_LeadL1Jet_Et   =new TH1D("SJ28_LeadL1Jet_Et"   , "L1_SingleJet28==1, Lead L1 Jet E_{T}"        , 1000,0.,1000.);
    //all L1 Jets
    hpp_SJ28_L1Jet_E  =new TH1D("SJ28_L1Jet_E"   , "L1_SingleJet28==1, L1 Jet E"    , 1000, 0., 1000.);
    hpp_SJ28_L1Jet_Et =new TH1D("SJ28_L1Jet_Et"  , "L1_SingleJet28==1, L1 Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ28_L1Jet_Eta=new TH1D("SJ28_L1Jet_Eta" , "L1_SingleJet28==1, L1 Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ28_L1Jet_Phi=new TH1D("SJ28_L1Jet_Phi" , "L1_SingleJet28==1, L1 Jet #phi" , 64  , -3.2,3.2 );
    //Forward L1 Jets
    hpp_SJ28_L1ForJet_E  =new TH1D("SJ28_L1ForJet_E"   , "L1_SingleJet28==1, L1 Forward Jet E"    , 1000, 0., 1000.);
    hpp_SJ28_L1ForJet_Et =new TH1D("SJ28_L1ForJet_Et"  , "L1_SingleJet28==1, L1 Forward Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ28_L1ForJet_Eta=new TH1D("SJ28_L1ForJet_Eta" , "L1_SingleJet28==1, L1 Forward Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ28_L1ForJet_Phi=new TH1D("SJ28_L1ForJet_Phi" , "L1_SingleJet28==1, L1 Forward Jet #phi" , 64  , -3.2,3.2 );
    //Central Jets
    hpp_SJ28_L1CenJet_E  =new TH1D("SJ28_L1CenJet_E"   , "L1_SingleJet28==1, L1 Central Jet E"    , 1000, 0., 1000.);
    hpp_SJ28_L1CenJet_Et =new TH1D("SJ28_L1CenJet_Et"  , "L1_SingleJet28==1, L1 Central Jet E_{T}", 1000, 0., 1000.);
    hpp_SJ28_L1CenJet_Eta=new TH1D("SJ28_L1CenJet_Eta" , "L1_SingleJet28==1, L1 Central Jet #eta" , 102 , -5.1, 5.1);
    hpp_SJ28_L1CenJet_Phi=new TH1D("SJ28_L1CenJet_Phi" , "L1_SingleJet28==1, L1 Central Jet #phi" , 64  , -3.2,3.2 );
    
  }
  
  
  

  
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
  int pBeamScrapingFilter_I, pHBHENoiseFilterRun2Loose_I, pprimaryvertexFilter_I,puvertexFilter_I;
  int pHBHEIsoNoiseFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilterRun2Loose_I);//what i have been using
  jetpp[2]->SetBranchAddress("HBHEIsoNoiseFilterResult",&pHBHEIsoNoiseFilter_I);//only option
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
  
  
  int NL1ForJet_I;
  float L1ForJetE_F[4], L1ForJetEt_F[4], L1ForJetEta_F[4], L1ForJetPhi_F[4];
  jetpp[4]->SetBranchAddress("NL1ForJet"   , &NL1ForJet_I    );
  jetpp[4]->SetBranchAddress( "L1ForJetE"   , &L1ForJetE_F   );
  jetpp[4]->SetBranchAddress( "L1ForJetEt"  , &L1ForJetEt_F  );
  jetpp[4]->SetBranchAddress( "L1ForJetEta" , &L1ForJetEta_F );
  jetpp[4]->SetBranchAddress( "L1ForJetPhi" , &L1ForJetPhi_F );
  
  int NL1CenJet_I;
  float L1CenJetE_F[4], L1CenJetEt_F[4], L1CenJetEta_F[4], L1CenJetPhi_F[4];
  jetpp[4]->SetBranchAddress("NL1CenJet"   , &NL1CenJet_I   );
  jetpp[4]->SetBranchAddress( "L1CenJetE"   , &L1CenJetE_F  );
  jetpp[4]->SetBranchAddress( "L1CenJetEt"  , &L1CenJetEt_F );
  jetpp[4]->SetBranchAddress( "L1CenJetEta" , &L1CenJetEta_F);
  jetpp[4]->SetBranchAddress( "L1CenJetPhi" , &L1CenJetPhi_F);



 
  
  
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
  

  
  
  //pfcandAnalyzer/pfTree
  Int_t nPFpart_I;  
  std::vector<Float_t> *pfpt, *pfphi, *pfeta, *pfE;
  Float_t pfptArr[2000], pfphiArr[2000], pfetaArr[2000], pfEArr[2000];
  
  if(usePFCandAnalyzer){
    if(debugMode)std::cout<<"using PFCandAnalyzer. Setting branch addresses!"<<std::endl;
    jetpp[5]->SetBranchAddress( "nPFpart", &nPFpart_I);
    if(filelistIsJet80){
      if(debugMode)std::cout<<"filelistIsJet80, pointing branches to std::vector<Float_t>!"<<std::endl;
      jetpp[5]->SetBranchAddress( "pfPt"     , &pfpt );
      jetpp[5]->SetBranchAddress( "pfEta"    , &pfeta);
      jetpp[5]->SetBranchAddress( "pfPhi"    , &pfphi);    
      jetpp[5]->SetBranchAddress( "pfEnergy" , &pfE);    
    }
    else if(filelistIsLowerJets){
      if(debugMode)std::cout<<"filelistIsLowerJets, pointing branches to Float_t* (a.k.a. Float_t array)"<<std::endl;
      jetpp[5]->SetBranchAddress( "pfPt"   , &pfptArr );
      jetpp[5]->SetBranchAddress( "pfEta"  , &pfetaArr);
      jetpp[5]->SetBranchAddress( "pfPhi"  , &pfphiArr);    
      jetpp[5]->SetBranchAddress( "pfEnergy"  , &pfEArr);    
    }
    else assert(false);
  }


  


  
  //ONE HLT path ONE tree ONE trig obj pt branch
  std::vector<double>  *trgObjpt_40, *trgObjpt_60, *trgObjpt_80, *trgObjpt_100;
  jetpp[7] ->SetBranchAddress("pt",&trgObjpt_40);
  jetpp[8] ->SetBranchAddress("pt",&trgObjpt_60);  
  jetpp[9] ->SetBranchAddress("pt",&trgObjpt_80);  
  jetpp[10]->SetBranchAddress("pt",&trgObjpt_100);
  
  std::vector<double>  *trgObjeta_40, *trgObjeta_60, *trgObjeta_80, *trgObjeta_100;
  jetpp[7]->SetBranchAddress("eta",&trgObjeta_40);
  jetpp[8]->SetBranchAddress("eta",&trgObjeta_60);  
  jetpp[9]->SetBranchAddress("eta",&trgObjeta_80);  
  jetpp[10]->SetBranchAddress("eta",&trgObjeta_100);
  

  
  
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






  
  //begin event loop
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {



    
    if (nEvt%100000==0 || deepDebug) {
      //if (true) {
      float nEvt_fracdone=((float)nEvt)/((float)NEvents_toRead);
      float nEvt_percdone=nEvt_fracdone*100.;
      int   nEvt_intpercdone=(int)nEvt_percdone;                
      std::cout<<"from trees, grabbing Evt # = "<<nEvt<<", ~"<<nEvt_intpercdone<<"% of job completed"<<std::endl;
    }
    
    //grab an entry    
    jetpp[0]->GetEntry(nEvt);
    
    
    h_NEvents_read->Fill(0.);    //h_NEvents_read->Fill(1.,weight_eS);        
    

    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilterRun2Loose_I    ==0    || 
	pHBHEIsoNoiseFilter_I    ==0    || 
	pBeamScrapingFilter_I ==0 || 
	pprimaryvertexFilter_I==0  )  continue;  
    h_NEvents_skimCut->Fill(0.); 

    
    if( fabs(vz_F)>24. )     continue;
    h_NEvents_vzCut->Fill(0.); //h_NEvents_vzCut->Fill(1.,weight_eS);    








    if(usePFCandAnalyzer){
      //PFCANDLOOP
      //std::cout<<std::endl<<"nPFpart_I="<<nPFpart_I<<std::endl;
      
      Float_t sumpfpt_x=0.,sumpfpt_y=0.;
      Float_t sumpfEt=0.;//scalar sum of PF Et, not a vector sum.
      //h_nPFpart->Fill(nPFpart_I,weight_eS);
      
      for(Int_t part=0;part<nPFpart_I;part++){
	
	Float_t pfpt_part,pfeta_part,pfphi_part, pfE_part;
	if(filelistIsJet80){	  
	  pfpt_part =pfpt ->at(part);	  
	  pfeta_part=pfeta->at(part);	  
	  pfphi_part=pfphi->at(part);	  
	  pfE_part=pfE->at(part);	}
	else if (filelistIsLowerJets){
	  pfpt_part  = pfptArr[part];	  
	  pfeta_part =pfetaArr[part];	  
	  pfphi_part =pfphiArr[part];	  	
	  pfE_part  =pfEArr[part];	}
	else assert(false);
	
	//h_pfpt ->Fill(pfpt_part ,weight_eS);
	//h_pfeta->Fill(pfeta_part,weight_eS);
	//h_pfphi->Fill(pfphi_part,weight_eS);
	//h_pfE->Fill(pfE_part,weight_eS);
	
	Float_t pftheta_part=2.*atan(exp(-1.*pfeta_part));
	//h_pftheta->Fill(pftheta_part,weight_eS);
	
	Float_t pfEt_part=pfE_part*TMath::Sin(pftheta_part);		
	//h_pfEt->Fill(pfEt_part,weight_eS);
	
	//std::cout<<"pfpt->at("<<part<<")=" << pfpt_part<<std::endl;
	//std::cout<<"pfphi->at("<<part<<")="<<pfphi_part<<std::endl;
	//std::cout<<"pfeta->at("<<part<<")="<<pfeta_part<<std::endl;
	sumpfpt_x+=pfpt_part*TMath::Cos(pfphi_part);
	sumpfpt_y+=pfpt_part*TMath::Sin(pfphi_part);
	sumpfEt+=pfEt_part;
      }
      
      Float_t sumpfpt=sqrt(sumpfpt_x*sumpfpt_x + sumpfpt_y*sumpfpt_y );
      //Float_t sqrtsumpfEt=sqrt(sumpfEt);
      Float_t PFMETfrac=sumpfpt/sumpfEt;
      //Float_t PFMETsig=sumpfpt/sqrtsumpfEt;
      
      if(PFMETfrac>0.3)continue;//giving this a shot. see how it goes...
      //h_NEvents_PFMETfracCut->Fill(0.);
      
      //h_sumpfpt_x->Fill(sumpfpt_x,weight_eS);
      //h_sumpfpt_y->Fill(sumpfpt_y,weight_eS);
      //h_sumpfEt->Fill(sumpfEt,weight_eS);

      //h_sumpfpt->Fill(sumpfpt,weight_eS);
      //h_sqrtsumpfEt->Fill(sqrtsumpfEt,weight_eS);
      //h_PFMETfrac->Fill(PFMETfrac,weight_eS);
      //h_PFMETsig->Fill(PFMETsig,weight_eS);
      //assert(false);
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


    double maxTrgPt=-1.;//,maxTrgEta=-999.; //globle leading object for an event, regardless of trigger exclusion
    double HLT40maxTrgPt=-1.,  HLT60maxTrgPt=-1.,  HLT80maxTrgPt=-1.,  HLT100maxTrgPt=-1.; 
    double HLT40maxTrgEta=-1.,  HLT60maxTrgEta=-1.,  HLT80maxTrgEta=-1., HLT100maxTrgEta=-1.; 
    unsigned int trgObj40_size=trgObjpt_40->size(), trgObj60_size=trgObjpt_60->size();
    unsigned int trgObj80_size=trgObjpt_80->size(), trgObj100_size=trgObjpt_100->size();   
    if(trgDec[3]&&useHLT100){
      for(unsigned int itt=0; itt<trgObj100_size; ++itt){
	double trgpt=trgObjpt_100->at(itt);
	double trgeta=trgObjeta_100->at(itt);
    	if(trgpt > maxTrgPt) { 
	  //maxTrgEta = trgeta;
	  maxTrgPt = trgpt;	}
	if(trgpt > HLT100maxTrgPt){
	  HLT100maxTrgPt=trgpt;
	  HLT100maxTrgEta=trgeta;	}
      }
    }
    

    if(trgDec[2]){
      for(unsigned int itt=0; itt<trgObj80_size; ++itt){	
	double trgpt=trgObjpt_80->at(itt);
	double trgeta=trgObjeta_80->at(itt);
    	if(trgpt > maxTrgPt) { 
	  //maxTrgEta = trgeta;
	  maxTrgPt  = trgpt;	}	  	
	if(trgpt > HLT80maxTrgPt){
	  HLT80maxTrgPt=trgpt;
	  HLT80maxTrgEta=trgeta;	}
      }
    }
    
    if(trgDec[1]) {
      for(unsigned int itt=0; itt<trgObj60_size; ++itt){
	double trgpt=trgObjpt_60->at(itt);
	double trgeta=trgObjeta_60->at(itt);
    	if(trgpt > maxTrgPt) { 
	  //maxTrgEta =trgeta;
	  maxTrgPt  =trgpt;	}
	if(trgpt > HLT60maxTrgPt){	  
	  HLT60maxTrgPt=trgpt;
	  HLT60maxTrgEta=trgeta;	}    
      }
    }
    
    if(trgDec[0])  {
      for(unsigned int itt=0; itt<trgObj40_size; ++itt)	{
	double trgpt=trgObjpt_40->at(itt);
	double trgeta=trgObjeta_40->at(itt);
    	if(trgpt > maxTrgPt) { 
	  //maxTrgEta = trgeta;
	  maxTrgPt  = trgpt;	} 
	if(trgpt > HLT40maxTrgPt){
	  HLT40maxTrgPt=trgpt;
	  HLT40maxTrgEta=trgeta;	}
      }
    }



    
    double trgPt=-1., trgEta=-999.;//for lead obj from final trigger decision path
    double weight_eS=0.;
    bool is40  = false, is60  = false, is80  = false, is100 = false, isMB=false;
    if(useHLT100){
      if(      trgDec[3] && !(HLT100maxTrgPt<HLTthresh[3])               ) 
	{ is100 = true;  weight_eS=(double)trgPscl[3]; trgPt=HLT100maxTrgPt; trgEta=HLT100maxTrgEta;}
      else if( trgDec[2] && !(HLT80maxTrgPt<HLTthresh[2])  && HLT80maxTrgPt<HLTthresh[3] ) 
	{ is80  = true;  weight_eS=(double)trgPscl[2]; trgPt=HLT80maxTrgPt;trgEta=HLT80maxTrgEta;}
      else if( trgDec[1] && !(HLT60maxTrgPt<HLTthresh[1])  && HLT60maxTrgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=(double)trgPscl[1]; trgPt=HLT60maxTrgPt;trgEta=HLT60maxTrgEta;}
      else if( trgDec[0] && !(HLT40maxTrgPt<HLTthresh[0])  && HLT40maxTrgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=(double)trgPscl[0]; trgPt=HLT40maxTrgPt;trgEta=HLT40maxTrgEta;}      
      else if( MBtrgDec ) { 
	  { isMB=true; weight_eS=(double)MBtrgPscl; }
      }
    }
    else {
      if( trgDec[2] && !(HLT80maxTrgPt<HLTthresh[2])  ) 
	{ is80  = true;  weight_eS=(double)trgPscl[2];  trgPt=HLT80maxTrgPt;trgEta=HLT80maxTrgEta;}
      else if( trgDec[1] && !(HLT60maxTrgPt<HLTthresh[1])  && HLT60maxTrgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=(double)trgPscl[1];  trgPt=HLT60maxTrgPt;trgEta=HLT60maxTrgEta;}
      else if( trgDec[0] && !(HLT40maxTrgPt<HLTthresh[0])  && HLT40maxTrgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=(double)trgPscl[0];  trgPt=HLT40maxTrgPt;trgEta=HLT40maxTrgEta;}
      else if( MBtrgDec ) { 
	isMB=true; weight_eS=(double)MBtrgPscl;	    	
      }
    }
    
    
    
    
    //if(debugMode) std::cout<<"weight_eS="<<weight_eS<<std::endl;    
    bool firedJetTrigger=(is40 || is60 || is80 || is100 );
    bool firedMBTrigger=(isMB);
    bool firedTrigger=firedJetTrigger || firedMBTrigger;
    
    //event counts, pre-duplicate skip, pre-evt cuts, pre-trigger reqs.
    h_NEvents->Fill(0.);   h_NEvents->Fill(1.,weight_eS);      



    
    //event counts, post-duplicate skip, pre-evt cuts, having fired a trigger
    if(firedTrigger)
      h_NEvents_trigd_1->Fill(0.);   h_NEvents_trigd_1->Fill(1.,weight_eS );    

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
    
    
    float LeadL1CenJetEt=0., LeadL1ForJetEt=0., LeadL1JetEt=0.;
    if(fillDataJetTrigQAHists){ //only want to fill these trigger jet plots if they pass all our quality criteria      
      

      
      //L1 JETS SECTION
      int N_L1CenJets=0, N_L1ForJets=0;
      for(int i=0; i<4; i++){
	if(L1ForJetE_F[i]>0.)N_L1ForJets++;
	else continue;	
	if(L1ForJetEt_F[i]>LeadL1ForJetEt)LeadL1ForJetEt=L1ForJetEt_F[i];
	if(L1ForJetEt_F[i]>LeadL1JetEt)LeadL1JetEt=L1ForJetEt_F[i];
	hpp_L1ForJet_E->Fill(L1ForJetE_F[i],1.);
	hpp_L1ForJet_Et->Fill(L1ForJetEt_F[i],1.);
	hpp_L1ForJet_Eta->Fill(L1ForJetEta_F[i],1.);
	hpp_L1ForJet_Phi->Fill(L1ForJetPhi_F[i],1.);	
	hpp_L1Jet_E->Fill(L1ForJetE_F[i],1.);
	hpp_L1Jet_Et->Fill(L1ForJetEt_F[i],1.);
	hpp_L1Jet_Eta->Fill(L1ForJetEta_F[i],1.);
	hpp_L1Jet_Phi->Fill(L1ForJetPhi_F[i],1.);	
	
	if(jet80_l1s_I==1){//L1SingleJet48
	  hpp_SJ48_L1ForJet_E->Fill(L1ForJetE_F[i]      ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1ForJet_Et->Fill(L1ForJetEt_F[i]    ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1ForJet_Eta->Fill(L1ForJetEta_F[i]  ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1ForJet_Phi->Fill(L1ForJetPhi_F[i]  ,(float)jet80_l1s_ps_I);	
	  hpp_SJ48_L1Jet_E->Fill(L1ForJetE_F[i]         ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1Jet_Et->Fill(L1ForJetEt_F[i]       ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1Jet_Eta->Fill(L1ForJetEta_F[i]     ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1Jet_Phi->Fill(L1ForJetPhi_F[i]     ,(float)jet80_l1s_ps_I);	}
	if(jet60_l1s_I==1){//L1SingleJet40		                        
	  hpp_SJ40_L1ForJet_E->Fill(L1ForJetE_F[i]	,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1ForJet_Et->Fill(L1ForJetEt_F[i]	,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1ForJet_Eta->Fill(L1ForJetEta_F[i]	,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1ForJet_Phi->Fill(L1ForJetPhi_F[i]	,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1Jet_E->Fill(L1ForJetE_F[i]		,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1Jet_Et->Fill(L1ForJetEt_F[i]	,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1Jet_Eta->Fill(L1ForJetEta_F[i]	,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1Jet_Phi->Fill(L1ForJetPhi_F[i]	,(float)jet60_l1s_ps_I);	}
	if(jet40_l1s_I==1){//L1SingleJet28		                        
	  hpp_SJ28_L1ForJet_E->Fill(L1ForJetE_F[i]	,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1ForJet_Et->Fill(L1ForJetEt_F[i]	,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1ForJet_Eta->Fill(L1ForJetEta_F[i]	,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1ForJet_Phi->Fill(L1ForJetPhi_F[i]	,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1Jet_E->Fill(L1ForJetE_F[i]		,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1Jet_Et->Fill(L1ForJetEt_F[i]	,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1Jet_Eta->Fill(L1ForJetEta_F[i]	,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1Jet_Phi->Fill(L1ForJetPhi_F[i]	,(float)jet40_l1s_ps_I);	}
      }
      hpp_NL1ForJets->Fill(N_L1ForJets,1.);
      if(N_L1ForJets>0)
	hpp_LeadL1ForJet_Et->Fill(LeadL1ForJetEt,1.);	
      
      
      for(int i=0; i<4; i++){
	if(L1CenJetE_F[i]>0.)N_L1CenJets++;
	else continue;	
	if(L1CenJetEt_F[i]>LeadL1CenJetEt)LeadL1CenJetEt=L1CenJetEt_F[i];
	if(L1CenJetEt_F[i]>LeadL1JetEt)LeadL1JetEt=L1CenJetEt_F[i];
	hpp_L1CenJet_E->Fill(L1CenJetE_F[i],1.);
	hpp_L1CenJet_Et->Fill(L1CenJetEt_F[i],1.);
	hpp_L1CenJet_Eta->Fill(L1CenJetEta_F[i],1.);
	hpp_L1CenJet_Phi->Fill(L1CenJetPhi_F[i],1.);	
	hpp_L1Jet_E->Fill(L1CenJetE_F[i],1.);
	hpp_L1Jet_Et->Fill(L1CenJetEt_F[i],1.);
	hpp_L1Jet_Eta->Fill(L1CenJetEta_F[i],1.);
	hpp_L1Jet_Phi->Fill(L1CenJetPhi_F[i],1.);	
	if(jet80_l1s_I==1){//L1SingleJet48
	  hpp_SJ48_L1CenJet_E->Fill(L1CenJetE_F[i]     ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1CenJet_Et->Fill(L1CenJetEt_F[i]   ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1CenJet_Eta->Fill(L1CenJetEta_F[i] ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1CenJet_Phi->Fill(L1CenJetPhi_F[i] ,(float)jet80_l1s_ps_I);	
	  hpp_SJ48_L1Jet_E->Fill(L1CenJetE_F[i]	       ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1Jet_Et->Fill(L1CenJetEt_F[i]      ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1Jet_Eta->Fill(L1CenJetEta_F[i]    ,(float)jet80_l1s_ps_I);
	  hpp_SJ48_L1Jet_Phi->Fill(L1CenJetPhi_F[i]    ,(float)jet80_l1s_ps_I);		}	
	if(jet60_l1s_I==1){//L1SingleJet40	                               
	  hpp_SJ40_L1CenJet_E->Fill(L1CenJetE_F[i]     ,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1CenJet_Et->Fill(L1CenJetEt_F[i]   ,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1CenJet_Eta->Fill(L1CenJetEta_F[i] ,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1CenJet_Phi->Fill(L1CenJetPhi_F[i] ,(float)jet60_l1s_ps_I);	
	  hpp_SJ40_L1Jet_E->Fill(L1CenJetE_F[i]	       ,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1Jet_Et->Fill(L1CenJetEt_F[i]      ,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1Jet_Eta->Fill(L1CenJetEta_F[i]    ,(float)jet60_l1s_ps_I);
	  hpp_SJ40_L1Jet_Phi->Fill(L1CenJetPhi_F[i]    ,(float)jet60_l1s_ps_I);		}	
	if(jet40_l1s_I==1){//L1SingleJet28	                               
	  hpp_SJ28_L1CenJet_E->Fill(L1CenJetE_F[i]     ,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1CenJet_Et->Fill(L1CenJetEt_F[i]   ,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1CenJet_Eta->Fill(L1CenJetEta_F[i] ,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1CenJet_Phi->Fill(L1CenJetPhi_F[i] ,(float)jet40_l1s_ps_I);	
	  hpp_SJ28_L1Jet_E->Fill(L1CenJetE_F[i]	       ,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1Jet_Et->Fill(L1CenJetEt_F[i]      ,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1Jet_Eta->Fill(L1CenJetEta_F[i]    ,(float)jet40_l1s_ps_I);
	  hpp_SJ28_L1Jet_Phi->Fill(L1CenJetPhi_F[i]    ,(float)jet40_l1s_ps_I);		}	
      }
      hpp_NL1CenJets->Fill(N_L1CenJets,1.);      
      if(N_L1CenJets>0)
	hpp_LeadL1CenJet_Et->Fill(LeadL1CenJetEt,1.);
      
      
      hpp_NL1Jets->Fill(N_L1CenJets+N_L1ForJets,1.);      
      hpp_LeadL1Jet_Et->Fill(LeadL1JetEt,1.);
      
      if(jet80_l1s_I==1){//L1SingleJet48
	hpp_SJ48_NL1ForJets->Fill(N_L1ForJets                          ,(float)jet80_l1s_ps_I);
	hpp_SJ48_NL1CenJets->Fill(N_L1CenJets      		       ,(float)jet80_l1s_ps_I);
	hpp_SJ48_NL1Jets->Fill((N_L1ForJets+N_L1CenJets)               ,(float)jet80_l1s_ps_I);
	if(N_L1ForJets>0)hpp_SJ48_LeadL1ForJet_Et->Fill(LeadL1ForJetEt ,(float)jet80_l1s_ps_I);
	if(N_L1CenJets>0)hpp_SJ48_LeadL1CenJet_Et->Fill(LeadL1CenJetEt ,(float)jet80_l1s_ps_I);     
	hpp_SJ48_LeadL1Jet_Et->Fill(LeadL1JetEt      	               ,(float)jet80_l1s_ps_I); }
      
      if(jet60_l1s_I==1){//L1SingleJet40
	hpp_SJ40_NL1ForJets->Fill(N_L1ForJets                          ,(float)jet60_l1s_ps_I);
	hpp_SJ40_NL1CenJets->Fill(N_L1CenJets      		       ,(float)jet60_l1s_ps_I);
	hpp_SJ40_NL1Jets->Fill((N_L1ForJets+N_L1CenJets) 	       ,(float)jet60_l1s_ps_I);
	if(N_L1ForJets>0)hpp_SJ40_LeadL1ForJet_Et->Fill(LeadL1ForJetEt ,(float)jet60_l1s_ps_I);
	if(N_L1CenJets>0)hpp_SJ40_LeadL1CenJet_Et->Fill(LeadL1CenJetEt ,(float)jet60_l1s_ps_I);     
	hpp_SJ40_LeadL1Jet_Et->Fill(LeadL1JetEt      	               ,(float)jet60_l1s_ps_I); }
      if(jet40_l1s_I==1){//L1SingleJet28
	hpp_SJ28_NL1ForJets->Fill(N_L1ForJets                          ,(float)jet40_l1s_ps_I);
	hpp_SJ28_NL1CenJets->Fill(N_L1CenJets      		       ,(float)jet40_l1s_ps_I);
	hpp_SJ28_NL1Jets->Fill((N_L1ForJets+N_L1CenJets) 	       ,(float)jet40_l1s_ps_I);
	if(N_L1ForJets>0)hpp_SJ28_LeadL1ForJet_Et->Fill(LeadL1ForJetEt ,(float)jet40_l1s_ps_I);
	if(N_L1CenJets>0)hpp_SJ28_LeadL1CenJet_Et->Fill(LeadL1CenJetEt ,(float)jet40_l1s_ps_I);      
	hpp_SJ28_LeadL1Jet_Et->Fill(LeadL1JetEt      	               ,(float)jet40_l1s_ps_I); }
      //END L1 JETS SECTION





      // --- BEGIN TRIGGER OBJECT(S) HIST SECTION --- //

      //EXCLUSIVE TRIGER OBJECT HISTS
      if(is40)hpp_is40pscl ->Fill(  weight_eS, 1.);
      if(is60)hpp_is60pscl ->Fill(  weight_eS, 1.);
      if(is80)hpp_is80pscl ->Fill(  weight_eS, 1.);

      //all trigger jet object pts
      if(is40)
	for(unsigned int itt=0; itt<trgObj40_size; ++itt)hpp_is40trgPt_allobj->Fill(trgObjpt_40->at(itt),weight_eS);
      if(is60)
	for(unsigned int itt=0; itt<trgObj60_size; ++itt)hpp_is60trgPt_allobj->Fill(trgObjpt_60->at(itt),weight_eS);
      if(is80)
	for(unsigned int itt=0; itt<trgObj80_size; ++itt)hpp_is80trgPt_allobj->Fill(trgObjpt_80->at(itt),weight_eS);

      if( is40  )    hpp_is40trgPt  ->Fill(  HLT40maxTrgPt  , weight_eS );      
      else if( is60  )    hpp_is60trgPt  ->Fill(  HLT60maxTrgPt  , weight_eS );
      else if( is80  )    hpp_is80trgPt  ->Fill(  HLT80maxTrgPt  , weight_eS );
      if(firedJetTrigger) hpp_HLTCombtrgPt->Fill(  trgPt          , weight_eS );   

      if     ( is40  )    hpp_is40trgEta  ->Fill(  HLT40maxTrgEta  , weight_eS );      
      else if( is60  )    hpp_is60trgEta  ->Fill(  HLT60maxTrgEta  , weight_eS );
      else if( is80  )    hpp_is80trgEta  ->Fill(  HLT80maxTrgEta  , weight_eS );
      if(firedJetTrigger) hpp_HLTCombtrgEta->Fill(  trgEta          , weight_eS );                         

      if(      is40  && fabs(HLT40maxTrgEta) <= 2.0)   hpp_is40trgPt_00eta20  ->Fill(  HLT40maxTrgPt  , weight_eS );      
      else if( is60  && fabs(HLT60maxTrgEta) <= 2.0)   hpp_is60trgPt_00eta20  ->Fill(  HLT60maxTrgPt  , weight_eS );
      else if( is80  && fabs(HLT80maxTrgEta) <= 2.0)   hpp_is80trgPt_00eta20  ->Fill(  HLT80maxTrgPt  , weight_eS );
      if(firedJetTrigger && fabs(trgEta)  <= 2.0)   hpp_HLTCombtrgPt_00eta20->Fill(  trgPt          , weight_eS );   

      if(      is40  && fabs(HLT40maxTrgEta) > 2.0)   hpp_is40trgPt_20eta51  ->Fill(     HLT40maxTrgPt , weight_eS );      
      else if( is60  && fabs(HLT60maxTrgEta) > 2.0)   hpp_is60trgPt_20eta51  ->Fill(     HLT60maxTrgPt , weight_eS );
      else if( is80  && fabs(HLT80maxTrgEta) > 2.0)   hpp_is80trgPt_20eta51  ->Fill(     HLT80maxTrgPt , weight_eS );
      if(firedJetTrigger && fabs(trgEta)  > 2.0)   hpp_HLTCombtrgPt_20eta51->Fill(     trgPt         , weight_eS );   
      
      if( is40 ) hpp_is40trgPt_leadobj_globleadobj->Fill( HLT40maxTrgPt , maxTrgPt,  weight_eS);      	
      else if( is60 ) hpp_is60trgPt_leadobj_globleadobj->Fill( HLT60maxTrgPt , maxTrgPt,  weight_eS);
      else if( is80 ) hpp_is80trgPt_leadobj_globleadobj->Fill( HLT80maxTrgPt , maxTrgPt,  weight_eS);
      
      if( is40 ) hpp_is40trgPt_trgEta_leadobj->Fill( HLT40maxTrgPt , HLT40maxTrgEta,  weight_eS);      	
      else if( is60 ) hpp_is60trgPt_trgEta_leadobj->Fill( HLT60maxTrgPt , HLT60maxTrgEta,  weight_eS);
      else if( is80 ) hpp_is80trgPt_trgEta_leadobj->Fill( HLT80maxTrgPt , HLT80maxTrgEta,  weight_eS);
      //EXCLUSIVE TRIGGER OBJECT HISTS END


      //INCLUSIVE TRIGER OBJECT HISTS
      if(trgDec[0])hpp_IncHLT40pscl ->Fill( (double)trgPscl[0], 1.);
      if(trgDec[1])hpp_IncHLT60pscl ->Fill( (double)trgPscl[1], 1.);
      if(trgDec[2])hpp_IncHLT80pscl ->Fill( (double)trgPscl[2], 1.);

      //all trigger jet object pts
      if(trgDec[0])
	for(unsigned int itt=0; itt<trgObj40_size; ++itt)hpp_IncHLT40trgPt_allobj->Fill(trgObjpt_40->at(itt),(double)trgPscl[0]);
      if(trgDec[1])
	for(unsigned int itt=0; itt<trgObj60_size; ++itt)hpp_IncHLT60trgPt_allobj->Fill(trgObjpt_60->at(itt),(double)trgPscl[1]);
      if(trgDec[2])
	for(unsigned int itt=0; itt<trgObj80_size; ++itt)hpp_IncHLT80trgPt_allobj->Fill(trgObjpt_80->at(itt),(double)trgPscl[2]);

      if( trgDec[0] )  hpp_IncHLT40trgEta-> Fill(HLT40maxTrgEta , (double)trgPscl[0]);      	
      if( trgDec[1] )  hpp_IncHLT60trgEta-> Fill(HLT60maxTrgEta , (double)trgPscl[1]);
      if( trgDec[2] )  hpp_IncHLT80trgEta-> Fill(HLT80maxTrgEta , (double)trgPscl[2]);


      if( trgDec[0] )  hpp_IncHLT40trgPt_leadobj-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);      	
      if( trgDec[1] )  hpp_IncHLT60trgPt_leadobj-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);
      if( trgDec[2] )  hpp_IncHLT80trgPt_leadobj-> Fill(HLT80maxTrgPt , (double)trgPscl[2]);


      if( trgDec[0] && (fabs(HLT40maxTrgEta<=2.0)))  hpp_IncHLT40trgPt_leadobj_00eta20-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);      	
      if( trgDec[1] && (fabs(HLT60maxTrgEta<=2.0)))  hpp_IncHLT60trgPt_leadobj_00eta20-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);
      if( trgDec[2] && (fabs(HLT80maxTrgEta<=2.0)))  hpp_IncHLT80trgPt_leadobj_00eta20-> Fill(HLT80maxTrgPt , (double)trgPscl[2]);


      if( trgDec[0] && (fabs(HLT40maxTrgEta >2.0)))  hpp_IncHLT40trgPt_leadobj_20eta51-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);      	
      if( trgDec[1] && (fabs(HLT60maxTrgEta >2.0)))  hpp_IncHLT60trgPt_leadobj_20eta51-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);
      if( trgDec[2] && (fabs(HLT80maxTrgEta >2.0)))  hpp_IncHLT80trgPt_leadobj_20eta51-> Fill(HLT80maxTrgPt , (double)trgPscl[2]);

      
      if( trgDec[0] ) hpp_IncHLT40trgPt_leadobj_globleadobj->Fill( HLT40maxTrgPt , maxTrgPt,  (double)trgPscl[0]);      	
      if( trgDec[1] ) hpp_IncHLT60trgPt_leadobj_globleadobj->Fill( HLT60maxTrgPt , maxTrgPt,  (double)trgPscl[1]);
      if( trgDec[2] ) hpp_IncHLT80trgPt_leadobj_globleadobj->Fill( HLT80maxTrgPt , maxTrgPt,  (double)trgPscl[2]);
      
      if( trgDec[0] ) hpp_IncHLT40trgPt_trgEta_leadobj->Fill( HLT40maxTrgPt , HLT40maxTrgEta, (double)trgPscl[0] );      	
      if( trgDec[1] ) hpp_IncHLT60trgPt_trgEta_leadobj->Fill( HLT60maxTrgPt , HLT60maxTrgEta, (double)trgPscl[1] );
      if( trgDec[2] ) hpp_IncHLT80trgPt_trgEta_leadobj->Fill( HLT80maxTrgPt , HLT80maxTrgEta, (double)trgPscl[2] );
      //INCLUSIVE TRIGER OBJECT HISTS END
 
      // --- END TRIGGER OBJECT(S) HIST SECTION --- //

      




      // --- BEGIN TRIGGER EFFICIENCY HIST SECTION --- //
      
      //for measuring trigger eff. of HLT80 --------
      if(trgDec[1] && HLT60maxTrgPt>=80. ){
	for(unsigned int itt=0; itt<trgObj60_size; ++itt)hpp_HLT80ref_IncHLT60trgPt_allobj->Fill(trgObjpt_60->at(itt),(double)trgPscl[1]);
	hpp_HLT80ref_IncHLT60trgPt_leadobj-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);		
	if(fabs(HLT60maxTrgEta)<=2.0) hpp_HLT80ref_IncHLT60trgPt_leadobj_00eta20-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);	
	if(fabs(HLT60maxTrgEta)>2.0)  hpp_HLT80ref_IncHLT60trgPt_leadobj_20eta51-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);	
	if(LeadL1JetEt>=48.){
	  hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);	
	  if(fabs(HLT60maxTrgEta)<=2.0) hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_00eta20-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);	
	  if(fabs(HLT60maxTrgEta)>2.0)  hpp_HLT80ref_IncHLT60trgPt_wL1Thresh_leadobj_20eta51-> Fill(HLT60maxTrgPt , (double)trgPscl[1]);	
	}
      }  

      //for measuring trigger eff. of HLT60 --------
      if(trgDec[0] && HLT40maxTrgPt>=60. ){
	for(unsigned int itt=0; itt<trgObj40_size; ++itt)hpp_HLT60ref_IncHLT40trgPt_allobj->Fill(trgObjpt_40->at(itt),(double)trgPscl[0]);
	hpp_HLT60ref_IncHLT40trgPt_leadobj-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);	
	if(fabs(HLT40maxTrgEta)<=2.0) hpp_HLT60ref_IncHLT40trgPt_leadobj_00eta20-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);	
	if(fabs(HLT40maxTrgEta)>2.0)  hpp_HLT60ref_IncHLT40trgPt_leadobj_20eta51-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);	
	if(LeadL1JetEt>=40.){
	  hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);	
	  if(fabs(HLT40maxTrgEta)<=2.0) hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_00eta20-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);	
	  if(fabs(HLT40maxTrgEta)>2.0)  hpp_HLT60ref_IncHLT40trgPt_wL1Thresh_leadobj_20eta51-> Fill(HLT40maxTrgPt , (double)trgPscl[0]);	
	}
      }
      //for measuring trigger eff. of HLT40 --------
      //NULL

      // --- END TRIGGER EFFICIENCY HIST SECTION --- //      

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
    //if(!firedTrigger)continue;//doesn't matter if a trigger is fired exclusively for this macro; i want something very specific 
    
    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;

    // for event counting + avoiding duplicate fills in dijet hists
    bool dijetHistsFilled=false;    
    float evt_leadJetPt=-1., evt_leadJetPt_wCuts=-1.;
    
    int jetsPerEvent=-1;//=tupel_njets;    
    

    bool evtHasJetInFinalSpectra=false;

    if(useIncJetAnalyzer){
      
      jetsPerEvent=nref_I;
      for(int jet = 0; jet<nref_I; ++jet){
	if(deepDebug)std::cout<<"jet="<<jet<<std::endl;
	
	float jtpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
	//float rawpt=rawpt_F[jet];      
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
	
	

	

	// 13 TeV JetID criterion, loose or tight
	bool passesJetID=false;
	//float jetIDpt=jtpt;//ala HIN jetID, jtpt is corrected w/ L2/L3 residuals
	
	if (!(absreceta > 2.4)) 
	  passesJetID=(bool)jetID_00eta24( jtpt, 
					   neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet], muSum_F[jet],
					   numConst,  chMult, useTightJetID);
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  passesJetID=(bool) jetID_24eta27( jtpt,
					    neSum_F[jet],  phSum_F[jet], muSum_F[jet],
					    numConst, useTightJetID);
	else if( !(absreceta>3.0) && absreceta>2.7 )
	  passesJetID=(bool) jetID_27eta30( jtpt,
					    neSum_F[jet],  phSum_F[jet], 
					    neuMult, useTightJetID);
	else  
	  passesJetID=(bool)jetID_32eta47( jtpt, 
					   phSum_F[jet], useTightJetID);
	
	
	// get repidity bin
	int theRapBin=-1;	
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin)
	  if( absreceta>=absetabins[rapbin]  && 		
	      absreceta<absetabins[rapbin+1]    	      ) {	    
	    theRapBin=rapbin;
	    break;
	  }

	if(!evtHasJetInFinalSpectra)
	  if(!(absreceta>2.0) && passesJetID){
	    evtHasJetInFinalSpectra=true;
	    if(trgDec[0]) hpp_IncHLT40trgPt_leadobj_finalSpectra->Fill(HLT40maxTrgPt,trgPscl[0]);
	    if(trgDec[1]) hpp_IncHLT60trgPt_leadobj_finalSpectra->Fill(HLT60maxTrgPt,trgPscl[1]);
	    if(trgDec[2]) hpp_IncHLT80trgPt_leadobj_finalSpectra->Fill(HLT80maxTrgPt,trgPscl[2]);
	    if(trgDec[0]) hpp_IncHLT40trgEta_leadobj_finalSpectra->Fill(HLT40maxTrgEta,trgPscl[0]);
	    if(trgDec[1]) hpp_IncHLT60trgEta_leadobj_finalSpectra->Fill(HLT60maxTrgEta,trgPscl[1]);
	    if(trgDec[2]) hpp_IncHLT80trgEta_leadobj_finalSpectra->Fill(HLT80maxTrgEta,trgPscl[2]);

	    if(is40) hpp_is40trgPt_leadobj_finalSpectra->Fill(HLT40maxTrgPt,weight_eS);
	    if(is60) hpp_is60trgPt_leadobj_finalSpectra->Fill(HLT60maxTrgPt,weight_eS);
	    if(is80) hpp_is80trgPt_leadobj_finalSpectra->Fill(HLT80maxTrgPt,weight_eS);
	    if(is40) hpp_is40trgEta_leadobj_finalSpectra->Fill(HLT40maxTrgEta,weight_eS);
	    if(is60) hpp_is60trgEta_leadobj_finalSpectra->Fill(HLT60maxTrgEta,weight_eS);
	    if(is80) hpp_is80trgEta_leadobj_finalSpectra->Fill(HLT80maxTrgEta,weight_eS);

	  }
	
	
      


      
	
      
	if(!(jtpt>jetQAPtCut)){ jetsPerEvent--; continue;}//need JEC syst. stuff before this cut; so jets from below 56 GeV can migrate above	
	h_NJets_jetQAPtCut->Fill(0.);	  h_NJets_jetQAPtCut->Fill(1., weight_eS);


	

	//STUFF I WANT FILLED FOR 56 GEV AND ABOVE ONLY
	for(int jtid=0; jtid<2; jtid++){
	  
	  if(jtid==1 ){
	    
	    if(!passesJetID){ jetsPerEvent--; continue;	  }
	    else {  h_NJets_JetIDCut->Fill(0.);  h_NJets_JetIDCut->Fill(1., weight_eS); }
	    
	    if(jtpt>evt_leadJetPt_wCuts)
	      evt_leadJetPt_wCuts=jtpt;      
	    
	    
	  }//end jtid==1 specific. 	  
	  if(deepDebug)std::cout<<"jtid="<<jtid<<std::endl;
	  
	  
	  
	  
	  
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
	    

	    if(jtid==1){
	      //TRIGGER EFFICIENCY HISTS
	      //for measuring trigger eff. of HLT80 --------
	      if(trgDec[1] && HLT60maxTrgPt>80. ){
		hpp_HLT80ref_HLT60[theRapBin]->Fill(jtpt, (double)trgPscl[1]);  
	      }
	      if(trgDec[1] && HLT60maxTrgPt>80. && LeadL1JetEt>48. ){
		hpp_HLT80ref_HLT60_wL1Thresh[theRapBin]->Fill(jtpt, (double)trgPscl[1]);  
	      }
	      //if(trgDec[0] && HLT40maxTrgPt>80. ){
	      //	hpp_HLT80ref_HLT40[theRapBin]->Fill(jtpt, (double)trgPscl[0]);  
	      //}
	      //if(trgDec[0] && HLT40maxTrgPt>80. && LeadL1JetEt>48.){
	      //	hpp_HLT80ref_HLT40_wL1Thresh[theRapBin]->Fill(jtpt, (double)trgPscl[0]);  
	      //}
	      //for measuring trigger eff. of HLT60 --------
	      if(trgDec[0] && HLT40maxTrgPt>60. ){
		hpp_HLT60ref_HLT40[theRapBin]->Fill(jtpt, (double)trgPscl[0]);  
	      }
	      if(trgDec[0] && HLT40maxTrgPt>60. && LeadL1JetEt>40.){
		hpp_HLT60ref_HLT40_wL1Thresh[theRapBin]->Fill(jtpt, (double)trgPscl[0]);  
	      }
	      //for measuring trigger eff. of HLT40 --------
	      //NULL? USE MIN BIAS?
	      
	    }
	    
	    
	    //with exclusion
	    if(isMB)   {
	      hpp_excHLTMB[theRapBin][jtid]->Fill( jtpt , weight_eS);	//minbias	
	      if( fabs(eta_F[jet])<0.5 && jtid==1 && passesJetID ){
		
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
	      hpp_excHLT100[theRapBin][jtid]->Fill( jtpt , weight_eS); 	    //
	    }
	    
	  }//end fill data jettrig qa hists	
	    
	  
	  
	  if(fillDataJetQAHists){
	    
	    int ind=0;
	    //if(jtpt>jetQAPtCut){//second jet pt cut
	    
	    if(deepDebug){
	      std::cout<<"  theRapBin="<<  theRapBin<<std::endl;
	      std::cout<<"       jtid="<<       jtid<<std::endl;
	      std::cout<<"       jtpt="<<       jtpt<<std::endl;
	      std::cout<<" eta_F[jet]="<< eta_F[jet]<<std::endl;
	      assert(false);	    	    
	    }
	    

	    
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
    

    
    
    if(fillDataEvtQAHists){
      hjetsPEvt ->Fill(jetsPerEvent,1.0); 
      hWjetsPEvt->Fill(jetsPerEvent,weight_eS);     
      if(evt_leadJetPt>0.) 	  hLeadJetPt->Fill(evt_leadJetPt ,weight_eS);
      if(evt_leadJetPt_wCuts>0.)  hLeadJetPt_wCuts->Fill(evt_leadJetPt_wCuts ,weight_eS);
    }
    

    //assert(false);
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
  
  std::cout<<"readForests_ppData_jetTrigEff finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  std::cout<<"output file written is  "<<outfile<<std::endl<<std::endl;
  return 0;
}    // end readForests_ppData_jetTrigEff


////// main //////
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1 && argc!=10){
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppData_jetTrigEff.exe ";
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
  if(argc==1) rStatus = readForests_ppData_jetTrigEff();
  else if (argc==8){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppData_jetTrigEff( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
				  outputFileName);
  }
  else if (argc==10){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
    rStatus = readForests_ppData_jetTrigEff( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
					   outputFileName, jtEtaCutLo, jtEtaCutHi);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}



