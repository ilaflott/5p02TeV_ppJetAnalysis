// Custom header
#include "readForests_jetPlots.h"


// ppMC switches
const bool fillMCEvtQAHists=true;

const bool useIncJetAnalyzer=true;
const bool useTupel=false;

const bool fillMCJetQAHists=true;
const bool fillgenJetQA=true&fillMCJetQAHists;
const bool fillMCTupelJetQAHists=false&&useTupel;

const bool fillMCJetJECQAHists=false;
const bool fillMCDijetHists=false;//hardly used options




//// readForests_ppMC_jetPlots
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_jetPlots(std::string inFilelist , int startfile , int endfile ,  
			      int radius , std::string jetType , bool debugMode ,
			      std::string outfile, float jtEtaCutLo, float jtEtaCutHi){

  assert(fillMCJetQAHists!=fillMCTupelJetQAHists);//just dont use them at the same time. just. dont.
  
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
  std::cout<<", genjtPtCut = "<<genJetPtCut<<std::endl;
  std::cout<<", genjtPtCut_Hi = "<<genJetPtCut_Hi<<std::endl;  
  std::cout<<", gendrCut = "<<gendrCut<<std::endl;  
  std::cout<<"///////////////////"<<std::endl<<std::endl;
  
  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  
  // form jetTreeName
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;
  //std::string CaloJetTreeName="ak"+std::to_string(radius)+"CaloJetAnalyzer/t";
  //if(debugMode)std::cout<<"looking at CalojetTree "<<CaloJetTreeName<<std::endl;
  
  // initialize tree name array
  std::string trees[N_MCTrees];
  trees[0]=jetTreeName; //trees[1]=CaloJetTreeName;
  for(int i=1;i<N_MCTrees;++i)trees[i]=MCTreeNames[i];
  
  // declare TChains for each tree + friend  
  TChain* jetpp[N_MCTrees]={};
  for(int t = 0;t<N_MCTrees;++t)  {     
    if(!useTupel && t==4)continue;
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0)jetpp[0]->AddFriend( jetpp[t] );  
  }
  
  
  // open filelist
  std::cout<<"opening filelist: "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  std::string lastFileAdded=""; bool filesAdded=false; //Int_t filesCount=0;
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
    for(int t = 0;t<N_MCTrees;++t) {
      if(!useTupel &&t==4)continue;
      filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
    }
    lastFileAdded=filename_Forest;
  }//end file loop
  
  if(debugMode)std::cout<<"filesAdded="<<filesAdded<<std::endl;
  assert(filesAdded);//avoid segfault later
  
  
  // Declare the output file, declare hists after for easy write
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();
  
  
  // declare hists
  /////   EVT+JET COUNTS AND METADATA   ///// 
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"     ,(std::to_string(jtPtCut)).c_str()   ,    100, 0.,100.); hJetPtCut->Fill(jtPtCut);           
  TH1F *hJetPtCut_Hi        =new TH1F("hJetPtCut_Hi"     ,(std::to_string(jtPtCut_Hi)).c_str()   ,    1000, 1500.,2500.); hJetPtCut_Hi->Fill(jtPtCut_Hi);           
  TH1F *hGenJetPtCut     =new TH1F("hGenJetPtCut"  ,(std::to_string(genJetPtCut)).c_str()   ,    100, 0.,100.); hGenJetPtCut->Fill(genJetPtCut);           
  TH1F *hGenJetPtCut_Hi     =new TH1F("hGenJetPtCut_Hi"  ,(std::to_string(genJetPtCut_Hi)).c_str()   ,    1000, 1500.,2500.); hGenJetPtCut_Hi->Fill(genJetPtCut_Hi);             
  TH1F *hJetEtaCutHi     =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6  ); hJetEtaCutHi->Fill(jtEtaCutHi);	        
  TH1F *hJetEtaCutLo     =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6  ); hJetEtaCutLo->Fill(jtEtaCutLo);	        
  
  TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);  
  TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  
  
  
  //evt conts
  TH1D *h_NEvents          = new TH1D("NEvents"         , "NEvents",         2,0.,2.); //amount read directly from tree entries
  TH1D *h_NEvents_read     = new TH1D("NEvents_read"    , "NEvents read",    2,0.,2.);
  TH1D *h_NEvents_skimCut  = new TH1D("NEvents_skimCut" , "NEvents read post skimCut",  2,0.,2.);
  TH1D *h_NEvents_vzCut    = new TH1D("NEvents_vzCut"   , "NEvents read post vzCut AND skimCut", 2,0.,2.);
  
  /////   EVENT QA HISTS   ///// 
  TH1D *hVr=NULL, *hWVr=NULL;
  TH1D *hVx=NULL, *hWVx=NULL;
  TH1D *hVy=NULL, *hWVy=NULL;

  TH1D *hVz=NULL,  *hvzWVz=NULL,*hpthatWVz=NULL, *hWVz=NULL ;
  
  TH1D *hpthat=NULL, *hWpthat=NULL;  

  TH1D* hNref=NULL, *hWNref=NULL;
  TH1D* hjetsPEvt=NULL,* hWjetsPEvt=NULL;
  TH1D* hLeadJetPt=NULL,* hLeadJetPt_wCuts=NULL;
  
  if(fillMCEvtQAHists){
    
    hpthat    = new TH1D("hpthat","",1000,0,1000);//evt pthat, unweighted and weighted
    hWpthat   = new TH1D("hWeightedpthat","",1000,0,1000);  
    
    hVr = new TH1D("hVr","vr, no trig, no weights",         2000,0.,0.60); 
    hWVr = new TH1D("hWeightedVr","vr, trigd, with weights",2000,0.,0.60); 
    
    hVx = new TH1D("hVx","vx, no trig, no weights",  2000, -0.30,0.30);
    hWVx = new TH1D("hWeightedVx","vx, trigd, with weights", 2000, -0.30,0.30); 
    
    hVy = new TH1D("hVy","vy, no trig, no weights",   2000, -0.30,0.30); 
    hWVy = new TH1D("hWeightedVy","vy, trigd, with weights",   2000,-0.30,0.30);
    
    hVz       = new TH1D("hVz","vz, no trig, no weights", 1000,-25.,25.);//evtvz, weight=1 always
    hvzWVz    = new TH1D("hvzWeightedVz","vz, no trig, vz weights", 1000,-25.,25.);//vz-weighted only evtvz
    hpthatWVz = new TH1D("hpthatWeightedVz","vz, no trig, pthat weights", 1000,-25.,25.);//pthat-weighted only evtvz
    hWVz      = new TH1D("hWeightedVz","vz, no trig, with weights", 1000,-25.,25.);//pthat*vz-weighted evt vz
    
    hNref = new TH1D("hNref","numJets each evt",30,0,30);
    hWNref = new TH1D("hWNref","weighted numJets each evt",30,0,30);
    hjetsPEvt  = new TH1D("hjetsPEvt","numJets each evt",30,0,30);
    hWjetsPEvt = new TH1D("hWjetsPEvt","weighted numJets each evt",30,0,30);
    
    hLeadJetPt      = new TH1D("hLeadJetPt","lead jet pt for each evt",2500,0,2500);  
    hLeadJetPt_wCuts      = new TH1D("hLeadJetPt_wCuts","lead jet pt for each evt",2500,0,2500);
    
  }
  
  
  
  //jet counts
  TH1D *h_NJets          = new TH1D("NJets","NJets read", 2,0.,2.);
  
  TH1D *h_NJets_matched  = new TH1D("NJets_matched","NJets gen-reco matched", 2,0.,2.);//do truth/gen-level cuts first.
  TH1D *h_NJets_deltaRCut  = new TH1D("NJets_deltaRCut","NJets gen-reco deltaRcut", 2,0.,2.);  
  TH1D *h_NJets_genjtptCut  = new TH1D("NJets_genjtptCut ","NJets read post genjtptCut ", 2,0.,2.);
  TH1D *h_NJets_genjtptCut_Hi  = new TH1D("NJets_genjtptCut_Hi ","NJets read post genjtptCut_Hi ", 2,0.,2.);
  
  TH1D *h_NJets_jtptCut  = new TH1D("NJets_jtptCut ","NJets read post jtptCut ", 2,0.,2.); //reco level
  TH1D *h_NJets_jtptCut_Hi  = new TH1D("NJets_jtptCut_Hi ","NJets read post jtptCut_Hi ", 2,0.,2.);  
  TH1D *h_NJets_jtetaCut1  = new TH1D("NJets_jtetaCut1 ","NJets read post jtetaCut1 ", 2,0.,2.);//reco level
  TH1D *h_NJets_jtetaCut2  = new TH1D("NJets_jtetaCut2 ","NJets read post jtetaCut2 ", 2,0.,2.);  
  TH1D *h_NJets_jetQAPtCut  = new TH1D("NJets_jetQAPtCut ","NJets read post jetQAPtCut ", 2,0.,2.); //reco level
  TH1D *h_NJets_JetIDCut = new TH1D("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 2,0.,2.);
  


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


  if(fillMCJetJECQAHists){
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


  /////   JET QA HISTS   ///// 
  TH1D *hJetQA_jtptEntries[nbins_abseta][2]={};
  if(fillMCJetQAHists || fillMCTupelJetQAHists){
    for(int k = 0; k<nbins_abseta; ++k){      
      for(int i=0; i<2; i++){
	hJetQA_jtptEntries[k][i] = new TH1D( Form("hJetQA_%dwJetID_jtptEntries_etabin%d", i,k) , "N_{Jets}, no Pscl. Weights;N_{Jets};Jet p_{T} [GeV]" , 2500,0,2500);                 
      }
    }
  }
  
  
  
  //jet QA
  TH1D *hJetQA[nbins_abseta][N_vars][2] ={};
  
  //MC JetQA
  TH1D *hMCJetQA_genpt[nbins_abseta][2]={};
  TH1D *hMCJetQA_geneta[nbins_abseta][2]={};
  TH1D *hMCJetQA_genphi[nbins_abseta][2]={};
  TH1D *hMCJetQA_genrecpt[nbins_abseta][2]={};
  TH1D *hMCJetQA_genreceta[nbins_abseta][2]={};
  TH1D *hMCJetQA_genrecphi[nbins_abseta][2]={};
  TH1D *hMCJetQA_gendrjt[nbins_abseta][2]={};
  
  if(fillMCJetQAHists){

    for(int j = 0; j<N_vars; ++j){	
      for(int k = 0; k<nbins_abseta; ++k){
	for(int i=0; i<2; i++){
	  //binnings+hist declarations, etc.
	  
	  //jets
	  if(var[j]=="jtpt"||var[j]=="rawpt")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;", var[j].c_str()) , 2500,0,2500);
	  else if(var[j]=="jteta")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;",var[j].c_str()) , 100,-5,+5);
	  else if(var[j]=="jtphi")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;",var[j].c_str()) , 100,-4,+4);
	  
	  //jetconst. counts
	  else if(var[j]=="trkN"|| var[j]=="phN"|| var[j]=="chN"|| var[j]=="neN"|| var[j]=="eN"|| var[j]=="muN")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 60,0,60);
	  else if(var[j]=="trkHardN"|| var[j]=="phHardN"|| var[j]=="chHardN")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 60,0,60);
	  else if(var[j]=="neuMult"|| var[j]=="chMult"|| var[j]=="numConst")
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 100,0,100);
	  
	  //dijets
	  else if (fillMCDijetHists){
	    if(var[j]=="dphi") 
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k) , Form(";%s;",var[j].c_str()) , 50,0,+4);
	    else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt")
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;", var[j].c_str()), 2500,0,2500);
	    else //xj and Aj binnings
	      hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 200,0,2);         
	  }
	  
	  //consituent binnings
	  else 
	    hJetQA[k][j][i] = new TH1D( Form("hJetQA_%dwJetID_%s_etabin%d", i,var[j].c_str(),k), Form(";%s;",var[j].c_str()), 200,0,2);         
	}//jtid
      }//nbins abseta
    }//nvars
    
    if(fillgenJetQA){	
      for(int k = 0; k<nbins_abseta; ++k){
	for(int i=0; i<2; i++){
	  
	  hMCJetQA_genpt[k][i]     = new TH1D(    Form("hJetQA_%dwJetID_genpt_etabin%d",i,k)     , "genpt",    2500, 0,2500 );
	  hMCJetQA_geneta[k][i]    = new TH1D(   Form("hJetQA_%dwJetID_geneta_etabin%d",i,k)    , "geneta", 200,  -6.,6. );
	  hMCJetQA_genphi[k][i]    = new TH1D(   Form("hJetQA_%dwJetID_genphi_etabin%d",i,k)    , "genphi", 1280,  -6.4,6.4 );
	  hMCJetQA_genrecpt[k][i]  = new TH1D( Form("hJetQA_%dwJetID_genrecpt_etabin%d",i,k)  , "genrecpt",    200, 0., 2. );
	  hMCJetQA_genreceta[k][i] = new TH1D(Form("hJetQA_%dwJetID_genreceta_etabin%d",i,k) , "genreceta", 200,  -10, 10 ); 
	  hMCJetQA_genrecphi[k][i] = new TH1D(Form("hJetQA_%dwJetID_genrecphi_etabin%d",i,k) , "genrecphi", 256,  -12.8, 12.8 ); 
	  hMCJetQA_gendrjt[k][i] = new TH1D(    Form("hJetQA_%dwJetID_gendrjt_etabin%d",i,k) , "gendrjt", 200,  -1., 1. ); 		
	}//jtid 
      }//nbins abseta
    }//genjet qa
  }//fill mc jet qa
    


  


  //jet QA, i.e. hists filled once per jet
  TH1D *hTupelJetQA[nbins_abseta][N_tupelvars][2]={};    
  TH1D* hTupelGJetQA_genpt[nbins_abseta][2]={};
  TH1D* hTupelGJetQA_geneta[nbins_abseta][2]={};
  TH1D* hTupelGJetQA_genphi[nbins_abseta][2]={};
  TH1D* hTupelGJetQA_genenergy[nbins_abseta][2]={};
  TH1D* hTupelGJetQA_gendrjt[nbins_abseta][2]={};
  
  if(fillMCTupelJetQAHists){
    for(int k = 0; k<N_tupelvars; ++k){		      
      for(int j = 0; j<nbins_abseta; ++j){	
	for (int i=0; i<2; i++)	  {
	  hTupelJetQA[j][k][i]= new TH1D( ("hTupelJetQA_"+std::to_string(i)+"jetID_etabin"+std::to_string(j)+"_"+tupelvar[k]).c_str(),
					  (absetabins_str[j]+ "< #||{y} <"+absetabins_str[j+1]+";"+tupelvar[k]+";N_{Jets} x Pscl. Weight").c_str(), 
					  tupel_nbins[k], tupel_xlo[k], tupel_xhi[k] );	
	}      
      }    
    }
    
    if(fillgenJetQA){
      for(int j = 0; j<nbins_abseta; ++j){	
	for (int i=0; i<2; i++)	  {
	  hTupelGJetQA_genpt[j][i]=new TH1D ( Form("hTupelGJetQA_%djetID_etabin%d_genpt" , i,j ), "",  2000, 0, 2000);        
	  hTupelGJetQA_geneta[j][i]=new TH1D ( Form("hTupelGJetQA_%djetID_etabin%d_geneta" , i,j ), "",  1032, -5.16, 5.16);        
	  hTupelGJetQA_genphi[j][i]=new TH1D ( Form("hTupelGJetQA_%djetID_etabin%d_genphi" , i,j ), "",  730, -3.15, 3.15);        
	  hTupelGJetQA_genenergy[j][i]=new TH1D ( Form("hTupelGJetQA_%djetID_etabin%d_genenergy" , i,j ), "",  2000, 0., 2000);        
	  hTupelGJetQA_gendrjt[j][i]=new TH1D ( Form("hTupelGJetQA_%djetID_etabin%d_gendrjt" , i,j ), "",  60, 0, 0.6);        
	}	
      }    
    }    
  }
  //assert(false);
  
  
  
  // EVENT LOOP PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  float eta_F[1000];   float phi_F[1000];  
  float rawpt_F[1000]; //float y_F[1000], float jtpu_F[1000];

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
  
  //MC jet variables
  float pthat_F;  int subid_F[1000];    int refparton_F[1000];
  float refpt_F[1000];  float refeta_F[1000];  float refphi_F[1000];  
  float refdrjt_F[1000];

  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  //jetpp[0]->SetBranchAddress("jty",&y_F);
  //jetpp[0]->SetBranchAddress("jtpu",&jtpu_F);

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

  //MC jet variables
  jetpp[0]->SetBranchAddress("pthat",&pthat_F);
  jetpp[0]->SetBranchAddress( "subid"	, &subid_F     );
  jetpp[0]->SetBranchAddress( "refparton_flavor" , &refparton_F );
  jetpp[0]->SetBranchAddress( "refpt"	, &refpt_F     );
  jetpp[0]->SetBranchAddress( "refeta"  , &refeta_F    );
  jetpp[0]->SetBranchAddress( "refphi"  , &refphi_F    );
  jetpp[0]->SetBranchAddress( "refdrjt"	, &refdrjt_F   );  
  
  // HiEvtAnalyzer
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I; 
  float vz_F,vx_F,vy_F;
  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);
  jetpp[1]->SetBranchAddress("vz",&vz_F);
  jetpp[1]->SetBranchAddress("vx",&vx_F);
  jetpp[1]->SetBranchAddress("vy",&vy_F);

  
  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilter_I, pprimaryvertexFilter_I, puvertexFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);
  jetpp[2]->SetBranchAddress("pVertexFilterCutGtight",&puvertexFilter_I);
  

  // ppTrack/trackTree
  int nVtx_I;
  int nTrkVtx_I[1000];
  float xVtx_F[1000];
  float yVtx_F[1000];
  float zVtx_F[1000];

  //int nTrk_I[1000];
  //float trkPt_F[10000];
  //float trkEta_F[10000];
  //float trkPhi_F[10000];

  jetpp[3]->SetBranchAddress("nVtx",&nVtx_I);
  jetpp[3]->SetBranchAddress("nTrkVtx",&nTrkVtx_I);
  jetpp[3]->SetBranchAddress("xVtx",&xVtx_F);
  jetpp[3]->SetBranchAddress("yVtx",&yVtx_F);
  jetpp[3]->SetBranchAddress("zVtx",&zVtx_F);  




  //tupel/EventTree
  //RECO jet vars
  std::vector<float> *jetPt, *jetEta, *jetPhi, *jetE, *jetID, *jetRawPt, *jetRawE;
  std::vector<float> *jetHfHadE, *jetHfEmE, *jetChHadFrac, *jetNeutralHadAndHfFrac, *jetChEmFrac, *jetNeutralEmFrac, *jetChMult, *jetConstCnt;  
  if(useTupel){
    jetpp[4]->SetBranchAddress( "JetAk04Pt",                  &jetPt                  );
    jetpp[4]->SetBranchAddress( "JetAk04Eta",                 &jetEta                 );
    jetpp[4]->SetBranchAddress( "JetAk04Phi",                 &jetPhi                 );
    jetpp[4]->SetBranchAddress( "JetAk04E",                   &jetE                   );
    jetpp[4]->SetBranchAddress( "JetAk04Id",                  &jetID                  );
    jetpp[4]->SetBranchAddress( "JetAk04RawPt",               &jetRawPt               );
    jetpp[4]->SetBranchAddress( "JetAk04RawE",                &jetRawE                );
    jetpp[4]->SetBranchAddress( "JetAk04HfHadE",              &jetHfHadE              );
    jetpp[4]->SetBranchAddress( "JetAk04HfEmE",               &jetHfEmE               );
    jetpp[4]->SetBranchAddress( "JetAk04ChHadFrac",           &jetChHadFrac           );
    jetpp[4]->SetBranchAddress( "JetAk04NeutralHadAndHfFrac", &jetNeutralHadAndHfFrac );
    jetpp[4]->SetBranchAddress( "JetAk04ChEmFrac",            &jetChEmFrac            );
    jetpp[4]->SetBranchAddress( "JetAk04NeutralEmFrac",       &jetNeutralEmFrac       );
    jetpp[4]->SetBranchAddress( "JetAk04ChMult",              &jetChMult              );
    jetpp[4]->SetBranchAddress( "JetAk04ConstCnt",            &jetConstCnt            );
  }
  //MC/GEN jet vars in tupel tree only
  std::vector<int> *GjetIdx;
  std::vector<float> *GjetPt, *GjetEta, *GjetPhi, *GjetE; 
  if(useTupel){
    jetpp[4]->SetBranchAddress( "JetAk04GenJet",   &GjetIdx );
    jetpp[4]->SetBranchAddress( "GJetAk04Pt",      &GjetPt );
    jetpp[4]->SetBranchAddress( "GJetAk04Eta",     &GjetEta );
    jetpp[4]->SetBranchAddress( "GJetAk04Phi",     &GjetPhi );
    jetpp[4]->SetBranchAddress( "GJetAk04E",       &GjetE );
  }
  




  // event count from files
  UInt_t NEvents_jetAnalyzr=jetpp[0]->GetEntries();   // preskim event count from files
  UInt_t NEvents_skimAnalyzr=jetpp[2]->GetEntries();   // preskim event count from files
  bool jetEventsExist= (bool)NEvents_jetAnalyzr, skimEventsExist=(bool)NEvents_skimAnalyzr;
  bool EventCountsEqual= (NEvents_jetAnalyzr==NEvents_skimAnalyzr);
  if( !EventCountsEqual ||
      !jetEventsExist  || !skimEventsExist ){//if this fires, we're gonna segfault
    std::cout<<"jetAnalyzr # evts  = " << NEvents_jetAnalyzr  << std::endl;
    std::cout<<"skimAnalyzr # evts = " << NEvents_skimAnalyzr << std::endl; }
  assert( jetEventsExist && skimEventsExist && EventCountsEqual );//exit if the counts don't make sense, help avoid garbage results
  
  UInt_t NEvents_allFiles=NEvents_jetAnalyzr;   // preskim event count from files
  h_NEvents->SetBinContent(1, NEvents_allFiles);
  
  UInt_t NEvents_toRead=0;
  NEvents_toRead = NEvents_allFiles;
  
  std::cout<<"reading "<<NEvents_toRead<<" events"<<std::endl;   




   //  //readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0
  TF1* L2RelL3AbsJEC_arr[n_L2RelL3AbsJEC_etabins]={};
  if(fillMCJetJECQAHists){    
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
  }
  
  //  //readForests_JEC_v12_ak4PF_75X_data_Run2_v13
  TF1* L2RelL3AbsJEC_2_arr[n_L2RelL3AbsJEC_2_etabins]={};
  if(fillMCJetJECQAHists){    
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
  }
  
  int NJets_largeJECdiff=0;
  
  //begin event loop
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {
    
    if( debugMode     &&
        nEvt%1000==0     )std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if(nEvt%10000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;  
    
    
    jetpp[0]->GetEntry(nEvt);
    
    // grab vzweight
    double vzWeight=0.;
    if(doVzWeights){
      if(!(fabs(vz_F)>24.)) vzWeight=cpuVzWeight_poly(vz_F);
      //vzWeight=cpuVzWeight_gauss(vz_F);
      //vzWeight=cpuVzWeight_bins(vz_F);
    }
    
    double evtPthatWeight=-1.;  
    //int evtPthatBin=-1; float evtPthat=0.;  
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){       
      //evtPthatBin=i;
      //evtPthat=pthatbins[i];
      evtPthatWeight=pthatWeights[i];//in millibarns
      evtPthatWeight*=1e+06;//in nanobarns
    } //every event regardless of reco quality should have a pthat weight. no exceptions.
    if((evtPthatWeight<0.)){std::cout<<"something wrong. pthatweight. debug now."<<std::endl;assert(false);}
    
    double weight_eS=evtPthatWeight*vzWeight;                //std::cout<<"weight_eS is now:"<<weight_eS<<std::endl<<std::endl;    
    h_NEvents_read->Fill(0.); h_NEvents_read->Fill(1.,weight_eS);    
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0     || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0 )continue;
    //|| //puvertexFilter_I==0  	) continue;            
    h_NEvents_skimCut->Fill(0.);    h_NEvents_skimCut->Fill(1.,weight_eS);    
    
    if( fabs(vz_F)>24.              ) continue;   
    h_NEvents_vzCut->Fill(0.);    h_NEvents_vzCut->Fill(1.,weight_eS);
    
    
    
    //vz
    if(fillMCEvtQAHists){
      
      hpthat->Fill(pthat_F, 1.);
      hWpthat->Fill(pthat_F, weight_eS);
      
      float vr_F=std::sqrt(vx_F*vx_F + vy_F*vy_F);
      hVr->Fill(vr_F, 1.0);   
      hWVr->Fill(vr_F, weight_eS);  
      
      hVx->Fill(vx_F, 1.0);   
      hWVx->Fill(vx_F, weight_eS);  
      
      hVy->Fill(vy_F, 1.0);   
      hWVy->Fill(vy_F, weight_eS);

      hVz->Fill(vz_F, 1.);
      hvzWVz->Fill(vz_F, vzWeight);
      hpthatWVz->Fill(vz_F, evtPthatWeight);
      hWVz->Fill(vz_F, weight_eS);
      
      hNref->Fill(nref_I, 1.);
      hWNref->Fill(nref_I, weight_eS);
    }
    
    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;
    
    // for event counting + avoiding duplicate fills in dijet hists
    bool dijetHistsFilled=false;//, dijetHistsFilled_wJetID=false;
    
    int jetsPerEvent;//=nref_I;//, jetsPerEventJetID=nref_I;
    float evt_leadJetPt=-1., evt_leadJetPt_wCuts=-1.;
    
    if(useIncJetAnalyzer){
      jetsPerEvent=nref_I;
      
      for(int jet = 0; jet<nref_I; ++jet)      {
	
	float jtpt     = pt_F[jet];
	float rawpt=rawpt_F[jet];      
	float receta    = eta_F[jet];
	float absreceta = fabs(eta_F[jet]);
	float recphi    = phi_F[jet];
	
	float genpt   = refpt_F[jet];
	float genphi  = refphi_F[jet];
	float geneta  = refeta_F[jet];
	float gendrjt = refdrjt_F[jet];

	int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
	int neuMult = neN_I[jet] + phN_I[jet] ;
	int numConst  = chMult + neuMult;

	// event+jet counting
	h_NJets->Fill(0.);      h_NJets->Fill(1.,weight_eS);

	// largest jet pt in each event
	if(jtpt>evt_leadJetPt && absreceta<4.7)	evt_leadJetPt=jtpt;	     
	
	//if reco jet w/o matched gen jet, skip.
	if( subid_F[jet] != 0 ) { jetsPerEvent--;   continue; }
	h_NJets_matched->Fill(0.);      h_NJets_matched->Fill(1.,weight_eS);
	
	if( gendrjt > gendrCut ){	jetsPerEvent--;	continue;  }
	h_NJets_deltaRCut->Fill(0.);      h_NJets_deltaRCut->Fill(1.,weight_eS);
	
	if(!(absreceta<4.7)) {jetsPerEvent--;continue;}
	h_NJets_jtetaCut1->Fill(0.);      h_NJets_jtetaCut1->Fill(1.,weight_eS);
	
	if( absreceta  < jtEtaCutLo ) { jetsPerEvent--;	continue;}
	if( !(absreceta < jtEtaCutHi) ) { jetsPerEvent--;	continue;}
	h_NJets_jtetaCut2->Fill(0.);      h_NJets_jtetaCut2->Fill(1.,weight_eS);
	
	if( !(jtpt > jtPtCut)   )      { jetsPerEvent--;     continue; } //paranoia about float comparison
	h_NJets_jtptCut->Fill(0.);      h_NJets_jtptCut->Fill(1.,weight_eS);
	
	if( !(jtpt < jtPtCut_Hi)   )      { jetsPerEvent--;  continue; }           
	h_NJets_jtptCut_Hi->Fill(0.);      h_NJets_jtptCut_Hi->Fill(1.,weight_eS);	
	
	if( !(genpt > genJetPtCut)   )  { jetsPerEvent--;      continue; } //paranoia about float comparison
	h_NJets_genjtptCut->Fill(0.);      h_NJets_genjtptCut->Fill(1.,weight_eS);

	if( !(genpt < genJetPtCut_Hi)   )  { jetsPerEvent--;   continue; }
	h_NJets_genjtptCut_Hi->Fill(0.);      h_NJets_genjtptCut_Hi->Fill(1.,weight_eS);
				

	// readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0
	float L2RelL3Abs_rawpt=-1.;//L2RelL3Abs JEC
	int L2RelL3AbsJECetabin=-1;
	float L2L3Res_rawpt=-1.;//L2L3Res JEC WONT BE APPLIED (MC)
	//int L2L3ResJECetabin=-1;
	
	// readForests_JEC_v12_ak4PF_75X_data_Run2_v13
	float L2RelL3Abs_2_rawpt=-1.;//L2RelL3Abs JEC
	int L2RelL3Abs_2_JECetabin=-1;
	float L2L3Res_2_rawpt=-1.;//L2L3Res JEC WONT BE APPLIED (MC)
	//int L2L3Res_2_JECetabin=-1;
	
	if(fillMCJetJECQAHists){

	  //readForests_JEC_v12_ak4PF_74X_dataRun2_HLT_ppAt5TeV_v0
	  for(int uu=0; uu<1; uu++){//loop so that if the jet doesn't meet some criteria i can just quit while i'm ahead
	    for(int i=0; i<n_L2RelL3AbsJEC_etabins;i++){	    //L2RelL3Abs
	      if(receta>L2RelL3AbsJEC_etabins[i] && receta < L2RelL3AbsJEC_etabins[i+1]) { 
		L2RelL3AbsJECetabin=i; break;}
	      else continue;	    }
	    if(L2RelL3AbsJECetabin==-1)
	      continue;
	    
	    L2RelL3Abs_rawpt=rawpt*(L2RelL3AbsJEC_arr[L2RelL3AbsJECetabin]->Eval(rawpt));			    
	    
	    L2L3Res_rawpt=L2RelL3Abs_rawpt;
	  }//end loop for easy continuing
	  
	  //readForests_JEC_v12_ak4PF_75X_data_Run2_v13
	  for(int uu=0; uu<1; uu++){//loop so that if the jet doesn't meet some criteria i can just quit while i'm ahead
	    for(int i=0; i<n_L2RelL3AbsJEC_2_etabins;i++){	    //L2RelL3Abs
	      if(receta>L2RelL3AbsJEC_2_etabins[i] && receta < L2RelL3AbsJEC_2_etabins[i+1]) { 
		L2RelL3Abs_2_JECetabin=i; break;}
	      else continue;	    }
	    if(L2RelL3Abs_2_JECetabin==-1)
	      continue;
	    
	    L2RelL3Abs_2_rawpt=rawpt*(L2RelL3AbsJEC_2_arr[L2RelL3Abs_2_JECetabin]->Eval(rawpt));			    
	    L2L3Res_2_rawpt=L2RelL3Abs_2_rawpt;
	  }//end loop for easy continuing
	  
	}//end JEC QA
	
	// 13 TeV JetID criterion, loose or tight
	bool passesJetID=false;
	
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
	
	int theRapBin=-1;
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin){
	  if( absetabins[rapbin]<=absreceta  && 		
	      absreceta<absetabins[rapbin+1]    ) {	    
	    theRapBin=rapbin;	    
	    break;
	  }
	}//rap bin loop


	if(fillMCJetJECQAHists && passesJetID ){
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
	


	if(!(jtpt>jetQAPtCut)){ jetsPerEvent--; continue;}//need JEC syst. stuff before this cut; so jets from below 56 GeV can migrate above	
	h_NJets_jetQAPtCut->Fill(0.);	  h_NJets_jetQAPtCut->Fill(1., weight_eS);
	
	for(int jtid=0; jtid<2; jtid++){

	  if(jtid==1 ){
	    
	    if(!passesJetID){ jetsPerEvent--; continue;	  }
	    else {  h_NJets_JetIDCut->Fill(0.);  h_NJets_JetIDCut->Fill(1., weight_eS); }
	    
	    if(jtpt>evt_leadJetPt_wCuts)
	      evt_leadJetPt_wCuts=jtpt;      
	    
	  }
	  
	  /////   JETQA   ///// 
	  if(fillMCJetQAHists	 ){
	    
	    int ind=0;
	    //jets
	    hJetQA_jtptEntries[theRapBin][jtid]->Fill(jtpt,1.);	  	
	    hJetQA[theRapBin][ind][jtid]->Fill(jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(rawpt_F[jet], weight_eS); ind++;
	    //if(jtpt>jetQAPtCut){//second jet pt cut
	    
	    hJetQA[theRapBin][ind][jtid]->Fill(eta_F[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(phi_F[jet], weight_eS); ind++;
	    
	    //tracks
	    hJetQA[theRapBin][ind][jtid]->Fill(trkN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(trkSum_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(trkMax_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(trkHardN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(trkHardSum_F[jet]/jtpt, weight_eS); ind++;
	    
	    //PF photons
	    hJetQA[theRapBin][ind][jtid]->Fill(phN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(phSum_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(phMax_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(phHardN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(phHardSum_F[jet]/jtpt, weight_eS); ind++;
	    
	    //PF charged hadrons
	    hJetQA[theRapBin][ind][jtid]->Fill(chN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(chSum_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(chMax_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(chHardN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(chHardSum_F[jet]/jtpt, weight_eS); ind++;
	    
	    //PF neutral hadons
	    hJetQA[theRapBin][ind][jtid]->Fill(neN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(neSum_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(neMax_F[jet]/jtpt, weight_eS); ind++;
	    
	    //PF electrons
	    hJetQA[theRapBin][ind][jtid]->Fill(eN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(eSum_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(eMax_F[jet]/jtpt, weight_eS); ind++;
	    
	    //PF muons
	    hJetQA[theRapBin][ind][jtid]->Fill(muN_I[jet], weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(muSum_F[jet]/jtpt, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(muMax_F[jet]/jtpt, weight_eS); ind++;
	    
	    //PF particle sums
	    hJetQA[theRapBin][ind][jtid]->Fill(neuMult, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(chMult, weight_eS); ind++;
	    hJetQA[theRapBin][ind][jtid]->Fill(numConst, weight_eS); ind++;
	    
	    if(fillgenJetQA){
	      hMCJetQA_genpt[theRapBin][jtid]    ->Fill(genpt,weight_eS);
	      hMCJetQA_geneta[theRapBin][jtid]   ->Fill(geneta,weight_eS);
	      hMCJetQA_genphi[theRapBin][jtid]   ->Fill(genphi,weight_eS);
	      hMCJetQA_genrecpt[theRapBin][jtid] ->Fill(jtpt/genpt,weight_eS);
	      
	      hMCJetQA_genreceta[theRapBin][jtid]->Fill(fabs(receta-geneta), weight_eS);	
	      hMCJetQA_genrecphi[theRapBin][jtid]->Fill(fabs(recphi-genphi), weight_eS);	
	      
	      hMCJetQA_gendrjt[theRapBin][jtid]->Fill(gendrjt, weight_eS);
	    }//genjet qa
	    //}//jetQAcut
	    
	    //looking for the first two good jets that meet the criteria specified
	    if(fillMCDijetHists&&passesJetID){
	      if ( !firstGoodJetFound ){
		if(jtpt>ldJetPtCut) { firstGoodJetFound=true;
		  firstGoodJetPt =jtpt; 
		  firstGoodJetPhi=recphi; } }
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
		      secondGoodJetPhi = recphi; }}} }
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
		hJetQA[theRapBin][ind][jtid]->Fill( secondGoodJetPt , weight_eS ); ind++;
	      } }//mcdijet hists
	  }      //end fillMCJetQAHists	      
	}//end jetid loop
      }//end jet loop
    } //end use inc jet analyzer
      
    if(useTupel){
      int tupel_njets=jetPt->size();
      jetsPerEvent=tupel_njets;
      
      for(int jet=0; jet<(tupel_njets) ; jet++){
	
	float jtpt  = jetPt->at(jet); //this will be the jetpt pulled from the forest
	float receta = jetEta->at(jet);
	float recphi    = jetPhi->at(jet);
	float absreceta = fabs(receta);

	float jetNeuMult=jetConstCnt->at(jet)-jetChMult->at(jet);


	//all reco jets in reco jet branches are matched
	int gjet=GjetIdx->at(jet)-1;
	//std::cout<<"gjet="<<gjet<<std::endl;
	//std::cout<<"GjetPt->size()="<<GjetPt->size()<<std::endl;
	//assert(false);
	if(gjet<0)continue;
	h_NJets_matched->Fill(0.);      h_NJets_matched->Fill(1.,weight_eS);
	
	float genpt   = GjetPt->at(gjet);
	float genphi  = GjetPhi->at(gjet);
	float geneta  = GjetEta->at(gjet);
	float genenergy= GjetE->at(gjet);
	
	float gendrjt = sqrt(
			     pow(receta - geneta, 2) + 
			     pow(recphi - genphi, 2)
			     );
	
	if( gendrjt > gendrCut ){	jetsPerEvent--;	continue;  }
	h_NJets_deltaRCut->Fill(0.);      h_NJets_deltaRCut->Fill(1.,weight_eS);
	
	if(!(absreceta<4.7)) {jetsPerEvent--;continue;}
	h_NJets_jtetaCut1->Fill(0.);      h_NJets_jtetaCut1->Fill(1.,weight_eS);
	
	if( absreceta  < jtEtaCutLo ) { jetsPerEvent--;	continue;}
	if( !(absreceta < jtEtaCutHi) ) { jetsPerEvent--;	continue;}
	h_NJets_jtetaCut2->Fill(0.);      h_NJets_jtetaCut2->Fill(1.,weight_eS);
	
	if( !(jtpt > jtPtCut)   )      { jetsPerEvent--;     continue; } //paranoia about float comparison
	h_NJets_jtptCut->Fill(0.);      h_NJets_jtptCut->Fill(1.,weight_eS);
	
	if( !(jtpt < jtPtCut_Hi)   )      { jetsPerEvent--;  continue; }           
	h_NJets_jtptCut_Hi->Fill(0.);      h_NJets_jtptCut_Hi->Fill(1.,weight_eS);	
	
	if( !(genpt > genJetPtCut)   )  { jetsPerEvent--;      continue; } //paranoia about float comparison
	h_NJets_genjtptCut->Fill(0.);      h_NJets_genjtptCut->Fill(1.,weight_eS);

	if( !(genpt < genJetPtCut_Hi)   )  { jetsPerEvent--;   continue; }
	h_NJets_genjtptCut_Hi->Fill(0.);      h_NJets_genjtptCut_Hi->Fill(1.,weight_eS);
	

	if(jtpt>evt_leadJetPt && absreceta<4.7)	evt_leadJetPt=jtpt;	
	
	
	
	bool passesJetID=false;
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
	
	
	//if(!passesJetID)continue;
	if(passesJetID){
	  h_NJets_JetIDCut->Fill(0.);		h_NJets_JetIDCut->Fill(1.,weight_eS);	
	}

	int theRapBin=-1;
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin){
	  if( absetabins[rapbin]<=absreceta  && 		
	      absreceta<absetabins[rapbin+1]    ) {	    
	    theRapBin=rapbin;	    
	    break;
	  }
	}//rap bin loop
	
	for(int jtid=0; jtid<2; jtid++){

	  if(jtid==1 && !passesJetID)continue;
	  
	  
	  
	  
	  if(fillMCTupelJetQAHists){
	    
	    
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
	    
	    if(fillgenJetQA){
	      hTupelGJetQA_genpt[theRapBin][jtid]    ->Fill(genpt ,weight_eS);
	      hTupelGJetQA_geneta[theRapBin][jtid]   ->Fill(geneta ,weight_eS);
	      hTupelGJetQA_genphi[theRapBin][jtid]   ->Fill(genphi ,weight_eS);
	      hTupelGJetQA_genenergy[theRapBin][jtid]     ->Fill(genenergy ,weight_eS);
	      hTupelGJetQA_gendrjt[theRapBin][jtid] ->Fill(gendrjt,weight_eS);
	    }
	    
	  }//end fill data tupel jet qa hists	  
	}//end jtid loop
      }//end tupel jet loop
    }//end use Tupel
    

      
      
    if(fillMCEvtQAHists){
      hjetsPEvt ->Fill(jetsPerEvent,1.0);   
      hWjetsPEvt->Fill(jetsPerEvent,weight_eS); 
      if(evt_leadJetPt>0.)hLeadJetPt->Fill(evt_leadJetPt ,weight_eS);
      if(evt_leadJetPt_wCuts>0.)  hLeadJetPt_wCuts->Fill(evt_leadJetPt_wCuts ,weight_eS);
    } 
    
  }//end event loop
  
  

  
  std::cout<<std::endl;
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"/// Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"# of Events (criteria/description) = count / pthat*vz weighted count" <<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"# evts passing w cut applied = count / xsec weighted count" <<std::endl;
  std::cout<<"in file(s) opened            = " << NEvents_allFiles<<std::endl; 
  std::cout<<"looped over                  = " <<h_NEvents->GetBinContent(1)          << " / " << h_NEvents->GetBinContent(2)<<std::endl; 
  std::cout<<"read                         = " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;
  std::cout<<"skimCuts  = read - !skim     = " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;//same
  std::cout<<"vzCuts    = skimCuts - !vz   = " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl; //same
  std::cout << std::endl;


  std::cout<<"//////////////////////////////"<<std::endl;
  std::cout<<"/// Jet-Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////////"<<std::endl;

  std::cout<<std::endl;
  std::cout<<"# of Jets w/ large corrected jet p_T diff = "<< NJets_largeJECdiff << std::endl;
  std::cout<<"# of Jets read                     = " << h_NJets->GetBinContent(1) << " / " << h_NJets->GetBinContent(2) <<std::endl;
  std::cout<<"# of GEN-RECO Matched Jets         = " << h_NJets_matched->GetBinContent(1) << " / " << h_NJets_matched->GetBinContent(2) <<std::endl;
  std::cout<<"# of Matched Jets meeting dR cut   = " << h_NJets_deltaRCut->GetBinContent(1) << " / " << h_NJets_deltaRCut->GetBinContent(2) <<std::endl;  
  std::cout<<"# of Jets w/ jtpt > jtptCut        = " << h_NJets_jtptCut->GetBinContent(1) << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ jtpt > jtptCut_Hi     = " << h_NJets_jtptCut_Hi->GetBinContent(1) << " / " << h_NJets_jtptCut_Hi->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ |jteta|<4.7           = " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ passing jtEtaCut(s)   = " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ genjtpt>genjtptCut    = " << h_NJets_genjtptCut->GetBinContent(1) << " / " << h_NJets_genjtptCut->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ genjtpt<genjtptCut_Hi = " << h_NJets_genjtptCut_Hi->GetBinContent(1) << " / " << h_NJets_genjtptCut_Hi->GetBinContent(2) <<std::endl;  
  std::cout<<"# of Jets w/ jetIDCut              = " << h_NJets_JetIDCut->GetBinContent(1) << " / " << h_NJets_JetIDCut->GetBinContent(2) <<std::endl;   
  
  std::cout<<"writing output file... ";
  fout->Write(); 
  
  std::cout<<std::endl<<"readForests_ppMC_jetPlots finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real tyime (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  std::cout<<"output file written is  "<<outfile<<std::endl<<std::endl;
  
  return 0;
}// end readForests_ppMC_jetPlots
								 

////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1 && argc!=10){
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppMC_jetPlots.exe ";
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
  if(argc==1) rStatus = readForests_ppMC_jetPlots();
  else if (argc==8){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppMC_jetPlots( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
				  outputFileName);
  }
  else if (argc==10){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
    rStatus = readForests_ppMC_jetPlots( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
					   outputFileName, jtEtaCutLo, jtEtaCutHi);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
