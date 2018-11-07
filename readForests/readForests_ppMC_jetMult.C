// Custom header
#include "readForests_jetMult.h"

// ppMC switches
const bool fillMCEvtQAHists=true;
const bool fillMCJetQAHists=true;
const bool fillgenJetQA=true&fillMCJetQAHists;
const bool fillMCJetIDHists=true;//, tightJetID=false;
const bool fillMCJetMultHists=true;


const bool fillMCDijetHists=false;

const int jetIDint=(int)fillMCJetIDHists; //might not use this?

//// readForests_ppMC_jetMult
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_jetMult(std::string inFilelist , int startfile , int endfile ,  
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
  std::cout<<", debugMode = "<<debugMode<<std::endl;
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
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0)jetpp[0]->AddFriend( jetpp[t] );  }
  
  
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
    for(int t = 0;t<N_MCTrees;++t) filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
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
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"     ,(std::to_string(jtPtCut)).c_str()   ,    100, 0,100); hJetPtCut->Fill(jtPtCut);           
  TH1F *hJetPtCut_Hi        =new TH1F("hJetPtCut_Hi"     ,(std::to_string(jtPtCut_Hi)).c_str()   ,    1000, 50,1500); hJetPtCut_Hi->Fill(jtPtCut_Hi);           
  TH1F *hGenJetPtCut     =new TH1F("hGenJetPtCut"  ,(std::to_string(genJetPtCut)).c_str()   ,    100, 0,100); hGenJetPtCut->Fill(genJetPtCut);           
  TH1F *hGenJetPtCut_Hi     =new TH1F("hGenJetPtCut_Hi"  ,(std::to_string(genJetPtCut_Hi)).c_str()   ,    1000, 1000,2000); hGenJetPtCut_Hi->Fill(genJetPtCut_Hi);           
  TH1F *hJetEtaCutHi     =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6  ); hJetEtaCutHi->Fill(jtEtaCutHi);	        
  TH1F *hJetEtaCutLo     =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6  ); hJetEtaCutLo->Fill(jtEtaCutLo);	        
  TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);  
  TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  


  //evt conts
  TH1F *h_NEvents          = new TH1F("NEvents"         , "NEvents",         1,0,2);
  TH1F *h_NEvents_read     = new TH1F("NEvents_read"    , "NEvents read",    1,0,2);
  //TH1F *h_NEvents_skipped  = new TH1F("NEvents_skipped" , "NEvents skipped", 1,0,2);
  //TH1F *h_NEvents_trigd    = new TH1F("NEvents_trigd"   , "NEvents trigd",   1,0,2);
  //TH1F *h_NEvents_jet40    = new TH1F("NEvents_jet40Trigd"   , "NEvents jet40Trigd",    1,0,2);
  //TH1F *h_NEvents_jet60    = new TH1F("NEvents_jet60Trigd"   , "NEvents jet60Trigd",    1,0,2);
  //TH1F *h_NEvents_jet80    = new TH1F("NEvents_jet80Trigd"   , "NEvents jet80Trigd",    1,0,2);
  //TH1F *h_NEvents_jet100   = new TH1F("NEvents_jet100Trigd"  , "NEvents jet100Trigd",   1,0,2);
  TH1F *h_NEvents_skimCut  = new TH1F("NEvents_skimCut" , "NEvents read post skimCut",  1,0,2);
  TH1F *h_NEvents_vzCut    = new TH1F("NEvents_vzCut"   , "NEvents read post vzCut AND skimCut", 1,0,2);

  TH1F *h_NEvents_trkCuts_1    = new TH1F("NEvents_trkCuts_1"   , "NEvents read post prev cuts AND trkHiPur cut", 1,0,2);
  TH1F *h_NEvents_trkCuts_2    = new TH1F("NEvents_trkCuts_2"   , "NEvents read post prev cuts AND Vtxs spaced out cut", 1,0,2);
  TH1F *h_NEvents_trkCuts_3    = new TH1F("NEvents_trkCuts_3"   , "NEvents read post prev cuts AND oneVtxInRange cut", 1,0,2);    


  //counts post all evt cuts, that are thrown out because of no good jets pre and post jet cuts
  TH1F *h_NEvents_withJets           = new TH1F("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut   = new TH1F("NEvents_withJets_kmatCut" , 
						"NEvents read post evt cuts, w/ jets post kmatCut", 1,0,2);
  TH1F *h_NEvents_withJets_JetIDCut = NULL;
  if(fillMCJetIDHists)
    h_NEvents_withJets_JetIDCut = new TH1F("NEvents_withJets_JetIDCut" , 
					   "NEvents read post evt cuts, w/ jets post kmatCut AND JetID Cut", 1,0,2);
  
  /////   EVENT QA HISTS   ///// 
  TH1F *hVr=NULL, *hWVr=NULL;
  TH1F *hVx=NULL, *hWVx=NULL;
  TH1F *hVy=NULL, *hWVy=NULL;

  TH1F *hVz=NULL, *hWVz=NULL , *hvzWVz=NULL,*hpthatWVz=NULL;
  
  TH1F *hpthat=NULL, *hWpthat=NULL;  

  TH1F* hNref=NULL, *hWNref=NULL;
  TH1F* hjetsPEvt=NULL,* hWjetsPEvt=NULL;
  //TH1F* hjetsPEvt_wJetID=NULL,* hWjetsPEvt_wJetID=NULL;
  TH1F* hLeadJetPt=NULL;//,* hLeadJetPt_wJetID=NULL;

  if(fillMCEvtQAHists){

    hpthat    = new TH1F("hpthat","",1000,0,1000);//evt pthat, unweighted and weighted
    hWpthat   = new TH1F("hWeightedpthat","",1000,0,1000);  

    hVr = new TH1F("hVr","vr, no trig, no weights",         2000,0.,0.60); 
    hWVr = new TH1F("hWeightedVr","vr, trigd, with weights",2000,0.,0.60); 
    
    hVx = new TH1F("hVx","vx, no trig, no weights",  2000, -0.30,0.30);
    hWVx = new TH1F("hWeightedVx","vx, trigd, with weights", 2000, -0.30,0.30); 

    hVy = new TH1F("hVy","vy, no trig, no weights",   2000, -0.30,0.30); 
    hWVy = new TH1F("hWeightedVy","vy, trigd, with weights",   2000,-0.30,0.30);
    
    hVz       = new TH1F("hVz","", 1000,-25.,25.);//evtvz
    hWVz      = new TH1F("hWeightedVz","", 1000,-25.,25.);//pthat*vz-weighted evt vz

    hvzWVz    = new TH1F("hvzWeightedVz","", 100,-25.,25.);//vz-weighted evtvz
    hpthatWVz = new TH1F("hpthatWeightedVz","", 100,-25.,25.);//pthat-weighted evtvz
    
    hNref = new TH1F("hNref","numJets each evt",30,0,30);
    hWNref = new TH1F("hWNref","weighted numJets each evt",30,0,30);
    hjetsPEvt  = new TH1F("hjetsPEvt","numJets each evt",30,0,30);
    hWjetsPEvt = new TH1F("hWjetsPEvt","weighted numJets each evt",30,0,30);
    
//    if(fillMCJetIDHists){
//      hjetsPEvt_wJetID = new TH1F("hjetsPEvt_wJetID","num JetID Jets each evt",30,0,30);
//      hWjetsPEvt_wJetID = new TH1F("hWjetsPEvt_wJetID","weighted num JetID Jets each evt",30,0,30);}
    
    hLeadJetPt      = new TH1F("hLeadJetPt","lead jet pt for each evt",2000,0,2000);  
    //if(fillMCJetIDHists)
    //  hLeadJetPt_wJetID      = new TH1F("hLeadJetPt_wJetId","lead jet pt w Jet ID for each evt",1000,0,1000);
    
  }
  
  
  
  //jet counts
  TH1F *h_NJets          = new TH1F("NJets","NJets read", 1,0,2);
  TH1F *h_NJets_kmatCut  = new TH1F("NJets_kmatCut ","NJets read post kmatCut ", 1,0,2);
  TH1F *h_NJets_JetIDCut = NULL;
  if(fillMCJetIDHists)
    h_NJets_JetIDCut = new TH1F("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 1,0,2);


  //weighted jet counts
  TH1F *h_WNJets          = new TH1F("WNJets","WNJets read", 1,0,2);
  TH1F *h_WNJets_kmatCut  = new TH1F("WNJets_kmatCut ","WNJets read post kmatCut ", 1,0,2);
  TH1F *h_WNJets_JetIDCut=NULL;  
  if(fillMCJetIDHists)
    h_WNJets_JetIDCut = new TH1F("WNJets_JetIDCut","WNJets read post JetIDCut AND kmatCut", 1,0,2);
  


  /////   JET QA HISTS   ///// 
  //jet QA
  TH1F *hJetQA[2][N_vars]={};
  //MC JetQA
  TH1F *hMCJetQA_genpt[2]={}, *hMCJetQA_geneta[2]={};
  TH1F *hMCJetQA_genrecpt[2]={}, *hMCJetQA_genreceta[2]={};
  if(fillMCJetQAHists){
    for(int k = 0; k<2; ++k){

      if(!fillMCJetIDHists && k==1)continue;
      if( fillMCJetIDHists && k==0)continue;	
      
      for(int j = 0; j<N_vars; ++j){	
	//binnings+hist declarations, etc.
	
        //jets
        if(var[j]=="jtpt"||var[j]=="rawpt")
          hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;", var[j].c_str()) , 2000,0,2000);
        else if(var[j]=="jteta")
          hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-5,+5);
        else if(var[j]=="jtphi")
          hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-4,+4);
	
        //jetconst. counts
        else if(var[j]=="trkN"|| var[j]=="phN"|| var[j]=="chN"|| var[j]=="neN"|| var[j]=="eN"|| var[j]=="muN")
          hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 60,0,60);
        else if(var[j]=="trkHardN"|| var[j]=="phHardN"|| var[j]=="chHardN")
          hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 60,0,60);
        else if(var[j]=="neuMult"|| var[j]=="chMult"|| var[j]=="numConst")
          hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 100,0,100);

	//dijets
	else if (fillMCDijetHists){
	  if(var[j]=="dphi") 
	    hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 50,0,+4);
	  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt")
	    hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;", var[j].c_str()), 2000,0,2000);
	  else //xj and Aj binnings
	    hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);         
	}
	
	//consituent binnings
	else 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);         
      }
      
      if(fillgenJetQA){	
	hMCJetQA_genpt[k]    = new TH1F(Form("hJetQA_%dwJetID_genpt",k)     , "genpt",    2000, 0,2000 );
	hMCJetQA_geneta[k]   = new TH1F(Form("hJetQA_%dwJetID_geneta",k)    , "geneta", 200,  -6.,6. );
	hMCJetQA_genrecpt[k] = new TH1F(Form("hJetQA_%dwJetID_genrecpt",k)  , "genrecpt",    100, 0, 2. );
	hMCJetQA_genreceta[k]= new TH1F(Form("hJetQA_%dwJetID_genreceta",k) , "gen-receta", 20000,  -100, 100 ); }
      
    }}
  
  // jetMult hists, explore how yields/etc change w/ multiplicity of collision(s)

  // hists for counts
  TH1F *h_NTrks        =NULL;
  TH1F *h_NTrks_kmatCut=NULL;

  // hists filled once per event
  TH1F *Nvtxs=NULL;
  TH1F *Ntrks=NULL,  *Ntrks_passCuts=NULL;  //# trks in event
  //TH1F *Ntrks_vtxs=NULL,  *Ntrks_vtxs_passCuts=NULL;  // # trks assoc w/ vtx in event

  TH1F *vtx_maxPt=NULL, *vtx_maxMult=NULL, *vtx_minNormChi2=NULL;
  TH2F *vtx_maxPt_v_maxMult=NULL;

  TH2F *Nvtxs_v_Nref=NULL, *Ntrks_v_Nref=NULL, *Ntrks_v_Nvtxs=NULL;

  // vtx hists

  // x/y/z/cylindrical r
  TH1F *vtx_z=NULL,*vtx_y=NULL,*vtx_x=NULL;
  TH1F *vtx_r=NULL;
  TH2F *vtx_rVz=NULL  , *vtx_xVy=NULL  ;

  // ratios/steps used in eta calculation
  TH1F *vtx_yOverx=NULL, *vtx_rOverz=NULL,*vtx_zOverrad=NULL;
  TH1F *vtx_theta=NULL, *vtx_tantheta2, *vtx_lntantheta2;

  // spherical r/phi/eta
  TH1F *vtx_rad=NULL,*vtx_phi=NULL, *vtx_eta=NULL;


  TH1F *Nvtxtrks=NULL;
  TH2F *NvtxtrksVz=NULL, *NvtxtrksVr=NULL;

  TH1F *normChi2_vtx=NULL,*sumPtVtx_vtx=NULL;

  // vtxgt0 hists
  TH1F *Nvtxtrks_vtxgt0=NULL;
  TH2F *NvtxtrksVz_vtxgt0=NULL, *NvtxtrksVr_vtxgt0=NULL;

  TH1F *normChi2_vtxgt0=NULL,*sumPtVtx_vtxgt0=NULL;

  // vtx0 hists
  TH1F *Nvtxtrks_vtx0=NULL;
  TH2F *NvtxtrksVz_vtx0=NULL, *NvtxtrksVr_vtx0=NULL;

  TH1F *normChi2_vtx0=NULL,*sumPtVtx_vtx0=NULL;

  //trk hists
  TH1F *trkPt=NULL,  *trkPtError=NULL, *trkPtPercError=NULL;
  TH1F *trkEta=NULL,*trkPhi=NULL;

  TH1F *trkDxy=NULL, *trkDxyError=NULL, *trkDxySig=NULL;
  TH1F *trkDz=NULL , *trkDzError=NULL , *trkDzSig=NULL;

  TH1F *trkCharge=NULL;
  TH1F *trkNVtx=NULL,*trkChi2=NULL,*trkNdof=NULL;

  TH1F *trkNHit=NULL, *trkNLayer=NULL;

  TH1F *trkHighPur=NULL,*trkTight=NULL, *trkLoose=NULL,*trkFake=NULL;

  TH2F *trkPtvEta=NULL, *trkEtavPhi=NULL;


  if(fillMCJetMultHists){

    std::string jetID=std::to_string((int)fillMCJetIDHists);

    h_NTrks = new TH1F("NTrks","NTrks read", 1,0,2);
    h_NTrks_kmatCut = new TH1F("NTrks_kmatCut ","NTrks read post kmatCut ", 1,0,2);

    // hists filled once per event
    Nvtxs    = new TH1F(Form("hJetMultQA_%swJetID_Nvtx"    ,jetID.c_str()), "Nvtxs", 20,0.,20.);

    Ntrks    = new TH1F(Form("hJetMultQA_%swJetID_Ntrks"   ,jetID.c_str()), "Ntrks", 300,0.,300.);
    Ntrks_passCuts    = new TH1F(Form("hJetMultQA_%swJetID_Ntrks_passCuts"   ,jetID.c_str()), "Ntrks_passCuts", 300,0.,300.);

    //Ntrks_vtxs    = new TH1F(Form("hJetMultQA_%swJetID_Ntrks_vtxs"   ,jetID.c_str()), "Ntrks_vtxs", 300,0.,300.);
    //Ntrks_vtxs_passCuts    = new TH1F(Form("hJetMultQA_%swJetID_Ntrks_vtxs_passCuts"   ,jetID.c_str()), "Ntrks_vtxs_passCuts", 300,0.,300.);

    vtx_maxPt   = new TH1F(Form("hJetMultQA_%swJetID_vtx_maxPt"  , jetID.c_str()) , "vtx_maxPt" , 10 , 0 , 10  );
    vtx_maxMult   = new TH1F(Form("hJetMultQA_%swJetID_vtx_maxMult"  , jetID.c_str()) , "vtx_maxMult" , 10 , 0 , 10  );
    vtx_minNormChi2   = new TH1F(Form("hJetMultQA_%swJetID_vtx_minNormChi2"  , jetID.c_str()) , "vtx_minNormChi2" , 10, 0, 10  );

    vtx_maxPt_v_maxMult   = new TH2F ( Form("hJetMultQA_%swJetID_maxPt_v_maxMult", jetID.c_str()),  "maxPt_v_maxMult",  10,0,10,     10,0,10 );

    Nvtxs_v_Nref = new TH2F ( Form("hJetMultQA_%swJetID_Nvtxs_v_Nref", jetID.c_str()),  "nvtxs_v_nref",  20,0,20,     20,0,20 );
    Ntrks_v_Nref = new TH2F ( Form("hJetMultQA_%swJetID_Ntrks_v_Nref", jetID.c_str()),  "ntrks_v_nref",  300,0,300,   20,0,20 );
    Ntrks_v_Nvtxs= new TH2F ( Form("hJetMultQA_%swJetID_Ntrks_v_Nvtxs", jetID.c_str()), "ntrks_v_nvtxs",  300,0,300,   20,0,20 );

    // vtx hists
    // x/y/z/cylindrical r
    vtx_z   = new TH1F(Form("hJetMultQA_%swJetID_vtx_z"  , jetID.c_str()) , "vtx_z" , 2000 , -25. , 25.  );
    vtx_y   = new TH1F(Form("hJetMultQA_%swJetID_vtx_y"  , jetID.c_str()) , "vtx_y" , 2000 , -0.4 , 0.4  );
    vtx_x   = new TH1F(Form("hJetMultQA_%swJetID_vtx_x"  , jetID.c_str()) , "vtx_x" , 2000 , -0.4 , 0.4  );
    vtx_r   = new TH1F(Form("hJetMultQA_%swJetID_vtx_r"  , jetID.c_str()) , "vtx_r" , 2000 ,   0. , 0.60 );
  
    vtx_rVz   = new TH2F(Form("hJetMultQA_%swJetID_vtx_rVz"    ,jetID.c_str()), "vtx_rVz", 2000,   0., 0.60, 2000 , -25. , 25. );
    vtx_xVy   = new TH2F(Form("hJetMultQA_%swJetID_vtx_xVy"    ,jetID.c_str()), "vtx_xVy", 2000, -0.4, 0.4 , 2000 , -0.4 , 0.4 );

    // ratios/steps used in eta calculation of vtxs
    vtx_yOverx = new TH1F(Form("hJetMultQA_%swJetID_vtx_yOverx"      ,jetID.c_str()), "vtx_yOverx"  , 2000 ,  -100., 100.);
    vtx_rOverz = new TH1F(Form("hJetMultQA_%swJetID_vtx_rOverz"      ,jetID.c_str()), "vtx_rOverz"    , 2000,  -100., 100.);
    vtx_zOverrad = new TH1F(Form("hJetMultQA_%swJetID_vtx_zOverrad"  ,jetID.c_str()), "vtx_zOverrad"  , 200 ,  -10., 10.);

    vtx_theta = new TH1F(Form("hJetMultQA_%swJetID_vtx_theta"      ,jetID.c_str()), "vtx_theta", 80, 0., 4.);
    vtx_tantheta2 = new TH1F(Form("hJetMultQA_%swJetID_vtx_tantheta2"      ,jetID.c_str()), "vtx_tantheta2", 600, 0., 600.);
    vtx_lntantheta2 = new TH1F(Form("hJetMultQA_%swJetID_vtx_lntantheta2"      ,jetID.c_str()), "vtx_lntantheta2", 200, -10., 10.);

    // spherical r (aka rad)/theta/phi, with pseudorapidity
    vtx_rad = new TH1F(Form("hJetMultQA_%swJetID_vtx_rad"      ,jetID.c_str()), "vtx_rad", 300, 0., 30.);
    vtx_phi   = new TH1F(Form("hJetMultQA_%swJetID_vtx_phi"      ,jetID.c_str()), "vtx_phi", 160, -4., 4.);
    vtx_eta = new TH1F(Form("hJetMultQA_%swJetID_vtx_eta"      ,jetID.c_str()), "vtx_eta", 200, -10., 10.);
  

    // other vtx propoerties
    Nvtxtrks = new TH1F(Form("hJetMultQA_%swJetID_Nvtxtrks",jetID.c_str()), "NVtxTrks", 300,0.,300.);
    NvtxtrksVz = new TH2F(Form("hJetMultQA_%swJetID_NvtxtrksVz",jetID.c_str()), "NVtxTrksVz", 300,0,300, 2000,-25.,25.);
    NvtxtrksVr = new TH2F(Form("hJetMultQA_%swJetID_NvtxtrksVr",jetID.c_str()), "NVtxTrksVr", 300,0,300, 2000,0.,0.60);

    normChi2_vtx = new TH1F(Form("hJetMultQA_%swJetID_normChi2_vtx",jetID.c_str()), "normChi2_vtx", 100, 0., 10. );
    sumPtVtx_vtx = new TH1F(Form("hJetMultQA_%swJetID_sumPtVtx_vtx",jetID.c_str()), "sumPtVtx_vtx", 700 , 0., 700.);

    // vtxgt0 hists
    Nvtxtrks_vtxgt0 = new TH1F(Form("hJetMultQA_%swJetID_Nvtxtrks_vtxgt0",jetID.c_str()), "NVtxTrks_vtxgt0", 300,0.,300.);
    NvtxtrksVz_vtxgt0 = new TH2F(Form("hJetMultQA_%swJetID_NvtxtrksVz_vtxgt0",jetID.c_str()), "NVtxTrksVz_vtxgt0", 300,0,300, 2000,-25.,25.);
    NvtxtrksVr_vtxgt0 = new TH2F(Form("hJetMultQA_%swJetID_NvtxtrksVr_vtxgt0",jetID.c_str()), "NVtxTrksVr_vtxgt0", 300,0,300, 2000,0.,0.60);

    normChi2_vtxgt0 = new TH1F(Form("hJetMultQA_%swJetID_normChi2_vtxgt0",jetID.c_str()), "normChi2_vtxgt0", 100, 0., 10. );
    sumPtVtx_vtxgt0 = new TH1F(Form("hJetMultQA_%swJetID_sumPtVtx_vtxgt0",jetID.c_str()), "sumPtVtx_vtxgt0", 700, 0., 700.);

    // vtx0 hists
    Nvtxtrks_vtx0 = new TH1F(Form("hJetMultQA_%swJetID_Nvtxtrks_vtx0",jetID.c_str()), "NVtxTrks_vtx0", 300,0.,300.);
    NvtxtrksVz_vtx0 = new TH2F(Form("hJetMultQA_%swJetID_NvtxtrksVz_vtx0",jetID.c_str()), "NVtxTrksVz_vtx0", 300,0,300, 2000,-25.,25.);
    NvtxtrksVr_vtx0 = new TH2F(Form("hJetMultQA_%swJetID_NvtxtrksVr_vtx0",jetID.c_str()), "NVtxTrksVr_vtx0", 300,0,300, 2000,0.,0.60);

    normChi2_vtx0 = new TH1F(Form("hJetMultQA_%swJetID_normChi2_vtx0",jetID.c_str()), "normChi2_vtx0", 100, 0., 10. );
    sumPtVtx_vtx0 = new TH1F(Form("hJetMultQA_%swJetID_sumPtVtx_vtx0",jetID.c_str()), "sumPtVtx_vtx0", 700, 0., 700.);

    // trk hists
    trkPt    = new TH1F(Form("hJetMultQA_%swJetID_trkPt"   ,jetID.c_str()), "trkPt", 1000,0.,100.);
    trkPtError    = new TH1F(Form("hJetMultQA_%swJetID_trkPtError"   ,jetID.c_str()), "trkPtError", 1000,0.,100.);
    trkPtPercError    = new TH1F(Form("hJetMultQA_%swJetID_trkPtPercError"   ,jetID.c_str()), "trkPtPercError", 200,0.,200.);

    trkEta    = new TH1F(Form("hJetMultQA_%swJetID_trkEta"   ,jetID.c_str()), "trkEta", 100,-5.,5.);
    trkPhi    = new TH1F(Form("hJetMultQA_%swJetID_trkPhi"   ,jetID.c_str()), "trkPhi", 100,-5.,5.);

    trkDz    = new TH1F(Form("hJetMultQA_%swJetID_trkDz"   ,jetID.c_str()), "trkDz", 400,-20.,20.);
    trkDzError    = new TH1F(Form("hJetMultQA_%swJetID_trkDzError"   ,jetID.c_str()), "trkDzError", 200,0.,20.);
    trkDzSig    = new TH1F(Form("hJetMultQA_%swJetID_trkDzSig"   ,jetID.c_str()), "trkDzSig", 100, -5., 5.);

    trkDxy    = new TH1F(Form("hJetMultQA_%swJetID_trkDxy"   ,jetID.c_str()), "trkDxy", 400,-20.,20.);
    trkDxyError    = new TH1F(Form("hJetMultQA_%swJetID_trkDxyError"   ,jetID.c_str()), "trkDxyError", 200,0.,20.);
    trkDxySig    = new TH1F(Form("hJetMultQA_%swJetID_trkDxySig"   ,jetID.c_str()), "trkDxySig", 100, -5., 5.);

    trkCharge    = new TH1F(Form("hJetMultQA_%swJetID_trkCharge"   ,jetID.c_str()), "trkCharge", 5,-2.5,2.5);

    trkNVtx    = new TH1F(Form("hJetMultQA_%swJetID_trkNVtx"   ,jetID.c_str()), "trkNVtx", 11,-1.,10.);
    trkChi2    = new TH1F(Form("hJetMultQA_%swJetID_trkChi2"   ,jetID.c_str()), "trkChi2", 2000,0.,2000.);
    trkNdof    = new TH1F(Form("hJetMultQA_%swJetID_trkNdof"   ,jetID.c_str()), "trkNdof", 100,0.,100.);

    trkNHit    = new TH1F(Form("hJetMultQA_%swJetID_trkNHit"   ,jetID.c_str()), "trkNHit", 40,0.,40.);
    trkNLayer    = new TH1F(Form("hJetMultQA_%swJetID_trkNLayer"   ,jetID.c_str()), "trkNLayer", 40,0.,40.);

    trkHighPur    = new TH1F(Form("hJetMultQA_%swJetID_trkHighPur"   ,jetID.c_str()), "trkHighPur", 5,-2.,3.);

    trkTight    = new TH1F(Form("hJetMultQA_%swJetID_trkTight"   ,jetID.c_str()), "trkTight", 5,-2.,3.);
    trkLoose    = new TH1F(Form("hJetMultQA_%swJetID_trkLoose"   ,jetID.c_str()), "trkLoose", 5,-2.,3.);

    trkFake    = new TH1F(Form("hJetMultQA_%swJetID_trkFake"   ,jetID.c_str()), "trkFake", 5,-2.,3.);

    trkPtvEta    = new TH2F(Form("hJetMultQA_%swJetID_trkPtvEta"   ,jetID.c_str()), "trkPtvEta", 1000,0,100, 100,-5.,5.);
    trkEtavPhi    = new TH2F(Form("hJetMultQA_%swJetID_trkEtaVPhi"   ,jetID.c_str()), "trkEtavPhi", 100,-5.,5.,100,-5.,5.);
  }
  
  
  
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
  // data written a few different ways in this tree;
  // 1) by-event info (i.e. # of tracks, # of vertices)
  // 2) by-vertex info (i.e. vertex positions, # of tracks assoc to vertex)
  // 3) by-track into (i.e. track pt, eta phi)
  // 4) by-track-vertex combination (i.e. is-track-associated, transverse/longitudinal minimum approach, etc.)

  // evt //
  int nVtx_I;
  int nTrk_I;
  int nTrkTimesnVtx_I;
  int maxPtVtx_I, maxMultVtx_I;

  jetpp[3]->SetBranchAddress("nVtx",&nVtx_I);
  jetpp[3]->SetBranchAddress("nTrk"   ,   &nTrk_I  );
  jetpp[3]->SetBranchAddress("nTrkTimesnVtx"   ,   &nTrkTimesnVtx_I  );
  jetpp[3]->SetBranchAddress("maxPtVtx",&maxPtVtx_I);
  jetpp[3]->SetBranchAddress("maxMultVtx",&maxMultVtx_I);

  // vtxs //
  int nTrkVtx_I[1000];
  float normChi2Vtx_F[1000];
  float sumPtVtx_F[1000];

  float xVtx_F[1000];
  float yVtx_F[1000];
  float zVtx_F[1000];

  jetpp[3]->SetBranchAddress("nTrkVtx",&nTrkVtx_I);
  jetpp[3]->SetBranchAddress("normChi2Vtx",&normChi2Vtx_F);
  jetpp[3]->SetBranchAddress("sumPtVtx",&sumPtVtx_F);

  jetpp[3]->SetBranchAddress("xVtx",&xVtx_F);
  jetpp[3]->SetBranchAddress("yVtx",&yVtx_F);
  jetpp[3]->SetBranchAddress("zVtx",&zVtx_F);
  
  // trks //
  float trkPt_F[10000];
  float trkPtError_F[10000];
  float trkEta_F[10000];
  float trkPhi_F[10000];

  float trkDxy_F     [10000];
  float trkDxyError_F[10000];
  float trkDz_F      [10000];
  float trkDzError_F [10000];

  int   trkCharge_I [10000];

  UChar_t   trkNVtx_I   [10000];

  float trkChi2_F   [10000];
  UChar_t          trkNdof_I   [10000];
  UChar_t          trkNHit_I   [10000];
  UChar_t          trkNLayer_I [10000];

  bool   trkHighPur_O[10000];
  bool   trkTight_O     [10000];
  bool   trkLoose_O     [10000];
  bool   trkFake_O   [10000];
  
  jetpp[3]->SetBranchAddress("trkPt"   ,   &trkPt_F  );
  jetpp[3]->SetBranchAddress("trkPtError"   ,   &trkPtError_F  );
  jetpp[3]->SetBranchAddress("trkEta"   ,   &trkEta_F  );
  jetpp[3]->SetBranchAddress("trkPhi"   ,   &trkPhi_F  );

  jetpp[3]->SetBranchAddress("trkDxy1"          , &trkDxy_F            );
  jetpp[3]->SetBranchAddress("trkDxyError1"     , &trkDxyError_F       );
  jetpp[3]->SetBranchAddress("trkDz1"          , &trkDz_F            );
  jetpp[3]->SetBranchAddress("trkDzError1"     , &trkDzError_F       );


  jetpp[3]->SetBranchAddress("trkCharge"   ,   &trkCharge_I  );
  jetpp[3]->SetBranchAddress("trkNVtx"   ,   &trkNVtx_I  );

  jetpp[3]->SetBranchAddress("trkChi2"   ,   &trkChi2_F  );
  jetpp[3]->SetBranchAddress("trkNdof"   ,   &trkNdof_I  );
  jetpp[3]->SetBranchAddress("trkNHit"   ,   &trkNHit_I  );
  jetpp[3]->SetBranchAddress("trkNlayer"   ,   &trkNLayer_I  );

  jetpp[3]->SetBranchAddress("highPurity", &trkHighPur_O );
  jetpp[3]->SetBranchAddress("tight", &trkTight_O );
  jetpp[3]->SetBranchAddress("loose", &trkLoose_O );
  jetpp[3]->SetBranchAddress("trkFake"   ,   &trkFake_O  );
  
  // MC only, tracks matched to reco/PF particles are then matched to up to 5 gen particles.
  int   matchedGenID_I[10000][5];
  jetpp[3]->SetBranchAddress("matchedGenID"   ,   &matchedGenID_I  );

  




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
  for(UInt_t i=0;i < NEvents_allFiles; ++i) h_NEvents->Fill(1);
  
  UInt_t NEvents_toRead=0;
  NEvents_toRead = NEvents_allFiles;
  
  std::cout<<"reading "<<NEvents_toRead<<" events"<<std::endl;   
  
  //float jetIDCut_neSum, jetIDCut_phSum;
  //if(tightJetID){     jetIDCut_neSum=0.90;  jetIDCut_phSum=0.90;}
  //else{     jetIDCut_neSum=0.99;  jetIDCut_phSum=0.99;}
  
  //begin event loop
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {
    
    if( debugMode     &&
        nEvt%1000==0     )std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if(nEvt%10000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;  
    
    jetpp[0]->GetEntry(nEvt);
    h_NEvents_read->Fill(1);
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0     || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0 )continue;
    //||
    //puvertexFilter_I==0  	) continue;    
    
    h_NEvents_skimCut->Fill(1);
    
    
    if( fabs(vz_F)>24.              ) continue;
    h_NEvents_vzCut->Fill(1);
    
    // grab vzweight
    float vzWeight=1.;
    float vzStart=minbinValue_vzWeights, vzBinLeftSide=vzStart, vzBinRightSide=vzBinLeftSide+binsize_vzWeights;
    if(doVzWeights){
      for( int i=0; i<nbins_vzWeights ; i++ ) { 
	if(vzBinLeftSide<vz_F && vz_F<=vzBinRightSide) {
	  vzWeight=vzWeights[i];  
	  break; } 
	else {
	  vzBinLeftSide+=binsize_vzWeights;
	  vzBinRightSide+=binsize_vzWeights; } }
    }

    float evtPthatWeight=0.;    
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){ evtPthatWeight=pthatWeights[i]; } 
    
    float trigWeight=1.;
    //trigWeight = trigComb(trgDec, treePrescl, triggerPt);    
    
    float weight_eS=evtPthatWeight*trigWeight*vzWeight;              
    
    //vz
    if(fillMCEvtQAHists){

      hpthat->Fill(pthat_F, 1.);
      hWpthat->Fill(pthat_F, weight_eS);
      
      float vr_F=std::sqrt(vx_F*vx_F + vy_F*vy_F);
      hVr->Fill(vr_F, 1.0);   
      hWVr->Fill(vr_F, weight_eS);  
      //if(weight_eS!=0.)hTrgVr_noW->Fill(vr_F,1.0);

      hVx->Fill(vx_F, 1.0);   
      hWVx->Fill(vx_F, weight_eS);  
      //if(weight_eS!=0.)hTrgVx_noW->Fill(vx_F,1.0);

      hVy->Fill(vy_F, 1.0);   
      hWVy->Fill(vy_F, weight_eS);
      //if(weight_eS!=0.)hTrgVy_noW->Fill(vy_F,1.0);

      hVz->Fill(vz_F, 1.);
      hvzWVz->Fill(vz_F, vzWeight);
      hpthatWVz->Fill(vz_F, evtPthatWeight);
      hWVz->Fill(vz_F, weight_eS);
      
      hNref->Fill(nref_I, 1.);
      hWNref->Fill(nref_I, weight_eS);
    }
    

    if(fillMCJetMultHists){



      //// high Purity Fraction Event Cut
      bool has25PercHighPurTrks=false;
      if(nTrk_I<=10)has25PercHighPurTrks=true;
      else{
        int highPurCount=0;
        for(int itrk=0;itrk<nTrk_I;itrk++)
          if(trkHighPur_O[itrk]) highPurCount++;

        float highPurFrac=((float) highPurCount)/ ((float)nTrk_I);
	if(highPurFrac>0.25)
          has25PercHighPurTrks=true;
      }
      if(!has25PercHighPurTrks) continue;
      h_NEvents_trkCuts_1->Fill(0);

      //// cut on vtx_z spacings
      bool vtxsSpacedOutEnough=true;
      float vtx_z_1, vtx_z_2;
      for(int ivtx=0; ivtx<(nVtx_I-1);ivtx++){
        vtx_z_1=zVtx_F[ivtx];
        for(int jvtx=(ivtx+1); jvtx<nVtx_I;jvtx++){
          vtx_z_2=zVtx_F[jvtx];
          float deltaZ_vtx_1_2=fabs(vtx_z_1 - vtx_z_2);
          if(deltaZ_vtx_1_2 > 1.)
            continue;
          else{
            vtxsSpacedOutEnough=false; break;  }
	}//end 2nd vtx loop
        if(!vtxsSpacedOutEnough)break;
      }//end 1st vtx loop
      if(!vtxsSpacedOutEnough)continue;
      h_NEvents_trkCuts_2->Fill(0);

      //// cut on vtx positions w.r.t. beamspot/avgvtx positions
      // avg vtx positions consistent with beamspot plots made from payloads on conddb
      // tag: BeamSpotObjects_PCL_prompt_v0
      // gtag pointing to this beamspot tag: 75X_dataRun2_Prompt_ppAt5TeV_v0
      // no prescales, first vtx only,   0.0760011  / 0.0969312  / 0.285931
      // with prescales, first vtx only, 0.0759278  / 0.0969829  / 0.287621
      // with prescales, all vtxs,       0.07587616 / 0.09695976 / 0.2812252;
      float avgVtx_x=0.07587616 , avgVtx_y=0.09695976, avgVtx_z= 0.2812252;
      bool atLeastOnePVinRange=false;
      bool onePVinvzRange=false, onePVinvrRange=false;
      for(int ivtx=0; ivtx<nVtx_I;ivtx++){
        // check # of tracks, if good, check delta z
        if(nTrkVtx_I[ivtx]<2)
          continue;

        // check delta z, if good, check delta r
        float vtx_deltaZ=zVtx_F[ivtx]-avgVtx_z;
        if( fabs(vtx_deltaZ)<15. )
          onePVinvzRange=true;
        else {
          onePVinvzRange=false;
          continue;     }
	// check delta r, if good, stop checking for a vtx in range
        float vtx_deltaY=yVtx_F[ivtx]-avgVtx_y;
        float vtx_deltaX=xVtx_F[ivtx]-avgVtx_x;
        float vtx_deltaR=sqrt( vtx_deltaX*vtx_deltaX + vtx_deltaY*vtx_deltaY );
        if(vtx_deltaR<0.15)
          onePVinvrRange=true;
        else{
          onePVinvrRange=false;
          onePVinvzRange=false;
          continue;     }

        // if everything checks out and we make it here, we're good, we can exit the loop
        if( onePVinvzRange &&
            onePVinvrRange ) {
          atLeastOnePVinRange=true;
          break;        }

      }
      if(!atLeastOnePVinRange) continue;
      h_NEvents_trkCuts_3->Fill(0);


      // trk counting, # of tracks passing cuts
      int nTrk_passCuts_I=0;
      //int nTrk_vtxs_I=0, nTrk_vtxs_passCuts_I=0;
      for(int itrk=0;itrk<nTrk_I;itrk++){

        if (debugMode&&nEvt%10000==0){
	  std::cout<< " trkCharge =" << trkCharge_I[itrk]<< std::endl;
	  std::cout<< " trkNVtx   =" << (int)trkNVtx_I[itrk]<< std::endl;
	  std::cout<< " trkNdof   =" << (int)trkNdof_I[itrk]<< std::endl;
	  std::cout<< " trkNHit   =" << (int)trkNHit_I[itrk]<< std::endl;
	  std::cout<< " trkNLayer =" << (int)trkNLayer_I[itrk]<< std::endl;
	  std::cout<< " trkHighPur=" << (int)trkHighPur_O[itrk]<< std::endl;
	  std::cout<< " trkTight  =" << (int)trkTight_O[itrk]<< std::endl;
	  std::cout<< " trkLoose  =" << (int)trkLoose_O[itrk]<< std::endl;
	  std::cout<< " trkFake   =" << (int)trkFake_O[itrk]<< std::endl;
        }

        h_NTrks->Fill(0);

        if(!(trkHighPur_O[itrk]))continue;

	else if(trkFake_O[itrk])continue; //tracks not matched to PF particles are pinged as fake
	else if(matchedGenID_I[itrk][0] == 0)continue;//tracks matched to PF particles, but not subsequently matche to GEN particles have a matchedGenID of 0
        else if(trkPt_F[itrk]<0.4)continue;
        else if(fabs(trkEta_F[itrk])>2.4)continue;
        else if ((trkPtError_F[itrk]/trkPt_F[itrk])>0.1)continue;
        else if(fabs(trkDz_F[itrk]/trkDzError_F[itrk]) > 3.0) continue;
        else if(fabs(trkDxy_F[itrk]/trkDxyError_F[itrk]) > 3.0) continue;

        h_NTrks_kmatCut->Fill(0);

        nTrk_passCuts_I++;
	
      }
      
      
      
      ////// multiplicity cuts go here
      //bool inTrkMultClass=false;
      ////int minTrks=34, maxTrks=44; // peripheral, ~22-43% as of 12/12/17
      //int minTrks=80, maxTrks=300; // central, ~ 0-11% as of 12/12/17
      //if (nTrk_passCuts_I < minTrks || nTrk_passCuts_I > maxTrks)
      //        inTrkMultClass=false;
      //else
      //        inTrkMultClass=true;
      //
      //if(!inTrkMultClass)continue;



      // now loop over tracks for plots (kept separate for multiplicity class reasons)
      for(int itrk=0;itrk<nTrk_I;itrk++){
	
        if(!(trkHighPur_O[itrk]))continue;

	else if(trkFake_O[itrk])continue; //tracks not matched to PF particles are pinged as fake
	else if(matchedGenID_I[itrk][0] == 0)continue;//tracks matched to PF particles, but not subsequently matche to GEN particles have a matchedGenID of 0
        else if(trkPt_F[itrk]<0.4)continue;
        else if(fabs(trkEta_F[itrk])>2.4)continue;
        else if ((trkPtError_F[itrk]/trkPt_F[itrk])>0.1)continue;
        else if(fabs(trkDz_F[itrk]/trkDzError_F[itrk]) > 3.0) continue;
        else if(fabs(trkDxy_F[itrk]/trkDxyError_F[itrk]) > 3.0) continue;

        // TH1 QA
        trkPt->Fill( trkPt_F[itrk]  , weight_eS);
        trkPtError->Fill( trkPtError_F[itrk]  , weight_eS);
        trkPtPercError->Fill( (trkPtError_F[itrk]/trkPt_F[itrk])*(100.)  , weight_eS);

        trkEta->Fill(trkEta_F[itrk] , weight_eS);
        trkPhi->Fill(trkPhi_F[itrk] , weight_eS);

        trkDxy->Fill(trkDxy_F[itrk] , weight_eS);
        trkDxyError->Fill(trkDxyError_F[itrk] , weight_eS);
        trkDxySig->Fill( trkDxy_F[itrk]/trkDxyError_F[itrk] , weight_eS);

        trkDz->Fill(trkDz_F[itrk] , weight_eS);
        trkDzError->Fill(trkDzError_F[itrk] , weight_eS);
        trkDzSig->Fill( trkDz_F[itrk]/trkDzError_F[itrk] , weight_eS);

        trkCharge->Fill(trkCharge_I[itrk], weight_eS);

        trkNVtx->Fill((int)trkNVtx_I[itrk], weight_eS);
        trkChi2->Fill(trkChi2_F[itrk], weight_eS);
        trkNdof->Fill((int)trkNdof_I[itrk], weight_eS);

        trkNHit->Fill((int)trkNHit_I[itrk], weight_eS);
        trkNLayer->Fill((int)trkNLayer_I[itrk], weight_eS);

        trkHighPur->Fill(((int)trkHighPur_O[itrk]), weight_eS);
        trkTight->Fill(  ((int)trkTight_O[itrk]), weight_eS);
        trkLoose->Fill(  ((int)trkLoose_O[itrk]), weight_eS);
        trkFake->Fill(   ((int)trkFake_O[itrk]), weight_eS);

        // TH2 QA
        trkPtvEta->Fill( trkPt_F[itrk]  , trkEta_F[itrk],weight_eS);
        trkEtavPhi->Fill(trkEta_F[itrk], trkPhi_F[itrk],weight_eS);
      }//end trk loop


      // hists filled once-per-vtx
      //int nTrk_vtxs_I=0, nTrk_vtxs_passCuts_I;
      float vtx_minNormChi2_F=999.;
      int vtx_minNormChi2vtx_I=-1;
      for(int ivtx=0; ivtx<nVtx_I;ivtx++){

        if(nTrkVtx_I[ivtx]==0)continue;

        if(debugMode){
          if(ivtx!=0 && nTrkVtx_I[ivtx]==0)
	    std::cout<<"ivtx=="<<ivtx<<" nTrkVtx_I[ivtx]="<<nTrkVtx_I[ivtx]<<std::endl;
        }

        if(normChi2Vtx_F[ivtx]<vtx_minNormChi2_F) {
          vtx_minNormChi2_F=normChi2Vtx_F[ivtx];
          vtx_minNormChi2vtx_I=ivtx;
        }

        vtx_z->Fill(   zVtx_F[ivtx]          , weight_eS);
        vtx_y->Fill(   yVtx_F[ivtx]          , weight_eS);
        vtx_x->Fill(   xVtx_F[ivtx]          , weight_eS);

        float rVtx_F=sqrt(xVtx_F[ivtx]*xVtx_F[ivtx] + yVtx_F[ivtx]*yVtx_F[ivtx]);
        vtx_r->Fill(   rVtx_F                , weight_eS);
        vtx_rVz->Fill( rVtx_F,zVtx_F[ivtx]   , weight_eS);
        vtx_xVy->Fill( xVtx_F[ivtx],yVtx_F[ivtx]   , weight_eS);

        vtx_yOverx->Fill( yVtx_F[ivtx]/xVtx_F[ivtx], weight_eS);
	vtx_rOverz->Fill( rVtx_F/zVtx_F[ivtx], weight_eS);

        float radVtx_F=sqrt(rVtx_F*rVtx_F + zVtx_F[ivtx]*zVtx_F[ivtx]);
        vtx_zOverrad->Fill( zVtx_F[ivtx]/radVtx_F, weight_eS);

        //future; try computing theta two ways... or something?
        float theta=acos(zVtx_F[ivtx]/radVtx_F);
        vtx_theta->Fill( theta , weight_eS);

        float tantheta2=tan(theta/2.);
        vtx_tantheta2->Fill( tantheta2 , weight_eS);

        float lntantheta2=log(tantheta2);
        vtx_lntantheta2->Fill( lntantheta2 , weight_eS);

        float eta=-1.*lntantheta2;
        vtx_eta->Fill(eta,weight_eS);

        float phi=atan2(yVtx_F[ivtx],xVtx_F[ivtx]);
        vtx_phi->Fill(phi, weight_eS);

        vtx_rad->Fill(radVtx_F, weight_eS);



        Nvtxtrks->Fill(nTrkVtx_I[ivtx],weight_eS);
        NvtxtrksVz->Fill(nTrkVtx_I[ivtx], zVtx_F[ivtx], weight_eS);
        NvtxtrksVr->Fill(nTrkVtx_I[ivtx], rVtx_F      , weight_eS);

        normChi2_vtx->Fill(normChi2Vtx_F[ivtx], weight_eS);
        sumPtVtx_vtx->Fill(sumPtVtx_F[ivtx], weight_eS);
        //nTrk_vtxs_I+=nTrkVtx_I[ivtx];

	if(ivtx!=0){
          Nvtxtrks_vtxgt0->Fill(nTrkVtx_I[ivtx],weight_eS);
          NvtxtrksVz_vtxgt0->Fill(nTrkVtx_I[ivtx], zVtx_F[ivtx], weight_eS);
          NvtxtrksVr_vtxgt0->Fill(nTrkVtx_I[ivtx], rVtx_F      , weight_eS);

          normChi2_vtxgt0->Fill(normChi2Vtx_F[ivtx], weight_eS);
          sumPtVtx_vtxgt0->Fill(sumPtVtx_F[ivtx], weight_eS);

        }
        else{// ivtx==0, vtx hists filled once-per-event
          //Ntrks_vtx0->Fill(nTrkVtx_I[0], weight_eS);
          Nvtxtrks_vtx0->Fill(nTrkVtx_I[0], weight_eS);
          NvtxtrksVz_vtx0->Fill(nTrkVtx_I[0], zVtx_F[0], weight_eS);
          NvtxtrksVr_vtx0->Fill(nTrkVtx_I[0], rVtx_F   , weight_eS);

          normChi2_vtx0->Fill(normChi2Vtx_F[0], weight_eS);
          sumPtVtx_vtx0->Fill(sumPtVtx_F[0], weight_eS);
          //Nvtxtrks_vtx0_passCuts->Fill(nTrk_passCuts_vtx0_I, weight_eS);
          //Ntrks_passCuts_vtx0->Fill(nTrk_passCuts_vtx0_I, weight_eS);

        }

      }//end vtx loop

      // hists filled once-per-event
      Nvtxs->Fill(nVtx_I,weight_eS);

      Ntrks->Fill(nTrk_I,weight_eS);
      Ntrks_passCuts->Fill(nTrk_passCuts_I,weight_eS);// pass cuts, not necessarily assoc w/ vertices

      //Ntrks_vtxs->Fill(nTrk_vtxs_I, weight_eS);       //assoc w/ vertecies, not passing cuts
      //Ntrks_vtxs_passCuts->Fill(nTrk_vtxs_passCuts_I, weight_eS);       //assoc w/ vertecies, not passing cuts

      vtx_maxPt->Fill(maxPtVtx_I,weight_eS);
      vtx_maxMult->Fill(maxMultVtx_I,weight_eS);
      vtx_minNormChi2->Fill(vtx_minNormChi2vtx_I,weight_eS);

      vtx_maxPt_v_maxMult->Fill(maxPtVtx_I,maxMultVtx_I,weight_eS);

      Nvtxs_v_Nref->Fill(  nVtx_I   , nref_I , weight_eS);
      Ntrks_v_Nref->Fill(  nTrk_passCuts_I , nref_I , weight_eS);
      Ntrks_v_Nvtxs->Fill( nTrk_passCuts_I , nVtx_I , weight_eS);




    }//end jetMult portion
    
    
    

    
    
    
    
    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;
    
    //bool firstGoodJetFound_wJetID=false, secondGoodJetFound_wJetID=false; 
    //float firstGoodJetPt_wJetID=-1., secondGoodJetPt_wJetID=-1.;
    //float firstGoodJetPhi_wJetID=-1., secondGoodJetPhi_wJetID=-1.;
    
    // for event counting + avoiding duplicate fills in dijet hists
    bool dijetHistsFilled=false;//, dijetHistsFilled_wJetID=false;
    bool hNEvts_withJets_Filled=false, hNEvts_withJets_kmatCut_Filled=false, hNEvts_withJets_JetIDCut_Filled=false; 
    
    int jetsPerEvent=nref_I;//, jetsPerEventJetID=nref_I;
    float evt_leadJetPt=-1.;//, evt_leadJetPt_wJetID=-1.;
    
    for(int jet = 0; jet<nref_I; ++jet)      {

      // event+jet counting
      h_NJets->Fill(1);
      h_WNJets->Fill(1,weight_eS);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      
      
      float recpt  = pt_F[jet];
      float receta = eta_F[jet];
      float absreceta = fabs(eta_F[jet]);
      float genpt  = refpt_F[jet];
      float geneta = refeta_F[jet];

      //if reco jet w/o matched gen jet, skip.
      if( subid_F[jet] != 0 )             { jetsPerEvent--; /*jetsPerEventJetID--;*/ continue; }
      else if( genpt <= genJetPtCut   )   { jetsPerEvent--; /*jetsPerEventJetID--;*/ continue; }
      else if( genpt > genJetPtCut_Hi   ) { jetsPerEvent--; /*jetsPerEventJetID--;*/ continue; }
      else if( recpt <= jtPtCut   )       { jetsPerEvent--; /*jetsPerEventJetID--;*/ continue; }           
      else if( recpt > jtPtCut_Hi   )     { jetsPerEvent--; /*jetsPerEventJetID--;*/ continue; }           
      
      
      if( absreceta >= jtEtaCutHi ) { 
	//jetsPerEvent--; jetsPerEventJetID--;
	continue; }
      else if( absreceta < jtEtaCutLo )  { 
	//jetsPerEvent--; jetsPerEventJetID--; 
	continue; }
      

      if(absreceta>=4.7) jetsPerEvent--;

      // largest jet pt in each event
      if(recpt>evt_leadJetPt && absreceta<4.7)
	evt_leadJetPt=recpt;
      //if(recpt>evt_leadJetPt_wJetID && fillMCJetIDHists && passesJetID)
      //evt_leadJetPt_wJetID=recpt;
      
      // jet/event counts
      h_NJets_kmatCut->Fill(1);
      h_WNJets_kmatCut->Fill(1,weight_eS);
      if(!hNEvts_withJets_kmatCut_Filled){
	h_NEvents_withJets_kmatCut->Fill(1);
	hNEvts_withJets_kmatCut_Filled=true;  }           
      
      //float rawpt  = rawpt_F[jet];      
      float recphi = phi_F[jet];
      float jetIDpt=recpt;
      
      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;
      
      // 13 TeV JetID criterion, loose or tight
      bool passesJetID=false;
      if(fillMCJetIDHists) 	{
	if (!(absreceta > 2.4)) 
	  passesJetID=(bool)jetID_00eta24( jetIDpt, 
					   neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet],
					   numConst,  chMult);
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  passesJetID=(bool) jetID_24eta27( jetIDpt,
					    neSum_F[jet],  phSum_F[jet], 
					    numConst);
	else if( !(absreceta>3.0) && absreceta>2.7 )
	  passesJetID=(bool) jetID_27eta30( jetIDpt,
					    neSum_F[jet],  phSum_F[jet], 
					    numConst);
	else  
	  passesJetID=(bool)jetID_32eta47( jetIDpt, 
					   phSum_F[jet]);
	if(fillMCJetIDHists && !passesJetID)continue;
      }
      
      
      
      /////   JETQA   ///// 
      if(fillMCJetQAHists	 ){
	
	
        int ind=0;
        //jets
        hJetQA[jetIDint][ind]->Fill(recpt, weight_eS); ind++;
        hJetQA[jetIDint][ind]->Fill(rawpt_F[jet], weight_eS); ind++;
	if(recpt>jetQAPtCut){//second jet pt cut
          hJetQA[jetIDint][ind]->Fill(eta_F[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(phi_F[jet], weight_eS); ind++;
	  
          //tracks
          hJetQA[jetIDint][ind]->Fill(trkN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(trkSum_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(trkMax_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(trkHardN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(trkHardSum_F[jet]/jetIDpt, weight_eS); ind++;

          //PF photons
          hJetQA[jetIDint][ind]->Fill(phN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(phSum_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(phMax_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(phHardN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(phHardSum_F[jet]/jetIDpt, weight_eS); ind++;

          //PF charged hadrons
          hJetQA[jetIDint][ind]->Fill(chN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(chSum_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(chMax_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(chHardN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(chHardSum_F[jet]/jetIDpt, weight_eS); ind++;
	  
          //PF neutral hadons
          hJetQA[jetIDint][ind]->Fill(neN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(neSum_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(neMax_F[jet]/jetIDpt, weight_eS); ind++;
	  
          //PF electrons
          hJetQA[jetIDint][ind]->Fill(eN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(eSum_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(eMax_F[jet]/jetIDpt, weight_eS); ind++;
	  
          //PF muons
          hJetQA[jetIDint][ind]->Fill(muN_I[jet], weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(muSum_F[jet]/jetIDpt, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(muMax_F[jet]/jetIDpt, weight_eS); ind++;

          //PF particle sums
          hJetQA[jetIDint][ind]->Fill(neuMult, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(chMult, weight_eS); ind++;
          hJetQA[jetIDint][ind]->Fill(numConst, weight_eS); ind++;

	  if(fillgenJetQA){
	    hMCJetQA_genpt[jetIDint]    ->Fill(genpt,weight_eS);
	    hMCJetQA_geneta[jetIDint]   ->Fill(geneta,weight_eS);
	    hMCJetQA_genrecpt[jetIDint] ->Fill(recpt/genpt,weight_eS);
	    hMCJetQA_genreceta[jetIDint]->Fill(receta-geneta,weight_eS);	}

	}
	
	//looking for the first two good jets that meet the criteria specified
	if(fillMCDijetHists){

	  if ( !firstGoodJetFound ){
	    if(recpt>ldJetPtCut) { firstGoodJetFound=true;
	      firstGoodJetPt =recpt; 
	      firstGoodJetPhi=recphi; } }
	  
	  else if ( !secondGoodJetFound && 
		    firstGoodJetFound      ) {
	    if (recpt>subldJetPtCut) { 
	      float checkdPhi=deltaphi(firstGoodJetPhi,recphi);
	      float ptAve=(firstGoodJetPt+recpt)/2.;
	      if( checkdPhi>dPhiCut && ptAve>ptAveCut) {
		if(recpt>firstGoodJetPt){secondGoodJetFound=true;
		  std::cout<<std::endl<<
		    "WARNING: picked wrong jet for lead jet! Swapping..."<<std::endl<<std::endl;
		  secondGoodJetPt  = firstGoodJetPt;          
		  secondGoodJetPhi = firstGoodJetPhi;
		  firstGoodJetPt   = recpt;          
		  firstGoodJetPhi  = recphi; }
		else { secondGoodJetFound=true;
		  secondGoodJetPt  = recpt;
		  secondGoodJetPhi = recphi; }}} }
	  
	  if (!dijetHistsFilled &&
	      firstGoodJetFound && 
	      secondGoodJetFound ) { dijetHistsFilled=true;
	    float A_j=(firstGoodJetPt-secondGoodJetPt)
	      /(firstGoodJetPt+secondGoodJetPt);
	    float x_j=(secondGoodJetPt/firstGoodJetPt); 	
	    float dphi=deltaphi(firstGoodJetPhi,secondGoodJetPhi);
	    int ind=dijtInd;
            hJetQA[jetIDint][ind]->Fill( A_j , weight_eS ); ind++;
            hJetQA[jetIDint][ind]->Fill( x_j , weight_eS ); ind++;
            hJetQA[jetIDint][ind]->Fill( dphi , weight_eS ); ind++;
            hJetQA[jetIDint][ind]->Fill( firstGoodJetPt , weight_eS ); ind++;
            hJetQA[jetIDint][ind]->Fill( secondGoodJetPt , weight_eS ); ind++;
	  } }
      }//end fillMCJetQAHists	      
      
      
      if(fillMCJetIDHists){	
	// jet/event counts
	h_NJets_JetIDCut->Fill(1);
	h_WNJets_JetIDCut->Fill(1,weight_eS);
	if(!hNEvts_withJets_JetIDCut_Filled){
	  h_NEvents_withJets_JetIDCut->Fill(1);
	  hNEvts_withJets_JetIDCut_Filled=true;      }	
      }//end fillMCJetIDHists
      
    }//end jet loop
    
    if(fillMCEvtQAHists){
      //if(fillMCJetIDHists){
	//hLeadJetPt_wJetID->Fill(evt_leadJetPt_wJetID ,weight_eS);
	//hjetsPEvt_wJetID ->Fill(jetsPerEventJetID,1.0);   
	//hWjetsPEvt_wJetID->Fill(jetsPerEventJetID,weight_eS); 
	//}    
      if(evt_leadJetPt>0.)hLeadJetPt->Fill(evt_leadJetPt ,weight_eS);
      hjetsPEvt ->Fill(jetsPerEvent,1.0);   
      hWjetsPEvt->Fill(jetsPerEvent,weight_eS); }
    
    //assert(false);
  }//end event loop
  
  
  
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Event-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl<<std::endl;
  
  std::cout<<"Total Num of Events passing skimCuts  = " << h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events passing vzCuts    = " << h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
  
  std::cout<<"Total Num of Events passing 25% highpur track cut = " <<h_NEvents_trkCuts_1->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events passing spaced out vtxs cut = " <<h_NEvents_trkCuts_2->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events passing avg vtx position cut = " <<h_NEvents_trkCuts_3->GetEntries()<<std::endl;
  
  std::cout<<std::endl<<"/// Job Track-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Tracks read from good events                          = " <<    h_NTrks->GetEntries()<<std::endl;
  std::cout<<"Total Num of Tracks read from good events kmatCuts              = " <<    h_NTrks_kmatCut->GetEntries()<<std::endl;
  std::cout<<std::endl;


  std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
  
  std::cout<<"Total Num of good Events w/ jets = " << h_NEvents_withJets->GetEntries()<<std::endl; //note, approximate, some events make it through with nref=0 for some reason
  std::cout<<"Total Num of good Events, w/ jets, post kmatCuts = " << h_NEvents_withJets_kmatCut->GetEntries() << std::endl;   
  if(fillMCJetIDHists)
    std::cout << "Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " << h_NEvents_withJets_JetIDCut->GetEntries() << std::endl ;
  
  std::cout<<"Total Num of Jets read from good events                          = " <<    h_NJets->GetEntries()<<std::endl;  
  std::cout<<"Total Num of Jets read from good events post kinCut              = " <<    h_NJets_kmatCut->GetEntries()<<std::endl;
  if(fillMCJetIDHists)
    std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " << h_NJets_JetIDCut->GetEntries()<<std::endl;
  std::cout<<std::endl;
  
  


  std::cout<<"writing output file... ";
  fout->Write(); 
  
  std::cout<<std::endl<<"readForests_ppMC_jetMult finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real tyime (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  std::cout<<"output file written is  "<<outfile<<std::endl<<std::endl;
  
  return 0;
}// end readForests_ppMC_jetMult
								 

////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1 && argc!=10){
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppMC_jetMult.exe ";
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
  if(argc==1) rStatus = readForests_ppMC_jetMult();
  else if (argc==8){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppMC_jetMult( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
				  outputFileName);
  }
  else if (argc==10){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
    rStatus = readForests_ppMC_jetMult( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
					   outputFileName, jtEtaCutLo, jtEtaCutHi);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
