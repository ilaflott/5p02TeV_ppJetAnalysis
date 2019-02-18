// Custom header
#include "readForests_jetPlots.h"

// ppMC switches
const bool fillMCEvtQAHists=true;
const bool fillMCJetQAHists=true;
const bool fillgenJetQA=true&fillMCJetQAHists;
const bool fillMCJetIDHists=true;//, tightJetID=false;
const bool fillMCJetSpectraRapHists=false; //other
const bool fillgenJetRapHists=false&&fillMCJetSpectraRapHists;  //other switches

const bool fillMCDijetHists=false;//hardly used options

const int jetIDint=(int)fillMCJetIDHists;

//// readForests_ppMC_jetPlots
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_jetPlots(std::string inFilelist , int startfile , int endfile ,  
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
  TH1D* hLeadJetPt=NULL;

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
  TH1D *h_NJets_JetIDCut = NULL;
  if(fillMCJetIDHists)
    h_NJets_JetIDCut = new TH1D("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 2,0.,2.);
  


  /////   JET QA HISTS   ///// 
  //jet QA
  TH1D *hJetQA[2][N_vars]={};

  //MC JetQA
  TH1D *hMCJetQA_genpt[2]={};
  TH1D *hMCJetQA_geneta[2]={};
  TH1D *hMCJetQA_genphi[2]={};
  TH1D *hMCJetQA_genrecpt[2]={};
  TH1D *hMCJetQA_genreceta[2]={};
  TH1D *hMCJetQA_genrecphi[2]={};
  TH1D *hMCJetQA_gendrjt[2]={};
  
  if(fillMCJetQAHists){
    for(int k = 0; k<2; ++k){

      if(!fillMCJetIDHists && k==1)continue;
      if(fillMCJetIDHists && k==0)continue;
      //if( fillMCJetIDHists && k==0)continue;	
      for(int j = 0; j<N_vars; ++j){	
	//binnings+hist declarations, etc.

        //jets
        if(var[j]=="jtpt"||var[j]=="rawpt")
          hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;", var[j].c_str()) , 2500,0,2500);
        else if(var[j]=="jteta")
          hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-5,+5);
        else if(var[j]=="jtphi")
          hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 100,-4,+4);

        //jetconst. counts
        else if(var[j]=="trkN"|| var[j]=="phN"|| var[j]=="chN"|| var[j]=="neN"|| var[j]=="eN"|| var[j]=="muN")
          hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 60,0,60);
        else if(var[j]=="trkHardN"|| var[j]=="phHardN"|| var[j]=="chHardN")
          hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 60,0,60);
        else if(var[j]=="neuMult"|| var[j]=="chMult"|| var[j]=="numConst")
          hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 100,0,100);

	//dijets
	else if (fillMCDijetHists){
	  if(var[j]=="dphi") 
	    hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 50,0,+4);
	  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt")
	    hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;", var[j].c_str()), 2500,0,2500);
	  else //xj and Aj binnings
	    hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);         
	}
	
	//consituent binnings
	else 
	  hJetQA[k][j] = new TH1D( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);         
      }
      
      if(fillgenJetQA){	
	hMCJetQA_genpt[k]     = new TH1D(Form("hJetQA_%dwJetID_genpt",k)     , "genpt",    2500, 0,2500 );
	hMCJetQA_geneta[k]    = new TH1D(Form("hJetQA_%dwJetID_geneta",k)    , "geneta", 200,  -6.,6. );
	hMCJetQA_genphi[k]    = new TH1D(Form("hJetQA_%dwJetID_genphi",k)    , "genphi", 1280,  -6.4,6.4 );
	hMCJetQA_genrecpt[k]  = new TH1D(Form("hJetQA_%dwJetID_genrecpt",k)  , "genrecpt",    200, 0., 2. );
	hMCJetQA_genreceta[k] = new TH1D(Form("hJetQA_%dwJetID_genreceta",k) , "genreceta", 200,  -10, 10 ); 
	hMCJetQA_genrecphi[k] = new TH1D(Form("hJetQA_%dwJetID_genrecphi",k) , "genrecphi", 256,  -12.8, 12.8 ); 
	hMCJetQA_gendrjt[k] = new TH1D(Form("hJetQA_%dwJetID_gendrjt",k) , "gendrjt", 200,  -1., 1. ); 
	
	
      }
      
    }}
  
  TH1D *hJetSpectraRap[2][nbins_abseta]={};
  TH1D *hMCJetQA_rapBins_genpt[2][nbins_abseta]={};
  TH1D *hMCJetQA_rapBins_geneta[2][nbins_abseta]={};
  TH1D *hMCJetQA_rapBins_genrecpt[2][nbins_abseta]={};
  TH1D *hMCJetQA_rapBins_genreceta[2][nbins_abseta]={};
  if(fillMCJetSpectraRapHists){
    for(int k = 0; k<2; ++k){
      if(!fillMCJetIDHists && k==1)continue;	
      if(fillMCJetIDHists && k==0)continue;

      for(int j = 0; j<nbins_abseta; ++j){
	std::string h_Title="hJetSpectraRap_";
	if(k==1)h_Title+="wJetID_";
	h_Title+="bin"+std::to_string(j);
	std::stringstream stream1, stream2;	
	stream1.precision(1); stream1 << std::fixed << absetabins[j];
	stream2.precision(1); stream2 << std::fixed << absetabins[j+1];
	std::string h_Desc="JetPt Spectra for "+stream1.str()+"<abs(y)<"+ stream2.str();	
	hJetSpectraRap[k][j]=new TH1D(h_Title.c_str(),h_Desc.c_str(), 2500,0,2500);  
	
	if(fillgenJetRapHists){
	  hMCJetQA_rapBins_genpt[k][j]     = new TH1D ( Form("hJetSpectraRap_%dwJetID_bin%d_genpt",k,j), "genpt", 2500,0,2500);
	  hMCJetQA_rapBins_geneta[k][j]    = new TH1D ( Form("hJetSpectraRap_%dwJetID_bin%d_geneta",k,j), "geneta", 200,-6.,6.);
	  hMCJetQA_rapBins_genrecpt[k][j]  = new TH1D ( Form("hJetSpectraRap_%dwJetID_bin%d_genrecpt",k,j),  "rec/gen,pt", 100, 0.,2.);
	  hMCJetQA_rapBins_genreceta[k][j] = new TH1D ( Form("hJetSpectraRap_%dwJetID_bin%d_genreceta",k,j), "rec-gen,eta", 2000, -100,100);
	}
	
      }}}
  
  
  
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
  
  //float jetIDCut_neSum, jetIDCut_phSum;
  //if(tightJetID){     jetIDCut_neSum=0.90;  jetIDCut_phSum=0.90;}
  //else{     jetIDCut_neSum=0.99;  jetIDCut_phSum=0.99;}
  
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
    
    int jetsPerEvent=nref_I;//, jetsPerEventJetID=nref_I;
    float evt_leadJetPt=-1.;//, evt_leadJetPt_wJetID=-1.;
    
    for(int jet = 0; jet<nref_I; ++jet)      {
      
      // event+jet counting
      h_NJets->Fill(0.);      h_NJets->Fill(1.,weight_eS);
      
      float recpt     = pt_F[jet];
      float receta    = eta_F[jet];
      float recphi    = phi_F[jet];
      float absreceta = fabs(eta_F[jet]);
      //float rawpt  = rawpt_F[jet];      

      float genpt   = refpt_F[jet];
      float genphi  = refphi_F[jet];
      float geneta  = refeta_F[jet];
      float gendrjt = refdrjt_F[jet];
      
      //if reco jet w/o matched gen jet, skip.
      if( subid_F[jet] != 0 ) { jetsPerEvent--;   continue; }
      h_NJets_matched->Fill(0.);      h_NJets_matched->Fill(1.,weight_eS);

      if( gendrjt > gendrCut ){	jetsPerEvent--;	continue;  }
      h_NJets_deltaRCut->Fill(0.);      h_NJets_deltaRCut->Fill(1.,weight_eS);
      
      if( !(recpt > jtPtCut)   )      { jetsPerEvent--;     continue; } //paranoia about float comparison
      h_NJets_jtptCut->Fill(0.);      h_NJets_jtptCut->Fill(1.,weight_eS);

      if( !(recpt < jtPtCut_Hi)   )      { jetsPerEvent--;  continue; }           
      h_NJets_jtptCut_Hi->Fill(0.);      h_NJets_jtptCut_Hi->Fill(1.,weight_eS);
      
      if(!(absreceta<4.7)) {jetsPerEvent--;continue;}
      h_NJets_jtetaCut1->Fill(0.);      h_NJets_jtetaCut1->Fill(1.,weight_eS);
      
      if( absreceta  < jtEtaCutLo ) { jetsPerEvent--;	continue;}
      if( !(absreceta < jtEtaCutHi) ) { jetsPerEvent--;	continue;}
      h_NJets_jtetaCut2->Fill(0.);      h_NJets_jtetaCut2->Fill(1.,weight_eS);
      
      if( !(genpt > genJetPtCut)   )  { jetsPerEvent--;      continue; } //paranoia about float comparison
      h_NJets_genjtptCut->Fill(0.);      h_NJets_genjtptCut->Fill(1.,weight_eS);
      if( !(genpt < genJetPtCut_Hi)   )  { jetsPerEvent--;   continue; }
      h_NJets_genjtptCut_Hi->Fill(0.);      h_NJets_genjtptCut_Hi->Fill(1.,weight_eS);
      
      // largest jet pt in each event
      if(recpt>evt_leadJetPt && absreceta<4.7)	evt_leadJetPt=recpt;
      
      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;
      
      float jetIDpt=recpt;//ala HIN jetID, recpt is corrected w/ L2/L3 residuals
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
					    neuMult);
	else  
	  passesJetID=(bool)jetID_32eta47( jetIDpt, 
					   phSum_F[jet]);
	if(!passesJetID)continue;
	h_NJets_JetIDCut->Fill(0.);		h_NJets_JetIDCut->Fill(1.,weight_eS);	
      }
      
      
      //fill jetspectraRapHists w/ passing jetID criterion
      if(  fillMCJetSpectraRapHists ) { 
	int theRapBin=-1;
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin){
	  if( absetabins[rapbin]<=absreceta  && 		
	      absreceta<absetabins[rapbin+1]    ) {	    
	    theRapBin=rapbin;	    
	    hJetSpectraRap[jetIDint][theRapBin]->Fill(recpt,weight_eS);  
	    if(fillgenJetRapHists){
	      hMCJetQA_rapBins_genpt[jetIDint][theRapBin]     ->Fill(genpt,weight_eS);
	      hMCJetQA_rapBins_geneta[jetIDint][theRapBin]    ->Fill(geneta,weight_eS);
	      hMCJetQA_rapBins_genrecpt[jetIDint][theRapBin]  ->Fill(recpt/genpt,weight_eS);
	      hMCJetQA_rapBins_genreceta[jetIDint][theRapBin] ->Fill(receta-geneta,weight_eS);	 	    }	    
	    break;  
	  }
	}
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
	    hMCJetQA_genphi[jetIDint]   ->Fill(genphi,weight_eS);
	    hMCJetQA_genrecpt[jetIDint] ->Fill(recpt/genpt,weight_eS);

	    hMCJetQA_genreceta[jetIDint]->Fill(fabs(receta-geneta), weight_eS);	
	    hMCJetQA_genrecphi[jetIDint]->Fill(fabs(recphi-genphi), weight_eS);	
	    
	    hMCJetQA_gendrjt[jetIDint]->Fill(gendrjt, weight_eS);
	  }
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
      

    }//end jet loop
    
    if(fillMCEvtQAHists){
      hjetsPEvt ->Fill(jetsPerEvent,1.0);   
      hWjetsPEvt->Fill(jetsPerEvent,weight_eS); 
      if(evt_leadJetPt>0.)hLeadJetPt->Fill(evt_leadJetPt ,weight_eS);
    }    
  }//end event loop

  

  
  std::cout<<std::endl;
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"/// Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"# of Events (criteria/description) = count / pthat*vz weighted count" <<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"# evts passing w cut applied = count / xsec weighted count" <<std::endl;
  std::cout<<"read                         = " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;
  std::cout<<"skimCuts  = read - !skim     = " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;//same
  std::cout<<"vzCuts    = skimCuts - !vz   = " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl; //same
  std::cout << std::endl;


  std::cout<<"//////////////////////////////"<<std::endl;
  std::cout<<"/// Jet-Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////////"<<std::endl;

  std::cout<<std::endl;
  std::cout<<"# of Jets read                     = " << h_NJets->GetBinContent(1) << " / " << h_NJets->GetBinContent(2) <<std::endl;
  std::cout<<"# of GEN-RECO Matched Jets         = " << h_NJets_matched->GetBinContent(1) << " / " << h_NJets_matched->GetBinContent(2) <<std::endl;
  std::cout<<"# of Matched Jets meeting dR cut   = " << h_NJets_deltaRCut->GetBinContent(1) << " / " << h_NJets_deltaRCut->GetBinContent(2) <<std::endl;  
  std::cout<<"# of Jets w/ jtpt > jtptCut        = " << h_NJets_jtptCut->GetBinContent(1) << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ jtpt > jtptCut_Hi     = " << h_NJets_jtptCut_Hi->GetBinContent(1) << " / " << h_NJets_jtptCut_Hi->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ |jteta|<4.7           = " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ passing jtEtaCut(s)   = " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ genjtpt>genjtptCut    = " << h_NJets_genjtptCut->GetBinContent(1) << " / " << h_NJets_genjtptCut->GetBinContent(2) <<std::endl;
  std::cout<<"# of Jets w/ genjtpt<genjtptCut_Hi = " << h_NJets_genjtptCut_Hi->GetBinContent(1) << " / " << h_NJets_genjtptCut_Hi->GetBinContent(2) <<std::endl;
  if(fillMCJetIDHists)
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
