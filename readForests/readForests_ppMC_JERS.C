// custom header
#include "readForests_JERS.h"

// ppMC switches
const bool fillMCEvtQAHists=true;
const bool fillJERSHists=true;
const bool fillMCUnfoldingHists=false;
const bool fillMCEffHists=true;
const bool fillMCJetIDHists=true;//, tightJetID=false;

//// readForests_ppMC_JERS
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_JERS(std::string inFilelist , int startfile , int endfile ,  
			  int radius , std::string jetType , bool debugMode ,
			  std::string outfile  ,
			  float jtEtaCutLo, float jtEtaCutHi ){ 
  
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
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"      ,(std::to_string(jtPtCut)).c_str()   ,    100, 0,100); hJetPtCut->Fill(jtPtCut);           
  TH1F *hGenJetPtCut     =new TH1F("hGenJetPtCut"  ,(std::to_string(genJetPtCut)).c_str()   ,    100, 0,100); hGenJetPtCut->Fill(genJetPtCut);           
  TH1F *hJetEtaCutHi       =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6  ); hJetEtaCutHi->Fill(jtEtaCutHi);
  TH1F *hJetEtaCutLo       =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6  ); hJetEtaCutLo->Fill(jtEtaCutLo);	


  //TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  //TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  //TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);  
  //TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  

  TH1F *h_NEvents         = new TH1F("NEvents","NEvents", 1,0,2);
  TH1F *h_NEvents_read    = new TH1F("NEvents_read","NEvents read", 1,0,2);
  TH1F *h_NEvents_skimCut = new TH1F("NEvents_skimCut"      ,"NEvents read post skimCut", 1,0,2);
  TH1F *h_NEvents_vzCut   = new TH1F("NEvents_vzCut"        ,"NEvents read post vzCut AND skimCut", 1,0,2);
  
  //counts post all evt cuts, that are thrown out because of no good jets pre and post jet cuts
  TH1F *h_NEvents_withJets           = new TH1F("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut1   = new TH1F("NEvents_withJets_kmatCut1" , 
						"NEvents read post evt cuts, w/ jets post kmatCut1", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut2   = new TH1F("NEvents_withJets_kmatCut2" , 
						"NEvents read post evt cuts, w/ jets post kmatCut2", 1,0,2);
  
  //jet counts
  TH1F *h_NJets          = new TH1F("NJets","NJets read", 1,0,2);
  TH1F *h_NJets_kmatCut1  = new TH1F("NJets_kmatCut1 ","NJets read post kmatCut1 ", 1,0,2);
  TH1F *h_NJets_kmatCut2  = new TH1F("NJets_kmatCut2 ","NJets read post kmatCut2 ", 1,0,2);
  
  
  /////   EVENT QA   ///// 
  TH1F *hVz=NULL, *hpthatWVz=NULL, *hWVz=NULL, *hvzWVz=NULL; 
  TH1F *hpthat=NULL, *hWpthat=NULL;  
  if(fillMCEvtQAHists){
    hVz       = new TH1F("hVz","", 60,-15.,15.);//evtvz
    hpthatWVz = new TH1F("hpthatWeightedVz","", 60,-15.,15.);//pthat-weighted evtvz
    hvzWVz    = new TH1F("hvzWeightedVz","", 60,-15.,15.);//vz-weighted evtvz
    hWVz      = new TH1F("hWeightedVz","", 60,-15.,15.);//pthat*vz-weighted evt vz
    hpthat    = new TH1F("hpthat","",1000,0,1000);//evt pthat, unweighted and weighted
    hWpthat   = new TH1F("hWeightedpthat","",1000,0,1000);  }
  
  
  /////  UNFOLDING   /////
  //to unfold ppData 
  TH1F *hpp_gen[2]={};    
  TH1F *hpp_reco[2]={}; 
  TH2F *hpp_matrix[2]={}; 
  //to test MC sample consistency in unfolding
  TH1F *hpp_mcclosure_gen[2]={};      //the first three are for the "truth" response matrix
  TH1F *hpp_mcclosure_reco[2]={};         
  TH2F *hpp_mcclosure_matrix[2]={};       
  TH1F *hpp_mcclosure_reco_test[2]={};   //set of MC "data" to test "truth" response matrix on 
  
  if(fillMCUnfoldingHists){

    //int jtID=0;
    //if(fillMCJetIDHists)jtID=1;
    
    std::string hUnfTitleArray[]={ "gen", "reco", "matrix", 
				   "mcclosure_gen" , "mcclosure_reco" , "mcclosure_matrix" , "mcclosure_reco_test" };
    const int nUnfTitles=sizeof(hUnfTitleArray)/sizeof(std::string);    
    
    for(int jtID=0;jtID<2;jtID++){
      if(!fillMCJetIDHists&&jtID==1)continue;
      
      for(int k=0; k<nUnfTitles; k++){     

	std::string hTitle="hpp_"+hUnfTitleArray[k];
	if(jtID==1)hTitle+="_wJetID";
	hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;      
	
	if(hUnfTitleArray[k]=="gen")	
	  hpp_gen[jtID]    = new TH1F( hTitle.c_str(), "refpt", 1000,0,1000);
	else if(hUnfTitleArray[k]=="reco")	
	  hpp_reco[jtID]    = new TH1F( hTitle.c_str(), "recopt", 1000,0,1000);
	else if(hUnfTitleArray[k]=="matrix")	
	  hpp_matrix[jtID]    = new TH2F( hTitle.c_str(), "reftpt v. recopt", 1000, 0,1000, 1000, 0,1000);
	else if(hUnfTitleArray[k]=="mcclosure_gen")	
	  hpp_mcclosure_gen[jtID]    = new TH1F( hTitle.c_str(), "refpt for mcclosure truth", 1000,0,1000);
	else if(hUnfTitleArray[k]=="mcclosure_reco")	
	  hpp_mcclosure_reco[jtID]    = new TH1F( hTitle.c_str(), "recopt for mcclosure truth", 1000,0,1000);
	else if(hUnfTitleArray[k]=="mcclosure_matrix")	
	  hpp_mcclosure_matrix[jtID]    = new TH2F( hTitle.c_str(), "reftpt v. recopt for mcclosure truth", 1000, 0,1000, 1000, 0,1000);
	else if(hUnfTitleArray[k]=="mcclosure_reco_test")	
	  hpp_mcclosure_reco_test[jtID]    = new TH1F( hTitle.c_str(), "recopt for mcclosure test", 1000,0,1000);
      }
    }  
  }
  
  
  /////   GEN/RECO MATCHING   ///// 
  TH2F *hpp_mceff_pt[2]={}, *hpp_mceff_eta[2]={}, *hpp_mceff_phi[2]={};//, *hpp_mceff_drjt[2]={};
  TH2F *hpp_mceff_ptrat_drjt[2]={}, *hpp_mceff_ptrat_eta[2]={}, *hpp_mceff_ptrat_phi[2]={};// all v. genpt/recpt  

  if(fillMCEffHists){

    //int jtID=0;
    //if(fillMCJetIDHists)jtID=1;    

    std::string hMCEffTitleArray[]={ "pt", "eta", "phi", "ptrat_drjt", "ptrat_eta", "ptrat_phi" };
    const int nMCEffTitles=sizeof(hMCEffTitleArray)/sizeof(std::string);    
   
    for(int jtID=0;jtID<2;jtID++){

      if(!fillMCJetIDHists&&jtID==1)continue;

      for(int k=0; k<nMCEffTitles; k++){      	
	std::string hTitle="hpp_mceff_"+hMCEffTitleArray[k];
	if(jtID==1)hTitle+="_wJetID";
	hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;
	if(debugMode)std::cout<<"k="<<k<<std::endl;
	if(debugMode)std::cout<<"hMCEffTitleArray="<<hMCEffTitleArray[k]<<std::endl;
	
	if(hMCEffTitleArray[k]=="pt")	
	  hpp_mceff_pt[jtID] =  new TH2F( hTitle.c_str(), "refpt v recopt/refpt"          , 100, 0.,5., 1000, 0.,1000. );
	else if(hMCEffTitleArray[k]=="eta")						          
	  hpp_mceff_eta[jtID] = new TH2F( hTitle.c_str(), "refeta v recoeta"       , 100, -5.,5., 100, -5.,5.    );
	else if(hMCEffTitleArray[k]=="phi")						          
	  hpp_mceff_phi[jtID] = new TH2F( hTitle.c_str(), "refphi v recophi"       , 100, -5.,5., 100, -4.,4.    );
	else if(hMCEffTitleArray[k]=="ptrat_drjt")					          
	  hpp_mceff_ptrat_drjt[jtID] = new TH2F( hTitle.c_str(), "refdrjt v recopt/refpt" , 100, 0.,5., 30, 0. ,0.6    );   
	else if(hMCEffTitleArray[k]=="ptrat_eta")						          
	  hpp_mceff_ptrat_eta[jtID]  = new TH2F( hTitle.c_str(), "refeta v recopt/refpt"  , 100, 0.,5., 40, -5.,5.    );
	else if(hMCEffTitleArray[k]=="ptrat_phi")						          
	  hpp_mceff_ptrat_phi[jtID]  = new TH2F( hTitle.c_str(), "refphi v recopt/refpt"  , 100, 0.,5., 40, -4.,4.    );     
      }}  
  }
  
  
  /////   JEC AND JER   /////

  // TH1F * hJEC_check[2][nbins_JEC_ptbins][nbins_eta]={};
  // TH1F * hJER[2][nbins_pt]={};  
  // TH3F * hJEC[2]={};

  //sanity jteta check
  TH1F * hJER_jtetaQA;

  //all eta in bins of genpt
  TH1F * hJER[2][nbins_pt]={};  

  //two specific genpt ranges across eta bins
  TH1F * hJER_eta_30pt50[2][nbins_eta]={};
  TH1F * hJER_eta_150pt200[2][nbins_eta]={};

  //all genpt, and two specifically set genpt bins, across abs eta
  TH1F *hJER_absetabins[2][nbins_abseta];

  TH1F *hJER_absetabins_genptBin1[2][nbins_abseta]; 
  TH1F *hGenPtCutHi1   = new TH1F("hGenPtCutHi1"  , (std::to_string(genptBin1High)).c_str() ,  100,  0., 100.  ); 
  TH1F *hGenPtCutLo1   = new TH1F("hGenPtCutLo1"  , (std::to_string(genptBin1Low)).c_str()  ,  100,  0., 100.  ); 	
  hGenPtCutHi1->Fill((float)genptBin1High);  hGenPtCutLo1->Fill((float)genptBin1Low); 

  TH1F *hJER_absetabins_genptBin2[2][nbins_abseta]; 
  TH1F *hGenPtCutHi2   = new TH1F("hGenPtCutHi2"  , (std::to_string(genptBin2High)).c_str() ,  100,  0., 100.  ); 
  TH1F *hGenPtCutLo2   = new TH1F("hGenPtCutLo2"  , (std::to_string(genptBin2Low)).c_str()  ,  100,  0., 100.  ); 
  hGenPtCutHi2->Fill((float)genptBin2High);  hGenPtCutLo2->Fill((float)genptBin2Low); 

  if(fillJERSHists){        
    hJER_jtetaQA=new TH1F("jtEtaSanityCheckJER", "jetEta Dist Check" , 100, -6., 6.);

    for(int jtID=0;jtID<2;jtID++){
      if(!fillMCJetIDHists&&jtID==1)continue;
      
      

      //int jtID=0;
      //if(fillMCJetIDHists)jtID=1;
      //hJEC[0]= new TH3F( ("hJEC_"+std::to_string(0)+"wJetID").c_str(), ";raw p_{T};#eta;JEC", 500,0.,500., 200,-5.,+5., 300,0.,5.);    
      
      for(int y = 0; y<nbins_eta; ++y){	               
	std::string hTitleJER_30to50="hJER_"+std::to_string(jtID)+"wJetID_"+
	  "etabin"+std::to_string(y)+"_30_pt_50";//"etabin"+std::to_string(etabins[y])+"_30_pt_50";
	hJER_eta_30pt50[jtID][y] = new TH1F(hTitleJER_30to50.c_str(), "recpt/genpt, 30<genpt<50",100, 0, 2.);      
	std::string hTitleJER_150to200="hJER_"+std::to_string(jtID)+"wJetID_"+
	  "etabin"+std::to_string(y)+"_150_pt_200";
	hJER_eta_150pt200[jtID][y] = new TH1F(hTitleJER_150to200.c_str(), "recpt/genpt, 150<genpt<200",100, 0, 2.);          
	
	/*for(int x = 0; x<nbins_JEC_ptbins; ++x) {
	  std::string hTitleJEC_check="hJEC_"+std::to_string(jtID)+"wJetID_check_"+
	  "ptbin"+std::to_string(x)+"_etabin"+std::to_string(y);
	  hJEC_check[jtID][x][y] = new TH1F(hTitleJEC_check.c_str(),
	  Form("rawpt/genpt, %2.0f<genpt<%2.0f, %2.4f<geneta<%2.4f", 
	  JEC_ptbins[x],JEC_ptbins[x+1], etabins[y],etabins[y+1]),
	  100,0.,3.);           }*/ 
      }    
      
      for(int bin = 0; bin<nbins_pt; ++bin) {      
	std::string hTitleJER="hJER_"+std::to_string(jtID)+"wJetID_ptbin"+std::to_string(bin);
	hJER[jtID][bin] = new TH1F(hTitleJER.c_str()," jets, recpt/genpt, genptbins", 100,0.,2.); 
      } 
      
      for(int bin = 0; bin<nbins_abseta; ++bin) {      
	std::string hTitleJER="hJER_"+std::to_string(jtID)+"wJetID_jty_bin"+std::to_string(bin);      
	hJER_absetabins[jtID][bin] = new TH1F(hTitleJER.c_str(),("recpt/genpt, absetaBin"+std::to_string(bin)).c_str(), 100,0.,2.); 
	
	std::string hTitleJER_1="hJER_"+std::to_string(jtID)+"wJetID_jty_etabin"+std::to_string(bin)+
	  "_genptBin1";
	hJER_absetabins_genptBin1[jtID][bin] = new TH1F(hTitleJER_1.c_str(),("absetaBin "+std::to_string(bin)+", "						   
									     "genpt "+std::to_string(genptBin1Low)+" to "+std::to_string(genptBin1High)).c_str(), 100,0.,2.); 
	
	std::string hTitleJER_2="hJER_"+std::to_string(jtID)+"wJetID_jty_etabin"+std::to_string(bin)+
	  "_genptBin2";
	hJER_absetabins_genptBin2[jtID][bin] = new TH1F(hTitleJER_2.c_str(),("absetaBin "+std::to_string(bin)+", "						   
									     "genpt "+std::to_string(genptBin2Low)+" to "+std::to_string(genptBin2High)).c_str(), 100,0.,2.); 
      }
    }
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
  //MC jet variable
  float pthat_F;
  int subid_F[1000];    int refparton_F[1000];
  float refpt_F[1000];  float refeta_F[1000];  float refphi_F[1000];  
  float refdrjt_F[1000];

  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
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
  jetpp[0]->SetBranchAddress( "pthat"   , &pthat_F);
  jetpp[0]->SetBranchAddress( "subid"	, &subid_F     );
  jetpp[0]->SetBranchAddress( "refparton_flavor" , &refparton_F );
  jetpp[0]->SetBranchAddress( "refpt"	, &refpt_F     );
  jetpp[0]->SetBranchAddress( "refeta"  , &refeta_F    );
  jetpp[0]->SetBranchAddress( "refphi"  , &refphi_F    );
  jetpp[0]->SetBranchAddress( "refdrjt"	, &refdrjt_F   );  
  
  // HiEvtAnalyzer
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I; float vz_F;
  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);
  jetpp[1]->SetBranchAddress("vz",&vz_F);
  
  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilter_I, pprimaryvertexFilter_I, puvertexFilter_I;
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);
  jetpp[2]->SetBranchAddress("pVertexFilterCutGtight",&puvertexFilter_I);  
  
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
  
  UInt_t NEvents_read=0;
  NEvents_read = NEvents_allFiles;
  
  std::cout<<"reading "<<NEvents_read<<" events"<<std::endl;   

  //MCTruthResidual* MCResidual = new MCTruthResidual( "pp5");
  int mcclosureInt=0;
  //float jetIDCut_neSum, jetIDCut_phSum;
  //if(tightJetID){     jetIDCut_neSum=0.90;  jetIDCut_phSum=0.90;}
  //else{     jetIDCut_neSum=0.99;  jetIDCut_phSum=0.99;}

  for(UInt_t nEvt = 0; nEvt < NEvents_read; ++nEvt) {//event loop   
    
    if(nEvt%1000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;  
    jetpp[0]->GetEntry(nEvt);    
    h_NEvents_read->Fill(1);
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0     || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0  ) continue;    
    h_NEvents_skimCut->Fill(1);    

    if( fabs(vz_F)>15.              ) continue;
    h_NEvents_vzCut->Fill(1);
        
    // compute weights, this way only works because of the previous vz>15. cut
    float vzWeight=1.;
    float vzStart=minbinValue_vzWeights, vzBinLeftSide=vzStart, vzBinRightSide=vzBinLeftSide+binsize_vzWeights;
    if(doVzWeights){
      for( int i=0; i<nbins_vzWeights ; i++ ) { 
	if(vzBinLeftSide<vz_F && vz_F<=vzBinRightSide) {
	  vzWeight=vzWeights[i];  
	  break; }       
	vzBinLeftSide+=binsize_vzWeights;
	vzBinRightSide+=binsize_vzWeights;  }    
    }

    float evtPthatWeight=1.;    
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){ evtPthatWeight=pthatWeights[i]; }     
    float trigWeight=1.;
    float weight_eS=evtPthatWeight*trigWeight*vzWeight;              
    
    //vz hists
    if(fillMCEvtQAHists){

      hVz->Fill(vz_F, 1.);
      hvzWVz->Fill(vz_F, vzWeight);
      hpthatWVz->Fill(vz_F, evtPthatWeight);
      hWVz->Fill(vz_F, weight_eS);      

      hpthat->Fill(pthat_F, 1.);
      hWpthat->Fill(pthat_F, weight_eS);}
        
    // for event counting + avoiding duplicate fills in dijet hists
    bool hNEvts_withJets_Filled=false;
    bool hNEvts_withJets_kmatCut1_Filled=false, hNEvts_withJets_kmatCut2_Filled=false; 

    for(int jet = 0; jet<nref_I; ++jet){

      // event+jet counting
      h_NJets->Fill(1);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      
      
      //if reco jet w/o matched gen jet, skip.
      float genpt= refpt_F[jet];
      float recpt  = pt_F[jet];
      float receta = eta_F[jet];
      float absreceta=fabs(receta);
      float geneta = refeta_F[jet];
      float absgeneta=fabs(geneta);
      
      if( subid_F[jet]!=0 ) continue;
      else if ( recpt <= jtPtCut    ) continue;                 
      else if ( genpt <= genJetPtCut ) continue;
      else if (absgeneta >= jtEtaCutHi)continue;
      else if (absgeneta < jtEtaCutLo) continue;
      //else if ( absreceta > 4.7    ) continue;

      //bool passesEtaCut;
      //passesEtaCut=false;
      //else passesEtaCut=true;
      //if(!passesEtaCut)continue; 

      float genphi = refphi_F[jet];
      float gendrjt = refdrjt_F[jet];     
      float rawpt  = rawpt_F[jet];
      float recphi = phi_F[jet];
      
      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;

      // 13 TeV JetID criterion
      bool passesJetID=false; //int jtID=0;
      if(fillMCJetIDHists) 
	{
	  if (absreceta<=2.4) 
	    { 
	      if( neSum_F[jet]/rawpt    < 0.99 &&
		  phSum_F[jet]/rawpt    < 0.99 &&
		  numConst              > 1    &&      
		  chSum_F[jet]/rawpt    > 0.00 && 
		  chMult                > 0    &&
		  eSum_F[jet]/rawpt     < 0.99    ) passesJetID=true;	      
	    }
	  else if ( absreceta<=2.7 && absreceta>2.4 ) 
	    {	  
	      if( neSum_F[jet]/rawpt    < 0.99 &&
		  phSum_F[jet]/rawpt    < 0.99 &&
		  numConst              > 1       ) passesJetID=true;	      
	    }		  
	  else if( absreceta<=3.0 && absreceta>2.7 ) 
	    {                                                         // CMSSW [76,80]X criterion
	      if(  phSum_F[jet]/rawpt > 0.00 &&                       // else if(  phSum_F[jet]/rawpt [< 0.90 ] / [ > 0.01 &&]		     
		   neSum_F[jet]/rawpt > 0.00 &&                       //           neSum_F[jet]/rawpt [null   ] / [ < 0.98 &&]		     
		   neuMult            > 1       ) passesJetID=true;   //           neuMult            [> 2    ] / [ > 2      ] ) passesJetID=true;
	    }							      
	  else //( absreceta>3.0) 
	    {                                                          // CMSSW 76X criterion
	      if( (phSum_F[jet]/rawpt > 0.00 ||                         // else if( phSum_F[jet]/rawpt < 0.90 &&
		   neSum_F[jet]/rawpt  > 0.00 ) &&                         //          neSum_F[jet]/rawpt < null &&
	      	  neuMult             > 0       ) passesJetID=true;     //          neuMult            > 10
	    }	  	  
	}
      
      
      // jet/event counts
      h_NJets_kmatCut1->Fill(1);
      if(!hNEvts_withJets_kmatCut1_Filled){
	h_NEvents_withJets_kmatCut1->Fill(1);
	hNEvts_withJets_kmatCut1_Filled=true;      }      
      
      /////   JERS HISTS   ///// 
      if(fillJERSHists ){
	
	// dual-diff xsec abseta bins for two spec genpt bins and all genpt bins
	int absetabin=-1;
	if(absgeneta>absetabins[nbins_abseta]) absetabin = -1;//check end of array
	else{
	  for(int bin=0; bin < nbins_abseta; bin++)
	    if( absgeneta>=absetabins[bin] ) absetabin = bin; 				
	}
	
	if(absetabin != -1 ){

	  for(int jtID=0; jtID<2;jtID++){
	    if(jtID==1 && !passesJetID) continue;
	    else if (jtID==1 && !fillMCJetIDHists)continue;
	    
	    hJER_absetabins[jtID][absetabin]->Fill( (float)(recpt/genpt), weight_eS);     
	    
	    if( genpt >= ( (float) genptBin1Low  )  && 
		genpt <  ( (float) genptBin1High )    ) 
	      hJER_absetabins_genptBin1[jtID][absetabin]->Fill( (float)(recpt/genpt), weight_eS);     
	    
	    if( genpt >= ( (float) genptBin2Low  )  && 
		genpt <  ( (float) genptBin2High )    ) 
	      hJER_absetabins_genptBin2[jtID][absetabin]->Fill( (float)(recpt/genpt), weight_eS);     	
	  }
	  
	}
	
	// hJER, genpt30-50, 150-200, across geneta bins array (note, NOT abs-value geneta bins)
	int etabin = -1;
	if(fabs(geneta) >= etabins[nbins_eta]) etabin=-1; //check end of array
	else {
	  for(int bin = 0; bin<nbins_eta; ++bin)  //notice that the loop does not hit end of this array
	    if(geneta >= etabins[bin]) etabin = bin;	  }
	
	if(etabin != -1) {
	  
	  for(int jtID=0; jtID<2;jtID++){
	    if(jtID==1 && !passesJetID)continue;
	    else if (jtID==1 && !fillMCJetIDHists)continue;
	    
	    if(genpt >= 30.  && genpt<50. ) 
	      hJER_eta_30pt50[jtID][etabin]->Fill((float)recpt/genpt, weight_eS);
	    else if(genpt >= 150. && genpt<200.) 
	      hJER_eta_150pt200[jtID][etabin]->Fill((float)recpt/genpt, weight_eS);	
	  }
	}
	
	//JER in genptbins array, according to input absetacuts
	int ptbin = -1; 
	if( genpt >= ptbins[nbins_pt] ) ptbin=-1;//check end of binning array    //else if (!passesEtaCut) ptbin=-1;
	else {
	  for(int bin = 0; bin<nbins_pt; ++bin)      
	    if(genpt >= ptbins[bin]) ptbin = bin;  	}
	
	if(ptbin != -1) {
	  for(int jtID=0; jtID<2;jtID++){
	    if(jtID==1 && !passesJetID)continue;
	    else if (jtID==1 && !fillMCJetIDHists)continue;
	    if(jtID==0)hJER_jtetaQA->Fill(receta,weight_eS);      
	    hJER[jtID][ptbin]->Fill( (float)(recpt/genpt), weight_eS); 
	  }	}
	
      }
      
      
      /////   GEN/RECO MATCHING   /////
      
      if(fillMCEffHists){	
	for(int jtID=0; jtID<2;jtID++){
	  if(jtID==1 && !passesJetID)continue;
	  else if (jtID==1 && !fillMCJetIDHists)continue;
	  hpp_mceff_pt[jtID]->Fill(         recpt/genpt   , genpt   ,  weight_eS);
	  hpp_mceff_eta[jtID]->Fill(        receta , geneta  ,  weight_eS);
	  hpp_mceff_phi[jtID]->Fill(        recphi , genphi  ,  weight_eS);	
	  hpp_mceff_ptrat_drjt[jtID]->Fill( recpt/genpt   , gendrjt ,  weight_eS);	
	  hpp_mceff_ptrat_eta[jtID]->Fill(  recpt/genpt   , geneta  ,  weight_eS);
	  hpp_mceff_ptrat_phi[jtID]->Fill(  recpt/genpt   , genphi  ,  weight_eS);       
	}
      }
      
      /////   UNFOLDING   ///// 
      if(fillMCUnfoldingHists){
	for(int jtID=0; jtID<2;jtID++){

	  if(jtID==1 && !passesJetID)continue;
	  else if (jtID==1 && !fillMCJetIDHists)continue;

	  hpp_gen[jtID]->Fill(genpt, weight_eS);
	  hpp_reco[jtID]->Fill(recpt, weight_eS);
	  hpp_matrix[jtID]->Fill(recpt, genpt, weight_eS);

	  if(mcclosureInt%2 == 0){
	    hpp_mcclosure_reco_test[jtID]->Fill(recpt, weight_eS);      }
	  else {
	    hpp_mcclosure_matrix[jtID]->Fill(recpt, genpt, weight_eS);
	    hpp_mcclosure_gen[jtID]->Fill(genpt, weight_eS);
	    hpp_mcclosure_reco[jtID]->Fill(recpt, weight_eS); 	} 
	}
      }
      
      // jet/event counts
      if(passesJetID){
	h_NJets_kmatCut2->Fill(1);
	if(!hNEvts_withJets_kmatCut2_Filled){
	  h_NEvents_withJets_kmatCut2->Fill(1);
	  hNEvts_withJets_kmatCut2_Filled=true;      }}
      
      
    }//end jet loop
    mcclosureInt++;    
  }//end event loop
  
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Event-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;
  
  
  std::cout<<"Total Num of Events read passing skimCuts                              = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Total Num of good Events w/ jets                             = " <<
    h_NEvents_withJets->GetEntries()<<std::endl;	   
  std::cout<<"Total Num of good Events, w/ jets, post kmatCut1             = " <<
    h_NEvents_withJets_kmatCut1->GetEntries()<<std::endl; 
  std::cout<<"Total Num of good Events, w/ jets, post kmatCut2             = " <<
    h_NEvents_withJets_kmatCut2->GetEntries()<<std::endl ;
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Jets read from good events                          = " <<
    h_NJets->GetEntries()<<std::endl;  
  std::cout<<"Total Num of Jets read from good events post kmatCut1              = " <<
    h_NJets_kmatCut1->GetEntries()<<std::endl;
  std::cout<<"Total Num of Jets read from good events post kmatCut2              = " <<
    h_NJets_kmatCut2->GetEntries()<<std::endl;
  std::cout<<std::endl;
  

  std::cout<<"writing output file "<<outfile<<std::endl;
  fout->Write(); 


  std::cout<<std::endl<<"readForests_ppMC_JERS finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;

return 0;
}// end readForests_ppMC_JERS
								 

////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=10 && argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForests_ppMC_JERS.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppMC_JERS.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> ";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppMC_JERS();
  else{//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    if(argc==minArgs){
      rStatus = readForests_ppMC_JERS( inputFileList, startfile, endfile, 
				       jetRadius, jetType, debug,
				       outputFileName);                              }    
    else if(argc==(minArgs+2)){
      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
      rStatus = readForests_ppMC_JERS( inputFileList, startfile, endfile,
				       jetRadius, jetType, debug,
				       outputFileName, jtEtaCutLo, jtEtaCutHi);    }
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
