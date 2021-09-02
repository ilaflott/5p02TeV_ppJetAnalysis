// custom header
#include "readForests_JERS.h"

// ppMC switches
const bool fillMCEvtQAHists=true;
const bool fillJERSHists=true;
const bool fillMCJetIDHists=true;//, tightJetID=false;
const int jetIDint=(int)fillMCJetIDHists;

const bool verbose=false;
const bool useTightJetID=false;

const bool usePFCandAnalyzer=true;
//const bool useTupel=true;

const bool REDOjetMatching=true;
const float deltaR_req4Match=0.2;

//// readForests_ppMC_JERS
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_JERS(std::string inFilelist , int startfile , int endfile ,  
			  int radius , std::string jetType , bool debugMode ,
			  std::string outfile  ,
			  float jtEtaCutLo, float jtEtaCutHi ){ 
  
  // for monitoring performance + debugging
  TStopwatch timer;  timer.Start();
  const bool debugMode_verbose=debugMode&&verbose; 
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  if(debugMode_verbose)std::cout<<std::endl<<"WARNING debugMode_verbose is ON, LARGE AMOUNTS OF TEXT INCOMING"<<std::endl;
  
  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading filelist "<< inFilelist<<std::endl;
  std::cout<<"reading files #'s "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius = " << radius;
  std::cout<<", jetType = " << jetType<<std::endl;
  
  
  std::cout<<"///////////////////"<<std::endl<<std::endl;
  
  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  
  // form jetTreeName
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;
  
  // initialize tree name array
  std::string trees[N_MCTrees];
  trees[0]=jetTreeName; //trees[1]=CaloJetTreeName;
  for(int i=1;i<N_MCTrees;++i)trees[i]=MCTreeNames[i];
  
  // declare TChains for each tree + friend  
  TChain* jetpp[N_MCTrees]={};
  for(int t = 0;t<N_MCTrees;++t)  {     
    if(!usePFCandAnalyzer && t==4)continue;
    //if(!useTupel && t==5)continue;
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
      if(!usePFCandAnalyzer && t==4)continue;
      //if(!useTupel &&t==5)continue;
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

  TH1F *hGenDRCut     =new TH1F("hGenDRCut"  ,(std::to_string(gendrCut)).c_str()   ,    100, 0.,100.); hGenDRCut->Fill(gendrCut);           
  
  TH1F *hJetEtaCutHi       =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6  ); hJetEtaCutHi->Fill(jtEtaCutHi);
  TH1F *hJetEtaCutLo       =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6  ); hJetEtaCutLo->Fill(jtEtaCutLo);	
  
  TH1D *h_NEvents         = new TH1D("NEvents","NEvents", 1,0,2);
  TH1D *h_NEvents_read    = new TH1D("NEvents_read","NEvents read", 1,0,2);
  TH1D *h_NEvents_skimCut = new TH1D("NEvents_skimCut"      ,"NEvents read post skimCut", 1,0,2);
  TH1D *h_NEvents_vzCut   = new TH1D("NEvents_vzCut"        ,"NEvents read post vzCut AND skimCut", 1,0,2);
  
  //counts post all evt cuts, that are thrown out because of no good jets pre and post jet cuts
  TH1D *h_NEvents_withJets           = new TH1D("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1D *h_NEvents_withJets_kmatCut1   = new TH1D("NEvents_withJets_kmatCut1" , 
						 "NEvents read post evt cuts, w/ jets post kmatCut1", 1,0,2);
  TH1D *h_NEvents_withJets_kmatCut2   = new TH1D("NEvents_withJets_kmatCut2" , 
						"NEvents read post evt cuts, w/ jets post kmatCut2", 1,0,2);
  
  //jet counts
  TH1D *h_NJets          = new TH1D("NJets","NJets read", 1,0,2);
  TH1D *h_NJets_kmatCut1  = new TH1D("NJets_kmatCut1 ","NJets read post kmatCut1 ", 1,0,2);
  TH1D *h_NJets_kmatCut2  = new TH1D("NJets_kmatCut2 ","NJets read post kmatCut2 ", 1,0,2);
  
  
  /////   EVENT QA   ///// 
  TH1D *hVr=NULL, *hWVr=NULL;
  TH1D *hVx=NULL, *hWVx=NULL;
  TH1D *hVy=NULL, *hWVy=NULL;
  TH1D *hVz=NULL, *hWVz=NULL, *hpthatWVz=NULL, *hvzWVz=NULL; 

  TH1D *hpthat=NULL, *hWpthat=NULL;  
  if(fillMCEvtQAHists){
    hpthatWVz = new TH1D("hpthatWeightedVz","", 500,-25.,25.);//pthat-weighted evtvz
    hvzWVz    = new TH1D("hvzWeightedVz","", 500,-25.,25.);//vz-weighted evtvz    
    hVr = new TH1D("hVr","vr, no trig, no weights",         2000,0.,0.60); 
    hWVr = new TH1D("hWeightedVr","vr, trigd, with weights",2000,0.,0.60);     
    hVx = new TH1D("hVx","vx, no trig, no weights",  2000, -0.30,0.30);
    hWVx = new TH1D("hWeightedVx","vx, trigd, with weights", 2000, -0.30,0.30); 
    hVy = new TH1D("hVy","vy, no trig, no weights",   2000, -0.30,0.30); 
    hWVy = new TH1D("hWeightedVy","vy, trigd, with weights",   2000,-0.30,0.30);    
    hVz       = new TH1D("hVz","", 1000,-25.,25.);//evtvz
    hWVz      = new TH1D("hWeightedVz","", 1000,-25.,25.);//pthat*vz-weighted evt vz
    hpthat    = new TH1D("hpthat","",1000,0,1000);//evt pthat, unweighted and weighted
    hWpthat   = new TH1D("hWeightedpthat","",1000,0,1000);    }
  
  
  
  
  /////   JEC AND JER   /////
  
  // TH1D * hJEC_check[nbins_JEC_ptbins][nbins_eta]={};
  // TH1D * hJER[nbins_pt]={};  
  // TH3F * hJEC={};
  
  //sanity jteta check
  TH1D * hJER_jtetaQA;
  TH1D * hJER_jtyQA;
  TH1D * hJER_genetaQA;
  TH1D * hJER_jtptQA[nbins_abseta]={};
  TH1D * hJER_genptQA[nbins_abseta]={};
  
  //all eta in bins of genpt
  TH1D * hJER[nbins_pt]={};       // rec/gen , bins of gen pt
  //TH1D * hJER_abseta[nbins_pt][nbins_abseta]={};       // rec/gen , bins of gen pt, receta
  TH1D * hJER_absgeneta[nbins_pt][nbins_abseta]={};       // rec/gen , bins of gen pt, geneta
  //TH1D * hJER_absy[nbins_pt][nbins_absy]={};       // rec/gen , bins of gen pt, recy
  TH1D * hJER_absgeny[nbins_pt][nbins_absy]={};       // rec/gen , bins of gen pt, geny

  if(fillJERSHists){   
    
    hJER_jtetaQA=new TH1D("jtEtaSanityCheckJER", "jetEta Dist Check" , 100, -6., 6.);
    hJER_jtyQA=new TH1D("jtYSanityCheckJER", "jetY Dist Check" , 100, -6., 6.);
    hJER_genetaQA=new TH1D("genEtaSanityCheckJER", "genEta Dist Check" , 100, -6., 6.);


    for(int j=0;j<nbins_abseta;++j){
      std::string absetabins_i=absetabins_str[j];	std::string absetabins_ip1=absetabins_str[j+1];	  
      std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID"+"_jtptQA_etabin"+std::to_string(j);
      if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
      std::string hDescJER="rec pt, "+absetabins_i+"<|Eta|<"+absetabins_ip1;
      if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
      hJER_jtptQA[j] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 1000,0.,1000.); 	  
    }//etabin loop
    
    for(int j=0;j<nbins_abseta;++j){
      std::string absetabins_i=absetabins_str[j];	std::string absetabins_ip1=absetabins_str[j+1];	  
      std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID"+"_genptQA_etabin"+std::to_string(j);
      if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
      std::string hDescJER="gen pt, "+absetabins_i+"<|#eta|<"+absetabins_ip1;
      if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
      hJER_genptQA[j] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 1000,0.,1000.); 	  
    }//etabin loop
    

    for(int bin = 0; bin<nbins_pt; ++bin) {      
      std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );
      std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
      std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID_"+ptbins_i+"ptbin"+ptbins_ip1;
      if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
      std::string hDescJER="rec/gen pt, "+ptbins_i+"<genpt<"+ptbins_ip1;
      if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
      hJER[bin] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 200,0.,2.); 
    }//ptbin loop 
    
    //for(int j=0;j<nbins_abseta;++j){
    //  for(int bin = 0; bin<nbins_pt; ++bin) {      	  
    //	std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );	std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
    //	std::string absetabins_i=absetabins_str[j];	std::string absetabins_ip1=absetabins_str[j+1];	  
    //	std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID"+"_etabin"+std::to_string(j)+"_"+ptbins_i+"ptbin"+ptbins_ip1;
    //	if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
    //	std::string hDescJER="rec/gen pt, "+absetabins_i+"<|#eta|<"+absetabins_ip1+", "+ptbins_i+"<genpt<"+ptbins_ip1;	  
    //	if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
    //	hJER_abseta[bin][j] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 200,0.,2.); 
    //  }//ptbin loop
    //}//etabin loop

    for(int j=0;j<nbins_abseta;++j){
      for(int bin = 0; bin<nbins_pt; ++bin) {      	  
	std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );	std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
	std::string absetabins_i=absetabins_str[j];	std::string absetabins_ip1=absetabins_str[j+1];	  
	std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID"+"_genetabin"+std::to_string(j)+"_"+ptbins_i+"ptbin"+ptbins_ip1;
	if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
	std::string hDescJER="rec/gen pt, "+absetabins_i+"<|#eta^{GEN}|<"+absetabins_ip1+", "+ptbins_i+"<genpt<"+ptbins_ip1;	  
	if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
	hJER_absgeneta[bin][j] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 200,0.,2.); 
      }//ptbin loop
    }//etabin loop
    
    //for(int j=0;j<nbins_absy;++j){
    //  for(int bin = 0; bin<nbins_pt; ++bin) {      	  
    //	std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );	std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
    //	std::string absybins_i=absybins_str[j];	std::string absybins_ip1=absybins_str[j+1];	  
    //	std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID"+"_ybin"+std::to_string(j)+"_"+ptbins_i+"ptbin"+ptbins_ip1;
    //	if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
    //	std::string hDescJER="rec/gen pt, "+absybins_i+"<|y|<"+absybins_ip1+", "+ptbins_i+"<genpt<"+ptbins_ip1;	  
    //	if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
    //	hJER_absy[bin][j] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 200,0.,2.); 
    //  }//ptbin loop
    //}//etabin loop        
    
    for(int j=0;j<nbins_absy;++j){
      for(int bin = 0; bin<nbins_pt; ++bin) {      	  
	std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );	std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
	std::string absybins_i=absybins_str[j];	std::string absybins_ip1=absybins_str[j+1];	  
	std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID"+"_genybin"+std::to_string(j)+"_"+ptbins_i+"ptbin"+ptbins_ip1;
	if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
	std::string hDescJER="rec/gen pt, "+absybins_i+"<|y^{GEN}|<"+absybins_ip1+", "+ptbins_i+"<genpt<"+ptbins_ip1;	  
	if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
	hJER_absgeny[bin][j] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 200,0.,2.); 
      }//ptbin loop
    }//etabin loop
    
    
  }//fill JERShists

  
  
  // EVENT LOOP PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  float eta_F[1000];   float phi_F[1000];
  float rawpt_F[1000]; float m_F[1000];
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
  //MC reco jet variables
  float pthat_F; 
  int subid_I[1000];    int refparton_I[1000];
  float refpt_F[1000];  float refeta_F[1000];  float refphi_F[1000];  
  float refdrjt_F[1000];
  //MC gen jet variables
  int ngen_I;  int genmatchindex_I[1000];
  float genpt_F[1000];  float geneta_F[1000];  float genphi_F[1000];
  //evt pthat used to generate jets, MC specific
  jetpp[0]->SetBranchAddress( "pthat"   , &pthat_F);
  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jtm",&m_F);
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
  //MC specific matched jet variables
  jetpp[0]->SetBranchAddress( "subid"	, &subid_I     );
  jetpp[0]->SetBranchAddress( "refparton_flavor" , &refparton_I );
  jetpp[0]->SetBranchAddress( "refpt"	, &refpt_F     );
  jetpp[0]->SetBranchAddress( "refeta"  , &refeta_F    );
  jetpp[0]->SetBranchAddress( "refphi"  , &refphi_F    );
  jetpp[0]->SetBranchAddress( "refdrjt"	, &refdrjt_F   );  
  //gen jet stuff; for making "reco match index" array
  jetpp[0]->SetBranchAddress( "ngen" , &ngen_I );
  jetpp[0]->SetBranchAddress( "genmatchindex", &genmatchindex_I);
  jetpp[0]->SetBranchAddress( "genpt" , &genpt_F);
  jetpp[0]->SetBranchAddress( "geneta" , &geneta_F);
  jetpp[0]->SetBranchAddress( "genphi" , &genphi_F);
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
  int pBeamScrapingFilter_I, pHBHENoiseFilterRun2Loose_I, pprimaryvertexFilter_I, puvertexFilter_I;
  int pHBHEIsoNoiseFilter_I;
  jetpp[2]->SetBranchAddress("HBHEIsoNoiseFilterResult",&pHBHEIsoNoiseFilter_I);//only option
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilterRun2Loose_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);
  jetpp[2]->SetBranchAddress("pVertexFilterCutGtight",&puvertexFilter_I);  


  //pfcandAnalyzer/pfTree
  Int_t nPFpart_I;  
  //std::vector<Float_t> *pfpt, *pfphi, *pfeta, *pfE;
  Float_t pfptArr[2000], pfphiArr[2000], pfetaArr[2000], pfEArr[2000];
  
  if(usePFCandAnalyzer){
    if(debugMode)std::cout<<"using PFCandAnalyzer. Setting branch addresses!"<<std::endl;
    jetpp[4]->SetBranchAddress( "nPFpart", &nPFpart_I);
    if(debugMode)std::cout<<"filelistIsLowerJets, pointing branches to Float_t* (a.k.a. Float_t array)"<<std::endl;
    jetpp[4]->SetBranchAddress( "pfPt"   , &pfptArr );
    jetpp[4]->SetBranchAddress( "pfEta"  , &pfetaArr);
    jetpp[4]->SetBranchAddress( "pfPhi"  , &pfphiArr);    
    jetpp[4]->SetBranchAddress( "pfEnergy"  , &pfEArr);    
  }

  //tupel/EventTree (for calculating genjet rapidity for AK4 only [R=0.3 jet info absent])
  std::vector<float> *jetPt, *jetEta, *jetPhi, *jetE;// *jetID, *jetRawPt, *jetRawE;
  std::vector<float> *GjetPt, *GjetEta, *GjetPhi, *GjetE;// *jetID, *jetRawPt, *jetRawE;
  int evtvtxCnt_I;
  //std::vector<float> *jetHfHadE, *jetHfEmE, *jetChHadFrac, *jetNeutralHadAndHfFrac, *jetChEmFrac, *jetNeutralEmFrac, *jetChMult, *jetConstCnt;  
  if(debugMode)std::cout<<"setting branch addresses for tupel/EventTree"<<std::endl;
  
  jetpp[5]->SetBranchAddress( "JetAk04Pt",                  &jetPt                  );
  jetpp[5]->SetBranchAddress( "JetAk04Eta",                 &jetEta                 );
  jetpp[5]->SetBranchAddress( "JetAk04Phi",                 &jetPhi                 );
  jetpp[5]->SetBranchAddress( "JetAk04E",                   &jetE                   );
  jetpp[5]->SetBranchAddress( "GJetAk04Pt",                 &GjetPt                  );
  jetpp[5]->SetBranchAddress( "GJetAk04Eta",                &GjetEta                 );
  jetpp[5]->SetBranchAddress( "GJetAk04Phi",                &GjetPhi                 );
  jetpp[5]->SetBranchAddress( "GJetAk04E",                  &GjetE                   );
  jetpp[5]->SetBranchAddress("EvtVtxCnt",      &evtvtxCnt_I);
  
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
  int goodgenjetcount=0,badgenjetcount=0;  
  std::cout<<"reading "<<NEvents_read<<" events"<<std::endl;   

  for(UInt_t nEvt = 0; nEvt < NEvents_read; ++nEvt) {//event loop   
    
    if(nEvt%1000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;  
    jetpp[0]->GetEntry(nEvt);    
    h_NEvents_read->Fill(1);
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilterRun2Loose_I==0     || 
	pHBHEIsoNoiseFilter_I    ==0    || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0  ) continue;    
    h_NEvents_skimCut->Fill(1);    

    if( fabs(vz_F)>24.              ) continue;
    h_NEvents_vzCut->Fill(1);

    if( sqrt( pow(vx_F,2) + pow(vy_F,2) ) > 0.1) continue;//DONT FORGET I AM HERE
    if(evtvtxCnt_I>1)continue;//DONT FORGET I AM HERE

    if(usePFCandAnalyzer){
      //PFCANDLOOP
      //std::cout<<std::endl<<"nPFpart_I="<<nPFpart_I<<std::endl;
      
      Float_t sumpfpt_x=0.,sumpfpt_y=0.;
      Float_t sumpfEt=0.;//scalar sum of PF Et, not a vector sum.
      //h_nPFpart->Fill(nPFpart_I,weight_eS);
      
      for(Int_t part=0;part<nPFpart_I;part++){
	
	Float_t pfpt_part,pfeta_part,pfphi_part, pfE_part;
	pfpt_part  = pfptArr[part];	  
	pfeta_part =pfetaArr[part];	  
	pfphi_part =pfphiArr[part];	  	
	pfE_part  =pfEArr[part];	
	
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
      
      //h_PFMETfrac->Fill(PFMETfrac,weight_eS);
      if(PFMETfrac>0.3)continue;//giving this a shot. see how it goes...
      //if(PFMETfrac>0.5)continue;//giving this a shot. see how it goes...
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
    
    // grab vzweight
    double vzWeight=1.;
    if(doVzWeights&&false){
      //std::cout<<"vzWeight was:"<<vzWeight<<std::endl;
      vzWeight=cpuVzWeight_poly(vz_F);
    }
    
    double evtPthatWeight=1.;    
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){ 
      evtPthatWeight=pthatWeights[i]; 
      evtPthatWeight*=1e+06;//in nanobarns
    }     
    
    double trigWeight=1.;
    //no trig selection, trig weight always 1

    double weight_eS=evtPthatWeight*trigWeight*vzWeight;              
    float vr_F=std::sqrt( vx_F*vx_F + vy_F*vy_F);
    //vz hists
    if(fillMCEvtQAHists){

      hVr->Fill( vr_F, 1.);
      hWVr->Fill(vr_F, weight_eS);      
      hVx->Fill( vx_F, 1.);
      hWVx->Fill(vx_F, weight_eS);      
      hVy->Fill( vy_F, 1.);
      hWVy->Fill(vy_F, weight_eS);      

      hVz->Fill(vz_F, 1.);
      hWVz->Fill(vz_F, weight_eS);      
      hvzWVz->Fill(vz_F, vzWeight);
      hpthatWVz->Fill(vz_F, evtPthatWeight);

      hpthat->Fill(pthat_F, 1.);
      hWpthat->Fill(pthat_F, weight_eS); }
    
    // for event counting + avoiding duplicate fills in dijet hists
    bool hNEvts_withJets_Filled=false;
    bool hNEvts_withJets_kmatCut1_Filled=false, hNEvts_withJets_kmatCut2_Filled=false; 

    std::vector<int>recomatchindex_I;
    recomatchindex_I.assign(nref_I,-1);//set nref_I elements of vector to -1, indicating that reco jet N doesn't yet have matching genjet

    bool tempdebug=false;//debugMode&&(nEvt%200==0);    
    if(tempdebug)std::cout<<std::endl<<"_______ nEvt="<<nEvt<<" ________"<<std::endl;
    if(tempdebug)std::cout<<"ngen_I="<<ngen_I<<std::endl;
    for(int gjet=0; gjet<ngen_I; ++gjet){
      if(tempdebug)std::cout<<"// genmatchindex_I["<<gjet<<"]="<<genmatchindex_I[gjet]<<std::endl;
      if(tempdebug)std::cout<<"genpt_F["<<gjet<<"]="<<genpt_F[gjet]<<std::endl;
      if((genmatchindex_I[gjet]!=-1) &&
	 (genmatchindex_I[gjet]<nref_I)) recomatchindex_I[genmatchindex_I[gjet]]=gjet;
    }
    //if(tempdebug)std::cout<<std::endl;
    //if(tempdebug)std::cout<<"nref_I="<<nref_I<<std::endl;
    //for(int jet=0;jet<nref_I;++jet){
    //  if(tempdebug)std::cout<<"// recomatchindex_I["<<jet<<"]="<<recomatchindex_I.at(jet)<<std::endl;
    //  if(tempdebug)std::cout<<"subid_I["<<jet<<"]="<<subid_I[jet]<<std::endl;
    //  if(tempdebug)std::cout<<"refpt["<<jet<<"]="<<refpt_F[jet]<<std::endl;
    //  if(recomatchindex_I[jet]!=-1)
    //	if(tempdebug){
    //	  std::cout<<"genpt["<<recomatchindex_I[jet]<<"]="<<genpt_F[recomatchindex_I[jet]]<<std::endl;
    //	  std::cout<<"GJetAk04JetPt["<<recomatchindex_I[jet]<<"]="<<GjetPt->at(recomatchindex_I[jet])<<std::endl;
    //	}
    //}





    









    /*

    std::vector<bool> refjetIsSorted;  //will have length of nref_I    
    std::vector<double> refpt_sorted, refeta_sorted, refphi_sorted; //these have length of #of matched GEN jet
    std::vector<double> refy_sorted;
    
    std::vector<bool> refjetIsMatched; std::vector<int> refjetIsMatchedTo; 
    std::vector<bool> recjetIsMatched  ;       //length of nref_I
    std::vector<int> recjetIsMatchedTo; 
    if(REDOjetMatching){
      //if(nEvt>10000)assert(false);


      //I SHOULD BE USING INFO FROM THE TUPEL TREE
      //PROBABLY NEED TO DO SOME REWRITING OF ALL OF THIS

      //first, copy+reorder the refpt/refeta/refphi array from highest --> lowest in pT
      refjetIsSorted.assign(nref_I, false);
      int Nsorted=0, NsortAttempts=0;
      //std::cout<<"about to start do-while loop"<<std::endl;
      do{//sort
	//std::cout<<std::endl<<"now in do-while loop"<<std::endl;
	float maxjtpt=-1.;
	int maxjtpt_ind=-1;
	for(int jet=0; jet<nref_I; ++jet){//jet loop for 1 sort iteration
	  //std::cout<<"(in sort loop) checking jet#"<<jet<<" with pt="<<refpt_F[jet]<<std::endl;

	  if(refjetIsSorted.at(jet) || refpt_F[jet]<0.){
	    //std::cout<<"refjet is already sorted into array. moving on."<<std::endl;
	    continue;
	  }
	  else if(refpt_F[jet]>maxjtpt){
	    //std::cout<<"higher refpt found! was jet #"<< maxjtpt_ind<<" and pt was "<<maxjtpt;
	    maxjtpt=refpt_F[jet];
	    maxjtpt_ind=jet;
	    //std::cout<<", now is jet #"<<maxjtpt_ind<<" and pt is " <<maxjtpt<<std::endl;
	  }

	}//end jet loop for 1 sort iteration
	//std::cout<<"done with sort-loop in do-while statement. maxjtpt_ind="<<maxjtpt_ind<<", maxjtpt="<<maxjtpt<<std::endl;
	
	if(maxjtpt_ind>=0){//we found a highest jet pT!
	  //std::cout<<"highest pT jet found!"<<std::endl;
	  refpt_sorted.push_back(refpt_F[maxjtpt_ind]);
	  //std::cout<<"checking std::vector element assignment: refpt_sorted.at(0)="<<refpt_sorted.at(0)<<std::endl;
	  refeta_sorted.push_back(refeta_F[maxjtpt_ind]);
	  refphi_sorted.push_back(refphi_F[maxjtpt_ind]);
	  refjetIsSorted.at(maxjtpt_ind)=true;
	  Nsorted++;
	}
	
	NsortAttempts++;//every iteration finds the highest pt jet.
      }//end do sort
      while(NsortAttempts<nref_I);

      //TODO




      
      //here; i need to look for gen jets in the genjet arrays (GEN jets that went unmatched)
      //i expect the impact from this to be largest as *very* low pT, decreasing to negligible elsewhere.
      //so i will finish and debug the other steps first, and double back to this next.
      bool insertedagenjet=false; bool foundunmatchedgenjet=false;
      for(int gjet=0;gjet<ngen_I;++gjet){
	if(genmatchindex_I[gjet]>=0)continue;//already matched, so already written to the refjet array.
	else{
	  if(!foundunmatchedgenjet)
	    std::cout<<"nEvt="<<nEvt<<", FOUND UNMATCHED GENJET(S)"<<std::endl;//assert(false);
	  foundunmatchedgenjet=true;
	  std::cout<<"genpt/eta/phi_F["<<gjet<<"]="<<genpt_F[gjet]<<"/"<<geneta_F[gjet]<<"/"<<genphi_F[gjet]<<std::endl;
	  for(int refjet=0;refjet<Nsorted;++refjet){
	    std::cout<<"refjet="<<refjet<<std::endl;
  
	    if(refpt_sorted.at(refjet)<genpt_F[gjet]){//if it's not the lowest pTjet, and IT IS BIGGER than this ref pT, insert. 
	      if(!insertedagenjet){
		std::cout<<std::endl<<"------------BEFORE LOOKING FOR MISSED GEN JETS--------------------------"<<std::endl;
		for(int jet=0; jet<Nsorted; ++jet){
		  std::cout<<"sorted refpt/eta/phi["<<jet<<"]="<< refpt_sorted.at(jet)  << " / " << 
		    refeta_sorted.at(jet) << " / " << 
		    refphi_sorted.at(jet) << std::endl;
		}
		std::cout<<std::endl;
	      }
	      refpt_sorted .insert( refpt_sorted.begin()+refjet, genpt_F[gjet]);
	      refeta_sorted.insert(refeta_sorted.begin()+refjet, geneta_F[gjet]);
	      refphi_sorted.insert(refphi_sorted.begin()+refjet, genphi_F[gjet]);
	      Nsorted++;
	      insertedagenjet=true;
	      break;	    }//end finding a valid unmatched genjet to redo the matching with	    

	    else if(refjet==(Nsorted-1)){//boundary case; if the unmatched genjet is the lowest pT jet of any put in array so far.
	      if(!insertedagenjet){
		std::cout<<std::endl<<"------------BEFORE LOOKING FOR MISSED GEN JETS--------------------------"<<std::endl;
		for(int jet=0; jet<Nsorted; ++jet){
		  std::cout<<"sorted refpt/eta/phi["<<jet<<"]="<< refpt_sorted.at(jet)  << " / " << 
		    refeta_sorted.at(jet) << " / " << 
		    refphi_sorted.at(jet) << std::endl;
		}
		std::cout<<std::endl;
	      }
	      refpt_sorted .push_back( genpt_F[gjet]);
	      refeta_sorted.push_back( geneta_F[gjet]);
	      refphi_sorted.push_back( genphi_F[gjet]);
	      Nsorted++;
	      insertedagenjet=true;
	      break;//necessary because the size of the array increases by 1, so the for-loop condition won't fire.
	    }//end boundary case
	    //else if(refpt_sorted.at(refjet)>genpt_F[gjet]){//if it's not the lowest pT jet, and it's not bigger than this ref pT, move on.
	    //  continue;
	    //}	

	    

	  }//end refjet loop for inserting a valid genjet into the array.
	  
	  


	}
      }//end genjet loop looking for unmatched genjets
      

	if(insertedagenjet){
	std::cout<<std::endl<<"------------AFTER FINDING AT LEAST 1 MISSED GEN JET--------------------------"<<std::endl;
	for(int jet=0; jet<Nsorted; ++jet){
	  std::cout<<"sorted refpt/eta/phi["<<jet<<"]="<< refpt_sorted.at(jet)  << " / " << 
	  refeta_sorted.at(jet) << " / " << 
	  refphi_sorted.at(jet) << std::endl;
	  }
	  std::cout<<std::endl;
	  }

      //TODO END INSERTING MISSED GEN JETS




      ////BEGIN DEBUG this test was fine.
      //std::cout<<std::endl<<"--------------------------------------"<<std::endl;
      //for(int jet=0; jet<Nsorted; ++jet){
      //	std::cout<<"sorted refpt/eta/phi["<<jet<<"]="<< refpt_sorted.at(jet)  << " / " << 
      //	                                                refeta_sorted.at(jet) << " / " << 
      //                                         	        refphi_sorted.at(jet) << std::endl;
      //}
      //std::cout<<std::endl;
      //for(int jet=0; jet<nref_I; ++jet){
      //	std::cout<<"original order: refpt/eta/phi["<<jet<<"]="<< refpt_F[jet]  << " / " << 
      //	                                                                    refeta_F[jet] << " / " << 
      //                                         	                            refphi_F[jet] << std::endl;
      //}
      //if(nEvt>100)assert(false);
      //END DEBUG
 

      //if(nEvt>0)assert(false);
      
      //std::cout<<std::endl;
      //std::cout<<"-------------------------------------"<<std::endl;
      //std::cout<<"redoing jet matching for nEvt="<<nEvt<<std::endl;
      //std::cout<<"-------------------------------------"<<std::endl;
      

      
      //third, create a boolean array for the gen AND reco jets, each entry initialized to false (not matched yet)
      refjetIsMatched.assign(Nsorted, false);
      refjetIsMatchedTo.assign(Nsorted, -1);
      

      recjetIsMatched.assign(nref_I, false);
      recjetIsMatchedTo.assign(nref_I, -1);
      
      for(int refjet=0; refjet<Nsorted; ++refjet){	//GEN jet loop
	float mindeltaR=10;
	int mindeltaR_ind=-1;
	//std::cout<<"trying to match refjet #"<<refjet<<", refpt/refeta/refphi="<<refpt_sorted.at(refjet)<<"/"<<refeta_sorted.at(refjet)<<"/"<<refphi_sorted.at(refjet)<<std::endl;
	for(int jet=0; jet<nref_I;++jet){ //RECO jet loop
	  
	  //std::cout<<"recjetIsMatched.at("<< jet<<")="<<recjetIsMatched.at(jet)<<std::endl;		  
	  if(recjetIsMatched.at(jet))//if already matched, move on
	    continue;

	  float absdeltaeta=fabs(eta_F[jet]-refeta_sorted.at(refjet));	
	  //std::cout<<"absdeltaeta="<<absdeltaeta <<std::endl;
	  if(absdeltaeta>0.4)//if deltaeta bigger than 0.4, won't work, move on
	    continue;

	  float absdeltaphi=fabs(deltaphi(phi_F[jet],refphi_sorted.at(refjet)));
	  //std::cout<<"absdeltaphi="<<absdeltaphi <<std::endl;
	  if(absdeltaphi>0.4){//if deltaphi bigger than 0.4, won't work, move on
	    //std::cout<<"phi_F["<<jet<<"]="<<phi_F[jet]<<" AND refphi_sorted.at("<<refjet<<")"<<refphi_sorted.at(refjet)<<std::endl;
	    continue;
	  }

	  float deltaR=sqrt(pow(absdeltaeta,2.)+pow(absdeltaphi,2.));
	  //std::cout<<"deltaR="<<deltaR <<std::endl;
	  if(deltaR>0.4)//obvious...
	    continue;
	  

	  if(deltaR<mindeltaR){//if it's the smallest deltaR found yet; note that down
	    mindeltaR_ind=jet;//index of the RECO jet matched to the GEN jet
	    mindeltaR=deltaR;
	  }
	  else if(deltaR==mindeltaR){
	    //std::cout<<"can't believe i happened! make sure you account for this case somehow."<<std::endl;
	    assert(false);
	  }

	}//end loop over RECO jets

	
	if(mindeltaR_ind>=0){
	  refjetIsMatched.at(refjet)=true;
	  refjetIsMatchedTo.at(refjet)=mindeltaR_ind;
	  recjetIsMatched.at(mindeltaR_ind)=true;
	  recjetIsMatchedTo.at(mindeltaR_ind)=refjet;
	  //std::cout<<"BEST MATCH IS JET#"<<refjetIsMatchedTo.at(refjet)<<"! recpt/receta/recphi="<<pt_F[refjetIsMatchedTo.at(refjet)]
	  //	   <<"/"<<eta_F[refjetIsMatchedTo.at(refjet)]<<"/"<<phi_F[refjetIsMatchedTo.at(refjet)]<<std::endl;
	}
	else{
	  //std::cout<<"THIS REF JET HAS NO MATCH"<<std::endl;
	}
	
      }//end loop over GEN jets


      if(insertedagenjet){
	//here i am going to check the jet collection + matches
	std::cout<<std::endl<<"ORIGINAL GENRECO MATCHING"<<std::endl;
	for(int jet=0; jet<nref_I;++jet){
	  if(subid_I[jet]!=0)continue;
	  std::cout<<"                        pt/eta/phi["<<jet<<"]="<<pt_F[jet]<<"/"<<eta_F[jet]<<"/"<<phi_F[jet]<<std::endl;
	  std::cout<<"WAS MATCHED TO refpt/refeta/refphi["<<jet<<"]="<<refpt_F[jet]<<"/"<<refeta_F[jet]<<"/"<<refphi_F[jet]<<std::endl;
	}
	
	//std::cout<<std::endl<<"MY GENRECO MATCHING"<<std::endl;
	for(int jet=0; jet<nref_I;++jet){
	  if(!recjetIsMatched.at(jet))continue;
	  int matchind=recjetIsMatchedTo.at(jet);
	  std::cout<<"                                    pt/eta/phi["<<jet<<"]="<<pt_F[jet]<<"/"<<eta_F[jet]<<"/"<<phi_F[jet]<<std::endl;
	  if(matchind>=0)
	    std::cout<<"IS NOW MATCHED TO (sorted) refpt/refeta/refphi["<< matchind<<"]="<< refpt_sorted.at(matchind)<<"/"<< refeta_sorted.at(matchind)<<"/"<< refphi_sorted.at(matchind)<<std::endl;
	  else
	    std::cout<<"IS UNMATCHED."<<std::endl;
	}
      }
    

    }//end redoing the jet matching
      */
      


    std::vector<bool> GENjetIsMatched; std::vector<int> GENjetIsMatchedTo;
    std::vector<bool> RECOjetIsMatched; std::vector<int> RECOjetIsMatchedTo;
    unsigned int nGENjets =GjetPt->size();
    unsigned int nRECOjets= jetPt->size();
    bool printout=false; int printoutcount=0;
    if(REDOjetMatching){
            //third, create a boolean array for the gen AND reco jets, each entry initialized to false (not matched yet)
      GENjetIsMatched.assign(nGENjets, false);
      GENjetIsMatchedTo.assign(nGENjets, -1);      

      RECOjetIsMatched.assign(nRECOjets, false);
      RECOjetIsMatchedTo.assign(nRECOjets, -1);
      
      for(unsigned int GENjet=0; GENjet<nGENjets; ++GENjet){	//GEN jet loop

	float mindeltaR=10;
	int mindeltaR_ind=-1;
	//std::cout<<"trying to match GENjet #"<<GENjet<<", GENpt/GENeta/GENphi="<<GENpt_sorted.at(GENjet)<<"/"<<GENeta_sorted.at(GENjet)<<"/"<<GENphi_sorted.at(GENjet)<<std::endl;
	for(unsigned int jet=0; jet<nRECOjets;++jet){ //RECO jet loop
	  
	  //std::cout<<"recjetIsMatched.at("<< jet<<")="<<recjetIsMatched.at(jet)<<std::endl;		  
	  if(RECOjetIsMatched.at(jet))//if already matched, move on
	    continue;

	  float absdeltaeta=fabs(jetEta->at(jet)-GjetEta->at(GENjet));	
	  //std::cout<<"absdeltaeta="<<absdeltaeta <<std::endl;
	  if(absdeltaeta>deltaR_req4Match)//if deltaeta bigger than 0.4, won't work, move on
	    continue;
	  
	  float absdeltaphi=fabs(deltaphi(jetPhi->at(jet), GjetPhi->at(GENjet)));
	  //std::cout<<"absdeltaphi="<<absdeltaphi <<std::endl;
	  if(absdeltaphi>deltaR_req4Match){//if deltaphi bigger than 0.4, won't work, move on
	    //std::cout<<"phi_F["<<jet<<"]="<<phi_F[jet]<<" AND refphi_sorted.at("<<refjet<<")"<<refphi_sorted.at(refjet)<<std::endl;
	    continue;
	  }

	  float deltaR=sqrt(pow(absdeltaeta,2.)+pow(absdeltaphi,2.));
	  //std::cout<<"deltaR="<<deltaR <<std::endl;
	  if(deltaR>deltaR_req4Match)//obvious...
	    continue;
	  

	  if(deltaR<mindeltaR){//if it's the smallest deltaR found yet; note that down
	    mindeltaR_ind=jet;//index of the RECO jet matched to the GEN jet
	    mindeltaR=deltaR;
	  }
	  else if(deltaR==mindeltaR){
	    //std::cout<<"can't believe i happened! make sure you account for this case somehow."<<std::endl;
	    assert(false);
	  }

	}//end loop over RECO jets

	
	if(mindeltaR_ind>=0){
	  GENjetIsMatched.at(GENjet)=true;
	  GENjetIsMatchedTo.at(GENjet)=mindeltaR_ind;
	  RECOjetIsMatched.at(mindeltaR_ind)=true;
	  RECOjetIsMatchedTo.at(mindeltaR_ind)=GENjet;
	  if(GjetPt->at(GENjet)>15.){
	    if(genmatchindex_I[GENjet]!=mindeltaR_ind)
	      if(mindeltaR_ind<nref_I)//so the difference in the jet matching will actually get printed out
		{       	  
		  if(false)std::cout<<"genmatchindex_I["<<GENjet<<"]="<<genmatchindex_I[GENjet]<<" BUT mindeltaR_ind = "<< mindeltaR_ind<<std::endl;
		  printout=true;
		}
	  }

	  //std::cout<<"BEST MATCH IS JET#"<<refjetIsMatchedTo.at(refjet)<<"! recpt/receta/recphi="<<pt_F[refjetIsMatchedTo.at(refjet)]
	  //	   <<"/"<<eta_F[refjetIsMatchedTo.at(refjet)]<<"/"<<phi_F[refjetIsMatchedTo.at(refjet)]<<std::endl;
	}
	else{
	  //std::cout<<"THIS REF JET HAS NO MATCH"<<std::endl;
	}
	
      }//end loop over GEN jets


      
      
      if(printout&&false){
	
	std::cout<<"nEvt="<<nEvt<<", GENRECO MATCHING COMPARISON PRINTOUT"<<std::endl;
	for(int jet=0; jet<nref_I;++jet){

	  int matchind=RECOjetIsMatchedTo.at(jet);//from new matching
	  bool resultIsDiff=false;
	  if(subid_I[jet]<0&&matchind<0)continue;
	  else if(subid_I[jet]<0 && matchind>=0)resultIsDiff=false;
	  else if(subid_I[jet]>=0 && matchind<0)resultIsDiff=false;//ok this is actually true but i dont want to print right now
	  else{
	    float refpt=refpt_F[jet], refeta=refeta_F[jet];//, refphi=refphi_F[jet];
	    float Gjetpt=GjetPt->at(matchind), Gjeteta=GjetEta->at(matchind);//, GjetPhi=Gjetphi->at(matchind);
	    if((Gjetpt-refpt)<0.00001 && (Gjeteta-refeta)<0.00001)continue;
	    else {
	      if(Gjetpt>43. || refpt >43. || pt_F[jet]>43.)//cause i dont care at SUPER low pT...
		resultIsDiff=true;
	      else continue;
	    }
	  }
	  
	  if(!resultIsDiff)continue;

	  std::cout<<"                        pt/eta/phi["<<jet<<"] = "<<pt_F[jet]<<"/"<<eta_F[jet]<<"/"<<phi_F[jet]<<std::endl;

	    
	  if(subid_I[jet]==0){
	    std::cout<<"   WAS MATCHED TO    refpt/refeta/refphi["<<jet<<"] = "<<refpt_F[jet]<<"/"<<refeta_F[jet]<<"/"<<refphi_F[jet]<<std::endl;
	    if(matchind>=0)
	      std::cout<<"IS NOW MATCHED TO GjetPt/GjetEta/GjetPhi["<< matchind<<"] = "<< GjetPt->at(matchind)<<"/"<< GjetEta->at(matchind)<<"/"<< GjetPhi->at(matchind)<<std::endl;
	    else
	      std::cout<<"IS UNMATCHED."<<std::endl;	    
	  }
	  else {
	    std::cout<<"WAS UNMATCHED."<<std::endl;
	    if(matchind>=0)
	      std::cout<<"IS NOW MATCHED TO GjetPt/GjetEta/GjetPhi["<< matchind<<"] = "<< GjetPt->at(matchind)<<"/"<< GjetEta->at(matchind)<<"/"<< GjetPhi->at(matchind)<<std::endl;
	    else
	      std::cout<<"IS UNMATCHED."<<std::endl;	    
	  }
	  
	  std::cout<<std::endl;
	  printoutcount++;

	}
	std::cout<<std::endl;
	
	
	
	////std::cout<<std::endl<<"MY GENRECO MATCHING"<<std::endl;
	//for(unsigned int jet=0; jet<nRECOjets;++jet){
	//  if(!RECOjetIsMatched.at(jet))continue;
	//  int matchind=RECOjetIsMatchedTo.at(jet);
	//  std::cout<<"                           pt/eta/phi["<<jet<<"]="<<jetPt->at(jet)<<"/"<<jetEta->at(jet)<<"/"<<jetPhi->at(jet)<<std::endl;
	//}
	//if(printoutcount>100)assert(false);
      }
      



    }//end redoing jet matching








    

    for(int jet = 0; jet<nref_I; ++jet){ //main jet loop
      
      //float rawpt=rawpt_F[jet];      	
      float jtpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
      float recpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
      float receta = eta_F[jet];
      //float recphi = phi_F[jet];
      float jtm=m_F[jet];
      float jetIDpt=recpt;//ala HIN jetID, recpt is corrected w/ L2/L3 residuals
      
      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;
      
      float absreceta = fabs(receta);	
      float jttheta=2.*atan(exp(-1.*receta));
      float jtpz=jtpt/tan(jttheta);
      float jtp=sqrt(jtpt*jtpt + jtpz*jtpz);
      float jtE=sqrt(jtp*jtp + jtm*jtm);
      //std::cout<<"jtE="<<jtE<< ", jtenergy="<<jtenergy<<", jtp="<<jtpt<<std::endl;

      float jty=0.5*log((jtE+jtpz)/(jtE-jtpz));//experimentalist version
      float absjty=fabs(jty);
      //float jtEt=jtE*sin(jttheta);
      //float jtmt=jtm*sin(jttheta);	
      
      float genpt;//   = refpt_F[jet];    
      //float genphi;//  = refphi_F[jet];   
      float geneta;//  = refeta_F[jet];   
      float absgeneta;//  = fabs(geneta); 
      float gendrjt;// = refdrjt_F[jet];  


      if(REDOjetMatching){
	if(!RECOjetIsMatched.at(jet)) continue;//redone gen-reco matching
	else{
	  genpt   = GjetPt->at(RECOjetIsMatchedTo.at(jet));   
	  //genphi  = GjetPhi->at(RECOjetIsMatchedTo.at(jet));   
	  geneta  = GjetEta->at(RECOjetIsMatchedTo.at(jet));   
	  absgeneta  = fabs(geneta); 
	  gendrjt    = -1.;//refdrjt_F[jet];  	
	}
      }
      else{
	if( subid_I[jet]!=0 ) continue;//matching gen-reco
	else{
	  genpt   = refpt_F[jet];    
	  //genphi  = refphi_F[jet];   
	  geneta  = refeta_F[jet];   
	  absgeneta  = fabs(geneta); 
	  gendrjt = refdrjt_F[jet];  
	}
      }
      
      // event+jet counting
      h_NJets->Fill(1);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      
      
      
      else if ( !(recpt > jtPtCut)    ) continue;        //low recopt cut          
      else if ( !(recpt < jtPtCut_Hi)    ) continue;     //high recopt cut
      else if ( !(genpt > genJetPtCut) ) continue;       //low genpt cut
      else if ( !(genpt < genJetPtCut_Hi) ) continue;   //high genpt cut
      else if ( absreceta < jtEtaCutLo ) continue; // lower abseta cut 
      else if (!(absreceta < jtEtaCutHi))continue; // higher abseta cut
      else if ( gendrjt > gendrCut ) continue;       //delta-r cut, proxy for gen-reco matching quality
      
      
      // jet/event counts
      h_NJets_kmatCut1->Fill(1);
      if(!hNEvts_withJets_kmatCut1_Filled){
	h_NEvents_withJets_kmatCut1->Fill(1);
	hNEvts_withJets_kmatCut1_Filled=true;      }      
      
      

      
      // 13 TeV JetID criterion
      bool passesJetID=false; //int jtID=0;
      if(fillMCJetIDHists) 	{
	if (!(absreceta > 2.4)) 
	  passesJetID=(bool)jetID_00eta24( jetIDpt, 
					   neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet], muSum_F[jet],
					   numConst,  chMult, useTightJetID);
//	if (!(absreceta > 2.4)) 
//	  passesJetID=(bool)jetID_00eta24( jetIDpt, 
//					   neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet], muSum_F[jet],
//					   numConst,  chMult, true);
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  passesJetID=(bool) jetID_24eta27( jetIDpt,
					    neSum_F[jet],  phSum_F[jet], muSum_F[jet],
					    numConst, useTightJetID);
	else if( !(absreceta>3.0) && absreceta>2.7 )
	  passesJetID=(bool) jetID_27eta30( jetIDpt,
					    neSum_F[jet],  phSum_F[jet], 
					    neuMult, useTightJetID);
	else  
	  passesJetID=(bool)jetID_32eta47( jetIDpt, 
					   phSum_F[jet], useTightJetID);
	if(!passesJetID) continue;
      }
      
      

      ///// this code to check the way i read info from the tupel tree given the info in the jetanalyzer tree
      //float genE_tup   =-999.;
      //float genpt_tup  =-999.;
      //float geneta_tup =-999.;
      //float genphi_tup =-999.;
      ////float recE_tup   =-999.;
      ////float recpt_tup  =-999.;
      //float receta_tup =-999.;
      //float recphi_tup =-999.;
      //if(recomatchindex_I.at(jet)>=0 ){
      //	//genE_tup  =GjetE->at(genmatchidx); 
      //	//genpt_tup =GjetPt->at(genmatchidx); 
      //	geneta_tup=GjetEta->at(recomatchindex_I.at(jet));
      //	genphi_tup=GjetPhi->at(recomatchindex_I.at(jet));
      //	float genDR=sqrt(pow((geneta_tup-geneta),2) + pow((genphi_tup-genphi),2));
      //	float radius_F=((float)radius)/10.;
      //	if(genDR>radius_F){
      //	  if(debugMode){
      //	    std::cout<<"WARNING: "<< radius_F << " < genDR ="<<genDR<<std::endl;
      //	    std::cout<<"(geneta_tup,geneta) = ("<< geneta_tup <<","<< geneta<<")"<<std::endl;
      //	    std::cout<<"(genphi_tup,genphi) = ("<< genphi_tup <<","<< genphi<<")"<<std::endl<<std::endl;
      //	  }
      //	  badgenjetcount++;
      //	}
      //	else goodgenjetcount++;
      //	//recE_tup  =jetE  ->at(jet); 
      //	//recpt_tup =jetPt ->at(jet); 
      //	receta_tup=jetEta->at(jet);
      //	recphi_tup=jetPhi->at(jet);
      //	float recDR=sqrt(pow((receta_tup-receta),2) + pow((recphi_tup-recphi),2));
      //	if(recDR>radius_F)
      //	  if(debugMode){
      //	    std::cout<<"WARNING: "<< radius_F << " < recDR ="<<recDR<<std::endl;
      //	  }
      //}
      
      //this only ever happens because the gen jets in the RECO tree have a low pT cut of 10 GeV but in the GEN tree have a low pT cut of 15 GeV.  
      //so sometimes, the matched GEN jet can be found in the RECO tree, but not in the GEN tree.

      if(refpt_F[jet]>15.)
	if(recomatchindex_I.at(jet)<0){
	  
	  std::cerr<<std::endl;
	  std::cerr<<"//////////////////////////////////////////////"<<std::endl;
	  std::cerr<<"////////////// !!! WARNING !!! ///////////////"<<std::endl;
	  std::cerr<<"reco jet w/ nonnegative subid has a negative recomatchindex!!!"<<std::endl;
	  std::cerr<<"---------DETAILS--------"<<std::endl;
	  std::cerr<<std::endl<<"_______ nEvt="<<nEvt<<" ________"<<std::endl;
	  std::cerr<<"ngen_I="<<ngen_I<<std::endl;
	  for(int gjet=0; gjet<ngen_I; ++gjet){
	    std::cerr<<"// genmatchindex_I["<<gjet<<"]="<<genmatchindex_I[gjet]<<std::endl;
	    std::cerr<<"genpt_F["<<gjet<<"]="<<genpt_F[gjet]<<std::endl;
	  }
	  std::cerr<<std::endl;
	  std::cerr<<"nref_I="<<nref_I<<std::endl;
	  for(int kjet=0;kjet<nref_I;++kjet){
	    std::cerr<<"// recomatchindex_I["<<kjet<<"]="<<recomatchindex_I.at(kjet)<<std::endl;
	    std::cerr<<"subid_I["<<kjet<<"]="<<subid_I[kjet]<<std::endl;
	    std::cerr<<"refpt["<<kjet<<"]="<<refpt_F[kjet]<<std::endl;
	    if(recomatchindex_I[kjet]!=-1){
	      std::cerr<<"genpt["<<recomatchindex_I[kjet]<<"]="<<genpt_F[recomatchindex_I[kjet]]<<std::endl;
	      std::cerr<<"GJetAk04JetPt["<<recomatchindex_I[kjet]<<"]="<<GjetPt->at(recomatchindex_I[kjet])<<std::endl;
	    }
	  }
	  badgenjetcount++;
	  continue;
	}

      float genjttheta=2.*atan(exp(-1.*geneta));
      float genjtpz=genpt/tan(genjttheta);
      float genjtE;
      if(REDOjetMatching)
	genjtE=GjetE->at(RECOjetIsMatchedTo.at(jet));
      else 
	genjtE=GjetE->at(recomatchindex_I.at(jet));
      float genjty=0.5*log((genjtE+genjtpz)/(genjtE-genjtpz));//experimentalist version
      float absgenjty=fabs(genjty);
      

      h_NJets_kmatCut2->Fill(1);
      if(!hNEvts_withJets_kmatCut2_Filled){
	h_NEvents_withJets_kmatCut2->Fill(1);
	hNEvts_withJets_kmatCut2_Filled=true;      }

      
      int theYBin=-1;
      for(int ybin=0;ybin<nbins_absy;ybin++)
	if( absybins[ybin]<=absjty  && 		
	    absjty<absybins[ybin+1]    	      ) {	    
	  theYBin=ybin;
	  break;	  }       	
      
      int theGENYBin=-1;
      for(int ybin=0;ybin<nbins_absy;ybin++)
	if( absybins[ybin]<=absgenjty  && 		
	    absgenjty<absybins[ybin+1]    	      ) {	    
	  theGENYBin=ybin;
	  break;	  }       	
      
      //int theEtaBin=-1;
      //for(int etabin=0;etabin<nbins_abseta;etabin++)
      //if( absetabins[etabin]<=absreceta  && 		
      //    absreceta<absetabins[etabin+1]    	      ) {	    
      //  theEtaBin=etabin;
      //  break;	  }       	
      
      int theGENEtaBin=-1;
      for(int etabin=0;etabin<nbins_abseta;etabin++)
	if( absetabins[etabin]<=absgeneta  && 		
	    absgeneta<absetabins[etabin+1]    	      ) {	    
	  theGENEtaBin=etabin;
	  break;	  }       	
      
      
      
      
      
      /////   JERS HISTS   ///// 
      if(fillJERSHists ){
	
	//// dual-diff xsec abseta bins for two spec genpt bins and all genpt bins
	//int absetabin=-1;
	//if(absreceta>absetabins[nbins_abseta]) absetabin = -1;//check end of array
	//else{
	//  for(int bin=0; bin < nbins_abseta; bin++)
	//    if( absreceta>=absetabins[bin] ) absetabin = bin; 				
	//}
	//
	//if(absetabin != -1 ){
	//  
	//  for(int jtID=0; jtID<2;jtID++){
	//    if(jtID==1 && !passesJetID) continue;
	//    else if (jtID==1 && !fillMCJetIDHists)continue;
	//    
	//    hJER_absetabins[absetabin]->Fill( (float)(recpt/genpt), weight_eS);     
	//    
	//    if( recpt >= ( (float) genptBin1Low  )  && 
	//	recpt <  ( (float) genptBin1High )    ) 
	//      hJER_absetabins_genptBin1[jtID][absetabin]->Fill( (float)(recpt/genpt), weight_eS);     
	//    
	//    if( recpt >= ( (float) genptBin2Low  )  && 
	//	recpt <  ( (float) genptBin2High )    ) 
	//      hJER_absetabins_genptBin2[jtID][absetabin]->Fill( (float)(recpt/genpt), weight_eS);     	
	//  }
	//}
	
	//// hJER, genpt30-50, 150-200, across receta bins array (note, NOT abs-value receta bins)
	//int etabin = -1;
	//if(fabs(receta) >= etabins[nbins_eta]) etabin=-1; //check end of array
	//else {
	//  for(int bin = 0; bin<nbins_eta; ++bin)  //notice that the loop does not hit end of this array
	//    if(receta >= etabins[bin]) etabin = bin;	  }
	//
	//if(etabin != -1) {
	//  
	//  for(int jtID=0; jtID<2;jtID++){
	//    if(jtID==1 && !passesJetID)continue;
	//    else if (jtID==1 && !fillMCJetIDHists)continue;
	//    
	//    if(recpt >= 30.  && recpt<50. ) 
	//      hJER_eta_30pt50[jtID][etabin]->Fill((float)recpt/genpt, weight_eS);
	//    else if(recpt >= 150. && recpt<200.) 
	//      hJER_eta_150pt200[jtID][etabin]->Fill((float)recpt/genpt, weight_eS);	
	//  }
	//}
	
	//JER in genptbins array, according to input absetacuts
	//	int ptbin = -1; 
	//	if( !(genpt < ptbins[nbins_pt]) ) ptbin=-1;//check end of binning array    //else if (!passesEtaCut) ptbin=-1;
	//	else {
	//	  for(int bin = 0; bin<nbins_pt; ++bin)      
	//	    if(! (genpt < ptbins[bin])) ptbin = bin;  	
	//	}
	
	
	int ptbin=-1;
	for(int bin=0;bin<nbins_pt;bin++)
	  if( ptbins[bin]<=genpt  && 		
	      genpt<ptbins[bin+1]    	      ) {	    
	    ptbin=bin;
	    if(debugMode_verbose)std::cout<<"genpt="<<genpt<<" should be between"<<std::endl;
	    if(debugMode_verbose)std::cout<<"ptbins["<<bin<<"]="<<ptbins[ptbin]<<" and ptbins["<<bin<<"]="<<ptbins[bin+1]<<std::endl;	    
	    break;	  }       	
	
	//	int recptbin=-1;
	//	for(int bin=0;bin<nbins_pt;bin++)
	//	  if( ptbins[bin]<=recpt  && 		
	//	      recpt<ptbins[bin+1]    	      ) {	    
	//	    recptbin=bin;
	//	    if(debugMode_verbose)std::cout<<"recpt="<<recpt<<" should be between"<<std::endl;
	//	    if(debugMode_verbose)std::cout<<"ptbins["<<bin<<"]="<<ptbins[ptbin]<<" and ptbins["<<bin<<"]="<<ptbins[bin+1]<<std::endl;	    
	//	    break;	  }       	
	
	if( (ptbin!=-1) ){	
	  hJER[ptbin]->Fill( (float)(recpt/genpt), weight_eS); 	  	  
	  if(theYBin!=-1)hJER_jtptQA[theYBin]->Fill(recpt,weight_eS);      
	  if(theGENYBin!=-1)hJER_genptQA[theGENYBin]->Fill(genpt,weight_eS);      
	  hJER_genetaQA->Fill(geneta,weight_eS);      
	  hJER_jtetaQA->Fill(receta,weight_eS);      
	  hJER_jtyQA->Fill(jty,weight_eS);      
	  //if( (theEtaBin!=-1) ) hJER_abseta[ptbin][theEtaBin]->Fill( (float)(recpt/genpt), weight_eS);         // genpt bin, receta bin
	  if( (theGENEtaBin!=-1) ) hJER_absgeneta[ptbin][theGENEtaBin]->Fill( (float)(recpt/genpt), weight_eS);// genpt bin, geneta bin
	  //if(   (theYBin!=-1) ) hJER_absy[ptbin][theYBin]->Fill(     (float)(recpt/genpt), weight_eS);         // genpt bin, recy bin
	  if( (theGENYBin!=-1) ) hJER_absgeny[ptbin][theGENYBin]->Fill((float)(recpt/genpt), weight_eS);       // genpt bin, geny bin
	}
	else continue;
	
      }//fill JERS hists      
    }//end jet loop
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
  std::cout<<"Total Num of good gen jets = "<<goodgenjetcount<<std::endl;
  std::cout<<"Total Num of bad gen jets  = "<<badgenjetcount<<std::endl;
  std::cout<<std::endl;
  

  std::cout<<"writing output file "<<outfile<<std::endl;
  fout->Write(); 


  std::cout<<std::endl<<"readForests_ppMC_JERS finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  
  return 0;
  }

// end readForests_ppMC_JERS
								 

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













  //TH1D * hJER_absgeneta[nbins_pt][nbins_abseta]={};       // rec/gen , bins of gen pt
  //  TH1D * hJES[nbins_pt]={};       // rec/raw , bins of raw pt
  //  TH1D * hJER_raw[nbins_pt]={};   // raw/gen , bins of gen pt
  
  ////two specific genpt ranges across eta bins
  //TH1D * hJER_eta_30pt50[nbins_eta]={};
  //TH1D * hJER_eta_150pt200[nbins_eta]={};
  
  ////all genpt, and two specifically set genpt bins, across abs eta
  //TH1D *hJER_absetabins[nbins_abseta];
  //
  //TH1D *hJER_absetabins_genptBin1[nbins_abseta]; 
  //TH1D *hGenPtCutHi1   = new TH1D("hGenPtCutHi1"  , (std::to_string(genptBin1High)).c_str() ,  100,  0., 100.  ); 
  //TH1D *hGenPtCutLo1   = new TH1D("hGenPtCutLo1"  , (std::to_string(genptBin1Low)).c_str()  ,  100,  0., 100.  ); 	
  //hGenPtCutHi1->Fill((float)genptBin1High);  hGenPtCutLo1->Fill((float)genptBin1Low); 
  //
  //TH1D *hJER_absetabins_genptBin2[nbins_abseta]; 
  //TH1D *hGenPtCutHi2   = new TH1D("hGenPtCutHi2"  , (std::to_string(genptBin2High)).c_str() ,  100,  0., 100.  ); 
  //TH1D *hGenPtCutLo2   = new TH1D("hGenPtCutLo2"  , (std::to_string(genptBin2Low)).c_str()  ,  100,  0., 100.  ); 
  //hGenPtCutHi2->Fill((float)genptBin2High);  hGenPtCutLo2->Fill((float)genptBin2Low); 














    //int jtID=0;
    //if(fillMCJetIDHists)jtID=1;
    //hJEC[0]= new TH3F( ("hJEC_"+std::to_string(0)+"wJetID").c_str(), ";raw p_{T};#eta;JEC", 500,0.,500., 200,-5.,+5., 300,0.,5.);    
    
    //for(int y = 0; y<nbins_eta; ++y){	               
    //
    //	std::string hTitleJER_30to50="hJER_"+std::to_string(jtID)+"wJetID_"+
    //	  "etabin"+std::to_string(y)+"_30_pt_50";//"etabin"+std::to_string(etabins[y])+"_30_pt_50";
    //	hJER_eta_30pt50[y] = new TH1D(hTitleJER_30to50.c_str(), "recpt/genpt, 30<genpt<50",100, 0, 2.);      
    //	std::string hTitleJER_150to200="hJER_"+std::to_string(jtID)+"wJetID_"+
    //	  "etabin"+std::to_string(y)+"_150_pt_200";
    //	hJER_eta_150pt200[y] = new TH1D(hTitleJER_150to200.c_str(), "recpt/genpt, 150<genpt<200",100, 0, 2.);          
    //	
    //	/*for(int x = 0; x<nbins_JEC_ptbins; ++x) {
    //	  std::string hTitleJEC_check="hJEC_"+std::to_string(jtID)+"wJetID_check_"+
    //	  "ptbin"+std::to_string(x)+"_etabin"+std::to_string(y);
    //	  hJEC_check[x][y] = new TH1D(hTitleJEC_check.c_str(),
    //	  Form("rawpt/genpt, %2.0f<genpt<%2.0f, %2.4f<geneta<%2.4f", 
    //	  JEC_ptbins[x],JEC_ptbins[x+1], etabins[y],etabins[y+1]),
    //	  100,0.,3.);           }*/ 
    //}    









    //for(int j=0;j<nbins_abseta;++j){
    //	for(int bin = 0; bin<nbins_pt; ++bin) {      	  
    //	  std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );	std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
    //	  std::string absetabins_i=absetabins_str[j];	std::string absetabins_ip1=absetabins_str[j+1];
    //	  std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID"+"_genetabin"+std::to_string(j)+"_"+ptbins_i+"ptbin"+ptbins_ip1;
    //	  if(debugMode&&verbose)std::cout<<"hTitleJER="<<hTitleJER<<std::endl;
    //	  std::string hDescJER="rec/gen pt, "+absetabins_i+"<|GEN y|<"+absetabins_ip1+", "+ptbins_i+"<genpt<"+ptbins_ip1;	  
    //	  if(debugMode&&verbose)std::cout<<"hDescJER="<<hDescJER<<std::endl;
      //	  hJER_absgeneta[bin][j] = new TH1D(hTitleJER.c_str(),hDescJER.c_str(), 200,0.,2.); 
      //	}//ptbin loop
      //}//etabin loop
      //assert(false);
      //      //TH1D * hJES[nbins_pt]={};       // rec/raw , bins of raw pt
      //      for(int bin = 0; bin<nbins_pt; ++bin) {      
      //	std::string hTitleJES="hJES_"+std::to_string(jetIDint)+"wJetID_ptbin"+std::to_string(bin);
      //	std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );
      //	std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
      //	std::string hDescJES="rec/raw pt, "+ptbins_i+"<genpt<"+ptbins_ip1;
      //	std::cout<<"hDescJES="<<hDescJES<<std::endl;
      //	hJES[bin] = new TH1D(hTitleJES.c_str(),hDescJES.c_str(), 200,0.,2.); 
      //       }
      //
      //      //TH1D * hJER_raw[nbins_pt]={};   // raw/gen , bins of gen pt
      //      for(int bin = 0; bin<nbins_pt; ++bin) {      
      //	std::string hTitleJER_raw="hJER_raw_"+std::to_string(jetIDint)+"wJetID_ptbin"+std::to_string(bin);
      //	std::string ptbins_i=std::to_string( ( (int)ptbins[bin]) );
      //	std::string ptbins_ip1=std::to_string( ( (int)ptbins[bin+1]) );
      //	std::string hDescJER_raw="raw/gen pt, "+ptbins_i+"<genpt<"+ptbins_ip1;
      //	std::cout<<"hDescJER_raw="<<hDescJER_raw<<std::endl;
      //	hJER_raw[bin] = new TH1D(hTitleJER_raw.c_str(),hDescJER_raw.c_str(), 200,0.,2.); 
      //      } 
      
      //for(int bin = 0; bin<nbins_abseta; ++bin) {      
      //	std::string hTitleJER="hJER_"+std::to_string(jetIDint)+"wJetID_jty_bin"+std::to_string(bin);      
      //	hJER_absetabins[bin] = new TH1D(hTitleJER.c_str(),("recpt/genpt, absetaBin"+std::to_string(bin)).c_str(), 100,0.,2.); 
      //	
      //	std::string hTitleJER_1="hJER_"+std::to_string(jetIDint)+"wJetID_jty_etabin"+std::to_string(bin)+
      //	  "_genptBin1";
      //	hJER_absetabins_genptBin1[bin] = new TH1D(hTitleJER_1.c_str(),("absetaBin "+std::to_string(bin)+", "						   
      //									     "genpt "+std::to_string(genptBin1Low)+" to "+std::to_string(genptBin1High)).c_str(), 100,0.,2.); 
      //	
      //	std::string hTitleJER_2="hJER_"+std::to_string(jetIDint)+"wJetID_jty_etabin"+std::to_string(bin)+
      //	  "_genptBin2";
      //	hJER_absetabins_genptBin2[bin] = new TH1D(hTitleJER_2.c_str(),("absetaBin "+std::to_string(bin)+", "						   
      //									     "genpt "+std::to_string(genptBin2Low)+" to "+std::to_string(genptBin2High)).c_str(), 100,0.,2.); 
      //}
