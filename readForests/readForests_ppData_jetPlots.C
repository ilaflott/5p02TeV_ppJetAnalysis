// custom header
#include "readForests_jetPlots.h"

// ppData switches
const bool fillDataEvtQAHists=true;
//const bool fillDataVtxTrkQAHists=true; //in the works

const bool fillDataJetQAHists=true;
const bool fillDataJetIDHists=false;//, tightJetID=false;
const bool fillDataJetTrigQAHists=true; //data-specific

const bool fillDataDijetHists=false;
const bool fillDataJetSpectraRapHists=false; //other

const bool useHLT100=false;

const int jetIDint=(int)fillDataJetIDHists;

const std::string trgCombType="Calo";
//const std::string trgCombType="PF";


//// readForests_ppData_jetPlots
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData_jetPlots( std::string inFilelist , int startfile , int endfile , 
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
  std::cout<<"///////////////////"<<std::endl<<std::endl;

  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  

  //check out what kind of filelist we're using
  bool filelistIsLowerJets=false; //to know if i should apply the duplicate check or not
  bool filelistIsJet80=false;  
  if( inFilelist.find("LowerJets")!=std::string::npos ) {
    filelistIsLowerJets=true; 
    std::cout<<"running LowerJets filelist!"<<std::endl; }
  else if( inFilelist.find("Jet80")    !=std::string::npos ) {
    filelistIsJet80=true;
    std::cout<<"running Jet80 filelist!"<<std::endl; }
  else assert(false);
  assert( !(filelistIsJet80&&filelistIsLowerJets) );
  
  // jet tree name+directory formation
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;

  // initialize tree name array
  std::string trees[N_dataTrees];
  trees[0]=jetTreeName;
  for(int i=1;i<N_dataTrees;++i){
    if(i<5)trees[i]=dataTreeNames[i];
    else{
      if(trgCombType=="PF")
	trees[i]=dataTreeNames[i]+PF_HLTBitStrings[i-5]+"_v";
      else if(trgCombType=="Calo")
	trees[i]=dataTreeNames[i]+Calo_HLTBitStrings[i-5]+"_v";
      else assert(false);
    }
  }
  if(debugMode)for(int i=0;i<N_dataTrees;++i)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  

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
    for(int t = 0;t<N_dataTrees;++t) filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
    lastFileAdded=filename_Forest;  }//end file loop
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
  
  //evt counts
  TH1D *h_NEvents          = new TH1D("NEvents"         , "NEvents",         2,0,2); //0-1 bin --> fill w/ unweighted entry, 1-2 bin --> fill w/ weighted entry
  TH1D *h_NEvents_trigd_1    = new TH1D("NEvents_trigd_1"   , "NEvents trigd",   2,0,2);  
  TH1D *h_NEvents_nopscl   = new TH1D("NEvents_nopscl"  , "NEvents triggered w/ pscl==1",   2,0,2);
  TH1D *h_NEvents_notTrigd    = new TH1D("NEvents_notTrigd"   , "NEvents not Trigd",   2,0,2);  
  TH1D *h_NEvents_jet40    = new TH1D("NEvents_jet40Trigd"   , "NEvents jet40Trigd",    2,0,2);
  TH1D *h_NEvents_jet60    = new TH1D("NEvents_jet60Trigd"   , "NEvents jet60Trigd",    2,0,2);
  TH1D *h_NEvents_jet80    = new TH1D("NEvents_jet80Trigd"   , "NEvents jet80Trigd",    2,0,2);
  TH1D *h_NEvents_jet100   = new TH1D("NEvents_jet100Trigd"  , "NEvents jet100Trigd",   2,0,2);
  
  TH1D *h_NEvents_skipped  = new TH1D("NEvents_skipped" , "NEvents skipped", 2,0,2);
  TH1D *h_NEvents_read     = new TH1D("NEvents_read"    , "NEvents read",    2,0,2);
  TH1D *h_NEvents_trigd_2    = new TH1D("NEvents_trigd_2"   , "NEvents trigd after duplicate skip",   2,0,2);  
  TH1D *h_NEvents_skimCut  = new TH1D("NEvents_skimCut" , "NEvents read post skimCut",  2,0,2);
  TH1D *h_NEvents_vzCut    = new TH1D("NEvents_vzCut"   , "NEvents read post vzCut AND skimCut", 2,0,2);
  
  //evt counts post all evt cuts. thrown out because no jets left in the event pre or post cuts
  TH1D *h_NEvents_withJets           = new TH1D("NEvents_withJets" , 
		  				"NEvents post evt cuts/trigger, w/ nref_I>0", 2,0,2);
  TH1D *h_NEvents_withJets_jtptCut   = new TH1D("NEvents_withJets_jtptCut" , 
						"NEvents w/ jets post jtptCut", 2,0,2);
  TH1D *h_NEvents_withJets_jtetaCut1   = new TH1D("NEvents_withJets_jtetaCut1" , 
						  "NEvents w/ jets post |jteta|<4.7", 2,0,2);
  TH1D *h_NEvents_withJets_jtetaCut2   = new TH1D("NEvents_withJets_jtetaCut2" , 
						  "NEvents w/ jets post jteta Hi/Lo cuts", 2,0,2);
  TH1D *h_NEvents_withJets_JetIDCut = NULL;
  if(fillDataJetIDHists)
    h_NEvents_withJets_JetIDCut= new TH1D("NEvents_withJets_JetIDCut" , 
					  "NEvents read post evt cuts, w/ jets post kmatCut AND JetID Cut", 2,0,2);
  
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

  TH1D *h_NJets          = new TH1D("NJets","NJets read", 2,0,2);
  TH1D *h_NJets_jtptCut  = new TH1D("NJets_jtptCut ","NJets read post jtptCut ", 2,0,2);
  TH1D *h_NJets_jtetaCut1  = new TH1D("NJets_jtetaCut1 ","NJets read post jtetaCut1 ", 2,0,2);
  TH1D *h_NJets_jtetaCut2  = new TH1D("NJets_jtetaCut2 ","NJets read post jtetaCut2 ", 2,0,2);
  TH1D *h_NJets_JetIDCut=NULL;
  if(fillDataJetIDHists)
    h_NJets_JetIDCut = new TH1D("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 2,0,2);
  
  //jet QA, i.e. hists filled once per jet
  TH1D *hJetQA[2][N_vars]={};
  TH1D *hJetQA_jtptEntries=NULL;
  TH1D *hJetQA_jtpt=NULL, *hJetQA_jtpt_L2Res=NULL, *hJetQA_jtpt_L3Res=NULL;
  TH2D *hJetQA_jtpt_v_L2Res=NULL,*hJetQA_jtpt_v_L3Res=NULL,*hJetQA_jtpt_L2Res_v_L3Res=NULL;
  if(fillDataJetQAHists)
    for(int k = 0; k<2; ++k){
      if(!fillDataJetIDHists && k==1)continue;	
      if(fillDataJetIDHists && k==0)continue;	
      
      hJetQA_jtptEntries = new TH1D( Form("hJetQA_%dwJetID_jtptEntries", k) , Form("jtpt hist, weight=1 always") , 2500,0,2500);       

      hJetQA_jtpt = new TH1D( Form("hJetQA_%dwJetID_jtpt_forRes", k) , Form(";%s;", "jtpt_forRes") , 2500,0,2500);       
      hJetQA_jtpt_L2Res = new TH1D( Form("hJetQA_%dwJetID_jtpt_L2Res", k) , Form(";%s;", "jtpt_L2Res") , 2500,0,2500);       
      hJetQA_jtpt_L3Res = new TH1D( Form("hJetQA_%dwJetID_jtpt_L3Res", k) , Form(";%s;", "jtpt_L3Res") , 2500,0,2500);       

      hJetQA_jtpt_v_L2Res = new TH2D( Form("hJetQA_%dwJetID_jtpt_v_L2Res", k) , Form("jtpt v. jetpt_L2Res" ) , 2500,0,2500 , 2500,0,2500);
      hJetQA_jtpt_v_L3Res = new TH2D( Form("hJetQA_%dwJetID_jtpt_v_L3Res", k) , Form("jtpt v. jetpt_L3Res" ) , 2500,0,2500 , 2500,0,2500);
      hJetQA_jtpt_L2Res_v_L3Res = new TH2D( Form("hJetQA_%dwJetID_jtpt_L2Res_v_L3Res", k) , Form("jtpt_L2Res v. jetpt_L3Res" ) , 2500,0,2500 , 2500,0,2500);
      
      for(int j = 0; j<N_vars; ++j){	
	
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
	else if (fillDataDijetHists){
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
      } }
  
  TH1D *hJetSpectraRap[2][nbins_abseta]={};
  if(fillDataJetSpectraRapHists)
    for(int k = 0; k<2; ++k) {
      if(!fillDataJetIDHists && k==1)continue;	
      if(fillDataJetIDHists && k==0)continue;	
      for(int j = 0; j<nbins_abseta; ++j){
	std::string h_Title="hJetSpectraRap_";
	if(k==1)h_Title+="wJetID_";
	h_Title+="bin"+std::to_string(j);
	std::string h_Desc;
	std::stringstream stream1, stream2;	
	stream1.precision(1); stream1 << std::fixed << absetabins[j];
	stream2.precision(1); stream2 << std::fixed << absetabins[j+1];
	h_Desc="JetPt Spectra for "+stream1.str()+"<abs(y)<"+ stream2.str();	
	hJetSpectraRap[k][j]=new TH1D(h_Title.c_str(),h_Desc.c_str(), 2000,0,2000);  
      } }

  // trigger (also jet) level
  TH1D  *hpp_TrgObj40[2]={}, *hpp_TrgObj60[2]={}, *hpp_TrgObj80[2]={}, *hpp_TrgObj100[2]={};// plots for incl spectra for each individual trigger
  TH1D  *hpp_ExcTrgObj40[2]={}, *hpp_ExcTrgObj60[2]={}, *hpp_ExcTrgObj80[2]={}, *hpp_ExcTrgObj100[2]={};// plots for trigger-exclusive inclusive jet spectra
  TH1D *hpp_TrgObjComb[2]={};//combination of triggers w/ exclusion
  
  TH1D  *hpp_IncHLT40trgPt =NULL ,   *hpp_IncHLT40trgEta =NULL;
  TH1D  *hpp_IncHLT60trgPt =NULL ,   *hpp_IncHLT60trgEta =NULL;
  TH1D  *hpp_IncHLT80trgPt =NULL ,   *hpp_IncHLT80trgEta =NULL;
  TH1D *hpp_IncHLT100trgPt =NULL ,  *hpp_IncHLT100trgEta =NULL;


  TH1D  *hpp_HLT40trgPt =NULL ,   *hpp_HLT40trgEta =NULL;
  TH1D  *hpp_HLT60trgPt =NULL ,   *hpp_HLT60trgEta =NULL;
  TH1D  *hpp_HLT80trgPt =NULL ,   *hpp_HLT80trgEta =NULL;
  TH1D *hpp_HLT100trgPt =NULL ,  *hpp_HLT100trgEta =NULL;
  
  TH1D *hpp_HLTCombtrgPt=NULL, *hpp_HLTCombtrgEta=NULL;
  
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

    std::string hJetTrigQATitleArray[]={ "hpp_HLT40", "hpp_HLT60",  "hpp_HLT80", "hpp_HLT100", "hpp_HLTComb"};
    std::string hJetTrigQADescArray[]={ "HLT40Jet","HLT60Jet",  "HLT80Jet", "HLT100Jet", "HLTCombo"};
    const int nJetTrigQATitles=sizeof(hJetTrigQATitleArray)/sizeof(std::string);
    
    std::string hJetExcTrigQATitleArray[]={ "hpp_excHLT40", "hpp_excHLT60",  "hpp_excHLT80", "hpp_excHLT100"};
    std::string hJetExcTrigQADescArray[]={ "ExcHLT40Jet","ExcHLT60Jet",  "ExcHLT80Jet", "ExcHLT100Jet"};
    const int nJetExcTrigQATitles=sizeof(hJetExcTrigQATitleArray)/sizeof(std::string);
    
    for(int j=0;j<2;j++){      
      
      if(!fillDataJetIDHists && j==1) continue;    
      if(fillDataJetIDHists && j==0)continue;	

      if(debugMode) std::cout<<"j=="<<j<<std::endl;
      
      for(int k=0; k<nJetTrigQATitles; k++){
	
	if(debugMode) std::cout<<"k=="<<k<<std::endl;	
	
	std::string hTitle=hJetTrigQATitleArray[k];
	std::string hDesc=hJetTrigQADescArray[k];
	if(j==1)hTitle+="_wJetID";
	if(j==1)hDesc+=" w/ JetID";
	hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;
	
	if(debugMode)std::cout<<"inclusive HLT plot hTitle="<<hTitle<<std::endl;
	if(debugMode)std::cout<<"inclusive HLT plot hDesc="<<hDesc<<std::endl;
	
	if     (hJetTrigQATitleArray[k]=="hpp_HLT40"  )  hpp_TrgObj40[j]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	else if(hJetTrigQATitleArray[k]=="hpp_HLT60"  )  hpp_TrgObj60[j]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	else if(hJetTrigQATitleArray[k]=="hpp_HLT80"  )  hpp_TrgObj80[j]     = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	else if(hJetTrigQATitleArray[k]=="hpp_HLT100" )  hpp_TrgObj100[j]    = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);	
	else if(hJetTrigQATitleArray[k]=="hpp_HLTComb")  hpp_TrgObjComb[j]   = new TH1D( hTitle.c_str(), (hDesc).c_str(), 2500,0,2500);
	
	if( k<nJetExcTrigQATitles ){
	  
	  std::string hExcTitle=hJetExcTrigQATitleArray[k];
	  std::string hExcDesc=hJetExcTrigQADescArray[k];
	  if(j==1)hExcTitle+="_wJetID";
	  if(j==1)hExcDesc+=" w/ JetID";
	  hExcTitle+="_R"+std::to_string(radius)+"_"+etaWidth;
	  if(debugMode)std::cout<<" exclusive HLT plothExcTitle="<<hExcTitle<<std::endl;
	  if(debugMode)std::cout<<" exclusive HLT plothExcDesc="<<hExcDesc<<std::endl<<std::endl;
	  
	  if     (hJetExcTrigQATitleArray[k]=="hpp_excHLT40"  )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj40[j]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	  else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT60"  )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj60[j]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	  else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT80"  )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj80[j]     = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}
	  else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT100" )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj100[j]    = new TH1D( hExcTitle.c_str(), (hExcDesc).c_str(), 2500,0,2500);}	  
	}
	
      } } }

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
    if(debugMode){std::cout<<std::endl;
      std::cout<<"PF_HLTBranches        [i="<<i<<"]= "<<PF_HLTBranches[i]<<std::endl;
      std::cout<<"PF_HLTPresclBranches  [i="<<i<<"]= "<<PF_HLTPresclBranches[i]<<std::endl;
      std::cout<<"Calo_HLTBranches      [i="<<i<<"]= "<<Calo_HLTBranches[i]<<std::endl;
      std::cout<<"Calo_HLTPresclBranches[i="<<i<<"]= "<<Calo_HLTPresclBranches[i]<<std::endl<<std::endl; }    }
  
  int Jet40_I, Jet60_I, Jet80_I, Jet100_I;   //primary HLT bits, either Calo or PF
  int Jet40_p_I, Jet60_p_I, Jet80_p_I, Jet100_p_I;    //prescales
  int Jet40_2_I, Jet60_2_I; //the other 40/60 bits from whichever one wasn't picked (for duplicate skip)
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
    jetpp[4]->SetBranchAddress( Calo_HLTBranches[1].c_str() , &Jet60_2_I);  }
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
    jetpp[4]->SetBranchAddress( PF_HLTBranches[1].c_str() , &Jet60_2_I);}
  else assert(false);
  
  //ONE HLT path ONE tree ONE trig obj pt branch
  std::vector<double>  *trgObjpt_40, *trgObjpt_60, *trgObjpt_80, *trgObjpt_100;
  jetpp[5]->SetBranchAddress("pt",&trgObjpt_40);
  jetpp[6]->SetBranchAddress("pt",&trgObjpt_60);  
  jetpp[7]->SetBranchAddress("pt",&trgObjpt_80);  
  jetpp[8]->SetBranchAddress("pt",&trgObjpt_100);
  
  std::vector<double>  *trgObjeta_40, *trgObjeta_60, *trgObjeta_80, *trgObjeta_100;
  jetpp[5]->SetBranchAddress("eta",&trgObjeta_40);
  jetpp[6]->SetBranchAddress("eta",&trgObjeta_60);  
  jetpp[7]->SetBranchAddress("eta",&trgObjeta_80);  
  jetpp[8]->SetBranchAddress("eta",&trgObjeta_100);
  
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
  std::cout<<"reading "<<NEvents_toRead<<" events"<<std::endl;  
  
  UInt_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;  //trigger event counts

  //JetEnergy Residuals
  int etacutForResid=3;
  if(radius==3)etacutForResid=4;
  
  L2ResidualJES* L2JES =NULL;//= new L2ResidualJES(radius, etacutForResid, "pp5");
  L3ResidualJES* L3JES =NULL;//= new L3ResidualJES("pp5");
  bool doResidualCorr=false; 
  if(radius==3 || radius==4){
    doResidualCorr=true;
    L2JES = new L2ResidualJES(radius, etacutForResid, "pp5");
    L3JES = new L3ResidualJES("pp5");  }
    
  //float jetIDCut_neSum, jetIDCut_phSum;
  //if(tightJetID){     jetIDCut_neSum=0.90;  jetIDCut_phSum=0.90;}
  //else{     jetIDCut_neSum=0.99;  jetIDCut_phSum=0.99;}
      
  //begin event loop
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {
    
    //grab an entry
    if( debugMode && (nEvt%10000==0) ) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if (nEvt%100000==0) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    jetpp[0]->GetEntry(nEvt);
    
    int trgPscl[N_HLTBits]={ Jet40_p_I * jet40_l1s_ps_I ,  Jet60_p_I *  jet60_l1s_ps_I ,   
                             Jet80_p_I * jet80_l1s_ps_I , Jet100_p_I * jet100_l1s_ps_I };    
    bool trgDec[N_HLTBits]={ (bool) Jet40_I , (bool)  Jet60_I ,
			     (bool) Jet80_I , (bool) Jet100_I };    
    
    // grab largest triggerPt among HLT paths' trigger object collectiosn
    double maxTrgPt=0.,maxTrgEta=0.; 
    unsigned int trgObj40_size=trgObjpt_40->size(), trgObj60_size=trgObjpt_60->size();
    unsigned int trgObj80_size=trgObjpt_80->size(), trgObj100_size=trgObjpt_100->size();    
    
    if(trgDec[3]&&useHLT100)
      for(unsigned int itt=0; itt<trgObj100_size; ++itt)
    	if(trgObjpt_100->at(itt) > maxTrgPt) { 
	  maxTrgEta = trgObjeta_100->at(itt);
	  maxTrgPt = trgObjpt_100->at(itt); 	} 
    
    if(trgDec[2])
      for(unsigned int itt=0; itt<trgObj80_size; ++itt)
    	if(trgObjpt_80->at(itt) > maxTrgPt)  { 
	  maxTrgEta = trgObjeta_80->at(itt);
	  maxTrgPt = trgObjpt_80->at(itt); 	} 
    
    if(trgDec[1]) 
      for(unsigned int itt=0; itt<trgObj60_size; ++itt)
	if(trgObjpt_60->at(itt) > maxTrgPt)  { 
	  maxTrgEta = trgObjeta_60->at(itt);
	  maxTrgPt = trgObjpt_60->at(itt); 	} 
    
    if(trgDec[0])  
      for(unsigned int itt=0; itt<trgObj40_size; ++itt)	
	if(trgObjpt_40->at(itt) > maxTrgPt)  { 
	  maxTrgEta = trgObjeta_40->at(itt);
	  maxTrgPt = trgObjpt_40->at(itt); 	} 
    
    double trgPt=maxTrgPt; 
    double trgEta=maxTrgEta;
    
    // check trigger decisions for events + exclusivity between them, count events, assign prescale weight
    double weight_eS=0.;
    //float weight_eS = trigComb(trgDec, trgPscl, trgPt); // trig comb function replicates the procedure below
    
    bool is40  = false, is60  = false, is80  = false, is100 = false;
    if(useHLT100){
      if(      trgDec[3] && !(trgPt<HLTthresh[3])               ) 
	{ is100 = true;  weight_eS=(double)trgPscl[3]; }
      else if( trgDec[2] && !(trgPt<HLTthresh[2])  && trgPt<HLTthresh[3] ) 
	{ is80  = true;  weight_eS=(double)trgPscl[2]; }
      else if( trgDec[1] && !(trgPt<HLTthresh[1])  && trgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=(double)trgPscl[1]; }
      else if( trgDec[0] && !(trgPt<HLTthresh[0])  && trgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=(double)trgPscl[0]; }                             }
    else {
      if( trgDec[2] && !(trgPt<HLTthresh[2])  ) 
	{ is80  = true;  weight_eS=(double)trgPscl[2]; }
      else if( trgDec[1] && !(trgPt<HLTthresh[1])  && trgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=(double)trgPscl[1]; }
      else if( trgDec[0] && !(trgPt<HLTthresh[0])  && trgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=(double)trgPscl[0]; }                             }
    
    bool firedTrigger=(!((bool)weight_eS<1.) );//including exclusivity req. 
    
    //event counts, pre-duplicate skip, pre-evt cuts, pre-trigger reqs.
    h_NEvents->Fill(0.);        
    if (firedTrigger)
      h_NEvents->Fill(1.,weight_eS);  
    else //!firedTrigger
      h_NEvents->Fill(1.,1.);  //i don't want this # to represent evts not firing a jet trigger. evts not firing a jet trigger get weight 1 here.
    
    //event counts, pre-duplicate skip, pre-evt cuts, having fired a trigger
    if(firedTrigger){ 
      h_NEvents_trigd_1->Fill(0.);   h_NEvents_trigd_1->Fill(1.,weight_eS );    
      if(weight_eS==1.) { 
	h_NEvents_nopscl->Fill(0.);   	h_NEvents_nopscl->Fill(1.,weight_eS );    }    
    }
    else
      { h_NEvents_notTrigd->Fill(0.); h_NEvents_notTrigd->Fill(1.,1.);}
    
    //event counts by individual trigger, pre-duplicate skip, pre-evt cuts, having fired a trigger
    if     ( is100 )  { NEvents_100++ ; h_NEvents_jet100->Fill(0.)  ; h_NEvents_jet100->Fill(1.,weight_eS ) ; }
    else if( is80  )  { NEvents_80++  ;  h_NEvents_jet80->Fill(0.)  ;  h_NEvents_jet80->Fill(1.,weight_eS ) ; }
    else if( is60  )  { NEvents_60++  ;  h_NEvents_jet60->Fill(0.)  ;  h_NEvents_jet60->Fill(1.,weight_eS ) ; }
    else if( is40  )  { NEvents_40++  ;  h_NEvents_jet40->Fill(0.)  ;  h_NEvents_jet40->Fill(1.,weight_eS ) ; }                    
    
    //duplicate skipping between LowerJets and Jet80
    //only evts that actually fire one of the triggers can be skipped as duplicates
    if(filelistIsJet80)
      if( (bool)Jet40_I || (bool)Jet60_I || (bool)Jet40_2_I || (bool)Jet60_2_I ) {
	if(debugMode&&nEvt%1000==0)std::cout<<"this event is in Jet80 AND LowerJets dataset.!"<<std::endl;
	if(debugMode&&nEvt%1000==0)std::cout<<"Skipping event, will read it in LowerJets instead"<<std::endl;	
	h_NEvents_skipped->Fill(0.);	
	if(firedTrigger)
	  h_NEvents_skipped->Fill(1.,weight_eS);  
	else //!firedTrigger
	  h_NEvents_skipped->Fill(1., 1.);  
	continue;      }        
    
    h_NEvents_read->Fill(0.);    h_NEvents_read->Fill(1.,weight_eS);        
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilter_I==0    || 
	pBeamScrapingFilter_I==0 || 
	pprimaryvertexFilter_I==0  ) continue;    //	//	puvertexFilter_I==0  	) continue;    
    h_NEvents_skimCut->Fill(0.);  h_NEvents_skimCut->Fill(1.,weight_eS);
    
    if( fabs(vz_F)>24. )     continue;
    h_NEvents_vzCut->Fill(0.); h_NEvents_vzCut->Fill(1.,weight_eS);    
    
    if(fillDataEvtQAHists){
      
      hVz->Fill(vz_F);                        //this will be # of events post cuts, no trigger, no weights, post skim/vz Cut 
      if(weight_eS>0.)
	hTrgVz_noW->Fill(vz_F); //         ''       trigger events w/o weights post-skim/vz Cut      
      hWVz->Fill(vz_F, weight_eS);            //         ''       trigger events w/ weights, post skim/vz Cut
      if(weight_eS>0.)
	hTrgVz_wW->Fill(vz_F, weight_eS); //         ''       trigger events w/o weights post-skim/vz Cut
    }


    if( !firedTrigger ) continue;  //we didn't fire one of our triggers of interest. skip.
    else{ //firedTrigger
      h_NEvents_trigd_2->Fill(0.);   
      h_NEvents_trigd_2->Fill(1.,weight_eS);    }
    
    if(fillDataJetTrigQAHists){ //only want to fill these trigger jet plots if they pass all our quality criteria
      
      //make these trigger plots an array?
      // trgPt
      if(trgDec[0]&&!(trgPt<HLTthresh[0]))   hpp_IncHLT40trgPt->Fill(  trgPt, (double)trgPscl[0] );
      if(trgDec[1]&&!(trgPt<HLTthresh[1]))   hpp_IncHLT60trgPt->Fill(  trgPt, (double)trgPscl[1] );
      if(trgDec[2]&&!(trgPt<HLTthresh[2]))   hpp_IncHLT80trgPt->Fill(  trgPt, (double)trgPscl[2] );
      if(trgDec[3]&&!(trgPt<HLTthresh[3]))   hpp_IncHLT100trgPt->Fill( trgPt, (double)trgPscl[3] );
      
      if(     is100 )  hpp_HLT100trgPt->Fill(   trgPt , weight_eS );   
      else if(is80  )   hpp_HLT80trgPt->Fill(   trgPt , weight_eS );
      else if(is60  )   hpp_HLT60trgPt->Fill(   trgPt , weight_eS );
      else if(is40  )   hpp_HLT40trgPt->Fill(   trgPt , weight_eS );      
      if(true)         hpp_HLTCombtrgPt->Fill(  trgPt , weight_eS );   
      
      // trgEta
      if(trgDec[0]&&!(trgPt<HLTthresh[0]))   hpp_IncHLT40trgEta->Fill(  trgEta, (double)trgPscl[0] );
      if(trgDec[1]&&!(trgPt<HLTthresh[1]))   hpp_IncHLT60trgEta->Fill(  trgEta, (double)trgPscl[1] );
      if(trgDec[2]&&!(trgPt<HLTthresh[2]))   hpp_IncHLT80trgEta->Fill(  trgEta, (double)trgPscl[2] );
      if(trgDec[3]&&!(trgPt<HLTthresh[3]))   hpp_IncHLT100trgEta->Fill( trgEta, (double)trgPscl[3] );
      
      if(     is100 )  hpp_HLT100trgEta->Fill(   trgEta , weight_eS );   
      else if(is80  )  hpp_HLT80trgEta->Fill(    trgEta , weight_eS );
      else if(is60  )  hpp_HLT60trgEta->Fill(    trgEta , weight_eS );
      else if(is40  )  hpp_HLT40trgEta->Fill(    trgEta , weight_eS );      
      if(true)         hpp_HLTCombtrgEta->Fill(  trgEta , weight_eS );                   }
    
    
    // fill evt vz histo
    if(fillDataEvtQAHists){
      
      //      hVz->Fill(vz_F);                        //this will be # of events post cuts, no trigger, no weights, post skim/vz Cut 
      //      hWVz->Fill(vz_F, weight_eS);            //         ''       trigger events w/ weights, post skim/vz Cut
      //      if(weight_eS>0.)hTrgVz_noW->Fill(vz_F); //         ''       trigger events w/o weights post-skim/vz Cut
      
      hVz_BSshift->Fill((vz_F-pp5TeV2k15_BS));                        
      hWVz_BSshift->Fill((vz_F-pp5TeV2k15_BS), weight_eS);            
      if(weight_eS>0.)hTrgVz_noW_BSshift->Fill((vz_F-pp5TeV2k15_BS)); 
      
      hVz_baryshift_pos->Fill((vz_F+pp5TeV2k15_baryshift));                        
      hWVz_baryshift_pos->Fill((vz_F+pp5TeV2k15_baryshift), weight_eS);            
      if(weight_eS>0.)hTrgVz_noW_baryshift_pos->Fill((vz_F+pp5TeV2k15_baryshift)); 
      
      hVz_baryshift_neg->Fill((vz_F-pp5TeV2k15_baryshift));                        
      hWVz_baryshift_neg->Fill((vz_F-pp5TeV2k15_baryshift), weight_eS);            
      if(weight_eS>0.)hTrgVz_noW_baryshift_neg->Fill((vz_F-pp5TeV2k15_baryshift)); 

      float vr_F=std::sqrt(vx_F*vx_F + vy_F*vy_F);
      hVr->Fill(vr_F);   
      hWVr->Fill(vr_F, weight_eS);  
      if(weight_eS>0.)hTrgVr_noW->Fill(vr_F);
      
      hVx->Fill(vx_F);   
      hWVx->Fill(vx_F, weight_eS);  
      if(weight_eS>0.)hTrgVx_noW->Fill(vx_F);
      
      hVy->Fill(vy_F);   
      hWVy->Fill(vy_F, weight_eS);
      if(weight_eS>0.)hTrgVy_noW->Fill(vy_F);
      
      hNref->Fill(nref_I);
      hWNref->Fill(nref_I,weight_eS);
      //if(weight_eS>0.) hTrgNref_noW->Fill(nref_I);
    }
    
    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;

    // for event counting + avoiding duplicate fills in dijet hists
    bool dijetHistsFilled=false;
    
    bool hNEvts_withJets_jtptCut_Filled=false;
    bool hNEvts_withJets_jtetaCut1_Filled=false, hNEvts_withJets_jtetaCut2_Filled=false;
    bool hNEvts_withJets_JetIDCut_Filled=false; 
    
    int jetsPerEvent=nref_I;
    float evt_leadJetPt=-1., evt_leadJetPt_wCuts=-1.;
    
    if(nref_I>0){ h_NEvents_withJets->Fill(0.);	h_NEvents_withJets->Fill(1.,weight_eS); }    
    for(int jet = 0; jet<nref_I; ++jet){
      
      // event+jet counting
      h_NJets->Fill(0.);            h_NJets->Fill(1.,weight_eS);      
      
      float recpt  = pt_F[jet]; //this will L3 pt after corrections
      float jtpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
      float jtpt_L2Res=-1., jtpt_L3Res=-1.; //jtpt_L3Res==recpt if done right
      
      float receta = eta_F[jet];
      float absreceta = fabs(receta);
      
      if(doResidualCorr){
	
	recpt=L2JES->getCorrectedPt(recpt,receta);      
	jtpt_L2Res=recpt;
	recpt=L3JES->getCorrectedPt(recpt);        
	jtpt_L3Res=recpt;
	
	//hJetQA_jtpt->Fill(jtpt,weight_eS);
	//hJetQA_jtpt_L2Res->Fill(jtpt_L2Res,weight_eS);
	//hJetQA_jtpt_L3Res->Fill(jtpt_L3Res,weight_eS);
	//hJetQA_jtpt_v_L2Res->Fill(jtpt,jtpt_L2Res,weight_eS);
	//hJetQA_jtpt_v_L3Res->Fill(jtpt,jtpt_L3Res,weight_eS);
	//hJetQA_jtpt_L2Res_v_L3Res->Fill(jtpt_L2Res,jtpt_L3Res,weight_eS);
      }      
      
      if(recpt>evt_leadJetPt){
	evt_leadJetPt=recpt;      }
      
//// kmatCuts      
//if( recpt < jtPtCut ||
//	  recpt > jtPtCut_Hi ){ 
//	jetsPerEvent--;
//	continue;}
      if( !(recpt > jtPtCut) ){ jetsPerEvent--; continue;} //paranoia about float comparisons
      else if( !(recpt < jtPtCut_Hi) ){ jetsPerEvent--; continue;}
      
      if(!hNEvts_withJets_jtptCut_Filled){ 
	h_NEvents_withJets_jtptCut->Fill(0.);  	h_NEvents_withJets_jtptCut->Fill(1.,weight_eS);  
	hNEvts_withJets_jtptCut_Filled=true;}    
      h_NJets_jtptCut->Fill(0.); h_NJets_jtptCut->Fill(1.,weight_eS);      
      
      if( !(absreceta < 4.7) ) {       
	jetsPerEvent--; 
	continue;}      
      if(!hNEvts_withJets_jtetaCut1_Filled){ 
	h_NEvents_withJets_jtetaCut1->Fill(0.);  	h_NEvents_withJets_jtetaCut1->Fill(1.,weight_eS);  
	hNEvts_withJets_jtetaCut1_Filled=true;}    
      h_NJets_jtetaCut1->Fill(0.); h_NJets_jtetaCut1->Fill(1.,weight_eS);      
      
      if( absreceta < jtEtaCutLo ) { jetsPerEvent--;	continue;}
      if( !(absreceta < jtEtaCutHi) ) { jetsPerEvent--;	continue;}
      
      if(!hNEvts_withJets_jtetaCut2_Filled){ 
	h_NEvents_withJets_jtetaCut2->Fill(0.);  	h_NEvents_withJets_jtetaCut2->Fill(1.,weight_eS);  
	hNEvts_withJets_jtetaCut2_Filled=true;}    
      h_NJets_jtetaCut2->Fill(0.); h_NJets_jtetaCut2->Fill(1.,weight_eS);      
      
      // largest jet pt in each event
      if(recpt>evt_leadJetPt_wCuts){
	evt_leadJetPt_wCuts=recpt;      }
      
      //float rawpt  = rawpt_F[jet];
      float recphi = phi_F[jet];

      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;
      
      // 13 TeV JetID criterion, loose or tight
      bool passesJetID=false;
      //float jetIDpt=recpt;//ala HIN jetID, recpt is corrected w/ L2/L3 residuals
      float jetIDpt=jtpt;//ala HIN jetID, jtpt is what is read from the HiForest
      //float jetIDpt=rawpt_F[jet];//ala SMP jetID (which should really use energy fractions, not pt)


      if(fillDataJetIDHists) 	{
	//if (absreceta<2.4) 
	if (!(absreceta > 2.4)) 
	  { 
	    if( neSum_F[jet]/jetIDpt    < 0.99 &&
		phSum_F[jet]/jetIDpt    < 0.99 &&
		numConst              > 1    &&      
		chSum_F[jet]/jetIDpt    > 0.00 && 
		chMult                > 0    &&
		eSum_F[jet]/jetIDpt     < 0.99    ) passesJetID=true;	      
	  }
	//else if ( absreceta<2.7 && absreceta>2.4 ) 
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  {	  
	    if( neSum_F[jet]/jetIDpt    < 0.99 &&
		phSum_F[jet]/jetIDpt    < 0.99 &&
		numConst              > 1       ) passesJetID=true;	      
	  }		  
	//else if( absreceta<3.0 && absreceta>2.7 ) 
	else if( !(absreceta>3.0) && absreceta>2.7 ) 
	  {                                                         // CMSSW [76,80]X criterion
	    if(  true && //phSum_F[jet]/jetIDpt > 0.00 &&                       // else if(  phSum_F[jet]/jetIDpt [< 0.90 ] / [ > 0.01 &&]		     
		 true && //neSum_F[jet]/jetIDpt < 1.00 &&                       //           neSum_F[jet]/jetIDpt [null   ] / [ < 0.98 &&]		     
		 numConst            > 0       ) passesJetID=true;   //           neuMult            [> 2    ] / [ > 2      ] ) passesJetID=true;
	  }							      
	else //( absreceta>3.0) 
	  {                                                                                        // CMSSW 76X criterion
	    if( phSum_F[jet] < 0.4 &&         //( phSum_F[jet]/jetIDpt > 0.                      // else if( phSum_F[jet]/jetIDpt < 0.90 &&
		true && //trkSum_F[jet] < 0.4 && //  neSum_F[jet]/jetIDpt > 0. ) &&                         //          neSum_F[jet]/jetIDpt < null &&
		true       ) passesJetID=true;     //          neuMult            > 10
	  }	  	  
	
	if(!passesJetID) {
	  jetsPerEvent--;
	  continue; }
	if(!hNEvts_withJets_JetIDCut_Filled){ 
	  h_NEvents_withJets_JetIDCut->Fill(0.);  	h_NEvents_withJets_JetIDCut->Fill(1.,weight_eS);  
	  hNEvts_withJets_JetIDCut_Filled=true;}    			
	h_NJets_JetIDCut->Fill(0.);	h_NJets_JetIDCut->Fill(1.,weight_eS);
      }//end fillDataJetIDHists
      
      
      
      //fill jetspectraRapHists w/ passing jetID criterion
      if( fillDataJetSpectraRapHists ) { 
	int theRapBin=-1;
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin)
	  if( absetabins[rapbin]<=absreceta  && 		
	      absreceta<absetabins[rapbin+1]    	      ) {	    
	    theRapBin=rapbin;
	    hJetSpectraRap[jetIDint][theRapBin]->Fill(recpt,weight_eS);  	    
	    break;	  }       
      }
      
      
      // trig plots
      if(fillDataJetTrigQAHists){
	
	hpp_TrgObjComb[jetIDint]->Fill(recpt, weight_eS); //kurts method
	
	//no exclusion between triggers
	if(trgDec[0]&&!(trgPt<HLTthresh[0]))  hpp_TrgObj40[jetIDint]->Fill(recpt, (double)trgPscl[0]);
	if(trgDec[1]&&!(trgPt<HLTthresh[1]))  hpp_TrgObj60[jetIDint]->Fill(recpt, (double)trgPscl[1]);
	if(trgDec[2]&&!(trgPt<HLTthresh[2]))  hpp_TrgObj80[jetIDint]->Fill(recpt, (double)trgPscl[2]);
	if(trgDec[3]&&!(trgPt<HLTthresh[3]))  hpp_TrgObj100[jetIDint]->Fill(recpt, (double)trgPscl[3]);  
	
	//with exclusion
	if(is40 )       hpp_ExcTrgObj40[jetIDint]->Fill( recpt , weight_eS);
	else if(is60 )  hpp_ExcTrgObj60[jetIDint]->Fill( recpt , weight_eS);
	else if(is80 )  hpp_ExcTrgObj80[jetIDint]->Fill( recpt , weight_eS);
	else if(is100)  hpp_ExcTrgObj100[jetIDint]->Fill( recpt , weight_eS); 
	
	//assert(false);
      }	
            
      
      // jetQA noJetID
      if(fillDataJetQAHists){
	
	hJetQA_jtpt->Fill(jtpt,weight_eS);
	hJetQA_jtpt_L2Res->Fill(jtpt_L2Res,weight_eS);
	hJetQA_jtpt_L3Res->Fill(jtpt_L3Res,weight_eS);
	hJetQA_jtpt_v_L2Res->Fill(jtpt,jtpt_L2Res,weight_eS);
	hJetQA_jtpt_v_L3Res->Fill(jtpt,jtpt_L3Res,weight_eS);
	hJetQA_jtpt_L2Res_v_L3Res->Fill(jtpt_L2Res,jtpt_L3Res,weight_eS);
	
	int ind=0;
	//jets
	hJetQA_jtptEntries->Fill(recpt,1.);
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

	}
	
	//looking for the first two good jets that meet the criteria specified
	if(fillDataDijetHists){
	  if ( !firstGoodJetFound ){
	    if(recpt>ldJetPtCut) { firstGoodJetFound=true;
	      firstGoodJetPt =recpt; 
	      firstGoodJetPhi=recphi; 
	    } 	  }
	  
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
	    hJetQA[jetIDint][ind]->Fill( A_j , weight_eS ); ind++; 
	    hJetQA[jetIDint][ind]->Fill( x_j , weight_eS ); ind++;       
	    hJetQA[jetIDint][ind]->Fill( dphi , weight_eS ); ind++;       
	    hJetQA[jetIDint][ind]->Fill( firstGoodJetPt , weight_eS ); ind++;       
	    hJetQA[jetIDint][ind]->Fill( secondGoodJetPt , weight_eS ); ind++; }
	}
      }
      
      //if(fillDataJetIDHists){
      //	// jet/event counts
      //	h_NJets_JetIDCut->Fill(0.);	h_NJets_JetIDCut->Fill(1.,weight_eS);
      //	if(!hNEvts_withJets_JetIDCut_Filled){
      //	  h_NEvents_withJets_JetIDCut->Fill(0.);
      //	  h_NEvents_withJets_JetIDCut->Fill(1.,weight_eS);
      //	  hNEvts_withJets_JetIDCut_Filled=true;   }      
      //}//end fillDataAJetIDhists      
      
    }//end jet loop
    
    if(fillDataEvtQAHists){
      hjetsPEvt ->Fill(jetsPerEvent,1.0); 
      hWjetsPEvt->Fill(jetsPerEvent,weight_eS);     
      if(evt_leadJetPt>0.) 	  hLeadJetPt->Fill(evt_leadJetPt ,weight_eS);
      if(evt_leadJetPt_wCuts>0.)  hLeadJetPt_wCuts->Fill(evt_leadJetPt_wCuts ,weight_eS);
    }
        
  }//end event loop

  std::cout<<std::endl;  
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"/// Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////"<<std::endl;
  std::cout<<"# of Events (criteria/description) = count / pscl weighted count" <<std::endl;

  std::cout<<std::endl;
  std::cout<<"pre-duplicate skip, pre-trigger reqs, pre-evt cuts"<<std::endl;
  std::cout<<"# evts in file(s) opened = "    <<h_NEvents->GetBinContent(1)        << " / "<< h_NEvents->GetBinContent(2)<<std::endl;   
  std::cout<<"# evts trigd             = "    <<h_NEvents_trigd_1->GetBinContent(1) << " / "<< h_NEvents_trigd_1->GetBinContent(2)<<std::endl;
  std::cout<<"# evts w/ pscl==1        = " << h_NEvents_nopscl->GetBinContent(1) << " / " << h_NEvents_nopscl->GetBinContent(2) <<std::endl ;
  std::cout<<"# evts NOT trigd         = "    <<h_NEvents_notTrigd->GetBinContent(1) << " / "<< h_NEvents_notTrigd->GetBinContent(2)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Trigger Counts summary, pre-duplicate skip, pre-evt cuts"<<std::endl;
  std::cout<<"Trigger Combination type = " << hTrigComb->GetTitle() << std::endl;
  std::cout<<"Trigger Thresholds for 40/60/80/100 = " << HLTthresh[0] <<"/ "<< HLTthresh[1]<<"/ "<< HLTthresh[2]<<"/ "<< HLTthresh[3]<<std::endl;
  std::cout<<std::endl;    
  std::cout<<"# evts passing is40      = " << h_NEvents_jet40->GetBinContent(1)  << " / " << h_NEvents_jet40->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is60      = " << h_NEvents_jet60->GetBinContent(1)  << " / " << h_NEvents_jet60->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is80      = " << h_NEvents_jet80->GetBinContent(1)  << " / " << h_NEvents_jet80->GetBinContent(2)  <<std::endl;
  std::cout<<"# evts passing is100     = " << h_NEvents_jet100->GetBinContent(1) << " / " << h_NEvents_jet100->GetBinContent(2) <<std::endl;
  std::cout<<std::endl;
  
  std::cout<<"duplicate-skip and trigger req.(s)s applied, pre-evt cuts"<<std::endl;
  std::cout<<"# evts in file(s) opened                     = "    <<h_NEvents->GetBinContent(1)        << " / "<< h_NEvents->GetBinContent(2)<<std::endl; 
  std::cout<<"# evts actually read from those file(s)      = " << h_NEvents_read->GetBinContent(1)    << " / " << h_NEvents_read->GetBinContent(2)<<std::endl;  
  std::cout<<"# evts skipped from those file(s)            = " << h_NEvents_skipped->GetBinContent(1) << " / " << h_NEvents_skipped->GetBinContent(2)<<std::endl;
  std::cout<<"# evts passing a trigger post duplicate skip = " << h_NEvents_trigd_2->GetBinContent(1) << " / " << h_NEvents_trigd_2->GetBinContent(2) <<std::endl;
  std::cout << std::endl;
  
  
  std::cout<<"evt-cuts applied, post duplicate-skip, post trigger req.(s)"<<std::endl;
  std::cout<<"# evts passing skimCuts  = " << h_NEvents_skimCut->GetBinContent(1) << " / " << h_NEvents_skimCut->GetBinContent(2) << std::endl;
  std::cout<<"# evts passing vzCuts    = " << h_NEvents_vzCut->GetBinContent(1)   << " / " << h_NEvents_vzCut->GetBinContent(2)   << std::endl; 
  std::cout << std::endl;
  
  std::cout<<"post evt-cuts, post duplicate-skip, post trigger req.(s)"<<std::endl;
  std::cout<<"# evts passing trig  = " << hpp_HLTCombtrgPt ->GetEntries() << " / " <<  hpp_HLTCombtrgPt->Integral()<<std::endl;
  std::cout<<"# evts passing is40  = " << hpp_HLT40trgEta ->GetEntries()  << " / " <<  hpp_HLT40trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is60  = " << hpp_HLT60trgEta ->GetEntries()  << " / " <<  hpp_HLT60trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is80  = " << hpp_HLT80trgEta ->GetEntries()  << " / " <<  hpp_HLT80trgEta->Integral()<<std::endl;
  std::cout<<"# evts passing is100 = " << hpp_HLT100trgEta->GetEntries()  << " / " <<  hpp_HLT100trgEta->Integral()<<std::endl<<std::endl;  

  std::cout<<"# evts passing HLT40  = " << hpp_IncHLT40trgEta ->GetEntries()  << " / " <<  hpp_IncHLT40trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT60  = " << hpp_IncHLT60trgEta ->GetEntries()  << " / " <<  hpp_IncHLT60trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT80  = " << hpp_IncHLT80trgEta ->GetEntries()  << " / " <<  hpp_IncHLT80trgEta ->Integral()<<std::endl;
  std::cout<<"# evts passing HLT100 = " << hpp_IncHLT100trgEta->GetEntries()  << " / " <<  hpp_IncHLT100trgEta->Integral()<<std::endl<<std::endl;  
  
  
  std::cout<<"//////////////////////////////"<<std::endl;
  std::cout<<"/// Jet-Event Loop Summary ///"<<std::endl;
  std::cout<<"//////////////////////////////"<<std::endl;

  std::cout<<"# of Evts w/ min. 1 jet                   = " << h_NEvents_withJets->GetBinContent(1) << " / " << h_NEvents_withJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Evts w/ min. 1 jet w/ jtpt > jtPtCut = " << h_NEvents_withJets_jtptCut->GetBinContent(1) << " / " << h_NEvents_withJets_jtptCut->GetBinContent(2) << std::endl;
  std::cout<<"# of Evts w/ min. 1 jet w/ |jteta| < 4.7  = "                 << h_NEvents_withJets_jtetaCut1->GetBinContent(1) << " / " << h_NEvents_withJets_jtetaCut1->GetBinContent(2) << std::endl; 
  std::cout<<"# of Evts w/ min. 1 jet w/ jtEtaCut(s)    = " << h_NEvents_withJets_jtetaCut2->GetBinContent(1) << " / " << h_NEvents_withJets_jtetaCut2->GetBinContent(2) << std::endl;  
  if(fillDataJetIDHists) 
    std::cout<<"# of Evts w/ min. 1 jet w/ JetIDCut       = " << h_NEvents_withJets_JetIDCut->GetBinContent(1) << " / " << h_NEvents_withJets_JetIDCut->GetBinContent(2) << std::endl;        

  std::cout<<std::endl;    
  std::cout<<"# of Jets read                   = " << h_NJets->GetBinContent(1) << " / " << h_NJets->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Jets w/ jtpt > jtptCut      = " << h_NJets_jtptCut->GetBinContent(1) << " / " << h_NJets_jtptCut->GetBinContent(2) <<std::endl;   
  std::cout<<"# of Jets w/ |jteta|<4.7         = " << h_NJets_jtetaCut1->GetBinContent(1) << " / " << h_NJets_jtetaCut1->GetBinContent(2) <<std::endl; 
  std::cout<<"# of Jets w/ passing jtEtaCut(s) = " << h_NJets_jtetaCut2->GetBinContent(1) << " / " << h_NJets_jtetaCut2->GetBinContent(2) <<std::endl; 
  if(fillDataJetIDHists)
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
  std::cout<<"Real tyime (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
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
  
