// custom header
#include "readForests_findEvt.h"

// ppData switches
const bool fillDataEvtQAHists=true;
const bool fillDataJetQAHists=false;
const bool fillDataDijetHists=false;
const bool fillDataJetIDHists=false;//, tightJetID=false;

const bool fillDataJetTrigQAHists=false; //data-specific
const bool fillDataJetSpectraRapHists=false; //other

const bool useHLT100=false;

//const std::string trgCombType="Calo";
const std::string trgCombType="PF";

//// LowerJets
//const int targetEvt=113080452;
//const int targetRun=262252;
//const int targetEvt=112211851;
//const int targetRun=262252;
//const int targetEvt=6608245;
//const int targetRun=262163;
//const int targetEvt=197246756;
//const int targetRun=262271;
//const int targetEvt=1000276728;
//const int targetRun=262274;
const int targetEvt=190748136;
const int targetRun=262274;
const int maxEvtsToPrint=100;
const bool onelinemode=true;


//Run:262252, event: 112211851,
//const int targetLumi=;

//// readForests_ppData_findEvt
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppData_findEvt( std::string inFilelist , int startfile , int endfile , 
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
  

  //check out what kind of filelist we're using
  bool filelistIsLowerJets=false; //to know if i should apply the duplicate check or not
  if( inFilelist.find("LowerJets")!=std::string::npos ) {
    filelistIsLowerJets=true; 
    std::cout<<"running LowerJets filelist!"<<std::endl; }
  bool filelistIsJet80=false;
  if( inFilelist.find("Jet80")    !=std::string::npos ) {
    filelistIsJet80=true;
    std::cout<<"running Jet80 filelist!"<<std::endl; }
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
  if(debugMode&&false)for(int i=0;i<N_dataTrees;++i)std::cout<<"trees[i="<<i<<"]="<<trees[i]<<std::endl;  

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
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"      ,(std::to_string(jtPtCut)).c_str()   ,    100, 0,100); hJetPtCut->Fill(jtPtCut);           
  TH1F *hJetEtaCutHi       =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6.  ); hJetEtaCutHi->Fill(jtEtaCutHi); 
  TH1F *hJetEtaCutLo       =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6.  ); hJetEtaCutLo->Fill(jtEtaCutLo);     
  TH1F *hJetQAPtCut      =new TH1F("hJetQAPtCut"    ,(std::to_string(jetQAPtCut)).c_str(),    100, 0,100); hJetQAPtCut->Fill(jetQAPtCut);     
  TH1F *hLeadJetPtCut    =new TH1F("hLdJetPtCut"    ,(std::to_string(ldJetPtCut)).c_str(),    100, 0,100); hLeadJetPtCut->Fill(ldJetPtCut);     
  TH1F *hSubLeadJetPtCut =new TH1F("hSubldJetPtCut" ,(std::to_string(subldJetPtCut)).c_str(), 100, 0,100); hSubLeadJetPtCut->Fill(subldJetPtCut);
  TH1F *hPtAveCut        =new TH1F("hPtAveCut"      ,(std::to_string(ptAveCut)).c_str(),      100, 0,100); hPtAveCut->Fill(ptAveCut);  
  

  //evt counts
  TH1F *h_NEvents          = new TH1F("NEvents"         , "NEvents",         1,0,2);
  TH1F *h_NEvents_read     = new TH1F("NEvents_read"    , "NEvents read",    1,0,2);
  TH1F *h_NEvents_skipped  = new TH1F("NEvents_skipped" , "NEvents skipped", 1,0,2);
  TH1F *h_NEvents_skimCut  = new TH1F("NEvents_skimCut" , "NEvents read post skimCut",  1,0,2);
  TH1F *h_NEvents_vzCut    = new TH1F("NEvents_vzCut"   , "NEvents read post vzCut AND skimCut", 1,0,2);

  TH1F *h_NEvents_trigd    = new TH1F("NEvents_trigd"   , "NEvents trigd",   1,0,2);
  TH1F *h_NEvents_jet40    = new TH1F("NEvents_jet40Trigd"   , "NEvents jet40Trigd",    1,0,2);
  TH1F *h_NEvents_jet60    = new TH1F("NEvents_jet60Trigd"   , "NEvents jet60Trigd",    1,0,2);
  TH1F *h_NEvents_jet80    = new TH1F("NEvents_jet80Trigd"   , "NEvents jet80Trigd",    1,0,2);
  TH1F *h_NEvents_jet100   = new TH1F("NEvents_jet100Trigd"  , "NEvents jet100Trigd",   1,0,2);

  //evt counts post all evt cuts. thrown out because no jets left in the event pre or post cuts
  TH1F *h_NEvents_withJets           = new TH1F("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1F *h_NEvents_withJets_kmatCut   = new TH1F("NEvents_withJets_kmatCut" , 
						"NEvents read post evt cuts, w/ jets post kmatCut", 1,0,2);
  TH1F *h_NEvents_withJets_JetIDCut = NULL;
  if(fillDataJetIDHists)
    h_NEvents_withJets_JetIDCut= new TH1F("NEvents_withJets_JetIDCut" , 
					  "NEvents read post evt cuts, w/ jets post kmatCut AND JetID Cut", 1,0,2);
  
  //EvtQA, i.e. hists filled once per event
  TH1F *hVr=NULL, *hWVr=NULL, *hTrgVr_noW=NULL;
  TH1F *hVx=NULL, *hWVx=NULL, *hTrgVx_noW=NULL;
  TH1F *hVy=NULL, *hWVy=NULL, *hTrgVy_noW=NULL;
  TH1F *hVz=NULL, *hWVz=NULL, *hTrgVz_noW=NULL;
  
  TH1F* hNref=NULL, *hWNref=NULL;
  TH1F* hjetsPEvt=NULL,* hWjetsPEvt=NULL;
  TH1F* hjetsPEvt_wJetID=NULL,* hWjetsPEvt_wJetID=NULL;
  TH1F* hLeadJetPt=NULL;//,* hLeadJetPt_wJetID=NULL;

  if(fillDataEvtQAHists){

    hVr = new TH1F("hVr","vr, no trig, no weights",  400,-0.25,0.25); 
    hWVr = new TH1F("hWeightedVr","vr, trigd, with weights",  400,-0.25,0.25);    
    hTrgVr_noW = new TH1F("hTriggerVr_noWeights","vr, trigd, no weights",  400,-0.25,0.25);    

    hVx = new TH1F("hVx","vx, no trig, no weights",  400,-0.25,0.25); //bin size of .0025
    hWVx = new TH1F("hWeightedVx","vx, trigd, with weights",  400,-0.25,0.25);    
    hTrgVx_noW = new TH1F("hTriggerVx_noWeights","vx, trigd, no weights",  400,-0.25,0.25);    

    hVy = new TH1F("hVy","vy, no trig, no weights",  400,-0.25,0.25); 
    hWVy = new TH1F("hWeightedVy","vy, trigd, with weights",  400,-0.25,0.25);    
    hTrgVy_noW = new TH1F("hTriggerVy_noWeights","vy, trigd, no weights",  400,-0.25,0.25);    

    hVz = new TH1F("hVz","vz, no trig, no weights",  100,-25.,25.); //bin size of .5
    hWVz = new TH1F("hWeightedVz","vz, trigd, with weights",  100,-25.,25.);    
    hTrgVz_noW = new TH1F("hTriggerVz_noWeights","vz, trigd, no weights",  100,-25.,25.);  

    hNref      = new TH1F("hNref","nref each evt",30,0,30);
    hWNref     = new TH1F("hWNref","weighted nref each evt",30,0,30);
    hjetsPEvt  = new TH1F("hjetsPEvt","numJets each evt",30,0,30);
    hWjetsPEvt = new TH1F("hWjetsPEvt","weighted numJets each evt",30,0,30);

    if(fillDataJetIDHists){
      hjetsPEvt_wJetID = new TH1F("hjetsPEvt_wJetID","num JetID Jets each evt",30,0,30);
      hWjetsPEvt_wJetID = new TH1F("hWjetsPEvt_wJetID","weighted num JetID Jets each evt",30,0,30);}
    
    hLeadJetPt      = new TH1F("hLeadJetPt","lead jet pt for each evt",1000,0,1000);
    //if(fillDataJetIDHists)
    //  hLeadJetPt_wJetID      = new TH1F("hLeadJetPt_wJetId","lead jet pt w Jet ID for each evt",1000,0,1000);

  }
  
  

  //jet counts
  TH1F *h_NJets          = new TH1F("NJets","NJets read", 1,0,2);
  TH1F *h_NJets_kmatCut  = new TH1F("NJets_kmatCut ","NJets read post kmatCut ", 1,0,2);
  TH1F *h_NJets_JetIDCut=NULL;
  if(fillDataJetIDHists)
    h_NJets_JetIDCut = new TH1F("NJets_JetIDCut","NJets read post JetIDCut AND kmatCut", 1,0,2);

  //weighted jet counts
  TH1F *h_WNJets          = new TH1F("WNJets","WNJets read", 1,0,2);
  TH1F *h_WNJets_kmatCut  = new TH1F("WNJets_kmatCut ","WNJets read post kmatCut ", 1,0,2);
  TH1F *h_WNJets_JetIDCut=NULL;
  if(fillDataJetIDHists)
    h_WNJets_JetIDCut = new TH1F("WNJets_JetIDCut","WNJets read post JetIDCut AND kmatCut", 1,0,2);
  
  //jet QA, i.e. hists filled once per jet
  TH1F *hJetQA[2][N_vars]={};
  if(fillDataJetQAHists)
    for(int k = 0; k<2; ++k){
      if(!fillDataJetIDHists && k==1)continue;	
      for(int j = 0; j<N_vars; ++j){	

	//jets
	if(var[j]=="jtpt"||var[j]=="rawpt")
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;", var[j].c_str()) , 1000,0,1000);       
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
	else if (fillDataDijetHists){
	  if(var[j]=="dphi") 
	    hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()) , Form(";%s;",var[j].c_str()) , 50,0,+4);
	  else if(var[j]=="leadJetPt"||var[j]=="subleadJetPt")
	    hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;", var[j].c_str()), 1000,0,1000);
	  else //xj and Aj binnings
	    hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);         
	}

	//consituent binnings
	else 
	  hJetQA[k][j] = new TH1F( Form("hJetQA_%dwJetID_%s", k,var[j].c_str()), Form(";%s;",var[j].c_str()), 200,0,2);         
      } }
  
  TH1F *hJetSpectraRap[2][nbins_abseta]={};
  if(fillDataJetSpectraRapHists)
    for(int k = 0; k<2; ++k) {
      if(!fillDataJetIDHists && k==1)continue;	
      for(int j = 0; j<nbins_abseta; ++j){
	std::string h_Title="hJetSpectraRap_";
	if(k==1)h_Title+="wJetID_";
	h_Title+="bin"+std::to_string(j);
	std::string h_Desc;
	std::stringstream stream1, stream2;	
	stream1.precision(1); stream1 << std::fixed << absetabins[j];
	stream2.precision(1); stream2 << std::fixed << absetabins[j+1];
	h_Desc="JetPt Spectra for "+stream1.str()+"<abs(y)<"+ stream2.str();	
	hJetSpectraRap[k][j]=new TH1F(h_Title.c_str(),h_Desc.c_str(), 1000,0,1000);  
      } }

  // trigger (also jet) level
  TH1F  *hpp_TrgObj40[2]={}, *hpp_TrgObj60[2]={}, *hpp_TrgObj80[2]={}, *hpp_TrgObj100[2]={};// plots for incl spectra for each individual trigger
  TH1F  *hpp_ExcTrgObj40[2]={}, *hpp_ExcTrgObj60[2]={}, *hpp_ExcTrgObj80[2]={}, *hpp_ExcTrgObj100[2]={};// plots for trigger-exclusive inclusive jet spectra
  TH1F *hpp_TrgObjComb[2]={};//combination of triggers w/ exclusion
  
  TH1F  *hpp_IncHLT40trgPt =NULL ,   *hpp_IncHLT40trgEta =NULL;
  TH1F  *hpp_IncHLT60trgPt =NULL ,   *hpp_IncHLT60trgEta =NULL;
  TH1F  *hpp_IncHLT80trgPt =NULL ,   *hpp_IncHLT80trgEta =NULL;
  TH1F *hpp_IncHLT100trgPt =NULL ,  *hpp_IncHLT100trgEta =NULL;


  TH1F  *hpp_HLT40trgPt =NULL ,   *hpp_HLT40trgEta =NULL;
  TH1F  *hpp_HLT60trgPt =NULL ,   *hpp_HLT60trgEta =NULL;
  TH1F  *hpp_HLT80trgPt =NULL ,   *hpp_HLT80trgEta =NULL;
  TH1F *hpp_HLT100trgPt =NULL ,  *hpp_HLT100trgEta =NULL;
  
  TH1F *hpp_HLTCombtrgPt=NULL, *hpp_HLTCombtrgEta=NULL;
  
  //TH2F *hpp_pscl_HLTPFvCalo=NULL;
  //TH2F *hpp_totpscl_HLTPFvCalo=NULL;
  //TH2F *hpp_pscl_HLTvL1=NULL;

  if(fillDataJetTrigQAHists){

    ////trigger prescale TH2s
    //if(false)hpp_totpscl_HLTPFvCalo = new TH2F ( "totpsclTH2_HLTPFvCalo" , "pscls for PForCalo trigd evts"    , 101, -1,100, 101, -1,100);
    //if(false)hpp_pscl_HLTPFvCalo    = new TH2F ( "psclTH2_HLTPFvCalo"    , "totpscls for PForCalo trigd evts" , 101, -1,100, 101, -1,100);
    //if(false)hpp_pscl_HLTvL1        = new TH2F ( "psclTH2_HLTvL1"        , "pscls for HLTorL1 trigd evts"     , 101, -1,100, 101, -1,100);

    //no exclusion trigger pt
    hpp_IncHLT40trgPt   = new TH1F( "IncHLT40_trgPt"  , "inc trgPt for HLT40   ", 1000, 0, 1000 ); 
    hpp_IncHLT60trgPt   = new TH1F( "IncHLT60_trgPt"  , "inc trgPt for HLT60   ", 1000, 0, 1000 ); 
    hpp_IncHLT80trgPt   = new TH1F( "IncHLT80_trgPt"  , "inc trgPt for HLT80   ", 1000, 0, 1000 ); 
    hpp_IncHLT100trgPt  = new TH1F( "IncHLT100_trgPt" , "inc trgPt for HLT100  ", 1000, 0, 1000 );       

    //with exclusion trigger pt
    hpp_HLT40trgPt   = new TH1F( "isHLT40_trgPt"  , "exc trgPt for HLT40   ", 1000, 0, 1000 ); 
    hpp_HLT60trgPt   = new TH1F( "isHLT60_trgPt"  , "exc trgPt for HLT60   ", 1000, 0, 1000 ); 
    hpp_HLT80trgPt   = new TH1F( "isHLT80_trgPt"  , "exc trgPt for HLT80   ", 1000, 0, 1000 ); 
    hpp_HLT100trgPt  = new TH1F( "isHLT100_trgPt" , "exc trgPt for HLT100  ", 1000, 0, 1000 );       

    hpp_HLTCombtrgPt = new TH1F( "HLTComb_trgPt"  , "trgPt for HLTComb ", 1000, 0, 1000 );       

    //no exclusion trigger eta
    hpp_IncHLT40trgEta   = new TH1F( "IncHLT40_trgEta"  , "inc trgEta for HLT40   ", 100, -5., 5. ); 
    hpp_IncHLT60trgEta   = new TH1F( "IncHLT60_trgEta"  , "inc trgEta for HLT60   ", 100, -5., 5. ); 
    hpp_IncHLT80trgEta   = new TH1F( "IncHLT80_trgEta"  , "inc trgEta for HLT80   ", 100, -5., 5. ); 
    hpp_IncHLT100trgEta  = new TH1F( "IncHLT100_trgEta" , "inc trgEta for HLT100  ", 100, -5., 5. );       

    //with exclusion trigger eta
    hpp_HLT40trgEta   = new TH1F( "isHLT40_trgEta"  , "exc trgEta for HLT40   ", 100, -5., 5. ); 
    hpp_HLT60trgEta   = new TH1F( "isHLT60_trgEta"  , "exc trgEta for HLT60   ", 100, -5., 5. ); 
    hpp_HLT80trgEta   = new TH1F( "isHLT80_trgEta"  , "exc trgEta for HLT80   ", 100, -5., 5. ); 
    hpp_HLT100trgEta  = new TH1F( "isHLT100_trgEta" , "exc trgEta for HLT100  ", 100, -5., 5. );       

    hpp_HLTCombtrgEta = new TH1F( "HLTComb_trgEta"  , "trgEta for HLTComb ", 100, -5., 5. );       

    std::string hJetTrigQATitleArray[]={ "hpp_HLT40", "hpp_HLT60",  "hpp_HLT80", "hpp_HLT100", "hpp_HLTComb"};//,"hpp_TrgCombTest" };
    std::string hJetTrigQADescArray[]={ "HLT40Jet","HLT60Jet",  "HLT80Jet", "HLT100Jet", "HLTCombo"};//,"HLTAltComb" };
    const int nJetTrigQATitles=sizeof(hJetTrigQATitleArray)/sizeof(std::string);
    
    std::string hJetExcTrigQATitleArray[]={ "hpp_excHLT40", "hpp_excHLT60",  "hpp_excHLT80", "hpp_excHLT100"};//, "hpp_excHLTComb"};
    std::string hJetExcTrigQADescArray[]={ "ExcHLT40Jet","ExcHLT60Jet",  "ExcHLT80Jet", "ExcHLT100Jet"};
    const int nJetExcTrigQATitles=sizeof(hJetExcTrigQATitleArray)/sizeof(std::string);
    
    for(int j=0;j<2;j++){      
      
      if(!fillDataJetIDHists && j==1) continue;    
      if(debugMode&&true) std::cout<<"j=="<<j<<std::endl;
      
      for(int k=0; k<nJetTrigQATitles; k++){
	
	if(debugMode&&true) std::cout<<"k=="<<k<<std::endl;	
	
	std::string hTitle=hJetTrigQATitleArray[k];
	std::string hDesc=hJetTrigQADescArray[k];
	if(j==1)hTitle+="_wJetID";
	if(j==1)hDesc+=" w/ JetID";
	hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;
	
	if(debugMode&&true)std::cout<<"inclusive HLT plot hTitle="<<hTitle<<std::endl;
	if(debugMode&&true)std::cout<<"inclusive HLT plot hDesc="<<hDesc<<std::endl;
	
	if     (hJetTrigQATitleArray[k]=="hpp_HLT40"  )  hpp_TrgObj40[j]     = new TH1F( hTitle.c_str(), (hDesc).c_str(), 1000,0,1000);
	else if(hJetTrigQATitleArray[k]=="hpp_HLT60"  )  hpp_TrgObj60[j]     = new TH1F( hTitle.c_str(), (hDesc).c_str(), 1000,0,1000);
	else if(hJetTrigQATitleArray[k]=="hpp_HLT80"  )  hpp_TrgObj80[j]     = new TH1F( hTitle.c_str(), (hDesc).c_str(), 1000,0,1000);
	else if(hJetTrigQATitleArray[k]=="hpp_HLT100" )  hpp_TrgObj100[j]    = new TH1F( hTitle.c_str(), (hDesc).c_str(), 1000,0,1000);	
	else if(hJetTrigQATitleArray[k]=="hpp_HLTComb")  hpp_TrgObjComb[j]   = new TH1F( hTitle.c_str(), (hDesc).c_str(), 1000,0,1000);
	
	if( k<nJetExcTrigQATitles ){
	  
	  std::string hExcTitle=hJetExcTrigQATitleArray[k];
	  std::string hExcDesc=hJetExcTrigQADescArray[k];
	  if(j==1)hExcTitle+="_wJetID";
	  if(j==1)hExcDesc+=" w/ JetID";
	  hExcTitle+="_R"+std::to_string(radius)+"_"+etaWidth;
	  if(debugMode&&true)std::cout<<" exclusive HLT plothExcTitle="<<hExcTitle<<std::endl;
	  if(debugMode&&true)std::cout<<" exclusive HLT plothExcDesc="<<hExcDesc<<std::endl<<std::endl;
	  
	  if     (hJetExcTrigQATitleArray[k]=="hpp_excHLT40"  )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj40[j]     = new TH1F( hExcTitle.c_str(), (hExcDesc).c_str(), 1000,0,1000);}
	  else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT60"  )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj60[j]     = new TH1F( hExcTitle.c_str(), (hExcDesc).c_str(), 1000,0,1000);}
	  else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT80"  )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj80[j]     = new TH1F( hExcTitle.c_str(), (hExcDesc).c_str(), 1000,0,1000);}
	  else if(hJetExcTrigQATitleArray[k]=="hpp_excHLT100" )  {
	    std::cout<<"booking hist "<<hJetExcTrigQATitleArray[k]<<std::endl;
	    hpp_ExcTrgObj100[j]    = new TH1F( hExcTitle.c_str(), (hExcDesc).c_str(), 1000,0,1000);}	  
	}
	
      } } }

  // EVENT LOOP + PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  float rawpt_F[1000];
  float eta_F[1000];   float phi_F[1000];  
  float y_F[1000]; //float jtpu_F[1000]; 
  
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

  ////other
  //float hcalSum_F[1000], ecalSum_F[1000];

  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("jty",&y_F);

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

  ////ecal/hcal sums from PF?
  //jetpp[0]->SetBranchAddress("hcalSum",&hcalSum_F);
  //jetpp[0]->SetBranchAddress("ecalSum",&ecalSum_F);

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

  //jetpp[8]->SetBranchAddress("trkEta",&trkEta_F);
  //jetpp[8]->SetBranchAddress("nTrk",&nTrk_I);


  // hltanalysis
  assert(N_HLTBits==N_L1Bits);

  //string array formation for HLT branches
  std::string L1Branches[N_L1Bits], L1PresclBranches[N_L1Bits];
  std::string PF_HLTBranches[N_HLTBits],PF_HLTPresclBranches[N_HLTBits];
  std::string Calo_HLTBranches[N_HLTBits], Calo_HLTPresclBranches[N_HLTBits];

  for(int i=0; i<N_HLTBits; ++i) {
    L1Branches[i] = L1BitStrings[i];
    L1PresclBranches[i] = L1BitStrings[i] + "_Prescl";  
    PF_HLTBranches[i] = PF_HLTBitStrings[i] + "_v1"; 
    PF_HLTPresclBranches[i] = PF_HLTBitStrings[i] + "_v1_Prescl";   
    Calo_HLTBranches[i] = Calo_HLTBitStrings[i] + "_v1"; 
    Calo_HLTPresclBranches[i] = Calo_HLTBitStrings[i] + "_v1_Prescl"; 
    
    if(debugMode&&false){std::cout<<std::endl;
      std::cout<<"L1Branches      [i="<<i<<"]= "<<L1Branches[i]<<std::endl;
      std::cout<<"L1PresclBranches      [i="<<i<<"]= "<<L1PresclBranches[i]<<std::endl;
      std::cout<<"PF_HLTBranches        [i="<<i<<"]= "<<PF_HLTBranches[i]<<std::endl;
      std::cout<<"PF_HLTPresclBranches  [i="<<i<<"]= "<<PF_HLTPresclBranches[i]<<std::endl;
      std::cout<<"Calo_HLTBranches      [i="<<i<<"]= "<<Calo_HLTBranches[i]<<std::endl;
      std::cout<<"Calo_HLTPresclBranches[i="<<i<<"]= "<<Calo_HLTPresclBranches[i]<<std::endl<<std::endl; }  }
  
  int jet40_l1s_I, jet60_l1s_I,  jet80_l1s_I,  jet100_l1s_I;   //L1 
  int jet40_l1s_ps_I, jet60_l1s_ps_I,  jet80_l1s_ps_I,  jet100_l1s_ps_I;   
  int PFJet40_I, PFJet60_I, PFJet80_I, PFJet100_I;   //PF HLT
  int PFJet40_p_I, PFJet60_p_I, PFJet80_p_I, PFJet100_p_I;
  int CaloJet40_I, CaloJet60_I, CaloJet80_I, CaloJet100_I;   //Calo HLT
  int CaloJet40_p_I, CaloJet60_p_I, CaloJet80_p_I, CaloJet100_p_I;

  jetpp[4]->SetBranchAddress( L1Branches[0].c_str()  , &jet40_l1s_I);
  jetpp[4]->SetBranchAddress( L1Branches[1].c_str()  , &jet60_l1s_I);
  jetpp[4]->SetBranchAddress( L1Branches[2].c_str()  , &jet80_l1s_I);
  jetpp[4]->SetBranchAddress( L1Branches[3].c_str()  , &jet100_l1s_I);

  jetpp[4]->SetBranchAddress( L1PresclBranches[0].c_str()  , &jet40_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[1].c_str()  , &jet60_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[2].c_str()  , &jet80_l1s_ps_I);
  jetpp[4]->SetBranchAddress( L1PresclBranches[3].c_str()  , &jet100_l1s_ps_I);

  jetpp[4]->SetBranchAddress( PF_HLTBranches[0].c_str() , &PFJet40_I);
  jetpp[4]->SetBranchAddress( PF_HLTBranches[1].c_str() , &PFJet60_I);
  jetpp[4]->SetBranchAddress( PF_HLTBranches[2].c_str() , &PFJet80_I);
  jetpp[4]->SetBranchAddress( PF_HLTBranches[3].c_str() , &PFJet100_I);

  jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[0].c_str() , &PFJet40_p_I);
  jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[1].c_str() , &PFJet60_p_I);
  jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[2].c_str() , &PFJet80_p_I);
  jetpp[4]->SetBranchAddress( PF_HLTPresclBranches[3].c_str() , &PFJet100_p_I);

  jetpp[4]->SetBranchAddress( Calo_HLTBranches[0].c_str() , &CaloJet40_I);
  jetpp[4]->SetBranchAddress( Calo_HLTBranches[1].c_str() , &CaloJet60_I);
  jetpp[4]->SetBranchAddress( Calo_HLTBranches[2].c_str() , &CaloJet80_I);
  jetpp[4]->SetBranchAddress( Calo_HLTBranches[3].c_str() , &CaloJet100_I);

  jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[0].c_str() , &CaloJet40_p_I);
  jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[1].c_str() , &CaloJet60_p_I);
  jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[2].c_str() , &CaloJet80_p_I);
  jetpp[4]->SetBranchAddress( Calo_HLTPresclBranches[3].c_str() , &CaloJet100_p_I);

  double HLTthresh[4]={55.,75.,105.,135.};//hard coded temporarily to answer john quickly
  




  //ONE HLT path ONE tree ONE trig obj pt branch
  //e.g. trgObjpt_40 is filled with jet pt from the specific jet40 HLT tree/branch 
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
  for(UInt_t i=0;i < NEvents_allFiles; ++i) h_NEvents->Fill(1);
  
  UInt_t NEvents_toRead=0;
  NEvents_toRead = NEvents_allFiles;
  
  UInt_t NEvents_40=0, NEvents_60=0, NEvents_80=0, NEvents_100=0;  //trigger event counts
  std::cout<<"reading "<<NEvents_toRead<<" events"<<std::endl;  
  
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
  int numevtsprinted=0;
  bool targetEvtFound=false;
  for(UInt_t nEvt = 0; nEvt < NEvents_toRead; ++nEvt) {
    //if(targetEvtFound)break;
    //if(targetEvtFound&&nEvt%100==0)continue;
    if(numevtsprinted>maxEvtsToPrint){
      std::cout<<"max # evts printed reached; maxEvts="<<maxEvtsToPrint<<std::endl;
      break;
    }
    //grab an entry
    if( debugMode && (nEvt%1000==0) ) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    else if (nEvt%10000==0) std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;
    jetpp[0]->GetEntry(nEvt);
    




    //duplicate skipping between LowerJets and Jet80
    if(filelistIsJet80)
      if( (bool)CaloJet40_I || (bool)CaloJet60_I || (bool)PFJet40_I || (bool)PFJet60_I ) {
	if(debugMode&&nEvt%100==0)std::cout<<"this event is in Jet80 AND LowerJets dataset.!"<<std::endl;
	if(debugMode&&nEvt%100==0)std::cout<<"Skipping event, will read it in LowerJets instead"<<std::endl;	
	h_NEvents_skipped->Fill(1);
	continue; }
    h_NEvents_read->Fill(1);
    



    //if(true) { std::cout<<std::endl<<std::endl<<"target evt/run/lumi found"<<std::endl<<std::endl;
    //    if(evt_I==targetEvt && run_I==targetRun && lumi_I==targetLumi) std::cout<<std::endl<<std::endl<<"target evt/run/lumi found"<<std::endl<<std::endl;
    if((evt_I==targetEvt && run_I==targetRun) || targetEvtFound) { 
      
      if(evt_I==targetEvt && run_I==targetRun){
	std::cout<<std::endl<<std::endl<<"target evt/run/lumi found"<<std::endl<<std::endl;
	if(onelinemode)std::cout<<"format for one-line is:"<<std::endl;
	//if(onelinemode)std::cout<<"run lumi evt ";
	if(onelinemode)std::cout<<"run evt ";
	//if(onelinemode)std::cout<<"ak4PF40 ak4PF60 ak4PF80 ";
	//	if(onelinemode)std::cout<<"ak4PF40_PS ak4PF60_PS ak4PF80_PS ";
	//	if(onelinemode)std::cout<<"L1_SingleJet28_PS L1_SingleJet40_PS L1_SingleJet48_PS ";//<<std::endl;
	//if(onelinemode)std::cout<<"trgObj_maxTrgPt trgObj_maxTrgEta ";//<<std::endl;
	//if(onelinemode)std::cout<<"is40 is60 is80 ";//<<std::endl;
	if(onelinemode)std::cout<<"nJets jet1 p_T jet2 p_T ... jet(nJets) p_T ";
	if(onelinemode)std::cout<<"FinalPrescaleWeight ";

	std::cout<<std::endl;
      }
      
      if(!onelinemode){
	std::cout<<"run=="<<run_I<<std::endl;
	std::cout<<"lumi=="<<lumi_I<<std::endl;
	std::cout<<"evt=="<<evt_I<<std::endl;
	std::cout<<std::endl<<"HiEvtAnalyzer/EvtFilers:"<<std::endl;
	std::cout<<"pHBHENoiseFilter_I     =="<< pHBHENoiseFilter_I     <<std::endl;
	std::cout<<"pBeamScrapingFilter_I  =="<< pBeamScrapingFilter_I  <<std::endl;
	std::cout<<"pprimaryvertexFilter_I =="<< pprimaryvertexFilter_I <<std::endl;
	std::cout<<"vx_F =="<< vx_F <<std::endl;
	std::cout<<"vy_F =="<< vy_F <<std::endl;
	std::cout<<"vz_F =="<< vz_F <<std::endl;      
	std::cout<<std::endl<<"trackTree vtx info:"<<std::endl;
	std::cout<<"nVtx="<<nVtx_I << std::endl;
	std::cout<<"nTrkVtx="<<nTrkVtx_I[0] << std::endl;
	std::cout<<"xVtx="<<xVtx_F[0] << std::endl;
	std::cout<<"yVtx="<<yVtx_F[0] << std::endl;
	std::cout<<"zVtx="<<zVtx_F[0] << std::endl;
      }
      else if (onelinemode){
	std::cout<<run_I<<" ";
	//std::cout<<lumi_I<<" ";
	std::cout<<evt_I<<" ";
      }
      
      targetEvtFound=true;
      numevtsprinted++;
    }
    else{ 
      if(nEvt%10000==0){
	std::cout<<"evt=="<<evt_I<<std::endl;
	std::cout<<"run=="<<run_I<<std::endl;
	std::cout<<"lumi=="<<lumi_I<<std::endl;
	//std::cout<<"nvtxs="<<nvtx_I << std::endl;
      }
      continue;
    }





    //prescale/decision arrays, total prescale=L1_ps*HLT_ps
    bool PFtrgDec[N_HLTBits]   ={ (bool)PFJet40_I, (bool)PFJet60_I, (bool)PFJet80_I, (bool)PFJet100_I };
    int PFtrgPrescl[N_HLTBits] ={ PFJet40_p_I*jet40_l1s_ps_I , PFJet60_p_I*jet60_l1s_ps_I , 
				  PFJet80_p_I*jet80_l1s_ps_I , PFJet100_p_I*jet100_l1s_ps_I };
    
    bool CalotrgDec[N_HLTBits]   ={ (bool)CaloJet40_I, (bool)CaloJet60_I, (bool)CaloJet80_I, (bool)CaloJet100_I};
    int CalotrgPrescl[N_HLTBits] ={ CaloJet40_p_I*jet40_l1s_ps_I , CaloJet60_p_I*jet60_l1s_ps_I,
				    CaloJet80_p_I*jet80_l1s_ps_I , CaloJet100_p_I*jet100_l1s_ps_I };     
    
    if(targetEvtFound ) {// std::cout<<std::endl<<std::endl<<"target evt/run/lumi found"<<std::endl<<std::endl;
      if(!onelinemode){
	std::cout<<std::endl;
	std::cout<<"HLT PF 40/60/80     == "<<PFJet40_I << "/"<<PFJet60_I << "/"<<PFJet80_I  <<std::endl;
	std::cout<<"HLT PF 40/60/80 PSCL== "<< (PFJet40_p_I*jet40_l1s_ps_I) << "/"<< (PFJet60_p_I*jet60_l1s_ps_I)<< "/"<< (PFJet80_p_I*jet80_l1s_ps_I)  <<std::endl;
	
	std::cout<<std::endl;
	std::cout<<"HLT Calo 40/60/80     == "<<CaloJet40_I << "/"<<CaloJet60_I << "/"<<CaloJet80_I  <<std::endl;
	std::cout<<"HLT Calo 40/60/80 PSCL== "<< (CaloJet40_p_I*jet40_l1s_ps_I) << "/"<< (CaloJet60_p_I*jet60_l1s_ps_I)<< "/"<< (CaloJet80_p_I*jet80_l1s_ps_I)  <<std::endl;
      }
      else if(onelinemode){
	//std::cout<< PFJet40_I <<" ";
	//std::cout<< PFJet60_I <<" ";
	//std::cout<< PFJet80_I <<" ";
//	std::cout<< PFJet40_p_I <<" ";
//	std::cout<< PFJet60_p_I <<" ";
//	std::cout<< PFJet80_p_I <<" ";
//	std::cout<< jet40_l1s_ps_I <<" ";
//	std::cout<< jet60_l1s_ps_I <<" ";
//	std::cout<< jet80_l1s_ps_I <<" ";
	//assert(false);
      }
    }
    
    bool *trgDec=NULL; int *trgPscl=NULL;
    if(trgCombType=="Calo")   { 
      trgDec =CalotrgDec     ;   
      trgPscl=CalotrgPrescl ;     }
    else if(trgCombType=="PF"){ 
      trgDec=PFtrgDec   ;   
      trgPscl=PFtrgPrescl   ;     }
    else assert(false);//should never fire if this is working right
    
    
    // grab largest triggerPt among HLT paths' trigger object collectiosn
    double maxTrgPt=0.,maxTrgEta=0.; 
    unsigned int trgObj40_size=trgObjpt_40->size(), trgObj60_size=trgObjpt_60->size();
    unsigned int trgObj80_size=trgObjpt_80->size(), trgObj100_size=trgObjpt_100->size();
    
    
    if(trgDec[3]&&useHLT100)
      for(unsigned int itt=0; itt<trgObj100_size; ++itt)
    	if(trgObjpt_100->at(itt) > maxTrgPt) { 
	  maxTrgEta = trgObjeta_100->at(itt);
	  maxTrgPt = trgObjpt_100->at(itt); 
	} 
    
    if(trgDec[2])
      for(unsigned int itt=0; itt<trgObj80_size; ++itt)
    	if(trgObjpt_80->at(itt) > maxTrgPt)  { 
	  maxTrgEta = trgObjeta_80->at(itt);
	  maxTrgPt = trgObjpt_80->at(itt); 
	} 
    
    if(trgDec[1]) 
      for(unsigned int itt=0; itt<trgObj60_size; ++itt)
	if(trgObjpt_60->at(itt) > maxTrgPt)  { 
	  maxTrgEta = trgObjeta_60->at(itt);
	  maxTrgPt = trgObjpt_60->at(itt); 
	} 
    
    if(trgDec[0])  
      for(unsigned int itt=0; itt<trgObj40_size; ++itt)	
	if(trgObjpt_40->at(itt) > maxTrgPt)  { 
	  maxTrgEta = trgObjeta_40->at(itt);
	  maxTrgPt = trgObjpt_40->at(itt); 
	} 
    
    
    double trgPt=maxTrgPt; 
    double trgEta=maxTrgEta;
    
    // check trigger decisions for events + exclusivity between them, count events, assign prescale weight
    float weight_eS=0.;
    //float weight_eS = trigComb(trgDec, trgPscl, trgPt); // trig comb function replicates the procedure below
    
    bool is40  = false, is60  = false, is80  = false, is100 = false;
    if(useHLT100){
      if(      trgDec[3] && trgPt>=HLTthresh[3]               ) 
	{ is100 = true;  weight_eS=trgPscl[3]; }
      else if( trgDec[2] && trgPt>=HLTthresh[2]  && trgPt<HLTthresh[3] ) 
	{ is80  = true;  weight_eS=trgPscl[2]; }
      else if( trgDec[1] && trgPt>=HLTthresh[1]  && trgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=trgPscl[1]; }
      else if( trgDec[0] && trgPt>=HLTthresh[0]  && trgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=trgPscl[0]; }            
    }
    else {
      if( trgDec[2] && trgPt>=HLTthresh[2]  ) 
	{ is80  = true;  weight_eS=trgPscl[2]; }
      else if( trgDec[1] && trgPt>=HLTthresh[1]  && trgPt<HLTthresh[2]  ) 
	{ is60  = true;  weight_eS=trgPscl[1]; }
      else if( trgDec[0] && trgPt>=HLTthresh[0]  && trgPt<HLTthresh[1]  ) 
	{ is40  = true;  weight_eS=trgPscl[0]; }            
    }
    
     if(targetEvtFound ) {// std::cout<<std::endl<<std::endl<<"target evt/run/lumi found"<<std::endl<<std::endl;
       if(!onelinemode){
	 std::cout<<std::endl;
	 std::cout<<"trigComboType="<< trgCombType<<std::endl;
	 std::cout<<"HLT 40/60/80 thresholds== "<< HLTthresh[0] << "/"<< HLTthresh[1]<< "/"<< HLTthresh[2]  <<std::endl;
	 std::cout<<"TrgObj 40/60/80 pt/eta == "<< trgPt << "/"<< trgEta<< std::endl;
	 std::cout<<"is 40/60/80     == "<<is40 << "/"<<is60 << "/"<<is80  <<std::endl;
	 std::cout<<std::endl;
       }
       else if(onelinemode){
	 //	 std::cout << trgPt << " ";
	 //	 std::cout << trgEta << " ";	 
	 //std::cout << is40 <<" ";
	 //std::cout << is60 <<" ";
	 //std::cout << is80 <<" ";
	 //	 std::cout << weight_eS <<" ";
	 //assert(false);
       }
     }
     
     
     
    if     ( is100 )  { NEvents_100++ ; h_NEvents_jet100->Fill(1)  ; }
    else if( is80  )  { NEvents_80++  ;  h_NEvents_jet80->Fill(1)  ; }
    else if( is60  )  { NEvents_60++  ;  h_NEvents_jet60->Fill(1)  ; }
    else if( is40  )  { NEvents_40++  ;  h_NEvents_jet40->Fill(1)  ; }

    if(fillDataJetTrigQAHists){
      
      //make these trigger plots an array?
      if(trgDec[0]&&trgPt>=HLTthresh[0])   hpp_IncHLT40trgPt->Fill(  (float)trgPt, trgPscl[0] );
      if(trgDec[1]&&trgPt>=HLTthresh[1])   hpp_IncHLT60trgPt->Fill(  (float)trgPt, trgPscl[1] );
      if(trgDec[2]&&trgPt>=HLTthresh[2])   hpp_IncHLT80trgPt->Fill(  (float)trgPt, trgPscl[2] );
      if(trgDec[3]&&trgPt>=HLTthresh[3])   hpp_IncHLT100trgPt->Fill( (float)trgPt, trgPscl[3] );
      
      if(     is100 )  hpp_HLT100trgPt->Fill(   (float)trgPt , weight_eS );   
      else if(is80  )  hpp_HLT80trgPt->Fill(    (float)trgPt , weight_eS );
      else if(is60  )  hpp_HLT60trgPt->Fill(    (float)trgPt , weight_eS );
      else if(is40  )  hpp_HLT40trgPt->Fill(    (float)trgPt , weight_eS );
      
      if(true)         hpp_HLTCombtrgPt->Fill(  (float)trgPt , weight_eS );   

      
      if(trgDec[0]&&trgPt>=HLTthresh[0])   hpp_IncHLT40trgEta->Fill(  (float)trgEta, trgPscl[0] );
      if(trgDec[1]&&trgPt>=HLTthresh[1])   hpp_IncHLT60trgEta->Fill(  (float)trgEta, trgPscl[1] );
      if(trgDec[2]&&trgPt>=HLTthresh[2])   hpp_IncHLT80trgEta->Fill(  (float)trgEta, trgPscl[2] );
      if(trgDec[3]&&trgPt>=HLTthresh[3])   hpp_IncHLT100trgEta->Fill( (float)trgEta, trgPscl[3] );
      
      if(     is100 )  hpp_HLT100trgEta->Fill(   (float)trgEta , weight_eS );   
      else if(is80  )  hpp_HLT80trgEta->Fill(    (float)trgEta , weight_eS );
      else if(is60  )  hpp_HLT60trgEta->Fill(    (float)trgEta , weight_eS );
      else if(is40  )  hpp_HLT40trgEta->Fill(    (float)trgEta , weight_eS );
      
      if(true)         hpp_HLTCombtrgEta->Fill(  (float)trgEta , weight_eS );               
      
    }


    //if(weight_eS!=0.) h_NEvents_trigd->Fill(1); //this will be # of triggered events w/o weights pre vz/skim cuts
    //else continue;

    if(targetEvtFound&&!onelinemode) std::cout<<"prescale =="<<weight_eS<<std::endl;
    if(targetEvtFound&&!onelinemode) std::cout<<"nref     =="<<nref_I<<std::endl;
    
    if(targetEvtFound&&onelinemode) std::cout<<nref_I<<" ";//<<std::endl;
    
    // fill evt vz histo
    if(fillDataEvtQAHists){

      float vr_F=std::sqrt(vx_F*vx_F + vy_F*vy_F);
      hVr->Fill(vr_F, 1.0);   
      hWVr->Fill(vr_F, weight_eS);  
      if(weight_eS!=0.)hTrgVr_noW->Fill(vr_F,1.0);
      
      hVx->Fill(vx_F, 1.0);   
      hWVx->Fill(vx_F, weight_eS);  
      if(weight_eS!=0.)hTrgVx_noW->Fill(vx_F,1.0);

      hVy->Fill(vy_F, 1.0);   
      hWVy->Fill(vy_F, weight_eS);
      if(weight_eS!=0.)hTrgVy_noW->Fill(vy_F,1.0);

      hVz->Fill(vz_F, 1.0);   //this will be # of events post cuts, no trigger, no weights, post skim/vz Cut 
      hWVz->Fill(vz_F, weight_eS);    //this will be # of trigger events w/ weights, post skim/vz Cut
      if(weight_eS!=0.)hTrgVz_noW->Fill(vz_F,1.0); //this will be # of trigger events w/o weights post-skim/vz Cut

      hNref->Fill(nref_I,1.0);   //# of jets per event, no trigger, no weights, post cuts
      hWNref->Fill(nref_I,weight_eS); //# of jets per event, triggered, weighted, post cuts
    }

    //if(run_I==262252 && evt_I==113080452) std::cout<<std::endl<<std::endl<<"target evt/run/lumi found"<<std::endl<<std::endl;

    // JET LOOP 
    // for Aj/xj computation
    bool firstGoodJetFound=false, secondGoodJetFound=false; 
    float firstGoodJetPt=-1., secondGoodJetPt=-1.;
    float firstGoodJetPhi=-1., secondGoodJetPhi=-1.;
    
    bool firstGoodJetFound_wJetID=false, secondGoodJetFound_wJetID=false; 
    float firstGoodJetPt_wJetID=-1., secondGoodJetPt_wJetID=-1.;
    float firstGoodJetPhi_wJetID=-1., secondGoodJetPhi_wJetID=-1.;
    
    // for event counting + avoiding duplicate fills in dijet hists
    bool dijetHistsFilled=false, dijetHistsFilled_wJetID=false;
    bool hNEvts_withJets_Filled=false, hNEvts_withJets_kmatCut_Filled=false, hNEvts_withJets_JetIDCut_Filled=false; 

    int jetsPerEvent=nref_I, jetsPerEventJetID=nref_I;
    float evt_leadJetPt=-1.;//, evt_leadJetPt_wJetID=-1.;

    for(int jet = 0; jet<nref_I; ++jet){
      

      
      // event+jet counting
      h_NJets->Fill(1);
      h_WNJets->Fill(1,weight_eS);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      
      
      
      float recpt  = pt_F[jet];
      float receta = eta_F[jet];
      float absreceta = fabs(receta);
      float rawpt  = rawpt_F[jet];
      float recy   = y_F[jet];
      float recphi = phi_F[jet];

      
      //if(run_I==262252 && evt_I==113080452) std::cout<<std::endl<<std::endl<<"target evt/run/lumi found"<<std::endl<<std::endl;

      if(doResidualCorr){
	
	if(targetEvtFound&&!onelinemode)std::cout<<"jet #"<<jet<<" , eta = "<<receta<<" , y = "<<recy<<std::endl;
	if(targetEvtFound&&!onelinemode)std::cout<<"jet #"<<jet<<" , phi = "<<recphi<<std::endl;	
	if(targetEvtFound&&!onelinemode) std::cout<<"jet #"<<jet<< " , rawpt = "<<rawpt<<std::endl;	
	if(targetEvtFound&&!onelinemode) std::cout<<"jet #"<<jet<< " , recpt = "<<recpt<<std::endl;	

      	recpt=L2JES->getCorrectedPt(recpt,receta);  
	if(targetEvtFound&&!onelinemode) std::cout<<"jet #"<<jet<<" : recptL2 = "<<recpt<<std::endl;

	recpt=L3JES->getCorrectedPt(recpt);        
	if(targetEvtFound&&!onelinemode) std::cout<<"jet #"<<jet<<" : recptL3 = "<<recpt<<std::endl;
      }
      else{
	if(targetEvtFound&&!onelinemode) std::cout<<"jet #"<<jet<< " , y = "<<receta<<" , phi = "<<recphi<<std::endl;	
	if(targetEvtFound&&!onelinemode) std::cout<<"jet #"<<jet<< " , rawpt = "<<rawpt<<std::endl;	
	if(targetEvtFound&&!onelinemode) std::cout<<"jet #"<<jet<< " , recpt = "<<recpt<<std::endl;	
      }
      //std::cout<<std::endl;
      
      
      if(targetEvtFound&&onelinemode) std::cout<< recpt << " ";
      
      
      // kmatCuts      
      if( recpt <= jtPtCut ){ 
	jetsPerEvent--;	jetsPerEventJetID--;
	continue;}
      
      if( absreceta >= jtEtaCutHi ) {
	//jetsPerEvent--;	jetsPerEventJetID--;
	continue;}
      else if( absreceta < jtEtaCutLo ){ 
	//jetsPerEvent--;	jetsPerEventJetID--;
	continue;}

      if(absreceta >= 4.7 ) jetsPerEvent--; 

      // largest jet pt in each event
      if(recpt>evt_leadJetPt && absreceta<4.7)
	evt_leadJetPt=recpt;
      //if(recpt>evt_leadJetPt_wJetID && fillDataJetIDHists && passesJetID)
      //evt_leadJetPt_wJetID=recpt;
      
      h_NJets_kmatCut->Fill(1);
      h_WNJets_kmatCut->Fill(1,weight_eS);
      if(!hNEvts_withJets_kmatCut_Filled){
	h_NEvents_withJets_kmatCut->Fill(1);
	hNEvts_withJets_kmatCut_Filled=true;  }           
      
      
      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;
      float jetIDpt=recpt;
      
      
      
      // 13 TeV JetID criterion, loose or tight
      bool passesJetID=false;
      if(fillDataJetIDHists) 	{
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
      }
      
      
      //fill jetspectraRapHists w/ passing jetID criterion
      if( fillDataJetSpectraRapHists ) { 
	int theRapBin=-1;
	for(int rapbin=0;rapbin<nbins_abseta;++rapbin)
	  if( absetabins[rapbin]<=absreceta  && 		
	      absreceta<absetabins[rapbin+1]    	      ) {
	    
	    theRapBin=rapbin;
	    hJetSpectraRap[0][theRapBin]->Fill(recpt,weight_eS);  
	    if( passesJetID ) 
	      hJetSpectraRap[1][theRapBin]->Fill(recpt,weight_eS);    
	    
	    break;
	  } 
      }
      

      // trig plots
      if(fillDataJetTrigQAHists){
	
	hpp_TrgObjComb[0]->Fill(recpt, weight_eS); //kurts method
	
	//no exclusion between triggers
	if(trgDec[0]&&trgPt>=HLTthresh[0])  hpp_TrgObj40[0]->Fill(recpt, trgPscl[0]);
	if(trgDec[1]&&trgPt>=HLTthresh[1])  hpp_TrgObj60[0]->Fill(recpt, trgPscl[1]);
	if(trgDec[2]&&trgPt>=HLTthresh[2])  hpp_TrgObj80[0]->Fill(recpt, trgPscl[2]);
	if(trgDec[3]&&trgPt>=HLTthresh[3])  hpp_TrgObj100[0]->Fill(recpt, trgPscl[3]);  
	
	//with exclusion
	if(is40 )       hpp_ExcTrgObj40[0]->Fill( recpt , weight_eS);
	else if(is60 )  hpp_ExcTrgObj60[0]->Fill( recpt , weight_eS);
	else if(is80 )  hpp_ExcTrgObj80[0]->Fill( recpt , weight_eS);
	else if(is100)  hpp_ExcTrgObj100[0]->Fill( recpt , weight_eS); 
	
	//assert(false);
      }	
      
      
      
      // jetQA noJetID
      if(fillDataJetQAHists){
	int ind=0;
	//jets
	hJetQA[0][ind]->Fill(recpt, weight_eS); ind++;
	hJetQA[0][ind]->Fill(rawpt_F[jet], weight_eS); ind++;
	if(recpt>jetQAPtCut){//second jet pt cut
	  hJetQA[0][ind]->Fill(eta_F[jet], weight_eS); ind++;
	  hJetQA[0][ind]->Fill(phi_F[jet], weight_eS); ind++;

	  //tracks
	  hJetQA[0][ind]->Fill(trkN_I[jet], weight_eS); ind++;
	  hJetQA[0][ind]->Fill(trkSum_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(trkMax_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(trkHardN_I[jet], weight_eS); ind++;
	  hJetQA[0][ind]->Fill(trkHardSum_F[jet]/jetIDpt, weight_eS); ind++;

	  //PF photons
	  hJetQA[0][ind]->Fill(phN_I[jet], weight_eS); ind++;
	  hJetQA[0][ind]->Fill(phSum_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(phMax_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(phHardN_I[jet], weight_eS); ind++;
	  hJetQA[0][ind]->Fill(phHardSum_F[jet]/jetIDpt, weight_eS); ind++;

	  //PF charged hadrons
	  hJetQA[0][ind]->Fill(chN_I[jet], weight_eS); ind++;
	  hJetQA[0][ind]->Fill(chSum_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(chMax_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(chHardN_I[jet], weight_eS); ind++;
	  hJetQA[0][ind]->Fill(chHardSum_F[jet]/jetIDpt, weight_eS); ind++;

	  //PF neutral hadons
	  hJetQA[0][ind]->Fill(neN_I[jet], weight_eS); ind++;   
	  hJetQA[0][ind]->Fill(neSum_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(neMax_F[jet]/jetIDpt, weight_eS); ind++;

	  //PF electrons
	  hJetQA[0][ind]->Fill(eN_I[jet], weight_eS); ind++;   
	  hJetQA[0][ind]->Fill(eSum_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(eMax_F[jet]/jetIDpt, weight_eS); ind++;

	  //PF muons
	  hJetQA[0][ind]->Fill(muN_I[jet], weight_eS); ind++;   
	  hJetQA[0][ind]->Fill(muSum_F[jet]/jetIDpt, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(muMax_F[jet]/jetIDpt, weight_eS); ind++;

	  //PF particle sums
	  hJetQA[0][ind]->Fill(neuMult, weight_eS); ind++;	  
	  hJetQA[0][ind]->Fill(chMult, weight_eS); ind++;
	  hJetQA[0][ind]->Fill(numConst, weight_eS); ind++;

	  //hJetQA[0][ind]->Fill(hcalSum_F[jet]/jetIDpt, weight_eS); ind++;
	  //hJetQA[0][ind]->Fill(ecalSum_F[jet]/jetIDpt, weight_eS); ind++;	
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
	    hJetQA[0][ind]->Fill( A_j , weight_eS ); ind++; 
	    hJetQA[0][ind]->Fill( x_j , weight_eS ); ind++;       
	    hJetQA[0][ind]->Fill( dphi , weight_eS ); ind++;       
	    hJetQA[0][ind]->Fill( firstGoodJetPt , weight_eS ); ind++;       
	    hJetQA[0][ind]->Fill( secondGoodJetPt , weight_eS ); ind++; }
	}
      }
	
      
      // apply JetID
      if(fillDataJetIDHists){
	
	//jetsPerEventJetID=jetsPerEvent;
	if(!passesJetID) {
	  //jetsPerEventJetID--;
	  continue;}
	

	// trig plots
	if(fillDataJetTrigQAHists){

	  hpp_TrgObjComb[1]->Fill(recpt, weight_eS); //kurts method
	  
	  //no exclusion between triggers
	  if(trgDec[0]&&trgPt>=HLTthresh[0])   hpp_TrgObj40[1]->Fill(recpt, trgPscl[0]);
	  if(trgDec[1]&&trgPt>=HLTthresh[1])   hpp_TrgObj60[1]->Fill(recpt, trgPscl[1]);
	  if(trgDec[2]&&trgPt>=HLTthresh[2])   hpp_TrgObj80[1]->Fill(recpt, trgPscl[2]);
	  if(trgDec[3]&&trgPt>=HLTthresh[3])  hpp_TrgObj100[1]->Fill(recpt, trgPscl[3]);  
	  
	  //with exclusion
	  if(is40 )        hpp_ExcTrgObj40[1]->Fill( recpt , weight_eS);
	  else if(is60 )   hpp_ExcTrgObj60[1]->Fill( recpt , weight_eS);
	  else if(is80 )   hpp_ExcTrgObj80[1]->Fill( recpt , weight_eS);
	  else if(is100)   hpp_ExcTrgObj100[1]->Fill( recpt , weight_eS); 
	
	
	}

	
	// jetQA wJetID
	if(fillDataJetQAHists){
	  int ind=0;
	  
	  //jets
	  hJetQA[1][ind]->Fill(recpt, weight_eS); ind++;
	  hJetQA[1][ind]->Fill(rawpt_F[jet], weight_eS); ind++;
	  if(recpt>jetQAPtCut){ // second jetpt cut
	    hJetQA[1][ind]->Fill(eta_F[jet], weight_eS); ind++;
	    hJetQA[1][ind]->Fill(phi_F[jet], weight_eS); ind++;

	    //tracks
	    hJetQA[1][ind]->Fill(trkN_I[jet], weight_eS); ind++;
	    hJetQA[1][ind]->Fill(trkSum_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(trkMax_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(trkHardN_I[jet], weight_eS); ind++;
	    hJetQA[1][ind]->Fill(trkHardSum_F[jet]/jetIDpt, weight_eS); ind++;

	    //PF photons
	    hJetQA[1][ind]->Fill(phN_I[jet], weight_eS); ind++;
	    hJetQA[1][ind]->Fill(phSum_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(phMax_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(phHardN_I[jet], weight_eS); ind++;
	    hJetQA[1][ind]->Fill(phHardSum_F[jet]/jetIDpt, weight_eS); ind++;

	    //PF charged hadrons
	    hJetQA[1][ind]->Fill(chN_I[jet], weight_eS); ind++;
	    hJetQA[1][ind]->Fill(chSum_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(chMax_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(chHardN_I[jet], weight_eS); ind++;
	    hJetQA[1][ind]->Fill(chHardSum_F[jet]/jetIDpt, weight_eS); ind++;

	    //PF neutral hadons
	    hJetQA[1][ind]->Fill(neN_I[jet], weight_eS); ind++;   
	    hJetQA[1][ind]->Fill(neSum_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(neMax_F[jet]/jetIDpt, weight_eS); ind++;

	    //PF electrons
	    hJetQA[1][ind]->Fill(eN_I[jet], weight_eS); ind++;   
	    hJetQA[1][ind]->Fill(eSum_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(eMax_F[jet]/jetIDpt, weight_eS); ind++;

	    //PF muons
	    hJetQA[1][ind]->Fill(muN_I[jet], weight_eS); ind++;   
	    hJetQA[1][ind]->Fill(muSum_F[jet]/jetIDpt, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(muMax_F[jet]/jetIDpt, weight_eS); ind++;

	    //PF particle sums
	    hJetQA[1][ind]->Fill(neuMult, weight_eS); ind++;	  
	    hJetQA[1][ind]->Fill(chMult, weight_eS); ind++;
	    hJetQA[1][ind]->Fill(numConst, weight_eS); ind++;
	  }
	  
	  //looking for the first two good jets that meet the criteria specified
	  if(fillDataDijetHists){
	    if ( !firstGoodJetFound_wJetID ){
	      if(recpt>ldJetPtCut) { firstGoodJetFound_wJetID=true;
		firstGoodJetPt_wJetID =recpt; 
		firstGoodJetPhi_wJetID=recphi; } }
	    
	    else if ( !secondGoodJetFound_wJetID && 
		      firstGoodJetFound_wJetID      ) {
	      if (recpt>subldJetPtCut) { 
		float checkdPhi=deltaphi(firstGoodJetPhi_wJetID,recphi);
		float ptAve=(firstGoodJetPt_wJetID+recpt)/2.;
		if( checkdPhi>dPhiCut && ptAve>ptAveCut) {
		  if(recpt>firstGoodJetPt_wJetID){secondGoodJetFound_wJetID=true;
		    std::cout<<std::endl<<
		      "WARNING: picked wrong jet for lead jet! Swapping..."<<std::endl<<std::endl;
		    secondGoodJetPt_wJetID  = firstGoodJetPt_wJetID;          
		    secondGoodJetPhi_wJetID = firstGoodJetPhi_wJetID;
		    firstGoodJetPt_wJetID   = recpt;          
		    firstGoodJetPhi_wJetID  = recphi; }
		  else { secondGoodJetFound_wJetID=true;
		    secondGoodJetPt_wJetID  = recpt;
		    secondGoodJetPhi_wJetID = recphi; }}} }
	    
	    if (!dijetHistsFilled_wJetID &&
		firstGoodJetFound_wJetID && 
		secondGoodJetFound_wJetID ) { dijetHistsFilled_wJetID=true;
	      float A_j=(firstGoodJetPt_wJetID-secondGoodJetPt_wJetID)
		/(firstGoodJetPt_wJetID+secondGoodJetPt_wJetID);
	      float x_j=(secondGoodJetPt_wJetID/firstGoodJetPt_wJetID); 	
	      float dphi=deltaphi(firstGoodJetPhi_wJetID,secondGoodJetPhi_wJetID);
	      int ind=dijtInd;
	      hJetQA[1][ind]->Fill( A_j , weight_eS ); ind++; 
	      hJetQA[1][ind]->Fill( x_j , weight_eS ); ind++;       
	      hJetQA[1][ind]->Fill( dphi , weight_eS ); ind++;       
	      hJetQA[1][ind]->Fill( firstGoodJetPt_wJetID , weight_eS ); ind++;       
	      hJetQA[1][ind]->Fill( secondGoodJetPt_wJetID , weight_eS ); ind++;       } 
	  }
	}//end fillDataJetQAhists
	
	// jet/event counts
	h_NJets_JetIDCut->Fill(1);
	h_WNJets_JetIDCut->Fill(1,weight_eS);
	if(!hNEvts_withJets_JetIDCut_Filled){
	  h_NEvents_withJets_JetIDCut->Fill(1);
	  hNEvts_withJets_JetIDCut_Filled=true;   }
	
      }//end fillDataAJetIDhists
      
      
    }//end jet loop
    std::cout << weight_eS <<" ";
    //if(onelinemode){ 
    //std::cout<<endl;
    //std::cout<<endl;
    //assert(false);
    //}
    if(fillDataEvtQAHists){
      if(fillDataJetIDHists){
	//hLeadJetPt_wJetID->Fill(evt_leadJetPt_wJetID ,weight_eS);
	hjetsPEvt_wJetID ->Fill(jetsPerEventJetID,1.0);   
	hWjetsPEvt_wJetID->Fill(jetsPerEventJetID,weight_eS); }   
      if(evt_leadJetPt>0.) hLeadJetPt->Fill(evt_leadJetPt ,weight_eS);
      hjetsPEvt ->Fill(jetsPerEvent,1.0);   
      hWjetsPEvt->Fill(jetsPerEvent,weight_eS); 
    }
    std::cout<<std::endl;
  }//end event loop


  if(!onelinemode){
    std::cout<<std::endl;
    
    std::cout<<std::endl<<"/// Job Event-Loop Summary ///"<<std::endl<<std::endl;
    std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events skipped from those file(s) = " <<h_NEvents_skipped->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl<<std::endl;
    
    std::cout<<"Total Num of Events passing jet40 = "  << h_NEvents_jet40->GetEntries()  << std::endl;
    std::cout<<"Total Num of Events passing jet60 = "  << h_NEvents_jet60->GetEntries()  << std::endl;
    std::cout<<"Total Num of Events passing jet80 = "  << h_NEvents_jet80->GetEntries()  << std::endl;
    std::cout<<"Total Num of Events passing jet100 = " << h_NEvents_jet100->GetEntries() << std::endl;
    std::cout<<"Total Num of Events passing a trigger = " << h_NEvents_trigd->GetEntries() << std::endl << std::endl;
    
    std::cout<<"Total Num of Events passing skimCuts  = " << h_NEvents_skimCut->GetEntries()<<std::endl;
    std::cout<<"Total Num of Events passing vzCuts    = " << h_NEvents_vzCut->GetEntries()<<std::endl<<std::endl;
    
    std::cout<<"Total Num of good Events w/ jets = " << h_NEvents_withJets->GetEntries()<<std::endl; //note, approximate, some events make it through with nref=0 for some reason
    std::cout<<"Total Num of good Events, w/ jets, post kmatCuts = " << h_NEvents_withJets_kmatCut->GetEntries() << std::endl; 
    if(fillDataJetIDHists) 
      std::cout<<"Total Num of good Events, w/ jets, post kmatCut AND JetIDCut = " << h_NEvents_withJets_JetIDCut->GetEntries()<<std::endl;        
    
    std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
    std::cout<<"Total Num of Jets read from good events                          = " <<    h_NJets->GetEntries()<<std::endl;  
    std::cout<<"Total Num of Jets read from good events post kinCut              = " <<    h_NJets_kmatCut->GetEntries()<<std::endl;
    if(fillDataJetIDHists) 
      std::cout<<"Total Num of Jets read from good events post kinCut AND JetIDCut = " << h_NJets_JetIDCut->GetEntries()<<std::endl;
    std::cout<<std::endl;
    
    std::cout<<"writing output file... ";
  }
  fout->Write(); 
  
  trgObjpt_40->clear();
  trgObjpt_60->clear();
  trgObjpt_80->clear();
  trgObjpt_100->clear();
  
  trgObjeta_40->clear();
  trgObjeta_60->clear();
  trgObjeta_80->clear();
  trgObjeta_100->clear();
  
  if(!onelinemode)std::cout<<std::endl<<"readForests_ppData_findEvt finished."<<std::endl;  
  timer.Stop();
  if(!onelinemode)std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  if(!onelinemode)std::cout<<"Real tyime (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  if(!onelinemode)std::cout<<"output file written is  "<<outfile<<std::endl<<std::endl;
  return 0;
}    // end readForests_ppData_findEvt


////// main //////
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=8 && argc!=1 && argc!=10){
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppData_findEvt.exe ";
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
  if(argc==1) rStatus = readForests_ppData_findEvt();
  else if (argc==8){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    rStatus = readForests_ppData_findEvt( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
				  outputFileName);
  }
  else if (argc==10){//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
    rStatus = readForests_ppData_findEvt( inputFileList, startfile, endfile, 
				  jetRadius, jetType, debug,
					   outputFileName, jtEtaCutLo, jtEtaCutHi);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}










//TH1F *hpp_HLT40InEff =NULL;
//TH1F *hpp_HLT60InEff =NULL;
//TH1F *hpp_HLT80InEff =NULL;
//TH1F *hpp_HLT100InEff=NULL;
//TH1F *hpp_HLTAllInEff=NULL;      

////ineff
//hpp_HLT40InEff  = new TH1F( "isHLT40_inEff"  , "isHLT40  inEff Check", 1000, 0,1000 ); 
//hpp_HLT60InEff  = new TH1F( "isHLT60_inEff"  , "isHLT60  inEff Check", 1000, 0,1000 ); 
//hpp_HLT80InEff  = new TH1F( "isHLT80_inEff"  , "isHLT80  inEff Check", 1000, 0,1000 ); 
//hpp_HLT100InEff = new TH1F( "isHLT100_inEff" , "isHLT100 inEff Check", 1000, 0,1000 );       
//
//hpp_HLTAllInEff = new TH1F( "isHLTAll_inEff" , "isHLTAll inEff Check", 1000, 0,1000 );                 

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
  
