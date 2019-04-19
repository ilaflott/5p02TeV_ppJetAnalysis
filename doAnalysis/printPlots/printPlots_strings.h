// C++, C, etc.  
#include <cstdlib>  
#include <cstdio>   
#include <cmath>   

const std::string CMSSW_BASE=
  "  /home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/";
// "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/";
const std::string SCRATCH_BASE=
  "/cms/se/store/user/ilaflott/jobOutput/T2_US_MIT_SCRATCH/5p02TeV_ppJetAnalysis/readForests/";
//  "/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/";
const std::string HOME_BASE=
  "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/printPlots/";
  //  "/Users/ilaflott/Working/CERNBox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/";
//"/Users/ilaflott/Working/Dropbox/localAnalysis/src/5p02TeV_ppJetAnalysis-master/doAnalysis/printPlots/";

//const std::string inputDir=HOME_BASE+"outputCondor/";
//const std::string inputDir=SCRATCH_BASE+"outputCondor/";
const std::string inputDir=SCRATCH_BASE;

const std::string outputDir=HOME_BASE+"output/";

// misc text ------------------------
const std::string crossSectionAxTitle  ="#sigma (nb)";
//const std::string ddcrossSectionAxTitle="#frac{d^{2}#sigma}{dp_{T}d#||{y}} (nb/GeV)";
const std::string ddcrossSectionAxTitle="#frac{d^{2}#sigma}{dp_{T}dy} (nb/GeV)";
const std::string dcrossSectionAxTitle="#frac{d#sigma^{}}{dp_{T}} (nb/GeV)";
const std::string AUAxTitle  ="A.U.";
const std::string ratioTitle ="MC/Data";
//const std::string jetIDratioTitle ="JetID/no JetID";
//const std::string jetIDratioTitle ="1 - (w/ Jet ID)/(w/o JetID)";
const std::string jetIDratioTitle ="% Removed w/ Jet ID";

const std::string intLumiString="L_{int} = 27.4  ^{}pb^{ -1}";
const std::string sqrts2k15ppString="2015 pp, ^{}#sqrt{s} = 5.02 TeV";
const std::string Py8TuneString="PYTHIA8 Tune CUETP8M1";

const std::string promptRecoPDString="pp 2015 promptReco, Jet80+LowerJets";
const std::string PDStatsString="L_{int} = 27.4 pb^{-1}, #sqrt{s} = 5.02 TeV";
//const std::string PDStatsString_2="27.4 pb^{-1}, #sqrt{s} = 5.02 TeV, #it{CMS Work In Progress}";
const std::string PDStatsString_2="#bf{CMS} #it{Work In Progress}                                             27.4 pb^{-1} (5.02 TeV)";

// variable names for QA Plots
const std::string var[]={
  //jets, 0-6, 7 vars tot
  "jtpt"  , "rawpt" ,
  //"jtpt_forRes", "jtpt_L2Res", "jtpt_L3Res",
  "jteta" , "jtphi" ,

  //jet constituents, 7-33, 27 vars tot
  "trkN", "trkSum", "trkMax", //trk
  "trkHardN", "trkHardSum",

  "phN", "phSum", "phMax",    //ph
  "phHardN", "phHardSum",

  "chN", "chSum", "chMax",    //ch
  "chHardN", "chHardSum",

  "neN" ,  "neSum" , "neMax" , //ne
  "eN"  ,  "eSum"  , "eMax"  ,    //e
  "muN" ,  "muSum" , "muMax" , //mu

  "neuMult", "chMult", "numConst",

  // dijets, 34-38=5 vars
  "Aj" , "xj" , "dphi",    //JERS might not like me
  "leadJetPt", "subleadJetPt"  
};
const int N_vars=sizeof(var)/sizeof(std::string);
const int jetQA_varStart=0;
//const int jetQA_Nvars=7;
const int jetQA_Nvars=4;
const int jetConstits_varStart=jetQA_varStart+jetQA_Nvars;
const int jetConstits_Nvars=27;
const int dijet_varStart=jetConstits_varStart+jetConstits_Nvars;

// variable names, string array
const std::string tupelvar[]={
  "Pt",
  "Eta",
  "Phi",
  "E",
  "Id",
  "RawPt",
  "RawE",
  "HfHadE",
  "HfEmE",
  "ChHadFrac",
  "NeutralHadAnHfFrac",
  "ChEmFrac",
  "NeutralEmFrac",
  "ChMult",
  "ConstCnt",
  "NeuMult"
};
const int N_tupelvars=sizeof(tupelvar)/sizeof(std::string);

const std::string tupelvar_xAx_Titles[N_tupelvars]={
  "Jet p_{T}^{RECO} [GeV]",//"Pt",
  "Jet y^{RECO}",//"Eta",
  "Jet #phi^{RECO}",//"Phi",
  "Jet E^{RECO} [GeV]",//"E",
  "SMP Jet ID [0/1/2/3]=[!Pass/PassTight/PassMed/PassLoose]",//"Id",
  "Jet p_{T}^{Raw} [GeV]",//"RawPt",
  "Jet E^{RAW} [GeV]",//"RawE",
  "HF Calo Hadronic Energy [GeV]",//"HfHadE",
  "HF Calo Electromagnetic Energy[GeV]",//"HfEmE",
  "Charged Hadron Energy/Jet E^{RAW} (PF Charged Hadrons)",//"ChHadFrac",
  "Neutral Hadron Energy/Jet E^{RAW} (PF Neutral Hadrons)",//"NeutralHadAnHfFrac",
  "Charged EM Energy/Jet E^{RAW} (PF Electrons)",//"ChEmFrac",
  "Neutral EM Energy/Jet E^{RAW} (PF Photons)",//"NeutralEmFrac",
  "Charged Multiplicity (PF Charged Hadrons + Electrons + Muons)",//"ChMult",
  "Constituent Count (All PF Particles)",//"ConstCnt",
  "Neutral Multiplicity (PF Neutral Hadrons + Photons)"//"NeuMult"
};

const std::string tupelvar_incjetanaequiv[]={
"jtpt",//  "Pt",
"jteta",//  "Eta",
"jtphi",//  "Phi",
"NULL",//  "E",
"NULL",//  "Id",
"rawpt",//  "RawPt",
"NULL",//  "RawE",
"NULL",//  "HfHadE",
"NULL",//  "HfEmE",
"chSum",//  "ChHadFrac",
"neSum",//  "NeutralHadAnHfFrac",
"eSum",//  "ChEmFrac",
"phSum",//  "NeutralEmFrac",
"chMult",//  "ChMult",
"numConst",//  "ConstCnt",
"neuMult"//  "NeuMult"
};


const std::string var_xAx_Titles[] = {   
  "Jet p_{T}^{RECO} (GeV)" ,  "Jet p_{T}^{RAW} (GeV)",  
//  "Jet p_{T}^{RECO} (GeV) no Residuals" , 
//  "Jet p_{T}^{RECO} (GeV) L2 Residuals" , 
//  "Jet p_{T}^{RECO} (GeV) L2+L3 Residuals"  ,
  "#eta_{Jet}", "#phi_{Jet}", //"y_{jet}" ,
  
  "N^{Trks}",   "#Sigma p_{T}^{Trk}/p_{T}^{RECO}", "Max p_{T}^{Trk}/p_{T}^{RECO}", 
  "N_{Hard}^{Trks}",   "#Sigma p_{T}^{Hard Trk}/p_{T}^{RECO}",
  
  "phN", "phSum/p_{T}^{RECO}", "phMax/p_{T}^{RECO}", 
  "phHardN", "phHardSum/p_{T}^{RECO}",
  
  "chN", "chSum/p_{T}^{RECO}", "chMax/p_{T}^{RECO}", 
  "chHardN", "chHardSum/p_{T}^{RECO}",
  
  "neN", "neSum/p_{T}^{RECO}",  "neMax/p_{T}^{RECO}", 
  "eN",   "eSum/p_{T}^{RECO}",   "eMax/p_{T}^{RECO}", 
  "muN", "muSum/p_{T}^{RECO}",  "muMax/p_{T}^{RECO}",
  
  "N_{neu}^{Jet}",   "N_{chg}^{Jet}",   "N_{sum}^{Jet}",
  
  "A_{j}",     "x_{j}" , "#Delta #phi", 
  "Jet p_{T,1}^{RECO} (GeV)", "Jet p_{T,2}^{RECO} (GeV)"
};

// gen variable names for Th2 profiles//MCEff Ploits
const std::string genVars[] = {   
  "pt","eta","phi" //,"drjt"
};
const int N_genVars = sizeof(genVars)/sizeof(std::string);

const std::string genVars_ptrat[] = {   
  "eta","phi" ,"drjt"
};
const int N_genVars_ptrat = sizeof(genVars_ptrat)/sizeof(std::string);


// trigs
const std::string HLTName[] = {
  "HLTComb","HLT40","HLT60","HLT80","HLT100", "HLTMB"
};
// trigs
const std::string HLTCaloName_Leg[] = {
  "HLT ak4CaloJet Combo","HLT ak4CaloJet40","HLT ak4CaloJet60","HLT ak4CaloJet80","HLT ak4CaloJet100", "HLT MinimumBias HF1OR"
};
// trigs
const std::string HLTPFName_Leg[] = {
  "HLT ak4PFJet Combo","HLT ak4PFJet40","HLT ak4PFJet60","HLT ak4PFJet80","HLT ak4PFJet100", "HLT MinimumBias HF1OR"
};
const int N_trigs=sizeof(HLTName)/sizeof(std::string);

// string arrays
//L1
const std::string L1BitStrings[]={
  "L1_SingleJet28_BptxAND",
  "L1_SingleJet40_BptxAND",
  "L1_SingleJet48_BptxAND",
  "L1_SingleJet52_BptxAND"
};
const int N_L1Bits=sizeof(L1BitStrings)/sizeof(std::string);

//HLT
const std::string HLTBitStrings[N_L1Bits]={
  "HLT_AK4CaloJet40_Eta5p1",
  "HLT_AK4CaloJet60_Eta5p1",
  "HLT_AK4CaloJet80_Eta5p1",
  "HLT_AK4CaloJet100_Eta5p1"    
};
const int N_HLTBits=sizeof(HLTBitStrings)/sizeof(std::string);

//// combos
const std::string etaWidth = "20_eta_20";//why do i still have this

bool str_replace(std::string& str, const std::string& from, const std::string& to) {
//bool str_replace(std::string str, const std::string from, const std::string to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}

bool findsubstr(std::string& str, const std::string& substr) {
  //bool str_replace(std::string str, const std::string from, const std::string to) {
  bool substrfound=false;
  size_t substr_pos= str.find(substr);
  //size_t start_pos = str.find(from);
  if(substr_pos != std::string::npos)
    substrfound=true;  
  return substrfound;
}



















// string arrays, etc.
const char *cdir[]  = {"PP"};
const char *ccent[] = {"PP"};
const int ncen=1;

const std::string algo ="ak";// strings of jet algo
//const std::string algo[] ={"ak","akPu"};// strings of jet algo

const std::string jetType="PF";// strings of jet type
//const std::string jetType[]={"PF","Calo","Cs"};// strings of jet type

const std::string srad[] ={"4"};// strings of jet radii
const int Nrad =sizeof(srad)/sizeof(std::string);//wasn't in original code, 10 arbitrarily chosen

const std::string PFCandType[] = {"unknown",
				  "chargedHadron",
				  "electron",  "muon", "photon",
				  "neutralHadron",
				  "HadEnergyinHF", "EMEnergyinHF"};  
const int PFType = sizeof(PFCandType)/sizeof(std::string);


//// tree names+directories
//const std::string ppDataTrees[]={ 
//  "GARBAGE ENTRY" , //use jet ana of choice later
//  "hiEvtAnalyzer/HiTree" ,
//  "skimanalysis/HltTree" ,
//  "hltanalysis/HltTree" ,
//  "hltobject/"+HLTBitStrings[0]+"_v" , 
//  "hltobject/"+HLTBitStrings[1]+"_v" , 
//  "hltobject/"+HLTBitStrings[2]+"_v" , 
//  "hltobject/"+HLTBitStrings[3]+"_v"   
//}; 
//const int N_trees_ppData = sizeof(treeNames)/sizeof(std::string);



//// other useful things
//const std::string HLThNames[]={"HLT40_","HLT60_","HLT80_","HLT100_","HLTComb_","TrgCombTest_"};
//const std::string HLThTitle[]={"HLT40 ","HLT60 ","HLT80 ","HLT100 ","HLTComb ","Kurt Meth. HLTComb "};
//const int N_HLThNames=sizeof(HLThNames)/sizeof(std::string);

void makeFullJetTypeString(std::string* fullJetType, std::string* radius, std::string input_ppData_condorDir){
  bool funcDebug=false;
  if(funcDebug)std::cout<<"in makeFullJetTypeString"<<std::endl;
  
  std::size_t radPos=input_ppData_condorDir.find("_ak")+3;  
  const std::string radiusInt= input_ppData_condorDir.substr( radPos,1 );
  
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=input_ppData_condorDir.find("Jets");;
  
  const std::string jetType=input_ppData_condorDir.substr( jetTypePos,(jetsPos-jetTypePos) );
  if(funcDebug)std::cout<<"jetType string is = "<<jetType<<std::endl;
  
  (*radius)="R"+radiusInt+"_";
  if(funcDebug)std::cout<<"radius string is = "<<*radius<<std::endl;
  
  (*fullJetType)=(std::string)("ak"+radiusInt+jetType);
  if(funcDebug)std::cout<<"fullJetType string is = "<<*fullJetType<<std::endl;

  return;
}
//put this in diff header TO DO
void makeJetCutStrings( std::string* jetCutString, 
		        std::string* jetEtaCutString, 
		        std::string* jtptQACut_str, 
			TFile* fin){
  bool funcDebug=false;
  if(funcDebug)std::cout<<"in makeJetCutStrings"<<std::endl;
  
  //grab cut hists + relevant titles+#'s
  TH1F* jtptCut_h= (TH1F*)fin->Get( "hJetPtCut" );
  std::string jtptCut_str = std::to_string( (int) jtptCut_h->GetMean() );
  if(funcDebug)std::cout<<"jtptCut_str = "<<jtptCut_str<<std::endl;  
  
  TH1F* jtptQACut_h= (TH1F*)fin->Get( "hJetQAPtCut" );
  
  TH1F* jtetaLoCut_h= (TH1F*)fin->Get( "hJetEtaCutLo" );
  float jtetaLoCut_F=jtetaLoCut_h->GetMean();
  std::stringstream etaLo; etaLo.precision(1);
  etaLo << std::fixed << jtetaLoCut_F;
  
  TH1F* jtetaHiCut_h= (TH1F*)fin->Get( "hJetEtaCutHi" );
  float jtetaHiCut_F=jtetaHiCut_h->GetMean();
  std::stringstream etaHi; etaHi.precision(1);
  etaHi << std::fixed << jtetaHiCut_F;
  
  //form the strings as desired
  (*jetCutString)=(std::string)("p_{T}>"+jtptCut_str+" GeV");//, "+jtetaLoCut_str+"<|#eta|<"+jtetaHiCut_str;
  if(funcDebug)  std::cout<<"jetCutString="<<*jetCutString<<std::endl;  
  
  (*jtptQACut_str) =(std::string)( std::to_string( (int) jtptQACut_h->GetMean() ));
  if(funcDebug)  std::cout<<"jtptQACut_str = "<<*jtptQACut_str<<std::endl;  
  
  std::string jtetaLoCut_str = etaLo.str();
  if(funcDebug)std::cout<<"jtetaLoCut_str = "<<jtetaLoCut_str<<std::endl;    
  
  std::string jtetaHiCut_str = etaHi.str();
  if(funcDebug)std::cout<<"jtetaHiCut_str = "<<jtetaHiCut_str<<std::endl;  
  
  (*jetEtaCutString)=(std::string)(jtetaLoCut_str+"<|y|<"+jtetaHiCut_str);
  if(funcDebug)std::cout<<"jetEtaCutString="<<*jetEtaCutString<<std::endl;  
  
  return;
}
