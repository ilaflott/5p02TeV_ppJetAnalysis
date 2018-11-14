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
  "jtpt_forRes", "jtpt_L2Res", "jtpt_L3Res",
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
const int jetQA_Nvars=7;
const int jetConstits_varStart=jetQA_varStart+jetQA_Nvars;
const int jetConstits_Nvars=27;
const int dijet_varStart=jetConstits_varStart+jetConstits_Nvars;

// variable names, string array
const std::string var_xAx_Titles[] = {   
  "Jet p_{T}^{RECO} (GeV)" ,  "Jet p_{T}^{RAW} (GeV)",  
  "Jet p_{T}^{RECO} (GeV) no Residuals" , 
  "Jet p_{T}^{RECO} (GeV) L2 Residuals" , 
  "Jet p_{T}^{RECO} (GeV) L2+L3 Residuals"  ,
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
  "HLTComb","HLT40","HLT60","HLT80","HLT100"
};
// trigs
const std::string HLTCaloName_Leg[] = {
  "HLTak4CaloJetsCombo","HLTak4CaloJets40","HLTak4CaloJets60","HLTak4CaloJets80","HLTak4CaloJets100"
};
// trigs
const std::string HLTPFName_Leg[] = {
  "HLTPFJets Combo","HLTak4PFJets40","HLTak4PFJets60","HLTak4PFJets80","HLTak4PFJets100"
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
const std::string etaWidth = "20_eta_20";

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

