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
//const int N_vars_xAx_Titles = sizeof(var_xAx_Titles)/sizeof(std::string);


