#include <cstdlib>
#include <cstring>

const std::string fNLOJetsSpectra_dir="/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/fNLOJetsSpectra";

//std::string fNLO_filename_start="fnl5020_LO2_R0";
//std::string fNLO_filename_mid="Jets_modify_";

//raghav's files
const std::string fNLOFile_R04_CT10nlo  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_LO2_R04Jets_modify_CT10nlo_HS.root";
//const std::string fNLOFile_R04_CT14nlo  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_LO2_R04Jets_modify_CT14nlo_HS.root";
const std::string fNLOFile_R04_CT14nlo  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQpt_v4.root";//mur=muf=jetpt, for each jet
const std::string fNLOFile_R04_CT14nlo2  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQpt1_v3.root";//mur=muf=leadjet pt, for all jets in an event
const std::string fNLOFile_R04_CT14nlo3  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQHTp_v4.root";//mur=muf=leadjet pt, for all jets in an event
const std::string fNLOFile_R04_HERAPDF  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_LO2_R04Jets_modify_HERAPDF15NLO_ALPHAS_HS.root";
const std::string fNLOFile_R04_MMHTnlo  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_LO2_R04Jets_modify_MMHT2014nlo68cl_HS.root";
const std::string fNLOFile_R04_NNPDFnnlo        = fNLOJetsSpectra_dir+"/R04/fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC.root";

//joao's files
//const std::string fNLOFile_R04_CT14nnlo  	  = fNLOJetsSpectra_dir+"/R04/fnl5020_R04Jets_ybins_CT14nnlo_murEQmufEQpt_v2.root";


//const std::string NPCorr_filename               = "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV.root";
//const std::string NPCorr_filename               = "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV_fromJohn.root";
const std::string NPCorr_filename               = "/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/doAnalysis/doUnfolding/smearTheory/NPCorr5TeV/NLOpNP_InclusiveJets5TeV_fromJohn_refit.root";
