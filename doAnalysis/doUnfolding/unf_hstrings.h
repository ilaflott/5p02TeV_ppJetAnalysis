#include <string>
#include <vector>
const std::string yaxtitle_xsec="#frac{d^{2}#sigma}{d#eta dp_{T}} [nb/GeV]";


const int NPLTS=4;
//const std::string UNF_DIR="/Users/ilaflott/Working/CERNbox/L2L3recoJetID_HLTPF_w_Py8LOMC_NLOMC_NLONPMC/";
const std::string UNF_DIR="/cms/heavyion/ilaflott/temp/4scp/L2L3recoJetID_HLTPF_w_Py8LOMC_NLOMC_NLONPMC/";

std::vector<std::string> HNAMES=
  {"MC Meas.", "MC Prior",
   "Data Meas.", "Data Unf."
  };

//BAYESIAN FILES
std::string     FILE_PY8_bayes = UNF_DIR+"ak4PF_wjtID_anabins_unf_HPtJetTrig_w_Py8_Bayes_07.02.18_0.0eta2.0_gendr0p1_recoJetID.root";
//std::string PY8_bayes_NAMES[NPLTS]=
std::vector<std::string> PY8_bayes_NAMES=
  {"hpp_reco_wJetID_R4_20_eta_20_clone_rebins", "hpp_gen_wJetID_R4_20_eta_20_clone_rebins", //MC MEAS AND PRIOR
     "hJetQA_1wJetID_jtpt_divBylumietabin_clone_rebins","ppData_BayesUnf_Spectra"//DATA MEAS AND DATA UNF
  };

std::string     FILE_NLO_bayes = UNF_DIR+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_Bayes_07.02.18_NNPDF30nnlo_NLOMC.root";
//std::string NLO_bayes_NAMES[NPLTS]=
std::vector<std::string> NLO_bayes_NAMES=
  {"smeared_rnd_ynew_clone_rebins", "theory_rnd_ynew_clone_rebins",//MC MEAS AND PRIOR
   "hJetQA_1wJetID_jtpt_divBylumietabin_clone_rebins","ppData_BayesUnf_Spectra"//DATA MEAS AND DATA UNF
};

std::string  FILE_NLOwNP_bayes = UNF_DIR+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_Bayes_07.02.18_NNPDF30nnlo_NLOMC_wNP.root";
//std::string NLOwNP_bayes_NAMES[NPLTS]=
std::vector<std::string> NLOwNP_bayes_NAMES=
  {"smeared_rnd_NPynew_clone_rebins", "theory_rnd_NPynew_clone_rebins",//MC MEAS AND PRIOR
   "hJetQA_1wJetID_jtpt_divBylumietabin_clone_rebins","ppData_BayesUnf_Spectra"//DATA MEAS AND DATA UNF
};

//SVD FILES            
std::string     FILE_PY8_SVD = UNF_DIR+"ak4PF_wjtID_anabins_unf_HPtJetTrig_w_Py8_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID.root";
//std::string PY8_SVD_NAMES[NPLTS]=
std::vector<std::string> PY8_SVD_NAMES=
  {"hpp_reco_wJetID_R4_20_eta_20_clone_rebins", "hpp_gen_wJetID_R4_20_eta_20_clone_rebins",//MC MEAS AND PRIOR
   "hJetQA_1wJetID_jtpt_divBylumi_clone_rebins", "hunf_svd_kReg8"//DATA MEAS AND DATA UNF
};

std::string     FILE_NLO_SVD = UNF_DIR+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID_NNPDF30nnlo_NLOMC.root";
//std::string NLO_SVD_NAMES[NPLTS]=
std::vector<std::string> NLO_SVD_NAMES=
  {"smeared_rnd_ynew_clone_rebins", "theory_rnd_ynew_clone_rebins",//MC MEAS AND PRIOR
   "hJetQA_1wJetID_jtpt_divBylumi_clone_rebins", "hunf_svd_kReg16"//DATA MEAS AND DATA UNF
};

std::string  FILE_NLOwNP_SVD = UNF_DIR+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID_NNPDF30nnlo_NLOMC_wNP.root";
//std::string NLOwNP_SVD_NAMES[NPLTS]=
std::vector<std::string> NLOwNP_SVD_NAMES=
  {"smeared_rnd_NPynew_clone_rebins", "theory_rnd_NPynew_clone_rebins",//MC MEAS AND PRIOR
   "hJetQA_1wJetID_jtpt_divBylumi_clone_rebins","hunf_svd_kReg16"//DATA MEAS AND DATA UNF
};
