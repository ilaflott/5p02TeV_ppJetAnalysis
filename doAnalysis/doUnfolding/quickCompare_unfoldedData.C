#include <cstdlib>
#include <iostream>
#include <cassert>
#include <string>

#include "TH2.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"

const int NPLTS=3;
const std::string unf_dir="/Users/ilaflott/Working/CERNbox/L2L3recoJetID_HLTPF_w_Py8LOMC_NLOMC_NLONPMC/";

//BAYESIAN FILES
std::string     file_PY8_bayes = unf_dir+"ak4PF_wjtID_anabins_unf_HPtJetTrig_w_Py8_Bayes_07.02.18_0.0eta2.0_gendr0p1_recoJetID.root";
std::string PY8_bayes_names[NPLTS]={
  "hpp_gen_wJetID_R4_20_eta_20_clone_rebins",
  "hJetQA_1wJetID_jtpt_divBylumietabin_clone_rebins","ppData_BayesUnf_Spectra"};
//std::string PY8_bayes_hmat_name         ="", PY8_bayes_hrec_sameside_name="", PY8_bayes_hgen_name         ="", PY8_bayes_htruthfold_name   ="";//MC
//std::string PY8_bayes_hrec_name ="", PY8_bayes_hunf_name ="", PY8_bayes_hfold_name="";//Data

std::string     file_NLO_bayes = unf_dir+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_Bayes_07.02.18_NNPDF30nnlo_NLOMC.root";
std::string NLO_bayes_names[NPLTS]=
  {"theory_rnd_ynew_clone_rebins",
   "hJetQA_1wJetID_jtpt_divBylumietabin_clone_rebins","ppData_BayesUnf_Spectra"};

std::string  file_NLOwNP_bayes = unf_dir+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_Bayes_07.02.18_NNPDF30nnlo_NLOMC_wNP.root";
std::string NLOwNP_bayes_names[NPLTS]=
  {"theory_rnd_ynew_clone_rebins",
   "hJetQA_1wJetID_jtpt_divBylumietabin_clone_rebins","ppData_BayesUnf_Spectra"};

//SVD FILES            
std::string     file_PY8_SVD = unf_dir+"ak4PF_wjtID_anabins_unf_HPtJetTrig_w_Py8_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID.root";
std::string PY8_SVD_names[NPLTS]=
  {"hpp_gen_wJetID_R4_20_eta_20_clone_rebins",
   "hJetQA_1wJetID_jtpt_divBylumi_clone_rebins", "hunf_svd_kReg8"};

std::string     file_NLO_SVD = unf_dir+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID_NNPDF30nnlo_NLOMC.root";
std::string NLO_SVD_names[NPLTS]=
  {"theory_rnd_ynew_clone_rebins",
   "hJetQA_1wJetID_jtpt_divBylumi_clone_rebins", "hunf_svd_kReg8"};

std::string  file_NLOwNP_SVD = unf_dir+"ak4PFJets_wjtID_anabins_unf_HPtJetTrig_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID_NNPDF30nnlo_NLOMC_wNP.root";
std::string NLOwNP_SVD_names[NPLTS]=
  {"theory_rnd_NPynew_clone_rebins",
   "hJetQA_1wJetID_jtpt_divBylumi_clone_rebins","hunf_svd_kReg8"};

//i.e. match the priors, but make sure diff methods
void print_compareMethods(TFile* fin_p8byes, TFile* fin_NLObyes, TFile* fin_NLOwNPbyes,	 
			  TFile* fin_p8svd, TFile* fin_NLOsvd,   TFile* fin_NLOwNPsvd   ) 
{
  if( !((bool) fin_p8byes) || !((bool) fin_NLObyes) || !((bool) fin_NLOwNPbyes) ||	 
      !((bool) fin_p8svd ) || !((bool) fin_NLOsvd ) || !((bool) fin_NLOwNPsvd )  ){
    std::cout<<"file not found. exit."<<endl;
    assert(false);  }      


  //compare py8byes file w/ p8svd file
  // // compare p8byes unf data w/ meas data
  // // compare p8svd unf data w/ meas data
  // // compare p8byes unf data w/ p8byes prior
  // // compare p8svd unf data w/ p8svd prior
  // // compare p8byes unf data w/ p8svd unf data

  //compare NLObyes w/ NLOsvd
  // // compare NLObyes unf data w/ meas data
  // // compare NLOsvd unf data w/ meas data
  // // compare NLObyes unf data w/ NLObyes prior
  // // compare NLOsvd unf data w/ NLOsvd prior
  // // compare NLObyes unf data w/ NLOsvd unf data

  //compare NLOwNPbyes w/ NLOwNPsvd
  // // compare NLOwNPbyes unf data w/ meas data
  // // compare NLOwNPsvd unf data w/ meas data
  // // compare NLOwNPbyes unf data w/ NLOwNPbyes prior
  // // compare NLOwNPsvd unf data w/ NLOwNPsvd prior
  // // compare NLOwNPbyes unf data w/ NLOwNPsvd unf data


  return;
}

void print_comparePriors(TFile* fin_p8byes, TFile* fin_NLObyes, TFile* fin_NLOwNPbyes,
			 TFile* fin_p8svd, TFile* fin_NLOsvd, TFile* fin_NLOwNPsvd){
  if( !((bool) fin_p8byes) || !((bool) fin_NLObyes) || !((bool) fin_NLOwNPbyes) ||	 
      !((bool) fin_p8svd ) || !((bool) fin_NLOsvd ) || !((bool) fin_NLOwNPsvd )  ){
    std::cout<<"file not found. exit."<<endl;
    assert(false);   }
  
  //compare py8byes w/ NLObyes
  //compare py8byes w/ NLOwNPbyes
  //compare NLObyes w/ NLOwNPbyes

  //compare py8svd w/ NLOsvd
  //compare py8svd w/ NLOwNPsvd
  //compare NLOsvd w/ NLOwNPsvd

  return;
}


void quickCompare_unfoldedData(bool compareMethods=false, bool comparePriors=false){
  TFile* fin_PY8_bayes   =TFile::Open(file_PY8_bayes.c_str(),"READ");
  TFile* fin_NLO_bayes   =TFile::Open(file_NLO_bayes.c_str(),"READ");
  TFile* fin_NLOwNP_bayes=TFile::Open(file_NLOwNP_bayes.c_str(),"READ");
  
  TFile* fin_PY8_SVD   =TFile::Open(file_PY8_SVD.c_str(),"READ");
  TFile* fin_NLO_SVD   =TFile::Open(file_NLO_SVD.c_str(),"READ");
    TFile* fin_NLOwNP_SVD=TFile::Open(file_NLOwNP_SVD.c_str(),"READ");
    if(compareMethods){
      print_compareMethods( fin_PY8_bayes, fin_NLO_bayes, fin_NLOwNP_bayes,
			    fin_PY8_SVD,   fin_NLO_SVD,   fin_NLOwNP_SVD	  );
    }
    
    if(comparePriors){
      print_comparePriors( fin_PY8_bayes, fin_NLO_bayes, fin_NLOwNP_bayes,
			   fin_PY8_SVD,   fin_NLO_SVD,   fin_NLOwNP_SVD	 );
      
    }
    return;
}
