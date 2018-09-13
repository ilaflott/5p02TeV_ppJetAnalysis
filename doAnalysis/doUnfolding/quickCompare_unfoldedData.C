#include <cstdlib>
#include <iostream>
#include "TH1.h"


//BAYESIAN FILES
std::string     file_PY8_bayes = "ak4PF_wjtID_anabins_unf_HPtJetTrig_w_Py8_Bayes_07.02.18_0.0eta2.0_gendr0p1_recoJetID.root";
std::string     file_NLO_bayes = "ak4PFJets_wjtID_anabins_unf_HPtJetTrig_Bayes_07.02.18_NNPDF30nnlo_NLOMC.root";
std::string  file_NLOwNP_bayes = "ak4PFJets_wjtID_anabins_unf_HPtJetTrig_Bayes_07.02.18_NNPDF30nnlo_NLOMC_wNP.root";

//SVD FILES            
std:string     file_PY8_SVD = "ak4PF_wjtID_anabins_unf_HPtJetTrig_w_Py8_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID.root"; 
std:string     file_NLO_SVD = "ak4PFJets_wjtID_anabins_unf_HPtJetTrig_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID_NNPDF30nnlo_NLOMC.root";
std:string  file_NLOwNP_SVD = "ak4PFJets_wjtID_anabins_unf_HPtJetTrig_SVD_07.02.18_0.0eta2.0_gendr0p1_recoJetID_NNPDF30nnlo_NLOMC_wNP.root"; 

void print_compareMethods(){
  return;
}

void print_comparePriors(){
  return;
}


void quickCompare_unfoldedData(bool compareMethods=false, bool comparePriors=false){
  if(compareMethods){
    print_compareMethods();
  }
  
  if(comparePriors){
    print_comparePriors();
  }
  return;
}
