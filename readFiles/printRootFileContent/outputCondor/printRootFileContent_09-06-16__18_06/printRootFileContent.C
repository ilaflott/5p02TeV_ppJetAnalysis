#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TProfile.h>




int main(int argc, char *argv[]){

  std::string Datafile1="/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtJet80/HiForestAOD_1.root";
  std::cout<<"opening + listing contents of Datafile1, HighPtJet80"<<std::endl;
  TFile* file1= (TFile*) TFile::Open(Datafile1.c_str());
  file1->ls();

  std::string Datafile2="/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtLowerJets/HiForestAOD_1.root";
  std::cout<<"opening + listing contents of Datafile2, HighPtLowerJets"<<std::endl;
  TFile* file2= (TFile*) TFile::Open(Datafile2.c_str());
  file2->ls();

  std::string Datafile3="/mnt/hadoop/cms/store/user/dgulhan/MinimumBias1/HiForest_MinimumBias1_Run2015E_PromptReco_pp_5020GeV/151201_113508/0000/HiForest_1.root";
  std::cout<<"opening + listing contents of Datafile3, MinBias"<<std::endl;
  TFile* file3= (TFile*) TFile::Open(Datafile3.c_str());
  file3->ls();


  std::string MCfile1="/mnt/hadoop/cms/store/user/krajczar/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV/Pythia8_Dijet120_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/160119_133414/0000/HiForestAOD_1.root";
  std::cout<<"opening + listing contents of MCFile1, PY8 QCD DijetMC"<<std::endl;
  TFile* file1_MC= (TFile*) TFile::Open(MCfile1.c_str());
  file1_MC->ls();

  std::string MCfile2="/mnt/hadoop/cms/store/user/rbi/merged/MinBias_TuneCUETP8M1_5p02TeV-pythia8-HINppWinter16DR-NoPU_75X_mcRun2_asymptotic_ppAt5TeV_forest_v2/0.root";
  std::cout<<"opening + listing contents of MCfile2, PY8 MinBiasMC"<<std::endl;
  TFile* file2_MC= (TFile*) TFile::Open(MCfile2.c_str());
  file2_MC->ls();
  return 0;
  }
