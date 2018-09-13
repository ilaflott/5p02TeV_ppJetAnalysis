#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>

#include <TH1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TH2D.h>






int main(int argc, char* argv[]){  int rStatus = -1;


  
  TFile* NPCorrFile=(TFile*)TFile::Open("NLOpNP_InclusiveJets5TeV.root");
  TFile* outputFile=new TFile("NLOpNP_InclusiveJets5TeV_customNP.root");
  //avg herwig EE4C and herwig EE5C together + save to root file... wait how to do this? misguided perhaps..
  
  return rStatus;
}
