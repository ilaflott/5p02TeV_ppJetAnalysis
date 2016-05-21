#include <cstdlib>
#include <vector>
#include <string>
#include <TLorentzVector.h>

///////////////////////////////////////////////
////////// (initializa/declara)tions //////////
///////////////////////////////////////////////

// main macros
int deriveDijetResponse(int startfile, int endfile, std::string infile_Forest,
                        int radius, bool isMC,
                        std::string outfile, bool debugMode);
const int deriveRespArgCount=7+2;
int sumDijetResponse(std::string filename, bool isMC, bool doDraw, bool debugMode);
const int sumRespArgCount=4+2;

// helper functions
int findBin(float value, int nbins, const double* xbins);
TLorentzVector findMissEt(int nPFpart_, std::vector<int>* pfId_, std::vector<float>* pfPt_,
                          std::vector<float>* pfEta_, std::vector<float>* pfPhi_,
                          float* eSum, float *phoSum,
                          int nref,    float* pt_F, float* rawpt_F,
                          float* eta_F, float* phi_F, float* m_F);

//for convenience during testing, coding, etc.
const bool atMIT  = true;
const std::string hadoopDir    = "/mnt/hadoop/cms";
const std::string xrootdDirect = "root://cmsxrootd.fnal.gov/";

//initializations
const int colors[] = {1,2,4,8,20};
const int N_colors= sizeof(colors)/sizeof(int);

const double xbins_pt[] = {40, 60, 80, 110, 200, 1000};
const int nbins_pt = sizeof(xbins_pt)/sizeof(double)-1;

const double xbins_eta[] = {
  -5.191, -4.716, -4.363, -4.013,
  -3.664, -3.314, -2.964,
  -2.650, -2.322, -2.043, -1.830,
  -1.653, -1.479, -1.305,
  -1.131, -0.957, -0.783, -0.609,
  -0.435, -0.261, -0.087,
  0.000,
  0.087,  0.261,  0.435, 0.609,
  0.783,  0.957,  1.131,
  1.305,  1.479,  1.653,  1.830,
  2.043,  2.322,  2.650,
  2.964,  3.314,   3.664, 4.013,
  4.363,  4.716,  5.191
};
const int nbins_eta = sizeof(xbins_eta)/sizeof(double)-1;

//const double minJetPt = 40;

////////////////
// findBin //
////////////////
int findBin(float value,
            int nbins, const double* xbins){
  for(int i=0; i<nbins; i++) if(value<xbins[i+1]) return i;
  return nbins-1;
}

////////////////
// findMissEt //
////////////////
TLorentzVector findMissEt(int nPFpart_, std::vector<int> * pfId_, std::vector<float> * pfPt_,
                          std::vector<float> * pfEta_, std::vector<float> * pfPhi_,
                          float* eSum, float *phoSum,
                          int nref, float* pt_F, float* rawpt_F,
                          float* eta_F, float* phi_F, float* m_F){

  // missEt = - Sum(Et you can see)
  TLorentzVector missEt(0,0,0,0);

//missing Et from particles
for(int i=0; i<nPFpart_; i++){//PF Loop
  TLorentzVector pfp;
  double pfMass=0;
  if(pfId_->at(i)==1 || pfId_->at(i)==4) pfMass = 0.1395702;
  //pfp.SetPtEtaPhiM(pfPt[i],pfEta[i],pfPhi[i],pfMass);
  pfp.SetPtEtaPhiM(pfPt_->at(i),pfEta_->at(i),pfPhi_->at(i),pfMass);
  pfp.SetPz(0.); //2d projection in x-y plane
  missEt += pfp;
 }//end PF loop


//missing Et from jets
for(int i=0; i<nref; i++){//begin jet loop
  if(pt_F[i]<=15) continue;//quick cut before expensive computations

  TLorentzVector jtTmp(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
  double jetEnergy=jtTmp.E();
  if( (eSum[i]+phoSum[i])/jetEnergy<0.9 ){
    TLorentzVector jt;
    jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
    jt.SetPz(0.); //project to x-y plane
    missEt += jt;
  }
 }//end jet loop
 
//  //missing Et from jets
//  for(int i=0; i<nref; i++){//begin jet loop
//    if(pf_F[i]<=15) continue;//quick cut before expensive computations
//    TLorentzVector jt(pt_F[i],eta_F[i],phi_F[i],m_F[i]);
//    if( (eSum[i]+phoSum[i])/jt.E()<0.9 ){
//      jt.SetPtEtaPhiM(pt_F[i]-rawpt_F[i],eta_F[i],phi_F[i],m_F[i]);
//      jt.SetPz(0.); //project to x-y plane
//      missEt += jt;
//    }
//  }//end jet loop
 
 missEt*=-1; // rotate 180 deg, to get missing Et,
 missEt.SetE(0.);  // only care about Et, not E
 //std::cout << "missEt:  Mag: " << missEt.Mag() << std::endl;
 return missEt;
}
