#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TF1.h"

using namespace std;
class L2ResidualJES
{
 private:
  int neta;
  int etacut;
  double lowerPtCut;
  double higherPtCut;
  bool dopPb;
  int radius;
  TFile *correctionFile;
  TF1 * fits[100], *fitInitial;
  TString algoCorr;
  TString mode;
  double etaMin[100];
  double etaMax[100];
  
  public:
  
   void reset()
   { 
    dopPb = false;
    for(int ieta=0;ieta<100;ieta++){
     fits[ieta] = NULL;
    }

    correctionFile = NULL;
    
   }
  
   L2ResidualJES(int radius=3, int etacut=3, TString mode="pp5")
   {
     reset();
   
     this->radius = radius;
     this->etacut = etacut;
     this->mode = mode;
     if(mode == "pPb5" || mode == "Pbp5") this->dopPb = true;
   
     algoCorr = Form("ak%dPF",radius);   
   
     lowerPtCut = 30;
     higherPtCut = 400;
   
     if(etacut==3){
       if(dopPb){
         neta=22;
	
         double etabins_hcalbins[]= {-3,-2.650,-2.500,  -2.172,-1.740, -1.479, -1.392, -1.218,-1.044,  -0.696,  -0.348,  0.000,  0.348,   0.696,   1.044,1.218, 1.392, 1.479,  1.740,   2.172,   2.500,   2.650, 3};
	     for(int ieta = 0; ieta < neta; ieta++){
	       etaMin[ieta] = etabins_hcalbins[ieta];
	       etaMax[ieta] = etabins_hcalbins[ieta+1];
	     }
	   }else{
         neta=58;
         double etabins_hcalbins[]= {-3, -2.853,
                         -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                         -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                         -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                         -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                          0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                          0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                          1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                          1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                          3};
					   
         for(int ieta = 0; ieta < neta; ieta++){
	         etaMin[ieta] = etabins_hcalbins[ieta];
	         etaMax[ieta] = etabins_hcalbins[ieta+1];
	     }
	   }
     }
     if(etacut==4){
	   neta=64;
	   double etabins_hcalbins4[]= {-4,      -3.664,  -3.314,  -2.964, -2.853,
                        -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
                        -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
                        -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
                        -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,  0.000,
                         0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
                         0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
                         1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
                         1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
                         2.964,  3.314,  3.664, 4};
       for(int ieta = 0; ieta < neta; ieta++){
	     etaMin[ieta] = etabins_hcalbins4[ieta];
	     etaMax[ieta] = etabins_hcalbins4[ieta+1];
       }
     }

     if( mode == "pPb5" ) correctionFile = new TFile(Form("JECDataDriven/L2L3VsPtEtaBinned_alphacut_high2_%s_etacut%d_dopPb%d_July4_pPb_double_const.root",algoCorr.Data(),etacut,dopPb));
     else if( mode == "Pbp5" ) correctionFile = new TFile(Form("JECDataDriven/L2L3VsPtEtaBinned_alphacut_high2_%s_etacut%d_dopPb%d_July4_Pbp_double_const.root",algoCorr.Data(),etacut,dopPb));
     else correctionFile = new TFile(Form("JECDataDriven/L2L3VsPtEtaBinned_alphacut_high2_%s_etacut%d_dopPb%d.root",algoCorr.Data(),etacut,dopPb));

     //cout << Form("JECDataDriven/L2L3VsPtEtaBinned_alphacut_high2_%s_etacut%d_dopPb%d.root",algoCorr.Data(),etacut,dopPb) << endl;
     
     if(mode == "pPb5") fitInitial = new TF1("fitInitial","9.87174e-01+(-5.45335e-03)*x+(1.31839e-02)*x*x",-3,3);
     if(mode == "Pbp5") fitInitial = new TF1("fitInitial","9.87174e-01+(-5.45335e-03)*x+(1.31839e-02)*x*x",-3,3);
     
     //cout << "neta = "<< neta <<  Form(" JECDataDriven/L2L3VsPtEtaBinned_alphacut_high2_%s_etacut%d_dopPb%d.root",algoCorr.Data(),etacut,dopPb) << endl;

     for(int ieta = 0;ieta<neta;ieta++){
       fits[ieta] = (TF1*)correctionFile->Get(Form("fit%d",ieta));
     } 
   }
   
   double getCorrectedPt(double jetpt, double jeteta)
   {
     double correction = 1;
     if( fabs(jeteta)> ((double)etacut)) return correction*jetpt;
     int etaindex = 0;
     for(int ieta = 0; ieta < neta; ieta++){
       if(etaMin[ieta] > jeteta ) continue;
	   else etaindex = ieta;
     }
      
   
    if(jetpt < lowerPtCut) return fits[etaindex]->Eval(lowerPtCut)*jetpt;
    if(jetpt > higherPtCut) return fits[etaindex]->Eval(higherPtCut)*jetpt;
 
    return fits[etaindex]->Eval(jetpt)*jetpt;
   }
  
   double getInitialCorr(float eta){
     return fitInitial->Eval(eta);
   }
   void test(){
   }
};

class L3ResidualJES
{
  private: 
    double constCorrection;
  public:
    L3ResidualJES(TString mode){
	  if( mode == "pp5" ) constCorrection = 9.81038e-01;
	  else constCorrection = 1;
    }	
	double getCorrectedPt(double pt){
	  return pt/constCorrection;
	}
};

class L2ResidualJER
{
  private:
    TF1 *fitSF[2];
    TF1 *fitGaus;
	
    const double ptLow = 50;
    const double ptHigh = 400;
	const double etaBound = 1.3;
	const double etaMax = 3.;
	int neta = 2;
	TString mode; //pp or pPb
  public:
  L2ResidualJER(TString mode){
    
	double parsPPb[2][2] = {{5.61139e+00,9.40866e-01},
                       {1.19898e+00,6.28925e-01}};
	  
	double parsPP[2][2] = {{7.57848e-01,6.23244e-01},
	                  {6.08142e+00,1.05733e+00}};
	  
	for(int ieta = 0; ieta < neta; ieta++){
     fitSF[ieta] = new TF1(Form("fitSF_eta%d",ieta),"[0]/pow(x,[1])",50,400);
    }
	
	this->mode = mode;
	for(int ieta = 0; ieta < neta; ieta++){
	  for(int ipar = 0; ipar < 2; ipar++){
  	    if( mode == "pPb5" ){
		  fitSF[ieta]->SetParameter(ipar, parsPPb[ieta][ipar]);
		}else{
		  fitSF[ieta]->SetParameter(ipar, parsPP[ieta][ipar]);
		}
	  }
	}
    fitGaus = new TF1("fitGaus","gaus(0)",-20,20);
    fitGaus->SetParameters(1,0,1);
  }
 
  int getIEta(double eta){
   if( fabs(eta) < 1.3 ) return 0;
   else return 1;
  }
 
 double getSmearedPt(double pt, double eta){
   if( pt < ptLow || pt > ptHigh || fabs(eta) > etaMax ) return pt;
   return pt*(1+(fitSF[getIEta(eta)]->Eval(pt))*fitGaus->GetRandom());
 }
};
