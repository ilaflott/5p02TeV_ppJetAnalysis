#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cassert>

#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>

#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <TSpline.h>
#include <TLine.h>

#include "RooUnfoldResponse.h"

using namespace std;
//TH1D *cross_section(TH1D *Spectra_);

TH1D* applyNPtoxsec(TH1D* xsec, TF1* fNP){
  const bool funcDebug=true;
  if(funcDebug)xsec->Print("base");
  if(funcDebug){
    std::cout<<"pre NP Corr, xsec mean along x axis is = " << xsec->GetMean(1) << std::endl;
    std::cout<<"pre NP Corr, xsec mean along y axis is = " << xsec->GetMean(2) << std::endl;
  }
  //if(funcDebug)fNP->Print("base");  
  TH1D* xsec_wNP=(TH1D*)xsec->Clone( ( 
				      ((std::string)xsec->GetName()) + "_wNP"
				       ).c_str() 
				     );
  xsec_wNP->Reset();
  xsec_wNP->SetTitle( (
		       "NP+"+((std::string)xsec->GetTitle())
		       ).c_str()
		      );
  xsec_wNP->Sumw2(true);
  int nbinsx=xsec->GetNbinsX();
  for(int i=1; i<=nbinsx; ++i){
    Double_t bincent = xsec->GetBinCenter(i);
    Double_t NPval=fNP->Eval(bincent);    
    Double_t binval  = xsec->GetBinContent(i);
    Double_t binerr = xsec->GetBinError(i);
    binval*= NPval;
    binerr*= NPval;
    xsec_wNP->SetBinContent(i,binval);
    xsec_wNP->SetBinError(i,binerr);
  }
  if(funcDebug)xsec_wNP->Print("base");
  cout<<"done making theory cross section hist"<<endl;
  if(funcDebug){
    std::cout<<"post NP Corr, xsec mean along x axis is = " << xsec_wNP->GetMean(1) << std::endl;
    std::cout<<"post NP Corr, xsec mean along y axis is = " << xsec_wNP->GetMean(2) << std::endl;
  }
  return (TH1D*)xsec_wNP;
}

TH1D* cross_section(TH1D* Spectra_){  
  /// Evaluate xSections
  TH1D *Xsection_=(TH1D*)Spectra_->Clone();
  for(int i=0; i<Spectra_->GetNbinsX(); ++i){
    Double_t bin_value = Spectra_->GetBinContent(1+i);
    Double_t bin_error = Spectra_->GetBinError(1+i);
    Double_t bin_width = Spectra_->GetBinWidth(1+i);
    Double_t new_bin_value = bin_value/bin_width;
    Double_t new_bin_error = bin_error/bin_width;
    Xsection_->SetBinContent(1+i,new_bin_value);
    Xsection_->SetBinError(1+i,new_bin_error);
  }
  //   std::cout<<Integral_<<std::endl;
  
  return Xsection_;
}
void divideBinWidth(TH1D* h){  
  /// Evaluate xSections
  //TH1D *Xsection_=(TH1D*)Spectra_->Clone();
  Int_t nbins=h->GetNbinsX();
  for(int i=1; i<=nbins; ++i){
    Double_t bin_value = h->GetBinContent(i);
    Double_t bin_error = h->GetBinError(i);
    Double_t bin_width = h->GetBinWidth(i);
    Double_t new_bin_value = bin_value/bin_width;
    Double_t new_bin_error = bin_error/bin_width;
    h->SetBinContent(i,new_bin_value);
    h->SetBinError(i,new_bin_error);
  }
  //   std::cout<<Integral_<<std::endl;
  return;
}
void multiplyBinWidth(TH1D* h){  
  /// Evaluate xSections
  //TH1D *Xsection_=(TH1D*)Spectra_->Clone();
  Int_t nbins=h->GetNbinsX();
  for(int i=1; i<=nbins; ++i){
    Double_t bin_value = h->GetBinContent(i);
    Double_t bin_error = h->GetBinError(i);
    Double_t bin_width = h->GetBinWidth(i);
    Double_t new_bin_value = bin_value*bin_width;
    Double_t new_bin_error = bin_error*bin_width;
    h->SetBinContent(i,new_bin_value);
    h->SetBinError(i,new_bin_error);
  }
  //   std::cout<<Integral_<<std::endl;
  return;
}

/////////////// Inclusive jets Binning form Ksenia ////////////////////////////////////////////////////////
double thyBins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //  32., 37., //junk bins above 
  //43., 49., 
  //37.,
  43., 49.,
  56., 64., 74., 84., 97., 
  114., 133., 153., 174., 196., 
  220., 245., 272., 
  300., 330., 362., 395., 
  430., 468.,
  507., 548., 592., 
  638., 
  686.//, 
  //1000.//,//967.//1032. 
  //1500.
  //  1032., 1101. //junk from here down
}; 
const int n_thybins_incl=sizeof(thyBins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1
double thyBins_incl2[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //  32., 37., //junk bins above 
  //43., 49., 
  //37.,
  43., 49.,
  56., 64., 74., 84., 97., 
  114., 133., 153., 174., 196., 
  220., 245., 272., 
  300., 330., 362., 395., 
  430., 468.,
  507., 548., 592., 
  638., 
  686.//, 
  //1000.//,//967.//1032. 
  //1500.
  //  1032., 1101. //junk from here down
}; 
const int n_thybins_incl2=sizeof(thyBins_incl2)/sizeof(double)-1; //this is # of bins = # of entries in array - 1


double smearedBins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //37.,
  43., 49.,
  56., 64., 74., 84., 97., 
  114., 133., 153., 174., 196., 
  220., 245., 272., 
  300., 330., 362., 395., 
  430., 468.,
  507., 548., 592., 
  638., 
  686.//, 
  //1000.//,
  //1500.//967.//1032.
  //,1500.
  //  1032., 1101. //junk from here down
}; 
const int n_smearedbins_incl=sizeof(smearedBins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1


////////// NPs from Paolo 4/9/2015
// THESE ARE LIKELY FOR R>0.4
// RAGHAV SAYS NP CORR FACTORS > 1 FOR R=0.5 AND ABOVE
/// y = [0] + [1] * pow (x , [2] )
///   The values of the parameters are for each bin:
///                      [0]            [1]             [2]
/// rapidity bin 1:    0.998262       314.329        -1.60473
/// rapidity bin 2:    1.00199        593.25         -1.73981
/// rapidity bin 3:    1.00209        325.315        -1.61465
/// rapidity bin 4:    1.01009        474.01         -1.72148
/// rapidity bin 5:    1.02156        1854.2         -2.05028
/// rapidity bin 6:    1.01768        277.946        -1.66097
/// rapidity bin 7:   -65.6541        67.1291        -0.00113689

double NP_y0(double x){
  return (0.998262 + 314.329 * pow (x ,-1.60473));}
string NP_y0_str="0.998262 + 314.329 * pow (x ,-1.60473)";

double NP_y1(double x){
  return (1.00199  + 593.25  * pow (x ,-1.73981));}
string NP_y1_str="1.00199  + 593.25  * pow (x ,-1.73981)";

double NP_y2(double x){
  return (1.00209  + 325.315 * pow (x ,-1.61465));}
string NP_y2_str="1.00209  + 325.315 * pow (x ,-1.61465)";

double NP_y3(double x){
  return (1.01009  + 474.01  * pow (x ,-1.72148));}
string NP_y3_str="1.01009  + 474.01  * pow (x ,-1.72148)";

double NP_y4(double x){
  return (1.02156  + 1854.2  * pow (x ,-2.05028));}
string NP_y4_str="1.02156  + 1854.2  * pow (x ,-2.05028)";

double NP_y5(double x){
  return (1.01768  + 277.946 * pow (x ,-1.66097));}
string NP_y5_str="1.01768  + 277.946 * pow (x ,-1.66097)";

double NP_y6(double x){
  return (-65.6541 + 67.1291 * pow (x ,-0.00113689));}
string NP_y6_str="-65.6541 + 67.1291 * pow (x ,-0.00113689)";

double NP_ynew(double x){  
  return ( ( NP_y0(x) + NP_y1(x) + NP_y2(x) + NP_y3(x) ) / 4. ); 
}
string NP_ynew_str= 
  "( (" + 
  NP_y0_str + " + " + 
  NP_y1_str + " + " + 
  NP_y2_str + " + " + 
  NP_y3_str + 
  " )/4. )";


////////// Gaussian Core Distributions from Ian L 10/3/17
// from JER fits using, sigma/mu = [0] + [1] /( pow (gen p_T , [2] ) + [3] * gen p_t)
double gJER_ynew(double x){
  return ( 4.80963e-02 + 1.51797e-02/(pow(x , -5.78569e-01)+(2.17736e-03*x)) );}
string gJER_ynew_str="4.80963e-02+1.51797e-02/(pow(x,-5.78569e-01)+(2.17736e-03*x))";



TH1D* make00eta20Hist( TH1D * h_1, 
		       TH1D * h_2, 
		       TH1D * h_3, 
		       TH1D * h_4  ){
  
  bool funcDebug=true;
  
  //TH1D* h00eta20 = new TH1D( "h00eta20","thyHist 00eta20", numBins, (const double*) thyBins_incl_00eta20 );
  std::string h00eta20_title="h00eta20";
  if( ((std::string)h_1->GetName()).find("_wNP") )
    h00eta20_title+="_wNP";
  TH1D* h00eta20 = (TH1D*)h_1->Clone(h00eta20_title.c_str());
  h00eta20->Reset("ICES");
  h00eta20->Reset("M");
  h00eta20->Sumw2(true);
  h00eta20->Print("base");
  
  //this part assumes all hists start at same lower limit; which they better!
  if(funcDebug)cout<<endl;
  if(funcDebug)cout<<"filling new hist!"<<endl;
  if(funcDebug)cout<<endl;
  int n_thybins=h00eta20->GetNbinsX();
  
  for(int i=1; i<=n_thybins; i++) {
    
    if(funcDebug)cout<<"i="<<i<<endl;
    if(funcDebug)cout<<"lower edge of bin i for h00eta20="<<h00eta20->GetXaxis()->GetBinLowEdge(i)<<endl;
    if(funcDebug)cout<<"lower edge of bin i for h_1="<<h_1->GetXaxis()->GetBinLowEdge(i)<<endl;
    float v=0., v_err=0.;
    
    float v_1=0., v_1err=0.;
    if(i <= h_1->GetNbinsX() ) {
      v_1=h_1->GetBinContent(i);
      v_1err=h_1->GetBinError(i); 
      v+=v_1;
      v_err+=v_1err*v_1err;    }
    else 
      cout<<"no more bins for hist h_1!"<<endl;
    
    float v_2=0., v_2err=0.;
    if(i <= h_2->GetNbinsX() ) {
      v_2=h_2->GetBinContent(i);
      v_2err=h_2->GetBinError(i); 
      v+=v_2;
      v_err+=v_2err*v_2err;      }
    else
      cout<<"no more bins for hist h_2!"<<endl;
    
    float v_3=0., v_3err=0.;
    if(i <= h_3->GetNbinsX() ) {
      v_3=h_3->GetBinContent(i);
      v_3err=h_3->GetBinError(i); 
      v+=v_3;
      v_err+=v_3err*v_3err;        }
    else 
      cout<<"no more bins for hist h_3!"<<endl;
    
    float v_4=0., v_4err=0.;
    if(i <= h_4->GetNbinsX() ) {
      v_4=h_4->GetBinContent(i);
      v_4err=h_4->GetBinError(i); 
      v+=v_4;
      v_err+=v_4err*v_4err;        }
    else 
      cout<<"no more bins for hist h_4!"<<endl;
    
    v/=4.;//for etabin width    
    //v/=1000.;//picobarns to nanobarns
    
    v_err=sqrt(v_err); //for err
    v_err/=4.;
    //v_err/=1000.;//picobarns to nanobarns
    
    if(funcDebug)cout<<"v="<<v<<endl;
    if(funcDebug)cout<<"v_err="<<v_err<<endl;
    if(funcDebug)cout<<"setting bin content..."<<endl;
    h00eta20->SetBinContent(i, v);
    h00eta20->SetBinError(i, v_err);
    if(funcDebug)cout<<"histogram contents..."<<endl;
    if(funcDebug)cout<<"bincontent="<<h00eta20->GetBinContent(i)<<endl;
    if(funcDebug)cout<<"binerr="<<h00eta20->GetBinError(i)<<endl;
    if(funcDebug)cout<<endl;
  }
  
  return h00eta20;
}



//int getRandomSeed(){
//  bool funcDebug=true;  
//  int
//}





void makeToySpectra(TH1D* hthy, 
		    TSpline3* hthy_spline, 
		    TF1* fJER, 
		    int nevts, 
		    TH1D* hthy_toyMC, TH1D* hsmeared_toyMC, TH1D* hsmeared_toyMC_test, 
		    TH2D* resp ){
  bool funcDebug=true;
  int tenth_nEvents=nevts/10;
  
  srand((unsigned)time(0));  
  UInt_t rnd_seed=rand();
  if(funcDebug)cout<<"rnd_seed="<<rnd_seed<<endl;
  TRandom3 *rnd = new TRandom3(rnd_seed);
  
  //srand((unsigned)time(0));
  UInt_t rnd_test_seed=rand();
  if(funcDebug)cout<<"rnd_test_seed="<<rnd_test_seed<<endl;
  TRandom3 *rnd_test = new TRandom3(rnd_test_seed);
 
  assert(false);          
  double ptmin_thy=hthy->GetBinLowEdge(1);
  double ptmax_thy=hthy->GetBinLowEdge(hthy->GetNbinsX()) + hthy->GetBinWidth(hthy->GetNbinsX());
  double ptmin_smeared=ptmin_thy, ptmax_smeared=ptmax_thy;
  
  int respcount=0, misscount=0, fakecount=0;

  for(int i=0;i<nevts;++i){      
    


    if(i%tenth_nEvents==0)
      cout<<"throwing random #'s for event # "<<i<<endl;

    double ptTrue  = rnd->Uniform(ptmin_thy,ptmax_thy);

    double sigma   = fJER->Eval(ptTrue);    //sigma*=1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment

    double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);

    double pt_w      =  hthy_spline->Eval(ptTrue);

    bool in_smearpt_range=(  ( ptSmeared>ptmin_smeared ) && ( ptSmeared<ptmax_smeared )  );
    bool in_trupt_range=(  ( ptTrue>ptmin_thy )        && ( ptTrue<ptmax_thy )  )  ;
    bool fillresp=in_smearpt_range&&in_trupt_range;


    if(fillresp){
      resp->Fill(ptSmeared,ptTrue,pt_w);
      hthy_toyMC    -> Fill(ptTrue,pt_w);
      hsmeared_toyMC -> Fill(ptSmeared,pt_w);      
      respcount++;
    }
    else if(in_trupt_range){//not in smear pt range but in truept range
      //      resp->Miss(ptTrue,pt_w);                                        
      hthy_toyMC    -> Fill(ptTrue,pt_w);
      misscount++;
    }
    else if(in_smearpt_range){//shouldnt get filled by construction. ptTrue always pulled s.t. in_trupt_range is true. if in_smearpt_range is true, so it in_trupt_range
      //      resp->Fake(ptSmeared,pt_w);                                        
      hsmeared_toyMC    -> Fill(ptSmeared,pt_w);
      fakecount++;
    }
    
    double    ptTrue_test = rnd_test->Uniform(ptmin_thy,ptmax_thy);
    double     sigma_test = fJER->Eval(ptTrue_test);    //sigma*=1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment
    double ptSmeared_test =  rnd_test->Gaus(ptTrue_test,ptTrue_test*sigma_test);
    double      pt_w_test =  hthy_spline->Eval(ptTrue_test);
    
    bool in_test_smearpt_range=(  ( ptSmeared_test>ptmin_smeared ) && ( ptSmeared_test<ptmax_smeared )  );
    if(in_test_smearpt_range)
      hsmeared_toyMC_test->Fill(ptSmeared_test,pt_w_test);
    

  }//end loop throwing random no's

  cout<<"done smearing ynew."<<endl;  
  cout<<"smear summary;"<<endl;
  cout<<"nEvents="<<nevts<<endl;
  cout<<"response_count="<<respcount<<endl;
  cout<<"miss_count="<<misscount<<endl;
  cout<<"fake_count="<<fakecount<<endl;
  
  return;
}



//void smear_theory_spectra(TH1D* hthy, TH1D* hsmrd, RooUnfoldResponse respmat, TSpline3* spline, TF1* fJER){
// 
//  bool funcDebug=true;
//  
//  return;
//}











void setTH2_ZAxisRange(TH2* h){
  bool funcDebug=false;
  int nbinsx=h->GetNbinsX();
  int nbinsy=h->GetNbinsY();
  if(funcDebug){
    std::cout<<"nbinsx="<<nbinsx<<std::endl;
    std::cout<<"nbinsy="<<nbinsy<<std::endl;  }
  double maxVal=1.e-30, minVal=1.e+30;
  if(funcDebug){
    std::cout<<"maxVal="<<maxVal<<std::endl;
    std::cout<<"minVal="<<minVal<<std::endl;  }
  
  for (int i=1;i<=nbinsx;i++){
    for (int j=1;j<=nbinsy;j++){
      if(funcDebug)std::cout<<"(i, j) = ("<<i<<", "<<j<<")"<<std::endl;
      double val=h->GetBinContent(i,j);

      if(val>maxVal){    
	maxVal=val;
	if(funcDebug)std::cout<<"maxVal="<<maxVal<<std::endl;
      }
      if(val<minVal){
	if(val>0.){
	  minVal=val;
	  if(funcDebug)
	    std::cout<<"minVal="<<minVal<<std::endl;		}
      }       
    }
  }
  
  
  h->SetAxisRange(minVal/2.,maxVal*2.,"Z");
  
  return;  
}
