#include <cstdlib>
#include <cassert>
#include <cmath>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>

//#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>

#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TSpline.h>
#include <TLine.h>
#include <TVirtualFitter.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>
#include <TStopwatch.h>

 //#include "RooUnfoldResponse.h"
#include "smearTheorySpectra_strs.h"

 //TH1D *cross_section(TH1D *Spectra_);

void debugcout(int debugInt=0){
  std::cout<<"location #"<<debugInt<<std::endl;
  debugInt++;
  return; 
}


TH1D* applyNPtoxsec(TH1D* xsec, TF1* fNP){
  const bool funcDebug=true;
  
  if(funcDebug)std::cout<<"pre NP Corr, xsec mean along x axis is = " << xsec->GetMean(1) << std::endl;
  xsec->Print("base");

  TH1D* xsec_wNP=(TH1D*)xsec->Clone( ( 
				      ((std::string)xsec->GetName()) + "_wNP"
				       ).c_str() 
				     );
  xsec_wNP->Reset("M ICES");
  xsec_wNP->SetTitle( (
		       "NP+"+((std::string)xsec->GetTitle())
		       ).c_str()
		      );
  xsec_wNP->Sumw2(true);
  
  if(funcDebug)std::cout<<"NP Fit Obj Name is: " << fNP->GetName()<<std::endl;
  if(funcDebug)std::cout<<"xsec Name is: " << xsec->GetName()<<std::endl;
  if(funcDebug)std::cout<<"xsec_wNP Name is: " << xsec_wNP->GetName()<<std::endl;
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
  
  if(funcDebug)std::cout<<"post NP Corr, xsec mean along x axis is = " << xsec_wNP->GetMean(1) << std::endl;
  xsec_wNP->Print("base");
  
  std::cout<<"done making theory cross section hist"<<std::endl;
  return (TH1D*)xsec_wNP;
}

//TH1D* cross_section(TH1D* Spectra_){  
//  /// Evaluate xSections
//  TH1D *Xsection_=(TH1D*)Spectra_->Clone();
//  for(int i=0; i<Spectra_->GetNbinsX(); ++i){
//    Double_t bin_value = Spectra_->GetBinContent(1+i);
//    Double_t bin_error = Spectra_->GetBinError(1+i);
//    Double_t bin_width = Spectra_->GetBinWidth(1+i);
//    Double_t new_bin_value = bin_value/bin_width;
//    Double_t new_bin_error = bin_error/bin_width;
//    Xsection_->SetBinContent(1+i,new_bin_value);
//    Xsection_->SetBinError(1+i,new_bin_error);
//  }
//  //   std::cout<<Integral_<<std::endl;
//  
//  return Xsection_;
//}

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
double py8Bins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //  32., 37., //junk bins above 
  //43., 49., 
  //37.,
  //43., 
  49.,
  56., 
  64., 
  74., 
  84., 
  97., 
  114.,
  133., 
  153., 
  174., 
  196., 
  220., 
  245.,
  272., 
  300.,
  330.,
  362., 
  395., 
  430., 
  468.,
  507., 
  548., 
  592., 
  638., 
  686.//,
  //1000.//,//967.//1032. 
  //1500.
  //  1032., 1101. //junk from here down
}; 
const int n_py8bins_incl=sizeof(py8Bins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1

double py8smearedBins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //37.,
  //43.,
  49.,
  56.,
  64., 
  74., 
  84., 
  97., 
  114., 
  133., 
  153., 
  174.,
  196., 
  220., 
  245., 
  272., 
  300., 
  330., 
  362.,
  395., 
  430., 
  468.,
  507.,
  548.,
  592., 
  638., 
  686.//, 
  //  1000.//,
  //1500.//967.//1032.
  //,1500.
  //  1032., 1101. //junk from here down
}; 
const int n_py8smearedbins_incl=sizeof(py8smearedBins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1


double thyBins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //  32., 37., //junk bins above 
  //43., 49., 
  //37.,
  43., 
  49.,
  56., 
  64., 
  74., 
  84., 
  97., 
  114.,
  133., 
  153., 
  174., 
  196., 
  220., 
  245.,
  272., 
  300.,
  330.,
  362., 
  395., 
  430., 
  468.,
  507., 
  548., 
  592., 
  638., 
  686.//,
  //1000.//,//967.//1032. 
  //1500.
  //  1032., 1101. //junk from here down
}; 
const int n_thybins_incl=sizeof(thyBins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1

double smearedBins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //37.,
  43.,
  49.,
  56.,
  64., 
  74., 
  84., 
  97., 
  114., 
  133., 
  153., 
  174.,
  196., 
  220., 
  245., 
  272., 
  300., 
  330., 
  362.,
  395., 
  430., 
  468.,
  507.,
  548.,
  592., 
  638., 
  686.//, 
  //  1000.//,
  //1500.//967.//1032.
  //,1500.
  //  1032., 1101. //junk from here down
}; 
const int n_smearedbins_incl=sizeof(smearedBins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1

const double absetabins[]={
  0.0,
  0.5,
  1.0,
  1.5,
  2.0,
  2.5,
  3.0//,
  //  3.2,
  //  4.7
};
const int n_absetabins=sizeof(absetabins)/sizeof(double)-1;
const std::string absetabins_str[]={
  "0.0",
  "0.5",
  "1.0",
  "1.5",
  "2.0",
  "2.5",
  "3.0"//,
  //  "3.2",
  //  "4.7"
};
const std::string absetabins_tags[]={
  "00eta05",
  "05eta10",
  "10eta15",
  "15eta20",
  "20eta25",
  "25eta30",
};



// 1.079 etc. --- > from 8TeV/2.76 TeV Incl Jets. analysis, AN-14-160, table 5
// 1.01 etc.  --- > from 5 TeV HIN studies showing excellent MC/Data JER closure
//one scale factor per |y| bin. #'s only available up to |y| < 3.0
const double JERscaleFacts[]={
1.01,//  1.079,
1.02,//  1.099,
1.03,//  1.121,
1.04,//  1.208,
1.05,//  1.254,
1.06//  1.395  
};
const double JERscaleFactErrs[]={
0.01,//  0.026,
0.02,//  0.028,
0.03,//  0.029,
0.04,//  0.046,
0.05,//  0.062,
0.06//  0.063
};


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
std::string NP_y0_str="0.998262 + 314.329 * pow (x ,-1.60473)";

double NP_y1(double x){
  return (1.00199  + 593.25  * pow (x ,-1.73981));}
std::string NP_y1_str="1.00199  + 593.25  * pow (x ,-1.73981)";

double NP_y2(double x){
  return (1.00209  + 325.315 * pow (x ,-1.61465));}
std::string NP_y2_str="1.00209  + 325.315 * pow (x ,-1.61465)";

double NP_y3(double x){
  return (1.01009  + 474.01  * pow (x ,-1.72148));}
std::string NP_y3_str="1.01009  + 474.01  * pow (x ,-1.72148)";

double NP_y4(double x){
  return (1.02156  + 1854.2  * pow (x ,-2.05028));}
std::string NP_y4_str="1.02156  + 1854.2  * pow (x ,-2.05028)";

double NP_y5(double x){
  return (1.01768  + 277.946 * pow (x ,-1.66097));}
std::string NP_y5_str="1.01768  + 277.946 * pow (x ,-1.66097)";

double NP_y6(double x){
  return (-65.6541 + 67.1291 * pow (x ,-0.00113689));}
std::string NP_y6_str="-65.6541 + 67.1291 * pow (x ,-0.00113689)";

double NP_ynew(double x){  
  return ( ( NP_y0(x) + NP_y1(x) + NP_y2(x) + NP_y3(x) ) / 4. ); 
}
std::string NP_ynew_str= 
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
std::string gJER_ynew_str="4.80963e-02+1.51797e-02/(pow(x,-5.78569e-01)+(2.17736e-03*x))";



//adds a parabolic extension of the spline, smoothly reducing weights to 0 (avoid neg weight in spline3)
//void make_spline3_ext(TH1D* hthy, TSpline3* hthy_spline3, TF1* hthy_spline3_ext){
void make_spline3_ext(TSpline3* hthy_spline3, TF1* hthy_spline3_ext , double x1, double x2, double y2){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"calculating spline3 extension for between x1="<<x1<<" and x2="<<x2<<std::endl;
  
  double y1=hthy_spline3->Eval(x1);
  //double y2=0.;//don't delete.
  
  //note; f=Ax^2+Bx+C
  //hthy_spline3_ext=new TF1( (hthy_spline3->GetName() + "_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1,x2);  
  //a few boundary conditions are assumed here;
  // 1; fext(x1)=spl3(x1)
  // 2; fext(x2)=0.
  // 3; fext'(x2)=0. (first deriv. condition) 
  
  //note, A/B/C here calculated assuming y2=0. more general formula so i can mess w/ end points should be put in here at some point for more tests
  if(funcDebug){
    std::cout<<"y1="<<y1<<std::endl;
    std::cout<<"y2="<<y2<<std::endl;
    //assert(false);
  }
  double A= (y1-y2)/( (x1*x1) + (x2*x2) - (2.*x1*x2) );
  double B= (-2.)*A*x2;
  double C= y2 - (1.*(A*pow(x2,2) + B*x2));
  hthy_spline3_ext->SetParameter(0,A);
  hthy_spline3_ext->SetParameter(1,B);
  hthy_spline3_ext->SetParameter(2,C);    
  hthy_spline3_ext->SetLineColor(kRed);
  
  if(funcDebug){
    std::cout<<"spline3 calculation done!"<<std::endl;
    std::cout<<"A="<<hthy_spline3_ext->GetParameter(0)<<std::endl;
    std::cout<<"B="<<hthy_spline3_ext->GetParameter(1)<<std::endl;
    std::cout<<"C="<<hthy_spline3_ext->GetParameter(2)<<std::endl;
    
    std::cout<<"compare spline3 at halfway point between x1="<<x1<<" and x2="<<x2<<std::endl;
    double xtest=x1+(x2-x1)/2.;
    std::cout<<"xtest=x1+(x2-x1)/2.="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;

    xtest+=10.;
    std::cout<<"xtest="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;

    xtest+=10.;
    std::cout<<"xtest="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;

    xtest+=10.;
    std::cout<<"xtest="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;
  }
  
  return;
}


void make_spline3_extv2(TSpline3* hthy_spline3, TF1* hthy_spline3_ext , double x1, double x2, double y2, double xprime, double yprime){
  bool funcDebug=true;
  if(funcDebug)std::cout<<"calculating spline3 extension for between x1="<<x1<<" and x2="<<x2<<std::endl;
  
  double y1=hthy_spline3->Eval(x1);
  //double y2=0.;//don't delete.
  
  //note; f=Ax^2+Bx+C
  //hthy_spline3_ext=new TF1( (hthy_spline3->GetName() + "_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1,x2);  
  //a few boundary conditions are assumed here;
  // 1; fext(x1)=spl3(x1)
  // 2; fext(xprime)=yprime
  // 3; fext'(x2)=0. (first deriv. condition) 
  //comparing w/ the original make_spline3_ext, this one makes sure that the parabola passes through the final bin's height at the bin center.
  
  //note, A/B/C here calculated assuming y2=0. more general formula so i can mess w/ end points should be put in here at some point for more tests
  if(funcDebug){
    std::cout<<"x1="<<x1<<std::endl;
    std::cout<<"y1="<<y1<<std::endl;
    std::cout<<"x2="<<x2<<std::endl;
    std::cout<<"y2="<<y2<<std::endl;
    std::cout<<"xprime="<<xprime<<std::endl;
    std::cout<<"yprime="<<yprime<<std::endl;
    //assert(false);
  }
  double A= (y1-yprime)/( (x1*x1) - (xprime*xprime) - 2.*(x2*x1 - x2*xprime) );
  double B= (-2.)*A*x2;
  double C= y2 - (1.*((A*x2*x2) + B*x2));
  //double C= yprime - (A*xprime*xprime) - (B*xprime);
  hthy_spline3_ext->SetParameter(0,A);
  hthy_spline3_ext->SetParameter(1,B);
  hthy_spline3_ext->SetParameter(2,C);    
  hthy_spline3_ext->SetLineColor(kRed);
  
  if(funcDebug){
    std::cout<<"spline3 calculation done!"<<std::endl;
    std::cout<<"A="<<hthy_spline3_ext->GetParameter(0)<<std::endl;
    std::cout<<"B="<<hthy_spline3_ext->GetParameter(1)<<std::endl;
    std::cout<<"C="<<hthy_spline3_ext->GetParameter(2)<<std::endl;
    
    std::cout<<"compare spline3 at halfway point between x1="<<x1<<" and x2="<<x2<<std::endl;
    double xtest=x1+(x2-x1)/2.;
    std::cout<<"xtest=x1+(x2-x1)/2.="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;

    xtest+=10.;
    std::cout<<"xtest="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;

    xtest+=10.;
    std::cout<<"xtest="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;

    xtest+=10.;
    std::cout<<"xtest="<<xtest<<std::endl;
    std::cout<<"spl3(xtest)="<<hthy_spline3->Eval(xtest)<<std::endl;
    std::cout<<"spl3ext(xtest)="<<hthy_spline3_ext->Eval(xtest)<<std::endl;
  }
  //  assert(false);
  return;
}

TH1D* make00eta20Hist( TH1D * h_1, TH1D * h_2, TH1D * h_3, TH1D * h_4,
		       std::string histTitle){
  
  bool funcDebug=false;  
  
  TH1D* hfinal_00eta20 = (TH1D*)h_1->Clone(histTitle.c_str());
  hfinal_00eta20->Reset("M ICES");
  hfinal_00eta20->Sumw2(true);
  hfinal_00eta20->Print("base");
  
  //this part assumes all hists start at same lower limit; which they better!
  if(funcDebug)std::cout<<std::endl;
  if(funcDebug)std::cout<<"filling new hist!"<<std::endl;
  if(funcDebug)std::cout<<std::endl;
  int n_thybins=hfinal_00eta20->GetNbinsX();
  
  for(int i=1; i<=n_thybins; i++) {
    
    if(funcDebug)std::cout<<"i="<<i<<std::endl;
    if(funcDebug)std::cout<<"lower edge of bin i for hfinal_00eta20="<<hfinal_00eta20->GetXaxis()->GetBinLowEdge(i)<<std::endl;
    if(funcDebug)std::cout<<"lower edge of bin i for h_1="<<h_1->GetXaxis()->GetBinLowEdge(i)<<std::endl;
    float v=0., v_err=0.;
    
    float v_1=0., v_1err=0.;
    if(i <= h_1->GetNbinsX() ) {
      v_1=h_1->GetBinContent(i);
      v_1err=h_1->GetBinError(i); 
      v+=v_1;
      v_err+=v_1err*v_1err;    }
    else 
      std::cout<<"no more bins for hist h_1!"<<std::endl;
    
    float v_2=0., v_2err=0.;
    if(i <= h_2->GetNbinsX() ) {
      v_2=h_2->GetBinContent(i);
      v_2err=h_2->GetBinError(i); 
      v+=v_2;
      v_err+=v_2err*v_2err;      }
    else
      std::cout<<"no more bins for hist h_2!"<<std::endl;
    
    float v_3=0., v_3err=0.;
    if(i <= h_3->GetNbinsX() ) {
      v_3=h_3->GetBinContent(i);
      v_3err=h_3->GetBinError(i); 
      v+=v_3;
      v_err+=v_3err*v_3err;        }
    else 
      std::cout<<"no more bins for hist h_3!"<<std::endl;
    
    float v_4=0., v_4err=0.;
    if(i <= h_4->GetNbinsX() ) {
      v_4=h_4->GetBinContent(i);
      v_4err=h_4->GetBinError(i); 
      v+=v_4;
      v_err+=v_4err*v_4err;        }
    else 
      std::cout<<"no more bins for hist h_4!"<<std::endl;
    
    v/=4.;//for etabin width    
    //v/=1000.;//picobarns to nanobarns
    
    v_err=sqrt(v_err); //for err
    v_err/=4.;
    //v_err/=1000.;//picobarns to nanobarns
    
    if(funcDebug)std::cout<<"v="<<v<<std::endl;
    if(funcDebug)std::cout<<"v_err="<<v_err<<std::endl;
    if(funcDebug)std::cout<<"setting bin content..."<<std::endl;
    hfinal_00eta20->SetBinContent(i, v);
    hfinal_00eta20->SetBinError(i, v_err);
    if(funcDebug)std::cout<<"histogram contents..."<<std::endl;
    if(funcDebug)std::cout<<"bincontent="<<hfinal_00eta20->GetBinContent(i)<<std::endl;
    if(funcDebug)std::cout<<"binerr="<<hfinal_00eta20->GetBinError(i)<<std::endl;
    if(funcDebug)std::cout<<std::endl;
  }
  
  return hfinal_00eta20;
}





void makeToySpectra(TH1D* hthy, 
		    //		    TSpline3* hthy_spline, 
		    TObject* weights, 
		    TF1* fJER, 
		    int nevts, 
		    TH1D* hthy_toyMC, TH1D* hsmeared_toyMC, TH1D* hsmeared_toyMC_test, 
		    TH2D* resp , TF1* spline3ext=NULL
		    //		    int etabin=-1, int JERsysdir=0
		    ){
  
  bool funcDebug=false;
  int tenth_nEvents=nevts/10;
  TSpline3* spl3weights=NULL;
  TF1* fitweights=NULL;
  bool useFitWeights=false, useSplineWeights=false, useSplineExt=false;
  if( weights->InheritsFrom("TSpline3") ) {
    std::cout<<"using spline weights"<<std::endl;
    spl3weights=(TSpline3*)weights->Clone(((std::string)weights->GetName()+"_funcclone").c_str() );
    std::cout<<"weights object name = "<< spl3weights->GetName()<<std::endl;
    useSplineWeights=true;
    if((bool)spline3ext) useSplineExt=true;    
  }
  else if(weights->InheritsFrom("TF1") ) {
    std::cout<<"using fit weights"<<std::endl;
    fitweights=(TF1*)weights->Clone(((std::string)weights->GetName()+"_funcclone").c_str() );    
    std::cout<<"weights object name = "<< fitweights->GetName()<<std::endl;
    useFitWeights=true;
  }
  
  //if both sets of weights are there, or if both sets are not there, shut everything down, we have a weird scenario on our hands...
  if( useFitWeights == useSplineWeights ) assert(false);
  
  srand((unsigned)time(0));  
  UInt_t rnd_seed=rand();
  if(funcDebug)std::cout<<"rnd_seed="<<rnd_seed<<std::endl;
  TRandom3 *rnd = new TRandom3(rnd_seed);
  
  //srand((unsigned)time(0));
  UInt_t rnd_test_seed=rand();
  if(funcDebug)std::cout<<"rnd_test_seed="<<rnd_test_seed<<std::endl;
  TRandom3 *rnd_test = new TRandom3(rnd_test_seed);
 
  //assert(false);          
  int nbins=hthy->GetNbinsX();
  double ptmin_thy=hthy->GetBinLowEdge(1);
  double ptmax_thy=hthy->GetBinLowEdge(nbins) + hthy->GetBinWidth(nbins);
  double ptmin_smeared=ptmin_thy, ptmax_smeared=ptmax_thy;  
  int respcount=0, misscount=0, fakecount=0;
  
  //  double spline3_maxptval=hthy->GetBinLowEdge(nbins)+0.5*hthy->GetBinWidth(nbins);
  double spline3_ext_ptmin=-1., spline3_ext_ptmax=-1.;
  if(useSplineExt){
    spline3ext->GetRange( spline3_ext_ptmin, spline3_ext_ptmax);
    std::cout<<"using spline3 extension for p_T range " << (spline3_ext_ptmin) << " GeV - " << (spline3_ext_ptmax) << std::endl;
    if(spline3_ext_ptmin<0. || spline3_ext_ptmax<0.)assert(false);
  }
  
  //bool doJERscalefact=(etabin>=0), doJERsys=(JERsysdir!=0);
  //double JERscaleFactor=1.;
  //if(doJERscalefact){    
  //  JERscaleFactor=JERscaleFacts[etabin];    
  //  if(doJERsys){
  //    if(abs(JERsysdir)!=1){std::cout<<"err, JERsysdir = -1, 0, or 1 only."<<std::endl; assert(false);}
  //    JERscaleFactor+=JERsysdir*JERscaleFactErrs[etabin];      }  }
  
  
  for(int i=0;i<nevts;++i){      

    if(i%tenth_nEvents==0)
      std::cout<<"throwing random #'s for event # "<<i<<std::endl;
    
    double ptTrue  = rnd->Uniform(ptmin_thy,ptmax_thy);
    double sigma   = fJER->Eval(ptTrue);    //sigma*=1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment
    //sigma*=JERscaleFactor;
    double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
    
    double pt_w =0.;//     =  hthy_spline->Eval(ptTrue);
    if(useSplineWeights){
      pt_w=spl3weights->Eval(ptTrue);
      if(useSplineExt)
	{
	  if(ptTrue>spline3_ext_ptmin &&
	     ptTrue<spline3_ext_ptmax ){
	    //std::cout<<"pt_w="<<pt_w<<std::endl;
	    pt_w=spline3ext->Eval(ptTrue);      
	    //std::cout<<"pt_w="<<pt_w<<std::endl;
	  }	  
	}
    }
    else if(useFitWeights){
      pt_w=fitweights->Eval(ptTrue);    }
    
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
    //sigma_test*=JERscaleFactor;
    double ptSmeared_test =  rnd_test->Gaus(ptTrue_test,ptTrue_test*sigma_test);
    
    double      pt_w_test =  0.;//hthy_spline->Eval(ptTrue_test);
    if(useSplineWeights){
      pt_w_test=spl3weights->Eval(ptTrue_test);
      if(useSplineExt){
	if(ptTrue_test>spline3_ext_ptmin && 
	   ptTrue_test<spline3_ext_ptmax )
	  pt_w_test=spline3ext->Eval(ptTrue_test);}      }
    else if((bool)fitweights){
      pt_w_test=fitweights->Eval(ptTrue_test);    }
    
    bool in_test_smearpt_range=(  ( ptSmeared_test>ptmin_smeared ) && ( ptSmeared_test<ptmax_smeared )  );
    if(in_test_smearpt_range)
      hsmeared_toyMC_test->Fill(ptSmeared_test,pt_w_test);
    
    
  }//end loop throwing random no's

  std::cout<<"done smearing ynew."<<std::endl;  
  if(funcDebug)std::cout<<"smear summary;"<<std::endl;
  if(funcDebug)std::cout<<"nEvents="<<nevts<<std::endl;
  if(funcDebug)std::cout<<"response_count="<<respcount<<std::endl;
  if(funcDebug)std::cout<<"miss_count="<<misscount<<std::endl;
  if(funcDebug)std::cout<<"fake_count="<<fakecount<<std::endl;
  
  return;
}

double calc_spline3ext_y2(TH1D* hthy){
  bool funcDebug=true;
  int nbins=hthy->GetNbinsX();
  if(funcDebug)std::cout<<"nbins="<< nbins<< std::endl;
  double y_lastbin=hthy->GetBinContent(nbins);
  double y_nxt2lastbin=hthy->GetBinContent(nbins-1);
  if(funcDebug)std::cout<<"y_lastbin     ="<< y_lastbin     << std::endl;
  if(funcDebug)std::cout<<"y_nxt2lastbin ="<< y_nxt2lastbin << std::endl;
  
  double log_y_lastbin=TMath::Log10(y_lastbin);
  double log_y_nxt2lastbin=TMath::Log10(y_nxt2lastbin);
  if(funcDebug)std::cout<<"log_y_lastbin     ="<< log_y_lastbin     << std::endl;
  if(funcDebug)std::cout<<"log_y_nxt2lastbin ="<< log_y_nxt2lastbin << std::endl;

  double binwidth_scalefact=hthy->GetBinWidth(nbins)/hthy->GetBinWidth(nbins-1);
  double logdiff=log_y_nxt2lastbin-log_y_lastbin;
  //double logdiff_ov2=logdiff/2.;
  if(funcDebug)std::cout<<"logdiff     ="<< logdiff     << std::endl;
  if(funcDebug)std::cout<<"binwidth_scalefact     ="<< binwidth_scalefact     << std::endl;
  //if(funcDebug)std::cout<<"logdiff_ov2 ="<< logdiff_ov2 << std::endl;
  
  
  //double y2_log=log_y_lastbin-logdiff_ov2;  
  //double y2_log=log_y_lastbin-16.0*logdiff;  
  double y2_log=log_y_lastbin-8.0*logdiff;  
  //double y2_log=log_y_lastbin-(binwidth_scalefact*logdiff);  
  double y2=TMath::Power(10., y2_log);
  if(funcDebug)std::cout<<"y2_log ="<< y2_log<< std::endl;
  if(funcDebug)std::cout<<"y2     ="<< y2<< std::endl;  
  
  return y2;
}








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

void makeModLogFits(TH1D* hthy, 
		    TF1* logFit    =NULL , 
		    TF1* modLogFit =NULL , 
		    TF1* modLog2Fit=NULL , 
		    TF1* modLog3Fit=NULL ,
		    TF1* modLog4Fit=NULL ,
		    TF1* modLog5Fit=NULL 
		    ){
  if( !((bool)logFit) || !((bool)modLogFit) || !((bool)modLog2Fit)){
    std::cout<<"don't do this fit unless you plan on going to at least third order, bare minimum! exit."<<std::endl;
    assert(false);
  }

  
  bool funcDebug=false;
  std::string fitOpt="MRE";
  if(funcDebug)fitOpt+="V";
  //void makeModLogFits(TH1D* hthy, std::vector<TF1*> logfitarr){  
  // fit of real interest is; A(B/pt)^[C+ D Log(pt/B) + E(Log(pt/B))^2]
  
  const int nbins=hthy->GetNbinsX();
  const int parambin=nbins/2;
  //const int startbin=1, endbin=hthy->GetNbinsX();
  //const float xlo=hthy->GetBinLowEdge(startbin);
  //const float xhi=hthy->GetBinLowEdge(endbin) + hthy->GetBinWidth(endbin);
  
  TH1D* hthy_log    =(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_logfit").c_str());
  TH1D* hthy_modlog =(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_modlogfit").c_str());   
  TH1D* hthy_modlog2=(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_modlog2fit").c_str());
  TH1D* hthy_modlog3=(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_modlog3fit").c_str());
  TH1D* hthy_modlog4=(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_modlog4fit").c_str());
  TH1D* hthy_modlog5=(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_modlog5fit").c_str());
  
  float logfit_p0=hthy->GetBinContent(parambin); //std::cout<<"fit param 0 being set to="<<logFitp0<<std::endl;  //parameter 0 val of xsec
  float logfit_p1=hthy->GetBinCenter(parambin);  //std::cout<<"fit param 1 being set to="<<logFitp1<<std::endl;  //parameter 1 val bin center of param 0's bin
  //basically, using my hist, in logx/logy coordinates, this is a linear fit
  //TF1* logFit=new TF1("logfit",    "[0]*( TMath::Power( ( [1]/x ) , [2] ) )",xlo,xhi);
  logFit->SetParameter(0,logfit_p0);
  logFit->SetParameter(1,logfit_p1);
  hthy_log->Fit(logFit->GetName(),fitOpt.c_str());  
  logFit->SetLineColor(kRed);
  
  float modlogfit_p0=logFit->GetParameter(0);
  float modlogfit_p1=logFit->GetParameter(1);
  float modlogfit_p2=logFit->GetParameter(2);
  //thus, this is a parabolic fit in logx/logy
  //TF1* modLogFit=new TF1("modlogfit","[0]*( TMath::Power( ( [1]/x ) , [2]+[3]*TMath::Log10( x/[1] ) ) )",xlo,xhi);
  modLogFit->SetParameter(0,modlogfit_p0);
  modLogFit->SetParameter(1,modlogfit_p1);
  modLogFit->SetParameter(2,modlogfit_p2);
  hthy_modlog->Fit(modLogFit->GetName(),fitOpt.c_str());
  modLogFit->SetLineColor(kBlue);
  
  float modlog2fit_p0=modLogFit->GetParameter(0);
  float modlog2fit_p1=modLogFit->GetParameter(1);
  float modlog2fit_p2=modLogFit->GetParameter(2);
  float modlog2fit_p3=modLogFit->GetParameter(3);
  //thus, this is a cubic fit fit in logx/logy
  //TF1* modLog2Fit=new TF1("modlog2fit", "[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
  modLog2Fit->SetParameter(0,modlog2fit_p0);
  modLog2Fit->SetParameter(1,modlog2fit_p1);
  modLog2Fit->SetParameter(2,modlog2fit_p2);
  modLog2Fit->SetParameter(3,modlog2fit_p3);
  hthy_modlog2->Fit(modLog2Fit->GetName(),fitOpt.c_str());
  modLog2Fit->SetLineColor(kMagenta);


  if((bool)modLog3Fit){
    float modlog3fit_p0=modLog2Fit->GetParameter(0);
    float modlog3fit_p1=modLog2Fit->GetParameter(1);
    float modlog3fit_p2=modLog2Fit->GetParameter(2);
    float modlog3fit_p3=modLog2Fit->GetParameter(3);
    float modlog3fit_p4=modLog2Fit->GetParameter(4);
    //thus, this is a cubic fit fit in logx/logy
    //TF1* modLog3Fit=new TF1("modlog2fit", "[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
    modLog3Fit->SetParameter(0,modlog3fit_p0);
    modLog3Fit->SetParameter(1,modlog3fit_p1);
    modLog3Fit->SetParameter(2,modlog3fit_p2);
    modLog3Fit->SetParameter(3,modlog3fit_p3);
    modLog3Fit->SetParameter(4,modlog3fit_p4);
    hthy_modlog3->Fit(modLog3Fit->GetName(),fitOpt.c_str());
    modLog3Fit->SetLineColor(kGreen+2);
  }

  if((bool)modLog4Fit){
    float modlog4fit_p0=modLog3Fit->GetParameter(0);
    float modlog4fit_p1=modLog3Fit->GetParameter(1);
    float modlog4fit_p2=modLog3Fit->GetParameter(2);
    float modlog4fit_p3=modLog3Fit->GetParameter(3);
    float modlog4fit_p4=modLog3Fit->GetParameter(4);
    float modlog4fit_p5=modLog3Fit->GetParameter(5);
    //thus, this is a cubic fit fit in logx/logy
    //TF1* modLog3Fit=new TF1("modlog2fit", "[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
    modLog4Fit->SetParameter(0,modlog4fit_p0);
    modLog4Fit->SetParameter(1,modlog4fit_p1);
    modLog4Fit->SetParameter(2,modlog4fit_p2);
    modLog4Fit->SetParameter(3,modlog4fit_p3);
    modLog4Fit->SetParameter(4,modlog4fit_p4);
    modLog4Fit->SetParameter(5,modlog4fit_p5);
    hthy_modlog4->Fit(modLog4Fit->GetName(),fitOpt.c_str());
    modLog4Fit->SetLineColor(kTeal);
  }

  if((bool)modLog5Fit){
    float modlog5fit_p0=modLog4Fit->GetParameter(0);
    float modlog5fit_p1=modLog4Fit->GetParameter(1);
    float modlog5fit_p2=modLog4Fit->GetParameter(2);
    float modlog5fit_p3=modLog4Fit->GetParameter(3);
    float modlog5fit_p4=modLog4Fit->GetParameter(4);
    float modlog5fit_p5=modLog4Fit->GetParameter(5);
    float modlog5fit_p6=modLog4Fit->GetParameter(6);
    //thus, this is a cubic fit fit in logx/logy
    //TF1* modLog3Fit=new TF1("modlog2fit", "[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
    modLog5Fit->SetParameter(0,modlog5fit_p0);
    modLog5Fit->SetParameter(1,modlog5fit_p1);
    modLog5Fit->SetParameter(2,modlog5fit_p2);
    modLog5Fit->SetParameter(3,modlog5fit_p3);
    modLog5Fit->SetParameter(4,modlog5fit_p4);
    modLog5Fit->SetParameter(5,modlog5fit_p5);
    modLog5Fit->SetParameter(6,modlog5fit_p6);
    hthy_modlog5->Fit(modLog5Fit->GetName(),fitOpt.c_str());
    modLog5Fit->SetLineColor(kOrange);
  }
  std::cout<<"    (bool)logFit = "<<    (bool)logFit<<std::endl;
  std::cout<<" (bool)modLogFit = "<< (bool)modLogFit<<std::endl;
  std::cout<<"(bool)modLog2Fit = "<<(bool)modLog2Fit<<std::endl;
  std::cout<<"(bool)modLog3Fit = "<<(bool)modLog3Fit<<std::endl;
  std::cout<<"(bool)modLog4Fit = "<<(bool)modLog4Fit<<std::endl;
  std::cout<<"(bool)modLog5Fit = "<<(bool)modLog5Fit<<std::endl;
  


  return;

}






void make7TeVFits(TH1D* hthy=NULL, 
		  TF1* _FitA =NULL , 
		  TF1* _FitB =NULL ,
		  float norm=-1.
		  ){
  
  if( !((bool)_FitA) || !((bool)_FitB) || !((bool)hthy)){
    std::cout<<"Fits or hist not passed correclty to make7TeVFits, exit"<<std::endl;
    assert(false);
  }  
  bool funcDebug=false;
  std::string fitOpt="MRE";
  if(funcDebug)fitOpt+="V";
  
  //const int nbins=hthy->GetNbinsX();
 //const int parambin=nbins/2;
  //const int startbin=1, endbin=hthy->GetNbinsX();
  //const float xlo=hthy->GetBinLowEdge(startbin);
  //const float xhi=hthy->GetBinLowEdge(endbin) + hthy->GetBinWidth(endbin);
  
  
  TH1D* hthy_fitB=(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_7tevfitB").c_str());  
  hthy_fitB->Fit(_FitB->GetName(),fitOpt.c_str());//do fitB first, it has 
  _FitB->SetLineColor(kBlue);
  

  TH1D* hthy_fitA=(TH1D*)hthy->Clone( ( ((std::string)hthy->GetName())+"_7tevfitA").c_str());   
  
  float fitA_p0=_FitB->GetParameter(0);//take in the three param from fitB [norm isnt free param, preset] as suggestions
  float fitA_p1=_FitB->GetParameter(1);
  float fitA_p2=_FitB->GetParameter(2);
  
  _FitA->SetParameter(0,fitA_p0);//set the free params
  _FitA->SetParameter(1,fitA_p1);
  _FitA->SetParameter(2,fitA_p2);
  if(norm>0.){    
    _FitA->SetParameter(3,norm);
    _FitA->SetParLimits(3,.5*norm,2.0*norm);
  }
  
  hthy_fitA->Fit(_FitA->GetName(),fitOpt.c_str());
  _FitA->SetLineColor(kMagenta);
  
  std::cout << "(bool) _FitA = " << (bool)_FitA            <<std::endl;
  std::cout<<"_FitA par0 (alpha) ="<<_FitA->GetParameter(0)<<std::endl;
  std::cout<<"_FitA par1 (beta)  ="<<_FitA->GetParameter(1)<<std::endl;
  std::cout<<"_FitA par2 (gamma) ="<<_FitA->GetParameter(2)<<std::endl;
  std::cout<<"_FitA par3 (N0)    ="<<_FitA->GetParameter(3)<<std::endl;
  
  std::cout << "(bool) _FitB = " << (bool)_FitB            <<std::endl;
  std::cout<<"_FitA par0 (alpha) ="<<_FitB->GetParameter(0)<<std::endl;
  std::cout<<"_FitA par1 (beta)  ="<<_FitB->GetParameter(1)<<std::endl;
  std::cout<<"_FitA par2 (gamma) ="<<_FitB->GetParameter(2)<<std::endl;
  std::cout<<"norm set to        ="<<norm                  <<std::endl;


  //assert(false);
  return;

}


//shift PDF spectra up/down by 6or2 pt scale only
TH1D* makePDFsys(TH1D* thy, TH1D* thyerrs, std::string updown_str="sysup", double errfact=1.0){
  TH1D* thy_PDFsys=(TH1D*)thy->Clone(("thy_PDFsyshist_NPterr_"+updown_str).c_str()); 
  int nbinsx=thy->GetNbinsX();
  int nbinsx_err=thyerrs->GetNbinsX();
  if(nbinsx!=nbinsx_err){
    std::cout<<"error, # of bins mismatched between thy and thyerrs hist. exit."<<std::endl;    assert(false);  }
  
  if(updown_str=="sysup" || updown_str=="sysdown"){
    for(int i=1; i<=nbinsx;i++){
      double err=thy->GetBinContent(i)*(thyerrs->GetBinContent(i)/100.);
      //thy_PDFsys->SetBinContent(i, thy->GetBinContent(i)+thyerrs->GetBinContent(i));//WRONG; errors are RELATIVE in fastNLO files           
      //thy_PDFsys->SetBinContent(i, thy->GetBinContent(i)+errfact*thy->GetBinContent(i)*(thyerrs->GetBinContent(i)/100.) );//scaled to PERCENT in my code now
      thy_PDFsys->SetBinContent( i, thy->GetBinContent(i) + errfact*err );//scaled to PERCENT in my code now
    }      
  }
  else { 
    std::cout<<"error, updown_str not recognized. exit." <<std::endl; assert(false);}
  return (TH1D*)thy_PDFsys;
}


//shit PDF spectra up/down by 6or2 pt scale unc + PDF unc added in quad.
TH1D* makePDFsys_wPDFerr(TH1D* thy, TH1D* thyerrs, std::string updown_str="sysup", double errfact=1.0){
  TH1D* thy_PDFsys=(TH1D*)thy->Clone(("thy_PDFsyshist_NPterr_wPDFerr_"+updown_str).c_str()); // PDF errs on thy hist set CORRECTLY now
  int nbinsx=thy->GetNbinsX();
  int nbinsx_err=thyerrs->GetNbinsX();
  if(nbinsx!=nbinsx_err){
    std::cout<<"error, # of bins mismatched between thy and thyerrs hist. exit."<<std::endl;    assert(false);  }
  
  if(updown_str=="sysup" || updown_str=="sysdown"){
    for(int i=1; i<=nbinsx;i++){
      double err=sqrt( 
		      thy->GetBinError(i)*thy->GetBinError(i) +  ( 
								  (thy->GetBinContent(i)*
								   thyerrs->GetBinContent(i)/
								   1.e+02) *
								  (thy->GetBinContent(i)*								   
								   thyerrs->GetBinContent(i)/
								   1.e+02  )
								   )
		       );
      if(      updown_str== "sysup" ) thy_PDFsys->SetBinContent(i, thy->GetBinContent(i)+errfact*err);//cause error will be positive due to squaring
      else if( updown_str=="sysdown") thy_PDFsys->SetBinContent(i, thy->GetBinContent(i)-errfact*err);
      else assert(false);
    }
  }
  else { 
    std::cout<<"error, updown_str not recognized. exit." <<std::endl; assert(false);}
  return (TH1D*)thy_PDFsys;
}








//int getRandomSeed(){
//  bool funcDebug=true;  
//  int
//}




//BACKUP
//void makeToySpectra(TH1D* hthy, 
//		    TSpline3* hthy_spline, 
//		    TF1* fJER, 
//		    int nevts, 
//		    TH1D* hthy_toyMC, TH1D* hsmeared_toyMC, TH1D* hsmeared_toyMC_test, 
//		    TH2D* resp ){
//  bool funcDebug=true;
//  int tenth_nEvents=nevts/10;
//  
//  srand((unsigned)time(0));  
//  UInt_t rnd_seed=rand();
//  if(funcDebug)std::cout<<"rnd_seed="<<rnd_seed<<std::endl;
//  TRandom3 *rnd = new TRandom3(rnd_seed);
//  
//  //srand((unsigned)time(0));
//  UInt_t rnd_test_seed=rand();
//  if(funcDebug)std::cout<<"rnd_test_seed="<<rnd_test_seed<<std::endl;
//  TRandom3 *rnd_test = new TRandom3(rnd_test_seed);
// 
//  //assert(false);          
//  double ptmin_thy=hthy->GetBinLowEdge(1);
//  double ptmax_thy=hthy->GetBinLowEdge(hthy->GetNbinsX()) + hthy->GetBinWidth(hthy->GetNbinsX());
//  double ptmin_smeared=ptmin_thy, ptmax_smeared=ptmax_thy;
//  
//  int respcount=0, misscount=0, fakecount=0;
//
//  for(int i=0;i<nevts;++i){      
//    
//
//
//    if(i%tenth_nEvents==0)
//      std::cout<<"throwing random #'s for event # "<<i<<std::endl;
//
//    double ptTrue  = rnd->Uniform(ptmin_thy,ptmax_thy);
//
//    double sigma   = fJER->Eval(ptTrue);    //sigma*=1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment
//
//    double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
//
//    double pt_w      =  hthy_spline->Eval(ptTrue);
//
//    bool in_smearpt_range=(  ( ptSmeared>ptmin_smeared ) && ( ptSmeared<ptmax_smeared )  );
//    bool in_trupt_range=(  ( ptTrue>ptmin_thy )        && ( ptTrue<ptmax_thy )  )  ;
//    bool fillresp=in_smearpt_range&&in_trupt_range;
//
//
//    if(fillresp){
//      resp->Fill(ptSmeared,ptTrue,pt_w);
//      hthy_toyMC    -> Fill(ptTrue,pt_w);
//      hsmeared_toyMC -> Fill(ptSmeared,pt_w);      
//      respcount++;
//    }
//    else if(in_trupt_range){//not in smear pt range but in truept range
//      //      resp->Miss(ptTrue,pt_w);                                        
//      hthy_toyMC    -> Fill(ptTrue,pt_w);
//      misscount++;
//    }
//    else if(in_smearpt_range){//shouldnt get filled by construction. ptTrue always pulled s.t. in_trupt_range is true. if in_smearpt_range is true, so it in_trupt_range
//      //      resp->Fake(ptSmeared,pt_w);                                        
//      hsmeared_toyMC    -> Fill(ptSmeared,pt_w);
//      fakecount++;
//    }
//    
//    double    ptTrue_test = rnd_test->Uniform(ptmin_thy,ptmax_thy);
//    double     sigma_test = fJER->Eval(ptTrue_test);    //sigma*=1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment
//    double ptSmeared_test =  rnd_test->Gaus(ptTrue_test,ptTrue_test*sigma_test);
//    double      pt_w_test =  hthy_spline->Eval(ptTrue_test);
//    
//    bool in_test_smearpt_range=(  ( ptSmeared_test>ptmin_smeared ) && ( ptSmeared_test<ptmax_smeared )  );
//    if(in_test_smearpt_range)
//      hsmeared_toyMC_test->Fill(ptSmeared_test,pt_w_test);
//    
//
//  }//end loop throwing random no's
//
//  std::cout<<"done smearing ynew."<<std::endl;  
//  std::cout<<"smear summary;"<<std::endl;
//  std::cout<<"nEvents="<<nevts<<std::endl;
//  std::cout<<"response_count="<<respcount<<std::endl;
//  std::cout<<"miss_count="<<misscount<<std::endl;
//  std::cout<<"fake_count="<<fakecount<<std::endl;
//  
//  return;
//}
//BACKUP


//void smear_theory_spectra(TH1D* hthy, TH1D* hsmrd, RooUnfoldResponse respmat, TSpline3* spline, TF1* fJER){
// 
//  bool funcDebug=true;
//  
//  return;
//}

