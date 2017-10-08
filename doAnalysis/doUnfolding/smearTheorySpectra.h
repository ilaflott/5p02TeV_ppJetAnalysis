#include <iostream>
#include <cmath>
#include <cassert>

#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>

#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <TSpline.h>

#include "RooUnfoldResponse.h"

using namespace std;

//TH1D *cross_section(TH1D *Spectra_);

TH1D* cross_section(TH1D* Spectra_){
  
  /// Evaluate xSections
  TH1D *Xsection_=(TH1D*)Spectra_->Clone();
  for(int i=0; i<Spectra_->GetNbinsX(); ++i){
    float bin_value = Spectra_->GetBinContent(1+i);
    float bin_error = Spectra_->GetBinError(1+i);
    float bin_width = Spectra_->GetBinWidth(1+i);
    float new_bin_value = bin_value/bin_width;
    float new_bin_error = bin_error/bin_width;
    Xsection_->SetBinContent(1+i,new_bin_value);
    Xsection_->SetBinError(1+i,new_bin_error);
  }
  //   std::cout<<Integral_<<std::endl;
  
  return Xsection_;
}

/////////////// Inclusive jets Binning form Ksenia ////////////////////////////////////////////////////////
double thyBins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  //  32., 37., //junk bins above 
  43., 49., 56., 64., 74., 84.,
  97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362., 395., 430., 468.,
  507., 548., 592., 638., 686., 
  //737., 790., 846., 905., 967 ., 1032. //this line; old bins pre 1000
  1000., 
  1500.
  //  1032., 1101. //junk from here down
}; 
const int n_thybins_incl=sizeof(thyBins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1


double smearedBins_incl[]={
  //0., 1., 5., 6., 8., 10., 12., 15., 18., 21., 24.,   // junk
  56., 64., 74., 84.,
  97., 114., 133., 153., 174., 196., 220., 245., 272., 300., 330., 362., 395., 430., 468.,
  507., 548., 592., 638., 686., 
  //737., 790., 846., 905., 967 ., 1032. //this line; old bins pre 1000
  1000.
  //  1032., 1101. //junk from here down
}; 
const int n_smearedbins_incl=sizeof(smearedBins_incl)/sizeof(double)-1; //this is # of bins = # of entries in array - 1


////////// NPs from Paolo 4/9/2015
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
  return (  (  NP_y0(x) + NP_y1(x) + NP_y2(x) + NP_y3(x) )/4. ); 
}
string NP_ynew_str=NP_y0_str+" + "+NP_y1_str+" + "+NP_y2_str+" + "+NP_y3_str;




////////// Gaussian Core Distributions from Ian L 10/3/17
// from JER fits using, sigma/mu = [0] + [1] /( pow (gen p_T , [2] ) + [3] * gen p_t)
double gres_ynew(double x){
  return ( 4.80963e-02 + 1.51797e-02/(pow(x , -5.78569e-01)+(2.17736e-03*x)) );}
string gres_ynew_str="4.80963e-02+1.51797e-02/(pow(x,-5.78569e-01)+(2.17736e-03*x))";



TH1D* make00eta20Hist( TH1D * h_1, 
		       TH1D * h_2, 
		       TH1D * h_3, 
		       TH1D * h_4  ){
  
  bool funcDebug=true;
  
  //double lowestBinEdge[4]={};
  //lowestBinEdge[0]=h_1->GetXaxis()->GetBinLowEdge(1);
  //lowestBinEdge[1]=h_2->GetXaxis()->GetBinLowEdge(1);
  //lowestBinEdge[2]=h_3->GetXaxis()->GetBinLowEdge(1);
  //lowestBinEdge[3]=h_4->GetXaxis()->GetBinLowEdge(1);
  //
  //double highestBinEdge[4]={};
  //highestBinEdge[0]=h_1->GetXaxis()->GetBinLowEdge(h_1->GetNbinsX())    +h_1->GetBinWidth(h_1->GetNbinsX());
  //highestBinEdge[1]=h_2->GetXaxis()->GetBinLowEdge(h_2->GetNbinsX())    +h_2->GetBinWidth(h_2->GetNbinsX());
  //highestBinEdge[2]=h_3->GetXaxis()->GetBinLowEdge(h_3->GetNbinsX())    +h_3->GetBinWidth(h_3->GetNbinsX());
  //highestBinEdge[3]=h_4->GetXaxis()->GetBinLowEdge(h_4->GetNbinsX())    +h_4->GetBinWidth(h_4->GetNbinsX());
  //
  //double min_lowestBinEdge   =+9999.;
  //for(int i=0; i<3 ; i++) if( lowestBinEdge[i] < min_lowestBinEdge ) min_lowestBinEdge=lowestBinEdge[i];
  //if(funcDebug)cout<<"min lowest bin Edge = "<< min_lowestBinEdge << endl;
  //
  //double max_highestBinEdge  = -1. ;
  //for(int i=0; i< 3; i++) if( highestBinEdge[i] > max_highestBinEdge ) max_highestBinEdge=highestBinEdge[i];
  //if(funcDebug)cout<<"max highest bin Edge = "<< max_highestBinEdge << endl;
  //cout<<endl;
  //
  //int numEntries00eta20=0; //#of entries in new ptbinning array
  //int startingEntry=-1, endingEntry=-1;
  //for(int i=0; i<=n_thybins_incl; i++) {
  //  if(thyBins_incl[i] == min_lowestBinEdge) {
  //    if(funcDebug)cout<<"lowest bin edge found!"<<endl;
  //    startingEntry=i;
  //    if(funcDebug)cout<<"starting entry="<<startingEntry<<endl;
  //    numEntries00eta20++;
  //    if(funcDebug)cout<<"thyBins_incl["<<i<<"]= "<<thyBins_incl[i]<<endl;
  //    if(funcDebug)cout<<"numEntries00eta20 now = "<<numEntries00eta20<<endl;
  //  }
  //  //if(thyBins_incl[i] == max_highestBinEdge) endingEntry=i;
  //  else if( thyBins_incl[i] > min_lowestBinEdge && 
  //	thyBins_incl[i] <= max_highestBinEdge ) {
  //
  //    if(funcDebug)cout<<"pt val inbetween min/max bin edges!"<<endl;
  //    if(funcDebug)cout<<"thyBins_incl["<<i<<"]= "<<thyBins_incl[i]<<endl;
  //    numEntries00eta20++;       
  //    if(funcDebug)cout<<"numEntries00eta20 now = "<<numEntries00eta20<<endl;
  //  }
  //}
  //endingEntry=startingEntry+numEntries00eta20;
  ////numEntries00eta20++;//need +1 here because for some reason <= max_highestBinEdge in above loop fucks up
  //
  //if(funcDebug)cout<<endl;
  //if(funcDebug)cout<<"numEntries 00eta20 = "<< numEntries00eta20 << endl;
  //if(funcDebug)cout<<"startingEntry = "<< startingEntry << endl;
  //if(funcDebug)cout<<"ptbins val =" << thyBins_incl[startingEntry] << endl;
  //
  //if(funcDebug)cout<<"endingEntry = "<< endingEntry << endl;
  //if(funcDebug)cout<<"ptbins val =" << thyBins_incl[endingEntry] << endl;
  //
  //if(funcDebug)cout<<endl<<"creating new array..."<<endl<<endl;
  //const int numEntries=numEntries00eta20;
  //double thyBins_incl_00eta20[numEntries];
  ////const int numBins=numEntries-1;
  //const int numBins=numEntries-1;
  //int newBinsEntry=0;
  ////for(int i=startingEntry; i<=endingEntry; i++) {
  //for(int i=0; i<n_thybins_incl; i++) {
  //  
  //
  //  if(newBinsEntry>=(numEntries) ){
  //    if(funcDebug)cout<<"warning, about to go off the end of the new binning array..."<<endl;
  //    break;
  //  }
  //  
  //  
  //  
  //  if(funcDebug)cout << "i=" << i << endl;
  //  if(funcDebug)cout << "thyBins_incl["<<i<<"] = "<< thyBins_incl[i] << endl;    
  //  
  //  thyBins_incl_00eta20[newBinsEntry]=thyBins_incl[i];    
  //  if(funcDebug)cout << "newBinsEntry=" << newBinsEntry << endl;
  //  if(funcDebug)cout << "thyBins_incl_00eta20["<<newBinsEntry<<"] = "<< thyBins_incl_00eta20[newBinsEntry] << endl;    
  //  
  //  newBinsEntry++;
  //}
  
  //TH1D* h00eta20 = new TH1D( "h00eta20","thyHist 00eta20", numBins, (const double*) thyBins_incl_00eta20 );
  TH1D* h00eta20 = new TH1D( "h00eta20","thyHist 00eta20", n_thybins_incl, (const double*) thyBins_incl );
  h00eta20->Print("base");
  
  //this part assumes all hists start at same lower limit; which they better!
  if(funcDebug)cout<<endl;
  if(funcDebug)cout<<"filling new hist!"<<endl;
  if(funcDebug)cout<<endl;
  
  for(int i=1; i<=n_thybins_incl; i++) {
    
    if(funcDebug)cout<<"i="<<i<<endl;
    cout<<"lower edge of bin i for h00eta20="<<h00eta20->GetXaxis()->GetBinLowEdge(i)<<endl;
    cout<<"lower edge of bin i for h_1="<<h_1->GetXaxis()->GetBinLowEdge(i)<<endl;
    float v=0., v_err=0.;
    
    float v_1=0., v_1err=0.;
    if(i <= h_1->GetNbinsX() ) {
      v_1=h_1->GetBinContent(i);
      v_1err=h_1->GetBinError(i); 
      v+=v_1;
      v_err+=v_1err*v_1err;    }
    else{ cout<<"no more bins for hist h_1!"<<endl;}

    float v_2=0., v_2err=0.;
    if(i <= h_2->GetNbinsX() ) {
      v_2=h_2->GetBinContent(i);
      v_2err=h_2->GetBinError(i); 
      v+=v_2;
      v_err+=v_2err*v_2err;  
    }
    else{ cout<<"no more bins for hist h_2!"<<endl;}

    float v_3=0., v_3err=0.;
    if(i <= h_3->GetNbinsX() ) {
      v_3=h_3->GetBinContent(i);
      v_3err=h_3->GetBinError(i); 
      v+=v_3;
      v_err+=v_3err*v_3err;    
    }
    else{ cout<<"no more bins for hist h_3!"<<endl;}

    float v_4=0., v_4err=0.;
    if(i <= h_4->GetNbinsX() ) {
      v_4=h_4->GetBinContent(i);
      v_4err=h_4->GetBinError(i); 
      v+=v_4;
      v_err+=v_4err*v_4err;    
    }
    else{ cout<<"no more bins for hist h_4!"<<endl;}
    
    v/=4.;//for etabin width    
    v_err=sqrt(v_err); //for err
    v_err/=4.;
    
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

















//#include <math.h>
//#include <Riostream.h>
//#include <Strlen.h>
//#include <TSystem.h>
//#include <TDirectory.h>
//#include <TApplication.h>
//#include <TInterpreter.h>
//#include <THashList.h>
//#include <TRegexp.h>
//#include <TClassTable.h>
//#include <TClass.h>
//#include <TVirtualMutex.h>
//#include <TError.h>




//////////// Gaussian Core Distributions from Panos 11/9/2015
//double gres_y0(double x){
//  return ( 2.56956e-02 + 1.09091e+00/( pow( x,5.74784e-01 )+( -2.82630e-03*x ) )  ); }
//
//double gres_y1(double x){
//  return ( 2.82669e-02 + 1.17381e+00/(pow(x ,5.95671e-01)+(-4.44816e-03*x)) );}
//
//double gres_y2(double x){
//  return ( 5.52461e-02 + 4.62664e+00/(pow(x ,9.19063e-01)+(2.86983e-02*x)) + (-5.90944e-06*x) );}
//
//double gres_y3(double x){
//  return ( 4.38816e-02 + 7.84121e+00/(pow(x ,9.88655e-01)+(2.42292e-01*x)) + (-6.79034e-06*x) );}
//
//double gres_y4(double x){
//  return ( 6.72175e-02 + 3.40566e+04/(pow(x ,2.76679e+00)+(5.72477e+03*x)) + (-2.95586e-05*x) );}
//
//double gres_y5(double x){
//  return ( 7.72268e-02 + 2.24271e+03/(pow(x ,2.25985e+00)+(2.25981e+02*x)) + (-5.08042e-05*x) );}
//
//double gres_y6(double x){
//  if(x<68.) {
//    return (-4.35645e-01 + 1.39255e+00/(pow(x ,2.51204e-01)+(-6.87795e-03*x)) );
//  }
//  else {   //if(x>=68){
//    return ( 1.18352e-01 + 1.99019e+00/(pow(x ,1.38662e+00)+(-2.46112e+02)) );
//  }
//}



  //TF1 *fNP_y0 = new TF1("fNP_y0","NP_y0(x)",32,6000);
  //TF1 *fNP_y1 = new TF1("fNP_y1","NP_y1(x)",32,5200);
  //TF1 *fNP_y2 = new TF1("fNP_y2","NP_y2(x)",32,4000);
  //TF1 *fNP_y3 = new TF1("fNP_y3","NP_y3(x)",32,2500);
  //TF1 *fNP_y4 = new TF1("fNP_y4","NP_y4(x)",32,1588);
  //TF1 *fNP_y5 = new TF1("fNP_y5","NP_y5(x)",32,967);
  //TF1 *fNP_y6 = new TF1("fNP_y6","NP_y6(x)",32,468);
  //plotNPs->cd(1);
  //fNP_y0->Draw();
  //
  //plotNPs->cd(2);
  //fNP_y1->Draw();
  //
  //plotNPs->cd(3);
  //fNP_y2->Draw();
  //
  //plotNPs->cd(4);
  //fNP_y3->Draw();
  //
  //plotNPs->cd(5);
  //fNP_y4->Draw();
  //
  //plotNPs->cd(6);
  //fNP_y5->Draw();
  //
  //plotNPs->cd(7);
  //fNP_y6->Draw();
  






  //TF1 *fres_y0 = new TF1("fres_y0","gres_y0(x)",56,3000);
  //TF1 *fres_y1 = new TF1("fres_y1","gres_y1(x)",56,3000);
  //TF1 *fres_y2 = new TF1("fres_y2","gres_y2(x)",56,2500);
  //TF1 *fres_y3 = new TF1("fres_y3","gres_y3(x)",56,1600);
  //TF1 *fres_y4 = new TF1("fres_y4","gres_y4(x)",56,1100);
  //TF1 *fres_y5 = new TF1("fres_y5","gres_y5(x)",56,600);
  //TF1 *fres_y6 = new TF1("fres_y6","gres_y6(x)",56,180);

  //fres_y0->SetLineColor(2);
  //fres_y0->SetMinimum(0); fres_y0->SetMaximum(0.25);
  //
  //fres_y1->SetLineColor(3);
  //fres_y2->SetLineColor(4);
  //fres_y3->SetLineColor(6);
  //fres_y4->SetLineColor(kOrange+1);
  //fres_y5->SetLineColor(kCyan-3);
  //fres_y6->SetLineColor(kYellow-5);
  //fres_y0->Draw();
  //fres_y1->Draw("same");
  //fres_y2->Draw("same");
  //fres_y3->Draw("same");
  //fres_y4->Draw("same");
  //fres_y5->Draw("same");
  //fres_y6->Draw("same");

  //  leg->AddEntry(fres_y0,"PYTHIA8 CUETM1 @ 13TeV", "");
  //  leg->AddEntry(fres_y0,"antikt7 and PFchs Jets", "");
  //  leg->AddEntry(fres_y0,"0<|y|<0.5", "l"); 
  //  leg->AddEntry(fres_y1,"0.5<|y|<1.0", "l"); 
  //  leg->AddEntry(fres_y2,"1.0<|y|<1.5", "l"); 
  //  leg->AddEntry(fres_y3,"1.5<|y|<2.0", "l"); 
  //  leg->AddEntry(fres_y4,"2.0<|y|<2.5", "l"); 
  //  leg->AddEntry(fres_y5,"2.5<|y|<3.0", "l"); 
  //  leg->AddEntry(fres_y6,"3.2<|y|<3.7", "l"); 
  


  ///// Cross section y0
  //TH1D *Xsection_theory_y0=(TH1D*)theory_y0->Clone();
  //for(int i=0; i<theory_y0->GetNbinsX(); ++i){
  //   double bin_value  = theory_y0->GetBinContent(1+i);
  //   double bin_center = theory_y0->GetBinCenter(1+i);
  //   double new_bin_value = bin_value * NP_y0(bin_center);
  //   Xsection_theory_y0->SetBinContent(1+i,new_bin_value);
  //}
  //
  ///// Cross section y1
  //TH1D *Xsection_theory_y1=(TH1D*)theory_y1->Clone();
  //for(int i=0; i<theory_y1->GetNbinsX(); ++i){
  //   double bin_value  = theory_y1->GetBinContent(1+i);
  //   double bin_center = theory_y1->GetBinCenter(1+i);
  //   double new_bin_value = bin_value * NP_y1(bin_center);
  //   Xsection_theory_y1->SetBinContent(1+i,new_bin_value);
  //}
  //
  ///// Cross section y2
  //TH1D *Xsection_theory_y2=(TH1D*)theory_y2->Clone();
  //for(int i=0; i<theory_y2->GetNbinsX(); ++i){
  //   double bin_value  = theory_y2->GetBinContent(1+i);
  //   double bin_center = theory_y2->GetBinCenter(1+i);
  //   double new_bin_value = bin_value * NP_y2(bin_center);
  //   Xsection_theory_y2->SetBinContent(1+i,new_bin_value);
  //}
  //
  ///// Cross section y3
  //TH1D *Xsection_theory_y3=(TH1D*)theory_y3->Clone();
  //for(int i=0; i<theory_y3->GetNbinsX(); ++i){
  //   double bin_value  = theory_y3->GetBinContent(1+i);
  //   double bin_center = theory_y3->GetBinCenter(1+i);
  //   double new_bin_value = bin_value * NP_y3(bin_center);
  //   Xsection_theory_y3->SetBinContent(1+i,new_bin_value);
  //}
  //
  ///// Cross section y4
  //TH1D *Xsection_theory_y4=(TH1D*)theory_y4->Clone();
  //for(int i=0; i<theory_y4->GetNbinsX(); ++i){
  //  double bin_value  = theory_y4->GetBinContent(1+i);
  //   double bin_center = theory_y4->GetBinCenter(1+i);
  //   double new_bin_value = bin_value * NP_y4(bin_center);
  //   Xsection_theory_y4->SetBinContent(1+i,new_bin_value);
  //}
  //
  ///// Cross section y5
  //TH1D *Xsection_theory_y5=(TH1D*)theory_y5->Clone();
  //for(int i=0; i<theory_y5->GetNbinsX(); ++i){
  //   double bin_value  = theory_y5->GetBinContent(1+i);
  //   double bin_center = theory_y5->GetBinCenter(1+i);
  //   double new_bin_value = bin_value * NP_y5(bin_center);
  //   Xsection_theory_y5->SetBinContent(1+i,new_bin_value);
  //}
  //
  ///// Cross section y6
  //TH1D *Xsection_theory_y6=(TH1D*)theory_y6->Clone();
  //for(int i=0; i<theory_y6->GetNbinsX(); ++i){
  //   double bin_value  = theory_y6->GetBinContent(1+i);
  //   double bin_center = theory_y6->GetBinCenter(1+i);
  //   double new_bin_value = bin_value * NP_y6(bin_center);
  //   Xsection_theory_y6->SetBinContent(1+i,new_bin_value);
  //}




  
  //TSpline3 *spline3_y0 = new TSpline3(Xsection_theory_y0);
  //TSpline3 *spline3_y1 = new TSpline3(Xsection_theory_y1);
  //TSpline3 *spline3_y2 = new TSpline3(Xsection_theory_y2);
  //TSpline3 *spline3_y3 = new TSpline3(Xsection_theory_y3);
  //TSpline3 *spline3_y4 = new TSpline3(Xsection_theory_y4);
  //TSpline3 *spline3_y5 = new TSpline3(Xsection_theory_y5);
  //TSpline3 *spline3_y6 = new TSpline3(Xsection_theory_y6);
  
  //  spline3_y0->SetLineColor(6);
  //  spline3_y1->SetLineColor(6);
  //  spline3_y2->SetLineColor(6);
  //  spline3_y3->SetLineColor(6);
  //  spline3_y4->SetLineColor(6);
  //  spline3_y5->SetLineColor(6);
  //  spline3_y6->SetLineColor(6);
  
  
  //Xsection_theory_y0->SetAxisRange(56,6000,"X"); 
  //Xsection_theory_y1->SetAxisRange(56,5200,"X"); 
  //Xsection_theory_y2->SetAxisRange(56,4000,"X"); 
  //Xsection_theory_y3->SetAxisRange(56,2500,"X"); 
  //Xsection_theory_y4->SetAxisRange(56,1588,"X"); 
  //Xsection_theory_y5->SetAxisRange(56,967,"X"); 
  //Xsection_theory_y6->SetAxisRange(56,448,"X"); 
  
  
  //TheorySpectra->Divide(4,2);
  //TheorySpectra->cd(1);
  //Xsection_theory_y0->Draw();
  //spline3_y0->Draw("same");
  //
  //TheorySpectra->cd(2);gPad->SetLogy();
  //Xsection_theory_y1->Draw();
  //spline3_y1->Draw("same");
  //
  //TheorySpectra->cd(3);gPad->SetLogy();
  //Xsection_theory_y2->Draw();
  //spline3_y2->Draw("same");  
  //
  //TheorySpectra->cd(4);gPad->SetLogy();
  //Xsection_theory_y3->Draw();
  //spline3_y3->Draw("same");  
  //
  //TheorySpectra->cd(5);gPad->SetLogy();
  //Xsection_theory_y4->Draw();
  //spline3_y4->Draw("same");  
  //
  //TheorySpectra->cd(6);gPad->SetLogy();
  //Xsection_theory_y5->Draw();
  //spline3_y5->Draw("same");  
  //
  //TheorySpectra->cd(7);gPad->SetLogy();
  //Xsection_theory_y6->Draw();
  //spline3_y6->Draw("same");  
  
  
  //TH1D *True_incl_jet_y0    = new TH1D("True_incl_jet_y0","True_incl_jet_y0", n_thybins_incl, thyBins_incl); True_incl_jet_y0->Sumw2();
  //TH1D *Smeared_incl_jet_y0 = new TH1D("Smeared_incl_jet_y0","Smeared_incl_jet_y0", n_thybins_incl, thyBins_incl); Smeared_incl_jet_y0->Sumw2();
  //TH1D *NormBinWidth_True_incl_jet_y0    = new TH1D("NormBinWidth_True_incl_jet_y0","NormBinWidth_True_incl_jet_y0", n_thybins_incl, thyBins_incl); NormBinWidth_True_incl_jet_y0->Sumw2();
  //TH1D *NormBinWidth_Smeared_incl_jet_y0 = new TH1D("NormBinWidth_Smeared_incl_jet_y0","NormBinWidth_Smeared_incl_jet_y0", n_thybins_incl, thyBins_incl); NormBinWidth_Smeared_incl_jet_y0->Sumw2();
  //RooUnfoldResponse response_y0(Smeared_incl_jet_y0,True_incl_jet_y0);
  //
  //TH1D *True_incl_jet_y1    = new TH1D("True_incl_jet_y1","True_incl_jet_y1", n_thybins_incl, thyBins_incl); True_incl_jet_y1->Sumw2();
  //TH1D *Smeared_incl_jet_y1 = new TH1D("Smeared_incl_jet_y1","Smeared_incl_jet_y1", n_thybins_incl, thyBins_incl); Smeared_incl_jet_y1->Sumw2();
  //TH1D *NormBinWidth_True_incl_jet_y1    = new TH1D("NormBinWidth_True_incl_jet_y1","NormBinWidth_True_incl_jet_y1", n_thybins_incl, thyBins_incl); NormBinWidth_True_incl_jet_y1->Sumw2();
  //TH1D *NormBinWidth_Smeared_incl_jet_y1 = new TH1D("NormBinWidth_Smeared_incl_jet_y1","NormBinWidth_Smeared_incl_jet_y1", n_thybins_incl, thyBins_incl); NormBinWidth_Smeared_incl_jet_y1->Sumw2();
  //RooUnfoldResponse response_y1(Smeared_incl_jet_y1,True_incl_jet_y1);
  //
  //TH1D *True_incl_jet_y2    = new TH1D("True_incl_jet_y2","True_incl_jet_y2", n_thybins_incl, thyBins_incl); True_incl_jet_y2->Sumw2();
  //TH1D *Smeared_incl_jet_y2 = new TH1D("Smeared_incl_jet_y2","Smeared_incl_jet_y2", n_thybins_incl, thyBins_incl); Smeared_incl_jet_y2->Sumw2();
  //TH1D *NormBinWidth_True_incl_jet_y2    = new TH1D("NormBinWidth_True_incl_jet_y2","NormBinWidth_True_incl_jet_y2", n_thybins_incl, thyBins_incl); NormBinWidth_True_incl_jet_y2->Sumw2();
  //TH1D *NormBinWidth_Smeared_incl_jet_y2 = new TH1D("NormBinWidth_Smeared_incl_jet_y2","NormBinWidth_Smeared_incl_jet_y2", n_thybins_incl, thyBins_incl); NormBinWidth_Smeared_incl_jet_y2->Sumw2();
  //RooUnfoldResponse response_y2(Smeared_incl_jet_y2,True_incl_jet_y2);
  //
  //TH1D *True_incl_jet_y3    = new TH1D("True_incl_jet_y3","True_incl_jet_y3", n_thybins_incl, thyBins_incl); True_incl_jet_y3->Sumw2();
  //TH1D *Smeared_incl_jet_y3 = new TH1D("Smeared_incl_jet_y3","Smeared_incl_jet_y3", n_thybins_incl, thyBins_incl); Smeared_incl_jet_y3->Sumw2();
  //TH1D *NormBinWidth_True_incl_jet_y3    = new TH1D("NormBinWidth_True_incl_jet_y3","NormBinWidth_True_incl_jet_y3", n_thybins_incl, thyBins_incl); NormBinWidth_True_incl_jet_y3->Sumw2();
  //TH1D *NormBinWidth_Smeared_incl_jet_y3 = new TH1D("NormBinWidth_Smeared_incl_jet_y3","NormBinWidth_Smeared_incl_jet_y3", n_thybins_incl, thyBins_incl); NormBinWidth_Smeared_incl_jet_y3->Sumw2();
  //RooUnfoldResponse response_y3(Smeared_incl_jet_y3,True_incl_jet_y3);
  //
  //TH1D *True_incl_jet_y4    = new TH1D("True_incl_jet_y4","True_incl_jet_y4", n_thybins_incl, thyBins_incl); True_incl_jet_y4->Sumw2();
  //TH1D *Smeared_incl_jet_y4 = new TH1D("Smeared_incl_jet_y4","Smeared_incl_jet_y4", n_thybins_incl, thyBins_incl); Smeared_incl_jet_y4->Sumw2();
  //TH1D *NormBinWidth_True_incl_jet_y4    = new TH1D("NormBinWidth_True_incl_jet_y4","NormBinWidth_True_incl_jet_y4", n_thybins_incl, thyBins_incl); NormBinWidth_True_incl_jet_y4->Sumw2();
  //TH1D *NormBinWidth_Smeared_incl_jet_y4 = new TH1D("NormBinWidth_Smeared_incl_jet_y4","NormBinWidth_Smeared_incl_jet_y4", n_thybins_incl, thyBins_incl); NormBinWidth_Smeared_incl_jet_y4->Sumw2();
  //RooUnfoldResponse response_y4(Smeared_incl_jet_y4,True_incl_jet_y4);
  //
  //TH1D *True_incl_jet_y5    = new TH1D("True_incl_jet_y5","True_incl_jet_y5", n_thybins_incl, thyBins_incl); True_incl_jet_y5->Sumw2();
  //TH1D *Smeared_incl_jet_y5 = new TH1D("Smeared_incl_jet_y5","Smeared_incl_jet_y5", n_thybins_incl, thyBins_incl); Smeared_incl_jet_y5->Sumw2();
  //TH1D *NormBinWidth_True_incl_jet_y5    = new TH1D("NormBinWidth_True_incl_jet_y5","NormBinWidth_True_incl_jet_y5", n_thybins_incl, thyBins_incl); NormBinWidth_True_incl_jet_y5->Sumw2();
  //TH1D *NormBinWidth_Smeared_incl_jet_y5 = new TH1D("NormBinWidth_Smeared_incl_jet_y5","NormBinWidth_Smeared_incl_jet_y5", n_thybins_incl, thyBins_incl); NormBinWidth_Smeared_incl_jet_y5->Sumw2();
  //RooUnfoldResponse response_y5(Smeared_incl_jet_y5,True_incl_jet_y5);
  //
  //TH1D *True_incl_jet_y6    = new TH1D("True_incl_jet_y6","True_incl_jet_y6", n_thybins_incl, thyBins_incl); True_incl_jet_y6->Sumw2();
  //TH1D *Smeared_incl_jet_y6 = new TH1D("Smeared_incl_jet_y6","Smeared_incl_jet_y6", n_thybins_incl, thyBins_incl); Smeared_incl_jet_y6->Sumw2();
  //TH1D *NormBinWidth_True_incl_jet_y6    = new TH1D("NormBinWidth_True_incl_jet_y6","NormBinWidth_True_incl_jet_y6", n_thybins_incl, thyBins_incl); NormBinWidth_True_incl_jet_y6->Sumw2();
  //TH1D *NormBinWidth_Smeared_incl_jet_y6 = new TH1D("NormBinWidth_Smeared_incl_jet_y6","NormBinWidth_Smeared_incl_jet_y6", n_thybins_incl, thyBins_incl); NormBinWidth_Smeared_incl_jet_y6->Sumw2();
  //RooUnfoldResponse response_y6(Smeared_incl_jet_y6,True_incl_jet_y6);
  
  









  /*
    printf("Evaluating y0\n");
    
    //// y0
    for(int i=0;i<nEvents;++i){  
    double ptmin = 32.;
    double ptmax = 6000.;
    
    double ptTrue  = rnd->Uniform(ptmin,ptmax);
    double sigma   = gres_y0(ptTrue);
    
    //sigma=sigma*1.079; //JER Scaling factor 8 TeV  Acoording to Mikko no scaling for the moment
    
    double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
    double pt_w      =  spline3_y0->Eval(ptTrue);

    True_incl_jet_y0    -> Fill(ptTrue,pt_w);
    Smeared_incl_jet_y0 -> Fill(ptSmeared,pt_w);
    
    //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
    
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
       response_y0.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
     response_y0.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
     response_y0.Fake(ptSmeared,pt_w);
     }
     
     
  }
  
  NormBinWidth_True_incl_jet_y0    = cross_section(True_incl_jet_y0);
  NormBinWidth_Smeared_incl_jet_y0 = cross_section(Smeared_incl_jet_y0);
  NormBinWidth_True_incl_jet_y0->SetLineColor(2);
  NormBinWidth_Smeared_incl_jet_y0->SetLineColor(2);  
  
  printf("Evaluating y1\n");
  
  //// y1
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32.;
     double ptmax = 5200.;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y1(ptTrue);

     //sigma=sigma*1.099; //JER Scaling factor 8 TeV

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y1->Eval(ptTrue);

     True_incl_jet_y1    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y1 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
     
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
       response_y1.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
       response_y1.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
       response_y1.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y1    = cross_section(True_incl_jet_y1);
   NormBinWidth_Smeared_incl_jet_y1 = cross_section(Smeared_incl_jet_y1);
   NormBinWidth_True_incl_jet_y1->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y1->SetLineColor(2);


  printf("Evaluating y2\n");

  //// y2
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 4000;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y2(ptTrue);

     //sigma=sigma*1.121; //JER Scaling factor 8 TeV

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y2->Eval(ptTrue);
     
     True_incl_jet_y2    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y2 -> Fill(ptSmeared,pt_w);
     
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
     
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y2.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y2.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y2.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y2    = cross_section(True_incl_jet_y2);
   NormBinWidth_Smeared_incl_jet_y2 = cross_section(Smeared_incl_jet_y2);
   NormBinWidth_True_incl_jet_y2->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y2->SetLineColor(2);



  printf("Evaluating y3\n");

  //// y3
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 2500;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y3(ptTrue);

     //sigma=sigma*1.208; //JER Scaling factor 8 TeV

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y3->Eval(ptTrue);

     True_incl_jet_y3    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y3 -> Fill(ptSmeared,pt_w);
     
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
 
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y3.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y3.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y3.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y3    = cross_section(True_incl_jet_y3);
   NormBinWidth_Smeared_incl_jet_y3 = cross_section(Smeared_incl_jet_y3);
   NormBinWidth_True_incl_jet_y3->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y3->SetLineColor(2);


  printf("Evaluating y4\n");


  //// y4
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 1588;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y4(ptTrue);     
          
     //sigma=sigma*1.208; //JER Scaling factor 8 TeV

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y4->Eval(ptTrue);

     True_incl_jet_y4    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y4 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
 
     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y4.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y4.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y4.Fake(ptSmeared,pt_w);
     }


  }

   NormBinWidth_True_incl_jet_y4    = cross_section(True_incl_jet_y4);
   NormBinWidth_Smeared_incl_jet_y4 = cross_section(Smeared_incl_jet_y4);
   NormBinWidth_True_incl_jet_y4->SetLineColor(2);
   NormBinWidth_Smeared_incl_jet_y4->SetLineColor(2);


  printf("Evaluating y5\n");

  //// y5
  for(int i=0;i<nEvents;++i){  
     double ptmin = 32;
     double ptmax = 967;
    
     double ptTrue  = rnd->Uniform(ptmin,ptmax);
     double sigma   = gres_y5(ptTrue);     
  
     //sigma=sigma*1.254; //JER Scaling factor 8 TeV

     double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
     double pt_w      =  spline3_y5->Eval(ptTrue);

     True_incl_jet_y5    -> Fill(ptTrue,pt_w);
     Smeared_incl_jet_y5 -> Fill(ptSmeared,pt_w);
  
     //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);

     if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
          response_y5.Fill(ptSmeared,ptTrue,pt_w);
     }
     else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
          response_y5.Miss(ptTrue,pt_w);
     }
     else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
          response_y5.Fake(ptSmeared,pt_w);
     }


  }

  NormBinWidth_True_incl_jet_y5    = cross_section(True_incl_jet_y5);
  NormBinWidth_Smeared_incl_jet_y5 = cross_section(Smeared_incl_jet_y5);
  NormBinWidth_True_incl_jet_y5->SetLineColor(2);
  NormBinWidth_Smeared_incl_jet_y5->SetLineColor(2);
  
  
  printf("Evaluating y6\n");
  
  //// y6
  for(int i=0;i<nEvents;++i){  
    double ptmin = 32;
    double ptmax = 468;   //468;
    
    double ptTrue  = rnd->Uniform(ptmin,ptmax);
    double sigma   = gres_y6(ptTrue);     
    
    //sigma=sigma*1.056; //JER Scaling factor 8 TeV
    
    double ptSmeared =  rnd->Gaus(ptTrue,ptTrue*sigma);
    double pt_w      =  spline3_y6->Eval(ptTrue);
    
    True_incl_jet_y6    -> Fill(ptTrue,pt_w);
    Smeared_incl_jet_y6 -> Fill(ptSmeared,pt_w);
    
    //printf("ptTrue=%f  ptSmeared=%f sigma=%f  ptTrue/ptSmeared=%f  pt_w=%f\n", ptTrue, ptSmeared, sigma, ptTrue/ptSmeared ,pt_w);
    
    if (((ptSmeared>=ptmin)&&(ptSmeared<=ptmax))&&((ptTrue>=ptmin)&&(ptTrue<=ptmax))){
      response_y6.Fill(ptSmeared,ptTrue,pt_w);
    }
    else if ((ptTrue>=ptmin)&&(ptTrue<=ptmax)){
      response_y6.Miss(ptTrue,pt_w);
    }
    else if ((ptSmeared>=ptmin)&&(ptSmeared<=ptmax)){
      response_y6.Fake(ptSmeared,pt_w);
    }        
  }
  
  NormBinWidth_True_incl_jet_y6    = cross_section(True_incl_jet_y6);
  NormBinWidth_Smeared_incl_jet_y6 = cross_section(Smeared_incl_jet_y6);
  NormBinWidth_True_incl_jet_y6->SetLineColor(2);
  NormBinWidth_Smeared_incl_jet_y6->SetLineColor(2);
  */








//  TH1F *rat_y0_true_smeared=(TH1F*)True_incl_jet_y0->Clone(); rat_y0_true_smeared->SetTitle("y0 True/Smeared;p_T (GeV); y0 True/Smeared          ");
//  rat_y0_true_smeared->Divide(True_incl_jet_y0,Smeared_incl_jet_y0,1.,1.,"B");
//  
//  TH1F *rat_y1_true_smeared=(TH1F*)True_incl_jet_y1->Clone(); rat_y1_true_smeared->SetTitle("y1 True/Smeared;p_T (GeV); y1 True/Smeared          ");
//  rat_y1_true_smeared->Divide(True_incl_jet_y1,Smeared_incl_jet_y1,1.,1.,"B");
//  
//  TH1F *rat_y2_true_smeared=(TH1F*)True_incl_jet_y2->Clone(); rat_y2_true_smeared->SetTitle("y2 True/Smeared;p_T (GeV); y2 True/Smeared          ");
//  rat_y2_true_smeared->Divide(True_incl_jet_y2,Smeared_incl_jet_y2,1.,1.,"B");
//  
//  TH1F *rat_y3_true_smeared=(TH1F*)True_incl_jet_y3->Clone(); rat_y3_true_smeared->SetTitle("y3 True/Smeared;p_T (GeV); y3 True/Smeared          ");
//  rat_y3_true_smeared->Divide(True_incl_jet_y3,Smeared_incl_jet_y3,1.,1.,"B");
//  
//  TH1F *rat_y4_true_smeared=(TH1F*)True_incl_jet_y4->Clone(); rat_y4_true_smeared->SetTitle("y4 True/Smeared;p_T (GeV); y4 True/Smeared          ");
//  rat_y4_true_smeared->Divide(True_incl_jet_y4,Smeared_incl_jet_y4,1.,1.,"B");
//  
//  TH1F *rat_y5_true_smeared=(TH1F*)True_incl_jet_y5->Clone(); rat_y5_true_smeared->SetTitle("y5 True/Smeared;p_T (GeV); y5 True/Smeared          ");
//  rat_y5_true_smeared->Divide(True_incl_jet_y5,Smeared_incl_jet_y5,1.,1.,"B");
//  
//  TH1F *rat_y6_true_smeared=(TH1F*)True_incl_jet_y6->Clone(); rat_y6_true_smeared->SetTitle("y6 True/Smeared;p_T (GeV); y6 True/Smeared          ");
//  rat_y6_true_smeared->Divide(True_incl_jet_y6,Smeared_incl_jet_y6,1.,1.,"B");
    
  
  //  rat_y0_true_smeared->SetMaximum(1.1); rat_y0_true_smeared->SetAxisRange(97,6000,"X"); 
  //  rat_y1_true_smeared->SetMaximum(1.1); rat_y1_true_smeared->SetAxisRange(97,5200,"X"); 
  //  rat_y2_true_smeared->SetMaximum(1.1); rat_y2_true_smeared->SetAxisRange(97,4000,"X"); 
  //  rat_y3_true_smeared->SetMaximum(1.1); rat_y3_true_smeared->SetAxisRange(97,2500,"X"); 
  //  rat_y4_true_smeared->SetMaximum(1.1); rat_y4_true_smeared->SetAxisRange(97,1588,"X"); 
  //  rat_y5_true_smeared->SetMaximum(1.1); rat_y5_true_smeared->SetAxisRange(97,967,"X"); 
  //  rat_y6_true_smeared->SetMaximum(1.1); rat_y6_true_smeared->SetAxisRange(97,448,"X"); 
  
  



  //trueOverSmeared->Divide(4,2);
  //trueOverSmeared->cd(1);
  //gPad->SetGrid(); 
  //rat_y0_true_smeared->Draw();     
  //trueOverSmeared->cd(2);gPad->SetGrid(); 
  //rat_y1_true_smeared->Draw();     
  //trueOverSmeared->cd(3);gPad->SetGrid(); 
  //rat_y2_true_smeared->Draw();     
  //trueOverSmeared->cd(4);gPad->SetGrid(); 
  //rat_y3_true_smeared->Draw();     
  //trueOverSmeared->cd(5);gPad->SetGrid(); 
  //rat_y4_true_smeared->Draw();     
  //trueOverSmeared->cd(6); gPad->SetGrid();
  //rat_y5_true_smeared->Draw();     
  //trueOverSmeared->cd(7);gPad->SetGrid(); 
  //rat_y6_true_smeared->Draw();     

  //Xsection_theory_y0->Write();
  //True_incl_jet_y0->Write();
  //NormBinWidth_True_incl_jet_y0->Write("NormBinWidth_True_incl_jet_y0");    
  //Smeared_incl_jet_y0->Write(); 
  //NormBinWidth_Smeared_incl_jet_y0->Write("NormBinWidth_Smeared_incl_jet_y0"); 
  //response_y0.Write("response_y0");
  //
  //Xsection_theory_y1->Write();
  //True_incl_jet_y1->Write();
  //NormBinWidth_True_incl_jet_y1->Write("NormBinWidth_True_incl_jet_y1");    
  //Smeared_incl_jet_y1->Write(); 
  //NormBinWidth_Smeared_incl_jet_y1->Write("NormBinWidth_Smeared_incl_jet_y1"); 
  //response_y1.Write("response_y1");
  //
  //Xsection_theory_y2->Write();
  //True_incl_jet_y2->Write();
  //NormBinWidth_True_incl_jet_y2->Write("NormBinWidth_True_incl_jet_y2");    
  //Smeared_incl_jet_y2->Write(); 
  //NormBinWidth_Smeared_incl_jet_y2->Write("NormBinWidth_Smeared_incl_jet_y2"); 
  //response_y2.Write("response_y2");
  //
  //Xsection_theory_y3->Write();
  //True_incl_jet_y3->Write();
  //NormBinWidth_True_incl_jet_y3->Write("NormBinWidth_True_incl_jet_y3");    
  //Smeared_incl_jet_y3->Write(); 
  //NormBinWidth_Smeared_incl_jet_y3->Write("NormBinWidth_Smeared_incl_jet_y3"); 
  //response_y3.Write("response_y3");
  //
  //Xsection_theory_y4->Write();
  //True_incl_jet_y4->Write();
  //NormBinWidth_True_incl_jet_y4->Write("NormBinWidth_True_incl_jet_y4");    
  //Smeared_incl_jet_y4->Write(); 
  //NormBinWidth_Smeared_incl_jet_y4->Write("NormBinWidth_Smeared_incl_jet_y4"); 
  //response_y4.Write("response_y4");
  //
  //Xsection_theory_y5->Write();
  //True_incl_jet_y5->Write();
  //NormBinWidth_True_incl_jet_y5->Write("NormBinWidth_True_incl_jet_y5");    
  //Smeared_incl_jet_y5->Write(); 
  //NormBinWidth_Smeared_incl_jet_y5->Write("NormBinWidth_Smeared_incl_jet_y5"); 
  //response_y5.Write("response_y5");
  //
  //Xsection_theory_y6->Write();
  //True_incl_jet_y6->Write();
  //NormBinWidth_True_incl_jet_y6->Write("NormBinWidth_True_incl_jet_y6");    
  //Smeared_incl_jet_y6->Write(); 
  //NormBinWidth_Smeared_incl_jet_y6->Write("NormBinWidth_Smeared_incl_jet_y6"); 
  //response_y6.Write("response_y6");
  
  
  
  
