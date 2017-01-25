// unfoldSpectra.C and unfoldSpectra.h written by Ian W. Laflotte
// June 28 2016
// heavily based on unfold.C and utilities.h bu Raghav Kunnawalkam Elayavalli
// Macro to perform Bayes and SVD Unfolding direct from the official RooUnfold classes
// compile with...
// g++ -L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold -I"${ROOUNFOLD_PATH}/src" unfoldSpectra.C $(root-config --cflags --libs) -Werror -Wall -O2 -o unfoldSpectra.exe

// C++, C, etc. 
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <vector>

// ROOTSYS
#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>

// I/O
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TNtuple.h>

// utilities
#include <TEventList.h>
#include <TCut.h>
#include <TMath.h>
#include <TError.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSpline.h>
#include <TKey.h>
#include <TMatrix.h>
#include <TMatrixT.h>
#include <TMatrixD.h>
#include <TSVDUnfold.h>

// plotting
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

// histos
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>

// RooUnfold
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"


//-----------------------------------------------------------------------------------------------------------------------
const double boundaries_pt[] = {
3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 
174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 1000};
const int nbins_pt = sizeof(boundaries_pt)/sizeof(double)-1;

const double ptbins[] = {
0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 
200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 
380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 
560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 
740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 
920, 930, 940, 950, 960, 970, 980, 990, 1000};
const int nbins = sizeof(ptbins)/sizeof(double)-1;

//-----------------------------------------------------------------------------------------------------------------------
// color scheme for different color
const Int_t color[13]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,kRed+2,kMagenta+1,kRed+1,kGreen+3};
// color scheme for fill area - fcolor with lcolor as a boundary
const Int_t fcolor[5]={kRed-10,kBlue-10,kGreen-8,kOrange-4, kRed-6};
const Int_t lcolor[5]={kRed+1,kBlue-3,kGreen+3,kOrange+3, kRed+2};
// marker scheme 
const Int_t fmstyle[6] = {20,21,22,23,29,3};
const Int_t emstyle[6] = {24,25,26,27,30,28};


//-----------------------------------------------------------------------------------------------------------------------
void SetUnfoldBins1D(TH1F *&h, float minpt, float maxpt){
  int lbin = h->GetXaxis()->FindBin(minpt);
  int hbin = h->GetXaxis()->FindBin(maxpt);
  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    if(ix > lbin && ix < hbin)continue;
    h->SetBinContent(ix,0);
    h->SetBinError(ix,0);
  }
  //return;
}


void SetUnfoldBins2D(TH2F *&h, float minGenPt, float maxGenPt, float minRecPt, float maxRecPt){
  int lgen_bin = h->GetYaxis()->FindBin(minGenPt);
  int hgen_bin = h->GetYaxis()->FindBin(maxGenPt);
  int lrec_bin = h->GetXaxis()->FindBin(minRecPt);
  int hrec_bin = h->GetXaxis()->FindBin(maxRecPt);
  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    for(int iy=1; iy<=h->GetNbinsY(); iy++){
      if( (ix>lrec_bin && ix<hrec_bin) && 
	  (iy>lgen_bin && iy<hgen_bin) ) continue;
      h->SetBinContent(ix, iy, 0);
      h->SetBinError(ix, iy, 0);
    }
  }
  //return;
}


void Truncate1D(TH1 * hInput, TH1 * hSetup){
  for(int k = 1; k<=hSetup->GetNbinsX(); ++k){
    int binx = hInput->FindBin(hSetup->GetBinCenter(k));
    Double_t val = hInput->GetBinContent(binx);
    Double_t val_Err = hInput->GetBinError(binx);
    hSetup->SetBinContent(k, val);
    hSetup->SetBinError(k, val_Err);
  }
  //hInput->Reset();
  //hInput = (TH1F*)hSetup->Clone(Form("%s",newname));
}


void doCorr(TH1F *h1, TH1F *h2, int mode = 0){
  for(int ix = 1; ix<=h2->GetNbinsX(); ++ix){
    int xbin = h1->FindBin(h2->GetBinCenter(ix));
    double val = h1->GetBinContent(xbin);
    double err = h1->GetBinError(xbin);
    switch (mode){
    case 0:
      val = val * h2->GetBinContent(ix);
      err = err * h2->GetBinContent(ix);
      break;
    case 1:
      val = val * (double)1./h2->GetBinContent(ix);
      err = err * (double)1./h2->GetBinContent(ix);
      break;
    case 2:
      val = val * 1./(1 - h2->GetBinContent(ix));
      err = err * 1./(1 - h2->GetBinContent(ix));
    }
    h1->SetBinContent(xbin, val);    
    h1->SetBinError(xbin, err);    
  }
  //return;
}


void doCorr2D(TH2F *h1, TH1F * h2, int mode = 0){
  for(int ix = 1; ix<=h2->GetNbinsX(); ++ix){
    int xbin = h1->GetXaxis()->FindBin(h2->GetBinCenter(ix));
    for(int iy = 1; iy<=h1->GetNbinsY(); iy++){
      double val = h1->GetBinContent(xbin, iy);
      //std::cout<<"bin = "<<h1->GetYaxis()->GetBinCenter(xbin)<<"  initial value = "<<val<<std::endl;
      switch (mode){
      case 0:	val = val * (1.0 + h2->GetBinContent(ix));	break;
      case 1:	val = val * (double)1./h2->GetBinContent(ix);	break;
      case 2:	val = val * 1./(1 - h2->GetBinContent(ix));
      }
      h1->SetBinContent(xbin, iy, val);
      //std::cout<<"final value = "<<val<<std::endl;
    }
  }
  //return;
}


void Truncate2D(TH2 * hInput, TH2 * hSetup){
  for(int k = 1; k<=hSetup->GetNbinsX(); ++k){
    for(int l = 1; l<=hSetup->GetNbinsY(); ++l){      
      int binx = hInput->GetXaxis()->FindBin(hSetup->GetXaxis()->GetBinCenter(k));
      int biny = hInput->GetYaxis()->FindBin(hSetup->GetYaxis()->GetBinCenter(l));
      Double_t val = hInput->GetBinContent(binx, biny);
      Double_t val_Err = hInput->GetBinError(binx, biny);
      hSetup->SetBinContent(k, l, val);
      hSetup->SetBinError(k, l, val_Err);
    }
  }
  //return;
} // end Truncate2D
//  //hInput->Reset();
//  //hInput = (TH2F*)hSetup->Clone(Form("%s",newname));
//} // end Truncate2D


void doTrigCorr(TH1F *h1, TH1F *h2){
  for(int ix = 1; ix<=h2->GetNbinsX(); ++ix){
    int xbin = h1->FindBin(h2->GetBinCenter(ix));
    double val = h1->GetBinContent(xbin);
    double err = h1->GetBinError(xbin);
    val = val * 1./h2->GetBinContent(ix);
    h1->SetBinContent(xbin, val);    
    err = err * 1./h2->GetBinContent(ix);
    h1->SetBinError(xbin, err);
  }
  //return;
}


TH1D* GetDummyHist1(Int_t bins,Float_t xmin, Float_t xmax, Float_t min, Float_t max,const char *ttl,const char *xttl,const char *yttl, Bool_t stat){
   TH1D *dum;
   dum = new TH1D(ttl,"",bins,xmin,xmax);
   dum->SetMinimum(min);
   dum->SetMaximum(max);
   dum->SetStats(stat);
   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();
   return dum;
}


TH1D* GetDummyHist(Float_t xmin, Float_t xmax, Float_t min, Float_t max,const char *ttl,const char *xttl,const char *yttl, Bool_t stat) {
   TH1D *dum;
   dum = new TH1D(ttl,"",100,xmin,xmax);
   dum->SetMinimum(min);
   dum->SetMaximum(max);
   dum->SetStats(stat);
   dum->GetYaxis()->SetTitle(yttl);
   dum->GetYaxis()->CenterTitle();
   dum->GetXaxis()->SetTitle(xttl);
   dum->GetXaxis()->CenterTitle();
   return dum;
}


// convert hist to tgraph
TGraphErrors* TgraphItWithRange(TH1D* hist, double min, double max){

  int nbins = 0;
  for(int i=0;i<hist->GetNbinsX();i++){
    double xvalue = hist->GetBinCenter(i+1);
    if(xvalue>min && xvalue<max) nbins++;
  }
  const int nlines = nbins;

  TGraphErrors *tg;
  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];
  int counter = 0;
  for(int i = 0; i<hist->GetNbinsX(); i++ ){
    double xvalue = hist->GetBinCenter(i+1);
    if(xvalue>min && xvalue<max){
      pt[counter] = hist->GetBinCenter(i+1);
      xsec[counter] = hist->GetBinContent(i+1);
      xsecerr[counter] = hist->GetBinError(i+1);
      pterr[counter] = 0;
      counter++;
    }
  }
  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}


TGraphErrors* TgraphIt(TH1D* hist){

  int nbins = hist->GetNbinsX();
  const int nlines = nbins;

  TGraphErrors *tg;
  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];
  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    pterr[i] = 0;
  }
  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}


TGraphErrors* TgraphIt(TH1F* hist){

  int nbins = hist->GetNbinsX();
  const int nlines = nbins;

  TGraphErrors *tg;
  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];
  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    pterr[i] = 0;
  }
  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}


TH1D* HistIt(TGraph* tg){

  int nbins = tg->GetN();
  const int nlines = nbins;
  double pt[nlines], xsec[nlines], xsecerr[nlines];//, pterr[nlines];

  for(int i = 0; i<nbins; i++ ){
    tg->GetPoint(i,pt[i],xsec[i]);
    //pterr[i]   = tg->GetErrorX(i);
    xsecerr[i] = tg->GetErrorY(i);
  }

  //init. of array here makes me nervous
  double bins[nlines];
  bins[0] = 0.0;
  for(int i=0;i<nlines;i++){
    double dbin = pt[i] - bins[i];
    bins[i+1]   = pt[i] + dbin;
  }

  TH1D *hist = new TH1D("","",nbins,bins);
  for(int i = 0; i<nbins; i++ ){
    hist->SetBinContent(i+1,xsec[i]);
    hist->SetBinError  (i+1,xsecerr[i]);
  }  
  return hist;
}


TH1D* ImprovedHistIt(TGraph* tg){

  int nbins = tg->GetN();
  const int nlines = nbins;

  double pt[nlines], xsec[nlines], xsecerr[nlines];//,double pterr[nlines];

  for(int i = 0; i<nbins; i++ ){
    tg->GetPoint(i,pt[i],xsec[i]);
    //pterr[i]   = tg->GetErrorX(i);
    xsecerr[i] = tg->GetErrorY(i);
  }

  double bins[nlines];
  if(pt[0]!=0) bins[0] = pt[0] - 0.5*(pt[1]-pt[0]);
  else bins[0] = 0.0;

  int ncount=0;
  for(int i=0;i<nlines;i++){
    double dbin = pt[i] - bins[i];
    bins[i+1] = pt[i] + dbin;
    std::cout<<"pt = "<<pt[i]<<" dbin = "<<dbin<<" bin boundary = "<<bins[i]<<std::endl;
    ncount++;
  }
  bins[ncount+1] = bins[ncount] + (bins[ncount]-bins[ncount-1]); // final bin

  TH1D *hist = new TH1D("","",nbins,bins);
  for(int i = 0; i<nbins; i++ ){
    hist->SetBinContent(i+1,xsec[i]);
    hist->SetBinError  (i+1,xsecerr[i]);
  }
  return hist;
}


TH1D* DivideByNearBinRAA(TH1D* num, TH1D* den){

  int nNum = num->GetXaxis()->GetNbins();
  TH1D* ratio = (TH1D*) num->Clone("ratio");
  ratio->Reset();

  for(int i=0;i<nNum;i++){
    double xvalue = num->GetXaxis()->GetBinCenter(i+1);
    double denvalue = den->GetBinContent(den->GetXaxis()->FindBin(xvalue));
    if (denvalue==0) return 0;
    double rat = num->GetBinContent(i+1)/denvalue;
    ratio->SetBinContent(i+1,rat);
    ratio->SetBinError(i+1,0.0001*rat); // need to be calculated correctly                                                                                  
  }
  return ratio;
}


TH1D* ratio_hist_to_tspline(TH1D* num, TSpline3* ts, double minx, double maxx){

  std::cout<<"[Ratio to fit with range used]"<<std::endl;

  TH1D *hRatio = (TH1D*) num->Clone("hRatio");
  int nbin = num->GetNbinsX();

  for(int i=0;i<nbin;i++){
    double cms_value = (double) ts->Eval(hRatio->GetBinCenter(i+1));
    double ratio = -999;
    double ratio_err = 0;
    if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
      ratio = hRatio->GetBinContent(i+1)/cms_value;
      ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
    }
    else{
      ratio = -999;
      ratio_err = 0.0;
    }
    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,ratio_err);
  }
  return hRatio;
}


TGraphErrors* ratio_func_to_tspline(TGraphErrors* num, TSpline3* ts, double minx, double maxx){

  std::cout<<"[Ratio to fit used]"<<std::endl;

  TGraphErrors *tg;
  int nbin = num->GetN();
  const int nlines = nbin;
  double pt[nlines], pterr[nlines];
  double xsec[nlines], xsecerr[nlines];
  double ratio[nlines], ratioerr[nlines];
  for(int i=0;i<nbin;i++){
    num->GetPoint(i,pt[i],xsec[i]);
    xsecerr[i] = num->GetErrorY(i);
    double cms_value = (double) ts->Eval(pt[i]);
    if(pt[i]>minx && pt[i]<maxx){
      ratio[i] = xsec[i]/cms_value;
      ratioerr[i] = (xsecerr[i]/xsec[i])*ratio[i];
    }else{
      ratio[i] = -999;
      ratioerr[i] = 0.0;
    }
    pterr[i] = 0;
  }
  tg = new TGraphErrors(nlines,pt,ratio,pterr,ratioerr);
  return tg;
}


TGraphErrors* ratio_tg_to_hist(TGraphErrors* num, TH1D *hist, double minx, double maxx){

  std::cout<<"[Ratio to fit used]"<<std::endl;

  TGraphErrors *tg;

  int nbin = num->GetN();
  const int nlines = nbin;
  double pt[nlines], pterr[nlines];
  double xsec[nlines], xsecerr[nlines];
  double ratio[nlines], ratioerr[nlines];

  for(int i=0;i<nbin;i++){
    num->GetPoint(i,pt[i],xsec[i]);
    xsecerr[i] = num->GetErrorY(i);
    if(pt[i]>minx && pt[i]<maxx){
      double ptx= pt[i]+0.0001;
      int ptbin = hist->FindBin(ptx);
      double cms_value = (double) hist->GetBinContent(ptbin);
      ratio[i] = xsec[i]/cms_value;
      ratioerr[i] = (xsecerr[i]/xsec[i])*ratio[i];
      std::cout<<"pt = "<<pt[i]<<" pt bin = "<<hist->GetBinCenter(ptbin)<<std::endl;
      //std::cout<<"cms value = "<<cms_value<<" xsec[i] = "<<xsec[i]<<" and ratio = "<<ratio[i]<<std::endl;
    }
    else{
      ratio[i] = -999;
      ratioerr[i] = 0.0;
    }
    //std::cout<<"pt = "<<pt[i]<<" ratio = "<<ratio[i]<<std::endl;
    pterr[i] = 0;
  }
  tg = new TGraphErrors(nlines,pt,ratio,pterr,ratioerr);
  return tg;
}


TH1D* ratio_hist_to_func(TH1D* num, TF1* f3, double minx, double maxx){

  std::cout<<"[Ratio to fit used]"<<std::endl;

  TH1D *hRatio = (TH1D*) num->Clone("hRatio");

  int nbin = hRatio->GetNbinsX();
  for(int i=0;i<nbin;i++){
    double cms_value = (double) f3->Eval(hRatio->GetBinCenter(i+1));
    double ratio = -999;
    double ratio_err = 0;
    if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
      ratio = hRatio->GetBinContent(i+1)/cms_value;
      ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
    }
    else{
      ratio = -999;
      ratio_err = 0.0;
    }
    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,ratio_err);
  }
  return hRatio;
}


TH1F* ratio_hist_to_func(TH1F* num, TF1* f3, double minx, double maxx){

  std::cout<<"[Ratio to fit used]"<<std::endl;

  TH1F *hRatio = (TH1F*) num->Clone("hRatio");

  int nbin = hRatio->GetNbinsX();
  for(int i=0;i<nbin;i++){
    double cms_value = (double) f3->Eval(hRatio->GetBinCenter(i+1));
    double ratio = -999;
    double ratio_err = 0;
    if(hRatio->GetBinCenter(i+1)>minx && hRatio->GetBinCenter(i+1)<maxx){
      ratio = hRatio->GetBinContent(i+1)/cms_value;
      ratio_err =  (hRatio->GetBinError(i+1)/hRatio->GetBinContent(i+1))*ratio;
    }
    else{
      ratio = -999;
      ratio_err = 0.0;
    }
    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,ratio_err);
  }
  return hRatio;
}

TGraphErrors* ratio_tg_to_tg(TGraphErrors* num1, TGraphErrors *num2){

  std::cout<<"[Ratio to tg used]"<<std::endl;

  TGraphErrors *tg;

  int nbin = num1->GetN();
  const int nlines = nbin;
  double x1[nlines], x1err[nlines];
  double y1[nlines], y1err[nlines];
  double y2[nlines], y2err[nlines];
  double ratio[nlines], ratioerr[nlines];
  for(int i=0;i<nbin;i++){
    num1->GetPoint(i,x1[i],y1[i]);
    y1err[i] = num1->GetErrorY(i);
    num2->GetPoint(i,x1[i],y2[i]);
    y2err[i] = num2->GetErrorY(i);
    if(y2[i]!=0){
      ratio[i]    = y1[i]/y2[i];
      ratioerr[i] = ratio[i]*sqrt(pow(y1err[i]/y1[i],2) + pow(y2err[i]/y2[i],2));
    }
    else{
      ratio[i] = -999;
      ratioerr[i] = 0.0;
    }
  }
  tg = new TGraphErrors(nlines,x1,ratio,x1err,ratioerr);
  return tg;
}


void makeMultiPanelCanvas(TCanvas*& canv,  const Int_t columns,   const Int_t rows,
                          const Float_t leftOffset, const Float_t bottomOffset, const Float_t leftMargin, const Float_t bottomMargin, 
			  const Float_t edge,const Float_t asyoffset) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];
  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];

  Float_t PadWidth = (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
				       (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight = (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
					  (1.0/(1.0-edge))+(Float_t)rows-2.0);
  
   Xlow[0] = leftOffset - asyoffset;
   Xup[0]  = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
     Xlow[i] = Xup[0] + (i-1)*PadWidth;
     Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
     Ylow[i] = Yup[rows-1] + ct*PadHeight;
     Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
     ct++;
   }

   TString padName;
   for(Int_t i=0; i<columns; i++) {
     for(Int_t j=0; j<rows; j++) {
       canv->cd();
       padName = Form("p_%d_%d",i,j);
       pad[i][j] = new TPad(padName.Data(),padName.Data(),
                            Xlow[i],Ylow[j],Xup[i],Yup[j]);
       // this is hacked version to create aysmmetric pads around low 
       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0);
       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0);
       if(j==0) pad[i][j]->SetTopMargin(edge);
       else pad[i][j]->SetTopMargin(0);
       //if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
       //else pad[i][j]->SetTopMargin(0);
       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0);
       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
     }
   }
   //return;
} // end makeMultiPanelCanvas


void makeMultiPanelCanvasAsy(TCanvas*& canv,  const Int_t columns,  const Int_t rows,
			     const Float_t leftOffset,  const Float_t bottomOffset,  const Float_t leftMargin,  const Float_t bottomMargin,
			     const Float_t edge, const Float_t asyoffset ){
  if( canv==0 ){
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
   Float_t PadWidth =
     (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                       (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
     (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                         (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
     Xlow[i] = Xup[0] + (i-1)*PadWidth;
     Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;

   for(Int_t i=rows-2;i>0;i--) {
     Ylow[i] = Yup[rows-1] + ct*PadHeight;
     Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
     ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
     for(Int_t j=0;j<rows;j++) {
       canv->cd();
       padName = Form("p_%d_%d",i,j);

       // this is hacked version to create aysmmetric pads around low 
       if(j==0){
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
       }else{
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
       }

       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0);

       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0);

       if(j==0) pad[i][j]->SetTopMargin(edge);
       else pad[i][j]->SetTopMargin(0);

       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0);

       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
     }
   }
} // end makeMultiPanelCanvasAsy


void makeMultiPanelCanvasWithGap(TCanvas* &canv,  const Int_t columns,  const Int_t rows,
				 const Float_t leftOffset,  const Float_t bottomOffset,  const Float_t leftMargin,  const Float_t bottomMargin,
				 const Float_t edge, const Float_t asyoffset) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];

  Float_t PadWidth =
     (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                       (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
     (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                         (1.0/(1.0-edge))+(Float_t)rows-2.0);

   //PadHeight = 0.5*PadHeight;

   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
     Xlow[i] = Xup[0] + (i-1)*PadWidth;
     Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
     if(i==rows-2){
       Ylow[i] = Yup[rows-1] + ct*PadHeight;
       Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
     }else{
       Ylow[i] = Yup[rows-1] + ct*PadHeight;
       Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
       //Yup[i] = 0.2*Yup[i];
     }
     ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
     for(Int_t j=0;j<rows;j++) {
       canv->cd();
       padName = Form("p_%d_%d",i,j);
       //pad[i][j] = new TPad(padName.Data(),padName.Data(),
       //Xlow[i],Ylow[j],Xup[i],Yup[j]);
       // this is hacked version to create aysmmetric pads around low 
       if(j==0){
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
       }else{
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
                              Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
       }


       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0);

       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0);

       if(j==0) pad[i][j]->SetTopMargin(edge);
       //else pad[i][j]->SetTopMargin(0.01);
       else pad[i][j]->SetTopMargin(0.02);

       //if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
       //else pad[i][j]->SetTopMargin(0.0);

       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0.15);

       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
     }
   }
} // end makeMultiPanelCanvasWithGap


void makeMultiPanelCanvasXMerged(TCanvas*& canv, const Int_t columns, const Int_t rows,
				 const Float_t leftOffset, const Float_t rightOffset, const Float_t bottomOffset, 
				 const Float_t leftMargin, const Float_t rightMargin, const Float_t topMargin, const Float_t bottomMargin,
				 const Float_t edge, const Float_t asyoffset){
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];
  
  Float_t Xlow[columns];
  Float_t Xup [columns];
  Float_t Ylow[rows];
  Float_t Yup [rows];
  Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
		      (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			(1.0/(1.0-edge))+(Float_t)rows-2.0);

  std::cout<<"PadWidth : "<<PadWidth<<"\t PadHeight : "<<PadHeight<<std::endl;
  Xlow[0] = leftOffset;
  Xup [0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
  
  Yup [0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup [rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
  
  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i]  = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i]  = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);


      // this is hacked version to create aysmmetric pads around low 
      if(j==0){
	pad[i][j] = new TPad(padName.Data(),padName.Data(),
			     Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
      }else{
	pad[i][j] = new TPad(padName.Data(),padName.Data(),
			     Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
      }

      pad[i][j]->SetLeftMargin(leftMargin);
      pad[i][j]->SetRightMargin(rightMargin);
      pad[i][j]->SetBottomMargin(bottomMargin);
      pad[i][j]->SetTopMargin(topMargin);
      /*
      if(i==0 || i%columns==0){
      pad[i][j]->SetLeftMargin(leftMargin);
	pad[i][j]->SetRightMargin(rightMargin);
      }
      else pad[i][j]->SetRightMargin(edge+0.01);

      if(j==(rows-1)){
	pad[i][j]->SetBottomMargin(bottomMargin);
	pad[i][j]->SetTopMargin(edge);
      }else{
	pad[i][j]->SetBottomMargin(edge);
	pad[i][j]->SetTopMargin(topMargin);
      }
      */

      /*
       if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
       else pad[i][j]->SetLeftMargin(0.17);

       if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
       else pad[i][j]->SetRightMargin(0.1);

       if(j==0) pad[i][j]->SetTopMargin(edge);
       else pad[i][j]->SetTopMargin(0);

       if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
       else pad[i][j]->SetBottomMargin(0);
      */

       pad[i][j]->Draw();
       pad[i][j]->cd();
       pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
} // end makeMultiPanelCanvasXMerged


void putAorB(const char *text){
  TLatex * texB;
  texB = new TLatex(0.2,0.88,text);
  texB->SetTextSize(0.045);
  texB->SetTextFont(42);
  texB->SetLineWidth(2);
  texB->SetNDC();
  texB->Draw();
}


void printCanvases(TCanvas * MyCanvas, const char * name, int log=0)
{
   MyCanvas->cd();
   MyCanvas->SetLogx(log);
   printf("canvas name: %s\n",name);

   // add some text labels  
   //double ndcX = 0.25;
   //double ndcY = 0.9;
   //TLatex * tex;

   //MyCanvas->Print(Form("./fig/%s.eps",name));
   MyCanvas->Print(Form("./fig/%s.gif",name));
   MyCanvas->Print(Form("./fig/%s.pdf",name));
}


void th1Style1( TGraph *g=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1){
   g->SetMarkerStyle(marker);
   if(mcol>95 || mcol<14){
      g->SetMarkerColor(mcol);
   }else{
      g->SetMarkerColor(color[mcol-14]);
   }
   g->SetMarkerSize(msize);
   if(lcol>95 || lcol<14){
      g->SetLineColor(lcol);
   }else{
      g->SetLineColor(color[lcol-14]);
   }
   g->SetLineWidth(lsize);
   g->SetLineStyle(lstyle);

   if(dStyle==1) g->Draw("PZsame");
   if(dStyle==2) g->Draw("Xlsame");
   if(dStyle==3) g->Draw("histsame");
   if(dStyle==4) g->Draw("plxsame");
}


void th1Style1( TH1 *h=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1){
   h->SetMarkerStyle(marker);
   if(mcol>95 || mcol<14){
      h->SetMarkerColor(mcol);
   }else{
      h->SetMarkerColor(color[mcol-14]);
   }
   h->SetMarkerSize(msize);
   if(lcol>95 || lcol<14){
      h->SetLineColor(lcol);
   }else{
      h->SetLineColor(color[lcol-14]);
   }
   h->SetLineWidth(lsize);
   h->SetLineStyle(lstyle);

   if(dStyle==1) h->Draw("PZsame");
   if(dStyle==2) h->Draw("Csame");
   if(dStyle==3) h->Draw("histsame");
   if(dStyle==4) h->Draw("PCsame");
}


void th1Style2( TH1 *h=0, int mcol = 1, int marker = 20, float msize = 1, int lcol = 1, float lsize = 1, int lstyle = 1, int dStyle=1, bool drawstat=false){
  h->SetMarkerStyle(marker);
  if(mcol>95 || mcol<14){
    h->SetMarkerColor(mcol);
  }else{
    h->SetMarkerColor(color[mcol-14]);
  }
  h->SetMarkerSize(msize);
  if(lcol>95 || lcol<14){
    h->SetLineColor(lcol);
  }else{
    h->SetLineColor(color[lcol-14]);
  }
  h->SetLineWidth(lsize);
  h->SetLineStyle(lstyle);

  if(drawstat){
    if(dStyle==1) h->Draw("PZsames");
    if(dStyle==2) h->Draw("Csames");
    if(dStyle==3) h->Draw("histsames");
    if(dStyle==4) h->Draw("PCsames");
  }else{
    if(dStyle==1) h->Draw("PZsame");
    if(dStyle==2) h->Draw("Csame");
    if(dStyle==3) h->Draw("histsame");
    if(dStyle==4) h->Draw("PCsame");
  }
}


void setPad1by2(TCanvas *can){
  can->cd();

  TPad *pp1=0, *pp2=0;
  pp1 = new TPad("p1","",0,0.34,1,1,0,0,0);
  pp1->SetBottomMargin(0.0);
  pp1->SetTopMargin(0.05*(1/0.72));
  pp1->Draw();
  pp1->cd();
  pp1->SetNumber(1);
  can->cd();

  pp2 = new TPad("pp2","",0,0.0,1,0.34,0,0,0);
  pp2->SetTopMargin(0);
  pp2->SetBottomMargin(0.14*(1/0.34));
  pp2->Draw();
  pp2->cd();
  pp2->SetNumber(2);
}


void setLowerPad1by2(TH1D *hlow){
  hlow->GetYaxis()->SetLabelFont(43);
  hlow->GetYaxis()->SetTitleFont(43);
  hlow->GetYaxis()->SetLabelSize(20);
  hlow->GetYaxis()->SetTitleSize(22);
  hlow->GetYaxis()->SetTitleOffset(2.1);
  hlow->GetXaxis()->SetLabelFont(43);
  hlow->GetXaxis()->SetTitleFont(43);
  hlow->GetXaxis()->SetLabelSize(20);
  hlow->GetXaxis()->SetTitleSize(22);
  hlow->GetXaxis()->SetTitleOffset(3.8);
  hlow->GetYaxis()->SetNdivisions(906);
}


void MakeHist(TH1 *histo,const char *xtitle,const char *ytitle){
  histo->SetStats(0);
  histo->SetTitle("");

  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->CenterTitle(true);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetXaxis()->SetTitleSize(0.07);
  histo->GetXaxis()->SetTitleOffset(1.15);
  histo->GetXaxis()->SetLabelSize(0.07);
  histo->GetXaxis()->SetLabelOffset(0.005);
  histo->GetXaxis()->SetNdivisions(507);

  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetYaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetTitleOffset(1.50);
  histo->GetYaxis()->SetLabelSize(0.07);
  histo->GetYaxis()->SetNdivisions(507);
}


void MakeHist2D(TH2 *histo,const char *xtitle,const char *ytitle){
  histo->SetStats(0);
  histo->SetTitle("");

  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->CenterTitle(true);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetXaxis()->SetTitleOffset(1.15);
  histo->GetXaxis()->SetLabelSize(0.06);
  histo->GetXaxis()->SetLabelOffset(0.005);
  histo->GetXaxis()->SetNdivisions(507);

  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetYaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitleSize(0.06);
  histo->GetYaxis()->SetTitleOffset(1.23);
  histo->GetYaxis()->SetLabelSize(0.06);
  histo->GetYaxis()->SetNdivisions(507);

  histo->GetZaxis()->SetTitleSize(0.06);
  histo->GetZaxis()->SetTitleOffset(1.50);
  histo->GetZaxis()->SetLabelSize(0.06);
  histo->GetZaxis()->SetNdivisions(507);
}


void LoadStyle(){
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetErrorX(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetPadBorderSize(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetPalette(1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
}


void putCMSPrel(double x, double y, double size){
  TLatex *tex=0;
  tex = new TLatex(x,y,"CMS Preliminary");
  tex->SetTextSize(size);
  tex->SetLineWidth(2);
  tex->SetNDC();
  tex->Draw();
}


void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}


TGraph* GetTGraphErrorBand(TH1D *hist, Double_t xoffset=0) {
   //returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   //xoffset: points at the edges can be moved such that band does
   //  not pass through the marker. It is assumed that points are 
   //  in increasing order in x.
   TGraph tmp(hist);
   TGraph *ing = &tmp;
   if(!ing)  return 0;
   Int_t n=ing->GetN();
   if(!n) return 0;
   TGraph *outg=new TGraph(2*n+1);
   Double_t *x=ing->GetX();
   Double_t *y=ing->GetY();
   Double_t *yerr=ing->GetEY();
   Double_t *yerrlow=ing->GetEYlow();
//   Double_t *yerrhi=ing->GetEYhigh();
   if(yerr && yerrlow) {
      Error("GetTGraphErrorBand","both yerr and yerrlow are here,"
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for(Int_t i=0;i<n;i++) {
      Double_t shiftedx=x[i];
      std::cout <<x[i]<<std::endl;
      if(i==0) shiftedx-=xoffset;
      if(i==n-1) shiftedx+=xoffset;
      outg->SetPoint(i,shiftedx,y[i]-hist->GetBinError(i+1));
      outg->SetPoint(2*n-1-i,shiftedx,y[i]+hist->GetBinError(i+1));
      
      outg->SetPoint(2*n,x[0],y[i]-hist->GetBinError(1)); 
      
   }
   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(ing->GetFillColor());
   outg->SetFillStyle(ing->GetFillStyle());
   return outg;
}


TGraph* GetErrorBand(TH1D *hist, Double_t ratio1, Double_t ratio2, Double_t xoffset) {
   //returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   //xoffset: points at the edges can be moved such that band does
   //  not pass through the marker. It is assumed that points are 
   //  in increasing order in x.
   TGraph tmp(hist);
   TGraph *ing = &tmp;
   if(!ing)  return 0;
   Int_t n= ing->GetN();
   std::cout <<"n = "<<n<<std::endl;
   if(!n) return 0;
   TGraph *outg=new TGraph(2*n+1-4);
   Double_t *x=ing->GetX();
   Double_t *y=ing->GetY();
   Double_t *yerr=ing->GetEY();
   Double_t *yerrlow=ing->GetEYlow();
//   Double_t *yerrhi=ing->GetEYhigh();
   if(yerr && yerrlow) {
      Error("GetTGraphErrorBand","both yerr and yerrlow are here,"
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for(Int_t i=1;i<n-1;i++) {
      Double_t shiftedx=x[i];
      std::cout <<x[i]<<std::endl;
      if(i==1) shiftedx-=xoffset;
      if(i==n-2) shiftedx+=xoffset;
      outg->SetPoint(i-1,shiftedx,y[i]*(1-ratio2));
      outg->SetPoint(2*n-4-i,shiftedx,y[i]*(1+ratio1));
      outg->SetPoint(2*n-4,x[1]-xoffset,y[i]*(1-ratio2)); 
      
   }
   outg->Print();
   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(17);
   outg->SetFillStyle(1002);
   return outg;
}


TGraph* GetErrorBand(TH1D *hist, Double_t *ratio1, Double_t *ratio2, Double_t xoffset,Int_t nBin=0) {
   //returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   //xoffset: points at the edges can be moved such that band does
   //  not pass through the marker. It is assumed that points are 
   //  in increasing order in x.
   TGraph tmp(hist);
   TGraph *ing = &tmp;
   if(!ing)  return 0;
   Int_t n=nBin;
   if (nBin==0) n=ing->GetN();
   std::cout <<"n = "<<n<<std::endl;
   if(!n) return 0;
   TGraph *outg=new TGraph(2*n);
   Double_t *x=ing->GetX();
   Double_t *y=ing->GetY();
   Double_t *yerr=ing->GetEY();
   Double_t *yerrlow=ing->GetEYlow();
//   Double_t *yerrhi=ing->GetEYhigh();
   if(yerr && yerrlow) {
      Error("GetTGraphErrorBand","both yerr and yerrlow are here,"
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for(Int_t i=0;i<n;i++) {
      Double_t shiftedx=x[i];
      std::cout <<x[i]<<std::endl;
      if(i==0) shiftedx-=xoffset;
      if(i==n-1) shiftedx+=xoffset;
      if(i==n) shiftedx+=xoffset;
      outg->SetPoint(i,shiftedx,y[i]*(1-ratio2[i]));
      outg->SetPoint(2*n-i-1,shiftedx,y[i]*(1+ratio1[i]));
//      outg->SetPoint(2*n,x[1]-xoffset,y[i]*(1-ratio2[i])); 
      
   }
   outg->Print();
   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(17);
   outg->SetFillStyle(1002);
   return outg;
}


TGraph* GetErrorBand2(TH1D *hist, Double_t *ratio1, Double_t *ratio2, Double_t xoffset) {
   //returns a band surrounding a TGraphErrors or TGraphAsymmErrors
   //xoffset: points at the edges can be moved such that band does
   //  not pass through the marker. It is assumed that points are 
   //  in increasing order in x.
   TGraph tmp(hist);
   TGraph *ing = &tmp;
   if(!ing)  return 0;
   Int_t n=ing->GetN();
   if(!n) return 0;
   TGraph *outg=new TGraph(2*n+1-8);
   std::cout <<n<<std::endl;
   Double_t *x=ing->GetX();
   Double_t *y=ing->GetY();
   Double_t *yerr=ing->GetEY();
   Double_t *yerrlow=ing->GetEYlow();
//   Double_t *yerrhi=ing->GetEYhigh();
   if(yerr && yerrlow) {
      Error("GetTGraphErrorBand","both yerr and yerrlow are here,"
	    " GetTGraphErrorsBand is written with wrong assumptions");
      return 0;
   }
   for(Int_t i=2;i<n-2;i++) {
      Double_t shiftedx=x[i];
      std::cout <<x[i]<<std::endl;
      if(i==2) shiftedx-=xoffset;
      if(i==n-3) shiftedx+=xoffset;
      outg->SetPoint(i-1-1,shiftedx,y[i]*(1-ratio2[i]));
      outg->SetPoint(2*n-6-i-1,shiftedx,y[i]*(1+ratio1[i]));
      outg->SetPoint(2*n-6-1,x[1]-xoffset,y[i]*(1-ratio2[i])); 
      
   }
   outg->Print();
   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(17);
   outg->SetFillStyle(1002);
   return outg;
}


TGraph* GetTGraphNoErrors(TGraph* ing) {
   //returns a TGraph with no errors
   if(!ing) return 0;

   Int_t n=ing->GetN();
   if(!n) return 0;

   TGraph *outg=new TGraph(n);
   Double_t *x=ing->GetX();
   Double_t *y=ing->GetY();
   
   for(Int_t i=0;i<n;i++) {
      outg->SetPoint(i,x[i],y[i]);
   }
   outg->SetLineColor(ing->GetLineColor());
   outg->SetLineWidth(ing->GetLineWidth());
   outg->SetLineStyle(ing->GetLineStyle());
   outg->SetFillColor(ing->GetFillColor());
   outg->SetFillStyle(ing->GetFillStyle());
   outg->SetMarkerColor(ing->GetMarkerColor());
   outg->SetMarkerSize(ing->GetMarkerSize());
   outg->SetMarkerStyle(ing->GetMarkerStyle());
   return outg;
}


void DrawTick(Double_t Y, Double_t Erru, Double_t Errd, Double_t X, Double_t Tick_x=0.012, Double_t Tick_y=1.8, Int_t TickColor=1){
  TLine *line;
  line = new TLine(X-Tick_y,Y+Erru,X-Tick_y,Y+Erru-Tick_x);
  line->SetLineColor(TickColor);
  line->Draw();
  line = new TLine(X+Tick_y,Y+Erru,X+Tick_y,Y+Erru-Tick_x);
  line->SetLineColor(TickColor);
  line->Draw();
  line = new TLine(X-Tick_y,Y+Erru,X+Tick_y,Y+Erru);
  line->SetLineColor(TickColor);
  line->Draw();
  
  line = new TLine(X-Tick_y,Y-Errd,X-Tick_y,Y-Errd+Tick_x);
  line->SetLineColor(TickColor);
  line->Draw();
  line = new TLine(X+Tick_y,Y-Errd,X+Tick_y,Y-Errd+Tick_x);
  line->SetLineColor(TickColor);
  line->Draw();
  line = new TLine(X-Tick_y,Y-Errd,X+Tick_y,Y-Errd);
  line->SetLineColor(TickColor);
  line->Draw();
}


void getLogBins(const Float_t lower, const Float_t higher, const Int_t nBins, Float_t bins[]){
  Float_t logBins[nBins+1];
  bins[0] = lower;
  bins[nBins] = higher;
  logBins[0] = TMath::Log10(lower);
  logBins[nBins] = TMath::Log10(higher);
  Float_t interval = (logBins[nBins] - logBins[0])/nBins;
  for(Int_t iter = 1; iter < nBins; iter++){
    logBins[iter] = logBins[0] + iter*interval;
    bins[iter] = TMath::Power(10, logBins[iter]);
  }
  //return;
}


void CorrectBinWidth(TH1F *&h1){
  for(int ix=1; ix<=h1->GetNbinsX(); ix++){
    float val    = h1->GetBinContent(ix);
    float valErr = h1->GetBinError(ix);
    float width  = h1->GetBinWidth(ix);
    if(val!=0){
      h1->SetBinContent(ix, val/width);
      h1->SetBinError(ix, valErr/width);
    }
  }
}


void CorrectBinWidth(TH1 *h1){
  for(int ix=1; ix<=h1->GetNbinsX(); ix++){
    float val    = h1->GetBinContent(ix);
    float valErr = h1->GetBinError(ix);
    float width  = h1->GetBinWidth(ix);
    if(val!=0){
      h1->SetBinContent(ix, val/width);
      h1->SetBinError(ix, valErr/width);
    }
  }
}


void setBinomialErrors(TGraphAsymmErrors* grEff,const TH1F* hnum, const TH1F* hden){
  for (int i=0;i<grEff->GetN();i++) {
    float nenum =hnum ->GetBinContent(i+1);
    float ndenom=hden->GetBinContent(i+1);
    float eeff=(ndenom>0.0) ? sqrt(nenum/(ndenom*ndenom)*(1-nenum/ndenom)):0.0;
    grEff->SetPointError(i,0,0,eeff,eeff);
  }
}


void MakeGraph(TGraphErrors *&gr,const char *xtitle, const char *ytitle){
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle(xtitle);
  gr->GetXaxis()->CenterTitle(true);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(1.15);
  gr->GetXaxis()->SetLabelSize(0.06);
  gr->GetXaxis()->SetLabelOffset(0.005);
  gr->GetXaxis()->SetNdivisions(507);

  gr->GetYaxis()->SetTitle(ytitle);
  gr->GetYaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetTitleOffset(1.35);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetNdivisions(507);
}


void set_range_truncatedRMS(TH1F *&hist,float frac){
  if (0==hist) return;

  const float nevts = hist->Integral(); if (0==nevts) return;
  const int   nbins = hist->GetNbinsX();

  if (frac<=0.0 || frac==1.) return;

  for (int ibin=1;ibin<nbins;++ibin) {
    int binx1   = ibin;
    int binx2   = nbins+1-ibin;
    float ievts = hist->Integral(binx1,binx2);

    if ( (ievts/nevts)>frac ) continue;
    else { hist->GetXaxis()->SetRange(binx1,binx2); break; }
  }
  return;
}


void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++)
    {
      Float_t val = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);
      val/=h->GetBinWidth(i);
      valErr/=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}


void multiplyBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++)
    {
      Float_t val = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);
      val*=h->GetBinWidth(i);
      valErr*=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

void adjust_fitrange(TH1 *&h,double& min,double& max){
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}


void MakeZero(TH1F *&h1){
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
}


void MakeHist(TH1F *&histo,int istat,const char *xname, const char *yname){
  histo->SetStats(istat);
  histo->SetMarkerStyle(24);
  histo->SetMarkerColor(1);
  histo->SetLineColor(1);
  histo->SetLineStyle(1);
  histo->GetXaxis()->SetTitle(xname);
  histo->GetXaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitle(yname);
  histo->GetYaxis()->CenterTitle(true);
}


void MakeHistMean(TH1F *&h1,float xmin, float xmax, float ymin,float ymax){
  h1->SetMinimum(ymin);
  h1->SetMaximum(ymax);
  h1->SetTitle("");
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetLabelOffset(0.005);
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitle("#mu");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
}


void MakeHistRMS(TH1F *&h1,float xmin, float xmax, float ymin,float ymax){
  h1->SetTitle("");
  h1->SetMinimum(ymin);
  h1->SetMaximum(ymax);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelOffset(0.01);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitle("#sigma / #mu");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelOffset(0.01);
  h1->GetYaxis()->SetLabelSize(0.09);
  h1->GetYaxis()->SetTitleSize(0.09);
  h1->GetYaxis()->SetTitleOffset(1.12);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetDecimals(true);
}


void drawText2(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}


TLegend *getLegend(double x1, double y1, double x2, double y2){
  TLegend *leg = new TLegend(x1,y1,x2,y2,NULL,"BRNDC");
  leg->SetHeader("");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(1001);
  return leg;
}


double roundoff(double val, int digi){
  return ceil( ( val * pow( 10, digi ) ) - 0.49 ) / pow( 10, digi );
}

