#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include <TVirtualFitter.h>

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"


using namespace std;

const int digi=3;

const int ncen=7;
const char *cdir[ncen]  = {"05","510","1030","3050","5070","7090","pp"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","pp"};


double ptbins[] = {60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

double xmin=ptbins[0];
double xmax=ptbins[nbins];

void  LoadStyle();
void drawText2(const char */*text*/, float /*xp*/, float /*yp*/, int /*size*/);
void makeMultiPanelCanvas(TCanvas*& /*canv*/,
                          const Int_t /*columns*/,
                          const Int_t /*rows*/,
                          const Float_t /*leftOffset*/,
                          const Float_t /*bottomOffset*/,
                          const Float_t /*leftMargin*/,
                          const Float_t /*bottomMargin*/,
                          const Float_t /*edge*/,const Float_t /*asyoffset*/); 
void MakeHist(TH1F *&/*hist*/,int /*istat*/, const char */*xname*/, const char */*yname*/);
void MakeHistRMS(TH1F *&/*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1F *&/*Mean*/,float /*max*/,float /*min*/);
void MakeZero(TH1F *&/*hist*/);
double roundoff(double /*val*/);
TLegend *getLegend(double /*x1*/, double /*y1*/, double /*x2*/, double /*y2*/);

void FitGauss(TH1F* inHist_p, double &mean, double &meanError, double &res, double &resError);


void fit_double_gaussian(TH1F *&hrsp);
int fit_dscb(TH1F *&hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg);

/// double sided crystal ball function definition  
double fnc_dscb(double*xx,double*pp);
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}


void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter);
void adjust_fitrange(TH1 *h,double& min,double& max);

void set_range_truncatedRMS(TH1F *&hist,float frac);
double fracRMS = 1.00;

void FitDist(TH1F *&/*hrsp*/, double &/*mean*/, double &/*emean*/, double &/*sig*/, double &/*esig*/);

int SmearFac(int wJetID=1, int radius = 3)
{
  int iSave=1;

  string  plotAlgo="PF";
  int id = wJetID;
  // int id=0;
  // if(plotAlgo == "Calo"){
  //   id=1;
  //   wJetID=0;
  // }

  string algname="Pu";

  const int knj = 1;
  std::string srad[knj]={Form("%d",radius)};

  cout <<" npt : " << nbins << endl;

  LoadStyle();

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  //! PbPb Reco
  //TFile *fin = new TFile("OutputHist_pawan_ntuples_pp_pbpb_ak3.root","r");
  //TFile *fin = new TFile("OutputHist_raghav_ntuples_pp_pbpb_ak3.root","r");

  TFile *finPbPb = TFile::Open(Form("PbPb_MC_new_akPu%d_20_eta_20.root",radius));
  TFile *finPP = TFile::Open(Form("../PP_MC/pp_MC_new_ak%d_20_eta_20.root",radius));

  TH1F *hrsp [knj][ncen][nbins];
  TH1F *hMean[knj][ncen], *hSigma[knj][ncen];
  TH1F *hJER[knj][ncen];

  TF1 *fReso = new TF1("fReso","sqrt(-1.0*pow([0],2)/x/x + (pow([1],2)/x)*pow(x,[2]) + pow([3],2))",0,800);
  fReso->SetParNames("N^{2}","S^{2}","m","C^{2}");
  fReso->SetParameters(0.01,0.01,0.01,0.01);
  fReso->SetLineColor(2);
  fReso->SetLineStyle(2);

  TH1F *hN = new TH1F("hN","N traditional noise",ncen,-0.5,ncen-0.5);
  hN->SetMarkerStyle(20);
  hN->SetMarkerColor(1);
  hN->SetLineColor(1);
  
  TH1F *hS  = new TH1F("hS","S : stochastic term",ncen,-0.5,ncen-0.5);
  hS->SetMarkerStyle(20);
  hS->SetMarkerColor(1);
  hS->SetLineColor(1);

  TH1F *hm  = new TH1F("hm","m : models the pT dependence of S and C",ncen,-0.5,ncen-0.5);
  hm->SetMarkerStyle(20);
  hm->SetMarkerColor(1);
  hm->SetLineColor(1);

  TH1F *hC  = new TH1F("hC","C : constant term",ncen,-0.5,ncen-0.5);
  hC->SetMarkerStyle(20);
  hC->SetMarkerColor(1);
  hC->SetLineColor(1);
  for(int ic=0;ic<ncen;ic++){
    cout<<"ic : " <<ic<<"\t label "<<ccent[ic]<<endl;
    hN->GetXaxis()->SetBinLabel(ic+1,ccent[(ncen-ic)-1]);
    hS ->GetXaxis()->SetBinLabel(ic+1,ccent[(ncen-ic)-1]);
    hm ->GetXaxis()->SetBinLabel(ic+1,ccent[(ncen-ic)-1]);
    hC ->GetXaxis()->SetBinLabel(ic+1,ccent[(ncen-ic)-1]);
  }


  TF1 *fgaus=0;

  for(int nj=0; nj<knj; nj++){
    for(int ic=0; ic<ncen; ic++){

      hMean[nj][ic] = new TH1F(Form("hMean%d_%d",nj,ic),Form("Mean %s %s",algname.c_str(),ccent[ic]),nbins,ptbins);
      hMean[nj][ic]->SetMarkerColor(1);
      hMean[nj][ic]->SetMarkerStyle(20);
      hMean[nj][ic]->SetLineColor(1);
      hMean[nj][ic]->SetMarkerSize(1.3);
      //MakeHistMean(hMean[nj][ic],1.022,0.978); 
      //MakeHistMean(hMean[nj][ic],1.082,0.898); 
      //MakeHistMean(hMean[nj][ic],1.148,0.848); 

      MakeHistMean(hMean[nj][ic],1.082,0.858); 
      //MakeHistMean(hMean[nj][ic],1.052,0.934); 
 
      hSigma[nj][ic] = new TH1F(Form("hSigma%d_%d",nj,ic),Form("Sigma %s %s",algname.c_str(),ccent[ic]),nbins,ptbins);
      hSigma[nj][ic]->SetMarkerColor(1);
      hSigma[nj][ic]->SetMarkerStyle(20);
      hSigma[nj][ic]->SetLineColor(1);
      hSigma[nj][ic]->SetMarkerSize(1.3);
      MakeHistRMS(hSigma[nj][ic],0.563,0.001); 

      hJER[nj][ic] = new TH1F(Form("hJER%d_%d",nj,ic),Form("JER %s %s",algname.c_str(),ccent[ic]),nbins,ptbins);
      hJER[nj][ic]->SetMarkerColor(1);
      hJER[nj][ic]->SetMarkerStyle(20);
      hJER[nj][ic]->SetLineColor(1);
      hJER[nj][ic]->SetMarkerSize(1.3);
      MakeHistRMS(hJER[nj][ic],0.563,0.001); 

      for(int ip=0; ip<nbins; ip++){

	if(ic<ncen-1)hrsp[nj][ic][ip]  = (TH1F*)finPbPb->Get(Form("hJER_%d_pt_%d_cent%d",(int)ptbins[ip], (int)ptbins[ip+1], ic));
        else hrsp[nj][ic][ip]  = (TH1F*)finPP->Get(Form("hJER_%d_pt_%d", (int)ptbins[ip], (int)ptbins[ip+1]));
	hrsp[nj][ic][ip]->Scale(1./hrsp[nj][ic][ip]->Integral());

	//hrsp[nj][ic][ip]->Rebin(3);
	double norm  = hrsp[nj][ic][ip]->GetMaximumStored();
	double err   = hrsp[nj][ic][ip]->GetStdDevError();
	double mean  = hrsp[nj][ic][ip]->GetMean();
	double emean = hrsp[nj][ic][ip]->GetMeanError();
	double sig   = hrsp[nj][ic][ip]->GetStdDev()/mean;
	double esig = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));

	double jer=hrsp[nj][ic][ip]->GetStdDev();
	double ejer=hrsp[nj][ic][ip]->GetStdDevError();

	fgaus = new TF1("fgaus","gaus", 0.50,1.50);
	fgaus->SetParameters(norm, 0.9999, 0.5);

	int fitstatus = 0;
	fitstatus = hrsp[nj][ic][ip]->Fit(fgaus,"RQ");
	// if(nj==0 && ic==ncen-1){
	//   cout << "*******START*************** "<< ccent[ic] << " " << ptbins[ip]  << " to " << ptbins[ip+1] << endl;
	//   fitstatus = hrsp[nj][ic][ip]->Fit(fgaus,"RLQ");
	//   cout <<" \t Fit Mean : " << fgaus->GetParameter(1) << " Hist Mean : " << mean << endl;
	//   cout << "*******END***************** "<< ccent[ic] << " " << ptbins[ip]  << " to " << ptbins[ip+1] << endl;
	//   cout << endl;
	// }
	// else fitstatus = hrsp[nj][ic][ip]->Fit(fgaus,"RLQ");
	//cout << " fitstatus : "<< fitstatus << " Fit Error : " << fgaus->GetParError(1) << " Hist Error " << hrsp[nj][ic][ip]->GetMeanError() << endl;
	//if(ic == ncen-1)fitstatus=-1;

	// mean  = (fitstatus!=0) ? hrsp[nj][ic][ip]->GetMean()     :  0.5*(fgaus->GetParameter(1) + hrsp[nj][ic][ip]->GetMean());
	// emean = (fitstatus!=0) ? hrsp[nj][ic][ip]->GetMeanError(): sqrt(pow(fgaus->GetParError(1),2) + pow(hrsp[nj][ic][ip]->GetMeanError(),2));
	mean  = (fitstatus!=0) ? hrsp[nj][ic][ip]->GetMean()     : fgaus->GetParameter(1);
	emean = (fitstatus!=0) ? hrsp[nj][ic][ip]->GetMeanError(): fgaus->GetParError(1);//hrsp[nj][ic][ip]->GetMeanError();
	sig   = (fitstatus!=0) ? hrsp[nj][ic][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
	esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ic][ip]->GetStdDevError(),2))+(pow(-hrsp[nj][ic][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1/mean,2)*pow(hrsp[nj][ic][ip]->GetStdDevError(),2))+(pow(-hrsp[nj][ic][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))); 
	//esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ic][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ic][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1/mean,2)*pow(hrsp[nj][ic][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ic][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 


	jer=(fitstatus!=0) ? hrsp[nj][ic][ip]->GetStdDev() : fgaus->GetParameter(2);
	ejer=esig;

	hMean [nj][ic]->SetBinContent (ip+1, mean);
	hMean [nj][ic]->SetBinError   (ip+1, emean);
	hSigma[nj][ic]->SetBinContent (ip+1, sig);
	hSigma[nj][ic]->SetBinError   (ip+1, esig);

	hJER[nj][ic]->SetBinContent (ip+1, jer);
	hJER[nj][ic]->SetBinError   (ip+1, ejer);
      }
    }
  }
  //return 0;
  gPad->Close();

  int maxc=6;
  int maxr=2;
  int ipad=0;

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  TLegend *leg0 = getLegend(0.4418803,0.3887016,0.8692251,0.6763714);  
  if( wJetID )leg0->SetHeader("w JetID");
  else leg0->SetHeader("w/o JetID");
  leg0->SetTextSize(0.10);
  leg0->AddEntry(hMean[0][ncen-1],Form("pp ak%sPF",srad[0].c_str()),"p");
  leg0->AddEntry(hMean[0][0],Form("PbPb akPu%sPF",srad[0].c_str()),"p");
  
  TLine *l0 = new TLine(xmin,1.00,xmax,1.0);
  l0->SetLineStyle(2);
  l0->SetLineWidth(2);
  TLine *l1 = new TLine(xmin,0.99,xmax,0.99);
  l1->SetLineStyle(2);
  TLine *l2 = new TLine(xmin,1.01,xmax,1.01);
  l2->SetLineStyle(2);

  TCanvas *c11[knj];
  for(int nj=0; nj<knj; nj++){
    c11[nj] = new TCanvas(Form("c11_%d",nj),Form("%s",(srad[nj]+plotAlgo).c_str()),1695,530);
    makeMultiPanelCanvas(c11[nj],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
    ipad=0;

    for(int ic=ncen-3; ic>=0; ic--){

      c11[nj]->cd(++ipad);
      gPad->SetLogx();

      hSigma[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hSigma[nj][ic]->Draw("p");

      hSigma[nj][ncen-1]->SetMarkerColor(2);
      hSigma[nj][ncen-1]->SetLineColor(2);
      hSigma[nj][ncen-1]->Draw("psame");

      if(ipad==1 && ic==5){
	drawText2("PYTHIA+HYDJET",0.35,0.85,21);      
	leg0->Draw();
      }
      if(ipad==2 && ic==4){
	drawText2("#sqrt{s_{NN}} = 2.76 TeV, |#eta| < 2",0.10,0.85,21);
      }
      drawText2(ccent[ic],0.28,0.70,21);

      c11[nj]->cd(ipad+maxc);
      gPad->SetLogx();
    
      if(ipad!=1)
      hMean[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hMean[nj][ic]->Draw("p");

      hMean[nj][ncen-1]->SetMarkerColor(2);
      hMean[nj][ncen-1]->SetLineColor(2);
      hMean[nj][ncen-1]->Draw("psame");

      l1->Draw();
      l0->Draw();
      l2->Draw();
    }
    if(iSave){
      c11[nj]->SaveAs(Form("JER_%d_PbPb_%s.gif",wJetID,(srad[nj]+plotAlgo).c_str()));
      c11[nj]->SaveAs(Form("JER_%d_PbPb_%s.pdf",wJetID,(srad[nj]+plotAlgo).c_str()));
    }
  }
  //return 0;

  ipad=0;
  TCanvas *c99[knj][ncen];
  for(int nj=0;nj<knj;nj++){
    for(int ic=ncen-3;ic>=0;ic--){
      c99[nj][ic] = new TCanvas(Form("c99%d_%d",nj,ic),Form("%s %s ",srad[nj].c_str(),ccent[ic]),100,102,1399,942);
      c99[nj][ic]->Divide(4,4,0,0);
      ipad=0;
      for(int ip=0;ip<nbins;ip++){      
	c99[nj][ic]->cd(++ipad);
	
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	
	// hrsp[nj][ic][ip]->Scale(1./hrsp[nj][ic][ip]->Integral());
	// if(ic==ncen-3)hrsp[nj][ncen-1][ip]->Scale(1./hrsp[nj][ncen-1][ip]->Integral());	

	double ymax = TMath::Max(hrsp[nj][ic][ip]->GetMaximum(),hrsp[nj][ncen-1][ip]->GetMaximum());
	//double ymax = hrsp[nj][ic][ip]->GetMaximum();
	hrsp[nj][ic][ip]->SetMaximum(1.45*ymax);
	hrsp[nj][ic][ip]->SetMinimum(0);
	hrsp[nj][ic][ip]->SetTitle(0);
	hrsp[nj][ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
	hrsp[nj][ic][ip]->GetXaxis()->SetTitleFont(42);
	hrsp[nj][ic][ip]->GetXaxis()->SetLabelFont(42);
	hrsp[nj][ic][ip]->GetXaxis()->SetLabelSize(0.08);
	hrsp[nj][ic][ip]->GetXaxis()->SetTitleSize(0.07);
	hrsp[nj][ic][ip]->GetYaxis()->SetTitle("");
	hrsp[nj][ic][ip]->GetYaxis()->SetTitleFont(42);
	hrsp[nj][ic][ip]->GetYaxis()->SetLabelFont(42);
	hrsp[nj][ic][ip]->GetYaxis()->SetLabelSize(0.08);
      
	hrsp[nj][ic][ip]->SetMarkerStyle(20);
	hrsp[nj][ic][ip]->SetMarkerColor(1);
	hrsp[nj][ic][ip]->SetLineColor(1);
	hrsp[nj][ic][ip]->SetMarkerSize(1.3);
	hrsp[nj][ic][ip]->Draw("p");  

	hrsp[nj][ncen-1][ip]->SetMarkerStyle(20);
	hrsp[nj][ncen-1][ip]->SetMarkerColor(2);
	hrsp[nj][ncen-1][ip]->SetLineColor(2);
	hrsp[nj][ncen-1][ip]->SetMarkerSize(1.3);
	hrsp[nj][ncen-1][ip]->Draw("psame");  
      
	std::ostringstream strs; 
	strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
	std::string spt = strs.str();

	if(ipad==1){
	  //drawText2(Form("%s",(srad[nj]+plotAlgo).c_str()),0.28,0.90,21);      

	  drawText2(ccent[ic],0.75,0.87,21);	  
	  drawText2(spt.c_str(),0.22,0.80,21);		

	  if(wJetID)drawText2("w JetID",0.22,0.65,21);
	  else drawText2("w/o JetID",0.22,0.66,21);
		
	} else drawText2(spt.c_str(),0.17,0.80,21);		
      }
      if(iSave){
	c99[nj][ic]->SaveAs(Form("IndDist_JER_PbPb_%s_%s.gif",ccent[ic],(srad[nj]+plotAlgo).c_str()));
	c99[nj][ic]->SaveAs(Form("IndDist_JER_PbPb_%s_%s.pdf",ccent[ic],(srad[nj]+plotAlgo).c_str()));
      }
    }
  }


  //! Fit Sigma
  ipad=0;
  int plotNJ=0;
  double min=40.;
  double max=340.;
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  TCanvas *c97 = new TCanvas("c97",Form("Sigma Fitting  %s",(srad[plotNJ]+plotAlgo).c_str()),121,112,1506,637);
  c97->Divide(4,2);
  for(int ic=ncen-1;ic>=0;ic--){
    
    if( ic == ncen-2 )continue; //! 90-100

    c97->cd(++ipad);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.22);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.2);

    hSigma[plotNJ][ic]->SetMaximum(0.43);
    hSigma[plotNJ][ic]->SetMinimum(0.053);
    hSigma[plotNJ][ic]->SetMarkerStyle(20);
    hSigma[plotNJ][ic]->SetMarkerColor(1);
    hSigma[plotNJ][ic]->SetLineColor(1);
    hSigma[plotNJ][ic]->SetMarkerSize(1.0);

    fReso->SetParameters(0.1,1,0.001,0.01);

    //fReso->SetParameters(0.1,1,0.001,0.01);
    // fReso->SetParLimits(0,0,0.5);
    // fReso->SetParLimits(1,0,15.0);
    // fReso->SetParLimits(2,0.0,0.008);
    // fReso->SetParLimits(3,0.,0.1);

    // fReso->SetParLimits(0,0,0.5);
    // fReso->SetParLimits(1,0,5.0);
    // fReso->SetParLimits(2,0.0,0.004);
    //fReso->SetParLimits(3,0.,0.1);

    hSigma[plotNJ][ic]->Draw("p");
    hSigma[plotNJ][ic]->Fit(fReso,"R+Q","",min,max+50);
    gPad->Update();

    cout<<"cent : "<<ccent[ic]<<"\t"<<endl
	<<"  \t p0   : "<<fReso->GetParameter(0)<<"\t +- \t "<<fReso->GetParError(0)<<endl
	<<"  \t p1   : "<<fReso->GetParameter(1)<<"\t +- \t "<<fReso->GetParError(1)<<endl
	<<"  \t p2   : "<<fReso->GetParameter(2)<<"\t +- \t "<<fReso->GetParError(2)<<endl
	<<"  \t p3   : "<<fReso->GetParameter(3)<<"\t +- \t "<<fReso->GetParError(3)<<endl
	<<"  \t chi2 : "<<fReso->GetChisquare() <<"\t NDF \t"<<fReso->GetNDF()<<endl
	<<"\t"<<endl;

    double val = sqrt(fReso->GetParameter(0));
    double err = val*0.5*(fReso->GetParError(0)/fReso->GetParameter(0));
    hN->SetBinContent(ic+1,val);
    hN->SetBinError  (ic+1,err);

    val = sqrt(fReso->GetParameter(1));
    err = val*0.5*(fReso->GetParError(1)/fReso->GetParameter(1));
    hS->SetBinContent(ic+1,val);
    hS->SetBinError  (ic+1,err);

    val = fReso->GetParameter(2);
    err = fReso->GetParError(2);
    hm->SetBinContent(ic+1,val);
    hm->SetBinError  (ic+1,err);

    val = sqrt(fReso->GetParameter(3));
    err = val*0.5*(fReso->GetParError(3)/fReso->GetParameter(3));
    hC->SetBinContent(ic+1,val);
    hC->SetBinError  (ic+1,err);

    TPaveStats *ps = (TPaveStats*) hSigma[plotNJ][ic]->GetListOfFunctions()->FindObject("stats");
    ps->SetX1NDC(0.37);
    ps->SetY1NDC(0.48);       
    ps->SetX2NDC(0.98);
    ps->SetY2NDC(0.99);
    ps->SetTextSize(0.045);
    ps->SetTextFont(42);
    ps->Draw();

    TPaveText *pt1 = new TPaveText(0.53546,0.3913862,0.8336132,0.4598043,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(10);
    pt1->SetTextFont(42);
    TText *tx1 = pt1->AddText(Form("%s",ccent[ic]));
    tx1->SetTextSize(0.08);
    pt1->Draw();
  }



  // TFile *fout = new TFile(Form("JERJES_R%s_PbPb_pp.root",srad[0].c_str()),"RECREATE");
  // fout->cd();
  // fout->mkdir(Form("ak%sJetAnalyzer",srad[0].c_str()));
  // fout->cd(Form("ak%sJetAnalyzer",srad[0].c_str()));
  // for(int nj=0; nj<knj; nj++){
  //   for(int ic=0; ic<ncen; ic++){
  //     if(ic == ncen-1 ){
  // 	hMean[nj][ic]->SetName(Form("hMean_pp_R%s_PF_%s",srad[nj].c_str(),cdir[ic]));
  // 	hMean[nj][ic]->SetTitle(Form("hMean_pp_R%s_PF_%s",srad[nj].c_str(),cdir[ic]));
  // 	hSigma[nj][ic]->SetName(Form("hSigma_pp_R%s_PF_%s",srad[nj].c_str(),cdir[ic]));
  // 	hSigma[nj][ic]->SetTitle(Form("hSigma_pp_R%s_PF_%s",srad[nj].c_str(),cdir[ic]));
  //     }else{
  // 	hMean[nj][ic]->SetName(Form("hMean_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
  // 	hMean[nj][ic]->SetTitle(Form("hMean_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
  // 	hSigma[nj][ic]->SetName(Form("hSigma_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
  // 	hSigma[nj][ic]->SetTitle(Form("hSigma_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
  //     }
  //     hMean[nj][ic]->Write();
  //     hSigma[nj][ic]->Write();
  //   }
  // }
  // fout->Close();

  return 0;
}
void LoadStyle()
{
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
void drawText2(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge,const Float_t asyoffset) {
  if (canv==0) {
    cout<<"makeMultiPanelCanvas :  Got null canvas."<<endl;
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
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
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
      
      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);
      
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}
void MakeHist(TH1F *&histo,int istat,const char *xname, const char *yname)
{
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
void MakeHistMean(TH1F *&h1,float ymax,float ymin)
{
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
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
void MakeHistRMS(TH1F *&h1,float ymax,float ymin)
{

  h1->SetTitle("");
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
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

void fit_gaussian(TH1F *&hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();
  double ptRefMax(1.0),rspMax(0.0);

  //cout << " Mean : "  << mean << " \t  RMS  : " << rms << endl;

  double norm  = hrsp->GetMaximumStored();
  double peak  = mean;
  double sigma = rms;

  rspMax =jtptmin/ptRefMax;

  double xmin  = hrsp->GetXaxis()->GetXmin();
  double xmax  = hrsp->GetXaxis()->GetXmax();
  TF1* fitfnc(0); int fitstatus(-1);
  for (int iiter=0;iiter<niter;iiter++) {
    vector<double> vv;
    vv.push_back(rspMax);
    vv.push_back(xmin);
    vv.push_back(peak-nsigma*sigma);
    double fitrange_min = *std::max_element(vv.begin(),vv.end());
    double fitrange_max = std::min(xmax,peak+nsigma*sigma);

    adjust_fitrange(hrsp,fitrange_min,fitrange_max);
    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    fitstatus = hrsp->Fit(fitfnc,"RQ0");
    delete fitfnc;
    fitfnc = hrsp->GetFunction("fgaus");
    //fitfnc->ResetBit(TF1::kNotDraw);
    norm  = fitfnc->GetParameter(0);
    peak  = fitfnc->GetParameter(1);
    sigma = fitfnc->GetParameter(2);
  }
  if(hrsp->GetFunction("fgaus")==0){
      cout << "No function recorded in histogram " << hrsp->GetName() << endl;
    }
  //   if (0==fitstatus){
  //     cout<<"fit_gaussian() to "<<hrsp->GetName()    <<"  sucessful : " <<endl;
  //   }          



  if (0!=fitstatus){
    cout<<"fit_gaussian() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}
void adjust_fitrange(TH1 *h,double& min,double& max)
{
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}
void MakeZero(TH1F *&h1)
{
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
}
int fit_dscb(TH1F *&hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_dscb()"<<endl;return -1;
  }


  // first use a gaussian to constrain crystal ball gaussian core
  fit_gaussian(hrsp, nsigma, jtptmin, niter);
  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold 
  string histname = hrsp->GetName();
  //double ptRefMax(1.0),rspMax(0.0);

  double fitrange_min(0.4);
  double fitrange_max(1.6);

  //cout <<" \t \t  xmin : "  << fitrange_min << "\t" << fitrange_max << endl;


  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);
  fdscb->SetLineWidth(2);
  fdscb->SetLineStyle(2);

  double norm = fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);

  double aone(2.0),atwo(2.0),pone(10.0),ptwo(10.0);
  TVirtualFitter::SetDefaultFitter("Minuit");
  
  int fitstatus(0);
  for (int i=0;i<niter;i++) {
    fdscb->SetParameter(0,norm); // N
    fdscb->SetParameter(1,mean); // mean
    fdscb->SetParameter(2,sigma);// sigma
    fdscb->SetParameter(3,aone); // a1
    fdscb->SetParameter(4,pone); // p1
    fdscb->SetParameter(5,atwo); // a2
    fdscb->SetParameter(6,ptwo); // p2                

    fdscb->FixParameter(1,mean);
    fdscb->FixParameter(2,sigma);
    
    if (i>0) fdscb->FixParameter(3,aone);
    else fdscb->SetParLimits(3,1.,50.);
    
    if (i>1) fdscb->FixParameter(5,atwo);
    else fdscb->SetParLimits(5,1.,50.);

    fdscb->SetParLimits(4,0.,100.);
    fdscb->SetParLimits(6,0.,100.);

    fitstatus = hrsp->Fit(fdscb,"RQB+");
    if (0==fitstatus) i=999;
    delete fdscb;
    fdscb = hrsp->GetFunction("fdscb");
    if (0==fdscb) return -1;

    norm  = fdscb->GetParameter(0);
    aone  = fdscb->GetParameter(3);
    pone  = fdscb->GetParameter(4);
    atwo  = fdscb->GetParameter(5);
    ptwo  = fdscb->GetParameter(6);

    // reset sigma and mean to gauss values...
    fdscb->SetParameter(1,fgaus->GetParameter(1));
    fdscb->SetParError (1,fgaus->GetParError(1));
    fdscb->SetParameter(2,fgaus->GetParameter(2));
    fdscb->SetParError (2,fgaus->GetParError(2));
  }
  if (0!=fitstatus){
    cout<<"fit_fdscb() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus<<endl;
    hrsp->GetFunction("fdscb")->Delete();
  }
  return fitstatus;
}
double fnc_dscb(double*xx,double*pp)
{
  double x   = xx[0];
  // gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance

  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];

  double u   = (x-mu)/sig;
  double A1  = pow(p1/fabs(a1),p1)*exp(-a1*a1/2);
  double A2  = pow(p2/abs(a2),p2)*exp(-a2*a2/2);
  double B1  = p1/fabs(a1) - fabs(a1);
  double B2  = p2/fabs(a2) - fabs(a2);

  double result(N);
  if      (u<-a1) result *= A1*pow(B1-u,-p1);
  else if (u<a2)  result *= exp(-u*u/2);
  else            result *= A2*pow(B2+u,-p2);
  return result;
}
void fit_double_gaussian(TH1F *&hrsp)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_double_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  hrsp->Scale(1./hrsp->Integral());
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();

  int maxbin    = hrsp->GetMaximumBin();
  double norm1  = hrsp->GetBinContent(maxbin);
  double peak1  = hrsp->GetBinCenter(maxbin);
  double sigma1 = 0.04;

  double norm2  = norm1/2.0;
  double peak2  = mean;
  double sigma2 = 2*rms;

  //cout << " Mean  : "  << mean  << " \t  RMS  : " << rms    << endl;
  //cout << " norm1 : "  << norm1 << " \t  norm2 : " << norm2 << endl;
  //cout << " peak1 : "  << peak1 << " \t  sig1 : " << sigma1 << endl;
  //cout << " peak2 : "  << peak2 << " \t  sig2 : " << sigma2 << endl;

  double fitrange_min = 0.4;
  double fitrange_max = 1.4;

  TF1* fitfnc(0); int fitstatus(-1);
  TF1 *fitg1(0), *fitg2(0);
  fitfnc = new TF1("fdgaus","gaus(0)+gaus(3)",fitrange_min,fitrange_max);
  fitfnc->SetLineColor(1);
  fitfnc->SetLineStyle(2);
  fitfnc->SetLineWidth(2);

  fitfnc->SetParNames("N_{1}", "#mu_{1}", "#sigma_{1}",
		      "N_{2}", "#mu_{2}", "#sigma_{2}");
  fitfnc->SetParameters(norm1, peak1, sigma1, 
   			norm2, peak2, sigma2); 
  fitstatus = hrsp->Fit(fitfnc,"RQ");

  fitfnc->SetParLimits(0,0.01,50*norm1);
  // fitfnc->SetParLimits(1,0.7,1.2);
  // fitfnc->SetParLimits(2,0.01,5.0);

  //fitfnc->SetParLimits(3,0.0,2*norm2);
  // fitfnc->SetParLimits(4,0.2,1.7);
  // fitfnc->SetParLimits(5,1.0,10.0);

  //fitfnc->SetParLimits(4,peak2-3.0*sigma2,peak2+3.0*sigma2);
  //fitfnc->SetParLimits(5,0.10,2.0*sigma2);


  // if (0!=fitstatus){
  //   fitfnc->SetParLimits(4,0.2,1.7);
  //   fitfnc->SetParLimits(5,2.5,20.0);
  //   //cout <<" Not able to Fit this pt bin " << hrsp->GetName() << endl;
  // }
  fitstatus = hrsp->Fit(fitfnc,"RQ");
  hrsp->SetMaximum(norm1+0.2*norm1);
  fitg1 = new TF1("fg1","gaus(0)",fitrange_min,fitrange_max);
  fitg1->SetParameters(fitfnc->GetParameter(0),
		       fitfnc->GetParameter(1),
		       fitfnc->GetParameter(2));
  fitg1->SetLineColor(2);
  fitg1->SetLineStyle(2);
  hrsp->GetListOfFunctions()->Add(fitg1);

  fitg2 = new TF1("fg2","gaus(0)",fitrange_min,fitrange_max);
  fitg2->SetParameters(fitfnc->GetParameter(3),
		       fitfnc->GetParameter(4),
		       fitfnc->GetParameter(5));
  fitg2->SetLineColor(4);
  fitg2->SetLineStyle(4);
  hrsp->GetListOfFunctions()->Add(fitg2);

  // if(hrsp->GetFunction("fdgaus")==0){
  //   cout << "No function recorded in histogram " << hrsp->GetName() << endl;
  // }
  // if (0!=fitstatus){
  //   cout<<"fit_double_gaussian() to "<<hrsp->GetName()
  //       <<" failed. Fitstatus: "<<fitstatus
  //       <<" - FNC deleted."<<endl;
  //   hrsp->GetListOfFunctions()->Delete();
  // }
}
void set_range_truncatedRMS(TH1F *&hist,float frac)
{
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
double roundoff(double val)
{
  return ceil( ( val * pow( 10, digi ) ) - 0.49 ) / pow( 10, digi );
}
TLegend *getLegend(double x1, double y1, double x2, double y2)
{
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

void FitDist(TH1F *&hrsp, double &mean, double &emean, double &sig, double &esig)
{
  //hrsp->Rebin(2);
  set_range_truncatedRMS(hrsp, fracRMS);

  double norm  = hrsp->GetMaximumStored();
  mean  = hrsp->GetMean();
  emean = hrsp->GetMeanError();
  double rms   = hrsp->GetRMS();
  double err   = hrsp->GetRMSError();
  sig   = hrsp->GetRMS()/mean;
  esig = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));
  esig = sqrt(esig);
      
  //TF1 *fgaus = new TF1("fgaus","gaus", mean - 1.50*rms, mean + 1.50*rms);
  TF1 *fgaus = new TF1("fgaus","gaus", 0.65, 1.30);
  fgaus->SetParameters(norm, 1.0, 0.5);
  fgaus->SetParLimits(1,0.38,2.00);
  fgaus->SetParLimits(2,0.2,1.00);
  int fitstatus  = hrsp->Fit(fgaus,"RMLQ");
  //fitstatus=-1;
  // mean  = (fitstatus!=4000) ? hrsp->GetMean()     :  0.5*(fgaus->GetParameter(1) + hrsp->GetMean());
  // emean = (fitstatus!=4000) ? hrsp->GetMeanError(): sqrt(pow(fgaus->GetParError(1),2) + pow(hrsp->GetMeanError(),2));
  mean  = (fitstatus!=4000) ? hrsp->GetMean()     : fgaus->GetParameter(1);
  emean = (fitstatus!=4000) ? hrsp->GetMeanError(): fgaus->GetParError(1);
  sig   = (fitstatus!=4000) ? hrsp->GetRMS()/mean :  fgaus->GetParameter(2)/fgaus->GetParameter(1);
  esig  = (fitstatus!=4000) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));

  hrsp->GetFunction("fgaus")->SetBit(TF1::kNotDraw);
  delete fgaus;

  // mean  = (fitstatus!=0) ? hrsp->GetMean()     :  ((fgaus->GetParameter(1)+hrsp->GetMean())*0.5);
  // emean = (fitstatus!=0) ? hrsp->GetMeanError():   0.5*sqrt(pow(hrsp->GetMeanError(),2)+pow(fgaus->GetParError(1),2));//fgaus->GetParError(1);
  // sig   = (fitstatus!=0) ? hrsp->GetRMS()/mean :  fgaus->GetParameter(2)/fgaus->GetParameter(1);
  // esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));

  //cout <<" fit status : " << fitstatus << endl;

  // fit_gaussian(hrsp, 1.0, 1.0, 5);      
  // TF1*  fgaus = (TF1*) hrsp->GetListOfFunctions()->Last();
  // mean  = (fgaus==0) ? hrsp->GetMean()     :  fgaus->GetParameter(1);      
  // emean = (fgaus==0) ? hrsp->GetMeanError():  fgaus->GetParError(1);
  // sig   = (fgaus==0) ? hrsp->GetRMS()/mean :  fgaus->GetParameter(2)/mean;
  // esig  = (fgaus==0) ? hrsp->GetRMSError() :  fgaus->GetParError(2);
	
  //mean  = (fgaus==0) ? hrsp->GetMean()     :  ((fgaus->GetParameter(1)+hrsp->GetMean())*0.5);
  //emean = (fgaus==0) ? hrsp->GetMeanError():   0.5*sqrt(pow(hrsp->GetMeanError(),2)+pow(fgaus->GetParError(1),2));//fgaus->GetParError(1);
  //esig  = (fgaus==0) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));
	
	
  //! Double sided crystal ball
  fit_dscb(hrsp, 1.0, 1.0, 15, Form("ak%d%s",nj+1,algname[nj].c_str()));
  TF1*  frelrsp = (TF1*) hrsp->GetListOfFunctions()->Last();
  //mean  = (frelrsp==0) ? hrsp->GetMean()     :  frelrsp->GetParameter(1);
  //emean = (frelrsp==0) ? hrsp->GetMeanError():  frelrsp->GetParError(1);
  mean  = (frelrsp==0) ? hrsp->GetMean()     :   ((frelrsp->GetParameter(1)+hrsp->GetMean())*0.5);
  emean = (frelrsp==0) ? hrsp->GetMeanError():   0.5*sqrt(pow(hrsp->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
  sig   = (frelrsp==0) ? hrsp->GetRMS()/mean :  frelrsp->GetParameter(2)/frelrsp->GetParameter(1);
  esig  = (frelrsp==0) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
  	: sqrt((pow(1./frelrsp->GetParameter(1),2)*pow(frelrsp->GetParError(2),2))+pow(-frelrsp->GetParameter(2)/pow(frelrsp->GetParameter(1),2),2)*pow(frelrsp->GetParError(1),2));	
  //sig   = (frelrsp==0) ? hrsp->GetRMS()/mean :  frelrsp->GetParameter(2);
  //esig  = (frelrsp==0) ? hrsp->GetRMSError() :  frelrsp->GetParError(2);
	
  //peak    =(frelrsp==0) ? hrsp[i][j]->GetMean()     : ((frelrsp->GetParameter(1)+hrsp[i][j]->GetMean())*0.5);
  //epeak   =(frelrsp==0) ? hrsp[i][j]->GetMeanError(): 0.5*sqrt(pow(hrsp[i][j]->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
	
	
  //fit_double_gaussian(hrsp);
  // TF1*  frelrsp = (TF1*) hrsp->GetFunction("fgaus");
  // mean         = (frelrsp==0) ? hrsp->GetMean()     :  frelrsp->GetParameter(1);
  // double emean = (frelrsp==0) ? hrsp->GetMeanError():  frelrsp->GetParError(1);
  // double sig   = (frelrsp==0) ? hrsp->GetRMS()      :  frelrsp->GetParameter(2);
  // double esig  = (frelrsp==0) ? hrsp->GetRMSError() :  frelrsp->GetParError(2);


}
void FitGauss(TH1F* inHist_p, double &mean, double &meanError, double &res, double &resError)
{
  inHist_p->Fit("gaus", "Q L M", "");
  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);
  Float_t prob = inHist_p->GetFunction("gaus")->GetProb();
  // if(TMath::Abs(1.00 - mean) < .01) return;
  Int_t meanBin = inHist_p->FindBin(mean);
  Float_t meanRMS = -1;
  Float_t total = inHist_p->Integral();
  for(Int_t iter = 0; iter < inHist_p->GetNbinsX(); iter++){
    Int_t lowBound = 0;
    if(meanBin - iter > 0) lowBound = meanBin - iter;
    if(inHist_p->Integral(lowBound, meanBin + iter)/total > .95 || lowBound == 0 || inHist_p->GetBinContent(lowBound) < .01){
      meanRMS = inHist_p->GetBinCenter(meanBin + iter) - inHist_p->GetBinCenter(meanBin);
      break;
    }
  }
  double minPt = inHist_p->GetBinCenter(meanBin) - meanRMS;
  double maxPt = inHist_p->GetBinCenter(meanBin) + meanRMS;
  minPt = std::max(std::max(minPt, 0.0), inHist_p->GetXaxis()->GetXmin());
  maxPt = std::min(maxPt, inHist_p->GetXaxis()->GetXmax());
  inHist_p->Fit("gaus", "Q L M", "", minPt, maxPt);
  if(TMath::Abs(1.00 - mean) < TMath::Abs(1.00 - inHist_p->GetFunction("gaus")->GetParameter(1)) && prob > 0.0001){
    inHist_p->Fit("gaus", "Q L M", "");
    return;
  }
  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);
  return;
}
