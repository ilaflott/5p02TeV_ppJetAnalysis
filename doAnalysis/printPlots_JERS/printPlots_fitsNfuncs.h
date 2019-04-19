// C/C++
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>

// file IO/settings
#include <TROOT.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TFile.h>

// draw/print
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
//#include <TPaveText.h>
//#include <TPaveStats.h>
#include <TLatex.h>
//#include <TString.h>

// histsNfits
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

const int digi = 3;
const double fracRMS = 1.00;
const double pi=3.14159265359;

// fit misc ----------------------------------------------------------------------------
void FitDist(TH1F *&, //hrsp
	     double &, double &, double &, double &);//mean,emean,sig,esig
void set_range_truncatedRMS(TH1F *&hist,float frac);
void adjust_fitrange(TH1 *h, 
		     double& min,double& max);

// the methods 
void FitDist(TH1F *&hrsp, 
	     double &mean, double &emean, double &sig, double &esig){

  //hrsp->Rebin(2);
  set_range_truncatedRMS(hrsp, fracRMS);

  double norm  = hrsp->GetMaximumStored();
  mean  = hrsp->GetMean();
  emean = hrsp->GetMeanError();

  //double rms   = hrsp->GetRMS();
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

  //std::cout <<" fit status : " << fitstatus << std::endl;

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
  // fit_dscb(hrsp, 1.0, 1.0, 15, Form("ak%d%s",nj+1,algname[nj].c_str()));
  // TF1*  frelrsp = (TF1*) hrsp->GetListOfFunctions()->Last();
  // //mean  = (frelrsp==0) ? hrsp->GetMean()     :  frelrsp->GetParameter(1);
  // //emean = (frelrsp==0) ? hrsp->GetMeanError():  frelrsp->GetParError(1);
  // mean  = (frelrsp==0) ? hrsp->GetMean()     :   ((frelrsp->GetParameter(1)+hrsp->GetMean())*0.5);
  // emean = (frelrsp==0) ? hrsp->GetMeanError():   0.5*sqrt(pow(hrsp->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
  // sig   = (frelrsp==0) ? hrsp->GetRMS()/mean :  frelrsp->GetParameter(2)/frelrsp->GetParameter(1);
  // esig  = (frelrsp==0) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
  // 	: sqrt((pow(1./frelrsp->GetParameter(1),2)*pow(frelrsp->GetParError(2),2))+pow(-frelrsp->GetParameter(2)/pow(frelrsp->GetParameter(1),2),2)*pow(frelrsp->GetParError(1),2));
	
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


void set_range_truncatedRMS(TH1F *&hist,float frac){

  if (0==hist) return;

  const float nevts = hist->Integral(); 
  if (0.==nevts) return;
  
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


void adjust_fitrange(TH1 *h,double& min,double& max){
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}


// styledraw/hist/visuals methods -----------------------------------------------------------------------
void LoadStyle();
void makeMultiPanelCanvas(TCanvas*& , const Int_t , const Int_t , // canv, columns, rows
                          const Float_t , const Float_t , // left/bottom Offset
                          const Float_t , const Float_t, // left/bottom Margin
                          const Float_t , const Float_t );// edge, asyoffset

//void MakeHist(TH1F *&, int,//hist,istat 
	      //	      const char*, const char*);//xname,yname

void MakeHist( TH1F* ,
	       float, float);


//void MakeHistMean(TH1F * , // mean 
//		 float , float);// max, min 
void MakeHistMu(TH1F * , // mean 
		 float , float);// max, min 
//void MakeHistRMS(TH1F *  , // hRMS 
//		 float , float);// max, min 
void MakeHistSigma(TH1F *  , // hRMS 
		   float , float);// max, min 
void MakeHistSigmaMu(TH1F *  , // hRMS 
		     float , float);// max, min 
void MakeHistChi2NDF(TH1F *  , // hRMS 
		  float , float);// max, min 



void MakeHistRMSEta(TH1F *&  , // hRMS 
		 float , float);// max, min 
void MakeHistRMSAbsEta(TH1F *&  , // hRMS 
		 float , float);// max, min 
void MakeHistMeanEta(TH1F *& , // mean 
		 float , float);// max, min 




void MakeZero(TH1F *& ); //hist 
void divideBinWidth(TH1*);
void drawText2(const char *,  //text 
	       float, float , int ); //xp,yp,size

TLegend* getLegend(double , double ,//x1,y1
		   double , double );//x2,y2

// the methods
void LoadStyle(){
  
  std::cout<<std::endl<<"loading style..." <<std::endl;
  
  gStyle->SetOptStat(0);
  //gStyle->SetOptStat("emr");
  //gStyle->SetOptFit(1);
  //gStyle->SetOptFit(101);
  //gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);
  
  
  //color
  //gStyle->SetPalette(1);
  //gStyle->SetCanvasColor(10);
  //gStyle->SetFillColor(10);
  //gStyle->SetFrameFillColor(10);
  //gStyle->SetPadColor(10);
  //gStyle->SetStatColor(10);
  //gStyle->SetTitleFillColor(10);
  
  //pad
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);  
  
  //gStyle->SetPadBorderSize(0);
  //gStyle->SetPadBottomMargin(0.1);
  //gStyle->SetPadLeftMargin(0.1);

  //other
  //gStyle->SetLegendBorderSize(1);
  //gStyle->SetTitleBorderSize(1);
  //gStyle->SetErrorX(0);

  gROOT->ForceStyle();

  return;
}


void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, const Int_t rows,
                          const Float_t leftOffset, const Float_t bottomOffset,
                          const Float_t leftMargin, const Float_t bottomMargin,
                          const Float_t edge, const Float_t asyoffset  ){
  if (canv==0) {
    std::cout<<"makeMultiPanelCanvas :  Got null canvas."<<std::endl;
    return;
  }
  canv->Clear();
  
  Float_t PadWidth = (1.0-leftOffset)/
    ( (1.0/(1.0-leftMargin)) + (1.0/(1.0-edge)) + ((Float_t)columns-2.0) );
  Float_t PadHeight = (1.0-bottomOffset)/
    ( (1.0/(1.0-bottomMargin)) + (1.0/(1.0-edge)) + ((Float_t)rows-2.0) );

  Float_t Xlow[columns], Xup[columns];  
  Xlow[0] = leftOffset - asyoffset;
  Xup[0]  = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }

  Float_t Ylow[rows], Yup[rows];
  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  
  TPad* pad[columns][rows];  
  //TString padName;
  std::string padName="";
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName +="p_"+std::to_string(i)+"_"+std::to_string(j); //Form("p_%d_%d",i,j);
      pad[i][j] = new TPad( padName.c_str(), padName.c_str(),
			    Xlow[i],Ylow[j], Xup[i], Yup[j]);

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
  return;
}

//void MakeHist(TH1F *&histo,int istat,const char *xname, const char *yname){
void MakeHist(TH1F* h1, float xmin, float xmax){
  //h1->SetLineColor(kBlue);  

  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);    
  h1->SetMarkerSize(0.8);    

  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.05);//default is 0.08
  h1->GetYaxis()->SetTitleOffset(0.60);//no numerical default; auto calculated
  h1->GetYaxis()->SetLabelSize(.045);//default is 0.4

  h1->SetAxisRange(xmin,xmax,"X");
  h1->GetXaxis()->SetTitle("GEN Jet p_{T} (GeV/c)");
  h1->GetXaxis()->SetMoreLogLabels(true);
  h1->GetXaxis()->SetNoExponent(true);
  h1->GetXaxis()->SetLabelSize(.045);//default is 0.4
  //h1->GetXaxis()->SetLabelOffset(.01); //default is 0.005
  h1->GetXaxis()->SetTitleSize(0.045);  

  return;
}



void MakeHistMu(TH1F * h1, float xmin, float xmax){
  MakeHist((TH1F*)h1,xmin,xmax);
  //h1->SetTitle("mean, resolution, chi-square, v gen p_{T}");
  h1->SetTitle("");
  h1->GetYaxis()->SetTitle("#mu");
  h1->GetYaxis()->SetDecimals(true);
  h1->SetAxisRange(0.98,1.03, "Y");    
  return;
}

void MakeHistSigma(TH1F * h1,float xmin,float xmax){
  MakeHist((TH1F*)h1,xmin,xmax);
  h1->SetTitle("");
  h1->GetYaxis()->SetTitle("#sigma");
  h1->GetYaxis()->SetDecimals(true);
  h1->SetAxisRange(0.04,0.16, "Y");    

  return;
}
void MakeHistSigmaMu(TH1F * h1,float xmin,float xmax){
  MakeHist((TH1F*)h1,xmin,xmax);
  h1->SetTitle("");
  h1->GetYaxis()->SetTitle("#sigma / #mu");
  h1->GetYaxis()->SetDecimals(true);
  h1->SetAxisRange(0.04,0.16, "Y");    

  return;
}
void MakeHistChi2NDF(TH1F * h1, float xmin, float xmax){
  MakeHist((TH1F*)h1,xmin,xmax);
  h1->SetTitle("");
  //  h1->GetYaxis()->SetTitle("Gaussian JER #chi^{2}/NDF");
  h1->GetYaxis()->SetTitle("Gaussian Fit #chi^{2}/NDF");
  h1->SetAxisRange(1.e-1,1.e+2,"Y");
  h1->SetMarkerStyle(kOpenStar);
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kMagenta);
  return;
}


void MakeHistChi2Prob(TH1F * h1, float xmin, float xmax){
  MakeHist((TH1F*)h1,xmin,xmax);
  h1->SetTitle("");
  h1->GetYaxis()->SetTitle("Gaussian JER #chi^{2} Probability");
  //h1->SetAxisRange(1e-10 , 1., "Y");
  h1->SetMarkerStyle(kOpenStar);
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kMagenta);
  return;
}





void MakeHistMeanEta(TH1F *&h1,float ymax,float ymin){
  // hist
  //h1->SetMaximum(ymax);
  //h1->SetMinimum(ymin);
  //h1->SetTitle("");
  h1->SetLineColor(kBlue);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerSize(.8);

  // x axis
  //h1->GetXaxis()->SetRangeUser(xmin,xmax);
  //h1->GetXaxis()->SetTitle("gen #eta");
  //h1->GetXaxis()->CenterTitle(true);
  //h1->GetXaxis()->SetMoreLogLabels();
  //h1->GetXaxis()->SetNoExponent();

  //title/label

  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleSize(0.04);
  //h1->GetXaxis()->SetTitleOffset(1.0);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelSize(0.03);
  //h1->GetXaxis()->SetLabelOffset(0.005);

  //h1->GetXaxis()->SetNdivisions(507);

  // y axis
  h1->GetYaxis()->SetTitle("#mu");
  //h1->GetYaxis()->CenterTitle(true);
  //title/label
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.04);
  //h1->GetYaxis()->SetTitleOffset(1.00);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.03);

  //h1->GetYaxis()->SetNdivisions(507);
  //h1->GetYaxis()->SetDecimals(true);

}

void MakeHistRMSEta(TH1F *&h1,float ymax,float ymin){
  // hist
  h1->SetTitle("");
  //h1->SetMaximum(ymax);
  //h1->SetMinimum(ymin);

  h1->SetLineColor(kBlue);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerSize(.8);

  // x axis
  //h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("#eta");
  //h1->GetXaxis()->CenterTitle(true);
  //h1->GetXaxis()->SetMoreLogLabels();
  //h1->GetXaxis()->SetNoExponent();
  //title/label
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleSize(0.04);
  //h1->GetXaxis()->SetTitleOffset(1.0);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelSize(0.03);
  //h1->GetXaxis()->SetLabelOffset(0.01);
  
  //h1->GetXaxis()->SetNdivisions(507);

  // y axis
  h1->GetYaxis()->SetTitle("#sigma / #mu");
  //h1->GetYaxis()->CenterTitle(true);
  //title/label
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.04);
  //h1->GetYaxis()->SetTitleOffset(1.0);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.03);
  //h1->GetYaxis()->SetLabelOffset(0.01);

  //h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
}


void MakeHistRMSAbsEta(TH1F *&h1,float ymax,float ymin){
  // hist
  h1->SetTitle("");
  //h1->SetMaximum(ymax);
  //h1->SetMinimum(ymin);

  h1->SetLineColor(kBlue);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerSize(.8);

  // x axis
  //h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("|#eta|");
  //h1->GetXaxis()->CenterTitle(true);
  //h1->GetXaxis()->SetMoreLogLabels();
  //h1->GetXaxis()->SetNoExponent();
  //title/label
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleSize(0.04);
  //h1->GetXaxis()->SetTitleOffset(1.0);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelSize(0.03);
  //h1->GetXaxis()->SetLabelOffset(0.01);
  
  //h1->GetXaxis()->SetNdivisions(507);

  // y axis
  h1->GetYaxis()->SetTitle("#sigma / #mu");
  //h1->GetYaxis()->CenterTitle(true);
  //title/label
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.04);
  //h1->GetYaxis()->SetTitleOffset(1.0);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.03);
  //h1->GetYaxis()->SetLabelOffset(0.01);

  //h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
}

















void MakeZero(TH1F *&h1){
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
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


TLegend *getLegend(double x1, double y1, 
		   double x2, double y2){
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

// define gaussian fits
void FitGauss(TH1F* inHist_p, 
	      double &mean, double &meanError, double &res, double &resError);
void fit_gaussian(TH1F*& hrsp,
                  const double nsigma, const double jtptmin, const int niter);
void fit_double_gaussian(TH1F *&hrsp);

// the methods
void FitGauss(TH1F* inHist_p, 
	      double &mean, double &meanError, 
	      double &res, double &resError){

  inHist_p->Fit("gaus", "Q L M", "");

  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);

  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);

  Float_t prob = inHist_p->GetFunction("gaus")->GetProb();
  // if(TMath::Abs(1.00 - mean) < .01) return;

  Float_t total = inHist_p->Integral();
  Int_t meanBin = inHist_p->FindBin(mean);
  Float_t meanRMS = -1;
  for(Int_t iter = 0; iter < inHist_p->GetNbinsX(); iter++){
    Int_t lowBound = 0;
    if(meanBin - iter > 0) lowBound = meanBin - iter;
    if(inHist_p->Integral(lowBound, meanBin + iter)/total > .95 || 
       lowBound == 0 || 
       inHist_p->GetBinContent(lowBound) < .01){
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


void fit_gaussian(TH1F *&hrsp,
                  const double nsigma, const double jtptmin, const int niter){
  if (0==hrsp) {
    std::cout<<"ERROR: Empty pointer to fit_gaussian()"<<std::endl;return;
  }
  
  std::string histname = hrsp->GetName();
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();
  double ptRefMax(1.0),rspMax(0.0);
  rspMax =jtptmin/ptRefMax;
  //std::cout << " Mean : "  << mean << " \t  RMS  : " << rms << std::endl;

  double norm  = hrsp->GetMaximumStored();
  double peak  = mean;
  double sigma = rms;

  double min_x  = hrsp->GetXaxis()->GetXmin(); //=xmin;
  double max_x  = hrsp->GetXaxis()->GetXmax(); //=xmax;
  
  TF1* fitfnc(0); int fitstatus(-1);
  for (int iiter=0;iiter<niter;iiter++) {
    std::vector<double> vv;
    vv.push_back(rspMax);
    vv.push_back(min_x);
    vv.push_back(peak-nsigma*sigma);
    double fitrange_min = *std::max_element(vv.begin(),vv.end());
    double fitrange_max = std::min(max_x,peak+nsigma*sigma);
    adjust_fitrange(hrsp,fitrange_min,fitrange_max);
    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    fitstatus = hrsp->Fit(fitfnc,"RQ0");
    delete fitfnc;
    fitfnc = hrsp->GetFunction("fgaus");
    //fitfnc->ResetBit(TF1::kNotDraw;
    norm  = fitfnc->GetParameter(0);
    peak  = fitfnc->GetParameter(1);
    sigma = fitfnc->GetParameter(2);
  }
  
  if(hrsp->GetFunction("fgaus")==0) std::cout << "No function recorded in histogram " << hrsp->GetName() << std::endl;
  //   if (0==fitstatus) std::cout<<"fit_gaussian() to "<<hrsp->GetName()    <<"  sucessful : " <<std::endl;

  if (0!=fitstatus){
    std::cout<<"fit_gaussian() to " <<hrsp->GetName() <<" failed. " <<std::endl
	     <<" Fitstatus: " <<fitstatus <<" - FNC deleted." <<std::endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}


void fit_double_gaussian(TH1F *&hrsp){

  if (0==hrsp) {
    std::cout<<"ERROR: Empty pointer to fit_double_gaussian()"<<std::endl;return;
  }
  
  std::string histname = hrsp->GetName();
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

  //std::cout << " Mean  : "  << mean  << " \t  RMS  : " << rms    << std::endl;
  //std::cout << " norm1 : "  << norm1 << " \t  norm2 : " << norm2 << std::endl;
  //std::cout << " peak1 : "  << peak1 << " \t  sig1 : " << sigma1 << std::endl;
  //std::cout << " peak2 : "  << peak2 << " \t  sig2 : " << sigma2 << std::endl;

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
  std::cout << "fitstatus="<<fitstatus<<std::endl;
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
  //   //std::cout <<" Not able to Fit this pt bin " << hrsp->GetName() << std::endl;
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
  //   std::cout << "No function recorded in histogram " << hrsp->GetName() << std::endl;
  // }
  // if (0!=fitstatus){
  //   std::cout<<"fit_double_gaussian() to "<<hrsp->GetName()
  //       <<" failed. Fitstatus: "<<fitstatus
  //       <<" - FNC deleted."<<std::endl;
  //   hrsp->GetListOfFunctions()->Delete();
  // }
}

// double sided crystal ball fit+function ----------------------------------------------------
int fit_dscb(TH1F *&hrsp,
             const double nsigma, const double jtptmin, const int niter,
             const std::string alg);
double fnc_dscb(double*xx,double*pp);

// the methods
int fit_dscb(TH1F *&hrsp,
             const double nsigma, const double jtptmin, const int niter,
             const std::string alg){

  if (0==hrsp) {
    std::cout<<"ERROR: Empty pointer to fit_dscb()"<<std::endl;return -1;
  }

  // first use a gaussian to constrain crystal ball gaussian core
  fit_gaussian(hrsp, nsigma, jtptmin, niter);
  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold 
  std::string histname = hrsp->GetName();
  //double ptRefMax(1.0),rspMax(0.0);

  double fitrange_min(0.4);
  double fitrange_max(1.6);

  //std::cout <<" \t \t  xmin : "  << fitrange_min << "\t" << fitrange_max << std::endl;

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
    std::cout<<"fit_fdscb() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus<<std::endl;
    hrsp->GetFunction("fdscb")->Delete();
  }
  return fitstatus;
}


double fnc_dscb(double*xx,double*pp) {

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

// maths ----------------------------------------------------------------------------------------------------
int findBin(int );
float deltaR(float, float, float, float);
float deltaphi(float, float);
double roundoff(double    );//val  

// the methods
const int findBinArray[]={20,60,100,200};//const int findBinArray[]={20,60,80,100,120,140,160,180,200};
const int nbins_findBinArray=sizeof(findBinArray)/sizeof(int);
int findBin(int bin){
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins
  int ibin=-1;
  for(int i=0; i<nbins_findBinArray; ++i){
    if(bin<findBinArray[i])ibin=i;
    else continue;     }
  return ibin;
}

float deltaphi(float phi1, float phi2){
  //float pi=TMath::Pi();

  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}

float deltaR(float eta1, float phi1, float eta2, float phi2){
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}

double roundoff(double val){
  return ceil( ( val * pow( 10, digi ) ) - 0.49 ) / pow( 10, digi );
}


//// class T template defintion ?
//template <class T>
//bool from_string( T& t, const std::string& s,  std::ios_base& (*f)(std::ios_base&)){
//  std::istringstream iss(s);
//  return !(iss >> f >> t).fail();
//}

