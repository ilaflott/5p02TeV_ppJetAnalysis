// unfoldSpectra.h written by Ian W. Laflotte
// June 28 2016
// heavily based on unfold.C and utilities.h by Raghav Kunnawalkam Elayavalli
// provides functions and binning common to doMCClosureTests and unfoldDataSpectra
// compile doMCClosureTests or unfoldDataSpectra with...
// g++ -L"${ROOUNFOLD_PATH}" -I"${ROOUNFOLD_PATH}" -lRooUnfold -I"${ROOUNFOLD_PATH}/src" <target>.C $(root-config --cflags --libs) -Werror -Wall -O2 -o <target>.exe
// OR, if rooUnfoldCompile function is set in .bashrc
// rooUnfoldCompile <target>.C

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
const std::string CMSSW_BASE="/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/";
const std::string inFile_MC_dir="readFiles/ppMC/saved_outputCondor/readFiles_ppMC_5p02TeV_Py8_CUETP8M1_QCDjetAllPtBins_2016-07-09/";
const std::string inFile_MC_name="QCDjetAllPtBins_ak4PF-allFiles.root";//const std::string inFile_MC_name="readFiles_ppMC_QCDJets_forUnfoldTesting_101Files.root";
const std::string inFile_Data_dir="readFiles/ppData/saved_outputCondor/readFiles_ppData_5p02TeV_HighPtJetTrig_2016-06-10_allFiles/";
const std::string inFile_Data_name="HighPtJetTrig_ak4PF-allFiles.root";

// may make these inputs to unfoldDataSpectra ----------------------------------------------------------------------------
const std::string inFile_MC=CMSSW_BASE+inFile_MC_dir+inFile_MC_name;
const std::string inFile_Data=CMSSW_BASE+inFile_Data_dir+inFile_Data_name;
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


// CalculatePearsonCoefficients
//-----------------------------------------------------------------------------------------------------------------------
TMatrixD* CalculatePearsonCoefficients(TMatrixD* covmat) {

  //std::cout<<"calculating pearson coefficients!"<<std::endl;
  //std::cout<<"covmat="<< covmat<<std::endl;
  //std::cout<<"&covmat="<< &covmat<<std::endl;
  Int_t nrows = covmat->GetNrows();
  Int_t ncols = covmat->GetNcols();
  std::cout<<"nrows,ncols="<<nrows<<" , "<<ncols<<std::endl;
  TMatrixD* pearsonCoefs=new TMatrixD( nrows,ncols );

  //pearsonCoefs->Clear();
  //std::cout<<"pearsonCoefs="<< pearsonCoefs<<std::endl;
  //std::cout<<"&pearsonCoefs="<< &pearsonCoefs<<std::endl;

  // loop over columns+rows of the intput covariance matrix's pearson coefficients 
  for(int row = 0; row<nrows; row++){
    Double_t covmatElement_rowrow = (TMatrixD (*covmat))(row,row);
    for(int col = 0; col<ncols; col++){
      Double_t covmatElement_rowcol = (TMatrixD (*covmat))(row,col);
      Double_t covmatElement_colcol = (TMatrixD (*covmat))(col,col);
      Double_t pearson = covmatElement_rowcol/TMath::Sqrt(covmatElement_rowrow*covmatElement_colcol);
      (TMatrixD(*pearsonCoefs))(row,col) = pearson;
      if(pearson!=pearson)pearson=-1.;//nan
      if((row%10==0&&col%10==0)&&pearson!=-1)std::cout<<"for (row,col)= "<<row <<" , "<<col <<std::endl;
      if((row%10==0&&col%10==0)&&pearson!=-1)std::cout<<"pearson="<<pearson <<std::endl;
    }
  }
  return pearsonCoefs;
}


// PrintMatrix
//-----------------------------------------------------------------------------------------------------------------------
// Print the matrix as a table of elements.
// Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
// By default, format ="%11.4g", is used to print one element.
// One can specify an alternative format, e.g. format ="%6.2f  " 
int PrintMatrix( const TMatrixD& m, 
		  const char* format, const char* name, Int_t cols_per_sheet ){

  if (!m.IsValid()) {
    m.Error("PrintMatrix","%s is invalid",name);
    return -1;
  }

  const Int_t ncols  = m.GetNcols();
  const Int_t nrows  = m.GetNrows();
  const Int_t collwb = m.GetColLwb();
  const Int_t rowlwb = m.GetRowLwb();

  if (!(format && format[0])) format= "%11.4g ";
  char topbar[1000];
  snprintf(topbar,1000,format,123.456789);
  Int_t nch = strlen(topbar)+1;
  if (nch > 18) nch = 18;
  char ftopbar[20];
  for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
  Int_t nk = 1 + Int_t(log10(ncols));
  snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
  Int_t nch2 = strlen(ftopbar);
  for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
  ftopbar[nch] = '|';
  ftopbar[nch+1] = 0;

  printf("\n%dx%d %s is as follows",nrows,ncols,name);

  if (cols_per_sheet <= 0) {
    cols_per_sheet = 5;
    if (nch <= 8) cols_per_sheet =10;
  }
  nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
  for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
  topbar[nk] = 0;
  for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
    printf("\n\n     |");
    for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
      printf(ftopbar,j+collwb-1);
    printf("\n%s\n",topbar);
    if (m.GetNoElements() <= 0) continue;
    for (Int_t i=1 ; i<=nrows ; i++) {
      printf("%4d |",i+rowlwb-1);
      for( Int_t j=sheet_counter ; j<sheet_counter+cols_per_sheet && j<=ncols ; j++) printf(format,m(i+rowlwb-1,j+collwb-1));
      printf("\n");
    }
  }
  printf("\n");
  return 0;
} // end PrintMatrix


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
