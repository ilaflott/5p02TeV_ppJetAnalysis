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
#include <cassert>
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

// custom
#include "unfoldSpectra_simpbins.h"
#include "unfoldSpectra_anabins.h"



// -----------------------------------------------------------------------------------------------------------------------
const std::string CMSSW_BASE=
  "/net/hisrv0001/home/ilaflott/5p02TeV_ppJetAnalysis/CMSSW_7_5_8/src/readForests/outputCondor/";
const std::string SCRATCH_BASE=
  "/export/d00/scratch/ilaflott/5p02TeV_ppJetAnalysis/readForests/";
//const std::string unfoldSpectra_outdir="output/";
const std::string unfoldDataSpectra_outdir="output/unfoldDataSpectra/";
const std::string unfoldMCSpectra_outdir="output/unfoldMCSpectra/";





//Bayes setting that don't change too much
const bool doBayes=true; 
const int kIter = 4; //,kIterRange=4, kIterDraw = 3, kIterCenter=21;

//SVD setting that don't change too much
const bool doSVD=true; //!(doBayes); 
const int nKregMax  = 9 , kRegRange=(nKregMax-1)/2 ;//max num of diff kregs to do

//other options
const bool doOverUnderflows=false;
const bool clearOverUnderflows=false;

const bool normalizedMCMatrix=false;
const bool fillRespHists=false;
const bool useSimplePtBinning=false;//bin by ten everywhere instead of custom binning

//useful strings, numbers
const double integratedLuminosity=27.4*pow(10.,9.);//+/-2.4%
const std::string MCdesc= "QCD PY8 Tune CUETP8M1"; 
const std::string Datadesc1= "pp promptReco, #sqrt{s}=5.02 TeV"; 
const std::string Datadesc2= "L_{int}=27.4 pb^{-1} +/- 2.4%"; 



// color/styles
const Int_t color[13]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,kRed+2,kMagenta+1,kRed+1,kGreen+3};
// fill/fill-border color/style
const Int_t fcolor[5]={kRed-10,kBlue-10,kGreen-8,kOrange-4, kRed-6};
const Int_t lcolor[5]={kRed+1,kBlue-3,kGreen+3,kOrange+3, kRed+2};
// marker color/style
const Int_t fmstyle[6] = {20,21,22,23,29,3};
const Int_t emstyle[6] = {24,25,26,27,30,28};





void drawText(const char *text, float xp, float yp, int size){
  std::cout<<std::endl<<"in drawText"<<std::endl<<std::endl;

  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();

  std::cout<<std::endl<<"drawText done"<<std::endl<<std::endl;
  return;
}



TMatrixD* CalculatePearsonCoefficients(TMatrixD* covmat, bool debugPearson=false) {
  std::cout<<std::endl<<"in CalculatePearsonCoefficients"<<std::endl<<std::endl;

  int startRowColInt=0, skipRowColInt=1;   
  int NanCount=0, gt1Count=0, EntryCount=0;
  int nrows = covmat->GetNrows(), ncols = covmat->GetNcols();

  TMatrixD* pearsonCoefs= (TMatrixD*)covmat->Clone("pearsonCoefs");

  if(debugPearson)std::cout<<"looping over (nrows,ncols)=("<<nrows<<", "<<ncols<<")"<<std::endl;
  for(int row = startRowColInt; row<nrows; row++){ double covmatElement_rowrow = (*covmat)(row,row);
    for(int col = startRowColInt; col<ncols; col++){ double covmatElement_colcol = (*covmat)(col,col);
      double covmatElement_rowcol = (*covmat)(row,col);
      double pearson = covmatElement_rowcol/TMath::Sqrt(covmatElement_rowrow*covmatElement_colcol);

      // NaN protection
      bool notNan=true;
      if(pearson!=pearson){ 	pearson=-10.;
	if(debugPearson&&(row%skipRowColInt==0||col%skipRowColInt==0)) { 
	  //std::cout<<"for (row,col)= "<<row <<" , "<<col <<std::endl;
	  std::cout<<"NaN! pearson("<<row<<","<<col<<")=-10."<<std::endl; 	}
	NanCount++;	notNan=false;    gt1Count++;  }
      
      // assign the value to the TMatrix
      (*pearsonCoefs)(row,col) = pearson;            

      // debug
      
      if(notNan&&debugPearson&&abs(pearson)>1.) {
	std::cout<<"warning, abs(pearson)="<<abs(pearson)<<std::endl; 
	gt1Count++;      } 
      if( debugPearson&&(row%skipRowColInt==0&&col%skipRowColInt==0&&notNan) ) {
	//std::cout<<"for (row,col)= "<<row <<" , "<<col <<std::endl;	
	//std::cout<<"covmat(row,row)="<< covmatElement_rowrow<<std::endl;       
	//std::cout<<"covmat(col,col)="<< covmatElement_colcol<<std::endl;       
	//std::cout<<"covmat(row,col)="<< covmatElement_rowcol<<std::endl;       
	//std::cout<<"pearson="<<pearson <<std::endl;       
	std::cout<<"pearsonCoefs("<<row<<","<<col<<")="<<( TMatrixD (*pearsonCoefs) )(row,col) <<std::endl;       }
      EntryCount++;      
    } // end col loop
  } // end row loop

  if(debugPearson)std::cout<<"Entries in Pearson Matrix="<<EntryCount<<std::endl;  
  if(debugPearson)std::cout<<"pearson>1 Entries="<<gt1Count<<std::endl;  
  if(debugPearson)std::cout<<"Nan Entries="<<NanCount<<std::endl;  

  std::cout<<std::endl<<"calculatePearsonCoefficients done"<<std::endl<<std::endl;
  return pearsonCoefs;
}


//-----------------------------------------------------------------------------------------------------------------------


void divideBinWidth(TH1 *h){
  std::cout<<std::endl<<"in divideBinWidth"<<std::endl<<std::endl;

  h->Sumw2();
  double nEntries=h->GetEntries();
  int numbins=h->GetNbinsX();

  //std::cout<<"#of bins is "<<numbins<<std::endl;
  //std::cout<<"according to consts array... "<<nbins_pt<<std::endl;
  std::cout<<"integral before:  "<<h->Integral()<<std::endl;
  //std::cout<<"entries before:  "<<nEntries<<std::endl;
  
  for (int i=1;i<=numbins;i++)
    {
      
      //std::cout<<"before.. i = numbins = "<< numbins<<std::endl;
      //std::cout<<"bin width = "<<   h->GetBinWidth(i)   << std::endl;
      //std::cout<<"bin center = "<<  h->GetBinCenter(i)  << std::endl;
      //std::cout<<"val = "   <<      h->GetBinContent(i) << std::endl;
      //std::cout<<"valErr = "<<      h->GetBinError(i)   << std::endl;
      
      
      Float_t val    = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);
      
      val    /= h->GetBinWidth(i);
      valErr /= h->GetBinWidth(i);
      
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
      
      //std::cout<<"after... i = numbins = "<< numbins<<std::endl;
      //std::cout<<"bin width = "<<   h->GetBinWidth(i)   << std::endl;
      //std::cout<<"bin center = "<<  h->GetBinCenter(i)  << std::endl;
      //std::cout<<"val = "   <<      h->GetBinContent(i) << std::endl;
      //std::cout<<"valErr = "<<      h->GetBinError(i)   << std::endl;
    }

  std::cout<<"integral after:  "<<h->Integral()<<std::endl;
  //std::cout<<"entries after:  "<<h->GetEntries()<<std::endl;

  h->SetEntries(nEntries);

  //std::cout<<"exiting function"<<std::endl;
  std::cout<<std::endl<<"divideBinWidth done"<<std::endl<<std::endl;
  return;
}



void divideBinWidth_TH2(TH2F *h){
  std::cout<<std::endl<<"in divideBinWidth_TH2"<<std::endl<<std::endl;

  h->Sumw2();
  h->Print("base");      std::cout<<std::endl;

  double numEntries=h->GetEntries();

  for (int i=1;i<=h->GetNbinsX();i++)    {    
    for (int j=1;j<=h->GetNbinsY();j++)    {
      
      Float_t xWidth=h->GetXaxis()->GetBinWidth(i);
      Float_t yWidth=h->GetYaxis()->GetBinWidth(j);
      
      Float_t val = (h->GetBinContent(i,j))/(xWidth*yWidth);
      Float_t valErr = (h->GetBinError(i,j))/(xWidth*yWidth);
      
      h->SetBinContent(i,j,val);
      h->SetBinError(i,j,valErr);
    }
  }

  h->SetEntries(numEntries);
  h->Print("base");      std::cout<<std::endl;
  
  std::cout<<std::endl<<"divideBinWidth_TH2 done"<<std::endl<<std::endl;
  return;
}

//normalizes MC response TH2 to the rowSum of the gen-axis
void normalizeMC_TH2(TH2F* inputTH2){
  std::cout<<std::endl<<"in normalizeMC_TH2"<<std::endl<<std::endl;
  
  inputTH2->Sumw2();
  inputTH2->Print("base");  std::cout<<std::endl;

  double numEntries=inputTH2->GetEntries();
  
  int colStart=1;
  TAxis *xaxis = inputTH2->GetXaxis();
  int nbins_x= xaxis->GetNbins();
  
  int rowStart=1;
  TAxis *yaxis = inputTH2->GetYaxis();
  int nbins_y = yaxis->GetNbins(); 
  
  if(doOverUnderflows){
    nbins_x++;    nbins_y++;
    rowStart--;   colStart--; }
  
  const int rowSums_len=nbins_x;
  float rowSums[rowSums_len]={0.};
  float rowSumErrs[rowSums_len]={0.};
  
  //calculate sums
  std::cout<<"computing row sums..."<<std::endl;
  for(int rowNum=rowStart;rowNum<=nbins_y;rowNum++){
    float theSum=0.; 
    float theSumErr=0.; 
    for(int colNum=colStart;colNum<=nbins_x;colNum++){
      float theVal=inputTH2->GetBinContent(colNum,rowNum);
      float theValErr=inputTH2->GetBinError(colNum,rowNum);
      theSum+=theVal;
      theSumErr+=(theValErr*theValErr);
    }
    rowSums[rowNum]=theSum;
    rowSumErrs[rowNum]=std::sqrt(theSumErr);
  }
  std::cout<<std::endl;

  //scale each entry in each row down by that row's sum
  std::cout<<"scaling rows..."<<std::endl;
  for(int  rowNum=rowStart; rowNum<=nbins_y; rowNum++){

    float theSum=rowSums[rowNum]; 
    float theSumErr=rowSumErrs[rowNum]; 

    for( int colNum=colStart; colNum<=nbins_x; colNum++){

      float theVal=(inputTH2->GetBinContent(colNum,rowNum));
      float theValErr=(inputTH2->GetBinError(colNum,rowNum));

      float finalVal=theVal/theSum;
      float finalErr=(theVal/theSum)* std::sqrt( 
						(theSumErr/theSum)*(theSumErr/theSum) + 
						(theValErr/theVal)*(theValErr/theVal) )  ;


      inputTH2->SetBinContent(colNum,rowNum,finalVal);
      inputTH2->SetBinError(colNum,rowNum,finalErr);      
    }
  }
  std::cout<<std::endl;

  inputTH2->SetEntries(numEntries);
  inputTH2->Print("base");

  std::cout<<std::endl<<"exiting normalizeMC_TH2"<<std::endl<<std::endl;
  return ;
  

}


//TH2F* reBinTH2(TH2F* inputTH2, std::string inputTH2_title, const double* boundaries_pt, const int nbins_pt){
TH2F* reBinTH2(TH2F* inputTH2, std::string inputTH2_title, 
	       double* boundaries_pt_reco, int nbins_pt_reco,
	       double* boundaries_pt_gen  , int nbins_pt_gen	       ){

  std::cout<<std::endl<<"in reBinTH2"<<std::endl<<std::endl;


  double numEntries=inputTH2->GetEntries();
  inputTH2->Sumw2();
  inputTH2->Print("base");  std::cout<<std::endl;
  
  
  int colStart=1;
  TAxis *xaxis = inputTH2->GetXaxis(); //reco pt axis
  int nbins_x= xaxis->GetNbins();

  int rowStart=1;
  TAxis *yaxis = inputTH2->GetYaxis(); //gen pt axis
  int nbins_y = yaxis->GetNbins();

  int finalxbin=nbins_x, finalybin=nbins_y;
  if(doOverUnderflows){
    finalxbin=nbins_x+1;    finalybin=nbins_y+1;
    rowStart--;   colStart--; }

  TH2F *reBinnedTH2 = new TH2F(inputTH2_title.c_str(), inputTH2->GetTitle(), 
			       nbins_pt_reco, boundaries_pt_reco , 
			       nbins_pt_gen, boundaries_pt_gen );
  reBinnedTH2->Sumw2();
  
  for (  int rowNum=rowStart ; rowNum <= (finalybin) ; rowNum++ ) {
    for (int colNum=colStart ; colNum <= (finalxbin) ; colNum++ ) {
      
      //std::cout<<"i="<<i<<" and j="<<j<<std::endl; 
      //std::cout<<"bin x center = "<<xaxis->GetBinCenter(i)<< std::endl ;
      //std::cout<<"bin y center = "<<yaxis->GetBinCenter(j)<< std::endl ;
      
      //if(j==0 || i==0 || j==nbins_y || i==nbins_x){	
      //	if( inputTH2->GetBinContent(i,j) != 0 )	  {
      //	  //std::cout<<"i="<<i<<" and j="<<j<<std::endl; 
      //	  //std::cout<<"under/overflow bin in TH2 nonzero!!"<<std::endl;
      //	  //std::cout<<"under/overflow bin content="<<inputTH2->GetBinContent(i,j)<<std::endl;
      //	  continue;	}      }      
      //else 
      //if(rowNum!=0 && colNum !=0 && rowNum != (nbinsy+1) && colNum !=(ybins+1) )
      reBinnedTH2->Fill( xaxis->GetBinCenter(colNum) , yaxis->GetBinCenter(rowNum) ,
			 inputTH2->GetBinContent(colNum,rowNum) );            
      //else if(doOverUnderFlows)	{
      //
      //reBinnedTH1->Fill
      
      
    }  }//end x-y loop
  
  reBinnedTH2->SetEntries(numEntries);
  reBinnedTH2->Print("base");  std::cout<<std::endl;
  
  //std::cout<<" done rebinning TH2" <<std::endl<<std::endl;    
  std::cout<<std::endl<<"reBinTH2 done"<<std::endl<<std::endl;  
  return reBinnedTH2;
}



TH2F* reBinPearsonTH2(TMatrixD* pearson, const double* boundaries_pt, const int nbins_pt){
  std::cout<<std::endl<<"in reBinPearsonTH2"<<std::endl<<std::endl;

  bool debugPearson=false;
  //const double* the_ptBins=boundaries_pt;//ptbins
  //const int numbins=nbins_pt;//nbins

  TH2F *reBinnedTH2 = new TH2F(  "pearsonCustomBins" , "pearsonMatrix custom bins"   , 
				 nbins_pt,boundaries_pt,
				 nbins_pt,boundaries_pt);
  std::cout<<std::endl<<"in reBinPearsonTH2 function..."<<std::endl;  
  reBinnedTH2->Print("base");
  
  for (int j=1;   j<=nbins_pt;j++) {
    for (int i=1; i<=nbins_pt;i++) {

      double pearson_i_j=(*pearson)(i-1,j-1);
      
      if(fabs(pearson_i_j) > 1.0) {
	if(debugPearson)std::cout <<std::endl<<"WARNING:: pearson["<<i-1<<"]["<<j-1<<"]="<<pearson_i_j<<std::endl<<std::endl;
	//	pearson_i_j=10.;
	//std::cout<< "pearson_i_j=10"<<std::endl
	continue; }
      else {
	if(debugPearson)std::cout <<"pearson["<<i-1<<"]["<<j-1<<"]="<<pearson_i_j<<std::endl;	}

      reBinnedTH2->Fill( i, j, fabs( (float) pearson_i_j ) );
    }
  }


  reBinnedTH2->Print("base");

  std::cout<<std::endl<<"reBinPearsonTH2 done"<<std::endl<<std::endl;
  return reBinnedTH2;
}


void TH1clearOverUnderflows(TH1* h)
{
  std::cout<<std::endl<<"WARNING!!!!"<<std::endl<<std::endl;
  std::cout<<"doOverUnderflows="<<doOverUnderflows<<std::endl;
  std::cout<<"clearing Over/Underflow Bins...."<<std::endl;
  h->TH1::ClearUnderflowAndOverflow();
  return;  
}











// PrintMatrix
//-----------------------------------------------------------------------------------------------------------------------
// Print the matrix as a table of elements.
// Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
// By default, format ="%11.4g", is used to print one element.
// One can specify an alternative format, e.g. format ="%6.2f  " 
int PrintMatrix( const TMatrixD& m, 
		  const char* format, const char* name, Int_t cols_per_sheet ){
  std::cout<<std::endl<<"in PrintMatrix"<<std::endl<<std::endl;
  
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
  std::cout<<std::endl<<"PrintMatrix done"<<std::endl<<std::endl;
  return 0;
} // end PrintMatrix



