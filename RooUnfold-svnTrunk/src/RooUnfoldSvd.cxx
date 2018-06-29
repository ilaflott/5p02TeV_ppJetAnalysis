//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldSvd.cxx 346 2013-07-30 13:37:06Z T.J.Adye $
//
// Description:
//      SVD unfolding. Just an interface to TSVDUnfold.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>Links to TSVDUnfold class which unfolds using Singular Value Decomposition (SVD).</p>
<p>Regularisation parameter defines the level at which values are deemed to be due to statistical fluctuations and are cut out. (Default= number of bins/2)
<p>Returns errors as a full matrix of covariances
<p>Can only handle 1 dimensional distributions
<p>Can account for both smearing and biasing
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldSvd.h"

#include <iostream>
#include <iomanip>
//#include <string>
#include <cassert>
#include "TString.h"
#include "TClass.h"
#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#if defined(HAVE_TSVDUNFOLD) || ROOT_VERSION_CODE < ROOT_VERSION(5,34,99)
#include "TSVDUnfold_local.h"  /* Use local copy of TSVDUnfold.h */
#else
#include "TSVDUnfold.h"
#if ROOT_VERSION_CODE < ROOT_VERSION(5,34,0)
#define TSVDUNFOLD_LEAK 1
#endif
#endif

#include "RooUnfoldResponse.h"

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldSvd);

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldSvd& rhs)
  : RooUnfold (rhs)
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  // The regularisation parameter is kreg.
  Init();
}

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg, Int_t ntoyssvd,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  // Constructor with old ntoyssvd argument. No longer required.
  Init();
  _NToys = ntoyssvd;
}

RooUnfoldSvd*
RooUnfoldSvd::Clone (const char* newname) const
{
  RooUnfoldSvd* unfold= new RooUnfoldSvd(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void
RooUnfoldSvd::Reset()
{
  Destroy();
  Init();
  RooUnfold::Reset();
}

void
RooUnfoldSvd::Destroy()
{
  delete _svd;
  delete _meas1d;
#ifdef TSVDUNFOLD_LEAK
  delete _meascov;
#endif
  delete _train1d;
  delete _truth1d;
  delete _reshist;
}

void
RooUnfoldSvd::Init()
{
  _svd= 0;
  _meas1d= _train1d= _truth1d= 0;
  _reshist= _meascov= 0;
  GetSettings();
}

void
RooUnfoldSvd::Assign (const RooUnfoldSvd& rhs)
{
  RooUnfold::Assign (rhs);
  CopyData (rhs);
}

void
RooUnfoldSvd::CopyData (const RooUnfoldSvd& rhs)
{
  _kreg= rhs._kreg;
}

TSVDUnfold*
RooUnfoldSvd::Impl()
{
  return _svd;
}

void
RooUnfoldSvd::Unfold()
{

  // Subtract fakes from measured distribution if they're there...
  cout<<endl;
  bool fakesExist=(bool)_res->FakeEntries();
  if(fakesExist)cout<<"//---------------------- \"fakes\" in response matrix detected!!! --------------------------- //"<<endl;
  else cout<<"no fakes in response matrix."<<endl;
  
  cout<<endl;
  bool subtractFakes=true&&fakesExist; 
  if(subtractFakes)cout<<"Will subtract fakes from measured data distribution."<<endl;
  else if (fakesExist&&!subtractFakes)cout<<"WARNING: fakes exist, but not subtracting fakes."<<endl;
  
  cout<<endl;
  bool resetBinErrs=true&&subtractFakes;
  if(resetBinErrs)cout<<"resetting measured data's errors"<<endl;
  if(subtractFakes&&!resetBinErrs)cout<<"WARNING: subtracting fakes but not resetting bin errors in measured data. Error bars may not accurately reflect real error!"<<endl;
  
  




  if (_res->GetDimensionTruth() != 1 || _res->GetDimensionMeasured() != 1) {
    cerr << "RooUnfoldSvd may not work very well for multi-dimensional distributions" << endl;
  }
  if (_kreg < 0) {
    cerr << "RooUnfoldSvd invalid kreg: " << _kreg << endl;
    return;
  }
  
  _nb= _nm > _nt ? _nm : _nt;
  
  if (_kreg > _nb) {
    cerr << "RooUnfoldSvd invalid kreg=" << _kreg << " with " << _nb << " bins" << endl;
    return;
  }


  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _meas1d=  HistNoOverflow (_meas,             _overflow); // data measured
  _train1d= HistNoOverflow (_res->Hmeasured(), _overflow); // MC meas.
  _truth1d= HistNoOverflow (_res->Htruth(),    _overflow); // MC truth
  _reshist= _res->HresponseNoOverflow() ; // Proj(MC Response, xaxis) = MC meas - fakes



  
  //TH1D* hfakes=(TH1D*)_res->Hfakes();    //for the error-resetting of the histo.
  TH1D* hfakes=(TH1D*) HistNoOverflow(_res->Hfakes(), _overflow);   //for the error-resetting of the histo.
  
  
  Resize ( _meas1d, _nb );
  Resize (_train1d, _nb );
  Resize (_truth1d, _nb );
  Resize (_reshist, _nb, _nb);
  Resize (  hfakes, _nb );
  
  //cout<<endl;
  


  
  //if ( (_res->FakeEntries()) && true ) {
  if ( fakesExist && subtractFakes) {
    
    
    Double_t fac_array[_nm]={0.}; //fac as a funcn of pt bin 
    Double_t fac_err_array[_nm]={0.}; //fac as a funcn of pt bin 
    
    Double_t fac=0.; //= data meas. norm/MC meas. norm, s.t. # fakes subtracted is proportional (what if i just leave the fakes in?)
    Double_t fac_err=0.;
    
    Int_t sm=hfakes->GetSumw2N(); //>1 if there's sum sq weights. 0 otherwise
    bool errsExist= (bool)sm;
    
    bool fill_fac_err_array=( ((bool)_meas1d->GetSumw2N()) && ((bool)_train1d->GetSumw2N()) );
    
    if(_verbose>=1){
      cout<<endl;
      cout << "# of meas bins=_nm="<<_nm<<endl;
      cout << "# of fake bins="<<hfakes->GetNbinsX()<<endl;        
      
      cout<<endl;
      cout <<"computing fac=data meas. norm/ MC meas. norm."<<endl;
      if(_verbose>=2)cout <<"computing fac_array(s). fac_array[ptbin -1]=data meas.[ptbin]/MC meas.[ptbin]"<<endl;
    }
    
    for (int i = 0; i< _nm; i++){
      
      if(_verbose>=2){
	cout<<endl;
	//cout <<"_reshist["<< i+1 << "] =" << _reshist->GetBinContent(i+1) <<endl;//DOESN'T WORK FO RSOME WEIRD REASON
	cout <<" _meas1d["<< i+1 << "] =" <<  _meas1d->GetBinContent(i+1) <<endl;
	cout <<"_train1d["<< i+1 << "] =" << _train1d->GetBinContent(i+1) <<endl;      }
      
      //make sure we don't divide by zero
      if(_train1d->GetBinContent(i+1) > 0.){
	
	fac_array[i]=( _meas1d->GetBinContent(i+1) / _train1d->GetBinContent(i+1) );
	
	if(fill_fac_err_array){
	  
	  if(_meas1d->GetBinContent(i+1)>0.){
	    Double_t meas_err_i= (_meas1d->GetBinError(i+1)/_meas1d->GetBinContent(i+1));
	    meas_err_i*=meas_err_i;
	    Double_t train_err_i=(_train1d->GetBinError(i+1)/_train1d->GetBinContent(i+1));
	    train_err_i*=train_err_i;
	    fac_err_array[i] = fac_array[i]*sqrt( train_err_i + meas_err_i );
	  }//end bin content check of _meas1d
	  else {
	    cout<<"WARNING bin i="<<i<<" for either meas or train 1d has no bin contents. fac_err_array entry not computed (div 0 issue). "<<endl;
	  }
	  
	  if(_verbose>=2) cout<<"fac_array["<<i<<"] = " << fac_array[i] << "+/- "<< fac_err_array[i] << endl;
	  
	  
	}//end fill err array
	
      }//end bin content check of _train1d 
      
      
    }//end fill fac_array
    
    
    Double_t  _hfakes_integ_err = 0.;
    Double_t  _meas1d_integ_err = 0.;
    Double_t _train1d_integ_err = 0.;
    Double_t _truth1d_integ_err = 0.;
    Double_t _reshist_integ_err = 0.;
    
    Double_t  _hfakes_integ =  hfakes->TH1::IntegralAndError(1, _nb,  _hfakes_integ_err , "");
    Double_t  _meas1d_integ = _meas1d->TH1::IntegralAndError(1, _nb,  _meas1d_integ_err , ""); //data
    Double_t _train1d_integ =_train1d->TH1::IntegralAndError(1, _nb, _train1d_integ_err , ""); //mc meas
    Double_t _truth1d_integ =_truth1d->TH1::IntegralAndError(1, _nb, _truth1d_integ_err , "");
    Double_t _reshist_integ =_reshist->TH1::IntegralAndError(1, _nb, _reshist_integ_err , "");
    
    if(_verbose>=1){
      cout<<endl;
      cout<<"  _hfakes_integ = " <<  _hfakes_integ << " +/- " <<  _hfakes_integ_err << endl;
      cout<<"  _meas1d_integ = " <<  _meas1d_integ << " +/- " <<  _meas1d_integ_err << endl;
      cout<<" _train1d_integ = " << _train1d_integ << " +/- " << _train1d_integ_err << endl;
      cout<<" _truth1d_integ = " << _truth1d_integ << " +/- " << _truth1d_integ_err << endl;
      cout<<" _reshist_integ = " << _reshist_integ << " +/- " << _reshist_integ_err << endl;
    }
    
    //assert(false);
    //TVectorD fakes= _res->Vfakes();             // diff bten. meas MC spectra + response proj onto x/meas axis
    //TVectorD _res_meas1d=_res->Vmeasured();
    //TVectorD _meas1d_vect=Vmeasured();    
    //Double_t fac_meas1d= Vmeasured().Sum(); //i think this is the meas data spectra; how can i test?
    //Double_t fac_resp= _res->Vmeasured().Sum(); //meas MC spectra i fed RooUnfold
    //_res->Vmeasured == vector version of (_res->Hmeasured)
    
    
    
    
    
    
    
    //cout<<"fakes.Sum()="<<fakes.Sum()<<endl;        
    //Double_t fakesHistSumTest=0.;      
    //for (int i = 0; i< _nm; i++)
    //  fakesHistSumTest+=hfakes->GetBinContent(i+1);
    //cout<<"fakesHistSumTest="<<fakesHistSumTest<<endl;    
    //trying tout fac as a func of pt bin
    //for (int i = 0; i< _nm; i++){

    
    
    
    //overall normalization factor between data and MC
    if (_train1d_integ>0.0) {
      
      //fac= fakes.Sum()/_res_meas1d.Sum();
      fac=  _meas1d_integ / _train1d_integ;
      if(errsExist){
	fac_err= sqrt( ((_meas1d_integ_err*_meas1d_integ_err)/(_meas1d_integ*_meas1d_integ)) 
		       + ((_train1d_integ_err*_train1d_integ_err)/(_train1d_integ*_train1d_integ)) );
	fac_err*=fac;
      }
      
    }
    else {
      cout<<"WARNING _train1d_integ is NOT greater than 0... check the file and/or code!"<<endl;
    }
    
    //if(true) fac=0.; //don't subtract fakes even if the fakes are there    
    //    cout << "from meas data hist, subtracting fac * fakes.Sum = " 
    //	 << fac << " * " << fakes.Sum() << " = "
    //	 << fac*fakes.Sum() << endl;
    
    if(_verbose>=1){ 
      cout<<endl;    
      cout<<" fac = _meas1d_integ/_train1d_integ = "<<fac <<" +/- "<<fac_err<<endl;
    }
    
    cout<<endl;
    cout << "from meas data hist, subtracting fac * _hfakes_integ = " 
      	 << fac << " * " << _hfakes_integ << " = "
      	 << fac*(_hfakes_integ) << endl;
    //cout<<endl;
    
    for (Int_t i= 1; i<=_nm; i++){
      
      if(_verbose>=2){
	cout <<"before fake subtraction..."<<std::endl;	
	cout <<"_meas1d->GetBinContent("<<i<<")="<<_meas1d->GetBinContent(i)<< " +/- "<< _meas1d->GetBinError(i)<<endl;
	cout <<"hfakes->GetBinContent("<<i<<")="<<hfakes->GetBinContent(i)<< " +/- "<< hfakes->GetBinError(i)<<endl;
      }
      
      

      //_meas1d->SetBinContent (i, _meas1d->GetBinContent(i)-(fac*fakes[i-1])); //ORIGINAL LINE //is okay 
      
      TString fakesubtype="subtract fakes[ptbin] * (data meas integ/MC meas integ)";       //is okay
      _meas1d->SetBinContent (i, _meas1d->GetBinContent(i)-(fac*hfakes->GetBinContent(i)));
      
      //TString fakesubtype="uniform sub across pt; subtract _meas1d[ptbin] * total # fakes/total # MC [no error reset]"; // is bad
      //_meas1d->SetBinContent (i, _meas1d->GetBinContent(i)-(fac*_meas1d->GetBinContent(i)) );
      
      //TString fakesubtype="not subtracting fakes"; //is bad
      //fac=0.;
      //_meas1d->SetBinContent (i, _meas1d->GetBinContent(i)-(fac*fakes[i-1]));
      
      //TString fakesubtype="uniform sub across pt; subtract fakes * (1.)";   // is too good (it's a lie :-()
      //fac=1.;
      //_meas1d->SetBinContent (i, _meas1d->GetBinContent(i)-(fac*fakes[i-1]));
      
      //TString fakesubtype="non-uniform sub across pt; subtract fakes[ptbin] * data norm[ptbin]/mc norm[ptbin]";      // is okay
      //_meas1d->SetBinContent (i, _meas1d->GetBinContent(i)-(fac_array[i-1]*hfakes->GetBinContent(i)));
      
      
      
      if(errsExist && resetBinErrs){
	
      	//Double_t fakesErr = fac*hfakes->GetBinError(i)*fac*hfakes->GetBinError(i);
	Double_t fakesErr = fac_array[i-1]*hfakes->GetBinError(i)*fac_array[i-1]*hfakes->GetBinError(i);
      	Double_t _meas1dErr= _meas1d->GetBinError(i)*_meas1d->GetBinError(i);
      	Double_t _meas1dErrNew= _meas1dErr+fakesErr;
      	
	cout<<fakesubtype<<endl;
	cout<<endl;
      	cout<<"fakesErr     ="<<sqrt(fakesErr)<<endl;
      	cout<<"_meas1dErr   ="<<sqrt(_meas1dErr)<<endl;
      	cout<<"_meas1dErrNew="<<sqrt(_meas1dErrNew)<<endl;
      	cout<<endl;
      	
      	//_meas1d->SetBinError(  i, sqrt( _meas1dErrNew )  );
	
	cout<<endl;
	cout<<"final err = "<<_meas1d->GetBinError(i)<<endl;
	cout<<endl;
	
      }//end resetting bin errs in data
      
      if(_verbose>=2){
	cout<<"fakesubtype="<<endl<<fakesubtype<<endl<<endl;
	cout <<"after fake subtraction + re-setting bin content..."<<std::endl;
	cout <<"_meas1d->GetBinContent("<<i<<")="<<_meas1d->GetBinContent(i)<<endl;
	cout <<"  _meas1d->GetBinError("<<i<<")="<<_meas1d->GetBinError(i)<<endl<<endl;    	      }
      
    }// end loop over bins
  }//end if fake entries
  
  //assert(false);
  _meascov= new TH2D ("meascov", "meascov", _nb, 0.0, 1.0, _nb, 0.0, 1.0);
  const TMatrixD& cov= GetMeasuredCov();
  
  for (Int_t i= 0; i<_nm; i++)
    for (Int_t j= 0; j<_nm; j++)
      _meascov->SetBinContent (i+1, j+1, cov(i,j));
  
  if (_verbose>=1) cout << "SVD init " << _reshist->GetNbinsX() << " x " << _reshist->GetNbinsY()
                        << " bins, kreg=" << _kreg << endl;
  
  //_meas1d --> data meas/RECO
  //_train 1d --> MC meas/RECO 
  //_truth 1d --> MC truth/GEN
  //_meascov --> ??? not exactly sure
  //_reshist=_res->HresponseNoOverflow() --> ??? not exactly sure, response matrix projection onto x axis?
  _svd= new TSVDUnfold (_meas1d, _meascov, _train1d, _truth1d, _reshist);

  TH1D* rechist= _svd->Unfold (_kreg);                 
  
  _rec.ResizeTo (_nt);
  for (Int_t i= 0; i<_nt; i++) {
    _rec[i]= rechist->GetBinContent(i+1);
  }

  if (_verbose>=2) {
    PrintTable (cout, _truth1d, _train1d, 0, _meas1d, rechist, _nb, _nb, kFALSE, kErrors);
    TMatrixD* resmat= RooUnfoldResponse::H2M (_reshist, _nb, _nb);
    RooUnfoldResponse::PrintMatrix(*resmat,"TSVDUnfold response matrix");
    delete resmat;
  }

  delete rechist;
  TH1::AddDirectory (oldstat);

  _unfolded= true;
  _haveCov=  false;
}

void
RooUnfoldSvd::GetCov()
{
  if (!_svd) return;
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);

  TH2D *unfoldedCov= 0, *adetCov= 0;
  //Get the covariance matrix for statistical uncertainties on the measured distribution
  if (_dosys!=2) unfoldedCov= _svd->GetXtau();
  //Get the covariance matrix for statistical uncertainties on the response matrix
  //Uses Poisson or Gaussian-distributed toys, depending on response matrix histogram's Sumw2 setting.
  if (_dosys)        adetCov= _svd->GetAdetCovMatrix (_NToys);

  _cov.ResizeTo (_nt, _nt);
  for (Int_t i= 0; i<_nt; i++) {
    for (Int_t j= 0; j<_nt; j++) {
      Double_t v  = 0;
      if (unfoldedCov) v  = unfoldedCov->GetBinContent(i+1,j+1);
      if (adetCov)     v += adetCov    ->GetBinContent(i+1,j+1);
      _cov(i,j)= v;
    }
  }

  delete adetCov;
#ifdef TSVDUNFOLD_LEAK
  delete unfoldedCov;
#endif
  TH1::AddDirectory (oldstat);

  _haveCov= true;
}

void RooUnfoldSvd::GetWgt()
{
  // Get weight matrix
  if (_dosys) RooUnfold::GetWgt();   // can't add sys errors to weight, so calculate weight from covariance
  if (!_svd) return;
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);

  //Get the covariance matrix for statistical uncertainties on the measured distribution
  TH2D* unfoldedWgt= _svd->GetXinv();

  _wgt.ResizeTo (_nt, _nt);
  for (Int_t i= 0; i<_nt; i++) {
    for (Int_t j= 0; j<_nt; j++) {
      _wgt(i,j)= unfoldedWgt->GetBinContent(i+1,j+1);
    }
  }

#ifdef TSVDUNFOLD_LEAK
  delete unfoldedWgt;
#endif
  TH1::AddDirectory (oldstat);

  _haveWgt= true;
}

void RooUnfoldSvd::GetSettings(){
    _minparm=0;
    _maxparm= _meas ? _meas->GetNbinsX() : 0;
    _stepsizeparm=1;
    _defaultparm=_maxparm/2;
}

void RooUnfoldSvd::Streamer (TBuffer &R__b)
{
  // Stream an object of class RooUnfoldSvd.
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldSvd::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldSvd::Class()->WriteBuffer (R__b, this);
  }
}
