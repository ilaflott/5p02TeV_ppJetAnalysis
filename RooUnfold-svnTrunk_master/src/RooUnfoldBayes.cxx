//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBayes.cxx 335 2012-06-08 18:18:04Z T.J.Adye $
//
// Description:
//      Bayesian unfolding.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>Links to the RooUnfoldBayesImpl class which uses Bayesian unfolding to reconstruct the truth distribution.</p>
<p>Works for 2 and 3 dimensional distributions
<p>Returned errors can be either as a diagonal matrix or as a full matrix of covariances
<p>Regularisation parameter sets the number of iterations used in the unfolding (default=4)
<p>Is able to account for bin migration and smearing
<p>Can unfold if test and measured distributions have different binning.
<p>Returns covariance matrices with conditions approximately that of the machine precision. This occasionally leads to very large chi squared values
END_HTML */

/////////////////////////////////////////////////////////////

//#define OLDERRS   // restore old (incorrect) error calculation
//#define OLDERRS2  // restore old (incorrect) systematic error calculation
//#define OLDMULT   // restore old (slower) matrix multiplications

#include "RooUnfoldBayes.h"

#include <iostream>
#include <iomanip>
#include <math.h>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"

#include "RooUnfoldResponse.h"

using std::min;
using std::cerr;
using std::endl;
using std::cout;
using std::setw;
using std::left;
using std::right;

ClassImp (RooUnfoldBayes);

RooUnfoldBayes::RooUnfoldBayes (const RooUnfoldBayes& rhs)
  : RooUnfold (rhs)
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

RooUnfoldBayes::RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas, Int_t niter, Bool_t smoothit,
                                const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _niter(niter), _smoothit(smoothit)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  // The regularisation parameter is niter (number of iterations).
  Init();
}

RooUnfoldBayes* RooUnfoldBayes::Clone (const char* newname) const
{
  // Creates a copy of the RooUnfoldBayes object
  RooUnfoldBayes* unfold= new RooUnfoldBayes(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void RooUnfoldBayes::Init()
{
  _nc= _ne= 0;
  _nbartrue= _N0C= 0.0;
  GetSettings();
}

void RooUnfoldBayes::Reset()
{
  Init();
  RooUnfold::Reset();
}

void RooUnfoldBayes::Assign (const RooUnfoldBayes& rhs)
{
  RooUnfold::Assign (rhs);
  CopyData (rhs);
}

void RooUnfoldBayes::CopyData (const RooUnfoldBayes& rhs)
{
  _niter=    rhs._niter;
  _smoothit= rhs._smoothit;
}

void RooUnfoldBayes::Unfold()
{
  setup();
  if (verbose() >= 2) {
    Print();
    RooUnfoldResponse::PrintMatrix(_Nji,"RooUnfoldBayes response matrix (Nji)");
  }
  if (verbose() >= 1) cout << "Now unfolding..." << endl;
  unfold();
  if (verbose() >= 2) Print();
  _rec.ResizeTo(_nc);
  _rec = _nbarCi;
  _rec.ResizeTo(_nt);  // drop fakes in final bin
  _unfolded= true;
  _haveCov=  false;
}

void RooUnfoldBayes::GetCov()
{
  getCovariance();
  _cov.ResizeTo (_nt, _nt);  // drop fakes in final bin
  _haveCov= true;
}

void RooUnfoldBayes::GetSettings()
{
    _minparm=1;
    _maxparm=15;
    _stepsizeparm=1;
    _defaultparm=4;
}

TMatrixD& RooUnfoldBayes::H2M (const TH2* h, TMatrixD& m, Bool_t overflow)
{
  // TH2 -> TMatrixD
  if (!h) return m;
  Int_t first= overflow ? 0 : 1;
  Int_t nm= m.GetNrows(), nt= m.GetNcols();
  for (Int_t j= 0; j < nt; j++)
    for (Int_t i= 0; i < nm; i++)
      m(i,j)= h->GetBinContent(i+first,j+first);
  return m;
}

//-------------------------------------------------------------------------
void RooUnfoldBayes::setup()
{
  _nc = _nt; //NUMBER OF GEN pT BINS //number of truth [_nt] bins (aka number of cause [_nc] bins)
  _ne = _nm; //NUMBER OF RECO pT BINS //number of measured [_nm] bins (aka number of effect [_ne] bins)

  _nEstj.ResizeTo(_ne); //number of measured events from effect E_j (aka effect/measured bin j)
  _nEstj= Vmeasured(); //RECO DATA HISTOGRAM//this is the measured events histogram/vector

  _nCi.ResizeTo(_nt);
  _nCi= _res->Vtruth();// GEN pT HISTOGRAM I GIVE TO ROOUNFOLD [NOT A PROJECTION]
                       

  _Nji.ResizeTo(_ne,_nt);
  H2M (_res->Hresponse(), _Nji, _overflow);   //TURN RESPONSE MATRIX INTO TMATRIX CALLED _Nji
  // don't normalise, which is what _res->Mresponse() would give us

  if (_res->FakeEntries()) {// IF THERE EXIST FAKES, DETERMINE THE FAKES VECTOR
    TVectorD fakes= _res->Vfakes();
    Double_t nfakes= fakes.Sum();
    if (verbose()>=0) cout << "Add truth bin for " << nfakes << " fakes" << endl;
    _nc++;
    _nCi.ResizeTo(_nc);
    _nCi[_nc-1]= nfakes;//truth bin for fakes added at end
    _Nji.ResizeTo(_ne,_nc);
    for (Int_t i= 0; i<_nm; i++) _Nji(i,_nc-1)= fakes[i];
  }

  _nbarCi.ResizeTo(_nc);//estimated number of true events from cause C_i [aka cause/truth bin i] ??
  _efficiencyCi.ResizeTo(_nc); //efficiency vector (accounts for misses/inefficiency)
  _Mij.ResizeTo(_nc,_ne); //unfolding matrix
  _P0C.ResizeTo(_nc);//prior before last iteration  // how does this get used?
  _UjInv.ResizeTo(_ne);// (1 / folded prior) from last iteration ??
#ifndef OLDERRS
  if (_dosys!=2) _dnCidnEj.ResizeTo(_nc,_ne);
#endif
  if (_dosys)    _dnCidPjk.ResizeTo(_nc,_ne*_nc);

  // Initial distribution //WHERE _P0C GETS SET INITIALLY
                          //_P0C IS _nCi initially, "number of true events from cause C_i"
  _N0C= _nCi.Sum(); //INTEGRAL OF GEN pT HISTOGRAM I GIVE ROO UNFOLD
  if (_N0C!=0.0) {
    _P0C= _nCi;
    _P0C *= 1.0/_N0C;//_P0C == GEN pT HISTOGRAM NORMALIZED TO 1
  }
}

//-------------------------------------------------------------------------
void RooUnfoldBayes::unfold()
{
  // Calculate the unfolding matrix.
  // _niter = number of iterations to perform (3 by default).
  // _smoothit = smooth the matrix in between iterations (default false).
  

  TMatrixD PEjCi(_ne,_nc);    //this is response matrix w/ all entries normalized by the gen pT hist gentries in their respective gen pT bin, each reco row doesn't sum to 1, it sums to the efficiency
  TMatrixD PEjCiEff(_ne,_nc); //this should be response matrix accounted for misses, i.e. each reco row should sum to 1
  for (Int_t i = 0 ; i < _nc ; i++) {
    if (_nCi[i] <= 0.0) { _efficiencyCi[i] = 0.0; continue; }
    Double_t eff = 0.0;
    for (Int_t j = 0 ; j < _ne ; j++) {
      Double_t response = _Nji(j,i) / _nCi[i];
      PEjCi(j,i) = PEjCiEff(j,i) = response;  //response = response matrix entry j i / truth spectra i 
      eff += response;//efficiency = N1i/nCi + N2i/nCi + ... N_nei/nCi
    }
    _efficiencyCi[i] = eff;
    Double_t effinv = eff > 0.0 ? 1.0/eff : 0.0;   // reset PEjCiEff if eff=0
                                                   // effinv = eff^-1 as long as it's not 0 (which would signify no truth events in that truth bin)
    for (Int_t j = 0 ; j < _ne ; j++) PEjCiEff(j,i) *= effinv;
  }

  //for the first iteration, PbarCi empty. it will take in the new numbers after the unfolding is applied to the initial dist/previous iteration's dist
  TVectorD PbarCi(_nc); //PbarCi created here with _nc bins.
  _chi2iter= new TH1D("_chi2_iter", "chi2 between iterations", _niter, 0.5, _niter+0.5);
  
  for (Int_t kiter = 0 ; kiter < _niter; kiter++) {
    
    if (verbose()>=1) cout << "Iteration : " << kiter << endl;
    
    // update prior, P0C, from previous iteration, PbarCi
    if (kiter>0) {
      _P0C = PbarCi;//for iterations past the first iteration
                    // for _P0C's first iteration, search 'WHERE _P0C GETS SET INITIALLY'
      _N0C = _nbartrue;
    }
    
    //calculated Uj, UjInv. Uj = resp matrix norm'd to eff * truth. UjInv=1/Uj for each entry in vector
    for (Int_t j = 0 ; j < _ne ; j++) {
      Double_t Uj = 0.0;
      for (Int_t i = 0 ; i < _nc ; i++)
        Uj += PEjCi(j,i) * _P0C[i];
      _UjInv[j] = Uj > 0.0 ? 1.0/Uj : 0.0;
    }

    // calculate Unfolding matrix Mij = resp matrix norm'd to eff * truth * resp matrix norm'd to 1 * truth
    _nbartrue = 0.0;
    for (Int_t i = 0 ; i < _nc ; i++) {
      Double_t nbarC = 0.0;
      for (Int_t j = 0 ; j < _ne ; j++) {
        Double_t Mij = _UjInv[j] * PEjCiEff(j,i) * _P0C[i];//compute response matrix element
        _Mij(i,j) = Mij;
        nbarC += Mij * _nEstj[j]; //apply response matrix to RECO data
      }
      _nbarCi[i] = nbarC; //unfolded RECO data value for bin i
      _nbartrue += nbarC;  // best estimate of true number of events
    }
    
    // new estimate of true distribution
    PbarCi= _nbarCi; 
    PbarCi *= 1.0/_nbartrue; //PbarCi == 1/ unfolded RECO data value for bin i

#ifndef OLDERRS
    if (_dosys!=2) {
      if (kiter <= 0) {
        _dnCidnEj= _Mij;
      } else {
#ifndef OLDMULT
        TVectorD en(_nc);//-1 * 1/gen Pt Hist corr for efficiency
	TVectorD nr(_nc);
        for (Int_t i = 0 ; i < _nc ; i++) {
          if (_P0C[i]<=0.0) continue;
          Double_t ni= 1.0/(_N0C*_P0C[i]);// should just be 1 / the gen Pt histogram value for gen pT bin i
	  // !!!! THIS IS WHERE THE EFFICIENCY/INEFFICIENCY/MISSES MAGIC HAPPENS !!!! //
	  cout<<"LOCATION 1: _efficiencyCi["<<i<<"]="<<_efficiencyCi[i]<<endl;
          en[i]= -ni*_efficiencyCi[i];
          nr[i]=  ni*_nbarCi[i];
        }
        TMatrixD M1= _dnCidnEj;
        M1.NormByColumn(nr,"M");
        TMatrixD M2 (TMatrixD::kTransposed, _Mij);
        M2.NormByColumn(_nEstj,"M");
        M2.NormByRow(en,"M");
        TMatrixD M3 (M2, TMatrixD::kMult, _dnCidnEj);
        _dnCidnEj.Mult (_Mij, M3);
        _dnCidnEj += _Mij;
        _dnCidnEj += M1;
#else /* OLDMULT */
        TVectorD ksum(_ne);
        for (Int_t j = 0 ; j < _ne ; j++) {
          for (Int_t k = 0 ; k < _ne ; k++) {
            Double_t sum = 0.0;
            for (Int_t l = 0 ; l < _nc ; l++) {
              if (_P0C[l]>0.0) {
		sum += _efficiencyCi[l]*_Mij(l,k)*_dnCidnEj(l,j)/_P0C[l];
		//cout<<"LOCATION 2: _efficiencyCi["<<i<<"]="<<_efficiencyCi[i]<<endl;
	      }
            }
            ksum[k]= sum;
          }
          for (Int_t i = 0 ; i < _nc ; i++) {
            Double_t dsum = _P0C[i]>0 ? _dnCidnEj(i,j)*_nbarCi[i]/_P0C[i] : 0.0;
            for (Int_t k = 0 ; k < _ne ; k++) {
              dsum -= _Mij(i,k)*_nEstj[k]*ksum[k];
            }
            // update dnCidnEj. Note that we can do this in-place due to the ordering of the accesses.
            _dnCidnEj(i,j) = _Mij(i,j) + dsum/_N0C;
          }
        }
#endif /*end OLDMULT */
      }
    }
#endif/* end OLDERRS */

    if (_dosys) {
#ifndef OLDERRS2
      if (kiter > 0) {
        TVectorD mbyu(_ne);
        for (Int_t j = 0 ; j < _ne ; j++) {
          mbyu[j]= _UjInv[j]*_nEstj[j]/_N0C;
        }
        TMatrixD A= _Mij;
        A.NormByRow (mbyu, "M");
        TMatrixD B(A, TMatrixD::kMult, PEjCi);
        TMatrixD dnCidPjkUpd (B, TMatrixD::kMult, _dnCidPjk);
        Int_t nec= _ne*_nc;
        for (Int_t i = 0 ; i < _nc ; i++) {
          if (_P0C[i]<=0.0) continue;  // skip loop: dnCidPjkUpd(i,jk) will also be 0 because _Mij(i,j) will be 0
          Double_t r= PbarCi[i]/_P0C[i];
          for (Int_t jk= 0; jk<nec; jk++)
            _dnCidPjk(i,jk)= r*_dnCidPjk(i,jk) - dnCidPjkUpd(i,jk);
        }
      }
#else  /* OLDERRS2 */
      if (kiter == _niter-1)   // used to only calculate _dnCidPjk for the final iteration
#endif /*end OLDERRS2 */
      for (Int_t j = 0 ; j < _ne ; j++) {
        if (_UjInv[j]==0.0) continue;
        Double_t mbyu= _UjInv[j]*_nEstj[j];
        Int_t j0= j*_nc;
        for (Int_t i = 0 ; i < _nc ; i++) {
          Double_t b= -mbyu * _Mij(i,j);
          for (Int_t k = 0 ; k < _nc ; k++) _dnCidPjk(i,j0+k) += b*_P0C[k];
          if (_efficiencyCi[i]!=0.0){	    
            _dnCidPjk(i,j0+i) += (_P0C[i]*mbyu - _nbarCi[i]) / _efficiencyCi[i];
	    //cout<<"LOCATION 3: _efficiencyCi["<<i<<"]="<<_efficiencyCi[i]<<endl;
	  }
        }
      }
    }

    // no need to smooth the last iteraction
    if (_smoothit && kiter < (_niter-1)) smooth(PbarCi);

    // Chi2 based on Poisson errors
    //    Double_t chi2 = getChi2(PbarCi, _P0C, _nbartrue);
    // PbarCi == latest iteration of unfolded spectra, normalized to #of events. 
    // _P0C == prev iteration, if kiter=0, this is simply the truth histogram from roounfold response object, normlized to 1.
    // _nbartrue == sum of #of events in PbarCi
    Double_t chi2 = getChi2(PbarCi, _P0C, _nbartrue, (TH1D*) _chi2iter, kiter);
    if (verbose()>=1) cout << "Chi^2 of change " << chi2 << endl;
    if(verbose()>=1)cout<< "sanity check in RooUnfoldBayes::unfold(); bin content of _chi2iter hist is "<<  _chi2iter->GetBinContent(kiter+1) << endl;
    
    // and repeat
  }//end kIter loop
}//end unfold function

//-------------------------------------------------------------------------
void RooUnfoldBayes::getCovariance()
{
  if (_dosys!=2) {
    if (verbose()>=1) cout << "Calculating covariances due to number of measured events" << endl;

    // Create the covariance matrix of result from that of the measured distribution
    _cov.ResizeTo (_nc, _nc);
#ifdef OLDERRS
    const TMatrixD& Dprop= _Mij;
#else
    const TMatrixD& Dprop= _dnCidnEj;
#endif
    if (_haveCovMes) {
      ABAT (Dprop, GetMeasuredCov(), _cov);
    } else {
      TVectorD v= Emeasured();
      v.Sqr();
      ABAT (Dprop, v, _cov);
    }
  }

  if (_dosys) {
    if (verbose()>=1) cout << "Calculating covariance due to unfolding matrix..." << endl;

    const TMatrixD& Eres= _res->Eresponse();
    TVectorD Vjk(_ne*_nc);           // vec(Var(j,k))
    for (Int_t j = 0 ; j < _ne ; j++) {
      Int_t j0= j*_nc;
      for (Int_t i = 0 ; i < _nc ; i++) {
        Double_t e= Eres(j,i);
        Vjk[j0+i]= e*e;
      }
    }

    if (_dosys!=2) {
      TMatrixD covres(_nc,_nc);
      ABAT (_dnCidPjk, Vjk, covres);
      _cov += covres;
    } else {
      _cov.ResizeTo (_nc, _nc);
      ABAT (_dnCidPjk, Vjk, _cov);
    }
  }
}

//-------------------------------------------------------------------------
void RooUnfoldBayes::smooth(TVectorD& PbarCi) const
{
  // Smooth unfolding distribution. PbarCi is the array of proababilities
  // to be smoothed PbarCi; nevts is the numbers of events
  // (needed to calculate suitable errors for the smearing).
  // PbarCi is returned with the smoothed distribution.

  if (_res->GetDimensionTruth() != 1) {
    cerr << "Smoothing only implemented for 1-D distributions" << endl;
    return;
  }
  if (verbose()>=1) cout << "Smoothing." << endl;
  TH1::SmoothArray (_nc, PbarCi.GetMatrixArray(), 1);
  return;
}

//-------------------------------------------------------------------------
Double_t RooUnfoldBayes::getChi2(const TVectorD& prob1,
                                 const TVectorD& prob2,
                                 Double_t nevents) const
{
  // calculate the chi^2. prob1 and prob2 are the probabilities
  // and nevents is the number of events used to calculate the probabilities
  Double_t chi2= 0.0;
  Int_t n= prob1.GetNrows();
  if (verbose()>=2) cout << "chi2 " << n << " " << nevents << endl;
  for (Int_t i = 0 ; i < n ; i++) {
    Double_t psum  = (prob1[i] + prob2[i])*nevents;
    Double_t pdiff = (prob1[i] - prob2[i])*nevents;
    if (psum > 1.0) {
      chi2 = chi2 + (pdiff*pdiff)/psum;
    } else {
      chi2 = chi2 + (pdiff*pdiff);
    }
  }
  return(chi2);
}



Double_t RooUnfoldBayes::getChi2(const TVectorD& prob1,
                                 const TVectorD& prob2,
                                 Double_t nevents, TH1D* h_chi2iter, int curr_kiter) const
{
  // calculate the chi^2. prob1 and prob2 are the probabilities
  // and nevents is the number of events used to calculate the probabilities
  Double_t chi2= 0.0;
  Int_t n= prob1.GetNrows();
  if (verbose()>=2) cout << "chi2 " << n << " " << nevents << endl;
  for (Int_t i = 0 ; i < n ; i++) {
    Double_t psum  = (prob1[i] + prob2[i])*nevents;
    Double_t pdiff = (prob1[i] - prob2[i])*nevents;
    if (psum > 1.0) {
      chi2 = chi2 + (pdiff*pdiff)/psum;
    } else {
      chi2 = chi2 + (pdiff*pdiff);
    }
  }
  h_chi2iter->SetBinContent(curr_kiter+1, chi2);
  if(verbose()>=2)cout<<"sanity check in new getChi2, bin content h_chi2iter is " 
		      << h_chi2iter->GetBinContent(curr_kiter+1)<<endl;
  return(chi2);
}

//-------------------------------------------------------------------------
void RooUnfoldBayes::Print(Option_t* option) const
{
  RooUnfold::Print (option);
  if (_nc<=0 || _ne<=0) return;

  // Print out some useful info of progress so far

  cout << "-------------------------------------------" << endl;
  cout << "Unfolding Algorithm" << endl;
  cout << "Generated (Training):" << endl;
  cout << "  Total Number of bins   : " << _nc << endl;
  cout << "  Total Number of events : " << _nCi.Sum() << endl;

  cout << "Measured (Training):" << endl;
  cout << "  Total Number of bins   : " << _ne << endl;

  cout << "Input (for unfolding):" << endl;
  cout << "  Total Number of events : " << _nEstj.Sum() << endl;

  cout << "Output (unfolded):" << endl;
  cout << "  Total Number of events : " << _nbarCi.Sum() <<endl;

  cout << "-------------------------------------------\n" << endl;

  if (((_nEstj.Sum())!=0) || ((_nCi.Sum())!=0)) {
    Int_t iend = min(_nCi.GetNrows(),_nEstj.GetNrows());
    cout << "    \tTrain \tTest\tUnfolded"<< endl;
    cout << "Bin \tTruth \tInput\tOutput"<< endl;
    for (Int_t i=0; i < iend ; i++) {
      if ((_nCi[i] == 0) && (_nEstj[i] == 0) &&
          (_nEstj[i] == 0) && (_nbarCi[i]==0)) continue;
      cout << i << "\t" << _nCi[i]                                      \
           << "\t " << _nEstj[i] << "\t " << _nbarCi[i] << endl;
    }

    // if the number of bins is different
    if (_nCi.GetNrows() > _nEstj.GetNrows() ) {
      for (Int_t i=iend; i < _nCi.GetNrows() ; i++) {
        cout << i << "\t " << _nCi[i] << endl;
      }
    }

    cout << "--------------------------------------------------------" << endl;
    cout << " \t" << (_nCi.Sum())
         << "\t " << (_nEstj.Sum()) << "\t " << (_nbarCi.Sum()) << endl;
    cout << "--------------------------------------------------------\n" << endl;
  }
}


TH1D* RooUnfoldBayes::getChi2iter ()
{
  TH1D* h_chi2iter=(TH1D*)(_chi2iter->Clone(
					    "h_chi2iter_clone"
					    )
			   );
  h_chi2iter->SetTitle("#chi^{2} between iteration i and i-1");
  h_chi2iter->GetXaxis()->SetTitle("iteration #i");
  h_chi2iter->GetYaxis()->SetTitle("#chi^{2}_{i}");
  return h_chi2iter;
}
