#include "JERandJES_plotNfit.h"
#include "JERandJES_binsNconsts.h"
#include "JERandJES_fitsNfuncs.h"

const bool debugMode=true;
const std::string defOutFilename="plotNfit_ppMC_QCDDijetMC_defOut";
const std::string outRootFilename=defOutFilename+".root";
const std::string outPdfFilename=defOutFilename+".pdf";

int JERandJES_plotNfit(){

  // root style settings
  std::cout<<std::endl<<"loading style..." <<std::endl;
  LoadStyle();  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);  

  // fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 

  // input file + input hist + input jet analyzer
  std::cout<<std::endl<<"opening input file in dir "<< inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  TFile *finPP=new TFile(inFile_MC.c_str());
  TH1F *hrsp[Nrad][nbins_pt];  
  double array_mean[Nrad][nbins_pt];
  
  // output hists+fits
  TH1F *hMean[Nrad], *hSigma[Nrad];   // mean/sigma hists
  TF1 *fgaus=0;   // gaussian fit

  // first loop, pp loop over available radii
  for(int nj=0; nj<Nrad; nj++){
    const std::string algname=algo+srad[nj]+jetType;
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;

    // declare Mean hist
    hMean[nj] = new TH1F(Form( "hMean%d", nj),
			 Form( "Mean %s %s", algname.c_str(), ccent[nj]),
			  nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);
    //Error in <TAxis::TAxis::Set>: bins must be in increasing order    //std::cout<<"checkDebug"<<std::endl;

    hMean[nj]->SetMarkerColor(1);
    hMean[nj]->SetMarkerStyle(20);
    hMean[nj]->SetLineColor(1);
    hMean[nj]->SetMarkerSize(1.3);
    MakeHistMean(hMean[nj],1.082,0.858); 
    //MakeHistMean(hMean[nj],1.022,0.978); 
    //MakeHistMean(hMean[nj],1.082,0.898); 
    //MakeHistMean(hMean[nj],1.148,0.848); 
    //MakeHistMean(hMean[nj],1.052,0.934); 
    
    // declare RMS/Sigma hist
    hSigma[nj] = new TH1F(Form( "hSigma%d", nj),
			  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),
			  nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);
    //Error in <TAxis::TAxis::Set>: bins must be in increasing order //std::cout<<"checkDebug"<<std::endl;    

    hSigma[nj]->SetMarkerColor(1);
    hSigma[nj]->SetMarkerStyle(20);
    hSigma[nj]->SetLineColor(1);
    hSigma[nj]->SetMarkerSize(1.3);
    MakeHistRMS(hSigma[nj],0.563,0.001); 

    // pt bin hJER-fit loop
    std::cout<<"looping over nbins_pt= "<<nbins_pt<<" input JER hists + fitting..." <<std::endl;
    for(int ip=0; ip<nbins_pt; ip++){

      std::string inputHistName="hJER_"+std::to_string(ptbins[ip])+"_pt_"+std::to_string(ptbins[ip+1]);
      hrsp[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );    std::cout<<std::endl; 
      hrsp[nj][ip]->Print("base");
      hrsp[nj][ip]->Scale( 1./ hrsp[nj][ip]->Integral() );
      hrsp[nj][ip]->Rebin( 4 );
      
      // get some quick stats
      double norm  = hrsp[nj][ip]->GetMaximumStored();
      double err   = hrsp[nj][ip]->GetStdDevError();
      double mean  = hrsp[nj][ip]->GetMean();
      double emean = hrsp[nj][ip]->GetMeanError();
      double sig   = hrsp[nj][ip]->GetStdDev()/mean;
      double esig  = pow(1/mean,2)*pow(err,2) + pow( -sig/pow(mean,2), 2)*pow(emean, 2);
      
      // i dont know what this part is doing or what is has to do with the radius
      //// R 0.3
      //fgaus = new TF1("fgaus","gaus", 0.50,1.50);
      //fgaus->SetParameters(norm, 1.00, 0.3);            

      // R 0.2
      fgaus = new TF1("fgaus","gaus", 0.50,1.50);
      fgaus->SetParameters(norm, 0.9999, 0.5);
            
      //fgaus = new TF1("fgaus","gaus", mean - 2.15*hrsp[nj][ip]->GetRMS(), mean + 2.15*hrsp[nj][ip]->GetRMS());
      //fgaus->SetParameters(norm, 1.00, 0.6);
      //fgaus->SetParLimits(1,0.60,1.30);
      //fgaus->SetParLimits(2,0.20,0.80);
      
      //fgaus = new TF1("fgaus","gaus", 0.80, 1.20);
      //fgaus->SetParameters(norm, 0.9999, 0.5);
      //fgaus->SetParLimits(1,0.92,1.05);
      //fgaus->SetParLimits(2,0.20,0.80);
      
      int fitstatus = 0;
      fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");
      
      std::cout<< "Fit Status: "<< fitstatus<< ", Fit Error: "<< fgaus->GetParError(1)<< std::endl;
      std::cout<< "Hist Error: "<< hrsp[nj][ip]->GetMeanError()<< std::endl;
      
      // grab metrics
      mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     : fgaus->GetParameter(1);
      emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): fgaus->GetParError(1) ; // hrsp[nj][ip]->GetMeanError(); //debug
      sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
      esig  = (fitstatus!=0) 
	? sqrt( (pow(1/mean,2) * pow(hrsp[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ) 
	: sqrt( (pow(1/mean,2) * pow(hrsp[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ); 	
      
      //// ALT METRICS
      // mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     :  0.5*(fgaus->GetParameter(1) + hrsp[nj][ip]->GetMean());
      // mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     : fgaus->GetParameter(1);
      // emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): hrsp[nj][ip]->GetMeanError(1);
      // emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): sqrt(pow(fgaus->GetParError(1),2) + pow(hrsp[nj][ip]->GetMeanError(),2));      
      // sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetRMS()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
      // esig  = (fitstatus!=0) 
      // 	? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
      // 	: sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
      // esig  = (fitstatus!=0) 
      // 	? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
      // 	: sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
      //       
      
      //hrsp[nj][ip]->Write();

      //array_norm[nj][ip]=norm ; 
      //array_err[nj][ip]=err  ;  
      array_mean[nj][ip]=mean ; 
      //array_emean[nj][ip]=emean;
      //array_sig[nj][ip]=sig  ;  
      //array_esig[nj][ip]=esig ; 



      // set contents+errors in specific bin
      if(debugMode)std::cout<<"setting bin content for [nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;
      if(debugMode)std::cout<<"mean ="<< mean<<std::endl;
      if(debugMode)std::cout<<"emean ="<< emean<<std::endl;
      hMean[nj]->SetBinContent (ip+1, mean);
      hMean[nj]->SetBinError   (ip+1, emean);
      // hMean[nj]->SetBinContent (ip+1, roundoff(mean));
      // hMean[nj]->SetBinError   (ip+1, roundoff(emean));

      if(debugMode)std::cout<<"sig ="<< sig<<std::endl;
      if(debugMode)std::cout<<"esig ="<< esig<<std::endl;
      hSigma[nj]->SetBinContent (ip+1, sig);
      hSigma[nj]->SetBinError   (ip+1, esig);
    }// end fit-loop over ptbins
  }// end loop over available jet radii

  std::cout<<std::endl<<"closing gPad..."<<std::endl;
  gPad->Close();

  {
    std::cout<<std::endl<<"opening output pdf file "<< outPdfFilename<<std::endl;
    TCanvas* pdfoutCanv=new TCanvas("outputPdf","outputPdf", 800, 800);
    std::string open_outPdfFilename=outPdfFilename+"[";
    std::string close_outPdfFilename=outPdfFilename+"]";
    pdfoutCanv->Print( open_outPdfFilename.c_str() );
    pdfoutCanv->cd();
    std::cout<<std::endl<<"printing hists to pdf file..."<<std::endl;
    for(int i=0;i<Nrad;++i){
      hMean[i] ->Draw(); pdfoutCanv->Print(outPdfFilename.c_str());
      hSigma[i]->Draw(); pdfoutCanv->Print(outPdfFilename.c_str());    
      for(int j=0;j<nbins_pt;++j){
	std::string hrspTitle="hJER_"+std::to_string(ptbins[j])+"_pt_"+std::to_string(ptbins[j+1]);	
	hrsp[i][j]->SetTitle(hrspTitle.c_str());
	hrsp[i][j]->SetMarkerStyle(8);
	hrsp[i][j]->SetMarkerSize(1.1);
	hrsp[i][j]->GetXaxis()->SetTitle("recpt/genpt");
	hrsp[i][j]->Draw("e");
	TLine* histMeanLine=new TLine( array_mean[i][j], 0., // x1,y1 
				       array_mean[i][j], hrsp[i][j]->GetYaxis()->GetXmax() );// x2,y2
	histMeanLine->Draw("same");
	
	pdfoutCanv->Print(outPdfFilename.c_str());
      }
    }
    std::cout<<"closing output pdf file "<<outPdfFilename<<std::endl<<std::endl;
    pdfoutCanv->Print( close_outPdfFilename.c_str() );
  }
  

  {
    std::cout<<std::endl<<"opening output root file "<< outRootFilename<<std::endl;
    TFile *rootfout = new TFile(outRootFilename.c_str(),"RECREATE"); 
    rootfout->cd();  
    std::cout<<"writing output to root file..."<<std::endl;
    for(int i=0;i<Nrad;++i){
      hMean[i]->Write();
      hSigma[i]->Write();    
      for(int j=0;j<nbins_pt;++j) hrsp[i][j]->Write();
    }
    std::cout<<"closing output root file"<<std::endl<<std::endl;
    rootfout->Close();
  }
  
    return 0;
}


int main(int argc, char*argv[]){
  int rStatus=-1;
  if(argc!=1){
    std::cout<<"no arguments, just do "<<std::endl<<"./JERandJES_plotNfit.exe"<<std::endl<<std::endl;
    return rStatus;
  }
  rStatus=1;
  rStatus=JERandJES_plotNfit();
  std::cout<<"done, rStatus="<<rStatus <<std::endl<<std::endl;
  return rStatus;
}
