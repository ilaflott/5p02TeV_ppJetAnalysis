//#include "plotNfit.h"
#include "printPlots.h"
#include "fitsNfuncs.h"
#include "binsNConsts.h"

const bool debugMode=true;
//const std::string defOutFilename="printPlots_jtCalibHists_ppMC_defOut";

int printPlots_jtCalibHists(const std::string inFile_MC_dir,const std::string baseName){

  // root style settings.
  std::cout<<std::endl<<"loading style..." <<std::endl;
  LoadStyle();  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);  
  
  
  // input stuff. //

  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+3;
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=inFile_MC_dir.find("Jets");;

  //make final jetType strings
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  const std::string fullJetType="ak"+radiusInt+jetType;//"ak3PFJets";//find in dir name
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  //std::string doJetID=inFile_MC_dir.substr( inFile_MC_dir.find("MCJEC_JetID")+11 );
  //if(debugMode)std::cout<<"doJetID="<<doJetID<<std::endl;
  std::string doJetID="1";

  // input file
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  std::cout<<std::endl<<"opening input file in dir "<< inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  //TFile *finPP=new TFile( (SCRATCH_BASE+inFile_MC_dir+inFile_MC_name).c_str() );
  TFile *finPP=new TFile( (CMSSW_BASE+inFile_MC_dir+inFile_MC_name).c_str() );

  // input hists
  TH1F *hrsp[Nrad][nbins_pt];  
  double array_mean[Nrad][nbins_pt];
  double array_sig[Nrad][nbins_pt];
  

  // output stuff. //

  // output hists
  TH1F *hMean[Nrad], *hSigma[Nrad];   

  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  TF1 *fgaus=NULL;   


  // FOR JER HISTS ACROSS ALL GENPT BINS
  // ABS(RECOETA)<2, RECOPT>15, REFID!=0, JET ID OPTIONAL

  // first loop, pp loop over available radii
  for(int nj=0; nj<Nrad; nj++){

    const std::string algname=algo+srad[nj]+jetType;
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;


    // initialize mean hist
    hMean[nj] = new TH1F(Form( "hMean%d", nj),
			 Form( "Mean %s %s", algname.c_str(), ccent[nj]),
			 nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);
    hMean[nj]->SetMarkerColor(1);
    hMean[nj]->SetMarkerStyle(20);
    hMean[nj]->SetLineColor(1);
    hMean[nj]->SetMarkerSize(1.1);
    MakeHistMean(hMean[nj],1.082,0.858); 
    hMean[nj]->SetAxisRange(0.9,1.1, "Y");
    hMean[nj]->SetAxisRange(10,500, "X");

    // initialize RMS/Sigma hist
    hSigma[nj] = new TH1F(Form( "hSigma%d", nj),
			  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),
			  nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);
    hSigma[nj]->SetMarkerColor(1);
    hSigma[nj]->SetMarkerStyle(20);
    hSigma[nj]->SetLineColor(1);
    hSigma[nj]->SetMarkerSize(1.1);
    MakeHistRMS(hSigma[nj],0.563,0.001); 
    hSigma[nj]->SetAxisRange(0.0,0.5, "Y");
    hSigma[nj]->SetAxisRange(10,500, "X");


    // pt bin hJER-fit loop
    std::cout<<"looping over nbins_pt="<<nbins_pt<<std::endl;
    for(int ip=0; ip<nbins_pt; ip++){

      // input hist title string
      if(debugMode)std::cout<<"for pt range "<<ptbins[ip]<<" to "<<ptbins[ip+1]<<std::endl;
      std::string inputHistName="hJER_"+doJetID+"wJetID_"+std::to_string(ptbins[ip])+
	"_pt_"+std::to_string(ptbins[ip+1]);//"hJER_"+std::to_string(ptbins[ip])+"_pt_"+std::to_string(ptbins[ip+1]);      

      // open the input hist
      hrsp[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );    
      if(!hrsp[nj][ip]){std::cout<<"no input hist, exiting..."<<std::endl;assert(false);}      

      // print the details...
      std::cout<<std::endl;  hrsp[nj][ip]->Print("base");
      hrsp[nj][ip]->Scale( 1./ hrsp[nj][ip]->Integral() );
      //hrsp[nj][ip]->Rebin( 4 );
      
      // get some quick stats
      double norm  = hrsp[nj][ip]->GetMaximumStored();
      double err   = hrsp[nj][ip]->GetStdDevError();
      double mean  = hrsp[nj][ip]->GetMean();
      double emean = hrsp[nj][ip]->GetMeanError();
      double sig   = hrsp[nj][ip]->GetStdDev()/mean;
      double esig  = pow(1/mean,2)*pow(err,2) + pow( -sig/pow(mean,2), 2)*pow(emean, 2);
      int fitstatus = 0;
      
      //fgaus = new TF1("fgaus","gaus", mean - 2.15*hrsp[nj][ip]->GetRMS(), mean + 2.15*hrsp[nj][ip]->GetRMS());
      //fgaus = new TF1("fgaus","gaus", 0.80, 1.20);
      fgaus = new TF1("fgaus","gaus", 0.50,1.50);
      fgaus->SetParameters(norm, 0.9999, 0.5);
      //fgaus->SetParameters(norm, 1.00, 0.6);
      //fgaus->SetParameters(norm, 0.9999, 0.5);
      //fgaus->SetParLimits(1,0.60,1.30);
      //fgaus->SetParLimits(2,0.20,0.80);
      //fgaus->SetParLimits(1,0.92,1.05);
      //fgaus->SetParLimits(2,0.20,0.80);
      
      std::cout<<"for inputHist "<<inputHistName<<", we have.."<<std::endl;
      std::cout<< "Mean= "<< hrsp[nj][ip]->GetMean()<<" +/- " <<hrsp[nj][ip]->GetMeanError()<< std::endl;

      std::cout<<"fitting..."<<std::endl;
      fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");
      std::cout<< "Fit Status: "<< fitstatus<< ", Fit Error: "<< fgaus->GetParError(1)<< std::endl;

      std::cout<<"grabbing metrics (whatever this means)...."<<std::endl;
      mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     : fgaus->GetParameter(1);
      emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): fgaus->GetParError(1) ; // hrsp[nj][ip]->GetMeanError(); //debug
      sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
      esig  = (fitstatus!=0) 
	? sqrt( (pow(1/mean,2) * pow(hrsp[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ) 
	: sqrt( (pow(1/mean,2) * pow(hrsp[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ); 
      
      // set contents+errors in specific bin
      if(debugMode)std::cout<<"setting bin content for hMean[nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;
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


      array_mean[nj][ip]=mean ; //for a striaght line for plots later.
      //array_emean[nj][ip]=emean;
      array_sig[nj][ip]=sig  ;  
      //array_esig[nj][ip]=esig ;       
    }// end fit-loop over ptbins
  }// end loop over available jet radii
  if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;
  gPad->Close();

  // DRAW PDFS. //
  std::string jobType="_MCJEC_JetID"+doJetID;
  if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;
  { 
    std::string thePDFFileName=CMSSW_BASE+JERandJES_outdir+fullJetType+jobType+"_"+baseName+".pdf";
    std::string open_thePDFFileName=thePDFFileName+"[";
    std::string close_thePDFFileName=thePDFFileName+"]";

    std::cout<<"opening output pdf file "<< thePDFFileName<<std::endl<<std::endl;
    TCanvas* pdfoutCanv_wLog=new TCanvas("outputPdfwLog","outputPdfwLog", 800, 800);
    pdfoutCanv_wLog->Print( open_thePDFFileName.c_str() );
    pdfoutCanv_wLog->cd();
    pdfoutCanv_wLog->SetLogy(1);

    TCanvas* pdfoutCanv=new TCanvas("outputPdf","outputPdf", 800, 600);
    pdfoutCanv->cd();
    pdfoutCanv->SetLogy(0);
    pdfoutCanv->Divide(1,2);


    // draw and print pdfs. //
    std::cout<<std::endl<<"printing hists to pdf file..."<<std::endl;
    for(int i=0;i<Nrad;++i){
      pdfoutCanv->cd(1);
      hMean[i] ->Draw("e"); //pdfoutCanv->Print(thePDFFileName.c_str());
      TLine* meanLine=new TLine(15,1.,500,1.);
      meanLine->Draw();

      pdfoutCanv->cd(2);
      hSigma[i]->Draw("e"); 
      TLine* sigmaLine=new TLine(15,0.1,500,0.1);
      sigmaLine->Draw();

      pdfoutCanv->Print(thePDFFileName.c_str());    


      pdfoutCanv_wLog->cd();
      for(int j=0;j<nbins_pt;++j){
	std::string hrspTitle="hJER_"+std::to_string(ptbins[j])+"_pt_"+std::to_string(ptbins[j+1]);
	hrsp[i][j]->SetTitle(hrspTitle.c_str());
	hrsp[i][j]->SetMarkerStyle(8);
	hrsp[i][j]->SetMarkerSize(1.1);
	hrsp[i][j]->GetXaxis()->SetTitle("recpt/genpt");
	//hrsp[i][j]->SetAxisRange(0.,0.130,"Y");
	hrsp[i][j]->Draw("e");

	float textx=0.5,texty=0.25;
	TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean[i][j])).c_str() );
	meanText->SetNDC();
	meanText->Draw();
	TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig[i][j])).c_str() );
	sigmaText->SetNDC();
	sigmaText->Draw();
	//TLine* histMeanLine=new TLine( array_mean[i][j], 0., // x1,y1 
	//array_mean[i][j], 0.130 );// x2,y2
	//array_mean[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2
	//histMeanLine->Draw("same");
	
	pdfoutCanv_wLog->Print(thePDFFileName.c_str());
      }
    }
    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;
    pdfoutCanv->Print( close_thePDFFileName.c_str() );
  }

  // save output root file. //

  {

    std::string outFileName=CMSSW_BASE+JERandJES_outdir+fullJetType+jobType+"_"+baseName+".root";
    std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
    TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE"); 
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
  if(argc!=3){
    std::cout<<"no defaults. Do...."<<std::endl;
    std::cout<<"./printPlots_jtCalibHists.exe "<<
      "<target_ppMC_dir> <outputNameBase>"<<std::endl<<std::endl;
    return rStatus;  }
  
  rStatus=1;
  rStatus=printPlots_jtCalibHists((const std::string)argv[1], (const std::string)argv[2] );
  return rStatus;
}


//// ALT METRICS
// mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     :  0.5*(fgaus->GetParameter(1) + hrsp[nj][ip]->GetMean());
// mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     : fgaus->GetParameter(1);
// emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): hrsp[nj][ip]->GetMeanError(1);
// emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): sqrt(pow(fgaus->GetParError(1),2) + pow(hrsp[nj][ip]->GetMeanError(),2));      
// sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetRMS()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
// esig  = (fitstatus!=0) 
// ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
// : sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
// esig  = (fitstatus!=0) 
// ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
// : sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
//       
      
//hrsp[nj][ip]->Write();
