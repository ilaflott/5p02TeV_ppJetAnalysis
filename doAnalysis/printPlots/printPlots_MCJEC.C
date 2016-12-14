//#include "plotNfit.h"
#include "printPlots.h"
#include "printPlots_fitsNfuncs.h"
#include "printPlots_binsNconsts.h"

const bool debugMode=true;
//const std::string defOutFilename="printPlots_MCJEC_ppMC_defOut";

int printPlots_MCJEC(const std::string inFile_MC_dir,const std::string baseName){

  // root style settings.
  std::cout<<std::endl<<"loading style..." <<std::endl;
  LoadStyle();  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);  
  
  
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+strlen("_ak");
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
  std::string didJetID=inFile_MC_dir.substr( inFile_MC_dir.find("MCJEC_jtID")+strlen("MCJEC_jtID"), 1 );
  //if(debugMode)std::cout<<"doJetID="<<doJetID<<std::endl;
  
  std::string doJetID;
  if(didJetID=="0") doJetID=didJetID;
  else if (didJetID=="L"||didJetID=="T") doJetID="1";
  
  std::string jobType="_MCJEC_jtID"+didJetID;
  if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;

  // input file
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  std::cout<<std::endl<<"opening input file in dir "<< inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  TFile *finPP=new TFile( (inputDir+inFile_MC_dir+inFile_MC_name).c_str() );
  
  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  TF1 *fgaus=NULL;   
  
  

  // input hists
  TH1F *hrsp[Nrad][nbins_pt];  
  double array_mean[Nrad][nbins_pt];
  double array_sig[Nrad][nbins_pt];
  
  // output hists
  TH1F *hMean[Nrad], *hSigma[Nrad];   
  
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
      if(!hrsp[nj][ip]){std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;assert(false);}      

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


  {     
    // DRAW THOSE PDFS //
    std::string thePDFFileName=outputDir+fullJetType+jobType+"_"+baseName+".pdf";
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
  std::string outFileName=outputDir+fullJetType+jobType+"_"+baseName+".root";
  std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
  TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE"); 

    {
      rootfout->cd();  
      std::cout<<"writing output to root file..."<<std::endl;
      for(int i=0;i<Nrad;++i){
	for(int j=0;j<nbins_pt;++j) hrsp[i][j]->Write();
	hMean[i]->Write();
	hSigma[i]->Write();          
      }
      //std::cout<<"closing output root file"<<std::endl<<std::endl;
      //rootfout->Close();
    }
  



  ////     JER/JES PLOTS, 150 < GENPT < 200        ////
  
  // input hists
  TH1F *hrsp_gen150[Nrad][nbins_eta]={};  
  double array_mean_gen150[Nrad][nbins_eta];
  double array_sig_gen150[Nrad][nbins_eta];
  

  // output hists
  TH1F *hMean_gen150[Nrad]={}, *hSigma_gen150[Nrad]={};   

  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  TF1 *fgaus_gen150=NULL;   


  // FOR JER 30to50 and 150to200 HISTS ACROSS ALL GENETA BINS
  // SUBJECT TO RECO CUTS: 
  // ABS(RECOETA)<2, RECOPT>15, REFID!=0, JET ID OPTIONAL

  // first loop, pp loop over available radii
  for(int nj=0; nj<Nrad; nj++){

    const std::string algname=algo+srad[nj]+jetType;
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;


    // initialize mean hist
    hMean_gen150[nj] = new TH1F(Form( "hMean_gen150%d", nj),
			 Form( "Mean %s %s", algname.c_str(), ccent[nj]),
			 nbins_eta, etabins); //nbins_pt, ptbins_bound);
    hMean_gen150[nj]->SetMarkerColor(1);
    hMean_gen150[nj]->SetMarkerStyle(20);
    hMean_gen150[nj]->SetLineColor(1);
    hMean_gen150[nj]->SetMarkerSize(1.1);
    MakeHistMean(hMean_gen150[nj],1.082,0.858); 
    hMean_gen150[nj]->SetAxisRange(0.8,1.2, "Y");
    hMean_gen150[nj]->SetAxisRange(-5.3,5.3, "X");

    // initialize RMS/Sigma hist
    hSigma_gen150[nj] = new TH1F(Form( "hSigma_gen150%d", nj),
			  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),
			 nbins_eta, etabins); //nbins_pt, ptbins_bound);
    hSigma_gen150[nj]->SetMarkerColor(1);
    hSigma_gen150[nj]->SetMarkerStyle(20);
    hSigma_gen150[nj]->SetLineColor(1);
    hSigma_gen150[nj]->SetMarkerSize(1.1);
    MakeHistRMS(hSigma_gen150[nj],0.563,0.001); 
    hSigma_gen150[nj]->SetAxisRange(0.0,0.5, "Y");
    hSigma_gen150[nj]->SetAxisRange(-5.3,5.3, "X");


    // pt bin hJER-fit loop
    std::cout<<"looping over nbins_eta="<<nbins_eta<<std::endl;
    std::cout<<"for pt range 150 to 200"<<std::endl;
    for(int ip=0; ip<nbins_eta; ip++){

      //geneta cut for these hists
      //if( abs(etabins[ip])>2.1 )continue;

      // input hist title string
      //std::string inputHistName="hJER_"+doJetID+"wJetID_etabin"+std::to_string(etabins[ip])+"_150_pt_200";
      std::string inputHistName="hJER_"+doJetID+"wJetID_etabin"+std::to_string(ip)+"_150_pt_200";

      // open the input hist
      hrsp_gen150[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );    
      if(!hrsp_gen150[nj][ip]){std::cout<<"no input hist, exiting..."<<std::endl;assert(false);}      
      if(hrsp_gen150[nj][ip]->GetEntries()<5){//not enough entries
	std::cout<<"hist has less than five entries... continuing"<<std::endl;
	hrsp_gen150[nj][ip]=NULL;
	continue;      }

      // print the details...
      std::cout<<std::endl;  
      std::cout<<"for eta range "<< etabins[ip]<<" to "<<etabins[ip+1] <<std::endl;
      hrsp_gen150[nj][ip]->Print("base");
      hrsp_gen150[nj][ip]->Scale( 1./ hrsp_gen150[nj][ip]->Integral() );
      //hrsp_gen150[nj][ip]->Rebin( 4 );
      
      // get some quick stats
      double norm  = hrsp_gen150[nj][ip]->GetMaximumStored();
      double err   = hrsp_gen150[nj][ip]->GetStdDevError();
      double mean  = hrsp_gen150[nj][ip]->GetMean();
      double emean = hrsp_gen150[nj][ip]->GetMeanError();
      double sig   = hrsp_gen150[nj][ip]->GetStdDev()/mean;
      double esig  = pow(1/mean,2)*pow(err,2) + pow( -sig/pow(mean,2), 2)*pow(emean, 2);
      int fitstatus = 0;
      
      //fgaus_gen150 = new TF1("fgaus_gen150","gaus", mean - 2.15*hrsp_gen150[nj][ip]->GetRMS(), mean + 2.15*hrsp_gen150[nj][ip]->GetRMS());
      //fgaus_gen150 = new TF1("fgaus_gen150","gaus", 0.80, 1.20);
      fgaus_gen150 = new TF1("fgaus_gen150","gaus", 0.50,1.50);
      fgaus_gen150->SetParameters(norm, 0.9999, 0.5);
      //fgaus_gen150->SetParameters(norm, 1.00, 0.6);
      //fgaus_gen150->SetParameters(norm, 0.9999, 0.5);
      //fgaus_gen150->SetParLimits(1,0.60,1.30);
      //fgaus_gen150->SetParLimits(2,0.20,0.80);
      //fgaus_gen150->SetParLimits(1,0.92,1.05);
      //fgaus_gen150->SetParLimits(2,0.20,0.80);
      
      std::cout<<"for inputHist "<<inputHistName<<", we have.."<<std::endl;
      std::cout<< "Mean= "<< hrsp_gen150[nj][ip]->GetMean()<<" +/- " <<hrsp_gen150[nj][ip]->GetMeanError()<< std::endl;

      std::cout<<"fitting..."<<std::endl;
      fitstatus = hrsp_gen150[nj][ip]->Fit(fgaus_gen150,"RQ");
      std::cout<< "Fit Status: "<< fitstatus<< ", Fit Error: "<< fgaus_gen150->GetParError(1)<< std::endl;

      std::cout<<"grabbing metrics (whatever this means)...."<<std::endl;
      mean  = (fitstatus!=0) ? hrsp_gen150[nj][ip]->GetMean()     : fgaus_gen150->GetParameter(1);
      emean = (fitstatus!=0) ? hrsp_gen150[nj][ip]->GetMeanError(): fgaus_gen150->GetParError(1) ; // hrsp_gen150[nj][ip]->GetMeanError(); //debug
      sig   = (fitstatus!=0) ? hrsp_gen150[nj][ip]->GetStdDev()/mean : fgaus_gen150->GetParameter(2)/fgaus_gen150->GetParameter(1);
      esig  = (fitstatus!=0) 
	? sqrt( (pow(1/mean,2) * pow(hrsp_gen150[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp_gen150[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ) 
	: sqrt( (pow(1/mean,2) * pow(hrsp_gen150[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp_gen150[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ); 
      
      // set contents+errors in specific bin
      if(debugMode)std::cout<<"setting bin content for hMean_gen150[nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;
      if(debugMode)std::cout<<"mean ="<< mean<<std::endl;
      if(debugMode)std::cout<<"emean ="<< emean<<std::endl;
      hMean_gen150[nj]->SetBinContent (ip+1, mean);
      hMean_gen150[nj]->SetBinError   (ip+1, emean);
      // hMean_gen150[nj]->SetBinContent (ip+1, roundoff(mean));
      // hMean_gen150[nj]->SetBinError   (ip+1, roundoff(emean));

      if(debugMode)std::cout<<"sig ="<< sig<<std::endl;
      if(debugMode)std::cout<<"esig ="<< esig<<std::endl;
      hSigma_gen150[nj]->SetBinContent (ip+1, sig);
      hSigma_gen150[nj]->SetBinError   (ip+1, esig);


      array_mean_gen150[nj][ip]=mean ; //for a striaght line for plots later.
      //array_emean[nj][ip]=emean;
      array_sig_gen150[nj][ip]=sig  ;  
      //array_esig[nj][ip]=esig ;       
    }// end fit-loop over ptbins
  }// end loop over available jet radii
  if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;
  gPad->Close();

  // DRAW PDFS. //
  //std::string jobType="_MCJEC_JetID"+doJetID;
  //if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;
  { 
    std::string thePDFFileName=outputDir+fullJetType+jobType+"_150to200_"+baseName+".pdf";
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
      hMean_gen150[i] ->Draw("e"); //pdfoutCanv->Print(thePDFFileName.c_str());
      TLine* meanLine=new TLine(-5.3,1.,5.3,1.);
      meanLine->Draw();

      pdfoutCanv->cd(2);
      hSigma_gen150[i]->Draw("e"); 
      TLine* sigmaLine=new TLine(-5.3,0.1,5.3,0.1);
      sigmaLine->Draw();
      pdfoutCanv->Print(thePDFFileName.c_str());    
      
      pdfoutCanv_wLog->cd();
      for(int j=0;j<nbins_eta;++j){
	if( !(hrsp_gen150[i][j]) )continue;
	std::string hrsp_gen150Title="JER, 150<genpt<200, "+std::to_string(etabins[j])+"<geneta<"+std::to_string(etabins[j+1]);
	hrsp_gen150[i][j]->SetTitle(hrsp_gen150Title.c_str());
	hrsp_gen150[i][j]->SetMarkerStyle(8);
	hrsp_gen150[i][j]->SetMarkerSize(1.1);
	hrsp_gen150[i][j]->GetXaxis()->SetTitle("recpt/genpt");
	//hrsp_gen150[i][j]->SetAxisRange(0.,0.130,"Y");
	hrsp_gen150[i][j]->Draw("e");

	float textx=0.5,texty=0.25;
	TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean_gen150[i][j])).c_str() );
	meanText->SetNDC();
	meanText->Draw();
	TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig_gen150[i][j])).c_str() );
	sigmaText->SetNDC();
	sigmaText->Draw();
//TLine* histMeanLine=new TLine( array_mean_gen150[i][j], 0., // x1,y1 
	//array_mean_gen150[i][j], 0.130 );// x2,y2
	//array_mean_gen150[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2
	//histMeanLine->Draw("same");
	
	pdfoutCanv_wLog->Print(thePDFFileName.c_str());
      }
    }
    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;
    pdfoutCanv->Print( close_thePDFFileName.c_str() );
  }

  // save output root file. //

  {

    //std::string outFileName=outputDir+fullJetType+jobType+"_150to200_"+baseName+".root";
    //std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
    //TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE"); 
    rootfout->cd();  
    std::cout<<"writing output to root file..."<<std::endl;
    for(int i=0;i<Nrad;++i){
      hMean_gen150[i]->Write();
      hSigma_gen150[i]->Write();    
      for(int j=0;j<nbins_eta;++j) {
	if( !(hrsp_gen150[i][j]) ) continue;
	else hrsp_gen150[i][j]->Write();	
      }
    }

  }


 
  ////     JER/JES PLOTS, 30 < GENPT < 50        ////

  // input hists
  TH1F *hrsp_gen30[Nrad][nbins_eta]={};  
  double array_mean_gen30[Nrad][nbins_eta];
  double array_sig_gen30[Nrad][nbins_eta];
  

  // output hists
  TH1F *hMean_gen30[Nrad]={}, *hSigma_gen30[Nrad]={};   

  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  TF1 *fgaus_gen30=NULL;   


  // first loop, pp loop over available radii
  for(int nj=0; nj<Nrad; nj++){

    const std::string algname=algo+srad[nj]+jetType;
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;


    // initialize mean hist
    hMean_gen30[nj] = new TH1F(Form( "hMean_gen30%d", nj),
			 Form( "Mean %s %s", algname.c_str(), ccent[nj]),
			 nbins_eta, etabins); //nbins_pt, ptbins_bound);
    hMean_gen30[nj]->SetMarkerColor(1);
    hMean_gen30[nj]->SetMarkerStyle(20);
    hMean_gen30[nj]->SetLineColor(1);
    hMean_gen30[nj]->SetMarkerSize(1.1);
    MakeHistMean(hMean_gen30[nj],1.082,0.858); 
    hMean_gen30[nj]->SetAxisRange(0.8,1.2, "Y");
    hMean_gen30[nj]->SetAxisRange(-5.3,5.3, "X");

    // initialize RMS/Sigma hist
    hSigma_gen30[nj] = new TH1F(Form( "hSigma_gen30%d", nj),
			  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),
			 nbins_eta, etabins); //nbins_pt, ptbins_bound);
    hSigma_gen30[nj]->SetMarkerColor(1);
    hSigma_gen30[nj]->SetMarkerStyle(20);
    hSigma_gen30[nj]->SetLineColor(1);
    hSigma_gen30[nj]->SetMarkerSize(1.1);
    MakeHistRMS(hSigma_gen30[nj],0.563,0.001); 
    hSigma_gen30[nj]->SetAxisRange(0.0,0.5, "Y");
    hSigma_gen30[nj]->SetAxisRange(-5.3,5.3, "X");


    // pt bin hJER-fit loop
    std::cout<<"looping over nbins_eta="<<nbins_eta<<std::endl;
    std::cout<<"for pt range 30 to 50"<<std::endl;
    for(int ip=0; ip<nbins_eta; ip++){

      //geneta cut for these hists
      //if( abs(etabins[ip])>2.1 )continue;

      // input hist title string
      //std::string inputHistName="hJER_"+doJetID+"wJetID_etabin"+std::to_string(etabins[ip])+"_30_pt_50";
      std::string inputHistName="hJER_"+doJetID+"wJetID_etabin"+std::to_string(ip)+"_30_pt_50";

      // open the input hist
      hrsp_gen30[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );    
      if(!hrsp_gen30[nj][ip]){std::cout<<"no input hist, exiting..."<<std::endl;assert(false);}      
      if(hrsp_gen30[nj][ip]->GetEntries()<5){//not enough entries
	std::cout<<"hist has less than five entries... continuing"<<std::endl;
	hrsp_gen30[nj][ip]=NULL;
	continue;      }

      // print the details...
      std::cout<<std::endl;  
      std::cout<<"for eta range "<< etabins[ip]<<" to "<<etabins[ip+1] <<std::endl;
      hrsp_gen30[nj][ip]->Print("base");
      hrsp_gen30[nj][ip]->Scale( 1./ hrsp_gen30[nj][ip]->Integral() );
      //hrsp_gen30[nj][ip]->Rebin( 4 );
      
      // get some quick stats
      double norm  = hrsp_gen30[nj][ip]->GetMaximumStored();
      double err   = hrsp_gen30[nj][ip]->GetStdDevError();
      double mean  = hrsp_gen30[nj][ip]->GetMean();
      double emean = hrsp_gen30[nj][ip]->GetMeanError();
      double sig   = hrsp_gen30[nj][ip]->GetStdDev()/mean;
      double esig  = pow(1/mean,2)*pow(err,2) + pow( -sig/pow(mean,2), 2)*pow(emean, 2);
      int fitstatus = 0;
      
      //fgaus_gen30 = new TF1("fgaus_gen30","gaus", mean - 2.15*hrsp_gen30[nj][ip]->GetRMS(), mean + 2.15*hrsp_gen30[nj][ip]->GetRMS());
      //fgaus_gen30 = new TF1("fgaus_gen30","gaus", 0.80, 1.20);
      fgaus_gen30 = new TF1("fgaus_gen30","gaus", 0.50,1.50);
      fgaus_gen30->SetParameters(norm, 0.9999, 0.5);
      //fgaus_gen30->SetParameters(norm, 1.00, 0.6);
      //fgaus_gen30->SetParameters(norm, 0.9999, 0.5);
      //fgaus_gen30->SetParLimits(1,0.60,1.30);
      //fgaus_gen30->SetParLimits(2,0.20,0.80);
      //fgaus_gen30->SetParLimits(1,0.92,1.05);
      //fgaus_gen30->SetParLimits(2,0.20,0.80);
      
      std::cout<<"for inputHist "<<inputHistName<<", we have.."<<std::endl;
      std::cout<< "Mean= "<< hrsp_gen30[nj][ip]->GetMean()<<" +/- " <<hrsp_gen30[nj][ip]->GetMeanError()<< std::endl;

      std::cout<<"fitting..."<<std::endl;
      fitstatus = hrsp_gen30[nj][ip]->Fit(fgaus_gen30,"RQ");
      std::cout<< "Fit Status: "<< fitstatus<< ", Fit Error: "<< fgaus_gen30->GetParError(1)<< std::endl;

      std::cout<<"grabbing metrics (whatever this means)...."<<std::endl;
      mean  = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetMean()     : fgaus_gen30->GetParameter(1);
      emean = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetMeanError(): fgaus_gen30->GetParError(1) ; // hrsp_gen30[nj][ip]->GetMeanError(); //debug
      sig   = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetStdDev()/mean : fgaus_gen30->GetParameter(2)/fgaus_gen30->GetParameter(1);
      esig  = (fitstatus!=0) 
	? sqrt( (pow(1/mean,2) * pow(hrsp_gen30[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp_gen30[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ) 
	: sqrt( (pow(1/mean,2) * pow(hrsp_gen30[nj][ip]->GetStdDevError(),2)) + 
		(pow(-hrsp_gen30[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) ); 
      
      // set contents+errors in specific bin
      if(debugMode)std::cout<<"setting bin content for hMean_gen30[nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;
      if(debugMode)std::cout<<"mean ="<< mean<<std::endl;
      if(debugMode)std::cout<<"emean ="<< emean<<std::endl;
      hMean_gen30[nj]->SetBinContent (ip+1, mean);
      hMean_gen30[nj]->SetBinError   (ip+1, emean);
      // hMean_gen30[nj]->SetBinContent (ip+1, roundoff(mean));
      // hMean_gen30[nj]->SetBinError   (ip+1, roundoff(emean));

      if(debugMode)std::cout<<"sig ="<< sig<<std::endl;
      if(debugMode)std::cout<<"esig ="<< esig<<std::endl;
      hSigma_gen30[nj]->SetBinContent (ip+1, sig);
      hSigma_gen30[nj]->SetBinError   (ip+1, esig);


      array_mean_gen30[nj][ip]=mean ; //for a striaght line for plots later.
      //array_emean[nj][ip]=emean;
      array_sig_gen30[nj][ip]=sig  ;  
      //array_esig[nj][ip]=esig ;       
    }// end fit-loop over ptbins
  }// end loop over available jet radii
  if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;
  gPad->Close();

  // DRAW PDFS. //
  //std::string jobType="_MCJEC_JetID"+doJetID;
  //if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;
  { 
    std::string thePDFFileName=outputDir+fullJetType+jobType+"_30to50_"+baseName+".pdf";
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
      hMean_gen30[i] ->Draw("e"); //pdfoutCanv->Print(thePDFFileName.c_str());
      TLine* meanLine=new TLine(-5.3,1.,5.3,1.);
      meanLine->Draw();

      pdfoutCanv->cd(2);
      hSigma_gen30[i]->Draw("e"); 
      TLine* sigmaLine=new TLine(-5.3,0.1,5.3,0.1);
      sigmaLine->Draw();
      pdfoutCanv->Print(thePDFFileName.c_str());    

      pdfoutCanv_wLog->cd();
      for(int j=0;j<nbins_eta;++j){
	if( !(hrsp_gen30[i][j]) )continue;
	std::string hrsp_gen30Title="JER, 30<genpt<50, "+std::to_string(etabins[j])+"<geneta<"+std::to_string(etabins[j+1]);
	hrsp_gen30[i][j]->SetTitle(hrsp_gen30Title.c_str());
	hrsp_gen30[i][j]->SetMarkerStyle(8);
	hrsp_gen30[i][j]->SetMarkerSize(1.1);
	hrsp_gen30[i][j]->GetXaxis()->SetTitle("recpt/genpt");
	//hrsp_gen30[i][j]->SetAxisRange(0.,0.130,"Y");
	hrsp_gen30[i][j]->Draw("e");

	float textx=0.5,texty=0.25;
	TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean_gen30[i][j])).c_str() );
	meanText->SetNDC();
	meanText->Draw();
	TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig_gen30[i][j])).c_str() );
	sigmaText->SetNDC();
	sigmaText->Draw();
//TLine* histMeanLine=new TLine( array_mean_gen30[i][j], 0., // x1,y1 
	//array_mean_gen30[i][j], 0.130 );// x2,y2
	//array_mean_gen30[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2
	//histMeanLine->Draw("same");
	
	pdfoutCanv_wLog->Print(thePDFFileName.c_str());
      }
    }
    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;
    pdfoutCanv->Print( close_thePDFFileName.c_str() );
  }

  // save output root file. //

  {

    //std::string outFileName=outputDir+fullJetType+jobType+"_30to50_"+baseName+".root";
    //std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
    //TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE"); 
    rootfout->cd();  
    std::cout<<"writing output to root file..."<<std::endl;
    for(int i=0;i<Nrad;++i){
      hMean_gen30[i]->Write();
      hSigma_gen30[i]->Write();    
      for(int j=0;j<nbins_eta;++j) {
	if( !(hrsp_gen30[i][j]) ) continue;
	else hrsp_gen30[i][j]->Write();	
      }
    }
  }


  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  { 
    std::cout<<" drawing MC Eff. QA Plots..."<<std::endl;
    
    std::string thePDFFileName=outputDir+fullJetType+jobType+"_MCEff_"+baseName+".pdf";
    std::string open_thePDFFileName=thePDFFileName+"[";
    std::string close_thePDFFileName=thePDFFileName+"]";
    std::cout<<std::endl<<"creating temporary canvas for printing MCEff plots..."<<std::endl;

    TCanvas *temp_canvMCEff = new TCanvas("tempMCEff", "tempMCEff", 1200, 600);
    temp_canvMCEff->Print(open_thePDFFileName.c_str());
    temp_canvMCEff->cd();
  
    //if(drawProfiles){     
      for(int j=0; j<(N_genVars); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;
  
  	  
  	std::string inHistName="hpp_mceff_"+genVars[j];
  	if(doJetID=="1")inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	std::cout<<"opening TH2F "<<inHistName << std::endl;
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	if(!the2DMCEffQAHist) {std::cout<<"no MCEff plot, continuing..."<<std::endl; continue;}
  	the2DMCEffQAHist->Print("base");
  	  
  	std::string h_Title   ="MC Eff. QA, TH2 profile";
  	if(doJetID=="1")h_Title+=", w/ JetIDCut";      
  	  
  	//std::string h_ZAx_Title="Entries"; 
  	std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
  	std::string h_YAx_Title="avg "+stringRatio+"/bin";
  	  
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Eta"       ;
  	if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
  	  
  	std::cout<<"taking the profile of "<<inHistName << std::endl;
  	TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX("MCEffProfile",1,-1, "o");
  	  
  	theMCEffQAHist->SetTitle (    h_Title.c_str() );
  	theMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	theMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	  
  	theMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)theMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
  	if(j==1)theMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
  	if(j==2)theMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
  	  
  	theMCEffQAHist->Draw("E HIST"); 
  	 
  	temp_canvMCEff->Print(thePDFFileName.c_str());   
      }// gen genvar loop
      
      for(int j=0; j<(N_genVars_ptrat); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;      
  	
  	std::string inHistName="hpp_mceff_ptrat_"+genVars_ptrat[j];
  	if(doJetID=="1")inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	
  	std::string h_Title   ="MC Eff. QA, TH2 profile";
  	if(doJetID=="1")h_Title+=", w/ JetIDCut";      
  	    
  	std::string h_YAx_Title="avg (recpt)/(genpt)/bin";	    
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Eta (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Phi (rad)"  ;
  	if(j==2)h_XAx_Title="gen Jet dRJet"      ; 
  	
  	TH1F* theMCEffQAHist= (TH1F*)the2DMCEffQAHist->TH2::ProfileX("MCEffProfile",1,-1, "o");
  	
  	theMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)theMCEffQAHist->SetAxisRange( -3.0, 3.0, "X");
  	if(j==1)theMCEffQAHist->SetAxisRange( -4.0, 4.0, "X");
  	if(j==2)theMCEffQAHist->SetAxisRange(   0., 0.5, "X");
  	
  	theMCEffQAHist->SetTitle (    h_Title.c_str() );
  	theMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	theMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	
  	theMCEffQAHist->Draw("HIST E"); 
  	
  	temp_canvMCEff->Print(thePDFFileName.c_str());   	
      }//end genvar ptrat loops      
      //}//end drawProfles



  //if(drawTH2s){     
      for(int j=0; j<(N_genVars); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;
  
  	  
  	std::string inHistName="hpp_mceff_"+genVars[j];
  	if(doJetID=="1")inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	std::cout<<"opening TH2F "<<inHistName << std::endl;
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	if(!the2DMCEffQAHist) {std::cout<<"no MCEff plot, continuing..."<<std::endl; continue;}
  	  
  	std::string h_Title   ="MC Eff. QA, TH2";
  	if(doJetID=="1")h_Title+=", w/ JetIDCut";      
  	  
  	//std::string h_ZAx_Title="Entries"; 
  	std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
  	std::string h_YAx_Title=stringRatio;
  	  
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Eta"       ;
  	if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
  	  
  	  
  	the2DMCEffQAHist->SetTitle (    h_Title.c_str() );
  	the2DMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	the2DMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	  
  	the2DMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)the2DMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
  	if(j==1)the2DMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
  	if(j==2)the2DMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
  	
  	the2DMCEffQAHist->Draw("COLZ"); 
  	 
  	temp_canvMCEff->Print(thePDFFileName.c_str());   
      }// gen genvar loop
      
      for(int j=0; j<(N_genVars_ptrat); j++){     
  	if(debugMode)std::cout<<std::endl<<"j= "<<j<<std::endl;      
  	
  	std::string inHistName="hpp_mceff_ptrat_"+genVars_ptrat[j];
  	if(doJetID=="1")inHistName+="_wJetID";
  	inHistName+="_"+radius+etaWidth;
  	  
  	TH2F* the2DMCEffQAHist= (TH2F*)finPP->Get( inHistName.c_str() );
  	
  	std::string h_Title   ="MC Eff. QA TH2";
  	if(doJetID=="1")h_Title+=", w/ JetIDCut";      
  	    
  	std::string h_YAx_Title="(recpt)/(genpt)";	    
  	std::string h_XAx_Title="";	    
  	if(j==0)h_XAx_Title="gen Jet Eta (GeV)"  ;
  	if(j==1)h_XAx_Title="gen Jet Phi (rad)"  ;
  	if(j==2)h_XAx_Title="gen Jet dRJet"      ; 
  	
  	the2DMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
  	if(j==0)the2DMCEffQAHist->SetAxisRange( -3.0, 3.0, "X");
  	if(j==1)the2DMCEffQAHist->SetAxisRange( -4.0, 4.0, "X");
  	if(j==2)the2DMCEffQAHist->SetAxisRange(   0., 0.5, "X");
  	
  	the2DMCEffQAHist->SetTitle (    h_Title.c_str() );
  	the2DMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
  	the2DMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
  	
  	the2DMCEffQAHist->Draw("COLZ"); 
  	
  	temp_canvMCEff->Print(thePDFFileName.c_str());   	
      }//end genvar ptrat loops      
      //}//end draw TH2s
    temp_canvMCEff->Print(close_thePDFFileName.c_str());   	
    temp_canvMCEff->Close();        
  }//end  MCEff Plots else


  std::cout<<"closing output root file"<<std::endl<<std::endl;
  rootfout->Close();
  return 0;
}


int main(int argc, char*argv[]){
  int rStatus=-1;
  if(argc!=3){
    std::cout<<"no defaults. Do...."<<std::endl;
    std::cout<<"./printPlots_MCJEC.exe "<<
      "<target_ppMC_dir> <outputNameBase>"<<std::endl<<std::endl;
    return rStatus;  }
  
  rStatus=1;
  rStatus=printPlots_MCJEC((const std::string)argv[1], (const std::string)argv[2] );
  return rStatus;
}


//// ALT METRICS
// mean  = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetMean()     :  0.5*(fgaus_gen30->GetParameter(1) + hrsp_gen30[nj][ip]->GetMean());
// mean  = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetMean()     : fgaus_gen30->GetParameter(1);
// emean = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetMeanError(): hrsp_gen30[nj][ip]->GetMeanError(1);
// emean = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetMeanError(): sqrt(pow(fgaus_gen30->GetParError(1),2) + pow(hrsp_gen30[nj][ip]->GetMeanError(),2));      
// sig   = (fitstatus!=0) ? hrsp_gen30[nj][ip]->GetRMS()/mean : fgaus_gen30->GetParameter(2)/fgaus_gen30->GetParameter(1);
// esig  = (fitstatus!=0) 
// ? sqrt((pow(1/mean,2)*pow(hrsp_gen30[nj][ip]->GetRMSError(),2))+(pow(-hrsp_gen30[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
// : sqrt((pow(1/mean,2)*pow(hrsp_gen30[nj][ip]->GetRMSError(),2))+(pow(-hrsp_gen30[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
// esig  = (fitstatus!=0) 
// ? sqrt((pow(1/mean,2)*pow(hrsp_gen30[nj][ip]->GetRMSError(),2))+(pow(-hrsp_gen30[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
// : sqrt((pow(1/mean,2)*pow(hrsp_gen30[nj][ip]->GetRMSError(),2))+(pow(-hrsp_gen30[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
//       
      
//hrsp_gen30[nj][ip]->Write();
