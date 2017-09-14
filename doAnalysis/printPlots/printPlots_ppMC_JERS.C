#include "printPlots.h"

const bool debugMode=true;

const bool draw_hJER=true;



//other options
const bool draw_MCEff=true;
const bool draw_hJERRapBins=false, doGenBinsToo=false;//RapBins -> dual-diff xsec bins, GenBins -> variable, depends on readForests
const bool draw_JERgen150to200=false, draw_JERgen30to50=false;

//const std::string defOutFilename="printPlots_MCJEC_ppMC_defOut";

int printPlots_ppMC_JERS(std::string inFile_MC_dir,const std::string outputTag){
  
  // root style settings.
  LoadStyle();  
    
  //figure out what radius/jetcollection we are looking at using the ppData filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+strlen("_ak");
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt+"_";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;



  //make final jetType strings
  std::size_t jetTypePos=radPos+1;
  std::size_t jetsPos=inFile_MC_dir.find("Jets");;  
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  
  const std::string fullJetType="ak"+radiusInt+jetType;//"ak3PFJets";//find in dir name
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  std::string didJetID="1";
  std::string doJetID="1";
  
  
  
  // input file
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";  
  std::cout<<std::endl<<"input file dir "<< inputDir+inFile_MC_dir <<std::endl;
  std::cout<<"input file name in dir "<< inFile_MC_name <<std::endl;
  TFile *finPP=new TFile( (inputDir+inFile_MC_dir+inFile_MC_name).c_str() );
  

  // output file
  std::string outFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+".root";
  std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
  TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE");   
  rootfout->cd();    


  // fit + fit settings
  std::cout<<"setting up fitters..."<< std::endl;
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  TF1 *fgaus=NULL;   
  
  // input hists
  if(draw_hJER){

    // array of pointers for input/output hists
    TH1F *hrsp[Nrad][nbins_pt_debug]={};  //input
    TH1F *hMean[Nrad]={}, *hSigma[Nrad]={}; //output    
    
    // loop over available radii (one radii for now)
    for(int nj=0; nj<Nrad; nj++){
      
      const std::string algname=algo+srad[nj]+jetType; 
      std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;

      // initialize mean hist    
      hMean[nj] = new TH1F(Form( "hMean%d", nj),    
			   Form( "Mean %s %s", algname.c_str(), ccent[nj]),    
			   nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);    

      MakeHistMean( (TH1F*)hMean[nj],ptbins_debug[0],ptbins_debug[nbins_pt_debug]);     
      
      // initialize RMS/Sigma hist    
      hSigma[nj] = new TH1F(Form( "hSigma%d", nj),    
			    Form( "Sigma %s %s", algname.c_str(), ccent[nj]),    
			    nbins_pt_debug, ptbins_debug); //nbins_pt, ptbins_bound);    
      MakeHistRMS( (TH1F*) hSigma[nj], ptbins_debug[0],ptbins_debug[nbins_pt_debug]);     
      
      // pt bin hJER-fit loop    
      if(debugMode)std::cout<<"nbins_pt="<<nbins_pt_debug<<std::endl;    
      for(int ip=0; ip<nbins_pt_debug; ip++){    
	
	int ptbin_ip=(int)ptbins_debug[ip], ptbin_ip1=(int)ptbins_debug[ip+1];
	if(debugMode)std::cout<<"pt range for bin: "<<ptbin_ip<<" - "<<ptbin_ip1<< " GeV "<<std::endl;    
	
	// open the input hist    
	std::string inputHistName="hJER_"+doJetID+"wJetID_ptbin"+std::to_string(ip);
	hrsp[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );        
	if(!hrsp[nj][ip]){ 
	  std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;
	  assert(false);}          
	else if(hrsp[nj][ip]->GetEntries()<3)
	  continue;
	
	if(debugMode)std::cout<<std::endl;  
	if(debugMode)hrsp[nj][ip]->Print("base");    

	hrsp[nj][ip]->Scale( 1./ hrsp[nj][ip]->Integral() );    
	
	// get some quick stats    
	//float norm  = hrsp[nj][ip]->GetMaximumStored();    
	float mean  = hrsp[nj][ip]->GetMean();    
	float emean = hrsp[nj][ip]->GetMeanError();    
	float sig   = hrsp[nj][ip]->GetStdDev()/mean;    
	float esig  = sig* 
	  sqrt (  pow( emean/mean,2) + pow( hrsp[nj][ip]->GetStdDevError()/hrsp[nj][ip]->GetStdDev(),2) ) ; 
	
	fgaus = new TF1("fgaus","gaus", 0.80,1.20);    
	//fgaus->SetParameters(norm, 1.0, 0.1);    
	//fgaus->SetParLimits(1,0.60,1.30);    
	
	
	std::cout<<"inputHist "<<inputHistName<<std::endl;    
	if(debugMode)std::cout<< "Mean= "<< hrsp[nj][ip]->GetMean()<<" +/- " <<hrsp[nj][ip]->GetMeanError()<< std::endl;    
	
	std::cout<<"fitting..."<<std::endl;    
	
	int fitstatus = 1;    	
	//fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");    
	fitstatus = hrsp[nj][ip]->Fit(fgaus,"R");    
	std::cout<< "Fit Status: "<< fitstatus<< std::endl;    
	
	mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()        : fgaus->GetParameter(1);    
	emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError()   : fgaus->GetParError(1) ; // hrsp[nj][ip]->GetMeanError(); //debug    
	sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/mean;    
	esig  = (fitstatus!=0) 
	  ? sig * sqrt( pow( (emean/mean),2) + 
			pow( ( hrsp[nj][ip]->GetStdDevError()/hrsp[nj][ip]->GetStdDev() ),2) 
			) 
	  : sig * sqrt( pow( (emean/mean), 2) + 
			pow( ( fgaus->GetParError(2)/fgaus->GetParameter(2) ), 2) 
			)   ;
	//? sqrt( (pow(1/mean,2) * pow(hrsp[nj][ip]->GetStdDevError(),2)) +     
	//	  (pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) )     
	//: sqrt( (pow(1/mean,2) * pow(hrsp[nj][ip]->GetStdDevError(),2)) +     
	//	  (pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) );     
	
	// set contents+errors in specific bin    
	std::cout<<"setting bin content for hMean[nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;    
	if(debugMode)std::cout<<"mean ="<< mean<<std::endl;    
	if(debugMode)std::cout<<"emean ="<< emean<<std::endl;    
	if(debugMode)std::cout<<"sig ="<< sig<<std::endl;    
	if(debugMode)std::cout<<"esig ="<< esig<<std::endl;    
	
	hMean[nj]->SetBinContent (ip+1, mean);    
	hMean[nj]->SetBinError   (ip+1, emean);    
	hSigma[nj]->SetBinContent (ip+1, sig);    
	hSigma[nj]->SetBinError   (ip+1, esig);    
	
      }// end fit-loop over ptbins    
     
    }// end loop over available jet radii    
    
    
    //if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;    
    //gPad->Close();    
    
    
    
    // DRAW THOSE PDFS //    
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+".pdf";    
    std::string open_thePDFFileName=thePDFFileName+"[";    
    std::string close_thePDFFileName=thePDFFileName+"]";    
    
    std::cout<<"opening output pdf file "<< thePDFFileName<<std::endl<<std::endl;    
    
    TCanvas* pdfoutCanv=new TCanvas("outputPdf","outputPdf", 800, 800);    
    pdfoutCanv->cd();    
    pdfoutCanv->Print( open_thePDFFileName.c_str() );    
    
    TCanvas* pdfoutCanv_wLogy=new TCanvas("outputPdfwLogx","output Pdf wLogx", 800, 800);    
    pdfoutCanv_wLogy->SetLogy(1);    
    pdfoutCanv_wLogy->cd();        

    
    
    // draw and print pdfs. //    
    std::cout<<std::endl<<"printing hists to pdf file..."<<std::endl;    
    for(int i=0;i<Nrad;++i){    
            
      pdfoutCanv_wLogy->cd();    
      for(int j=0;j<nbins_pt_debug;++j){    
	int ptbin_j=(int)ptbins_debug[j];
	int ptbin_j1=(int)ptbins_debug[j+1];
	std::string hrspTitle=std::to_string(ptbin_j)+" GeV < gen jet p_{T} < "+std::to_string(ptbin_j1)+" GeV";    
	hrsp[i][j]->SetTitle(hrspTitle.c_str());    
	hrsp[i][j]->SetTitleSize(0.02);    
	hrsp[i][j]->SetMarkerStyle(8);    
	hrsp[i][j]->SetMarkerSize(1.0);    
	hrsp[i][j]->GetXaxis()->SetTitle("recpt/genpt");    
	hrsp[i][j]->GetYaxis()->SetTitle("A.U.");    

	hrsp[i][j]->Draw("E1");    
	hrsp[i][j]->Write();    
	
	//float textx=0.5,texty=0.25;    
	//TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean[i][j])).c_str() );    
	//meanText->SetNDC();    
	//meanText->Draw();    
	//TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig[i][j])).c_str() );    
	//sigmaText->SetNDC();    
	//sigmaText->Draw();    
	//TLine* histMeanLine=new TLine( array_mean[i][j], 0., // x1,y1     
	//array_mean[i][j], 0.130 );// x2,y2    
	//array_mean[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2    
	//histMeanLine->Draw("same");    
	
	pdfoutCanv_wLogy->Print(thePDFFileName.c_str());    
	pdfoutCanv_wLogy->Write();
	
	
	
      }    
      
      TCanvas* pdfoutCanv_muSigma=new TCanvas("outputPdfwLog_muSigma","outputPdfwLog", 800, 800);    
      pdfoutCanv_muSigma->Divide(1,2);
      pdfoutCanv_muSigma->cd();
      
      TPad* p1=(TPad*)pdfoutCanv_muSigma->cd(1);        
      p1->SetLogx(1);    
      p1->cd();
      
      hMean[i]->Draw("HIST E1"); 
      //hMean[i]->Write();
      
      TLine* meanLine=new TLine(ptbins_debug[0],1.,ptbins_debug[nbins_pt_debug],1.);    
      meanLine->SetLineStyle(2);  
      meanLine->SetLineColor(kBlue);    
      meanLine->Draw();    
      
      
      
      
      TPad* p2=(TPad*)pdfoutCanv_muSigma->cd(2);        
      p2->SetLogx(1);    	
      p2->cd();    
      
      hSigma[i]->Draw("HIST E1");               

      //TVirtualFitter::SetDefaultFitter("TFitter");
      //std::string fitFunc_str="[0]+[1]*pow((x-"+std::to_string( ((int)ptbins_debug[0]) )+".),[2])";
      //std::cout<<"fitFunc_str="<<fitFunc_str<<std::endl;
      //TF1 *hSigmaFit=new TF1("hSigmaFit", fitFunc_str.c_str() ,

      TF1 *hSigmaFit=new TF1("hSigmaFit","[0]+[1]*pow(x,[2])", 
			     ptbins_debug[0],
			     ptbins_debug[nbins_pt_debug]);
      //hSigmaFit->SetParameter(0,1.0);
      //hSigmaFit->SetParameter(1,300. );
      //hSigmaFit->SetParameter(2,-1.5);
      
      float fitParam_0    =hSigmaFit->GetParameter(0);
      float fitParam_0_err=hSigmaFit->GetParError(0);
      float fitParam_1    =hSigmaFit->GetParameter(1);
      float fitParam_1_err=hSigmaFit->GetParError(1);
      float fitParam_2    =hSigmaFit->GetParameter(2);
      float fitParam_2_err=hSigmaFit->GetParError(2); 
      
      std::cout << "fitParam_0     = " << fitParam_0    << std::endl;
      std::cout << "fitParam_0_err = " << fitParam_0_err<< std::endl;
      std::cout << "%err fitParam0 = " << (fitParam_0_err/fitParam_0)*100.    << std::endl;
      std::cout << "fitParam_1     = " << fitParam_1    << std::endl;
      std::cout << "fitParam_1_err = " << fitParam_1_err<< std::endl;
      std::cout << "%err fitParam1 = " << (fitParam_1_err/fitParam_1)*100.    << std::endl;
      std::cout << "fitParam_2     = " << fitParam_2    << std::endl;
      std::cout << "fitParam_2_err =   " << fitParam_2_err<< std::endl;
      std::cout << "%err fitParam2 = " << (fitParam_2_err/fitParam_2)*100.    << std::endl;
      
      
      bool sigmaFitStatus= true;
      std::cout<<std::endl<<"fitting JER width v. pt... init. sigmaFitStatus="<<sigmaFitStatus<<std::endl<<std::endl; 
      //sigmaFitStatus=hSigma[i]->Fit(hSigmaFit,"MEVR");
      //sigmaFitStatus=hSigma[i]->Fit(hSigmaFit,"MR");
      sigmaFitStatus=hSigma[i]->Fit(hSigmaFit,"R");	
      fitParam_0    =hSigmaFit->GetParameter(0);
      fitParam_0_err=hSigmaFit->GetParError(0);
      fitParam_1    =hSigmaFit->GetParameter(1);
      fitParam_1_err=hSigmaFit->GetParError(1);
      fitParam_2    =hSigmaFit->GetParameter(2);
      fitParam_2_err=hSigmaFit->GetParError(2); 
      
      if(!sigmaFitStatus){	
	std::cout<<std::endl<<"fit success!, sigmaFitStatus="<<sigmaFitStatus<<std::endl<<std::endl; 
	std::cout << "fitParam_0     = " << fitParam_0    << std::endl;
	std::cout << "fitParam_0_err = " << fitParam_0_err<< std::endl;
	std::cout << "%err fitParam0 = " << (fitParam_0_err/fitParam_0)*100.    << std::endl;
	std::cout << "fitParam_1     = " << fitParam_1    << std::endl;
	std::cout << "fitParam_1_err = " << fitParam_1_err<< std::endl;
	std::cout << "%err fitParam1 = " << (fitParam_1_err/fitParam_1)*100.    << std::endl;
	std::cout << "fitParam_2     = " << fitParam_2    << std::endl;
	std::cout << "fitParam_2_err =   " << fitParam_2_err<< std::endl;
	std::cout << "%err fitParam2 = " << (fitParam_2_err/fitParam_2)*100.    << std::endl;
	//	assert(false);
      }
      else {
	std::cout<<std::endl<<"fit failed, sigmaFitStatus="<<sigmaFitStatus<<std::endl<<std::endl; 
	//	assert(false);
      }
      
      
      hSigmaFit->Draw("SAME");
      pdfoutCanv_muSigma->Print(thePDFFileName.c_str());        
      

      hSigma[i]->Write();               
      hSigmaFit->Write();

      pdfoutCanv_muSigma->Write();
      
      
      
      
      
      
      
      
      std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;                
      
      pdfoutCanv->Print( close_thePDFFileName.c_str() );    
      
      pdfoutCanv->Close();
      pdfoutCanv_muSigma->Close();
      pdfoutCanv_wLogy->Close();
      

    }              
    
  }
  else std::cout<<std::endl<<std::endl<<"skipping hJER hists"<<std::endl<<std::endl;
                       
  
  
  
  
  
  
  
  
  
  
  
  
  //////////////////////////////////////////////////
  ////////////// DUAL DIFF ETA BINS ////////////////
  //////////////////////////////////////////////////
  
  // input hists    
  if(!draw_hJERRapBins)std::cout<<std::endl<<std::endl<<"skipping hJER RapBins hists"<<std::endl<<std::endl;   
  else{ //draw h_JER    
  TH1F *hrsp[Nrad][nbins_abseta]={};      
  //double array_mean[Nrad][nbins_abseta]={};    
  //double array_sig[Nrad][nbins_abseta]={};    
      
  // output hists    
  TH1F *hMean[Nrad]={}, *hSigma[Nrad]={};       
      
  // first loop, pp loop over available radii    
  for(int nj=0; nj<Nrad; nj++){    
      
    const std::string algname=algo+srad[nj]+jetType;    
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;    
      
      
    // initialize mean hist    
    hMean[nj] = new TH1F(Form( "hMean_abseta%d", nj),    
			 Form( "abseta bins %s %s", algname.c_str(), ccent[nj]),    
			 nbins_abseta, absetabins); //nbins_abseta, absetabins_bound);    
    hMean[nj]->SetMarkerColor(1);    
    hMean[nj]->SetMarkerStyle(20);    
    hMean[nj]->SetLineColor(1);    
    hMean[nj]->SetMarkerSize(1.0);    
    MakeHistMeanEta(hMean[nj],1.082,0.858);     
    hMean[nj]->SetAxisRange(0.9, 1.4, "Y");    
    //hMean[nj]->SetAxisRange(0. , 4.7, "X");    
    
    // initialize RMS/Sigma hist    
    hSigma[nj] = new TH1F(Form( "hSigma_abseta%d", nj),    
			  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),    
			  nbins_abseta, absetabins); //nbins_abseta, absetabins_bound);    
    hSigma[nj]->SetMarkerColor(1);    
    hSigma[nj]->SetMarkerStyle(20);    
    hSigma[nj]->SetLineColor(1);    
    hSigma[nj]->SetMarkerSize(1.0);    
    MakeHistRMSAbsEta(hSigma[nj],0.563,0.001);     
    hSigma[nj]->SetAxisRange(0.0, 0.3 , "Y");    
    //hSigma[nj]->SetAxisRange(0.0, 4.7 , "X");    
      
      
    // pt bin hJER-fit loop    
    std::cout<<"looping over nbins_abseta="<<nbins_abseta<<std::endl;    
    for(int ip=0; ip<nbins_abseta; ip++){    
      
      // input hist title string    
      if(debugMode)std::cout<<"for abseta range "<<absetabins[ip]<<" to "<<absetabins[ip+1]<<std::endl;    
      std::string inputHistName="hJER_"+doJetID+"wJetID_jty_bin"+std::to_string(ip);//+

      // open the input hist    
      hrsp[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );        
      if(!hrsp[nj][ip]){std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;assert(false);}          
      
      // print the details...    
      //std::cout<<std::endl;  hrsp[nj][ip]->Print("base");    
      float binSize=1.0;
      if(absetabins[ip]==3.0) binSize=0.2*2.;//3.0-3.2
      else if(absetabins[ip]==3.2) binSize=1.5*2.;//3.2-4.7
      hrsp[nj][ip]->Scale( 1./ (binSize * hrsp[nj][ip]->Integral()) );    
          
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
      fgaus = new TF1("fgaus","gaus", 0.80,1.20);    
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
      
      
      //array_mean[nj][ip]=mean ; //for a striaght line for plots later.    
      //array_emean[nj][ip]=emean;    
      //array_sig[nj][ip]=sig  ;      
      //array_esig[nj][ip]=esig ;           
    }// end fit-loop over absetabins    
  }// end loop over available jet radii    
      
      
  if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;    
  gPad->Close();    
      
  {         
    // DRAW THOSE PDFS //    
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+"_RapBins.pdf";    
    std::string open_thePDFFileName=thePDFFileName+"[";    
    std::string close_thePDFFileName=thePDFFileName+"]";    
        
    std::cout<<"opening output pdf file "<< thePDFFileName<<std::endl<<std::endl;    
    TCanvas* pdfoutCanv=new TCanvas("outputPdf","outputPdf", 800, 600);    
    pdfoutCanv->Print( open_thePDFFileName.c_str() );    
    pdfoutCanv->cd();    
  
    TCanvas* pdfoutCanv_wLogy=new TCanvas("outputPdfwLogx","output Pdf wLogx", 800, 800);    
    //pdfoutCanv_wLogx->Print( open_thePDFFileName.c_str() );    
    pdfoutCanv_wLogy->cd();    
    pdfoutCanv_wLogy->SetLogy(1);    
      
    TCanvas* pdfoutCanv_muSigma=new TCanvas("outputPdfwLog_muSigma","outputPdfwLog", 800, 800);    
    //pdfoutCanv_muSigma->cd();    
    pdfoutCanv_muSigma->Divide(1,2);
    TPad* p2=(TPad*)pdfoutCanv_muSigma->cd(2);        
    //p2->SetLogx(1);    
    TPad* p1=(TPad*)pdfoutCanv_muSigma->cd(1);        
    //p1->SetLogx(1);    
        
    // draw and print pdfs. //    
    std::cout<<std::endl<<"printing hists to pdf file..."<<std::endl;    
    for(int i=0;i<Nrad;++i){    
      p1->cd();      
      hMean[i] ->Draw("HIST E1"); //pdfoutCanv->Print(thePDFFileName.c_str());    

      //TLine* meanLine=new TLine(0.,1.,4.69,1.);    
      TLine* meanLine=new TLine(absetabins[0],1.,absetabins[nbins_abseta],1.);    
      meanLine->SetLineStyle(2);       meanLine->SetLineColor(kBlue);    
      meanLine->Draw();    
          
      p2->cd();
      hSigma[i]->Draw("HIST E1");     
          
      //TLine* sigmaLine=new TLine(15,0.1,500,0.1);    
      //sigmaLine->SetLineStyle(2);       sigmaLine->SetLineColor(kBlue);    
      //sigmaLine->Draw();    
          
      pdfoutCanv_muSigma->Print(thePDFFileName.c_str());        
      
      pdfoutCanv_wLogy->cd();    
      for(int j=0;j<nbins_abseta;++j){    
	std::string hrspTitle="hJER "+std::to_string((int)(10*absetabins[j]))+"< abseta <"+std::to_string((int)(10*absetabins[j+1]));    
	hrsp[i][j]->SetTitle(hrspTitle.c_str());    
	hrsp[i][j]->SetMarkerStyle(8);    
	hrsp[i][j]->SetMarkerSize(1.0);    
	hrsp[i][j]->GetXaxis()->SetTitle("recpt/genpt");    
	//hrsp[i][j]->SetAxisRange(0.,0.130,"Y");    
	hrsp[i][j]->Draw("E1");    
	
	//float textx=0.5,texty=0.25;    
	//TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean[i][j])).c_str() );    
	//meanText->SetNDC();    
	//meanText->Draw();    
	//TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig[i][j])).c_str() );    
	//sigmaText->SetNDC();    
	//////////sigmaText->Draw();    
	////////////TLine* histMeanLine=new TLine( array_mean[i][j], 0., // x1,y1     
	////////////array_mean[i][j], 0.130 );// x2,y2    
	//array_mean[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2    
	//histMeanLine->Draw("same");    
	
	pdfoutCanv_wLogy->Print(thePDFFileName.c_str());    
      }    
    }    
    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;    
    pdfoutCanv->Print( close_thePDFFileName.c_str() );    
    pdfoutCanv->Close();
    pdfoutCanv_muSigma->Close();
    pdfoutCanv_wLogy->Close();
    
  }              
      
//  // save output root file. //    
//  rootfout->cd();      
//  std::cout<<"writing output to root file..."<<std::endl;    
//  for(int i=0;i<Nrad;++i){    
//    for(int j=0;j<nbins_abseta;++j) hrsp[i][j]->Write();    
//    hMean[i]->Write();    
//    hSigma[i]->Write();              
//  }    
//  //std::cout<<"closing output root file"<<std::endl<<std::endl;    
//  //rootfout->Close();    
  
  }                       







  ///////////////////////////////////////////////////////////////////////
  ////////////// DUAL DIFF ETA BINS IN SELECT GENPT BINS ////////////////
  ///////////////////////////////////////////////////////////////////////

  // input hists    
  if(!draw_hJERRapBins||!doGenBinsToo)std::cout<<std::endl<<std::endl<<"skipping hJER RapBins hists"<<std::endl<<std::endl;   
  else{ //draw_hJERRapBins && doGenBinsToo

  TH1F *hrsp_genBin1[Nrad][nbins_abseta]={};      
  //double array_mean_genBin1[Nrad][nbins_abseta]={};    
  //double array_sig_genBin1[Nrad][nbins_abseta]={};    
      
  // output hists    
  TH1F *hMean_genBin1[Nrad]={}, *hSigma_genBin1[Nrad]={};       

  // first loop, pp loop over available radii    
  for(int nj=0; nj<Nrad; nj++){    
      
    const std::string algname=algo+srad[nj]+jetType;    
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;    
      
      
    // initialize mean hist    
    hMean_genBin1[nj] = new TH1F(Form( "hMean_genBin1_abseta%d", nj),    
			 Form( "12<genpt<15 %s %s", algname.c_str(), ccent[nj]),    
			 nbins_abseta, absetabins); //nbins_abseta, absetabins_bound);    
    //hMean_genBin1[nj]->SetTitle("12<genpt<15");
    hMean_genBin1[nj]->SetMarkerColor(1);    
    hMean_genBin1[nj]->SetMarkerStyle(20);    
    hMean_genBin1[nj]->SetLineColor(1);    
    hMean_genBin1[nj]->SetMarkerSize(1.0);    
    MakeHistMeanEta(hMean_genBin1[nj],1.082,0.858);     
    hMean_genBin1[nj]->SetAxisRange(0.9, 1.4, "Y");    
    //hMean_genBin1[nj]->SetAxisRange(0. , 4.7, "X");    
    
    // initialize RMS/Sigma hist    
    hSigma_genBin1[nj] = new TH1F(Form( "hSigma_genBin1_abseta%d", nj),    
			  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),    
			  nbins_abseta, absetabins); //nbins_abseta, absetabins_bound);    
    hSigma_genBin1[nj]->SetMarkerColor(1);    
    hSigma_genBin1[nj]->SetMarkerStyle(20);    
    hSigma_genBin1[nj]->SetLineColor(1);    
    hSigma_genBin1[nj]->SetMarkerSize(1.0);    
    MakeHistRMSAbsEta(hSigma_genBin1[nj],0.563,0.001);     
    hSigma_genBin1[nj]->SetAxisRange(0.0, 0.3 , "Y");    
    //hSigma_genBin1[nj]->SetAxisRange(0.0, 4.7 , "X");    




    // pt bin hJER-fit loop    
    std::cout<<"looping over genBin1 with nbins_abseta="<<nbins_abseta<<std::endl;    
    for(int ip=0; ip<nbins_abseta; ip++){    
      
      // input hist title string    
      if(debugMode)std::cout<<"for abseta range "<<absetabins[ip]<<" to "<<absetabins[ip+1]<<std::endl;    

     std::string inputHistName="hJER_"+doJetID+"wJetID_jty_etabin"+std::to_string(ip)+"_genptBin1";//+

      // open the input hist    
      hrsp_genBin1[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );        
      if(!hrsp_genBin1[nj][ip]){std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;assert(false);}          
      
      // print the details...    
      //std::cout<<std::endl;  hrsp_genBin1[nj][ip]->Print("base");    
      float binSize=1.0;
      if(absetabins[ip]==3.0) binSize=0.2*2.;//3.0-3.2
      else if(absetabins[ip]==3.2) binSize=1.5*2.;//3.2-4.7
      hrsp_genBin1[nj][ip]->Scale( 1./ (binSize * hrsp_genBin1[nj][ip]->Integral()) );    
          
      // get some quick stats    
      double norm  = hrsp_genBin1[nj][ip]->GetMaximumStored();    
      double err   = hrsp_genBin1[nj][ip]->GetStdDevError();    
      double mean  = hrsp_genBin1[nj][ip]->GetMean();    
      double emean = hrsp_genBin1[nj][ip]->GetMeanError();    
      double sig   = hrsp_genBin1[nj][ip]->GetStdDev()/mean;    
      double esig  = pow(1/mean,2)*pow(err,2) + pow( -sig/pow(mean,2), 2)*pow(emean, 2);    
      int fitstatus = 0;    
          
      //fgaus = new TF1("fgaus","gaus", mean - 2.15*hrsp_genBin1[nj][ip]->GetRMS(), mean + 2.15*hrsp_genBin1[nj][ip]->GetRMS());    
      //fgaus = new TF1("fgaus","gaus", 0.80, 1.20);    
      fgaus = new TF1("fgaus","gaus", 0.80,1.20);    
      fgaus->SetParameters(norm, 0.9999, 0.5);    
      //fgaus->SetParameters(norm, 1.00, 0.6);    
      //fgaus->SetParameters(norm, 0.9999, 0.5);    
      //fgaus->SetParLimits(1,0.60,1.30);    
      //fgaus->SetParLimits(2,0.20,0.80);    
      //fgaus->SetParLimits(1,0.92,1.05);    
      //fgaus->SetParLimits(2,0.20,0.80);    
          
      std::cout<<"for inputHist "<<inputHistName<<", we have.."<<std::endl;    
      std::cout<< "Mean= "<< hrsp_genBin1[nj][ip]->GetMean()<<" +/- " <<hrsp_genBin1[nj][ip]->GetMeanError()<< std::endl;    
      
      std::cout<<"fitting..."<<std::endl;    
      fitstatus = hrsp_genBin1[nj][ip]->Fit(fgaus,"RQ");    
      std::cout<< "Fit Status: "<< fitstatus<< ", Fit Error: "<< fgaus->GetParError(1)<< std::endl;    
      
      std::cout<<"grabbing metrics (whatever this means)...."<<std::endl;    
      mean  = (fitstatus!=0) ? hrsp_genBin1[nj][ip]->GetMean()     : fgaus->GetParameter(1);    
      emean = (fitstatus!=0) ? hrsp_genBin1[nj][ip]->GetMeanError(): fgaus->GetParError(1) ; // hrsp_genBin1[nj][ip]->GetMeanError(); //debug    
      sig   = (fitstatus!=0) ? hrsp_genBin1[nj][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);    
      esig  = (fitstatus!=0)     
	? sqrt( (pow(1/mean,2) * pow(hrsp_genBin1[nj][ip]->GetStdDevError(),2)) +     
		(pow(-hrsp_genBin1[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) )     
	: sqrt( (pow(1/mean,2) * pow(hrsp_genBin1[nj][ip]->GetStdDevError(),2)) +     
		(pow(-hrsp_genBin1[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) );     
      
      // set contents+errors in specific bin    
      if(debugMode)std::cout<<"setting bin content for hMean_genBin1[nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;    
      if(debugMode)std::cout<<"mean ="<< mean<<std::endl;    
      if(debugMode)std::cout<<"emean ="<< emean<<std::endl;    
      hMean_genBin1[nj]->SetBinContent (ip+1, mean);    
      hMean_genBin1[nj]->SetBinError   (ip+1, emean);    
      // hMean_genBin1[nj]->SetBinContent (ip+1, roundoff(mean));    
      // hMean_genBin1[nj]->SetBinError   (ip+1, roundoff(emean));    
      
      if(debugMode)std::cout<<"sig ="<< sig<<std::endl;    
      if(debugMode)std::cout<<"esig ="<< esig<<std::endl;    
      hSigma_genBin1[nj]->SetBinContent (ip+1, sig);    
      hSigma_genBin1[nj]->SetBinError   (ip+1, esig);    
      
      
      //array_mean_genBin1[nj][ip]=mean ; //for a striaght line for plots later.    
      //array_emean_genBin1[nj][ip]=emean;    
      //array_sig_genBin1[nj][ip]=sig  ;      
      //array_esig_genBin1[nj][ip]=esig ;           
    }// end fit-loop over absetabins    
  }// end loop over available jet radii          
  
  // initialize mean hist    
  TH1F *hrsp_genBin2[Nrad][nbins_abseta]={};      
  //double array_mean_genBin2[Nrad][nbins_abseta]={};    
  //double array_sig_genBin2[Nrad][nbins_abseta]={};    
  
  // output hists    
  TH1F *hMean_genBin2[Nrad]={}, *hSigma_genBin2[Nrad]={};       
  
      
  // first loop, pp loop over available radii    
  for(int nj=0; nj<Nrad; nj++){    
      
    const std::string algname=algo+srad[nj]+jetType;    
    std::cout<<"making hists for algorithm "<< algname<<"..."<<std::endl;    
      
    // initialize mean hist    
    hMean_genBin2[nj] = new TH1F(Form( "hMean_genBin2_abseta%d", nj),    
			 Form( "84<genpt<97 %s %s", algname.c_str(), ccent[nj]),    
			 nbins_abseta, absetabins); //nbins_abseta, absetabins_bound);    
    //hMean_genBin2[nj]->SetTitle("84<genpt<97");
    hMean_genBin2[nj]->SetMarkerColor(1);    
    hMean_genBin2[nj]->SetMarkerStyle(20);    
    hMean_genBin2[nj]->SetLineColor(1);    
    hMean_genBin2[nj]->SetMarkerSize(1.0);    
    MakeHistMeanEta(hMean_genBin2[nj],1.082,0.858);     
    hMean_genBin2[nj]->SetAxisRange(0.9, 1.4, "Y");    
    //hMean_genBin2[nj]->SetAxisRange(0. , 4.7, "X");    
    
    // initialize RMS/Sigma hist    
    hSigma_genBin2[nj] = new TH1F(Form( "hSigma_genBin2_abseta%d", nj),    
			  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),    
			  nbins_abseta, absetabins); //nbins_abseta, absetabins_bound);    
    hSigma_genBin2[nj]->SetMarkerColor(1);    
    hSigma_genBin2[nj]->SetMarkerStyle(20);    
    hSigma_genBin2[nj]->SetLineColor(1);    
    hSigma_genBin2[nj]->SetMarkerSize(1.0);    
    MakeHistRMSAbsEta(hSigma_genBin2[nj],0.563,0.001);     
    hSigma_genBin2[nj]->SetAxisRange(0.0, 0.3 , "Y");    
    //hSigma_genBin2[nj]->SetAxisRange(0.0, 4.7 , "X");    




    // pt bin hJER-fit loop    
    std::cout<<"looping over genBin2 with nbins_abseta="<<nbins_abseta<<std::endl;    
    for(int ip=0; ip<nbins_abseta; ip++){    
      
      // input hist title string    
      if(debugMode)std::cout<<"for abseta range "<<absetabins[ip]<<" to "<<absetabins[ip+1]<<std::endl;          
      std::string inputHistName="hJER_"+doJetID+"wJetID_jty_etabin"+std::to_string(ip)+"_genptBin2";

      // open the input hist    
      hrsp_genBin2[nj][ip] = (TH1F*)finPP->Get( inputHistName.c_str() );        
      if(!hrsp_genBin2[nj][ip]){std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;assert(false);}          
      //if(hrsp_genBin2[nj][ip]->GetEntries()<100){std::cout<<"no input hist named " <<  inputHistName<< ", exiting..."<<std::endl;assert(false);}          
      
      // print the details...    
      //std::cout<<std::endl;  hrsp_genBin2[nj][ip]->Print("base");    
      //float binSize=1.0;
      //if(absetabins[ip]==3.0) binSize=0.2*2.;//3.0-3.2
      //else if(absetabins[ip]==3.2) binSize=1.5*2.;//3.2-4.7
      hrsp_genBin2[nj][ip]->Scale( 1./ ( hrsp_genBin2[nj][ip]->Integral()) );    
      
      // get some quick stats    
      //double norm  = hrsp_genBin2[nj][ip]->GetMaximumStored();    
      double err   = hrsp_genBin2[nj][ip]->GetStdDevError();    
      double mean  = hrsp_genBin2[nj][ip]->GetMean();    
      double emean = hrsp_genBin2[nj][ip]->GetMeanError();    
      double sig   = hrsp_genBin2[nj][ip]->GetStdDev()/mean;    
      double esig  = pow(1/mean,2)*pow(err,2) + pow( -sig/pow(mean,2), 2)*pow(emean, 2);    
      int fitstatus = 0;    
          
      //fgaus = new TF1("fgaus","gaus", mean - 2.15*hrsp_genBin2[nj][ip]->GetRMS(), mean + 2.15*hrsp_genBin2[nj][ip]->GetRMS());    
      //fgaus = new TF1("fgaus","gaus", 0.80, 1.20);    
      //if(ip==0)fgaus->Clear();
      fgaus = new TF1("fgaus","gaus", 0.80,1.20);    
      //fgaus->SetParameters(norm, 1.05, 0.05);    
      //fgaus->SetParameters(norm, 1.00, 0.6);    
      //fgaus->SetParameters(norm, 0.9999, 0.5);    
      //fgaus->SetParLimits(1,0.60,1.30);    
      //fgaus->SetParLimits(2,0.20,0.80);    
      //fgaus->SetParLimits(1,0.92,1.05);    
      //fgaus->SetParLimits(2,0.20,0.80);    
          
      std::cout<<"for inputHist "<<inputHistName<<", we have.."<<std::endl;    
      std::cout<< "Mean= "<< hrsp_genBin2[nj][ip]->GetMean()<<" +/- " <<hrsp_genBin2[nj][ip]->GetMeanError()<< std::endl;    
      
      std::cout<<"fitting..."<<std::endl;    
      std::cout<< "Fit Status: "<< fitstatus<< ", Fit Error: "<< fgaus->GetParError(1)<< std::endl;    

      //if(fitstatus!=0)continue;
      std::cout<<"grabbing metrics (whatever this means)...."<<std::endl;    
      mean  = (fitstatus!=0) ? hrsp_genBin2[nj][ip]->GetMean()     : fgaus->GetParameter(1);    
      emean = (fitstatus!=0) ? hrsp_genBin2[nj][ip]->GetMeanError(): fgaus->GetParError(1) ; // hrsp_genBin2[nj][ip]->GetMeanError(); //debug    
      sig   = (fitstatus!=0) ? hrsp_genBin2[nj][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);    
      esig  = (fitstatus!=0)     
	? sqrt( (pow(1/mean,2) * pow(hrsp_genBin2[nj][ip]->GetStdDevError(),2)) +     
		(pow(-hrsp_genBin2[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) )     
	: sqrt( (pow(1/mean,2) * pow(hrsp_genBin2[nj][ip]->GetStdDevError(),2)) +     
		(pow(-hrsp_genBin2[nj][ip]->GetStdDev()/pow(mean,2), 2) * pow(emean,2)) );     
      
      if(fitstatus!=0){std::cout<<"fit did not converge on first try; continuing...."<<std::endl;continue;}


      // set contents+errors in specific bin    
      if(debugMode)std::cout<<"setting bin content for hMean_genBin2[nj][ip]=["<< nj<<"]["<< ip<<"] "<<std::endl;    
      if(debugMode)std::cout<<"mean ="<< mean<<std::endl;    
      if(debugMode)std::cout<<"emean ="<< emean<<std::endl;    
      hMean_genBin2[nj]->SetBinContent (ip+1, mean);    
      hMean_genBin2[nj]->SetBinError   (ip+1, emean);    
      // hMean_genBin2[nj]->SetBinContent (ip+1, roundoff(mean));    
      // hMean_genBin2[nj]->SetBinError   (ip+1, roundoff(emean));    
      
      if(debugMode)std::cout<<"sig ="<< sig<<std::endl;    
      if(debugMode)std::cout<<"esig ="<< esig<<std::endl;    
      hSigma_genBin2[nj]->SetBinContent (ip+1, sig);    
      hSigma_genBin2[nj]->SetBinError   (ip+1, esig);    
      
      
      //array_mean_genBin2[nj][ip]=mean ; //for a striaght line for plots later.    
      //array_emean_genBin2[nj][ip]=emean;    
      //array_sig_genBin2[nj][ip]=sig  ;      
      //array_esig_genBin2[nj][ip]=esig ;           
    }// end fit-loop over absetabins    
  }// end loop over available jet radii    
      








if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;    
  gPad->Close();    
      
  {         
    // DRAW THOSE PDFS //    
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+"_RapBins_genBin1.pdf";    
    std::string open_thePDFFileName=thePDFFileName+"[";    
    std::string close_thePDFFileName=thePDFFileName+"]";    
        
    std::cout<<"opening output pdf file "<< thePDFFileName<<std::endl<<std::endl;    
    TCanvas* pdfoutCanv=new TCanvas("outputPdf","outputPdf", 800, 600);    
    pdfoutCanv->Print( open_thePDFFileName.c_str() );    
    pdfoutCanv->cd();    
  
    TCanvas* pdfoutCanv_wLogy=new TCanvas("outputPdfwLogx","output Pdf wLogx", 800, 800);    
    pdfoutCanv_wLogy->cd();    
    pdfoutCanv_wLogy->SetLogy(1);    
      
    TCanvas* pdfoutCanv_muSigma=new TCanvas("outputPdfwLog_muSigma","outputPdfwLog", 800, 800);    
    pdfoutCanv_muSigma->Divide(1,2);

    TPad* p2=(TPad*)pdfoutCanv_muSigma->cd(2);        
    //p2->SetLogx(1);    
    TPad* p1=(TPad*)pdfoutCanv_muSigma->cd(1);        
    //p1->SetLogx(1);    
        
    // draw and print pdfs. //    
    std::cout<<std::endl<<"printing hists to pdf file..."<<std::endl;    
    for(int i=0;i<Nrad;++i){    
      p1->cd();      
      hMean_genBin1[i] ->Draw("HIST E1"); //pdfoutCanv->Print(thePDFFileName.c_str());    

      //TLine* meanLine=new TLine(0.,1.,4.69,1.);    
      TLine* meanLine=new TLine(absetabins[0],1.,absetabins[nbins_abseta],1.);    
      meanLine->SetLineStyle(2);       meanLine->SetLineColor(kBlue);    
      meanLine->Draw();    
          
      p2->cd();
      hSigma_genBin1[i]->Draw("HIST E1");     
          
      //TLine* sigmaLine=new TLine(15,0.1,500,0.1);    
      //sigmaLine->SetLineStyle(2);       sigmaLine->SetLineColor(kBlue);    
      //sigmaLine->Draw();    
          
      pdfoutCanv_muSigma->Print(thePDFFileName.c_str());        
      
      pdfoutCanv_wLogy->cd();    
      for(int j=0;j<nbins_abseta;++j){    
	std::string hrspTitle="hJER "+std::to_string((int)(10*absetabins[j]))+"< abseta <"+std::to_string((int)(10*absetabins[j+1]));    
	hrsp_genBin1[i][j]->SetTitle(hrspTitle.c_str());    
	hrsp_genBin1[i][j]->SetMarkerStyle(8);    
	hrsp_genBin1[i][j]->SetMarkerSize(1.0);    
	hrsp_genBin1[i][j]->GetXaxis()->SetTitle("recpt/genpt");    
	//hrsp_genBin1[i][j]->SetAxisRange(0.,0.130,"Y");    
	hrsp_genBin1[i][j]->Draw("E1");    
	
	//float textx=0.5,texty=0.25;    
	//TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean_genBin1[i][j])).c_str() );    
	//meanText->SetNDC();    
	//meanText->Draw();    
	//TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig_genBin1[i][j])).c_str() );    
	//sigmaText->SetNDC();    
	//sigmaText->Draw();    
	//TLine* histMeanLine=new TLine( array_mean[i][j], 0., // x1,y1     
	//array_mean[i][j], 0.130 );// x2,y2    
	//array_mean[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2    
	//histMeanLine->Draw("same");    
	
	pdfoutCanv_wLogy->Print(thePDFFileName.c_str());    
      }    
    }    
    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;    
    pdfoutCanv->Print( close_thePDFFileName.c_str() );    






    thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+"_RapBins_genBin2.pdf";    
    open_thePDFFileName=thePDFFileName+"[";    
    close_thePDFFileName=thePDFFileName+"]";    
    pdfoutCanv->Print( open_thePDFFileName.c_str() );


    // draw and print pdfs. //    
    std::cout<<std::endl<<"printing hists to pdf file..."<<std::endl;    
    for(int i=0;i<Nrad;++i){    
      p1->cd();      
      hMean_genBin2[i] ->Draw("HIST E1"); //pdfoutCanv->Print(thePDFFileName.c_str());    

      //TLine* meanLine=new TLine(0.,1.,4.69,1.);    
      TLine* meanLine=new TLine(absetabins[0],1.,absetabins[nbins_abseta],1.);    
      meanLine->SetLineStyle(2);       meanLine->SetLineColor(kBlue);    
      meanLine->Draw();    
          
      p2->cd();
      hSigma_genBin2[i]->Draw("HIST E1");     
          
      //TLine* sigmaLine=new TLine(15,0.1,500,0.1);    
      //sigmaLine->SetLineStyle(2);       sigmaLine->SetLineColor(kBlue);    
      //sigmaLine->Draw();    
          
      pdfoutCanv_muSigma->Print(thePDFFileName.c_str());        
      
      pdfoutCanv_wLogy->cd();    
      for(int j=0;j<nbins_abseta;++j){    
	std::string hrspTitle="hJER "+std::to_string((int)(10*absetabins[j]))+"< abseta <"+std::to_string((int)(10*absetabins[j+1]));    
	hrsp_genBin2[i][j]->SetTitle(hrspTitle.c_str());    
	hrsp_genBin2[i][j]->SetMarkerStyle(8);    
	hrsp_genBin2[i][j]->SetMarkerSize(1.0);    
	hrsp_genBin2[i][j]->GetXaxis()->SetTitle("recpt/genpt");    
	//hrsp_genBin2[i][j]->SetAxisRange(0.,0.130,"Y");    
	hrsp_genBin2[i][j]->Draw("E1");    
	
	//float textx=0.5,texty=0.25;    
	//TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean_genBin2[i][j])).c_str() );    
	//meanText->SetNDC();    
	//meanText->Draw();    
	//TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig_genBin2[i][j])).c_str() );    
	//sigmaText->SetNDC();    
	//sigmaText->Draw();    
	//TLine* histMeanLine=new TLine( array_mean[i][j], 0., // x1,y1     
	//array_mean[i][j], 0.130 );// x2,y2    
	//array_mean[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2    
	//histMeanLine->Draw("same");    
	
	pdfoutCanv_wLogy->Print(thePDFFileName.c_str());    
      }    
    }    

    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;    
    pdfoutCanv->Print( close_thePDFFileName.c_str() );    


    pdfoutCanv->Close();
    pdfoutCanv_muSigma->Close();
    pdfoutCanv_wLogy->Close();
  }              
  
//  // save output root file. //    
//  
//  rootfout->cd();      
//  std::cout<<"writing output to root file..."<<std::endl;    
//  for(int i=0;i<Nrad;++i){    
//    for(int j=0;j<nbins_abseta;++j) hrsp_genBin2[i][j]->Write();    
//    hMean_genBin2[i]->Write();    
//    hSigma_genBin2[i]->Write();              
//  }    
//  
//  //std::cout<<"closing output root file"<<std::endl<<std::endl;    
//  //rootfout->Close();    
  
  }                       
  




















































































































  
  
  ////     JER/JES PLOTS, 150 < GENPT < 200        ////
  
  
  if(!draw_JERgen150to200)std::cout<<std::endl<<std::endl<<"skipping JER gen150-200 hists"<<std::endl<<std::endl;
  else{
    
  // input hists
  TH1F *hrsp_gen150[Nrad][nbins_eta]={};  
  //double array_mean_gen150[Nrad][nbins_eta]={};
  //double array_sig_gen150[Nrad][nbins_eta]={};
  
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
			       Form( "150<genpt<50 , etabins, %s %s", algname.c_str(), ccent[nj]),
				nbins_eta, etabins); //nbins_pt, ptbins_bound);
    hMean_gen150[nj]->SetMarkerColor(1);
    hMean_gen150[nj]->SetMarkerStyle(20);
    hMean_gen150[nj]->SetLineColor(1);
    hMean_gen150[nj]->SetMarkerSize(1.0);
    MakeHistMeanEta(hMean_gen150[nj],1.082,0.858); 
    hMean_gen150[nj]->SetAxisRange(0.8,1.2, "Y");
    //hMean_gen150[nj]->SetAxisRange(-5.3,5.3, "X");
    
    // initialize RMS/Sigma hist
    hSigma_gen150[nj] = new TH1F(Form( "hSigma_gen150%d", nj),
				 Form( "Sigma %s %s", algname.c_str(), ccent[nj]),
				 nbins_eta, etabins); //nbins_pt, ptbins_bound);
    hSigma_gen150[nj]->SetMarkerColor(1);
    hSigma_gen150[nj]->SetMarkerStyle(20);
    hSigma_gen150[nj]->SetLineColor(1);
    hSigma_gen150[nj]->SetMarkerSize(1.0);
    MakeHistRMSEta(hSigma_gen150[nj],0.563,0.001); 
    hSigma_gen150[nj]->SetAxisRange(0.0,0.5, "Y");
    //hSigma_gen150[nj]->SetAxisRange(-5.3,5.3, "X");


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
      fgaus_gen150 = new TF1("fgaus_gen150","gaus", 0.80,1.20);
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


      //array_mean_gen150[nj][ip]=mean ; //for a striaght line for plots later.
      //array_emean[nj][ip]=emean;
      //array_sig_gen150[nj][ip]=sig  ;  
      //array_esig[nj][ip]=esig ;       
    }// end fit-loop over ptbins
  }// end loop over available jet radii
  if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;
  gPad->Close();

  // DRAW PDFS. //
  //std::string jobType="_MCJEC_JetID"+doJetID;
  //if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;
  { 
    std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+"_gen150to200.pdf";
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
      hMean_gen150[i] ->Draw("HIST E1"); //pdfoutCanv->Print(thePDFFileName.c_str());
      TLine* meanLine=new TLine(etabins[0],1., etabins[nbins_eta],1.);
      meanLine->SetLineStyle(2); meanLine->SetLineColor(kBlue);
      meanLine->Draw();

      pdfoutCanv->cd(2);
      hSigma_gen150[i]->Draw("HIST E1"); 
      //TLine* sigmaLine=new TLine(-5.3,0.1,5.3,0.1);
      //sigmaLine->Draw();
      pdfoutCanv->Print(thePDFFileName.c_str());    
      
      pdfoutCanv_wLog->cd();
      for(int j=0;j<nbins_eta;++j){
	if( !(hrsp_gen150[i][j]) )continue;
	std::string hrsp_gen150Title="JER, 150<genpt<200, "+std::to_string(etabins[j])+"<geneta<"+std::to_string(etabins[j+1]);
	hrsp_gen150[i][j]->SetTitle(hrsp_gen150Title.c_str());
	hrsp_gen150[i][j]->SetMarkerStyle(8);
	hrsp_gen150[i][j]->SetMarkerSize(1.0);
	hrsp_gen150[i][j]->GetXaxis()->SetTitle("recpt/genpt");
	//hrsp_gen150[i][j]->SetAxisRange(0.,0.130,"Y");
	hrsp_gen150[i][j]->Draw("E1");
	
	//float textx=0.5,texty=0.25;
	//TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean_gen150[i][j])).c_str() );
	//meanText->SetNDC();
	//meanText->Draw();
	//TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig_gen150[i][j])).c_str() );
	//sigmaText->SetNDC();
	//sigmaText->Draw();
	//TLine* histMeanLine=new TLine( array_mean_gen150[i][j], 0., // x1,y1 
	//array_mean_gen150[i][j], 0.130 );// x2,y2
	//array_mean_gen150[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2
	//histMeanLine->Draw("same");
	
	pdfoutCanv_wLog->Print(thePDFFileName.c_str());
      }
    }
    std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;
    pdfoutCanv->Print( close_thePDFFileName.c_str() );
    pdfoutCanv->Close();
    pdfoutCanv_wLog->Close();
    //pdfoutCanv_wLogx->Close();
  }

//  // save output root file. //
//
//  {
//
//    //std::string outFileName=outputDir+fullJetType+jobType+"_150to200_"+outputTag+".root";
//    //std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
//    //TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE"); 
//    rootfout->cd();  
//    std::cout<<"writing output to root file..."<<std::endl;
//    for(int i=0;i<Nrad;++i){
//      hMean_gen150[i]->Write();
//      hSigma_gen150[i]->Write();    
//      for(int j=0;j<nbins_eta;++j) {
//	if( !(hrsp_gen150[i][j]) ) continue;
//	else hrsp_gen150[i][j]->Write();
//      }
//    }
//    
//  }
  
  }

 
  ////     JER/JES PLOTS, 30 < GENPT < 50        ////
  if(!draw_JERgen30to50)std::cout<<std::endl<<std::endl<<"skipping JER gen30-50 hists"<<std::endl<<std::endl;
  else{
    
    // input hists
    TH1F *hrsp_gen30[Nrad][nbins_eta]={};  
    //double array_mean_gen30[Nrad][nbins_eta];
  //double array_sig_gen30[Nrad][nbins_eta];
    
    
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
				 Form( "30<genpt<50, etabins, %s %s", algname.c_str(), ccent[nj]),
				 nbins_eta, etabins); //nbins_pt, ptbins_bound);
      hMean_gen30[nj]->SetMarkerColor(1);
      hMean_gen30[nj]->SetMarkerStyle(20);
      hMean_gen30[nj]->SetLineColor(1);
      hMean_gen30[nj]->SetMarkerSize(1.0);
      MakeHistMeanEta(hMean_gen30[nj],1.082,0.858); 
      hMean_gen30[nj]->SetAxisRange(0.8,1.2, "Y");
      //hMean_gen30[nj]->SetAxisRange(-5.3,5.3, "X");
      
      // initialize RMS/Sigma hist
      hSigma_gen30[nj] = new TH1F(Form( "hSigma_gen30%d", nj),
				  Form( "Sigma %s %s", algname.c_str(), ccent[nj]),
				  nbins_eta, etabins); //nbins_pt, ptbins_bound);
      hSigma_gen30[nj]->SetMarkerColor(1);
      hSigma_gen30[nj]->SetMarkerStyle(20);
      hSigma_gen30[nj]->SetLineColor(1);
      hSigma_gen30[nj]->SetMarkerSize(1.0);
      MakeHistRMSEta(hSigma_gen30[nj],0.563,0.001); 
      hSigma_gen30[nj]->SetAxisRange(0.0,0.5, "Y");
      //hSigma_gen30[nj]->SetAxisRange(-5.3,5.3, "X");
      
      
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
	fgaus_gen30 = new TF1("fgaus_gen30","gaus", 0.80,1.20);
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
	
	
	//array_mean_gen30[nj][ip]=mean ; //for a striaght line for plots later.
	//array_emean[nj][ip]=emean;
	//array_sig_gen30[nj][ip]=sig  ;  
	//array_esig[nj][ip]=esig ;       
      }// end fit-loop over ptbins
    }// end loop over available jet radii
    if(debugMode)std::cout<<"closing gPad..."<<std::endl<<std::endl;
    gPad->Close();
    
    // DRAW PDFS. //
    //std::string jobType="_MCJEC_JetID"+doJetID;
    //if(debugMode)std::cout<<"jobType="<<jobType<<std::endl;
    { 
      std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+"_gen30to50.pdf";
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
	hMean_gen30[i]->Draw("HIST E1"); //pdfoutCanv->Print(thePDFFileName.c_str());
	
	TLine* meanLine=new TLine(etabins[0],1., etabins[nbins_eta],1.);     
	meanLine->SetLineStyle(2);      meanLine->SetLineColor(kBlue);
	meanLine->Draw();

	pdfoutCanv->cd(2);
	hSigma_gen30[i]->Draw("HIST E1"); 
	//TLine* sigmaLine=new TLine(-5.3,0.1,5.3,0.1);
	//sigmaLine->Draw();
	pdfoutCanv->Print(thePDFFileName.c_str());    
	
	pdfoutCanv_wLog->cd();
	for(int j=0;j<nbins_eta;++j){
	  if( !(hrsp_gen30[i][j]) )continue;
	  std::string hrsp_gen30Title="JER, 30<genpt<50, "+std::to_string(etabins[j])+"<geneta<"+std::to_string(etabins[j+1]);
	  hrsp_gen30[i][j]->SetTitle(hrsp_gen30Title.c_str());
	  hrsp_gen30[i][j]->SetMarkerStyle(8);
	  hrsp_gen30[i][j]->SetMarkerSize(1.0);
	  hrsp_gen30[i][j]->GetXaxis()->SetTitle("recpt/genpt");
	  //hrsp_gen30[i][j]->SetAxisRange(0.,0.130,"Y");
	  hrsp_gen30[i][j]->Draw("E1");
	  
	  //float textx=0.5,texty=0.25;
	  //TLatex* meanText=new TLatex(textx,texty, ("mean="+std::to_string(array_mean_gen30[i][j])).c_str() );
	  //meanText->SetNDC();
	  //meanText->Draw();
	  //TLatex* sigmaText=new TLatex(textx,texty-.05, ("sigma="+std::to_string(array_sig_gen30[i][j])).c_str() );
	  //sigmaText->SetNDC();
	  //sigmaText->Draw();
	  //TLine* histMeanLine=new TLine( array_mean_gen30[i][j], 0., // x1,y1 
	  //array_mean_gen30[i][j], 0.130 );// x2,y2
	  //array_mean_gen30[i][j], pdfoutCanv->GetFrame()->GetY2() );// x2,y2
	  //histMeanLine->Draw("same");
	  
	  pdfoutCanv_wLog->Print(thePDFFileName.c_str());
	}
      }
      std::cout<<"closing output pdf file "<<thePDFFileName<<std::endl<<std::endl;
      pdfoutCanv->Print( close_thePDFFileName.c_str() );
      pdfoutCanv->Close();
      pdfoutCanv_wLog->Close();
      //pdfoutCanv_wLogx->Close();
    }
    
//  // save output root file. //
//    
//    {
//      
//      //std::string outFileName=outputDir+fullJetType+jobType+"_30to50_"+outputTag+".root";
//      //std::cout<<std::endl<<"opening output root file "<< outFileName<<std::endl;
//      //TFile *rootfout = new TFile(outFileName.c_str(),"RECREATE"); 
//      rootfout->cd();  
//      std::cout<<"writing output to root file..."<<std::endl;
//      for(int i=0;i<Nrad;++i){
//	hMean_gen30[i]->Write();
//	hSigma_gen30[i]->Write();    
//	for(int j=0;j<nbins_eta;++j) {
//	  if( !(hrsp_gen30[i][j]) ) continue;
//	  else hrsp_gen30[i][j]->Write();
//	}
//      }
//    }
    
  }
  
  
  
  if(!draw_MCEff)  std::cout<<std::endl<<std::endl<<"skipping MCEff gen30-50 hists"<<std::endl<<std::endl;
  else{
    
    //gStyle->SetOptStat(0);
    //gROOT->ForceStyle();
    { 
      std::cout<<" drawing MC Eff. QA Plots..."<<std::endl;
      
      //std::string thePDFFileName=outputDir+fullJetType+jobType+"_MCEff_"+outputTag+".pdf";
      std::string thePDFFileName=outputDir+fullJetType+"_MCJEC_"+outputTag+"_MCEff.pdf";
      std::string open_thePDFFileName=thePDFFileName+"[";
      std::string close_thePDFFileName=thePDFFileName+"]";
      std::cout<<std::endl<<"creating temporary canvas for printing MCEff plots..."<<std::endl;
      
      TCanvas *temp_canvMCEff = new TCanvas("tempMCEff", "tempMCEff", 1200, 600);
      temp_canvMCEff->Print(open_thePDFFileName.c_str());
      temp_canvMCEff->SetLogz(1);
      temp_canvMCEff->cd();
      
      
      bool drawProfiles=false;
      if(drawProfiles){     
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
	  
	  theMCEffQAHist->Draw("E HIST COLZ"); 
	  
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
	  
	  theMCEffQAHist->Draw("HIST E COLZ"); 
	  
	  temp_canvMCEff->Print(thePDFFileName.c_str());   
	}//end genvar ptrat loops      
      }//end drawProfles
      
      
      
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
	//std::string stringRatio="(rec"+genVars[j]+"/gen"+genVars[j]+")";
	//std::string h_YAx_Title=stringRatio;
	
	//std::string h_XAx_Title="";    
	//if(j==0)h_XAx_Title="gen Jet Pt (GeV)"  ;
	//if(j==1)h_XAx_Title="gen Jet Eta"       ;
	//if(j==2)h_XAx_Title="gen Jet Phi (rad)" ;
	
	
	the2DMCEffQAHist->SetTitle (    h_Title.c_str() );
	//the2DMCEffQAHist->SetXTitle( h_XAx_Title.c_str() );
	//the2DMCEffQAHist->SetYTitle( h_YAx_Title.c_str() );
	
	//the2DMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
	//if(j==0)the2DMCEffQAHist->SetAxisRange(  0., 200., "X");//genpt
	//if(j==1)the2DMCEffQAHist->SetAxisRange(-3.0, 3.0 , "X");//eta
	//if(j==2)the2DMCEffQAHist->SetAxisRange(-4.0, 4.0 , "X");//phi
	
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
	
	//the2DMCEffQAHist->SetAxisRange(0.5,2.0,"Y");
	//if(j==0)the2DMCEffQAHist->SetAxisRange( -3.0, 3.0, "X");
	//if(j==1)the2DMCEffQAHist->SetAxisRange( -4.0, 4.0, "X");
	//if(j==2)the2DMCEffQAHist->SetAxisRange(   0., 0.5, "X");
	
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
  }
  
  if(debugMode)std::cout<<"closing output root file"<<std::endl<<std::endl;
  rootfout->Close();
  return 0;
}





int main(int argc, char*argv[]){
  int rStatus=-1;
  if(argc!=3){
    std::cout<<"no defaults. Do...."<<std::endl;
    std::cout<<"./printPlots_ppMC_JERS.exe "<<
      "<target_ppMC_dir> <output tag>"<<std::endl<<std::endl;
    return rStatus;  }
  
  rStatus=1;
  rStatus=printPlots_ppMC_JERS((const std::string)argv[1], (const std::string)argv[2] );
  return rStatus;
}


