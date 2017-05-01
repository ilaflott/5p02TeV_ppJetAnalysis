#include "unfoldSpectra.h"

// procedural settings
const bool doBayes=true; 
const int kIter = 4;//,kIterRange=4, kIterDraw = 3, kIterCenter=21;

const bool doSVD=true; //!(doBayes); 
const int  nKregMax = 9, kRegRange=4, kRegDraw = 4, kRegCenter=5;

const bool drawPDFs=true; 
const bool drawPDFs_BayesInputHistos= doBayes && drawPDFs;
const bool drawPDFs_SVDInputHistos  = doSVD   && drawPDFs;

//other settings
const bool doToyErrs=true;
const bool debugMode=true;
const bool debugPearson=(false && debugMode) ;
const bool useSimplePtBinning=false;//bin by ten everywhere instead of custom binning

//  the code --------------------------------------------------
int unfoldDataSpectra( std::string inFile_Data_dir , std::string inFile_MC_dir , 
		       const std::string baseName ,    const bool doJetID=true ){
  
  inFile_Data_dir=SCRATCH_BASE+inFile_Data_dir;
  inFile_MC_dir=SCRATCH_BASE+inFile_MC_dir;
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl; 

  //figure out what radius/jetcollection we are looking at using the MC filename
  std::size_t radPos=inFile_MC_dir.find("_ak")+3;
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt;
  const std::string RandEtaRange="_"+radius+"_20_eta_20";
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  
  const std::string Rstring="_"+radius; 
  const std::string Rstring_plotTitle=" "+radius; 
  const std::string RandEtaRange_plotTitle=Rstring_plotTitle+" 20eta20";

  std::size_t jetTypePos=radPos+1, jetsPos=inFile_MC_dir.find("Jets");
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;

  const std::string fullJetType="ak"+radiusInt+jetType;//"ak3PFJets";//find in dir name
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  const std::string inFile_Data_name="/HighPtJetTrig_"+fullJetType+"-allFiles.root";

  const std::string outFileName=unfoldSpectra_outdir+fullJetType+"_"+baseName;//+".root";  
  const std::string outRootFile=unfoldSpectra_outdir+fullJetType+"_"+baseName+".root";  
  const std::string outBayesPdfFile=baseName+"_Bayes";
  const std::string outSVDPdfFile=baseName+"_SVD"; // see drawPDFs part for rest of string



  // set error handing, stat info, other settings  // fix me
  RooUnfold::ErrorTreatment errorTreatment;
  if(!doToyErrs) errorTreatment = RooUnfold::kCovariance;
  else errorTreatment = RooUnfold::kCovToy; 
  if(debugMode)std::cout<<"doToyErrs="<<doToyErrs<<std::endl; 

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  //// script performance tracking if needed
  //TStopwatch timer; 
  //timer.Start();


  //  ppData input histos --------------------------------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from DATA file"<<std::endl;

  if(debugMode){
    //std::cout<<"BASE : "<< (CMSSW_BASE) <<std::endl; 
    std::cout<<"input data dir : "<< (inFile_Data_dir)  <<std::endl; 
    std::cout<<"data file name : "<< (inFile_Data_name)<<std::endl; 
    std::cout<<std::endl<<std::endl;  }

  //TFile *fpp_Data = TFile::Open( (CMSSW_BASE+inFile_Data_dir+inFile_Data_name).c_str());
  TFile *fpp_Data = TFile::Open( (inFile_Data_dir+inFile_Data_name).c_str());

  //reco jetpt spectra to unfold
  std::string histTitle="hJetQA_";
  if(doJetID)histTitle+="1wJetID_jtpt";
  else histTitle+="0wJetID_jtpt";
  TH1F *hrec, *hrec_anabin;
  hrec = (TH1F*)fpp_Data->Get( histTitle.c_str() );
  hrec->Print("base");

  std::string newHistTitle="hpp_anaBin_"+histTitle;
  hrec_anabin = (TH1F*)hrec->Clone( newHistTitle.c_str() );

  //std::cout<<std::endl<<"rebinning; i'll likely crash!!"<<std::endl<<std::endl;  
  //rebin/scaling and binwidth
  //hrec_anabin = (TH1F*)hrec_anabin->Rebin(10);
  //hrec_anabin->Scale(1./10.);
  //hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins, newHistTitle.c_str() , ptbins);

  hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins_pt, newHistTitle.c_str() , boundaries_pt);
  hrec_anabin->Print("base");    

  divideBinWidth(hrec_anabin);
  hrec_anabin->Print("base");    

  //assert(false);

  //intlumi efficiency/scaling
  TH1F *h_NEvents_vzCut   = (TH1F*)fpp_Data->Get("NEvents_vzCut");
  TH1F *h_NEvents         = (TH1F*)fpp_Data->Get("NEvents");
  long double effIntgrtdLumi_vz=integratedLuminosity, LumiEff_vz=1.;
  LumiEff_vz = h_NEvents_vzCut->GetEntries()/h_NEvents->GetEntries();
  effIntgrtdLumi_vz*=LumiEff_vz;
  std::cout<<std::endl<<"integrated luminosity of jet80+LowerJets datasets: ";
  std::cout<< std::endl << integratedLuminosity << std::endl;
  std::cout<<"luminosity efficiency: " << std::endl << LumiEff_vz << std::endl;
  std::cout<<"effective integrated luminosity: " << std::endl ;
  std::cout<< effIntgrtdLumi_vz << std::endl<< std::endl;

  //hrec_anabin = (TH1F*)hrec_anabin->Rebin(10);
  //hrec_anabin->Scale(1./10.);
  hrec_anabin->Scale(1/effIntgrtdLumi_vz);
  hrec_anabin->Print("base");    
  

  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  if(debugMode){
    //std::cout<<"BASE : "<< (CMSSW_BASE) <<std::endl; 
    std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
    std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl; 
    std::cout<<std::endl<<std::endl;  }
  //TFile *fpp_MC = TFile::Open( (CMSSW_BASE+inFile_MC_dir+inFile_MC_name).c_str());
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());

  //mat, NOT NORMALIZED
  std::string TH2_title="hpp_matrix"+RandEtaRange;
  std::string TH2_newTitle="hpp_anaBin_Trans_matrix_HLT"+RandEtaRange;
  TH2F *hmat,*hmat_anabin;
  hmat = (TH2F*)fpp_MC->Get( TH2_title.c_str() );
  hmat->Print("base");

  std::cout<<std::endl<<"trying to rebin a TH2; brace for impact!!!"<<std::endl<<std::endl;

  //hmat_anabin = (TH2F*)hmat->Clone( TH2_newTitle.c_str() );
  hmat_anabin=reBinTH2(hmat, TH2_newTitle);
  hmat_anabin->Print("base");    

  //assert(false);
  //hmat_anabin->Rebin2D(10, 10);
  
  //gen
  std::string genHistTitle="hpp_gen"+RandEtaRange;
  std::string newGenHistTitle="hpp_anaBin_gen"+RandEtaRange;//"hpp_anaBin_gen_JetComb"+RandEtaRange;
  TH1F*hgen,*hgen_anabin;
  hgen= (TH1F*)fpp_MC->Get( genHistTitle.c_str() );
  hgen->Print("base");

  hgen_anabin = (TH1F*)hgen->Clone( newGenHistTitle.c_str() );
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins_pt, newGenHistTitle.c_str() , boundaries_pt);
  hgen_anabin->Print("base");    

  divideBinWidth(hgen_anabin);
  hgen_anabin->Print("base");

  //assert(false);
  //hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins, ("hpp_anaBin_gen_JetComb"+RandEtaRange).c_str(), ptbins);
  //hgen_anabin = (TH1F*)hgen_anabin->Rebin(10);
  //hgen_anabin->Scale(1/10.);
  //hgen_anabin->Scale(hrec_anabin->Integral()/hgen_anabin->Integral());



  // output histos -------------------------
  std::cout<<std::endl<<std::endl<<"creating new OUTPUT histos..."<<std::endl;

  //gen
  TH1F *hgen_resp, *hgen_resp_anabin;
  hgen_resp = new TH1F( ("hpp_gen_response"+RandEtaRange).c_str(),"", 
			hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  hgen_resp->Sumw2();
  hgen_resp->Print("base");
  //  hgen_resp_anabin = new TH1F( ("hpp_gen_response_anabin"+RandEtaRange).c_str() ,"", nbins, ptbins);
  hgen_resp_anabin = new TH1F( ("hpp_gen_response_anabin"+RandEtaRange).c_str() ,"", nbins_pt, boundaries_pt);
  hgen_resp_anabin->Sumw2();
  hgen_resp_anabin->Print("base");
  
  //rec
  if(debugMode)std::cout<<"opening output rec/recresponse histograms..."<<std::endl;
  TH1F *hrec_resp, *hrec_resp_anabin; 
  hrec_resp = new TH1F( ("hpp_rec_response"+RandEtaRange).c_str(),"", 
		       hrec->GetNbinsX(), hrec->GetXaxis()->GetXmin(), hrec->GetXaxis()->GetXmax());
  hrec_resp->Sumw2();
  hrec_resp->Print("base");
  hrec_resp_anabin = new TH1F( ("hpp_rec_response_anabin"+RandEtaRange).c_str(),"", nbins_pt, boundaries_pt);//nbins, ptbins);
  hrec_resp_anabin->Sumw2();
  hrec_resp_anabin->Print("base");
  
  // open output file for writing output to a root file as we go along
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE"); 

  
  // Bayesian unfolding ------------------------- 
  if(!doBayes)  std::cout<<std::endl<<"   skipping Bayesian Unfolding..."<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<std::endl<<"   beginning Bayesian Unfolding..."<<std::endl;

    std::cout<<"calling RooUnfoldResponse "<<std::endl;
    RooUnfoldResponse roo_resp( hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix"+RandEtaRange).c_str());

    std::cout<<"calling RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );

    TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);
    hunf->SetName("ppData_BayesianUnfolding_UnfoldSpectra");
    hunf->SetTitle("ppData, Bayesian Unfolded Spectra");
    std::cout<<std::endl; hunf->Print("base");

    TH1F *hratio = (TH1F*)hunf->Clone( "ppData_BayesianUnfolding_UnfoldRatio" );
    hratio->SetTitle( "ppData, Bayesian Unfolded/ppMC Gen" );
    hratio->SetMarkerStyle(24);
    hratio->SetMarkerColor(kRed);
    hratio->Divide(hgen_anabin);
    hratio->Print("base");

    std::cout<<std::endl<<"writing bayesian unfolding output to file..."<<std::endl;
    fout->cd();

    hunf->SetMarkerStyle(24);
    hunf->SetMarkerColor(kRed);
    hunf->Write();

    hgen_resp->Write();      hgen_resp_anabin->Write();
    hrec_resp->Write();	     hrec_resp_anabin->Write();

    hratio->SetMarkerStyle(24);    
    hratio->SetMarkerColor(kRed);  
    hratio->Write();               

    std::cout<<"writing input histos to output file for easy access later..."<<std::endl;
    hrec->Write(); hrec_anabin->Write(); 
    hmat->Write(); hmat_anabin->Write(); 
    hgen->Write(); hgen_anabin->Write(); 

    //  drawPDFS -------------------------------------------------- 
    if(!drawPDFs)std::cout<<std::endl<<"NOT drawing PDFs for Bayesian Unfolding!"<<std::endl<<std::endl;
    else{ std::cout<<std::endl<<"drawing output PDFs for Bayesian Unfolding..."<<std::endl;
      std::string outPdfFile=unfoldSpectra_outdir+outBayesPdfFile+".pdf";
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";

      // temp canvas for printing
      TCanvas* tempCanvForPdfPrint_wLogy=new TCanvas("tempCanv_Logy","",1000,800);    
      tempCanvForPdfPrint_wLogy->cd();
      tempCanvForPdfPrint_wLogy->SetLogy(1);

      TCanvas* tempCanvForPdfPrint_wLogz=new TCanvas("tempCanv_Logz","",1200,800);    
      tempCanvForPdfPrint_wLogz->cd();
      tempCanvForPdfPrint_wLogz->SetLogz(1);
      
      TCanvas* tempCanvForPdfPrint=new TCanvas("tempCanv","",1000,800);    
      tempCanvForPdfPrint->cd();
      tempCanvForPdfPrint->SetLogy(0);
      tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); 
      
      if(drawPDFs_BayesInputHistos){

	//---------------

	tempCanvForPdfPrint_wLogz->cd();

	std::cout<<std::endl<<"drawing input histos to Bayesian Unfolding..."<<std::endl;
	hmat_anabin->SetAxisRange(0.0000000000001,.001,"Z");
        hmat_anabin->SetTitle("ppMC jet input");
        hmat_anabin->GetZaxis()->SetLabelSize(0.028);

        hmat_anabin->GetYaxis()->SetLabelSize(0.02);
        hmat_anabin->GetYaxis()->SetTitleSize(0.03);
        hmat_anabin->GetYaxis()->SetTitle("gen p_{t}");

        hmat_anabin->GetXaxis()->SetLabelSize(0.02);
        hmat_anabin->GetXaxis()->SetTitleSize(0.03);
        hmat_anabin->GetXaxis()->SetTitle("reco p_{t}   ");
	hmat_anabin->Draw("COLZ");           

	tempCanvForPdfPrint_wLogz->Print(outPdfFile.c_str());

	//---------------

	tempCanvForPdfPrint_wLogy->cd();

        hrec_anabin->SetMarkerStyle(kOpenTriangleUp);
        hrec_anabin->SetMarkerColor(kBlue);     
        hrec_anabin->SetMarkerSize(0.90);     
        hrec_anabin->SetTitle("input hists and unfolded data");
        hrec_anabin->Draw();           

        hgen_anabin->SetMarkerStyle(kOpenSquare);
        hgen_anabin->SetMarkerColor(kGreen+2);
        hgen_anabin->SetMarkerSize(0.90);     
	hgen_anabin->Draw("SAME");           

	hunf->SetMarkerStyle(kCircle);
        hunf->SetMarkerColor(kRed);
        hunf->SetMarkerSize(0.92);     
	hunf->Draw("SAME");           

	
	TLegend* legend = new TLegend( 0.7,0.8,0.9,0.9 );
	legend->AddEntry(hgen_anabin, "MC gen jet", "p");
	legend->AddEntry(hrec_anabin, "data reco jet", "p");
	legend->AddEntry(hunf, "unfolded data", "p");
	legend->Draw();

	tempCanvForPdfPrint_wLogy->Print(outPdfFile.c_str());

	//---------------
	

      }

      // draw and print unfold/ratio plots
      //tempCanvForPdfPrint_wLogy->cd();
      //hunf->Draw();         tempCanvForPdfPrint_wLogy->Print(outPdfFile.c_str());

      tempCanvForPdfPrint->cd();
      hratio->Draw();       tempCanvForPdfPrint->Print(outPdfFile.c_str());

      // draw and print response output, all of these plots draw empty though...
      //hgen_resp->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());              //empty
      //hgen_resp_anabin->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());   //empty     
      //hrec_resp->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());	      //empty
      //hrec_resp_anabin->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());   //empty

      tempCanvForPdfPrint->Print(close_outPdfFile.c_str());      
    }// end draw pdfs
  }  // end doBayes
  std::cout<<std::endl<<"Bayesian Unfolding done, all output written!"<<std::endl;


  // SVD unfolding ------------------------------------------------------------------------- 
  if(!doSVD) std::cout<<std::endl<<"   skipping SVD unfolding..."<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<std::endl<<"   beginning SVD unfolding..."<<std::endl;
    std::cout<<std::endl<<"initializing kReg parameter array w/ "<<nKregMax<<" elements"<<std::endl;
    std::cout<<kRegCenter-kRegRange<<" <= kReg[i] <= "<<kRegCenter+kRegRange<<std::endl;
    
    int kReg[nKregMax];
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) kReg[i+kRegRange]=kRegCenter+i;

    if(debugMode)for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) 
		   std::cout<<"kReg["<<i+kRegRange<<"]="<<kReg[i+kRegRange]<<std::endl;
    
    if(debugMode)std::cout<<std::endl<<"creating histo-arrays..."<<std::endl;
    TH1 *hrec_folded_ratio[nKregMax];     TH1 *hrec_unfolded_ratio[nKregMax];

    TH2 *hPearsonSVDPriorMeas[nKregMax];      //TH2D *hPearsonSVDPriorMeas[nKregMax];  
    TH1 *hFoldedSVDPriorMeas[nKregMax];     
    TH1F *hunf_svd[nKregMax];        
    TH1F *hratio_svd[nKregMax];      


    std::cout<<"calling RooUnfoldResponse"<<std::endl;
    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix_anabin"+Rstring).c_str());

    if(debugMode)std::cout<<"opening TCanvases..."<<std::endl;
    TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","", 1200, 1000);      cPearsonMatrixIter->Divide(3,3);	
    TCanvas *cRatio             = new TCanvas("cRatio","",             1200, 1000);      cRatio->Divide(3,3);      		
    TCanvas *cSpectra           = new TCanvas("cSpectra","",           1200, 1000);      cSpectra->Divide(3,3);	   
    TCanvas *cRatioCheck      = new TCanvas("cRatioCheck","",      1200, 1000);    
    TCanvas *c11 = new TCanvas("c11"," Singular Values and divectors", 1200, 1000);      c11->Divide(2);                    

    TLegend *leg[nKregMax],*leg1[nKregMax];      

    std::cout<<"entering SVD Unfolding Loop..."<<std::endl;
    for(int kr = 0; kr < nKregMax; ++kr){
      std::string kRegRandEtaRange_plotTitle=" kReg"+std::to_string(kReg[kr])+RandEtaRange_plotTitle;
      std::string kRegRandEtaRange="_kReg"+std::to_string(kReg[kr])+RandEtaRange;
      
      std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl;
   
      // unfold the MC dist, then create a ratio with the gen-level MC dist
      if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd and Hreco..."<<std::endl<<std::endl;
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      
      std::cout<<std::endl;   
      hunf_svd[kr]->SetName( ("ppData_SVDUnfolding_UnfoldSpectra"+kRegRandEtaRange).c_str() );
      hunf_svd[kr]->SetLineStyle(33);
      hunf_svd[kr]->SetLineColor(kBlue);
      hunf_svd[kr]->Print("base");
  
      hratio_svd[kr] = (TH1F*)hunf_svd[kr]->Clone( ("ppData_SVDUnfolding_UnfoldRatio"+kRegRandEtaRange).c_str());
      hratio_svd[kr]->Divide(hgen_anabin);
      hratio_svd[kr]->SetMarkerStyle(33);
      hratio_svd[kr]->SetMarkerColor(kRed);
      hratio_svd[kr]->Print("base");

      // get covariance matrix and calculate pearson coefficients
      if(debugMode)std::cout<<std::endl<<"calling Ereco..."<<std::endl;
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);


      if(debugPearson)std::cout<<"debugPearson is ON!"<<std::endl;

      if(debugPearson)std::cout<<std::endl;
      if(debugMode)std::cout<<"calling CalculatePearsonCoefficients..."<<std::endl;
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat, debugPearson );
      if(debugPearson)std::cout<<std::endl;

      if(debugMode)std::cout<<"creating TH2 for pearson matrix..."<<std::endl;
      hPearsonSVDPriorMeas[kr] = new TH2D (*pearson);
      //hPearsonSVDPriorMeas[kr] = reBinPearsonTH2( pearson );
      
      
      // drawing/visuals
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter canvas..."<<std::endl;
      cPearsonMatrixIter->cd(kr+1);
      //cPearsonMatrixIter->cd(kr+1)->SetLogz(1);
      //int kr_add1=kr+1;
      //i'm sorry, programming gods
      //TPad * thePad;
      //if       (kr<3)   thePad = new TPad( ("pad"+std::to_string(kr)).c_str(), "thepad" , 0.33*kr    , 0.67 , 0.33*(kr+1)     , 1.00 );
      //else if  (kr<6)   thePad = new TPad( ("pad"+std::to_string(kr)).c_str(), "thepad" , 0.33*(kr%3), 0.34 , 0.33*((kr%3)+1) , 0.67 );
      //else /*if(kr<9)*/ thePad = new TPad( ("pad"+std::to_string(kr)).c_str(), "thepad" , 0.33*(kr%3), 0.00 , 0.33*((kr%3)+1) , 0.34 );
      //cPearsonMatrixIter->SetLogz(1);
      //thePad->SetLogz(1);
      //thePad->Draw();
      //thePad->cd();
      
      //float ZminNmax=1.0;             //int XYmin=1; int XYmax=99;
      //int XYmin=1; int XYmax=1000;
      //hPearsonSVDPriorMeas[kr]->SetMinimum(-1*ZminNmax);  
      //hPearsonSVDPriorMeas[kr]->SetMaximum(1.0);
      hPearsonSVDPriorMeas[kr]->GetZaxis()->SetLabelSize(0.035);
      hPearsonSVDPriorMeas[kr]->GetXaxis()->SetTitle("reco bin num");
      hPearsonSVDPriorMeas[kr]->GetYaxis()->SetTitle("gen bin num");
      //hPearsonSVDPriorMeas[kr]->SetAxisRange(0, nbins_pt, "X");
      //hPearsonSVDPriorMeas[kr]->SetAxisRange(0, nbins_pt, "Y");
      //hPearsonSVDPriorMeas[kr]->SetAxisRange(.0000000001, 1, "Z");
      hPearsonSVDPriorMeas[kr]->SetAxisRange(-1., 1., "Z");
      hPearsonSVDPriorMeas[kr]->SetName( ("hPearsonSVDPriorMeas"+kRegRandEtaRange).c_str() );
      hPearsonSVDPriorMeas[kr]->SetTitle( ("pearsonMatrix SVDPriorMeas"+kRegRandEtaRange_plotTitle).c_str()  );
      hPearsonSVDPriorMeas[kr]->Draw("COLZ");

      if(debugMode)std::cout<<"applying roo_resp to \"truth\" histo hunf_svd[kr="<<kr<<"]..."<<std::endl;
      hFoldedSVDPriorMeas[kr] = roo_resp.ApplyToTruth(hunf_svd[kr]);
      hFoldedSVDPriorMeas[kr]->SetName( ("hFoldedSVDPriorMeas"+kRegRandEtaRange).c_str() );
      hFoldedSVDPriorMeas[kr]->SetLineStyle(33);
      hFoldedSVDPriorMeas[kr]->SetLineColor(kRed);
      



      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
      cSpectra->cd(kr+1);
      cSpectra->cd(kr+1)->SetLogy();

      // what's happening here exactly?
      //std::cout <<"CHECK: kr="<<kr<<"  and kReg[kr]="<<kReg[kr]<<std::endl<<std::endl;
      TH1F* hrec_anabin_clone=(TH1F*)hrec_anabin->Clone("hrec_anabin_clone");
      hrec_anabin_clone->SetTitle( ("Reco Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_anabin_clone->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_anabin_clone->SetMarkerStyle(24);
      hrec_anabin_clone->SetMarkerColor(kBlack);
      hrec_anabin_clone->SetAxisRange(45, 1000, "X");
      hrec_anabin_clone->Print("base");

      hrec_anabin_clone->Draw();
      hunf_svd[kr]->Draw("same");
      hFoldedSVDPriorMeas[kr]->Draw("same");

      leg[kr] = new TLegend(0.1, 0.1, 0.40, 0.30, NULL,"NBNDC");//x1,y1,x2,y2,header,option
      leg[kr]->AddEntry(hrec_anabin_clone,"Measured","pl");
      leg[kr]->AddEntry(hunf_svd[kr],"Unfolded","pl");
      leg[kr]->AddEntry(hFoldedSVDPriorMeas[kr],"Folded","pl");
      leg[kr]->SetTextSize(0.04);
      leg[kr]->Draw();

      // draw on cRatio canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio canvas..."<<std::endl;
      cRatio->cd(kr+1);

      hrec_folded_ratio[kr] = (TH1F*)hFoldedSVDPriorMeas[kr]->Clone( ("ppData_SVDUnfolding_FoldedRatio"+kRegRandEtaRange).c_str() );
      hrec_folded_ratio[kr]->SetTitle( ("Ratios w/ Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_folded_ratio[kr]->SetMarkerColor(kRed);
      hrec_folded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_folded_ratio[kr]->SetAxisRange(0., 1000., "X");
      hrec_folded_ratio[kr]->SetAxisRange(0.1, 1.9, "Y");
      hrec_folded_ratio[kr]->Divide(hrec_anabin);
      hrec_folded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Draw();

      hrec_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("ppData_SVDUnfolding_UnfoldedRatio"+kRegRandEtaRange).c_str());
      hrec_unfolded_ratio[kr]->SetTitle( ("Ratio, Unf. Over Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_unfolded_ratio[kr]->SetMarkerColor(kBlue);
      hrec_unfolded_ratio[kr]->Divide(hrec_anabin);
      hrec_unfolded_ratio[kr]->Print("base");
      hrec_unfolded_ratio[kr]->Draw("same");

      leg1[kr] = new TLegend(0.1, 0.1, 0.30, 0.3, NULL,"BRNDC");
      leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],"Unfold/Meas","pl");
      leg1[kr]->AddEntry(hrec_folded_ratio[kr],"Fold/Meas","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg1[kr]->SetTextSize(0.04); 
      leg1[kr]->Draw();
  	
      //  singular values and d_i vector ---------------------------
      //Note that these do not depend on the regularization.
      //The opposite: they tell you which regularization to use! (ian note: how?)
      if(kr == kRegDraw){      TSVDUnfold *svdUnfold = unf_svd.Impl();
        //if(debugMode) std::cout<<std::endl<<"    !!! kr==0 !!!"<<std::endl<<std::endl;

        // singular values
        if(debugMode)std::cout << "  getting singular values... " << std::endl;
        TH1 *hSVal = (TH1*)svdUnfold->GetSV();
        if(debugMode)for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
		       std::cout<<"bin: "<<bin<<",  SV: "<<hSVal->GetBinContent(bin)<< std::endl;

        if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
        c11->cd(1);
        gPad->SetLogy();        //c11->SetLogy();
        hSVal->SetTitle(" singular values ");
        hSVal->SetXTitle(" singular values ");
        hSVal->SetAxisRange(0,35,"X");
        hSVal->DrawCopy();
        drawText( "5.02 TeV pp, ak4PF Jets",0.358173, 0.8459761, 19);
	drawText( "Prompt-Reco, Jet80+LowerJets",0.358173, 0.8159761, 19);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);

        // di vector values
        if(debugMode)std::cout << "  getting di vector... " <<  std::endl;
        TH1D *hdi = (TH1D*)svdUnfold->GetD();
        if(debugMode)for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
		       std::cout<<"i: "<<bin<<",  di: "<<hdi->GetBinContent(bin)<<std::endl;

        if(debugMode)std::cout<<std::endl<<"drawing di vectors on c11 canvas.."<<std::endl;
        c11->cd(2);
        gPad->SetLogy();        //cdi->SetLogy();
        hdi->SetTitle(" di vectors ");
        hdi->SetXTitle(" |d_{i}^{kreg}| ");
        hdi->SetAxisRange(0,35,"X");
        hdi->DrawCopy();
        //drawText( "5.02 TeV ak4PFJets",0.358173, 0.8459761, 19);
	//drawText( "ppData, Prompt Reco HighPtJets",0.358173, 0.8159761, 19);
	//drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	//drawText( ("kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);
	const std:: string MCdesc= "QCD Dijets, PY8 Tune CUETP8M1"; 
        if(debugMode)std::cout<<std::endl<<"done with kr==0 specifics"<<std::endl<<std::endl;
      }

      hrec_unfolded_ratio[kr]->Print("base");
      hrec_folded_ratio[kr]->Print("base");
    }// kReg loop ends
    
    // loop over histos pointers in arrays to write
    fout->cd();

    std::cout<<std::endl<<"writing SVD histo arrays to file..."<<std::endl;    fout->cd();
    for(int kr = 0; kr<nKregMax; ++kr){
      hunf_svd[kr]->Write();            
      hratio_svd[kr]->Write();          
      hPearsonSVDPriorMeas[kr]->Write();
      hFoldedSVDPriorMeas[kr]->Write(); 
    }
    
    if(!drawPDFs){ if(debugMode)std::cout<<std::endl<<"NOT drawing PDFs for SVD Unfolding"<<std::endl<<std::endl; }
    else{ std::cout<<std::endl<<"drawing PDFs for SVD Unfolding..."<<std::endl;
      
      // form filename string, open pdf file and draw the canvases we have so far
      std::string outPdfFile=unfoldSpectra_outdir+outSVDPdfFile+".pdf";
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
      cPearsonMatrixIter->Print(open_outPdfFile.c_str());


      //  already-drawn canvases --------------------------------------------------
      cSpectra->cd()           ;   cSpectra->Print(outPdfFile.c_str());
      cRatio->cd()             ;   cRatio->Print(outPdfFile.c_str());
      cPearsonMatrixIter->cd() ;   cPearsonMatrixIter->Print(outPdfFile.c_str());
      c11->cd()                ;   c11->Print(outPdfFile.c_str());


      //  cRatio Check --------------------------------------------------
      cRatioCheck->cd();

      hrec_folded_ratio[kRegDraw]->SetAxisRange(0.1, 1.9, "Y");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(0., 1000., "X");
      hrec_folded_ratio[kRegDraw]->SetTitle("(Un)Fold/Meas");
      hrec_folded_ratio[kRegDraw]->Draw();

      TLegend * leg2 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
      leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Unfolded Data/Data","pl");
      leg2->AddEntry(hrec_folded_ratio[kRegDraw],"reFolded Data/Data","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2->SetTextSize(0.02);
      leg2->Draw();

      hrec_unfolded_ratio[kRegDraw]->Draw("same");
      drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
      drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
      drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);

      cRatioCheck->Print(outPdfFile.c_str());
      cRatioCheck->Print(close_outPdfFile.c_str());

      std::cout<<std::endl<<"done drawing SVD PDFs!"<<std::endl<<std::endl;
    }// end drawPDFs
  }// end SVD specific

  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();

  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_MC->Close();
  fpp_Data->Close();
  return 0;

} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=5 ){
    std::cout<<"do ./unfoldDataSpectra.exe "<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;
    return rStatus;  }
  
  rStatus=1; // runtime error
  rStatus=unfoldDataSpectra( (const std::string)argv[1], (const std::string)argv[2], (const std::string)argv[3], (int)std::atoi(argv[4]) ); 
  
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}

