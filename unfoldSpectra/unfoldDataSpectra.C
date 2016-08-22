#include "unfoldSpectra.h"


//  global settings -------------------------------------------------- 
const int radius=4;
const bool doToyErrs=false;

// Bayesian settings 
const int kIter = 4;

// SVD settings 
const int  nKregMax = 9, kRegRange=4, kRegDraw = 3, kRegCenter=21;

// procedural settings
const bool doBayes=true; 
const bool doSVD=true; 

const bool drawPDFs=true; 
const bool drawPDFs_BayesInputHistos= true&&drawPDFs;
const bool drawPDFs_SVDInputHistos  = true&&drawPDFs;

// I/O
const std::string baseName="unfoldDataSpectra_defOut"; // output name
const std::string outRootFile=baseName+".root";
const std::string outBayesPdfFile=baseName+"_BayesianUnfoldingPlots.pdf";
const std::string outSVDPdfFile=baseName+"_SVDUnfoldingPlots"; // see drawPDFs part for rest of string

const bool defDebugMode=false;
const bool defDebugPearson=false;   // for easy debugging of pearson matrix stuff

const std::string Rstring="_R"+std::to_string(radius); 
const std::string RandEtaRange=Rstring+"_20_eta_20";

const std::string Rstring_plotTitle=" R"+std::to_string(radius); 
const std::string RandEtaRange_plotTitle=Rstring_plotTitle+" 20eta20";

const double integratedLuminosity=25.8*pow(10,9)*0.99;
const std::string MCdesc= "MC, QCD PY8 Tune CUETP8M1"; 

//  the code --------------------------------------------------
int unfoldDataSpectra( bool debugMode=defDebugMode){
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl; 

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
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from DATA file"; 
  if(debugMode)std::cout<<": "<<inFile_Data_name; 
  std::cout<<std::endl<<std::endl;
  TFile *fpp_Data = TFile::Open(inFile_Data.c_str());

  //rec
  TH1F *hrec, *hrec_anabin;
  hrec = (TH1F*)fpp_Data->Get( ("hpp_TrgCombTest_JetID"+RandEtaRange).c_str() );
  hrec->Scale(1/integratedLuminosity);
  hrec->Print("base");

  hrec_anabin = (TH1F*)hrec->Clone( ("hpp_anaBin_HLTComb"+RandEtaRange).c_str() );
  hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins, ("hpp_anaBin_HLTComb"+RandEtaRange).c_str(), ptbins);
  divideBinWidth(hrec_anabin);
  //hrec_anabin->Scale(1/integratedLuminosity);
  hrec_anabin->Print("base");    
  

  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"; 
  if(debugMode)std::cout<<": "<<inFile_MC_name; 
  std::cout<<std::endl<<std::endl;
  TFile *fpp_MC = TFile::Open(inFile_MC.c_str());

  //mat, NOT NORMALIZED
  TH2F *hmat,*hmat_anabin;
  hmat = (TH2F*)fpp_MC->Get( ("hpp_matrix"+RandEtaRange).c_str() );
  hmat->Print("base");

  hmat_anabin = (TH2F*)hmat->Clone( ("hpp_anaBin_Trans_matrix_HLT"+RandEtaRange).c_str() );
  hmat_anabin->Rebin2D(10, 10);
  hmat_anabin->Print("base");    
  
  //gen
  TH1F*hgen,*hgen_anabin;
  hgen= (TH1F*)fpp_MC->Get( ("hpp_gen"+RandEtaRange).c_str() );
  hgen->Print("base");

  hgen_anabin = (TH1F*)hgen->Clone( ("hpp_anaBin_gen_JetComb"+RandEtaRange).c_str() );
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins, ("hpp_anaBin_gen_JetComb"+RandEtaRange).c_str(), ptbins);  //hgen_anabin = (TH1F*)hgen_anabin->Rebin(10);

  divideBinWidth(hgen_anabin);
  hgen_anabin->Scale(hrec_anabin->Integral()/hgen_anabin->Integral());
  hgen_anabin->Print("base");


  // output histos -------------------------
  std::cout<<std::endl<<std::endl<<"creating new OUTPUT histos..."<<std::endl;

  //gen
  TH1F *hgen_resp, *hgen_resp_anabin;
  hgen_resp = new TH1F( ("hpp_gen_response"+RandEtaRange).c_str(),"", 
			hgen->GetNbinsX(), hgen->GetXaxis()->GetXmin(), hgen->GetXaxis()->GetXmax());
  hgen_resp->Sumw2();
  hgen_resp->Print("base");
  hgen_resp_anabin = new TH1F( ("hpp_gen_response_anabin"+RandEtaRange).c_str() ,"", nbins, ptbins);
  hgen_resp_anabin->Sumw2();
  hgen_resp_anabin->Print("base");
  
  //rec
  if(debugMode)std::cout<<"opening output rec/recresponse histograms..."<<std::endl;
  TH1F *hrec_resp, *hrec_resp_anabin; 
  hrec_resp = new TH1F( ("hpp_rec_response"+RandEtaRange).c_str(),"", 
		       hrec->GetNbinsX(), hrec->GetXaxis()->GetXmin(), hrec->GetXaxis()->GetXmax());
  hrec_resp->Sumw2();
  hrec_resp->Print("base");
  hrec_resp_anabin = new TH1F( ("hpp_rec_response_anabin"+RandEtaRange).c_str(),"", nbins, ptbins);
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
      std::string outPdfFile=outBayesPdfFile;
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";

      // temp canvas for printing
      TCanvas* tempCanvForPdfPrint=new TCanvas("","",800,800);    
      tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); 
      tempCanvForPdfPrint->cd();

      // draw and print unfold/ratio plots
      tempCanvForPdfPrint->SetLogy(1);
      hunf->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());

      tempCanvForPdfPrint->SetLogy(0);
      hratio->Draw();       tempCanvForPdfPrint->Print(outPdfFile.c_str());

      // draw and print response output, all of these plots draw empty though...
      //hgen_resp->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());              //empty
      //hgen_resp_anabin->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());   //empty     
      //hrec_resp->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());	      //empty
      //hrec_resp_anabin->Draw();         tempCanvForPdfPrint->Print(outPdfFile.c_str());   //empty

      if(drawPDFs_BayesInputHistos){
	std::cout<<std::endl<<"drawing input histos to Bayesian Unfolding..."<<std::endl;
        hmat->SetTitle("ppMC jet input, genpt v. recopt");
        hmat->GetYaxis()->SetLabelSize(0.02);
        hmat->GetYaxis()->SetTitleSize(0.04);
        hmat->GetXaxis()->SetLabelSize(0.02);
        hmat->GetXaxis()->SetTitleSize(0.04);
        tempCanvForPdfPrint->SetLogy(0);
        hmat->Draw();           tempCanvForPdfPrint->Print(outPdfFile.c_str());

        hgen->SetMarkerStyle(24);
        hgen->SetMarkerColor(kRed);
        hgen->SetTitle("ppMC genjet input, genpt");
        tempCanvForPdfPrint->SetLogy(1);
        hgen->Draw();           tempCanvForPdfPrint->Print(outPdfFile.c_str());

	TH1F* hrec_clone=(TH1F*)hrec->Clone("hrec_clone");
	hrec_clone->SetAxisRange(0.,1000.,"X");
        hrec_clone->SetMarkerStyle(24);
        hrec_clone->SetMarkerColor(kRed);
        hrec_clone->SetTitle("ppData recojet input, recopt");
        tempCanvForPdfPrint->SetLogy(1);
        hrec_clone->Draw();           tempCanvForPdfPrint->Print(outPdfFile.c_str());
      }
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

      bool debugPearson=(bool)(defDebugPearson && debugMode) ;
      if(debugPearson)std::cout<<"debugPearson is ON!"<<std::endl;

      if(debugPearson)std::cout<<std::endl;
      if(debugMode)std::cout<<"calling CalculatePearsonCoefficients..."<<std::endl;
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat, debugPearson );
      if(debugPearson)std::cout<<std::endl;

      if(debugMode)std::cout<<"creating TH2 for pearson matrix..."<<std::endl;
      hPearsonSVDPriorMeas[kr] = new TH2D (*pearson);
      
      // drawing/visuals
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter canvas..."<<std::endl;
      cPearsonMatrixIter->cd(kr+1);

      float ZminNmax=1.0;             int XYmin=1; int XYmax=99;
      hPearsonSVDPriorMeas[kr]->SetMinimum(-1*ZminNmax);  hPearsonSVDPriorMeas[kr]->SetMaximum(ZminNmax);
      hPearsonSVDPriorMeas[kr]->GetZaxis()->SetLabelSize(0.035);
      hPearsonSVDPriorMeas[kr]->SetAxisRange(XYmin, XYmax, "X");
      hPearsonSVDPriorMeas[kr]->SetAxisRange(XYmin, XYmax, "Y");
      hPearsonSVDPriorMeas[kr]->SetName( ("hPearsonSVDPriorMeas"+kRegRandEtaRange).c_str() );
      hPearsonSVDPriorMeas[kr]->SetTitle( ("pearsonMatrix SVDPriorMeas"+kRegRandEtaRange_plotTitle).c_str()  );
      hPearsonSVDPriorMeas[kr]->Draw("colz");

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
      hrec_folded_ratio[kr]->SetAxisRange(45, 1000, "X");
      hrec_folded_ratio[kr]->SetAxisRange(0.5, 1.5, "Y");
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
      if(kr == 0){      TSVDUnfold *svdUnfold = unf_svd.Impl();
        if(debugMode) std::cout<<std::endl<<"    !!! kr==0 !!!"<<std::endl<<std::endl;

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
        drawText( "5.02 TeV pp, ak4PFJets",0.358173, 0.8459761, 19);
	drawText( "Data, Prompt-Reco, HighPtJets",0.358173, 0.8159761, 19);
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
      std::string outPdfFile=outSVDPdfFile+".pdf";
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
      cPearsonMatrixIter->Print(open_outPdfFile.c_str());


      //  already-drawn canvases --------------------------------------------------
      cSpectra->cd()           ;   cSpectra->Print(outPdfFile.c_str());
      cRatio->cd()             ;   cRatio->Print(outPdfFile.c_str());
      cPearsonMatrixIter->cd() ;   cPearsonMatrixIter->Print(outPdfFile.c_str());
      c11->cd()                ;   c11->Print(outPdfFile.c_str());


      //  cRatio Check --------------------------------------------------
      cRatioCheck->cd();

      hrec_folded_ratio[kRegDraw]->SetAxisRange(0.4, 1.2, "Y");
      hrec_folded_ratio[kRegDraw]->SetAxisRange(45, 1000, "X");
      hrec_folded_ratio[kRegDraw]->SetTitle("(Un)Fold/Meas");
      hrec_folded_ratio[kRegDraw]->Draw();

      TLegend * leg2 = new TLegend(0.1, 0.1, 0.30, 0.20, NULL,"NBNDC");
      leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Unfold/Meas","pl");
      leg2->AddEntry(hrec_folded_ratio[kRegDraw],"Fold/Meas","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2->SetTextSize(0.03);
      leg2->Draw();

      hrec_unfolded_ratio[kRegDraw]->Draw("same");
      drawText( "5.02 TeV pp, ak4PFJets",                               0.13, 0.46, 22);
      drawText( "Data, Prompt-Reco HighPtJets",                  0.13, 0.43, 22);
      drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.13, 0.40, 22);

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

  if( argc!=1 && argc!=2 ) { // input error
    std::cout<<"do ./unfoldDataSpectra.exe <debugMode= 0or1, def="<<defDebugMode<<" to run"<<std::endl;
    return rStatus;  }

  rStatus=1; // runtime error
  if(argc==1)   rStatus=unfoldDataSpectra();
  if(argc==2)   rStatus=unfoldDataSpectra( (const bool) atoi(argv[1]) ); 

  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}

