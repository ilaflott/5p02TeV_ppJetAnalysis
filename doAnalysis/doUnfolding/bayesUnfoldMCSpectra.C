#include "unfoldSpectra.h"

//other settings
const bool drawPDFs=true; 
const bool debugMode=false;



// CODE --------------------------------------------------
int bayesUnfoldMCSpectra( std::string inFile_MC_dir , const std::string baseName , 
		     const bool doJetID=true , const bool useSimpBins=false){
  
  
  
  // BINNING -----------
  
  if(!useSimpBins)std::cout<<"using analysis pt bins"<<std::endl<<std::endl;
  else std::cout<<"using simple pt bins"<<std::endl<<std::endl;
  
  double* boundaries_pt_gen=setBinning(  useSimpBins , "gen"  );
  int     nbins_pt_gen=setNBins(  useSimpBins , "gen"  );
  
  double* boundaries_pt_reco      = setBinning( useSimpBins , "reco" );
  int     nbins_pt_reco =   setNBins(  useSimpBins , "reco"  );
  
  double* boundaries_pt_gen_mat   = setBinning( useSimpBins , "gen"  );
  int     nbins_pt_gen_mat =   setNBins(  useSimpBins , "gen"  );
  
  double* boundaries_pt_reco_mat  = setBinning( useSimpBins , "reco" );
  int     nbins_pt_reco_mat =   setNBins(  useSimpBins , "reco"  );
  
  

  // STRINGS -----------
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl; 
  //inFile_Data_dir=SCRATCH_BASE+inFile_Data_dir;
  inFile_MC_dir  =SCRATCH_BASE+inFile_MC_dir;
  
  std::size_t radPos=inFile_MC_dir.find("_ak")+3;  
  const std::string radius="R"+inFile_MC_dir.substr( radPos,1 );
  const std::string RandEtaRange="_"+radius+"_20_eta_20";
  const std::string RandEtaRange_plotTitle=" "+radius+" 20eta20";// R4 20eta20
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  if(debugMode)std::cout<<"RandEtaRange string is = "<<RandEtaRange<<std::endl;  
  
  std::size_t jetTypePos=radPos+1, jetsPos=inFile_MC_dir.find("Jets");
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );//PFJets
  const std::string fullJetType="ak"+inFile_MC_dir.substr( radPos,1 )+jetType;//"ak4PFJets"
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  
  
  
  // INFILE NAME(S) -----------
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  //const std::string inFile_Data_name="/HighPtJetTrig_"+fullJetType+"-allFiles.root";
  
  
  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldMCSpectra_outdir+fullJetType;
  if(doJetID)outFileName+="_wjtID";
  if(!useSimpBins)outFileName+="_anabins";
  else outFileName+="_simpbins";
  outFileName+="_"+baseName;
  
  std::string outRespMatPdfFile =  outFileName+"_respMat.pdf";
  std::string outBayesPdfFile =  outFileName+".pdf";
  std::string outRootFile     =  outFileName+".root";  
  

  
  // ERROR/WEIGHTS/STATS HANDLING ------------------
  RooUnfold::ErrorTreatment errorTreatment;
  if(!doToyErrs) errorTreatment = RooUnfold::kCovariance;
  else errorTreatment = RooUnfold::kCovToy; 
  if(debugMode)std::cout<<"doToyErrs="<<doToyErrs<<std::endl; 
  
  if(debugMode)std::cout<<"TH2 GetDefaultSumw2="<<TH2::GetDefaultSumw2()<<std::endl;
  if(debugMode)std::cout<<"setting TH2DefSumw2 to true..."<<std::endl; 
  TH2::SetDefaultSumw2(true);
  
  if(debugMode)std::cout<<"TH1 GetDefaultSumw2="<<TH1::GetDefaultSumw2()<<std::endl;
  if(debugMode)std::cout<<"setting TH1DefSumw2 to true..."<<std::endl; 
  TH1::SetDefaultSumw2(true);
  
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();



  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
  std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());



  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE");   
  fout->cd();  
  
  
  // ---------- reco, measured spectra to unfold
  std::string histTitle="hpp_mcclosure_reco_test";
  if(doJetID)histTitle+="_wJetID";
  histTitle+=RandEtaRange;
  
  TH1F*  hrec = (TH1F*)fpp_MC->Get( histTitle.c_str() ); 
  hrec->Write(histTitle.c_str());
  if(debugMode)hrec->Print("base");

  histTitle+="_divBylumietabin";
  hrec->Scale(1./4.); // lumi
  hrec->Write(histTitle.c_str());
  if(debugMode)hrec->Print("base");
  
  histTitle+="_clone";
  TH1F *hrec_anabin = (TH1F*)hrec->Clone( (histTitle).c_str() );
  hrec_anabin->Write();
  if(debugMode)hrec_anabin->Print("base");
  
  std::cout<<"rebinning hrec..."<<std::endl;
  histTitle+="_anabins";
  hrec_anabin = (TH1F*)hrec_anabin->Rebin( nbins_pt_reco, (histTitle).c_str() , boundaries_pt_reco);
  hrec_anabin->Write(histTitle.c_str());   
  if(debugMode)hrec_anabin->Print("base");  
  
  histTitle+="_normbinwidth";
  divideBinWidth(hrec_anabin); 
  hrec_anabin->Write(histTitle.c_str());
  if(debugMode)hrec_anabin->Print("base");  
  
  if(clearOverUnderflows){
    histTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_anabin);
    hrec_anabin->Write(histTitle.c_str());
    if(debugMode)hrec_anabin->Print("base");  
  }
  
  
  
  
  
  //// response hist, for output? what is this for if it's empty?
  //TH1F* hrec_resp_anabin;
  //if(fillRespHists) hrec_resp_anabin = (TH1F*)hrec_anabin->Clone("recanabinClone4unf");
  //else{    hrec_resp_anabin = new TH1F( ("hpp_rec_response_anabin"+RandEtaRange).c_str(),"", 
  //					nbins_pt_reco, boundaries_pt_reco); }
  //hrec_resp_anabin->Write();
  //if(debugMode)hrec_resp_anabin->Print(" base");  
  
  
  
  
  // ---------- reco, measured same-side spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="hpp_mcclosure_reco";
  if(doJetID)histTitle2+="_wJetID";
  histTitle2+=RandEtaRange;
  
  TH1F*  hrec_sameside = (TH1F*)fpp_MC->Get( histTitle2.c_str() ); 
  hrec_sameside->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  histTitle2+="_divByetabin";
  hrec_sameside->Scale(1./4.); // eta bin width for 0.<|y|<2.
  hrec_sameside->Write( histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  histTitle2+="_clone";
  TH1F *hrec_sameside_anabin = (TH1F*)hrec_sameside->Clone( (histTitle2).c_str() );
  hrec_sameside_anabin->Write();
  if(debugMode)hrec_sameside_anabin->Print("base");
  
  std::cout<<"rebinning hrec_sameside..."<<std::endl;
  histTitle2+="_anabins";
  hrec_sameside_anabin = (TH1F*)hrec_sameside_anabin->Rebin( nbins_pt_reco, (histTitle2).c_str() , boundaries_pt_reco);
  hrec_sameside_anabin->Write();   
  if(debugMode)hrec_sameside_anabin->Print("base");  
  
  histTitle2+="_normbinwidth";
  divideBinWidth(hrec_sameside_anabin); 
  hrec_sameside_anabin->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside_anabin->Print("base");  
  
  if(clearOverUnderflows){
    histTitle2+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_sameside_anabin);
    hrec_sameside_anabin->Write(histTitle2.c_str());
    if(debugMode)hrec_sameside_anabin->Print("base");    }
  
  // response hist, for output? what is this for if it's empty?
  TH1F* hrec_sameside_resp_anabin;
  //if(fillRespHists) hrec_sameside_resp_anabin = (TH1F*)hrec_sameside_anabin->Clone("recanabinsamesideClone4unf");
  if(!useSimpBins) hrec_sameside_resp_anabin = (TH1F*)hrec_sameside_anabin->Clone("recanabinsamesideClone4unf");
  else{    hrec_sameside_resp_anabin = new TH1F( ("hpp_rec_sameside_response_anabin"+RandEtaRange).c_str(),"", 
						 nbins_pt_reco, boundaries_pt_reco); }
  hrec_sameside_resp_anabin->Write();
  if(debugMode)hrec_sameside_resp_anabin->Print(" base");  



  if(doMCIntegralScaling){
    std::cout<<"scaling hrec_anabin to hrec_sameside_anabin"<<std::endl;
    histTitle+="_integNormd";
    float integral_rec_anabin=hrec_anabin->Integral();
    float integral_rec_ss_anabin=hrec_sameside_anabin->Integral();
    std::cout<<"scale factor="<<integral_rec_ss_anabin/integral_rec_anabin<<std::endl;
    hrec_anabin->Scale(integral_rec_ss_anabin/integral_rec_anabin);
    hrec_anabin->Write(histTitle.c_str());
    if(debugMode)hrec_anabin->Print("base");  }
  
  
  
  
  
  
  
  
  // ---------- gen, MC truth spectra
  std::string genHistTitle="hpp_mcclosure_gen";
  if(doJetID)genHistTitle+="_wJetID";
  genHistTitle+=RandEtaRange;
  
  TH1F* hgen = (TH1F*)fpp_MC->Get( genHistTitle.c_str() );
  hgen->Write();
  if(debugMode)hgen->Print("base");    
  
  genHistTitle+="_divByetabin";
  hgen->Scale(1./4.); // eta bin width for 0.<|y|<2.
  hgen->Write( genHistTitle.c_str() );
  if(debugMode)hgen->Print("base");


  genHistTitle+="_clone";
  TH1F* hgen_anabin = (TH1F*)hgen->Clone( (genHistTitle).c_str() );
  hgen_anabin->Write(genHistTitle.c_str());
  if(debugMode)hgen_anabin->Print("base");
  
  genHistTitle+="_anabins";
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins_pt_gen, (genHistTitle).c_str() , boundaries_pt_gen);
  hgen_anabin->Write(genHistTitle.c_str());
  if(debugMode)hgen_anabin->Print("base"); 
  
  genHistTitle+="_normbinwidth";
  divideBinWidth(hgen_anabin);
  hgen_anabin->Write(genHistTitle.c_str());
  if(debugMode)hgen_anabin->Print("base");  
  
  if(clearOverUnderflows){
    genHistTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hgen_anabin);
    hgen_anabin->Write(genHistTitle.c_str());
    if(debugMode)hgen_anabin->Print("base");    }
  
  TH1F* hgen_resp_anabin;
  //if(fillRespHists) hgen_resp_anabin = (TH1F*)hgen_anabin->Clone("genanabinClone4unf");
  if(!useSimpBins) hgen_resp_anabin = (TH1F*)hgen_anabin->Clone("genanabinClone4unf");
  else{    hgen_resp_anabin = new TH1F( ("hpp_gen_response_anabin"+RandEtaRange).c_str() ,"", 
  					nbins_pt_gen, boundaries_pt_gen);  }
  hgen_resp_anabin->Write();
  if(debugMode)hgen_resp_anabin->Print("base");  
  
  
  
  
  //std::cout<<std::endl<<"writing input hists to file..."<<std::endl;
  //
  //fout->cd();
  //
  //hgen->Write(); 
  //hrec->Write();
  //hgen_anabin->Write(); 
  //hrec_anabin->Write(); 
  //
  //hrec_sameside->Write();
  //hrec_sameside_anabin->Write(); 
  //
  //if(fillRespHists) hgen_resp_anabin->Write();
  //if(fillRespHists) hrec_resp_anabin->Write();
  
  
  
  
  
  // ---------- open MC "response" matrix
  std::string TH2_title="hpp_mcclosure_matrix";//+RandEtaRange;
  if(doJetID)TH2_title+="_wJetID";
  TH2_title+=RandEtaRange;
  
  //get the response matrix made by readforests
  TH2F* hmat = (TH2F*)fpp_MC->Get( TH2_title.c_str() );
  hmat->Write();
  if(debugMode)hmat->Print("base");
  
  TH2_title+="_divByetabin";
  hmat->Scale(1./4.); // eta bin width for 0.<|y|<2.
  hmat->Write( TH2_title.c_str());
  if(debugMode)hmat->Print("base");
  
  // rebinned matrix ---------------
  TH2_title+="_clone";
  TH2F* hmat_anabin = (TH2F*)hmat->Clone( (TH2_title).c_str() );
  hmat_anabin->Write(TH2_title.c_str());
  if(debugMode)hmat_anabin->Print("base"); 

  TH2_title+="_anabins";
  hmat_anabin=(TH2F*) reBinTH2(hmat_anabin, (TH2_title).c_str(), 
			       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );  
  hmat_anabin->Write(TH2_title.c_str());
  if(debugMode)hmat_anabin->Print("base"); 
  
  TH2_title+="_normbinwidth";
  divideBinWidth_TH2(hmat_anabin);
  hmat_anabin->Write(TH2_title.c_str());
  if(debugMode)hmat_anabin->Print("base"); 
  
  if(clearOverUnderflows){
    TH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2F*)hmat_anabin);
    hmat_anabin->Write(TH2_title.c_str());
    if(debugMode)hmat_anabin->Print("base");  }
  
  
  // error and %error for response matrix ---------------
  std::string errTH2_title="hmat_errors_anabins";
  TH2F* hmat_errors=makeRespMatrixErrors( (TH2F*) hmat,
					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  hmat_errors->Write(errTH2_title.c_str());
  if(debugMode)hmat_errors->Print("base");
  
  errTH2_title+="_normbinwidth";
  divideBinWidth_TH2(hmat_errors);
  hmat_errors->Write(errTH2_title.c_str());
  if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    errTH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2F*)hmat_errors);
    hmat_errors->Write(errTH2_title.c_str());
    if(debugMode)hmat_errors->Print("base");  }
  
  
  // give response matrix the correct errors
  setRespMatrixErrs( (TH2F*)hmat_anabin, (TH2F*) hmat_errors , (bool)zeroBins);

  TH2_title+="_wseterrs";
  hmat_anabin->Write(TH2_title.c_str());
  if(debugMode)hmat_anabin->Print("base");
  
  
  TH2F* hmat_percenterrs= makeRespMatrixPercentErrs( (TH2F*) hmat_errors, (TH2F*) hmat_anabin,
						     (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						     (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
  hmat_percenterrs->Write();
  if(debugMode) hmat_percenterrs->Print("base");
  
  
  
  
  // row/col normd matrix, has correct errors since hmat_anabin has correct errors ---------------
  // POTENTIAL ISSUE: using the resp matrix post rebinning/divide bin width/clearing overflows... should i be using the original matrix?
  TH2F* hmat_anabin_colnormd = normalizeCol_RespMatrix( (TH2F*)  hmat_anabin,
							(double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
							(double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  hmat_anabin_colnormd->Write();
  if(debugMode)  hmat_anabin_colnormd->Print("base");

  
  TH2F*  hmat_anabin_rownormd = normalizeRow_RespMatrix( (TH2F*)  hmat_anabin,
							 (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
							 (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  hmat_anabin_rownormd->Write();
  if(debugMode)  hmat_anabin_rownormd->Print("base");
  




  bool drawRespMatrix=false;
  if(drawPDFs && drawRespMatrix){    
    
    
    std::cout<<std::endl<<"drawing input response matrices..."<<std::endl;
    
    std::string outPdfFile=outRespMatPdfFile;
    std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
    
    TCanvas* tempCanvForPdfPrint=new TCanvas("tempCanv_respMat","",1200,1200);    
    tempCanvForPdfPrint->cd();
    
    if(useSimpBins){	tempCanvForPdfPrint->SetLogx(0);
      tempCanvForPdfPrint->SetLogy(0);       
      tempCanvForPdfPrint->SetLogz(1);          }
    else {	tempCanvForPdfPrint->SetLogx(1);
      tempCanvForPdfPrint->SetLogy(1);       
      tempCanvForPdfPrint->SetLogz(1);         }       
    
    // open file
    
    tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); 
    
    // orig matrix ---------------
    
    matStylePrint(hmat, "ppMC Resp Matrix, original", tempCanvForPdfPrint, outPdfFile, useSimpBins);      
    
    // orig matrix w/ used pt range ---------------
    
    tempCanvForPdfPrint->cd();
    
    hmat->SetTitle("ppMC Resp Matrix, used pt range");
    hmat->SetAxisRange(boundaries_pt_reco_mat[0],boundaries_pt_reco_mat[nbins_pt_reco_mat],"X");
    hmat->SetAxisRange(boundaries_pt_gen_mat[0],boundaries_pt_gen_mat[nbins_pt_gen_mat],"Y");      
    hmat->Draw("COLZ");           
    
    tempCanvForPdfPrint->Print(outPdfFile.c_str());      
    
    // matrix rebinned ---------------
    
    matStylePrint(hmat_anabin, "ppMC Resp Matrix, rebinned", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // error matrix in binning of interest ---------------
    
    matStylePrint(hmat_errors, "ppMC Resp Matrix errors", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // percent error matrix in binning of interest ---------------
    
    matStylePrint(hmat_percenterrs, "ppMC Resp Matrix % Errors", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // col normd matrix in binning of interest  ---------------
    
    matStylePrint(hmat_anabin_colnormd, "ppMC Resp Matrix, Columns Normalized", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // row normd matrix in binning of interest  ---------------
    
    matStylePrint(hmat_anabin_rownormd, "ppMC Resp Matrix, Rows Normalized", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // close file 
    
    tempCanvForPdfPrint->Print(close_outPdfFile.c_str());   }
  













  // Bayesian unfolding ------------------------- 

    // RooUnfoldResponse, Bayes --------------------
    std::cout<<"calling RooUnfoldResponse "<<std::endl;
    
    //for inversion of response matrix
    RooUnfoldResponse roo_resp_nullresp( 0, 0, hmat_anabin, ("Response_matrix_nullresp"+RandEtaRange).c_str()) ;
    roo_resp_nullresp.UseOverflow(doOverUnderflows);    
    roo_resp_nullresp.Write();
    
    RooUnfoldResponse roo_resp( hrec_sameside_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix"+RandEtaRange).c_str()) ;
    roo_resp.UseOverflow(doOverUnderflows);    
    roo_resp.Write();
    
    std::cout<<"calling RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );
    
    TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);     std::cout<<std::endl; 
    hunf->SetName("ppMC_BayesUnf_Spectra");
    hunf->SetTitle("ppMC BayesUnf Spectra");
    if(debugMode)hunf->Print("base");
    
    
    std::cout<<"reFolding unfolded histogram!!"<<std::endl;
    TH1F* hfold=(TH1F*)roo_resp_nullresp.ApplyToTruth(hunf);        
    //TH1F* hfold=(TH1F*)roo_resp.ApplyToTruth(hunf);        
    hfold->SetName("ppMC_BayesFold_Spectra");
    hfold->SetTitle("ppMC BayesFold Spectra");
    if(debugMode)hfold->Print("base");

    //if(doMCIntegralScaling)hfold->Scale(hrec_anabin->Integral()/hfold->Integral());
    if(debugMode)hfold->Print("base");

    // RATIOS WITH OPP SIDE GEN
    TH1F *hgen_anabin_ratiobin=(TH1F*)hgen_anabin->Clone("ppMC_Gen_Ratio_denom");
    hgen_anabin_ratiobin=(TH1F*)hgen_anabin_ratiobin->Rebin(nbins_pt_reco,"ppMC_Gen_Ratio_denom_rebin",boundaries_pt_reco);
    if(debugMode)hgen_anabin_ratiobin->Print("base");
    
    TH1F *h_genratio_oppunf = (TH1F*)hunf->Clone( "ppMC_Gen_Ratio_OppUnf" );
    h_genratio_oppunf->SetTitle( "Bayes Unf./Gen. Truth" );
    h_genratio_oppunf->SetMarkerStyle(24);
    h_genratio_oppunf->SetMarkerColor(2);
    h_genratio_oppunf->Divide(hgen_anabin);
    if(debugMode)h_genratio_oppunf->Print("base");
    
    TH1F *h_genratio_oppfold = (TH1F*)hfold->Clone( "ppMC_Gen_Ratio_OppFold" );
    h_genratio_oppfold->SetTitle( "Bayes Fold/Gen. Truth" );
    h_genratio_oppfold->SetMarkerStyle(24);
    h_genratio_oppfold->SetMarkerColor(3);
    h_genratio_oppfold->Divide(hgen_anabin_ratiobin);
    if(debugMode)h_genratio_oppfold->Print("base");
    
    TH1F *h_genratio_oppmeas = (TH1F*)hrec_anabin->Clone( "ppMC_Gen_Ratio_Meas" );
    h_genratio_oppmeas->SetTitle( "Meas./Gen. Truth" );
    h_genratio_oppmeas->SetMarkerStyle(24);
    h_genratio_oppmeas->SetMarkerColor(4);
    h_genratio_oppmeas->Divide(hgen_anabin_ratiobin);
    if(debugMode)h_genratio_oppmeas->Print("base");
    
    
    // RATIOS WITH OPP SIDE MEAS
    TH1F *h_recratio_oppunf = (TH1F*)hunf->Clone( "ppMC_Meas_Ratio_OppUnf" );
    h_recratio_oppunf=(TH1F*)h_recratio_oppunf->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_OppUnf_rebin" , boundaries_pt_reco);
    h_recratio_oppunf->SetTitle( "Bayes Unf./Meas." );
    h_recratio_oppunf->SetMarkerStyle(24);
    h_recratio_oppunf->SetMarkerColor(2);
    h_recratio_oppunf->Divide(hrec_anabin);
    if(debugMode)h_recratio_oppunf->Print("base");
    
    TH1F *h_recratio_oppfold = (TH1F*)hfold->Clone( "ppMC_Meas_Ratio_OppFold" );
    //h_recratio_oppfold=(TH1F*)h_recratio_oppfold->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_OppFold_rebin" , boundaries_pt_reco);
    h_recratio_oppfold->SetTitle( "Bayes Fold./Meas." );
    h_recratio_oppfold->SetMarkerStyle(24);
    h_recratio_oppfold->SetMarkerColor(3);
    h_recratio_oppfold->Divide(hrec_anabin);
    if(debugMode)h_recratio_oppfold->Print("base");
    
    


    // SAMESIDE UNFOLDING
    RooUnfoldResponse roo_ss_resp_nullresp( 0,0, hmat_anabin, ("ss_Response_matrix_nullresp"+RandEtaRange).c_str()) ;
    roo_ss_resp_nullresp.UseOverflow(doOverUnderflows);    
    roo_ss_resp_nullresp.Write();
    
    RooUnfoldResponse roo_ss_resp( hrec_sameside_resp_anabin, hgen_resp_anabin, hmat_anabin, ("ss_Response_matrix"+RandEtaRange).c_str()) ;
    roo_ss_resp.UseOverflow(doOverUnderflows);    
    roo_ss_resp.Write();
    
    std::cout<<"calling sameside RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_sameside_bayes( &roo_ss_resp, hrec_sameside_anabin, kIter );
    
    TH1F *hunf_sameside = (TH1F*)unf_sameside_bayes.Hreco(errorTreatment);    std::cout<<std::endl; 
    hunf_sameside->SetName("ppMC_BayesUnf_sameSideSpectra");
    hunf_sameside->SetTitle("ppMC BayesUnf sameSideSpectra");
    if(debugMode)hunf_sameside->Print("base");
    
    std::cout<<"reFolding unfolded histogram!!"<<std::endl;
    TH1F* hfold_sameside=(TH1F*)roo_ss_resp_nullresp.ApplyToTruth(hunf_sameside);    
    //TH1F* hfold_sameside=(TH1F*)roo_ss_resp.ApplyToTruth(hunf_sameside);    
    hfold_sameside->SetName("ppMC_BayesFold_sameSide_Spectra");
    hfold_sameside->SetTitle("ppMC BayesFold sameSide Spectra");
    if(debugMode)hfold_sameside->Print("base");

    //if(doMCIntegralScaling)hfold_sameside->Scale(hrec_sameside_anabin->Integral()/hfold_sameside->Integral());
    if(debugMode)hfold_sameside->Print("base");

    
    
    // SAMESIDE RATIO W/ OPPSIDE MEAS
    TH1F *h_recratio_ssunf = (TH1F*)hunf_sameside->Clone( "ppMC_Meas_Ratio_SSUnf" );
    h_recratio_ssunf=(TH1F*)h_recratio_ssunf->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSUnf_rebin" , boundaries_pt_reco);
    h_recratio_ssunf->SetTitle( "SameSide Bayes Unf./Meas." );
    h_recratio_ssunf->SetMarkerStyle(25);
    h_recratio_ssunf->SetMarkerColor(2);
    h_recratio_ssunf->Divide(hrec_anabin);
    if(debugMode)h_recratio_ssunf->Print("base");
    
    
    TH1F *h_recratio_ssfold = (TH1F*)hfold_sameside->Clone( "ppMC_Meas_Ratio_SSFold" );
    //h_recratio_ssfold=(TH1F*)h_recratio_ssfold->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSFold_rebin" , boundaries_pt_reco);
    h_recratio_ssfold->SetTitle( "SameSide Bayes Fold./Meas." );
    h_recratio_ssfold->SetMarkerStyle(25);
    h_recratio_ssfold->SetMarkerColor(3);
    h_recratio_ssfold->Divide(hrec_anabin);
    if(debugMode)h_recratio_ssfold->Print("base");
    
    
    TH1F *h_recratio_ssgen = (TH1F*)hgen_anabin->Clone( "ppMC_Meas_Ratio_SSTruth" );    
    h_recratio_ssgen=(TH1F*)h_recratio_ssgen->Rebin(nbins_pt_reco, "ppMC_Meas_Ratio_SSTruth_rebin" , boundaries_pt_reco);
    h_recratio_ssgen->SetTitle( "Truth/Meas." );
    h_recratio_ssgen->SetMarkerStyle(24);
    h_recratio_ssgen->SetMarkerColor(6);
    h_recratio_ssgen->Divide(hrec_anabin);
    if(debugMode)h_recratio_ssgen->Print("base");
    
    
    // SAMESIDE RATIO W/ OPPSIDE GEN
    TH1F *h_genratio_ssunf = (TH1F*)hunf_sameside->Clone( "ppMC_Gen_Ratio_SSUnf" );
    h_genratio_ssunf->SetTitle( "SameSide Bayes Unf./Gen. Truth" );
    h_genratio_ssunf->SetMarkerStyle(25);
    h_genratio_ssunf->SetMarkerColor(2);
    h_genratio_ssunf->Divide(hgen_anabin);
    if(debugMode)h_genratio_ssunf->Print("base");
    
    
    TH1F *h_genratio_ssfold = (TH1F*)hfold_sameside->Clone( "ppMC_Gen_Ratio_SSFold" );
    h_genratio_ssfold->SetTitle( "SameSide Bayes Fold/Gen. Truth" );
    h_genratio_ssfold->SetMarkerStyle(25);
    h_genratio_ssfold->SetMarkerColor(3);
    h_genratio_ssfold->Divide(hgen_anabin_ratiobin);
    if(debugMode)h_genratio_ssfold->Print("base");
    
    TH1F *h_genratio_ssmeas = (TH1F*)hrec_sameside_anabin->Clone( "ppMC_Gen_Ratio4_SSMeas" );
    h_genratio_ssmeas->SetTitle( "Same Side Meas./Gen. Truth" );
    h_genratio_ssmeas->SetMarkerStyle(25);
    h_genratio_ssmeas->SetMarkerColor(4);
    h_genratio_ssmeas->Divide(hgen_anabin_ratiobin);
    if(debugMode)h_genratio_ssmeas->Print("base");
    
    
    
    std::cout<<"writing output hists to file... "<<std::endl;
    
    hunf->Write();    
    hunf_sameside->Write();

    hfold->Write();    
    hfold_sameside->Write();
    
    h_genratio_oppunf->Write();  
    h_genratio_oppfold->Write();  
    h_genratio_oppmeas->Write();  
    
    h_genratio_ssunf->Write();  
    h_genratio_ssfold->Write();  
    h_genratio_ssmeas->Write();  
    
    h_recratio_oppunf->Write();  
    h_recratio_oppfold->Write();  
    
    h_recratio_ssunf->Write();  
    h_recratio_ssfold->Write();  
    h_recratio_ssgen->Write();  
    
    





    //  drawPDFS -------------------------------------------------- 
    if(!drawPDFs)std::cout<<std::endl<<"NOT drawing PDFs for Bayesian Unfolding!"<<std::endl<<std::endl;
    else{ std::cout<<std::endl<<"drawing output PDFs for Bayesian Unfolding..."<<std::endl;
      std::string outPdfFile=outBayesPdfFile;
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
      
      // canvas for printing
      TCanvas* canvForPrint=new TCanvas("printCanvas","for printing",           1000,1000);    
      canvForPrint->cd();
      canvForPrint->Print(open_outPdfFile.c_str()); 
      
      canvForPrint->cd();
      if(!useSimpBins)canvForPrint->SetLogx(1);
      canvForPrint->SetLogy(1);
      
//      TLine* theLineAtOne_gen= new TLine( boundaries_pt_gen_mat[0]
//					  ,1.
//					  ,(boundaries_pt_gen_mat[nbins_pt_gen_mat] 
//						  + (boundaries_pt_gen_mat[nbins_pt_gen_mat] 
//						     - boundaries_pt_gen_mat[nbins_pt_gen_mat-1]))
//					  ,1.);

      TLine* theLineAtOne_gen= new TLine( 56.,   1.,
					  1000., 1.);
      theLineAtOne_gen->SetLineWidth(1);
      theLineAtOne_gen->SetLineStyle(2);
      theLineAtOne_gen->SetLineColor(36);
      
      TLine* theLineAtOne_reco= new TLine( boundaries_pt_reco_mat[0]  
					   ,1.   
					   ,(boundaries_pt_reco_mat[nbins_pt_reco_mat] 
					     + (boundaries_pt_reco_mat[nbins_pt_reco_mat] 
						- boundaries_pt_reco_mat[nbins_pt_reco_mat-1]))
					   ,1.);
      
      theLineAtOne_reco->SetLineWidth(1);
      theLineAtOne_reco->SetLineStyle(2);
      theLineAtOne_reco->SetLineColor(36);
      
      // hist styles ---------------
      
      hrec_anabin->SetMarkerStyle(kOpenCircle);
      hrec_anabin->SetMarkerColor(kBlue);     
      hrec_anabin->SetMarkerSize(1.02);     
      
      hunf->SetMarkerStyle(kOpenCircle);
      hunf->SetMarkerColor(kRed);
      hunf->SetMarkerSize(1.02);     

      TH1F* hfold_clone=(TH1F*)hfold->Clone("hfold_clone");
      hfold_clone->SetMarkerStyle(kOpenCircle);
      hfold_clone->SetMarkerColor(kGreen);
      hfold_clone->SetMarkerSize(1.02);     
      
      hrec_sameside_anabin->SetMarkerStyle(kOpenSquare);
      hrec_sameside_anabin->SetMarkerColor(kBlue-3);     
      hrec_sameside_anabin->SetMarkerSize(1.02);     
      
      hunf_sameside->SetMarkerStyle(kOpenSquare);
      hunf_sameside->SetMarkerColor(kRed-3);
      hunf_sameside->SetMarkerSize(1.02);     

      TH1F* hfold_sameside_clone=(TH1F*)hfold_sameside->Clone("hfold_sameside_clone");
      hfold_sameside_clone->SetMarkerStyle(kOpenSquare);
      hfold_sameside_clone->SetMarkerColor(kGreen-3);
      hfold_sameside_clone->SetMarkerSize(1.02);     
      
      hgen_anabin->SetMarkerStyle(kOpenStar);
      hgen_anabin->SetMarkerColor(kMagenta);
      hgen_anabin->SetMarkerSize(1.02);     	
      if(useSimpBins)hgen_anabin->SetAxisRange(1.0e-16, 3.0e-4, "Y");
      else hgen_anabin->SetAxisRange(9.0e-15, 3.0e-4, "Y");
      
      
      // gen/meas/unf spectra hists ---------------
      
      hgen_anabin->SetTitle("ppMC Spectra, input and output");
      setupSpectraHist(hgen_anabin, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      
      setupSpectraHist(hrec_anabin	  , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      setupSpectraHist(hfold_clone	  , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      setupSpectraHist(hunf		  , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      
      setupSpectraHist(hrec_sameside_anabin, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      setupSpectraHist(hfold_sameside_clone	  , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      setupSpectraHist(hunf_sameside	  , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      
      hgen_anabin->Draw("E");           
      
      hrec_anabin->Draw("E SAME");           
      hfold_clone->Draw("P E SAME");           
      hunf->Draw("E SAME");           
      
      hrec_sameside_anabin->Draw("E SAME");           
      hfold_sameside_clone->Draw("P E SAME");           
      hunf_sameside->Draw("E SAME");           
      
      TLegend* legend_in = new TLegend( 0.6,0.7,0.9,0.9 );
      
      legend_in->AddEntry(hrec_anabin,          "OppSide Meas.", "p");	
      legend_in->AddEntry(hunf,                    "OppSide Unf.", "p");
      legend_in->AddEntry(hfold_clone,                    "OppSide Fold", "p");
      
      legend_in->AddEntry(hrec_sameside_anabin, "SameSide Meas.", "p");
      legend_in->AddEntry(hunf_sameside,           "SameSide Unf.", "p");
      legend_in->AddEntry(hfold_sameside_clone,           "SameSide Fold", "p");
      
      legend_in->AddEntry(hgen_anabin,          "SameSide Truth", "p");
      
      legend_in->Draw();
      
      canvForPrint->Print(outPdfFile.c_str());      
      
      //---------------
      
      canvForPrint->cd();
      if(!useSimpBins)canvForPrint->SetLogx(1);
      canvForPrint->SetLogy(0);
      
      h_genratio_oppunf->SetTitle( "Ratios w/ Gen. Truth" );
      setupRatioHist(h_genratio_oppunf, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
      
      h_genratio_oppunf->SetAxisRange(56. ,1000.,"X") ;   h_genratio_oppunf->SetAxisRange(  0.7 ,1.6,"Y") ;
      h_genratio_oppmeas->SetAxisRange(56.,1000.,"X") ;	  h_genratio_oppmeas->SetAxisRange( 0.7 ,1.6,"Y") ;
      h_genratio_ssunf->SetAxisRange(56.  ,1000.,"X") ;	  h_genratio_ssunf->SetAxisRange(   0.7 ,1.6,"Y") ;
      h_genratio_ssmeas->SetAxisRange(56. ,1000.,"X");	  h_genratio_ssmeas->SetAxisRange(  0.7 ,1.6,"Y");
     
      h_genratio_oppunf->Draw();
      h_genratio_oppmeas->Draw("SAME");
      //h_genratio_oppfold->Draw("SAME");

      h_genratio_ssunf->Draw("SAME");
      h_genratio_ssmeas->Draw("SAME");
      //h_genratio_ssfold->Draw("SAME");

      TLegend* legend2 = new TLegend( 0.5,0.75,0.9,0.9 );
      legend2->AddEntry(h_genratio_oppunf, "Unf./Gen. Truth", "p");
      legend2->AddEntry(h_genratio_oppmeas, NULL, "p");
      //legend2->AddEntry(h_genratio_oppfold, NULL, "p");
      legend2->AddEntry(h_genratio_ssunf, "Same Side Unf./Gen. Truth", "p");
      legend2->AddEntry(h_genratio_ssmeas, NULL, "p");
      //legend2->AddEntry(h_genratio_ssfold, NULL, "p");
      
      legend2->Draw();
      
      theLineAtOne_gen->Draw();
      
      canvForPrint->Print(outPdfFile.c_str());
      
      //---------------
      
      canvForPrint->cd();
      if(!useSimpBins)canvForPrint->SetLogx(1);
      canvForPrint->SetLogy(0);      
      
      h_recratio_oppunf->SetTitle("Ratios w/ Meas.");
      setupRatioHist(h_recratio_oppunf, useSimpBins, boundaries_pt_reco_mat, nbins_pt_reco_mat);
      
      h_recratio_oppunf->Draw();
      h_recratio_oppfold->Draw("SAME");
      h_recratio_ssunf->Draw("SAME");
      h_recratio_ssfold->Draw("SAME");
      h_recratio_ssgen->Draw("SAME");
      
      TLegend* legend3 = new TLegend( 0.1,0.8,0.4,0.9 );
      legend3->AddEntry(h_recratio_oppunf, NULL, "p");
      legend3->AddEntry(h_recratio_oppfold, NULL, "p");
      legend3->AddEntry(h_recratio_ssunf, NULL, "p");
      legend3->AddEntry(h_recratio_ssfold, NULL, "p");
      legend3->AddEntry(h_recratio_ssgen, NULL, "p");
      legend3->Draw();
      
      theLineAtOne_reco->Draw();
      
      canvForPrint->Print(outPdfFile.c_str());
      
      //---------------      
      
      canvForPrint->Print(close_outPdfFile.c_str());      
      
    }// end draw pdfs
    std::cout<<std::endl<<"draw pdfs for Bayesian Unfolding done!"<<std::endl;

  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();
  
  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_MC->Close();
  return 0;

} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=5 ){
    std::cout<<"do ./bayesUnfoldMCSpectra.exe <targMCDir> <baseOutputName> <doJetID> <useSimpleBins>"<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;

    return rStatus;  }
  
  rStatus=1; // runtime error
  
  rStatus=bayesUnfoldMCSpectra(  (const std::string)argv[1], (const std::string)argv[2], 
			    (int)std::atoi(argv[3]) , (int)std::atoi(argv[4]) 			    ); 
  
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}














  //std::cout<<"underflow="<<hrec_anabin->GetBinContent(0)<<std::endl;
  //std::cout<<"overflow="<<hrec_anabin->GetBinContent(1001)<<std::endl;
  //float _56to64GeVSum = 0.;
  //for(int i=57 ; i < 65; i++) {
  //  std::cout<<"i="<<i<<std::endl;
  //  std::cout<<"bin center ="<< hrec_anabin->GetBinCenter(i) <<std::endl;
  //  std::cout<<"bin content ="<< hrec_anabin->GetBinContent(i) <<std::endl;
  //  _56to64GeVSum += hrec_anabin->GetBinContent(i);
  //  std::cout<<"Sum now ="<< _56to64GeVSum<<std::endl;
  //
  //}
  //std::cout<<"before rebinning: 56-64GeV's contents are:"<< _56to64GeVSum <<std::endl;//57 because 1 GeV bins-> 56-57 is 57'th bin because 0-1 is the first bin
  //std::cout<<"before rebinning: first bin's width is:"    << hrec_anabin->GetBinWidth(57)<<std::endl;
  //std::cout<<"after rebinning: first bin's center is    :"<< hrec_anabin->GetBinCenter(1)<<std::endl;
  //std::cout<<"after rebinning: first bin's contents are :"<< hrec_anabin->GetBinContent(1)<<std::endl;
  //std::cout<<"after rebinning: first bin's width is     :"<< hrec_anabin->GetBinWidth(1)<<std::endl;








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




    //std::cout<<"GetSumw2N="<<h_recratio_ssgen->TH1::GetSumw2N()<<std::endl<<std::endl; 



      //TH1F* hgen_anabin_forDiv=(TH1F*)hgen_anabin->Clone("hgen_anabin_clone_forSVDdiv");
      //hgen_anabin_forDiv=(TH1F*)hgen_anabin_forDiv->Rebin(nbins_pt_reco, "hgen_anabin_clone_forSVDdiv_recobins" , boundaries_pt_reco);
      //hratio_svd[kr]->Divide(hgen_anabin_forDiv);



      //hgen_unfolded_ratio[kr] = (TH1F*)hgen_unfolded_ratio[kr]->Rebin(nbins_pt_reco, ("ppMC_SVDUnf_Ratio_rebin4div"+kRegRandEtaRange).c_str() , boundaries_pt_reco);
      //if(clearOverUnderflows)TH1clearOverUnderflows((TH1*)hgen_unfolded_ratio[kr]);

//hrec_unfolded_ratio[kr]->Print("base");
      //hrec_unfolded_ratio[kr] = (TH1F*)hrec_unfolded_ratio[kr]->Rebin(nbins_pt_reco, ("ppMC_SVDUnf_Ratio_rebin4div"+kRegRandEtaRange).c_str() , boundaries_pt_reco);
      //if(clearOverUnderflows)TH1clearOverUnderflows((TH1*)hrec_unfolded_ratio[kr]);
      //hrec_folded_ratio[kr]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_reco[nbins_pt_reco]);      
      //hrec_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_reco[nbins_pt_reco])
      //hrec_unfolded_ratio[kr]->SetTitle( ("SVD Unf./Meas.,"+kRegRandEtaRange_plotTitle).c_str() );



      //hgen_folded_ratio[kr] = (TH1F*)hfold_svd[kr]->Clone( ("ppMC_SVDUnf_reFold_Ratio"+kRegRandEtaRange).c_str() );
      //hgen_folded_ratio[kr]->Print("base");
      //hgen_folded_ratio[kr]->Divide(hgen_anabin);

      //hgen_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("ppMC_SVDUnf_Ratio"+kRegRandEtaRange).c_str());
      //hgen_unfolded_ratio[kr]->Print("base");
      //hgen_unfolded_ratio[kr]->Divide(hgen_anabin);


//std::cout<<"GetSumw2N="<<hunf->TH1::GetSumw2N()<<std::endl<<std::endl ;



//  if(zeroBins){
//    TH1zeroBins((TH1*)hrec_anabin,Nbins2Clear_reco_lowpt, Nbins2Clear_reco_highpt);
//    hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned_divBinWidth_noOverUnderFlows_zeroBins");
//    hrec_anabin->Print("base");  
//    //std::cout<<"GetSumw2N="<<hrec_anabin->TH1::GetSumw2N()<<std::endl<<std::endl;    
//  }
  
  
 
 
