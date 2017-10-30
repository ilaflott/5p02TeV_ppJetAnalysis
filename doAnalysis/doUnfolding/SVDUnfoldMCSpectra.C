#include "unfoldSpectra.h"

//other settings
const int kRegDraw  = 1 ; // array entries w/ arguments 0-8. 4 -> middle hist on 3x3 SVDplot
const bool drawPDFs=true; 
const bool debugMode=false;



// CODE --------------------------------------------------
int SVDUnfoldMCSpectra( std::string inFile_MC_dir , const std::string baseName , 
			const bool doJetID=true , 
			const bool useSimpBins=false ,
			const int kRegCenter= 5 
			){
  
  
  
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
  std::string outSVDPdfFile   =  outFileName+".pdf"; 
  std::string outSSSVDPdfFile   =  outFileName+"_SS.pdf"; 
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
  
  hrec->Scale(1./4.); // |y| bin width
  hrec->Write( (histTitle).c_str() );
  if(debugMode)hrec->Print("base");
  
  
  histTitle+="_clone";
  TH1F *hrec_anabin = (TH1F*)hrec->Clone( (histTitle).c_str() );
  hrec_anabin->Write((histTitle).c_str());
  if(debugMode)hrec_anabin->Print("base");
  
  std::cout<<"rebinning hrec..."<<std::endl;
  histTitle+="_anabins";
  hrec_anabin = (TH1F*)hrec_anabin->Rebin( nbins_pt_reco, (histTitle).c_str() , boundaries_pt_reco);
  hrec_anabin->Write((histTitle).c_str());   
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
  hrec_sameside_anabin->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside_anabin->Print("base");
  
  std::cout<<"rebinning hrec_sameside..."<<std::endl;
  histTitle2+="_anabins";
  hrec_sameside_anabin = (TH1F*)hrec_sameside_anabin->Rebin( nbins_pt_reco, (histTitle2).c_str() , boundaries_pt_reco);
  hrec_sameside_anabin->Write(histTitle2.c_str());   
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
  
  
  TH1F* hrec_sameside_resp_anabin_empty= new TH1F( ("hpp_rec_sameside_response_anabin_empty"+RandEtaRange).c_str(),"", 
						    nbins_pt_reco, boundaries_pt_reco); 
  
  


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
  hgen->Write(genHistTitle.c_str());
  if(debugMode)hgen->Print("base");    

  genHistTitle+="_divByetabin";
  hgen->Scale(1./4.); // eta bin width for 0.<|y|<2.
  hgen->Write( genHistTitle.c_str());
  if(debugMode)hgen->Print("base");
  
  genHistTitle+="_clone";
  TH1F* hgen_anabin = (TH1F*)hgen->Clone( (genHistTitle).c_str() );
  hgen_anabin->Write(genHistTitle.c_str());
  if(debugMode)hgen_anabin->Print("base");
  
  genHistTitle+="_anabins";
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins_pt_gen, (genHistTitle).c_str() , boundaries_pt_gen);
  hgen_anabin->Write();
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
  
  
  TH1F* hgen_resp_anabin_empty= new TH1F( ("hpp_gen_response_anabin_empty"+RandEtaRange).c_str() ,"", 
					  nbins_pt_gen, boundaries_pt_gen);  
  
  
  
  
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
  hmat->Write(TH2_title.c_str());  
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
    
    matStylePrint(hmat_percenterrs, "ppMC Resp Matrix %errs", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // col normd matrix in binning of interest  ---------------
    
    matStylePrint(hmat_anabin_colnormd, "ppMC Resp Matrix, Columns Normalized", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // row normd matrix in binning of interest  ---------------
    
    matStylePrint(hmat_anabin_rownormd, "ppMC Resp Matrix, Rows Normalized", tempCanvForPdfPrint, outPdfFile, useSimpBins);
    
    // close file 
    
    tempCanvForPdfPrint->Print(close_outPdfFile.c_str());   }
 






  
  
    // setup RooUnfoldResponse, and histogram arrays
    std::cout<<"calling RooUnfoldResponse"<<std::endl;
    //RooUnfoldResponse roo_resp_inv(0, 0, hmat_anabin, ("Response_matrix_anabin_nullresp_"+radius).c_str());    
    RooUnfoldResponse roo_resp_inv(hrec_sameside_resp_anabin_empty, hgen_resp_anabin_empty, hmat_anabin, ("Response_matrix_anabin_nullresp_"+radius).c_str());    
    roo_resp_inv.UseOverflow(doOverUnderflows);

    RooUnfoldResponse roo_resp(hrec_sameside_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix_anabin_"+radius).c_str());    
    roo_resp.UseOverflow(doOverUnderflows);
    
    
    if(debugMode)std::cout<<"creating histos/arrays/canvases for SVD unfolding..."<<std::endl;
    TCanvas *cSpectra               = new TCanvas("cSpectra","",           1500, 1500);      cSpectra->Divide(3,3);	   
    TCanvas *cRatio_gen             = new TCanvas("cRatio_gen","",             1500, 1500);      cRatio_gen->Divide(3,3);      		
    TCanvas *cRatio_rec             = new TCanvas("cRatio_rec","",             1500, 1500);      cRatio_rec->Divide(3,3);      		
    TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","", 1500, 1500);      cPearsonMatrixIter->Divide(3,3);    
    TCanvas *c11                = new TCanvas("c11"," Singular Values and divectors", 1500, 1500);      c11->Divide(2);                    
    
    
    if(debugMode)std::cout<<"cloning input hists..."<<std::endl;
    TH1F* hrec_anabin_clone=(TH1F*)hrec_anabin->Clone("hrec_anabin_clone");
    setupSpectraHist(hrec_anabin_clone, useSimpBins, boundaries_pt_reco, nbins_pt_reco);
    hrec_anabin_clone->SetMarkerStyle(kOpenCircle);
    hrec_anabin_clone->SetMarkerColor(kBlue);
    if(debugMode)hrec_anabin_clone->Print("base");
    
    TH1F* hgen_anabin_clone=(TH1F*)hgen_anabin->Clone("hgen_anabin_clone");
    setupSpectraHist(hgen_anabin_clone, useSimpBins, boundaries_pt_gen, nbins_pt_gen);
    hgen_anabin_clone->SetMarkerStyle(kOpenStar);
    hgen_anabin_clone->SetMarkerColor(kMagenta);
    if(debugMode)hgen_anabin_clone->Print("base");
    
    
    // spectra
    TH1F *hunf_svd[nKregMax]={};        
    TH1F *hfold_svd[nKregMax]={};     
    
    // ratios 
    TH1F *hgen_unfolded_ratio[nKregMax]={};
    TH1F *hrec_unfolded_ratio[nKregMax]={};
    TH1F *hgen_folded_ratio[nKregMax]={};  
    TH1F *hrec_folded_ratio[nKregMax]={};     

    TH2D *hPearsonSVD[nKregMax];      
    


    //RooUnfoldResponse roo_resp_ss_inv(0, 0, hmat_anabin, ("Response_matrix_anabin_nullresp_ss"+radius).c_str());    
    RooUnfoldResponse roo_resp_ss_inv(hrec_sameside_resp_anabin_empty, hgen_resp_anabin_empty, hmat_anabin, ("Response_matrix_anabin_nullresp_ss"+radius).c_str());    
    roo_resp_inv.UseOverflow(doOverUnderflows);

    RooUnfoldResponse roo_resp_ss(hrec_sameside_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix_anabin_ss"+radius).c_str());    
    roo_resp.UseOverflow(doOverUnderflows);



    if(debugMode)std::cout<<"creating histos/arrays/canvases for same-side SVD unfolding..."<<std::endl;
    TCanvas *cSpectra_ss               = new TCanvas("cSpectra_ss","",           1500, 1500);      cSpectra_ss->Divide(3,3);	   
    TCanvas *cRatio_gen_ss            = new TCanvas("cRatio_gen_ss","",             1500, 1500);      cRatio_gen_ss->Divide(3,3);      	    
    TCanvas *cRatio_rec_ss             = new TCanvas("cRatio_rec_ss","",             1500, 1500);      cRatio_rec_ss->Divide(3,3);          	
    TCanvas *cPearsonMatrixIter_ss = new TCanvas("cPearsonMatrixIter_ss","", 1500, 1500);      cPearsonMatrixIter_ss->Divide(3,3);   
    TCanvas *c11_ss                = new TCanvas("c11_ss","S.S. Singular Values and divectors", 1500, 1500);      c11_ss->Divide(2);          
    
    TH1F* hrec_sameside_anabin_clone=(TH1F*)hrec_sameside_anabin->Clone("hrec_sameside_anabin_clone");
    setupSpectraHist(hrec_sameside_anabin_clone, useSimpBins, boundaries_pt_reco, nbins_pt_reco);
    hrec_sameside_anabin_clone->SetMarkerStyle(kOpenSquare);
    hrec_sameside_anabin_clone->SetMarkerColor(kBlue-3);
    if(debugMode)hrec_sameside_anabin_clone->Print("base");
    
    // spectra
    TH1F *hunf_ss_svd[nKregMax]={};        
    TH1F *hfold_ss_svd[nKregMax]={};     
    
    // ratios 
    TH1F *hgen_ss_unfolded_ratio[nKregMax]={};
    TH1F *hrec_ss_unfolded_ratio[nKregMax]={};
    TH1F *hgen_ss_folded_ratio[nKregMax]={};  
    TH1F *hrec_ss_folded_ratio[nKregMax]={};     
    
    TH2D *hPearsonSVD_SS[nKregMax];      
    
    
    
    
    
    // prep for svd unfolding loop
    if(debugMode) std::cout<<std::endl<<"initializing kReg parameter array w/ "<<nKregMax<<" elements"<<std::endl;
    if(debugMode) std::cout<<kRegCenter-kRegRange<<" <= kReg[i] <= "<<kRegCenter+kRegRange<<std::endl;
    
    int kReg[nKregMax]={0};
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) kReg[i+kRegRange]=kRegCenter+i;
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) 
      std::cout<<"kReg["<<i+kRegRange<<"]="<<kReg[i+kRegRange]<<std::endl;
    
    TLegend *leg[nKregMax],*leg1[nKregMax], *leg2[nKregMax];      
    TLegend *leg_ss[nKregMax],*leg_ss1[nKregMax], *leg_ss2[nKregMax];      
    
    TLine* theGenLine= new TLine(boundaries_pt_gen_mat[0]
				 ,1.
				 ,(boundaries_pt_gen_mat[nbins_pt_gen_mat] 
				   + (boundaries_pt_gen_mat[nbins_pt_gen_mat] 
				      - boundaries_pt_gen_mat[nbins_pt_gen_mat-1]))
				 ,1.);
    theGenLine->SetLineWidth(1);
    theGenLine->SetLineStyle(2);
    theGenLine->SetLineColor(36);    
    
    TLine* theRecoLine= new TLine(boundaries_pt_reco_mat[0]
				  ,1.   
				  ,(boundaries_pt_reco_mat[nbins_pt_reco_mat]) ,1.);
				    //				    + (boundaries_pt_reco_mat[nbins_pt_reco_mat] 
				    //				       - boundaries_pt_reco_mat[nbins_pt_reco_mat-1]))
				    //,1.);
    theRecoLine->SetLineWidth(1);
    theRecoLine->SetLineStyle(2);
    theRecoLine->SetLineColor(36);    


    // SVD unfolding loop
    for(int kr = 0; kr < nKregMax; ++kr){
      if(debugMode&&kr==0)    std::cout<<"beginning SVD Unfolding Loop..."<<std::endl;
      std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl;
      
      std::string kRegRandEtaRange_plotTitle=" kReg"+std::to_string(kReg[kr]);//+RandEtaRange_plotTitle;
      std::string kRegRandEtaRange="_kReg"+std::to_string(kReg[kr]);//+RandEtaRange;
         

      if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd..."<<std::endl<<std::endl;
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      
      
      if(debugMode)std::cout<<"applying roo_resp to histo hunf_svd[kr="<<kr<<"]..."<<std::endl;
      //hfold_svd[kr] = (TH1F*)roo_resp.ApplyToTruth(hunf_svd[kr]);
      hfold_svd[kr] = (TH1F*)roo_resp_inv.ApplyToTruth(hunf_svd[kr]);
      
      
      // get covariance matrix and calculate pearson coefficients
      if(debugMode)std::cout<<"calling Ereco..."<<std::endl;
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);
      
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat, false );      
      hPearsonSVD[kr] = new TH2D (*pearson);
      hPearsonSVD[kr]->SetName(("pearson_oppside_"+kRegRandEtaRange).c_str());
      if(debugMode)hPearsonSVD[kr]->Print("base");
      //if(debugMode)std::cout<<"creating \"rebinned\" pearson matrix..."<<std::endl;
      //hPearsonSVD[kr] = reBinPearsonTH2( pearson );
      
      // 3x3 spectra canvases  ---------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
      cSpectra->cd(kr+1);
      cSpectra->cd(kr+1)->SetLogy();
      if(!useSimpBins)cSpectra->cd(kr+1)->SetLogx(1);
      
      hunf_svd[kr]->SetTitle( ("Reco Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );
      hunf_svd[kr]->SetName( ("hunf_svd"+kRegRandEtaRange).c_str() );
      hunf_svd[kr]->SetMarkerStyle(kOpenCircle);
      hunf_svd[kr]->SetMarkerColor(kRed);     //hunf_svd[kr]->SetLineStyle(33);      //hunf_svd[kr]->SetLineColor(kBlue);
      if(useSimpBins)hunf_svd[kr]->SetAxisRange(1.0e-16, 3.0e-4, "Y");
      else hunf_svd[kr]->SetAxisRange(9.0e-15, 3.0e-4, "Y");

      if(debugMode)hunf_svd[kr]->Print("base");      
      
      hfold_svd[kr]->SetName( ("hfold_svd"+kRegRandEtaRange).c_str() );
      hfold_svd[kr]->SetMarkerStyle(kOpenCircle);
      hfold_svd[kr]->SetMarkerColor(kGreen);      //hfold_svd[kr]->SetLineColor(kRed);      
      if(debugMode)hfold_svd[kr]->Print("base");
      
      hunf_svd[kr]->Draw("P");
      hfold_svd[kr]->Draw("P SAME");      
      hgen_anabin_clone->Draw("P SAME");
      hrec_anabin_clone->Draw("P SAME");
      
      leg[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option 
      leg[kr]->AddEntry(hunf_svd[kr],"OppSide Unf.","p");
      leg[kr]->AddEntry(hfold_svd[kr],"Oppside Fold","p");
      leg[kr]->AddEntry(hgen_anabin_clone,"SameSide Truth","p");
      leg[kr]->AddEntry(hrec_anabin_clone,"OppSide Meas.","p");
      leg[kr]->SetTextSize(0.02);
      leg[kr]->Draw();


      // 3x3 unfORfold/meas canvases  ---------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_rec canvas..."<<std::endl;
      cRatio_rec->cd(kr+1);
      if(!useSimpBins)cRatio_rec->cd(kr+1)->SetLogx(1);
      
      hrec_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("hrec_unf_ratio"+kRegRandEtaRange).c_str());
      hrec_unfolded_ratio[kr] = (TH1F*)hrec_unfolded_ratio[kr]->Rebin( nbins_pt_reco, ("hrec_unf_ratio"+kRegRandEtaRange+"_rebin").c_str() , boundaries_pt_reco);      
      hrec_unfolded_ratio[kr]->Divide(hrec_anabin);
      hrec_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
      hrec_unfolded_ratio[kr]->SetTitle( ("SVD/Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
      if(debugMode)hrec_unfolded_ratio[kr]->Print("base");
      
      
      hrec_folded_ratio[kr] = (TH1F*)hfold_svd[kr]->Clone( ("hrec_fold_ratio"+kRegRandEtaRange).c_str() );
      hrec_folded_ratio[kr] = (TH1F*)hrec_folded_ratio[kr]->Rebin( nbins_pt_reco, ("hrec_fold_ratio"+kRegRandEtaRange+"_rebin").c_str() , boundaries_pt_reco);      
      hrec_folded_ratio[kr]->Divide(hrec_anabin);
      hrec_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_folded_ratio[kr]->SetMarkerColor(kRed);
      if(debugMode)hrec_folded_ratio[kr]->Print("base");
      
      hrec_unfolded_ratio[kr]->Draw();
      hrec_folded_ratio[kr]->Draw("SAME");
      
      leg1[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option 
      leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],"Unf./Meas.","pl");
      leg1[kr]->AddEntry(hrec_folded_ratio[kr],"reFold/Meas.","pl");
      leg1[kr]->SetTextSize(0.02); 
      leg1[kr]->Draw();
      
      theRecoLine->Draw();



      // draw on cRatio_gen canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_gen canvas..."<<std::endl;
      cRatio_gen->cd(kr+1);
      if(!useSimpBins)cRatio_gen->cd(kr+1)->SetLogx(1);

      hgen_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("hgen_unf_ratio"+kRegRandEtaRange).c_str());
      hgen_unfolded_ratio[kr]->Divide(hgen_anabin);      
      hgen_unfolded_ratio[kr]->SetMarkerStyle(33);
      hgen_unfolded_ratio[kr]->SetMarkerColor(kRed);
      if(debugMode)hgen_unfolded_ratio[kr]->Print("base");

      hgen_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hgen_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
      hgen_unfolded_ratio[kr]->SetTitle( ("SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      hgen_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");      
      
      TH1F* hgen_anabin_ratClone=(TH1F*)hgen_anabin->Clone("hgen_anabin_ratioClone");
      hgen_anabin_ratClone=(TH1F*)hgen_anabin_ratClone->Rebin( nbins_pt_reco, ("hgen_anabin_ratioClone_rebin") , boundaries_pt_reco);      
      
      hgen_folded_ratio[kr] = (TH1F*)hfold_svd[kr]->Clone( ("hgen_fold_ratio"+kRegRandEtaRange).c_str());            
      hgen_folded_ratio[kr]->Divide(hgen_anabin_ratClone);      
      hgen_folded_ratio[kr]->SetMarkerStyle(24);
      hgen_folded_ratio[kr]->SetMarkerColor(kBlue);
      if(debugMode)hgen_folded_ratio[kr]->Print("base");      
      hgen_folded_ratio[kr]->SetTitle( ("SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      
      hgen_unfolded_ratio[kr]->Draw();
      hgen_folded_ratio[kr]->Draw("same");
      
      
      leg2[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg2[kr]->AddEntry(hgen_unfolded_ratio[kr],"Unf./Gen","pl");
      leg2[kr]->AddEntry(hgen_folded_ratio[kr],"Folded/Gen","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2[kr]->SetTextSize(0.02); 
      leg2[kr]->Draw();
      
      theGenLine->Draw();








      // same-side unfolding/refolding --------------
      std::cout<<std::endl<<"calling RooUnfoldSvd and Hreco for sameside MC test..."<<std::endl<<std::endl;
      RooUnfoldSvd unf_ss_svd(&roo_resp_ss, hrec_sameside_anabin, kReg[kr]);
      hunf_ss_svd[kr] = (TH1F*)unf_ss_svd.Hreco(errorTreatment);
      std::cout<<std::endl;         
      
      std::cout<<"applying roo_resp to sameside histo hunf_ss_svd[kr="<<kr<<"]..."<<std::endl;
      //hfold_ss_svd[kr] = (TH1F*)roo_resp.ApplyToTruth(hunf_ss_svd[kr]);
      hfold_ss_svd[kr] = (TH1F*)roo_resp_ss_inv.ApplyToTruth(hunf_ss_svd[kr]);
      
      //////////////////////////////////////////////
      //sameside covariance/pearson matrices??? maybe, we will see.
      TMatrixD covmat_ss = unf_ss_svd.Ereco(errorTreatment);
      if(debugMode)std::cout<<std::endl;
      
      TMatrixD *pearson_ss = CalculatePearsonCoefficients(&covmat_ss, false );
      if(debugMode)std::cout<<std::endl;
      
      std::cout<<"creating TH2 for pearson matrix..."<<std::endl;
      hPearsonSVD_SS[kr] = new TH2D (*pearson_ss);
      hPearsonSVD_SS[kr]->SetName(("pearson_sameside_"+kRegRandEtaRange).c_str());
      if(debugMode)hPearsonSVD_SS[kr]->Print("base");
      if(debugMode)std::cout<<std::endl;
      //sameside covariance/pearson matrices???
      /////////////////////////////////////////////
      
      // 3x3 spectra canvases  ------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
      cSpectra_ss->cd(kr+1);
      cSpectra_ss->cd(kr+1)->SetLogy();
      if(!useSimpBins)cSpectra_ss->cd(kr+1)->SetLogx(1);            
      
      //hrec_sameside_anabin_clone->SetTitle( ("S.S. Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );

      hunf_ss_svd[kr]->SetTitle( ("S.S. Spectra,"+kRegRandEtaRange_plotTitle).c_str() );
      hunf_ss_svd[kr]->SetName( ("hunf_ss_svd"+kRegRandEtaRange).c_str() );
      hunf_ss_svd[kr]->SetLineStyle(33);
      hunf_ss_svd[kr]->SetLineColor(kBlue);
      if(debugMode)hunf_ss_svd[kr]->Print("base");
      
      hfold_ss_svd[kr]->SetName( ("hfold_ss_svd"+kRegRandEtaRange).c_str() );
      hfold_ss_svd[kr]->SetLineColor(kRed);      
      if(debugMode)hfold_ss_svd[kr]->Print("base");
      
      
      hunf_ss_svd[kr]->Draw();
      hrec_sameside_anabin_clone->Draw("same");
      hfold_ss_svd[kr]->Draw("same");
      
      
      leg_ss[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option
      leg_ss[kr]->AddEntry(hrec_sameside_anabin_clone,"S.S. Measured","pl");
      leg_ss[kr]->AddEntry(hunf_ss_svd[kr],"S.S. Unfolded","pl");
      leg_ss[kr]->AddEntry(hfold_ss_svd[kr],"S.S. Folded","pl");
      leg_ss[kr]->SetTextSize(0.02);
      leg_ss[kr]->Draw();





      // 3x3 unfORfold/meas canvases  ---------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_rec canvas..."<<std::endl;
      cRatio_rec_ss->cd(kr+1);
      if(!useSimpBins)cRatio_rec_ss->cd(kr+1)->SetLogx(1);
      
      
      
      hrec_ss_folded_ratio[kr] = (TH1F*)hfold_ss_svd[kr]->Clone( ("hrec_ss_fold_ratio"+kRegRandEtaRange).c_str() );
      hrec_ss_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_ss_folded_ratio[kr]->SetMarkerColor(kRed);
      hrec_ss_folded_ratio[kr]->Divide(hrec_anabin);
      if(debugMode)hrec_ss_folded_ratio[kr]->Print("base");
      
      
      
      //assert(false);
      
      hrec_ss_unfolded_ratio[kr] = (TH1F*)hunf_ss_svd[kr]->Clone( ("hrec_ss_unf_ratio"+kRegRandEtaRange).c_str());
      hrec_ss_unfolded_ratio[kr] = (TH1F*)hrec_ss_unfolded_ratio[kr]->Rebin (nbins_pt_reco, ("hrec_ss_unf_ratio"+kRegRandEtaRange+"_rebin").c_str(), boundaries_pt_reco);
      hrec_ss_unfolded_ratio[kr]->Divide(hrec_anabin);
      
      hrec_ss_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
      hrec_ss_unfolded_ratio[kr]->SetTitle( ("S.S. SVD/Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_ss_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");      
      hrec_ss_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_ss_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
      if(debugMode)hrec_ss_unfolded_ratio[kr]->Print("base");

      hrec_ss_unfolded_ratio[kr]->Draw();
      hrec_ss_folded_ratio[kr]->Draw("same");

      leg_ss1[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg_ss1[kr]->AddEntry(hrec_ss_unfolded_ratio[kr],"S.S. Unf./Meas.","pl");
      leg_ss1[kr]->AddEntry(hrec_ss_folded_ratio[kr],"S.S. Folded/Meas.","pl");
      leg_ss1[kr]->SetTextSize(0.02); 
      leg_ss1[kr]->Draw();
      
      theRecoLine->Draw();



      // draw on cRatio_gen canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_gen_ss canvas..."<<std::endl;
      cRatio_gen_ss->cd(kr+1);
      if(!useSimpBins)cRatio_gen_ss->cd(kr+1)->SetLogx(1);

      hgen_ss_unfolded_ratio[kr] = (TH1F*)hunf_ss_svd[kr]->Clone( ("hgen_ss_unf_ratio"+kRegRandEtaRange).c_str());
      hgen_ss_unfolded_ratio[kr]->Divide(hgen_anabin);      
      hgen_ss_unfolded_ratio[kr]->SetMarkerStyle(33);
      hgen_ss_unfolded_ratio[kr]->SetMarkerColor(kRed);
      if(debugMode)hgen_ss_unfolded_ratio[kr]->Print("base");

      hgen_ss_unfolded_ratio[kr]->SetTitle( ("S.S. SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      //hgen_unfolded_ratio[kr]->SetMarkerStyle(27);
      //hgen_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
      

      hgen_ss_folded_ratio[kr] = (TH1F*)hfold_ss_svd[kr]->Clone( ("hgen_ss_fold_ratio"+kRegRandEtaRange).c_str());
      hgen_ss_folded_ratio[kr]->Divide(hgen_anabin_ratClone);      
      hgen_ss_folded_ratio[kr]->SetMarkerStyle(24);
      hgen_ss_folded_ratio[kr]->SetMarkerColor(kBlue);
      if(debugMode)      hgen_ss_folded_ratio[kr]->Print("base");
      
      hgen_ss_unfolded_ratio[kr]->SetTitle( ("S.S. SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      hgen_ss_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hgen_ss_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hgen_ss_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
      
      hgen_ss_unfolded_ratio[kr]->Draw();
      hgen_ss_folded_ratio[kr]->Draw("same");
      
      leg_ss2[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg_ss2[kr]->AddEntry(hgen_ss_unfolded_ratio[kr],"S.S. Unf./Gen","pl");
      leg_ss2[kr]->AddEntry(hgen_ss_folded_ratio[kr],"S.S. Folded/Gen","pl");
      leg_ss2[kr]->SetTextSize(0.02); 
      leg_ss2[kr]->Draw();
      
      theGenLine->Draw();
      
      
      
      //  singular values and d_i vector ---------------------------
      //Note that these do not depend on the regularization.
      //The opposite: they tell you which regularization to use! (ian note: how?)
      if(kr == kRegDraw){      

	std::cout << "  getting divectors and singular values for oppSide ... " << std::endl;
	TSVDUnfold *svdUnfold = unf_svd.Impl();
	//svdUnfold->SetNormalize(true);
        //if(debugMode) std::cout<<std::endl<<"    !!! kr==0 !!!"<<std::endl<<std::endl;
	
        // singular values
        if(debugMode)std::cout << "  getting singular values... " << std::endl;
        TH1 *hSVal = (TH1*)svdUnfold->GetSV();
        if(debugMode)
	  for(int bin=1; bin<=hSVal->GetNbinsX(); bin++) std::cout<<"bin: "<<bin<<",  SV: "<<hSVal->GetBinContent(bin)<< std::endl;
	
        if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
        c11->cd(1);
	c11->cd(1)->SetLogy();        //gPad->SetLogy();     
	
        hSVal->SetTitle(" singular values ");
        hSVal->SetXTitle(" singular values ");        
        hSVal->DrawCopy("E");	//hSVal->Draw("E");
	
	// di vector values
        std::cout << "  getting di vector values... " <<  std::endl;
        TH1D *hdi = (TH1D*)svdUnfold->GetD();
        if(debugMode)
	  for(int bin=1; bin<=hdi->GetNbinsX(); bin++)std::cout<<"i: "<<bin<<",  di: "<<hdi->GetBinContent(bin)<<std::endl;
	
	c11->cd(2);
        c11->cd(2)->SetLogy();    //gPad->SetLogy();        //cdi->SetLogy();
        
	hdi->SetTitle(" di vectors ");
        hdi->SetXTitle(" |d_{i}^{kreg}| ");
	hdi->DrawCopy("E"); //hdi->Draw("E");

        //drawText( "5.02 TeV ak4PFJets",0.358173, 0.8459761, 19);
	//drawText( "ppData, Prompt Reco HighPtJets",0.358173, 0.8159761, 19);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);	
        if(debugMode)std::cout<<std::endl<<"done with kr==0 specifics"<<std::endl<<std::endl;      





	std::cout << "  getting divectors and singular values for sameSide ... " << std::endl;
	TSVDUnfold *svdUnfold_ss = unf_ss_svd.Impl();
	//svdUnfold->SetNormalize(true);
        //if(debugMode) std::cout<<std::endl<<"    !!! kr==0 !!!"<<std::endl<<std::endl;
	
        // singular values
        if(debugMode)std::cout << "  getting singular values... " << std::endl;
        TH1 *hSVal_ss = (TH1*)svdUnfold_ss->GetSV();
        if(debugMode)
	  for(int bin=1; bin<=hSVal_ss->GetNbinsX(); bin++) std::cout<<"bin: "<<bin<<",  SV: "<<hSVal_ss->GetBinContent(bin)<< std::endl;
	
        if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
        c11_ss->cd(1);
	c11_ss->cd(1)->SetLogy();        //gPad->SetLogy();     
	
        hSVal_ss->SetTitle(" s.s. singular values ");
        hSVal_ss->SetXTitle(" s.s singular values ");        
        hSVal_ss->DrawCopy("E");	//hSVal->Draw("E");
	
	// di vector values
        std::cout << "  getting di vector values... " <<  std::endl;
        TH1D *hdi_ss = (TH1D*)svdUnfold_ss->GetD();
        if(debugMode)
	  for(int bin=1; bin<=hdi_ss->GetNbinsX(); bin++)std::cout<<"i: "<<bin<<",  di: "<<hdi_ss->GetBinContent(bin)<<std::endl;
	
	c11_ss->cd(2);
        c11_ss->cd(2)->SetLogy();    //gPad->SetLogy();        //cdi->SetLogy();
        
	hdi_ss->SetTitle(" s.s. di vectors ");
        hdi_ss->SetXTitle(" |d_{i}^{kreg}| ");
	hdi_ss->DrawCopy("E"); //hdi->Draw("E");

        //drawText( "5.02 TeV ak4PFJets",0.358173, 0.8459761, 19);
	//drawText( "ppData, Prompt Reco HighPtJets",0.358173, 0.8159761, 19);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);	
        if(debugMode)std::cout<<std::endl<<"done with kr==0 specifics"<<std::endl<<std::endl;      
      }//end drawing singular/divector values
      


      // Pearson Drawing
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter canvas..."<<std::endl;
      cPearsonMatrixIter->cd(kr+1);
      
      hPearsonSVD[kr]->SetName( ("hPearsonSVD"+kRegRandEtaRange).c_str() );      
      hPearsonSVD[kr]->SetAxisRange(-1., 1., "Z");
      
      hPearsonSVD[kr]->GetZaxis()->SetLabelSize(0.035);
      hPearsonSVD[kr]->GetXaxis()->SetTitle("reco bin num");
      hPearsonSVD[kr]->GetYaxis()->SetTitle("gen bin num");
      hPearsonSVD[kr]->SetTitle( ("pearsonMatrix SVDPriorMeas"+kRegRandEtaRange_plotTitle).c_str()  );
      hPearsonSVD[kr]->Draw("COLZ");


      // Pearson Drawing
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter_SS canvas..."<<std::endl;
      cPearsonMatrixIter_ss->cd(kr+1);
      
      hPearsonSVD_SS[kr]->SetName( ("hPearsonSVD SS"+kRegRandEtaRange).c_str() );      
      hPearsonSVD_SS[kr]->SetAxisRange(-1., 1., "Z");
      
      hPearsonSVD_SS[kr]->GetZaxis()->SetLabelSize(0.035);
      hPearsonSVD_SS[kr]->GetXaxis()->SetTitle("reco bin num");
      hPearsonSVD_SS[kr]->GetYaxis()->SetTitle("gen bin num");
      hPearsonSVD_SS[kr]->SetTitle( ("pearsonMatrix SVD SS"+kRegRandEtaRange_plotTitle).c_str()  );
      hPearsonSVD_SS[kr]->Draw("COLZ");
            
    }// kReg loop ends
    
    
    std::cout<<std::endl<<"writing SVD histo arrays to file..."<<std::endl;    
    fout->cd();
    for(int kr = 0; kr<nKregMax; ++kr){
      
      hunf_svd[kr]->Write();            
      hfold_svd[kr]->Write();           
      hrec_unfolded_ratio[kr]->Write();          
      hrec_folded_ratio[kr]->Write();          
      hgen_unfolded_ratio[kr]->Write();          
      hgen_folded_ratio[kr]->Write();          

      hPearsonSVD[kr]->Write();



      hunf_ss_svd[kr]->Write();            
      hfold_ss_svd[kr]->Write();           

      hrec_ss_unfolded_ratio[kr]->Write();          
      hrec_ss_folded_ratio[kr]->Write();          
      hgen_ss_unfolded_ratio[kr]->Write();          
      hgen_ss_folded_ratio[kr]->Write();          

      hPearsonSVD_SS[kr]->Write();
    }
    
    if(!drawPDFs){ if(debugMode)std::cout<<std::endl<<"NOT drawing PDFs for SVD Unfolding"<<std::endl<<std::endl; }
    else{ std::cout<<std::endl<<"drawing PDFs for SVD Unfolding..."<<std::endl;
      {
	// form filename string, open pdf file and draw the canvases we have so far
	std::string outPdfFile=outSVDPdfFile;
	std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
	
	//  already-drawn oppSide canvases --------------------------------------------------
	c11->Print(open_outPdfFile.c_str());
	c11->cd()                ;   c11->Print(outPdfFile.c_str());


	//  cSpectra 3X3 and Check for meas --------------------------------------------------
	cPearsonMatrixIter->cd() ;   cPearsonMatrixIter->Print(outPdfFile.c_str());
	
	TCanvas *cPearsonMatCheck        = new TCanvas("cPearsonMatCheck","",        1500, 1500);    
	
	cPearsonMatCheck->cd();
	
	hPearsonSVD[kRegDraw]->Draw("COLZ");
	
	cPearsonMatCheck->Print(outPdfFile.c_str());
	
	//  cSpectra 3X3 and Check for meas --------------------------------------------------
	cSpectra->cd()           ;   cSpectra->Print(outPdfFile.c_str());
	
	TCanvas *cSpectraCheck_rec        = new TCanvas("cSpectraCheck_rec","",        1500, 1500);    
	cSpectraCheck_rec->cd();
	cSpectraCheck_rec->SetLogy(1);
	if(!useSimpBins)cSpectraCheck_rec->SetLogx(1);
	
	hunf_svd[kRegDraw]->SetTitle("SVD, Spectra");
	
	hunf_svd[kRegDraw]->Draw("P");
	hfold_svd[kRegDraw]->Draw("P SAME");
	hgen_anabin_clone->Draw("P SAME");	
	hrec_anabin_clone->Draw("P SAME");

	
	TLegend * leg1 = new TLegend(0.60, 0.79, 0.90, 0.87, NULL,"NBNDC");
	leg1->AddEntry(hunf_svd[kRegDraw],"OppSide Unf.","p");
	leg1->AddEntry(hfold_svd[kRegDraw],"OppSide Fold","p");
	leg1->AddEntry(hgen_anabin_clone,"SameSide Truth","p");	
	leg1->AddEntry(hrec_anabin_clone,"OppSide Meas","p");
	leg1->SetTextSize(0.02);
	leg1->Draw();
	
	//drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
	//drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
	
	cSpectraCheck_rec->Print(outPdfFile.c_str());
	
	
	//  cRatio 3X3 and Check for meas --------------------------------------------------
	cRatio_rec->cd()         ;   cRatio_rec->Print(outPdfFile.c_str());
	
	TCanvas *cRatioCheck_rec        = new TCanvas("cRatioCheck_rec","",        1500, 1500);    
	cRatioCheck_rec->cd();
	if(!useSimpBins)cRatioCheck_rec->SetLogx(1);
	
	hrec_folded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");
	//	hrec_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hrec_folded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ Meas.");
	
	hrec_unfolded_ratio[kRegDraw]->Draw();
	hrec_folded_ratio[kRegDraw]->Draw("same");
	
	TLegend * leg2 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
	leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Unf./Meas.","pl");
	leg2->AddEntry(hrec_folded_ratio[kRegDraw],"Refold/Meas.","pl");
	//leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
	leg2->SetTextSize(0.02);
	leg2->Draw();
	
	//drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
	//drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
	
	//TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
	//theLine->SetLineWidth(1);
	//theLine->SetLineStyle(2);
	//theLine->SetLineColor(36);

	theRecoLine->Draw();
	
	cRatioCheck_rec->Print(outPdfFile.c_str());
	
	
	
	
	//  cRatio 3X3 Check for gen --------------------------------------------------
	cRatio_gen->cd()         ;   cRatio_gen->Print(outPdfFile.c_str());      
	
	TCanvas *cRatioCheck_gen        = new TCanvas("cRatioCheck_gen","",        1500, 1500);     
	cRatioCheck_gen->cd();
	if(!useSimpBins)cRatioCheck_gen->SetLogx(1);
	
	hgen_folded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");
	//hgen_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hgen_folded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ Gen Truth");
	
	hgen_unfolded_ratio[kRegDraw]->Draw();
	hgen_folded_ratio[kRegDraw]->Draw("same");
	
	TLegend * leg3 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
	leg3->AddEntry(hgen_unfolded_ratio[kRegDraw],"Unf./Gen.","pl");
	leg3->AddEntry(hgen_folded_ratio[kRegDraw],"Refold/Gen.","pl");
	//leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
	leg3->SetTextSize(0.02);
	leg3->Draw();
	
	//drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
	//drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
	
	//TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
	//theLine->SetLineWidth(1);
	//theLine->SetLineStyle(2);
	//theLine->SetLineColor(36);
	theGenLine->Draw();
	
	cRatioCheck_gen->Print(outPdfFile.c_str());
	cRatioCheck_gen->Print(close_outPdfFile.c_str());
      }	
      
      
      {
	// form filename string, open pdf file and draw the canvases we have so far
	std::string outPdfFile=outSSSVDPdfFile;
	std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
	
	
	//  already-drawn sameside canvases --------------------------------------------------
	c11->Print(open_outPdfFile.c_str());
	c11_ss->cd()                ;   c11_ss->Print(outPdfFile.c_str());
	cPearsonMatrixIter_ss->cd() ;   cPearsonMatrixIter_ss->Print(outPdfFile.c_str());
	cSpectra_ss->cd()           ;   cSpectra_ss->Print(outPdfFile.c_str());
	
	
	//  cRatio Check for meas --------------------------------------------------
	cRatio_rec_ss->cd()         ;   cRatio_rec_ss->Print(outPdfFile.c_str());      
	
	TCanvas *cRatioCheck_rec_ss        = new TCanvas("cRatioCheck_rec_ss","",        1500, 1500);    
	cRatioCheck_rec_ss->cd();
	if(!useSimpBins)cRatioCheck_rec_ss->SetLogx(1);
	
	hrec_ss_folded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");
	//hrec_ss_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hrec_ss_folded_ratio[kRegDraw]->SetTitle("S.S. SVD Unf. Ratios w/ Meas.");
	
	hrec_ss_unfolded_ratio[kRegDraw]->Draw();
	hrec_ss_folded_ratio[kRegDraw]->Draw("same");
	
	TLegend * leg_ss2 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
	leg_ss2->AddEntry(hrec_ss_unfolded_ratio[kRegDraw],"Unf./Meas.","pl");
	leg_ss2->AddEntry(hrec_ss_folded_ratio[kRegDraw],"Refold/Meas.","pl");
	leg_ss2->SetTextSize(0.02);
	leg_ss2->Draw();
	
	//drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
	//drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
	
	//TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
	//theLine->SetLineWidth(1);
	//theLine->SetLineStyle(2);
	//theLine->SetLineColor(36);
	theRecoLine->Draw();
	
	cRatioCheck_rec_ss->Print(outPdfFile.c_str());
	
	
	//  cRatio Check for gen --------------------------------------------------
	cRatio_gen_ss->cd()         ;   cRatio_gen_ss->Print(outPdfFile.c_str());
	
	TCanvas *cRatioCheck_gen_ss        = new TCanvas("cRatioCheck_gen_ss","",        1500, 1500);    
	cRatioCheck_gen_ss->cd();
	if(!useSimpBins)cRatioCheck_gen_ss->SetLogx(1);
	
	hgen_ss_folded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");//SetAxisRange(0.5, 1.5, "Y");
	//hgen_ss_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hgen_ss_folded_ratio[kRegDraw]->SetTitle("S.S. SVD Unf. Ratios w/ Gen Truth");
	
	hgen_ss_unfolded_ratio[kRegDraw]->Draw();
	hgen_ss_folded_ratio[kRegDraw]->Draw("same");
	
	TLegend * leg_ss3 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
	leg_ss3->AddEntry(hgen_ss_unfolded_ratio[kRegDraw],"Unf./Gen.","pl");
	leg_ss3->AddEntry(hgen_ss_folded_ratio[kRegDraw],"Refold/Gen.","pl");
	leg_ss3->SetTextSize(0.02);
	leg_ss3->Draw();
	
	//drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
	//drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
	drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
	drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
	
	//TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
	//theLine->SetLineWidth(1);
	//theLine->SetLineStyle(2);
	//theLine->SetLineColor(36);
	theGenLine->Draw();
	
	cRatioCheck_gen_ss->Print(outPdfFile.c_str());
	
	cRatioCheck_gen_ss->Print(close_outPdfFile.c_str());
      }
      std::cout<<std::endl<<"done drawing SVD PDFs!"<<std::endl<<std::endl;
      
    }// end drawPDFs
    
  
  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();
  
  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_MC->Close();
  return 0;

} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=6 ){
    std::cout<<"do ./SVDUnfoldMCSpectra.exe <targMCDir> <baseOutputName> <doJetID> <useSimpleBins> <kRegCenter>"<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;

    return rStatus;  }
  
  rStatus=1; // runtime error
  
  rStatus=SVDUnfoldMCSpectra(  (const std::string)argv[1], (const std::string)argv[2], 
			       (int)std::atoi(argv[3]) , (int)std::atoi(argv[4]) ,  (int)std::atoi(argv[5]) 
			       ); 
  
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
  
  
 
 
