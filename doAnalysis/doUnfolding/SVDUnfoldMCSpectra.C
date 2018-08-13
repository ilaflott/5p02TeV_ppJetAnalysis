#include "unfoldSpectra.h"

//other settings
const int kRegDraw  = 4 ; // array entries w/ arguments 0-8. 4 -> middle hist on 3x3 SVDplot
const int kRegDrawSS = 1;

const bool drawPDFs=true; 
const bool debugMode=false;
const bool drawRespMatrix=false;

// CODE --------------------------------------------------
int SVDUnfoldMCSpectra( std::string inFile_MC_dir , const std::string baseName   ,
			const bool doJetID=true   , const bool useSimpBins=false , const int kRegCenter= 5 ){
  
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

  // ppMC input file -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
  std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  
  
  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldMCSpectra_outdir+fullJetType;
  if(doJetID)outFileName+="_wjtID";
  if(!useSimpBins)outFileName+="_anabins";
  else outFileName+="_simpbins";
  outFileName+="_"+baseName;
  
  std::string outRespMatPdfFile =  outFileName+"_respMat.pdf";//add SS to outfile name itself
  std::string outSVDPdfFile     =  outFileName+".pdf";  
  std::string outRootFile       =  outFileName+".root";  
  std::string out3x3SVDPdfFile  =  outFileName+"_3x3.pdf"; 
  
  checkNRenameFiles( outFileName, 
		     &outRespMatPdfFile, 
		     &outSVDPdfFile, 
		     &outRootFile ,
		     &out3x3SVDPdfFile);
  
  //for output
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE");   
  if(debugMode)fout->cd();  
  


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
  
  // ---------- reco, measured spectra to unfold
  //std::string histTitle="hJetQA";
  std::string histTitle="hpp_mcclosure_reco_test";
  if(doJetID)histTitle+="_wJetID";
  //else histTitle+="_0wJetID";
  //histTitle+="_jtpt";
  histTitle+=RandEtaRange;
  if(debugMode)std::cout<<"histTitle="<<histTitle<<std::endl;

  TH1D*  hrec = (TH1D*)fpp_MC->Get( histTitle.c_str() ); 
  if(debugMode)hrec->Print("base");
  if(debugMode)hrec->Write();
  
  
  histTitle+="_clone";
  TH1D *hrec_rebin = (TH1D*)hrec->Clone( (histTitle).c_str() );
  if(debugMode)hrec_rebin->Print("base");
  if(debugMode)hrec_rebin->Write(histTitle.c_str());
  
  std::cout<<"rebinning hrec..."<<std::endl;
  histTitle+="_rebins";
  hrec_rebin = (TH1D*)hrec_rebin->Rebin( nbins_pt_reco, (histTitle).c_str() , boundaries_pt_reco);
  if(debugMode)hrec_rebin->Print("base");  
  if(debugMode)hrec_rebin->Write(histTitle.c_str());   
  
  if(clearOverUnderflows){
    histTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_rebin);
    if(debugMode)hrec_rebin->Print("base");  
    if(debugMode)hrec_rebin->Write(histTitle.c_str());
  }
  
  //cosmetics
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetLineColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02);     
  
  
  
  // ---------- reco, measured same-side spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="hpp_mcclosure_reco";
  if(doJetID)histTitle2+="_wJetID";
  histTitle2+=RandEtaRange;
  
  TH1D*  hrec_sameside = (TH1D*)fpp_MC->Get( histTitle2.c_str() );   
  if(debugMode)hrec_sameside->Print("base");
  if(debugMode)hrec_sameside->Write(histTitle2.c_str());
  
  histTitle2+="_clone";
  TH1D *hrec_sameside_rebin = (TH1D*)hrec_sameside->Clone( (histTitle2).c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");
  if(debugMode)hrec_sameside_rebin->Write( histTitle2.c_str() );
  
  std::cout<<"rebinning hrec_sameside..."<<std::endl;
  histTitle2+="_rebins";
  hrec_sameside_rebin = (TH1D*)hrec_sameside_rebin->Rebin( nbins_pt_reco, (histTitle2).c_str() , boundaries_pt_reco);
  if(debugMode)hrec_sameside_rebin->Print("base");  
  if(debugMode)hrec_sameside_rebin->Write( histTitle2.c_str() );   
  
  if(clearOverUnderflows){
    histTitle2+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_sameside_rebin);
    if(debugMode)hrec_sameside_rebin->Print("base");  
    if(debugMode)hrec_sameside_rebin->Write( histTitle2.c_str() );  } 
  
  //cosmetics
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);
  hrec_sameside_rebin->SetLineColor(kBlue-3);
  hrec_sameside_rebin->SetMarkerSize(1.02);
  
  TH1D* hrec_sameside_resp_rebin=NULL;
  if(fillRespHists) 
    hrec_sameside_resp_rebin = (TH1D*)hrec_sameside_rebin->Clone("recrebinsamesideClone4unf");
  else
    hrec_sameside_resp_rebin = new TH1D( ("hpp_rec_sameside_response_rebin"+RandEtaRange).c_str(),"", 
					 nbins_pt_reco, boundaries_pt_reco); 
  
  if(debugMode)hrec_sameside_resp_rebin->Print(" base");  
  if(debugMode)hrec_sameside_resp_rebin->Write( histTitle2.c_str() );
  
  
  // ---------- gen, MC truth spectra
  std::string genHistTitle="hpp_mcclosure_gen";
  if(doJetID)genHistTitle+="_wJetID";
  genHistTitle+=RandEtaRange;
  
  TH1D* hgen = (TH1D*)fpp_MC->Get( genHistTitle.c_str() );
  if(debugMode)hgen->Print("base");    
  if(debugMode)hgen->Write();
  
  genHistTitle+="_clone";
  TH1D* hgen_rebin = (TH1D*)hgen->Clone( (genHistTitle).c_str() );
  if(debugMode)hgen_rebin->Print("base");
  if(debugMode)hgen_rebin->Write(genHistTitle.c_str());
  
  genHistTitle+="_rebins";
  hgen_rebin = (TH1D*)hgen_rebin->Rebin(nbins_pt_gen, (genHistTitle).c_str() , boundaries_pt_gen);
  if(debugMode)hgen_rebin->Print("base"); 
  if(debugMode)hgen_rebin->Write(genHistTitle.c_str());

  if(clearOverUnderflows){
    genHistTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hgen_rebin);
    if(debugMode)hgen_rebin->Print("base");    
    if(debugMode)hgen_rebin->Write(genHistTitle.c_str());}
  
  //cosmetics
  hgen_rebin->SetMarkerStyle(kOpenStar);
  hgen_rebin->SetMarkerColor(kMagenta);
  hgen_rebin->SetLineColor(kMagenta);
  hgen_rebin->SetMarkerSize(1.02);     	
  
  

  TH1D* hgen_resp_rebin;
  if(fillRespHists) 
    hgen_resp_rebin = (TH1D*)hgen_rebin->Clone("genrebinClone4unf");
  else
    hgen_resp_rebin = new TH1D( ("hpp_gen_response_rebin"+RandEtaRange).c_str() ,"", 
				nbins_pt_gen, boundaries_pt_gen);  
  if(debugMode)hgen_resp_rebin->Write();
  if(debugMode)hgen_resp_rebin->Print("base");  
  
  
  
  
  
  // ---------- open MC "response" matrix
  std::string TH2_title="hpp_mcclosure_matrix";//+RandEtaRange;
  if(doJetID)TH2_title+="_wJetID";
  TH2_title+=RandEtaRange;
  
  //get the response matrix made by readforests
  TH2D* hmat = (TH2D*)fpp_MC->Get( TH2_title.c_str() );
  if(debugMode)hmat->Write();
  if(debugMode)hmat->Print("base");
  
  // rebinned matrix ---------------
  TH2_title+="_clone";
  TH2D* hmat_rebin = (TH2D*)hmat->Clone( (TH2_title).c_str() );
  if(debugMode)hmat_rebin->Write( TH2_title.c_str() );
  if(debugMode)hmat_rebin->Print("base"); 
  
  TH2_title+="_rebins";
  hmat_rebin=(TH2D*) reBinTH2(hmat_rebin, (TH2_title).c_str(), 
			       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );  
  if(debugMode)hmat_rebin->Write( TH2_title.c_str() );
  if(debugMode)hmat_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    TH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2D*)hmat_rebin);
    if(debugMode)hmat_rebin->Write( TH2_title.c_str() );
    if(debugMode)hmat_rebin->Print("base");  }

    
  // error and %error for response matrix ---------------
  std::string errTH2_title="hmat_errors_rebins";
  TH2D* hmat_errors=makeRespMatrixErrors( (TH2D*) hmat,
					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  if(debugMode)hmat_errors->Write( errTH2_title.c_str() );
  if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    errTH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2D*)hmat_errors);
    if(debugMode)hmat_errors->Write(errTH2_title.c_str());
    if(debugMode)hmat_errors->Print("base");  }
  
  
  // give response matrix the correct errors
  setRespMatrixErrs( (TH2D*)hmat_rebin, (TH2D*) hmat_errors , (bool)zeroBins);
  TH2_title+="_wseterrs";
  if(debugMode)hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base");
  
  
  
  if(debugMode)std::cout<<"cloning input hists..."<<std::endl;
  TH1D* hrec_rebin_clone=(TH1D*)hrec_rebin->Clone("hrec_rebin_clone");
  setupSpectraHist(hrec_rebin_clone, useSimpBins);
  if(debugMode)hrec_rebin_clone->Print("base");

  TH1D* hgen_rebin_clone=(TH1D*)hgen_rebin->Clone("hgen_rebin_clone");
  setupSpectraHist(hgen_rebin_clone, useSimpBins);
  if(debugMode)hgen_rebin_clone->Print("base");
  
  TH1D* hrec_sameside_rebin_clone=(TH1D*)hrec_sameside_rebin->Clone("hrec_sameside_rebin_clone");
  setupSpectraHist(hrec_sameside_rebin_clone, useSimpBins);
  if(debugMode)hrec_sameside_rebin_clone->Print("base");  
  

  // SVD spectra
  TH1D *hunf_svd[nKregMax]={};        
  TH1D *hfold_svd[nKregMax]={};     
  TH1D *hunf_ss_svd[nKregMax]={};        
  TH1D *hfold_ss_svd[nKregMax]={};     
  TH1D *hfold_truth_svd[nKregMax]={};     
  
  // SVD ratios 
  TH1D* hgen_unfolded_ratio[nKregMax]={}; //data unf/gen truth
  TH1D* hgen_ssunfolded_ratio[nKregMax]={}; //data unf/gen truth
  TH1D* hgen_folded_ratio[nKregMax]={};   //not drawn anymore
  //  TH1D* hgen_ssfolded_ratio[nKregMax]={};   //not drawn anymore
  
  TH1D* hgen_meas_ratio=(TH1D*)hrec_rebin->Clone("hgen_meas_ratio_clone"); //data meas/gen truth
  hgen_meas_ratio->Divide(hgen_rebin);
  
  TH1D* hgen_ssmeas_ratio=(TH1D*)hrec_sameside_rebin->Clone("hgen_ssmeas_ratio_clone"); //mc meas/gen truth
  hgen_ssmeas_ratio->Divide(hgen_rebin);
  
  TH1D* hrec_unfolded_ratio[nKregMax]={}; //data unf/data meas
  TH1D* hrec_ssunfolded_ratio[nKregMax]={}; //data unf/data meas
  
  TH1D* hrec_fold_ratio[nKregMax]={};     // not really draw anymore
  TH1D* hrec_ssfold_ratio[nKregMax]={};     // not really draw anymore
  TH1D* hrec_truthfold_ratio[nKregMax]={};     // not really draw anymore
  
  TH1D* hrec_truth_ratio=(TH1D*)hgen_rebin->Clone("hrec_truth_ratio_clone");  //gen truth/data meas
  hrec_truth_ratio->Divide(hrec_rebin);

  TH1D* hrec_ssmeas_ratio=(TH1D*)hrec_sameside_rebin->Clone("hrec_ssmeas_ratio_clone"); //MC meas/data meas  
  hrec_ssmeas_ratio->Divide(hrec_rebin);
  
  TH1D* hgen_rebin_ratClone=(TH1D*)hgen_rebin->Clone("hgen_rebin_ratioClone");
  
  TH1D *hSVal=NULL;
  TH1D *hdi=NULL;
  

  // TH2s
  TH2D *hPearsonSVD[nKregMax];      
  TH2D *hCovmatSVD[nKregMax];
  TH2D *hCovmatAbsValSVD[nKregMax];
  //TH2D *hUnfmatSVD[nKregMax];
  
  TH2D *hRegCovMat   =NULL;
  TH2D *hRegCovMatInv=NULL;//i *think* it's the inv of the reg covmat, could be wrong and it's the data covmat instead...
  TH2D *hDataCovMat  =NULL;
  
  
  // TCanvases (do i really want this many? Do i really want this here in my code?)
  if(debugMode)std::cout<<"creating histos/arrays/canvases for SVD unfolding..."<<std::endl;
  TCanvas *cSpectra               = new TCanvas("cSpectra","",                      1800, 1500);      cSpectra->Divide(3,3);	   
  TCanvas *cRatio_gen             = new TCanvas("cRatio_gen","",                    1800, 1500);      cRatio_gen->Divide(3,3);      		
  TCanvas *cRatio_rec             = new TCanvas("cRatio_rec","",                    1800, 1500);      cRatio_rec->Divide(3,3);      		
  TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","",                1800, 1500);      cPearsonMatrixIter->Divide(3,3);    
  TCanvas *cCovMatrix         = new TCanvas("cCovMatrix","",                1800, 1500);      cCovMatrix->Divide(3,3);    
  TCanvas *cCovMatrixAbs         = new TCanvas("cCovMatrixAbs","",                1800, 1500);      cCovMatrixAbs->Divide(3,3);    
  TCanvas *c11                = new TCanvas("c11"," Singular Values and divectors", 1400, 1000);      c11->Divide(2);                    
  
  TLegend *leg[nKregMax],*leg1[nKregMax], *leg2[nKregMax];      
  
  
  
  TLine* theGenLine= new TLine( boundaries_pt_gen_mat[0]   	        , 1.,
				boundaries_pt_gen_mat[nbins_pt_gen_mat] , 1.  ); 
  theGenLine->SetLineWidth(1);
  theGenLine->SetLineStyle(2);
  theGenLine->SetLineColor(36);    
  

  TLine* theLineAtp9_gen= new TLine( boundaries_pt_gen_mat[0], 0.9, 
				     boundaries_pt_gen_mat[nbins_pt_gen_mat], 0.9);
  theLineAtp9_gen->SetLineWidth(1);
  theLineAtp9_gen->SetLineStyle(2);
  theLineAtp9_gen->SetLineColor(36);
  if(debugMode)std::cout<<"theLineAtp9_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
  if(debugMode)std::cout<<"theLineAtp9_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_gen_mat])  <<std::endl;
  
  TLine* theLineAt1p1_gen= new TLine( boundaries_pt_gen_mat[0], 1.1,
				      boundaries_pt_gen_mat[nbins_pt_gen_mat], 1.1);
  theLineAt1p1_gen->SetLineWidth(1);
  theLineAt1p1_gen->SetLineStyle(2);
  theLineAt1p1_gen->SetLineColor(36);

  TLine* theLineAtOne_gen= new TLine( boundaries_pt_gen_mat[0], 1.0,
				      boundaries_pt_gen_mat[nbins_pt_gen_mat], 1.0);
  theLineAtOne_gen->SetLineWidth(1);
  theLineAtOne_gen->SetLineStyle(2);
  theLineAtOne_gen->SetLineColor(36);
  
  
  TLine* theRecoLine= new TLine( boundaries_pt_reco_mat[0], 1.,
				 boundaries_pt_reco_mat[nbins_pt_reco_mat] , 1.  );
  theRecoLine->SetLineWidth(1);
  theRecoLine->SetLineStyle(2);
  theRecoLine->SetLineColor(36);      
  
  // prep for svd unfolding loop
  if(debugMode) std::cout<<std::endl<<"initializing kReg parameter array w/ "<<nKregMax<<" elements"<<std::endl;
  if(debugMode) std::cout<<kRegCenter-kRegRange<<" <= kReg[i] <= "<<kRegCenter+kRegRange<<std::endl;
  
  //for kreg = user spec
  int kReg[nKregMax]={0};
  for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) 
    kReg[i+kRegRange]=kRegCenter+i;  
  if(debugMode){
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i)
      std::cout<<"kReg["<<i+kRegRange<<"]="<<kReg[i+kRegRange]<<std::endl;}
  
  //for kregSS = 1 through 9
  int kRegSS[nKregMax]={0};
  for(int i=1; i<=9; ++i) 
    kRegSS[i-1]=i; 
  if(debugMode){
    for(int i=1; i<=9; ++i) 
      std::cout<<"kRegSS["<<i<<"]="<<kRegSS[i-1]<<std::endl;}
  
  std::cout<<std::endl;
  std::cout<<"creating instance of RooUnfoldResponse class"<<std::endl;    
  std::cout<<std::endl;
  std::string roo_resp_title = "Response_matrix_rebin_"+radius+std::to_string(kRegCenter);  
  RooUnfoldResponse roo_resp(hrec_sameside_resp_rebin, hgen_resp_rebin, hmat_rebin, (roo_resp_title).c_str());    
  roo_resp.UseOverflow(doOverUnderflows);
  if(debugMode)roo_resp.Write();
  
  TH1D* hfak=  (TH1D*) roo_resp.Hfakes() ;
  hfak->Print("base");
  hfak->SetMarkerColor(kGreen);
  hfak->SetLineColor(kGreen);
  hfak->SetMarkerSize(1.02);
  
  
  // thy spectra  
  TH1D* CT10nlo  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_CT10nlo  ).c_str());
  CT10nlo->SetMarkerSize(0);
  CT10nlo->SetLineColor(kBlack);  
  TH1D* CT14nlo  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_CT14nlo  ).c_str());
  CT14nlo->SetMarkerSize(0);
  CT14nlo->SetLineColor(kGreen);  
  TH1D* HERAPDF  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_HERAPDF  ).c_str());
  HERAPDF->SetMarkerSize(0);
  HERAPDF->SetLineColor(kViolet-5);  
  TH1D* MMHTnlo  =(TH1D*)makeThyHist_00eta20((fNLOFile_R04_MMHTnlo  ).c_str());
  MMHTnlo->SetMarkerSize(0);
  MMHTnlo->SetLineColor(kOrange+7);  
  TH1D* NNPDFnnlo=(TH1D*)makeThyHist_00eta20((fNLOFile_R04_NNPDFnnlo).c_str());
  NNPDFnnlo->SetMarkerSize(0);
  NNPDFnnlo->SetLineColor(kCyan-6);  
  
  // SVD unfolding loop
  for(int kr = 0; kr < nKregMax; ++kr){
    
    if(debugMode&&kr==0)    std::cout<<"beginning SVD Unfolding Loop..."<<std::endl;
    std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl;        
    
    std::string kRegRandEtaRange_plotTitle="kReg="+std::to_string(kReg[kr]);
    std::string kRegRandEtaRange="_kReg"+std::to_string(kReg[kr]);
    
    
    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd..."<<std::endl<<std::endl;
    RooUnfoldSvd unf_svd(&roo_resp, hrec_rebin, kReg[kr]);
    unf_svd.SetVerbose(0);    
    //if(kr==0) unf_svd.SetVerbose(2);
    //else  unf_svd.SetVerbose(0);
    //if(debugMode) unf_svd.SetVerbose(2);
    
    if(doToyErrs){
      std::cout<<"using toy errors, suppressing text output"<<std::endl;
      unf_svd.SetVerbose(1);    //trust me, so much output it slows the code down
      unf_svd.SetNToys(10000);     }
    if(debugMode) std::cout<<"RooUnfoldSvd Overflow Status: " << unf_svd.Overflow()<<std::endl;
    
    if(debugMode)std::cout<<"unfolding into hunf_svd[kr="<<kr<<"]..."<<std::endl;
    hunf_svd[kr] = (TH1D*)unf_svd.Hreco(errorTreatment);
    
    if(debugMode)std::cout<<"applying roo_resp to histo hunf_svd[kr="<<kr<<"]..."<<std::endl;
    hfold_svd[kr] = (TH1D*)roo_resp.ApplyToTruth(hunf_svd[kr]);
    
    
    
    
    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd..."<<std::endl<<std::endl;
    RooUnfoldSvd unf_ss_svd(&roo_resp, hrec_sameside_rebin, kReg[kr]);
    unf_ss_svd.SetVerbose(0);
    //if(kr==0) unf_ss_svd.SetVerbose(2);
    //else  unf_ss_svd.SetVerbose(0);
    //if(debugMode) unf_ss_svd.SetVerbose(2);
    
    if(doToyErrs){
      std::cout<<"using toy errors, suppressing text output"<<std::endl;
      unf_ss_svd.SetVerbose(1);    //trust me, so much output it slows the code down
      unf_ss_svd.SetNToys(10000);     }
    if(debugMode) std::cout<<"RooUnfoldSvd Overflow Status: " << unf_ss_svd.Overflow()<<std::endl;
    
    if(debugMode)std::cout<<"unfolding into hunf_svd[kr="<<kr<<"]..."<<std::endl;
    hunf_ss_svd[kr] = (TH1D*)unf_ss_svd.Hreco(errorTreatment);
        
    if(debugMode)std::cout<<"applying roo_resp to histo hunf_ss_svd[kr="<<kr<<"]..."<<std::endl;
    hfold_ss_svd[kr] = (TH1D*)roo_resp.ApplyToTruth(hunf_ss_svd[kr]);

    if(debugMode)std::cout<<"applying roo_resp to mc truth"<<std::endl;
    hfold_truth_svd[kr] = (TH1D*)roo_resp.ApplyToTruth(hgen_rebin);
    
    
    


    
    if(debugMode)std::cout<<"calling Ereco... getting Cov Mat"<<std::endl;
    TMatrixD covmat = unf_svd.Ereco(errorTreatment);
    hCovmatSVD[kr]=new TH2D(covmat);
    hCovmatSVD[kr]->SetName( ("hCovmatSVD"+kRegRandEtaRange).c_str() );      
    
    hCovmatAbsValSVD[kr]=(TH2D*)absVal_TH2Content((TH2D*)hCovmatSVD[kr]);
    
    if(debugMode)std::cout<<"calling CalculatePearsonCoefficients... getting pearson matrix"<<std::endl;
    hPearsonSVD[kr] = (TH2D*) CalculatePearsonCoefficients(&covmat, false, ("SVD_pearson_"+std::to_string(kr)) );      
    
    //hUnfmatSVD[kr]=new TH2D(unf_svd.UnfoldingMatrix());
    //hUnfmatSVD[kr]->SetName(("SVD_unfmatTH2_kr"+std::to_string(kr)).c_str());
    
    
    
    
    //  singular values and d_i vector ---------------------------
    //Note that these do not depend on the regularization.
    //The opposite: they tell you which regularization to use! (ian note: how?)
    if(kr == kRegDraw ){      
      
      std::cout << "  getting divectors and singular values for oppSide ... " << std::endl;
      TSVDUnfold *svdUnfold = unf_svd.Impl();
      //svdUnfold->SetNormalize(true);
      
      std::cout << "  getting singular values... " << std::endl;
      TH1D* temphSVal = (TH1D*)svdUnfold->GetSV();      
      
      std::cout << "  getting di vector values... " <<  std::endl;
      TH1D* temphdi = (TH1D*)svdUnfold->GetD();
      
      TH2D* temphRegCovMat    = (TH2D*)svdUnfold->GetXtau();
      TH2D* temphRegCovMatInv = (TH2D*)svdUnfold->GetXinv();
      TH2D* temphDataCovMat   = (TH2D*)svdUnfold->GetBCov();
      
      hSVal         = (TH1D*) temphSVal         ->Clone();      
      if(debugMode)
	for(int bin=1; bin<=hSVal->GetNbinsX(); bin++) std::cout<<"bin: "<<bin<<",  SV: "<<hSVal->GetBinContent(bin)<< std::endl;
      
      hdi           = (TH1D*) temphdi           ->Clone();
      if(debugMode)
	for(int bin=1; bin<=hdi->GetNbinsX(); bin++)	  std::cout<<"i: "<<bin<<",  di: "<<hdi->GetBinContent(bin)<<std::endl;
      
      hRegCovMat    = (TH2D*) temphRegCovMat    ->Clone();
      hRegCovMatInv = (TH2D*) temphRegCovMatInv ->Clone();
      hDataCovMat   = (TH2D*) temphDataCovMat   ->Clone();


      
      if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
      c11->cd(1);
      c11->cd(1)->SetLogy();  
      
      
      //divBinWidth_DiAndSVals( (double*)boundaries_pt_gen, (int)nbins_pt_gen, (TH1D*)hSVal);
      
      //hSVal->SetTitle(" singular values ");
      hSVal->SetAxisRange(1.,(double)(hSVal->GetNbinsX()-1),"X");
      //hSVal->SetAxisRange(0.,(double)(hSVal->GetNbinsX()),"X");
      hSVal->SetXTitle("i");        
      hSVal->SetYTitle("Singular Values s_{i}");        
      //hSVal->DrawCopy("HIST E");
      hSVal->Draw("HIST E");
      
      
      // di vector values
      c11->cd(2);
      c11->cd(2)->SetLogy();    

      //dont do this. misguided at best.
      //divBinWidth_DiAndSVals( (double*)boundaries_pt_gen, (int)nbins_pt_gen, (TH1D*)hdi);
      
      hdi->SetAxisRange(1.,(double)(hdi->GetNbinsX()-1),"X");
      //hdi->SetAxisRange(0.,(double)(hdi->GetNbinsX()),"X");
      hdi->SetXTitle("i");
      hdi->SetYTitle("Divector Values d_{i}");
      //hdi->DrawCopy("HIST E"); 
      hdi->Draw("HIST E"); 
      
      TLine* theLineAtOne_hdi=new TLine(1., 1., (double)(hdi->GetNbinsX()), 1.);
      theLineAtOne_hdi->SetLineWidth(1);
      theLineAtOne_hdi->SetLineStyle(2);
      theLineAtOne_hdi->SetLineColor(36);
      theLineAtOne_hdi->Draw();
      
      //drawText( "5.02 TeV ak4PFJets",0.358173, 0.8459761, 19);
      //drawText( "2015 Prompt Reco",0.358173, 0.8159761, 19);
      drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
      drawText( ("Current kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);	
      if(debugMode)std::cout<<std::endl<<"done with kr=="<< kRegDraw<<" specifics"<<std::endl<<std::endl;          }

    
  }// kReg loop ends
  
  //scale rec/gen clones to bin widths
  divideBinWidth((TH1D*)hgen_rebin);
  hgen_rebin->Scale(1./etaBinWidth);

  divideBinWidth((TH1D*)hgen_resp_rebin);
  hgen_resp_rebin->Scale(1./etaBinWidth);

  divideBinWidth((TH1D*)hgen_rebin_clone);
  hgen_rebin_clone->Scale(1./etaBinWidth);

  divideBinWidth((TH1D*)hgen_rebin_ratClone);
  hgen_rebin_ratClone->Scale(1./etaBinWidth);    
  
  divideBinWidth((TH1D*)hrec_rebin);
  hrec_rebin->Scale(1./etaBinWidth);
  
  divideBinWidth((TH1D*)hrec_rebin_clone);
  hrec_rebin_clone->Scale(1./etaBinWidth);
  
  divideBinWidth((TH1D*)hrec_sameside_rebin);
  hrec_sameside_rebin->Scale(1./etaBinWidth);
  
  divideBinWidth((TH1D*)hrec_sameside_rebin_clone);
  hrec_sameside_rebin_clone->Scale(1./etaBinWidth);
  
  divideBinWidth((TH1D*)hrec_sameside_resp_rebin);
  hrec_sameside_resp_rebin->Scale(1./etaBinWidth);
  
  divideBinWidth((TH1D*) hfak);
  hfak->Scale(1./etaBinWidth);
  
  //correct hrec hist here for fakes
  Double_t fac_min1=1.;//defined + set this to 1 because i'm too lazy to fix 
  Double_t fac=(hrec_rebin->Integral()/hrec_sameside_rebin->Integral());
  
  TH1D* hrec_rebin_fakecorr         = (TH1D*)hrec_rebin->Clone(          ( ( (std::string)hrec_rebin->GetName()         )+"_fakecorr").c_str() );
  TH1D* hrec_sameside_rebin_fakecorr= (TH1D*)hrec_sameside_rebin->Clone( ( ( (std::string)hrec_sameside_rebin->GetName())+"_fakecorr").c_str() );
  
  //fake correction
  int nbins_fakcorr=hrec_rebin_fakecorr->GetNbinsX();
  for(int i=1;i<=nbins_fakcorr; i++){
    Double_t hrec_i     = hrec_rebin_fakecorr->GetBinContent(i); Double_t hrec_i_err = hrec_rebin_fakecorr->GetBinError(i);
    Double_t hrecss_i     = hrec_sameside_rebin_fakecorr->GetBinContent(i); Double_t hrecss_i_err = hrec_sameside_rebin_fakecorr->GetBinError(i);
    Double_t hfak_i     = hfak->GetBinContent(i); Double_t hfak_i_err = hfak->GetBinError(i);
    hrec_rebin_fakecorr->SetBinContent( i, hrec_i-fac*hfak_i );     
    hrec_rebin_fakecorr->SetBinError(   i, sqrt(hrec_i_err*hrec_i_err + fac*hfak_i_err*fac*hfak_i_err));  
    hrec_sameside_rebin_fakecorr->SetBinContent( i, hrecss_i-fac_min1*hfak_i );     
    hrec_sameside_rebin_fakecorr->SetBinError(   i, sqrt(hrecss_i_err*hrecss_i_err + fac_min1*hfak_i_err*fac_min1*hfak_i_err));  
  }

  

  
  // -- MC OS UNF x2 (for thy comparison) -- //
  TH1D* hunf_x2=NULL;

  for(int kr = 0; kr < nKregMax; ++kr){
    
    //scale unf to bin widths
    divideBinWidth((TH1D*)hunf_svd[kr]);
    hunf_svd[kr]->Scale(1./etaBinWidth);    

    divideBinWidth((TH1D*)hunf_ss_svd[kr]);
    hunf_ss_svd[kr]->Scale(1./etaBinWidth);    
    
    //scale fold to bin widths
    divideBinWidth((TH1D*)hfold_svd[kr]);
    hfold_svd[kr]->Scale(1./etaBinWidth);

    divideBinWidth((TH1D*)hfold_ss_svd[kr]);
    hfold_ss_svd[kr]->Scale(1./etaBinWidth);

    divideBinWidth((TH1D*)hfold_truth_svd[kr]);
    hfold_truth_svd[kr]->Scale(1./etaBinWidth);
    
    std::string kRegRandEtaRange_plotTitle="kReg="+std::to_string(kReg[kr]);
    std::string kRegRandEtaRange="_kReg"+std::to_string(kReg[kr]);    
    

    //hunf cosmetics + other
    hunf_svd[kr]->SetMarkerStyle(kOpenCircle);
    hunf_svd[kr]->SetMarkerColor(kRed);     
    hunf_svd[kr]->SetLineColor(kRed);     
    hunf_svd[kr]->SetMarkerSize(1.01);     
    
    hunf_svd[kr]->SetName( ("hunf_svd"+kRegRandEtaRange).c_str() );
    hunf_svd[kr]->SetTitle( ("OS MC SVD Unf. "+kRegRandEtaRange_plotTitle).c_str() );
    if(debugMode)hunf_svd[kr]->Print("base");                
    
    

    //hunf cosmetics + other
    hunf_ss_svd[kr]->SetMarkerStyle(kOpenSquare);
    hunf_ss_svd[kr]->SetMarkerColor(kRed+1);     
    hunf_ss_svd[kr]->SetLineColor(kRed+1);     
    hunf_ss_svd[kr]->SetMarkerSize(1.01);     
    
    hunf_ss_svd[kr]->SetName( ("hunf_ss_svd"+kRegRandEtaRange).c_str() );
    hunf_ss_svd[kr]->SetTitle( ("SS MC SVD Unf. "+kRegRandEtaRange_plotTitle).c_str() );
    if(debugMode)hunf_svd[kr]->Print("base");                
    

    
    
    //hfold cosmetics
    hfold_svd[kr]->SetMarkerStyle(kOpenCircle);
    hfold_svd[kr]->SetMarkerColor(kGreen-5);      
    hfold_svd[kr]->SetLineColor(kGreen-5);      
    hfold_svd[kr]->SetMarkerSize(1.01);      
    
    hfold_svd[kr]->SetName( ("hfold_svd"+kRegRandEtaRange).c_str() );
    if(debugMode)hfold_svd[kr]->Print("base");    

    //hfold ss cosmetics
    hfold_ss_svd[kr]->SetMarkerStyle(kOpenSquare);
    hfold_ss_svd[kr]->SetMarkerColor(kGreen-5);      
    hfold_ss_svd[kr]->SetLineColor(kGreen-5);      
    hfold_ss_svd[kr]->SetMarkerSize(1.01);      
    
    hfold_ss_svd[kr]->SetName( ("hfold_ss_svd"+kRegRandEtaRange).c_str() );
    if(debugMode)hfold_ss_svd[kr]->Print("base");    
    
    //hfold truth cosmetics
    hfold_truth_svd[kr]->SetMarkerStyle(kOpenStar);
    hfold_truth_svd[kr]->SetMarkerColor(kGreen+3);      
    hfold_truth_svd[kr]->SetLineColor(kGreen+3);      
    hfold_truth_svd[kr]->SetMarkerSize(1.01);      
    
    hfold_truth_svd[kr]->SetName( ("hfold_truth_svd"+kRegRandEtaRange).c_str() );
    if(debugMode)hfold_truth_svd[kr]->Print("base");    
    

    
    
    
    // ----------------------- 3X3 CANVAS DRAWING  ----------------------- //
    // for data unfolding
    // CAN I DO THIS ELSEWHERE/DIFFERENTLY? //
    
    // 3x3 spectra canvases  ---------------------
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
    cSpectra->cd(kr+1);
    cSpectra->cd(kr+1)->SetLogy(1);
    if(!useSimpBins)
      cSpectra->cd(kr+1)->SetLogx(1);
    
    //setupSpectraHist(hrec_rebin_clone,useSimpBins);
    setupSpectraHist(hunf_svd[kr],useSimpBins);
    //hrec_rebin_clone->SetTitle( ("Jet Spectra, "+kRegRandEtaRange_plotTitle).c_str() ); //FOR TITLE OF PLOT
    //hrec_rebin_clone->GetYaxis()->SetTitle("N_{Jets}/L_{int}");
    hunf_svd[kr]->SetTitle( ("Jet Spectra,"+kRegRandEtaRange_plotTitle).c_str() );

    hunf_svd[kr]->Draw("P E");
    hrec_rebin_clone->Draw("P E SAME");
    hgen_rebin_clone->Draw("P E SAME");
    //hfold_svd[kr]->Draw("P E SAME");      
    
    leg[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option 
    leg[kr]->AddEntry(hrec_rebin_clone,"OS MC Meas.","p");
    leg[kr]->AddEntry(hgen_rebin_clone,"SS MC Truth","p");
    leg[kr]->AddEntry(hunf_svd[kr],"OS MC Unf.","p");
    //leg[kr]->AddEntry(hfold_svd[kr],"OS MC Fold(Unf.)","p");
    leg[kr]->SetTextSize(0.03);
    leg[kr]->Draw();
    
    
    // 3x3 unfORfold/meas ratio canvases  ---------------------
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_rec canvas..."<<std::endl;
    cRatio_rec->cd(kr+1);
    if(!useSimpBins)
      cRatio_rec->cd(kr+1)->SetLogx(1);
    
    hrec_unfolded_ratio[kr] = (TH1D*)hunf_svd[kr]->Clone( ("hrec_unf_ratio"+kRegRandEtaRange).c_str());
    hrec_unfolded_ratio[kr]->Divide(hrec_rebin);
    if(debugMode)hrec_unfolded_ratio[kr]->Print("base");

    hrec_ssunfolded_ratio[kr] = (TH1D*)hunf_ss_svd[kr]->Clone( ("hrec_ssunf_ratio"+kRegRandEtaRange).c_str());
    hrec_ssunfolded_ratio[kr]->Divide(hrec_rebin);
    if(debugMode)hrec_unfolded_ratio[kr]->Print("base");
    
    

    

    hrec_fold_ratio[kr] = (TH1D*)hfold_svd[kr]->Clone( ("hrec_fold_ratio"+kRegRandEtaRange).c_str() );
    hrec_fold_ratio[kr]->Divide(hrec_rebin_fakecorr);
    if(debugMode)hrec_fold_ratio[kr]->Print("base");
    
    hrec_ssfold_ratio[kr] = (TH1D*)hfold_ss_svd[kr]->Clone( ("hrec_ssfold_ratio"+kRegRandEtaRange).c_str() );
    hrec_ssfold_ratio[kr]->Divide(hrec_sameside_rebin_fakecorr);
    if(debugMode)hrec_ssfold_ratio[kr]->Print("base");
    
    hrec_truthfold_ratio[kr] = (TH1D*)hfold_truth_svd[kr]->Clone( ("hrec_truthfold_ratio"+kRegRandEtaRange).c_str() );
    hrec_truthfold_ratio[kr]->Divide(hrec_sameside_rebin_fakecorr);
    if(debugMode)hrec_truthfold_ratio[kr]->Print("base");
    
    
    hrec_unfolded_ratio[kr]->SetTitle( ("OS MC Unf./OS MC Meas., "+kRegRandEtaRange_plotTitle).c_str() );
    setupRatioHist(hrec_unfolded_ratio[kr],useSimpBins);
    //hrec_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
    hrec_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
    
    hrec_unfolded_ratio[kr]->Draw("P E");
    //hrec_fold_ratio[kr]->Draw("P E SAME");
    
    leg1[kr] = new TLegend(0.67, 0.80, 0.9, 0.9, NULL,"BRNDC");   
    leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],"OS MC Unf./OS MC Meas.","pl");
    //leg1[kr]->AddEntry(hrec_fold_ratio[kr],"OS MC Fold(Unf.)/OS MC Meas.","pl");
    leg1[kr]->SetTextSize(0.02); 
    leg1[kr]->Draw();
    
    theRecoLine->Draw();
    
    
    
    // draw on cRatio_gen ratio canvas -----------------------------------------
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_gen canvas..."<<std::endl;
    cRatio_gen->cd(kr+1);
    if(!useSimpBins)
      cRatio_gen->cd(kr+1)->SetLogx(1);
    
    hgen_unfolded_ratio[kr] = (TH1D*)hunf_svd[kr]->Clone( ("hgen_unf_ratio"+kRegRandEtaRange).c_str());
    hgen_unfolded_ratio[kr]->Divide(hgen_rebin);      
    if(debugMode)hgen_unfolded_ratio[kr]->Print("base");

    hgen_ssunfolded_ratio[kr] = (TH1D*)hunf_ss_svd[kr]->Clone( ("hgen_ssunf_ratio"+kRegRandEtaRange).c_str());
    hgen_ssunfolded_ratio[kr]->Divide(hgen_rebin);      
    if(debugMode)hgen_ssunfolded_ratio[kr]->Print("base");
    
    hgen_folded_ratio[kr] = (TH1D*)hfold_svd[kr]->Clone( ("hgen_fold_ratio"+kRegRandEtaRange).c_str());            
    hgen_folded_ratio[kr]->Divide(hgen_rebin_ratClone);      
    if(debugMode)hgen_folded_ratio[kr]->Print("base");      
    
    
    hgen_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");    
    hgen_unfolded_ratio[kr]->SetTitle( ("OS MC Unf./SS MC Truth, "+kRegRandEtaRange_plotTitle).c_str() );
    hgen_folded_ratio[kr]->SetTitle( ("SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
    
    setupRatioHist(hgen_unfolded_ratio[kr],useSimpBins);
    
    hgen_unfolded_ratio[kr]->Draw("P E");
    //hgen_folded_ratio[kr]->Draw("P E SAME");
    
    leg2[kr] = new TLegend(0.67, 0.80, 0.9, 0.9, NULL,"BRNDC");
    leg2[kr]->AddEntry(hgen_unfolded_ratio[kr],"OS MC Unf./SS MC Truth","pl");
    //leg2[kr]->AddEntry(hgen_folded_ratio[kr],"OS MC Fold(Unf.)/SS MC Trtuh","pl");
    leg2[kr]->SetTextSize(0.02); 

    leg2[kr]->Draw();
    
    theGenLine->Draw();
    
    
    
    
    // Covariance Matrix Drawing
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cCovMatrix canvas..."<<std::endl;
    cCovMatrix->cd(kr+1);
    cCovMatrix->cd(kr+1)->SetLogz(1);
    
    hCovmatSVD[kr]->SetTitle( ("SVD Covariance Matrix, "+kRegRandEtaRange_plotTitle).c_str()  );
    hCovmatSVD[kr]->GetXaxis()->SetTitle("RECO Jet p_{T} Bin #");
    hCovmatSVD[kr]->GetYaxis()->SetTitle("GEN Jet p_{T} Bin #");
    
    hCovmatSVD[kr]->SetAxisRange(1e-15,1e+02,"Z");
    hCovmatSVD[kr]->GetZaxis()->SetLabelSize(0.035);
    hCovmatSVD[kr]->Draw("COLZ");



    // AbsVal Covariance Matrix Drawing
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cCovMatrix canvas..."<<std::endl;
    cCovMatrixAbs->cd(kr+1);
    cCovMatrixAbs->cd(kr+1)->SetLogz(1);
    
    hCovmatAbsValSVD[kr]->SetTitle( ("SVD |Covariance Matrix|, "+kRegRandEtaRange_plotTitle).c_str()  );
    hCovmatAbsValSVD[kr]->GetXaxis()->SetTitle("RECO Jet p_{T} Bin #");
    hCovmatAbsValSVD[kr]->GetYaxis()->SetTitle("GEN Jet p_{T} Bin #");
    
    hCovmatAbsValSVD[kr]->SetAxisRange(1e-15,1e+02,"Z");
    hCovmatAbsValSVD[kr]->GetZaxis()->SetLabelSize(0.035);
    hCovmatAbsValSVD[kr]->Draw("COLZ");
    
    if(kReg[kr]==kRegDraw){
      hunf_x2=(TH1D*)hunf_svd[kr]->Clone( ( ((std::string)hunf_svd[kr]->GetName()) + "_2xScaled").c_str() );
      hunf_x2->Scale(2.);
      if(debugMode)hunf_x2->Write();
      if(debugMode)hunf_x2->Print("base");
    }
    
    // Pearson Drawing
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter canvas..."<<std::endl;
    cPearsonMatrixIter->cd(kr+1);
    
    hPearsonSVD[kr]->SetName( ("hPearsonSVD"+kRegRandEtaRange).c_str() );      
    hPearsonSVD[kr]->SetTitle( ("SVD Pearson Matrix, "+kRegRandEtaRange_plotTitle).c_str()  );
    hPearsonSVD[kr]->GetXaxis()->SetTitle("RECO Jet p_{T} Bin #");
    hPearsonSVD[kr]->GetYaxis()->SetTitle("GEN Jet p_{T} Bin #");
    
    hPearsonSVD[kr]->GetZaxis()->SetLabelSize(0.035);
    hPearsonSVD[kr]->SetAxisRange(-1., 1., "Z");    
    hPearsonSVD[kr]->Draw("COLZ");



    // ----------------------- END 3X3 CANVAS DRAWING for data unfolding ----------------------- //
  }
  
  

  

  
  // -- MC OS UNF x2 (for thy comparison) -- //
  TH1D* hgen_rebin_x2=(TH1D*)hgen_rebin->Clone( ( ((std::string)hgen_rebin->GetName()) + "_2xScaled").c_str() );
  hgen_rebin_x2->Scale(2.);
  if(debugMode)hgen_rebin_x2->Write();
  if(debugMode)hgen_rebin_x2->Print("base");
  
  // ---------------- THY RATIOS (THY / DATA) ----------------- //
  // must be after the binwidth divisions + normalization (the thy hists are made this way by default)
  
  TH1D* h_thyratio_CT10nlo  =(TH1D*)CT10nlo  ->Clone("");
  h_thyratio_CT10nlo = (TH1D*)h_thyratio_CT10nlo->Rebin(nbins_pt_gen,"pp_CT10Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT10nlo  ->SetTitle("CT10 NLO/OS MC Unf.");  
  h_thyratio_CT10nlo  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_CT10nlo  ->Print("base");
  
  TH1D* h_thyratio_CT14nlo  =(TH1D*)CT14nlo  ->Clone("");  
  h_thyratio_CT14nlo=(TH1D*)h_thyratio_CT14nlo->Rebin(nbins_pt_gen,"pp_CT14Thy_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_CT14nlo  ->SetTitle("CT14 NLO/OS MC Unf.");
  h_thyratio_CT14nlo  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_CT14nlo  ->Print("base");
  
  TH1D* h_thyratio_HERAPDF  =(TH1D*)HERAPDF  ->Clone("");  
  h_thyratio_HERAPDF=(TH1D*)h_thyratio_HERAPDF->Rebin(nbins_pt_gen,"pp_HERAPDF_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_HERAPDF  ->SetTitle("HERAPDF15 NLO/OS MC Unf.");
  h_thyratio_HERAPDF  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_HERAPDF  ->Print("base");
  
  TH1D* h_thyratio_MMHTnlo  =(TH1D*)MMHTnlo  ->Clone("");
  h_thyratio_MMHTnlo=(TH1D*)h_thyratio_MMHTnlo->Rebin(nbins_pt_gen,"pp_MMHT_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_MMHTnlo  ->SetTitle("MMHT14 NLO/OS MC Unf.");
  h_thyratio_MMHTnlo  ->Divide(hunf_x2);
  if(debugMode)h_thyratio_MMHTnlo  ->Print("base");
  
  TH1D* h_thyratio_NNPDFnnlo=(TH1D*)NNPDFnnlo->Clone("");
  h_thyratio_NNPDFnnlo=(TH1D*)h_thyratio_NNPDFnnlo->Rebin(nbins_pt_gen,"pp_NNPDFnlo_Ratio_rebin",boundaries_pt_gen);
  h_thyratio_NNPDFnnlo->SetTitle("NNPDF NLO/OS MC Unf.");
  h_thyratio_NNPDFnnlo->Divide(hunf_x2);
  if(debugMode)h_thyratio_NNPDFnnlo->Print("base");
  
  
  if(drawPDFs){
    
    std::cout<<std::endl<<"drawing PDFs for SVD Unfolding..."<<std::endl;
    
    
    // form filename string, open pdf file and draw the canvases we have so far
    std::string out3x3PdfFile=out3x3SVDPdfFile;
    std::string open_out3x3PdfFile=out3x3PdfFile+"[";      std::string close_out3x3PdfFile=out3x3PdfFile+"]";

    //open 3x3 pdf file
    cSpectra->Print(open_out3x3PdfFile.c_str());
    

    //print already-drawn 3x3 canvases first
    //  cSpectra 3X3 --------------------------------------------------
    cSpectra->cd()           ;   
    cSpectra->Print(out3x3PdfFile.c_str()); 
    //  cRatio 3X3 Check for gen --------------------------------------------------
    cRatio_gen->cd()         ;   
    cRatio_gen->Print(out3x3PdfFile.c_str());    
    //  cRatio 3X3 Check for meas --------------------------------------------------
    cRatio_rec->cd()         ;   
    cRatio_rec->Print(out3x3PdfFile.c_str());    
    //  already drawn data-side 3x3 pearson matrix --------------------------------------------------            
    cPearsonMatrixIter->cd() ; 
    cPearsonMatrixIter->Print(out3x3PdfFile.c_str());     
    //  already drawn data-side 3x3 pearson matrix --------------------------------------------------            
    cCovMatrix->cd() ; 
    cCovMatrix->Print(out3x3PdfFile.c_str());     
    //  already drawn data-side 3x3 pearson matrix --------------------------------------------------            
    cCovMatrixAbs->cd() ; 
    cCovMatrixAbs->Print(out3x3PdfFile.c_str());     
    
    
    //close 3x3 pdf file
    cCovMatrixAbs->Print(close_out3x3PdfFile.c_str());    
    
    
    
    
    // form filename string, open pdf file and draw the canvases we have so far
    std::string outPdfFile=outSVDPdfFile;
    std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
    
    //   --------------------------------------------------
    // now the specific kr=kRegDraw plots
    TCanvas *cCheck        = new TCanvas("cCheck","",        1400, 1000);          
    cCheck->Print(open_outPdfFile.c_str());


    //  data-side singular values+divectors --------------------------------------------------
    c11->cd()                ;   
    c11->Print(outPdfFile.c_str());
    



    //   --------------------------------------------------
    // spectra for best kreg choice
    cCheck->cd();
    cCheck->SetLogy(1);
    if(!useSimpBins)	
      cCheck->SetLogx(1);
    
    hrec_rebin_clone->SetTitle("SVD, OS and SS MC Spectra");
    setupSpectraHist(hrec_rebin_clone, useSimpBins, boundaries_pt_reco, nbins_pt_reco);
    
    hrec_rebin_clone->Draw("P E");
    hrec_sameside_rebin_clone->Draw("P E SAME");
    hgen_rebin_clone->Draw("P E SAME");	
    hunf_svd[kRegDraw]->Draw("P E SAME");
    hunf_ss_svd[kRegDraw]->Draw("P E SAME");
    //hfold_svd[kRegDraw]->Draw("P E SAME");
    
    
    TLegend * leg1 = new TLegend(0.7, 0.7, 0.9, 0.9, NULL,"NBNDC");
    leg1->AddEntry(hrec_rebin_clone,   "OS MC Meas","lp");
    leg1->AddEntry(hrec_sameside_rebin_clone,   "SS MC Meas","lp");
    leg1->AddEntry(hunf_svd[kRegDraw], "OS MC Unf.","lp");
    leg1->AddEntry(hunf_ss_svd[kRegDraw], "SS MC Unf.","lp");
    leg1->AddEntry(hgen_rebin_clone,   "SS MC Truth", "lp");	
    //leg1->AddEntry(hfold_svd[kRegDraw],"Data Fold","p");
    //leg1->SetTextSize(0.02);
    leg1->Draw();
    
    //drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
    //drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
    //drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    
    cCheck->Print(outPdfFile.c_str());
    
    
    //  ratios w/ MC truth --------------------------------------------------
    cCheck->cd();
    cCheck->SetLogy(0);
    if(!useSimpBins)
      cCheck->SetLogx(1);
    
    hgen_unfolded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ SS MC Truth Spectra");            
    hgen_unfolded_ratio[kRegDraw]->GetYaxis()->SetTitle("Ratio w/ SS MC Truth");    
    setupRatioHist(hgen_unfolded_ratio[kRegDraw], useSimpBins);
    
    
    hgen_unfolded_ratio[kRegDraw]->Draw("P E");
    hgen_ssunfolded_ratio[kRegDraw]->Draw("P E SAME");
    hgen_meas_ratio->Draw("P E SAME");
    hgen_ssmeas_ratio->Draw("P E SAME");
    //hgen_folded_ratio[kRegDraw]->Draw("P E SAME");
    
    TLegend * leg3 = new TLegend(0.1, 0.8, 0.2, 0.9, NULL,"NBNDC");
    leg3->AddEntry(hgen_meas_ratio,              "OS MC Meas.","pl");
    leg3->AddEntry(hgen_unfolded_ratio[kRegDraw],"OS MC Unf.","pl");
    leg3->AddEntry(hgen_ssmeas_ratio,            "SS MC Meas.","pl");
    leg3->AddEntry(hgen_ssunfolded_ratio[kRegDraw],"SS MC Unf.","pl");

    //leg3->AddEntry(hgen_folded_ratio[kRegDraw],"MC Fold(Unf.)/SS MC Truth","pl");
    //leg3->SetTextSize(0.02);
    leg3->Draw();
    
    //drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
    //drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
    //drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    
    theGenLine->Draw();
    theLineAtp9_gen->Draw();
    theLineAt1p1_gen->Draw();
    
    
    cCheck->Print(outPdfFile.c_str());
    
    
    


    
    // Ratios w/ Meas. --------------------------------------------------
    cCheck->cd();
    cCheck->SetLogy(0);
    if(!useSimpBins)
      cCheck->SetLogx(1);
    
    hrec_unfolded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ OS MC Meas. Spectra");    
    hrec_unfolded_ratio[kRegDraw]->GetYaxis()->SetTitle("Ratio w/ OS MC Meas.");
    setupRatioHist(hrec_unfolded_ratio[kRegDraw], useSimpBins);
    
    
    hrec_unfolded_ratio[kRegDraw]->Draw("P E");
    hrec_ssunfolded_ratio[kRegDraw]->Draw("P E SAME");
    hrec_truth_ratio->Draw("P E SAME");
    hrec_ssmeas_ratio->Draw("P E SAME");
    //hrec_fold_ratio[kRegDraw]->Draw("P E SAME");
    
    TLegend * leg2 = new TLegend(0.1, 0.8, 0.2, 0.9, NULL,"NBNDC");
    leg2->AddEntry(hrec_ssmeas_ratio,"SS MC Meas.","pl");
    leg2->AddEntry(hrec_ssunfolded_ratio[kRegDraw],"SS MC Unf.","pl");
    leg2->AddEntry(hrec_truth_ratio,"SS MC Truth","pl");
    leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"OS MC Unf.","pl");
    //leg2->AddEntry(hrec_fold_ratio[kRegDraw],"OS MC Fold(Unf.)/OS MC Meas.","pl");
    //leg2->SetTextSize(0.02);
    leg2->Draw();
    
    //drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
    //drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
    //drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    
    theRecoLine->Draw();
    theLineAtp9_gen->Draw();
    theLineAt1p1_gen->Draw();
    
    cCheck->Print(outPdfFile.c_str());



    // Folding Ratios --------------------------------------------------
    cCheck->cd();
    cCheck->SetLogy(0);
    if(!useSimpBins)
      cCheck->SetLogx(1);
    
    hrec_fold_ratio[kRegDraw]->SetTitle("Folded Ratios w/ (Corrected) MC Meas. Spectra");
    hrec_fold_ratio[kRegDraw]->GetYaxis()->SetTitle("Ratio w/ MC Meas. Spectra");
    setupRatioHist(hrec_fold_ratio[kRegDraw], useSimpBins);
    
    
    hrec_fold_ratio[kRegDraw]->Draw("P E");
    hrec_ssfold_ratio[kRegDraw]->Draw("P E SAME");
    hrec_truthfold_ratio[kRegDraw]->Draw("P E SAME");
    
    TLegend * legfold = new TLegend(0.1, 0.8, 0.4, 0.9, NULL,"NBNDC");
    legfold->AddEntry(hrec_fold_ratio[kRegDraw],"OS MC Fold(Unf.)/OS MC Meas.","pl");
    legfold->AddEntry(hrec_ssfold_ratio[kRegDraw],"SS MC Fold(Unf.)/SS MC Meas.","pl");
    legfold->AddEntry(hrec_truthfold_ratio[kRegDraw],"SS MC Fold(Truth)/SS MC Meas.","pl");
    legfold->Draw();
    
    //drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
    //drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
    //drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    
    theRecoLine->Draw();
    theLineAtp9_gen->Draw();
    theLineAt1p1_gen->Draw();
    
    cCheck->Print(outPdfFile.c_str());


    //   --------------------------------------------------
    //  MC response matrix fakes
    cCheck->cd();
    cCheck->SetLogy(1);
    if(!useSimpBins)
      cCheck->SetLogx(1);
    
    hfak->SetTitle("MC Only, Spectra and Kinematic Fakes");
    setupSpectraHist(hfak, useSimpBins, boundaries_pt_reco, nbins_pt_reco);
    
    hfak->Draw("P E");
    hrec_sameside_rebin->Draw("P E SAME");
    hrec_rebin->Draw("P E SAME");
    hgen_resp_rebin->Draw("P E SAME");
    
    TLegend * legfake = new TLegend(0.6, 0.7, 0.9, 0.9, NULL,"NBNDC");
    legfake->AddEntry(hfak,"SS MC Kinematic Fakes","lp");
    legfake->AddEntry(hgen_resp_rebin,"SS MC Truth" ,"lp");
    legfake->AddEntry(hrec_sameside_rebin,"SS MC Meas." ,"lp");
    legfake->AddEntry(hrec_rebin,"OS MC Meas." ,"lp");
    //legfake->SetTextSize(0.02);
    legfake->Draw();
    
//    drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
//    drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
//    drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    
    cCheck->Print(outPdfFile.c_str());
    
    
    // thy spectra CT10/14 NNPDF NLO---------------------------

    TLegend* legendThy1 =new TLegend( 0.7,0.7,0.9,0.9 );    
    
    //CT10nlo->SetLineColor(kRed);
    setupSpectraHist(CT10nlo  ,useSimpBins);
    legendThy1->AddEntry(CT10nlo  ,"CT10 NLO","l");
    if(debugMode){fout->cd(); CT10nlo->Write("CT10_NLO_R04_jtpt");}
    
    //CT14nlo->SetLineColor(kGreen);
    setupSpectraHist(CT14nlo  ,useSimpBins);
    legendThy1->AddEntry(CT14nlo  ,"CT14 NLO","l");
    if(debugMode){fout->cd(); CT14nlo->Write("CT14_NLO_R04_jtpt");}
    
    if(debugMode)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    cCheck->cd();
    if(!useSimpBins)cCheck->SetLogx(1);
    cCheck->SetLogy(1);                  
    
    legendThy1->AddEntry(hunf_x2,"2 x (OS MC Unf.)","lp");
    legendThy1->AddEntry(hgen_rebin_x2,"2 x (SS MC Truth)", "lp");
        
    CT10nlo->SetAxisRange(1e-07,1e+03,"Y");//for y axis in nanbarns
    CT10nlo->SetTitle("NLO Thy w/ MC, OS Unf., SS Truth");

    CT10nlo  ->Draw("][HIST E");
    CT14nlo  ->Draw("][HIST E SAME");
    
    hgen_rebin_x2->Draw("P E SAME");
    hunf_x2->Draw("P E SAME");   //just for axis range
    

    legendThy1->Draw();

    cCheck->Print(outPdfFile.c_str());


    //   --------------------------------------------------
    // thy spectra HERA/MMHT/NNPDF ---------------------------

    TLegend* legendThy =new TLegend( 0.7,0.7,0.9,0.9 );    

    //HERAPDF->SetLineColor(kBlue);
    setupSpectraHist(HERAPDF  ,useSimpBins);
    legendThy->AddEntry(HERAPDF  ,"HERAPDF 2015 NLO","l");
    if(debugMode){fout->cd(); HERAPDF->Write("HERAPDF105_NLO_R04_jtpt");}
    
    //MMHTnlo->SetLineColor(kMagenta);
    setupSpectraHist(MMHTnlo  ,useSimpBins);
    legendThy->AddEntry(MMHTnlo  ,"MMHT 2014 NLO   ","l");
    if(debugMode){fout->cd(); MMHTnlo->Write("MMHT2014_NLO_R04_jtpt"); }
    
    //NNPDFnnlo->SetLineColor(kBlack);
    setupSpectraHist(NNPDFnnlo,useSimpBins);
    legendThy->AddEntry(NNPDFnnlo,"NNPDF NNLO","l");
    if(debugMode){fout->cd(); NNPDFnnlo->Write("NNPDF_NLO_R04_jtpt");    }
    
    if(debugMode)fout->cd();//makeThyHist will make me cd to last thy file opened...
    
    cCheck->cd();
    if(!useSimpBins)cCheck->SetLogx(1);
    cCheck->SetLogy(1);                  
    
    legendThy->AddEntry(hunf_x2,"2 x (OS MC Unf.)","lp");
    legendThy->AddEntry(hgen_rebin_x2,"2 x (SS MC Truth)", "lp");
    
    HERAPDF->SetAxisRange(1e-07,1e+03,"Y");//for y axis in nanbarns
    HERAPDF->SetTitle("NLO Thy w/ MC, OS Unf., SS Truth");
    
    HERAPDF  ->Draw("][HIST E");
    MMHTnlo  ->Draw("][HIST E SAME");
    NNPDFnnlo->Draw("][HIST E SAME");
    
    hgen_rebin_x2->Draw("P E SAME");
    hunf_x2->Draw("P E SAME");   //just for axis range
    

    legendThy->Draw();

    cCheck->Print(outPdfFile.c_str());

    //   --------------------------------------------------
    // thy ratios w hunf------------    
    cCheck->cd();
    if(!useSimpBins)cCheck->SetLogx(1);
    cCheck->SetLogy(0);
    
    setupRatioHist(h_thyratio_CT10nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_CT14nlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_HERAPDF , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_MMHTnlo , useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    setupRatioHist(h_thyratio_NNPDFnnlo, useSimpBins, boundaries_pt_gen_mat, nbins_pt_gen_mat);
    
    h_thyratio_CT10nlo->SetTitle( "Thy Ratios w/ 2 x (SVD Unf. OS MC)" );
    h_thyratio_CT10nlo->GetYaxis()->SetTitle("Thy / 2 x (Unf. OS MC)");
    
    h_thyratio_CT10nlo ->Draw( "][HIST ");      
    h_thyratio_CT14nlo ->Draw( "][HIST SAME"); 
    h_thyratio_HERAPDF ->Draw( "][HIST SAME"); 
    h_thyratio_MMHTnlo ->Draw( "][HIST SAME"); 
    h_thyratio_NNPDFnnlo->Draw("][HIST SAME"); 

    
    TLegend* legendthyrat = new TLegend( 0.1,0.7,0.3,0.9 );
    legendthyrat->AddEntry(h_thyratio_CT10nlo ,  "CT10 PDF NLO" ,    "l");
    legendthyrat->AddEntry(h_thyratio_CT14nlo ,  "CT14 PDF NLO" ,    "l"); 
    legendthyrat->AddEntry(h_thyratio_HERAPDF ,  "HERAPDF 2015 NLO", "l");
    legendthyrat->AddEntry(h_thyratio_MMHTnlo ,  "MMHT 2014 NLO",    "l");
    legendthyrat->AddEntry(h_thyratio_NNPDFnnlo, "NNPDF NNLO",       "l");
    
    legendthyrat->Draw();
    
    theLineAtOne_gen->Draw();
    theLineAtp9_gen->Draw();
    theLineAt1p1_gen->Draw();
 


    cCheck->Print(outPdfFile.c_str());



    //   --------------------------------------------------
    // data-side cov matrix for best kreg choice      
    cCheck->cd();      
    cCheck->SetLogx(0);
    cCheck->SetLogy(0);
    cCheck->SetLogz(1);
    
    matStylePrint( (TH2D*)hCovmatSVD[kRegDraw], hCovmatSVD[kRegDraw]->GetTitle(), cCheck, outPdfFile, true);    

    //   --------------------------------------------------
    // data-side abs cov matrix for best kreg choice      
    cCheck->cd();      
    cCheck->SetLogx(0);
    cCheck->SetLogy(0);
    cCheck->SetLogz(1);
    
    matStylePrint( (TH2D*)hCovmatAbsValSVD[kRegDraw], hCovmatAbsValSVD[kRegDraw]->GetTitle(), cCheck, outPdfFile, true);    
    
    //   --------------------------------------------------
    // data-side pearson matrix for best kreg choice      
    cCheck->cd();      
    cCheck->SetLogx(0);
    cCheck->SetLogy(0);
    cCheck->SetLogz(0);
    
    matStylePrint( (TH2D*)hPearsonSVD[kRegDraw], hPearsonSVD[kRegDraw]->GetTitle(), cCheck, outPdfFile, true);        

    // ----------------------------
    cCheck->Print(close_outPdfFile.c_str());
    
  
    // ------ RESPONSE MATRIX DRAWING    
    if(drawRespMatrix)
      drawRespMatrixFile((TH2D*)hmat, (TH2D*)hmat_rebin, (TH2D*)hmat_errors,
			 //(TH2*)covmat_TH2, (TH2*)PearsonBayes, (TH2*)unfmat_TH2,
			 (double*)boundaries_pt_reco_mat, (int)nbins_pt_reco_mat,
			 (double*)boundaries_pt_gen_mat , (int)nbins_pt_gen_mat,
			 (std::string)outRespMatPdfFile , (bool)useSimpBins, (TFile*)fpp_MC);
    
  } // end drawPDFs
  if(debugMode)std::cout<<std::endl<<"done drawing SVD PDFs!"<<std::endl<<std::endl;
  





  std::cout<<std::endl<<"writing hists to file..."<<std::endl;    
  fout->cd();
  
  
  //output, 3x3 plots
  for(int kr = 0; kr<nKregMax; ++kr) hunf_svd[kr]->Write();            
  for(int kr = 0; kr<nKregMax; ++kr) hfold_svd[kr]->Write();           
  for(int kr = 0; kr<nKregMax; ++kr) hrec_unfolded_ratio[kr]->Write();          
  for(int kr = 0; kr<nKregMax; ++kr) hrec_fold_ratio[kr]->Write();          
  for(int kr = 0; kr<nKregMax; ++kr) hrec_ssfold_ratio[kr]->Write();          
  for(int kr = 0; kr<nKregMax; ++kr) hrec_truthfold_ratio[kr]->Write();          
  for(int kr = 0; kr<nKregMax; ++kr) hgen_unfolded_ratio[kr]->Write();          
  for(int kr = 0; kr<nKregMax; ++kr) hgen_folded_ratio[kr]->Write();              
  for(int kr = 0; kr<nKregMax; ++kr) hPearsonSVD[kr]->Write();
  for(int kr = 0; kr<nKregMax; ++kr) hCovmatSVD[kr]->Write();
  for(int kr = 0; kr<nKregMax; ++kr) hCovmatAbsValSVD[kr]->Write();
  

  //output, other
  hdi->Write();  
  hSVal->Write();
  //assert(false);
  hRegCovMat->Write();
  hRegCovMatInv->Write();
  hDataCovMat->Write();

  // input data ------------------
  hrec_rebin->Write();
  
  // input MC ---------------------
  hgen_rebin->Write();
  hrec_sameside_rebin->Write();
  //hJetQA_jtptEntries->Write();
  
  //response matrix
  hmat->Write();
  hmat_rebin->Write();
  hmat_errors->Write();    
  
  hfak->Write();








  
  if(debugMode)std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  //fpp_Data->Close();
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
			       (int)std::atoi(argv[3]) , 
			       (int)std::atoi(argv[4]) , 
			       (int)std::atoi(argv[5])   ); 
  
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}
