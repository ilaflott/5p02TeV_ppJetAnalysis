#include "unfoldSpectra.h"

//other settings
const int kRegDraw  = 4 ; // array entries w/ arguments 0-8. 4 -> middle hist on 3x3 SVDplot
const int kRegDrawSS = 1;

const bool drawPDFs=true; 
const bool debugMode=false;
const bool drawRespMatrix=true;
const float etaBinWidth=4.;


// CODE --------------------------------------------------
int SVDUnfoldDataSpectra( std::string inFile_Data_dir, std::string inFile_MC_dir , const std::string baseName , 
			  const bool doJetID=true ,  const bool useSimpBins=false, const int kRegCenter= 5 ){
  
  
  
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
  inFile_Data_dir=SCRATCH_BASE+inFile_Data_dir;
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
  const std::string inFile_Data_name="/HighPtJetTrig_"+fullJetType+"-allFiles.root";
  
  
  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldDataSpectra_outdir+fullJetType;
  if(doJetID)outFileName+="_wjtID";
  if(!useSimpBins)outFileName+="_anabins";
  else outFileName+="_simpbins";
  outFileName+="_"+baseName;
  
  std::string outRespMatPdfFile =  outFileName+"_respMat.pdf";
  std::string outSVDPdfFile   =  outFileName+".pdf";  
  std::string out3x3SVDPdfFile   =  outFileName+"_3x3.pdf"; 
  //std::string outSSSVDPdfFile   =  outFileName+"_SS.pdf"; 
  //std::string out3x3SSSVDPdfFile   =  outFileName+"_3x3_SS.pdf"; 
  
  std::string outRootFile     =  outFileName+".root";  
  
  checkNRenameFiles( outFileName, 
		     &outRespMatPdfFile, 
		     &outSVDPdfFile, 
		     &out3x3SVDPdfFile, 
		     &outRootFile );
  
  
  
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
  
  



  // ppData input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from Data file"<<std::endl; 
  std::cout<<"input Data dir : "<< (inFile_Data_dir)  <<std::endl; 
  std::cout<<"Data file name : "<< (inFile_Data_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fpp_Data = TFile::Open( (inFile_Data_dir+inFile_Data_name).c_str());


  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
  std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl;   std::cout<<std::endl<<std::endl;  
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  
  
  
  //for output
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE");   
  fout->cd();  
  
  
  
  // ---------- reco, measured spectra to unfold
  std::string histTitle="hJetQA";
  if(doJetID)histTitle+="_1wJetID";
  else histTitle+="_0wJetID";
  histTitle+="_jtpt";
  //histTitle+=RandEtaRange;
  
  TH1F*  hrec = (TH1F*)fpp_Data->Get( histTitle.c_str() ); 
  hrec->Write();
  if(debugMode)hrec->Print("base");

  histTitle+="_divBylumietabin";
  float effIntgrtdLumi=computeEffLumi(fpp_Data);
  hrec->Scale(1./effIntgrtdLumi); // lumi
  hrec->Scale(1./etaBinWidth); // |y|
  hrec->Write( (histTitle).c_str() );
  if(debugMode)hrec->Print("base");
  
  histTitle+="_clone";
  TH1F *hrec_rebin = (TH1F*)hrec->Clone( (histTitle).c_str() );
  hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");
  
  std::cout<<"rebinning hrec..."<<std::endl;
  histTitle+="_rebins";
  hrec_rebin = (TH1F*)hrec_rebin->Rebin( nbins_pt_reco, (histTitle).c_str() , boundaries_pt_reco);
  hrec_rebin->Write(histTitle.c_str());   
  if(debugMode)hrec_rebin->Print("base");  
  
  histTitle+="_normbinwidth";
  //divideBinWidth(hrec_rebin); 
  hrec_rebin->Write(histTitle.c_str());
  if(debugMode)hrec_rebin->Print("base");  
  
  if(clearOverUnderflows){
    histTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_rebin);
    hrec_rebin->Write(histTitle.c_str());
    if(debugMode)hrec_rebin->Print("base");  
  }
  
  //cosmetics
  hrec_rebin->SetMarkerStyle(kOpenCircle);
  hrec_rebin->SetMarkerColor(kBlue);     
  hrec_rebin->SetMarkerSize(1.02);     
  
  
  
  
  // ---------- reco, measured same-side spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="hpp_reco";
  if(doJetID)histTitle2+="_wJetID";
  histTitle2+=RandEtaRange;
  
  TH1F*  hrec_sameside = (TH1F*)fpp_MC->Get( histTitle2.c_str() ); 
  hrec_sameside->Write(histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  histTitle2+="_divByetabin";
  hrec_sameside->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  hrec_sameside->Write( histTitle2.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  
  histTitle2+="_clone";
  TH1F *hrec_sameside_rebin = (TH1F*)hrec_sameside->Clone( (histTitle2).c_str() );
  hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");
  
  std::cout<<"rebinning hrec_sameside..."<<std::endl;
  histTitle2+="_rebins";
  hrec_sameside_rebin = (TH1F*)hrec_sameside_rebin->Rebin( nbins_pt_reco, (histTitle2).c_str() , boundaries_pt_reco);
  hrec_sameside_rebin->Write( histTitle2.c_str() );   
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
  histTitle2+="_normbinwidth";
  //divideBinWidth(hrec_sameside_rebin); 
  hrec_sameside_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_rebin->Print("base");  
  
  if(clearOverUnderflows){
    histTitle2+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hrec_sameside_rebin);
    hrec_sameside_rebin->Write( histTitle2.c_str() );
    if(debugMode)hrec_sameside_rebin->Print("base");    } 

  //cosmetics
  hrec_sameside_rebin->SetMarkerStyle(kOpenSquare);
  hrec_sameside_rebin->SetMarkerColor(kBlue-3);
  hrec_sameside_rebin->SetMarkerSize(1.02);
  
  TH1F* hrec_sameside_resp_rebin=NULL;
  if(fillRespHists) 
    hrec_sameside_resp_rebin = (TH1F*)hrec_sameside_rebin->Clone("recrebinsamesideClone4unf");
  else
    hrec_sameside_resp_rebin = new TH1F( ("hpp_rec_sameside_response_rebin"+RandEtaRange).c_str(),"", 
					 nbins_pt_reco, boundaries_pt_reco); 
  
  hrec_sameside_resp_rebin->Write( histTitle2.c_str() );
  if(debugMode)hrec_sameside_resp_rebin->Print(" base");  
  
  
  
  
  
  
  
  
  
  
  
  
  // ---------- gen, MC truth spectra
  std::string genHistTitle="hpp_gen";
  if(doJetID)genHistTitle+="_wJetID";
  genHistTitle+=RandEtaRange;
  
  TH1F* hgen = (TH1F*)fpp_MC->Get( genHistTitle.c_str() );
  hgen->Write();
  if(debugMode)hgen->Print("base");    
  
  genHistTitle+="_divByetabin";
  hgen->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  hgen->Write( genHistTitle.c_str());
  if(debugMode)hrec_sameside->Print("base");
  
  genHistTitle+="_clone";
  TH1F* hgen_rebin = (TH1F*)hgen->Clone( (genHistTitle).c_str() );
  hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");
  
  genHistTitle+="_rebins";
  hgen_rebin = (TH1F*)hgen_rebin->Rebin(nbins_pt_gen, (genHistTitle).c_str() , boundaries_pt_gen);
  hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base"); 
  
  genHistTitle+="_normbinwidth";
  //divideBinWidth(hgen_rebin);
  hgen_rebin->Write(genHistTitle.c_str());
  if(debugMode)hgen_rebin->Print("base");  
  
  if(clearOverUnderflows){
    genHistTitle+="_noOverUnderFlows";
    TH1clearOverUnderflows((TH1*)hgen_rebin);
    hgen_rebin->Write(genHistTitle.c_str());
    if(debugMode)hgen_rebin->Print("base");    }

  //cosmetics
  hgen_rebin->SetMarkerStyle(kOpenStar);
  hgen_rebin->SetMarkerColor(kMagenta);
  hgen_rebin->SetMarkerSize(1.02);     	
  
  
  TH1F* hgen_resp_rebin;
  if(fillRespHists) 
    hgen_resp_rebin = (TH1F*)hgen_rebin->Clone("genrebinClone4unf");
  else
    hgen_resp_rebin = new TH1F( ("hpp_gen_response_rebin"+RandEtaRange).c_str() ,"", 
				nbins_pt_gen, boundaries_pt_gen);  

  hgen_resp_rebin->Write();
  if(debugMode)hgen_resp_rebin->Print("base");  
  
  
  
  
  
  // ---------- open MC "response" matrix
  std::string TH2_title="hpp_matrix";//+RandEtaRange;
  if(doJetID)TH2_title+="_wJetID";
  TH2_title+=RandEtaRange;
  
  //get the response matrix made by readforests
  TH2F* hmat = (TH2F*)fpp_MC->Get( TH2_title.c_str() );
  hmat->Write();
  if(debugMode)hmat->Print("base");
  
  TH2_title+="_divByetabin";
  hmat->Scale(1./etaBinWidth); // eta bin width for 0.<|y|<2.
  hmat->Write( TH2_title.c_str());
  if(debugMode)hmat->Print("base");
  
  
  // rebinned matrix ---------------
  TH2_title+="_clone";
  TH2F* hmat_rebin = (TH2F*)hmat->Clone( (TH2_title).c_str() );
  hmat_rebin->Write( TH2_title.c_str() );
  if(debugMode)hmat_rebin->Print("base"); 
  
  TH2_title+="_rebins";
  hmat_rebin=(TH2F*) reBinTH2(hmat_rebin, (TH2_title).c_str(), 
			       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );  
  hmat_rebin->Write( TH2_title.c_str() );
  if(debugMode)hmat_rebin->Print("base"); 
  
  TH2_title+="_normbinwidth";
  //divideBinWidth_TH2(hmat_rebin);
  hmat_rebin->Write( TH2_title.c_str() );
  if(debugMode)hmat_rebin->Print("base"); 
  
  if(clearOverUnderflows){
    TH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2F*)hmat_rebin);
    hmat_rebin->Write( TH2_title.c_str() );
    if(debugMode)hmat_rebin->Print("base");  }

    
  // error and %error for response matrix ---------------
  std::string errTH2_title="hmat_errors_rebins";
  TH2F* hmat_errors=makeRespMatrixErrors( (TH2F*) hmat,
					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  hmat_errors->Write( errTH2_title.c_str() );
  if(debugMode)hmat_errors->Print("base");
  
  errTH2_title+="_normbinwidth";
  //divideBinWidth_TH2(hmat_errors);
  hmat_errors->Write(errTH2_title.c_str());
  if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    errTH2_title+="_noOverUnderFlows";
    TH2clearOverUnderflows((TH2F*)hmat_errors);
    hmat_errors->Write(errTH2_title.c_str());
    if(debugMode)hmat_errors->Print("base");  }
  
  
  // give response matrix the correct errors
  setRespMatrixErrs( (TH2F*)hmat_rebin, (TH2F*) hmat_errors , (bool)zeroBins);
  TH2_title+="_wseterrs";
  hmat_rebin->Write(TH2_title.c_str());
  if(debugMode)hmat_rebin->Print("base");
  
  
  
  if(debugMode)std::cout<<"cloning input hists..."<<std::endl;
  TH1F* hrec_rebin_clone=(TH1F*)hrec_rebin->Clone("hrec_rebin_clone");
  setupHistXAxis(hrec_rebin_clone, useSimpBins);
  if(debugMode)hrec_rebin_clone->Print("base");
  
  TH1F* hgen_rebin_clone=(TH1F*)hgen_rebin->Clone("hgen_rebin_clone");
  setupHistXAxis(hgen_rebin_clone, useSimpBins);
  if(debugMode)hgen_rebin_clone->Print("base");
  
  TH1F* hrec_sameside_rebin_clone=(TH1F*)hrec_sameside_rebin->Clone("hrec_sameside_rebin_clone");
  setupHistXAxis(hrec_sameside_rebin_clone, useSimpBins);
  if(debugMode)hrec_sameside_rebin_clone->Print("base");  
  
  // SVD spectra
  TH1F *hunf_svd[nKregMax]={};        
  TH1F *hfold_svd[nKregMax]={};     
  
  // SVD ratios 
  TH1F *hgen_unfolded_ratio[nKregMax]={};
  TH1F *hrec_unfolded_ratio[nKregMax]={};
  TH1F *hgen_folded_ratio[nKregMax]={};  
  TH1F *hrec_folded_ratio[nKregMax]={};     
  
  TH2D *hPearsonSVD[nKregMax];      
 
  
  //  // MC-side SVD spectra
  //  TH1F *hunf_ss_svd[nKregMax]={};        
  //  TH1F *hfold_ss_svd[nKregMax]={};     
  //  
  //  // MC-side SVD ratios 
  //  TH1F *hgen_ss_unfolded_ratio[nKregMax]={};
  //  TH1F *hrec_ss_unfolded_ratio[nKregMax]={};
  //  TH1F *hgen_ss_folded_ratio[nKregMax]={};  
  //  TH1F *hrec_ss_folded_ratio[nKregMax]={};     
  
  //TH2D *hPearsonSVD_SS[nKregMax];      
  


  
  // TCanvases
  if(debugMode)std::cout<<"creating histos/arrays/canvases for SVD unfolding..."<<std::endl;
  TCanvas *cSpectra               = new TCanvas("cSpectra","",                      1800, 1500);      cSpectra->Divide(3,3);	   
  TCanvas *cRatio_gen             = new TCanvas("cRatio_gen","",                    1800, 1500);      cRatio_gen->Divide(3,3);      		
  TCanvas *cRatio_rec             = new TCanvas("cRatio_rec","",                    1800, 1500);      cRatio_rec->Divide(3,3);      		
  TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","",                1800, 1500);      cPearsonMatrixIter->Divide(3,3);    
  TCanvas *c11                = new TCanvas("c11"," Singular Values and divectors", 1400, 1000);      c11->Divide(2);                    
  
  TLegend *leg[nKregMax],*leg1[nKregMax], *leg2[nKregMax];      


//  // MC-side TCanvases
//  if(debugMode)std::cout<<"creating histos/arrays/canvases for same-side SVD unfolding..."<<std::endl;
//  TCanvas *cSpectra_ss               = new TCanvas("cSpectra_ss","",           1500, 1500);      cSpectra_ss->Divide(3,3);	   
//  TCanvas *cRatio_gen_ss            = new TCanvas("cRatio_gen_ss","",             1500, 1500);      cRatio_gen_ss->Divide(3,3);      	    
//  TCanvas *cRatio_rec_ss             = new TCanvas("cRatio_rec_ss","",             1500, 1500);      cRatio_rec_ss->Divide(3,3);          
//  TCanvas *cPearsonMatrixIter_ss = new TCanvas("cPearsonMatrixIter_ss","", 1500, 1500);      cPearsonMatrixIter_ss->Divide(3,3);   
//  TCanvas *c11_ss                = new TCanvas("c11_ss","S.S. Singular Values and divectors", 1500, 1500);      c11_ss->Divide(2);          
//
//  TLegend *leg_ss[nKregMax],*leg_ss1[nKregMax], *leg_ss2[nKregMax];      

    
  
  
  TLine* theGenLine= new TLine( boundaries_pt_gen_mat[0]   	        , 1.,
				boundaries_pt_gen_mat[nbins_pt_gen_mat] , 1.  ); 
  theGenLine->SetLineWidth(1);
  theGenLine->SetLineStyle(2);
  theGenLine->SetLineColor(36);    
  

  TLine* theLineAtp9_gen= new TLine( boundaries_pt_gen_mat[0]  
				     ,0.9   
				     ,(boundaries_pt_gen_mat[nbins_pt_reco_mat])
				     ,0.9);
  if(debugMode)std::cout<<"theLineAtp9_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
  if(debugMode)std::cout<<"theLineAtp9_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_reco_mat])  <<std::endl;
  theLineAtp9_gen->SetLineWidth(1);
  theLineAtp9_gen->SetLineStyle(2);
  theLineAtp9_gen->SetLineColor(36);
  
  TLine* theLineAt1p1_gen= new TLine( boundaries_pt_gen_mat[0]  
				      ,1.1   
				      ,(boundaries_pt_gen_mat[nbins_pt_reco_mat])
				      ,1.1);
  if(debugMode)std::cout<<"theLineAt1p1_gen starts at "<<boundaries_pt_gen_mat[0]  <<std::endl;
  if(debugMode)std::cout<<"theLineAt1p1_gen ends at   "<<(boundaries_pt_gen_mat[nbins_pt_reco_mat])  <<std::endl;
  theLineAt1p1_gen->SetLineWidth(1);
  theLineAt1p1_gen->SetLineStyle(2);
  theLineAt1p1_gen->SetLineColor(36);
  
  
  TLine* theRecoLine= new TLine( boundaries_pt_reco_mat[0]                 , 1.,
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
  

  //std::cout<<"calling RooUnfoldResponse"<<std::endl;    
  std::string roo_resp_title = "Response_matrix_rebin_"+radius+std::to_string(kRegCenter);
  
  RooUnfoldResponse roo_resp(hrec_sameside_resp_rebin, hgen_resp_rebin, hmat_rebin, (roo_resp_title).c_str());    
  roo_resp.UseOverflow(doOverUnderflows);
  roo_resp.Write();
  

  TH1F* hfak=  (TH1F*) roo_resp.Hfakes() ;
  hfak->Print("base");
  hfak->SetMarkerColor(kGreen);
  hfak->SetMarkerSize(1.02);

  
  
  
  // SVD unfolding loop
  for(int kr = 0; kr < nKregMax; ++kr){
    
    if(debugMode&&kr==0)    std::cout<<"beginning SVD Unfolding Loop..."<<std::endl;
    std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl;        
    
    std::string kRegRandEtaRange_plotTitle="kReg="+std::to_string(kReg[kr]);
    std::string kRegRandEtaRange="_kReg"+std::to_string(kReg[kr]);
    
    
    if(debugMode)std::cout<<std::endl<<"calling RooUnfoldSvd..."<<std::endl<<std::endl;
    RooUnfoldSvd unf_svd(&roo_resp, hrec_rebin, kReg[kr]);
    unf_svd.SetVerbose(2);
    std::cout<<"RooUnfoldSvd Overflow Status: " << unf_svd.Overflow()<<std::endl;
    
    //std::cout<<"calling Ereco... getting Cov Mat"<<std::endl;
    //TMatrixD covmat = unf_svd.Ereco(errorTreatment);
    //
    //std::cout<<"converting TMatrixD covmat to TH2D* covmat_TH2"<<std::endl;
    //TH2D* covmat_TH2= new TH2D(covmat);  //this gives the TH2 a name of "TMatrixDBase" by default.
    //covmat_TH2->SetName(("SVD_covarianceMatrix_TH2_kr"+std::to_string(kr)).c_str());
    //if(debugMode)covmat_TH2->Print("base");
    //covmat_TH2->Write();
    //
    //std::cout<<"getting unfolding matrix"<<std::endl;
    //TH2D* unfmat_TH2=new TH2D(unf_svd.UnfoldingMatrix());
    //unfmat_TH2->SetName( ("SVD_unfoldingMatrix_TH2_kr"+std::to_string(kr)).c_str());  
    //if(debugMode)unfmat_TH2->Print("base");
    //unfmat_TH2->Write();
    

    hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);

    //cosmetics
    hunf_svd[kr]->SetMarkerStyle(kOpenCircle);
    hunf_svd[kr]->SetMarkerColor(kRed);     
    hunf_svd[kr]->SetMarkerSize(1.01);     
          

      
    if(debugMode)std::cout<<"applying roo_resp to histo hunf_svd[kr="<<kr<<"]..."<<std::endl;
    hfold_svd[kr] = (TH1F*)roo_resp.ApplyToTruth(hunf_svd[kr]);
    
    
    //cosmetics
    hfold_svd[kr]->SetMarkerStyle(kOpenCircle);
    hfold_svd[kr]->SetMarkerColor(kGreen);      
    hfold_svd[kr]->SetMarkerSize(1.01);      
    
    
    if(debugMode)std::cout<<"calling Ereco... getting Cov Mat"<<std::endl;
    TMatrixD covmat = unf_svd.Ereco(errorTreatment);
    
    //if(debugMode)std::cout<<"converting TMatrixD covmat to TH2D* covmat_TH2"<<std::endl;
    //TH2D* covmat_TH2= new TH2D(covmat);
    
    hPearsonSVD[kr] = (TH2D*) CalculatePearsonCoefficients(&covmat, false, ("SVD_pearson_"+std::to_string(kr)) );      
    
    //if(debugMode)std::cout<<"creating \"rebinned\" pearson matrix..."<<std::endl;
    //hPearsonSVD[kr] = reBinPearsonTH2( pearson );
    
   





    // ----------------------- 3X3 CANVAS DRAWING  ----------------------- //
    // for data unfolding

    // CAN I DO THIS ELSEWHERE/DIFFERENTLY? //

    // 3x3 spectra canvases  ---------------------
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
    cSpectra->cd(kr+1);
    cSpectra->cd(kr+1)->SetLogy(1);
    if(!useSimpBins)
      cSpectra->cd(kr+1)->SetLogx(1);
    
    //hunf_svd[kr]->SetTitle( ("Reco Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );

    hunf_svd[kr]->SetName( ("hunf_svd"+kRegRandEtaRange).c_str() );
    if(debugMode)hunf_svd[kr]->Print("base");      
    
    hfold_svd[kr]->SetName( ("hfold_svd"+kRegRandEtaRange).c_str() );
    if(debugMode)hfold_svd[kr]->Print("base");
    
    hrec_rebin_clone->SetTitle( ("Jet Spectra, "+kRegRandEtaRange_plotTitle).c_str() );
    hrec_rebin_clone->GetYaxis()->SetTitle("A.U.");
    
    hrec_rebin_clone->Draw("P E");
    hgen_rebin_clone->Draw("P E SAME");
    hunf_svd[kr]->Draw("P E SAME");
    //hfold_svd[kr]->Draw("P E SAME");      
    
    leg[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option 
    leg[kr]->AddEntry(hrec_rebin_clone,"Data Meas.","p");
    leg[kr]->AddEntry(hgen_rebin_clone,"MC Truth","p");
    leg[kr]->AddEntry(hunf_svd[kr],"Data Unf.","p");
    //leg[kr]->AddEntry(hfold_svd[kr],"Data Fold(Unf.)","p");
    leg[kr]->SetTextSize(0.03);
    leg[kr]->Draw();
    
    
    // 3x3 unfORfold/meas canvases  ---------------------
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_rec canvas..."<<std::endl;
    cRatio_rec->cd(kr+1);
    if(!useSimpBins)
      cRatio_rec->cd(kr+1)->SetLogx(1);
    
    hrec_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("hrec_unf_ratio"+kRegRandEtaRange).c_str());
    hrec_unfolded_ratio[kr] = (TH1F*)hrec_unfolded_ratio[kr]->Rebin( nbins_pt_reco, ("hrec_unf_ratio"+kRegRandEtaRange+"_rebin").c_str() , boundaries_pt_reco);      
    hrec_unfolded_ratio[kr]->Divide(hrec_rebin);
    hrec_unfolded_ratio[kr]->SetTitle( ("Data Unf./Data Meas., "+kRegRandEtaRange_plotTitle).c_str() );
    hrec_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
    hrec_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
    if(debugMode)hrec_unfolded_ratio[kr]->Print("base");
    
    
    hrec_folded_ratio[kr] = (TH1F*)hfold_svd[kr]->Clone( ("hrec_fold_ratio"+kRegRandEtaRange).c_str() );
    hrec_folded_ratio[kr] = (TH1F*)hrec_folded_ratio[kr]->Rebin( nbins_pt_reco, ("hrec_fold_ratio"+kRegRandEtaRange+"_rebin").c_str() , boundaries_pt_reco);      
    hrec_folded_ratio[kr]->Divide(hrec_rebin);
    if(debugMode)hrec_folded_ratio[kr]->Print("base");
    
    hrec_unfolded_ratio[kr]->Draw("P E");
    //hrec_folded_ratio[kr]->Draw("P E SAME");
    
    leg1[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option 
    leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],"Data Unf./Data Meas.","pl");
    //leg1[kr]->AddEntry(hrec_folded_ratio[kr],"Data Fold(Unf.)/Data Meas.","pl");
    leg1[kr]->SetTextSize(0.02); 
    leg1[kr]->Draw();
    
    theRecoLine->Draw();
    
    
    
    // draw on cRatio_gen canvas
    if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_gen canvas..."<<std::endl;
    cRatio_gen->cd(kr+1);
    if(!useSimpBins)
      cRatio_gen->cd(kr+1)->SetLogx(1);
    
    hgen_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("hgen_unf_ratio"+kRegRandEtaRange).c_str());
    //hgen_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
    hgen_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
    hgen_unfolded_ratio[kr]->SetTitle( ("Data Unf./MC Truth, "+kRegRandEtaRange_plotTitle).c_str() );
    hgen_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");      
    hgen_unfolded_ratio[kr]->Divide(hgen_rebin);      
    if(debugMode)hgen_unfolded_ratio[kr]->Print("base");
    
    
    TH1F* hgen_rebin_ratClone=(TH1F*)hgen_rebin->Clone("hgen_rebin_ratioClone");
    hgen_rebin_ratClone=(TH1F*)hgen_rebin_ratClone->Rebin( nbins_pt_reco, ("hgen_rebin_ratioClone_rebin") , boundaries_pt_reco);          

    hgen_folded_ratio[kr] = (TH1F*)hfold_svd[kr]->Clone( ("hgen_fold_ratio"+kRegRandEtaRange).c_str());            
    hgen_folded_ratio[kr]->SetTitle( ("SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
    hgen_folded_ratio[kr]->Divide(hgen_rebin_ratClone);      
    if(debugMode)hgen_folded_ratio[kr]->Print("base");      

    
    hgen_unfolded_ratio[kr]->Draw("P E");
    //hgen_folded_ratio[kr]->Draw("P E SAME");
    
    
    leg2[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"BRNDC");
    leg2[kr]->AddEntry(hgen_unfolded_ratio[kr],"Data Unf./MC Truth","pl");
    //leg2[kr]->AddEntry(hgen_folded_ratio[kr],"Data Fold(Unf.)/MC Trtuh","pl");
    //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
    leg2[kr]->SetTextSize(0.02); 
    leg2[kr]->Draw();
    
    theGenLine->Draw();

    // ----------------------- END 3X3 CANVAS DRAWING for data unfolding ----------------------- //
    

    //  singular values and d_i vector ---------------------------
    //Note that these do not depend on the regularization.
    //The opposite: they tell you which regularization to use! (ian note: how?)
    if(kr == kRegDraw){      
      
      std::cout << "  getting divectors and singular values for oppSide ... " << std::endl;
      TSVDUnfold *svdUnfold = unf_svd.Impl();
      //svdUnfold->SetNormalize(true);
      
      // singular values
      if(debugMode)std::cout << "  getting singular values... " << std::endl;
      TH1 *hSVal = (TH1*)svdUnfold->GetSV();
      if(debugMode)
	for(int bin=1; bin<=hSVal->GetNbinsX(); bin++) std::cout<<"bin: "<<bin<<",  SV: "<<hSVal->GetBinContent(bin)<< std::endl;
      
      if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
      c11->cd(1);
      c11->cd(1)->SetLogy();  
      
      //hSVal->SetTitle(" singular values ");
      hSVal->SetXTitle("i");        
      hSVal->SetYTitle("Singular Values s_{i}");        
      hSVal->DrawCopy("HIST E");
      
      // di vector values
      std::cout << "  getting di vector values... " <<  std::endl;
      TH1D *hdi = (TH1D*)svdUnfold->GetD();
      if(debugMode)
	for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
	  std::cout<<"i: "<<bin<<",  di: "<<hdi->GetBinContent(bin)<<std::endl;
      
      c11->cd(2);
      c11->cd(2)->SetLogy();    
      
      hdi->SetXTitle("i");
      hdi->SetYTitle("Divector Values d_{i}");
      hdi->DrawCopy("HIST E"); 
      
      drawText( "5.02 TeV ak4PFJets",0.358173, 0.8459761, 19);
      drawText( "2015 Prompt Reco",0.358173, 0.8159761, 19);
      drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
      drawText( ("Current Best kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);	
      if(debugMode)std::cout<<std::endl<<"done with kr=="<< kRegDraw<<" specifics"<<std::endl<<std::endl;          }


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
    


    //// SS MC unfolding aka sanity check ---------------------------------
    //std::string kRegRandEtaRangess_plotTitle=" kReg"+std::to_string(kRegSS[kr]);//+RandEtaRange_plotTitle;
    //std::string kRegRandEtaRangess="_kReg"+std::to_string(kRegSS[kr]);//+RandEtaRange;
    //
    //// same-side unfolding/refolding --------------
    //std::cout<<std::endl<<"calling RooUnfoldSvd and Hreco for sameside MC test..."<<std::endl<<std::endl;
    ////RooUnfoldSvd unf_ss_svd(&roo_resp, hrec_sameside_rebin, kRegSS[kr]);
    //hunf_ss_svd[kr] = (TH1F*)unf_ss_svd.Hreco(errorTreatment);
    //
    ////cosmetics
    //hunf_ss_svd[kr]->SetMarkerStyle(kOpenSquare);
    //hunf_ss_svd[kr]->SetMarkerColor(kRed);     
    //hunf_ss_svd[kr]->SetMarkerSize(1.01);     
    ////std::cout<<std::endl;         
    //
    //std::cout<<"applying roo_resp to sameside histo hunf_ss_svd[kr="<<kr<<"]..."<<std::endl;
    //hfold_ss_svd[kr] = (TH1F*)roo_resp.ApplyToTruth(hunf_ss_svd[kr]);      
    ////hfold_ss_svd[kr]->Scale(hrec_sameside_resp_rebin->Integral()/hfold_ss_svd[kr]->Integral());
    //
    ////cosmetics
    //hfold_ss_svd[kr]->SetMarkerStyle(kOpenSquare);
    //hfold_ss_svd[kr]->SetMarkerColor(kGreen);      
    //hfold_ss_svd[kr]->SetMarkerSize(1.01);      
    //
    //////////////////////////////////////////////
    //sameside covariance/pearson matrices??? maybe, we will see.
    //TMatrixD covmat_ss = unf_ss_svd.Ereco(errorTreatment);
    //if(debugMode)std::cout<<std::endl;
    
    //TMatrixD *pearson_ss = CalculatePearsonCoefficients(&covmat_ss, false );
    //TMatrixD pearson_ss = CalculatePearsonCoefficients(&covmat_ss, false );
    //hPearsonSVD_SS[kr] = (TH2D*) CalculatePearsonCoefficients(&covmat_ss, false, ("SVD_SS_pearson_"+std::to_string(kr)) );      
    //if(debugMode)std::cout<<std::endl;
    
    //std::cout<<"creating TH2 for pearson matrix..."<<std::endl;
    //hPearsonSVD_SS[kr] = new TH2D (*pearson_ss);
    //hPearsonSVD_SS[kr] = new TH2D (pearson_ss);
    //hPearsonSVD_SS[kr]->SetName(("pearson_sameside_"+kRegRandEtaRangess).c_str());
    //if(debugMode)hPearsonSVD_SS[kr]->Print("base");
    //if(debugMode)std::cout<<std::endl;
    //sameside covariance/pearson matrices???

    
    
    
    
    
    // ----------------------- 3X3 CANVAS DRAWING  ----------------------- //
    // for MC sanity-check unfolding
    
    // 3x3 spectra canvases  ------------------
    //    if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
    //    cSpectra_ss->cd(kr+1);
    //    cSpectra_ss->cd(kr+1)->SetLogy();
    //    if(!useSimpBins)
    //      cSpectra_ss->cd(kr+1)->SetLogx(1);            
    
    //hrec_sameside_rebin_clone->SetTitle( ("S.S. Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );
    
    //hunf_ss_svd[kr]->SetTitle( ("S.S. Unf. Spectra,"+kRegRandEtaRangess_plotTitle).c_str() );
    //hunf_ss_svd[kr]->SetName( ("hunf_ss_svd"+kRegRandEtaRangess).c_str() );
    //if(debugMode)hunf_ss_svd[kr]->Print("base");
    //
    //hfold_ss_svd[kr]->SetTitle( ("S.S. Fold(Unf.) Spectra,"+kRegRandEtaRangess_plotTitle).c_str() );
    //hfold_ss_svd[kr]->SetName( ("hfold_ss_svd"+kRegRandEtaRangess).c_str() );
    //if(debugMode)hfold_ss_svd[kr]->Print("base");
    
    
    //hunf_ss_svd[kr]->Draw("E");
    // hrec_sameside_rebin_clone->Draw("E SAME");
    //hfold_ss_svd[kr]->Draw("E SAME");
    
      
    //leg_ss[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option
    //leg_ss[kr]->AddEntry(hrec_sameside_rebin_clone,"MC Meas.","pl");
    //leg_ss[kr]->AddEntry(hunf_ss_svd[kr],"MC Unf.","pl");
    ////leg_ss[kr]->AddEntry(hfold_ss_svd[kr],"MC Fold(Unf.)","pl");
    //leg_ss[kr]->SetTextSize(0.02);
    //leg_ss[kr]->Draw();

    
    
    
    
    //if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_rec canvas..."<<std::endl;
    //cRatio_rec_ss->cd(kr+1);
    //if(!useSimpBins)
    //cRatio_rec_ss->cd(kr+1)->SetLogx(1);
    
    
    
      
    //hrec_ss_folded_ratio[kr] = (TH1F*)hfold_ss_svd[kr]->Clone( ("hrec_ss_fold_ratio"+kRegRandEtaRangess).c_str() );
    //hrec_ss_folded_ratio[kr]->Divide(hrec_sameside_rebin);
    //if(debugMode)hrec_ss_folded_ratio[kr]->Print("base");
    //
    //
    //hrec_ss_unfolded_ratio[kr] = (TH1F*)hunf_ss_svd[kr]->Clone( ("hrec_ss_unf_ratio"+kRegRandEtaRangess).c_str());
    //hrec_ss_unfolded_ratio[kr] = (TH1F*)hrec_ss_unfolded_ratio[kr]->Rebin (nbins_pt_reco, ("hrec_ss_unf_ratio"+kRegRandEtaRangess+"_rebin").c_str(), boundaries_pt_reco);
    //hrec_ss_unfolded_ratio[kr]->Divide(hrec_sameside_rebin);
    //
    ////hrec_ss_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
    //hrec_ss_unfolded_ratio[kr]->SetTitle( ("S.S. Unf./Meas.,"+kRegRandEtaRangess_plotTitle).c_str() );
    //hrec_ss_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");      
    ////hrec_ss_unfolded_ratio[kr]->SetMarkerStyle(27);
    ////hrec_ss_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
    //if(debugMode)hrec_ss_unfolded_ratio[kr]->Print("base");
    
    //hrec_ss_unfolded_ratio[kr]->Draw("E");
    //hrec_ss_folded_ratio[kr]->Draw("E SAME");
    
    //leg_ss1[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"BRNDC");
    //leg_ss1[kr]->AddEntry(hrec_ss_unfolded_ratio[kr],"MC Unf./MC Meas.","pl");
    //leg_ss1[kr]->AddEntry(hrec_ss_folded_ratio[kr],"MC Fold(Unf.)/MC Meas.","pl");
    //leg_ss1[kr]->SetTextSize(0.02); 
    //leg_ss1[kr]->Draw();
    
    //theRecoLine->Draw();
      
    
    
    //// draw on cRatio_gen canvas
    //if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_gen_ss canvas..."<<std::endl;
    //cRatio_gen_ss->cd(kr+1);
    //if(!useSimpBins)
    //  cRatio_gen_ss->cd(kr+1)->SetLogx(1);
    
    //hgen_ss_unfolded_ratio[kr] = (TH1F*)hunf_ss_svd[kr]->Clone( ("hgen_ss_unf_ratio"+kRegRandEtaRangess).c_str());
    //hgen_ss_unfolded_ratio[kr]->Divide(hgen_rebin);      
    //if(debugMode)hgen_ss_unfolded_ratio[kr]->Print("base");
    //hgen_ss_unfolded_ratio[kr]->SetTitle( ("S.S. Fold(Unf.)/Gen.,"+kRegRandEtaRangess_plotTitle).c_str() );
    //
    //
    //hgen_ss_folded_ratio[kr] = (TH1F*)hfold_ss_svd[kr]->Clone( ("hgen_ss_fold_ratio"+kRegRandEtaRangess).c_str());
    //hgen_ss_folded_ratio[kr]->Divide(hgen_rebin_ratClone);          
    //if(debugMode)      hgen_ss_folded_ratio[kr]->Print("base");
    //hgen_ss_unfolded_ratio[kr]->SetTitle( ("S.S. SVD/Gen.,"+kRegRandEtaRangess_plotTitle).c_str() );
    //
    //hgen_ss_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
    ////hgen_ss_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
    //hgen_ss_unfolded_ratio[kr]->SetAxisRange(0.2, 1.8, "Y");
    
    //hgen_ss_unfolded_ratio[kr]->Draw("P E");
    //hgen_ss_folded_ratio[kr]->Draw("P E SAME");
    
    //leg_ss2[kr] = new TLegend(0.62, 0.75, 0.9, 0.9, NULL,"BRNDC");
    //leg_ss2[kr]->AddEntry(hgen_ss_unfolded_ratio[kr],"MC Unf./MC Truth","pl");
    ////leg_ss2[kr]->AddEntry(hgen_ss_folded_ratio[kr],"MC Fold(Unf.)/MC Truth","pl");
    //leg_ss2[kr]->SetTextSize(0.02); 
    //leg_ss2[kr]->Draw();
    //
    //theGenLine->Draw();
    // ----------------------- END 3X3 CANVAS DRAWING FOR MC SANITY CHECK UNFOLDING  ----------------------- //

      
    
    //if(kr == kRegDrawSS){      
    //  std::cout << "  getting divectors and singular values for sameSide ... " << std::endl;
    //  TSVDUnfold *svdUnfold_ss = unf_ss_svd.Impl();
    //  //svdUnfold->SetNormalize(true);
    //  
    //  // singular values
    //  if(debugMode)std::cout << "  getting singular values... " << std::endl;
    //  TH1 *hSVal_ss = (TH1*)svdUnfold_ss->GetSV();
    //  if(debugMode)
    //	for(int bin=1; bin<=hSVal_ss->GetNbinsX(); bin++) 
    //	  std::cout<<"bin: "<<bin<<",  SV: "<<hSVal_ss->GetBinContent(bin)<< std::endl;
    //  
    //  if(debugMode)std::cout<<std::endl<<"drawing singular values on c11 canvas.."<<std::endl<<std::endl;
    //  c11_ss->cd(1);
    //  c11_ss->cd(1)->SetLogy();    
    //  
    //  //hSVal_ss->SetTitle(" s.s. singular values ");
    //  hSVal_ss->SetXTitle(" s.s singular values ");        
    //  hSVal_ss->DrawCopy("HIST E");
    //  
    //  // di vector values
    //  std::cout << "  getting di vector values... " <<  std::endl;
    //  TH1D *hdi_ss = (TH1D*)svdUnfold_ss->GetD();
    //  if(debugMode)
    //	for(int bin=1; bin<=hdi_ss->GetNbinsX(); bin++)
    //	  std::cout<<"i: "<<bin<<",  di: "<<hdi_ss->GetBinContent(bin)<<std::endl;
    //  
    //  c11_ss->cd(2);
    //  c11_ss->cd(2)->SetLogy();    
    //  
    //  //hdi_ss->SetTitle(" s.s. di vectors ");
    //  hdi_ss->SetXTitle(" |d_{i}^{kreg}| ");
    //  hdi_ss->DrawCopy("HIST E"); 
    //  
    //  drawText( "5.02 TeV ak4PFJets",0.358173, 0.8459761, 19);
    //  //drawText( "ppData, Prompt Reco HighPtJets",0.358173, 0.8159761, 19);
    //  drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    //  //drawText( ("Best kReg="+std::to_string(kReg[kr])).c_str() ,0.408173, 0.7359761, 19);	
    //  if(debugMode)std::cout<<std::endl<<"done with kr=="<< kRegDrawSS<<" specifics"<<std::endl<<std::endl;      	
    //}//end drawing singular/divector values
    
    
    
    
    
    //// Pearson Drawing
    //if(debugMode)std::cout<<std::endl<<"drawing stuff on cPearsonMatrixIter_SS canvas..."<<std::endl;
    //cPearsonMatrixIter_ss->cd(kr+1);
    //
    //hPearsonSVD_SS[kr]->SetName( ("hPearsonSVDSS"+kRegRandEtaRange).c_str() );      
    //hPearsonSVD_SS[kr]->SetTitle( ("pearsonMatrix SVD SS"+kRegRandEtaRange_plotTitle).c_str()  );
    //hPearsonSVD_SS[kr]->GetXaxis()->SetTitle("reco bin num");
    //hPearsonSVD_SS[kr]->GetYaxis()->SetTitle("gen bin num");
    //
    //hPearsonSVD_SS[kr]->GetZaxis()->SetLabelSize(0.035);
    //hPearsonSVD_SS[kr]->SetAxisRange(-1., 1., "Z");    
    //hPearsonSVD_SS[kr]->Draw("COLZ");
    
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
    
    
    //hunf_ss_svd[kr]->Write();            
    //hfold_ss_svd[kr]->Write();               
    //hrec_ss_unfolded_ratio[kr]->Write();          
    //hrec_ss_folded_ratio[kr]->Write();          
    //hgen_ss_unfolded_ratio[kr]->Write();          
    //hgen_ss_folded_ratio[kr]->Write();              
    //hPearsonSVD_SS[kr]->Write();

  }
  
  if(drawPDFs){
    
    std::cout<<std::endl<<"drawing PDFs for SVD Unfolding..."<<std::endl;
    
    
    // form filename string, open pdf file and draw the canvases we have so far
    std::string out3x3PdfFile=out3x3SVDPdfFile;
    std::string open_out3x3PdfFile=out3x3PdfFile+"[";      std::string close_out3x3PdfFile=out3x3PdfFile+"]";
    
    //open 3x3 pdf file
    cSpectra->Print(open_out3x3PdfFile.c_str());
    
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
    
    
    //close 3x3 pdf file
    cPearsonMatrixIter->Print(close_out3x3PdfFile.c_str());    
    
    
    
    
    // form filename string, open pdf file and draw the canvases we have so far
    std::string outPdfFile=outSVDPdfFile;
    std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
    
    //   --------------------------------------------------
    // now the specific kr=kRegDraw plots
    TCanvas *cCheck        = new TCanvas("cCheck","",        1400, 1000);          
    cCheck->Print(open_outPdfFile.c_str());

    //   --------------------------------------------------
    // spectra for best kreg choice
    cCheck->cd();
    cCheck->SetLogy(1);
    if(!useSimpBins)	
      cCheck->SetLogx(1);
    
    hrec_rebin_clone->SetTitle("SVD MC and Data Spectra");
    setupHistXAxis(hrec_rebin_clone, useSimpBins, boundaries_pt_reco, nbins_pt_reco);
    
    hrec_rebin_clone->Draw("P E");
    hunf_svd[kRegDraw]->Draw("P E SAME");
    hgen_rebin_clone->Draw("P E SAME");	
    //hfold_svd[kRegDraw]->Draw("P E SAME");
    
    
    TLegend * leg1 = new TLegend(0.7, 0.7, 0.9, 0.9, NULL,"NBNDC");
    leg1->AddEntry(hrec_rebin_clone,   "Data Meas","p");
    leg1->AddEntry(hunf_svd[kRegDraw], "Data Unf.","p");
    leg1->AddEntry(hgen_rebin_clone,   "MC Truth","p");	
    //leg1->AddEntry(hfold_svd[kRegDraw],"Data Fold","p");
    //leg1->SetTextSize(0.02);
    leg1->Draw();
    
    //drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
    //drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
    //drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    
    cCheck->Print(outPdfFile.c_str());
    
    
    //   --------------------------------------------------
    //  MC response matrix fakes
    cCheck->cd();
    cCheck->SetLogy(1);
    if(!useSimpBins)
      cCheck->SetLogx(1);
    
    hrec_sameside_rebin->SetTitle("MC Only, Response Spectra and Fakes");
    setupHistXAxis(hrec_sameside_rebin, useSimpBins, boundaries_pt_reco, nbins_pt_reco);
    hrec_sameside_rebin->GetYaxis()->SetTitle("A.U.");
    
    hrec_sameside_rebin->Draw("P E");
    hgen_resp_rebin->Draw("P E SAME");
    hfak->Draw("P E SAME");
    
    TLegend * legfake = new TLegend(0.6, 0.7, 0.9, 0.9, NULL,"NBNDC");
    legfake->AddEntry(hfak,"MC Response Matrix Fakes","lp");
    legfake->AddEntry(hgen_resp_rebin,"MC Truth Response" ,"lp");
    legfake->AddEntry(hrec_sameside_rebin,"MC Meas. Response" ,"lp");
    //legfake->SetTextSize(0.02);
    legfake->Draw();
    
//    drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
//    drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
//    drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    
    cCheck->Print(outPdfFile.c_str());
    
    
    
    
    //  ratios w/ MC truth --------------------------------------------------
    cCheck->cd();
    cCheck->SetLogy(0);
    if(!useSimpBins)
      cCheck->SetLogx(1);
    
    hgen_unfolded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");
    hgen_unfolded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ MC Truth Spectra");            
    hgen_unfolded_ratio[kRegDraw]->GetYaxis()->SetTitle("Ratio w/ MC Truth");

    setupHistXAxis(hgen_unfolded_ratio[kRegDraw], useSimpBins);
    
    
    hgen_unfolded_ratio[kRegDraw]->Draw("P E");
    //hgen_folded_ratio[kRegDraw]->Draw("P E SAME");
    
    TLegend * leg3 = new TLegend(0.1, 0.8, 0.2, 0.9, NULL,"NBNDC");
    leg3->AddEntry(hgen_unfolded_ratio[kRegDraw],"Data Unf.","pl");
    //leg3->AddEntry(hgen_folded_ratio[kRegDraw],"MC Fold(Unf.)/MC Truth","pl");
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
    
    hrec_unfolded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");
    hrec_unfolded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ Data Meas.");    
    hrec_unfolded_ratio[kRegDraw]->GetYaxis()->SetTitle("Ratio w/ Data Meas.");
    setupHistXAxis(hrec_unfolded_ratio[kRegDraw], useSimpBins);
    
    
    hrec_unfolded_ratio[kRegDraw]->Draw("P E");
    //hrec_folded_ratio[kRegDraw]->Draw("P E SAME");
    
    TLegend * leg2 = new TLegend(0.1, 0.8, 0.2, 0.9, NULL,"NBNDC");
    leg2->AddEntry(hrec_unfolded_ratio[kRegDraw],"Data Unf.","pl");
    //leg2->AddEntry(hrec_folded_ratio[kRegDraw],"Data Fold(Unf.)/Data Meas.","pl");
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
    
    
    
    
    
    
    
    
    
    
    //   --------------------------------------------------
    // data-side pearson matrix for best kreg choice      
    cCheck->cd();      
    cCheck->SetLogx(0);
    cCheck->SetLogy(0);
    
    hPearsonSVD[kRegDraw]->Draw("COLZ");            
    cCheck->Print(outPdfFile.c_str());
    

    //  data-side singular values+divectors --------------------------------------------------
    c11->cd()                ;   
    c11->Print(outPdfFile.c_str());
    
    
    
    
    
    
    // ----------------------------
    cCheck->Print(close_outPdfFile.c_str());
    
  
    // ------ RESPONSE MATRIX DRAWING    
    if(drawRespMatrix)
      drawRespMatrixFile((TH2*)hmat, (TH2*)hmat_rebin, (TH2*)hmat_errors,
			 //(TH2*)covmat_TH2, (TH2*)PearsonBayes, (TH2*)unfmat_TH2,
			 (double*)boundaries_pt_reco_mat, (int)nbins_pt_reco_mat,
			 (double*)boundaries_pt_gen_mat , (int)nbins_pt_gen_mat,
			 (std::string)outRespMatPdfFile , (bool)useSimpBins);
    
  } // end drawPDFs
  if(debugMode)std::cout<<std::endl<<"done drawing SVD PDFs!"<<std::endl<<std::endl;

  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();
  
  if(debugMode)std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_Data->Close();
  fpp_MC->Close();
  return 0;
    
} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=7 ){
    std::cout<<"do ./SVDUnfoldDataSpectra.exe <targDataDir> <targMCDir> <baseOutputName> <doJetID> <useSimpleBins> <kRegCenter>"<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;

    return rStatus;  }
  
  rStatus=1; // runtime error
  
  rStatus=SVDUnfoldDataSpectra(  (const std::string)argv[1], (const std::string)argv[2], (const std::string)argv[3], 
			       (int)std::atoi(argv[4]) , (int)std::atoi(argv[5]) ,  (int)std::atoi(argv[6]) 
			       ); 
  
  std::cout<<std::endl<<"done!"<<std::endl<<" return status: "<<rStatus<<std::endl<<std::endl;
  return rStatus;
}

























    //    //SS SVD UNFOLDING
    //    {
    //      // form filename string, open pdf file and draw the canvases we have so far
    //      std::string outPdfFile=outSSSVDPdfFile;
    //      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
    //      
    //      //std::string out3x3PdfFile=out3x3SSSVDPdfFile;
    //      //	std::string open_out3x3PdfFile=out3x3SSSVDPdfFile+"[";      std::string close_out3x3PdfFile=out3x3SSSVDPdfFile+"]";
    //      
    //      //  already-drawn sameside canvases --------------------------------------------------
    //      c11->Print(open_outPdfFile.c_str());
    //      c11_ss->cd()                ;   c11_ss->Print(outPdfFile.c_str());
    //      cPearsonMatrixIter_ss->cd() ;   cPearsonMatrixIter_ss->Print(outPdfFile.c_str());
    //      cSpectra_ss->cd()           ;   cSpectra_ss->Print(outPdfFile.c_str());
    //      
    //      
    //      //  cRatio Check for meas --------------------------------------------------
    //      cRatio_rec_ss->cd()         ;   cRatio_rec_ss->Print(outPdfFile.c_str());   //   cRatio_rec_ss->Print(out3x3PdfFile.c_str());      
    //      
    //      TCanvas *cRatioCheck_rec_ss        = new TCanvas("cRatioCheck_rec_ss","",        1800, 1500);    
    //      cRatioCheck_rec_ss->cd();
    //      if(!useSimpBins)cRatioCheck_rec_ss->SetLogx(1);
    //      //cRatioCheck_rec_ss->SetLogx(1);
    //      
    //      hrec_ss_folded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");
    //      //hrec_ss_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
    //      hrec_ss_folded_ratio[kRegDraw]->SetTitle("S.S. SVD Unf. Ratios w/ Meas.");
    //      
    //      hrec_ss_unfolded_ratio[kRegDraw]->Draw("P E");
    //      //hrec_ss_folded_ratio[kRegDraw]->Draw("P E SAME");
    //      
    //      TLegend * leg_ss2 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
    //      leg_ss2->AddEntry(hrec_ss_unfolded_ratio[kRegDraw],"Data Unf./Data Meas.","pl");
    //      //leg_ss2->AddEntry(hrec_ss_folded_ratio[kRegDraw],"Data Fold(Unf.)/Data Meas.","pl");
    //      leg_ss2->SetTextSize(0.02);
    //      leg_ss2->Draw();
    //      
    //      //drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
    //      //drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
    //      drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    //      drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    //      
    //      //TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
    //      //theLine->SetLineWidth(1);
    //      //theLine->SetLineStyle(2);
    //      //theLine->SetLineColor(36);
    //      theRecoLine->Draw();
    //      
    //      cRatioCheck_rec_ss->Print(outPdfFile.c_str());
    //      
    //      
    //      //  cRatio Check for gen --------------------------------------------------
    //      cRatio_gen_ss->cd()         ;   cRatio_gen_ss->Print(outPdfFile.c_str());	//  cRatio_gen_ss->Print(out3x3PdfFile.c_str());	
    //      TCanvas *cRatioCheck_gen_ss        = new TCanvas("cRatioCheck_gen_ss","",        1800, 1500);    
    //      cRatioCheck_gen_ss->cd();
    //      if(!useSimpBins)cRatioCheck_gen_ss->SetLogx(1);
    //      //cRatioCheck_gen_ss->SetLogx(1);
    //      
    //      hgen_ss_folded_ratio[kRegDraw]->SetAxisRange(0.2, 1.8, "Y");//SetAxisRange(0.5, 1.5, "Y");
    //	//hgen_ss_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
    //      hgen_ss_folded_ratio[kRegDraw]->SetTitle("S.S. SVD Unf. Ratios w/ Gen Truth");
    //      
    //      hgen_ss_unfolded_ratio[kRegDraw]->Draw("P E");
    //      //hgen_ss_folded_ratio[kRegDraw]->Draw("P E SAME");
    //      
    //      TLegend * leg_ss3 = new TLegend(0.14, 0.79, 0.34, 0.87, NULL,"NBNDC");
    //      leg_ss3->AddEntry(hgen_ss_unfolded_ratio[kRegDraw],"MC Unf./MC Truth","pl");
    //      //leg_ss3->AddEntry(hgen_ss_folded_ratio[kRegDraw],"MC Fold(Unf.)/MC Truth","pl");
    //      leg_ss3->SetTextSize(0.02);
    //      leg_ss3->Draw();
    //      
    //      //drawText( "5.02 TeV pp, ak4PF Jets",          0.14, 0.75, 22);
    //      //drawText( "Prompt-Reco, Jet80+LowerJets",     0.14, 0.72, 22);
    //      drawText( MCdesc.c_str(),0.358173, 0.7859761, 19);
    //      drawText( ("kReg="+std::to_string(kReg[kRegDraw])).c_str(), 0.14, 0.69, 22);
    //      
    //      //TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
    //      //theLine->SetLineWidth(1);
    //      //theLine->SetLineStyle(2);
    //      //theLine->SetLineColor(36);
    //      theGenLine->Draw();
    //      
    //      cRatioCheck_gen_ss->Print(outPdfFile.c_str());
    //      
    //      cRatioCheck_gen_ss->Print(close_outPdfFile.c_str());
    //      //cRatioCheck_gen_ss->Print(close_out3x3PdfFile.c_str());
    //    }
    //std::cout<<std::endl<<"done drawing SS SVD PDFs!"<<std::endl<<std::endl;
