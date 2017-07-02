#include "unfoldSpectra.h"

const int kRegDraw  = 4 ; // array entries w/ arguments 0-8. 4 -> middle hist on 3x3 SVDplot
const bool drawPDFs=true; 

//other settings
const bool doToyErrs=true;
const bool debugMode=false;



// CODE --------------------------------------------------
int unfoldMCSpectra( std::string inFile_MC_dir , const std::string baseName , 
		     const bool doJetID=true , const int kRegCenter= 5 ,
		     const bool doBayes=true, const bool doSVD=false, const bool useSimplePtBinning=false){
  
  
  
  // BINNING -----------
  double* boundaries_pt_reco ;   int nbins_pt_reco ;   // reco spectra bins
  double* boundaries_pt_gen  ;   int nbins_pt_gen  ;   // gen spectra bins
  
  double* boundaries_pt_reco_mat ;   int nbins_pt_reco_mat ;   // matrix bins
  double* boundaries_pt_gen_mat  ;   int nbins_pt_gen_mat  ;   // 

  //int Nbins2Clear_gen_lowpt=0;
  //int Nbins2Clear_gen_highpt=0;
  //int Nbins2Clear_reco_lowpt=0;
  //int Nbins2Clear_reco_highpt=0;


  if(useSimplePtBinning){    std::cout<<"using simple pt bins"<<std::endl<<std::endl;

    boundaries_pt_reco = (double*)simpbins_pt_reco ; //spectra
    nbins_pt_reco      = (int)n_simpbins_pt_reco ; 
    boundaries_pt_gen = (double*)simpbins_pt_gen ;
    nbins_pt_gen      = (int)n_simpbins_pt_gen ; 
    
    boundaries_pt_reco_mat = (double*)simpbins_pt_reco ; //matrix
    nbins_pt_reco_mat      = (int)n_simpbins_pt_reco ; 
    boundaries_pt_gen_mat = (double*)simpbins_pt_gen ;
    nbins_pt_gen_mat      = (int)n_simpbins_pt_gen ;                    

    //Nbins2Clear_gen_lowpt=0;
    //Nbins2Clear_gen_highpt=0;
    //Nbins2Clear_reco_lowpt=0;
    //Nbins2Clear_reco_highpt=0;          
  } 

    else{   std::cout<<"using analysis pt bins"<<std::endl<<std::endl;

    boundaries_pt_reco = (double*)anabins_pt_reco ; //spectra
    nbins_pt_reco      = (int)n_anabins_pt_reco ; 
    boundaries_pt_gen = (double*)anabins_pt_gen ;
    nbins_pt_gen      = (int)n_anabins_pt_gen ; 

    boundaries_pt_reco_mat = (double*)anabins_pt_reco ; //matrix
    nbins_pt_reco_mat      = (int)n_anabins_pt_reco ; 
    boundaries_pt_gen_mat = (double*)anabins_pt_gen ;
    nbins_pt_gen_mat      = (int)n_anabins_pt_gen ;                  

    //Nbins2Clear_gen_lowpt=0;
    //Nbins2Clear_gen_highpt=0;
    //Nbins2Clear_reco_lowpt=0;
    //Nbins2Clear_reco_highpt=0;      
    }



  // STRINGS -----------
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl; 
  //inFile_Data_dir=SCRATCH_BASE+inFile_Data_dir;
  inFile_MC_dir  =SCRATCH_BASE+inFile_MC_dir;
  
  std::size_t radPos=inFile_MC_dir.find("_ak")+3;
  const std::string radiusInt= inFile_MC_dir.substr( radPos,1 );
  const std::string radius="R"+radiusInt;//R4
  if(debugMode)std::cout<<"radius string is = "<<radius<<std::endl;
  const std::string RandEtaRange="_"+radius+"_20_eta_20";//_R4_20_eta_20
  if(debugMode)std::cout<<"RandEtaRange string is = "<<RandEtaRange<<std::endl;
  
  const std::string Rstring="_"+radius; //_R4
  const std::string Rstring_plotTitle=" "+radius; // R4
  const std::string RandEtaRange_plotTitle=Rstring_plotTitle+" 20eta20";// R4 20eta20
  
  std::size_t jetTypePos=radPos+1, jetsPos=inFile_MC_dir.find("Jets");
  const std::string jetType=inFile_MC_dir.substr( jetTypePos,(jetsPos-jetTypePos) );//PFJets
  if(debugMode)std::cout<<"jetType string is = "<<jetType<<std::endl;

  const std::string fullJetType="ak"+radiusInt+jetType;//"ak4PFJets"
  if(debugMode)std::cout<<"fullJetType string is = "<<fullJetType<<std::endl;
  


  // INFILE NAME(S) -----------
  const std::string inFile_MC_name="/Py8_CUETP8M1_QCDjetAllPtBins_"+fullJetType+"-allFiles.root";
  //const std::string inFile_Data_name="/HighPtJetTrig_"+fullJetType+"-allFiles.root";

  // OUTPUT FILE, NAME(S) -----------
  std::string outFileName=unfoldMCSpectra_outdir+fullJetType;//+"_"+baseName;//+".root";  
  if(doJetID)outFileName+="_wjtID";//+".root";  
  if(!useSimplePtBinning)outFileName+="_anabins";
  else outFileName+="_simpbins";
  outFileName+="_"+baseName;

  std::string outRespMatPdfFile =  outFileName+"_respMat.pdf";
  std::string outBayesPdfFile =  outFileName+"_Bayes.pdf";
  std::string outSVDPdfFile   =  outFileName+"_SVD.pdf"; // see drawPDFs part for rest of strin
  std::string outSSSVDPdfFile   =  outFileName+"_SSSVD.pdf"; // see drawPDFs part for rest of string
  std::string outRootFile     =  outFileName+"_File.root";  
  
  
  if(debugMode)std::cout<<"opening output file: "<<outRootFile<<std::endl;
  TFile* fout = new TFile(outRootFile.c_str(),"RECREATE"); 



  // ERROR/WEIGHTS/STATS HANDLING ------------------
  RooUnfold::ErrorTreatment errorTreatment;
  if(!doToyErrs) errorTreatment = RooUnfold::kCovariance;
  else errorTreatment = RooUnfold::kCovToy; 

  if(debugMode)std::cout<<"doToyErrs="<<doToyErrs<<std::endl; 
   
  std::cout<<"TH2 GetDefaultSumw2="<<TH2::GetDefaultSumw2()<<std::endl;
  std::cout<<"setting TH2DefSumw2 to true..."<<std::endl; 
  TH2::SetDefaultSumw2(true);
  
  std::cout<<"TH1 GetDefaultSumw2="<<TH1::GetDefaultSumw2()<<std::endl;
  std::cout<<"setting TH1DefSumw2 to true..."<<std::endl; 
  TH1::SetDefaultSumw2(true);
  
  gStyle->SetOptStat(0);



  // ppMC input histos -------------------------
  std::cout<<std::endl<<std::endl<<"opening INPUT histos from MC file"<<std::endl; 
  if(debugMode){
    std::cout<<"input MC dir : "<< (inFile_MC_dir)  <<std::endl; 
    std::cout<<"MC file name : "<< (inFile_MC_name)<<std::endl; 
    std::cout<<std::endl<<std::endl;  }
  TFile *fpp_MC = TFile::Open( (inFile_MC_dir+inFile_MC_name).c_str());
  


  // ---------- reco, measured spectra to unfold
  std::string histTitle="hpp_mcclosure_reco_test";
  if(doJetID)histTitle+="_wJetID";
  histTitle+=RandEtaRange;
  
  TH1F *hrec;
  hrec = (TH1F*)fpp_MC->Get( histTitle.c_str() ); 

  if(debugMode){fout->cd();
    hrec->Write("hpp_mcclosure_reco_test_forChecking");
    hrec->Print("base");}

  TH1F *hrec_anabin;
  hrec_anabin = (TH1F*)hrec->Clone( (histTitle+"_clone").c_str() );
  //hrec_anabin->Scale(1/effIntgrtdLumi_vz);
  if(debugMode)hrec_anabin->Print("base");

  if(debugMode)std::cout<<"rebinning hrec..."<<std::endl;
  hrec_anabin = (TH1F*)hrec_anabin->Rebin(nbins_pt_reco, (histTitle+"_anabins").c_str() , boundaries_pt_reco);

  if(debugMode)hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned");
  if(debugMode)hrec_anabin->Print("base");  
  
  divideBinWidth(hrec_anabin); 
  if(debugMode)hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned_divBinWidth");
  if(debugMode)hrec_anabin->Print("base");  
  
  if(clearOverUnderflows){
    TH1clearOverUnderflows((TH1*)hrec_anabin);
    if(debugMode)hrec_anabin->Write("hpp_mcclosure_reco_test_rebinned_divBinWidth_noOverUnderFlows");
    if(debugMode)hrec_anabin->Print("base");  
  }

  // response hist, for output? what is this for if it's empty?
  TH1F* hrec_resp_anabin;
  if(fillRespHists) hrec_resp_anabin = (TH1F*)hrec_anabin->Clone("recanabinClone4unf");
  else{
    hrec_resp_anabin = new TH1F( ("hpp_rec_response_anabin"+RandEtaRange).c_str(),"", 
				 nbins_pt_reco, boundaries_pt_reco);}
  if(debugMode)hrec_resp_anabin->Print(" base");  

  

  // ---------- reco, measured spectra used to create response matrix, and for "sameside" unfolding test
  std::string histTitle2="hpp_mcclosure_reco";
  if(doJetID)histTitle2+="_wJetID";
  histTitle2+=RandEtaRange;
  
  TH1F *hrec_sameside;
  hrec_sameside = (TH1F*)fpp_MC->Get( histTitle2.c_str() ); 
  if(debugMode)hrec_sameside->Print("base");  
  
  TH1F *hrec_sameside_anabin;
  hrec_sameside_anabin = (TH1F*)hrec_sameside->Clone( (histTitle2+"_clone").c_str() );
  //hrec_sameside_anabin->Scale(1/effIntgrtdLumi_vz);
  if(debugMode)hrec_sameside_anabin->Print("base");
  
  if(debugMode)std::cout<<"rebinning hrec_sameside..."<<std::endl;
  hrec_sameside_anabin = (TH1F*)hrec_sameside_anabin->Rebin(nbins_pt_reco, (histTitle2+"_anabins").c_str() , boundaries_pt_reco); 
  if(debugMode)hrec_sameside_anabin->Print("base");
  
  divideBinWidth(hrec_sameside_anabin); 
  if(debugMode)hrec_sameside_anabin->Print("base");
  
  if(clearOverUnderflows){
    TH1clearOverUnderflows((TH1*)hrec_sameside_anabin);
    if(debugMode)hrec_sameside_anabin->Print("base");
  }

  if(doMCIntegralScaling){
    std::cout<<"scaling hrec_anabin to hrec_sameside_anabin"<<std::endl;
    float integral_rec_anabin=hrec_anabin->Integral();
    float integral_rec_ss_anabin=hrec_sameside_anabin->Integral();
    std::cout<<"scale factor="<<integral_rec_ss_anabin/integral_rec_anabin<<std::endl;
    hrec_anabin->Scale(integral_rec_ss_anabin/integral_rec_anabin);
  }
  


  // ---------- gen, MC truth spectra
  std::string genHistTitle="hpp_mcclosure_gen";
  if(doJetID)genHistTitle+="_wJetID";
  genHistTitle+=RandEtaRange;
  
  TH1F* hgen = (TH1F*)fpp_MC->Get( genHistTitle.c_str() );
  if(debugMode)hgen->Print("base");    
  
  TH1F* hgen_anabin = (TH1F*)hgen->Clone( (genHistTitle+"_clone").c_str() );
  hgen_anabin = (TH1F*)hgen_anabin->Rebin(nbins_pt_gen, (genHistTitle+"_anabins").c_str() , boundaries_pt_gen);
  if(debugMode)hgen_anabin->Print("base"); 
  
  divideBinWidth(hgen_anabin);
  if(debugMode)hgen_anabin->Print("base");  
  
  if(clearOverUnderflows){
    TH1clearOverUnderflows((TH1*)hgen_anabin);
    if(debugMode)hgen_anabin->Print("base");  
  }

  TH1F* hgen_resp_anabin;
  if(fillRespHists) hgen_resp_anabin = (TH1F*)hgen_anabin->Clone("genanabinClone4unf");
  else{
    hgen_resp_anabin = new TH1F( ("hpp_gen_response_anabin"+RandEtaRange).c_str() ,"", 
				 nbins_pt_gen, boundaries_pt_gen);
  }
  if(debugMode)hgen_resp_anabin->Print("base");  
  



  std::cout<<std::endl<<"writing input hists to file..."<<std::endl;
  
  fout->cd();

  hgen->Write(); 
  hrec->Write();
  hgen_anabin->Write(); 
  hrec_anabin->Write(); 
  
  hrec_sameside->Write();
  hrec_sameside_anabin->Write(); 
  
  if(fillRespHists) hgen_resp_anabin->Write();
  if(fillRespHists) hrec_resp_anabin->Write();
  


  
  
  // ---------- open MC "response" matrix
  std::string TH2_title="hpp_mcclosure_matrix";//+RandEtaRange;
  if(doJetID)TH2_title+="_wJetID";
  TH2_title+=RandEtaRange;
  
  //get the response matrix made by readforests
  TH2F* hmat = (TH2F*)fpp_MC->Get( TH2_title.c_str() );

  // cd to fout so i can save the hist
  if(debugMode){fout->cd();
    hmat->Write();
    hmat->Print("base");}



  // rebinned matrix ---------------
  TH2F* hmat_anabin = (TH2F*)hmat->Clone( (TH2_title+"_clone").c_str() );
  if(debugMode)hmat_anabin->Print("base"); 

  hmat_anabin=(TH2F*) reBinTH2(hmat_anabin, (TH2_title+"_anabins").c_str(), 
			       (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
			       (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );  
  if(debugMode)hmat_anabin->Write("hmat_rebinned");
  if(debugMode)hmat_anabin->Print("base"); 
  
  divideBinWidth_TH2(hmat_anabin);
  if(debugMode)hmat_anabin->Write("hmat_rebinned_binWidthDiv");
  if(debugMode)hmat_anabin->Print("base"); 

  if(clearOverUnderflows){
    TH2clearOverUnderflows((TH2F*)hmat_anabin);
    if(debugMode)hmat_anabin->Write("hmat_rebinned_binWidthDiv_noOverUnderFlow");
    if(debugMode)hmat_anabin->Print("base");
  }
  
  
  // error and %error for response matrix ---------------
  TH2F* hmat_errors=makeRespMatrixErrors( (TH2F*) hmat,
					  (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
					  (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  if(debugMode)hmat_errors->Write("hmat_errors_inanabins");
  if(debugMode)hmat_errors->Print("base");
  
  divideBinWidth_TH2(hmat_errors);
  if(debugMode)hmat_errors->Write("hmat_errors_inanabins_divbinWidth");
  if(debugMode)hmat_errors->Print("base");
  
  if(clearOverUnderflows){
    TH2clearOverUnderflows((TH2F*)hmat_errors);
    if(debugMode)hmat_errors->Write("hmat_errors_inanabins_divbinWidth_noOverUnder");
    if(debugMode)hmat_errors->Print("base");
  }
  
  // give response matrix the correct errors
  setRespMatrixErrs( (TH2F*)hmat_anabin, (TH2F*) hmat_errors );
  
  TH2F* hmat_percenterrs= makeRespMatrixPercentErrs( (TH2F*) hmat_errors, (TH2F*) hmat_anabin,
						     (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
						     (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );		     
  if(debugMode)hmat_percenterrs->Print("base");
  if(debugMode)  hmat_percenterrs->Write("hmat_percentErrors_inanabins");
  

    
  // row/col normd matrix, has correct errors since hmat_anabin has correct errors ---------------
  // POTENTIAL ISSUE: using the resp matrix post rebinning/divide bin width/clearing overflows... should i be using the original matrix?
  TH2F* hmat_anabin_colnormd = normalizeCol_RespMatrix( (TH2F*)  hmat_anabin,
							(double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
							(double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  if(debugMode)  hmat_anabin_colnormd->Write("hmat_rebinned_binWidthDiv_noOverUnderFlow_colNormd");
  if(debugMode)  hmat_anabin_colnormd->Print("base");
  
  TH2F*  hmat_anabin_rownormd = normalizeRow_RespMatrix( (TH2F*)  hmat_anabin,
							 (double*) boundaries_pt_reco_mat, nbins_pt_reco_mat,
							 (double*) boundaries_pt_gen_mat, nbins_pt_gen_mat  );
  if(debugMode)  hmat_anabin_rownormd->Write("hmat_rebinned_binWidthDiv_noOverUnderFlow_rowNormd");
  if(debugMode)  hmat_anabin_rownormd->Print("base");


  if(drawPDFs){    
    
    bool drawRespMatrix=true;
    //bool drawInputSpectra=true;
    
    if(drawRespMatrix){
      std::cout<<std::endl<<"drawing input response matrices..."<<std::endl;
      
      std::string outPdfFile=outRespMatPdfFile;
      std::string open_outPdfFile=outPdfFile+"[";      std::string close_outPdfFile=outPdfFile+"]";
      
      TCanvas* tempCanvForPdfPrint=new TCanvas("tempCanv_respMat","",1200,1200);    
      tempCanvForPdfPrint->cd();

      if(useSimplePtBinning){
	tempCanvForPdfPrint->SetLogx(0);
	tempCanvForPdfPrint->SetLogy(0);       
	tempCanvForPdfPrint->SetLogz(1);          }
      else {
	tempCanvForPdfPrint->SetLogx(1);
	tempCanvForPdfPrint->SetLogy(1);       
	tempCanvForPdfPrint->SetLogz(1);         } 

  
      tempCanvForPdfPrint->Print(open_outPdfFile.c_str()); 


      // general for drawRespMatrix ---------------
      tempCanvForPdfPrint->cd();

      if(!useSimplePtBinning){
      	hmat->GetYaxis()->SetMoreLogLabels(true);
	hmat->GetYaxis()->SetNoExponent(true);	
	hmat->GetXaxis()->SetMoreLogLabels(true);
	hmat->GetXaxis()->SetNoExponent(true);      }
      
      hmat->GetZaxis()->SetLabelSize(0.025);      
      hmat->GetYaxis()->SetLabelSize(0.02);
      hmat->GetYaxis()->SetTitleSize(0.023);
      hmat->GetYaxis()->SetTitle("gen p_{t}");      
      hmat->GetXaxis()->SetLabelSize(0.02);
      hmat->GetXaxis()->SetTitleSize(0.025);
      hmat->GetXaxis()->SetTitle("reco p_{t}");

      // input resp matrix w/ full range ---------------
      tempCanvForPdfPrint->cd();
      
      hmat->SetTitle("ppMC Resp Matrix, original");
      
      hmat->SetAxisRange(0.,1000.,"X");
      hmat->SetAxisRange(0.,1000.,"Y");
      hmat->SetAxisRange(0.000000000000000001,.001,"Z");
      
      hmat->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());
      
      // matrix w/ specific range ---------------
      
      tempCanvForPdfPrint->cd();
      
      hmat->SetTitle("ppMC Resp Matrix, used pt range");
      hmat->SetAxisRange(boundaries_pt_reco_mat[0],boundaries_pt_reco_mat[nbins_pt_reco_mat],"X");
      hmat->SetAxisRange(boundaries_pt_gen_mat[0],boundaries_pt_gen_mat[nbins_pt_gen_mat],"Y");
      
      hmat->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());
      

      // matrix rebinned ---------------
      
      tempCanvForPdfPrint->cd();
      if(!useSimplePtBinning){
      	hmat_anabin->GetYaxis()->SetMoreLogLabels(true);
	hmat_anabin->GetYaxis()->SetNoExponent(true);	
	hmat_anabin->GetXaxis()->SetMoreLogLabels(true);
	hmat_anabin->GetXaxis()->SetNoExponent(true);      
      }      
      
      hmat_anabin->GetZaxis()->SetLabelSize(0.025);      
      hmat_anabin->GetYaxis()->SetLabelSize(0.02);
      hmat_anabin->GetYaxis()->SetTitleSize(0.023);
      hmat_anabin->GetYaxis()->SetTitle("gen p_{t}");      
      hmat_anabin->GetXaxis()->SetLabelSize(0.02);
      hmat_anabin->GetXaxis()->SetTitleSize(0.025);
      hmat_anabin->GetXaxis()->SetTitle("reco p_{t}");
      hmat_anabin->SetAxisRange(0.000000000000000001,.001,"Z");
      hmat_anabin->SetTitle("ppMC Resp Matrix rebinned");
      
      hmat_anabin->Draw("COLZ");           
      
      tempCanvForPdfPrint->Print(outPdfFile.c_str());

      // error matrix in binning of interest ---------------
      
      tempCanvForPdfPrint->cd();
      
      if(!useSimplePtBinning){
      	hmat_errors->GetYaxis()->SetMoreLogLabels(true);
	hmat_errors->GetYaxis()->SetNoExponent(true);	
	hmat_errors->GetXaxis()->SetMoreLogLabels(true);
	hmat_errors->GetXaxis()->SetNoExponent(true);      
      }      
      
      hmat_errors->GetZaxis()->SetLabelSize(0.025);      
      hmat_errors->GetYaxis()->SetLabelSize(0.02);
      hmat_errors->GetYaxis()->SetTitleSize(0.023);
      hmat_errors->GetYaxis()->SetTitle("gen p_{t}");      
      hmat_errors->GetXaxis()->SetLabelSize(0.02);
      hmat_errors->GetXaxis()->SetTitleSize(0.025);
      hmat_errors->GetXaxis()->SetTitle("reco p_{t}");

      
      hmat_errors->SetAxisRange(0.000000000000000001,.001,"Z");      
      hmat_errors->SetTitle("ppMC Resp Matrix Errors");
      
      hmat_errors->Draw("COLZ");           

      tempCanvForPdfPrint->Print(outPdfFile.c_str());
      
      // percent error matrix in binning of interest ---------------
      
      tempCanvForPdfPrint->cd();
      if(!useSimplePtBinning){
      	hmat_percenterrs->GetYaxis()->SetMoreLogLabels(true);
	hmat_percenterrs->GetYaxis()->SetNoExponent(true);	
	hmat_percenterrs->GetXaxis()->SetMoreLogLabels(true);
	hmat_percenterrs->GetXaxis()->SetNoExponent(true);      
      }      
      
      hmat_percenterrs->GetZaxis()->SetLabelSize(0.025);      
      hmat_percenterrs->GetYaxis()->SetLabelSize(0.02);
      hmat_percenterrs->GetYaxis()->SetTitleSize(0.023);
      hmat_percenterrs->GetYaxis()->SetTitle("gen p_{t}");      
      hmat_percenterrs->GetXaxis()->SetLabelSize(0.02);
      hmat_percenterrs->GetXaxis()->SetTitleSize(0.025);
      hmat_percenterrs->GetXaxis()->SetTitle("reco p_{t}");
      
      hmat_percenterrs->SetAxisRange(0.001,200.,"Z");      
      hmat_percenterrs->SetTitle("Resp Matrix % Errors");
      hmat_percenterrs->Draw("COLZ");                 

      tempCanvForPdfPrint->Print(outPdfFile.c_str());


      // col normd matrix in binning of interest  ---------------
  
      tempCanvForPdfPrint->cd();
      
      if(!useSimplePtBinning){
      	hmat_anabin_colnormd->GetYaxis()->SetMoreLogLabels(true);
	hmat_anabin_colnormd->GetYaxis()->SetNoExponent(true);	
	hmat_anabin_colnormd->GetXaxis()->SetMoreLogLabels(true);
	hmat_anabin_colnormd->GetXaxis()->SetNoExponent(true);      
      }      
      
      hmat_anabin_colnormd->GetZaxis()->SetLabelSize(0.025);      
      hmat_anabin_colnormd->GetYaxis()->SetLabelSize(0.02);
      hmat_anabin_colnormd->GetYaxis()->SetTitleSize(0.023);
      hmat_anabin_colnormd->GetYaxis()->SetTitle("gen p_{t}");      
      hmat_anabin_colnormd->GetXaxis()->SetLabelSize(0.02);
      hmat_anabin_colnormd->GetXaxis()->SetTitleSize(0.025);
      hmat_anabin_colnormd->GetXaxis()->SetTitle("reco p_{t}");


      hmat_anabin_colnormd->SetAxisRange(0.00001,10,"Z");            
      hmat_anabin_colnormd->SetTitle("Resp Matrix, Columns Normallized");
      hmat_anabin_colnormd->Draw("COLZ");                 

      tempCanvForPdfPrint->Print(outPdfFile.c_str());


      // row normd matrix in binning of interest  ---------------
  
      tempCanvForPdfPrint->cd();
      
      if(!useSimplePtBinning){
      	hmat_anabin_rownormd->GetYaxis()->SetMoreLogLabels(true);
	hmat_anabin_rownormd->GetYaxis()->SetNoExponent(true);	
	hmat_anabin_rownormd->GetXaxis()->SetMoreLogLabels(true);
	hmat_anabin_rownormd->GetXaxis()->SetNoExponent(true);      
      }      
      
      hmat_anabin_rownormd->GetZaxis()->SetLabelSize(0.025);      
      hmat_anabin_rownormd->GetYaxis()->SetLabelSize(0.02);
      hmat_anabin_rownormd->GetYaxis()->SetTitleSize(0.023);
      hmat_anabin_rownormd->GetYaxis()->SetTitle("gen p_{t}");      
      hmat_anabin_rownormd->GetXaxis()->SetLabelSize(0.02);
      hmat_anabin_rownormd->GetXaxis()->SetTitleSize(0.025);
      hmat_anabin_rownormd->GetXaxis()->SetTitle("reco p_{t}");


      hmat_anabin_rownormd->SetAxisRange(0.00001,10,"Z");            
      hmat_anabin_rownormd->SetTitle("Resp Matrix, Rows Normallized");
      hmat_anabin_rownormd->Draw("COLZ");                 

      tempCanvForPdfPrint->Print(outPdfFile.c_str());

      // close file ---------------

      tempCanvForPdfPrint->Print(close_outPdfFile.c_str()); }
  }
  













  // Bayesian unfolding ------------------------- 
  if(!doBayes)  std::cout<<std::endl<<"   skipping Bayesian Unfolding..."<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<std::endl<<"   beginning Bayesian Unfolding..."<<std::endl;
    

    // RooUnfoldResponse, Bayes --------------------
    std::cout<<"calling RooUnfoldResponse "<<std::endl;
    RooUnfoldResponse roo_resp( hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix"+RandEtaRange).c_str()) ;
    roo_resp.UseOverflow(doOverUnderflows);    
        

    std::cout<<"calling RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_bayes( &roo_resp, hrec_anabin, kIter );

    TH1F *hunf = (TH1F*)unf_bayes.Hreco(errorTreatment);     std::cout<<std::endl; 
    hunf->SetName("ppMC_BayesUnf_Spectra");
    hunf->SetTitle("ppMC BayesUnf Spectra");
    if(debugMode)hunf->Print("base");
         
    TH1F *hfold = (TH1F*)roo_resp.ApplyToTruth(hunf);    std::cout<<std::endl; 
    hfold->SetName("ppMC_BayesFold_Spectra");
    hfold->SetTitle("ppMC BayesFold Spectra");
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
    std::cout<<"calling sameside RooUnfoldBayes..."<<std::endl;
    RooUnfoldBayes unf_sameside_bayes( &roo_resp, hrec_sameside_anabin, kIter );
    
    TH1F *hunf_sameside = (TH1F*)unf_sameside_bayes.Hreco(errorTreatment);    std::cout<<std::endl; 
    hunf_sameside->SetName("ppMC_BayesUnf_sameSideSpectra");
    hunf_sameside->SetTitle("ppMC BayesUnf sameSideSpectra");
    if(debugMode)hunf_sameside->Print("base");

    TH1F *hfold_sameside = (TH1F*)roo_resp.ApplyToTruth(hunf_sameside);    std::cout<<std::endl; 
    hfold_sameside->SetName("ppMC_BayesFold_sameSide_Spectra");
    hfold_sameside->SetTitle("ppMC BayesFold sameSide Spectra");
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
      

      TLine* theLineAtOne= new TLine( boundaries_pt_gen_mat[0],1.,boundaries_pt_gen_mat[nbins_pt_gen_mat],1.);
      theLineAtOne->SetLineWidth(1);
      theLineAtOne->SetLineStyle(2);
      theLineAtOne->SetLineColor(36);

      
      // gen/meas/unf spectra ---------------
      
      canvForPrint->cd();
      canvForPrint->SetLogx(1);
      canvForPrint->SetLogy(1);
      
      hgen_anabin->SetTitle("ppMC Spectra, input and output");
      hgen_anabin->GetXaxis()->SetTitle("jet p_{T} (GeV)");
      hgen_anabin->GetXaxis()->SetMoreLogLabels(true);
      hgen_anabin->GetXaxis()->SetNoExponent(true);
      hgen_anabin->SetAxisRange(boundaries_pt_gen_mat[0], boundaries_pt_gen_mat[nbins_pt_gen_mat],"X");           

      hrec_anabin->SetMarkerStyle(kOpenTriangleUp);
      hrec_anabin->SetMarkerColor(kBlue);     
      hrec_anabin->SetMarkerSize(1.02);     
      
      hunf->SetMarkerStyle(kOpenCircle);
      hunf->SetMarkerColor(kRed);
      hunf->SetMarkerSize(1.02);     
      
      hrec_sameside_anabin->SetMarkerStyle(kOpenTriangleDown);
      hrec_sameside_anabin->SetMarkerColor(kBlue-3);     
      hrec_sameside_anabin->SetMarkerSize(1.02);     
      
      hunf_sameside->SetMarkerStyle(kOpenSquare);
      hunf_sameside->SetMarkerColor(kRed-3);
      hunf_sameside->SetMarkerSize(1.02);     
      
      hgen_anabin->SetMarkerStyle(kOpenCross);
      hgen_anabin->SetMarkerColor(kGreen-2);
      hgen_anabin->SetMarkerSize(1.02);     	
      
      hgen_anabin->Draw();           
      hrec_anabin->Draw("SAME");           
      hunf->Draw("SAME");           
      hunf_sameside->Draw("SAME");           
      hrec_sameside_anabin->Draw("SAME");           
      
      
      TLegend* legend_in = new TLegend( 0.6,0.8,0.9,0.9 );
      legend_in->AddEntry(hrec_anabin,          "Meas.", "p");	
      legend_in->AddEntry(hunf,                    "Unf.", "p");
      legend_in->AddEntry(hrec_sameside_anabin, "Meas., Same Side", "p");
      legend_in->AddEntry(hunf_sameside,           "Unf., Same Side", "p");
      legend_in->AddEntry(hgen_anabin,          "Gen. Truth", "p");
      legend_in->Draw();
      
      canvForPrint->Print(outPdfFile.c_str());      
     

      //---------------
      
      canvForPrint->cd();
      canvForPrint->SetLogx(1);
      canvForPrint->SetLogy(0);

      h_genratio_oppunf->SetAxisRange(0., 2., "Y");
      h_genratio_oppunf->GetXaxis()->SetMoreLogLabels(true);
      h_genratio_oppunf->GetXaxis()->SetNoExponent(true);
      h_genratio_oppunf->SetAxisRange(boundaries_pt_gen_mat[0], boundaries_pt_gen_mat[nbins_pt_gen_mat],"X");           

      h_genratio_oppunf->Draw();
      h_genratio_oppfold->Draw("SAME");
      h_genratio_oppmeas->Draw("SAME");
      h_genratio_ssunf->Draw("SAME");
      h_genratio_ssfold->Draw("SAME");
      h_genratio_ssmeas->Draw("SAME");

      TLegend* legend2 = new TLegend( 0.6,0.8,0.9,0.9 );
      legend2->AddEntry(h_genratio_oppunf, "Unf./Gen. Truth", "p");
      legend2->AddEntry(h_genratio_oppfold, NULL, "p");
      legend2->AddEntry(h_genratio_oppmeas, NULL, "p");
      legend2->AddEntry(h_genratio_ssunf, "Same Side Unf./Gen. Truth", "p");
      legend2->AddEntry(h_genratio_ssfold, NULL, "p");
      legend2->AddEntry(h_genratio_ssmeas, NULL, "p");

      legend2->Draw();

      theLineAtOne->Draw();

      canvForPrint->Print(outPdfFile.c_str());

      //---------------
      	
      canvForPrint->cd();
      canvForPrint->SetLogx(1);
      canvForPrint->SetLogy(0);

      h_recratio_oppunf->SetAxisRange(0., 2., "Y");
      h_recratio_oppunf->GetXaxis()->SetMoreLogLabels(true);
      h_recratio_oppunf->GetXaxis()->SetNoExponent(true);
      h_recratio_oppunf->SetAxisRange(boundaries_pt_gen_mat[0], boundaries_pt_gen_mat[nbins_pt_gen_mat],"X");           

      h_recratio_oppunf->Draw();
      h_recratio_oppfold->Draw("SAME");
      h_recratio_ssunf->Draw("SAME");
      h_recratio_ssfold->Draw("SAME");
      h_recratio_ssgen->Draw("SAME");
      
      TLegend* legend3 = new TLegend( 0.6,0.8,0.9,0.9 );
      legend3->AddEntry(h_recratio_oppunf, NULL, "p");
      legend3->AddEntry(h_recratio_oppfold, NULL, "p");
      legend3->AddEntry(h_recratio_ssunf, NULL, "p");
      legend3->AddEntry(h_recratio_ssfold, NULL, "p");
      legend3->AddEntry(h_recratio_ssgen, NULL, "p");
      legend3->Draw();
      
      theLineAtOne->Draw();
      
      canvForPrint->Print(outPdfFile.c_str());

      //---------------
	
    
      canvForPrint->Print(close_outPdfFile.c_str());      

    }// end draw pdfs
  }  // end doBayes

  std::cout<<std::endl<<"Bayesian Unfolding done, all output written!"<<std::endl;



  // SVD unfolding ------------------------------------------------------------------------- 
  if(!doSVD) std::cout<<std::endl<<"   skipping SVD unfolding..."<<std::endl<<std::endl;
  else{ std::cout<<std::endl<<std::endl<<std::endl<<"   beginning SVD unfolding..."<<std::endl;


    // setup RooUnfoldResponse, and histogram arrays
    std::cout<<"calling RooUnfoldResponse"<<std::endl;

    RooUnfoldResponse roo_resp(hrec_resp_anabin, hgen_resp_anabin, hmat_anabin, ("Response_matrix_anabin"+Rstring).c_str());    
    roo_resp.UseOverflow(doOverUnderflows);


    std::cout<<std::endl<<"creating histos+histo-arrays+canvases for SVD unfolding..."<<std::endl;

    // prep for drawing
    if(debugMode)std::cout<<"opening TCanvases..."<<std::endl;
    TCanvas *cSpectra               = new TCanvas("cSpectra","",           1500, 1500);      cSpectra->Divide(3,3);	   
    TCanvas *cRatio_gen             = new TCanvas("cRatio_gen","",             1500, 1500);      cRatio_gen->Divide(3,3);      		
    TCanvas *cRatio_rec             = new TCanvas("cRatio_rec","",             1500, 1500);      cRatio_rec->Divide(3,3);      		
    TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","", 1500, 1500);      cPearsonMatrixIter->Divide(3,3);    
    TCanvas *c11                = new TCanvas("c11"," Singular Values and divectors", 1500, 1500);      c11->Divide(2);                    
    
   
    if(debugMode)std::cout<<"cloning input hists..."<<std::endl;
    TH1F* hrec_anabin_clone=(TH1F*)hrec_anabin->Clone("hrec_anabin_clone");
    hrec_anabin_clone->SetXTitle("Jet p_{T} (GeV/c)");
    hrec_anabin_clone->SetMarkerStyle(24);
    hrec_anabin_clone->SetMarkerColor(kBlack);
    hrec_anabin_clone->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
    hrec_anabin_clone->Print("base");
    
    TH1F* hgen_anabin_clone=(TH1F*)hgen_anabin->Clone("hgen_anabin_clone");
    hgen_anabin_clone->SetXTitle("Jet p_{T} (GeV/c)");
    hgen_anabin_clone->SetMarkerStyle(24);
    hgen_anabin_clone->SetMarkerColor(kBlack);
    hgen_anabin_clone->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
    hgen_anabin_clone->Print("base");
    

    // unfolded/folded measured spectra
    TH1F *hunf_svd[nKregMax];        
    TH1F *hfold_svd[nKregMax];     
    
    // ratios of unf/folded meas. spectra w/ gen/meas
    TH1F *hgen_unfolded_ratio[nKregMax];
    TH1F *hgen_folded_ratio[nKregMax];  
    TH1F *hrec_folded_ratio[nKregMax];     
    TH1F *hrec_unfolded_ratio[nKregMax];

    //TH2D *hPearsonSVDPriorMeas[nKregMax];      
    TH2D *hPearsonSVD[nKregMax];      




    std::cout<<std::endl<<"creating histos+histo-arrays+canvases for same-side SVD unfolding..."<<std::endl;

    // for same side stuff
    TCanvas *cSpectra_ss               = new TCanvas("cSpectra_ss","",           1500, 1500);      cSpectra_ss->Divide(3,3);	   
    TCanvas *cRatio_gen_ss            = new TCanvas("cRatio_gen_ss","",             1500, 1500);      cRatio_gen_ss->Divide(3,3);      		
    TCanvas *cRatio_rec_ss             = new TCanvas("cRatio_rec_ss","",             1500, 1500);      cRatio_rec_ss->Divide(3,3);      		
    TCanvas *cPearsonMatrixIter_ss = new TCanvas("cPearsonMatrixIter_ss","", 1500, 1500);      cPearsonMatrixIter_ss->Divide(3,3);   
    TCanvas *c11_ss                = new TCanvas("c11_ss","S.S. Singular Values and divectors", 1500, 1500);      c11_ss->Divide(2);                    

    TH1F* hrec_sameside_anabin_clone=(TH1F*)hrec_sameside_anabin->Clone("hrec_sameside_anabin_clone");
    hrec_sameside_anabin_clone->SetXTitle("Jet p_{T} (GeV/c)");
    hrec_sameside_anabin_clone->SetMarkerStyle(24);
    hrec_sameside_anabin_clone->SetMarkerColor(kBlack);
    hrec_sameside_anabin_clone->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
    hrec_sameside_anabin_clone->Print("base");
    
    // unfolded/folded same-side measured spectra
    TH1F *hunf_ss_svd[nKregMax];        
    TH1F *hfold_ss_svd[nKregMax];     
    //
    //// ratios of ss unf/folded meas. spectra w/ gen/meas
    TH1F *hgen_ss_unfolded_ratio[nKregMax];
    TH1F *hgen_ss_folded_ratio[nKregMax];  
    TH1F *hrec_ss_folded_ratio[nKregMax];     
    TH1F *hrec_ss_unfolded_ratio[nKregMax];

    TH2D *hPearsonSVD_SS[nKregMax];      





    // prep for svd unfolding loop
    if(debugMode) std::cout<<std::endl<<"initializing kReg parameter array w/ "<<nKregMax<<" elements"<<std::endl;
    if(debugMode) std::cout<<kRegCenter-kRegRange<<" <= kReg[i] <= "<<kRegCenter+kRegRange<<std::endl;
    
    int kReg[nKregMax];
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) kReg[i+kRegRange]=kRegCenter+i;
    for(int i=(-1*kRegRange); (i+kRegRange)<nKregMax; ++i) 
      std::cout<<"kReg["<<i+kRegRange<<"]="<<kReg[i+kRegRange]<<std::endl;

    TLegend *leg[nKregMax],*leg1[nKregMax], *leg2[nKregMax];      
    TLegend *leg_ss[nKregMax],*leg_ss1[nKregMax], *leg_ss2[nKregMax];      

    TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
    theLine->SetLineWidth(1);
    theLine->SetLineStyle(2);
    theLine->SetLineColor(36);    


    // SVD unfolding loop
    for(int kr = 0; kr < nKregMax; ++kr){
      if(debugMode&&kr==0)    std::cout<<"beginning SVD Unfolding Loop..."<<std::endl;
      std::cout<<std::endl<<std::endl<<"kr="<<kr<<" , kReg = "<<kReg[kr]<<std::endl<<std::endl;
      
      std::string kRegRandEtaRange_plotTitle=" kReg"+std::to_string(kReg[kr]);//+RandEtaRange_plotTitle;
      std::string kRegRandEtaRange="_kReg"+std::to_string(kReg[kr]);//+RandEtaRange;
         

      // opp-side unfolding/refolding --------------

      // unfold the MC dist, then create a ratio with the gen-level MC dist
      std::cout<<std::endl<<"calling RooUnfoldSvd and Hreco..."<<std::endl<<std::endl;
      RooUnfoldSvd unf_svd(&roo_resp, hrec_anabin, kReg[kr]);
      hunf_svd[kr] = (TH1F*)unf_svd.Hreco(errorTreatment);
      std::cout<<std::endl;         
      
      std::cout<<"applying roo_resp to histo hunf_svd[kr="<<kr<<"]..."<<std::endl;
      hfold_svd[kr] = (TH1F*)roo_resp.ApplyToTruth(hunf_svd[kr]);
      

      // get covariance matrix and calculate pearson coefficients
      std::cout<<std::endl<<"calling Ereco..."<<std::endl;
      TMatrixD covmat = unf_svd.Ereco(errorTreatment);
      if(debugMode)std::cout<<std::endl;
      
      //const bool debugPearson=(false && debugMode) ;        TMatrixD *pearson = CalculatePearsonCoefficients(&covmat, debugPearson );
      TMatrixD *pearson = CalculatePearsonCoefficients(&covmat, false );
      if(debugMode)std::cout<<std::endl;
      
      std::cout<<"creating TH2 for pearson matrix..."<<std::endl;
      hPearsonSVD[kr] = new TH2D (*pearson);
      hPearsonSVD[kr]->SetName(("pearson_oppside_"+kRegRandEtaRange).c_str());
      hPearsonSVD[kr]->Print("base");
      //std::cout<<"creating \"rebinned\" pearson matrix..."<<std::endl;
      //hPearsonSVD[kr] = reBinPearsonTH2( pearson );
      

      // 3x3 spectra canvases  ---------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
      cSpectra->cd(kr+1);
      cSpectra->cd(kr+1)->SetLogy();
      cSpectra->cd(kr+1)->SetLogx();

      //hrec_anabin_clone->SetTitle( ("Reco Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );
     
      hunf_svd[kr]->SetTitle( ("Reco Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );
      hunf_svd[kr]->SetName( ("hunf_svd"+kRegRandEtaRange).c_str() );
      hunf_svd[kr]->SetLineStyle(33);
      hunf_svd[kr]->SetLineColor(kBlue);
      hunf_svd[kr]->Print("base");      

      hfold_svd[kr]->SetName( ("hfold_svd"+kRegRandEtaRange).c_str() );
      hfold_svd[kr]->SetLineColor(kRed);      
      hfold_svd[kr]->Print("base");
      
      hunf_svd[kr]->Draw();
      hrec_anabin_clone->Draw("same");
      hfold_svd[kr]->Draw("same");      
      
      leg[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option
      leg[kr]->AddEntry(hrec_anabin_clone,"Measured","pl");
      leg[kr]->AddEntry(hunf_svd[kr],"Unfolded","pl");
      leg[kr]->AddEntry(hfold_svd[kr],"Folded","pl");
      leg[kr]->SetTextSize(0.04);
      leg[kr]->Draw();


      // 3x3 unfORfold/meas canvases  ---------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_rec canvas..."<<std::endl;
      cRatio_rec->cd(kr+1);
      cRatio_rec->cd(kr+1)->SetLogx();
      
      hrec_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("hrec_unf_ratio"+kRegRandEtaRange).c_str());
      hrec_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
      hrec_unfolded_ratio[kr]->Divide(hrec_anabin);
      hrec_unfolded_ratio[kr]->Print("base");
      hrec_unfolded_ratio[kr]->SetTitle( ("SVD/Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hrec_unfolded_ratio[kr]->SetAxisRange(0.5, 1.5, "Y");


      hrec_folded_ratio[kr] = (TH1F*)hfold_svd[kr]->Clone( ("hrec_fold_ratio"+kRegRandEtaRange).c_str() );
      hrec_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_folded_ratio[kr]->SetMarkerColor(kRed);
      hrec_folded_ratio[kr]->Divide(hrec_anabin);
      hrec_folded_ratio[kr]->Print("base");


      hrec_unfolded_ratio[kr]->Draw();
      hrec_folded_ratio[kr]->Draw("SAME");

      leg1[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg1[kr]->AddEntry(hrec_unfolded_ratio[kr],"Unf./Meas.","pl");
      leg1[kr]->AddEntry(hrec_folded_ratio[kr],"reFold/Meas.","pl");
      leg1[kr]->SetTextSize(0.04); 
      leg1[kr]->Draw();
      
      theLine->Draw();



      // draw on cRatio_gen canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_gen canvas..."<<std::endl;
      cRatio_gen->cd(kr+1);
      cRatio_gen->cd(kr+1)->SetLogx();

      hgen_unfolded_ratio[kr] = (TH1F*)hunf_svd[kr]->Clone( ("hgen_unf_ratio"+kRegRandEtaRange).c_str());
      hgen_unfolded_ratio[kr]->Divide(hgen_anabin);      
      hgen_unfolded_ratio[kr]->SetMarkerStyle(33);
      hgen_unfolded_ratio[kr]->SetMarkerColor(kRed);
      hgen_unfolded_ratio[kr]->Print("base");
      hgen_unfolded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hgen_unfolded_ratio[kr]->SetAxisRange(0.5, 1.5, "Y");
      hgen_unfolded_ratio[kr]->SetTitle( ("SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      hgen_unfolded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");      

      hgen_folded_ratio[kr] = (TH1F*)hfold_svd[kr]->Clone( ("hgen_fold_ratio"+kRegRandEtaRange).c_str());
      hgen_folded_ratio[kr]->Divide(hgen_anabin);      
      hgen_folded_ratio[kr]->SetMarkerStyle(24);
      hgen_folded_ratio[kr]->SetMarkerColor(kBlue);
      hgen_folded_ratio[kr]->Print("base");      
      hgen_folded_ratio[kr]->SetTitle( ("SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      
      hgen_folded_ratio[kr]->Draw();
      hgen_unfolded_ratio[kr]->Draw("same");
      
      
      leg2[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg2[kr]->AddEntry(hgen_unfolded_ratio[kr],"Unf./Gen","pl");
      leg2[kr]->AddEntry(hgen_folded_ratio[kr],"Folded/Gen","pl");
      //leg->AddEntry(hSVD_prior,"Prior, normalized to data","pl");
      leg2[kr]->SetTextSize(0.04); 
      leg2[kr]->Draw();
      
      theLine->Draw();








      // same-side unfolding/refolding --------------
      std::cout<<std::endl<<"calling RooUnfoldSvd and Hreco for sameside MC test..."<<std::endl<<std::endl;
      RooUnfoldSvd unf_ss_svd(&roo_resp, hrec_sameside_anabin, kReg[kr]);
      hunf_ss_svd[kr] = (TH1F*)unf_ss_svd.Hreco(errorTreatment);
      std::cout<<std::endl;         
      
      std::cout<<"applying roo_resp to sameside histo hunf_ss_svd[kr="<<kr<<"]..."<<std::endl;
      hfold_ss_svd[kr] = (TH1F*)roo_resp.ApplyToTruth(hunf_ss_svd[kr]);
      
      //////////////////////////////////////////////
      //sameside covariance/pearson matrices??? maybe, we will see.
      TMatrixD covmat_ss = unf_ss_svd.Ereco(errorTreatment);
      if(debugMode)std::cout<<std::endl;
      
      TMatrixD *pearson_ss = CalculatePearsonCoefficients(&covmat_ss, false );
      if(debugMode)std::cout<<std::endl;
      
      std::cout<<"creating TH2 for pearson matrix..."<<std::endl;
      hPearsonSVD_SS[kr] = new TH2D (*pearson_ss);
      hPearsonSVD_SS[kr]->SetName(("pearson_sameside_"+kRegRandEtaRange).c_str());
      hPearsonSVD_SS[kr]->Print("base");
      if(debugMode)std::cout<<std::endl;
      //sameside covariance/pearson matrices???
      /////////////////////////////////////////////
      
      // 3x3 spectra canvases  ------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cSpectra canvas..."<<std::endl<<std::endl;
      cSpectra_ss->cd(kr+1);
      cSpectra_ss->cd(kr+1)->SetLogy();
      cSpectra_ss->cd(kr+1)->SetLogx();            
      
      //hrec_sameside_anabin_clone->SetTitle( ("S.S. Jet spectra,"+kRegRandEtaRange_plotTitle).c_str() );

      hunf_ss_svd[kr]->SetTitle( ("S.S. Spectra,"+kRegRandEtaRange_plotTitle).c_str() );
      hunf_ss_svd[kr]->SetName( ("hunf_ss_svd"+kRegRandEtaRange).c_str() );
      hunf_ss_svd[kr]->SetLineStyle(33);
      hunf_ss_svd[kr]->SetLineColor(kBlue);
      hunf_ss_svd[kr]->Print("base");
      
      hfold_ss_svd[kr]->SetName( ("hfold_ss_svd"+kRegRandEtaRange).c_str() );
      hfold_ss_svd[kr]->SetLineColor(kRed);      
      hfold_ss_svd[kr]->Print("base");

      
      hunf_ss_svd[kr]->Draw();
      hrec_sameside_anabin_clone->Draw("same");
      hfold_ss_svd[kr]->Draw("same");
      
      
      leg_ss[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"NBNDC");//x1,y1,x2,y2,header,option
      leg_ss[kr]->AddEntry(hrec_sameside_anabin_clone,"S.S. Measured","pl");
      leg_ss[kr]->AddEntry(hunf_ss_svd[kr],"S.S. Unfolded","pl");
      leg_ss[kr]->AddEntry(hfold_ss_svd[kr],"S.S. Folded","pl");
      leg_ss[kr]->SetTextSize(0.04);
      leg_ss[kr]->Draw();


      // 3x3 unfORfold/meas canvases  ---------------------
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_rec canvas..."<<std::endl;
      cRatio_rec_ss->cd(kr+1);
      cRatio_rec_ss->cd(kr+1)->SetLogx();
      


      hrec_ss_folded_ratio[kr] = (TH1F*)hfold_ss_svd[kr]->Clone( ("hrec_ss_fold_ratio"+kRegRandEtaRange).c_str() );
      hrec_ss_folded_ratio[kr]->SetMarkerStyle(27);
      hrec_ss_folded_ratio[kr]->SetMarkerColor(kRed);
      hrec_ss_folded_ratio[kr]->Divide(hrec_anabin);
      hrec_ss_folded_ratio[kr]->Print("base");


      hrec_ss_folded_ratio[kr]->SetTitle( ("S.S. SVD/Meas.,"+kRegRandEtaRange_plotTitle).c_str() );
      hrec_ss_folded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hrec_ss_folded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hrec_ss_folded_ratio[kr]->SetAxisRange(0.5, 1.5, "Y");

      hrec_ss_unfolded_ratio[kr] = (TH1F*)hunf_ss_svd[kr]->Clone( ("hrec_ss_unf_ratio"+kRegRandEtaRange).c_str());
      hrec_ss_unfolded_ratio[kr]->Divide(hrec_anabin);
      hrec_ss_unfolded_ratio[kr]->SetMarkerStyle(27);
      hrec_ss_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
      hrec_ss_unfolded_ratio[kr]->Print("base");


      hrec_ss_folded_ratio[kr]->Draw();
      hrec_ss_unfolded_ratio[kr]->Draw("same");

      leg_ss1[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg_ss1[kr]->AddEntry(hrec_ss_unfolded_ratio[kr],"S.S. Unf./Meas.","pl");
      leg_ss1[kr]->AddEntry(hrec_ss_folded_ratio[kr],"S.S. Folded/Meas.","pl");
      leg_ss1[kr]->SetTextSize(0.04); 
      leg_ss1[kr]->Draw();
      
      theLine->Draw();



      // draw on cRatio_gen canvas
      if(debugMode)std::cout<<std::endl<<"drawing stuff on cRatio_gen_ss canvas..."<<std::endl;
      cRatio_gen_ss->cd(kr+1);
      cRatio_gen_ss->cd(kr+1)->SetLogx();

      hgen_ss_unfolded_ratio[kr] = (TH1F*)hunf_ss_svd[kr]->Clone( ("hgen_ss_unf_ratio"+kRegRandEtaRange).c_str());
      hgen_ss_unfolded_ratio[kr]->Divide(hgen_anabin);      
      hgen_ss_unfolded_ratio[kr]->SetMarkerStyle(33);
      hgen_ss_unfolded_ratio[kr]->SetMarkerColor(kRed);
      hgen_ss_unfolded_ratio[kr]->Print("base");

      hgen_ss_unfolded_ratio[kr]->SetTitle( ("S.S. SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      //hgen_unfolded_ratio[kr]->SetMarkerStyle(27);
      //hgen_unfolded_ratio[kr]->SetMarkerColor(kBlue);      
      

      hgen_ss_folded_ratio[kr] = (TH1F*)hfold_ss_svd[kr]->Clone( ("hgen_ss_fold_ratio"+kRegRandEtaRange).c_str());
      hgen_ss_folded_ratio[kr]->Divide(hgen_anabin);      
      hgen_ss_folded_ratio[kr]->SetMarkerStyle(24);
      hgen_ss_folded_ratio[kr]->SetMarkerColor(kBlue);
      hgen_ss_folded_ratio[kr]->Print("base");
      
      hgen_ss_folded_ratio[kr]->SetTitle( ("S.S. SVD/Gen.,"+kRegRandEtaRange_plotTitle).c_str() );
      hgen_ss_folded_ratio[kr]->SetXTitle("Jet p_{T} (GeV/c)");
      hgen_ss_folded_ratio[kr]->SetAxisRange(boundaries_pt_gen[0], boundaries_pt_gen[nbins_pt_gen], "X");
      hgen_ss_folded_ratio[kr]->SetAxisRange(0.5, 1.5, "Y");
      
      hgen_ss_folded_ratio[kr]->Draw();
      hgen_ss_unfolded_ratio[kr]->Draw("same");
      
      leg_ss2[kr] = new TLegend(0.57, 0.75, 0.9, 0.9, NULL,"BRNDC");
      leg_ss2[kr]->AddEntry(hgen_ss_unfolded_ratio[kr],"S.S. Unf./Gen","pl");
      leg_ss2[kr]->AddEntry(hgen_ss_folded_ratio[kr],"S.S. Folded/Gen","pl");
      leg_ss2[kr]->SetTextSize(0.04); 
      leg_ss2[kr]->Draw();
      
      theLine->Draw();

      

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
	cPearsonMatrixIter->cd() ;   cPearsonMatrixIter->Print(outPdfFile.c_str());
	cSpectra->cd()           ;   cSpectra->Print(outPdfFile.c_str());
	
	//  cRatio 3X3 and Check for meas --------------------------------------------------
	cRatio_rec->cd()         ;   cRatio_rec->Print(outPdfFile.c_str());
	
	TCanvas *cRatioCheck_rec        = new TCanvas("cRatioCheck_rec","",        1500, 1500);    
	cRatioCheck_rec->cd();
	cRatioCheck_rec->SetLogx(1);
	
	hrec_folded_ratio[kRegDraw]->SetAxisRange(0.5, 1.5, "Y");
	hrec_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hrec_folded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ Meas.");
	
	hrec_folded_ratio[kRegDraw]->Draw();
	hrec_unfolded_ratio[kRegDraw]->Draw("same");
	
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
	
	TLine* theLine= new TLine(boundaries_pt_gen[0],1.,boundaries_pt_gen[nbins_pt_gen],1.);
	theLine->SetLineWidth(1);
	theLine->SetLineStyle(2);
	theLine->SetLineColor(36);
	theLine->Draw();
	
	cRatioCheck_rec->Print(outPdfFile.c_str());
	
	
	
	
	//  cRatio 3X3 Check for gen --------------------------------------------------
	cRatio_gen->cd()         ;   cRatio_gen->Print(outPdfFile.c_str());      
	
	TCanvas *cRatioCheck_gen        = new TCanvas("cRatioCheck_gen","",        1500, 1500);     
	cRatioCheck_gen->cd();
	cRatioCheck_gen->SetLogx(1);
	
	hgen_folded_ratio[kRegDraw]->SetAxisRange(0.5, 1.5, "Y");
	hgen_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hgen_folded_ratio[kRegDraw]->SetTitle("SVD, Ratios w/ Gen Truth");
	
	hgen_folded_ratio[kRegDraw]->Draw();
	hgen_unfolded_ratio[kRegDraw]->Draw("same");
	
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
	theLine->Draw();
	
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
	cRatioCheck_rec_ss->SetLogx(1);
	
	hrec_ss_folded_ratio[kRegDraw]->SetAxisRange(0.5, 1.5, "Y");
	hrec_ss_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hrec_ss_folded_ratio[kRegDraw]->SetTitle("S.S. SVD Unf. Ratios w/ Meas.");
	
	hrec_ss_folded_ratio[kRegDraw]->Draw();
	hrec_ss_unfolded_ratio[kRegDraw]->Draw("same");
	
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
	theLine->Draw();
	
	cRatioCheck_rec_ss->Print(outPdfFile.c_str());
	
	
	//  cRatio Check for gen --------------------------------------------------
	cRatio_gen_ss->cd()         ;   cRatio_gen_ss->Print(outPdfFile.c_str());
	
	TCanvas *cRatioCheck_gen_ss        = new TCanvas("cRatioCheck_gen_ss","",        1500, 1500);    
	cRatioCheck_gen_ss->cd();
	cRatioCheck_gen_ss->SetLogx(1);
	
	hgen_ss_folded_ratio[kRegDraw]->SetAxisRange(0.5, 1.5, "Y");
	hgen_ss_folded_ratio[kRegDraw]->SetAxisRange(boundaries_pt_reco[0], boundaries_pt_gen[nbins_pt_gen], "X");
	hgen_ss_folded_ratio[kRegDraw]->SetTitle("S.S. SVD Unf. Ratios w/ Gen Truth");
	
	hgen_ss_folded_ratio[kRegDraw]->Draw();
	hgen_ss_unfolded_ratio[kRegDraw]->Draw("same");
	
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
	theLine->Draw();
	
	cRatioCheck_gen_ss->Print(outPdfFile.c_str());
	
	cRatioCheck_gen_ss->Print(close_outPdfFile.c_str());
      }
      std::cout<<std::endl<<"done drawing SVD PDFs!"<<std::endl<<std::endl;
      
    }// end drawPDFs
  }// end SVD specific
  
  if(debugMode) std::cout<<"writing output root file..."<<std::endl;
  fout->Write();
  
  std::cout<<"closing input/output root files..."<<std::endl;
  fout->Close();
  fpp_MC->Close();
  return 0;

} // end unfoldDataSpectra



//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  int rStatus = -1;
  
  if( argc!=8 ){
    std::cout<<"do ./unfoldMCSpectra.exe <targMCDir> <baseOutputName> <doJetID> <kRegCenter> <doBayes> <doSVD> <useSimpleBins>"<<std::endl;
    std::cout<<"actually... just open the damn code and look"<<std::endl;

    return rStatus;  }
  
  rStatus=1; // runtime error
  
  rStatus=unfoldMCSpectra(  (const std::string)argv[1], (const std::string)argv[2], 
			    (int)std::atoi(argv[3]) , (int)std::atoi(argv[4]) , 
			    (int)std::atoi(argv[5]) , (int)std::atoi(argv[6]) , (int)std::atoi(argv[7]) 
			    
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
  
  
 
 
