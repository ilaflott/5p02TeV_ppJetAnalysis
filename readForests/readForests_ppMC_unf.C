// custom header
#include "readForests_unf.h"

// ppMC switches
const bool fillMCEvtQAHists=true;
const bool fillMCUnfoldingHists=true;
const bool fillMCUnfJetSpectraRapHists=true&&fillMCUnfoldingHists;
const bool fillMCJetCovMatrix=true; // leave me on almost always
const bool fillMCJetIDHists=true;//, tightJetID=false;
const int jetIDint=(int)fillMCJetIDHists;
const bool useJERScaleFactors=true;

const bool usegenetabins=true;//

const bool usePFCandAnalyzer=true;
const bool useTightJetID=false;
const bool verbose=false;

//// readForests_ppMC_unf
// ---------------------------------------------------------------------------------------------------------------
int readForests_ppMC_unf(std::string inFilelist , int startfile , int endfile ,  
			 int radius , std::string jetType , bool debugMode ,
			 std::string outfile  ,
			 float jtEtaCutLo, float jtEtaCutHi ){ 
  
  // for monitoring performance + debugging
  TStopwatch timer;  timer.Start();
  const bool debugMode_verbose=debugMode&&verbose; 
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;
  if(debugMode_verbose)std::cout<<std::endl<<"WARNING debugMode_verbose is ON, LARGE AMOUNTS OF TEXT INCOMING"<<std::endl;
  
  // basic info the screen
  std::cout<<std::endl<<"///////////////////"<<std::endl;
  std::cout<<"reading filelist "<< inFilelist<<std::endl;
  std::cout<<"reading files #'s "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"radius = " << radius;
  std::cout<<", jetType = " << jetType<<std::endl;

  
  std::cout<<"///////////////////"<<std::endl<<std::endl;
  
  // plot settings for later
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  
  // form jetTreeName
  std::string jetTreeName="ak"+std::to_string(radius)+jetType+"JetAnalyzer/t";
  if(debugMode)std::cout<<"looking at jetTree "<<jetTreeName<<std::endl;
  
  // initialize tree name array
  std::string trees[N_MCTrees];
  trees[0]=jetTreeName;
  for(int i=1;i<N_MCTrees;++i)trees[i]=MCTreeNames[i];
  
  // declare TChains for each tree + friend  
  TChain* jetpp[N_MCTrees]={};
  for(int t = 0;t<N_MCTrees;++t)  {     
    if(!usePFCandAnalyzer && t==4)continue;
    //if(!useTupel && t==5)continue;
    jetpp[t] = new TChain( trees[t].data() );
    if(t>0)jetpp[0]->AddFriend( jetpp[t] );  
  }

  
  // open filelist
  std::cout<<"opening filelist: "<<inFilelist<<std::endl;
  std::ifstream instr_Forest(inFilelist.c_str(),std::ifstream::in);
  std::string filename_Forest;  
  std::string lastFileAdded=""; bool filesAdded=false; //Int_t filesCount=0;
  for(int ifile = 0; ifile<=endfile; ++ifile){
    
    // grab a filename, check startfile and end of filelist condition
    instr_Forest>>filename_Forest; 
    if(ifile<startfile){ lastFileAdded=filename_Forest; 
      continue; }
    if(filename_Forest==lastFileAdded){ std::cout<<"end of filelist!"<<std::endl; 
      break; }
    
    std::cout<<"adding file #"<<ifile; 
    if(debugMode)std::cout<<", "<<filename_Forest; 
    
    std::cout<<", to each TChain in array"<<std::endl;
    for(int t = 0;t<N_MCTrees;++t) {
      if(!usePFCandAnalyzer && t==4)continue;
      //if(!useTupel &&t==5)continue;
      filesAdded=jetpp[t]->AddFile(filename_Forest.c_str()); 
    }
    lastFileAdded=filename_Forest;
  }//end file loop
  
  if(debugMode)std::cout<<"filesAdded="<<filesAdded<<std::endl;
  assert(filesAdded);//avoid segfault later

  // EXTRA JER SCALE FACTOR SMEARING FITS
  //[0]+[1]/(pow(x,[2])+[3]*x)
  TFile* JERsmearingfiles[4]={};
  TF1* JER_sigmu[4]={};
  double ptmin[4]={0.};
  double ptmax[4]={0.};
  if(useJERScaleFactors){
    for(int i=0; i<4; i++){
      JERsmearingfiles[i]=TFile::Open((JERsmearingfilenames[i]).c_str(),"READ");      
      JER_sigmu[i]=((TF1*)
		    ((TF1*)JERsmearingfiles[i]->Get("SigmaFit_f"))->Clone(
									  ("SigmaFit_f_"+absetabins_str[i]+"y"+absetabins_str[i+1]).c_str()
									  )
		    )
	;
      JERsmearingfiles[i]->Close();
      JER_sigmu[i]->GetRange(ptmin[i],ptmax[i]);
      //if(debugMode)JER_sigmu[i]->Print("base");
      if(debugMode)std::cout<<" this smearing fit is defined between "<<(int)ptmin[i]<<" GeV - "<< (int)ptmax[i]<<" GeV"<<std::endl;
    }
  }
  
  
  // Declare the output file, declare hists after for easy write
  std::cout<<"opening output file "<<outfile<<std::endl;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  fout->cd();  
  
  // declare hists
  /////   EVT+JET COUNTS AND METADATA   ///// 
  TH1F *hJetPtCut        =new TH1F("hJetPtCut"     ,(std::to_string(jtPtCut)).c_str()   ,    100, 0.,100.); hJetPtCut->Fill(jtPtCut);           
  TH1F *hJetPtCut_Hi        =new TH1F("hJetPtCut_Hi"     ,(std::to_string(jtPtCut_Hi)).c_str()   ,    1000, 1500.,2500.); hJetPtCut_Hi->Fill(jtPtCut_Hi);           
  TH1F *hGenJetPtCut     =new TH1F("hGenJetPtCut"  ,(std::to_string(genJetPtCut)).c_str()   ,    100, 0.,100.); hGenJetPtCut->Fill(genJetPtCut);           
  TH1F *hGenJetPtCut_Hi     =new TH1F("hGenJetPtCut_Hi"  ,(std::to_string(genJetPtCut_Hi)).c_str()   ,    1000, 1500.,2500.); hGenJetPtCut_Hi->Fill(genJetPtCut_Hi);           

  TH1F *hGenDRCut     =new TH1F("hGenDRCut"  ,(std::to_string(gendrCut)).c_str()   ,    100, 0.,100.); hGenDRCut->Fill(gendrCut);           
  
  TH1F *hJetEtaCutHi       =new TH1F("hJetEtaCutHi"     ,(std::to_string(jtEtaCutHi)).c_str()  ,    60,   0,6  ); hJetEtaCutHi->Fill(jtEtaCutHi);
  TH1F *hJetEtaCutLo       =new TH1F("hJetEtaCutLo"     ,(std::to_string(jtEtaCutLo)).c_str()  ,    60,   0,6  ); hJetEtaCutLo->Fill(jtEtaCutLo);	
  
  TH1D *h_NEvents         = new TH1D("NEvents","NEvents", 1,0,2);
  TH1D *h_NEvents_read    = new TH1D("NEvents_read","NEvents read", 1,0,2);
  TH1D *h_NEvents_skimCut = new TH1D("NEvents_skimCut"      ,"NEvents read post skimCut", 1,0,2);
  TH1D *h_NEvents_vzCut   = new TH1D("NEvents_vzCut"        ,"NEvents read post vzCut AND skimCut", 1,0,2);
  
  //counts post all evt cuts, that are thrown out because of no good jets pre and post jet cuts
  TH1D *h_NEvents_withJets           = new TH1D("NEvents_withJets" , 
		  				"NEvents read post evt cuts, w/ jets", 1,0,2);
  TH1D *h_NEvents_withJets_kmatCut1   = new TH1D("NEvents_withJets_kmatCut1" , 
						 "NEvents read post evt cuts, w/ jets post kmatCut1", 1,0,2);
  TH1D *h_NEvents_withJets_kmatCut2   = new TH1D("NEvents_withJets_kmatCut2" , 
						"NEvents read post evt cuts, w/ jets post kmatCut2", 1,0,2);
  
  //jet counts
  TH1D *h_NJets          = new TH1D("NJets","NJets read", 1,0,2);
  TH1D *h_NJets_kmatCut1  = new TH1D("NJets_kmatCut1 ","NJets read post kmatCut1 ", 1,0,2);
  TH1D *h_NJets_kmatCut2  = new TH1D("NJets_kmatCut2 ","NJets read post kmatCut2 ", 1,0,2);
  
  
  /////   EVENT QA   ///// 
  TH1D *hVr=NULL, *hWVr=NULL;
  TH1D *hVx=NULL, *hWVx=NULL;
  TH1D *hVy=NULL, *hWVy=NULL;
  TH1D *hVz=NULL, *hWVz=NULL, *hpthatWVz=NULL, *hvzWVz=NULL; 

  TH1D *hpthat=NULL, *hWpthat=NULL;  
  if(fillMCEvtQAHists){
    hpthatWVz = new TH1D("hpthatWeightedVz","", 500,-25.,25.);//pthat-weighted evtvz
    hvzWVz    = new TH1D("hvzWeightedVz","", 500,-25.,25.);//vz-weighted evtvz    
    hVr = new TH1D("hVr","vr, no trig, no weights",         2000,0.,0.60); 
    hWVr = new TH1D("hWeightedVr","vr, trigd, with weights",2000,0.,0.60);     
    hVx = new TH1D("hVx","vx, no trig, no weights",  2000, -0.30,0.30);
    hWVx = new TH1D("hWeightedVx","vx, trigd, with weights", 2000, -0.30,0.30); 
    hVy = new TH1D("hVy","vy, no trig, no weights",   2000, -0.30,0.30); 
    hWVy = new TH1D("hWeightedVy","vy, trigd, with weights",   2000,-0.30,0.30);    
    hVz       = new TH1D("hVz","", 1000,-25.,25.);//evtvz
    hWVz      = new TH1D("hWeightedVz","", 1000,-25.,25.);//pthat*vz-weighted evt vz
    hpthat    = new TH1D("hpthat","",1000,0,1000);//evt pthat, unweighted and weighted
    hWpthat   = new TH1D("hWeightedpthat","",1000,0,1000);    }



  /////  UNFOLDING   /////
  //to unfold ppData 
  TH1D *hpp_gen=NULL;    
  TH1D *hpp_reco=NULL; 
  TH2D *hpp_matrix=NULL; 
  
  //to test MC sample consistency in unfolding
  TH1D *hpp_mcclosure_gen=NULL;      //the first three are for the "truth" response matrix
  TH1D *hpp_mcclosure_reco=NULL;         
  TH2D *hpp_mcclosure_matrix=NULL;         
  TH1D *hpp_mcclosure_reco_test=NULL;   
  
  //to unfold ppData 
  TH1D *hpp_gen_eta[nbins_abseta]={};    
  TH1D *hpp_reco_eta[nbins_abseta]={}; 
  TH2D *hpp_matrix_eta[nbins_abseta]={};
 
  TH1D *hpp_gen_y[nbins_abseta]={};    
  TH1D *hpp_reco_y[nbins_abseta]={}; 
  TH2D *hpp_matrix_y[nbins_abseta]={}; 

  TH1D *hpp_reco_y_JERsysup[nbins_abseta]={}; 
  TH2D *hpp_matrix_y_JERsysup[nbins_abseta]={}; 

  TH1D *hpp_reco_y_JERsysdown[nbins_abseta]={}; 
  TH2D *hpp_matrix_y_JERsysdown[nbins_abseta]={}; 
  
  //TH1D *hpp_gen_geneta[nbins_abseta]={};    
  //TH1D *hpp_reco_geneta[nbins_abseta]={}; 
  //TH2D *hpp_matrix_geneta[nbins_abseta]={}; 
  
  //to test MC sample consistency in unfolding
  TH1D *hpp_mcclosure_gen_eta[nbins_abseta]={};      //the first three are for the "truth" response matrix
  TH1D *hpp_mcclosure_reco_eta[nbins_abseta]={};         
  TH2D *hpp_mcclosure_matrix_eta[nbins_abseta]={};         
  TH1D *hpp_mcclosure_reco_test_eta[nbins_abseta]={};   

  TH1D *hpp_mcclosure_gen_y[nbins_abseta]={};      //the first three are for the "truth" response matrix
  TH1D *hpp_mcclosure_reco_y[nbins_abseta]={};         
  TH2D *hpp_mcclosure_matrix_y[nbins_abseta]={};         
  TH1D *hpp_mcclosure_reco_test_y[nbins_abseta]={};   
  
  if(fillMCUnfoldingHists){
    
    std::string hUnfTitleArray[]={ "gen", "reco", "matrix", 
				   "mcclosure_gen" , "mcclosure_reco" , "mcclosure_matrix" , 
				   "mcclosure_gen_test" , "mcclosure_reco_test" , "mcclosure_matrix_test"     }; 
    const int nUnfTitles=sizeof(hUnfTitleArray)/sizeof(std::string);    
    
    //    for(int jtID=0;jtID<2;jtID++){
    
    //if(!fillMCJetIDHists&&jtID==1)continue;
    //if(fillMCJetIDHists&&jtID==0)continue;
    
    for(int k=0; k<nUnfTitles; k++){     
      
      std::string hTitle="hpp_"+hUnfTitleArray[k];
      if(fillMCJetIDHists)hTitle+="_wJetID";
      hTitle+="_R"+std::to_string(radius)+"_"+etaWidth;      
      
      if(hUnfTitleArray[k]=="gen")	{
	hpp_gen    = new TH1D( hTitle.c_str(), "MC genpt for unf data", 2500,0,2500);
	if(fillMCUnfJetSpectraRapHists){
	  
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_gen_eta[j] = new TH1D( (hTitle+"_bin"+std::to_string(j)).c_str(), 
				       ("MC gen pt for unf data bin"+std::to_string(j)).c_str(), 2500,0,2500);	  
	  
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_gen_y[j] = new TH1D( (hTitle+"_ybin"+std::to_string(j)).c_str(), 
				     ("MC gen pt for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500);	  
	  //if(useJERScaleFactors)
	  //  for(int j=0;j<nbins_abseta;j++)
	  //    hpp_gen_y_JERsysup[j] = new TH1D( (hTitle+"_ybin"+std::to_string(j)+"_JERsysup").c_str(), 
	  //			       ("MC gen pt for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500);	  
	  //if(useJERScaleFactors)
	  //  for(int j=0;j<nbins_abseta;j++)
	  //    hpp_gen_y_JERsysdown[j] = new TH1D( (hTitle+"_ybin"+std::to_string(j)+"_JERsysdown").c_str(), 
	  //			       ("MC gen pt for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500);	  
	}
      }
      else if(hUnfTitleArray[k]=="reco")	{	  
	hpp_reco    = new TH1D( hTitle.c_str(), 
				"MC recopt for unf data", 2500,0,2500);
	if(fillMCUnfJetSpectraRapHists){
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_reco_eta[j] = new TH1D( (hTitle+"_bin"+std::to_string(j)).c_str(), 
					("MC reco pt for unf data bin"+std::to_string(j)).c_str(), 2500,0,2500); 	  

	  for(int j=0;j<nbins_abseta;j++)
	    hpp_reco_y[j] = new TH1D( (hTitle+"_ybin"+std::to_string(j)).c_str(), 
				      ("MC reco pt for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500); 	  
	  
	  if(useJERScaleFactors){
	    for(int j=0;j<nbins_abseta;j++)
	      hpp_reco_y_JERsysup[j] = new TH1D( (hTitle+"_JERsysup_ybin"+std::to_string(j)).c_str(), 
						 ("MC reco pt JERsysup for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500); 	  
	    for(int j=0;j<nbins_abseta;j++)
	      hpp_reco_y_JERsysdown[j] = new TH1D( (hTitle+"_JERsysdown_ybin"+std::to_string(j)).c_str(), 
						   ("MC reco pt for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500); 	  
	  }
	}
      }
      else if(hUnfTitleArray[k]=="matrix")	{
	hpp_matrix    = new TH2D( hTitle.c_str(), 
				  "MC gentpt v. recopt for unf data", 2500, 0,2500, 2500, 0,2500);
	if(fillMCUnfJetSpectraRapHists){
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_matrix_eta[j] = new TH2D( (hTitle+"_bin"+std::to_string(j)).c_str(), 
					  ("MC genpt v. recopt for unf data bin"+std::to_string(j)).c_str(), 2500,0,2500,2500,0,2500);  
	  
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_matrix_y[j] = new TH2D( (hTitle+"_ybin"+std::to_string(j)).c_str(), 
					("MC genpt v. recopt for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500,2500,0,2500);  
	  if(useJERScaleFactors){
	    for(int j=0;j<nbins_abseta;j++)
	      hpp_matrix_y_JERsysup[j] = new TH2D( (hTitle+"_JERsysup_ybin"+std::to_string(j)).c_str(), 
						   ("MC genpt v. recopt JERsysup for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500,2500,0,2500);  
	    for(int j=0;j<nbins_abseta;j++)
	      hpp_matrix_y_JERsysdown[j] = new TH2D( (hTitle+"_JERsysdown_ybin"+std::to_string(j)).c_str(), 
						     ("MC genpt v. recopt JERsysdown for unf data ybin"+std::to_string(j)).c_str(), 2500,0,2500,2500,0,2500);  
	  }
	}
      }
      else if(hUnfTitleArray[k]=="mcclosure_gen") {
	hpp_mcclosure_gen    = new TH1D( hTitle.c_str(), 
					 "genpt for mcclosure same side", 2500,0,2500);
	if(fillMCUnfJetSpectraRapHists){
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_gen_eta[j] = new TH1D( (hTitle+"_bin"+std::to_string(j)).c_str(), 
						 ("MC gen pt for mcclosure same side bin"+std::to_string(j)).c_str(), 2500,0,2500);
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_gen_y[j] = new TH1D( (hTitle+"_ybin"+std::to_string(j)).c_str(), 
					       ("MC gen pt for mcclosure same side ybin"+std::to_string(j)).c_str(), 2500,0,2500);
	}
      }
      else if(hUnfTitleArray[k]=="mcclosure_reco")	{
	hpp_mcclosure_reco    = new TH1D( hTitle.c_str(), 
					  "recopt for mcclosure same side", 2500,0,2500);
	if(fillMCUnfJetSpectraRapHists){
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_reco_eta[j] = new TH1D( (hTitle+"_bin"+std::to_string(j)).c_str(), 
						  ("MC reco pt for mcclosure same side bin"+std::to_string(j)).c_str(), 2500,0,2500);   
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_reco_y[j] = new TH1D( (hTitle+"_ybin"+std::to_string(j)).c_str(), 
						("MC reco pt for mcclosure same side ybin"+std::to_string(j)).c_str(), 2500,0,2500);   
	}
      }
      else if(hUnfTitleArray[k]=="mcclosure_matrix")	{
	hpp_mcclosure_matrix    = new TH2D( hTitle.c_str(), 
					    "genpt v. recopt for mcclosure same side", 2500, 0,2500, 2500, 0,2500);
	if(fillMCUnfJetSpectraRapHists){
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_matrix_eta[j] = new TH2D( (hTitle+"_bin"+std::to_string(j)).c_str(), 
						    ("MC genpt v. recopt for mcclosure same side bin"+std::to_string(j)).c_str(), 2500,0,2500,2500,0,2500);  
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_matrix_y[j] = new TH2D( (hTitle+"_ybin"+std::to_string(j)).c_str(), 
						  ("MC genpt v. recopt for mcclosure same side ybin"+std::to_string(j)).c_str(), 2500,0,2500,2500,0,2500);  
	}
      }
      else if(hUnfTitleArray[k]=="mcclosure_reco_test")	{
	hpp_mcclosure_reco_test    = new TH1D( hTitle.c_str(), 
					       "recopt for mcclosure opp side test", 2500,0,2500);
	if(fillMCUnfJetSpectraRapHists){
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_reco_test_eta[j] = new TH1D( (hTitle+"_bin"+std::to_string(j)).c_str(), 
						       ("MC reco pt for mcclosure opp side test bin"+std::to_string(j)).c_str(), 2500,0,2500);   
	  for(int j=0;j<nbins_abseta;j++)
	    hpp_mcclosure_reco_test_y[j] = new TH1D( (hTitle+"_ybin"+std::to_string(j)).c_str(), 
						       ("MC reco pt for mcclosure opp side test ybin"+std::to_string(j)).c_str(), 2500,0,2500);   
	}
      }
      else continue;
    }//end loop over unf titles
  }//end fill mc unfolding hists
  
  

  TH2D* hpp_mcclosure_covmat_test=NULL;
  TH1D* hpp_mcclosure_covmat_test_helper=NULL;
  bool bool_helper=false;
  
  TH2D* hpp_mcclosure_covmat_test_eta_arr[nbins_abseta]={};
  TH1D* hpp_mcclosure_covmat_test_eta_arr_helpers[nbins_abseta]={};
  bool etabin_bool_helper[nbins_abseta]={0};

  TH2D* hpp_mcclosure_covmat_test_y_arr[nbins_abseta]={};
  TH1D* hpp_mcclosure_covmat_test_y_arr_helpers[nbins_abseta]={};
  bool ybin_bool_helper[nbins_abseta]={0};
  
  if(fillMCJetCovMatrix && fillMCJetIDHists){
    std::string covmat_title="MCClosure, Test PY8 Jet Covariance Matrix, "+absetabins_str[0]+" < #||{#eta} < "+absetabins_str[nbins_abseta];    
    if(debugMode)std::cout<<"covmat_title =" << covmat_title<< std::endl;      
    std::string covmat_name="hpp_mcclosure_covmat_test_";
    covmat_name+="wJetID";
    if(debugMode)std::cout<<"covmat_name  =" << covmat_name << std::endl;      
    hpp_mcclosure_covmat_test=new TH2D(covmat_name.c_str(), covmat_title.c_str(), nbins_pt, ptbins, nbins_pt, ptbins);     
    hpp_mcclosure_covmat_test_helper=new TH1D( "mcclosure_covmat_test_helper", "helper only", nbins_pt, ptbins);      
    
    if(fillMCUnfJetSpectraRapHists){
      for(int i=0;i<nbins_abseta;i++){
	std::string covmat_eta_title="MCClosure, Test PY8 Jet Covariance Matrix, "+absetabins_str[i]+" < #||{#eta} < "+absetabins_str[i+1];    
	if(debugMode)std::cout<<"covmat_eta_title =" << covmat_eta_title<< std::endl;      
	std::string covmat_eta_name="hpp_mcclosure_covmat_test_";
	//if(fillDataJetIDHists)covmat_eta_name+="wJetID_";
	covmat_eta_name+="wJetID_";
	covmat_eta_name+="etabin_"+std::to_string(i);      
	if(debugMode)std::cout<<"covmat_eta_name  =" << covmat_eta_name << std::endl;      
	hpp_mcclosure_covmat_test_eta_arr[i]=new TH2D(covmat_eta_name.c_str(), covmat_eta_title.c_str(), nbins_pt, ptbins, nbins_pt, ptbins);     
	hpp_mcclosure_covmat_test_eta_arr_helpers[i]=new TH1D(("mcclosure_covmat_test_helper_etabin"+std::to_string(i)).c_str(), "helper only", nbins_pt, ptbins);      
      }//end etabin loop
      
      for(int i=0;i<nbins_abseta;i++){
	std::string covmat_y_title="MCClosure, Test PY8 Jet Covariance Matrix, "+absetabins_str[i]+" < #||{y} < "+absetabins_str[i+1];    
	if(debugMode)std::cout<<"covmat_y_title =" << covmat_y_title<< std::endl;      
	std::string covmat_y_name="hpp_mcclosure_covmat_test_";
	//if(fillDataJetIDHists)covmat_eta_name+="wJetID_";
	covmat_y_name+="wJetID_";
	covmat_y_name+="ybin_"+std::to_string(i);      
	if(debugMode)std::cout<<"covmat_y_name  =" << covmat_y_name << std::endl;      
	hpp_mcclosure_covmat_test_y_arr[i]=new TH2D(covmat_y_name.c_str(), covmat_y_title.c_str(), nbins_pt, ptbins, nbins_pt, ptbins);     
	hpp_mcclosure_covmat_test_y_arr_helpers[i]=new TH1D(("mcclosure_covmat_test_helper_ybin"+std::to_string(i)).c_str(), "helper only", nbins_pt, ptbins);      
      }//end etabin loop
    }
  }//end fillMCJetCovMatrix


  // EVENT LOOP PREP
  // declare variables/arrays + set branch address for each input tree
  //JetAnalyzer, jets
  int nref_I;  
  float pt_F[1000];  float eta_F[1000];   float phi_F[1000];
  float rawpt_F[1000]; float m_F[1000];
  //tracks
  int trkN_I[1000], trkHardN_I[1000];
  float trkSum_F[1000],  trkMax_F[1000], trkHardSum_F[1000];
  //charged hadrons
  int chN_I[1000], chHardN_I[1000];
  float chSum_F[1000],  chMax_F[1000], chHardSum_F[1000];
  //photons
  int phN_I[1000], phHardN_I[1000];
  float phSum_F[1000],  phMax_F[1000], phHardSum_F[1000];
  //leptons
  int eN_I[1000], muN_I[1000];
  float eSum_F[1000];   float eMax_F[1000];
  float muSum_F[1000];  float muMax_F[1000];
  //neutrals
  int neN_I[1000];
  float neSum_F[1000];  float neMax_F[1000];
  //MC jet variable
  float pthat_F;
  int subid_I[1000];    int refparton_F[1000];
  float refpt_F[1000];  float refeta_F[1000];  float refphi_F[1000];  
  float refdrjt_F[1000];
  //MC gen jet variables
  int ngen_I;  int genmatchindex_I[1000];
  float genpt_F[1000];

  //jets
  jetpp[0]->SetBranchAddress("nref",&nref_I);
  jetpp[0]->SetBranchAddress("jtpt",&pt_F);
  jetpp[0]->SetBranchAddress("jteta",&eta_F);
  jetpp[0]->SetBranchAddress("jtphi",&phi_F);
  jetpp[0]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[0]->SetBranchAddress("jtm",&m_F);
  //tracks
  jetpp[0]->SetBranchAddress("trackN",&trkN_I);
  jetpp[0]->SetBranchAddress("trackSum",&trkSum_F);
  jetpp[0]->SetBranchAddress("trackMax",&trkMax_F);
  jetpp[0]->SetBranchAddress("trackHardN",&trkHardN_I);
  jetpp[0]->SetBranchAddress("trackHardSum",&trkHardSum_F);
  //charged hadrons from PF (no HF-Calo ctontributions)
  jetpp[0]->SetBranchAddress("chargedN",&chN_I);
  jetpp[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[0]->SetBranchAddress("chargedMax",&chMax_F);
  jetpp[0]->SetBranchAddress("chargedHardN",&chHardN_I);
  jetpp[0]->SetBranchAddress("chargedHardSum",&chHardSum_F);
  //photons from PF (no HF-Calo contributions)
  jetpp[0]->SetBranchAddress("photonN",&phN_I);
  jetpp[0]->SetBranchAddress("photonSum",&phSum_F);
  jetpp[0]->SetBranchAddress("photonMax",&phMax_F);
  jetpp[0]->SetBranchAddress("photonHardN",&phHardN_I);
  jetpp[0]->SetBranchAddress("photonHardSum",&phHardSum_F);
  //electrons from PF
  jetpp[0]->SetBranchAddress("eN",&eN_I);
  jetpp[0]->SetBranchAddress("eSum",&eSum_F);
  jetpp[0]->SetBranchAddress("eMax",&eMax_F);
  //muons from PF
  jetpp[0]->SetBranchAddress("muN",&muN_I);
  jetpp[0]->SetBranchAddress("muSum",&muSum_F);
  jetpp[0]->SetBranchAddress("muMax",&muMax_F);
  //neutral hadrons from PF (no HF-Calo contributions)
  jetpp[0]->SetBranchAddress("neutralN",&neN_I);
  jetpp[0]->SetBranchAddress("neutralSum",&neSum_F);
  jetpp[0]->SetBranchAddress("neutralMax",&neMax_F);
  //MC jet variables
  jetpp[0]->SetBranchAddress( "pthat"   , &pthat_F);
  jetpp[0]->SetBranchAddress( "subid"	, &subid_I     );
  jetpp[0]->SetBranchAddress( "refparton_flavor" , &refparton_F );
  jetpp[0]->SetBranchAddress( "refpt"	, &refpt_F     );
  jetpp[0]->SetBranchAddress( "refeta"  , &refeta_F    );
  jetpp[0]->SetBranchAddress( "refphi"  , &refphi_F    );
  jetpp[0]->SetBranchAddress( "refdrjt"	, &refdrjt_F   );  
  //gen jet stuff; for making "reco match index" array
  jetpp[0]->SetBranchAddress( "ngen" , &ngen_I );
  jetpp[0]->SetBranchAddress( "genmatchindex", &genmatchindex_I);
  jetpp[0]->SetBranchAddress( "genpt" , &genpt_F);
  
  // HiEvtAnalyzer
  ULong64_t evt_I;   UInt_t run_I;   UInt_t lumi_I; 
  float vz_F,vx_F,vy_F;
  jetpp[1]->SetBranchAddress("evt",&evt_I);
  jetpp[1]->SetBranchAddress("run",&run_I);
  jetpp[1]->SetBranchAddress("lumi",&lumi_I);
  jetpp[1]->SetBranchAddress("vz",&vz_F);
  jetpp[1]->SetBranchAddress("vx",&vx_F);
  jetpp[1]->SetBranchAddress("vy",&vy_F);
  

  // skimanalysis
  int pBeamScrapingFilter_I, pHBHENoiseFilterRun2Loose_I, pprimaryvertexFilter_I, puvertexFilter_I;
  int pHBHEIsoNoiseFilter_I;
  jetpp[2]->SetBranchAddress("HBHEIsoNoiseFilterResult",&pHBHEIsoNoiseFilter_I);//only option
  jetpp[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_I);
  jetpp[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilterRun2Loose_I);
  jetpp[2]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_I);
  jetpp[2]->SetBranchAddress("pVertexFilterCutGtight",&puvertexFilter_I);  


  //pfcandAnalyzer/pfTree
  Int_t nPFpart_I;  
  //std::vector<Float_t> *pfpt, *pfphi, *pfeta, *pfE;
  Float_t pfptArr[2000], pfphiArr[2000], pfetaArr[2000], pfEArr[2000];
  
  if(usePFCandAnalyzer){
    if(debugMode)std::cout<<"using PFCandAnalyzer. Setting branch addresses!"<<std::endl;
    jetpp[4]->SetBranchAddress( "nPFpart", &nPFpart_I);
    if(debugMode)std::cout<<"filelistIsLowerJets, pointing branches to Float_t* (a.k.a. Float_t array)"<<std::endl;
    jetpp[4]->SetBranchAddress( "pfPt"   , &pfptArr );
    jetpp[4]->SetBranchAddress( "pfEta"  , &pfetaArr);
    jetpp[4]->SetBranchAddress( "pfPhi"  , &pfphiArr);    
    jetpp[4]->SetBranchAddress( "pfEnergy"  , &pfEArr);    
  }



  //tupel/EventTree (for calculating genjet rapidity for AK4 only [R=0.3 jet info absent])
  std::vector<float> *jetPt, *jetEta, *jetPhi, *jetE;// *jetID, *jetRawPt, *jetRawE;
  std::vector<float> *GjetPt, *GjetEta, *GjetPhi, *GjetE;// *jetID, *jetRawPt, *jetRawE;
  //std::vector<float> *jetHfHadE, *jetHfEmE, *jetChHadFrac, *jetNeutralHadAndHfFrac, *jetChEmFrac, *jetNeutralEmFrac, *jetChMult, *jetConstCnt;  
  if(debugMode)std::cout<<"setting branch addresses for tupel/EventTree"<<std::endl;
  
  jetpp[5]->SetBranchAddress( "JetAk04Pt",                  &jetPt                  );
  jetpp[5]->SetBranchAddress( "JetAk04Eta",                 &jetEta                 );
  jetpp[5]->SetBranchAddress( "JetAk04Phi",                 &jetPhi                 );
  jetpp[5]->SetBranchAddress( "JetAk04E",                   &jetE                   );
  jetpp[5]->SetBranchAddress( "GJetAk04Pt",                 &GjetPt                  );
  jetpp[5]->SetBranchAddress( "GJetAk04Eta",                &GjetEta                 );
  jetpp[5]->SetBranchAddress( "GJetAk04Phi",                &GjetPhi                 );
  jetpp[5]->SetBranchAddress( "GJetAk04E",                  &GjetE                   );

  
  // event count from files
  UInt_t NEvents_jetAnalyzr=jetpp[0]->GetEntries();   // preskim event count from files
  UInt_t NEvents_skimAnalyzr=jetpp[2]->GetEntries();   // preskim event count from files
  bool jetEventsExist= (bool)NEvents_jetAnalyzr, skimEventsExist=(bool)NEvents_skimAnalyzr;
  bool EventCountsEqual= (NEvents_jetAnalyzr==NEvents_skimAnalyzr);
  if( !EventCountsEqual ||
      !jetEventsExist  || !skimEventsExist ){//if this fires, we're gonna segfault
    std::cout<<"jetAnalyzr # evts  = " << NEvents_jetAnalyzr  << std::endl;
    std::cout<<"skimAnalyzr # evts = " << NEvents_skimAnalyzr << std::endl; }
  assert( jetEventsExist && skimEventsExist && EventCountsEqual );//exit if the counts don't make sense, help avoid garbage results
  
  UInt_t NEvents_allFiles=NEvents_jetAnalyzr;   // preskim event count from files
  for(UInt_t i=0;i < NEvents_allFiles; ++i) h_NEvents->Fill(1);
  
  UInt_t NEvents_read=0;
  NEvents_read = NEvents_allFiles;
  //NEvents_read=0;//debug

  int goodgenjetcount=0,badgenjetcount=0;  

  std::cout<<"reading "<<NEvents_read<<" events"<<std::endl;   

  int mcclosureInt=0;
  bool incrementClosureInt=false;
  
  //int mcclosureInt_arr[nbins_abseta]={0,0,0,0};
  int mcclosureInt_arr[nbins_abseta]={0};
  bool incrementClosureInt_arr[nbins_abseta]={0};
  int mcclosureInt_y_arr[nbins_abseta]={0};
  bool incrementClosureInt_y_arr[nbins_abseta]={0};
  //bool eventAccepted=false;
  
  for(UInt_t nEvt = 0; nEvt < NEvents_read; ++nEvt) {//event loop   
    //eventAccepted=false;
    
    incrementClosureInt=false;
    for(int i=0; i<nbins_abseta; i++)
      incrementClosureInt_arr[i]=false;

    for(int i=0; i<nbins_abseta; i++)
      incrementClosureInt_y_arr[i]=false;
      
    
    if(nEvt%1000==0)std::cout<<"from trees, grabbing Evt # = "<<nEvt<<std::endl;  
    jetpp[0]->GetEntry(nEvt);    
    h_NEvents_read->Fill(1);
    
    // skim/HiEvtAnalysis criteria
    if( pHBHENoiseFilterRun2Loose_I==0     || 
	pHBHEIsoNoiseFilter_I    ==0    || 
        pBeamScrapingFilter_I==0  || 
        pprimaryvertexFilter_I==0  ) continue;    
    h_NEvents_skimCut->Fill(1);    

    if( fabs(vz_F)>24.              ) continue;
    h_NEvents_vzCut->Fill(1);





    if(usePFCandAnalyzer){
      //PFCANDLOOP
      //std::cout<<std::endl<<"nPFpart_I="<<nPFpart_I<<std::endl;
      
      Float_t sumpfpt_x=0.,sumpfpt_y=0.;
      Float_t sumpfEt=0.;//scalar sum of PF Et, not a vector sum.
      //h_nPFpart->Fill(nPFpart_I,weight_eS);
      
      for(Int_t part=0;part<nPFpart_I;part++){
	
	Float_t pfpt_part,pfeta_part,pfphi_part, pfE_part;
	pfpt_part  = pfptArr[part];	  
	pfeta_part =pfetaArr[part];	  
	pfphi_part =pfphiArr[part];	  	
	pfE_part  =pfEArr[part];	
	
	//h_pfpt ->Fill(pfpt_part ,weight_eS);
	//h_pfeta->Fill(pfeta_part,weight_eS);
	//h_pfphi->Fill(pfphi_part,weight_eS);
	//h_pfE->Fill(pfE_part,weight_eS);
	
	Float_t pftheta_part=2.*atan(exp(-1.*pfeta_part));
	//h_pftheta->Fill(pftheta_part,weight_eS);
	
	Float_t pfEt_part=pfE_part*TMath::Sin(pftheta_part);		
	//h_pfEt->Fill(pfEt_part,weight_eS);
	
	//std::cout<<"pfpt->at("<<part<<")=" << pfpt_part<<std::endl;
	//std::cout<<"pfphi->at("<<part<<")="<<pfphi_part<<std::endl;
	//std::cout<<"pfeta->at("<<part<<")="<<pfeta_part<<std::endl;
	sumpfpt_x+=pfpt_part*TMath::Cos(pfphi_part);
	sumpfpt_y+=pfpt_part*TMath::Sin(pfphi_part);
	sumpfEt+=pfEt_part;
      }
      
      Float_t sumpfpt=sqrt(sumpfpt_x*sumpfpt_x + sumpfpt_y*sumpfpt_y );
      //Float_t sqrtsumpfEt=sqrt(sumpfEt);
      Float_t PFMETfrac=sumpfpt/sumpfEt;
      //Float_t PFMETsig=sumpfpt/sqrtsumpfEt;
      
      //h_PFMETfrac->Fill(PFMETfrac,weight_eS);
      if(PFMETfrac>0.3)continue;//giving this a shot. see how it goes...
      //h_NEvents_PFMETfracCut->Fill(0.);
      
      //h_sumpfpt_x->Fill(sumpfpt_x,weight_eS);
      //h_sumpfpt_y->Fill(sumpfpt_y,weight_eS);
      //h_sumpfEt->Fill(sumpfEt,weight_eS);

      //h_sumpfpt->Fill(sumpfpt,weight_eS);
      //h_sqrtsumpfEt->Fill(sqrtsumpfEt,weight_eS);
      //h_PFMETfrac->Fill(PFMETfrac,weight_eS);
      //h_PFMETsig->Fill(PFMETsig,weight_eS);
      //assert(false);
    }

    
    //eventAccepted=true;
    //if(debugMode && eventAccepted)
    //  std::cout<<"nEvt="<<nEvt<<", event accepted"<<std::endl;
    
    // grab vzweight
    double vzWeight=1.;
    if(doVzWeights){
      //std::cout<<"vzWeight was:"<<vzWeight<<std::endl;
      vzWeight=cpuVzWeight_poly(vz_F);
    }
    
    double evtPthatWeight=1.;    
    for( int i=0; i<nbins_pthat && pthat_F>=pthatbins[i]; i++ ){ 
      evtPthatWeight=pthatWeights[i]; 
      evtPthatWeight*=1e+06;//in nanobarns
    }     
    
    double trigWeight=1.;
    //no trig selection, trig weight always 1

    double weight_eS=evtPthatWeight*trigWeight*vzWeight;              
    float vr_F=std::sqrt( vx_F*vx_F + vy_F*vy_F);
    //vz hists
    if(fillMCEvtQAHists){
      
      hVr->Fill( vr_F, 1.);
      hWVr->Fill(vr_F, weight_eS);      
      hVx->Fill( vx_F, 1.);
      hWVx->Fill(vx_F, weight_eS);      
      hVy->Fill( vy_F, 1.);
      hWVy->Fill(vy_F, weight_eS);      

      hVz->Fill(vz_F, 1.);
      hWVz->Fill(vz_F, weight_eS);      
      hvzWVz->Fill(vz_F, vzWeight);
      hpthatWVz->Fill(vz_F, evtPthatWeight);

      hpthat->Fill(pthat_F, 1.);
      hWpthat->Fill(pthat_F, weight_eS); }



    std::vector<int>recomatchindex_I;
    recomatchindex_I.assign(nref_I,-1);//set nref_I elements of vector to -1, indicating that reco jet N doesn't yet have matching genjet

    bool tempdebug=false;//debugMode&&(nEvt%200==0);    
    if(tempdebug)std::cout<<std::endl<<"_______ nEvt="<<nEvt<<" ________"<<std::endl;
    if(tempdebug)std::cout<<"ngen_I="<<ngen_I<<std::endl;
    for(int gjet=0; gjet<ngen_I; ++gjet){
      if(tempdebug)std::cout<<"// genmatchindex_I["<<gjet<<"]="<<genmatchindex_I[gjet]<<std::endl;
      if(tempdebug)std::cout<<"genpt_F["<<gjet<<"]="<<genpt_F[gjet]<<std::endl;
      if((genmatchindex_I[gjet]!=-1) &&
	 (genmatchindex_I[gjet]<nref_I)) recomatchindex_I[genmatchindex_I[gjet]]=gjet;
    }
    
    // for event counting + avoiding duplicate fills in dijet hists
    bool hNEvts_withJets_Filled=false;
    bool hNEvts_withJets_kmatCut1_Filled=false, hNEvts_withJets_kmatCut2_Filled=false; 
    //bool jetAccepted=false;
    for(int jet = 0; jet<nref_I; ++jet){
      //jetAccepted=false;
      // event+jet counting
      h_NJets->Fill(1);
      if(!hNEvts_withJets_Filled){
	h_NEvents_withJets->Fill(1);
	hNEvts_withJets_Filled=true;      }      

      float jtpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
      float recpt  = pt_F[jet]; //this will be the jetpt pulled from the forest
      float receta = eta_F[jet];
      float jtm=m_F[jet];
      
      int chMult  = chN_I[jet] + eN_I[jet] + muN_I[jet] ;
      int neuMult = neN_I[jet] + phN_I[jet] ;
      int numConst  = chMult + neuMult;
      
      float absreceta = fabs(receta);	
      float jttheta=2.*atan(exp(-1.*receta));
      float jtpz=jtpt/tan(jttheta);
      float jtp=sqrt(jtpt*jtpt + jtpz*jtpz);
      float jtE=sqrt(jtp*jtp + jtm*jtm);
      float jty=0.5*log((jtE+jtpz)/(jtE-jtpz));//experimentalist version
      float absjty=fabs(jty);
      
      float genpt   = refpt_F[jet];
      //float genphi  = refphi_F[jet];
      float geneta  = refeta_F[jet];
      float absgeneta  = fabs(geneta);
      float gendrjt = refdrjt_F[jet];

      
      float jetIDpt=recpt;//ala HIN jetID, recpt is corrected w/ L2/L3 residuals
      
      
      if( subid_I[jet]!=0 ) {
	//if(debugMode)std::cout<<"jet not matched!"<<std::endl;
	continue;//matching gen-reco
      }
      else if ( !(recpt > jtPtCut)    ) {
	//if(debugMode)std::cout<<"recpt="<<recpt<<", too low!"<<std::endl;
	continue;        //low recopt cut          
      }
      else if ( !(recpt < jtPtCut_Hi)    ){ 
	//if(debugMode)std::cout<<"recpt too high!"<<std::endl;
	continue;}     //high recopt cut
      else if ( !(genpt > genJetPtCut) ) {
	//if(debugMode)std::cout<<"genpt too low!"<<std::endl;
	continue;}       //low genpt cut
      else if ( !(genpt < genJetPtCut_Hi) ) {
	//if(debugMode)std::cout<<"genpt too high!"<<std::endl;
	continue;}   //high genpt cut
      else if ( absreceta < jtEtaCutLo ) {
	//if(debugMode)std::cout<<"absreceta too low!"<<std::endl;
	continue;} // lower abseta cut 
      else if (!(absreceta < jtEtaCutHi)){
	//if(debugMode)std::cout<<"absreceta too high!"<<std::endl;
	continue;} // higher abseta cut
      else if ( gendrjt > gendrCut ) {
	//if(debugMode)std::cout<<"gendrjt too high!"<<std::endl;
	continue;}       //delta-r cut, proxy for gen-reco matching quality
      
      // jet/event counts
      h_NJets_kmatCut1->Fill(1);
      if(!hNEvts_withJets_kmatCut1_Filled){
	h_NEvents_withJets_kmatCut1->Fill(1);
	hNEvts_withJets_kmatCut1_Filled=true;      }      
      
      
      // 13 TeV JetID criterion
      bool passesJetID=false; //int jtID=0;
      if(fillMCJetIDHists) 	{
	if (!(absreceta > 2.4)) 
	  passesJetID=(bool)jetID_00eta24( jetIDpt, 
					   neSum_F[jet],  phSum_F[jet],  chSum_F[jet],  eSum_F[jet], muSum_F[jet],
					   numConst,  chMult, useTightJetID);
	else if ( !(absreceta>2.7) && absreceta>2.4 ) 
	  passesJetID=(bool) jetID_24eta27( jetIDpt,
					    neSum_F[jet],  phSum_F[jet], muSum_F[jet],
					    numConst, useTightJetID);
	else if( !(absreceta>3.0) && absreceta>2.7 )
	  passesJetID=(bool) jetID_27eta30( jetIDpt,
					    neSum_F[jet],  phSum_F[jet], 
					    neuMult, useTightJetID);
	else  
	  passesJetID=(bool)jetID_32eta47( jetIDpt, 
					   phSum_F[jet], useTightJetID);
	if(!passesJetID) continue;
      }

      
      if(recomatchindex_I.at(jet)<0){
	std::cerr<<std::endl;
	std::cerr<<"//////////////////////////////////////////////"<<std::endl;
	std::cerr<<"////////////// !!! WARNING !!! ///////////////"<<std::endl;
	std::cerr<<"reco jet w/ nonnegative subid has a negative recomatchindex!!!"<<std::endl;
	std::cerr<<"---------DETAILS--------"<<std::endl;
	std::cerr<<std::endl<<"_______ nEvt="<<nEvt<<" ________"<<std::endl;
	std::cerr<<"ngen_I="<<ngen_I<<std::endl;
	for(int gjet=0; gjet<ngen_I; ++gjet){
	  std::cerr<<"// genmatchindex_I["<<gjet<<"]="<<genmatchindex_I[gjet]<<std::endl;
	  std::cerr<<"genpt_F["<<gjet<<"]="<<genpt_F[gjet]<<std::endl;
	}
	std::cerr<<std::endl;
	std::cerr<<"nref_I="<<nref_I<<std::endl;
	for(int kjet=0;kjet<nref_I;++kjet){
	  std::cerr<<"// recomatchindex_I["<<kjet<<"]="<<recomatchindex_I.at(kjet)<<std::endl;
	  std::cerr<<"subid_I["<<kjet<<"]="<<subid_I[kjet]<<std::endl;
	  std::cerr<<"refpt["<<kjet<<"]="<<refpt_F[kjet]<<std::endl;
	  if(recomatchindex_I[kjet]!=-1){
	    std::cerr<<"genpt["<<recomatchindex_I[kjet]<<"]="<<genpt_F[recomatchindex_I[kjet]]<<std::endl;
	    std::cerr<<"GJetAk04JetPt["<<recomatchindex_I[kjet]<<"]="<<GjetPt->at(recomatchindex_I[kjet])<<std::endl;
	  }
	}
	badgenjetcount++;
	continue;
      }
      
      float genjttheta=2.*atan(exp(-1.*geneta));
      float genjtpz=genpt/tan(genjttheta);
      float genjtE=GjetE->at(recomatchindex_I.at(jet));
      float genjty=0.5*log((genjtE+genjtpz)/(genjtE-genjtpz));//experimentalist version
      float absgenjty=fabs(genjty);
      
      //jetAccepted=true;
      //eventAccepted=true;
      //if(debugMode&&jetAccepted&&eventAccepted)std::cout<<"nEvt="<<nEvt<<", jet="<<jet<<", jet accepted. filling hists"<<std::endl;
      
      h_NJets_kmatCut2->Fill(1);
      if(!hNEvts_withJets_kmatCut2_Filled){
	h_NEvents_withJets_kmatCut2->Fill(1);
	hNEvts_withJets_kmatCut2_Filled=true;      }
      
      
      int theYBin=-1;
      for(int ybin=0;ybin<nbins_absy;ybin++)
	if( absybins[ybin]<=absjty  && 		
	    absjty<absybins[ybin+1]    	      ) {	    
	  theYBin=ybin;
	  break;	  }       
      
      int theGENYBin=-1;
      for(int ybin=0;ybin<nbins_absy;ybin++)
	if( absybins[ybin]<=absgenjty  && 		
	    absgenjty<absybins[ybin+1]    	      ) {	    
	  theGENYBin=ybin;
	  break;	  }       		
      
      int theEtaBin=-1;
      for(int rapbin=0;rapbin<nbins_abseta;rapbin++)
	if( absetabins[rapbin]<=absreceta  && 		
	    absreceta<absetabins[rapbin+1]    	      ) {	    
	  theEtaBin=rapbin;
	  break;	  }       	
      
      int theGENEtaBin=-1;
      for(int rapbin=0;rapbin<nbins_abseta;rapbin++)
	if( absetabins[rapbin]<=absgeneta  && 		
	    absgeneta<absetabins[rapbin+1]    	      ) {	    
	  theGENEtaBin=rapbin;
	  break;	  }     
      
      if(usegenetabins){
	theYBin=theGENYBin;
	theEtaBin=theGENEtaBin;
      }
      
      float recpt_SFdown=-1.;
      float recpt_SFup=-1.;
      if(useJERScaleFactors){
	if(theGENYBin<4&&theGENYBin>=0){
	  recpt_SFdown=recpt;//b.c. SF=1.1 +/- .1 -> (1.0 , 1.1, 1.2)
	  //std::cout<<std::endl<<"genpt="<<genpt<<std::endl;
	  //std::cout<<"theYBin="<<theYBin<<std::endl;
	  //std::cout<<"WAS: recpt="<<recpt<<std::endl;
	  //std::cout<<"WAS: recpt/genpt="<<recpt/genpt<<std::endl;
	  smearPY8wJERScaleFactor(JER_sigmu[theGENYBin],genpt, &recpt);
	  recpt_SFup=genpt*(((recpt/genpt-1.)/1.1*1.2)+1.);//b.c. algebra
	  //std::cout<<"NOW: recpt_SFdown="<<recpt_SFdown<<std::endl;
	  //std::cout<<"NOW: recpt="<<recpt<<std::endl;
	  //std::cout<<"NOW: recpt_SFup="<<recpt_SFup<<std::endl;
	  //if(recpt/genpt<1.)assert(false);
	  //if(nEvt>100)assert(false);
	}
	else {
	  //have to do this since i'm running w/ abseta<2.5 and i don't have a JER fit for that region
	  //h_check_JERScaleFactorSmearing->Fill(recpt,weight_eS);
	  recpt_SFdown=recpt;
	  recpt_SFup=recpt;
	}
      }
      
      /////   UNFOLDING   ///// 
      if(fillMCUnfoldingHists){
	
	hpp_gen->Fill(genpt, weight_eS);
	hpp_reco->Fill(recpt, weight_eS);
	hpp_matrix->Fill(recpt, genpt, weight_eS);	

	if(mcclosureInt%2 == 0){
	  incrementClosureInt=true;
	  hpp_mcclosure_reco_test->Fill(recpt, weight_eS);      	
	  if(fillMCJetCovMatrix){	    
	    bool_helper=true;
	    hpp_mcclosure_covmat_test_helper->Fill(recpt); }	  
	}
	else {
	  incrementClosureInt=true;
	  hpp_mcclosure_gen->Fill(genpt, weight_eS);
	  hpp_mcclosure_reco->Fill(recpt, weight_eS); 	
	  hpp_mcclosure_matrix->Fill(recpt, genpt, weight_eS);	
	}	  
	

	//fill jetspectraRapHists
	if( fillMCUnfJetSpectraRapHists ) { 	  //these should be filled in according to the RECO eta rap bin
	                                          //why; data jet eta is all we can look at, so the resp matrix for a given RECO eta rap bin should reflect that occasionally, a "true" jet in a given RECO eta rap bin actually belongs in a diff GEN eta rap bin
	  
	  if(theEtaBin!=-1){
	    hpp_gen_eta[theEtaBin]->Fill(genpt,weight_eS);
	    hpp_reco_eta[theEtaBin]->Fill(recpt,weight_eS);
	    hpp_matrix_eta[theEtaBin]->Fill(recpt, genpt, weight_eS);
	    
	    if(mcclosureInt_arr[theEtaBin]%2==0){
	      incrementClosureInt_arr[theEtaBin]=true;
	      hpp_mcclosure_reco_test_eta[theEtaBin]->Fill(recpt, weight_eS);      	  
	      if(fillMCJetCovMatrix){
		etabin_bool_helper[theEtaBin]=true;
		hpp_mcclosure_covmat_test_eta_arr_helpers[theEtaBin]->Fill(recpt);	      	    }
	    }
	    else {
	      incrementClosureInt_arr[theEtaBin]=true;
	      hpp_mcclosure_gen_eta[theEtaBin]->Fill(genpt, weight_eS);
	      hpp_mcclosure_reco_eta[theEtaBin]->Fill(recpt, weight_eS); 	
	      hpp_mcclosure_matrix_eta[theEtaBin]->Fill(recpt, genpt, weight_eS);  
	    }	  
	  }//end theEtaBin	  
	  
	  if(theYBin!=-1){
	    hpp_gen_y[theYBin]->Fill(genpt,weight_eS);
	    hpp_reco_y[theYBin]->Fill(recpt,weight_eS);
	    hpp_matrix_y[theYBin]->Fill(recpt, genpt, weight_eS);
	    if(useJERScaleFactors){
	      hpp_reco_y_JERsysup[theYBin]->Fill(recpt_SFup,weight_eS);
	      hpp_matrix_y_JERsysup[theYBin]->Fill(recpt_SFup, genpt, weight_eS);
	      hpp_reco_y_JERsysdown[theYBin]->Fill(recpt_SFdown,weight_eS);
	      hpp_matrix_y_JERsysdown[theYBin]->Fill(recpt_SFdown, genpt, weight_eS);
	    }

	    if(mcclosureInt_y_arr[theYBin]%2==0){
	      incrementClosureInt_y_arr[theYBin]=true;
	      hpp_mcclosure_reco_test_y[theYBin]->Fill(recpt, weight_eS);      	  
	      if(fillMCJetCovMatrix){
		ybin_bool_helper[theYBin]=true;
		hpp_mcclosure_covmat_test_y_arr_helpers[theYBin]->Fill(recpt);	      	    }
	    }
	    else {
	      incrementClosureInt_y_arr[theYBin]=true;
	      hpp_mcclosure_gen_y[theYBin]->Fill(genpt, weight_eS);
	      hpp_mcclosure_reco_y[theYBin]->Fill(recpt, weight_eS); 	
	      hpp_mcclosure_matrix_y[theYBin]->Fill(recpt, genpt, weight_eS);  
	    }	  
	  }//end theYBin	  
	}//end fillmcunfjetspectraraphists	  	  	
	
      }//end fillmcunfoldinghists

    }//end jet loop

    //increment the closure int IF A CLOSURE HIST WAS FILLED IN 0 < |y| < |y_max|
    //if(debugMode && eventAccepted)std::cout<<"mcclosureInt="<<mcclosureInt<<std::endl;
    if(incrementClosureInt)mcclosureInt++;    
    
    //increment the closure int IF A CLOSURE HIST WAS FILLED IN |y_min|[i] < |y| < |y_max|[i+1]
    for(int i=0; i<nbins_abseta; i++){
      if(incrementClosureInt_arr[i]){
	//if(debugMode&&eventAccepted)std::cout<<"mcclosureInt_arr["<<i<<"]="<<mcclosureInt_arr[i]<<std::endl;
	mcclosureInt_arr[i]++;
      }
      if(incrementClosureInt_y_arr[i]){
	//if(debugMode&&eventAccepted)std::cout<<"mcclosureInt_arr["<<i<<"]="<<mcclosureInt_arr[i]<<std::endl;
	mcclosureInt_y_arr[i]++;
      }
    }
    
    
    if(fillMCJetCovMatrix&&fillMCJetIDHists){
      
      if(bool_helper){//don't undergo expensive covmatrix calc if the etabin wasn't filled.
	fillCovMatrix( (TH2D*)hpp_mcclosure_covmat_test, (TH1D*)hpp_mcclosure_covmat_test_helper, nbins_pt , (double) weight_eS);
	hpp_mcclosure_covmat_test_helper->Reset();//reset contents of TH1 without resetting binning, min/max, etc....	  
      }//end bool_helper
      bool_helper=false;	
      
      if(fillMCUnfJetSpectraRapHists){
	for(int i=0; i<nbins_abseta; i++){
	  if(etabin_bool_helper[i]){//don't undergo expensive covmatrix calc if the etabin wasn't filled.
	    fillCovMatrix( (TH2D*)hpp_mcclosure_covmat_test_eta_arr[i], (TH1D*)hpp_mcclosure_covmat_test_eta_arr_helpers[i], nbins_pt , (double) weight_eS);
	    hpp_mcclosure_covmat_test_eta_arr_helpers[i]->Reset();//reset contents of TH1 without resetting binning, min/max, etc....	  
	  }//end etabin_bool_helper	  
	  etabin_bool_helper[i]=false;	
	  if(ybin_bool_helper[i]){//don't undergo expensive covmatrix calc if the etabin wasn't filled.
	    fillCovMatrix( (TH2D*)hpp_mcclosure_covmat_test_y_arr[i], (TH1D*)hpp_mcclosure_covmat_test_y_arr_helpers[i], nbins_pt , (double) weight_eS);
	    hpp_mcclosure_covmat_test_y_arr_helpers[i]->Reset();//reset contents of TH1 without resetting binning, min/max, etc....	  
	  }//end etabin_bool_helper	  
	  ybin_bool_helper[i]=false;	
	}//end etabin loop
      }//end fillMCJetSpectraRapHists
    }//end fillMCJetCovMatrix
    
    //if((nEvt%1000==0) &&
    //   (nEvt!=0)) break;
    
  }//end event loop
  
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Event-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Events in file(s) opened       = " <<h_NEvents->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read from those file(s) = " <<h_NEvents_read->GetEntries()<<std::endl;
  
  
  std::cout<<"Total Num of Events read passing skimCuts                              = " <<
    h_NEvents_skimCut->GetEntries()<<std::endl;
  std::cout<<"Total Num of Events read passing vzCuts and skimCuts (aka good events) = " <<
    h_NEvents_vzCut->GetEntries()<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Total Num of good Events w/ jets                             = " <<
    h_NEvents_withJets->GetEntries()<<std::endl;	   
  std::cout<<"Total Num of good Events, w/ jets, post kmatCut1             = " <<
    h_NEvents_withJets_kmatCut1->GetEntries()<<std::endl; 
  std::cout<<"Total Num of good Events, w/ jets, post kmatCut2             = " <<
    h_NEvents_withJets_kmatCut2->GetEntries()<<std::endl ;
  std::cout<<std::endl;
  
  std::cout<<std::endl<<"/// Job Jet-Loop Summary ///"<<std::endl<<std::endl;
  std::cout<<"Total Num of Jets read from good events                          = " <<
    h_NJets->GetEntries()<<std::endl;  
  std::cout<<"Total Num of Jets read from good events post kmatCut1              = " <<
    h_NJets_kmatCut1->GetEntries()<<std::endl;
  std::cout<<"Total Num of Jets read from good events post kmatCut2              = " <<
    h_NJets_kmatCut2->GetEntries()<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Total Num of good gen jets = "<<goodgenjetcount<<std::endl;
  std::cout<<"Total Num of bad gen jets  = "<<badgenjetcount<<std::endl;
  std::cout<<std::endl;  

  std::cout<<"writing output file "<<outfile<<std::endl;
  fout->Write(); 


  std::cout<<std::endl<<"readForests_ppMC_UNF finished."<<std::endl;  timer.Stop();
  std::cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<std::endl;
  std::cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<std::endl;
  
  return 0;
}

// end readForests_ppMC_JERS
								 

////// main //////
// acts as the frontend control for .exe file
int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=10 && argc!=8 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForests_ppMC_JERS.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForests_ppMC_JERS.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> ";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  // good input, run
  rStatus=1;
  if(argc==1) rStatus = readForests_ppMC_unf();
  else{//read input argument vector
    std::string inputFileList=argv[1]; int startfile= atoi(argv[2]); int endfile= atoi(argv[3]);  
    int jetRadius= atoi(argv[4]); std::string jetType=argv[5];     bool debug=(bool)atoi(argv[6]);
    std::string outputFileName=argv[7];      
    if(argc==minArgs){
      rStatus = readForests_ppMC_unf( inputFileList, startfile, endfile, 
				      jetRadius, jetType, debug,
				      outputFileName);                              }    
    else if(argc==(minArgs+2)){
      float jtEtaCutLo= atof( argv[8] ) ; float jtEtaCutHi=atof( argv[9] ) ;
      rStatus = readForests_ppMC_unf( inputFileList, startfile, endfile,
				      jetRadius, jetType, debug,
				      outputFileName, jtEtaCutLo, jtEtaCutHi);    }
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
