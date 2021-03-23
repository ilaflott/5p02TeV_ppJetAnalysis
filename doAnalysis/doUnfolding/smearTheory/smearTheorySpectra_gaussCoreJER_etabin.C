#include "smearTheorySpectra.h"

//const std::string SPLINE_STR="spline";
//const std::string FIT_STR="fit";
const int CANVX=1200, CANVY=1000;
const int xsecorder=1;//0--> order of calculation (i.e. NLO or NNLO depending on fastNLO setting), 1--> just LO part if available, 2--> just NLO correction to LO part if available... etc.

const bool printBaseDebug=true;
const bool useHistSigmaFit=false;
const bool useJERscaleFactors=true;//note; alters JER syst. as well.



const int nEvents=1e+09;  ///new typical (smaller rel stat unc)
//const int nEvents=1e+08;  ///typical
//const int nEvents=5e+07;  /// debug nevents
//const int nEvents=1e+07;  /// debug nevents
//const int nEvents=1e+06;  /// debug nevents
//const int nEvents=1e+05;  /// debug nevents

//// LO NPCs

//const std::string NPCorrFits_str =HERWGEE4_NPS;    //LO, reasonable par errs
//const std::string NPCorrFits_text=HERWGEE4_NPS_TXT;
//const std::string NPCorrFits_str =HERWGEE5_NPS;    //LO, reasonable par errs
//const std::string NPCorrFits_text=HERWGEE5_NPS_TXT;
//const std::string NPCorrFits_str =_PYTHIA8_NPS;    //LO, reasonable par errs
//const std::string NPCorrFits_text=_PYTHIA8_NPS_TXT;
//const std::string NPCorrFits_str =_POW_PY8_NPS;    //NLO, unusable par errs
//const std::string NPCorrFits_text=_POW_PY8_NPS_TXT;
//const std::string NPCorrFits_str =POWPY8CT_NPS;    //NLO, bad par errs
//const std::string NPCorrFits_text=POWPY8CT_NPS_TXT;
std::string NPCorrFits_str ="f";    //LO, reasonable par errs
//std::string NPCorrFits_text=_PYTHIA8_NPS_TXT;
std::string NPCorrFits_text=_AVERAGE_NPS_TXT;

//SYSTEMATICS
//JER
const bool doJERsys=true;

//NPCs
const bool doNPsys12=false;//involves making a new thy hist, therefore, also a sep spline fit
//NP sys v1, using HERWIG EE4C/PYTHIA8 NPs for unfolding
const std::string NPsys1_CorrFits_str =HERWGEE4_NPS;
const std::string NPsys1_CorrFits_text=HERWGEE4_NPS_TXT;
const std::string NPsys2_CorrFits_str =_PYTHIA8_NPS;
const std::string NPsys2_CorrFits_text=_PYTHIA8_NPS_TXT;
//NP sys v2, shifting the NP fits up/down by 1 sigma or using the upper/lower fits from john that he made himself (how???)
const bool doNPsysupdown=false;//involves making a new thy hist, therefore, also a sep spline fit
const float NPerrfact=1.0;//# sigma to shift NP fit params by. 

//NP systs v1, using POW+PY8/POW+PY8 CTEQ NPs for unfolding, shifting NP fit params up/down by 2 sigma
//const float NPerrfact=1.0;//# sigma to shift NP fit params by. 
//const std::string NPsys1_CorrFits_str =_POW_PY8_NPS;
//const std::string NPsys1_CorrFits_text=_POW_PY8_NPS_TXT;
//const std::string NPsys2_CorrFits_str =POWPY8CT_NPS;
//const std::string NPsys2_CorrFits_text=POWPY8CT_NPS_TXT;

//PDFs
const bool doPDFsys12=false; //involves making a new thy hist from alt PDF choices, therefore, also a sep spline fit
//PDF systs v1; using CT14/HERA pdfs for unfolding, using the 6 pt scale uncertainty + PDF unc w/ err fact == 1
const std::string in_NLOFile_PDFsys1=_CT14FILESTR;
std::string PDFsys1_text=_CT14DESCTXT;
const std::string in_NLOFile_PDFsys2=NNPDFFILESTR;
std::string PDFsys2_text=NNPDFDESCTXT;
const bool doPDFsysupdown=false; //involves making a new thy hist utilizing scale errs, therefore, also a sep spline fit
const float PDFerrfact=1.0;
//const std::string in_NLOFile_PDFsys2=_HERAFILESTR;
//std::string PDFsys2_text=_HERADESCTXT;
//PDF systs v2; using CT10/MMHT pdfs for unfolding, using the 6 pt scale uncertainty + PDF unc w/ err fact == 2
//const float PDFerrfact=1.0;
//const std::string in_NLOFile_PDFsys1=_CT10FILESTR;
//std::string PDFsys1_text=_CT10DESCTXT;
//const std::string in_NLOFile_PDFsys2=_MMHTFILESTR;
//std::string PDFsys2_text=_MMHTDESCTXT;



//valid fitTypes; modLog, 7TeV, spl3wgts, spl3wgts + ext string
int smearTheorySpectra_gaussCoreJER_etabin( std::string in_NLOfileString=in_NLOFile_default, std::string fitType="modLog" , 
					    std::string in_JERFile=in_JERFile_default,  int etabin=0,  std::string outputFile="smearTheory_gaussCoreJER_etabin_defOut.root" 
					    ){
  
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);
 
  TVirtualFitter::SetDefaultFitter("Minuit2");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-02);      //default i think
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-03);      
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04);      
  
  // Output File to write to
  bool NPsfromJohn=false;
  if(NPCorrFile.find("fromJohn")!=std::string::npos)
    NPsfromJohn=true;
  const bool useSplineWeights=((bool)(fitType.find("spl3wgts")!=std::string::npos));
  const bool useSplineExt=(useSplineWeights && ((bool)(fitType.find("ext")!=std::string::npos)));
  const bool useFitWeights=!(useSplineWeights);
  const bool useModLogFit=useFitWeights && ((bool)(fitType=="modLog"));
  const bool use7TeVFit=useFitWeights && ((bool)(fitType=="7TeV"));
  if(!useModLogFit && !use7TeVFit && !useSplineWeights){
    std::cout<<"fit type not recognized! Exit!"<<std::endl; assert(false);}
    
  std::string outputTag=(std::string)absetabins_tags[etabin];
  
  std::cout<<"opening output:"<<outputFile<<std::endl<<std::endl;
  TFile *outf    = new TFile(outputFile.c_str(), "RECREATE" );
  
  // Input, Theory spectrum
  std::string inputFile = "fNLOJetsSpectra/R04/"+in_NLOfileString+".root";
  std::cout<<"opening input file:"<<inputFile<<std::endl<<std::endl;
  TFile* fin_NLO=TFile::Open(inputFile.c_str());

  std::string PDF_text="";
  if     (in_NLOfileString==_CT10FILESTR)PDF_text=_CT10DESCTXT;
  else if(in_NLOfileString==_CT14FILESTR)PDF_text=_CT14DESCTXT;
  else if(in_NLOfileString==_HERAFILESTR)PDF_text=_HERADESCTXT;
  else if(in_NLOfileString==_MMHTFILESTR)PDF_text=_MMHTDESCTXT;
  else if(in_NLOfileString==NNPDFFILESTR)PDF_text=NNPDFDESCTXT;
  else if(in_NLOfileString==_CT14NNLOFILESTR)PDF_text=_CT14NNLODESCTXT;
  else if(in_NLOfileString==_CT14NNLOFILESTR2)PDF_text=_CT14NNLODESCTXT;
  else if(in_NLOfileString==_CT14NNLOFILESTR3)PDF_text=_CT14NNLODESCTXT;
  else assert(false);

  
  
  /////////////// NPCs to NLO jet spectra ////////////////////////////////////////////////////////  
  std::cout<<"opening NP corr file " << NPCorrFile<<std::endl;
  TFile* fin_NP=TFile::Open(NPCorrFile.c_str());  
  if(!fin_NP){
    std::cout<<"error, NP file not found. exit."<<std::endl; assert(false);  }
  
  /////////////// Get NP fit functions
  TF1 *fNP = NULL;
  if(!NPsfromJohn)
    fNP=(TF1*)fin_NP->Get(("f"+NPCorrFits_str+"_etabin"+std::to_string(etabin)).c_str());
  else
    fNP=(TF1*)fin_NP->Get(("fNPC_AVG_R4_etabin"+std::to_string(etabin)).c_str());  
  //  else
  //    fNP=(TF1*)fin_NP->Get(("f"+std::to_string(etabin)).c_str());  

  //fin_NP->ls();
  //std::cout<<"NP name is: "<<"fNPC_AVG_R4_etabin"+std::to_string(etabin)<<std::endl;
  if(!fNP){
    std::cout<<"error, NP fit(s) not found."<<std::endl; assert(false);}
  
  TH1F *NPDataPoints= NULL;
  if(!NPsfromJohn)
    NPDataPoints=(TH1F*)fin_NP->Get( (NPCorrFits_str+"_etabin"+std::to_string(etabin)).c_str() );  
  else 
    NPDataPoints=(TH1F*)fin_NP->Get( ("hNPC_HERWIG_R4_etabin"+std::to_string(etabin)).c_str() );  

  TF1 *fNP_sys1=NULL;
  TH1F* NPsys1_DataPoints=NULL; 
  TF1 *fNP_sys2=NULL;
  TH1F* NPsys2_DataPoints=NULL;
  if(doNPsys12){
    fNP_sys1 = (TF1*)fin_NP->Get(("f"+NPsys1_CorrFits_str+"_etabin"+std::to_string(etabin)).c_str());
    NPsys1_DataPoints= (TH1F*)fin_NP->Get( (NPsys1_CorrFits_str+"_etabin"+std::to_string(etabin)).c_str() );

    fNP_sys2 = (TF1*)fin_NP->Get(("f"+NPsys2_CorrFits_str+"_etabin"+std::to_string(etabin)).c_str());
    NPsys2_DataPoints= (TH1F*)fin_NP->Get( (NPsys2_CorrFits_str+"_etabin"+std::to_string(etabin)).c_str() );
  }
  
  TF1* fNP_sysup=NULL, *fNP_sysdown=NULL;
  if(doNPsysupdown){
    if(!NPsfromJohn){
      //old np systemaics
      std::cout<<"fNP parameters + errs are..."<<std::endl;
      std::cout<<"Par 0="<<fNP->GetParameter(0)<<" +/- "<<fNP->GetParError(0)<<std::endl;
      std::cout<<"Par 1="<<fNP->GetParameter(1)<<" +/- "<<fNP->GetParError(1)<<std::endl;
      std::cout<<"Par 2="<<fNP->GetParameter(2)<<" +/- "<<fNP->GetParError(2)<<std::endl;
      
      fNP_sysup=(TF1*)fNP->Clone( ((std::string)fNP->GetName()+"_sysup").c_str());
      fNP_sysup->SetParameter(0, fNP->GetParameter(0)+NPerrfact*fNP->GetParError(0));
      fNP_sysup->SetParameter(1, fNP->GetParameter(1)+NPerrfact*fNP->GetParError(1));
      fNP_sysup->SetParameter(2, fNP->GetParameter(2)-NPerrfact*fNP->GetParError(2));
      
      std::cout<<"fNP_sysup parameters are..."<<std::endl;
      std::cout<<"Par 0="<<fNP_sysup->GetParameter(0)<<std::endl;
      std::cout<<"Par 1="<<fNP_sysup->GetParameter(1)<<std::endl;
      std::cout<<"Par 2="<<fNP_sysup->GetParameter(2)<<std::endl;
      
      fNP_sysdown=(TF1*)fNP->Clone( ((std::string)fNP->GetName()+"_sysdown").c_str());
      fNP_sysdown->SetParameter(0, fNP->GetParameter(0)-NPerrfact*fNP->GetParError(0));
      fNP_sysdown->SetParameter(1, fNP->GetParameter(1)-NPerrfact*fNP->GetParError(1));
      fNP_sysdown->SetParameter(2, fNP->GetParameter(2)+NPerrfact*fNP->GetParError(2));  
      
      std::cout<<"fNP_sysdown parameters are..."<<std::endl;
      std::cout<<"Par 0="<<fNP_sysdown->GetParameter(0)<<std::endl;
      std::cout<<"Par 1="<<fNP_sysdown->GetParameter(1)<<std::endl;
      std::cout<<"Par 2="<<fNP_sysdown->GetParameter(2)<<std::endl;
    }
    else{
      fNP_sysup=(TF1*)(
		       (TF1*)fin_NP->Get(("fNPC_AVG_R4_etabin"+std::to_string(etabin)+"_up").c_str())
		       )->Clone(("fNPC_AVG_R4_etabin"+std::to_string(etabin)+"_up_clone").c_str());

      fNP_sysdown=(TF1*)(
			 (TF1*)fin_NP->Get(("fNPC_AVG_R4_etabin"+std::to_string(etabin)+"_down").c_str())
			 )->Clone(("fNPC_AVG_R4_etabin"+std::to_string(etabin)+"_down_clone").c_str());
      
    }
  }

  

  //assert(false);
  
  /////////////// Create/Get JER fit function(s)
  std::cout<<"opening JER file " << in_JERFile<<std::endl;
  TFile* fin_JER=TFile::Open(in_JERFile.c_str());
  if(!fin_JER){
    std::cout<<"error, JER file not found or not open."<<std::endl; assert(false); }
  
  std::cout<<"etabin="<<etabin<<std::endl;  
  std::string absetarange_str=absetabins_str[etabin]+" < #||{y} < "+absetabins_str[etabin+1];  
  std::string thyname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"00";
  
  std::string orderstring="LO";
  if(xsecorder==0)orderstring="NN"+orderstring;
  else if(xsecorder==1)orderstring="N"+orderstring;
  else if(xsecorder==2)orderstring=orderstring;
  PDF_text+=" "+orderstring;
  

  int n_thybins=-1;
  double* thybins=setBinning_etabin(etabin, &n_thybins);
  ////why do i do this again?
  //// oh yes, it's because all the |y| bins have the same lowest pt bin, but different high pt bins. 
  //// idea is; rebin in the bins of |y|=0 to get rid of last bin, then use the next few lines to assign the correct higher pt bins for the given |y| bin.
  //// TO DO
  //int n_thybins=theory_rebin->GetNbinsX();
  //std::cout<<"n_thybins="<<n_thybins<<std::endl;
  //std::cout<<"n_thybins+1="<<n_thybins+1<<std::endl;
  //std::cout<<std::endl;
  //double thybins[n_thybins+1]={0.};//+1 because this array represents the bin edges, not the bins themselves
  //for(int i=0; i<(n_thybins+1); i++){
  //  if(i==n_thybins)
  //    thybins[n_thybins]=theory->GetBinLowEdge(i)+theory->GetBinWidth(i);
  //  else{
  //    thybins[i]=theory->GetBinLowEdge(i+1);      
  //  }
  //  std::cout<<"thybins["<<i<<"]="<<thybins[i]<<std::endl;
  //}
  ////assert(false);


  //////////// NLO Thy calculation xsecions get from file ////////////////////////////////////////////////////////
  std::cout<<"plotting NLO xsecs!"<<  std::endl;  
  
  TH1D *theory = (TH1D*)fin_NLO->Get(thyname.c_str());
  //theory->SetTitle((PDF_text+" Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T};"+ddxsec_yax).c_str());
  theory->SetTitle((PDF_text+" #sigma_{jet}^{Incl.};Jet p_{T};"+ddxsec_yax).c_str());
  theory->Sumw2(true);
  theory->Scale(1.e-03);//picobarns to nanobarns
  theory->GetXaxis()->SetNoExponent(true);
  theory->GetXaxis()->SetMoreLogLabels(true);
  theory->SetLineColor(kCyan+3);
  if(printBaseDebug)theory->Print("base");  

  //RELATIVE errors, PDF systematic
  std::string thyPDFNegerrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"01";//hOSPxRx2 == -1. * hOSPxRx1 (symm)  
  TH1D *theoryPDFNegerr = (TH1D*)fin_NLO->Get(thyPDFNegerrname.c_str());
  theoryPDFNegerr=(TH1D*)theoryPDFNegerr->TH1::Rebin(n_thybins,
						     ( ( (std::string) theoryPDFNegerr->GetName() ) +"_rebin").c_str() ,   thybins );
  theoryPDFNegerr->SetTitle((PDF_text+" #sigma_{Jet} (-) PDF Errors;Jet p_{T} [GeV]; #||{(-) Errors} [%]").c_str());
  theoryPDFNegerr->Sumw2(true);
  theoryPDFNegerr->Scale(1.e+02);//rel. err. to % err.
  if(printBaseDebug)  theoryPDFNegerr->Print("Base"); 
  
  std::string thyPDFPoserrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"02";//hOSPxRx2 == -1. * hOSPxRx1 (symm)  
  TH1D *theoryPDFPoserr = (TH1D*)fin_NLO->Get(thyPDFPoserrname.c_str());
  theoryPDFPoserr=(TH1D*)theoryPDFPoserr->TH1::Rebin(n_thybins,
						     ( ( (std::string) theoryPDFPoserr->GetName() ) +"_rebin").c_str() ,   thybins );
  theoryPDFPoserr->SetTitle((PDF_text+" #sigma_{Jet} (+) PDF Errors;Jet p_{T} [GeV]; (+) Errors [%]").c_str());
  theoryPDFPoserr->Sumw2(true);
  theoryPDFPoserr->Scale(1.e+02);//rel. err. to % err.
  if(printBaseDebug)  theoryPDFNegerr->Print("Base"); 
  
  //std::string thy2PtScaleNegerrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"06";//2 PT SCALE UNCERTAINTY, ASYMMETRIC
  //TH1D *theory2PtScaleNegerr = (TH1D*)fin_NLO->Get(thy2PtScaleNegerrname.c_str());
  //theory2PtScaleNegerr->SetTitle((PDF_text+" #sigma_{Jet} (-) 2 Pt. Errors;Jet p_{T} [GeV]; #||{(-) Errors} [%]").c_str());
  //theory2PtScaleNegerr->Sumw2(true);
  //theory2PtScaleNegerr->Scale(1.e+02);//rel. err. to % err.
  //if(printBaseDebug)  theory2PtScaleNegerr->Print("Base"); 
  //
  //std::string thy2PtScalePoserrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"07";//07 IS BROKEN WARNING
  //TH1D *theory2PtScalePoserr = (TH1D*)fin_NLO->Get(thy2PtScalePoserrname.c_str());
  //theory2PtScalePoserr->SetTitle((PDF_text+" #sigma_{Jet} (+) 2 Pt. Errors;Jet p_{T} [GeV]; (+) Errors [%]").c_str());
  //theory2PtScalePoserr->Sumw2(true);
  //theory2PtScalePoserr->Scale(1.e+02);//rel. err. to % err.
  //if(printBaseDebug)  theory2PtScalePoserr->Print("Base"); 
  
  std::string thy6PtScaleNegerrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"08";//2 PT SCALE UNCERTAINTY, ASYMMETRIC
  TH1D *theory6PtScaleNegerr = (TH1D*)fin_NLO->Get(thy6PtScaleNegerrname.c_str());
  theory6PtScaleNegerr=(TH1D*)theory6PtScaleNegerr->TH1::Rebin(n_thybins,
							       ( ( (std::string) theory6PtScaleNegerr->GetName() ) +"_rebin").c_str() ,   thybins );
  theory6PtScaleNegerr->SetTitle((PDF_text+" #sigma_{Jet} (-) 6 Pt. Errors;Jet p_{T} [GeV]; #||{(-) Errors} [%]").c_str());
  theory6PtScaleNegerr->Sumw2(true);
  theory6PtScaleNegerr->Scale(1.e+02);//rel. err. to % err.
  if(printBaseDebug)  theory6PtScaleNegerr->Print("Base"); 
  
  std::string thy6PtScalePoserrname="h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"09";//07 IS BROKEN FOR LO HISTS IIRC
  TH1D *theory6PtScalePoserr = (TH1D*)fin_NLO->Get(thy6PtScalePoserrname.c_str());
  theory6PtScalePoserr=(TH1D*)theory6PtScalePoserr->TH1::Rebin(n_thybins,
						     ( ( (std::string) theory6PtScalePoserr->GetName() ) +"_rebin").c_str() ,   thybins );
  theory6PtScalePoserr->SetTitle((PDF_text+" #sigma_{Jet} (+) 6 Pt. Errors;Jet p_{T} [GeV]; (+) Errors [%]").c_str());
  theory6PtScalePoserr->Sumw2(true);
  theory6PtScalePoserr->Scale(1.e+02);//rel. err. to % err.
  if(printBaseDebug)  theory6PtScalePoserr->Print("Base"); 
  
  std::cout<<"done making "<<orderstring<<" hist."<<std::endl;
  //assert(false);
  //-----------------------------------
  
  //////////// Apply NPs to theory + Plot////////////////////////////////////////////////////////
  std::cout<<"plotting NP+NLO xsecs!"<<  std::endl;  
  TH1D *theory_NP  = (TH1D*) applyNPtoxsec(theory, fNP);
  theory_NP->SetName( ((std::string)theory->GetName()+"_wNP").c_str() );
  
  theory_NP->SetTitle((PDF_text+" + "+NPCorrFits_text+" NPCs, #sigma_{jet}^{Incl.};Jet p_{T} [GeV];"+ddxsec_yax).c_str());
  theory_NP->GetXaxis()->SetNoExponent(true);
  theory_NP->GetXaxis()->SetMoreLogLabels(true);
  theory_NP->SetLineColor(theory->GetLineColor()-1);
  if(printBaseDebug)theory_NP->Print("base");  
  
  std::cout<<"done making NLO+NPhist."<<std::endl;
  
  
  //rebinned theory hists
  TH1D* theory_rebin=(TH1D*)theory->TH1::Rebin(n_thybins, 
					       ( ( (std::string) theory->GetName() ) +"_rebin").c_str() ,   thybins );      
  theory_rebin->SetLineColor(kCyan+4);    
  TH1D* theory_NP_rebin=(TH1D*)theory_NP->TH1::Rebin(n_thybins, 
						     ( ( (std::string) theory_NP->GetName() ) +"_rebin").c_str() ,     thybins ); 
  theory_NP_rebin->SetLineColor(kCyan+4);    
  
  
  TH1D* theoryRelStaterr=(TH1D*)theory->Clone( ("h"+std::to_string(xsecorder)+"100"+std::to_string(etabin+1)+"03").c_str());
  for(int i=1; i<=theoryRelStaterr->GetNbinsX();i++){
    float val=theoryRelStaterr->GetBinContent(i);
    float valerr=theoryRelStaterr->GetBinError(i);
    theoryRelStaterr->SetBinContent(i, (valerr/val)*100.);//in units of % for easy interpretation.
    theoryRelStaterr->SetBinError(i, 0.);    
  }
  



  
  
  //----------------------------------- PDF SYSTEMATICS 
  //v1 HISTS - using different PDFs to evaluate unc due to PDF choice
  TFile* fin_PDFsys1=NULL, * fin_PDFsys2=NULL;
  TH1D* theory_PDFsys1_PDFerrs=NULL, *theory_PDFsys2_PDFerrs=NULL; //FOR SETTING ERRS ON HIST
  TH1D* theory_PDFsys1=NULL, *theory_PDFsys2=NULL; //NO NPCs
  TH1D* theory_PDFsys1_rebin=NULL, *theory_PDFsys2_rebin=NULL; //NO NPCs
  TH1D* theory_NP_PDFsys1=NULL, *theory_NP_PDFsys2=NULL; //WITH NPCs
  TH1D* theory_NP_PDFsys1_rebin=NULL, *theory_NP_PDFsys2_rebin=NULL; //WITH NPCs
  if(doPDFsys12){
    std::cout<<"NOW making theory PDF sys1/2 hists."<<std::endl;
    
    std::cout<<"opening PDFsys1 file"<<std::endl;
    fin_PDFsys1=TFile::Open( ("fNLOJetsSpectra/R04/"+in_NLOFile_PDFsys1+".root").c_str());
    if(!fin_PDFsys1)assert(false);
    PDFsys1_text+=" "+orderstring;
    
    theory_PDFsys1=(TH1D*)fin_PDFsys1->Get(thyname.c_str());
    theory_PDFsys1->SetName( ((std::string)theory->GetName()+"_PDFsys1").c_str() );    
    theory_PDFsys1->SetTitle((PDFsys1_text+" Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_PDFsys1->GetXaxis()->SetNoExponent(true);
    theory_PDFsys1->GetXaxis()->SetMoreLogLabels(true);
    theory_PDFsys1->SetLineColor(kTeal);
    if(printBaseDebug)theory_PDFsys1->Print("base");  
    
    //set the PDF errs, 
    theory_PDFsys1_PDFerrs=(TH1D*)fin_PDFsys1->Get(thyPDFPoserrname.c_str());
    theory_PDFsys1_PDFerrs->Scale(1.e+02);//rel. err --> % err
    theory_PDFsys1_PDFerrs->SetName( ((std::string)theory->GetName()+"_PDFsys1_Poserr").c_str() );   
     for(int i=1; i<=(theory_PDFsys1->GetNbinsX());i++ )
      theory_PDFsys1->SetBinError(i, theory_PDFsys1->GetBinContent(i)*(theory_PDFsys1_PDFerrs->GetBinContent(i))/100.);
    
    
    std::cout<<"opening PDFsys2 file"<<std::endl;
    fin_PDFsys2=TFile::Open( ("fNLOJetsSpectra/R04/"+in_NLOFile_PDFsys2+".root").c_str());
    if(!fin_PDFsys2)assert(false);
    PDFsys2_text+=" "+orderstring;
    
    theory_PDFsys2=(TH1D*)fin_PDFsys2->Get(thyname.c_str());
    theory_PDFsys2->SetName( ((std::string)theory->GetName()+"_PDFsys2").c_str() );
    
    theory_PDFsys2->SetTitle((PDFsys2_text+" Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_PDFsys2->GetXaxis()->SetNoExponent(true);
    theory_PDFsys2->GetXaxis()->SetMoreLogLabels(true);
    theory_PDFsys2->SetLineColor(kTeal);
    if(printBaseDebug)theory_PDFsys2->Print("base");  
    
    //set the PDF errs, 
    theory_PDFsys2_PDFerrs=(TH1D*)fin_PDFsys2->Get(thyPDFPoserrname.c_str());
    theory_PDFsys2_PDFerrs->Scale(1.e+02);//rel. err --> % err
    theory_PDFsys2_PDFerrs->SetName( ((std::string)theory->GetName()+"_PDFsys2_Poserr").c_str() );    
    for(int i=1; i<=(theory_PDFsys2->GetNbinsX());i++ )
      theory_PDFsys2->SetBinError(i, theory_PDFsys2->GetBinContent(i)*(theory_PDFsys2_PDFerrs->GetBinContent(i))/100.);
    //    for(int i=1; i<=(theory_PDFsys2->GetNbinsX());i++ )
    //theory_PDFsys2->SetBinError(i, theory_PDFsys2_PDFerrs->GetBinContent(i));

    theory_PDFsys1_rebin=(TH1D*)theory_PDFsys1->TH1::Rebin(n_thybins, 
							   ( ( (std::string) theory_PDFsys1->GetName() ) +"_rebin").c_str() ,     thybins ); 
    theory_PDFsys2_rebin=(TH1D*)theory_PDFsys2->TH1::Rebin(n_thybins, 
							   ( ( (std::string) theory_PDFsys2->GetName() ) +"_rebin").c_str() ,     thybins ); 
    

    //else assert(false);//can't do np sysup/down for more forward bins yet; TO DO
    
    std::cout<<"done making theory PDF sys1/2 hists."<<std::endl;

    // ------------------------------------------------------------- //

    std::cout<<"NOW making theory_NP PDF sys1/2 hists."<<std::endl;
    
    theory_NP_PDFsys1  = (TH1D*) applyNPtoxsec(theory_PDFsys1, fNP);
    theory_NP_PDFsys1->SetTitle((PDFsys1_text+"+"+NPCorrFits_text+" NPCs, Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NP_PDFsys1->GetXaxis()->SetNoExponent(true);
    theory_NP_PDFsys1->GetXaxis()->SetMoreLogLabels(true);
    theory_NP_PDFsys1->SetLineColor(kTeal);
    if(printBaseDebug)theory_NP_PDFsys1->Print("base");  
    
    theory_NP_PDFsys2  = (TH1D*) applyNPtoxsec(theory_PDFsys2, fNP);    
    theory_NP_PDFsys2->SetTitle((PDFsys2_text+"+"+NPCorrFits_text+" NPCs, Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NP_PDFsys2->GetXaxis()->SetNoExponent(true);
    theory_NP_PDFsys2->GetXaxis()->SetMoreLogLabels(true);
    theory_NP_PDFsys2->SetLineColor(kTeal);
    if(printBaseDebug)theory_NP_PDFsys2->Print("base");  
    
    theory_NP_PDFsys1_rebin=(TH1D*)theory_NP_PDFsys1->TH1::Rebin(n_thybins, 
								 ( ( (std::string) theory_NP_PDFsys1->GetName() ) +"_rebin").c_str() ,     thybins ); 
    theory_NP_PDFsys2_rebin=(TH1D*)theory_NP_PDFsys2->TH1::Rebin(n_thybins, 
								 ( ( (std::string) theory_NP_PDFsys2->GetName() ) +"_rebin").c_str() ,     thybins ); 
      
    std::cout<<"done making theory_NP PDF sys1/2 hists."<<std::endl; 
    

  }
  
  
  //v2 HISTS - using 6 pt scale errors to evaluate unc due to PDF calculation
  TH1D* theory_PDFsysup=NULL, *theory_PDFsysdown=NULL; //NO NPCs
  TH1D* theory_PDFsysup_rebin=NULL, *theory_PDFsysdown_rebin=NULL; //NO NPCs
  TH1D* theory_NP_PDFsysup=NULL, *theory_NP_PDFsysdown=NULL; //WITH NPCs
  TH1D* theory_NP_PDFsysup_rebin=NULL, *theory_NP_PDFsysdown_rebin=NULL; //WITH NPCs  
  if(doPDFsysupdown){
    
    // ------------------------------------------------------------- //
    
    std::cout<<"NOW making theory PDF sysup/down hists."<<std::endl;
    theory_PDFsysup=(TH1D*)theory_rebin->Clone(((std::string)theory_rebin->GetName()+"_PDFsysup").c_str());
    theory_PDFsysup= (TH1D*) makePDFsys(theory_rebin, theory6PtScalePoserr, "sysup", PDFerrfact);//add/subtract 6pt scale unc only
    //theory_PDFsysup= (TH1D*) makePDFsys_wPDFerr(theory, theory6PtScalePoserr, "sysup", PDFerrfact);//add/subtract 6pt scale unc w/ PDF unc in quad. "v2"
    theory_PDFsysup->SetTitle((PDF_text+" sysup Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_PDFsysup->GetXaxis()->SetNoExponent(true);
    theory_PDFsysup->GetXaxis()->SetMoreLogLabels(true);
    theory_PDFsysup->SetLineColor(kTeal);
    if(printBaseDebug)theory_PDFsysup->Print("base");  
    
    theory_PDFsysdown=(TH1D*)theory_rebin->Clone(((std::string)theory_rebin->GetName()+"_PDFsysdown").c_str());
    theory_PDFsysdown  = (TH1D*) makePDFsys(theory_rebin, theory6PtScaleNegerr, "sysdown", PDFerrfact);// add/subtract 6pt scale unc only
    //theory_PDFsysdown  = (TH1D*) makePDFsys_wPDFerr(theory, theory6PtScaleNegerr, "sysdown", PDFerrfact);//add/subtract 6pt scale unc w/ PDF unc in quad. "v2"
    theory_PDFsysdown->SetTitle((PDF_text+" sysdown Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_PDFsysdown->GetXaxis()->SetNoExponent(true);
    theory_PDFsysdown->GetXaxis()->SetMoreLogLabels(true);
    theory_PDFsysdown->SetLineColor(kTeal);
    if(printBaseDebug)theory_PDFsysdown->Print("base");  
    
    theory_PDFsysup_rebin=(TH1D*)theory_PDFsysup->TH1::Rebin(n_thybins, 
							     ( ( (std::string) theory_PDFsysup->GetName() ) +"_rebin").c_str() ,     thybins );
    theory_PDFsysdown_rebin=(TH1D*)theory_PDFsysdown->TH1::Rebin(n_thybins, 
								 ( ( (std::string) theory_PDFsysdown->GetName() ) +"_rebin").c_str() ,     thybins ); 
    
    std::cout<<"done making theory PDF sysup/down hists."<<std::endl;
    
    // ------------------------------------------------------------- //
    
    std::cout<<"NOW making theory_NP PDF sysup/down hists."<<std::endl;
    
    theory_NP_PDFsysup  = (TH1D*) applyNPtoxsec(theory_PDFsysup, fNP);
    theory_NP_PDFsysup->SetTitle((PDF_text+"+"+NPCorrFits_text+" NPCs, sysup Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NP_PDFsysup->GetXaxis()->SetNoExponent(true);
    theory_NP_PDFsysup->GetXaxis()->SetMoreLogLabels(true);
    theory_NP_PDFsysup->SetLineColor(kTeal);
    if(printBaseDebug)theory_NP_PDFsysup->Print("base");  
    
    theory_NP_PDFsysdown  = (TH1D*) applyNPtoxsec(theory_PDFsysdown, fNP);    
    theory_NP_PDFsysdown->SetTitle((PDF_text+"+"+NPCorrFits_text+" NPCs, sysdown Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NP_PDFsysdown->GetXaxis()->SetNoExponent(true);
    theory_NP_PDFsysdown->GetXaxis()->SetMoreLogLabels(true);
    theory_NP_PDFsysdown->SetLineColor(kTeal);
    if(printBaseDebug)theory_NP_PDFsysdown->Print("base");  
    
    theory_NP_PDFsysup_rebin=(TH1D*)theory_NP_PDFsysup->TH1::Rebin(n_thybins, 
								   ( ( (std::string) theory_NP_PDFsysup->GetName() ) +"_rebin").c_str() ,     thybins ); 
    theory_NP_PDFsysdown_rebin=(TH1D*)theory_NP_PDFsysdown->TH1::Rebin(n_thybins, 
								       ( ( (std::string) theory_NP_PDFsysdown->GetName() ) +"_rebin").c_str() ,     thybins ); 
    
    //else assert(false);//can't do np sysup/down for more forward bins yet; TO DO
    
    std::cout<<"done making theory_NP PDF sysup/down hists."<<std::endl; 
    
    // ------------------------------------------------------------- //
    
  }
  



  //----------------------------------- NP SYSTEMATICS HISTS
  TH1D *theory_NPsys1=NULL, *theory_NPsys2=NULL;// using alt NP fit choices to eval syst err
  TH1D *theory_NPsys1_rebin=NULL, *theory_NPsys2_rebin=NULL;
  TH1D *theory_NPsys1_ratio=NULL, *theory_NPsys2_ratio=NULL;//QA to ensure NPs are being applied correctly. 
  TH1D * theory_NPsys1_ratio_NPsys2;
  if(doNPsys12){
    
    //ALT NP CHOICE 1
    std::cout<<"making theory NP sys1/2 hists"<<std::endl;
    theory_NPsys1  = (TH1D*) applyNPtoxsec(theory, fNP_sys1);
    theory_NPsys1->SetName( ((std::string)theory->GetName()+"_wNPsys1").c_str() );
    
    theory_NPsys1->SetTitle((PDF_text+"+"+NPsys1_CorrFits_text+" NPCs, Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NPsys1->GetXaxis()->SetNoExponent(true);
    theory_NPsys1->GetXaxis()->SetMoreLogLabels(true);
    theory_NPsys1->SetLineColor(kTeal);
    if(printBaseDebug)theory_NPsys1->Print("base");  
    
    theory_NPsys1_ratio=(TH1D*)theory_NPsys1->Clone(((std::string)theory_NPsys1->GetName() + "_ratio").c_str());
    theory_NPsys1_ratio->Divide(theory);
    
    //ALT NP CHOICE 2
    theory_NPsys2  = (TH1D*) applyNPtoxsec(theory, fNP_sys2);
    theory_NPsys2->SetName( ((std::string)theory->GetName()+"_wNPsys2").c_str() );
    
    theory_NPsys2->SetTitle((PDF_text+"+"+NPsys2_CorrFits_text+" NPCs, Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NPsys2->GetXaxis()->SetNoExponent(true);
    theory_NPsys2->GetXaxis()->SetMoreLogLabels(true);
    theory_NPsys2->SetLineColor(kTeal);
    if(printBaseDebug)theory_NPsys2->Print("base");  
    
    theory_NPsys2_ratio=(TH1D*)theory_NPsys2->Clone(((std::string)theory_NPsys2->GetName() + "_ratio").c_str());
    theory_NPsys2_ratio->Divide(theory);
    
    theory_NPsys1_rebin=(TH1D*)theory_NPsys1->TH1::Rebin(n_thybins, 
							 ( ( (std::string) theory_NPsys1->GetName() ) +"_rebin").c_str() ,     thybins ); 
    theory_NPsys2_rebin=(TH1D*)theory_NPsys2->TH1::Rebin(n_thybins, 
							 ( ( (std::string) theory_NPsys2->GetName() ) +"_rebin").c_str() ,     thybins );     
    //SYS UP/DOWN QA -- still does alt 1/2 QA, TO DO
    theory_NPsys1_ratio_NPsys2=(TH1D*)theory_NPsys1->Clone("theory_NPsys1_ratio_NPsys2");
    theory_NPsys1_ratio_NPsys2->Divide(theory_NPsys2);
    std::cout<<"done making theory NP sys1/2 hists"<<std::endl;

  }

  TH1D *theory_NPsysup=NULL, *theory_NPsysdown=NULL;   // varying the NP fit up/down
  TH1D *theory_NPsysup_rebin=NULL, *theory_NPsysdown_rebin=NULL;
  if(doNPsysupdown){
    
    //NP SYSUP
    std::cout<<"making theory NP sysup/down hists"<<std::endl;
    theory_NPsysup  = (TH1D*) applyNPtoxsec(theory, fNP_sysup);
    theory_NPsysup->SetName( ((std::string)theory->GetName()+"_wNPsysup").c_str() );
    
    theory_NPsysup->SetTitle((PDF_text+"+"+NPCorrFits_text+" sysup NPCs, Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NPsysup->GetXaxis()->SetNoExponent(true);
    theory_NPsysup->GetXaxis()->SetMoreLogLabels(true);
    theory_NPsysup->SetLineColor(kTeal);
    if(printBaseDebug)theory_NPsysup->Print("base");  
    
    //theory_NPsysup_ratio=(TH1D*)theory_NPsysup->Clone(((std::string)theory_NPsysup->GetName() + "_ratio").c_str());
    //theory_NPsysup_ratio->Divide(theory);
    
    //NP SYSDOWN
    theory_NPsysdown  = (TH1D*) applyNPtoxsec(theory, fNP_sysdown);
    theory_NPsysdown->SetName( ((std::string)theory->GetName()+"_wNPsysdown").c_str() );
    
    theory_NPsysdown->SetTitle((PDF_text+"+"+NPCorrFits_text+" sysdown NPCs, Inclusive #sigma_{jet}, "+absetarange_str+";Jet p_{T} [GeV];"+ddxsec_yax).c_str());
    theory_NPsysdown->GetXaxis()->SetNoExponent(true);
    theory_NPsysdown->GetXaxis()->SetMoreLogLabels(true);
    theory_NPsysdown->SetLineColor(kTeal);
    if(printBaseDebug)theory_NPsysdown->Print("base");  
    
    //theory_NPsysdown_ratio=(TH1D*)theory_NPsysdown->Clone(((std::string)theory_NPsysdown->GetName() + "_ratio").c_str());
    //theory_NPsysdown_ratio->Divide(theory);
    

    theory_NPsysup_rebin=(TH1D*)theory_NPsysup->TH1::Rebin(n_thybins, 
							   ( ( (std::string) theory_NPsysup->GetName() ) +"_rebin").c_str() ,     thybins ); 
    theory_NPsysdown_rebin=(TH1D*)theory_NPsysdown->TH1::Rebin(n_thybins, 
							       ( ( (std::string) theory_NPsysdown->GetName() ) +"_rebin").c_str() ,     thybins ); 
    
    //else {
    //  std::cout<<"error; don't do systematics for |y| bins starting beyond 1.0 for now. It's on your to do list to add this Ian."<<std::endl;
    //  assert(false);//can't do np sysup/down for more forward bins yet; TO DO
    //}
    



  }
  
  //TF1 *fJER_default = new TF1("fJER_def",gJER_ynew_str.c_str(),thybins[0], thybins[n_thybins]);  //for comparison only
  //fJER_default->SetLineColor(kBlue);
  //fJER_default->SetMinimum(0.); //do i need to set min/max really?
  //fJER_default->SetMaximum(0.25);
  
  TF1 *fJER = NULL;
  if(useHistSigmaFit)
    fJER=(TF1*)fin_JER->Get("SigmaFit_h");
  else
    fJER=(TF1*)fin_JER->Get("SigmaFit_f");
  if(!fJER){
    std::cout<<"error, JER fit not open or not found!."<<std::endl; assert(false);}

  fJER->SetLineColor(kRed);
  fJER->SetMinimum(0.); //do i need to set min/max really?
  fJER->SetMaximum(0.25);
  std::cout<<std::endl;

  //JER systematics
  TF1* fJER_sysup=NULL, * fJER_sysdown=NULL;
  if(doJERsys){
    std::string jersysup_name="SigmaFit_sysup_";
    if(useHistSigmaFit)jersysup_name+="h";
    else jersysup_name+="f";    
    //fJER_sysup=(TF1*)fin_JER->Get(jersysup_name.c_str()); //old way of JERsyst
    fJER_sysup=(TF1*)fJER->Clone(jersysup_name.c_str());
    std::string jersysdown_name="SigmaFit_sysdown_";
    if(useHistSigmaFit)jersysdown_name+="h";
    else jersysdown_name+="f";    
    //fJER_sysdown=(TF1*)fin_JER->Get(jersysdown_name.c_str());    //old way of JERsyst
    fJER_sysdown=(TF1*)fJER->Clone(jersysdown_name.c_str());
  }
  
  if(useJERscaleFactors){
    std::cout<<"using JER scale factors!!"<<std::endl;
    //float par0scalefact_err= fJER->GetParameter(0)*JERscaleFacts[etabin]*sqrt( pow( ( fJER->GetParError(0)/fJER->GetParameter(0) ), 2) +
    //										 pow( ( JERscaleFactErrs[etabin]/JERscaleFacts[etabin] ), 2) +
    //										 );//i wanna try this out
    fJER->SetParameter(0, fJER->GetParameter(0)*JERscaleFacts[etabin]);
    fJER->SetParameter(1, fJER->GetParameter(1)*JERscaleFacts[etabin]);
    if(doJERsys){
      std::cout<<"doing JER systematics + using JER scale factors!!"<<std::endl;

      fJER_sysup->SetParameter(  0,   fJER_sysup->GetParameter(0)*( JERscaleFacts[etabin]+JERscaleFactErrs[etabin] ) );//orig            
      fJER_sysdown->SetParameter(0, fJER_sysdown->GetParameter(0)*( JERscaleFacts[etabin]-JERscaleFactErrs[etabin] ) );//orig
      fJER_sysup->SetParameter(  1,   fJER_sysup->GetParameter(1)*( JERscaleFacts[etabin]+JERscaleFactErrs[etabin] ) );//orig
      fJER_sysdown->SetParameter(1, fJER_sysdown->GetParameter(1)*( JERscaleFacts[etabin]-JERscaleFactErrs[etabin] ) );//orig
    }
  }
  



  /////////////// plot Core p_T Resolution ////////////////////////////////////////////////////////
  std::cout<<"plotting JER!"<<  std::endl;  
  //TF1 is irritating to set titles for+blah... so use this  
  TH1D * hJER = (TH1D*)( (TH1D*)fJER->GetHistogram() )->Clone(("hJER_"+std::to_string(etabin)).c_str());
  if(!useHistSigmaFit)hJER->SetTitle(("Gauss Core JER Fit, "+absetarange_str+"; Jet p_{T} ; #sigma / #mu from Fit").c_str());
  else if(useHistSigmaFit)hJER->SetTitle(("Gauss Core JER Fit, "+absetarange_str+"; Jet p_{T} ; #sigma from Fit").c_str());
  hJER->GetXaxis()->SetNoExponent(true);
  hJER->GetXaxis()->SetMoreLogLabels(true);
  if(useHistSigmaFit)hJER->SetLineColor(kBlue-7);
  else hJER->SetLineColor(kMagenta-2);
  
  TH1D* hJER_sysup=NULL , *hJER_sysdown=NULL;
  if(doJERsys) hJER_sysup   = (TH1D*)( (TH1D*)fJER_sysup->GetHistogram() )->Clone(("hJER_sysup_"+std::to_string(etabin)).c_str());
  if(doJERsys) hJER_sysdown = (TH1D*)( (TH1D*)fJER_sysdown->GetHistogram() )->Clone(("hJER_sysdown_"+std::to_string(etabin)).c_str());  
  
  //-----------------------------------  
  TCanvas *plot_JER = new TCanvas("plot_JER", "plot JER",900,600);
  plot_JER->cd()->SetLogx(1);
  plot_JER->cd();  
  hJER->SetAxisRange(0.02,0.20,"Y");

  hJER->DrawClone("HIST E ][");    
  if(doJERsys)hJER_sysup->DrawClone("HIST E SAME ][");
  if(doJERsys)hJER_sysdown->DrawClone("HIST E SAME ][");
  //-----------------------------------

  //TO DO REPLACE ME W/TEXT SAYING SAME THING
  TLegend *leg=new TLegend(0.53,0.50,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0.);
  leg->AddEntry(hJER, "PY8 CUETP8M1 JER", "l");  
  if(doJERsys)leg->AddEntry(hJER_sysup, "Upper/Lower Unc.", "l");  
  leg->AddEntry((TObject*)0,"ak4PF Jets", "");
  if(useJERscaleFactors)leg->AddEntry((TObject*)0, "8 TeV JER Scale Factors Applied", "");
  //leg->AddEntry((TObject*)0,"0<|y|<2.0", "");   
  leg->Draw();
  //-----------------------------------
  
  /////////////// plots to to check NPs 
  // 2x2 canv of |y| bins
  TLine* lineatone=new TLine(thybins[0],1.,thybins[n_thybins],1.); lineatone->SetLineStyle(7);
  TCanvas *plot_NP = new TCanvas("plot_NP", ("plot y"+std::to_string(etabin)+" NPs").c_str(),1200,1000);
  //plot_y0_to_y3_NPs->Divide(2,2);
  //-----------------------------------
  plot_NP->cd()->SetLogx(1);
  plot_NP->cd();
  TH1F *hNP = (TH1F*)( (TH1F*)fNP->GetHistogram()
		       )->Clone( ("hNP_y"+std::to_string(etabin)).c_str());
  hNP->SetLineColor(kBlue);
  
  TH1F *hNP_sys1=NULL, *hNP_sys2=NULL;
  TH1F *hNP_sys1_ratio_sysdown=NULL;
  if(doNPsys12){
    //ALT CHOICE NP1/2 sys hists
    hNP_sys1 = (TH1F*)( (TH1F*)fNP_sys1->GetHistogram()
			 )->Clone( ("hNP_sys1_y"+std::to_string(etabin)).c_str());
    hNP_sys1->SetLineColor(kRed);
    hNP_sys2 = (TH1F*)( (TH1F*)fNP_sys2->GetHistogram()
			   )->Clone( ("hNP_sys2_y"+std::to_string(etabin)).c_str());
    hNP_sys2->SetLineColor(kRed); 
    //hNP->SetMaximum(hNP_sys1->GetMaximum()*1.2);
    //hNP->SetMinimum(hNP_sys2->GetMinimum()*0.8);    

    //SYSUP+SYSDOWN QA -- still does sys1/2 QA, TO DO
    hNP_sys1_ratio_sysdown=(TH1F*)hNP_sys1->Clone("hNP_sys1_ratio_sysdown");
    hNP_sys1_ratio_sysdown->Divide(hNP_sys2);
  }

  TH1F *hNP_sysup=NULL, *hNP_sysdown=NULL;
  if(doNPsysupdown){
    //NP sysup/down hists
    hNP_sysup = (TH1F*)( (TH1F*)fNP_sysup->GetHistogram()
			 )->Clone( ("hNP_sysup_y"+std::to_string(etabin)).c_str());
    hNP_sysup->SetLineColor(kRed);
    hNP_sysdown = (TH1F*)( (TH1F*)fNP_sysdown->GetHistogram()
			   )->Clone( ("hNP_sysdown_y"+std::to_string(etabin)).c_str());
    hNP_sysdown->SetLineColor(kRed); 
    hNP->SetMaximum(hNP_sysup->GetMaximum()*1.2);
    hNP->SetMinimum(hNP_sysdown->GetMinimum()*0.8);        
  }
  
  
  hNP->SetTitle( (NPCorrFits_text+" NPCs, "+absetarange_str+";Jet p_{T};NP Corr. Factor" ).c_str() );
  hNP->GetXaxis()->SetNoExponent(true);
  hNP->GetXaxis()->SetMoreLogLabels(true);
  hNP->SetAxisRange(0.5, 1.5, "Y");
  
  hNP->DrawClone("HIST ][");  
  NPDataPoints->DrawClone("E SAME");
  if(doNPsysupdown){
    hNP_sysup->DrawClone("HIST ][ SAME");
    hNP_sysdown->DrawClone("HIST ][ SAME");
  }
  
  TLegend* NPleg=new TLegend(0.15,0.75,0.45,0.85, "", "NDC");
  NPleg->SetFillStyle(0); NPleg->SetBorderSize(0);
  NPleg->AddEntry(NPDataPoints,(NPCorrFits_text).c_str(),"lp");
  NPleg->AddEntry(hNP,(NPCorrFits_text+" NPC Fit").c_str(),"l");
  if(doNPsysupdown)NPleg->AddEntry(hNP_sysup,"Upper/Lower Fit Uncertainty","l");
  NPleg->Draw();
  lineatone->SetX2( hNP->GetBinLowEdge(hNP->GetNbinsX())+ hNP->GetBinWidth(hNP->GetNbinsX()) );
  lineatone->DrawClone();
  //-----------------------------------  
  
  

  

  
  
  
  
  /////////////// draw xsec hists from file
  
  //2X2 CANV THEORY NLO xsec
  TCanvas *plot_NLOxsec = new TCanvas("plot_NLOxsec", "plot NLOxsec",1200,1000);
  //--------------------------------------------
  plot_NLOxsec->cd()->SetLogx(1);
  plot_NLOxsec->cd()->SetLogy(1);
  plot_NLOxsec->cd();
  theory->DrawClone("HIST E ][");    
  if(doPDFsysupdown)theory_PDFsysup->DrawClone("HIST SAME ][");
  if(doPDFsysupdown)theory_PDFsysdown->DrawClone("HIST SAME ][");
  
  TLegend* nloxsecleg= new TLegend(0.6,0.7,0.9,0.9);  
  nloxsecleg->SetBorderSize(0.);
  nloxsecleg->SetFillStyle(0.);
  nloxsecleg->AddEntry(theory,(orderstring+" "+PDF_text).c_str(),"l");
  if(doPDFsysupdown)nloxsecleg->AddEntry(theory_PDFsysup, "Upper/Lower 6 Pt. Scale Unc.", "l");
  nloxsecleg->AddEntry((TObject*)0,"R=0.4 Jets","");
  nloxsecleg->AddEntry((TObject*)0, absetarange_str.c_str(), "");
  nloxsecleg->Draw();    
  //--------------------------------------------
  
  //--------------------------------------------
  //2X2 CANV THEORY NPNLO xsce
  TCanvas *plot_NPNLOxsec = new TCanvas("plot_NPNLOxsec", "plot NPNLOxsec",1200,1000);
  plot_NPNLOxsec->cd()->SetLogx(1);
  plot_NPNLOxsec->cd()->SetLogy(1);
  plot_NPNLOxsec->cd();
  theory_NP->DrawClone("HIST E ][");      
  if(doPDFsysupdown)theory_NP_PDFsysup->DrawClone("HIST SAME ][");
  if(doPDFsysupdown)theory_NP_PDFsysdown->DrawClone("HIST SAME ][");
  
  TLegend* npnloxsecleg= new TLegend(0.55,0.7,0.9,0.9);  
  npnloxsecleg->SetBorderSize(0.);
  npnloxsecleg->SetFillStyle(0.);
  npnloxsecleg->AddEntry(theory_NP,(PDF_text+" + "+NPCorrFits_text).c_str(),"l");
  if(doPDFsysupdown)npnloxsecleg->AddEntry(theory_NP_PDFsysup, "Upper/Lower 6 Pt. Scale Unc.", "l");
  npnloxsecleg->AddEntry((TObject*)0,"R=0.4 Jets","");
  npnloxsecleg->AddEntry((TObject*)0, absetarange_str.c_str(), "");
  npnloxsecleg->Draw();  
  
  //--------------------------------------------
  //PREP FOR BOTH POS AND NEG SCALE ERR PLOTS


  float scaleerrmin=10000000000.;
  float scaleerrmax=-1.;


  theory6PtScaleNegerr->Scale(-1.);
  //theory2PtScaleNegerr->Scale(-1.);
  theoryPDFNegerr     ->Scale(-1.);
  
  //if(theory6PtScaleNegerr->GetMinimum()<scaleerrmin)scaleerrmin=theory6PtScaleNegerr->GetMinimum();
  //if(theoryPDFNegerr->GetMinimum()<scaleerrmin)scaleerrmin=theoryPDFNegerr->GetMinimum();
  if(getnonzeromin(theory6PtScaleNegerr)<scaleerrmin)scaleerrmin=getnonzeromin(theory6PtScaleNegerr);
  if(getnonzeromin(theoryPDFNegerr     )<scaleerrmin)scaleerrmin=getnonzeromin(theoryPDFNegerr     );
  if(getnonzeromin(theoryRelStaterr    )<scaleerrmin)scaleerrmin=getnonzeromin(theoryRelStaterr    );
  if(getnonzeromin(theory6PtScalePoserr)<scaleerrmin)scaleerrmin=getnonzeromin(theory6PtScalePoserr);
  if(getnonzeromin(theoryPDFPoserr     )<scaleerrmin)scaleerrmin=getnonzeromin(theoryPDFPoserr     );
  scaleerrmin*=0.1;
  
  //if(theory6PtScalePoserr->GetMinimum()<scaleerrmin)scaleerrmin=theory6PtScalePoserr->GetMinimum();
  //if(theoryPDFPoserr->GetMinimum()<scaleerrmin)scaleerrmin=theoryPDFPoserr->GetMinimum();
  
  if(theory6PtScalePoserr->GetMaximum()>scaleerrmax)scaleerrmax=theory6PtScalePoserr->GetMaximum();
  if(theoryPDFPoserr     ->GetMaximum()>scaleerrmax)scaleerrmax=theoryPDFPoserr     ->GetMaximum();  
  if(theory6PtScaleNegerr->GetMaximum()>scaleerrmax)scaleerrmax=theory6PtScaleNegerr->GetMaximum();
  if(theoryPDFNegerr     ->GetMaximum()>scaleerrmax)scaleerrmax=theoryPDFNegerr     ->GetMaximum();
  if(theoryRelStaterr    ->GetMaximum()>scaleerrmax)scaleerrmax=theoryRelStaterr    ->GetMaximum();
  scaleerrmax*=10.;
  
  std::cout<<"scaleerrmin="<<scaleerrmin<<std::endl;
  std::cout<<"scaleerrmax="<<scaleerrmax<<std::endl;
  
  //--------------------------------------------
  //NEG SCALE ERRS
  TCanvas *plot_xsecNegErrs = new TCanvas("plot_xsecNegErrs", "plot xsec +Errors",1200,1000);
  plot_xsecNegErrs->cd()->SetLogx(1);
  plot_xsecNegErrs->cd()->SetLogy(1);
  
  theoryPDFNegerr->SetLineColor(kBlue);       
  //theory2PtScaleNegerr->SetLineColor(kRed);
  theory6PtScaleNegerr->SetLineColor(kMagenta);
  
  std::string _6ptnegerrs_oldtitle=(std::string)theory6PtScalePoserr->GetTitle();
  theory6PtScaleNegerr->SetTitle(("Theory "+orderstring+" #sigma_{Jet} (-) Errors;Jet p_{T}[GeV]; #||{(-) Errors} [%]").c_str());
  
  theory6PtScaleNegerr->SetMinimum(scaleerrmin);
  theory6PtScaleNegerr->SetMaximum(scaleerrmax);
  theory6PtScaleNegerr->GetXaxis()->SetNoExponent(true);
  theory6PtScaleNegerr->GetXaxis()->SetMoreLogLabels(true);
  
//  std::cout<<std::endl;
//  for(int i=1; i<=(theory6PtScaleNegerr->GetNbinsX());i++)
//    std::cout<<"theory6PtScaleNegerr->GetBinContent("<<i<<")="<<theory6PtScaleNegerr->GetBinContent(i)<<std::endl;
//  std::cout<<std::endl;
//  for(int i=1; i<=(theoryPDFNegerr->GetNbinsX());i++)
//    std::cout<<"theoryPDFNegerr->GetBinContent("<<i<<")="<<theoryPDFNegerr->GetBinContent(i)<<std::endl;

  TLegend* negerrsleg=new TLegend(0.1,0.7,0.4,0.9);
  negerrsleg->SetBorderSize(0.);
  negerrsleg->SetFillStyle(0);  
  negerrsleg->AddEntry( theory6PtScaleNegerr, "6 Pt. Scale Error (Asymm.)" , "l");
  //negerrsleg->AddEntry( theory2PtScaleNegerr, "2 Pt. Scale Error (Asymm.)" , "l");
  negerrsleg->AddEntry( theoryPDFNegerr     , "PDF Error (Asymm.)" , "l");
  negerrsleg->AddEntry( theoryRelStaterr    , "Stat. Unc. (Symm.)" , "l");
  
  plot_xsecNegErrs->cd();  

  theory6PtScaleNegerr->DrawClone("HIST");
  //theory2PtScaleNegerr->DrawClone("HIST SAME");
  theoryPDFNegerr     ->DrawClone("HIST SAME");
  theoryRelStaterr    ->DrawClone("HIST SAME");
  
  negerrsleg->Draw();
  
  theory6PtScaleNegerr->SetTitle(_6ptnegerrs_oldtitle.c_str());
  
  
  //assert(false);
  

  
  //--------------------------------------------
  //POSITIVE SCALE ERRS
  TCanvas *plot_xsecPosErrs = new TCanvas("plot_xsecPosErrs", "plot xsec +Errors",1200,1000);
  plot_xsecPosErrs->cd()->SetLogx(1);
  plot_xsecPosErrs->cd()->SetLogy(1);
  
  theoryPDFPoserr->SetLineColor(kBlue);       
  //theory2PtScalePoserr->SetLineColor(kRed);
  theory6PtScalePoserr->SetLineColor(kMagenta);
  
  std::string _6ptposerrs_oldtitle=(std::string)theory6PtScalePoserr->GetTitle();
  theory6PtScalePoserr->SetTitle(("Theory "+orderstring+" #sigma_{Jet} (+) Errors;Jet p_{T} [GeV]; (+) Errors [%]").c_str());
  
  theory6PtScalePoserr->SetMinimum(scaleerrmin);
  theory6PtScalePoserr->SetMaximum(scaleerrmax);
  theory6PtScalePoserr->GetXaxis()->SetNoExponent(true);
  theory6PtScalePoserr->GetXaxis()->SetMoreLogLabels(true);  


  //std::cout<<"theory6PtScalePoserr->GetMean()="<<theory6PtScalePoserr->GetMean()<<std::endl;
  //std::cout<<"theoryPDFPoserr     ->GetMean()="<<theoryPDFPoserr     ->GetMean()<<std::endl;
  //std::cout<<"theory6PtScalePoserr->GetBinContent(5)="<<theory6PtScalePoserr->GetBinContent(5)<<std::endl;
  //std::cout<<"theoryPDFPoserr     ->GetBinContent(5)="<<theoryPDFPoserr     ->GetBinContent(5)<<std::endl;
  //std::cout<<std::endl;
  for(int i=1; i<=(theory6PtScalePoserr->GetNbinsX());i++)
    std::cout<<"theory6PtScalePoserr->GetBinContent("<<i<<")="<<theory6PtScalePoserr->GetBinContent(i)<<std::endl;
  std::cout<<std::endl;
  for(int i=1; i<=(theoryPDFPoserr->GetNbinsX());i++)
    std::cout<<"theoryPDFPoserr->GetBinContent("<<i<<")="<<theoryPDFPoserr->GetBinContent(i)<<std::endl;
  
  TLegend* poserrsleg=new TLegend(0.1,0.7,0.4,0.9);
  poserrsleg->SetBorderSize(0.);
  poserrsleg->SetFillStyle(0);  
  poserrsleg->AddEntry( theory6PtScalePoserr, "6 Pt. Scale Error(Asymm.)" , "l");
  //poserrsleg->AddEntry( theory2PtScalePoserr, "2 Pt. Scale Error(Asymm.)" , "l");
  poserrsleg->AddEntry( theoryPDFPoserr     , "PDF Error (ASymm.)" , "l");  
  poserrsleg->AddEntry( theoryRelStaterr    , "Stat. Unc. (Symm.)" , "l");  
  plot_xsecPosErrs->cd();
  
  theory6PtScalePoserr->DrawClone("HIST");
  //theory2PtScalePoserr->DrawClone("HIST SAME");  
  theoryPDFPoserr     ->DrawClone("HIST SAME");
  theoryRelStaterr->DrawClone("HIST SAME");
  
  poserrsleg->Draw();

  theory6PtScalePoserr->SetTitle(_6ptposerrs_oldtitle.c_str());
  

 
  

  


            

  
  ///Cubic Spline Interpolations using Cross sections; has issues for p_T beyond first/last bin centers
  TSpline3 *spline3=NULL, *spline3_NP=NULL;
  TF1      *spline3_ext=NULL, *spline3_NP_ext=NULL;
  TCanvas *plot_splines=NULL;
  TLegend *leg_spline=NULL;

  TSpline3 *spline3_NPsys1=NULL, *spline3_NPsys2=NULL;   // NP syst
  TSpline3 *spline3_NPsysup=NULL, *spline3_NPsysdown=NULL; 

  TSpline3 *spline3_PDFsys1=NULL, *spline3_PDFsys2=NULL;//PDF syst
  TSpline3 *spline3_NP_PDFsys1=NULL, *spline3_NP_PDFsys2=NULL;
  TSpline3 *spline3_PDFsysup=NULL, *spline3_PDFsysdown=NULL;
  TSpline3 *spline3_NP_PDFsysup=NULL, *spline3_NP_PDFsysdown=NULL;
  if(useSplineWeights){
    std::cout<<"using spline weights!!!"<<std::endl;
    if(useSplineExt)std::cout<<"also adding spline weight extension!!!"<<std::endl;


    /// check spectra with splines
    std::cout<<"drawing Tspline and cross section on canvas"<<std::endl;
    plot_splines = new TCanvas("plot_splines", "plot splines",1200,800);
    plot_splines->cd()->SetLogx(1);
    plot_splines->cd()->SetLogy(1);
    plot_splines->cd();
    
    // NLO //
    std::cout<<"creating TSpline of cross section hist"<<std::endl;  
    spline3 = new TSpline3( theory_rebin);
    spline3->SetName( ( (std::string)theory_rebin->GetName() + "_spline3").c_str() );
    spline3->SetLineColor(kAzure);  

    
    if(doPDFsys12){
      spline3_PDFsys1 = new TSpline3( theory_PDFsys1_rebin);
      spline3_PDFsys1->SetName( ( (std::string)theory_PDFsys1_rebin->GetName() + "_spline3").c_str() );
      spline3_PDFsys1->SetLineColor(kAzure-8);  
      
      spline3_PDFsys2 = new TSpline3( theory_PDFsys2_rebin);
      spline3_PDFsys2->SetName( ( (std::string)theory_PDFsys2_rebin->GetName() + "_spline3").c_str() );
      spline3_PDFsys2->SetLineColor(kAzure-8);      
    }
    
    if(doPDFsysupdown){      
      spline3_PDFsysup = new TSpline3( theory_PDFsysup_rebin);
      spline3_PDFsysup->SetName( ( (std::string)theory_PDFsysup_rebin->GetName() + "_spline3").c_str() );
      spline3_PDFsysup->SetLineColor(kAzure-8);  
      
      spline3_PDFsysdown = new TSpline3( theory_PDFsysdown_rebin);
      spline3_PDFsysdown->SetName( ( (std::string)theory_PDFsysdown_rebin->GetName() + "_spline3").c_str() );
      spline3_PDFsysdown->SetLineColor(kAzure-8);      
    }

    double y2=calc_spline3ext_y2(theory_rebin);//this is also used to help set range for the plot ; so it's outside of the 'usesplineext' flag
    double x1=-1., x2=-1.;
    //double y2=theory_rebin->GetBinContent(theory_rebin->GetNbinsX());
    if(useSplineExt){
      //figure out starting values for parabolic extension of spline
      //x1=theory_rebin->GetBinCenter(theory_rebin->GetNbinsX());// bad idea, spline3's val at x1 is approx the val of the last bin's height. 
      if(fitType.find("v1")!=std::string::npos)x1=theory_rebin->GetBinLowEdge(theory_rebin->GetNbinsX()); //v1
      else if(fitType.find("v2")!=std::string::npos)x1=theory_rebin->GetBinCenter(theory_rebin->GetNbinsX()-1);//v2
      else if(fitType.find("v3")!=std::string::npos)x1=theory_rebin->GetBinLowEdge(theory_rebin->GetNbinsX()-1); //v3
      else assert(false);
      x2=theory_rebin->GetBinLowEdge(theory_rebin->GetNbinsX()) + theory_rebin->GetBinWidth(theory_rebin->GetNbinsX());                    
      
      //double xprime=theory_rebin->GetBinCenter(theory_rebin->GetNbinsX());
      //double xprime=theory_rebin->GetBinCenter(theory_rebin->GetNbinsX()) - 0.25*theory_rebin->GetBinWidth(theory_rebin->GetNbinsX());
      //double yprime=theory_rebin->GetBinContent(theory_rebin->GetNbinsX());
      
      //spline3_ext= new TF1( ((std::string)spline3->GetName()+"_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1, x2);    
      spline3_ext= new TF1( ((std::string)spline3->GetName()+"_tf1ext").c_str(), "[0]*x+[1]", x1, x2);    
      make_spline3_ext(  (TSpline3*)spline3, (TF1*)spline3_ext , (double)x1, (double)x2, (double)y2);
      //make_spline3_extv2(  (TSpline3*)spline3, (TF1*)spline3_ext , (double)x1, (double)x2, (double)y2, (double)xprime, (double)yprime);    
      spline3_ext->SetLineColor(kAzure-1);
    }

    
    // NLO+NP //      
    spline3_NP = new TSpline3( theory_NP_rebin);
    spline3_NP->SetName( ( (std::string)theory_NP_rebin->GetName() + "_spline3").c_str() );
    spline3_NP->SetLineColor(kAzure-8);  
    
    if(doNPsys12){
      spline3_NPsys1 = new TSpline3( theory_NPsys1_rebin);
      spline3_NPsys1->SetName( ( (std::string)theory_NPsys1_rebin->GetName() + "_spline3").c_str() );
      spline3_NPsys1->SetLineColor(kAzure-8);  

      spline3_NPsys2 = new TSpline3( theory_NPsys2_rebin);
      spline3_NPsys2->SetName( ( (std::string)theory_NPsys2_rebin->GetName() + "_spline3").c_str() );
      spline3_NPsys2->SetLineColor(kAzure-8);      
    }
    
    if(doNPsysupdown){
      spline3_NPsysup = new TSpline3( theory_NPsysup_rebin);
      spline3_NPsysup->SetName( ( (std::string)theory_NPsysup_rebin->GetName() + "_spline3").c_str() );
      spline3_NPsysup->SetLineColor(kAzure-8);  
      
      spline3_NPsysdown = new TSpline3( theory_NPsysdown_rebin);
      spline3_NPsysdown->SetName( ( (std::string)theory_NPsysdown_rebin->GetName() + "_spline3").c_str() );
      spline3_NPsysdown->SetLineColor(kAzure-8);      
    }

    if(doPDFsys12){
      spline3_NP_PDFsys1 = new TSpline3( theory_NP_PDFsys1_rebin);
      spline3_NP_PDFsys1->SetName( ( (std::string)theory_NP_PDFsys1_rebin->GetName() + "_spline3").c_str() );
      spline3_NP_PDFsys1->SetLineColor(kAzure-8);  

      spline3_NP_PDFsys2 = new TSpline3( theory_NP_PDFsys2_rebin);
      spline3_NP_PDFsys2->SetName( ( (std::string)theory_NP_PDFsys2_rebin->GetName() + "_spline3").c_str() );
      spline3_NP_PDFsys2->SetLineColor(kAzure-8);      
    }
    if(doPDFsysupdown){
      spline3_NP_PDFsysup = new TSpline3( theory_NP_PDFsysup_rebin);
      spline3_NP_PDFsysup->SetName( ( (std::string)theory_NP_PDFsysup_rebin->GetName() + "_spline3").c_str() );
      spline3_NP_PDFsysup->SetLineColor(kAzure-8);  
      
      spline3_NP_PDFsysdown = new TSpline3( theory_NP_PDFsysdown_rebin);
      spline3_NP_PDFsysdown->SetName( ( (std::string)theory_NP_PDFsysdown_rebin->GetName() + "_spline3").c_str() );
      spline3_NP_PDFsysdown->SetLineColor(kAzure-8);            
    }

    if(useSplineExt){
      //double NP_y2=theory_NP_rebin->GetBinContent(theory_NP_rebin->GetNbinsX());
      //double NP_y2=calc_spline3ext_y2(theory_NP_rebin);    
      
      spline3_NP_ext= new TF1( ((std::string)spline3_NP->GetName()+"_tf1ext").c_str(), "[0]*pow(x,2)+[1]*x+[2]", x1, x2);    
      //make_spline3_ext(  (TSpline3*)spline3_NP, (TF1*)spline3_NP_ext , (double)x1, (double)x2, (double)NP_y2);    
      spline3_NP_ext->SetLineColor(kAzure-9);
    }

    //draw
    theory_rebin->SetTitle(("Cubic Splines for "+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_rebin->SetAxisRange( y2/2., theory_rebin->GetMaximum()*2., "Y");
    
    theory_rebin->DrawClone("HIST E");
    spline3->Draw("SAME");    
    if(useSplineExt)spline3_ext->Draw("SAME");
    
    theory_NP_rebin->DrawClone("HIST E SAME");  
    spline3_NP->Draw("SAME");
    if(useSplineExt)spline3_NP_ext->Draw("SAME");
    
    leg_spline=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    leg_spline->AddEntry(theory_rebin ,    ("Weighted "+orderstring+" Jet Counts,"+absetarange_str           ).c_str() , "lp");
    leg_spline->AddEntry(spline3   ,       ("Cubic Spline "+orderstring+","+absetarange_str                  ).c_str() , "l");
    //
    if(useSplineExt)leg_spline->AddEntry(spline3_ext  ,    ("Parabolic Spline Extension "+orderstring+","+absetarange_str    ).c_str() , "l");
    //
    leg_spline->AddEntry(theory_NP_rebin , ("Weighted NP+"+orderstring+" Jet Counts,"+absetarange_str        ).c_str() , "lp");
    leg_spline->AddEntry(spline3_NP   ,    ("Cubic Spline NP+"+orderstring+","+absetarange_str               ).c_str() , "l");
    if(useSplineExt)leg_spline->AddEntry(spline3_NP_ext  , ("Parabolic Spline Extension NP+"+orderstring+","+absetarange_str ).c_str() , "l");
    leg_spline->Draw();    
  }

  // for fits
  const int startbin=1, endbin=theory_rebin->GetNbinsX();//NP and nonNP hist should always have the same binning
  const float xlo=theory_rebin->GetBinLowEdge(startbin);
  const float xhi=theory_rebin->GetBinLowEdge(endbin) + theory_rebin->GetBinWidth(endbin);  
  
  
  ///Log fit: A(B/pt)^[C+ D Log(pt/B) + E(Log(pt/B))^2] using Cross sections; currently no known issues
  TF1 *logFit  =NULL, *modLogFit  =NULL, *modLog2Fit  =NULL, *modLog3Fit  =NULL, *modLog4Fit  =NULL,*modLog5Fit  =NULL;
  TF1* logFit_4Toy=NULL;

  TF1 *logFit_NP=NULL, *modLogFit_NP=NULL, *modLog2Fit_NP=NULL, *modLog3Fit_NP=NULL, *modLog4Fit_NP=NULL,*modLog5Fit_NP=NULL;
  TF1* logFit_NP_4Toy=NULL;
  
  TCanvas *plot_logfits=NULL,*plot_logfits_NP=NULL;
  TLegend *leg_logfits=NULL,*leg_logfits_NP=NULL;  
  //  if(useFitWeights){
  if(useFitWeights && useModLogFit){
    std::cout<<"using modlogfit weights!!!"<<std::endl;

    // --- NLO --- //
    /// check spectra with fits
    plot_logfits = new TCanvas("plot_logfits", "plot logfits",1200,800);
    plot_logfits->cd()->SetLogx(1);
    plot_logfits->cd()->SetLogy(1);
    plot_logfits->cd();
    
    std::cout<<"fitting hist theory_rebin"<<std::endl;
    logFit    =new TF1("logfit    ","[0] * ( TMath::Power( ( [1]/x ) ,   [2] ) )",xlo,xhi);
    modLogFit =new TF1("modlogfit ","[0] * ( TMath::Power( ( [1]/x ) ,   [2] + [3]*TMath::Log10( x/[1] ) ) )",xlo,xhi);
    modLog2Fit=new TF1("modlog2fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  ) ) )",xlo,xhi);
    modLog3Fit=new TF1("modlog3fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 )  ) ) )",xlo,xhi);
    modLog4Fit=new TF1("modlog4fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) )               ) )",xlo,xhi);
    modLog5Fit=new TF1("modlog5fit","[0] * ( TMath::Power( ( [1]/x ) , ( [2] + [3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )  + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) + [7]*TMath::Power( TMath::Log10( x/[1] ), 5 ) )               ) )",xlo,xhi);

    
    makeModLogFits(theory_rebin, logFit, modLogFit, modLog2Fit, modLog3Fit, modLog4Fit, modLog5Fit);
    
    if(etabin==0)logFit_4Toy=modLog5Fit;
    else if(etabin==1)logFit_4Toy=modLog2Fit;
    else if(etabin==2)logFit_4Toy=modLog4Fit;
    else if(etabin==3)logFit_4Toy=modLog2Fit;
    else if(etabin==4)logFit_4Toy=modLog4Fit;
    else if(etabin==5)logFit_4Toy=modLog4Fit;
    


    //draw
    std::cout<<"drawing fits and cross section on canvas"<<std::endl;
    theory_rebin->SetTitle(("Modified Logarithm Fits for "+orderstring+", "+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_rebin->DrawClone("HIST E");    
    
    logFit->Draw("SAME");
    modLogFit->Draw("SAME");
    modLog2Fit->Draw("SAME");
    modLog3Fit->Draw("SAME");    
    modLog4Fit->Draw("SAME");    
    modLog5Fit->Draw("SAME");    
    
    leg_logfits=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    
    leg_logfits->AddEntry(theory_rebin , ("Weighted "+orderstring+" Jet Counts for "+absetarange_str).c_str() , "lp");//change desc TO DO
    leg_logfits->AddEntry(logFit      , "Log Fit" , "l");//change desc TO DO
    leg_logfits->AddEntry(modLogFit   , "Mod. Log Fit" , "l");
    leg_logfits->AddEntry(modLog2Fit   , "2nd Order Mod. Log Fit" , "l");    
    leg_logfits->AddEntry(modLog3Fit   , "3rd Order Mod. Log Fit" , "l");    
    leg_logfits->AddEntry(modLog4Fit   , "4th Order Mod. Log Fit" , "l");    
    leg_logfits->AddEntry(modLog5Fit   , "5th Order Mod. Log Fit" , "l");    
    leg_logfits->Draw();    
    
    
    // --- NLO+NP --- //
    /// check spectra with fits    
    plot_logfits_NP = new TCanvas("plot_logfits_NP", "plot logfits NP",1200,800);
    plot_logfits_NP->cd()->SetLogx(1);
    plot_logfits_NP->cd()->SetLogy(1);
    plot_logfits_NP->cd();
    
    std::cout<<"fitting hist theory_NP_rebin"<<std::endl;
    logFit_NP    =new TF1("logfit_NP    ","[0] * ( TMath::Power( ( [1]/x ) , ( [2]                                                                                                                         ) ) )",xlo,xhi);
    modLogFit_NP =new TF1("modlogfit_NP ","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] )                                                                                               ) ) )",xlo,xhi);
    modLog2Fit_NP=new TF1("modlog2fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 )                                                ) ) )",xlo,xhi);
    modLog3Fit_NP=new TF1("modlog3fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) ) ) )",xlo,xhi);
    modLog4Fit_NP=new TF1("modlog4fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) ) ) )",xlo,xhi);
    modLog5Fit_NP=new TF1("modlog5fit_NP","[0] * ( TMath::Power( ( [1]/x ) , ( [2]+[3]*TMath::Log10( x/[1] ) + [4]*TMath::Power( TMath::Log10( x/[1] ), 2 ) + [5]*TMath::Power( TMath::Log10( x/[1] ), 3 ) + [6]*TMath::Power( TMath::Log10( x/[1] ), 4 ) + [7]*TMath::Power( TMath::Log10( x/[1] ), 5 ) ) ) )",xlo,xhi); 
    
    
    makeModLogFits(theory_NP_rebin, logFit_NP, modLogFit_NP, modLog2Fit_NP, modLog3Fit_NP, modLog4Fit_NP, modLog5Fit_NP);

    if(etabin==0)logFit_NP_4Toy=modLog5Fit_NP;
    if(etabin==1)logFit_NP_4Toy=modLog5Fit_NP;
    if(etabin==2)logFit_NP_4Toy=modLog4Fit_NP;
    if(etabin==3)logFit_NP_4Toy=modLog2Fit_NP;
    if(etabin==4)logFit_NP_4Toy=modLog4Fit_NP;
    if(etabin==5)logFit_NP_4Toy=modLog4Fit_NP;

    
    
    //draw
    std::cout<<"drawing NP fits and cross section on canvas"<<std::endl;
    theory_NP_rebin->SetTitle( ("Modified Logarithm Fits for "+orderstring+"+NP,"+absetarange_str+" ;Jet p_{T};Smear Weight").c_str());        
    theory_NP_rebin->DrawClone("HIST E");        
    
    logFit_NP->Draw("SAME");
    modLogFit_NP->Draw("SAME");
    modLog2Fit_NP->Draw("SAME");
    modLog3Fit_NP->Draw("SAME");    
    modLog4Fit_NP->Draw("SAME");    
    modLog5Fit_NP->Draw("SAME");    
    
    leg_logfits_NP=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    leg_logfits_NP->AddEntry(theory_NP_rebin , ("Weighted "+orderstring+"+NP Jet Counts for "+absetarange_str).c_str() , "lp");//change desc TO DO
    //leg_logfits_NP->AddEntry(theory_NP , "Weighted NLO+NP Jet Counts for"+absetarange_str , "lp");//change desc TO DO
    leg_logfits_NP->AddEntry(logFit_NP      , "Log Fit" , "l");//change desc TO DO
    leg_logfits_NP->AddEntry(modLogFit_NP   , "Mod. Log Fit" , "l");
    leg_logfits_NP->AddEntry(modLog2Fit_NP   , "2nd Order Mod. Log Fit" , "l");    
    leg_logfits_NP->AddEntry(modLog3Fit_NP   , "3rd Order Mod. Log Fit" , "l");    
    leg_logfits_NP->AddEntry(modLog4Fit_NP   , "4th Order Mod. Log Fit" , "l");    
    leg_logfits_NP->AddEntry(modLog5Fit_NP   , "5th Order Mod. Log Fit" , "l");    
    leg_logfits_NP->Draw();    

  }


  ///7 TeV fit: two variations using Cross sections, each w/ 3 free parameters min. note; sqrt(s) should be in GeV
  ///cite: CERN-PH-EP/2011-053, 2018/10/22, CMS-QCD-10-011, 
  ///    : "Measurement of the Inclusive Jet Cross Section in pp Collisions at sqrt(s)=7 TeV"
  ///    : arixiv:1106.0208v1 [hep-ex] 1 Jun 2011
  ///    : Submitted to Physical Review Letters
  /// A) ddxsec[p_T] = D * [(p_T/GeV)^A] * [(1.-(2. *p_T/sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
  /// B) ddxsec[p_T] = [integral(ddxsec[p_T])] * [(p_T/GeV)^A] * [(1 - (2 * p_T / sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
  TF1 *_7TeVFitA  =NULL, *_7TeVFitB  =NULL;
  TF1 *_7TeVFitA_NP=NULL, *_7TeVFitB_NP=NULL;
  
  TCanvas *plot_7tevfits=NULL,*plot_7tevfits_NP=NULL;
  TLegend  *leg_7tevfits=NULL, *leg_7tevfits_NP=NULL;  
  
  //const float xlo_7TeV=theory_rebin->GetBinLowEdge(startbin);//trying out for the 7TeV style fit
  //const float xhi_7TeV=theory_rebin->GetBinCenter(endbin);// + theory_rebin->GetBinWidth(endbin);  
  
  

  if(useFitWeights && use7TeVFit){
    std::cout<<"using 7TeVfit weights!!!"<<std::endl;
    
    std::cout<<"making 7tev fit strings"<<std::endl;
    float _coshymin=2.*TMath::CosH(absetabins[etabin]);//putting the 2. in here for simplicity
    std::string _coshymin_str=std::to_string(_coshymin);
    std::cout<<"_coshymin_str="<<_coshymin_str<<std::endl;
    
    float norm=theory_rebin->Integral("width");
    std::cout<<"norm="<<norm<<std::endl;

    std::string norm_str=std::to_string(norm);
    std::cout<<"norm_str="<<norm_str<<std::endl;

    float norm_NP=theory_NP_rebin->Integral("width");
    std::cout<<"norm_NP="<<norm_NP<<std::endl;
    
    std::string norm_NP_str=std::to_string(norm_NP);
    std::cout<<"norm_NP_str="<<norm_NP_str<<std::endl;    
    
    /// A) ddxsec[p_T] = D * [(p_T/GeV)^A] * [(1.-(2.*p_T/sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
    //    std::string _7tevfit_str="(TMath::Power(x,-1.*[0]))*(TMath::Power((1.-(2.*x*"+_coshymin_str+"*/(5.02*1000.) )),[1]))* (TMath::Exp(-1.*([2]/x)))";
    std::string _7tevfit_str="(TMath::Power(x,-1.*[0]))*(TMath::Power((1.-("+_coshymin_str+"*x/(5.02*1000.) )),[1]))* (TMath::Exp(-1.*([2]/x)))";
    std::cout<<"_7tevfit_str="<<_7tevfit_str<<std::endl;
    
    std::string _7tevfitA_str="([3])* "+_7tevfit_str;    //this string is same for nonNP and NP
    std::cout<<"_7tevfitA_str="<<_7tevfitA_str<<std::endl;
    
    // --- NLO --- //    
    /// check spectra with fits
    plot_7tevfits = new TCanvas("plot_7tevfits", "plot 7tevfits",1200,800);
    plot_7tevfits->cd()->SetLogx(1);
    plot_7tevfits->cd()->SetLogy(1);
    plot_7tevfits->cd();
    


  /// B) ddxsec[p_T] = [integral(ddxsec[p_T])] * [(p_T/GeV)^A] * [(1 - (2 * p_T / sqrt(s)]*cosh(ymin))^B] * [exp(-1*C/(p_T/GeV))]
    std::string _7tevfitB_str="("+norm_str+")* "+_7tevfit_str;
    std::cout<<"_7tevfitB_str="<<_7tevfitB_str<<std::endl;    
    
    std::cout<<std::endl<<"making _7TeVFitA TF1"<<std::endl;
    _7TeVFitA =new TF1("_7tevfitA", (_7tevfitA_str     ).c_str() ,xlo,xhi);
    
    std::cout<<std::endl<<"making _7TeVFitB TF1"<<std::endl;
    _7TeVFitB =new TF1("_7tevfitB" , (_7tevfitB_str).c_str(),xlo,xhi);
    
    make7TeVFits(theory_rebin,  _7TeVFitA, _7TeVFitB, norm);

    //draw
    std::cout<<"drawing fits and cross section on canvas"<<std::endl;
    theory_rebin->SetTitle(("7 TeV Fits for "+orderstring+", "+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_rebin->DrawClone("HIST E");        
    _7TeVFitA->Draw("SAME");
    _7TeVFitB->Draw("SAME");
    
    leg_7tevfits=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");    
    leg_7tevfits->AddEntry(theory_rebin , ("Weighted "+orderstring+" Jet Counts for"+absetarange_str).c_str() , "lp");//change desc TO DO
    leg_7tevfits->AddEntry(_7TeVFitB   , ("4 Param. 7 TeV Fit for "+orderstring+"").c_str() , "l");
    leg_7tevfits->AddEntry(_7TeVFitA   , ("3 Param. 7 TeV Fit for "+orderstring+"").c_str() , "l");    
    leg_7tevfits->Draw();    
    
    
    // --- NP+NLO --- //
    /// check spectra with fits
    plot_7tevfits_NP = new TCanvas("plot_7tevfits_NP", "plot 7tevfits NP",1200,800);
    plot_7tevfits_NP->cd()->SetLogx(1);
    plot_7tevfits_NP->cd()->SetLogy(1);
    plot_7tevfits_NP->cd();
    

    std::string _7tevfitB_NP_str="("+norm_NP_str+")* "+_7tevfit_str;
    std::cout<<"_7tevfitB_NP_str="<<_7tevfitB_NP_str<<std::endl;        
    
    std::cout<<"fitting hist theory_NP_rebin"<<std::endl;
    _7TeVFitA_NP =new TF1("_7tevfitA_NP",(_7tevfitA_str       ).c_str(),xlo,xhi);
    _7TeVFitB_NP=new TF1( "_7tevfitB_NP" ,(_7tevfitB_NP_str).c_str(),xlo,xhi);
    make7TeVFits(theory_NP_rebin, _7TeVFitA_NP, _7TeVFitB_NP, norm_NP);
    
    //draw    
    std::cout<<"drawing NP fits and cross section on canvas"<<std::endl;
    theory_NP_rebin->SetTitle( ("7 TeV Fits for NP+"+orderstring+","+absetarange_str+";Jet p_{T};Smear Weight").c_str());    
    theory_NP_rebin->DrawClone("HIST E");    
    
    _7TeVFitA_NP->Draw("SAME");
    _7TeVFitB_NP->Draw("SAME");
    
    leg_7tevfits_NP=new TLegend(0.65, 0.70, 0.9, 0.9, NULL,"BRNDC");
    
    leg_7tevfits_NP->AddEntry(theory_NP_rebin , ("Weighted NP+"+orderstring+" Jet Counts for"+absetarange_str).c_str() , "lp");//change desc TO DO
    leg_7tevfits_NP->AddEntry(_7TeVFitA_NP   ,  ("3 Param. 7 TeV Fit for NP+"+orderstring).c_str() , "l");
    leg_7tevfits_NP->AddEntry(_7TeVFitB_NP   ,  ("4 Param. 7 TeV Fit for NP+"+orderstring).c_str() , "l");    
    leg_7tevfits_NP->Draw();    
  }
  



  
  //-----------------------------------------------------------------------------------------//
  //////////////////////  START production of Smeared NLO spectra  ////////////////////////////
  //-----------------------------------------------------------------------------------------//
  std::cout<<"creating TH1 for toy NLO spectra generation, RooUnfoldResponse class, etc."<<std::endl<<std::endl;
  TH1D *theory_rnd       = new TH1D("theory_rnd","theory_rnd",             n_thybins    , thybins);   
  TH1D *smeared_rnd      = new TH1D("smeared_rnd","smeared_rnd",           n_thybins, thybins);   
  TH1D *smeared_rnd_test = new TH1D("smeared_rnd_test","smeared_rnd_test", n_thybins, thybins);     
  TH2D* response_th2=new TH2D("response_th2","response_th2",
			      (Int_t)n_thybins, (Double_t*)thybins,
			      (Int_t)n_thybins, (Double_t*)thybins);  

  //--------------------------------------------------------------------------------------------------------//
  std::cout<<"making rnd toy MC hists"<<std::endl;        
  if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
    makeToySpectra( (TH1D*)theory_rebin, 
		    (TSpline3*)spline3, 
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd, 
		    smeared_rnd, 
		    smeared_rnd_test, 
		    (TH2D*)response_th2,
		    (TF1*) spline3_ext); }//		    NULL); }
  else if(useFitWeights&&useModLogFit){//looks like modlog4Fit wins!
    makeToySpectra( (TH1D*)theory_rebin, 
		    (TF1*)logFit_4Toy, //		    (TF1*)modLog4Fit, // (TF1*)modLog3Fit, // 		    
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd, 
		    smeared_rnd, 
		    smeared_rnd_test, 
		    (TH2D*)response_th2); }  
  else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    makeToySpectra( (TH1D*)theory_rebin, 
		    (TF1*)_7TeVFitA,// (TF1*)_7TeVFitB, //
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd, 
		    smeared_rnd, 
		    smeared_rnd_test, 
		    (TH2D*)response_th2); }  

  TH1D* theory_rnd_RelStaterr=(TH1D*)theory_rnd->Clone("theory_rnd_relStaterr");
  for(int i=1; i<=theory_rnd->GetNbinsX();i++){
    float val=theory_rnd->GetBinContent(i);
    float valerr=theory_rnd->GetBinError(i);
    theory_rnd_RelStaterr->SetBinContent(i, valerr/val*100.);//units of %
    theory_rnd_RelStaterr->SetBinError(i, 0.);
  }
  theory_rnd_RelStaterr->SetTitle("Rel Stat Unc. for Toy Theory");

  TH1D* smeared_rnd_RelStaterr=(TH1D*)smeared_rnd->Clone("smeared_rnd_relStaterr");
  for(int i=1; i<=smeared_rnd->GetNbinsX();i++){
    float val=smeared_rnd->GetBinContent(i);
    float valerr=smeared_rnd->GetBinError(i);
    smeared_rnd_RelStaterr->SetBinContent(i, valerr/val*100.);//units of %
    smeared_rnd_RelStaterr->SetBinError(i, 0.);
  }
  smeared_rnd_RelStaterr->SetTitle("Rel Stat Unc. for Smeared Toy Theory");
  
  divideBinWidth(theory_rnd);  //for normalization only.
  divideBinWidth(smeared_rnd); 
  divideBinWidth(smeared_rnd_test); 
  
  double norm_rnd=theory_rebin->Integral()/theory_rnd->Integral();
  response_th2->Scale(   norm_rnd);    
  theory_rnd->Scale( norm_rnd);
  smeared_rnd->Scale(norm_rnd);  
  smeared_rnd_test->Scale(norm_rnd);  
  
  //plotting smeared+theory toy mc xsection
  theory_rnd->SetLineColor(kRed); 
  smeared_rnd->SetLineColor(kBlue);    
  smeared_rnd_test->SetLineColor(kBlue-5);    


  //-----------------------------------------------------------------------------------------//
  //////////////////////  END production of Smeared NLO spectra  //////////////////////////////
  //-----------------------------------------------------------------------------------------//  

  //----------------------------START TOY NLO COSMETICS/PLOTS--------------------------------------------//  
  Float_t maxy=theory_rnd->GetMaximum();//    std::cout<<"maxy="<<maxy<<std::endl;
  maxy*=2.;  //std::cout<<"maxy="<<maxy<<std::endl;
  Float_t miny=smeared_rnd->GetMinimum(); //  std::cout<<"miny="<<miny<<std::endl;
  miny/=2.;// std::cout<<"miny="<<miny<<std::endl;
  theory_rnd->SetAxisRange(miny,maxy,"Y");  
  theory_rnd->SetTitle( ((std::string)("Toy MC Incl. "+orderstring+" #sigma;;"+ddxsec_yax)).c_str() );
  
  
  ///// check plots True/Smeared 
  TCanvas *plot_true_smeared_rat = new TCanvas("plot_true_smeared_rat", "plot true smeared rat",1600,1600);  
  plot_true_smeared_rat->cd();
  
  TPad* overlay_pad=new TPad("overlay_pad","pad for spectra",0.,0.3,1.0,1.0);
  overlay_pad->SetLogx(1);
  overlay_pad->SetLogy(1);
  overlay_pad->SetBottomMargin(0);
  overlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor
  
  TPad* ratio_pad=new TPad("ratio_pad","pad for ratio",0.,0.,1.0,0.3);
  ratio_pad->SetLogx(1);
  ratio_pad->SetLogy(0);
  ratio_pad->SetTopMargin(0);
  ratio_pad->SetBottomMargin(0.20);
  ratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor

  // spectra pad
  overlay_pad->cd();  
  theory_rnd->DrawClone("HIST E");  
  smeared_rnd->DrawClone("HIST E SAME");
  smeared_rnd_test->DrawClone("HIST E SAME");
  
  TLegend* leg_smear=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  leg_smear->SetBorderSize(0.);
  leg_smear->SetFillStyle(0.);
  leg_smear->AddEntry(theory_rnd  ,      ("Truth Incl. "+orderstring+" #sigma"            ).c_str(), "lp" );
  leg_smear->AddEntry(smeared_rnd ,      ("Smeared Incl. "+orderstring+" #sigma"          ).c_str(), "lp" );
  leg_smear->AddEntry(smeared_rnd_test , ("Opp Side Smeared Incl. "+orderstring+" #sigma" ).c_str(), "lp" );
  leg_smear->Draw();  

  // ratio pad
  ratio_pad->cd();
  
  //make ratios
  TH1D *true_smeared_rat=(TH1D*)theory_rnd->Clone("true_smeared_rat");   
  if(printBaseDebug)true_smeared_rat->Print("base");
  true_smeared_rat->Divide(true_smeared_rat,smeared_rnd,1.,1.,"B");   

  TH1D *true_smeared_test_rat=(TH1D*)theory_rnd->Clone("true_smeared_test_rat");   
  if(printBaseDebug)true_smeared_test_rat->Print("base");  
  true_smeared_test_rat->Divide(true_smeared_test_rat,smeared_rnd_test,1.,1.,"B");   

  //set stuff and draw
  true_smeared_rat->SetTitle("; Jet p_{T} [GeV]; True/Smeared");
  true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  true_smeared_rat->GetXaxis()->SetNoExponent(true);
  true_smeared_rat->GetXaxis()->SetLabelSize(0.05);
  true_smeared_rat->GetXaxis()->SetTitleSize(0.055);

  true_smeared_rat->GetYaxis()->SetLabelSize(0.05);
  true_smeared_rat->GetYaxis()->SetTitleSize(0.05);
  true_smeared_rat->GetYaxis()->SetTitleOffset(0.7);    
  true_smeared_rat->GetYaxis()->CenterTitle(true);    
  true_smeared_rat->SetAxisRange(0.85,1.15,"Y");  

  true_smeared_rat->SetLineColor(smeared_rnd->GetLineColor());
  
  //true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  //true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  //true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  //true_smeared_test_rat->SetAxisRange(0.8,1.2,"Y");
  true_smeared_test_rat->SetLineColor(smeared_rnd_test->GetLineColor());
  
  true_smeared_rat->DrawClone("HIST E");        
  true_smeared_test_rat->DrawClone("HIST E SAME");        
  
  lineatone->SetX2(thybins[n_thybins]);
  lineatone->DrawClone();


  // response matrix plot
  TCanvas* plot_response_th2=new TCanvas("plot_response_th2","plot response th2",1200, 1000);
  plot_response_th2->cd()->SetLogx(1);
  plot_response_th2->cd()->SetLogy(1);
  plot_response_th2->cd()->SetLogz(1);
  plot_response_th2->cd();
  
  response_th2->SetTitle(("Toy MC Response Matrix from JER Smeared "+orderstring+" #sigma_{jet}").c_str());  
  response_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_th2->GetXaxis()->SetNoExponent(true);
  response_th2->GetXaxis()->SetMoreLogLabels(true);
  response_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_th2->GetYaxis()->SetNoExponent(true);
  response_th2->GetYaxis()->SetMoreLogLabels(true);  
  setTH2_ZAxisRange(response_th2);  

  response_th2->Draw("COLZ");
  





  //--------------------------------------------------------------------------------------------------------//

  //-----------------------------------------------------------------------------------------//
  //////////////////////  BEGIN production of Smeared NLO spectra for JER sys  ////////////////
  //-----------------------------------------------------------------------------------------//  

  TH1D * theory_rnd_JERsysup       = NULL  , * theory_rnd_JERsysdown       = NULL;    
  TH1D *smeared_rnd_JERsysup      = NULL   , *smeared_rnd_JERsysdown      = NULL;
  TH1D *smeared_rnd_JERsysup_test = NULL   , *smeared_rnd_JERsysdown_test = NULL;
  TH2D *   response_JERsysup_th2     = NULL, *   response_JERsysdown_th2     = NULL;

  if(doJERsys){
    //shift JER up
    theory_rnd_JERsysup       = new TH1D("theory_rnd_JERsysup","theory_rnd_JERsysup",             n_thybins    , thybins);   
    smeared_rnd_JERsysup      = new TH1D("smeared_rnd_JERsysup","smeared_rnd_JERsysup",           n_thybins, thybins);   
    smeared_rnd_JERsysup_test = new TH1D("smeared_rnd_JERsysup_test","smeared_rnd_JERsysup_test", n_thybins, thybins);       
    response_JERsysup_th2=new TH2D("response_JERsysup_th2","response_JERsysup_th2",
				(Int_t)n_thybins, (Double_t*)thybins,
				(Int_t)n_thybins, (Double_t*)thybins);      
    
    std::cout<<"making rnd_JERsysup toy MC hists"<<std::endl;
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_rebin, 
		      (TSpline3*)spline3, 
		      (TF1*)fJER_sysup, 
		      nEvents, 
		      theory_rnd_JERsysup, 
		      smeared_rnd_JERsysup, 
		      smeared_rnd_JERsysup_test, 
		      (TH2D*)response_JERsysup_th2,
		      (TF1*) spline3_ext); }//		    NULL); }
    //    else if(useFitWeights&&useModLogFit){//looks like modlog4Fit wins!
    //      makeToySpectra( (TH1D*)theory_rebin, 
    //		      (TF1*)logFit_4Toy, //		    (TF1*)modLog4Fit, // (TF1*)modLog3Fit, // 		    
    //		      (TF1*)fJER_sysup, 
    //		      nEvents, 
    //		      theory_rnd_JERsysup, 
    //		      smeared_rnd_JERsysup, 
    //		      smeared_rnd_JERsysup_test, 
    //		      (TH2D*)response_JERsysup_th2); }  
    //    else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    //      makeToySpectra( (TH1D*)theory_rebin, 
    //		      (TF1*)_7TeVFitA,// (TF1*)_7TeVFitB, //
    //		      (TF1*)fJER_sysup, 
    //		      nEvents, 
    //		      theory_rnd_JERsysup, 
    //		      smeared_rnd_JERsysup, 
    //		      smeared_rnd_JERsysup_test, 
    //		      (TH2D*)response_JERsysup_th2); }  

    divideBinWidth(theory_rnd_JERsysup);  //for normalization only.
    divideBinWidth(smeared_rnd_JERsysup); 
    divideBinWidth(smeared_rnd_JERsysup_test); 
    
    double norm_rnd_JERsysup=theory_rebin->Integral()/theory_rnd_JERsysup->Integral();
    response_JERsysup_th2->Scale(   norm_rnd_JERsysup);      
    theory_rnd_JERsysup->Scale( norm_rnd_JERsysup);
    smeared_rnd_JERsysup->Scale(norm_rnd_JERsysup);  
    smeared_rnd_JERsysup_test->Scale(norm_rnd_JERsysup);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_JERsysup->SetLineColor(kRed); 
    smeared_rnd_JERsysup->SetLineColor(kBlue);    
    smeared_rnd_JERsysup_test->SetLineColor(kBlue-5);    



    //shift JER down
    theory_rnd_JERsysdown       = new TH1D("theory_rnd_JERsysdown","theory_rnd_JERsysdown",             n_thybins    , thybins);   
    smeared_rnd_JERsysdown      = new TH1D("smeared_rnd_JERsysdown","smeared_rnd_JERsysdown",           n_thybins, thybins);   
    smeared_rnd_JERsysdown_test = new TH1D("smeared_rnd_JERsysdown_test","smeared_rnd_JERsysdown_test", n_thybins, thybins);       
    response_JERsysdown_th2=new TH2D("response_JERsysdown_th2","response_JERsysdown_th2",
				(Int_t)n_thybins, (Double_t*)thybins,
				(Int_t)n_thybins, (Double_t*)thybins);      

    std::cout<<"making rnd_JERsysdown toy MC hists"<<std::endl;
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_rebin, 
		      (TSpline3*)spline3, 
		      (TF1*)fJER_sysdown, 
		      nEvents, 
		      theory_rnd_JERsysdown, 
		      smeared_rnd_JERsysdown, 
		      smeared_rnd_JERsysdown_test, 
		      (TH2D*)response_JERsysdown_th2,
		    (TF1*) spline3_ext); }//		    NULL); }
//else if(useFitWeights&&useModLogFit){//looks like modlog4Fit wins!
//  makeToySpectra( (TH1D*)theory_rebin, 
//		      (TF1*)logFit_4Toy, //		    (TF1*)modLog4Fit, // (TF1*)modLog3Fit, // 		    
//		      (TF1*)fJER_sysdown, 
//		      nEvents, 
//		      theory_rnd_JERsysdown, 
//		      smeared_rnd_JERsysdown, 
//		      smeared_rnd_JERsysdown_test, 
//		      (TH2D*)response_JERsysdown_th2); }  
//else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
//  makeToySpectra( (TH1D*)theory_rebin, 
//		      (TF1*)_7TeVFitA,// (TF1*)_7TeVFitB, //
//		      (TF1*)fJER_sysdown, 
//		      nEvents, 
//		      theory_rnd_JERsysdown, 
//		      smeared_rnd_JERsysdown, 
//		      smeared_rnd_JERsysdown_test, 
//		      (TH2D*)response_JERsysdown_th2); }  

    divideBinWidth(theory_rnd_JERsysdown);  //for normalization only.
    divideBinWidth(smeared_rnd_JERsysdown); 
    divideBinWidth(smeared_rnd_JERsysdown_test); 
    
    double norm_rnd_JERsysdown=theory_rebin->Integral()/theory_rnd_JERsysdown->Integral();
    response_JERsysdown_th2->Scale(   norm_rnd_JERsysdown);      
    theory_rnd_JERsysdown->Scale( norm_rnd_JERsysdown);
    smeared_rnd_JERsysdown->Scale(norm_rnd_JERsysdown);  
    smeared_rnd_JERsysdown_test->Scale(norm_rnd_JERsysdown);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_JERsysdown->SetLineColor(kRed); 
    smeared_rnd_JERsysdown->SetLineColor(kBlue);    
    smeared_rnd_JERsysdown_test->SetLineColor(kBlue-5);    
  }

  //-----------------------------------------------------------------------------------------//
  //////////////////////  END production of Smeared NLO spectra for JER sys////////////////////
  //-----------------------------------------------------------------------------------------//  

  //--------------------------------------------------------------------------------------------------------//

  //-----------------------------------------------------------------------------------------//
  //////////////////////  BEGIN production of Smeared NLO spectra for PDF sys  ////////////////
  //-----------------------------------------------------------------------------------------//  


  TH1D * theory_rnd_PDFsys1       = NULL  , * theory_rnd_PDFsys2       = NULL;    
  TH1D *smeared_rnd_PDFsys1      = NULL   , *smeared_rnd_PDFsys2      = NULL;
  TH1D *smeared_rnd_PDFsys1_test = NULL   , *smeared_rnd_PDFsys2_test = NULL;
  TH2D *   response_PDFsys1_th2     = NULL, *   response_PDFsys2_th2     = NULL;

  if(doPDFsys12){

    // PDF alt choice 1 ---------------
    std::cout<<"making rnd_PDFsys1 toy MC hists"<<std::endl;        
    theory_rnd_PDFsys1       = new TH1D("theory_rnd_PDFsys1","theory_rnd_PDFsys1",             n_thybins    , thybins);   
    smeared_rnd_PDFsys1      = new TH1D("smeared_rnd_PDFsys1","smeared_rnd_PDFsys1",           n_thybins, thybins);   
    smeared_rnd_PDFsys1_test = new TH1D("smeared_rnd_PDFsys1_test","smeared_rnd_PDFsys1_test", n_thybins, thybins);       
    response_PDFsys1_th2=new TH2D("response_PDFsys1_th2","response_PDFsys1_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				   (Int_t)n_thybins, (Double_t*)thybins);      

    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_PDFsys1_rebin, 
		      (TSpline3*)spline3_PDFsys1, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_PDFsys1, 
		      smeared_rnd_PDFsys1, 
		      smeared_rnd_PDFsys1_test, 
		      (TH2D*)response_PDFsys1_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_PDFsys1);  //for normalization only.
    divideBinWidth(smeared_rnd_PDFsys1); 
    divideBinWidth(smeared_rnd_PDFsys1_test); 
    
    double norm_rnd_PDFsys1=theory_PDFsys1_rebin->Integral()/theory_rnd_PDFsys1->Integral();
    response_PDFsys1_th2->Scale(   norm_rnd_PDFsys1);      
    theory_rnd_PDFsys1->Scale( norm_rnd_PDFsys1);
    smeared_rnd_PDFsys1->Scale(norm_rnd_PDFsys1);  
    smeared_rnd_PDFsys1_test->Scale(norm_rnd_PDFsys1);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_PDFsys1->SetLineColor(kRed); 
    smeared_rnd_PDFsys1->SetLineColor(kBlue);    
    smeared_rnd_PDFsys1_test->SetLineColor(kBlue-5);    
    
    std::cout<<"done making rnd_PDFsys1 toy MC hists"<<std::endl;        
    
    // PDF alt choice 2 ---------------
    std::cout<<"making rnd_PDFsys2 toy MC hists"<<std::endl;        
    theory_rnd_PDFsys2       = new TH1D("theory_rnd_PDFsys2","theory_rnd_PDFsys2",             n_thybins    , thybins);   
    smeared_rnd_PDFsys2      = new TH1D("smeared_rnd_PDFsys2","smeared_rnd_PDFsys2",           n_thybins, thybins);   
    smeared_rnd_PDFsys2_test = new TH1D("smeared_rnd_PDFsys2_test","smeared_rnd_PDFsys2_test", n_thybins, thybins);       
    response_PDFsys2_th2=new TH2D("response_PDFsys2_th2","response_PDFsys2_th2",
				(Int_t)n_thybins, (Double_t*)thybins,
				(Int_t)n_thybins, (Double_t*)thybins);      
    
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_PDFsys2_rebin, 
		      (TSpline3*)spline3_PDFsys2, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_PDFsys2, 
		      smeared_rnd_PDFsys2, 
		      smeared_rnd_PDFsys2_test, 
		      (TH2D*)response_PDFsys2_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_PDFsys2);  //for normalization only.
    divideBinWidth(smeared_rnd_PDFsys2); 
    divideBinWidth(smeared_rnd_PDFsys2_test); 
    
    double norm_rnd_PDFsys2=theory_PDFsys2_rebin->Integral()/theory_rnd_PDFsys2->Integral();
    response_PDFsys2_th2->Scale(   norm_rnd_PDFsys2);      
    theory_rnd_PDFsys2->Scale( norm_rnd_PDFsys2);
    smeared_rnd_PDFsys2->Scale(norm_rnd_PDFsys2);  
    smeared_rnd_PDFsys2_test->Scale(norm_rnd_PDFsys2);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_PDFsys2->SetLineColor(kRed); 
    smeared_rnd_PDFsys2->SetLineColor(kBlue);    
    smeared_rnd_PDFsys2_test->SetLineColor(kBlue-5);    
    std::cout<<"done making rnd_PDFsys2 toy MC hists"<<std::endl;        
  }

  TH1D * theory_rnd_PDFsysup       = NULL  , * theory_rnd_PDFsysdown       = NULL;    
  TH1D *smeared_rnd_PDFsysup      = NULL   , *smeared_rnd_PDFsysdown      = NULL;
  TH1D *smeared_rnd_PDFsysup_test = NULL   , *smeared_rnd_PDFsysdown_test = NULL;
  TH2D *   response_PDFsysup_th2     = NULL, *   response_PDFsysdown_th2     = NULL;
  
  if(doPDFsysupdown){
    // PDF 6/2 Pt Scale Var Unc., sysup --------------------
    std::cout<<"making rnd_PDFsysup toy MC hists"<<std::endl;        
    theory_rnd_PDFsysup       = new TH1D("theory_rnd_PDFsysup","theory_rnd_PDFsysup",             n_thybins    , thybins);   
    smeared_rnd_PDFsysup      = new TH1D("smeared_rnd_PDFsysup","smeared_rnd_PDFsysup",           n_thybins, thybins);   
    smeared_rnd_PDFsysup_test = new TH1D("smeared_rnd_PDFsysup_test","smeared_rnd_PDFsysup_test", n_thybins, thybins);       
    response_PDFsysup_th2=new TH2D("response_PDFsysup_th2","response_PDFsysup_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				   (Int_t)n_thybins, (Double_t*)thybins);      

    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_PDFsysup_rebin, 
		      (TSpline3*)spline3_PDFsysup, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_PDFsysup, 
		      smeared_rnd_PDFsysup, 
		      smeared_rnd_PDFsysup_test, 
		      (TH2D*)response_PDFsysup_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_PDFsysup);  //for normalization only.
    divideBinWidth(smeared_rnd_PDFsysup); 
    divideBinWidth(smeared_rnd_PDFsysup_test); 
    
    double norm_rnd_PDFsysup=theory_PDFsysup_rebin->Integral()/theory_rnd_PDFsysup->Integral();
    response_PDFsysup_th2->Scale(   norm_rnd_PDFsysup);      
    theory_rnd_PDFsysup->Scale( norm_rnd_PDFsysup);
    smeared_rnd_PDFsysup->Scale(norm_rnd_PDFsysup);  
    smeared_rnd_PDFsysup_test->Scale(norm_rnd_PDFsysup);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_PDFsysup->SetLineColor(kRed); 
    smeared_rnd_PDFsysup->SetLineColor(kBlue);    
    smeared_rnd_PDFsysup_test->SetLineColor(kBlue-5);    
    std::cout<<"done making rnd_PDFsysup toy MC hists"<<std::endl;        

    
    // PDF 6/2 Pt Scale Var Unc., sysdown -------------------
    std::cout<<"making rnd_PDFsysdown toy MC hists"<<std::endl;        
    theory_rnd_PDFsysdown       = new TH1D("theory_rnd_PDFsysdown","theory_rnd_PDFsysdown",             n_thybins    , thybins);   
    smeared_rnd_PDFsysdown      = new TH1D("smeared_rnd_PDFsysdown","smeared_rnd_PDFsysdown",           n_thybins, thybins);   
    smeared_rnd_PDFsysdown_test = new TH1D("smeared_rnd_PDFsysdown_test","smeared_rnd_PDFsysdown_test", n_thybins, thybins);       
    response_PDFsysdown_th2=new TH2D("response_PDFsysdown_th2","response_PDFsysdown_th2",
				(Int_t)n_thybins, (Double_t*)thybins,
				(Int_t)n_thybins, (Double_t*)thybins);      
    
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_PDFsysdown_rebin, 
		      (TSpline3*)spline3_PDFsysdown, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_PDFsysdown, 
		      smeared_rnd_PDFsysdown, 
		      smeared_rnd_PDFsysdown_test, 
		      (TH2D*)response_PDFsysdown_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_PDFsysdown);  //for normalization only.
    divideBinWidth(smeared_rnd_PDFsysdown); 
    divideBinWidth(smeared_rnd_PDFsysdown_test); 
    
    double norm_rnd_PDFsysdown=theory_PDFsysdown_rebin->Integral()/theory_rnd_PDFsysdown->Integral();
    response_PDFsysdown_th2->Scale(   norm_rnd_PDFsysdown);      
    theory_rnd_PDFsysdown->Scale( norm_rnd_PDFsysdown);
    smeared_rnd_PDFsysdown->Scale(norm_rnd_PDFsysdown);  
    smeared_rnd_PDFsysdown_test->Scale(norm_rnd_PDFsysdown);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_PDFsysdown->SetLineColor(kRed); 
    smeared_rnd_PDFsysdown->SetLineColor(kBlue);    
    smeared_rnd_PDFsysdown_test->SetLineColor(kBlue-5);    

    std::cout<<" done making rnd_PDFsysdown toy MC hists"<<std::endl;        

    std::cout<<std::endl<<"========== done making PDFsys1/2/up/down toy MC hists ============"<<std::endl<<std::endl;        
  }
  
  
  
  ////cosmetics for JER sysup response //PUT SOMEWHERE ELSE/REMOVE, TO DO
  //if(doJERsys){
  //
  //  response_JERsysup_th2->SetTitle(("Toy MC Response Matrix from sysup-JER Smeared "+orderstring+" #sigma_{jet}").c_str());  
  //  response_JERsysup_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  //  response_JERsysup_th2->GetXaxis()->SetNoExponent(true);
  //  response_JERsysup_th2->GetXaxis()->SetMoreLogLabels(true);
  //  response_JERsysup_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  //  response_JERsysup_th2->GetYaxis()->SetNoExponent(true);
  //  response_JERsysup_th2->GetYaxis()->SetMoreLogLabels(true);  
  //  setTH2_ZAxisRange(response_JERsysup_th2);  
  //
  //  //response_JERsysdown_th2->SetTitle(("Toy MC Response Matrix from sysdown-JER Smeared "+orderstring+" #sigma_{jet}").c_str());  
  //  //response_JERsysdown_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  //  //response_JERsysdown_th2->GetXaxis()->SetNoExponent(true);
  //  //response_JERsysdown_th2->GetXaxis()->SetMoreLogLabels(true);
  //  //response_JERsysdown_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  //  //response_JERsysdown_th2->GetYaxis()->SetNoExponent(true);
  //  //response_JERsysdown_th2->GetYaxis()->SetMoreLogLabels(true);  
  //  //setTH2_ZAxisRange(response_JERsysdown_th2);  
  //
  //  response_JERsysdown_th2->SetTitle(("Toy MC Response Matrix from sysdown-JER Smeared "+orderstring+" #sigma_{jet}").c_str());  
  //  response_JERsysdown_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  //  response_JERsysdown_th2->GetXaxis()->SetNoExponent(true);
  //  response_JERsysdown_th2->GetXaxis()->SetMoreLogLabels(true);
  //  response_JERsysdown_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  //  response_JERsysdown_th2->GetYaxis()->SetNoExponent(true);
  //  response_JERsysdown_th2->GetYaxis()->SetMoreLogLabels(true);  
  //  setTH2_ZAxisRange(response_JERsysdown_th2);  
  //}
  //----------------------------END TOY MC COSMETICS/PLOTS--------------------------------------------//    

  




  ///////////////////////////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////////////////////////// 







  
  
  //-----------------------------------------------------------------------------------------//
  //////////////////////  START production of Smeared NLO+NP spectra  /////////////////////////
  //-----------------------------------------------------------------------------------------//
  std::cout<<"creating TH1 for toy NP+"+orderstring+" spectra generation"<<std::endl<<std::endl;
  TH1D *theory_rnd_NP    = new TH1D("theory_rnd_NP","theory_rnd_NP", n_thybins, thybins);   
  TH1D *smeared_rnd_NP = new TH1D("smeared_rnd_NP","smeared_rnd_NP", n_thybins, thybins);     
  TH1D *smeared_rnd_NP_test = new TH1D("smeared_rnd_NP_test","smeared_rnd_NP_test", n_thybins, thybins);     
  TH2D* response_NP_th2=new TH2D("response_NP_th2","response_NP_th2",
				     (Int_t)n_thybins, (Double_t*)thybins,
				     (Int_t)n_thybins, (Double_t*)thybins);
  
  //----------------------------------------------------------------------------------------------------------------//
  std::cout<<"making rnd_NP toy MC hists"<<std::endl;        
  if(useSplineWeights){
    makeToySpectra(theory_NP_rebin, 
		   (TSpline3*)spline3_NP, 
		   (TF1*)fJER, 
		   nEvents, 
		   theory_rnd_NP, 
		   smeared_rnd_NP, 
		   smeared_rnd_NP_test, 
		   (TH2D*) response_NP_th2, 	
		   (TF1*) spline3_NP_ext); }		    //
  else if(useFitWeights&&useModLogFit){
    makeToySpectra(theory_NP_rebin, 
		   (TF1*) logFit_NP_4Toy,//(TF1*)modLog4Fit_NP,
		   (TF1*)fJER, 
		   nEvents, 
		   theory_rnd_NP, 
		   smeared_rnd_NP, 
		   smeared_rnd_NP_test, 
		   (TH2D*) response_NP_th2);}
  else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    makeToySpectra( (TH1D*)theory_NP_rebin, 
		    (TF1*)_7TeVFitA_NP,// (TF1*)_7TeVFitB_NP, //
		    (TF1*)fJER, 
		    nEvents, 
		    theory_rnd_NP, 
		    smeared_rnd_NP, 
		    smeared_rnd_NP_test, 
		    (TH2D*)response_NP_th2); }  

  //----------------------------------------------------------------------------------------------------------------//
  TH1D* theory_rnd_NP_prenorm_prebindiv=(TH1D*)theory_rnd_NP->Clone("theory_rnd_NP_prenorm_prebindiv");//for double checking normalization only, nothing else
  
  divideBinWidth(theory_rnd_NP);  //for normalization only.
  divideBinWidth(smeared_rnd_NP); 
  divideBinWidth(smeared_rnd_NP_test);  

  //std::cout<<"::INTEGRALS::"<<std::endl;
  //std::cout<<"theory_NP_rebin->Integral()="         <<theory_NP_rebin->Integral()<<std::endl;
  //std::cout<<"theory_NP_rebin->Integral(\"width\")="         <<theory_NP_rebin->Integral("width")<<std::endl;
  //std::cout<<"theory_rnd_NP->  Integral()="         <<theory_rnd_NP->  Integral()<<std::endl;  
  //std::cout<<"theory_rnd_NP->  Integral(\"width\")="         <<theory_rnd_NP->  Integral("width")<<std::endl; 
  //assert(false); 
 
  double norm_rnd_NP=theory_NP_rebin->Integral("width")/theory_rnd_NP->Integral("width");
  //std::cout<<"without accounting for bin width: norm_rnd_NP="<<norm_rnd_NP<<std::endl;

  //double norm_rnd_NP_withwidth=theory_NP_rebin->Integral("width")/theory_rnd_NP->Integral("width");
  //std::cout<<"accounting for bin width division: norm_rnd_NP_postdiv="<<norm_rnd_NP_withwidth<<std::endl;
  //std::cout<<"theory_NP_rebin->Integral(\"width\")="<<theory_NP_rebin->Integral("width")<<std::endl;
  //std::cout<<"theory_rnd_NP->  Integral(\"width\")="<<theory_rnd_NP->  Integral("width")<<std::endl;  
  //std::cout<<"norm_rnd_NP_withwidth/norm_rnd_NP="<<norm_rnd_NP_withwidth/norm_rnd_NP<<std::endl;


  response_NP_th2->Scale(   norm_rnd_NP);   
  theory_rnd_NP->Scale( norm_rnd_NP);
  smeared_rnd_NP->Scale(norm_rnd_NP);
  smeared_rnd_NP_test->Scale(norm_rnd_NP);
  
  theory_rnd_NP->SetLineColor(kRed);
  smeared_rnd_NP->SetLineColor(kBlue);    
  smeared_rnd_NP_test->SetLineColor(kBlue-5);    
  
  //  if(theory_rnd_NP->GetBinContent(theory_rnd_NP->GetNbinsX() )<0.)
  //    theory_rnd_NP->SetBinContent(theory_rnd_NP->GetNbinsX(),1.e-10);
  //  if(smeared_rnd_NP->GetBinContent(smeared_rnd_NP->GetNbinsX() )<0.)
  //    smeared_rnd_NP->SetBinContent(smeared_rnd_NP->GetNbinsX(),1.e-08);


  //----------------------------START TOY NLO+NP COSMETICS/PLOTS--------------------------------------------//    
  Float_t NPmaxy=theory_rnd_NP->GetMaximum();    std::cout<<"NPmaxy="<<NPmaxy<<std::endl;
  NPmaxy*=2.;  std::cout<<"NPmaxy="<<NPmaxy<<std::endl;
  Float_t NPminy=smeared_rnd_NP->GetMinimum();   std::cout<<"NPminy="<<NPminy<<std::endl;
  NPminy/=2.; std::cout<<"NPminy="<<NPminy<<std::endl;
  theory_rnd_NP->SetAxisRange(NPminy,NPmaxy,"Y");
  theory_rnd_NP->SetTitle( ((std::string)("Toy MC Incl. NP+"+orderstring+" #sigma;;"+ddxsec_yax)).c_str() );
  
  ///// check plots True/Smeared 
  TCanvas *plot_NP_true_smeared_rat = new TCanvas("plot_NP_true_smeared_rat", "plot NP true smeared rat",1600,1600);  
  plot_NP_true_smeared_rat->cd();
  
  TPad* NPoverlay_pad=new TPad("NPoverlay_pad","pad for spectra",0.,0.3,1.0,1.0);
  NPoverlay_pad->SetLogx(1);
  NPoverlay_pad->SetLogy(1);
  NPoverlay_pad->SetBottomMargin(0);
  NPoverlay_pad->Draw();//tells the canvas theres a pad at the coordinates specified in the new constructor  
  
  TPad* NPratio_pad=new TPad("NPratio_pad","pad for ratio",0.,0.,1.0,0.3);
  NPratio_pad->SetLogx(1);
  NPratio_pad->SetLogy(0);
  NPratio_pad->SetTopMargin(0);
  NPratio_pad->SetBottomMargin(0.20);
  NPratio_pad->Draw();//tell canvas there's a pad at the specified coordinates in the new constructor
  
  // spectra pad
  NPoverlay_pad->cd();    
  theory_rnd_NP->DrawClone("HIST E");  
  smeared_rnd_NP->DrawClone("HIST E SAME");
  smeared_rnd_NP_test->DrawClone("HIST E SAME");
  
  TLegend* NPleg_smear=new TLegend(0.60, 0.70, 0.9, 0.9, NULL,"BRNDC");
  NPleg_smear->SetFillStyle(0);
  NPleg_smear->SetBorderSize(0.);
  NPleg_smear->AddEntry(theory_rnd_NP  , ("Truth Incl. NP+"+orderstring+" #sigma"   ).c_str(), "lp" );
  NPleg_smear->AddEntry(smeared_rnd_NP , ("Smeared Incl. NP+"+orderstring+" #sigma" ).c_str(), "lp" );
  NPleg_smear->AddEntry(smeared_rnd_NP_test , ("Opp Side Smeared Incl. NP+"+orderstring+" #sigma" ).c_str(), "lp" );
  NPleg_smear->Draw();  
      
  // ratio pad
  NPratio_pad->cd();
  
  std::cout<<"creating ratio hists for NP+"+orderstring+""<<std::endl;
  TH1D *NP_true_smeared_rat=(TH1D*)theory_rnd_NP->Clone("NP_true_smeared_rat"); 
  if(printBaseDebug)NP_true_smeared_rat->Print("base");
  NP_true_smeared_rat->Divide(NP_true_smeared_rat,smeared_rnd_NP,1.,1.,"B");     
  
  TH1D *NP_true_smeared_test_rat=(TH1D*)theory_rnd_NP->Clone("NP_true_smeared_test_rat"); 
  if(printBaseDebug)NP_true_smeared_test_rat->Print("base");
  NP_true_smeared_test_rat->Divide(NP_true_smeared_test_rat,smeared_rnd_NP_test,1.,1.,"B");     
  
  NP_true_smeared_rat->SetTitle("; Jet p_{T} [GeV]; True/Smeared");
  NP_true_smeared_rat->GetXaxis()->SetMoreLogLabels(true);
  NP_true_smeared_rat->GetXaxis()->SetNoExponent(true);
  NP_true_smeared_rat->GetXaxis()->SetLabelSize(0.05);
  NP_true_smeared_rat->GetXaxis()->SetTitleSize(0.055);

  NP_true_smeared_rat->GetYaxis()->SetLabelSize(0.05);
  NP_true_smeared_rat->GetYaxis()->SetTitleSize(0.05);
  NP_true_smeared_rat->GetYaxis()->SetTitleOffset(0.7);
  NP_true_smeared_rat->GetYaxis()->CenterTitle(true);
  NP_true_smeared_rat->SetAxisRange(0.85,1.15,"Y");  

  NP_true_smeared_rat->SetLineColor(smeared_rnd_NP->GetLineColor());
  
  //NP_true_smeared_test_rat->SetTitle("; Jet p_T (GeV); True/Smeared");
  //NP_true_smeared_test_rat->GetXaxis()->SetMoreLogLabels(true);
  //NP_true_smeared_test_rat->GetXaxis()->SetNoExponent(true);
  //NP_true_smeared_test_rat->SetAxisRange(0.8,1.2,"Y");
  NP_true_smeared_test_rat->SetLineColor(smeared_rnd_NP_test->GetLineColor());
  
  NP_true_smeared_rat->DrawClone("HIST E");        
  NP_true_smeared_test_rat->DrawClone("HIST E SAME");        
  
  lineatone->DrawClone();
  
  



  // RESPONSE MATRIX PLOT
  TCanvas* plot_response_NP_th2=new TCanvas("plot_response_NP_th2","plot response NP th2",1200, 1000);
  plot_response_NP_th2->cd()->SetLogx(1);
  plot_response_NP_th2->cd()->SetLogy(1);
  plot_response_NP_th2->cd()->SetLogz(1);
  plot_response_NP_th2->cd();
  
  response_NP_th2->SetTitle(("Toy MC Response Matrix from JER Smeared NP+"+orderstring+" #sigma_{jet}").c_str());
  response_NP_th2->GetXaxis()->SetTitle("RECO Jet p_{T} [GeV]");
  response_NP_th2->GetXaxis()->SetNoExponent(true);
  response_NP_th2->GetXaxis()->SetMoreLogLabels(true);
  
  response_NP_th2->GetYaxis()->SetTitle("GEN Jet p_{T} [GeV]");
  response_NP_th2->GetYaxis()->SetNoExponent(true);
  response_NP_th2->GetYaxis()->SetMoreLogLabels(true);
  
  setTH2_ZAxisRange(response_NP_th2);
  
  response_NP_th2->Draw("COLZ");
  //----------------------------END NLO+NP COSMETICS/PLOTS--------------------------------------------//    


  
  //--------------------------------------------------------------------------------------------------------//
  TH1D * theory_rnd_NP_JERsysup       = NULL  ,* theory_rnd_NP_JERsysdown       = NULL  ;
  TH1D *smeared_rnd_NP_JERsysup      = NULL   ,*smeared_rnd_NP_JERsysdown      = NULL   ;
  TH1D *smeared_rnd_NP_JERsysup_test = NULL   ,*smeared_rnd_NP_JERsysdown_test = NULL   ;
  TH2D *   response_NP_JERsysup_th2     = NULL,*   response_NP_JERsysdown_th2     = NULL;
				 
  if(doJERsys){
    //shift JER up
    theory_rnd_NP_JERsysup       = new TH1D("theory_rnd_NP_JERsysup","theory_rnd_NP_JERsysup",             n_thybins    , thybins);   
    smeared_rnd_NP_JERsysup      = new TH1D("smeared_rnd_NP_JERsysup","smeared_rnd_NP_JERsysup",           n_thybins, thybins);   
    smeared_rnd_NP_JERsysup_test = new TH1D("smeared_rnd_NP_JERsysup_test","smeared_rnd_NP_JERsysup_test", n_thybins, thybins);       
    response_NP_JERsysup_th2=new TH2D("response_NP_JERsysup_th2","response_NP_JERsysup_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				   (Int_t)n_thybins, (Double_t*)thybins);  
    
    std::cout<<"making rnd_NP_JERsysup toy MC hists"<<std::endl;
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NP_rebin, 
		      (TSpline3*)spline3_NP, 
		      (TF1*)fJER_sysup, 
		      nEvents, 
		      theory_rnd_NP_JERsysup, 
		      smeared_rnd_NP_JERsysup, 
		      smeared_rnd_NP_JERsysup_test, 
		      (TH2D*)response_NP_JERsysup_th2,
		      NULL); }//		    NULL); }
    //    else if(useFitWeights&&useModLogFit){//looks like modlog4Fit wins!
    //      makeToySpectra( (TH1D*)theory_rebin, 
    //		      (TF1*)logFit_4Toy, //		    (TF1*)modLog4Fit, // (TF1*)modLog3Fit, // 		    
    //		      (TF1*)fJER_sysup, 
    //		      nEvents, 
    //		      theory_rnd_NP_JERsysup, 
    //		      smeared_rnd_NP_JERsysup, 
    //		      smeared_rnd_NP_JERsysup_test, 
    //		      (TH2D*)response_NP_JERsysup_th2); }  
    //    else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    //      makeToySpectra( (TH1D*)theory_rebin, 
    //		      (TF1*)_7TeVFitA,// (TF1*)_7TeVFitB, //
    //		      (TF1*)fJER_sysup, 
    //		      nEvents, 
    //		      theory_rnd_NP_JERsysup, 
    //		      smeared_rnd_NP_JERsysup, 
    //		      smeared_rnd_NP_JERsysup_test, 
    //		      (TH2D*)response_NP_JERsysup_th2); }  
  
    divideBinWidth(theory_rnd_NP_JERsysup);  //for normalization only.
    divideBinWidth(smeared_rnd_NP_JERsysup); 
    divideBinWidth(smeared_rnd_NP_JERsysup_test); 
    
    double norm_rnd_NP_JERsysup=theory_NP_rebin->Integral()/theory_rnd_NP_JERsysup->Integral();
    response_NP_JERsysup_th2->Scale(   norm_rnd_NP_JERsysup);      
    theory_rnd_NP_JERsysup->Scale( norm_rnd_NP_JERsysup);
    smeared_rnd_NP_JERsysup->Scale(norm_rnd_NP_JERsysup);  
    smeared_rnd_NP_JERsysup_test->Scale(norm_rnd_NP_JERsysup);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NP_JERsysup->SetLineColor(kRed); 
    smeared_rnd_NP_JERsysup->SetLineColor(kBlue);    
    smeared_rnd_NP_JERsysup_test->SetLineColor(kBlue-5);    


    //shift JER down
    theory_rnd_NP_JERsysdown       = new TH1D("theory_rnd_NP_JERsysdown","theory_rnd_NP_JERsysdown",             n_thybins    , thybins);   
    smeared_rnd_NP_JERsysdown      = new TH1D("smeared_rnd_NP_JERsysdown","smeared_rnd_NP_JERsysdown",           n_thybins, thybins);   
    smeared_rnd_NP_JERsysdown_test = new TH1D("smeared_rnd_NP_JERsysdown_test","smeared_rnd_NP_JERsysdown_test", n_thybins, thybins);       
    response_NP_JERsysdown_th2=new TH2D("response_NP_JERsysdown_th2","response_NP_JERsysdown_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				   (Int_t)n_thybins, (Double_t*)thybins);  
    
    std::cout<<"making rnd_NP_JERsysdown toy MC hists"<<std::endl;
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NP_rebin, 
		      (TSpline3*)spline3_NP, 
		      (TF1*)fJER_sysdown, 
		      nEvents, 
		      theory_rnd_NP_JERsysdown, 
		      smeared_rnd_NP_JERsysdown, 
		      smeared_rnd_NP_JERsysdown_test, 
		      (TH2D*)response_NP_JERsysdown_th2,
		      NULL); }//		    NULL); }
    //    else if(useFitWeights&&useModLogFit){//looks like modlog4Fit wins!
    //      makeToySpectra( (TH1D*)theory_NP_rebin, 
    //		      (TF1*)logFit_4Toy, //		    (TF1*)modLog4Fit, // (TF1*)modLog3Fit, // 		    
    //		      (TF1*)fJER_sysdown, 
    //		      nEvents, 
    //		      theory_rnd_NP_JERsysdown, 
    //		      smeared_rnd_NP_JERsysdown, 
    //		      smeared_rnd_NP_JERsysdown_test, 
    //		      (TH2D*)response_NP_JERsysdown_th2); }  
    //    else if(useFitWeights&&use7TeVFit){//must decide between A and B TO DO
    //      makeToySpectra( (TH1D*)theory_NP_rebin, 
    //		      (TF1*)_7TeVFitA,// (TF1*)_7TeVFitB, //
    //		      (TF1*)fJER_sysdown, 
    //		      nEvents, 
    //		      theory_rnd_NP_JERsysdown, 
    //		      smeared_rnd_NP_JERsysdown, 
    //		      smeared_rnd_NP_JERsysdown_test, 
    //		      (TH2D*)response_NP_JERsysdown_th2); }  
  
    divideBinWidth(theory_rnd_NP_JERsysdown);  //for normalization only.
    divideBinWidth(smeared_rnd_NP_JERsysdown); 
    divideBinWidth(smeared_rnd_NP_JERsysdown_test); 
    
    double norm_rnd_NP_JERsysdown=theory_NP_rebin->Integral()/theory_rnd_NP_JERsysdown->Integral();
    response_NP_JERsysdown_th2->Scale(   norm_rnd_NP_JERsysdown);      
    theory_rnd_NP_JERsysdown->Scale( norm_rnd_NP_JERsysdown);
    smeared_rnd_NP_JERsysdown->Scale(norm_rnd_NP_JERsysdown);  
    smeared_rnd_NP_JERsysdown_test->Scale(norm_rnd_NP_JERsysdown);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NP_JERsysdown->SetLineColor(kRed); 
    smeared_rnd_NP_JERsysdown->SetLineColor(kBlue);    
    smeared_rnd_NP_JERsysdown_test->SetLineColor(kBlue-5);    
  }



  //--------------------------------------------------------------------------------------------------------//
  TH1D * theory_rnd_NP_PDFsys1       = NULL  ,* theory_rnd_NP_PDFsys2       = NULL  ;
  TH1D *smeared_rnd_NP_PDFsys1      = NULL   ,*smeared_rnd_NP_PDFsys2      = NULL   ;
  TH1D *smeared_rnd_NP_PDFsys1_test = NULL   ,*smeared_rnd_NP_PDFsys2_test = NULL   ;
  TH2D *   response_NP_PDFsys1_th2     = NULL,*   response_NP_PDFsys2_th2     = NULL;

  
  if(doPDFsys12){

    // PDF w/ NP alt choice 1 ---------------
    std::cout<<"making rnd_NP_PDFsys1 toy MC hists"<<std::endl;        
    theory_rnd_NP_PDFsys1       = new TH1D("theory_rnd_NP_PDFsys1","theory_rnd_NP_PDFsys1",             n_thybins    , thybins);   
    smeared_rnd_NP_PDFsys1      = new TH1D("smeared_rnd_NP_PDFsys1","smeared_rnd_NP_PDFsys1",           n_thybins, thybins);   
    smeared_rnd_NP_PDFsys1_test = new TH1D("smeared_rnd_NP_PDFsys1_test","smeared_rnd_NP_PDFsys1_test", n_thybins, thybins);       
    response_NP_PDFsys1_th2=new TH2D("response_NP_PDFsys1_th2","response_NP_PDFsys1_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				   (Int_t)n_thybins, (Double_t*)thybins);  

    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NP_PDFsys1_rebin, 
		      (TSpline3*)spline3_NP_PDFsys1, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NP_PDFsys1, 
		      smeared_rnd_NP_PDFsys1, 
		      smeared_rnd_NP_PDFsys1_test, 
		      (TH2D*)response_NP_PDFsys1_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_NP_PDFsys1);  //for normalization only.
    divideBinWidth(smeared_rnd_NP_PDFsys1); 
    divideBinWidth(smeared_rnd_NP_PDFsys1_test); 
    
    double norm_rnd_NP_PDFsys1=theory_NP_PDFsys1_rebin->Integral()/theory_rnd_NP_PDFsys1->Integral();
    response_NP_PDFsys1_th2->Scale(   norm_rnd_NP_PDFsys1);      
    theory_rnd_NP_PDFsys1->Scale( norm_rnd_NP_PDFsys1);
    smeared_rnd_NP_PDFsys1->Scale(norm_rnd_NP_PDFsys1);  
    smeared_rnd_NP_PDFsys1_test->Scale(norm_rnd_NP_PDFsys1);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NP_PDFsys1->SetLineColor(kRed); 
    smeared_rnd_NP_PDFsys1->SetLineColor(kBlue);    
    smeared_rnd_NP_PDFsys1_test->SetLineColor(kBlue-5);    

    std::cout<<"done making rnd_NP_PDFsys1 toy MC hists"<<std::endl;        


    // PDF w/ NP alt choice 2 ---------------
    std::cout<<"making rnd_NP_PDFsys2 toy MC hists"<<std::endl;        
    theory_rnd_NP_PDFsys2       = new TH1D("theory_rnd_NP_PDFsys2","theory_rnd_NP_PDFsys2",             n_thybins    , thybins);   
    smeared_rnd_NP_PDFsys2      = new TH1D("smeared_rnd_NP_PDFsys2","smeared_rnd_NP_PDFsys2",           n_thybins, thybins);   
    smeared_rnd_NP_PDFsys2_test = new TH1D("smeared_rnd_NP_PDFsys2_test","smeared_rnd_NP_PDFsys2_test", n_thybins, thybins);       
    response_NP_PDFsys2_th2=new TH2D("response_NP_PDFsys2_th2","response_NP_PDFsys2_th2",
					(Int_t)n_thybins, (Double_t*)thybins,
					(Int_t)n_thybins, (Double_t*)thybins);  

    std::cout<<"making rnd_NP_PDFsys2 toy MC hists"<<std::endl;        
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NP_PDFsys2_rebin, 
		      (TSpline3*)spline3_NP_PDFsys2, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NP_PDFsys2, 
		      smeared_rnd_NP_PDFsys2, 
		      smeared_rnd_NP_PDFsys2_test, 
		      (TH2D*)response_NP_PDFsys2_th2,
		      NULL); }
  
    divideBinWidth(theory_rnd_NP_PDFsys2);  //for normalization only.
    divideBinWidth(smeared_rnd_NP_PDFsys2); 
    divideBinWidth(smeared_rnd_NP_PDFsys2_test); 
    
    double norm_rnd_NP_PDFsys2=theory_NP_PDFsys2_rebin->Integral()/theory_rnd_NP_PDFsys2->Integral();
    response_NP_PDFsys2_th2->Scale(   norm_rnd_NP_PDFsys2);      
    theory_rnd_NP_PDFsys2->Scale( norm_rnd_NP_PDFsys2);
    smeared_rnd_NP_PDFsys2->Scale(norm_rnd_NP_PDFsys2);  
    smeared_rnd_NP_PDFsys2_test->Scale(norm_rnd_NP_PDFsys2);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NP_PDFsys2->SetLineColor(kRed); 
    smeared_rnd_NP_PDFsys2->SetLineColor(kBlue);    
    smeared_rnd_NP_PDFsys2_test->SetLineColor(kBlue-5);    

    std::cout<<"done making rnd_NP_PDFsys2 toy MC hists"<<std::endl;        
  }
  
  TH1D * theory_rnd_NP_PDFsysup       = NULL  ,* theory_rnd_NP_PDFsysdown       = NULL  ;
  TH1D *smeared_rnd_NP_PDFsysup      = NULL   ,*smeared_rnd_NP_PDFsysdown      = NULL   ;
  TH1D *smeared_rnd_NP_PDFsysup_test = NULL   ,*smeared_rnd_NP_PDFsysdown_test = NULL   ;
  TH2D *   response_NP_PDFsysup_th2     = NULL,*   response_NP_PDFsysdown_th2     = NULL;
  if(doPDFsysupdown){
    // PDF w/ NP alt choice 1 ---------------
    std::cout<<"making rnd_NP_PDFsysup toy MC hists"<<std::endl;        
    theory_rnd_NP_PDFsysup       = new TH1D("theory_rnd_NP_PDFsysup","theory_rnd_NP_PDFsysup",             n_thybins    , thybins);   
    smeared_rnd_NP_PDFsysup      = new TH1D("smeared_rnd_NP_PDFsysup","smeared_rnd_NP_PDFsysup",           n_thybins, thybins);   
    smeared_rnd_NP_PDFsysup_test = new TH1D("smeared_rnd_NP_PDFsysup_test","smeared_rnd_NP_PDFsysup_test", n_thybins, thybins);       
    response_NP_PDFsysup_th2=new TH2D("response_NP_PDFsysup_th2","response_NP_PDFsysup_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				   (Int_t)n_thybins, (Double_t*)thybins);  

    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NP_PDFsysup_rebin, 
		      (TSpline3*)spline3_NP_PDFsysup, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NP_PDFsysup, 
		      smeared_rnd_NP_PDFsysup, 
		      smeared_rnd_NP_PDFsysup_test, 
		      (TH2D*)response_NP_PDFsysup_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_NP_PDFsysup);  //for normalization only.
    divideBinWidth(smeared_rnd_NP_PDFsysup); 
    divideBinWidth(smeared_rnd_NP_PDFsysup_test); 
    
    double norm_rnd_NP_PDFsysup=theory_NP_PDFsysup_rebin->Integral()/theory_rnd_NP_PDFsysup->Integral();
    response_NP_PDFsysup_th2->Scale(   norm_rnd_NP_PDFsysup);      
    theory_rnd_NP_PDFsysup->Scale( norm_rnd_NP_PDFsysup);
    smeared_rnd_NP_PDFsysup->Scale(norm_rnd_NP_PDFsysup);  
    smeared_rnd_NP_PDFsysup_test->Scale(norm_rnd_NP_PDFsysup);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NP_PDFsysup->SetLineColor(kRed); 
    smeared_rnd_NP_PDFsysup->SetLineColor(kBlue);    
    smeared_rnd_NP_PDFsysup_test->SetLineColor(kBlue-5);    

    std::cout<<"done making rnd_NP_PDFsysup toy MC hists"<<std::endl;        


    // PDF w/ NP alt choice 2 ---------------
    std::cout<<"making rnd_NP_PDFsysdown toy MC hists"<<std::endl;        
    theory_rnd_NP_PDFsysdown       = new TH1D("theory_rnd_NP_PDFsysdown","theory_rnd_NP_PDFsysdown",             n_thybins    , thybins);   
    smeared_rnd_NP_PDFsysdown      = new TH1D("smeared_rnd_NP_PDFsysdown","smeared_rnd_NP_PDFsysdown",           n_thybins, thybins);   
    smeared_rnd_NP_PDFsysdown_test = new TH1D("smeared_rnd_NP_PDFsysdown_test","smeared_rnd_NP_PDFsysdown_test", n_thybins, thybins);       
    response_NP_PDFsysdown_th2=new TH2D("response_NP_PDFsysdown_th2","response_NP_PDFsysdown_th2",
					(Int_t)n_thybins, (Double_t*)thybins,
					(Int_t)n_thybins, (Double_t*)thybins);  

    std::cout<<"making rnd_NP_PDFsysdown toy MC hists"<<std::endl;        
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NP_PDFsysdown_rebin, 
		      (TSpline3*)spline3_NP_PDFsysdown, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NP_PDFsysdown, 
		      smeared_rnd_NP_PDFsysdown, 
		      smeared_rnd_NP_PDFsysdown_test, 
		      (TH2D*)response_NP_PDFsysdown_th2,
		      NULL); }
  
    divideBinWidth(theory_rnd_NP_PDFsysdown);  //for normalization only.
    divideBinWidth(smeared_rnd_NP_PDFsysdown); 
    divideBinWidth(smeared_rnd_NP_PDFsysdown_test); 
    
    double norm_rnd_NP_PDFsysdown=theory_NP_PDFsysdown_rebin->Integral()/theory_rnd_NP_PDFsysdown->Integral();
    response_NP_PDFsysdown_th2->Scale(   norm_rnd_NP_PDFsysdown);      
    theory_rnd_NP_PDFsysdown->Scale( norm_rnd_NP_PDFsysdown);
    smeared_rnd_NP_PDFsysdown->Scale(norm_rnd_NP_PDFsysdown);  
    smeared_rnd_NP_PDFsysdown_test->Scale(norm_rnd_NP_PDFsysdown);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NP_PDFsysdown->SetLineColor(kRed); 
    smeared_rnd_NP_PDFsysdown->SetLineColor(kBlue);    
    smeared_rnd_NP_PDFsysdown_test->SetLineColor(kBlue-5);    

    std::cout<<"done making rnd_NP_PDFsysdown toy MC hists"<<std::endl;        

    std::cout<<std::endl<<"========== done making NP_PDFsys1/2/up/down toy MC hists ============"<<std::endl<<std::endl;        
  }



  TH1D * theory_rnd_NPsys1       = NULL  ,* theory_rnd_NPsys2       = NULL  ;
  TH1D *smeared_rnd_NPsys1      = NULL   ,*smeared_rnd_NPsys2      = NULL   ;
  TH1D *smeared_rnd_NPsys1_test = NULL   ,*smeared_rnd_NPsys2_test = NULL   ;
  TH2D *   response_NPsys1_th2     = NULL,*   response_NPsys2_th2     = NULL;

  TH1D * theory_rnd_NPsys1_origthyratio       = NULL  ,* theory_rnd_NPsys2_origthyratio       = NULL  ;//QA
  TH1D * theory_rnd_NPsys1_NPsys1ratio       = NULL  ,* theory_rnd_NPsys2_NPsys2ratio       = NULL  ;//QA
  TH1D * theory_rnd_NPsys1_ratio_rnd_NPsys2=NULL;

  if(doNPsys12){
    // NP alt choice 1 ----------------

    std::cout<<"making rnd_NPsys1 toy MC hists"<<std::endl;        
    theory_rnd_NPsys1       = new TH1D("theory_rnd_NPsys1","theory_rnd_NPsys1",             n_thybins    , thybins);   
    smeared_rnd_NPsys1      = new TH1D("smeared_rnd_NPsys1","smeared_rnd_NPsys1",           n_thybins, thybins);   
    smeared_rnd_NPsys1_test = new TH1D("smeared_rnd_NPsys1_test","smeared_rnd_NPsys1_test", n_thybins, thybins);       
    response_NPsys1_th2=new TH2D("response_NPsys1_th2","response_NPsys1_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				  (Int_t)n_thybins, (Double_t*)thybins);  
    
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NPsys1_rebin, 
		      (TSpline3*)spline3_NPsys1, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NPsys1, 
		      smeared_rnd_NPsys1, 
		      smeared_rnd_NPsys1_test, 
		      (TH2D*)response_NPsys1_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_NPsys1);  //for normalization only.
    divideBinWidth(smeared_rnd_NPsys1); 
    divideBinWidth(smeared_rnd_NPsys1_test); 
    
    double norm_rnd_NPsys1=theory_NPsys1_rebin->Integral()/theory_rnd_NPsys1->Integral();
    response_NPsys1_th2->Scale(   norm_rnd_NPsys1);      
    theory_rnd_NPsys1->Scale( norm_rnd_NPsys1);
    smeared_rnd_NPsys1->Scale(norm_rnd_NPsys1);  
    smeared_rnd_NPsys1_test->Scale(norm_rnd_NPsys1);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NPsys1->SetLineColor(kRed); 
    smeared_rnd_NPsys1->SetLineColor(kBlue);    
    smeared_rnd_NPsys1_test->SetLineColor(kBlue-5);    
    
    //NP SYS12 QA
    theory_rnd_NPsys1_origthyratio  = (TH1D*)theory_rnd_NPsys1->Clone(((std::string)theory_rnd_NPsys1->GetName()+"_origthyratio").c_str());
    theory_rnd_NPsys1_origthyratio->Divide(theory_rebin);
    theory_rnd_NPsys1_NPsys1ratio  = (TH1D*)theory_rnd_NPsys1->Clone(((std::string)theory_rnd_NPsys1->GetName()+"_NPsys1ratio").c_str());
    theory_rnd_NPsys1_NPsys1ratio->Divide(theory_NPsys1_rebin);
    
    std::cout<<"done making rnd_NPsys1 toy MC hists"<<std::endl;        
    
    // NP alt choice 2 ----------------

    std::cout<<"making rnd_NPsys2 toy MC hists"<<std::endl;        
    theory_rnd_NPsys2       = new TH1D("theory_rnd_NPsys2","theory_rnd_NPsys2",             n_thybins    , thybins);   
    smeared_rnd_NPsys2      = new TH1D("smeared_rnd_NPsys2","smeared_rnd_NPsys2",           n_thybins, thybins);   
    smeared_rnd_NPsys2_test = new TH1D("smeared_rnd_NPsys2_test","smeared_rnd_NPsys2_test", n_thybins, thybins);       
    response_NPsys2_th2=new TH2D("response_NPsys2_th2","response_NPsys2_th2",
				    (Int_t)n_thybins, (Double_t*)thybins,
				    (Int_t)n_thybins, (Double_t*)thybins);  

    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NPsys2_rebin, 
		      (TSpline3*)spline3_NPsys2, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NPsys2, 
		      smeared_rnd_NPsys2, 
		      smeared_rnd_NPsys2_test, 
		      (TH2D*)response_NPsys2_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_NPsys2);  //for normalization only.
    divideBinWidth(smeared_rnd_NPsys2); 
    divideBinWidth(smeared_rnd_NPsys2_test); 
    
    double norm_rnd_NPsys2=theory_NPsys2_rebin->Integral()/theory_rnd_NPsys2->Integral();
    response_NPsys2_th2->Scale(   norm_rnd_NPsys2);      
    theory_rnd_NPsys2->Scale( norm_rnd_NPsys2);
    smeared_rnd_NPsys2->Scale(norm_rnd_NPsys2);  
    smeared_rnd_NPsys2_test->Scale(norm_rnd_NPsys2);  
    

    //plotting smeared+theory toy mc xsection
    theory_rnd_NPsys2->SetLineColor(kRed); 
    smeared_rnd_NPsys2->SetLineColor(kBlue);    
    smeared_rnd_NPsys2_test->SetLineColor(kBlue-5);    

    //NP SYSDOWN QA
    theory_rnd_NPsys2_origthyratio  = (TH1D*)theory_rnd_NPsys2->Clone(((std::string)theory_rnd_NPsys2->GetName()+"_origthyratio").c_str());
    theory_rnd_NPsys2_origthyratio->Divide(theory_rebin);
    theory_rnd_NPsys2_NPsys2ratio  = (TH1D*)theory_rnd_NPsys2->Clone(((std::string)theory_rnd_NPsys2->GetName()+"_NPsys2ratio").c_str());
    theory_rnd_NPsys2_NPsys2ratio->Divide(theory_NPsys2_rebin);


    //SYS12 QA
    theory_rnd_NPsys1_ratio_rnd_NPsys2=(TH1D*)theory_rnd_NPsys1->Clone("theory_rnd_NPsys1_ratio_rnd_NPsys2");
    theory_rnd_NPsys1_ratio_rnd_NPsys2->Divide(theory_rnd_NPsys2);
    
    std::cout<<"done making rnd_NPsys2 toy MC hists"<<std::endl;        

    std::cout<<std::endl<<"========== done making NPsys1/2 toy MC hists ============"<<std::endl<<std::endl;            
    
  }


  TH1D * theory_rnd_NPsysup       = NULL  ,* theory_rnd_NPsysdown       = NULL  ;
  TH1D *smeared_rnd_NPsysup      = NULL   ,*smeared_rnd_NPsysdown      = NULL   ;
  TH1D *smeared_rnd_NPsysup_test = NULL   ,*smeared_rnd_NPsysdown_test = NULL   ;
  TH2D *   response_NPsysup_th2     = NULL,*   response_NPsysdown_th2     = NULL;
  if(doNPsysupdown){

    // NP sysup ----------------
    std::cout<<"making rnd_NPsysup toy MC hists"<<std::endl;
    theory_rnd_NPsysup       = new TH1D("theory_rnd_NPsysup","theory_rnd_NPsysup",             n_thybins    , thybins);   
    smeared_rnd_NPsysup      = new TH1D("smeared_rnd_NPsysup","smeared_rnd_NPsysup",           n_thybins, thybins);   
    smeared_rnd_NPsysup_test = new TH1D("smeared_rnd_NPsysup_test","smeared_rnd_NPsysup_test", n_thybins, thybins);       
    response_NPsysup_th2=new TH2D("response_NPsysup_th2","response_NPsysup_th2",
				  (Int_t)n_thybins, (Double_t*)thybins,
				  (Int_t)n_thybins, (Double_t*)thybins);  
    
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NPsysup_rebin, 
		      (TSpline3*)spline3_NPsysup, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NPsysup, 
		      smeared_rnd_NPsysup, 
		      smeared_rnd_NPsysup_test, 
		      (TH2D*)response_NPsysup_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_NPsysup);  //for normalization only.
    divideBinWidth(smeared_rnd_NPsysup); 
    divideBinWidth(smeared_rnd_NPsysup_test); 
    
    double norm_rnd_NPsysup=theory_NPsysup_rebin->Integral()/theory_rnd_NPsysup->Integral();
    response_NPsysup_th2->Scale(   norm_rnd_NPsysup);      
    theory_rnd_NPsysup->Scale( norm_rnd_NPsysup);
    smeared_rnd_NPsysup->Scale(norm_rnd_NPsysup);  
    smeared_rnd_NPsysup_test->Scale(norm_rnd_NPsysup);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NPsysup->SetLineColor(kRed); 
    smeared_rnd_NPsysup->SetLineColor(kBlue);    
    smeared_rnd_NPsysup_test->SetLineColor(kBlue-5);    
        
    std::cout<<"done making rnd_NPsysup toy MC hists"<<std::endl;        

    // NP sysdown ----------------
    std::cout<<"making rnd_NPsysdown toy MC hists"<<std::endl;        
    theory_rnd_NPsysdown       = new TH1D("theory_rnd_NPsysdown","theory_rnd_NPsysdown",             n_thybins    , thybins);   
    smeared_rnd_NPsysdown      = new TH1D("smeared_rnd_NPsysdown","smeared_rnd_NPsysdown",           n_thybins, thybins);   
    smeared_rnd_NPsysdown_test = new TH1D("smeared_rnd_NPsysdown_test","smeared_rnd_NPsysdown_test", n_thybins, thybins);       
    response_NPsysdown_th2=new TH2D("response_NPsysdown_th2","response_NPsysdown_th2",
				   (Int_t)n_thybins, (Double_t*)thybins,
				  (Int_t)n_thybins, (Double_t*)thybins);  
    
    if(useSplineWeights){//spline fit w/o ext up to ptmax of 686 works well, w/ ext up to 1000 works about same w/ more ad-hoc-ness...
      makeToySpectra( (TH1D*)theory_NPsysdown_rebin, 
		      (TSpline3*)spline3_NPsysdown, 
		      (TF1*)fJER, 
		      nEvents, 
		      theory_rnd_NPsysdown, 
		      smeared_rnd_NPsysdown, 
		      smeared_rnd_NPsysdown_test, 
		      (TH2D*)response_NPsysdown_th2,
		      NULL); }
    
    divideBinWidth(theory_rnd_NPsysdown);  //for normalization only.
    divideBinWidth(smeared_rnd_NPsysdown); 
    divideBinWidth(smeared_rnd_NPsysdown_test); 
    
    double norm_rnd_NPsysdown=theory_NPsysdown_rebin->Integral()/theory_rnd_NPsysdown->Integral();
    response_NPsysdown_th2->Scale(   norm_rnd_NPsysdown);      
    theory_rnd_NPsysdown->Scale( norm_rnd_NPsysdown);
    smeared_rnd_NPsysdown->Scale(norm_rnd_NPsysdown);  
    smeared_rnd_NPsysdown_test->Scale(norm_rnd_NPsysdown);  
    
    //plotting smeared+theory toy mc xsection
    theory_rnd_NPsysdown->SetLineColor(kRed); 
    smeared_rnd_NPsysdown->SetLineColor(kBlue);    
    smeared_rnd_NPsysdown_test->SetLineColor(kBlue-5);    
    std::cout<<"done making rnd_NPsysdown toy MC hists"<<std::endl;        

    std::cout<<std::endl<<"========== done making NPsysup/down toy MC hists ============"<<std::endl<<std::endl;            
  }


  //-----------------------------------------------------------------------------------------//
  //////////////////////  END production of Smeared NLO+NP spectra  ///////////////////////////
  //-----------------------------------------------------------------------------------------//
  
  


  //--------------------------------------------
  //POSITIVE SCALE ERRS
  TCanvas *plot_xsecStatErrs = new TCanvas("plot_xsecStatErrs", "plot xsec Stat Errors",1200,1000);
  plot_xsecStatErrs->cd()->SetLogx(1);
  plot_xsecStatErrs->cd()->SetLogy(1);
  
  //theoryRelStaterr->SetLineColor(kRed);       
  theory_rnd_RelStaterr->SetLineColor(kMagenta);
  smeared_rnd_RelStaterr->SetLineColor(kBlue);

  theoryRelStaterr->SetTitle(("Theory "+orderstring+" #sigma_{Jet} Stat Errors;Jet p_{T} [GeV]; Rel. Errors [%]").c_str());
  theoryRelStaterr->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  theoryRelStaterr->GetYaxis()->SetTitle("Error [%]");
  
  float staterrmin=100000000000000000.;
  if(getnonzeromin(theoryRelStaterr)<staterrmin)staterrmin=getnonzeromin(theoryRelStaterr);
  if(getnonzeromin(theory_rnd_RelStaterr)<staterrmin)staterrmin=getnonzeromin(theory_rnd_RelStaterr);
  if(getnonzeromin(smeared_rnd_RelStaterr)<staterrmin)staterrmin=getnonzeromin(smeared_rnd_RelStaterr);
  staterrmin*=0.1;
  
  
  float staterrmax=-1.; 
  if(theoryRelStaterr->GetMaximum()>staterrmax)staterrmax=theoryRelStaterr->GetMaximum();
  if(theory_rnd_RelStaterr->GetMaximum()>staterrmax)staterrmax=theory_rnd_RelStaterr->GetMaximum();  
  if(smeared_rnd_RelStaterr->GetMaximum()>staterrmax)staterrmax=smeared_rnd_RelStaterr->GetMaximum(); 
  staterrmax*=10.;
  
  theoryRelStaterr->SetMinimum(staterrmin);
  theoryRelStaterr->SetMaximum(staterrmax);
  theoryRelStaterr->GetXaxis()->SetNoExponent(true);
  theoryRelStaterr->GetXaxis()->SetMoreLogLabels(true);  

  TLegend* staterrsleg=new TLegend(0.1,0.7,0.4,0.9);
  staterrsleg->SetBorderSize(0.);
  staterrsleg->SetFillStyle(0);  
  staterrsleg->AddEntry( theoryRelStaterr      , "Theory Stat. Unc." , "l");
  staterrsleg->AddEntry( theory_rnd_RelStaterr , "Toy Theory Stat. Unc" , "l");  
  staterrsleg->AddEntry( smeared_rnd_RelStaterr, "Toy Smeared Stat. Unc." , "l");  
  plot_xsecStatErrs->cd();
  
  theoryRelStaterr->DrawClone("HIST");
  theory_rnd_RelStaterr->DrawClone("HIST SAME");
  smeared_rnd_RelStaterr->DrawClone("HIST SAME");
  
  staterrsleg->Draw();
  
  
  
  
  
  // compare the toy mc truth w/ the original thy hist
  TCanvas *plot_true_rnd_rat = new TCanvas("plot_true_rnd_rat", "plot true rnd rat",1400,1200);  
  plot_true_rnd_rat->cd()->SetLogx(1);
  plot_true_rnd_rat->cd()->SetLogy(0);
  plot_true_rnd_rat->cd();
  std::string thy_rnd_rat_title="Toy MC Truth v. Thy. Spectra; Jet p_{T} [GeV]; Thy./Toy Truth";
  if(useSplineWeights)
    thy_rnd_rat_title="Spline Weighted "+thy_rnd_rat_title;
  if(useFitWeights)
    thy_rnd_rat_title="Fit Weighted "+thy_rnd_rat_title;
  
//  TH1D* theory_ratio_rnd=(TH1D*)theory_rebin->Clone ( 
//						     ( (std::string)theory->GetName() +"_theory_rnd_rat") .c_str()
//						      );
  TH1D* theory_ratio_rnd=(TH1D*)theory_rebin->Clone ( 
						     ( (std::string)theory->GetName() +"_theory_rnd_rat") .c_str()
						      );
  theory_ratio_rnd->Divide(theory_rnd);
  theory_ratio_rnd->SetLineColor(kBlack);
  theory_ratio_rnd->SetTitle(thy_rnd_rat_title.c_str());
  
  
  TH1D* theory_NP_ratio_rnd_NP=(TH1D*)theory_NP_rebin->Clone ( 
							      ( (std::string)theory_NP->GetName() +"_theory_rnd_rat") .c_str()
							       );
  theory_NP_ratio_rnd_NP->Divide(theory_rnd_NP);
  theory_NP_ratio_rnd_NP->SetLineColor(kGray+1);
  theory_NP_ratio_rnd_NP->SetTitle(thy_rnd_rat_title.c_str());
  
  theory_ratio_rnd->SetAxisRange    (0.6,1.4,"Y");
  theory_NP_ratio_rnd_NP->SetAxisRange(0.6,1.4,"Y");
  

  std::cout<<""+orderstring+" THY/RND THY LAST BIN CONTENT = "<<theory_ratio_rnd->GetBinContent(theory_ratio_rnd->GetNbinsX()) <<std::endl;
  std::cout<<""+orderstring+" THY/RND THY LAST BIN CONTENT = "<<theory_NP_ratio_rnd_NP->GetBinContent(theory_NP_ratio_rnd_NP->GetNbinsX()) <<std::endl;
  
  
  theory_ratio_rnd->DrawClone("HIST E");
  theory_NP_ratio_rnd_NP->DrawClone("HIST E SAME");
  
  lineatone->Draw();
  
  
  
  TLegend* leg_tru_rnd_rat=new TLegend(0.10, 0.70, 0.35, 0.9, NULL,"BRNDC");
  leg_tru_rnd_rat->SetFillStyle(0);
  leg_tru_rnd_rat->SetBorderSize(0.);
  leg_tru_rnd_rat->AddEntry( theory_ratio_rnd,       (orderstring+" Spectra"   ).c_str(), "lp");
  leg_tru_rnd_rat->AddEntry( theory_NP_ratio_rnd_NP, (orderstring+"+NP Spectra").c_str(), "lp");
  leg_tru_rnd_rat->Draw();
  
  

  


  


  outf->cd();

  // WRITE CANVASES
  plot_NLOxsec->Write();
  plot_xsecStatErrs->Write();
  plot_xsecPosErrs->Write();
  plot_xsecNegErrs->Write();
  plot_NP->Write(("plot_"+NPCorrFits_str).c_str());
  plot_NPNLOxsec->Write();    
  plot_true_smeared_rat->Write();
  plot_response_th2->Write();  
  plot_NP_true_smeared_rat->Write( );
  plot_response_NP_th2->Write();    
  plot_true_rnd_rat->Write();
  plot_JER->Write();//CANV

  if(useSplineWeights)plot_splines->Write();  
  else if(useFitWeights && useModLogFit) { 
    plot_logfits->Write();  
    plot_logfits_NP->Write();            }
  else if (useFitWeights && use7TeVFit){
    plot_7tevfits->Write();  
    plot_7tevfits_NP->Write();      }    
  
  


  // THY SPECTRA FITS/SPLINES
  if(useSplineWeights){  // SPLINES
    spline3->Write();
    spline3_NP->Write();
    if(useSplineExt)spline3_ext->Write();
    if(useSplineExt)spline3_NP_ext->Write();  }
  else if(useFitWeights && useModLogFit){//FITS
    logFit->Write();
    modLogFit->Write();
    modLog2Fit->Write();
    modLog3Fit->Write();
    modLog4Fit->Write();
    modLog5Fit->Write();
    //w/ NP
    logFit_NP->Write();
    modLogFit_NP->Write();
    modLog2Fit_NP->Write();
    modLog3Fit_NP->Write();
    modLog4Fit_NP->Write();
    modLog5Fit_NP->Write();  }
  else if(useFitWeights && use7TeVFit){//FITS
    _7TeVFitA->Write();
    _7TeVFitB->Write();
    _7TeVFitA_NP->Write();
    _7TeVFitB_NP->Write();  }
  
  
  // JET ENERGY RESOLUTION FITs
  fJER->Write();//FIT
  hJER->Write();//FIT  

  // NP FITS 
  fNP->Write();  
  hNP->Write();
  NPDataPoints->Write();
  
  // THY NLO w/ + w/o NP HISTS
  theory->Write();  
  theory_NP->Write();      


  // PDF ERR HISTS
  theoryRelStaterr->Write();
  theory_rnd_RelStaterr->Write();
  smeared_rnd_RelStaterr->Write();
  theoryPDFPoserr     ->Write();

  //theory2PtScalePoserr->Write();
  theory6PtScalePoserr->Write();
  
  theoryPDFNegerr     ->Write();
  //theory2PtScaleNegerr->Write();
  theory6PtScaleNegerr->Write();      
  
  // REBINNED 
  theory_rebin->Write(); 
  theory_NP_rebin->Write();  
  
  // TOY MC, TRUE + SMEARED + RESPONSE, w/ RATIOS
  theory_rnd->Write(); //truth spectra
  smeared_rnd->Write(); //smeraed spectra
  smeared_rnd_test->Write(); // test-side smeared spectra
  response_th2->Write(); //response matrix
  //ratios
  theory_ratio_rnd->Write(); // thy/truth  
  true_smeared_rat->Write();  // true/smeared
  true_smeared_test_rat->Write();  // true/test-side smeared  
  
  // TOY MC, JER SMEARED NLO+NP SPECTRA
  theory_rnd_NP_prenorm_prebindiv->Write();//truth spectra
  theory_rnd_NP->Write();//truth spectra
  smeared_rnd_NP->Write(); //smeared spectra
  smeared_rnd_NP_test->Write();  //test-side smeared spectra
  response_NP_th2->Write(); //response matrix
  //ratios
  theory_NP_ratio_rnd_NP->Write(); // thy/truth  
  NP_true_smeared_rat->Write();  // true/smeared
  NP_true_smeared_test_rat->Write();  //true/test-side-smeared  

  //JER syst toy mc
  TDirectory* JERsysdir=outf->mkdir("JERsysdir");
  JERsysdir->cd();
  fJER->Write();//FIT
  hJER->Write();

  if(doJERsys){
    
    
   

    
    fJER_sysup->Write();
    hJER_sysup->Write();
    theory_rnd_JERsysup->Write();
    smeared_rnd_JERsysup->Write();
    smeared_rnd_JERsysup_test->Write();
    response_JERsysup_th2->Write();  
    theory_rnd_JERsysdown->Write();
    smeared_rnd_JERsysdown->Write();
    smeared_rnd_JERsysdown_test->Write();
    response_JERsysdown_th2->Write();  

    fJER_sysdown->Write(); 
    hJER_sysdown->Write(); 
    theory_rnd_NP_JERsysup->Write();
    smeared_rnd_NP_JERsysup->Write();
    smeared_rnd_NP_JERsysup_test->Write();
    response_NP_JERsysup_th2->Write();  
    theory_rnd_NP_JERsysdown->Write();
    smeared_rnd_NP_JERsysdown->Write();
    smeared_rnd_NP_JERsysdown_test->Write();
    response_NP_JERsysdown_th2->Write();  

    outf->cd();
  }
  
  //PDF syst toy mc
  if(doPDFsys12||doPDFsysupdown){
    TDirectory* PDFsysdir=outf->mkdir("PDFsysdir");
    PDFsysdir->cd();

    if(doPDFsys12){
      //alt PDF choice 1
      theory_PDFsys1_PDFerrs->Write();//PDF errs
      
      theory_PDFsys1->Write();//no NP
      theory_PDFsys1_rebin->Write();    
      
      spline3_PDFsys1->Write();
      theory_rnd_PDFsys1->Write();
      smeared_rnd_PDFsys1->Write();
      smeared_rnd_PDFsys1_test->Write();
      response_PDFsys1_th2->Write();  
      
      theory_NP_PDFsys1->Write();//w NP
      theory_NP_PDFsys1_rebin->Write();    
      
      spline3_NP_PDFsys1->Write();
      theory_rnd_NP_PDFsys1->Write();
      smeared_rnd_NP_PDFsys1->Write();
      smeared_rnd_NP_PDFsys1_test->Write();
      response_NP_PDFsys1_th2->Write();  
      
      
      //alt PDF choice 2
      theory_PDFsys2_PDFerrs->Write();//PDF errs
      
      theory_PDFsys2->Write();//no NP
      theory_PDFsys2_rebin->Write();    
      
      spline3_PDFsys2->Write();
      theory_rnd_PDFsys2->Write();
      smeared_rnd_PDFsys2->Write();
      smeared_rnd_PDFsys2_test->Write();
      response_PDFsys2_th2->Write();  
      
      theory_NP_PDFsys2->Write();//w NP
      theory_NP_PDFsys2_rebin->Write();    
      
      spline3_NP_PDFsys2->Write();
      theory_rnd_NP_PDFsys2->Write();
      smeared_rnd_NP_PDFsys2->Write();
      smeared_rnd_NP_PDFsys2_test->Write();
      response_NP_PDFsys2_th2->Write();  
    }
    
    if(doPDFsysupdown){
      
      // PDF 6/2 Pt scale var up
      theory6PtScalePoserr->Write();//+scale errs
      //theory2PtScalePoserr->Write();    
      
      theory_PDFsysup->Write("theory_PDFsysup");//no NP
      theory_PDFsysup_rebin->Write("theory_PDFsysup_rebin");    
      
      spline3_PDFsysup->Write("theory_PDFsysup_rebin_spline3");
      theory_rnd_PDFsysup->Write();
      smeared_rnd_PDFsysup->Write();
      smeared_rnd_PDFsysup_test->Write();
      response_PDFsysup_th2->Write();  
      
      theory_NP_PDFsysup->Write("theory_NP_PDFsysup");//w NP
      theory_NP_PDFsysup_rebin->Write("theory_NP_PDFsysup_rebin");    
      
      spline3_NP_PDFsysup->Write("theory_NP_PDFsysup_rebin_spline3");
      theory_rnd_NP_PDFsysup->Write();
      smeared_rnd_NP_PDFsysup->Write();
      smeared_rnd_NP_PDFsysup_test->Write();
      response_NP_PDFsysup_th2->Write();  
      
      
      // PDF 6/2 Pt scale var down
      theory6PtScaleNegerr->Write();//-scale errs
      //theory2PtScaleNegerr->Write();    
      
      theory_PDFsysdown->Write("theory_PDFsysdown");//no NP
      theory_PDFsysdown_rebin->Write("theory_PDFsysdown_rebin");    
      
      spline3_PDFsysdown->Write("theory_PDFsysdown_rebin_spline3");
      theory_rnd_PDFsysdown->Write();
      smeared_rnd_PDFsysdown->Write();
      smeared_rnd_PDFsysdown_test->Write();
      response_PDFsysdown_th2->Write();  
      
      theory_NP_PDFsysdown->Write("theory_NP_PDFsysdown");//w NP
      theory_NP_PDFsysdown_rebin->Write("theory_NP_PDFsysdown_rebin");    
      
      spline3_NP_PDFsysdown->Write("theory_NP_PDFsysdown_rebin_spline3");
      theory_rnd_NP_PDFsysdown->Write();
      smeared_rnd_NP_PDFsysdown->Write();
      smeared_rnd_NP_PDFsysdown_test->Write();
      response_NP_PDFsysdown_th2->Write();  
    }
    outf->cd();
  }
  
  //NP syst toy mc
  if(doNPsys12){
    TDirectory* NPsys12dir=outf->mkdir("NPsys12dir");
    NPsys12dir->cd();
    
    fNP->Write();    
    hNP->Write();
    NPDataPoints->Write();
    
    //altchoice 1
    fNP_sys1->Write();
    hNP_sys1->Write();
    NPsys1_DataPoints->Write();

    theory_NPsys1->Write();
    theory_NPsys1_rebin->Write();
    
    theory_NPsys1_ratio->Write();//QA
    
    spline3_NPsys1->Write();
    theory_rnd_NPsys1->Write();
    smeared_rnd_NPsys1->Write();
    smeared_rnd_NPsys1_test->Write();
    response_NPsys1_th2->Write();  

    theory_rnd_NPsys1_origthyratio->Write();//QA
    theory_rnd_NPsys1_NPsys1ratio->Write();//QA
    
    //altchoice 2
    fNP_sys2->Write();     
    hNP_sys2->Write();
    NPsys2_DataPoints->Write();

    theory_NPsys2->Write();
    theory_NPsys2_rebin->Write();
    
    theory_NPsys2_ratio->Write(); //QA
    
    spline3_NPsys2->Write();
    theory_rnd_NPsys2->Write();
    smeared_rnd_NPsys2->Write();
    smeared_rnd_NPsys2_test->Write();
    response_NPsys2_th2->Write();

    theory_rnd_NPsys2_origthyratio->Write();//QA
    theory_rnd_NPsys2_NPsys2ratio->Write();//QA

    //up+down QA
    theory_NPsys1_ratio_NPsys2->Write();
    hNP_sys1_ratio_sysdown->Write();
    theory_rnd_NPsys1_ratio_rnd_NPsys2->Write();

    outf->cd();
  }
  
  if(doNPsysupdown){
    TDirectory* NPsysupdowndir=outf->mkdir("NPsysupdowndir");
    NPsysupdowndir->cd();

    //sysup
    fNP_sysup->Write();
    hNP_sysup->Write();
    
    theory_NPsysup->Write();
    theory_NPsysup_rebin->Write();
    
    spline3_NPsysup->Write();
    theory_rnd_NPsysup->Write();
    smeared_rnd_NPsysup->Write();
    smeared_rnd_NPsysup_test->Write();
    response_NPsysup_th2->Write();  
    
    //sysdown
    fNP_sysdown->Write();
    hNP_sysdown->Write();
    
    theory_NPsysdown->Write();
    theory_NPsysdown_rebin->Write();
    
    spline3_NPsysdown->Write();
    theory_rnd_NPsysdown->Write();
    smeared_rnd_NPsysdown->Write();
    smeared_rnd_NPsysdown_test->Write();
    response_NPsysdown_th2->Write();  
  }
    


  
  outf->Write();
  std::cout<<"output file "<<outf->GetName()<< " written!"<<std::endl;
  return 1;
  //    assert(false);    
}


//  steering ---------------------------------------------------------------------------------
int main(int argc, char* argv[]){  
  int rStatus = -1;   
  
  if(argc==6){
    rStatus=smearTheorySpectra_gaussCoreJER_etabin( 
						   (std::string) argv[1] ,  
						   (std::string) argv[2], 
						   (std::string) argv[3] ,
						   (int)std::atoi(argv[4]) ,
						   (std::string)argv[5] 						    ) ;      
  }
  
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}
  
  
















