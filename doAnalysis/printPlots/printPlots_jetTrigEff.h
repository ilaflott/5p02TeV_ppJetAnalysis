TH1D* safeGetTH1(TFile*, std::string);
std::string cpuPtThreshold(TF1*,float);
std::string cpuPtThreshold(float,float,float,float);
void drawAndFitTrigEff( TFile*, TFile*, 
			std::string, std::string, int, int, 
			std::string, std::string, 
			std::string, std::string, 
			std::string, TFile* );
		       
void drawAndFitTrigEff_trigPt( TFile*, TFile*, 
			       std::string, std::string, int, 
			       std::string, std::string, 
			       std::string, std::string, 
			       std::string, TFile* , std::string);

void printTrigPtHist( TFile*, TFile*, TFile*, 
		      std::string, std::string, std::string, 
		      std::string, TFile* , std::string);


const int TRGEFF_CANVX=1800, TRGEFF_CANVY=1400;
std::string etabin_str_arr[5]={
  "0.0",
  "0.5",
  "1.0",
  "1.5",
  "2.0"//,
  //"2.5"
};

void stupidtwofillfix(TH1* h){
  for(int i=1; i<=h->GetNbinsX();i++){
    h->SetBinContent(i,
		     h->GetBinContent(i)/2.);
    h->SetBinError(i,
		   h->GetBinError(i)/sqrt(2.));
  }
  return;
}

TH1D* safeGetTH1(TFile* file, std::string th1name){
  return (TH1D*)(
		 ((TH1D*)file->Get(th1name.c_str())
		  )->Clone((th1name+"_clone").c_str())
		 );
};

bool fitVerbose=true;
//int trgeff_rebinfactor=1;

//this is to determine the point of 99/98/99.9% efficiencies + return as a string
std::string cpuPtThreshold(TF1* theFit, float eff){
  std::string ptThresh_str="";
  float ptThresh=(float)theFit->TF1::GetX(eff,50.,200.);
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << ptThresh;
  //ptThresh_str=std::to_string(ptThresh);//this will give too many decimals but fine for now.
  ptThresh_str=stream.str();
  return ptThresh_str;
}

std::string cpuPtThreshold(float par0, float par1, float par2, float eff){
  std::string ptThresh_str="";
  float erfinv_arg=(par0+eff-1.)/par0;
  //float erfinv_val=(float)erfinv(erfinv_arg);
  float erfinv_val=3.*erfinv_arg;  //temp for debugging
  float ptThresh=(erfinv_val-par2)/par1;    
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << ptThresh;
  //ptThresh_str=std::to_string(ptThresh);//this will give too many decimals but fine for now.
  ptThresh_str=stream.str();
  return ptThresh_str;
}

//This is to measure the trigger efficiency using a fit function, etc.
void drawAndFitTrigEff( TFile* finDataJet80=NULL,TFile* finDataLowJets=NULL, 			      
			std::string fullJetType="", std::string radius="", int etabin_lo=-1, int etabin_hi=-1,
			std::string trigToMakeEffFor="", std::string refTrig="", 
			std::string method="emulation", std::string option="", 
			std::string thePDFFileName="",  TFile* fout=NULL      ){
  bool funcDebug=true;
  const int rebinfact=1;
  //CHECK METHOD STRING INPUT
  if(method=="emulation"){
    std::cout<<"using trigger emulation method."<<std::endl;  }
  else if(method=="reference"){
    std::cout<<"using reference trigger method."<<std::endl;  }
  //else if(method=="tagAndProbe"){
  //  std::cout<<"using tagAndProbe method."<<std::endl;
  //}
  else{
    //std::cout<<"error, method type not found. Must be of type \"emulation\", \"reference\", or \"tagAndProbe\""<<std::endl;
    std::cout<<"error, method type not found. Must be of type \"emulation\" or \"reference\"."<<std::endl;
    return;
  }
  
  //CHECK ETABIN+TRIGGER STRING INPUTS
  if(etabin_hi<etabin_lo){
    std::cout<<"error, etabin_hi cannot be lower than etabin_lo. return."<<std::endl;
    return;  }
  else if(etabin_hi<0 || etabin_lo<0){
    std::cout<<"error, etabin_hi and/or etabin_lo cannot be lower than 0. return."<<std::endl;
    return;  }  
  else if(trigToMakeEffFor==refTrig){
    std::cout<<"error, reference trigger cannot be the same as the trigger you want to determine the efficiency for. return."<<std::endl;
    return;  }
  else if ((trigToMakeEffFor=="")||(refTrig=="")){
    std::cout<<"error, need reference trigger and/or target trigger to be specified to determine the efficiency, return."<<std::endl;   
    return;  }
  else{
    std::cout<<"running drawAndFitTrigEff."<<std::endl;
    std::cout<<"measuring Eff for "<<trigToMakeEffFor<<" using "<<refTrig<< " as reference trigger"<<std::endl;
    std::cout<<std::endl;
  }
  
  //GET INPUT HISTS //HOW DO I RIG THIS SO I GET WHAT I WANT?!
  TH1D* numerator  =NULL;
  std::string numerator_name_base;//"hpp_"+trigToMakeEffFor+"ref_"+refTrig+"_wL1Thresh_wJetID_"+radius+"etabin";

  //make sure the names are right first
  if(     method=="reference"){   numerator_name_base="hpp_"+trigToMakeEffFor+"_wJetID_"+radius+"etabin";  }
  else if(method=="emulation"){
    if(option=="wL1Thresh")     { numerator_name_base="hpp_"+trigToMakeEffFor+"ref_"+refTrig+"_wL1Thresh_wJetID_"+radius+"etabin";}
    else                        { numerator_name_base="hpp_"+trigToMakeEffFor+"ref_"+refTrig+"_wJetID_"+radius+"etabin";}  }
  else assert(false);
  std::string numerator_name=numerator_name_base+std::to_string(etabin_lo);//final histogram name
  if(funcDebug)std::cout<<"(input) numerator_name="<<numerator_name<<std::endl;  

  //now decide which file the histogram comes from
  if(method=="reference"&&trigToMakeEffFor=="HLT80")
    numerator=(TH1D*)finDataJet80->Get(numerator_name.c_str());     
  else
    numerator=(TH1D*)finDataLowJets->Get(numerator_name.c_str());     
  //stupidtwofillfix(numerator);

  
  
  TH1D* denominator=NULL;
  std::string denominator_name_base="hpp_"+refTrig+"_wJetID_"+radius+"etabin";//this doesn't change between emulation and reference trigger methods
  std::string denominator_name=denominator_name_base+std::to_string(etabin_lo); 
  if(funcDebug)std::cout<<"(input) denominator_name="<<denominator_name<<std::endl;  
  denominator=(TH1D*)finDataLowJets->Get(denominator_name.c_str());//only HLT40/HLT60 get to be reference triggers, HLT80 is no one's reference trigger except higher triggers we don't use
  std::cout<<std::endl;
  //stupidtwofillfix(denominator);


  if(method=="reference"&&trigToMakeEffFor=="HLT80"){
    numerator->Scale(1./effJet80Lumi);
    denominator->Scale(1./effLowJetsLumi);
  }
  if(rebinfact!=1)numerator->Rebin(rebinfact);  
  if(rebinfact!=1)  denominator->Rebin(rebinfact);

  
  //ADD OTHER ETABINS IF RELEVANT
  if(etabin_hi>etabin_lo)
    for(int i=etabin_lo+1; i<=etabin_hi;i++){
      std::string addToNum_name=numerator_name_base+std::to_string(i);
      TH1D* addToNum=NULL;
      if(method=="reference"&&trigToMakeEffFor=="HLT80"){
	addToNum=(TH1D*)finDataJet80->Get(addToNum_name.c_str());
	//stupidtwofillfix(addToNum);
      }
      else{
	addToNum=(TH1D*)finDataLowJets->Get(addToNum_name.c_str());
	//stupidtwofillfix(addToNum);
      }
      if(rebinfact!=1)      addToNum->Rebin(rebinfact);
      numerator->Add(addToNum);
      
      std::string addToDenom_name=denominator_name_base+std::to_string(i);
      TH1D* addToDenom=NULL;
      addToDenom=(TH1D*)finDataLowJets->Get(addToDenom_name.c_str());
      //stupidtwofillfix(addToDenom);
      if(rebinfact!=1)      addToDenom->Rebin(rebinfact);
      denominator->Add(addToDenom);      
    }
  
  //EDIT NAMES ACCORDINGLY FOR WRITING OUTPUT LATER
  if(etabin_lo!=etabin_hi){
    numerator_name+="_to_etabin"+std::to_string(etabin_hi);
    denominator_name+="_to_etabin"+std::to_string(etabin_hi);  }
  numerator_name+="_numerator";
  if(funcDebug)std::cout<<"(output) numerator_name="<<numerator_name<<std::endl;
  denominator_name+="_denominator";      
  if(funcDebug)std::cout<<"(output) denominator_name="<<denominator_name<<std::endl;
  std::cout<<std::endl;
  
  //MAKE TRIGGER EFF CURDVE HIST
  std::string hTriggerEff_name="h"+trigToMakeEffFor+"Eff_using"+refTrig;
  if(etabin_lo==etabin_hi)
    hTriggerEff_name+="_etabin"+std::to_string(etabin_lo);
  else
    hTriggerEff_name+="_etabin"+std::to_string(etabin_lo)+"_to_etabin"+std::to_string(etabin_hi);
  if(funcDebug)    std::cout<<"hTriggerEff_name="<<hTriggerEff_name<<std::endl;
  TH1D* hTriggerEff=(TH1D*)numerator->Clone(hTriggerEff_name.c_str());
  hTriggerEff->Divide(denominator);
  std::cout<<std::endl;
  
  
  //SET PT LIMITS FOR FIT
  float fitpt_lo=-1., fitpt_hi=-1.;
  float drawpt_lo=-1., drawpt_hi=-1.;
  if(trigToMakeEffFor=="HLT80"){ 
    fitpt_lo=85.; fitpt_hi=115.;  
    drawpt_lo=70.;      drawpt_hi=150.;        }
  else if (trigToMakeEffFor=="HLT60"){
    fitpt_lo=65.; fitpt_hi=95.;  
    drawpt_lo=50.; drawpt_hi=130.;  }
  
  //MAKE TRIGGER EFF FIT
  std::string fTriggerEff_name="f"+trigToMakeEffFor+"Eff_using"+refTrig;
  if(etabin_lo==etabin_hi)
    fTriggerEff_name+="_etabin"+std::to_string(etabin_lo);
  else
    fTriggerEff_name+="_etabin"+std::to_string(etabin_lo)+"_to_etabin"+std::to_string(etabin_hi);
  TF1* fTriggerEff=new TF1(fTriggerEff_name.c_str(), "1.-[0]*(std::math::erfc([1]*x+[2]))",
			   fitpt_lo, fitpt_hi ); 
  if(funcDebug)    std::cout<<"fTriggerEff_name="<<fTriggerEff_name<<std::endl;
  fTriggerEff->SetLineColor(kBlue);
  std::cout<<std::endl;
  
  //FIT EFF HIST WITH CURVE
  std::cout<<"// -----------START FITTING----------- //"<<std::endl;
  bool fitTriggerEffStatus=true;
  std::cout<<"pre-fitting fitStatus="<<fitTriggerEffStatus<<std::endl;
  if(fitVerbose){
    std::cout<<"first round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"ER");
    std::cout<<"first round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;
    std::cout<<"second round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"MER");  
    std::cout<<"second round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;}
  else{
    std::cout<<"first round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"ER");
    std::cout<<"first round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;
    std::cout<<"second round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"MER");  
    std::cout<<"second round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;}


  //FIT DETAILS
  float fitpar_0=fTriggerEff->GetParameter(0);
  float fitpar_1=fTriggerEff->GetParameter(1);
  float fitpar_2=fTriggerEff->GetParameter(2);
  float fitparerr_0=fTriggerEff->GetParError(0);
  float fitparerr_1=fTriggerEff->GetParError(1);
  float fitparerr_2=fTriggerEff->GetParError(2);  
  
  std::cout<<"post-fitting fitStatus="<<fitTriggerEffStatus<<std::endl;
  std::cout<<"fitpar_0="<<fitpar_0<<" +/- "<<fitparerr_0<<std::endl;
  std::cout<<"fitpar_1="<<fitpar_1<<" +/- "<<fitparerr_1<<std::endl;
  std::cout<<"fitpar_2="<<fitpar_2<<" +/- "<<fitparerr_2<<std::endl;
  std::cout<<"// -----------END FITTING----------- //"<<std::endl;
  
  TF1* fTriggerEff_errup=(TF1*)fTriggerEff->Clone(((std::string)fTriggerEff->GetName()+"_errup").c_str());
  fTriggerEff_errup->SetLineColor(kRed);
  fTriggerEff_errup->SetLineWidth(fTriggerEff_errup->GetLineWidth()-1);
  fTriggerEff_errup->SetParameter(0,fTriggerEff->GetParameter(0)+fTriggerEff->GetParError(0));
  fTriggerEff_errup->SetParameter(1,fTriggerEff->GetParameter(1)+fTriggerEff->GetParError(1));
  fTriggerEff_errup->SetParameter(2,fTriggerEff->GetParameter(2)+fTriggerEff->GetParError(2));
  TF1* fTriggerEff_errdown=(TF1*)fTriggerEff->Clone(((std::string)fTriggerEff->GetName()+"_errdown").c_str());
  fTriggerEff_errdown->SetLineColor(kRed);
  fTriggerEff_errdown->SetLineWidth(fTriggerEff_errdown->GetLineWidth()-1);
  fTriggerEff_errdown->SetParameter(0,fTriggerEff->GetParameter(0)-fTriggerEff->GetParError(0));
  fTriggerEff_errdown->SetParameter(1,fTriggerEff->GetParameter(1)-fTriggerEff->GetParError(1));
  fTriggerEff_errdown->SetParameter(2,fTriggerEff->GetParameter(2)-fTriggerEff->GetParError(2));

  //DRAW HISTS AND PRINT CANVASES
  std::string title=trigToMakeEffFor+" Eff., ";
  if(method=="reference")
    title+="Reference Method, using "+refTrig;
  else if (method=="emulation"){
    title+="Emulation Method, ";
    if(option=="wL1Thresh")
      title+="with L1 requirement, ";
    title+="using "+refTrig;
  }
  hTriggerEff->SetTitle(title.c_str());  
  hTriggerEff->SetAxisRange(drawpt_lo, drawpt_hi, "X");  
  hTriggerEff->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  hTriggerEff->GetYaxis()->CenterTitle(true);
  hTriggerEff->GetYaxis()->SetTitle((trigToMakeEffFor+" Eff., Using Inclusive Jet p_{T}").c_str());
  hTriggerEff->SetMinimum(0.2);
  hTriggerEff->SetMaximum(1.6);

  std::string canvname=trigToMakeEffFor+"EffCanv_";
  if(method=="reference")
    canvname+="refmethod_";
  else if(method=="emulation"){
    canvname+="emumethod_";
    if(option=="wL1Thresh")
      canvname+="wL1thresh_";
  }
  if(etabin_lo==etabin_hi)
    canvname+="_etabin"+std::to_string(etabin_lo);
  else
    canvname+="_etabin"+std::to_string(etabin_lo)+"_etabin"+std::to_string(etabin_hi);




  TCanvas* trgEff_canv=new TCanvas(canvname.c_str(),canvname.c_str(), TRGEFF_CANVX,TRGEFF_CANVY);
  trgEff_canv->cd();  



  TLine* lineatone=new TLine(drawpt_lo, 1., drawpt_hi, 1.);
  lineatone->SetLineStyle(9);
  lineatone->SetLineColor(kBlack);
  lineatone->SetLineWidth(lineatone->GetLineWidth()-1);

  hTriggerEff->Draw("HIST E][");
  lineatone->Draw();
  hTriggerEff->Draw("HIST E][ SAME");
  fTriggerEff->Draw("SAME");
  //fTriggerEff_errup->Draw("SAME");
  //fTriggerEff_errdown->Draw("SAME");


  //std::string ptthresh_98  =(std::string)cpuPtThreshold(fitpar_0, fitpar_1, fitpar_2, 0.98);  
  //std::string ptthresh_99  =(std::string)cpuPtThreshold(fitpar_0, fitpar_1, fitpar_2, 0.99);
  //std::string ptthresh_99p9=(std::string)cpuPtThreshold(fitpar_0, fitpar_1, fitpar_2, 0.999);
  std::string ptthresh_98  =(std::string)cpuPtThreshold((TF1*)fTriggerEff, 0.98);  ptthresh_98  ="98%   Eff, p_{T} = "+ptthresh_98  +" GeV";
  std::string ptthresh_99  =(std::string)cpuPtThreshold((TF1*)fTriggerEff, 0.99);  ptthresh_99  ="99%   Eff, p_{T} = "+ptthresh_99  +" GeV";
  std::string ptthresh_99p9=(std::string)cpuPtThreshold((TF1*)fTriggerEff, 0.999); ptthresh_99p9="99.9% Eff, p_{T} = "+ptthresh_99p9+" GeV";
  std::cout<<"turn on threshold for 98%   eff="<<ptthresh_98<<" GeV"<<std::endl;
  std::cout<<"turn on threshold for 99%   eff="<<ptthresh_99<<" GeV"<<std::endl;
  std::cout<<"turn on threshold for 99.9% eff="<<ptthresh_99p9<<" GeV"<<std::endl;

  std::string etarange_string=etabin_str_arr[etabin_lo]+ " < #||{y} < "+etabin_str_arr[etabin_hi+1];
  TPaveText* etabin_thresholds=new TPaveText(0.15,0.6,0.45,0.85,"NDC");//, etarange_string.c_str(), "NDC");
  etabin_thresholds->SetBorderSize(0);
  etabin_thresholds->SetFillStyle(0); 
  etabin_thresholds->AddText(etarange_string.c_str());
  etabin_thresholds->AddText(ptthresh_98.c_str());
  etabin_thresholds->AddText(ptthresh_99.c_str());
  etabin_thresholds->AddText(ptthresh_99p9.c_str());  
  etabin_thresholds->Draw();

  trgEff_canv->Print(thePDFFileName.c_str());


  //WRITE OUTPUT TO ROOT FILE
  fout->cd();
  //canvs
  trgEff_canv->Write();
  //hists/fits to write
  numerator->Write();
  denominator->Write();
  hTriggerEff->Write();
  fTriggerEff->Write();
  
  //need to do this because the drawAndFitTrigEff method creates these objects in the same scope
  //meaning that they get the same hists added to them over and over again when they loop over the etabins. so just delete + recreate
  numerator->Delete(); 
  denominator->Delete();
  hTriggerEff->Delete();
  fTriggerEff->Delete();
  
  std::cout<<"done running drawAndFitTrigEff, returning."<<std::endl;
  return;
}





//This is to measure the trigger efficiency using a fit function, etc.
void drawAndFitTrigEff_trigPt( TFile* finDataJet80=NULL,TFile* finDataLowJets=NULL, 			      
			       std::string fullJetType="", std::string radius="", int rebinfact=1,
			       std::string trigToMakeEffFor="", std::string refTrig="", 
			       std::string method="emulation", std::string option="", 
			       std::string thePDFFileName="",  TFile* fout=NULL , std::string etarng=""     ){
  bool funcDebug=true;
  
  //if(method=="emulation"){
  //std::cout<<"can't use emulation method for lead trigger jet pt efficiency yet. exit."<<std::endl;
  //return;  }
  if(etarng!="" && etarng!="_00eta20" && etarng!="_20eta51"){
    std::cout<<"error! etrng can only be empty, 00eta20, or 20eta51. exit"<<std::endl;
    return;
  }
    
    
  //CHECK METHOD STRING INPUT
  if(method=="emulation"){
    std::cout<<"using trigger emulation method."<<std::endl;  }
  else if(method=="reference"){
    std::cout<<"using reference trigger method."<<std::endl;  }
  //else if(method=="tagAndProbe"){
  //  std::cout<<"using tagAndProbe method."<<std::endl;
  //}
  else{
    //std::cout<<"error, method type not found. Must be of type \"emulation\", \"reference\", or \"tagAndProbe\""<<std::endl;
    std::cout<<"error, method type not found. Must be of type \"emulation\" or \"reference\"."<<std::endl;
    return;
  }

  //CHECK ETABIN+TRIGGER STRING INPUTS
  if(trigToMakeEffFor==refTrig){
    std::cout<<"error, reference trigger cannot be the same as the trigger you want to determine the efficiency for. return."<<std::endl;
    return;  }
  else if ((trigToMakeEffFor=="")||(refTrig=="")){
    std::cout<<"error, need reference trigger and/or target trigger to be specified to determine the efficiency, return."<<std::endl;   
    return;  }
  else{
    std::cout<<"running drawAndFitTrigEff_trigPt."<<std::endl;
    std::cout<<"measuring Eff for "<<trigToMakeEffFor<<" using "<<refTrig<< " as reference trigger"<<std::endl;
    std::cout<<std::endl;
  }


  //GET INPUT HISTS //HOW DO I RIG THIS SO I GET WHAT I WANT?!
  TH1D* numerator  =NULL;
  std::string numerator_name_base;//"hpp_"+trigToMakeEffFor+"ref_"+refTrig+"_wL1Thresh_wJetID_"+radius+"etabin";

  //make sure the names are right first
  if(     method=="reference"){       numerator_name_base="Inc"+trigToMakeEffFor+"trgPt_leadobj"+etarng; }
  else if(method=="emulation"){
    if(option=="wL1Thresh")     {     numerator_name_base=trigToMakeEffFor+"ref_Inc"+refTrig+"trgPt_wL1Thresh_leadobj"+etarng;}//too much auto-correlation to use emulation method w/ trg pt
    else                        {     numerator_name_base="Inc"+trigToMakeEffFor+"ref_"+refTrig+"trgPt_leadobj"+etarng; } }
  else assert(false);
  std::string numerator_name=numerator_name_base;
  if(funcDebug)std::cout<<"(input) numerator_name="<<numerator_name<<std::endl;  
  
  //now decide which file the histogram comes from
  if(method=="reference"&&trigToMakeEffFor=="HLT80")
    numerator=(TH1D*)finDataJet80->Get(numerator_name.c_str());     
  else
    numerator=(TH1D*)finDataLowJets->Get(numerator_name.c_str());     
  //stupidtwofillfix(numerator);
	
  
  TH1D* denominator=NULL;
  std::string denominator_name_base="Inc"+refTrig+"trgPt_leadobj"+etarng;//this doesn't change between emulation and reference trigger methods
  std::string denominator_name=denominator_name_base;
  if(funcDebug)std::cout<<"(input) denominator_name="<<denominator_name<<std::endl;  
  denominator=(TH1D*)finDataLowJets->Get(denominator_name.c_str());//only HLT40/HLT60 get to be reference triggers, HLT80 is no one's reference trigger except higher triggers we don't use
  std::cout<<std::endl;
  //stupidtwofillfix(denominator);


  if(method=="reference"&&trigToMakeEffFor=="HLT80"){
    numerator->Scale(  1./effJet80Lumi);//Jet80PD Eff-intLumi
    denominator->Scale(1./effLowJetsLumi);//LowerJetsPD Eff-intLumi
  }
    
  if(rebinfact!=1)  denominator->Rebin(rebinfact);  
  if(rebinfact!=1)  numerator->Rebin(rebinfact);
  
  numerator_name+="_numerator";
  if(funcDebug)std::cout<<"(output) numerator_name="<<numerator_name<<std::endl;
  denominator_name+="_denominator";      
  if(funcDebug)std::cout<<"(output) denominator_name="<<denominator_name<<std::endl;
  std::cout<<std::endl;
  
  //MAKE TRIGGER EFF CURDVE HIST
  std::string hTriggerEff_name="h"+trigToMakeEffFor+"Eff_using"+refTrig+"_trgPt_leadobj";
  if(funcDebug)    std::cout<<"hTriggerEff_name="<<hTriggerEff_name<<std::endl;
  TH1D* hTriggerEff=(TH1D*)numerator->Clone(hTriggerEff_name.c_str());
  hTriggerEff->Divide(denominator);
  std::cout<<std::endl;
  
  //SET PT LIMITS FOR FIT
  float fitpt_lo=-1., fitpt_hi=-1.;
  float drawpt_lo=-1., drawpt_hi=-1.;
  if(trigToMakeEffFor=="HLT80"){ 
    fitpt_lo=90.; fitpt_hi=130.;  
    drawpt_lo=80.;      drawpt_hi=150.;          }  
  else if (trigToMakeEffFor=="HLT60"){
    fitpt_lo=70.; fitpt_hi=110.;  
    drawpt_lo=60.; drawpt_hi=130.;  }
  
  //MAKE TRIGGER EFF FIT
  std::string fTriggerEff_name="f"+trigToMakeEffFor+"Eff_using"+refTrig;
  TF1* fTriggerEff=new TF1(fTriggerEff_name.c_str(), "1.-[0]*(std::math::erfc([1]*x+[2]))",
			   fitpt_lo, fitpt_hi ); 
  if(funcDebug)    std::cout<<"fTriggerEff_name="<<fTriggerEff_name<<std::endl;
  fTriggerEff->SetLineColor(kBlue);
  std::cout<<std::endl;
  
  //FIT EFF HIST WITH CURVE
  std::cout<<"// -----------START FITTING----------- //"<<std::endl;
  bool fitTriggerEffStatus=true;
  std::cout<<"pre-fitting fitStatus="<<fitTriggerEffStatus<<std::endl;
  if(fitVerbose){
    std::cout<<"first round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"ER");
    std::cout<<"first round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;
    std::cout<<"second round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"MER");  
    std::cout<<"second round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;}
  else{
    std::cout<<"first round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"ER");
    std::cout<<"first round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;
    std::cout<<"second round of fitting begins..."<<std::endl;
    fitTriggerEffStatus=hTriggerEff->Fit(fTriggerEff,"MER");  
    std::cout<<"second round of fitting ends... fitStaus="<<fitTriggerEffStatus<<std::endl<<std::endl;}


  //FIT DETAILS
  float fitpar_0=fTriggerEff->GetParameter(0);
  float fitpar_1=fTriggerEff->GetParameter(1);
  float fitpar_2=fTriggerEff->GetParameter(2);
  float fitparerr_0=fTriggerEff->GetParError(0);
  float fitparerr_1=fTriggerEff->GetParError(1);
  float fitparerr_2=fTriggerEff->GetParError(2);  
  
  std::cout<<"post-fitting fitStatus="<<fitTriggerEffStatus<<std::endl;
  std::cout<<"fitpar_0="<<fitpar_0<<" +/- "<<fitparerr_0<<std::endl;
  std::cout<<"fitpar_1="<<fitpar_1<<" +/- "<<fitparerr_1<<std::endl;
  std::cout<<"fitpar_2="<<fitpar_2<<" +/- "<<fitparerr_2<<std::endl;
  std::cout<<"// -----------END FITTING----------- //"<<std::endl;
  
  TF1* fTriggerEff_errup=(TF1*)fTriggerEff->Clone(((std::string)fTriggerEff->GetName()+"_errup").c_str());
  fTriggerEff_errup->SetLineColor(kRed);
  fTriggerEff_errup->SetLineWidth(fTriggerEff_errup->GetLineWidth()-1);
  fTriggerEff_errup->SetParameter(0,fTriggerEff->GetParameter(0)+fTriggerEff->GetParError(0));
  fTriggerEff_errup->SetParameter(1,fTriggerEff->GetParameter(1)+fTriggerEff->GetParError(1));
  fTriggerEff_errup->SetParameter(2,fTriggerEff->GetParameter(2)+fTriggerEff->GetParError(2));
  TF1* fTriggerEff_errdown=(TF1*)fTriggerEff->Clone(((std::string)fTriggerEff->GetName()+"_errdown").c_str());
  fTriggerEff_errdown->SetLineColor(kRed);
  fTriggerEff_errdown->SetLineWidth(fTriggerEff_errdown->GetLineWidth()-1);
  fTriggerEff_errdown->SetParameter(0,fTriggerEff->GetParameter(0)-fTriggerEff->GetParError(0));
  fTriggerEff_errdown->SetParameter(1,fTriggerEff->GetParameter(1)-fTriggerEff->GetParError(1));
  fTriggerEff_errdown->SetParameter(2,fTriggerEff->GetParameter(2)-fTriggerEff->GetParError(2));

  //DRAW HISTS AND PRINT CANVASES
  std::string title=trigToMakeEffFor+" Eff., ";
  if(method=="reference")
    title+="Reference Method, using "+refTrig;
  else if (method=="emulation"){
    title+="Emulation Method, ";
    if(option=="wL1Thresh")
      title+="with L1 requirement, ";
    title+="using "+refTrig;
  }
  //hTriggerEff->SetTitle((trigToMakeEffFor+" Efficiency using "+refTrig+" as Reference Trig").c_str());
  hTriggerEff->SetTitle(title.c_str());
  
  hTriggerEff->SetAxisRange(drawpt_lo, drawpt_hi, "X");  
  hTriggerEff->GetXaxis()->SetTitle("Jet p_{T} [GeV]");

  hTriggerEff->GetYaxis()->CenterTitle(true);
  hTriggerEff->GetYaxis()->SetTitle((trigToMakeEffFor+" Eff., Using Lead Trigger Jet p_{T}").c_str());
  hTriggerEff->SetMinimum(0.2);
  hTriggerEff->SetMaximum(1.6);

  std::string canvname=trigToMakeEffFor+"EffCanv_";
  if(method=="reference")
    canvname+="refmethod_";
  else if(method=="emulation"){
    canvname+="emumethod_";
    if(option=="wL1Thresh")
      canvname+="wL1thresh_";
  }

  TCanvas* trgEff_canv=new TCanvas(canvname.c_str(),canvname.c_str(),TRGEFF_CANVX,TRGEFF_CANVY);  
  trgEff_canv->cd();  

  TLine* lineatone=new TLine(drawpt_lo, 1., drawpt_hi, 1.);
  lineatone->SetLineStyle(9);
  lineatone->SetLineColor(kBlack);
  lineatone->SetLineWidth(lineatone->GetLineWidth()-1);

  hTriggerEff->Draw("HIST E][");
  lineatone->Draw();
  hTriggerEff->Draw("HIST E][ SAME");
  fTriggerEff->Draw("SAME");
  //fTriggerEff_errup->Draw("SAME");
  //fTriggerEff_errdown->Draw("SAME");


  //std::string ptthresh_98  =(std::string)cpuPtThreshold(fitpar_0, fitpar_1, fitpar_2, 0.98);  
  //std::string ptthresh_99  =(std::string)cpuPtThreshold(fitpar_0, fitpar_1, fitpar_2, 0.99);
  //std::string ptthresh_99p9=(std::string)cpuPtThreshold(fitpar_0, fitpar_1, fitpar_2, 0.999);
  std::string ptthresh_98  =(std::string)cpuPtThreshold((TF1*)fTriggerEff, 0.98);  ptthresh_98  ="98%   Eff, p_{T} = "+ptthresh_98  +" GeV";
  std::string ptthresh_99  =(std::string)cpuPtThreshold((TF1*)fTriggerEff, 0.99);  ptthresh_99  ="99%   Eff, p_{T} = "+ptthresh_99  +" GeV";
  std::string ptthresh_99p9=(std::string)cpuPtThreshold((TF1*)fTriggerEff, 0.999); ptthresh_99p9="99.9% Eff, p_{T} = "+ptthresh_99p9+" GeV";
  std::cout<<"turn on threshold for 98%   eff="<<ptthresh_98<<" GeV"<<std::endl;
  std::cout<<"turn on threshold for 99%   eff="<<ptthresh_99<<" GeV"<<std::endl;
  std::cout<<"turn on threshold for 99.9% eff="<<ptthresh_99p9<<" GeV"<<std::endl;

  //std::string etarange_string=etabin_str_arr[etabin_lo]+ " < #||{y} < "+etabin_str_arr[etabin_hi+1];
  std::string etarange_string;
  if(etarng=="")etarange_string="#||{y} < 5.1";
  if(etarng=="_00eta20")etarange_string="#||{y} < 2.0";
  if(etarng=="_20eta51")etarange_string="2.0 < #||{y} < 5.1";
  TPaveText* etabin_thresholds=new TPaveText(0.15,0.6,0.45,0.85,"NDC");//, etarange_string.c_str(), "NDC");
  etabin_thresholds->SetBorderSize(0);
  etabin_thresholds->SetFillStyle(0); 
  etabin_thresholds->AddText(etarange_string.c_str());
  etabin_thresholds->AddText(ptthresh_98.c_str());
  etabin_thresholds->AddText(ptthresh_99.c_str());
  etabin_thresholds->AddText(ptthresh_99p9.c_str());  
  etabin_thresholds->Draw();
  
  trgEff_canv->Print(thePDFFileName.c_str());

  //WRITE OUTPUT TO ROOT FILE
  fout->cd();
  //canvs
  trgEff_canv->Write();
  //hists/fits to write
  numerator->Write();
  denominator->Write();
  hTriggerEff->Write();
  fTriggerEff->Write();
  
  //need to do this because the drawAndFitTrigEff method creates these objects in the same scope
  //meaning that they get the same hists added to them over and over again when they loop over the etabins. so just delete + recreate
  numerator->Delete(); 
  denominator->Delete();
  hTriggerEff->Delete();
  fTriggerEff->Delete();
  
  std::cout<<"done running drawAndFitTrigEff_trigPt, returning."<<std::endl;
  return;
}





void printTrigPtHist( TFile* finLowJets , TFile* finJet80 , TFile* fin , 
		      std::string fullJetType , std::string radius , std::string trigType,    
		      std::string thePDFFileName , TFile* fout=NULL, std::string etarng=""){
  bool funcDebug=true;
  const int rebinfact=1;
  if(funcDebug)std::cout<<"in printTrigPtHist, trigType="<<trigType<<", hope i don't crash!"<<std::endl;
  std::string jetCutString, jetEtaCutString, jtptQACut_str;
  makeJetCutStrings(&jetCutString, &jetEtaCutString, &jtptQACut_str, fin);
  std::string sampleDescString=sqrts2k15ppString;
  
  
  TH1D* HLTCombo_trgPt=(TH1D*)safeGetTH1(fin, "HLTComb_trgPt"+etarng);
  TH1D* HLT40_trgPt    =NULL, *HLT60_trgPt    =NULL, *HLT80_trgPt    =NULL;
  if(funcDebug)std::cout<<"now retreiving hists from files"<<std::endl;
  if(     trigType=="incl"){
    HLT40_trgPt=(TH1D*)safeGetTH1(finLowJets,"IncHLT40trgPt_leadobj"+etarng);
    HLT60_trgPt=(TH1D*)safeGetTH1(finLowJets,"IncHLT60trgPt_leadobj"+etarng);
    HLT80_trgPt=(TH1D*)safeGetTH1(finJet80,  "IncHLT80trgPt_leadobj"+etarng);
    //stupidtwofillfix(HLT40_trgPt);
    //stupidtwofillfix(HLT60_trgPt);
    //stupidtwofillfix(HLT80_trgPt);
    HLT40_trgPt->Scale(1./effLowJetsLumi); 
    HLT60_trgPt->Scale(1./effLowJetsLumi); 
    HLT80_trgPt->Scale(1./effJet80Lumi);   
        
    
    ////remake HLTCombo if i haven't fixed that bug yet --> i have now! 11/26/19
    //TH1D* HLT40=(TH1D*)safeGetTH1(finLowJets,"is40_trgPt"  );    
    //TH1D* HLT60=(TH1D*)safeGetTH1(finLowJets,"is60_trgPt"  );    
    //TH1D* HLT80=(TH1D*)safeGetTH1(finJet80,"is80_trgPt"  );     
    //HLT40->Scale(1./effLowJetsLumi); 
    //HLT60->Scale(1./effLowJetsLumi); 
    //HLT80->Scale(1./effJet80Lumi);   
  }
  else if(trigType=="excl"){
    //HLT40_trgPt=(TH1D*)safeGetTH1(fin,"is40_trgPt"  );    
    //HLT60_trgPt=(TH1D*)safeGetTH1(fin,"is60_trgPt"  );    
    //HLT80_trgPt=(TH1D*)safeGetTH1(fin,"is80_trgPt"  );     
    HLT40_trgPt=(TH1D*)safeGetTH1(finLowJets,"is40_trgPt"+etarng  );    
    HLT60_trgPt=(TH1D*)safeGetTH1(finLowJets,"is60_trgPt"+etarng  );    
    HLT80_trgPt=(TH1D*)safeGetTH1(finJet80,"is80_trgPt"  +etarng  );     
    HLT40_trgPt->Scale(1./effLowJetsLumi); 
    HLT60_trgPt->Scale(1./effLowJetsLumi); 
    HLT80_trgPt->Scale(1./effJet80Lumi);   
  }  

  ////rebin
  HLTCombo_trgPt=(TH1D*)HLTCombo_trgPt->TH1::Rebin(rebinfact,((std::string)HLTCombo_trgPt->GetName()+"_"+trigType+"_rebin").c_str());   
  HLT40_trgPt   =(TH1D*)HLT40_trgPt   ->TH1::Rebin(rebinfact,((std::string)HLT40_trgPt->GetName()+"_"+trigType+"_rebin").c_str());      
  HLT60_trgPt   =(TH1D*)HLT60_trgPt   ->TH1::Rebin(rebinfact,((std::string)HLT60_trgPt->GetName()+"_"+trigType+"_rebin").c_str());      
  HLT80_trgPt   =(TH1D*)HLT80_trgPt   ->TH1::Rebin(rebinfact,((std::string)HLT80_trgPt->GetName()+"_"+trigType+"_rebin").c_str());      
  HLTCombo_trgPt->Scale(1./(float)rebinfact);
  HLT40_trgPt   ->Scale(1./(float)rebinfact);
  HLT60_trgPt   ->Scale(1./(float)rebinfact);
  HLT80_trgPt   ->Scale(1./(float)rebinfact);
  
  if(funcDebug)std::cout<<"trigPtHistStyle"<<std::endl;
  trigPtHistStyle(HLT40_trgPt  ,1);
  trigPtHistStyle(HLT60_trgPt  ,2);
  trigPtHistStyle(HLT80_trgPt  ,3);
  trigPtHistStyle(HLTCombo_trgPt,0);
  
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  std::string canvname="leadTrigJetPt_"+trigType+"SpectraRatio";
  TCanvas *temp_canvJetTrig = new TCanvas(canvname.c_str(),canvname.c_str() , TRGEFF_CANVX,TRGEFF_CANVY);
  temp_canvJetTrig->cd();

  TPad *jetpad_spectra = new TPad(("jetpad_spectra_"+trigType).c_str(), "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  TPad *jetpad_ratio = new TPad(  ("jetpad_ratio_"  +trigType).c_str(), "Ratio Pad", 0.0, 0.05, 1.0, 0.3);  
  if(funcDebug)std::cout<<"setupJetTrigSpectraRatioCanvas"<<std::endl;
  setupJetTrigSpectraRatioCanvas(temp_canvJetTrig, jetpad_spectra,  jetpad_ratio);
  


  //DRAW SPECTRA
  jetpad_spectra->cd();
  
  TH1D* stylehist = (TH1D*)HLTCombo_trgPt->Clone(("styleonly_"+trigType).c_str());
  stylehist->Reset("ICES");
  stylehist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          
  stylehist->GetXaxis()->SetNoExponent(true);
  stylehist->GetXaxis()->SetMoreLogLabels(true);
  stylehist->SetMaximum(2.e+06);
  stylehist->SetMinimum(1.e-02);

  std::string h_Title;
  if(trigType=="excl")h_Title="Lead Trigger Jet p_{T} for Exclusive HLT40/60/80";
  else if(trigType=="incl")h_Title="Lead Trigger Jet p_{T} for HLT40/60/80";      
  
  stylehist->SetTitle(h_Title.c_str());
  //stylehist->GetXaxis()->SetTitle("Lead Trigger Jet p_{T} [GeV]");
  //stylehist->GetYaxis()->SetTitle(AUAxTitle.c_str());
  stylehist->GetYaxis()->SetTitle(ddcrossSectionAxTitle.c_str());
  
  stylehist->DrawClone();

  HLTCombo_trgPt->DrawClone("SAME");
  HLT40_trgPt->DrawClone("SAME");
  HLT60_trgPt->DrawClone("SAME");
  HLT80_trgPt->DrawClone("SAME");

  
  float t1Loc1=0.50, t1Loc2=0.82;      float textsize=30.;
  TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str(),textsize);      t1->Draw();
  TLatex *t2=NULL;
  if(etarng=="")t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0 < #||{y} < 5.1"),textsize);   
  else if(etarng=="_00eta20")t2=makeTLatex((t1Loc1),(t1Loc2-.04),("0.0 < #||{y} < 2.0"),textsize); 
  else if(etarng=="_20eta51")t2=makeTLatex((t1Loc1),(t1Loc2-.04),("2.0 < #||{y} < 5.1"),textsize); 
  else return;
  t2->Draw();  
  float pp1x=t1Loc1,pp1y=t1Loc2-0.08;
  TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str(),textsize );      pp1->Draw();      
  TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() ,textsize);      pp2->Draw();        
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.65,0.90,0.87, NULL,"brNDC");
  JetTrigLegend->SetBorderSize(0);
  JetTrigLegend->SetFillStyle(0);
  JetTrigLegend->AddEntry(HLT40_trgPt,"HLT40" , "lp");
  JetTrigLegend->AddEntry(HLT60_trgPt,"HLT60" , "lp");
  JetTrigLegend->AddEntry(HLT80_trgPt,"HLT80" , "lp");
  JetTrigLegend->AddEntry(HLTCombo_trgPt,"HLTCombo" , "lp");
  JetTrigLegend->Draw(); //spectra legend

  //END DRAW SPECTRA
  
  TH1D* HLT40_trgPt_ratio=(TH1D*)HLT40_trgPt->Clone("HLT40_trgPt_ratio");
  TH1D* HLT60_trgPt_ratio=(TH1D*)HLT60_trgPt->Clone("HLT60_trgPt_ratio");
  TH1D* HLT80_trgPt_ratio=(TH1D*)HLT80_trgPt->Clone("HLT80_trgPt_ratio");
  TH1D* theDenominator =(TH1D*)HLTCombo_trgPt->Clone("trgPt_denom");
  
  HLT40_trgPt_ratio->Divide(theDenominator);
  HLT60_trgPt_ratio->Divide(theDenominator);
  HLT80_trgPt_ratio->Divide(theDenominator);
  
  trigRatioHistStyle(HLT40_trgPt_ratio, 1);            
  trigRatioHistStyle(HLT60_trgPt_ratio, 2);            
  trigRatioHistStyle(HLT80_trgPt_ratio, 3);            
  trigRatioHistStyle(theDenominator,0);            
  

  //DRAW RATIO
  jetpad_ratio->cd();

  TH1D* stylehist_ratio = (TH1D*)HLT80_trgPt_ratio->Clone(("styleonly_ratio_"+trigType).c_str());
  stylehist_ratio->Reset("ICES");
  stylehist_ratio->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");          
  stylehist_ratio->SetAxisRange(0.0,1.3,"Y");
  stylehist_ratio->GetXaxis()->SetNoExponent(true);
  stylehist_ratio->GetXaxis()->SetMoreLogLabels(true);  
  stylehist_ratio->GetXaxis()->SetTitle("Lead Trigger Jet p_{T} [GeV]");
  stylehist_ratio->GetYaxis()->SetTitle("Ratio to HLT Combo");
  stylehist_ratio->GetYaxis()->CenterTitle(true);
  stylehist_ratio->SetTitle("");

  //float min=stylehist_ratio->GetBinLowEdge(1);
  //float max=stylehist_ratio->GetBinLowEdge(stylehist_ratio->GetNbinsX())+stylehist_ratio->GetBinWidth(stylehist_ratio->GetNbinsX());
  //TLine* lineAtOne          = new TLine(min,1.0,max,1.0);
  TLine* lineAtOne          = new TLine(ptbins_debug[0],1.0,ptbins_debug[nbins_pt_debug],1.0);
  lineAtOne->SetLineColor(12);          
  lineAtOne->SetLineStyle(2);
  
  stylehist_ratio->DrawClone();
  HLT40_trgPt_ratio->DrawClone("HIST ][ SAME");
  HLT60_trgPt_ratio->DrawClone("HIST ][ SAME");
  HLT80_trgPt_ratio->DrawClone("HIST ][ SAME");  
  lineAtOne->Draw("same");  
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());
  fout->cd();
  temp_canvJetTrig->Write();
  temp_canvJetTrig->Close();
  
  return;
}












/*


void printJetTrigHist_wRatio( TFile* fin , bool usedHLT100, bool usedMinBias, bool analysisRebin, int etabin,
			      std::string thePDFFileName , std::string fullJetType , 
			      std::string trigType, std::string radius , bool usedHLTPF, bool didJetID , TFile* fout=NULL){
  bool funcDebug=true;
  if(funcDebug){
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------- "<<std::endl;
  std::cout<<" --------------- DRAWING JET SPECTRA BY TRIGGER --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------- "<<std::endl<<std::endl<<std::endl;  }
  
  if(!fin ){    std::cout<<"input file not found, cannot look at event counts"<<std::endl; 
    return; }
  
  if(trigType=="excl")std::cout<<"drawing exclusive trigger spectra + ratios"<<std::endl;
  else if(trigType=="incl")std::cout<<"drawing inclusive trigger spectra + ratios"<<std::endl;
  else {
    std::cout<<"trigger spectra + ratios not recognized, exit"<<std::endl; 
    return; }
  
  std::string jetCutString, jetEtaCutString, jtptQACut_str;
  makeJetCutStrings(&jetCutString, &jetEtaCutString, &jtptQACut_str, fin);  
  std::string sampleDescString=sqrts2k15ppString;  
  if(etabin==0)jetEtaCutString="0.0 < #||{y} < 0.5";
  if(etabin==1)jetEtaCutString="0.5 < #||{y} < 1.0";
  if(etabin==2)jetEtaCutString="1.0 < #||{y} < 1.5";
  if(etabin==3)jetEtaCutString="1.5 < #||{y} < 2.0";
  if(etabin==4)jetEtaCutString="2.0 < #||{y} < 2.5";
  if(etabin==5)jetEtaCutString="2.5 < #||{y} < 3.0";
  
  if(funcDebug)std::cout<<std::endl<<"creating temporary canvas for printing JetTrig plots..."<<std::endl;
  TCanvas *temp_canvJetTrig = new TCanvas(("tempJetTrigSpec"+trigType).c_str(), ("blahTrigPt"+trigType).c_str(), canvx_trig, canvy_trig);
  temp_canvJetTrig->cd();  
  
  TPad *jetpad_spectra = new TPad("jetpad_spectra", "Overlay Pad", 0.0, 0.30, 1.0, 1.0);
  TPad *jetpad_ratio = new TPad("jetpad_ratio", "Ratio Pad", 0.0, 0.05, 1.0, 0.3);  
  setupJetTrigSpectraRatioCanvas(temp_canvJetTrig, jetpad_spectra,  jetpad_ratio);
  jetpad_spectra->cd();
  
  TLegend* JetTrigLegend=new TLegend(0.70,0.66,0.85,0.86, NULL,"brNDC");

  TH1D* theDenominator=NULL;//for later ratio drawing
  
  for(int j=0; j<(N_trigs); j++){
    
    std::string inHistName;
    if(j==0)inHistName="hpp_"+HLTName[j];
    else{
      if(trigType=="excl")inHistName="hpp_exc"+HLTName[j];
      else if(trigType=="incl")inHistName="hpp_"+HLTName[j];    
    }    
    if(didJetID)inHistName+="_wJetID";
    //inHistName+="_"+radius+etaWidth;    
    inHistName+="_"+radius+"etabin"+std::to_string(etabin);    
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;
    if(HLTName[j]=="HLTMB"&&!usedMinBias)continue;
    
    if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl;
    //assert(false);
    
    //get hist
    TH1D* theJetTrigQAHist= (TH1D*) ( (TH1*)fin->Get(inHistName.c_str()) );
    if(funcDebug)theJetTrigQAHist->Print("base"); std::cout<<std::endl;
    theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
    theJetTrigQAHist->GetXaxis()->SetNoExponent(true);

    if(!usedMinBias && HLTName[j]=="HLTComb"&&trigType=="excl"){
      std::string inMBHistName;
      if(trigType=="excl")inMBHistName="hpp_excHLTMB";
      else if(trigType=="incl")inMBHistName="hpp_HLTMB";
      else assert(false);
      if(didJetID)inMBHistName+="_wJetID";
      inMBHistName+="_"+radius+"etabin"+std::to_string(etabin);
      TH1D* HLTMBexcl_toSubtract=(TH1D*)fin->Get(inMBHistName.c_str());
      if((bool)HLTMBexcl_toSubtract)
	TH1Sub_correrrs(theJetTrigQAHist, HLTMBexcl_toSubtract);
    }
    
    //if((HLTName[j]=="HLTMB" || HLTName[j]=="HLTComb") && trigType=="excl"&&etabin==0){
    //  //      TH1D* MB_toSubtract_l1jet68=(TH1D*)fin->Get("excHLTMB_AND_l1jet68"); //subtract just this one first
    //  //      TH1D* MB_toSubtract_l1jet60=(TH1D*)fin->Get("excHLTMB_AND_l1jet60"); //then also subtract this one...
    //  //      TH1D* MB_toSubtract_l1jet52=(TH1D*)fin->Get("excHLTMB_AND_l1jet52"); //then subtract this one too...
    //  //      //  TH1D* MB_toSubtract_l1jet48=(TH1D*)fin->Get("excHLTMB_AND_l1jet48"); //etc...
    //  //      //  TH1D* MB_toSubtract_l1jet48=(TH1D*)fin->Get("excHLTMB_AND_l1jet44"); //etc...
    //  //      //  TH1D* MB_toSubtract_l1jet40=(TH1D*)fin->Get("excHLTMB_AND_l1jet40");
    //  //      //  TH1D* MB_toSubtract_l1jet28=(TH1D*)fin->Get("excHLTMB_AND_l1jet28");
    //  //      //  TH1D* MB_toSubtract_l1jet24=(TH1D*)fin->Get("excHLTMB_AND_l1jet24");
    //  //        
    //  //      TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet68);
    //  //      TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet60);
    //  //      TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet52);
    //  //      //  TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet48);
    //  //      //  TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet44);
    //  //      
    //  //      
    //  //      //TH1D* MB_toSubtract_HLTjet80=(TH1D*)fin->Get("excHLTMB_AND_HLTjet80"); //subtract just this one first
    //  //      //TH1D* MB_toSubtract_HLTjet60=(TH1D*)fin->Get("excHLTMB_AND_HLTjet60"); //then also subtract this one...
    //  //      //TH1D* MB_toSubtract_HLTjet40=(TH1D*)fin->Get("excHLTMB_AND_HLTjet40"); //then subtract this one too...
    //  //      
    //  //      //TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_HLTjet80);
    //  //      //TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_HLTjet60);
    //  //      //TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_HLTjet40);
    //  TH1D* MB_toSubtract_NOTHFcoinc=(TH1D*)fin->Get("excHLTMB_NOT_HFcoincidence");
    //  TH1Sub_correrrs(theJetTrigQAHist, MB_toSubtract_NOTHFcoinc);      
    //}
	

    
    //rebin
    if(analysisRebin) {
      theJetTrigQAHist=(TH1D*)theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_spec_anabins").c_str(), ptbins_debug    );
      divideBinWidth(theJetTrigQAHist);    }
    else {
      theJetTrigQAHist=(TH1D*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_spec_rebin"+std::to_string(rebinfactor)).c_str() );
      theJetTrigQAHist->Scale(1./((float)rebinfactor));    
    }
    theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X");

    //TEMPORARY!?
    if(trigType=="incl"){
      if      (HLTName[j]=="HLT80")
	theJetTrigQAHist->Scale(1./effJet80Lumi);
      else if (HLTName[j]=="HLT60" || HLTName[j]=="HLT40")
	theJetTrigQAHist->Scale(1./effLowJetsLumi);
    }    
    
    //legend entry titles
    if(funcDebug)std::cout<<"setting legend entry titles"<<std::endl;
    if(usedHLTPF){//      if(funcDebug)std::cout<<"using HLTPF Trig titles"<<std::endl;
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTPFName_Leg[j]).c_str(),"lp");  }
    else {//   if(funcDebug)std::cout<<"using HLTCalo Trig titles"<<std::endl;
      JetTrigLegend->AddEntry(theJetTrigQAHist,(HLTCaloName_Leg[j]).c_str(),"lp");    }
        
    //    if( j==0 ){      //draws combo first
    if(HLTName[j]=="HLTComb"){//j==0

      //for drawing ratios late
      theDenominator=(TH1D*)theJetTrigQAHist->Clone(( ((std::string)theJetTrigQAHist->GetTitle())+"_denomclone").c_str());
      trigRatioHistStyle(theDenominator,j);      
      
      std::string h_Title;
      if(trigType=="excl")h_Title="Excl. Jet Spectra by Trigger";
      else if(trigType=="incl")h_Title="Incl. Jet Spectra by Trigger";      
      std::string h_YAx_Title= AUAxTitle;

      //set titles
      theJetTrigQAHist->SetTitle (    h_Title.c_str() );      
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );      
      
      //style
      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw("");      
      
      //descriptive text on plot
      float t1Loc1=0.54, t1Loc2=0.82;
      TLatex *t1=makeTLatex(t1Loc1,t1Loc2,(fullJetType+"Jets").c_str());      t1->Draw();
      TLatex *t2=makeTLatex((t1Loc1),(t1Loc2-.04),(jetCutString).c_str());      t2->Draw();
      TLatex *t3=makeTLatex((t1Loc1),(t1Loc2-.08),(jetEtaCutString).c_str());      t3->Draw();
      
      float pp1x=t1Loc1,pp1y=t1Loc2-0.12;
      TLatex *pp1=makeTLatex(pp1x,pp1y,sqrts2k15ppString.c_str() );      pp1->Draw();
      TLatex *pp2=makeTLatex(pp1x,(pp1y-0.04),intLumiString.c_str() );    pp2->Draw();      
      if(fout){
	fout->cd();
	theJetTrigQAHist->Write();
	jetpad_spectra->cd();
      }

    }
    else{//rest of HLT spectra


      trigSpectraHistStyle(theJetTrigQAHist, j);
      theJetTrigQAHist->Draw(" SAME");
      if(fout){
	fout->cd();
	theJetTrigQAHist->Write();
	jetpad_spectra->cd();
      }
    }
    
  }//end loop over trigs for spectra
  
  JetTrigLegend->Draw();
  //done drawing spectra pad
  
  if(funcDebug){
  std::cout<<std::endl<<std::endl<<" ------------------------------------------------------------------ "<<std::endl;
  std::cout<<" --------------- DRAWING RATIOS  --------------- "<<std::endl;
  std::cout<<" ------------------------------------------------------------------ "<<std::endl<<std::endl<<std::endl;}
  
  
  jetpad_ratio->cd();

  for(int j=0; j<N_trigs; j++){
    
    std::string inHistName;
    if(HLTName[j]=="HLT100"&&!usedHLT100)continue;    
    if(HLTName[j]=="HLTMB"&&!usedMinBias)continue;
    

    if(j>0){//individial hlt paths, numerators, combohist setup by now
      
      if(trigType=="excl")      inHistName="hpp_exc"+HLTName[j];//+"_"+radius+etaWidth;
      else if(trigType=="incl") inHistName="hpp_"+HLTName[j];//+"_"+radius+etaWidth;
      if(didJetID)inHistName+="_wJetID";
      inHistName+="_"+radius+"etabin"+std::to_string(etabin);//etaWidth;      
      if(funcDebug)std::cout<<"inHistName="<<inHistName<<std::endl<<std::endl;
      
      TH1D* theJetTrigQAHist= (TH1D*) ( (TH1*)fin->Get(inHistName.c_str()) );
      if(funcDebug)theJetTrigQAHist->Print("base");
      //if((HLTName[j]=="HLTMB" || HLTName[j]=="HLTComb") && trigType=="excl"&&etabin==0){
      //	//TH1D* MB_toSubtract_l1jet68=(TH1D*)fin->Get("excHLTMB_AND_l1jet68"); //subtract just this one first
      //	//TH1D* MB_toSubtract_l1jet60=(TH1D*)fin->Get("excHLTMB_AND_l1jet60"); //then also subtract this one...
      //	//TH1D* MB_toSubtract_l1jet52=(TH1D*)fin->Get("excHLTMB_AND_l1jet52"); //then subtract this one too...
      //	//  TH1D* MB_toSubtract_l1jet48=(TH1D*)fin->Get("excHLTMB_AND_l1jet48"); //etc...
      //	//  TH1D* MB_toSubtract_l1jet40=(TH1D*)fin->Get("excHLTMB_AND_l1jet40");
      //	//  TH1D* MB_toSubtract_l1jet28=(TH1D*)fin->Get("excHLTMB_AND_l1jet28");
      //	//  TH1D* MB_toSubtract_l1jet24=(TH1D*)fin->Get("excHLTMB_AND_l1jet24");
      //	
      //	//TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet68);
      //	//TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet60);
      //	//TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet52);
      //	//  TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_l1jet48);
      //	  
      //	  
      //	//TH1D* MB_toSubtract_HLTjet80=(TH1D*)fin->Get("excHLTMB_AND_HLTjet80"); //subtract just this one first
      //	//TH1D* MB_toSubtract_HLTjet60=(TH1D*)fin->Get("excHLTMB_AND_HLTjet60"); //then also subtract this one...
      //	//TH1D* MB_toSubtract_HLTjet40=(TH1D*)fin->Get("excHLTMB_AND_HLTjet40"); //then subtract this one too...
      //	  
      //	//TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_HLTjet80);
      //	//TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_HLTjet60);
      //	//TH1Sub_correrrs(theJetTrigQAHist,MB_toSubtract_HLTjet40);
      //
      //	TH1D* MB_toSubtract_NOTHFcoinc=(TH1D*)fin->Get("excHLTMB_NOT_HFcoincidence");
      //	TH1Sub_correrrs(theJetTrigQAHist, MB_toSubtract_NOTHFcoinc);
      //	
      //}      

      if(analysisRebin) 
	theJetTrigQAHist= (TH1D*) theJetTrigQAHist->TH1::Rebin(nbins_pt_debug,(inHistName+"_rat_anabins").c_str(), ptbins_debug    );
      else {
	theJetTrigQAHist=(TH1D*)theJetTrigQAHist->TH1::Rebin(rebinfactor, (inHistName+"_rat_rebin"+std::to_string(rebinfactor)).c_str() );
	theJetTrigQAHist->Scale(1./((float)rebinfactor));    
      }

      //TEMPORARY!?
      if(trigType=="incl"){
	if      (HLTName[j]=="HLT80")
	  theJetTrigQAHist->Scale(1./effJet80Lumi);
	else if (HLTName[j]=="HLT60" || HLTName[j]=="HLT40")
	  theJetTrigQAHist->Scale(1./effLowJetsLumi);
      }    

      theJetTrigQAHist->SetAxisRange(ptbins_debug[0],ptbins_debug[nbins_pt_debug],"X"); 
      
      trigRatioHistStyle(theJetTrigQAHist, j);
      
      


      theJetTrigQAHist->Divide(theDenominator);
      
      theJetTrigQAHist->SetAxisRange(0.0,1.5,"Y");
      
      std::string h_XAx_Title="Jet p_{T} (GeV)     ";
      theJetTrigQAHist->SetXTitle( h_XAx_Title.c_str() );
      theJetTrigQAHist->GetXaxis()->SetMoreLogLabels(true);
      theJetTrigQAHist->GetXaxis()->SetNoExponent(true);
      
      std::string h_YAx_Title="Combo/Indiv.";
      theJetTrigQAHist->SetYTitle( h_YAx_Title.c_str() );
      
      theJetTrigQAHist->SetTitle("");	
      
      if(j==1)theJetTrigQAHist->Draw("");
      else theJetTrigQAHist->Draw(" SAME");
      
      
      if( j==1 )  {
	
	float min = theJetTrigQAHist->GetBinLowEdge(1);
	float max = ptbins_debug[nbins_pt_debug];
	
	TLine* lineAtOne          = new TLine(min,1.0,max,1.0);
	lineAtOne->SetLineColor(12);          lineAtOne->SetLineStyle(2);
	
	lineAtOne->Draw("same");
	
      }//end j=1 specific
      
      if(fout){
	fout->cd();
	theJetTrigQAHist->Write();
	jetpad_ratio->cd();
      }
    }//end j>0 
    
  }//end loop over trigs for ratio
  
  temp_canvJetTrig->Print(thePDFFileName.c_str());
  
  if(fout){
    fout->cd();
    TCanvas* outcanv=(TCanvas*)temp_canvJetTrig->DrawClone();
    
    std::string outcanvtitle = "HLTak4PF"; 
    if(usedHLT100)outcanvtitle+="40 to 100 ";
    else outcanvtitle+="40 to 80 ";
    outcanvtitle+=fullJetType+" Jet p_{T}";    
    if(analysisRebin)outcanvtitle+=" Ana. Bins ";
    else             outcanvtitle+=" Simp. Bins ";
    outcanvtitle +=trigType;
    outcanv->SetTitle((outcanvtitle+" Canvas").c_str());
    
    std::string outcanvname  = "HLTak4PF";
    if(usedHLT100)outcanvname+="40to100_";
    else          outcanvname+="40to80_";
    outcanvname+=fullJetType+"_jtpt_";
    if(analysisRebin) outcanvname+="anabins_";
    else              outcanvname+="simpbins_";
    outcanvname +=trigType+"_etabin"+std::to_string(etabin);  
    outcanv->Write((outcanvname+"_canv").c_str());   
    
    
   }//end fout
  
  temp_canvJetTrig->Close();
  
  
  
  
  return;
}





*/
