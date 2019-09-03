std::string cpuPtThreshold(TF1*,float);
std::string cpuPtThreshold(float,float,float,float);
void drawAndFitTrigEff(TFile*, std::string, std::string, std::string, TFile*,
		       int, int, 
		       std::string, std::string,
		       std::string, std::string);


std::string etabin_str_arr[5]={
  "0.0",
  "0.5",
  "1.0",
  "1.5",
  "2.0"//,
  //"2.5"
};


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
bool fitVerbose=true;
int trgeff_rebinfactor=1;
void drawAndFitTrigEff(TFile* finData=NULL, std::string thePDFFileName="", std::string fullJetType="", std::string radius="", TFile* fout=NULL,
		       int  etabin_lo=-1, int etabin_hi=-1, 
		       std::string trigToMakeEffFor="", std::string refTrig="", 
		       std::string method="emulation", std::string option=""){
  bool funcDebug=true;

  //CHECK METHOD STRING INPUT
  if(method=="emulation"){
    std::cout<<"using trigger emulation method."<<std::endl;
  }
  else if(method=="reference"){
    std::cout<<"using reference trigger method."<<std::endl;
  }
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
    return;
  }
  else if(etabin_hi<0 || etabin_lo<0){
    std::cout<<"error, etabin_hi and/or etabin_lo cannot be lower than 0. return."<<std::endl;
    return;
  }  
  else if(trigToMakeEffFor==refTrig){
    std::cout<<"error, reference trigger cannot be the same as the trigger you want to determine the efficiency for. return."<<std::endl;
    return;
  }
  else if ((trigToMakeEffFor=="")||(refTrig=="")){
    std::cout<<"error, need reference trigger and/or target trigger to be specified to determine the efficiency, return."<<std::endl;   
    return;
  }
  else{
    std::cout<<"running drawAndFitTrigEff."<<std::endl;
    std::cout<<"measuring Eff for "<<trigToMakeEffFor<<" using "<<refTrig<< " as reference trigger"<<std::endl;
    std::cout<<std::endl;
  }
  
  //GET INPUT HISTS
  std::string numerator_name_base;//"hpp_"+trigToMakeEffFor+"ref_"+refTrig+"_wL1Thresh_wJetID_"+radius+"etabin";
  if(method=="emulation"){
    if(option=="wL1Thresh")
      numerator_name_base="hpp_"+trigToMakeEffFor+"ref_"+refTrig+"_wL1Thresh_wJetID_"+radius+"etabin";
    else 
      numerator_name_base="hpp_"+trigToMakeEffFor+"ref_"+refTrig+"_wJetID_"+radius+"etabin";
  }
  else if(method=="reference")
    numerator_name_base="hpp_"+trigToMakeEffFor+"_wJetID_"+radius+"etabin";
  std::string numerator_name=numerator_name_base+std::to_string(etabin_lo);
  if(funcDebug)std::cout<<"(input) numerator_name="<<numerator_name<<std::endl;
  TH1D* numerator  =(TH1D*)finData->Get(numerator_name.c_str());
  if(trgeff_rebinfactor!=1)numerator->Rebin(trgeff_rebinfactor);
  
  std::string denominator_name_base="hpp_"+refTrig+"_wJetID_"+radius+"etabin";//this doesn't change between emulation and reference trigger methods
  std::string denominator_name=denominator_name_base+std::to_string(etabin_lo);
  if(funcDebug)std::cout<<"(input) denominator_name="<<denominator_name<<std::endl;
  TH1D* denominator=(TH1D*)finData->Get(denominator_name.c_str());
  if(trgeff_rebinfactor!=1)  denominator->Rebin(trgeff_rebinfactor);
  std::cout<<std::endl;
  
  
  //ADD OTHER ETABINS IF RELEVANT
  if(etabin_hi>etabin_lo)
    for(int i=etabin_lo+1; i<=etabin_hi;i++){

      std::string addToNum_name=numerator_name_base+std::to_string(i);
      TH1D* addToNum=(TH1D*)finData->Get(addToNum_name.c_str());
      if(trgeff_rebinfactor!=1)      addToNum->Rebin(trgeff_rebinfactor);
      numerator->Add(addToNum);
      
      std::string addToDenom_name=denominator_name_base+std::to_string(i);
      TH1D* addToDenom=(TH1D*)finData->Get(addToDenom_name.c_str());
      if(trgeff_rebinfactor!=1)      addToDenom->Rebin(trgeff_rebinfactor);
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
  
  //if(trigToMakeEffFor=="HLT60"){
  //  float mean=0.;
  //  for(int i=90;i<=110;i++)
  //    mean+=hTriggerEff->GetBinContent(i);
  //  mean/=21.;
  //  for(int i=0; i<=hTriggerEff->GetNbinsX();i++){
  //    hTriggerEff->SetBinContent(i,hTriggerEff->GetBinContent(i)/mean);
  //    hTriggerEff->SetBinError(  i,hTriggerEff->GetBinError(i) /mean);
  //  }    
  //}
  
  
  //SET PT LIMITS FOR FIT
  float fitpt_lo=-1., fitpt_hi=-1.;
  float drawpt_lo=-1., drawpt_hi=-1.;
  if(trigToMakeEffFor=="HLT80"){ 
    fitpt_lo=83.; fitpt_hi=123.;  
    drawpt_lo=fitpt_lo-5; drawpt_hi=fitpt_hi+5;  }
  else if (trigToMakeEffFor=="HLT60"){
    fitpt_lo=63.; fitpt_hi=103.;  
    drawpt_lo=56.; drawpt_hi=fitpt_hi+100.;  }
  
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
  //hTriggerEff->SetTitle((trigToMakeEffFor+" Efficiency using "+refTrig+" as Reference Trig").c_str());
  hTriggerEff->SetTitle(title.c_str());
  
  hTriggerEff->SetAxisRange(drawpt_lo, drawpt_hi, "X");  
  hTriggerEff->GetXaxis()->SetTitle("Jet p_{T} [GeV]");

  hTriggerEff->GetYaxis()->CenterTitle(true);
  hTriggerEff->GetYaxis()->SetTitle((trigToMakeEffFor+" Efficiency").c_str());
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




  TCanvas* trgEff_canv=new TCanvas(canvname.c_str(),canvname.c_str(), 1800,1400);
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
