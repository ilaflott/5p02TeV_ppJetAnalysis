void divideBinWidth(TH1 *h){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in divideBinWidth"<<std::endl<<std::endl;

  //h->Sumw2();
  if(funcDebug)h->Print("base");
  double nEntries=h->GetEntries();
  int numbins=h->GetNbinsX();

  //std::cout<<"#of bin should be = "<<nbins_pt<<std::endl;
  if(funcDebug)std::cout<<"#of bins = "<<numbins<<std::endl;
  if(funcDebug)std::cout<<"integral before = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries before = "<<nEntries<<std::endl;

  for (int i=1;i<=numbins;i++)
    {

      if(funcDebug)std::cout<<std::endl<<"bin-i = "<< i <<std::endl;
      if(funcDebug)std::cout<<"bin-i width = "<<  h->GetBinWidth(i)  << std::endl;
      if(funcDebug)std::cout<<"bin-i pt range is "<<   h->GetBinLowEdge(i) << " to " << (h->GetBinLowEdge(i)+h->GetBinWidth(i)) << std::endl;
      if(funcDebug)std::cout<<std::endl<<"old val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"old valErr = "<<      h->GetBinError(i)   << std::endl;


      Float_t val    = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);

      val    /= h->GetBinWidth(i);
      valErr /= h->GetBinWidth(i);

      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);

      if(funcDebug)std::cout<<"new val = "   <<      h->GetBinContent(i) << std::endl;
      if(funcDebug)std::cout<<"new valErr = "<<      h->GetBinError(i)   << std::endl;

    }

  if(funcDebug)std::cout<<"integral after = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries after = "<<h->GetEntries()<<std::endl;

  //h->SetEntries(nEntries);

  //std::cout<<"exiting function"<<std::endl;
  if(funcDebug)std::cout<<std::endl<<"divideBinWidth done"<<std::endl<<std::endl;
  return;
}





float getIntegralJetSpectra(TH1* h){
  bool funcDebug=false;
  
  float integral=0.;
  int numbins=h->GetNbinsX();
  std::cout<<"numbins="<<numbins<<std::endl;
  for(int i=1;i<=numbins;i++){
    
    Float_t val    = h->GetBinContent(i);
    if(funcDebug)std::cout<<"bin i="<<i<<std::endl;
    if(funcDebug)std::cout<<"bin width = "<<h->GetBinWidth(i)<<std::endl;
    if(funcDebug)std::cout<<"pt range of "<<h->GetBinLowEdge(i)<<" to "<<h->GetBinLowEdge(i)+h->GetBinWidth(i)<<std::endl;
    if(funcDebug)std::cout<<"val = "<<val<<std::endl;
    val    *= h->GetBinWidth(i);
    integral+=val;
    if(funcDebug)std::cout<<"integral now = "<<integral<<std::endl<<std::endl;
  }

  if(funcDebug)std::cout<<"final integral is:" << integral<<std::endl;
  if(funcDebug)std::cout<<"normal integral says:"<<h->Integral()<<std::endl;
  
  return integral;
}




float getIntegralErrJetSpectra(TH1* h){
  float integralErr=0.;
  int numbins=h->GetNbinsX();
  for(int i=1;i<=numbins;i++){
    Float_t valErr = h->GetBinError(i);
    valErr *= h->GetBinWidth(i);
    integralErr+=valErr*valErr;
  }
  integralErr=std::sqrt(integralErr);
  return integralErr;
}



//void specialRatio(TH1F* theRatio){
void fracSubtracted(TH1F* theRatio){
  bool funcDebug=false;
  int nbins=theRatio->GetNbinsX();
  for(int i=1; i<=nbins;i++){
    
    float content=theRatio->GetBinContent(i);
    
    if(funcDebug)std::cout<<"before bin i="<<i<<" content="<<content<<std::endl;
    content*=-1.;
    content+=1.;
    content*=100.;
    if(content==100.)continue;//avoid weird scaling shit
    //if(funcDebug)std::cout<<"now    bin i="<<i<<" content="<<content<<std::endl;
    theRatio->SetBinContent( i,content);
    if(funcDebug)std::cout<<"after  bin i="<<i<<" content="<<theRatio->GetBinContent(i)<<std::endl;
  }
  //  assert(false);
  return;
}



//template<typename T>
//void showMinMax() {
//  std::cout << "min: " << numeric_limits<T>::min() << std::endl;
//  std::cout << "max: " << numeric_limits<T>::max() << std::endl;
//  std::cout << std::endl;
//  return;
//}




float getmaxcontent(TH1* th1){
  
  //std::cout << "float: " << std::endl;
  //showMinMax<float>;
  
  float maxcontent=-1.e+40;
  int nbinsx=th1->GetNbinsX();
  for(int i=1;i<=nbinsx;i++){
    float content=th1->GetBinContent(i);
    if(content > maxcontent)//check that new content is less than current minium
      maxcontent=content;    
    else continue;//if the content is greater than current minimum
  }
  return maxcontent;  
}



float getnonzeromin(TH1* th1){
  float nonzeromin=1.e+40;
  int nbinsx=th1->GetNbinsX();
  for(int i=1;i<=nbinsx;i++){
    float content=th1->GetBinContent(i);
    if(!(content>0.)&&!(content<0.))//check for zero
      continue;
    else if(content < nonzeromin)//check that new content is less than current minium
      nonzeromin=content;
    else continue;//if the content is greater than current minimum
  }
  return nonzeromin;

}
