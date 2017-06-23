#include <cmath>

TH2F* reBinTH2(TH2F* inputTH2, std::string rebinTH2_name,
               double* boundaries_pt_reco, int nbins_pt_reco,
               double* boundaries_pt_gen  , int nbins_pt_gen           ){
  std::cout<<std::endl<<"in reBinTH2"<<std::endl<<std::endl;
  bool funcDebug = true;

  //double numEntries=inputTH2->GetEntries();
  //inputTH2->Sumw2();
  if(funcDebug)inputTH2->Print("base");  std::cout<<std::endl;


  TAxis *xaxis = inputTH2->GetXaxis(); //reco pt axis
  int nbins_x= xaxis->GetNbins();
  int xbinstart=1;

  TAxis *yaxis = inputTH2->GetYaxis(); //gen pt axis
  int nbins_y = yaxis->GetNbins();
  int ybinstart=1;

  TH2F *reBinnedTH2 = new TH2F(rebinTH2_name.c_str(), inputTH2->GetTitle(),
                               nbins_pt_reco, boundaries_pt_reco ,
                               nbins_pt_gen, boundaries_pt_gen );

  //if(doOverUnderflows){
  //  xbinstart--;
  //  ybinstart--;
  //  nbins_x++;
  //  nbins_y++;
  //}
  //reBinnedTH2->Sumw2();

  for (  int xbin=xbinstart ; xbin <= nbins_x ; xbin++ ) {
    for (int ybin=ybinstart ; ybin <= nbins_y ; ybin++ ) {

      reBinnedTH2->Fill( xaxis->GetBinCenter(xbin) , yaxis->GetBinCenter(ybin) ,
                         inputTH2->GetBinContent(xbin,ybin) );

    }  }//end x-y loop

  //reBinnedTH2->SetEntries(numEntries);

  if(funcDebug) std::cout<<"getSumW2N="<<reBinnedTH2->GetSumw2N()<<std::endl;
  if(funcDebug) reBinnedTH2->Print("base");  std::cout<<std::endl;


  std::cout<<std::endl<<"reBinTH2 done"<<std::endl<<std::endl;
  return reBinnedTH2;
}





void divideBinWidth_TH2(TH2F *h){
  std::cout<<std::endl<<"in divideBinWidth_TH2"<<std::endl<<std::endl;
  bool funcDebug=false;

  //h->Sumw2();
  if(funcDebug) h->Print("base");      std::cout<<std::endl;

  double numEntries=h->GetEntries();
  int nbinsx=h->GetNbinsX();
  int nbinsy=h->GetNbinsY();

  if(funcDebug)std::cout<<"#of xbins = "<<nbinsx<<std::endl;
  if(funcDebug)std::cout<<"#of ybins = "<<nbinsy<<std::endl;
  if(funcDebug)std::cout<<"total #of bins = "<<nbinsy*nbinsx<<std::endl;
  if(funcDebug)std::cout<<"integral before = "<<h->Integral()<<std::endl;
  if(funcDebug)std::cout<<"entries before = "<<numEntries<<std::endl;

  TAxis* the_xAxis=h->GetXaxis();
  TAxis* the_yAxis=h->GetYaxis();
  for (int i=1;i<=nbinsx;i++)    {
    for (int j=1;j<=nbinsy;j++)    {

      Float_t xWidth=h->GetXaxis()->GetBinWidth(i);
      Float_t yWidth=h->GetYaxis()->GetBinWidth(j);



      if(funcDebug)std::cout<<std::endl<<"xbin-i = "<< i <<std::endl;
      if(funcDebug)std::cout<<"xbin-i width = "<<  the_xAxis->GetBinWidth(i)  << std::endl;
      if(funcDebug)std::cout<<"xbin-i pt range is "<<   the_xAxis->GetBinLowEdge(i) << " to "
                            << (the_xAxis->GetBinLowEdge(i)+the_xAxis->GetBinWidth(i)) << std::endl;

      if(funcDebug)std::cout<<std::endl<<"ybin-j = "<< j <<std::endl;
      if(funcDebug)std::cout<<"ybin-j width = "<<  the_yAxis->GetBinWidth(j)  << std::endl;
      if(funcDebug)std::cout<<"ybin-j pt range is "<<   the_yAxis->GetBinLowEdge(j) << " to "
                            << (the_yAxis->GetBinLowEdge(j)+the_yAxis->GetBinWidth(j)) << std::endl;


      if(funcDebug)std::cout<<std::endl<<"old val = "   <<      h->GetBinContent(i,j) << std::endl;
      if(funcDebug)std::cout<<"old valErr = "<<      h->GetBinError(i,j)   << std::endl;


      Float_t val = (h->GetBinContent(i,j))/(xWidth*yWidth);
      //Float_t valErr = (h->GetBinError(i,j))/(xWidth*yWidth);

      h->SetBinContent(i,j,val);
      //h->SetBinError(i,j,valErr);
      h->SetBinError(i,j,0.0);

      if(funcDebug)std::cout<<std::endl<<"new val = "   <<      h->GetBinContent(i,j) << std::endl;
      if(funcDebug)std::cout<<"new valErr = "<<      h->GetBinError(i,j)   << std::endl;

    }
  }

  //h->SetEntries(numEntries);
  if(funcDebug)h->Print("base");      std::cout<<std::endl;

  std::cout<<std::endl<<"divideBinWidth_TH2 done"<<std::endl<<std::endl;
  return;
}

//normalizes MC response TH2 to the rowSum of the gen-axis
void normalizeMC_TH2(TH2F* inputTH2){
  std::cout<<std::endl<<"in normalizeMC_TH2"<<std::endl<<std::endl;

  inputTH2->Sumw2();
  inputTH2->Print("base");  std::cout<<std::endl;

  double numEntries=inputTH2->GetEntries();

  int colStart=1;
  TAxis *xaxis = inputTH2->GetXaxis();
  int nbins_x= xaxis->GetNbins();

  int rowStart=1;
  TAxis *yaxis = inputTH2->GetYaxis();
  int nbins_y = yaxis->GetNbins();

  //if(doOverUnderflows){
  //  nbins_x++;    nbins_y++;
  //  rowStart--;   colStart--; }

  const int rowSums_len=nbins_x;
  float rowSums[rowSums_len]={0.};
  float rowSumErrs[rowSums_len]={0.};

  //calculate sums
  std::cout<<"computing row sums..."<<std::endl;
  for(int rowNum=rowStart;rowNum<=nbins_y;rowNum++){
    float theSum=0.;
    float theSumErr=0.;
    for(int colNum=colStart;colNum<=nbins_x;colNum++){
      float theVal=inputTH2->GetBinContent(colNum,rowNum);
      float theValErr=inputTH2->GetBinError(colNum,rowNum);
      theSum+=theVal;
      theSumErr+=(theValErr*theValErr);
    }
    rowSums[rowNum]=theSum;
    rowSumErrs[rowNum]=std::sqrt(theSumErr);
  }
  std::cout<<std::endl;
  //scale each entry in each row down by that row's sum
  std::cout<<"scaling rows..."<<std::endl;
  for(int  rowNum=rowStart; rowNum<=nbins_y; rowNum++){

    float theSum=rowSums[rowNum];
    float theSumErr=rowSumErrs[rowNum];

    for( int colNum=colStart; colNum<=nbins_x; colNum++){

      float theVal=(inputTH2->GetBinContent(colNum,rowNum));
      float theValErr=(inputTH2->GetBinError(colNum,rowNum));

      float finalVal=theVal/theSum;
      float finalErr=(theVal/theSum)* std::sqrt(
                                                (theSumErr/theSum)*(theSumErr/theSum) +
                                                (theValErr/theVal)*(theValErr/theVal) )  ;


      inputTH2->SetBinContent(colNum,rowNum,finalVal);
      inputTH2->SetBinError(colNum,rowNum,finalErr);
    }
  }
  std::cout<<std::endl;

  inputTH2->SetEntries(numEntries);
  inputTH2->Print("base");

  std::cout<<std::endl<<"exiting normalizeMC_TH2"<<std::endl<<std::endl;
  return ;


}




TH2F* makeRespMatrixErrors( TH2F* inputTH2,			    
			    double* boundaries_pt_reco, int nbins_pt_reco,
			    double* boundaries_pt_gen  , int nbins_pt_gen           ){

  std::cout<<std::endl<<"in makeRespMatrixErrors"<<std::endl<<std::endl;
  bool funcDebug = true;
  
  //double numEntries=inputTH2->GetEntries();
  //inputTH2->Sumw2();
  if(funcDebug)inputTH2->Print("base");  std::cout<<std::endl;
  

  TAxis *xaxis = inputTH2->GetXaxis(); //reco pt axis
  int nbins_x= xaxis->GetNbins();
  int xbinstart=1;

  TAxis *yaxis = inputTH2->GetYaxis(); //gen pt axis
  int nbins_y = yaxis->GetNbins();
  int ybinstart=1;

  TH2F *errorTH2 = new TH2F("hmat_anabin_errors", "response matrix errors",
                               nbins_pt_reco, boundaries_pt_reco ,
                               nbins_pt_gen, boundaries_pt_gen );

  //if(doOverUnderflows){
  //  xbinstart--;
  //  ybinstart--;
  //  nbins_x++;
  //  nbins_y++;
  //}
  //reBinnedTH2->Sumw2();

  //fill each bin w/ sum of squares of errors from original matrix
  for (  int xbin=xbinstart ; xbin <= nbins_x ; xbin++ ) {
    for (int ybin=ybinstart ; ybin <= nbins_y ; ybin++ ) {
      
      errorTH2->Fill( xaxis->GetBinCenter(xbin) , yaxis->GetBinCenter(ybin) ,
		      (inputTH2->GetBinError(xbin,ybin)*inputTH2->GetBinError(xbin,ybin)) );
      
    }  }//end x-y loop

  
  //loop over each of the new bins for hmat_anabin, take the square root
  for (  int xbin=1 ; xbin <= nbins_pt_reco ; xbin++ ) {
    for (int ybin=1 ; ybin <= nbins_pt_gen ; ybin++ ) {
      
      float theVal=errorTH2->GetBinContent(xbin,ybin);
      theVal=std::sqrt(theVal);
      errorTH2->SetBinContent(xbin,ybin,theVal);
      
    }  }//end x-y loop
  
  //errorTH2->SetEntries(numEntries);

  if(funcDebug) std::cout<<"getSumW2N="<<errorTH2->GetSumw2N()<<std::endl;
  if(funcDebug) errorTH2->Print("base");  std::cout<<std::endl;


  std::cout<<std::endl<<"makeRespMatrixErrors Done"<<std::endl<<std::endl;
  return errorTH2;
}



TH2F* makeRespMatrixPercentErrs( TH2F* hmat_errors, TH2F* hmat_anabin,
				 double* boundaries_pt_reco, int nbins_pt_reco,
				 double* boundaries_pt_gen  , int nbins_pt_gen           ){

  
  std::cout<<std::endl<<"in makeRespMatrixPercentErrs"<<std::endl<<std::endl;
  bool funcDebug = true;
  

  int xbinstart=1;
  int ybinstart=1;
  
  TH2F *fracErrorTH2 = new TH2F("hmat_anabin_percenterrs", "response matrix %error",
				nbins_pt_reco, boundaries_pt_reco ,
				nbins_pt_gen, boundaries_pt_gen );
  
  //if(doOverUnderflows){
  //  xbinstart--;
  //  ybinstart--;
  //  nbins_x++;
  //  nbins_y++;
  //}
  //reBinnedTH2->Sumw2();
  
  //fill each bin w/ binerror/bincontent*100
  for (  int xbin=xbinstart ; xbin <= nbins_pt_reco ; xbin++ ) {
    for (int ybin=ybinstart ; ybin <= nbins_pt_gen ; ybin++ ) {

      float theVal=hmat_anabin->GetBinContent(xbin,ybin);
      float theErr=hmat_errors->GetBinContent(xbin,ybin);
      if(funcDebug)std::cout<<"for xbin="<<xbin<<" ybin="<<ybin<<std::endl;
      if(funcDebug)std::cout<<"theVal="<<theVal<<" +/- theErrn="<<theErr<<std::endl;

      //check NaN
      if( theVal!=theVal || theErr != theErr || theErr==0.) 
	{
	  if(funcDebug)std::cout<<"NaN or zero error!! Skipping xbin="<<xbin<<" ybin="<<ybin<<std::endl;
	  continue;
	}
      else {
	fracErrorTH2->SetBinContent( xbin ,  ybin ,
				     theErr/theVal*100.);      }
      
    }  }//end x-y loop
  
  
  
  //fracErrorTH2->SetEntries(numEntries);

  if(funcDebug) std::cout<<"getSumW2N="<<fracErrorTH2->GetSumw2N()<<std::endl;
  if(funcDebug) fracErrorTH2->Print("base");  std::cout<<std::endl;
  
  
  std::cout<<std::endl<<"makeRespMatrixPercentErrs done."<<std::endl<<std::endl;

  return fracErrorTH2;
}

void setRespMatrixErrs( TH2F* hmat_anabin , TH2F* hmat_errors ){   

  std::cout<<std::endl<<"in setRespMatrixErrors"<<std::endl<<std::endl;
  bool funcDebug = true;


  if(funcDebug) hmat_errors->Print("base");
  if(funcDebug) hmat_anabin->Print("base");

  //,
  //				 double* boundaries_pt_reco, int nbins_pt_reco,
  //				 double* boundaries_pt_gen  , int nbins_pt_gen           ){
  
  TAxis *xaxis_hmat = hmat_anabin->GetXaxis(); //reco pt axis
  int nbins_x= xaxis_hmat->GetNbins();
  int xbinstart=1;

  TAxis *yaxis_hmat = hmat_anabin->GetYaxis(); //gen pt axis
  int nbins_y = yaxis_hmat->GetNbins();
  int ybinstart=1;  

 
  TAxis *xaxis_hmaterr = hmat_errors->GetXaxis(); //reco pt axis
  int nbins_hmaterr_x= xaxis_hmaterr->GetNbins();

  TAxis *yaxis_hmaterr = hmat_errors->GetYaxis(); //gen pt axis
  int nbins_hmaterr_y = yaxis_hmaterr->GetNbins();

  assert(nbins_x==nbins_hmaterr_x);
  assert(nbins_y==nbins_hmaterr_y);


//  TH2F *fracErrorTH2 = new TH2F("hmat_anabin_percenterrs", "response matrix %error",
//				nbins_pt_reco, boundaries_pt_reco ,
//				nbins_pt_gen, boundaries_pt_gen );
  
  //if(doOverUnderflows){
  //  xbinstart--;
  //  ybinstart--;
  //  nbins_x++;
  //  nbins_y++;
  //}
  //reBinnedTH2->Sumw2();
  
  //fill each bin w/ binerror/bincontent*100
  for (  int xbin=xbinstart ; xbin <= nbins_x ; xbin++ ) {
    for (int ybin=ybinstart ; ybin <= nbins_y ; ybin++ ) {
      
      float theErr=hmat_errors->GetBinContent(xbin,ybin);
      hmat_anabin->SetBinError( xbin ,  ybin ,
				theErr );
      
    }  }//end x-y loop
  
  
  
  //fracErrorTH2->SetEntries(numEntries);

  if(funcDebug) std::cout<<"getSumW2N="<<hmat_errors->GetSumw2N()<<std::endl;
  if(funcDebug) hmat_errors->Print("base");  std::cout<<std::endl;
  
  
  std::cout<<std::endl<<"setRespMatrixErrs done."<<std::endl<<std::endl;

  return;
}




TH2F* reBinPearsonTH2(TMatrixD* pearson, const double* boundaries_pt, const int nbins_pt){
  std::cout<<std::endl<<"in reBinPearsonTH2"<<std::endl<<std::endl;

  bool debugPearson=false;
  //const double* the_ptBins=boundaries_pt;//ptbins
  //const int numbins=nbins_pt;//nbins

  TH2F *reBinnedTH2 = new TH2F(  "pearsonCustomBins" , "pearsonMatrix custom bins"   ,
                                 nbins_pt,boundaries_pt,
                                 nbins_pt,boundaries_pt);
  std::cout<<std::endl<<"in reBinPearsonTH2 function..."<<std::endl;
  reBinnedTH2->Print("base");

  for (int j=1;   j<=nbins_pt;j++) {
    for (int i=1; i<=nbins_pt;i++) {

      double pearson_i_j=(*pearson)(i-1,j-1);

      if(fabs(pearson_i_j) > 1.0) {
        if(debugPearson)std::cout <<std::endl<<"WARNING:: pearson["<<i-1<<"]["<<j-1<<"]="<<pearson_i_j<<std::endl<<std::endl;
        //      pearson_i_j=10.;
        //std::cout<< "pearson_i_j=10"<<std::endl
        continue; }
      else {
        if(debugPearson)std::cout <<"pearson["<<i-1<<"]["<<j-1<<"]="<<pearson_i_j<<std::endl;   }

      reBinnedTH2->Fill( i, j, fabs( (float) pearson_i_j ) );
    }
  }


  reBinnedTH2->Print("base");

  std::cout<<std::endl<<"reBinPearsonTH2 done"<<std::endl<<std::endl;
  return reBinnedTH2;
}



void TH2clearOverUnderflows(TH2F* h)
{
  //std::cout<<std::endl<<"WARNING!!!!"<<std::endl<<std::endl;
  //std::cout<<"doOverUnderflows="<<doOverUnderflows<<std::endl;
  std::cout<<"clearing Over/Underflow Bins for TH2...."<<std::endl;
  
  int nbinsx = h->GetXaxis()->GetNbins();
  int nbinsy = h->GetYaxis()->GetNbins();
  
  for(int i=0;i<=(nbinsx+1);i++){
    for(int j=0;j<=(nbinsy+1);j++){

      if( i==0          || 
	  j==0          || 
	  i==(nbinsx+1) || 
	  j==(nbinsy+1)    )	{
	  h->SetBinContent(i,j,0);
	  h->SetBinError(i,j,0);	}
      else continue;
    }
  }
  std::cout<<"done clearing Over/Underflow Bins for TH2...."<<std::endl;
  return;
}
