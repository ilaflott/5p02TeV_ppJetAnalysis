#include <cmath>

TH2D* reBinTH2( TH2D* inputTH2, std::string rebinTH2_name,
		double* boundaries_pt_reco, int nbins_pt_reco,
		double* boundaries_pt_gen  , int nbins_pt_gen           ){
  
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in reBinTH2"<<std::endl<<std::endl;
  
  //double numEntries=inputTH2->GetEntries();
  //inputTH2->Sumw2();
  if(funcDebug)inputTH2->Print("base");  std::cout<<std::endl;
  
  
  TAxis *xaxis = inputTH2->GetXaxis(); //reco pt axis
  int nbins_x= xaxis->GetNbins();
  int xbinstart=1;
  
  TAxis *yaxis = inputTH2->GetYaxis(); //gen pt axis
  int nbins_y = yaxis->GetNbins();
  int ybinstart=1;

  TH2D *reBinnedTH2 = new TH2D(rebinTH2_name.c_str(), inputTH2->GetTitle(),
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


  if(funcDebug)std::cout<<std::endl<<"reBinTH2 done"<<std::endl<<std::endl;
  return reBinnedTH2;
}





void divideBinWidth_TH2(TH2D *h){
  
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in divideBinWidth_TH2"<<std::endl<<std::endl;
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

      Double_t xWidth=h->GetXaxis()->GetBinWidth(i);
      Double_t yWidth=h->GetYaxis()->GetBinWidth(j);



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


      Double_t val = (h->GetBinContent(i,j))/(xWidth*yWidth);
      //Double_t valErr = (h->GetBinError(i,j))/(xWidth*yWidth);

      h->SetBinContent(i,j,val);
      //h->SetBinError(i,j,valErr);
      h->SetBinError(i,j,0.0);

      if(funcDebug)std::cout<<std::endl<<"new val = "   <<      h->GetBinContent(i,j) << std::endl;
      if(funcDebug)std::cout<<"new valErr = "<<      h->GetBinError(i,j)   << std::endl;

    }
  }

  //h->SetEntries(numEntries);
  if(funcDebug)h->Print("base");      std::cout<<std::endl;

  if(funcDebug)std::cout<<std::endl<<"divideBinWidth_TH2 done"<<std::endl<<std::endl;
  return;
}

//normalizes MC response TH2 columns to the sum of each columen, w/ error appropriately calculated.

TH2D* normalizeCol_RespMatrix( TH2D* inputTH2,
			       double* boundaries_pt_reco_mat, int nbins_pt_reco_mat,
			       double* boundaries_pt_gen_mat, int nbins_pt_gen_mat  ){
  
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in normalizeCol_RespMatix"<<std::endl<<std::endl;
  
  //inputTH2->Sumw2();
  inputTH2->Print("base");  
  if(funcDebug)std::cout<<std::endl;
  TH2D* colnormd_TH2 = new TH2D( "hmat_anabin_colnorm","RespMatrix, columns normalized to column sum",
				 nbins_pt_reco_mat, boundaries_pt_reco_mat ,
				 nbins_pt_gen_mat, boundaries_pt_gen_mat );
  
  
  //double numEntries=inputTH2->GetEntries();
  
  int xstart=1;
  TAxis *xaxis = inputTH2->GetXaxis();
  int nbins_x= xaxis->GetNbins();
  
  int ystart=1;
  TAxis *yaxis = inputTH2->GetYaxis();
  int nbins_y = yaxis->GetNbins();
  
  const int colSums_len=nbins_x;
  double colSums[colSums_len]={0.};
  double colSumErrs[colSums_len]={0.};
  
  //calculate sums
  if(funcDebug)std::cout<<"computing column sums st. each col val sums to 1"<<std::endl;
  if(funcDebug)std::cout<<"so val in cell of a given column, gives prob. a reco jet in bin i was gen'd from a jet from bin j"<<std::endl;

  for(int colNum=xstart;colNum<=nbins_x;colNum++){
    double theSum=0.;
    double theSumErr=0.;

    for(int rowNum=ystart;rowNum<=nbins_y;rowNum++){
      double theVal=inputTH2->GetBinContent(colNum,rowNum);
      double theValErr=inputTH2->GetBinError(colNum,rowNum);
      theSum+=theVal;
      theSumErr+=(theValErr*theValErr);
    }
    
    colSums[colNum]=theSum;
    colSumErrs[colNum]=std::sqrt(theSumErr);
  }
  if(funcDebug)std::cout<<std::endl;
  
  //scale each entry in each row down by that row's sum
  if(funcDebug) std::cout<<"scaling entries in each bin according to respective column sum..."<<std::endl;
  for(int  colNum=xstart; colNum<=nbins_x; colNum++){
    
    double theSum=colSums[colNum];
    double theSumErr=colSumErrs[colNum];
    
    for( int rowNum=ystart; rowNum<=nbins_y; rowNum++){
      
      double theVal=(inputTH2->GetBinContent(colNum,rowNum));
      double theValErr=(inputTH2->GetBinError(colNum,rowNum));
      
      double finalVal=theVal/theSum;
      double finalErr=(theVal/theSum)* std::sqrt(
                                                (theSumErr/theSum)*(theSumErr/theSum) +
                                                (theValErr/theVal)*(theValErr/theVal) )  ;
      if(finalVal!=finalVal){	 
	if(funcDebug)std::cout<<"NaN! skipping setbincontent/err"<<std::endl;
	if(funcDebug)std::cout<<"(colNum,rowNum)=("<<colNum<<","<<rowNum<<")"<<std::endl;
	continue;      }
      
      colnormd_TH2->SetBinContent(colNum,rowNum,finalVal);
      colnormd_TH2->SetBinError(colNum,rowNum,finalErr);
    }
  }
  
  if(funcDebug)std::cout<<std::endl;
  
  if(funcDebug) inputTH2->Print("base");
  
  //std::cout<<std::endl<<"exiting normalizeMC_TH2"<<std::endl<<std::endl;
  if(funcDebug)std::cout<<std::endl<<"exiting normalizeCol_RespMatrix"<<std::endl<<std::endl;
  return colnormd_TH2;
  
  
}




//normalizes MC response TH2 columns to the sum of each columen, w/ error appropriately calculated.

TH2D* normalizeRow_RespMatrix( TH2D* inputTH2,
			       double* boundaries_pt_reco_mat, int nbins_pt_reco_mat,
			       double* boundaries_pt_gen_mat, int nbins_pt_gen_mat  ){
  
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in normalizeRow_RespMatrix"<<std::endl<<std::endl;
  
  //inputTH2->Sumw2();
  inputTH2->Print("base");  std::cout<<std::endl;
  
  TH2D* rownormd_TH2 = new TH2D( "hmat_anabin_rownorm","RespMatrix, rows normalized to row sum",
				 nbins_pt_reco_mat, boundaries_pt_reco_mat ,
				 nbins_pt_gen_mat, boundaries_pt_gen_mat );
  

  //double numEntries=inputTH2->GetEntries();

  int xstart=1;
  TAxis *xaxis = inputTH2->GetXaxis();
  int nbins_x= xaxis->GetNbins();

  int ystart=1;
  TAxis *yaxis = inputTH2->GetYaxis();
  int nbins_y = yaxis->GetNbins();
  
  const int rowSums_len=nbins_y;
  double rowSums[rowSums_len]={0.};
  double rowSumErrs[rowSums_len]={0.};

  //calculate sums
  if(funcDebug)std::cout<<"computing row sums st. each row sums to 1"<<std::endl;
  if(funcDebug)std::cout<<"so val in cell of a given row, gives prob. a gen jet in bin i get's reco'd in bin j"<<std::endl;

  for(int rowNum=ystart;rowNum<=nbins_y;rowNum++){
    double theSum=0.;
    double theSumErr=0.;
    for(int colNum=xstart;colNum<=nbins_x;colNum++){
      
      double theVal=inputTH2->GetBinContent(colNum,rowNum);
      double theValErr=inputTH2->GetBinError(colNum,rowNum);
      theSum+=theVal;
      theSumErr+=(theValErr*theValErr);
    }
    
    rowSums[rowNum]=theSum;
    rowSumErrs[rowNum]=std::sqrt(theSumErr);
  }
  if(funcDebug)std::cout<<std::endl;
  
  //scale each entry in each row down by that row's sum
  if(funcDebug) std::cout<<"scaling entries in each bin according to respective row sum..."<<std::endl;
  for( int rowNum=ystart; rowNum<=nbins_y; rowNum++){
    
    double theSum=rowSums[rowNum];
    double theSumErr=rowSumErrs[rowNum];
    
    for(int  colNum=xstart; colNum<=nbins_x; colNum++){
      
      double theVal=(inputTH2->GetBinContent(colNum,rowNum));
      double theValErr=(inputTH2->GetBinError(colNum,rowNum));
      
      double finalVal=theVal/theSum;
      double finalErr=(theVal/theSum)* std::sqrt(
                                                (theSumErr/theSum)*(theSumErr/theSum) +
                                                (theValErr/theVal)*(theValErr/theVal) )  ;

      if(finalVal!=finalVal){	 
	if(funcDebug)std::cout<<"NaN! skipping setbincontent/err"<<std::endl;
	if(funcDebug)std::cout<<"(colNum,rowNum)=("<<colNum<<","<<rowNum<<")"<<std::endl;
	continue;      }
      rownormd_TH2->SetBinContent(colNum,rowNum,finalVal);
      rownormd_TH2->SetBinError(colNum,rowNum,finalErr);
    }
  }
  
  if(funcDebug)std::cout<<std::endl;
  
  if(funcDebug) inputTH2->Print("base");
  
  
  if(funcDebug)std::cout<<std::endl<<"exiting normalizeRow_RespMatrix"<<std::endl<<std::endl;
  return rownormd_TH2;
  
}




TH2D* makeRespMatrixErrors( TH2D* inputTH2,			    
			    double* boundaries_pt_reco, int nbins_pt_reco,
			    double* boundaries_pt_gen  , int nbins_pt_gen , std::string respmat_errs_name="hmat_anabin_errors"          ){


  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in makeRespMatrixErrors"<<std::endl<<std::endl;  
  //double numEntries=inputTH2->GetEntries();
  //inputTH2->Sumw2();
  if(funcDebug)inputTH2->Print("base");  std::cout<<std::endl;
  

  TAxis *xaxis = inputTH2->GetXaxis(); //reco pt axis
  int nbins_x= xaxis->GetNbins();
  int xbinstart=1;

  TAxis *yaxis = inputTH2->GetYaxis(); //gen pt axis
  int nbins_y = yaxis->GetNbins();
  int ybinstart=1;

  TH2D *errorTH2 = new TH2D(respmat_errs_name.c_str(), "response matrix errors",
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
      
      double theVal=errorTH2->GetBinContent(xbin,ybin);
      theVal=std::sqrt(theVal);
      errorTH2->SetBinContent(xbin,ybin,theVal);
      
    }  }//end x-y loop
  
  //errorTH2->SetEntries(numEntries);

  if(funcDebug) std::cout<<"getSumW2N="<<errorTH2->GetSumw2N()<<std::endl;
  if(funcDebug) errorTH2->Print("base");  std::cout<<std::endl;


  if(funcDebug)std::cout<<std::endl<<"makeRespMatrixErrors Done"<<std::endl<<std::endl;
  return errorTH2;
}



TH2D* makeRespMatrixPercentErrs( TH2D* hmat_errors, TH2D* hmat_anabin,
				 double* boundaries_pt_reco, int nbins_pt_reco,
				 double* boundaries_pt_gen  , int nbins_pt_gen           ){

  
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in makeRespMatrixPercentErrs"<<std::endl<<std::endl;
  

  int xbinstart=1;
  int ybinstart=1;
  
  
  //std::cout<<"memleak here?!"<<std::endl;
  TH2D *fracErrorTH2 = new TH2D("hmat_anabin_percenterrs", "response matrix %error",
				nbins_pt_reco, boundaries_pt_reco ,
				nbins_pt_gen, boundaries_pt_gen );
  //std::cout<<"memleak here in makeRespMatrixPercentErrors"<<std::endl;
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

      double theVal=hmat_anabin->GetBinContent(xbin,ybin);
      double theErr=hmat_errors->GetBinContent(xbin,ybin);
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
  
  
  if(funcDebug)std::cout<<std::endl<<"makeRespMatrixPercentErrs done."<<std::endl<<std::endl;

  return fracErrorTH2;
}

void setRespMatrixErrs( TH2D* hmat_anabin , TH2D* hmat_errors , bool zeroBins){   
  
  
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in setRespMatrixErrors"<<std::endl<<std::endl;    
  //if(funcDebug) hmat_errors->Print("base");
  //if(funcDebug) hmat_anabin->Print("base");

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


  //  TH2D *fracErrorTH2 = new TH2D("hmat_anabin_percenterrs", "response matrix %error",
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
  int binsZerod=0;
  for (  int xbin=xbinstart ; xbin <= nbins_x ; xbin++ ) {
    for (int ybin=ybinstart ; ybin <= nbins_y ; ybin++ ) {
      
      double theErr=hmat_errors->GetBinContent(xbin,ybin);
      double theVal=hmat_anabin->GetBinContent(xbin,ybin);
      if(theErr>theVal){
	if(funcDebug)std::cout<<"WARNING error greater than bin value. setting error to bin value."<<std::endl;
	if(funcDebug)std::cout<<"theErr is = "<< theErr << std::endl;
	if(funcDebug)std::cout<<"theVal is = "<< theVal<< std::endl;
	if(funcDebug)std::cout<<"theDiff is = "<< fabs(theVal-theErr)<< std::endl;
	hmat_anabin->SetBinError( xbin ,  ybin , theVal );      }
      else if (theErr==theVal && theVal>0.){
	if(zeroBins){
	  binsZerod++;
	  if(funcDebug)std::cout<<"error equal to bin value. setting content/err to zero!"<<std::endl;
	  hmat_anabin->SetBinContent( xbin ,  ybin , 0. );
	  hmat_anabin->SetBinError( xbin ,  ybin , 0. );	}
	else{
	  if(funcDebug)std::cout<<"setting error."<<std::endl;
	  hmat_anabin->SetBinError( xbin ,  ybin , theErr );
	}
      }
      else{
	if(funcDebug)std::cout<<"setting error."<<std::endl;
	hmat_anabin->SetBinError( xbin ,  ybin , theErr );
      }
      
    }  }//end x-y loop
  
  
  
  //fracErrorTH2->SetEntries(numEntries);

  if(funcDebug) std::cout<<"getSumW2N="<<hmat_errors->GetSumw2N()<<std::endl;
  if(funcDebug) std::cout<<"binsZerod="<<binsZerod<<" out of "<< nbins_x*nbins_y <<" possible bins"<<std::endl;
  if(funcDebug) hmat_errors->Print("base");  std::cout<<std::endl;
  
  
  if(funcDebug)std::cout<<std::endl<<"setRespMatrixErrs done."<<std::endl<<std::endl;

  return;
}




TH2D* reBinPearsonTH2(TMatrixD* pearson, const double* boundaries_pt, const int nbins_pt){

  bool debugPearson=false;
  if(debugPearson)std::cout<<std::endl<<"in reBinPearsonTH2"<<std::endl<<std::endl;
  //const double* the_ptBins=boundaries_pt;//ptbins
  //const int numbins=nbins_pt;//nbins

  TH2D *reBinnedTH2 = new TH2D(  "pearsonCustomBins" , "pearsonMatrix custom bins"   ,
                                 nbins_pt,boundaries_pt,
                                 nbins_pt,boundaries_pt);
  //std::cout<<std::endl<<"in reBinPearsonTH2 function..."<<std::endl;
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

      reBinnedTH2->Fill( i, j, fabs( (double) pearson_i_j ) );
    }
  }


  reBinnedTH2->Print("base");

  if(debugPearson)std::cout<<std::endl<<"reBinPearsonTH2 done"<<std::endl<<std::endl;
  return reBinnedTH2;
}



void TH2clearOverUnderflows(TH2D* h)
{
  bool funcDebug=false;
 //std::cout<<std::endl<<"WARNING!!!!"<<std::endl<<std::endl;
  //std::cout<<"doOverUnderflows="<<doOverUnderflows<<std::endl;
  if(funcDebug)std::cout<<"clearing Over/Underflow Bins for TH2...."<<std::endl;
  
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
  if(funcDebug)std::cout<<"done clearing Over/Underflow Bins for TH2...."<<std::endl;
  return;
}
