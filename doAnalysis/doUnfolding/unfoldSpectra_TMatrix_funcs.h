TH2D* CalculatePearsonCoefficients(TMatrixD* covmat, bool debugPearson=false, std::string histName="") {
  
  if(debugPearson)std::cout<<std::endl<<"in CalculatePearsonCoefficients"<<std::endl<<std::endl;
  
  int startRowColInt=0, skipRowColInt=1;
  int NanCount=0, gt1Count=0, EntryCount=0;
  int nrows = covmat->GetNrows(), ncols = covmat->GetNcols();
  
  TMatrixD* pearsonCoefs= (TMatrixD*)covmat->Clone();
  TH2D * pearsonCoefs_TH2= new TH2D(histName.c_str(),"perason coefs", nrows, 0, nrows, ncols, 0, ncols);
  
  if(debugPearson)std::cout<<"looping over (nrows,ncols)=("<<nrows<<", "<<ncols<<")"<<std::endl;
  for(int row = startRowColInt; row<nrows; row++){ double covmatElement_rowrow = (*covmat)(row,row);
    for(int col = startRowColInt; col<ncols; col++){ double covmatElement_colcol = (*covmat)(col,col);
      double covmatElement_rowcol = (*covmat)(row,col);
      double pearson = covmatElement_rowcol/TMath::Sqrt(covmatElement_rowrow*covmatElement_colcol);
      
      // NaN protection
      bool notNan=true;
      if(pearson!=pearson){     pearson=-10.;
        if(debugPearson&&(row%skipRowColInt==0||col%skipRowColInt==0)) {
          //std::cout<<"for (row,col)= "<<row <<" , "<<col <<std::endl;
	  std::cout<<"NaN! pearson("<<row<<","<<col<<")=-10."<<std::endl;       }
        NanCount++;     notNan=false;    gt1Count++;  }
      
      // assign the value to the TMatrix
      (*pearsonCoefs)(row,col) = pearson;
      pearsonCoefs_TH2->SetBinContent(row+1, col+1, pearson);
      // debug
      
      if(notNan&&debugPearson&&abs(pearson)>1.) {
	std::cout<<"warning, abs(pearson)="<<abs(pearson)<<std::endl;
        gt1Count++;      }
      if( debugPearson&&(row%skipRowColInt==0&&col%skipRowColInt==0&&notNan) ) {
        //std::cout<<"for (row,col)= "<<row <<" , "<<col <<std::endl;
        //std::cout<<"covmat(row,row)="<< covmatElement_rowrow<<std::endl;
        //std::cout<<"covmat(col,col)="<< covmatElement_colcol<<std::endl;
        //std::cout<<"covmat(row,col)="<< covmatElement_rowcol<<std::endl;
        //std::cout<<"pearson="<<pearson <<std::endl;
	std::cout<<"pearsonCoefs("<<row<<","<<col<<")="<<( TMatrixD (*pearsonCoefs) )(row,col) <<std::endl;       }
      EntryCount++;
    } // end col loop
  } // end row loop

  if(debugPearson)std::cout<<"Entries in Pearson Matrix="<<EntryCount<<std::endl;
  if(debugPearson)std::cout<<"pearson>1 Entries="<<gt1Count<<std::endl;
  if(debugPearson)std::cout<<"Nan Entries="<<NanCount<<std::endl;
  
  if(debugPearson)std::cout<<std::endl<<"calculatePearsonCoefficients done"<<std::endl<<std::endl;
  return pearsonCoefs_TH2;
}



// PrintMatrix
//-----------------------------------------------------------------------------------------------------------------------
// Print the matrix as a table of elements.
// Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
// By default, format ="%11.4g", is used to print one element.
// One can specify an alternative format, e.g. format ="%6.2f  "
int PrintMatrix( const TMatrixD& m,
		 const char* format, const char* name, Int_t cols_per_sheet ){
  bool funcDebug=false;
  if(funcDebug)std::cout<<std::endl<<"in PrintMatrix"<<std::endl<<std::endl;

  if (!m.IsValid()) {
    m.Error("PrintMatrix","%s is invalid",name);
    return -1;
  }

  const Int_t ncols  = m.GetNcols();
  const Int_t nrows  = m.GetNrows();
  const Int_t collwb = m.GetColLwb();
  const Int_t rowlwb = m.GetRowLwb();

  if (!(format && format[0])) format= "%11.4g ";
  char topbar[1000];
  snprintf(topbar,1000,format,123.456789);
  Int_t nch = strlen(topbar)+1;
  if (nch > 18) nch = 18;
  char ftopbar[20];
  for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
  Int_t nk = 1 + Int_t(log10(ncols));
  snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
  Int_t nch2 = strlen(ftopbar);
  for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
  ftopbar[nch] = '|';
  ftopbar[nch+1] = 0;

  printf("\n%dx%d %s is as follows",nrows,ncols,name);

  if (cols_per_sheet <= 0) {
    cols_per_sheet = 5;
    if (nch <= 8) cols_per_sheet =10;
  }
  nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
  for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
  topbar[nk] = 0;
  for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
    printf("\n\n     |");
    for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
      printf(ftopbar,j+collwb-1);
    printf("\n%s\n",topbar);
    if (m.GetNoElements() <= 0) continue;
    for (Int_t i=1 ; i<=nrows ; i++) {
      printf("%4d |",i+rowlwb-1);
      for( Int_t j=sheet_counter ; j<sheet_counter+cols_per_sheet && j<=ncols ; j++) printf(format,m(i+rowlwb-1,j+collwb-1));
      printf("\n");
    }
  }
  printf("\n");
  if(funcDebug)std::cout<<std::endl<<"PrintMatrix done"<<std::endl<<std::endl;
  return 0;
} // end PrintMatrix








// old pearson matrix routine, don't delete, 5/15/2018
//TMatrixD* CalculatePearsonCoefficients(TMatrixD* covmat, bool debugPearson=false, std::string cloneName="") {
//  //TMatrixD CalculatePearsonCoefficients(TMatrixD* covmat, bool debugPearson=false, std::string cloneName="") {
//  
//  if(debugPearson)std::cout<<std::endl<<"in CalculatePearsonCoefficients"<<std::endl<<std::endl;
//  
//  int startRowColInt=0, skipRowColInt=1;
//  int NanCount=0, gt1Count=0, EntryCount=0;
//  int nrows = covmat->GetNrows(), ncols = covmat->GetNcols();
//  
//  TMatrixD* pearsonCoefs= (TMatrixD*)covmat->Clone(cloneName.c_str());
//  //TMatrixD* pearsonCoefs= (TMatrixD*)covmat->Clone();
//
//  
//  if(debugPearson)std::cout<<"looping over (nrows,ncols)=("<<nrows<<", "<<ncols<<")"<<std::endl;
//  for(int row = startRowColInt; row<nrows; row++){ double covmatElement_rowrow = (*covmat)(row,row);
//    for(int col = startRowColInt; col<ncols; col++){ double covmatElement_colcol = (*covmat)(col,col);
//      double covmatElement_rowcol = (*covmat)(row,col);
//      double pearson = covmatElement_rowcol/TMath::Sqrt(covmatElement_rowrow*covmatElement_colcol);
//      
//      // NaN protection
//      bool notNan=true;
//      if(pearson!=pearson){     pearson=-10.;
//        if(debugPearson&&(row%skipRowColInt==0||col%skipRowColInt==0)) {
//          //std::cout<<"for (row,col)= "<<row <<" , "<<col <<std::endl;
//	  std::cout<<"NaN! pearson("<<row<<","<<col<<")=-10."<<std::endl;       }
//        NanCount++;     notNan=false;    gt1Count++;  }
//
//      // assign the value to the TMatrix
//      (*pearsonCoefs)(row,col) = pearson;
//
//      // debug
//
//      if(notNan&&debugPearson&&abs(pearson)>1.) {
//	std::cout<<"warning, abs(pearson)="<<abs(pearson)<<std::endl;
//        gt1Count++;      }
//      if( debugPearson&&(row%skipRowColInt==0&&col%skipRowColInt==0&&notNan) ) {
//        //std::cout<<"for (row,col)= "<<row <<" , "<<col <<std::endl;
//        //std::cout<<"covmat(row,row)="<< covmatElement_rowrow<<std::endl;
//        //std::cout<<"covmat(col,col)="<< covmatElement_colcol<<std::endl;
//        //std::cout<<"covmat(row,col)="<< covmatElement_rowcol<<std::endl;
//        //std::cout<<"pearson="<<pearson <<std::endl;
//	std::cout<<"pearsonCoefs("<<row<<","<<col<<")="<<( TMatrixD (*pearsonCoefs) )(row,col) <<std::endl;       }
//      EntryCount++;
//    } // end col loop
//  } // end row loop
//
//  if(debugPearson)std::cout<<"Entries in Pearson Matrix="<<EntryCount<<std::endl;
//  if(debugPearson)std::cout<<"pearson>1 Entries="<<gt1Count<<std::endl;
//  if(debugPearson)std::cout<<"Nan Entries="<<NanCount<<std::endl;
//
//  if(debugPearson)std::cout<<std::endl<<"calculatePearsonCoefficients done"<<std::endl<<std::endl;
//  return pearsonCoefs;
//  //return ((TMatrixD)(*pearsonCoefs));
//}



