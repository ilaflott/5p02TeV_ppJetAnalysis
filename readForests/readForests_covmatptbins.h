

std::vector<std::vector<double>> merged_SMP_ptbins{
  {//00eta05

    ////if i'm allowed to omit some bins at high pt, then i like this one, because it works for 05y10 as well as this one.
    43., 49., //junk bins here+below
      56.,      
      64.,       74.,       84.,       97.,       114.,
      133.,      153.,      174.,      196.,      
      220.,      
      245.,
      272.,      
      300.,      
      330.,            
      //362.,   
      395.,    
      //430.,   
      //468.,      
      507.,   
      //548.,            
      //592.,     
      638.,    
      //686.,     
      //737.,     
      // 790.,
      846.,
      905.,//junk bins here+above 
      967.,      
      1032. ,
      1101.},

//    //what i get if i keep all bins w/ statistics and don't care about how it plays with other y bins
//    //con: don't like the way the boundaries contrast w/ 05y10
//    //con: plus the statistics in the last two bins are tricky
//    56.,      
//    64.,       74.,       84.,       97.,       114.,
//    133.,      153.,      174.,      196.,      
//      220.,     
//      245.,
//      272.,  
//      300., 
//      330.,   
//      362.,
//      //395.,    
//      430.,   
//      //468.,      
//      507.,  
//      //548.,
//      //592.,
//      638., 
//      //686.,
//      //737.,
//      790.,  
//      //846.,
//      //905.},
//      //967.},
//      1032.},

//    //another possibility what i get if i omit some bins w/ crappy statistics and don't care about how it plays with other y bins
////    //con: don't like the way the boundaries contrast w/ 05y10
////    //con: plus the statistics in the last two bins are tricky
//    56.,      
//    64.,       74.,       84.,       97.,       114.,
//    133.,      153.,      174.,      196.,      
//      220.,     
//      245.,
//      272.,  
//      300., 
//      330.,   
//      362.,
//      //395.,    
//      430.,   
//      //468.,      
//      507.,  
//      //548.,
//      592.,
//      //638., 
//      //686.,
//      737.,
//      //790.,  
//      //846.,
//      905.},
//    //967.},//last entry in data spectra
//    //1032.},





  {//05eta10



    //if i'm allowed to omit some bins at high pt, and i care about preserving the bin edges used in 00y05, this is what i get
    43., 49., //junk bins below here
      56.,      
      64.,       74.,       84.,       97.,       114.,
      133.,      153.,      174.,      196.,      
      220.,      
      245.,
      272.,      
      300.,      
      330.,            
      //362.,   
      395.,   
      //430.,   
      //468.,      
      507.,   
      //548.,            
      //592.,     
      638.,    
      //686.,     
      //737.,     
      // 790.,
      846.,//}, //WARNING empty bin in data 846. - 905.      
      905.,//junk bins here+above       //last entry in data spectra
      967.,      
      1032. ,
      1101.},
    






  {//10eta15
  
//    56.,      
//    64.,       74.,       84.,       97.,       114.,
//    133.,      153.,      174.,      196.,      
//      220.,      
//      245.,
//      272.,      
//      300.,     
//      330.,     
//      362.,     
//      395.,      
//      430.,
//      468.,      
//      507.,     
//      548., 
//      592.,
//      638.,
//      686.,  
//      737.},
    //,      790.,
    //846.}, 
//this is what i arrive at if i restrict myself to using bin edges present in 00y05 and 05y10
    43., 49., //junk bins here+below
    56.,      
    64.,       74.,       84.,       97.,       114.,
    133.,      153.,      174.,      196.,      
      220.,      
      245.,
      272.,      
      //300.,     
      330.,     
      //362.,     
      395.,      
      //430.,
      //468.,      
      507.,     
      //548., 
      //592.,
      638.,//},//,//last entry in data spectra      
      846.,//junk bins here+above 
      905.,
      967.,      
      1032. ,
      1101.},

  {//15eta20
    //43.,       49.,       
    43., 49., //junk bins here+below
    56.,      
    64.,       74.,       84.,       97.,       114.,
    133.,      153.,      174.,      196.,      220.,      245.,
      272., 
      //300.,   
      330.,     
      //362.,      
      395.,  
      //430.,
      //468.,      
      507.,//      
      548.,//junk bins here+above      
      592.,////empty, bin in data here     
      638.,//limit of NLO smearing            
      846.,
      905.,
      1032.,
      1101.}
      //686.}//last entry in data spectra
//  {//15eta20 JOHNS SUGGESTIION
//    //43.,       49.,       
//    //43., 49., 
//    56.,      
//    64.,       74.,       84.,       97.,       114.,
//    133.,      153.,      174.,      196.,      220.,      245.,
//      272., 
//      300.,   
//      //330.,     
//      362.,      
//      //395.,  
//      430.,
//      //468.,      
//      //507.}//,      
//      //548.,      
//      592.}//empty, bin in data here     
//      //638.}//limit of NLO smearing      
//      //686.}//last entry in data spectra

  //  {//20eta25
  //    43.,       49.,       56.,      
  //    64.,       74.,       84.,       97.,       114.,
  //    133.,      153.,      174.,      196.,      220.,      245.,
  //    272.,      300.,      330.,      362.,      395.,      430.,
  //    468.,      507.,      548.,      592.,      638.,      686.}
  
  
};









double* setBinning_etabin(int etabinint=0, int* nbins=NULL    ){
  bool funcDebug=true;
  if(funcDebug)std::cout<<std::endl<<"in setBinning_etabin"<<std::endl<<std::endl;  
  
  double* binarray=NULL;   
  binarray=(double*)(merged_SMP_ptbins[etabinint].data());
  *(nbins)=(int)(merged_SMP_ptbins[etabinint].size()-1);
  
  if(funcDebug)std::cout<<std::endl<<"setBinning_etabin done, exiting."<<std::endl<<std::endl;
  return binarray;  
  
}
