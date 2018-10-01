#include <cstdlib>
#include <cstring>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"


//void testIter() {
int main(int argc, char* argv[]){
  
  
  TFile *f1 = TFile::Open("fnl5020_LO2_R04Jets_modify_NNPDF30_nnlo_as_0121_MC_gaussSmear_00eta20.root");
  
  TIter next(f1->GetListOfKeys());
  TKey *key;
  
  while ( (key = (TKey*)next()) ) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1 *h = (TH1*)key->ReadObj();
    std::cout<<"found TH1 histogram! h->GetName()="<<h->GetName()<<std::endl;    
    std::cout<<"found TH1 histogram! h->Title()="<<h->GetTitle()<<std::endl;    
  }
  return 0;
}


////////////PSEUDOCODE FOR NEW MACRO////////////

// ---- terminal behavior ---- //
// compareTH1s(file1,file0)

// ---- compareTH1s pseudo code ---- //
// for(list of stuff in file 1)
//   get object from file 1
//   check it's a histo
//   h1=this object from file 1
//   get it's name and/or title
//   h2=NULL
//   for(list of stuff in file 2)
//     get object from file 2
//     check it's a histo
//     get it's name and/or title
//     if(file 1 object matches file 2 object) h2=this object from file 2 AND break
//     else continue
//   if(h2=NULL, still) cout<<"HIST NOT FOUND"
//   else compareANDprint(h1, h2)
// return

////////////PSEUDOCODE FOR NEW MACRO////////////
