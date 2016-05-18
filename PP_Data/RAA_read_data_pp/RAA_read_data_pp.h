#include <cstdlib>	
#include <TH1.h>	    

//////////////////////////////////////////////
////////// (initializa/declara)tions //////////
///////////////////////////////////////////////

// declare helper functions
double trigComb(bool *trg, int *pscl, double pt);
void divideBinWidth(TH1 *h);

// for convenience during testing, coding, etc.
const bool atMIT = true;
const std::string hadoopDir = "/mnt/hadoop/cms";
const std::string xrootdDirect = "root://cmsxrootd.fnal.gov/";

// pt binning
const int ptbins[] = {15, 30, 50, 80, 120, 170, 220, 300, 500};
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one

// eta binning
static const char *etaWidth = (char*)"20_eta_20";

// root file directories + tree names
const std::string treeNames[]={ 
  "GARBAGE ENTRY" , //"ak"+((std::string)radius)+jetType+"JetAnalyzer/t"  //Form("ak%d%sJetAnalyzer/t",radius, jetType.c_str()),
  "hiEvtAnalyzer/HiTree"  ,
  "skimanalysis/HltTree"  ,
  "hltanalysis/HltTree"   ,
  //"ppTrack"       ,
  "hltobject/HLT_AK4CaloJet40_Eta5p1_v"     , //hlt40
  "hltobject/HLT_AK4CaloJet60_Eta5p1_v"     , //hlt60
  "hltobject/HLT_AK4CaloJet80_Eta5p1_v"     , //hlt80
  "hltobject/HLT_AK4CaloJet100_Eta5p1_v"     //hlt100
}; 
const int N_trees = sizeof(treeNames)/sizeof(std::string);

// Jet variable names
const std::string var[] = {   
  "jtpt" ,  "rawpt",  "jteta", "jtphi", 
  "trkMax", "trkSum", "trkHardSum", 
  "chMax",  "chSum",  "chHardSum", 
  "phMax",  "phSum",  "phHardSum", 
  "neMax",  "neSum", 
  "eMax",   "eSum", 
  "muMax",  "muSum", 
  "Aj",     "xj" 
};

const int N_vars = sizeof(var)/sizeof(std::string);

//////////////////////////////////////
////////// helper functions //////////
//////////////////////////////////////

/////////////////////////
// divide by bin width //
/////////////////////////
void divideBinWidth(TH1 *h){

  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){//binsX loop 
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//end nbinsX loop

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  return;
}

/////////////////////////
// trigger combination //
/////////////////////////
double trigComb(bool *trg, int *pscl, double pt){

  double weight=0;
  if(trg[3] && pt>=100 )          weight = pscl[3];
  if(trg[2] && pt>=80  && pt<100) weight = pscl[2];
  if(trg[1] && pt>=60  && pt<80 ) weight = pscl[1];
  if(trg[0] && pt>=40  && pt<60 ) weight = pscl[0];
  
  return weight;
}
