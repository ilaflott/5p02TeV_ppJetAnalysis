#include <cstdlib>	
#include <TH1.h>	    

//////////////////////////////////////////////
////////// (initializa/declara)tions //////////
///////////////////////////////////////////////

// declare helper functions
double trigComb(bool *trg, int *pscl, double pt);
void divideBinWidth(TH1 *h);

// for convenience during testing, coding, etc.
const bool atMIT  = true;
const std::string hadoopDir    = "/mnt/hadoop/cms";
const std::string xrootdDirect = "root://cmsxrootd.fnal.gov/";

// pt binning
const int ptbins[] = {15, 30, 50, 80, 120, 170, 220, 300, 500};
const int nbins_pt = sizeof(ptbins)/sizeof(int)-1;//above values define edges of bins, not centers, so subtract one

// eta binning
static const char *etaWidth = (char*)"20_eta_20";

// root file directories + tree names
const std::string defaultDirNames[]={ "hltanalysis"   ,
				      "skimanalysis"  ,
				      "GARBAGE ENTRY" , //later replaced with Form("ak%d%sJetAnalyzer",radius, jetType.c_str()),
				      //"ppTrack"       ,
				      "hiEvtAnalyzer" ,
				      "hltobject"     , //hlt40
				      "hltobject"     , //hlt60
				      "hltobject"     , //hlt80
				      "hltobject"     }; //hlt100

const std::string trees[]={ "HltTree"                    ,
			    "HltTree"                    ,
			    "t"                          , 
			    //"trackTree"                  ,
			    "HiTree"                     ,
			    "HLT_AK4CaloJet40_Eta5p1_v"  ,
			    "HLT_AK4CaloJet60_Eta5p1_v"  ,
			    "HLT_AK4CaloJet80_Eta5p1_v"  ,
			    "HLT_AK4CaloJet100_Eta5p1_v" };
const int N_trees = sizeof(defaultDirNames)/sizeof(std::string);
                                                                
// Jet variable names
const std::string var[] = {"jtpt" , "rawpt", "jteta"    , "jtphi", "trkMax", "trkSum"   , "trkHardSum" , 
			   "chMax", "chSum", "chHardSum", "phMax", "phSum" , "phHardSum", "neMax"      , 
			   "neSum", "eMax" , "eSum"     , "muMax", "muSum" , "Aj"       , "xj"         };
const int N_vars = sizeof(var)/sizeof(std::string);



//////////////////////////////////////
////////// helper functions //////////
//////////////////////////////////////


/////////////////////////
// divide by bin width //
/////////////////////////
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
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
