// const int ncen=10; 
// const int centbins[ncen+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}; 
// const char *cdir[ncen]  = {"010","1020","2030","3040","4050","5060","6070","7080","8090","90100"}; 
// const char *ccent[ncen] = {"0-10%","10-20%","20-30%","30-50%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"}; 

// const int ncen=4; 
// const int centbins[ncen+1] = {0, 20, 60, 100, 200}; 
// const char *cdir[ncen+1]  = {"010","1030","3050","50100", "PP"}; 
// const char *ccent[ncen+1] = {"0-10%","10-30%","30-50%","50-100%", "PP"}; 

// Adding PF candidates plots ----------------------------------------------------------------------------------------------------
// 1) 2D histograms for eta vs pT for the candidate types ( for PbPb this needs to be before and after Vs subtraction)
// 2) 1D histograms for candidate pT/ jet pT for candidates inside the jet radius, for each candidate type
// Particle::pdgId_ PFCandidate::particleId_
// PFCandidate::ParticleType Particle
// 0           0  X          unknown, or dummy
// +211, -211  1  h          charged hadron
// +11, -11    2  e          electron
// +13, -13    3  mu         muon
// 22          4  gamma      photon
// 130         5  h0         neutral hadron
// 130         6  h_HF       hadronic energy in an HF tower
// 22          7  egamma_HF  electromagnetic energy in an HF tower


// int findBin(int bin)
// {
//   int ibin=-1;
//   //! centrality is defined as 0.5% bins of cross section
//   //! in 0-200 bins
//   if(bin<20)ibin=0; //! 0-10%
//   else if(bin>=20  && bin<40 )ibin=1; //! 10-20%
//   else if(bin>=40  && bin<60 )ibin=2; //! 20-30%
//   else if(bin>=60  && bin<80 )ibin=3; //! 30-40%
//   else if(bin>=80  && bin<100 )ibin=4; //! 40-50%
//   else if(bin>=100  && bin<120 )ibin=5; //! 50-60%
//   else if(bin>=120  && bin<140 )ibin=6; //! 60-70%
//   else if(bin>=140  && bin<160 )ibin=7; //! 70-80%
//   else if(bin>=160  && bin<180 )ibin=8; //! 80-90%
//   else if(bin>=180  && bin<200 )ibin=9; //! 90-100%
//   return ibin;
// }

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
// const int no_radius = 3;
// const int list_radius[no_radius] = {2,3,4};


//// Jet Structure!
//struct Jet{
//  int id;
//  float pt;
//};
//bool compare_pt(Jet jet1, Jet jet2);
//bool compare_pt(Jet jet1, Jet jet2){
//  return jet1.pt > jet2.pt;
//}
//float LjCut = 100.0;
//float SbjCut = 40.0;
