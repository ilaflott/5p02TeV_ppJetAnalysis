

//// ------ SMP SPECTRA BINS ------
//const double anabins_pt_SMP[] = {//analysis pt bins

////1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 
//32, 37, 43, 49, 
//  56, 64, 74, 84,97, 114, 
//  133, 153, 174, 196, 220, 
//  245, 272, 300,   330, 
//  362, 395, 430, 
//  468,  507, 
//  548, 592, 638, 
//  686, 
//737, 790, 846,  905,   967,  1032 //, 
//1101, 1172, 1248, 1327, 
//  //1410, 1497, 1588, 1684, 1784, 1890, 2000,
//  //2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 
//  //4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000
//  };




//  //  0, 
//  //  3, 4, 5, 7, 9, 12, 15, 18, 21, 
//  28, // 24, 
//  30, // 28, 
//  32, 
//  37, 
//  43, 49, 56, 64, 74, 84, 97, 114, 
//  133, 153, 174, 196, 220, 245, 272, 300, 
//  330, 362, 395, 430, 468, 507, 548, 592,
//  638, 686 , 
//  800, 1000
//};
//const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;


//// ------ RECO SPECTRA BINS ------
//const double anabins_pt_reco[] = {//analysis pt bins
//  //  0, 
//  //  3, 4, 5, 7, 9, 12, 15, 18, 21, 
//  24, 
//  28,
//  32,
//  37, 
//  43, 49, 56, 64, 74, 84, 97, 114, 
//  133, 153, 174, 196, 220, 245, 272, 300, 
//  330, 362, 395, 430, 468, 507, 548, 592,
//  638, 686 , 
//  800, 1000
//};
//const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;



// for readForests jobs w/ 30/30 rec/gen pt cuts
// same as old binning... shifted up by 2 to respect the cuts at readforests-level
// ------ GEN SPECTRA BINS ------
const double anabins_pt_gen[] = {//analysis pt bins
  37.,
  43., 49., 56., 64., 74., 84., 97., 114., 
  133., 153., 174., 196., 220., 245., 272., 
  300., 330., 362., 395., 430., 468., 507., 
  548., 592., 638., 686., 737., 790., 846.,  
  905.,   967.,  1032., 1101., 1172. , 1248., 1327.
};
const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;

// ------ RECO SPECTRA BINS ------
const double anabins_pt_reco[] = {//analysis pt bins
  56., 64., 74., 84., 97., 114., 
  133., 153., 174., 196., 220., 245., 272., 
  300., 330., 362., 395., 430., 468., 507., 
  548., 592., 638., 686., 737., 790., 846.,  
  905., 967.,  1032., 1101.
};
const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;






//// for readForests jobs w/ 40/40 rec/gen pt cuts
//// same as old binning... shifted up by 2 to respect the cuts at readforests-level
//// ------ GEN SPECTRA BINS ------
//const double anabins_pt_gen[] = {//analysis pt bins
//  //  24, 
//  //30,  34,  
//  40,   46, 
//  52, 59, 67, 77, 87, 
//  100, 117, 136, 156, 
//  177, 199, 223, 248, 
//  275, 303, 333, 365, 
//  398, 433, 471, 510, 
//  551, 595, 641, 689 , 
//  739, 811, 900
//  //802,     1000
//};
//const int n_anabins_pt_gen = sizeof(anabins_pt_gen)/sizeof(double)-1;
//
//// ------ RECO SPECTRA BINS ------
//const double anabins_pt_reco[] = {//analysis pt bins
//  //  24, 
//  //30,  34,  
//  40,   46, 
//  52, 59, 67, 77, 87, 
//  100, 117, 136, 156, 
//  177, 199, 223, 248, 
//  275, 303, 333, 365, 
//  398, 433, 471, 510, 
//  551, 595, 641, 689 , 
//  739, 811, 900
//  //802,     1000
//};
//const int n_anabins_pt_reco = sizeof(anabins_pt_reco)/sizeof(double)-1;

