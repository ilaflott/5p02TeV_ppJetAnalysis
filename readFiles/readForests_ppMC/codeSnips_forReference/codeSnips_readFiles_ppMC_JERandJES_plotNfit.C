// old file input

//if(ic <  ncen)  hrsp[nj][ip] = (TH1F*)fin->Get(Form("hJER_%d_pt_%d_cent%d",ptbins[ip], ptbins[ip+1],ic));// PbPb
//if(ic == ncen) hrsp[nj][ip] = (TH1F*)finPP->Get(Form("hJER_%d_pt_%d",ptbins[ip], ptbins[ip+1])); // pp



// for gaussian fit of inputJER hists in ptbins, partially related to PbPb

//fgaus = new TF1("fgaus","gaus", mean - 2.15*hrsp[nj][ip]->GetRMS(), mean + 2.15*hrsp[nj][ip]->GetRMS());
//fgaus->SetParameters(norm, 1.00, 0.6);
//fgaus->SetParLimits(1,0.60,1.30);
//fgaus->SetParLimits(2,0.20,0.80);

//fgaus = new TF1("fgaus","gaus", 0.80, 1.20);
// fgaus->SetParameters(norm, 0.9999, 0.5);
// fgaus->SetParLimits(1,0.92,1.05);
// fgaus->SetParLimits(2,0.20,0.80);

// if(ic==ncen-1){
//   fgaus = new TF1("fgaus","gaus", 0.90, 1.10);
//   fgaus->SetParameters(norm, 0.9999, 0.1);
//   fgaus->SetParLimits(1,0.95,1.05);
//   fgaus->SetParLimits(2,0.20,0.80);
// }else{
//   fgaus = new TF1("fgaus","gaus", 0.80, 1.20);
//   fgaus->SetParameters(norm, 0.9999, 0.5);
//   fgaus->SetParLimits(1,0.90,1.10);
//   fgaus->SetParLimits(2,0.20,0.80);
// }

// post gaussian fit of input JET hists, PbPb related i think

//if(ic == ncen-1)fitstatus=-1;
// if(nj==0 && ic==ncen-1){
//   std::cout << "*******START*************** "<< ccent << " " << ptbins[ip]  << " to " << ptbins[ip+1] << std::endl;
//   fitstatus = hrsp[nj][ip]->Fit(fgaus,"RLQ");
//   std::cout <<" \t Fit Mean : " << fgaus->GetParameter(1) << " Hist Mean : " << mean << std::endl;
//   std::cout << "*******END***************** "<< ccent << " " << ptbins[ip]  << " to " << ptbins[ip+1] << std::endl;
//   std::cout << std::endl;
// }
// else fitstatus = hrsp[nj][ip]->Fit(fgaus,"RLQ");



// older implementation of fit functions, no longer used?

//FitDist(hrsp[nj][ip], mean, emean, sig, esig);
//FitGauss(hrsp[nj][ip], mean, emean, sig, esig);
