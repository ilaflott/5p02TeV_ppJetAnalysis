#include "readFiles_ppMC_plotJESandJER.h"

// the macro
//int plot_JESandJERClosure(int wJetID=1, std::string jetType = "PF"){
int readFiles_ppMC_plotJESandJER(){

  int iSave=1;
  int wJetID=1;  int id = wJetID;
  std::string jetType = "PF";
  std::cout <<" npt : " << nbins_pt << std::endl;

  LoadStyle();  
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  TFile *finPP = new TFile(Form("PP_5p02TeV_MC_ak%s_20_eta_20.root",(srad[0]+jetType).c_str()),"r");
  TH1F *hrsp [knj][nbins_pt];
  TH1F *hMean[knj], *hSigma[knj];
  TF1 *fgaus=0;
  for(int nj=0; nj<knj; nj++){
    std::cout<<"nj = "<<nj<<std::endl;
    
    std::string dirname="ak"+srad[nj]+jetType+"JetAnalyzer";
    std::string algname="ak"+srad[nj]+jetType;
    
    hMean[nj] = new TH1F(Form("hMean%d",nj),Form("Mean %s %s",algname.c_str(),ccent),nbins_pt,ptbins_bound);
    hMean[nj]->SetMarkerColor(1);
    hMean[nj]->SetMarkerStyle(20);
    hMean[nj]->SetLineColor(1);
    hMean[nj]->SetMarkerSize(1.3);
    //MakeHistMean(hMean[nj],1.022,0.978); 
    //MakeHistMean(hMean[nj],1.082,0.898); 
    //MakeHistMean(hMean[nj],1.148,0.848); 
    MakeHistMean(hMean[nj],1.082,0.858); 
    //MakeHistMean(hMean[nj],1.052,0.934); 
    
    hSigma[nj] = new TH1F(Form("hSigma%d_%d",nj,ic),Form("Sigma %s %s",algname.c_str(),ccent),nbins_pt,ptbins_bound);
    hSigma[nj]->SetMarkerColor(1);
    hSigma[nj]->SetMarkerStyle(20);
    hSigma[nj]->SetLineColor(1);
    hSigma[nj]->SetMarkerSize(1.3);
    MakeHistRMS(hSigma[nj],0.563,0.001); 
    
    for(int ip=0; ip<nbins_pt; ip++){
      std::cout<<"ip = "<<ip<<std::endl;
      //if(ic <  ncen)  hrsp[nj][ip] = (TH1F*)fin->Get(Form("hJER_%d_pt_%d_cent%d",ptbins[ip], ptbins[ip+1],ic));
      if(ic == ncen) hrsp[nj][ip] = (TH1F*)finPP->Get(Form("hJER_%d_pt_%d",ptbins[ip], ptbins[ip+1]));
      
      hrsp[nj][ip]->Print("base");
      hrsp[nj][ip]->Scale(1./hrsp[nj][ip]->Integral());
      hrsp[nj][ip]->Rebin(3);
      
      double norm  = hrsp[nj][ip]->GetMaximumStored();
      double err   = hrsp[nj][ip]->GetStdDevError();
      double mean  = hrsp[nj][ip]->GetMean();
      double emean = hrsp[nj][ip]->GetMeanError();
      double sig   = hrsp[nj][ip]->GetStdDev()/mean;

      double esig = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));
            
      //! R 0.3
      //fgaus = new TF1("fgaus","gaus", 0.50,1.50);
      //fgaus->SetParameters(norm, 1.00, 0.3);            
      //! R 0.2
      fgaus = new TF1("fgaus","gaus", 0.50,1.50);
      fgaus->SetParameters(norm, 0.9999, 0.5);
            
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
      
      int fitstatus = 0;
      fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");
      
      // if(nj==0 && ic==ncen-1){
      //   std::cout << "*******START*************** "<< ccent << " " << ptbins[ip]  << " to " << ptbins[ip+1] << std::endl;
      //   fitstatus = hrsp[nj][ip]->Fit(fgaus,"RLQ");
      //   std::cout <<" \t Fit Mean : " << fgaus->GetParameter(1) << " Hist Mean : " << mean << std::endl;
      //   std::cout << "*******END***************** "<< ccent << " " << ptbins[ip]  << " to " << ptbins[ip+1] << std::endl;
      //   std::cout << std::endl;
      // }
      // else fitstatus = hrsp[nj][ip]->Fit(fgaus,"RLQ");
      
      std::cout << " fitstatus : "<< fitstatus << " Fit Error : " << fgaus->GetParError(1) << " Hist Error " << hrsp[nj][ip]->GetMeanError() << std::endl;
      
      //if(ic == ncen-1)fitstatus=-1;
      // mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     :  0.5*(fgaus->GetParameter(1) + hrsp[nj][ip]->GetMean());
      // emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): sqrt(pow(fgaus->GetParError(1),2) + pow(hrsp[nj][ip]->GetMeanError(),2));

      mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     : fgaus->GetParameter(1);
      emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): fgaus->GetParError(1);//hrsp[nj][ip]->GetMeanError();
      sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
      esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetStdDevError(),2))+(pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetStdDevError(),2))+(pow(-hrsp[nj][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))); 
      
      //esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
      
      // mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     : fgaus->GetParameter(1);
      // emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError(): hrsp[nj][ip]->GetMeanError(1);
      // sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetRMS()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
      // esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))); 
      
      //FitDist(hrsp[nj][ip], mean, emean, sig, esig);
      //FitGauss(hrsp[nj][ip], mean, emean, sig, esig);
      
      // hMean [nj]->SetBinContent (ip+1, roundoff(mean));
      // hMean [nj]->SetBinError   (ip+1, roundoff(emean));
      hMean [nj]->SetBinContent (ip+1, mean);
      hMean [nj]->SetBinError   (ip+1, emean);
      hSigma[nj]->SetBinContent (ip+1, sig);
      hSigma[nj]->SetBinError   (ip+1, esig);
    }
  }

  gPad->Close();
  
  // hMean[1][0]->Draw("p");
  // hMean[1][ncen-1]->Draw("psame");
  // return 0;

  int maxc=7;
  int maxr=2;
  int ipad=0;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLegend *leg0 = getLegend(0.1950635,0.4625641,0.8242194,0.7497436);
  //TLegend *leg0 = getLegend(0.4418803,0.3887016,0.8692251,0.6763714);  
  //if( wJetID )leg0->SetHeader("w JetID");
  //else leg0->SetHeader("w/o JetID");
  leg0->SetTextSize(0.10);
  leg0->AddEntry(hMean[0][ncen],Form("pp ak%s",(srad[0]+jetType).c_str()),"p");
  leg0->AddEntry(hMean[0][0],Form("PbPb ak%s",(algo+srad[0]+jetType).c_str()),"p");
  //leg0->AddEntry(hMean_resd[0][0],"Resd. Correc.","p");
  
  TLine *l0 = new TLine(xmin,1.00,xmax,1.0);
  l0->SetLineStyle(2);
  l0->SetLineWidth(2);
  TLine *l1 = new TLine(xmin,0.99,xmax,0.99);
  l1->SetLineStyle(2);
  TLine *l2 = new TLine(xmin,1.01,xmax,1.01);
  l2->SetLineStyle(2);

  TCanvas *c11[knj];
  for(int nj=0; nj<knj; nj++){
    c11[nj] = new TCanvas(Form("c11_%d",nj),Form("%s",(srad[nj]+jetType).c_str()),1695,530);
    makeMultiPanelCanvas(c11[nj],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);

    ipad=0;
    for(int ic=ncen; ic>=0; ic--){
      c11[nj]->cd(++ipad);
      gPad->SetLogx();

      hSigma[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hSigma[nj][ic]->Draw("p");

      hSigma[nj][ncen]->SetMarkerColor(2);
      hSigma[nj][ncen]->SetLineColor(2);
      hSigma[nj][ncen]->Draw("psame");

      if(ipad==2 && ic==ncen-1){
	drawText2("PYTHIA+HYDJET",0.19,0.85,21);
	leg0->Draw();
      }
      if(ipad==3 && ic==4){
	drawText2("#sqrt{s_{NN}} = 5.02 TeV,|#eta| < 2",0.10,0.85,19);
      }
      // if(ipad==1 && ic==5){
      // 	if(wJetID)drawText2(Form("%s, w JetID",(srad[nj]+jetType).c_str()),0.28,0.60,21);
      // 	else drawText2(Form("%s, w/o JetID",(srad[nj]+jetType).c_str()),0.28,0.60,21);
      // }
      //if(ipad==3)

      drawText2(ccent[ic],0.28,0.70,21);
      c11[nj]->cd(ipad+maxc);
      gPad->SetLogx();
    
      // if(ipad!=1)
      hMean[nj][ic]->GetXaxis()->SetRangeUser(xmin,xmax);
      hMean[nj][ic]->GetYaxis()->SetRangeUser(0.8,1.05);
      hMean[nj][ic]->Draw("p");
      
      hMean[nj][ncen]->SetMarkerColor(2);
      hMean[nj][ncen]->SetLineColor(2);
      hMean[nj][ncen]->Draw("psame");

      l1->Draw();
      l0->Draw();
      l2->Draw();
    }
    
    if(iSave) c11[nj]->SaveAs(Form("Plots/JEC_%d_PbPb_%s.pdf",wJetID,(algo+srad[nj]+jetType).c_str()));
  }
  //return 0;



  TCanvas *c99[knj][ncen];   ipad=0;
  for(int nj=0;nj<knj;nj++){
    for(int ic=ncen-1;ic>=0;ic--){
      c99[nj][ic] = new TCanvas(Form("c99%d_%d",nj,ic),Form("%s %s ",srad[nj].c_str(),ccent[ic]),100,102,1399,942);
      c99[nj][ic]->Divide(4,4,0,0);

      ipad=0;
      for(int ip=0;ip<nbins_pt;ip++){      
	c99[nj][ic]->cd(++ipad);
	
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	
	// hrsp[nj][ic][ip]->Scale(1./hrsp[nj][ic][ip]->Integral());
	// if(ic==ncen-3)hrsp[nj][ncen-1][ip]->Scale(1./hrsp[nj][ncen-1][ip]->Integral());	

	double ymax = TMath::Max(hrsp[nj][ic][ip]->GetMaximum(),hrsp[nj][ncen-1][ip]->GetMaximum());
	//double ymax = hrsp[nj][ic][ip]->GetMaximum();
	hrsp[nj][ic][ip]->SetMaximum(1.45*ymax);
	hrsp[nj][ic][ip]->SetMinimum(0);
	hrsp[nj][ic][ip]->SetTitle(0);
	hrsp[nj][ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
	hrsp[nj][ic][ip]->GetXaxis()->SetTitleFont(42);
	hrsp[nj][ic][ip]->GetXaxis()->SetLabelFont(42);
	hrsp[nj][ic][ip]->GetXaxis()->SetLabelSize(0.08);
	hrsp[nj][ic][ip]->GetXaxis()->SetTitleSize(0.07);
	hrsp[nj][ic][ip]->GetYaxis()->SetTitle("");
	hrsp[nj][ic][ip]->GetYaxis()->SetTitleFont(42);
	hrsp[nj][ic][ip]->GetYaxis()->SetLabelFont(42);
	hrsp[nj][ic][ip]->GetYaxis()->SetLabelSize(0.08);
      
	hrsp[nj][ic][ip]->SetMarkerStyle(20);
	hrsp[nj][ic][ip]->SetMarkerColor(1);
	hrsp[nj][ic][ip]->SetLineColor(1);
	hrsp[nj][ic][ip]->SetMarkerSize(1.3);
	hrsp[nj][ic][ip]->Draw("p");  

	hrsp[nj][ncen-1][ip]->SetMarkerStyle(20);
	hrsp[nj][ncen-1][ip]->SetMarkerColor(2);
	hrsp[nj][ncen-1][ip]->SetLineColor(2);
	hrsp[nj][ncen-1][ip]->SetMarkerSize(1.3);
	// hrsp[nj][ncen-1][ip]->Draw("psame");  
      
	std::ostringstream strs; 
	strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
	std::string spt = strs.str();

	if(ipad==1){
	  //drawText2(Form("%s",(srad[nj]+jetType).c_str()),0.28,0.90,21);      
	  drawText2(ccent[ic],0.75,0.87,21);	  
	  drawText2(spt.c_str(),0.22,0.80,21);		
	  if(wJetID)drawText2("w JetID",0.22,0.65,21);
	  else drawText2("w/o JetID",0.22,0.66,21);	  
	} 
	else drawText2(spt.c_str(),0.17,0.80,21);
      }
      if(iSave)	c99[nj][ic]->SaveAs(Form("Plots/IndDist_PbPb_%s_%s.pdf",ccent[ic],(algo+srad[nj]+jetType).c_str()));
      
    }
  }
  //}

  TFile *fout = new TFile(Form("JERJES_R%s_PbPb.root",srad[0].c_str()),"RECREATE");
  fout->cd( );
  fout->mkdir( Form("ak%sJetAnalyzer",srad[0].c_str()) );
  fout->cd( Form("ak%sJetAnalyzer",srad[0].c_str()) );
  for(int nj=0; nj<knj; nj++){
    for(int ic=0; ic<ncen; ic++){
      hMean[nj][ic]->SetName(Form("hMean_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
      hMean[nj][ic]->SetTitle(Form("hMean_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
      hMean[nj][ic]->Write();
      hSigma[nj][ic]->SetName(Form("hSigma_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
      hSigma[nj][ic]->SetTitle(Form("hSigma_PbPb_R%s_PuPF_%s",srad[nj].c_str(),cdir[ic]));
      hSigma[nj][ic]->Write();
    }
  }

  fout->Close();
  return 0;
}


int main(int argc, char*argv[])
{
  int rStatus=-1;
  if(argc!=1){
    std::cout<<"no arguments, just do "<<std::endl<<"./readFiles_ppMC_plotJESandJER.exe"<<std::endl<<std::endl;
    return rStatus;
  }

  rStatus=1;
  rStatus=plot_JESandJERClosure();
  return rStatus;
}
