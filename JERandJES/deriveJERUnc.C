#include "deriveJERUnc.h"

const bool debugMode=true;
const bool drawPDFs=true;
const int radius = 4;

//pp_JER_plot
int JERandJES_deriveJERUnc(){

  // input file
  std::cout<< "input file dir is " << inFile_MC_dir<< std::endl;
  std::cout<< "input file name is " << inFile_MC_name<< std::endl;
  TFile *fin = TFile::Open(inFile_MC.c_str());

  // from fin
  TH1F *hJER[nbins_ana];
  TF1 *F1[nbins_ana];

  // for output/later use
  TCanvas *cJER[nbins_ana];

  float sigma[nbins_ana];  
  TH1F * hPP_JER_Uncert = new TH1F("hPP_JER_Uncert","",nbins_ana, ptbins_ana);
  TF1 * fResidual = new TF1("fResidual","1 - [0]/pow(x,[1])");

  if(debugMode)std::cout<<"nbins_ana= "<<nbins_ana<<std::endl;
  if(debugMode)std::cout<<"nbins= "<<nbins<<std::endl;

  for(int bin = 0; bin < (nbins_ana); ++bin){   //for(int bin = 1; bin < nbins_ana; ++bin){
    cJER[bin]= new TCanvas( ("cJER_"+std::to_string(bin) ).c_str() ,"",1200,1000);
    cJER[bin]->SetLogy();

    // open input hJER histos
    if(debugMode)std::cout<<"for bin = "<< bin<< std::endl;
    if(debugMode)std::cout<<"ptbins["<<bin<<"]= "<<ptbins[bin];
    if(debugMode)std::cout<<" , and ptbins["<<bin+1<<"]="<<ptbins[bin+1]<<std::endl;

    std::string histTitle="hJER_"+std::to_string(ptbins[bin])+"_pt_"+std::to_string(ptbins[bin+1]);
    std::cout<<"opening hist: "<<histTitle<<std::endl;
    hJER[bin] = (TH1F*)fin->Get(histTitle.c_str());  

    hJER[bin]->SetTitle(Form("PY6 Tune Z2 QCD %d < Gen p_{T} < %d; |#eta|<2", ptbins[bin], ptbins[bin+1]));
    hJER[bin]->SetXTitle("p_{T}^{detector jet}/p_{T}^{particle jet}");
    hJER[bin]->Print("base");

    // now fit
    F1[bin] = new TF1(Form("F1_bin%d",bin),"gaus", 0.4, 1.6);
    hJER[bin]->Fit(F1[bin],"R");
    hJER[bin]->Draw();

    drawText(Form("Mean = %f", F1[bin]->GetParameter(1)),0.15,0.75,14);
    drawText(Form("Sigma = %f", F1[bin]->GetParameter(2)),0.15,0.7,14);

    sigma[bin] = F1[bin]->GetParameter(2);
  }

  //cJER[14]->cd();
  //drawText(Form("PYTHIA Z2, ak%dPF Jets",radius),0.3,0.5,18);

  // define the residual correction factors
  if(radius == 2){
    fResidual->SetParameter(0, -0.05756);
    fResidual->SetParameter(1,  0.42750);
  }
  if(radius == 3){
    fResidual->SetParameter(0, -0.66229);
    fResidual->SetParameter(1,  1.02119);
  }
  if(radius == 4){
    fResidual->SetParameter(0, -1.28178);
    fResidual->SetParameter(1,  1.17348);
  }
  

  // set JER Uncertainties
  for(int bin = 1; bin <= hPP_JER_Uncert->GetNbinsX(); ++bin){
    float pt = hPP_JER_Uncert->GetBinCenter(bin);
    float res_val = (float)1./fResidual->Eval(pt) + sigma[bin-1];
    hPP_JER_Uncert->SetBinContent(bin, res_val);
  }


  // new canvas for JER uncertainties
  TCanvas * cPP_JER = new TCanvas("cPP_JER","",800,600);
  cPP_JER->cd();
  hPP_JER_Uncert->SetLineColor(kBlue);
  hPP_JER_Uncert->SetAxisRange(0.5, 1.5, "Y");
  hPP_JER_Uncert->SetAxisRange(60, 300, "X");
  hPP_JER_Uncert->SetTitle(" ");
  hPP_JER_Uncert->SetXTitle("Jet p_{T} (GeV/c)");
  hPP_JER_Uncert->SetYTitle("JER Uncertainty (%)");
  hPP_JER_Uncert->Draw();
  drawText(Form("ak%dPF Jets, PYTHIA", radius), 0.2, 0.2, 16);
  
  {
    // DRAW PDFS
    std::string thePDFFileName=baseName+".pdf";
    std::string open_thePDFFileName=thePDFFileName+"[";
    std::string close_thePDFFileName=thePDFFileName+"]";
    if(drawPDFs){
      TCanvas *temp_canv = new TCanvas("temp", "temp", 800, 600);
      temp_canv->Print( open_thePDFFileName.c_str() ); 
      
      temp_canv->cd();
      hPP_JER_Uncert->Draw();
      temp_canv->Print(thePDFFileName.c_str());
      
      for(int i=0;i<nbins_ana;++i){
      cJER[i]->cd();      
      cJER[i]->Print(thePDFFileName.c_str());
      }

      temp_canv->cd();
      temp_canv->Print( close_thePDFFileName.c_str() );      
    }
  }

  {
    // SAVE OUTPUT ROOT FILE
    std::string outFileName=baseName+".root";
    TFile *fout=new TFile(outFileName.c_str(),"RECREATE");   
    fout->cd();
    hPP_JER_Uncert->Write();  cPP_JER->Write();  
    for(int i=0;i<nbins_ana;++i){
      cJER[i]->Write();
    }
    fout->Close();  
  }

  fin->Close();
  return 0;

}

int main(int argc, char*argv[])
{
  int rStatus=-1;
  if(argc!=1){
    std::cout<<"no arguments, just do "<<std::endl<<"./JERandJES_deriveJERUnc.exe"<<std::endl<<std::endl;
    return rStatus;
  }

  rStatus=1;
  rStatus=JERandJES_deriveJERUnc();
  return rStatus;
}
