#include "PEhit.cxx"
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"

double PulseShape(double x, double Amp, double tmax, double s_rise, double s_fall){
  double Amplitude = Amp;
  double mean = tmax;
  double sigma_rise = s_rise;
  double sigma_fall = s_fall;
  if(x < mean){return Amplitude*TMath::Gaus(x,mean,sigma_rise);}
  else return Amplitude*TMath::Gaus(x,mean,sigma_fall);
}

double WFshape(double *x, double *par){
  double Npulses = par[0];
  double baseline = par[1];
  double value = baseline;
  for(int i = 0; i< Npulses;i++){
    value+=PulseShape(x[0],par[4*i+2], par[4*i+3], par[4*i+4], par[4*i+5]);
  }
  if(value >= 4095.)value = 4095.;
  return value;
}

void SinglePE(){
  //TFile *fIN = new TFile("RESULTS_16_12_13.root","READ");
  TFile *fIN = new TFile("RESULTS.root","READ");
  TTree *T = (TTree*)fIN->Get("T");
  
  
  int     run, subrun, event;
  double  BNB_recoTime, MC_T0, recoT0;
  bool    isMC, is20PE, MC_isCCQE;
  double  threshold = 75.79;
  double  Nsigmas = 1.61;
  double  waveforms[32][1500];
  double  Baselines[32];
  std::vector< std::vector<PEhit> > *PElist = 0;
  
  T->SetBranchAddress("run", &run);
  T->SetBranchAddress("subrun", &subrun);
  T->SetBranchAddress("event",&event);
  T->SetBranchAddress("BNB_recoTime",&BNB_recoTime);
  T->SetBranchAddress("isMC",&isMC);
  T->SetBranchAddress("MC_isCCQE", &MC_isCCQE);
  T->SetBranchAddress("is20PE", &is20PE);
  T->SetBranchAddress("MC_T0",&MC_T0);
  T->SetBranchAddress("recoT0",&recoT0);
  T->SetBranchAddress("PElist",&PElist);
  T->SetBranchAddress("waveforms",&waveforms);
  T->SetBranchAddress("Baselines",&Baselines);
  
  TH1D *hSigmaRise = new TH1D("hSigmaRise","hSigmaRise",300,0,3);
  TH1D *hSigmaFall = new TH1D("hSigmaFall","hSigmaFall",300,2,5);
  TH1D *hAmplitude = new TH1D("hAmplitude","hAmplitude",3000,0,300);
  TH1D *hDTmin = new TH1D("hDTmin","hDTmin;[ns]",100,0,1500);
  TH2D *hSigmaRiseAmp = new TH2D("hSigmaRiseAmp","hSigmaRiseAmp",300,0,3,3000,0,300);
  TH2D *hSigmaFallAmp = new TH2D("hSigmaFallAmp","hSigmaFallAmp",300,2,5,3000,0,300);
  TH2D *hSigmas = new TH2D("hSigmas","hSigmas",300,0,3,300,2,5);
  
  double T0 = 0;
  double DTmin = 1000;
  double DTmax = 10;
  PEhit CurrentPE(20,1500,1,1);
  TCanvas *cWF = new TCanvas();
  TGraph *g = new TGraph();
  TF1 *f = new TF1("f",WFshape,0,1500,6);
  f->SetParameter(0,1);
  f->FixParameter(0,1);
  for(int entry = 0; entry < /*T->GetEntries()*/ 100;entry++){
    T->GetEntry(entry);    
    for(int ch = 0;ch < 32;ch++){
      for(size_t pe = 0;pe<PElist->at(ch).size();pe++){
	if(PElist->at(ch).at(pe).GetAmplitude()<10){continue;}
	if(PElist->at(ch).at(pe).GetAmplitude()>30){continue;}
	//if(PElist->at(ch).at(pe).GetTime() < 210 || PElist->at(ch).at(pe).GetTime() > 325){continue;}
	//if(PElist->at(ch).at(pe).GetSigmaRise() == 1.4 && PElist->at(ch).at(pe).GetSigmaFall() == 3.5){continue;}
	DTmin = 1000;
	for(size_t pe_comp = 0;pe_comp<PElist->at(ch).size();pe_comp++){
	  if(pe_comp == pe){continue;}
	  if(TMath::Abs(PElist->at(ch).at(pe).GetTime()-PElist->at(ch).at(pe_comp).GetTime()) < DTmin){DTmin = TMath::Abs(PElist->at(ch).at(pe).GetTime()-PElist->at(ch).at(pe_comp).GetTime());}
	}
	if(DTmin < DTmax){continue;}
	hSigmaRise->Fill(PElist->at(ch).at(pe).GetSigmaRise());
	hAmplitude->Fill(PElist->at(ch).at(pe).GetAmplitude());
	hDTmin->Fill(DTmin*15.625);
	hSigmaRiseAmp->Fill(PElist->at(ch).at(pe).GetSigmaRise(),PElist->at(ch).at(pe).GetAmplitude());
	hSigmaFall->Fill(PElist->at(ch).at(pe).GetSigmaFall());
	hSigmaFallAmp->Fill(PElist->at(ch).at(pe).GetSigmaFall(),PElist->at(ch).at(pe).GetAmplitude());
	hSigmas->Fill(PElist->at(ch).at(pe).GetSigmaRise(),PElist->at(ch).at(pe).GetSigmaFall());
	/*if(PElist->at(ch).at(pe).GetSigmaRise() == 1.4 && PElist->at(ch).at(pe).GetTime() > 200 && PElist->at(ch).at(pe).GetTime() < 400){
	  g->Clear();
	  for(int i = 0;i< 40;i++){//PElist->at(ch).at(pe).GetTime()-20;i<PElist->at(ch).at(pe).GetTime()+20;i++){
	    g->SetPoint(i,i+PElist->at(ch).at(pe).GetTime()-20,waveforms[ch][(int)(i+PElist->at(ch).at(pe).GetTime()-20)]);
	  }
	  cWF->cd();
	  g->SetMarkerStyle(20);
	  g->Draw("AP");
	  f->SetParameter(1,Baselines[ch]);
	  f->SetParameter(2,PElist->at(ch).at(pe).GetAmplitude());
	  f->SetParameter(3,PElist->at(ch).at(pe).GetTime());
	  f->SetParameter(4,PElist->at(ch).at(pe).GetSigmaRise());
	  f->SetParameter(5,PElist->at(ch).at(pe).GetSigmaFall());
	  f->Draw("same");
	  cWF->Modified();
	  cWF->Update();
	  cWF->SaveAs(Form("cWF_%d_%d.png",entry,ch));
	  }*/
      }
    }
  }

  TCanvas *cRise = new TCanvas();
  hSigmaRise->Draw();

  TCanvas *cFall = new TCanvas();
  hSigmaFall->Draw();

  TCanvas *cAmplitude = new TCanvas();
  hAmplitude->Draw();

  TCanvas *cDTmin = new TCanvas();
  cDTmin->SetLogy();
  hDTmin->Draw();
  
  TCanvas *cSigmaRiseAmp = new TCanvas();
  hSigmaRiseAmp->Draw("colz");

  TCanvas *cSigmaFallAmp = new TCanvas();
  hSigmaFallAmp->Draw("colz");

  TCanvas *cSigmas = new TCanvas();
  hSigmas->Draw("colz");
}
