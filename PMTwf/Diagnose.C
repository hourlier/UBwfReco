#include "PEhit.cxx"
#include <vector>
#include "TCanvas.h"

TH1D *hDT = new TH1D("hDT","hDT;T_{0}(reco - th) [ns];N(/0.2 ns)",500,-50,50);

double GetResolution(double Nsigmas,double threshold, TFile *fIN){
  TTree *T = (TTree*)fIN->Get("T");
  
  
  int     run, subrun, event;
  double  BNB_recoTime, MC_T0, recoT0;
  bool    isMC, is20PE, MC_isCCQE;
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
  
  hDT->Reset();
  //hDT->Clear();
  
  double T0 = 0;
  PEhit firstPE(20,1500,1,1);
  PEhit firstsmallPE(20,1500,1,1);
  for(int entry = 0; entry < T->GetEntries();entry++){
    T->GetEntry(entry);
    
    int NPEtot = 0;
    double DTfirstPE = 1/15.625;// max DT between the first big peak and the first small PE considered
    for(UInt_t ch = 0; ch < PElist->size();ch++){
      NPEtot+=PElist->at(ch).size();
    }
    
    if(NPEtot==32){continue;}
    firstPE.SetValues(20,1500,1,1);
    for(UInt_t ch = 0; ch < PElist->size();ch++){
      for(UInt_t pe = 0;pe < PElist->at(ch).size();pe++){
	if(PElist->at(ch).at(pe).GetAmplitude() > threshold){
	  
	  if(PElist->at(ch).at(pe).GetTime()-Nsigmas*PElist->at(ch).at(pe).GetSigmaRise() < firstPE.GetTime()-Nsigmas*firstPE.GetSigmaRise()){
	    firstPE.SetValues(PElist->at(ch).at(pe).GetAmplitude(),PElist->at(ch).at(pe).GetTime(),PElist->at(ch).at(pe).GetSigmaRise(),PElist->at(ch).at(pe).GetSigmaFall());
	  }
	}
      }
    }
    
    firstsmallPE.SetValues(20,1500,1,1);
    for(UInt_t ch = 0; ch < PElist->size();ch++){
      for(UInt_t pe = 0;pe < PElist->at(ch).size();pe++){
	if(PElist->at(ch).at(pe).GetTime()-Nsigmas*PElist->at(ch).at(pe).GetSigmaRise() < firstsmallPE.GetTime()-Nsigmas*firstsmallPE.GetSigmaRise() && PElist->at(ch).at(pe).GetTime()-Nsigmas*PElist->at(ch).at(pe).GetSigmaRise() > firstPE.GetTime()-Nsigmas*firstPE.GetSigmaRise()+DTfirstPE){
	  firstsmallPE.SetValues(PElist->at(ch).at(pe).GetAmplitude(),PElist->at(ch).at(pe).GetTime(),PElist->at(ch).at(pe).GetSigmaRise(),PElist->at(ch).at(pe).GetSigmaFall());
	}
      }
    }

    hDT->Fill((firstsmallPE.GetTime()-Nsigmas*firstsmallPE.GetSigmaRise()-MC_T0)*15.625);
  }
  std::cout << "resolution ("<< Form("%.2f",Nsigmas) << "σ, " << Form("%.2f",threshold) << "ADC )  : " << hDT->GetRMS() << " ns" << std::endl;
  return hDT->GetRMS();
}


void Diagnose(int Npoints_ns = 1,int Npoints_th = 1){
  TFile *fIN = new TFile("RESULTS.root","READ");
  if(!fIN){std::cout << "ERROR, could not open file" << std::endl; return;}
  else{std::string name = fIN->GetName();std::cout << "Opened " << name << " OK" << std::endl;}
  double res = 0;
  double th, ns;
  double nsMin = 1;
  double nsMax = 2;
  double thMin = 60;
  double thMax = 80;
  double resMin = 25;
  double nsMinimum = 1;
  double thMinimum = 1;
  TH2D *hRes = new TH2D("hRes","hRes",3*Npoints_ns,-1,1+nsMax,3*Npoints_th,thMin-1,thMax+1);

  TGraph *gRes = new TGraph();
  for(int i = 0;i<Npoints_ns;i++){
    for(int j = 0;j<Npoints_th;j++){

      if(Npoints_th == 1){th = 63.84;}else{th = thMin+j*(thMax-thMin)/(Npoints_th-1);}
      if(Npoints_ns == 1){ns = 2;}else{ns =  nsMin+i*(nsMax-nsMin)/(Npoints_ns-1);}

      res = GetResolution(ns,th, fIN);
      if(res<resMin){resMin = res; nsMinimum = ns;thMinimum = th;}
      if(Npoints_th == 1 && Npoints_ns != 1){gRes->SetPoint(i,ns,res);}
      else if(Npoints_ns == 1 && Npoints_th != 1){gRes->SetPoint(j,th,res);}
      else if(Npoints_ns != 1 && Npoints_th != 1){hRes->Fill(ns,th,res);}
      else{std::cout << "done" << std::endl;}
    }
  }

  std::cout << std::endl;
  std::cout << "Minimum resolution for " << nsMinimum << " σ and min " << thMinimum << " ADC : " << resMin << std::endl;
  std::cout << std::endl;

  if((Npoints_th == 1 && Npoints_ns != 1) || (Npoints_th != 1 && Npoints_ns == 1)){
    gRes->SetMarkerStyle(20);
    gRes->Draw("AP");
  }
  if((Npoints_th != 1 && Npoints_ns != 1)){hRes->Draw("colz");}
  if(Npoints_th == 1 && Npoints_ns == 1){
    //TF1 *f = new TF1("f","gaus(0)",hDT->GetXaxis()->GetXmin(),hDT->GetXaxis()->GetXmax());
    //f->SetParameters(hDT->GetMaximum(),hDT->GetMean(),0.25*hDT->GetRMS());
    //hDT->Fit("f","l");
    hDT->Draw();
  }
}
