#include "PEhit.cxx"
#include <vector>
#include "TCanvas.h"

void FirstPE(){
  TFile *fIN = new TFile("RESULTS_1.root","READ");
  if(!fIN){std::cout << "ERROR, could not open file" << std::endl; return;}
  else{std::string name = fIN->GetName();std::cout << "Opened " << name << " OK" << std::endl;}
  TTree *T = (TTree*)fIN->Get("T");
  
  
  int     run, subrun, event;
  double  BNB_recoTime, MC_T0, recoT0;
  bool    isMC, is20PE, MC_isCCQE;
  double  threshold = 75.79;
  double  Nsigmas = 1.61;
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
  
  TH1D *hSigmaRise = new TH1D("hSigmaRise","hSigmaRise",100,0,3);
  TH1D *hAmplitude = new TH1D("hAmplitude","hAmplitude",10000,0,1000);
  
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

    hSigmaRise->Fill(firstsmallPE.GetSigmaRise());
    hAmplitude->Fill(firstsmallPE.GetAmplitude());
  }
  TCanvas *c1 = new TCanvas();
  hSigmaRise->Draw();

  TCanvas *c2 =new TCanvas();
  hAmplitude->Draw();
  
}
