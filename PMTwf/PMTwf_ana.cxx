#ifndef LARLITE_PMTWF_ANA_CXX
#define LARLITE_PMTWF_ANA_CXX

#include "PMTwf_ana.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TF1.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TTree.h"
#include <vector>

#include "DataFormat/opdetwaveform.h"
#include "DataFormat/ophit.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcnu.h"
#include "DataFormat/simphotons.h"
#include "Base/MCConstants.h"
#include "LArUtil/Geometry.h"



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


namespace larlite {
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::initialize() {
    printOK = true;
    //printOK = false;
    drawOutput = false;
    //drawOutput = true;
    if(!drawOutput){printOK = false;
      //gROOT->SetBatch(kTRUE);
    }
    gROOT->SetBatch(kTRUE);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    SetDTHitReco(25);
    GetChannelGeo();
    for(int wf = 0;wf < 32;wf++){
      for(int tick = 0;tick< 1500;tick++){
	waveforms[wf][tick] = 0.;
      }
    }
    hEventWF = new TH2D("hEventWF","hEventWF;tick;Channel",1500,0,1500,32,0,32);
    hBNBwf = new TH1D("hBNBwf","hBNBwf;tick;ADC",1500,0,1500);
    hWFtmp = new TH1D("hWFtmp","hWFtmp",1500,0,1500);
    hDeltaT0 = new TH1D("hDeltaT0","hDeltaT0",300,-10,20);
    hDerivativeWFtmp = new TH1D("hDerivativeWFtmp","hDerivativeWFtmp",1500,0,1500);
    h2ndDerivativeWFtmp = new TH1D("h2ndDerivativeWFtmp","h2ndDerivativeWFtmp",1500,0,1500);
    hFlashFinding = new TH1D("hFlashFinding","hFlashFinding",150 ,0,1500);
    hFlashFindingDeriv = new TH1D("hFlashFindingDeriv","hFlashFindingDeriv",150 ,0,1500);
    hFlashT0 = new TH1D("hFlashT0","hFlashT0",130,210,340);
    hPreCut = new TH1D("hPreCut","hPreCut",250,0,1500);
    cEvent = new TCanvas("cEvent","cEvent",1000,600);
    cIndivWF = new TCanvas("cIndivWF","cIndivWF",1600,900);
    cEventT0 = new TCanvas("cEventT0","cEventT0",1000,700);
    cEventT0->Divide(1,2);
    cEvent->Divide(1,2);
    
    //if(printOK)cIndivWF->Print("waveforms.png[");
    flashStartTimeOff = 50;
    flashEndTimeOff = 150;
    beamWindowStart = 210;
    beamWindowEnd = 340;
    errorLevel = 0.3;
    FitRangeMin = 100;
    FitRangeMax = 600;
    newPEth = 0.1;
    //newPEth = 0.06;
    //RecoWFinfo = new std::vector<wfInfo>();

    kNPE = 0;
    kBaseline = 1;
    kAmplitude = 2;
    kTmax = 3;
    kSigmaRise = 4;
    kSigmaFall = 5;

    T = new TTree("T","T");
    T->Branch("run", &run);
    T->Branch("subrun", &subrun);
    T->Branch("event" , &event);
    T->Branch("BNB_recoTime", &BNB_recoTime);
    T->Branch("BNBparameters",BNBparameters,"BNBparameters[5]/D");
    T->Branch("isMC", &isMC);
    T->Branch("MC_T0",&MC_T0);
    T->Branch("MC_Neutrino_vertex",MC_Neutrino_vertex,"MC_Neutrino_vertex[4]/D");
    T->Branch("MCparticles","std::vector<int>", &MCparticles);
    T->Branch("recoT0",&recoT0);
    T->Branch("MC_isCCQE",&MC_isCCQE);
    T->Branch("is20PE",&is20PE);
    T->Branch("isFlashBeforBeam",&isFlashBeforBeam);
    T->Branch("T0wellReconstructed",&T0wellReconstructed);
    T->Branch("Baselines",Baselines,"Baselines[32]/D");
    T->Branch("BLrms",BLrms,"BLrms[32]/D");
    T->Branch("PElist", "std::vector< std::vector<PEhit> >",&PElist);

    return true;
  }
  
  //
  //****************************************************************
  //

  void PMTwf_ana::initializePElist(){
    PElist.clear();
    for(int ch = 0;ch<32;ch++){
      PEhit emptyPE;
      std::vector<PEhit> pelistwf = {emptyPE};
      PElist.push_back(pelistwf);
    }
  }

  //
  //****************************************************************
  //

  void PMTwf_ana::CleanPElist(){
    for(UInt_t ch = 0;ch<PElist.size();ch++){
      for(UInt_t pe = 0;pe<PElist.at(ch).size();pe++){
	if(PElist.at(ch).at(pe).GetTime() == -1){PElist.at(ch).erase(PElist.at(ch).begin()+pe);}
      }
    }
  }
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::analyze(storage_manager* storage) {
    //RecoWFinfo.clear();
    initializePElist();
    if(!LoadEventInfo(storage)){std::cout << "could not fully load event" << std::endl;return true;}
    GetBaselines();
    FindAllPulses();
    recoFlashTime = FindFlash();
    GetReferenceTime();

    FitRangeMin = recoFlashTime-25;
    FitRangeMax = recoFlashTime+100;
    if(!PMTpreCut()){std::cout << "rejected by PMprecut" << std::endl; return true;}
    if(isMC && !MC_isCCQE){std::cout << "MC but not CCQE" << std::endl; return true;}
    if(recoFlashTime == -1){std::cout << "no flash found" << std::endl; return true;}
    for(int ch = 0;ch < 32; ch++){
      std::cout << "\t" << ch << std::endl;
      if(!IsAmplitudeHighEnough(ch))continue;
      IterateFitProcess(ch);
      if(drawOutput)DrawWF(ch);
      //RecoWFinfo[ch].SetWF(Baselines[ch],BLrms[ch],PElist[ch]);
    }
    CleanPElist();
    T0wellReconstructed = FindFlashT0();
    //if(isMC){std::cout << MC_T0 << "\t" << recoT0 << std::endl;}
    //else{std::cout << BNB_recoTime << "\t" << recoT0 << std::endl;}
    T->Fill();
    return true;
  }

  //
  //****************************************************************
  //
  
  bool PMTwf_ana::finalize() {
    if(_fout){
      hDeltaT0->Write();
      T->Write();
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::PMTpreCut(){
    isFlashBeforBeam = false;
    is20PE = true;

    if(IsSummedWFHighbeforeBeam()){std::cout << "REJECTED : flash before beam window" << std::endl; isFlashBeforBeam = true; return false;}
    return true;
    hPreCut->Reset();
    for(int ch = 0; ch<32; ch++){
      for(UInt_t pulse = 0; pulse < RecoPulseTimesEvent[ch].size();pulse++){
	if(RecoPulseTimesEvent[ch][pulse] > beamWindowStart-10 && RecoPulseTimesEvent[ch][pulse] < beamWindowEnd+20)hPreCut->Fill(RecoPulseTimesEvent[ch][pulse]);
      }
    }
    int bin = 0;
    while(hPreCut->GetBinContent(bin+1) < 20 && bin < hPreCut->GetNbinsX()){bin++;}
    if(hPreCut->GetBinCenter(bin+1) > beamWindowEnd){std::cout << "REJECTED : 20PE limit in beam window not achieved" << std::endl; is20PE = false; return false;}

    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::IsSummedWFHighbeforeBeam(){

    double wfmin(10000), wfmax(0);
    double threshold = 10;
    for(int tick = 0; tick < beamWindowStart-10; tick++){
      if(SummedWF[tick] > wfmax){wfmax = SummedWF[tick];}
      if(SummedWF[tick] < wfmin){wfmin = SummedWF[tick];}
    }
    if(wfmax - wfmin > threshold)return true;
    else return false;
  }

  //
  //****************************************************************
  //
  bool PMTwf_ana::GetChannelGeo(){
    for(int i = 0;i<32;i++){
      larutil::Geometry::GetME()->GetOpChannelPosition(i,xyzCh[i]);
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::IsFlashInBeam(Int_t flash){
    if(FlashesTimes[flash]-flashStartTimeOff > beamWindowEnd || FlashesTimes[flash]+flashEndTimeOff < beamWindowStart)return false;
    else return true;
  }
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::IsFlashBeforeBeam(Int_t flash){
    if(FlashesTimes[flash] < beamWindowStart)return true;
    else return false;
  }

  //
  //****************************************************************
  //
  
  void  PMTwf_ana::DrawWF(Int_t ch){
    
    cIndivWF->Clear();
    cIndivWF->cd();
    pEvent2D = new TPad(Form("pEvent2D_%d_%d_%d_%d",run,subrun,event,ch),Form("pEvent2D_%d_%d_%d_%d",run,subrun,event,ch),0.01,0.255,0.475,0.99);
    pEvent2D->Draw();
    pEvent2D->SetTicks();
    pEvent2D->SetGrid();
    
    cIndivWF->cd();
    pFullWF = new TPad(Form("pFullWF_%d_%d_%d_%d",run,subrun,event,ch),Form("pFullWF_%d_%d_%d_%d",run,subrun,event,ch),0.01,0.01,0.475,0.245);
    pFullWF->Draw();
    pFullWF->SetTicks();
    pFullWF->SetGrid();
    
    cIndivWF->cd();
    pWaveform = new TPad(Form("pWaveform_%d_%d_%d_%d",run,subrun,event,ch),Form("pWaveform_%d_%d_%d_%d",run,subrun,event,ch),0.525,0.425,0.99,0.99);
    pWaveform->Draw();
    pWaveform->SetTicks();
    pWaveform->SetGrid();
    
    cIndivWF->cd();
    pDerivative = new TPad(Form("pDerivative_%d_%d_%d_%d",run,subrun,event,ch),Form("pDerivative_%d_%d_%d_%d",run,subrun,event,ch),0.525,0.01,0.99,0.375);
    pDerivative->Draw();
    pDerivative->SetTicks();
    pDerivative->SetGrid();
    
    pEvent2D->cd();
    hEventWF->Draw("colz");
    TBox *bBeam = new TBox(beamWindowStart,0,beamWindowEnd,32);
    bBeam->SetLineColor(6);
    bBeam->SetLineWidth(1);
    bBeam->SetFillStyle(3004);
    bBeam->SetFillColor(6);
    bBeam->Draw();
    TLine *lFlashEvent = new TLine(recoFlashTime,0,recoFlashTime,32);
    lFlashEvent->SetLineWidth(2);
    lFlashEvent->SetLineColor(2);
    lFlashEvent->Draw();

    TLine *lCh = new TLine(0,ch+0.5,1500,ch+0.5);
    lCh->SetLineColor(1);
    lCh->SetLineWidth(1);
    lCh->Draw();
    
    pDerivative->cd();
    ComputeResidual(ch);
    TH1D *hResidual = new TH1D(Form("hResidual_%d_%d_%d_%02d",run,subrun,event,ch),Form("hResidual_%d_%d_%d_%02d;tick;(WF-fit)/fit",run,subrun,event,ch),1500,0,1500);
    for(int i = 0;i<1500;i++){hResidual->SetBinContent(i+1,WFresidual[i]);}
    hResidual->GetXaxis()->SetRange(FitRangeMin,FitRangeMax);
    hResidual->Draw();
    
    pWaveform->cd();
    hWFtmp = GetChWF(ch);
    hWFtmp->GetXaxis()->SetRange(beamWindowStart-10,beamWindowEnd+150);
    hWFtmp->Draw();
    
    //FitFcn->SetNpx(1500);
    FitFcn->Draw("same");
    
    TLine *lBaseline = new TLine(beamWindowStart-10,Baselines[ch],beamWindowEnd+150,Baselines[ch]);
    lBaseline->SetLineColor(2);
    lBaseline->SetLineWidth(2);
    lBaseline->SetLineStyle(2);
    pWaveform->cd();
    lBaseline->Draw();
    
    pFullWF->cd();
    TH1D *hwfFull = (TH1D*)hWFtmp->Clone("hwfFull");
    hwfFull->GetXaxis()->SetRange(1,1500);
    hwfFull->Draw();
    FitFcn->Draw("same");
    //FilterWaveform(ch);
    //TH1D *hfiltered = GetfilteredWF(ch);
    //hfiltered->SetLineColor(2);
    //hfiltered->Draw("same");
    
    
    cIndivWF->Modified();
    cIndivWF->Update();
    cIndivWF->SetName(Form("cIndivWF_%d_%d_%d_%02d",run,subrun,event,ch));
    cIndivWF->SetTitle(Form("cIndivWF_%d_%d_%d_%02d",run,subrun,event,ch));


    cEvent->cd();
    hwfFull->GetXaxis()->SetRange(FitRangeMin,FitRangeMax);
    hwfFull->SetMarkerStyle(20);
    hwfFull->Draw();
    FitFcn->Draw("same");
    TLine *lFlash = new TLine(recoFlashTime,hwfFull->GetBinContent(hwfFull->GetMinimumBin()),recoFlashTime,hwfFull->GetBinContent(hwfFull->GetMaximumBin()));
    lFlash->SetLineColor(2);
    lFlash->SetLineWidth(2);
    lFlash->Draw();

    TLine *lPEtruth;
    for(UInt_t pe = 0;pe < MChitTimes[ch].size();pe++){
      lPEtruth = new TLine(MChitTimes[ch][pe],hwfFull->GetBinContent(hwfFull->GetMinimumBin()),MChitTimes[ch][pe],hwfFull->GetBinContent(hwfFull->GetMaximumBin()));
      lPEtruth->SetLineColor(4);
      lPEtruth->Draw();
    }

    cEvent->SetName(Form("cEvent_%d_%d_%d_%02d",run,subrun,event,ch));
    cEvent->SetTitle(Form("cEvent_%d_%d_%d_%02d",run,subrun,event,ch));
    cEvent->Modified();
    cEvent->Update();
    //if(printOK)cIndivWF->Print("waveforms.png");
    if(printOK){
      cIndivWF->SaveAs(Form("plot/waveform_%d_%d_%d_%02d.png",run,subrun,event,ch));
      cEvent->SaveAs(Form("plot/FullWF_%d_%d_%d_%02d.png",run,subrun,event,ch));
    }
    //if(_fout){
    //  _fout->cd();
    //  cIndivWF->Write();
    //}
  }

  
  //
  //****************************************************************
  //

  bool PMTwf_ana::DrawFlashWF(Int_t flash, Int_t ch){

    cIndivWF->Clear();
    cIndivWF->cd();
    pEvent2D = new TPad(Form("pEvent2D_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),Form("pEvent2D_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),0.01,0.255,0.475,0.99);
    pEvent2D->Draw();
    pEvent2D->SetTicks();
    pEvent2D->SetGrid();

    cIndivWF->cd();
    pFullWF = new TPad(Form("pFullWF_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),Form("pFullWF_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),0.01,0.01,0.475,0.245);
    pFullWF->Draw();
    pFullWF->SetTicks();
    pFullWF->SetGrid();

    cIndivWF->cd();
    pWaveform = new TPad(Form("pWaveform_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),Form("pWaveform_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),0.525,0.425,0.99,0.99);
    pWaveform->Draw();
    pWaveform->SetTicks();
    pWaveform->SetGrid();

    cIndivWF->cd();
    pDerivative = new TPad(Form("pDerivative_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),Form("pDerivative_%d_%d_%d_%d_%d",run,subrun,event,flash,ch),0.525,0.01,0.99,0.375);
    pDerivative->Draw();
    pDerivative->SetTicks();
    pDerivative->SetGrid();

    pEvent2D->cd();
    hEventWF->Draw("colz");
    TBox *bBeam = new TBox(beamWindowStart,0,beamWindowEnd,32);
    bBeam->SetLineColor(6);
    bBeam->SetLineWidth(1);
    bBeam->SetFillStyle(3004);
    bBeam->SetFillColor(6);
    bBeam->Draw();
    TBox *bFlash = new TBox(TMath::Max(0.,FlashesTimes[flash]-flashStartTimeOff),0,TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.),32);
    bFlash->SetLineColor(2);
    bFlash->SetLineWidth(3);
    bFlash->SetFillStyle(0);
    bFlash->Draw();
    TLine *lFlash = new TLine(FlashesTimes[flash],0,FlashesTimes[flash],32);
    lFlash->SetLineColor(2);
    lFlash->SetLineWidth(1);
    lFlash->Draw();
    TLine *lCh = new TLine(0,ch+0.5,1500,ch+0.5);
    lCh->SetLineColor(1);
    lCh->SetLineWidth(1);
    lCh->Draw();

    pDerivative->cd();
    TH1D *hResidual = new TH1D(Form("hResidual_%d_%d_%d_%02d_%02d",run,subrun,event,flash,ch),Form("hResidual_%d_%d_%d_%02d_%02d;tick;(WF-fit)/fit",run,subrun,event,flash,ch),1500,0,1500);
    for(int i = 0;i<1500;i++){hResidual->SetBinContent(i+1,WFresidual[i]);}
    hResidual->GetXaxis()->SetRange(TMath::Max(1.,FlashesTimes[flash]-flashStartTimeOff),TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.));
    hResidual->Draw();

    pWaveform->cd();
    hWFtmp = GetChWF(ch);
    hWFtmp->GetXaxis()->SetRange(TMath::Max(1.,FlashesTimes[flash]-flashStartTimeOff),TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.));
    hWFtmp->Draw();

    //FitFcn->SetNpx(1500);
    //FitFcn->Draw("same");
    
    TLine *lBaseline = new TLine(TMath::Max(1.,FlashesTimes[flash]-flashStartTimeOff),Baselines[ch],TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.),Baselines[ch]);
    lBaseline->SetLineColor(2);
    lBaseline->SetLineWidth(2);
    lBaseline->SetLineStyle(2);
    pWaveform->cd();
    lBaseline->Draw();
    
    pFullWF->cd();
    TH1D *hwfFull = (TH1D*)hWFtmp->Clone("hwfFull");
    hwfFull->GetXaxis()->SetRange(1,1500);
    hwfFull->Draw();
    //FilterWaveform(ch);
    //TH1D *hfiltered = GetfilteredWF(ch);
    //hfiltered->SetLineColor(2);
    //hfiltered->Draw("same");
    

    cIndivWF->Modified();
    cIndivWF->Update();
    cIndivWF->SetName(Form("cIndivWF_%d_%d_%d_%02d_%02d",run,subrun,event,flash,ch));
    cIndivWF->SetTitle(Form("cIndivWF_%d_%d_%d_%02d_%02d",run,subrun,event,flash,ch));
    //if(printOK)cIndivWF->Print("waveforms.png");
    if(printOK)cIndivWF->SaveAs(Form("plot/waveform_%d_%d_%d_%02d_%02d.png",run,subrun,event,flash,ch));
    if(_fout){
      _fout->cd();
      cIndivWF->Write();
    }
    return true;
  }

  //
  //****************************************************************
  //

  void PMTwf_ana::FilterWaveform(Int_t ch){
    filteredWF[ch][0] = waveforms[ch][0];
    filteredWF[ch][1] = waveforms[ch][1];
    filteredWF[ch][2] = waveforms[ch][2];
    filteredWF[ch][3] = waveforms[ch][3];
    filteredWF[ch][1496] = waveforms[ch][1496];
    filteredWF[ch][1497] = waveforms[ch][1497];
    filteredWF[ch][1498] = waveforms[ch][1498];
    filteredWF[ch][1499] = waveforms[ch][1499];
    
    int LastGoodTick = 0;
    int NextGoodTick = 1;
    
    double LastGoodVal = 1;
    double NextGoodVal = 1;
    
    for(int tick = 3;tick < 1497;tick++){//initialize the baseline array with values of the WF for flat enough parts
      if(WFderivs[ch][tick-3] == 0 && WFderivs[ch][tick-2] == 0 && WFderivs[ch][tick-1] == 0 && WFderivs[ch][tick+1] == 0 && WFderivs[ch][tick+2] == 0 && WFderivs[ch][tick+3] == 0 && WFderivs[ch][tick] == 0 && waveforms[ch][tick] < 4095){
	filteredWF[ch][tick] = waveforms[ch][tick];
      }
      else filteredWF[ch][tick] = 0;
    }
    
    for(int i = 1;i<1499;i++){//re-run
      if(filteredWF[ch][i] != 0){//non-zero values => good flat parts
	LastGoodTick = i;
	LastGoodVal = filteredWF[ch][i];
      }
      else{// zero values, need to search for the end, evaluate values at the begining and end of the zero part and fill it with a straight line
	int j = i;
	while(filteredWF[ch][j] == 0 && j<1499){j++;}//loop over j until we are no longer in a zero region
	NextGoodTick = j;
	NextGoodVal = filteredWF[ch][j];

	filteredWF[ch][i] = LastGoodVal+(i-LastGoodTick)*(NextGoodVal-LastGoodVal)/(NextGoodTick-LastGoodTick);
      }
    }
    
  }
  
  //
  //****************************************************************
  //
  
  void PMTwf_ana::FilterWaveformOld(Int_t ch){
    double lastgoodsample = waveforms[ch][0];
    for(int tick = 0;tick < 1500;tick++){
      filteredWF[ch][tick] = waveforms[ch][tick];
    }
    bool filtered = false;
    int iter = 0;
    while(filtered == false && iter < 4){
      filtered = true;
      iter++;
      lastgoodsample = waveforms[ch][0];
      for(int tick = 0;tick < 1500;tick++){
	if(tick >=4 && tick <= 1496){
	  //if(filteredWF[ch][tick] > filteredWF[ch][tick-3]+2 || filteredWF[ch][tick] < filteredWF[ch][tick-3]-2 || WFderivs[ch][tick] > 1.5 || WFderivs[ch][tick] < -1.5){
	  //if(0.25*(filteredWF[ch][tick-2]+filteredWF[ch][tick-1]+filteredWF[ch][tick+1]+filteredWF[ch][tick+2])-filteredWF[ch][tick] > 1 || 0.25*(filteredWF[ch][tick-2]+filteredWF[ch][tick-1]+filteredWF[ch][tick+1]+filteredWF[ch][tick+2])-filteredWF[ch][tick] < -1){  
	  if(0.5*(filteredWF[ch][tick-3]+filteredWF[ch][tick+3])-filteredWF[ch][tick] > 1 || 0.5*(filteredWF[ch][tick-3]+filteredWF[ch][tick+3])-filteredWF[ch][tick] < -1){
	    filteredWF[ch][tick] = lastgoodsample;
	    filtered = false;
	  }
	  else{
	    lastgoodsample = filteredWF[ch][tick];
	  }
	}
	else{
	  lastgoodsample = filteredWF[ch][tick];
	  filteredWF[ch][tick] = waveforms[ch][tick];
	}
      }
    }
    
  }

  //
  //****************************************************************
  //

  TH1D* PMTwf_ana::GetfilteredWF(Int_t ch){
    TH1D *hfiltered = new TH1D(Form("hfiltered_run%d_%d_evt%d_ch%d",run,subrun,event,ch),Form("hfiltered_run%d_%d_evt%d_ch%d",run,subrun,event,ch),1500,0,1500);
    for(int tick = 0; tick < 1500; tick++){
      hfiltered->SetBinContent(tick+1,filteredWF[ch][tick]);
    }
    return hfiltered;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::DrawEvent(){
    return true;
  }
  
  //
  //****************************************************************
  //
  
  TH1D* PMTwf_ana::GetBNBwf(){
    TH1D *hBNBwf = new TH1D(Form("hBNBwf_run%d_%d_evt%d",run,subrun,event),Form("hBNBwf_run%d_%d_evt%d",run,subrun,event),1500,0,1500);
    for(int tick = 0; tick < 1500; tick++){
      hBNBwf->SetBinContent(tick+1,BNBwaveform[tick]);
      hBNBwf->SetBinError(tick+1,0.5);
    }
    return hBNBwf;
  }
  
  //
  //****************************************************************
  //

  TH1D* PMTwf_ana::GetChWF(Int_t ch){
    //TH1D *hChWF = new TH1D(Form("hChwf_run%d_%d_evt%d_ch%d",run,subrun,event,ch),Form("hChwf_run%d_%d_evt%d_ch%d",run,subrun,event,ch),1500,0,1500);
    for(int tick = 0; tick < 1500; tick++){
      hWFtmp->SetBinContent(tick+1,waveforms[ch][tick]);
      //hWFtmp->SetBonError(tick+1, errorLevel);
    }
    hWFtmp->SetTitle(Form("hChwf_run%d_%d_evt%d_ch%d",run,subrun,event,ch));
    return hWFtmp;
  }
  
  //
  //****************************************************************
  //

  TH1D* PMTwf_ana::GetDerivative1(Int_t ch){
    //TH1D *hDerivative1 = new TH1D(Form("h1stDerivative_run%d_%d_evt%d_ch%d",run,subrun,event,ch),Form("h1stDerivative_run%d_%d_evt%d_ch%d",run,subrun,event,ch),1500,0,1500);
    for(int tick = 0; tick < 1500; tick++){
      hDerivativeWFtmp->SetBinContent(tick+1,WFderivs[ch][tick]);
    }
    hDerivativeWFtmp->SetTitle(Form("h1stDerivative_run%d_%d_evt%d_ch%d",run,subrun,event,ch));
    return hDerivativeWFtmp;
  }

  //
  //****************************************************************
  //

  TH1D* PMTwf_ana::GetDerivative2(Int_t ch){
    //TH1D *hDerivative2 = new TH1D(Form("h2ndDerivative_run%d_%d_evt%d_ch%d",run,subrun,event,ch),Form("h2ndDerivative_run%d_%d_evt%d_ch%d",run,subrun,event,ch),1500,0,1500);
    for(int tick = 0; tick < 1500; tick++){
      h2ndDerivativeWFtmp->SetBinContent(tick+1,WF2ndDerivs[ch][tick]);
    }
    h2ndDerivativeWFtmp->SetTitle(Form("h2ndDerivative_run%d_%d_evt%d_ch%d",run,subrun,event,ch));
    return h2ndDerivativeWFtmp;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::GetAllDerivatives(){
    for(int ch = 0;ch < 32; ch++){
      GetDerivativeWF(ch);
      Get2ndDerivativeWF(ch);
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::GetDerivativeWF(Int_t ch){
    WFderivs[ch][0] = 0;
    for(int tick = 1; tick < 1500; tick++){
      WFderivs[ch][tick] = 0.5*(waveforms[ch][tick]-waveforms[ch][tick-1]);
      if(WFderivs[ch][tick] == -1.*WFderivs[ch][tick-1]){
	WFderivs[ch][tick] = 0;
	WFderivs[ch][tick-1] = 0;
      }
    }
    for(int tick = 2; tick < 1500; tick++){
      if(WFderivs[ch][tick-2] == 0 && WFderivs[ch][tick-1] != 0 && WFderivs[ch][tick] == 0){WFderivs[ch][tick-1] = 0;}
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::Get2ndDerivativeWF(Int_t ch){
    WF2ndDerivs[ch][0] = 0;
    for(int tick = 1; tick < 1500; tick++){
      WF2ndDerivs[ch][tick] = 0.5*(WFderivs[ch][tick]-WFderivs[ch][tick-1]);
      if(WF2ndDerivs[ch][tick] == -1.*WF2ndDerivs[ch][tick-1]){
        WF2ndDerivs[ch][tick] = 0;
        WF2ndDerivs[ch][tick-1] = 0;
      }
    }
    for(int tick = 2; tick < 1500; tick++){
      if(WF2ndDerivs[ch][tick-2] == 0 && WF2ndDerivs[ch][tick-1] != 0 && WF2ndDerivs[ch][tick] == 0){WF2ndDerivs[ch][tick-1] = 0;}
    }
    return true;
  }
  
  //
  //****************************************************************
  //
  
  void PMTwf_ana::DrawChPosition(){
    TCanvas *cPMTposition = new TCanvas("cPMTposition","cPMTposition",1100,220);
    TGraph *gPMTPosition = new TGraph();
    for(int i = 0;i<32;i++){
      gPMTPosition->SetPoint(i,xyzCh[i][2],xyzCh[i][1]);
    }
    gPMTPosition->SetMarkerStyle(20);
    gPMTPosition->SetMarkerSize(2);
    cPMTposition->cd();
    gPMTPosition->Draw("AP");
  }
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::FindPulses(Int_t ch){
   
    bool peakStarted = false;
    RecoPulseTimes.clear();
    for(int i = 0;i< 1500;i++){
      if(peakStarted == false && WFderivs[ch][i] > 0){// found new peak
	peakStarted = true;
	if((waveforms[ch][i]+waveforms[ch][i+1]+waveforms[ch][i+2]+waveforms[ch][i+3]+waveforms[ch][i+4]-5*Baselines[ch]) > 6*errorLevel){
	  RecoPulseTimes.push_back(i);
	}
      }
      if(peakStarted == true && WFderivs[ch][i] <= 0){// end peak
	peakStarted = false;
      }
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::FindAllPulses(){
    for(int ch = 0;ch<32;ch++){
      RecoPulseTimesEvent[ch].clear();
      FindPulses(ch);
      RecoPulseTimesEvent[ch] = RecoPulseTimes;
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::LoadEventInfo(storage_manager* storage){
    auto ophits        = storage->get_data<event_ophit>("ophit");
    auto opflashes     = storage->get_data<event_opflash>("opflash");
    auto mcdata        = storage->get_data<event_mctrack>("mcreco");
    auto mc            = storage->get_data<event_mctruth>("generator");
    auto opdetwaveform = storage->get_data<event_opdetwaveform>("pmtreadout");
    auto photonInfo    = storage->get_data<event_simphotons>("largeant");
    
    UInt_t correspondingPMTnumber=-1;
    if(photonInfo){
      for(UInt_t pm = 0; pm < photonInfo->size();pm++){
	for(UInt_t ch = 0; ch < 32;ch++){
	  correspondingPMTnumber = larutil::Geometry::GetME()->OpDetFromOpChannel(ch);
	  if(correspondingPMTnumber == pm){
	    MChitTimes[ch].clear(); // careful, this index number is NOT the number of ch but the number of PMT 
	    for(UInt_t pe = 0;pe < photonInfo->at(pm).size();pe++){
	      MChitTimes[ch].push_back((photonInfo->at(pm).at(pe).Time)*64/1000.);//geant4 time in ns, need to convert back in ticks
	    }
	  }
	}
      }
    }
    if(mc){
      MCparticles.clear();
      if(mc->at(0).GetNeutrino().InteractionType() == simb::kCCQE){MC_isCCQE = true;}
      else{MC_isCCQE = false;}
      for(size_t mctruth_index = 0; mctruth_index<mc->size();mctruth_index++){
	auto const& mctruth = mc->at(mctruth_index);
	for(size_t partcl = 0; partcl < (size_t)(mctruth.NParticles());partcl++){
	  auto const& mcp = mctruth.GetParticle(partcl);
	  if(mcp.StatusCode() == 0 && (mcp.PdgCode() == 12 || mcp.PdgCode() == -12 || mcp.PdgCode() == 14 || mcp.PdgCode() == -14 || mcp.PdgCode() == 16 || mcp.PdgCode() == -16)){
	    auto const& pos = mcp.Position(0);
	    MC_Neutrino_vertex[0] = pos.X();
	    MC_Neutrino_vertex[1] = pos.Y();
	    MC_Neutrino_vertex[2] = pos.Z();
	    MC_Neutrino_vertex[3] = pos.T();
	    MCparticles.push_back(mcp.PdgCode());
	  }
	}
      }
    }

    if(mc || photonInfo){isMC = true;}else{isMC = false;}
    if(mcdata){
      for(UInt_t track = 0; track < mcdata->size(); track++){
	MCparticles.push_back(mcdata->at(track).PdgCode());
      }
    }
    
    std::cout << "Run #" << storage->run_id() << "_" << storage->subrun_id() << "\t event #" << storage->event_id() << std::endl;
    run = storage->run_id();
    subrun = storage->subrun_id();
    event = storage->event_id();

    //GetWaveforms
    for(Int_t i = 0;i<1500;i++){SummedWF[i] = 0;}
    if(!opdetwaveform){std::cout << "\t \t ERROR >>>>> NO OPDETWAVEFORM" << std::endl;return false;}
    for(UInt_t wf = 0;wf < opdetwaveform->size();wf++){
      if(opdetwaveform->at(wf).ChannelNumber() == 38){
        for(UInt_t tick = 0;tick < 1500;tick++){BNBwaveform[tick] = opdetwaveform->at(wf)[tick];}
      }
      if(opdetwaveform->at(wf).ChannelNumber() < 32){
	for(UInt_t tick = 0; tick < 1500; tick++){
	  waveforms[opdetwaveform->at(wf).ChannelNumber()][tick] = opdetwaveform->at(wf)[tick];
	  hEventWF->SetBinContent(tick+1,opdetwaveform->at(wf).ChannelNumber()+1,waveforms[opdetwaveform->at(wf).ChannelNumber()][tick]-2000.);
	  SummedWF[tick]+=waveforms[opdetwaveform->at(wf).ChannelNumber()][tick]/32.;
	}
      }
    }
    
    // Get Derivatives
    GetAllDerivatives();

    // GetFlashInfo
    FlashesTimes.clear();
    PEwf1flash.clear();
    PEflashWF.clear();
    if(opflashes){
      for(UInt_t flash = 0;flash<opflashes->size();flash++){
	FlashesTimes.push_back(opflashes->at(flash).Time()*64.);
	flashStartTime = opflashes->at(flash).Time()*64.-flashStartTimeOff;
	flashEndTime = opflashes->at(flash).Time()*64.+flashEndTimeOff;
	PEwf1flash.clear();
	for(int wf = 0; wf < 32; wf++){
	  PEwf1flash.push_back(opflashes->at(flash).PE(wf));
	}
	PEflashWF.push_back(PEwf1flash);
      }
    }
    
    // GetHitInfo
    if(ophits){
      for(int wf = 0; wf < 32; wf++){
	HitTimes.clear();
	for(UInt_t hit = 0; hit < ophits->size();hit++){
	  if(ophits->at(hit).OpChannel() == wf){
	    HitTimes.push_back(ophits->at(hit).PeakTime()*64-DTHitReco);
	  }
	}
	AllHitTimes[wf] = HitTimes;
      }
    }
    return true;
  }
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::FitFlash(Int_t flash){
    for(int ch= 0;ch<32;ch++){
      if(PEflashWF[flash][ch] > 2){
	FitFlashWF(flash,ch);
      }
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::FitFlashWF(Int_t flash, Int_t ch){
    TH1D *hWFtmp = GetChWF(ch);
    for(int tick = 0;tick < 1500;tick++){hWFtmp->SetBinError(tick+1,errorLevel);}
    hWFtmp->GetXaxis()->SetRange(FlashesTimes[flash]-flashStartTimeOff,FlashesTimes[flash]+flashEndTimeOff);
    FitFcn = InitializeFitFcn(flash,ch);
    FitFcn->FixParameter(0,FitFcn->GetParameter(0));
    
    if(drawOutput == true){hWFtmp->Fit(Form("%s",FitFcn->GetName()),"qrl","");}
    else{hWFtmp->Fit(Form("%s",FitFcn->GetName()),"qrlO","");}
    //ComputeResidual(flash,ch);
    return true;
  }

  //
  //****************************************************************
  //

  void PMTwf_ana::FitWF(Int_t ch){
    TH1D *hWFtmp = GetChWF(ch);
    for(int tick = 0;tick < 1500;tick++){hWFtmp->SetBinError(tick+1,errorLevel);}
    FitFcn = InitializeFitFcn(ch);
    FitFcn->FixParameter(0,FitFcn->GetParameter(0));
    /*for(UInt_t i = 2;i<WFparameters.size();i++){
      if((i-4)/4==0 || (i-5)/4==0)FitFcn->FixParameter(i,FitFcn->GetParameter(i));
      }*/
    hWFtmp->Fit(Form("%s",FitFcn->GetName()),"rlq","",FitRangeMin,FitRangeMax);
    if(drawOutput == false){
      hWFtmp->Fit(Form("%s",FitFcn->GetName()),"rlOq","",FitRangeMin,FitRangeMax);
    }
    Baselines[ch] = FitFcn->GetParameter(1);
    SetPElist(ch);
    
  }

  //
  //****************************************************************
  //

  void PMTwf_ana::IterateFitProcess(Int_t ch){
    int IterFit = 0;
    //InitializeFitFcn(ch);
    FitWF(ch);
    while(NeedNewPE(ch) == true && IterFit<10){
      FitWF(ch);
      IterFit++;
    }
  }

  //
  //****************************************************************
  //

  void PMTwf_ana::ComputeResidual(Int_t flash, Int_t ch){
    for(int tick = 0; tick < 1500;tick++){
      if(tick > FlashesTimes[flash]-flashStartTimeOff && tick < FlashesTimes[flash]+flashEndTimeOff){
	WFresidual[tick] = (waveforms[ch][tick]-FitFcn->Eval(tick+0.5))/(FitFcn->Eval(tick+0.5)-(Baselines[ch]-50));
      }
      else WFresidual[tick] = 0;
    }
  }

  //
  //****************************************************************
  //

  void PMTwf_ana::ComputeResidual(Int_t ch){
    for(int tick = 0; tick < 1500;tick++){
      WFresidual[tick] = (waveforms[ch][tick]-FitFcn->Eval(tick+0.5))/(FitFcn->Eval(tick+0.5)-(Baselines[ch]-50));
    }
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::GetBaselines(){
    for(int ch = 0;ch < 32;ch++){
      EvalBaseline(ch);
    }
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::EvalBaseline(Int_t ch){
    int LastHighTick = 0;
    int DT =50;
    double threshold = 1.25;
    int TimeSinceHigh = DT-5;
    bool InIntegrationWindow = false;

    TH1D *hBaseline = new TH1D("hBaseline","hBaseline",1000,1800,2200);
    
    for(int tick = 0;tick<beamWindowStart;tick++){// eval baseline on the pre-beam window => this is where I want it stable (i.e. before the neutrino events)
      if(TMath::Abs(WFderivs[ch][tick]) < threshold){
	if(TimeSinceHigh >= DT && InIntegrationWindow == false){
	  InIntegrationWindow = true;
	  TimeSinceHigh = tick-LastHighTick;
	}
	TimeSinceHigh = tick-LastHighTick;
      }
      else{
	if(InIntegrationWindow == true){
	  InIntegrationWindow = false;
	  TimeSinceHigh = 0;
	  LastHighTick = tick;
	}
      }
      if(tick == 1499 && InIntegrationWindow == true){
	InIntegrationWindow = false;
      }
      if(InIntegrationWindow){
	hBaseline->Fill(waveforms[ch][tick]);
      }
    }
    
    
    Baselines[ch] = hBaseline->GetMean();
    BLrms[ch] = hBaseline->GetRMS();
    hBaseline->Delete();
    return true;
  }

  //
  //****************************************************************
  //

  TF1*  PMTwf_ana::InitializeFitFcn(Int_t flash, Int_t ch){
    NpeaksFound = 0;
    WFparameters.clear();
    WFparLimitHigh.clear();
    WFparLimitLow.clear();

    WFparameters.push_back(0);
    WFparLimitHigh.push_back(0);
    WFparLimitLow.push_back(0);
    WFparameters.push_back(Baselines[ch]);
    WFparLimitHigh.push_back(Baselines[ch]+5);
    WFparLimitLow.push_back(Baselines[ch]-5);

    int ThisPulseTime;
    int ThisPulseMaxTick = 0;;
    double ThisPulseMaxAmp = 0;
    for(UInt_t pulse = 0;pulse < RecoPulseTimesEvent[ch].size(); pulse++){
      if(RecoPulseTimesEvent[ch][pulse] > FlashesTimes[flash] - flashStartTimeOff && RecoPulseTimesEvent[ch][pulse] < FlashesTimes[flash] + flashEndTimeOff){
	ThisPulseTime = RecoPulseTimesEvent[ch][pulse];
	if(pulse < RecoPulseTimesEvent[ch].size()-1){
	  ThisPulseMaxTick = ThisPulseTime;
	  ThisPulseMaxAmp = 0;
	  for(int i = ThisPulseTime; i<RecoPulseTimesEvent[ch][pulse+1] ;i++){
	    if(waveforms[ch][i] >= ThisPulseMaxAmp){ThisPulseMaxAmp = waveforms[ch][i];ThisPulseMaxTick = i+0.5;}
	  }
	}
	else{ThisPulseMaxTick = ThisPulseTime+3.5;}
	
	WFparameters.push_back(waveforms[ch][ThisPulseMaxTick]-Baselines[ch]); // amplitude
	WFparLimitHigh.push_back(10000);
	WFparLimitLow.push_back(0);

	WFparameters.push_back(ThisPulseMaxTick); // Tmax
	WFparLimitHigh.push_back(ThisPulseMaxTick+3);
	WFparLimitLow.push_back(ThisPulseMaxTick-3);

	WFparameters.push_back(1.7); // rise time
	WFparLimitHigh.push_back(0);
	WFparLimitLow.push_back(0);
	
	WFparameters.push_back(3.5); // fall time
	WFparLimitHigh.push_back(0);
	WFparLimitLow.push_back(0);

	NpeaksFound++;
      }
    }
    WFparameters[0] = NpeaksFound;
    WFparLimitHigh[0] = 0;
    WFparLimitLow[0] = 0;
    
    FitFcn = new TF1(Form("FitFcn_%d_%d_%d_%d_%d", run, subrun, event, flash, ch),WFshape,0,1500,4*NpeaksFound+2);
    for(UInt_t i = 0;i<WFparameters.size();i++){
      FitFcn->SetParameter(i,WFparameters[i]);
      FitFcn->SetParLimits(i,WFparLimitLow[i],WFparLimitHigh[i]);
      if(i == 0)FitFcn->SetParName(i,"N_{peaks}");
      else if(i == 1)FitFcn->SetParName(i,"baseline");
      else{
	if((i-2)%4==0)FitFcn->SetParName(i,Form("A_{%d}",(i-2)/4));
	else if((i-3)%4==0)FitFcn->SetParName(i,Form("#mu_{%d}",(i-3)/4));
	else if((i-4)/4==0)FitFcn->SetParName(i,Form("#sigma_{r,%d}",(i-4)/4));
	else FitFcn->SetParName(i,Form("sigma_{f,%d}",(i-5)/4));
      }
    }
    FitFcn->SetLineWidth(1);
    FitFcn->SetNpx(3000);
    //ComputeResidual(flash, ch);
    return FitFcn;

  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::IsAmplitudeHighEnough(Int_t flash, Int_t ch){
    double min = 10000;
    double max = 0;
    
    for(int tick = 0; tick < 1500;tick++){
      if(!(tick < FlashesTimes[flash]+flashEndTimeOff && tick > FlashesTimes[flash]-flashStartTimeOff))continue;
      if(waveforms[ch][tick] < min){min = waveforms[ch][tick];}
      if(waveforms[ch][tick] > max){max = waveforms[ch][tick];}
    }
    maxAmplitudeinFlash = max - min;

    if(maxAmplitudeinFlash > 15){return true;}
    else{return false;}
  }

  //
  //****************************************************************
  //
  
  bool PMTwf_ana::IsAmplitudeHighEnough(Int_t ch){
    double min = 10000;
    double max = 0;

    for(int tick = FitRangeMin; tick < FitRangeMax;tick++){
      if(tick < beamWindowStart-10 || tick > beamWindowEnd) continue;
      if(waveforms[ch][tick] < min){min = waveforms[ch][tick];}
      if(waveforms[ch][tick] > max){max = waveforms[ch][tick];}
    }
    maxAmplitudeinFlash = max - min;

    if(maxAmplitudeinFlash > 15){return true;}
    else{return false;}
  }
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::NeedNewPE(Int_t flash, Int_t ch){
    ComputeResidual(flash, ch);
    double residualMax = -1000;
    int tickMax = 0;
    for(int tick = 0;tick<1500;tick++){
      if(WFresidual[tick]>residualMax){residualMax = WFresidual[tick];tickMax = tick;}
    }
    if(WFresidual[tickMax] >= 0.15){
      RecoPulseTimesEvent[ch].push_back(tickMax);
      return true;
    }
    else return false;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::NeedNewPE(Int_t ch){
    ComputeResidual(ch);
    double residualMax = -1000;
    int tickMax = 0;
    for(int tick = FitRangeMin;tick<FitRangeMax;tick++){
      if(WFresidual[tick]>residualMax){residualMax = WFresidual[tick];tickMax = tick;}
    }
    if(WFresidual[tickMax] >= newPEth){
      RecoPulseTimesEvent[ch].push_back(tickMax);
      return true;
    }
    else return false;
  }
  
  //
  //****************************************************************
  //

  void PMTwf_ana::DrawFilter(Int_t ch){
    cEvent->cd(1);
    TH1D *hWF = GetChWF(ch);
    hWF->Draw();

    TH1D *hfiltered = GetfilteredWF(ch);
    hfiltered->SetLineColor(2);
    hfiltered->Draw("same");

    cEvent->cd(2);
    TH1D *htrackBL = new TH1D(Form("htrackBL_%d_%d_%d_%d",run,subrun,event,ch),Form("htrackBL_%d_%d_%d_%d",run,subrun,event,ch),1000,1900,2400);
    for(int i = 0;i<1500;i++){
      htrackBL->Fill(filteredWF[ch][i]);
    }

    htrackBL->Draw();
    TLatex *tex = new TLatex(2100,0.5*htrackBL->GetBinContent(htrackBL->GetMaximumBin()),Form("baseline RMS %f",htrackBL->GetRMS()));
    tex->SetTextAlign(12);
    tex->SetTextSize(0.05);
    tex->Draw();
    TLatex *tex2 = new TLatex(2100,0.4*htrackBL->GetBinContent(htrackBL->GetMaximumBin()),Form("Eval baseline RMS %f",BLrms[ch]));
    tex2->SetTextAlign(12);
    tex2->SetTextSize(0.05);
    tex2->Draw();

    cEvent->Modified();
    cEvent->Update();
    cEvent->SaveAs(Form("plot/FilterWF_%d_%d_%d_%d.png",run,subrun,event,ch));

  }

  //
  //****************************************************************
  //
  
  double PMTwf_ana::FindFlash(){
    hFlashFinding->Reset();
    hFlashFindingDeriv->Reset();
    int time;
    for(int ch = 0;ch < 32;ch++){
      for(UInt_t pulse = 0; pulse < RecoPulseTimesEvent[ch].size();pulse++){
	time = RecoPulseTimesEvent[ch][pulse]+3;
	hFlashFinding->Fill(RecoPulseTimesEvent[ch][pulse],waveforms[ch][time]-Baselines[ch]);
      }
    }

    for(int i = 0;i<hFlashFindingDeriv->GetNbinsX();i++){
      if(i==0)hFlashFindingDeriv->SetBinContent(i+1,0);
      else{hFlashFindingDeriv->SetBinContent(i+1,hFlashFinding->GetBinContent(i+1)-hFlashFinding->GetBinContent(i));}
    }

    double flashPoint = -1;
    hFlashFindingDeriv->GetXaxis()->SetRange((beamWindowStart-10)/10,(beamWindowEnd+50)/10);
    if(hFlashFinding->GetBinContent(hFlashFindingDeriv->GetMaximumBin())>200){
      flashPoint = (hFlashFindingDeriv->GetMaximumBin())*hFlashFindingDeriv->GetBinWidth(1);
    }
    
    return flashPoint;

    cEvent->cd(1);
    hEventWF->Draw("colz");
    TLine *lFlashEvent = new TLine(flashPoint,0,flashPoint,32);
    lFlashEvent->SetLineColor(2);
    lFlashEvent->SetLineWidth(2);
    lFlashEvent->Draw();
    
    cEvent->cd(2);
    hFlashFinding->Draw();
    hFlashFindingDeriv->SetLineColor(2);
    hFlashFindingDeriv->Draw("same");
    TLine *lFlashFind = new TLine(flashPoint,hFlashFinding->GetMinimum(),flashPoint,hFlashFinding->GetMaximum());
    lFlashFind->SetLineColor(2);
    lFlashFind->SetLineWidth(2);
    lFlashFind->Draw();
    cEvent->Modified();
    cEvent->Update();
    cEvent->SaveAs(Form("plot/FlashFinding_%d_%d_%d.png",run,subrun,event));
    
    return flashPoint;
  }

  //
  //****************************************************************
  //

  TF1*  PMTwf_ana::InitializeFitFcn(Int_t ch){
    NpeaksFound = 0;
    WFparameters.clear();
    WFparLimitHigh.clear();
    WFparLimitLow.clear();
    
    WFparameters.push_back(0);
    WFparLimitHigh.push_back(0);
    WFparLimitLow.push_back(0);
    WFparameters.push_back(Baselines[ch]);
    WFparLimitHigh.push_back(Baselines[ch]+10);
    WFparLimitLow.push_back(Baselines[ch]-10);

    int ThisPulseTime;
    int ThisPulseMaxTick = 0;;
    double ThisPulseMaxAmp = 0;

    for(UInt_t pulse = 0;pulse < RecoPulseTimesEvent[ch].size();pulse++){
      ThisPulseTime = RecoPulseTimesEvent[ch][pulse];
      if(ThisPulseTime < FitRangeMin || ThisPulseTime > FitRangeMax)continue;
      if(pulse < RecoPulseTimesEvent[ch].size()-1){
	  ThisPulseMaxTick = ThisPulseTime;
	  ThisPulseMaxAmp = 0;
	  for(int i = ThisPulseTime; i<RecoPulseTimesEvent[ch][pulse+1] ;i++){
	    if(waveforms[ch][i] >= ThisPulseMaxAmp){ThisPulseMaxAmp = waveforms[ch][i];ThisPulseMaxTick = i+0.5;}
	  }
      }
      else{ThisPulseMaxTick = ThisPulseTime+3.5;}
      if(waveforms[ch][ThisPulseMaxTick]-Baselines[ch] < 2)continue;
	
      WFparameters.push_back(0.6*(waveforms[ch][ThisPulseMaxTick]-Baselines[ch])); // amplitude
      //if(ThisPulseMaxTick < WFparameters.at(WFparameters.size()-4)){WFparameters.at(WFparameters.size()-1) = 0.5*(waveforms[ch][ThisPulseMaxTick]-Baselines[ch]);}
      WFparLimitHigh.push_back(1.1*(waveforms[ch][ThisPulseMaxTick]-Baselines[ch]));
      WFparLimitLow.push_back(0);

      WFparameters.push_back(ThisPulseMaxTick); // Tmax
      WFparLimitHigh.push_back(ThisPulseMaxTick+3);
      WFparLimitLow.push_back(ThisPulseMaxTick-3);

      WFparameters.push_back(1.7); // rise time
      WFparLimitLow.push_back(1);
      WFparLimitHigh.push_back(2);
	
      WFparameters.push_back(3.9); // fall time
      WFparLimitLow.push_back(3.7);
      WFparLimitHigh.push_back(4.2);

      NpeaksFound++;
    }
    
    WFparameters[0] = NpeaksFound;
    WFparLimitHigh[0] = 0;
    WFparLimitLow[0] = 0;
    
    FitFcn = new TF1(Form("FitFcn_%d_%d_%d_%d", run, subrun, event, ch),WFshape,0,1500,4*NpeaksFound+2);
    for(UInt_t i = 0;i<WFparameters.size();i++){
      FitFcn->SetParameter(i,WFparameters[i]);
      FitFcn->SetParLimits(i,WFparLimitLow[i],WFparLimitHigh[i]);
      if(i == 0)FitFcn->SetParName(i,"N_{peaks}");
      else if(i == 1)FitFcn->SetParName(i,"baseline");
      else{
	if((i-kAmplitude)%4==0)FitFcn->SetParName(i,Form("A_{%d}",(i-kAmplitude)/4));
	else if((i-kTmax)%4==0)FitFcn->SetParName(i,Form("#mu_{%d}",(i-kTmax)/4));
	else if((i-kSigmaRise)%4==0)FitFcn->SetParName(i,Form("#sigma_{r,%d}",(i-kSigmaRise)/4));
	else FitFcn->SetParName(i,Form("sigma_{f,%d}",(i-kSigmaFall)/4));
      }
    }
    FitFcn->SetLineWidth(1);
    FitFcn->SetNpx(3000);
    //ComputeResidual(flash, ch);
    return FitFcn;
    
  }
  
  //
  //****************************************************************
  //
  
  void PMTwf_ana::SetPElist(Int_t ch){
    std::vector<PEhit> PElistWF;
    PEhit newPE;
    for(Int_t pe = 0; pe < WFparameters[0];pe++){
      newPE.SetValues(FitFcn->GetParameter(4*pe+kAmplitude), FitFcn->GetParameter(4*pe+kTmax),FitFcn->GetParameter(4*pe+kSigmaRise),FitFcn->GetParameter(4*pe+kSigmaFall));
      PElistWF.push_back(newPE);
    }
    for(UInt_t pe= 0;pe<PElistWF.size();pe++){
      PElist[ch].push_back(PElistWF[pe]);
    }
    //PElist[ch] = PElistWF;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::FindFlashT0(){
    bool wellReconstructed = true;

    PEhit firstPE(20,1500,1,1);
    PEhit firstsmallPE(20,1500,1,1);int NPEtot = 0;
    double DTfirstPE = 1/15.625;  // max DT between the first big peak and the first small PE considered                 
    for(UInt_t ch = 0; ch < PElist.size();ch++){NPEtot+=PElist.at(ch).size();}
    if(NPEtot==32){wellReconstructed = false; std::cout << "NPEtot = 32..." << std::endl; return wellReconstructed;}// no more than EXACTLY 1PE per WF => the initialization PE
    
    double Nsigmas   = 1.61; // probably put that in the .h and initialize  
    double threshold = 75.79;// ------------------------------------------

    firstPE.SetValues(20,1500,1,1);
    for(UInt_t ch = 0; ch < PElist.size();ch++){
      for(UInt_t pe = 0;pe < PElist.at(ch).size();pe++){
        if(PElist.at(ch).at(pe).GetAmplitude() > threshold){

          if(PElist.at(ch).at(pe).GetTime()-Nsigmas*PElist.at(ch).at(pe).GetSigmaRise() < firstPE.GetTime()-Nsigmas*firstPE.GetSigmaRise()){
            firstPE.SetValues(PElist.at(ch).at(pe).GetAmplitude(),PElist.at(ch).at(pe).GetTime(),PElist.at(ch).at(pe).GetSigmaRise(),PElist.at(ch).at(pe).GetSigmaFall());
          }
        }
      }
    }

    firstsmallPE.SetValues(20,1500,1,1);
    for(UInt_t ch = 0; ch < PElist.size();ch++){
      for(UInt_t pe = 0;pe < PElist.at(ch).size();pe++){
	if(PElist.at(ch).at(pe).GetTime()-Nsigmas*PElist.at(ch).at(pe).GetSigmaRise() < firstsmallPE.GetTime()-Nsigmas*firstsmallPE.GetSigmaRise() && PElist.at(ch).at(pe).GetTime()-Nsigmas*PElist.at(ch).at(pe).GetSigmaRise() > firstPE.GetTime()-Nsigmas*firstPE.GetSigmaRise()+DTfirstPE){
          firstsmallPE.SetValues(PElist.at(ch).at(pe).GetAmplitude(),PElist.at(ch).at(pe).GetTime(),PElist.at(ch).at(pe).GetSigmaRise(),PElist.at(ch).at(pe).GetSigmaFall());
        }
      }
    }  
    recoT0 = firstsmallPE.GetTime()-Nsigmas*firstsmallPE.GetSigmaRise();
    return wellReconstructed;
  }
  
  //
  //****************************************************************
  //

  void PMTwf_ana::GetBNBTime(){
    TF1 *fBNB = new TF1(Form("fBNB_%d_%d_%02d", run, subrun, event),WFshape,0,1500,6);
    fBNB->SetParameters(1,2050,2048,5,1.5,3);
    fBNB->FixParameter(0,1);
    fBNB->SetNpx(3000);
    cEvent->cd();
    TH1D *hBNB = GetBNBwf();
    hBNB->GetXaxis()->SetRange(1,200);
    hBNB->Fit(Form("%s",fBNB->GetName()),"rl",0,200);
    BNB_recoTime = fBNB->GetParameter(3);
    BNBparameters[0] = fBNB->GetParameter(1);
    BNBparameters[1] = fBNB->GetParameter(2);
    BNBparameters[2] = fBNB->GetParameter(3);
    BNBparameters[3] = fBNB->GetParameter(4);
    BNBparameters[4] = fBNB->GetParameter(5);
    cEvent->Modified();
    cEvent->Update();
  }

  //
  //****************************************************************
  //
  
  void PMTwf_ana::GetMCT0(){
    MC_T0=1500;
    for(int ch = 0;ch<32;ch++){
      for(UInt_t pe = 0;pe < MChitTimes[ch].size();pe++){
	if(MChitTimes[ch][pe] < MC_T0){MC_T0 = MChitTimes[ch][pe];}
      }
    }
  }
  

  //
  //
  //
  void PMTwf_ana::GetReferenceTime(){
    if(isMC){GetMCT0();}
    else{GetBNBTime();} 
  }
  //
  //****************************************************************
  //
}

#endif
