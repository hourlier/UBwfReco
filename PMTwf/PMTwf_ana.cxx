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
#include <vector>

#include "DataFormat/opdetwaveform.h"
#include "DataFormat/ophit.h"
#include "DataFormat/opflash.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"
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
  value = TMath::Min(value, 4096.);
  return value;
}


namespace larlite {
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::initialize() {
    printOK = true;
    //drawOutput = false;
    drawOutput = true;
    if(!drawOutput){printOK = false;}
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
    hSigma_rise_BNB = new TH1D("hSigma_rize_BNB","hSigma_rize_BNB",400,0,5);
    hSigma_fall_BNB = new TH1D("hSigma_fall_BNB","hSigma_fall_BNB",400,0,5);
    hTime_BNB = new TH1D("hTime_BNB","hTime_BNB",100,0,10);
    hWFtmp = new TH1D("hWFtmp","hWFtmp",1500,0,1500);
    hDerivativeWFtmp = new TH1D("hDerivativeWFtmp","hDerivativeWFtmp",1500,0,1500);
    h2ndDerivativeWFtmp = new TH1D("h2ndDerivativeWFtmp","h2ndDerivativeWFtmp",1500,0,1500);
    //cEvent = new TCanvas("cEvent","cEvent",1000,700);
    cIndivWF = new TCanvas("cIndivWF","cIndivWF",1600,900);
    if(printOK)cIndivWF->Print("waveforms.pdf[");
    flashStartTimeOff = 50;
    flashEndTimeOff = 150;
    beamWindowStart = 210;
    beamWindowEnd = 340;
    return true;
  }
  
  //
  //****************************************************************
  //

  bool PMTwf_ana::analyze(storage_manager* storage) {
    LoadEventInfo(storage);

    for(UInt_t flash = 0;flash < FlashesTimes.size(); flash++){
      if(!IsFlashInBeam(flash))continue;        // process only flashes that are in the beam
      if(IsFlashBeforeBeam(flash)) return true; // process only events that have no flash before the beam window
      if(IsSummedWFHighbeforeBeam()) return true;
      FindAllPulses();
      GetBaselines();
      for(Int_t ch = 0;ch < 32;ch++){
	//if(PEflashWF[flash][ch] > 2){
	  InitializeFitFcn(flash, ch);
	  if(drawOutput)DrawFlashWF(flash,ch);
	  //}
      }
    }
    

    return true;
  }

  //
  //****************************************************************
  //
  
  bool PMTwf_ana::finalize() {
    //DrawChPosition();
    if(printOK)cIndivWF->Print("waveforms.pdf]");
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::IsSummedWFHighbeforeBeam(){

    double wfmin(10000), wfmax(0);
    for(int tick = 0; tick < beamWindowStart-10; tick++){
      if(SummedWF[tick] > wfmax){wfmax = SummedWF[tick];}
      if(SummedWF[tick] < wfmin){wfmin = SummedWF[tick];}
    }
    std::cout << wfmax - wfmin << std::endl;
    if(wfmax - wfmin > 100)return true;
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

    pFullWF->cd();
    TH1D *hwfFull = GetChWF(ch);
    hwfFull->GetXaxis()->SetRange(1,1500);
    hwfFull->DrawClone();
    FilterWaveform(ch);
    TH1D *hfiltered = GetfilteredWF(ch);
    hfiltered->SetLineColor(2);
    hfiltered->Draw("same");

    pWaveform->cd();
    hWFtmp = GetChWF(ch);
    hWFtmp->GetXaxis()->SetRange(TMath::Max(1.,FlashesTimes[flash]-flashStartTimeOff),TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.));
    hWFtmp->Draw();
    TLine *lHit;
    for(UInt_t hit = 0;hit < AllHitTimes[ch].size();hit++){
      if(AllHitTimes[ch][hit] > TMath::Max(0.,FlashesTimes[flash]-flashStartTimeOff) && AllHitTimes[ch][hit] < TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.)){
	lHit = new TLine(AllHitTimes[ch][hit],hWFtmp->GetMinimum(),AllHitTimes[ch][hit],hWFtmp->GetMaximum());
	lHit->SetLineColor(2);
	lHit->SetLineWidth(2);
	lHit->Draw();
      }
    }
    TLine *lReco;
    for(UInt_t peak = 0; peak < RecoPulseTimesEvent[ch].size();peak++){
      if(RecoPulseTimesEvent[ch][peak] > TMath::Max(0.,FlashesTimes[flash]-flashStartTimeOff) && RecoPulseTimesEvent[ch][peak] < TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.)){
	lReco = new TLine(RecoPulseTimesEvent[ch][peak],hWFtmp->GetMinimum(),RecoPulseTimesEvent[ch][peak],hWFtmp->GetMaximum());
	lReco->SetLineColor(4);
	lReco->SetLineWidth(1);
	lReco->Draw();
      }
    }
    FitFcn->SetNpx(1500);
    FitFcn->Draw("same");
    
    TLine *lBaseline = new TLine(TMath::Max(1.,FlashesTimes[flash]-flashStartTimeOff),Baselines[ch],TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.),Baselines[ch]);
    lBaseline->SetLineColor(2);
    lBaseline->SetLineWidth(2);
    lBaseline->SetLineStyle(2);
    pWaveform->cd();
    lBaseline->Draw();
    
    pDerivative->cd();
    hDerivativeWFtmp = GetDerivative1(ch);
    hDerivativeWFtmp->GetXaxis()->SetRange(TMath::Max(1.,FlashesTimes[flash]-flashStartTimeOff),TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.));
    hDerivativeWFtmp->Draw();
    for(UInt_t peak = 0; peak < RecoPulseTimesEvent[ch].size();peak++){
      if(RecoPulseTimesEvent[ch][peak] > TMath::Max(0.,FlashesTimes[flash]-flashStartTimeOff) && RecoPulseTimesEvent[ch][peak] < TMath::Min(FlashesTimes[flash]+flashEndTimeOff,1500.)){
	lReco= new TLine(RecoPulseTimesEvent[ch][peak],0.95*hDerivativeWFtmp->GetMinimum(),RecoPulseTimesEvent[ch][peak],1.05*hDerivativeWFtmp->GetMaximum());
	lReco->SetLineColor(4);
	lReco->SetLineWidth(1);
	lReco->Draw();
      }
    }

    cIndivWF->Modified();
    cIndivWF->Update();
    if(printOK)cIndivWF->Print("waveforms.pdf");
    return true;
  }

  //
  //****************************************************************
  //

  void PMTwf_ana::FilterWaveform(Int_t ch){
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
	RecoPulseTimes.push_back(i);
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
    //auto mcdata        = storage->get_data<event_mctrack>("mcreco");
    auto mc            = storage->get_data<event_mctruth>("generator");
    auto opdetwaveform = storage->get_data<event_opdetwaveform>("pmtreadout");
    
    if(mc){isMC = true;}else{isMC = false;}
    
    std::cout << "Run #" << storage->run_id() << "_" << storage->subrun_id() << "\t event #" << storage->event_id() << std::endl;
    run = storage->run_id();
    subrun = storage->subrun_id();
    event = storage->event_id();

    //GetWaveforms
    for(Int_t i = 0;i<1500;i++){SummedWF[i] = 0;}
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
    
    // GetHitInfo
    for(int wf = 0; wf < 32; wf++){
      HitTimes.clear();
      for(UInt_t hit = 0; hit < ophits->size();hit++){
	if(ophits->at(hit).OpChannel() == wf){
	  HitTimes.push_back(ophits->at(hit).PeakTime()*64-DTHitReco);
	}
      }
      AllHitTimes[wf] = HitTimes;
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
    TH1D *hWFfit = GetChWF(ch);
    for(int tick = 0;tick < 1500;tick++){hWFfit->SetBinError(tick+1,0.3);}
    hWFfit->GetXaxis()->SetRange(FlashesTimes[flash]);
    return true;
  }

  //
  //****************************************************************
  //

  bool PMTwf_ana::GetBaselines(){
    for(int ch = 0;ch<32;ch++){
      EvalBaseline(ch);
      //std::cout << ch << "\t" << Baselines[ch] << "\t" << BLrms[ch] << std::endl;
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

    TH1D *hBaseline = new TH1D("hBaseline","hBaseline",100,2000,2100);
    
    for(int tick = 0;tick<1500;tick++){
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

  bool  PMTwf_ana::InitializeFitFcn(Int_t flash, Int_t ch){
    NpeaksFound = 0;
    WFparameters.clear();
    WFparameters.push_back(0);
    WFparameters.push_back(Baselines[ch]);
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
	    if(waveforms[ch][i] >= ThisPulseMaxAmp){ThisPulseMaxAmp = waveforms[ch][i];ThisPulseMaxTick = i;}
	  }
	}
	else{ThisPulseMaxTick = ThisPulseTime+3;}
	
	WFparameters.push_back(waveforms[ch][ThisPulseMaxTick]-Baselines[ch]); // amplitude
	WFparameters.push_back(ThisPulseMaxTick); // Tmax
	WFparameters.push_back(1.5); // rise time
	WFparameters.push_back(3); // fall time

	NpeaksFound++;
      }
    }
    WFparameters[0] = NpeaksFound;
    
    FitFcn = new TF1(Form("FitFcn_%d_%d_%d_%d_%d", run, subrun, event, flash, ch),WFshape,0,1500,4*NpeaksFound+2);
    for(UInt_t i = 0;i<WFparameters.size();i++){
      FitFcn->SetParameter(i,WFparameters[i]);
    }
    FitFcn->SetLineWidth(1);
    return true;

  }

  //
  //****************************************************************
  //
}

#endif
