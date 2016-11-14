#ifndef PEHIT_CXX
#define PEHIT_CXX
#include "PEhit.h"

PEhit::PEhit(){
  _time = -1;
  _amplitude = 22;
  _sigma_r = 1.5;
  _sigma_f = 3;
}

PEhit::PEhit(double Amp, double t, double sigma_rise, double sigma_fall){
  _time= t;
  _amplitude = Amp;
  _sigma_r = sigma_rise;
  _sigma_f = sigma_fall;
}

void PEhit::SetTime(double t){_time = t;}
void PEhit::SetAmplitude(double Amp){_amplitude = Amp;}
void PEhit::SetSigmas(double sigma_rise, double sigma_fall){_sigma_r = sigma_rise;_sigma_f = sigma_fall;}
void PEhit::SetValues(double Amp, double t, double sigma_rise, double sigma_fall){
  SetAmplitude(Amp);
  SetTime(t);
  SetSigmas(sigma_rise,sigma_fall);
}


double PEhit::GetTime(){return _time;}
double PEhit::GetAmplitude(){return _amplitude;}
double PEhit::GetSigmaRise(){return _sigma_r;}
double PEhit::GetSigmaFall(){return _sigma_f;}


//*****************************************************

wfInfo::wfInfo(){
  _baseline = -1;
  _baselineRMS = -1;
}
wfInfo::wfInfo(double baseline, double baselineRMS, std::vector<PEhit> PElist){
  _baseline = baseline;
  _baselineRMS = baselineRMS;
  _PElist = PElist;
}
wfInfo::~wfInfo(){}
void   wfInfo::SetBaseline(double baseline, double baselineRMS){_baseline = baseline;_baselineRMS = baselineRMS;}
void   wfInfo::SetPElist(std::vector<PEhit> PElist){_PElist = PElist;}
void   wfInfo::AddPE(PEhit newPE){_PElist.push_back(newPE);}
void   wfInfo::RemovePE(int PEindex){_PElist.erase(_PElist.begin()+PEindex-1);}
int    wfInfo::GetNPE(){return _PElist.size();}
double wfInfo::GetBaseline(){return _baseline;}
double wfInfo::GetBaselineRMS(){return _baselineRMS;}
std::vector<PEhit> wfInfo::GetPElist(){return _PElist;}
PEhit  wfInfo::GetPE(int PEindex){
  if(PEindex < (int)_PElist.size())return _PElist.at(PEindex);
  else{
    PEhit emptyPE;
    return emptyPE;
  }
}
#endif
