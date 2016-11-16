#ifndef PEHIT_CXX
#define PEHIT_CXX

#include "PEhit.h"
#include "TMath.h"

PEhit::PEhit(){
  _time = -1;
  _amplitude = 22;
  _sigma_r = 1.5;
  _sigma_f = 3;
  _peCharge = 0.5*sqrt(2*TMath::Pi())*_amplitude*(_sigma_r+_sigma_f);
}

PEhit::PEhit(double Amp, double t, double sigma_rise, double sigma_fall){
  _time= t;
  _amplitude = Amp;
  _sigma_r = sigma_rise;
  _sigma_f = sigma_fall;
  _peCharge = 0.5*sqrt(2*TMath::Pi())*_amplitude*(_sigma_r+_sigma_f);
}

void PEhit::SetTime(double t){_time = t;}
void PEhit::SetAmplitude(double Amp){_amplitude = Amp;_peCharge = 0.5*sqrt(2*TMath::Pi())*_amplitude*(_sigma_r+_sigma_f);}
void PEhit::SetSigmas(double sigma_rise, double sigma_fall){_sigma_r = sigma_rise;_sigma_f = sigma_fall;_peCharge = 0.5*sqrt(2*TMath::Pi())*_amplitude*(_sigma_r+_sigma_f);}
void PEhit::SetValues(double Amp, double t, double sigma_rise, double sigma_fall){
  SetAmplitude(Amp);
  SetTime(t);
  SetSigmas(sigma_rise,sigma_fall);
  _peCharge = 0.5*sqrt(2*TMath::Pi())*_amplitude*(_sigma_r+_sigma_f);
}


double PEhit::GetTime(){return _time;}
double PEhit::GetAmplitude(){return _amplitude;}
double PEhit::GetSigmaRise(){return _sigma_r;}
double PEhit::GetSigmaFall(){return _sigma_f;}
double PEhit::GetPEcharge(){return _peCharge;}

//
//*****************************************************
//

wfInfo::wfInfo(){
  _baseline = -1;
  _baselineRMS = -1;
}
wfInfo::wfInfo(double baseline, double baselineRMS, std::vector<PEhit> PElist){
  SetWF(baseline,baselineRMS,PElist);
}
wfInfo::~wfInfo(){}
void   wfInfo::SetWF(double baseline, double baselineRMS, std::vector<PEhit> PElist){
  SetBaseline(baseline,baselineRMS);
  SetPElist(PElist);
}

void   wfInfo::SetBaseline(double baseline, double baselineRMS){_baseline = baseline;_baselineRMS = baselineRMS;}

void   wfInfo::SetPElist(std::vector<PEhit> PElist){
  _chargeReconstructed = 0;
  for(unsigned int pe = 0; pe<PElist.size();pe++){
    _PElist.push_back(PElist[pe]);
    _chargeReconstructed+= PElist[pe].GetPEcharge();
  }
}

void   wfInfo::AddPE(PEhit newPE){
  _PElist.push_back(newPE);
  _chargeReconstructed+= newPE.GetPEcharge();
}

void   wfInfo::RemovePE(int PEindex){
  _chargeReconstructed-= _PElist[PEindex].GetPEcharge();
  _PElist.erase(_PElist.begin()+PEindex-1);
}

int    wfInfo::GetNPE(){return _PElist.size();}
double wfInfo::GetBaseline(){return _baseline;}
double wfInfo::GetBaselineRMS(){return _baselineRMS;}
double wfInfo::GetRecoCharge(){return _chargeReconstructed;}
std::vector<PEhit> wfInfo::GetPElist(){return _PElist;}
PEhit  wfInfo::GetPE(unsigned int PEindex){
  if(PEindex < _PElist.size())return _PElist.at(PEindex);
  else{
    PEhit emptyPE;
    return emptyPE;
  }
}
#endif
