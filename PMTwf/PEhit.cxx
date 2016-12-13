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

void PEhit::operator=(PEhit pe1){
  SetValues(pe1.GetAmplitude(),pe1.GetTime(),pe1.GetSigmaRise(),pe1.GetSigmaFall());
}


//double PEhit::GetTime(){return _time;}
//double PEhit::GetAmplitude(){return _amplitude;}
//double PEhit::GetSigmaRise(){return _sigma_r;}
//double PEhit::GetSigmaFall(){return _sigma_f;}
//double PEhit::GetPEcharge(){return _peCharge;}

#endif
