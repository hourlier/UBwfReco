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


#endif
