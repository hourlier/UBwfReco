#ifndef PEHIT_H
#define PEHIT_H

#include <vector>
//#include "PEhit.cxx"


class PEhit{
 public:
  PEhit();
  PEhit(double Amp, double t, double sigma_rise, double sigma_fall);
  inline ~PEhit(){}
  void SetValues(double Amp, double t, double sigma_rise, double sigma_fall);
  void SetTime(double t);
  void SetAmplitude(double Amp);
  void SetSigmas(double sigma_rise, double sigma_fall);
  double GetTime(){return _time;};
  double GetAmplitude(){return _amplitude;};
  double GetSigmaRise(){return _sigma_r;};
  double GetSigmaFall(){return _sigma_f;};
  double GetPEcharge(){return _peCharge;};
  void  operator=(PEhit pe1);

 private:
  double _time;
  double _amplitude;
  double _sigma_r;
  double _sigma_f;
  double _peCharge;
};


#endif
