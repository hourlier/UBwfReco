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
  double GetTime();
  double GetAmplitude();
  double GetSigmaRise();
  double GetSigmaFall();

 private:
  double _time;
  double _amplitude;
  double _sigma_r;
  double _sigma_f;
};

class wfInfo : public std::vector<PEhit>{
 
 public:
  
  wfInfo();
  wfInfo(double baseline, double baselineRMS, std::vector<PEhit> PElist);
  ~wfInfo();
  void               SetBaseline(double baseline, double baselineRMS);
  void               SetPElist(std::vector<PEhit> PElist);
  void               AddPE(PEhit newPE);
  void               RemovePE(int PEindex);
  int                GetNPE();
  double             GetBaseline();
  double             GetBaselineRMS();
  std::vector<PEhit> GetPElist();
  PEhit              GetPE(int PEindex);

 private:
  
  double _baseline;
  double _baselineRMS;
  std::vector<PEhit> _PElist;
  
};

#endif
