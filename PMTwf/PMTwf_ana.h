/**
 * \file PMTwf_ana.h
 *
 * \ingroup PMTwf
 * 
 * \brief Class def header for a class PMTwf_ana
 *
 * @author hourlier
 */

/** \addtogroup PMTwf

    @{*/

#ifndef LARLITE_PMTWF_ANA_H
#define LARLITE_PMTWF_ANA_H

#include "Analysis/ana_base.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <vector>
#include <array>

namespace larlite {
  /**
     \class PMTwf_ana
     User custom analysis class made by SHELL_USER_NAME
   */
  class PMTwf_ana : public ana_base{
  
  public:
    PMTwf_ana(){ _name="PMTwf_ana"; _fout=0;}
    virtual ~PMTwf_ana(){}

    virtual bool initialize();
    virtual bool analyze(storage_manager* storage);
    virtual bool finalize();

    bool  GetChannelGeo();
    void  DrawChPosition();
    void  FilterWaveform(Int_t ch);
    inline void  SetDTHitReco(double dTreco){DTHitReco = dTreco;}
    inline double  GetDTHitReco(){return DTHitReco;}
    bool  DrawEvent();
    TH1D* GetBNBwf();
    TH1D* GetChWF(Int_t ch);
    TH1D* GetfilteredWF(Int_t ch);
    TH1D* GetDerivative1(Int_t ch);
    TH1D* GetDerivative2(Int_t ch);
    bool  GetAllDerivatives();
    bool  GetDerivativeWF(Int_t ch);
    bool  Get2ndDerivativeWF(Int_t ch);
    bool  FindPulses(Int_t ch);
    bool  FindAllPulses();
    bool  LoadEventInfo(storage_manager* storage);
    bool  DrawFlashWF(Int_t flash, Int_t ch);
    bool  FitFlashWF(Int_t flash, Int_t ch);
    bool  FitFlash(Int_t flash);
    bool  EvalBaseline(Int_t ch);
    bool  GetBaselines();
    bool  InitializeFitFcn(Int_t flash, Int_t ch);
    bool  IsFlashInBeam(Int_t flash);
    bool  IsFlashBeforeBeam(Int_t flash);
    bool  IsSummedWFHighbeforeBeam();

  protected:
    // ROOT TObject members
    TH1D               *hBNBwf;
    TH1D               *hWFtmp;
    TH1D               *hDerivativeWFtmp;
    TH1D               *h2ndDerivativeWFtmp;
    TH2D               *hEventWF;
    TCanvas            *cEvent;
    TCanvas            *cIndivWF;
    TPad               *pEvent2D;
    TPad               *pWaveform;
    TPad               *pDerivative;
    TPad               *pFullWF;
    TH1D               *hSigma_rise_BNB;
    TH1D               *hSigma_fall_BNB;
    TH1D               *hTime_BNB;
    TF1                *FitFcn;

    std::vector<double> FlashesTimes;             // times of flashes for an event
    std::vector<double> HitTimes;                 // ophit times for a given WF
    std::vector<double> AllHitTimes[32];          // array of 32 vectors of the hit times for each waveforms
    std::vector<double> RecoPulseTimes;           // times of the peaks found by FindPulses() for one WF
    std::vector<double> RecoPulseTimesEvent[32];  // times of the peaks found by FindPulses() for all WF
    std::vector<double> PEwf1flash;               // PE per WF for 1 flash
    std::vector< std::vector<double> > PEflashWF; // PE per waveform for each flash
    
    // Fit parameters
    Int_t               kBaseline;                // baseline index
    Int_t               kTmax;                    // Tmax index
    Int_t               kSigmaRise;               // rise time index
    Int_t               kSigmaFall;               // fall time index
    Int_t               NpeaksFound;              // Full number of pulses found
    std::vector<double> WFparameters;             //all parameters for a given WF
    std::vector< std::vector<double> > AllWFparameters[32];// all parameters for a given event
    

    double              xyzCh[32][3];             // Position of the optical chanels
    double              DTHitReco;                // Time shift for comparing Treco to Tophit
    double              waveform[1500];           // individual waveform
    double              SummedWF[1500];           // allWF summed into one
    double              filteredWF[32][1500];     // Filterout the HF(?)
    double              WFderiv[1500];            // individual WF 1st derivative
    double              WF2ndDeriv[1500];         // individual WF 2nd derivative
    double              waveforms[32][1500];      // array with all waveforms
    double              WFderivs[32][1500];       // array with all WF derivatives
    double              WF2ndDerivs[32][1500];    // array with all WF 2nd derivatives
    double              BNBwaveform[1500];        // BNB waveform
    double              Baselines[32];            // baselines
    double              BLrms[32];                // baselineRMS
    bool                isMC;                     // is this MC (e.g. is there a BNB data?)
    bool                drawOutput;               // should I draw each reconstructed WF?
    bool                printOK;                  // should I print the output in a pdf file?
    double              flashStartTime;           // start time of the current flash
    double              flashEndTime;             // end time of the current flash
    double              flashStartTimeOff;        // offset for flash window determination
    double              flashEndTimeOff;          // offset for flash window determination
    int                 run;                      // current run number
    int                 subrun;                   // current subrun number
    int                 event;                    // current event number
    double              beamWindowStart;          // start of beam window
    double              beamWindowEnd;            // end of beam window
    

  };
}

#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
