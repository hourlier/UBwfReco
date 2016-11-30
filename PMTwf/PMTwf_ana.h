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
#include "TTree.h"
#include <vector>
#include <array>
#include "PEhit.h"

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

    void    initializePElist();
    void    CleanPElist();
    bool    GetChannelGeo();
    void    DrawChPosition();
    void    FilterWaveform(Int_t ch);
    void    FilterWaveformOld(Int_t ch);
    inline void  SetDTHitReco(double dTreco){DTHitReco = dTreco;}
    inline double  GetDTHitReco(){return DTHitReco;}
    bool    DrawEvent();
    TH1D*   GetBNBwf();
    TH1D*   GetChWF(Int_t ch);
    TH1D*   GetfilteredWF(Int_t ch);
    TH1D*   GetDerivative1(Int_t ch);
    TH1D*   GetDerivative2(Int_t ch);
    bool    PMTpreCut();
    bool    GetAllDerivatives();
    bool    GetDerivativeWF(Int_t ch);
    bool    Get2ndDerivativeWF(Int_t ch);
    bool    FindPulses(Int_t ch);
    bool    FindAllPulses();
    bool    LoadEventInfo(storage_manager* storage);
    bool    DrawFlashWF(Int_t flash, Int_t ch);
    void    DrawWF(Int_t ch);
    bool    FitFlashWF(Int_t flash, Int_t ch);
    bool    FitFlash(Int_t flash);
    void    FitWF(Int_t ch);
    void    IterateFitProcess(Int_t ch);
    bool    EvalBaseline(Int_t ch);
    void    ComputeResidual(Int_t flash, Int_t ch); 
    void    ComputeResidual(Int_t ch);
    bool    NeedNewPE(Int_t flash, Int_t ch);
    bool    NeedNewPE(Int_t ch);
    bool    GetBaselines();
    void    SetPElist(Int_t ch);
    TF1*    InitializeFitFcn(Int_t flash, Int_t ch);
    TF1*    InitializeFitFcn(Int_t ch);
    bool    IsFlashInBeam(Int_t flash);
    bool    IsFlashBeforeBeam(Int_t flash);
    bool    IsSummedWFHighbeforeBeam();
    bool    IsAmplitudeHighEnough(Int_t flash, Int_t ch);
    bool    IsAmplitudeHighEnough(Int_t ch);
    void    DrawFilter(Int_t ch);
    double  FindFlash();
    bool    FindFlashT0();
    void    GetReferenceTime();
    void    GetBNBTime();
    void    GetMCT0();

  protected:
    // ROOT TObject members
    TH1D               *hBNBwf;
    TH1D               *hWFtmp;
    TH1D               *hFlashT0;
    TH1D               *hDerivativeWFtmp;
    TH1D               *h2ndDerivativeWFtmp;
    TH1D               *hDeltaT0;
    TH1D               *hPreCut;
    TH2D               *hEventWF;
    TCanvas            *cEvent;
    TCanvas            *cIndivWF;
    TCanvas            *cEventT0;
    TPad               *pEvent2D;
    TPad               *pWaveform;
    TPad               *pDerivative;
    TPad               *pFullWF;
    TH1D               *hFlashFinding;
    TH1D               *hFlashFindingDeriv;
    TF1                *FitFcn;

    TTree              *T;

    std::vector<double> FlashesTimes;             // times of flashes for an event from opflash
    std::vector<double> RecoFlashesTimes;         // times of flashes for an event that I reco
    std::vector<double> HitTimes;                 // ophit times for a given WF
    std::vector<double> AllHitTimes[32];          // array of 32 vectors of the hit times for each waveforms
    std::vector<double> RecoPulseTimes;           // times of the peaks found by FindPulses() for one WF
    std::vector<double> RecoPulseTimesEvent[32];  // times of the peaks found by FindPulses() for all WF
    std::vector<double> PEwf1flash;               // PE per WF for 1 flash
    std::vector< std::vector<double> > PEflashWF; // PE per waveform for each flash
    std::vector< std::vector<PEhit> >  PElist;
    std::vector<wfInfo> RecoWFinfo;               // full information on the waveforms
    
    
    // Fit parameters
    Int_t               kNPE;                     // index of the number of PE
    Int_t               kAmplitude;               // index fo amplitude
    Int_t               kBaseline;                // baseline index
    Int_t               kTmax;                    // Tmax index
    Int_t               kSigmaRise;               // rise time index
    Int_t               kSigmaFall;               // fall time index
    Int_t               NpeaksFound;              // Full number of pulses found
    std::vector<double> WFparameters;             // all parameters for a given WF
    std::vector<double> WFparLimitHigh;           // higher limits
    std::vector<double> WFparLimitLow;            // Lower limits
    std::vector< std::vector<double> > AllWFparameters[32];// all parameters for a given event
    std::vector<double> MChitTimes[32];           // time for each MC photons
    std::vector<int>    MCparticles;              // list of MC particles in event_mctrack
    

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
    double              WFresidual[1500];         // residual for evaluating the fit performance
    double              Baselines[32];            // baselines
    double              BLrms[32];                // baselineRMS
    double              BNB_recoTime;              // time reconstructed for the BNB pulse
    bool                isMC;                     // is this MC (e.g. is there a BNB data?)
    bool                drawOutput;               // should I draw each reconstructed WF?
    bool                printOK;                  // should I print the output in a pdf file?
    bool                T0wellReconstructed;      // is T0 well reconstruected?
    bool                is20PE;                   // does this event satisfies the 20PE limit?
    bool                isFlashBeforBeam;         // is there a flash before the beam?
    double              flashStartTime;           // start time of the current flash
    double              flashEndTime;             // end time of the current flash
    double              flashStartTimeOff;        // offset for flash window determination
    double              flashEndTimeOff;          // offset for flash window determination
    double              newPEth;
    double              recoFlashTime;
    double              recoT0;
    double              MC_T0;                    // T0 in the waveforms
    double              MC_Neutrino_vertex[4];    // Neutrino interaction vertex from MC
    bool                MC_isCCQE;                // is CCQE?
    int                 run;                      // current run number
    int                 subrun;                   // current subrun number
    int                 event;                    // current event number
    double              beamWindowStart;          // start of beam window
    double              beamWindowEnd;            // end of beam window
    double              errorLevel;               // error introduced to fit the WF
    double              maxAmplitudeinFlash;
    int                 FitRangeMin;
    int                 FitRangeMax;
    double              BNBparameters[5];         // baseline, amplitude, time, sigmas (rise/fall) 
    

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
