/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ------------------------------------------------------------
// -----              R3BOnlineSpillAnalysis              -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BONLINESPILLANALYSIS_H
#define R3BONLINESPILLANALYSIS_H

#include "FairTask.h"

#include "TMath.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2F;
class R3BEventHeader;
class TCanvas;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BOnlineSpillAnalysis : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpillAnalysis();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpillAnalysis(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpillAnalysis();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }
    inline void SetSpillLength(Double_t SpillLength) { fSpillLength = SpillLength; }

    /**
     * Methods for setting number of planes and paddles
     */
    /* Put Custom functions here*/
    /** Function for exponential distribution*/
    double ran_expo(double lambda)
    {
        double u;

        u = rand() / (RAND_MAX + 1.0);

        return -log(1 - u) / lambda;
    }

    void Reset_Histo();
    void Update_Histo();

  private:
    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    Int_t fSamp;
    Double_t fSpillLength;
    TClonesArray* fSamplerMappedItems;

    TCanvas* cSpill;
    TCanvas* cFFT;

    unsigned long long time_start = 0, time = 0;

    unsigned long long time_spill_start = 0, time_spill_end = 0;
    long samplerPrev = 0, samplerCurr = 0, samplerSpill = 0, samplerBegin = 0;
    UInt_t spillCounter = 0;

    Bool_t spill_on = false;
    Bool_t spill_off = false;
    Bool_t spill_on_sampler = false;
    Double_t time_clear = -1.;
    Double_t tdiff = 0.;
    Double_t fNorm = 1.;
    Int_t iclear_count = 1;
    UInt_t reset_time;  // time after which bmon spectra are reseted
    Double_t read_time; // step in which scalers are read, in sec
    Double_t fDAQ_dead_time;

    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    // spill monitor
    Double_t LastWindowTime = 0.;
    Double_t tLastInWindow = 0.;
    Double_t time_global = 0.;
    Int_t FFT_Counter;
    Int_t xHits = 20;
    Double_t Schwelle = 90.;
    Double_t xmax;
    Double_t TIME_previous_FL = 0.;
    Double_t dtime_min = 80.;
    Double_t dtime_max = 95.;
    Int_t nHits_prev = 0;
    Double_t time_tmp = 0.;
    Int_t nHit_ctr = 0;
    Int_t count_dead = 0;
    Double_t tLastPrev = 0.;
    Double_t spill_times = 0.;
    Double_t spill_times_first = 0.;
    Int_t spill_ctr = 0;
    Int_t spill_ctr_LOS = 0;
    Double_t dt_LOS = 0.;
    Double_t dt_LOS_first = 0.;
    Int_t hit_counter = 0;
    Int_t SPILLS = 0;
    int a_ctr = 0;
    Int_t packagectr = 0;
    Int_t pHits = 0;
    Double_t F_duty = 0.;
    Double_t Hans_mean = 0.;
    Bool_t MissedSpillEnd = false;
    Bool_t spill_off_calc = false;
    Int_t hans_ctr = 0;
    Double_t Hans_check = 0.;
    Int_t short_ctr = 0;
    Double_t dt_prev = 0.;
    Double_t dt_dead_mean = 0.;
    Int_t dead_ctr = 0;
    Double_t deadTime = fDAQ_dead_time;
    Double_t tt_prev = 0.;
    Bool_t skip_spill = false;
    Double_t testcounter = 0.;
    Int_t numHits1a = 0;
    Int_t numHits1b = 0;
    Bool_t first_in_spill = true;
    Double_t CTAN = 0.;
    Double_t CTANN = 0.;
    Double_t CTF = 0.;
    Int_t PPSF = 0;

    TH1F* fh_spill_times;
    TH1F* fh_spill_times_Fine;
    TH1F* fh_spill_times_FFT;
    TH1F* fh_spill_times_Coarse;
    TH1F* fh_spill_times_Fine_adj;
    TH1F* fh_SAMP_tDiff;
    TH1F* fh_SAMP_tDiff_long;
    TH1F* fh_SAMP_tDiff_pois;
    TH1F* fh_dt_hits;
    TH1F* fh_SAMP_freq;
    TH1F* fh_SAMP_freq_long;
    TH1F* fh_rate;
    TH1F* fh_rate_sum;
    TH1F* fh_DutyFactor;
    TH1F* fh_DutyFactor_pois;
    TH1F* fh_DutyFactor_SAMP_clean;
    TH1F* fh_DutyFactor_MtA;
    TH1F* fh_DutyFactor_Avg;
    TH1F* fh_DutyFactor_Max;
    TH1F* fh_DutyFactor_MaxToAvg;
    TH1F* fh_DutyFactor_MaxRun;
    TH1F* fh_DutyFactor_AvgRun;
    TH1F* fh_DutyFactor_PLD;
    TH1F* fh_MtA_sum;
    TH1* fh_FFT;
    TH1* fh_FFT_adj;
    TH1* fh_FFT_add;
    TH1F* fh_spill_hans;
    TH1F* fh_spill_hans_byMax;
    TH1F* fh_spill_hans_sum_byMax;
    TH1F* fh_hans_sum;
    TH1F* fh_hans_sum_byMax;

    TH1F* fh_spill_times_pois;

  public:
    ClassDef(R3BOnlineSpillAnalysis, 2)
};

#endif /* R3BONLINESPILLANALYSIS_H */
