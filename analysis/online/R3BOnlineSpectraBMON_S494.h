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
// -----              R3BOnlineSpectraBMON_S494           -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BONLINESPECTRABMON_S494
#define R3BONLINESPECTRABMON_S494

#include "FairTask.h"
#include "R3BGlobalAnalysis.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

class R3BOnlineSpectraBMON_S494 : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraBMON_S494();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraBMON_S494(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraBMON_S494();

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
    inline void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }

    /**
     * Methods for setting reset and readout times for Bmon
     */
    inline void SetBmon(Int_t time_range, Double_t time_step, Int_t sens_SEE, Int_t sens_IC)
    {
        reset_time = time_range;
        read_time = time_step;
        fsens_SEE = sens_SEE;
        fsens_IC = sens_IC;
    }

    void Reset_ROLU_Histo();
    void Reset_BMON_Histo();

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;

    enum DetectorInstances
    {
        DET_BMON,
        DET_ROLU,
        DET_TOFD,
        DET_MAX
    };

    const char* fDetectorNames[DET_MAX + 1] = { "BeamMonitor", "Rolu", "Tofd", NULL };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat1, fTpat2;
    Int_t fSamp;
    Double_t fClockFreq; /**< Clock cycle in [ns]. */
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */

    //   TClonesArray *fbmonMappedItems;

    unsigned long long time_start = -1, time_begin = -1, time = 0, time_mem = 0, time_end = 0, time_mem_long = 0;
    unsigned long long time_prev_read = 0, time_to_read = 0;
    unsigned long int ic_mem = 0, see_mem = 0, see_first = 0, tofdor_mem = 0, srolu1_mem = 0, srolu2_mem = 0,
                      nBeamParticle = 0;
    unsigned long int ic_start = 0, see_start = 0, tofdor_start = 0, srolu1_start = 0, srolu2_start = 0;
    unsigned long long time_spill_start = 0, time_spill_end = 0;
    Double_t see_spill = 0, ic_spill = 0, tofdor_spill = 0;
    unsigned long long time_previous_event = 0;
    unsigned long ySEE = 0, ySEE_mem = 0;
    int yIC_mem_mem = 0, ySEE_mem_mem = 0, yTOFDOR_mem_mem = 0, ySROLU1_mem_mem = 0, ySROLU2_mem_mem = 0;

    Int_t maxevent;

    Bool_t newSpill = false;
    Bool_t spill_on = false;
    Double_t time_clear = -1.;
    Double_t tdiff = 0.;
    Double_t tdiff_long = 0.;
    Int_t counterWrongTpat = 0;
    Int_t counterWrongTrigger = 0;
    Double_t fNorm = 1.;
    Int_t iclear_count = 1;
    UInt_t reset_time;         // time after which bmon spectra are reseted
    Double_t read_time;        // step in which scalers are read, in sec
    Int_t fsens_SEE, fsens_IC; // SEETRAM and IC sensitivity, between -4 and -10
    Double_t calib_SEE = 1.;   // SEETRAM calibration factor
    Double_t see_offset = 0.2; // SEETRAM offset in kHz
    Double_t calib_IC = 1.;    // SEETRAM calibration factor
    Double_t ic_offset = 0;    // SEETRAM offset in kHz
    Int_t in_spill_off = 0;
    unsigned long seeLastSpill = 0, time_lastSpill = 0;
    UInt_t num_spills = 0;
    UInt_t num_spills_prev = 0;
    Double_t spilltime = 0.;
    Double_t counts_SEE = 0;
    Double_t counts_IC = 0;
    Double_t counts_TofD = 0;
    unsigned long fNEvents = 0, fNEvents_start = 0, fNEventsRolu = 0, fNSpills = 0; /**< Event counter. */

    TH1F* fh_Tpat;
    TH1F* fh_spill_length;
    TH1F* fh_Trigger;
    TH1F* fh_IC_spill;
    TH1F* fh_SEE_spill;
    TH1F* fh_SEE_spill_long;
    TH1F* fh_SEE_spill_raw;
    TH1F* fh_TOFDOR_spill;
    TH1F* fh_IC;
    TH1F* fh_SEE;
    TH1F* fh_SEEPS;
    TH1F* fh_SEER;
    TH1F* fh_TOFDOR;
    TH1F* fh_SROLU1;
    TH1F* fh_SROLU2;
    TH2F* fh_IC_cal;
    TH2F* fh_SEE_cal;
    TH1F* fh_SROLU1_spill;
    TH1F* fh_SROLU2_spill;
    TH1F* fh_IC_TOFDOR;
    TH1F* fh_IC_SEE;
    TH1F* fh_SEE_TOFDOR;

    TH2F* fh_rolu_tot;
    TH2F* fh_rolu_tof;
    TH1F* fh_rolu_channels;

  public:
    ClassDef(R3BOnlineSpectraBMON_S494, 1)
};

#endif
