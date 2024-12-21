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
// -----                  R3BOnlineSpectra                -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BONLINESPECTRA
#define R3BONLINESPECTRA

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

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectra(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectra();

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
     * Methods for setting position offset and effective velocity of light
     */
    inline void SetLosXYMCFD(Double_t offsetX, Double_t offsetY, Double_t veffX, Double_t veffY)
    {
        flosOffsetXV = offsetX;
        flosOffsetYV = offsetY;
        flosVeffXV = veffX;
        flosVeffYV = veffY;
    }

    inline void SetLosXYTAMEX(Double_t offsetXT, Double_t offsetYT, Double_t veffXT, Double_t veffYT)
    {
        flosOffsetXT = offsetXT;
        flosOffsetYT = offsetYT;
        flosVeffXT = veffXT;
        flosVeffYT = veffYT;
    }

    inline void SetLosXYToT(Double_t offsetXQ, Double_t offsetYQ, Double_t veffXQ, Double_t veffYQ)
    {
        flosOffsetXQ = offsetXQ;
        flosOffsetYQ = offsetYQ;
        flosVeffXQ = veffXQ;
        flosVeffYQ = veffYQ;
    }
    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

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
    /**
     * Methods for setting number of planes and paddles
     */
    inline void SetNofModules(Int_t planes, Int_t ppp)
    {
        fNofPlanes = planes;
        fPaddlesPerPlane = ppp;
    }
    void Reset_ROLU_Histo();
    void Reset_LOS_Histo();
    void Reset_SCI8_Histo();
    void Reset_BMON_Histo();
    void Reset_TOFD_Histo();
    void Reset_FIBERS_Histo();

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;

    enum DetectorInstances
    {
        DET_AMS,
        DET_CALIFA,
        DET_BMON,
        DET_FI_FIRST,
        DET_FI1A = DET_FI_FIRST,
        DET_FI1B,
        DET_FI2A,
        DET_FI2B,
        DET_FI3A,
        DET_FI3B,
        DET_FI4,
        DET_FI5,
        DET_FI6,
        DET_FI7,
        DET_FI8,
        DET_FI9,
        DET_FI10,
        DET_FI11,
        DET_FI12,
        DET_FI13,
        DET_FI_LAST = DET_FI13,
        DET_L3T,
        DET_LOS,
        DET_MUSIC,
        DET_NEULAND,
        DET_PSPX,
        DET_PTOF,
        DET_ROLU,
        DET_SCI8,
        DET_STRAW,
        DET_TOFD,
        DET_MAX
    };

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    const char* fDetectorNames[DET_MAX + 1] = { "Ams",  "Califa", "BeamMonitor", "Fi1a",  "Fi1b",    "Fi2a",
                                                "Fi2b", "Fi3a",   "Fi3b",        "Fi4",   "Fi5",     "Fi6",
                                                "Fi7",  "Fi8",    "Fi9",         "Fi10",  "Fi11",    "Fi12",
                                                "Fi13", "L3t",    "Los",         "Music", "Neuland", "Pspx",
                                                "Ptof", "Rolu",   "Sci8",        "Straw", "Tofd",    NULL };

    // If FiberI is present or not:
    Int_t ifibdet;
    // Number of fibers per detector
    Double_t n_fiber[NOF_FIB_DET] = { 256.,  256.,  256., 256,  512., 512.,  2048.,
                                      2048., 1024., 512., 512., 512., 1024., 1024. };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    Double_t fClockFreq; /**< Clock cycle in [ns]. */
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */

    //   TClonesArray *fbmonMappedItems;

    Double_t flosVeffXV;
    Double_t flosVeffYV;
    Double_t flosOffsetXV;
    Double_t flosOffsetYV;
    Double_t flosVeffXT;
    Double_t flosVeffYT;
    Double_t flosOffsetXT;
    Double_t flosOffsetYT;
    Double_t flosVeffXQ;
    Double_t flosVeffYQ;
    Double_t flosOffsetXQ;
    Double_t flosOffsetYQ;

    unsigned long long time_V_mem = 0, time_start = 0, time = 0, time_mem = 0;
    unsigned long long time_prev_read = 0, time_to_read = 0;
    unsigned long ic_mem = 0, see_mem = 0, tofdor_mem = 0;
    unsigned long ic_start = 0, see_start = 0, tofdor_start = 0;
    unsigned long long time_spill_start = 0, time_spill_end = 0;

    unsigned long long time_previous_event = 0;

    Double_t time_clear = -1.;
    Double_t tdiff = 0.;
    Double_t fNorm = 1.;
    Int_t iclear_count = 1;
    UInt_t reset_time;         // time after which bmon spectra are reseted
    Double_t read_time;        // step in which scalers are read, in sec
    Int_t fsens_SEE, fsens_IC; // SEETRAM and IC sensitivity, between -4 and -10
    Double_t calib_SEE = 1.;   // SEETRAM calibration factor
    Double_t see_offset = 7.1; // SEETRAM offset in kHz

    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    TH1F* fhTpat;
    TH1F* fh_spill_length;
    TH1F* fhTrigger;
    TH1F* fh_IC_spill;
    TH1F* fh_SEE_spill;
    TH1F* fh_TOFDOR_spill;
    TH1F* fh_IC;
    TH1F* fh_SEE;
    TH1F* fh_TOFDOR;
    TH2F* fh_IC_cal;
    TH2F* fh_SEE_cal;

    TH2F* fh_rolu_tot;
    TH1F* fh_rolu_channels;

    TH1F* fh_sci8_channels;
    TH1F* fh_sci8_tres_MCFD;
    TH1F* fh_sci8_tres_TAMEX;
    TH1F* fh_sci8_tot_mean;
    TH1F* fh_tof_sci8;
    TH2F* fh_sci8_tot;
    TH1F* fh_sci8_dt_hits;
    TH1F* fh_sci8_dt_hits_l;
    TH1F* fh_sci8_dt_hits_t;
    TH1F* fh_sci8_multihit;
    TH2F* fh_sci8_multihitVFTX;
    TH2F* fh_sci8_multihitLEAD;
    TH2F* fh_sci8_multihitTRAI;

    TH1F* fh_los_channels;
    TH1F* fh_los_tres_MCFD;
    TH1F* fh_los_tres_TAMEX;
    TH1F* fh_los_tres_MTDC;
    TH1F* fh_los_tot_mean;
    TH2F* fh_los_tot;
    TH1F* fh_los_dt_hits;
    TH1F* fh_los_multihit;
    TH2F* fh_los_multihitVFTX;
    TH2F* fh_los_multihitLEAD;
    TH2F* fh_los_multihitTRAI;
    TH2F* fh_los_pos_MCFD;
    TH2F* fh_los_pos_TAMEX;
    TH2F* fh_los_pos_ToT;
    TH2F* fh_los_ihit_ToT;
    TH2F* fh_los_dt_hits_ToT;
    TH2F* fh_los_dt_first_ToT;
    TH2F* fh_los_dt_hits_VT;

    TH1F* fh_channels_Fib[NOF_FIB_DET];
    TH1F* fh_channels_single_Fib[NOF_FIB_DET];
    TH1F* fh_fibers_Fib[NOF_FIB_DET];
    TH1F* fh_mult_Fib[NOF_FIB_DET];
    TH2F* fh_Fib_ToF[NOF_FIB_DET];
    TH1F* fh_Fib_pos[NOF_FIB_DET];
    TH2F* fh_time_Fib[NOF_FIB_DET];
    TH2F* fh_multihit_m_Fib[NOF_FIB_DET];
    TH2F* fh_multihit_s_Fib[NOF_FIB_DET];
    TH2F* fh_ToT_m_Fib[NOF_FIB_DET];
    TH2F* fh_ToT_s_Fib[NOF_FIB_DET];
    TH2F* fh_Fib_vs_Events[NOF_FIB_DET];
    TH2F* fh_ToT_single_Fib[NOF_FIB_DET];

    TH1F* fh_tofd_channels[N_PLANE_MAX_TOFD];
    TH2F* fh_tofd_multihit[N_PLANE_MAX_TOFD];
    TH2F* fh_tofd_ToF[N_PLANE_MAX_TOFD];
    TH2F* fh_tofd_TotPm[N_PLANE_MAX_TOFD];
    TH2F* fh_tofd_dt[3];
    TH1F* fh_TimePreviousEvent;

    TH1F* fh_ptof_channels;
    TH1F* fh_ptof_channels_cut;
    TH1F* fh_ptof_test1;
    TH1F* fh_ptof_test2;
    TH1F* fh_ptof_TotPm1[N_PADDLE_MAX_PTOF];
    TH1F* fh_ptof_TotPm2[N_PADDLE_MAX_PTOF];

    TH1F* fh_pspx_channel_x[N_PSPX];      /**< PSPX x channel profile on mapped level */
    TH1F* fh_pspx_channel_y[N_PSPX];      /**< PSPX y channel profile on mapped level */
    TH1F* fh_pspx_multiplicity_x[N_PSPX]; /**< PSPX x multiplicity on mapped level */
    TH1F* fh_pspx_multiplicity_y[N_PSPX]; /**< PSPX y multiplicity on mapped level */
    TH1F* fh_pspx_mapped_energy[N_PSPX];  /**< PSPX energy on mapped level as first check */

    TH2F* fh_pspx_strips_position[N_PSPX];  /**< PSPX 2D position on mapped level */
    TH2F* fh_pspx_hit_position[N_PSPX];     /**< PSPX 2D position on hit level */
    TH2F* fh_pspx12_strips_position_corr_x; // correlation between x1 vs x2
    TH2F* fh_pspx12_strips_position_corr_y; // and y1 vs y2, ONLY valid for 2 detectors at the moment
    TH2F* fh_pspx12_cal_energy_corr;        // correlation between en1 und en2

    TH1F* fh_pspx_hit_energy[N_PSPX]; /**< PSPX energy on hit level */

    TH2F* fh_pspx_cal_energy_frontback[N_PSPX]; /**< PSPX energy front vs back on cal level */
    TH2F* fh_pspx_hit_multi[(N_PSPX + 1) / 2];

  public:
    ClassDef(R3BOnlineSpectra, 2)
};

#endif
