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
// -----              R3BOnlineSpectraDec2019             -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BONLINESPECTRADEC2019
#define R3BONLINESPECTRADEC2019

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
class R3BOnlineSpectraDec2019 : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraDec2019();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraDec2019(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraDec2019();

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
    inline void SetLosXYMCFD(Double_t offsetX1,
                             Double_t offsetY1,
                             Double_t veffX1,
                             Double_t veffY1,
                             Double_t offsetX2,
                             Double_t offsetY2,
                             Double_t veffX2,
                             Double_t veffY2)
    {
        flosOffsetXV[0] = offsetX1;
        flosOffsetYV[0] = offsetY1;
        flosVeffXV[0] = veffX1;
        flosVeffYV[0] = veffY1;
        flosOffsetXV[1] = offsetX2;
        flosOffsetYV[1] = offsetY2;
        flosVeffXV[1] = veffX2;
        flosVeffYV[1] = veffY2;
    }

    inline void SetLosXYTAMEX(Double_t offsetXT1,
                              Double_t offsetYT1,
                              Double_t veffXT1,
                              Double_t veffYT1,
                              Double_t offsetXT2,
                              Double_t offsetYT2,
                              Double_t veffXT2,
                              Double_t veffYT2)
    {
        flosOffsetXT[0] = offsetXT1;
        flosOffsetYT[0] = offsetYT1;
        flosVeffXT[0] = veffXT1;
        flosVeffYT[0] = veffYT1;
        flosOffsetXT[1] = offsetXT2;
        flosOffsetYT[1] = offsetYT2;
        flosVeffXT[1] = veffXT2;
        flosVeffYT[1] = veffYT2;
    }

    inline void SetLosXYToT(Double_t offsetXQ1,
                            Double_t offsetYQ1,
                            Double_t veffXQ1,
                            Double_t veffYQ1,
                            Double_t offsetXQ2,
                            Double_t offsetYQ2,
                            Double_t veffXQ2,
                            Double_t veffYQ2)
    {
        flosOffsetXQ[0] = offsetXQ1;
        flosOffsetYQ[0] = offsetYQ1;
        flosVeffXQ[0] = veffXQ1;
        flosVeffYQ[0] = veffYQ1;
        flosOffsetXQ[1] = offsetXQ2;
        flosOffsetYQ[1] = offsetYQ2;
        flosVeffXQ[1] = veffXQ2;
        flosVeffYQ[1] = veffYQ2;
    }
    /* Method for setting number of LOS detectors */
    inline void SetNofLosModules(Int_t nDets) { fNofLosDetectors = nDets; }
    /**
     * Methods for setting condition on filling xs from LOS
     * 0 = no condition, 1 = condition from fiber
     */
    inline void SetLosPosition(Int_t optcond) { foptcond = optcond; }
    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Methods for setting reset and readout times for Bmon
     */
    inline void SetBmon(Int_t time_range, Double_t time_step, Double_t SpillLength, Double_t DAQ_dead_time)
    {
        reset_time = time_range;
        read_time = time_step;
        fDAQ_dead_time = DAQ_dead_time; // dead time in micro seconds
        fSpillLength = SpillLength;
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
        DET_SAMPLER,
        DET_MAX
    };

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    const char* fDetectorNames[DET_MAX + 1] = { "Ams",   "Califa",  "BeamMonitor", "Fi1a", "Fi1b", "Fi2a", "Fi2b",
                                                "Fi3a",  "Fi3b",    "Fi4",         "Fi5",  "Fi6",  "Fi7",  "Fi8",
                                                "Fi9",   "Fi10",    "Fi11",        "Fi12", "Fi13", "L3t",  "Los",
                                                "Music", "Neuland", "Pspx",        "Ptof", "Rolu", "Sci8", "Straw",
                                                "Tofd",  "Sampler", NULL };

    // If FiberI is present or not:
    Int_t ifibdet;
    // Number of fibers per detector
    Double_t n_fiber[NOF_FIB_DET] = { 256.,  256.,  256., 256,  512., 512.,  2048.,
                                      2048., 1024., 512., 512., 512., 1024., 1024. };
    Double_t tof_fiber[NOF_FIB_DET] = { 2580., 2580., 256., 256,  512., 512.,  2048.,
                                        2048., 1024., 512., 512., 512., 1024., 1024. };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    Int_t fSamp;
    Double_t fClockFreq; /**< Clock cycle in [ns]. */
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */

    //   TClonesArray *fbmonMappedItems;

    Int_t fNofLosDetectors; /**< Number of LOS detectors. */

#define NOF_LOS_DET 2

    Double_t flosVeffXV[2];
    Double_t flosVeffYV[2];
    Double_t flosOffsetXV[2];
    Double_t flosOffsetYV[2];
    Double_t flosVeffXT[2];
    Double_t flosVeffYT[2];
    Double_t flosOffsetXT[2];
    Double_t flosOffsetYT[2];
    Double_t flosVeffXQ[2];
    Double_t flosVeffYQ[2];
    Double_t flosOffsetXQ[2];
    Double_t flosOffsetYQ[2];
    Int_t foptcond;

    unsigned long long time_V_mem = 0, time_start = 0, time_mem = 0, time = 0, time_begin = 0;
    unsigned long long time_prev[NOF_LOS_DET] = { 0 };

    unsigned long long time_prev_read = 0, time_to_read = 0;
    unsigned long ic_mem = 0, see_mem = 0, tofdor_mem = 0;
    unsigned long ic_start = 0, see_start = 0, tofdor_start = 0;
    unsigned long long time_spill_start = 0, time_spill_end = 0;
    unsigned long long time_previous_event = 0;
    long samplerPrev = 0, samplerCurr = 0, samplerSpill = 0, samplerBegin = 0;
    UInt_t spillCounter = 0;

    Bool_t spill_on = false;
    Bool_t spill_off = false;
    Bool_t spill_on_sampler = false;
    Double_t time_clear = -1.;
    Double_t tdiff = 0.;
    Double_t fNorm = 1.;
    Int_t iclear_count = 1;
    UInt_t reset_time;         // time after which bmon spectra are reseted
    Double_t read_time;        // step in which scalers are read, in sec
    Double_t calib_SEE = 1.;   // SEETRAM calibration factor
    Double_t see_offset = 7.1; // SEETRAM offset in kHz
    Double_t fSpillLength;
    Double_t fDAQ_dead_time;
    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    TH1F* fhTpat;
    TH1F* fh_spill_length;
    TH1F* fhTrigger;
    TH1F* fh_particle;
    TH1F* fh_particle_pois;
    TH1F* fh_TOFDOR_spill;
    TH1F* fh_spill_times;
    TH1F* fh_spill_times_zoom;
    TH1F* fh_spill_times_pois;
    TH1F* fh_spill_times_daq;
    TH1F* fh_spill_times_daq_first;
    TH1F* fh_frequencies;
    TH1F* fh_TOFDOR;
    TH1F* fh_rate;
    // TH1F *h3;
    TH1F* fh_TimePreviousEvent;

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

    TH2F* fh_rolu_tot;
    TH1F* fh_rolu_channels;

    TH1F* fh_los_channels[NOF_LOS_DET];
    TH1F* fh_los_tres_MCFD[NOF_LOS_DET];
    TH1F* fh_los_tres_TAMEX[NOF_LOS_DET];
    TH1F* fh_los_tres_MTDC;
    TH1F* fh_los_tot_mean[NOF_LOS_DET];
    TH2F* fh_los_tot[NOF_LOS_DET];
    TH1F* fh_los_multihit[NOF_LOS_DET];
    TH2F* fh_los_multihitVFTX[NOF_LOS_DET];
    TH2F* fh_los_multihitLEAD[NOF_LOS_DET];
    TH2F* fh_los_multihitTRAI[NOF_LOS_DET];
    TH2F* fh_los_pos_MCFD[NOF_LOS_DET];
    TH2F* fh_los_pos_TAMEX[NOF_LOS_DET];
    TH2F* fh_los_pos_ToT[NOF_LOS_DET];
    TH2F* fh_los_ihit_ToT[NOF_LOS_DET];
    TH2F* fh_los_dt_hits_ToT[NOF_LOS_DET];
    TH2F* fh_los_ToT_vs_ToT;

    TH1F* fh_los_TOF[NOF_LOS_DET];
    TH1F* fh_los_dt_hits[NOF_LOS_DET];
    TH2F* fh_los_dt_first_ToT[NOF_LOS_DET];
    TH1F* fh_dt_hits[NOF_LOS_DET];

    TH1F* fh_channels_Fib[NOF_FIB_DET];
    TH1F* fh_channels_single_Fib[NOF_FIB_DET];
    TH1F* fh_fibers_Fib[NOF_FIB_DET];
    TH1F* fh_mult_Fib[NOF_FIB_DET];
    TH2F* fh_Fib_ToF[NOF_FIB_DET];
    TH2F* fh_Fib_pos_xy[NOF_FIB_DET];
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

    TH2F* fh_pspx_strips_position[N_PSPX];        /**< PSPX 2D position on mapped level */
    TH2F* fh_pspx_hit_position[(N_PSPX + 1) / 2]; /**< PSPX 2D position on hit level */

    TH1F* fh_pspx_hit_energy[(N_PSPX + 1) / 2]; /**< PSPX energy on hit level */
    TH2F* fh_pspx_cal_energy_frontback[N_PSPX]; /**< PSPX energy front vs back on cal level */

  public:
    ClassDef(R3BOnlineSpectraDec2019, 2)
};

#endif
