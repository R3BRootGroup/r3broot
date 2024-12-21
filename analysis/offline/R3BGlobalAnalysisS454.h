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
// -----                  R3BGlobalAnalysisS454                -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BGLOBALANALYSISS454
#define R3BGLOBALANALYSISS454

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

#include "TClonesArray.h"
#include "TMath.h"
#include <TLorentzVector.h>
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BGlobalAnalysisS454 : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BGlobalAnalysisS454();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BGlobalAnalysisS454(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BGlobalAnalysisS454();

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
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Methods for setting reset and readout times for Bmon
     */
    inline void SetBmon(Int_t sens_SEE, Int_t sens_IC)
    {
        fsens_SEE = sens_SEE;
        fsens_IC = sens_IC;
    }

    /**
     * Methods for setting cuts
     */
    inline void SetCuts(Bool_t cuts) { fCuts = cuts; }

    inline void SetGhost(Bool_t ghost) { fGhost = ghost; }

    inline void SetPairs(Bool_t p) { fPairs = p; }

    inline void SetBfield(Int_t B) { fB = B; }

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;
    TClonesArray* fMCTrack;
    TClonesArray* fTrack;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    Bool_t fCuts;
    Bool_t fGhost;
    Bool_t fPairs;
    Int_t fB;
    Bool_t tracker = true;

    unsigned long long time_start = 0, time = 0;
    unsigned long ic_start = 0, see_start = 0, tofdor_start = 0;
    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    Int_t maxevent;

    Int_t fsens_SEE, fsens_IC; // SEETRAM and IC sensitivity, between -4 and -10
    Double_t calib_SEE = 1.;   // SEETRAM calibration factor
    Double_t see_offset = 7.1; // SEETRAM offset in kHz
    Double_t counts_SEE = 0;
    Double_t counts_IC = 0;
    Double_t counts_TofD = 0;

    Double_t XHe, YHe, ZHe, XC, YC, ZC, THe, TC;
    Double_t pHex, pHey, pHez, pCx, pCy, pCz;
    Double_t Pxf, Pyf, Pzf, Xf, Yf, Zf, Pf_tot;

    Double_t amu = 931.49410242;
    //	Double_t mHe = 4.00260325413*amu;
    //	Double_t mC = 12. * amu;
    Double_t mHe = 3727.409;
    Double_t mC = 11174.950;
    Double_t mO = 15.99065084 * amu;

    TLorentzVector alpha, carbon, oxygen;
    TLorentzVector alpha_cm, carbon_cm, oxygen_cm;
    TLorentzVector photon;

    Int_t Q = 0;
    Double_t tPrev[10];
    Int_t detPrev[10];

    Int_t counter1 = 0;
    Int_t counter2 = 0;
    Int_t counter3 = 0;
    Int_t counter4 = 0;
    Int_t countdet;

    UInt_t num_spills = 0;

    TH1F* fh_Tpat;
    TH1F* fh_Trigger;
    TH1F* fh_IC;
    TH1F* fh_SEE;
    TH1F* fh_TOFDOR;
    TH2F* fh_Cave_position;
    TH2F* fh_target_xy;
    TH1F* fh_target_px;
    TH1F* fh_target_py;
    TH1F* fh_target_pz;
    TH1F* fh_target_p;
    TH1F* fh_chi2;
    TH1F* fh_px_He;
    TH1F* fh_py_He;
    TH1F* fh_pz_He;
    TH1F* fh_p_He;
    TH1F* fh_px_C;
    TH1F* fh_py_C;
    TH1F* fh_pz_C;
    TH1F* fh_p_C;

    TH2F* fh_chiy_vs_chix_He;
    TH2F* fh_chiy_vs_chix_C;

    TH1F* fh_dx;
    TH1F* fh_dy;
    TH1F* fh_dz;
    TH1F* fh_dpx;
    TH1F* fh_dpy;
    TH1F* fh_dpz;
    TH1F* fh_dp;
    TH2F* fh_thetax_dpx;
    TH2F* fh_thetax_dpx_abs;
    TH2F* fh_thetay_dpy;
    TH2F* fh_x_dpx;
    TH2F* fh_y_dpy;
    TH2F* fh_thetax_dpy;
    TH2F* fh_thetay_dpx;
    TH2F* fh_dpy_dpx;
    TH2F* fh_px_px_cm;
    TH2F* fh_py_py_cm;
    TH2F* fh_pz_pz_cm;
    TH2F* fh_px_px;
    TH2F* fh_py_py;
    TH2F* fh_pz_pz;

    TH1F* fh_theta26_simu;
    TH1F* fh_Erel_simu;
    TH1F* fh_theta26;
    TH1F* fh_theta_16O;
    TH1F* fh_phi_16O;
    TH1F* fh_theta26_cm;
    TH1F* fh_phi26_cm;
    TH1F* fh_theta_4He_cm;
    TH1F* fh_phi_4He_cm;
    TH1F* fh_phi_bc_cm;
    TH1F* fh_theta_bc_cm;
    TH1F* fh_theta_12C_cm;
    TH1F* fh_phi_12C_cm;
    TH1F* fh_Erel;
    TH1F* fh_ErelL;
    TH1F* fh_ErelR;
    TH1F* fh_ErelB;
    TH2F* fh_phi_bc_cm_polar;

    TH2F* fh_dErel_vs_x;
    TH2F* fh_dErel_vs_y;

  public:
    ClassDef(R3BGlobalAnalysisS454, 1)
};

#endif
