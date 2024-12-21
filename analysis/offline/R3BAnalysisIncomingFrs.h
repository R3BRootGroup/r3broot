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

#ifndef R3BAnalysisIncomingFrs_H
#define R3BAnalysisIncomingFrs_H

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "R3BFrsData.h"
#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BAnalysisIncomingFrs : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BAnalysisIncomingFrs();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BAnalysisIncomingFrs(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BAnalysisIncomingFrs();

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
    void SetLosXYMCFD(Double_t offsetX1,
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

    void SetLosXYTAMEX(Double_t offsetXT1,
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

    void SetLosXYToT(Double_t offsetXQ1,
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
    /* Method for pile-up */
    void SetEpileup(Double_t Epileup) { fEpileup = Epileup; }

    void SetZcut(Double_t zcut, Double_t zcutw)
    {
        ZMUSIC_cut = zcut;
        ZMUSIC_wcut = zcutw;
    }
    void SetAcut(Double_t acut, Double_t acutw)
    {
        AoQ_cut = acut;
        AoQ_wcut = acutw;
    }

    /* Method for setting number of LOS detectors */
    void SetNofLosModules(Int_t nDets)
    {
        fNofLosDetectors = nDets;

        time_prev.resize(nDets);
    }

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    void SetTpat(Int_t tpat) { fTpat = tpat; }

    // --- hard coded calibration
    void SetToFoffset(Double_t offset) { fToFoffset = offset; }
    void SetToFmin(Double_t t) { fToFmin = t; }
    void SetToFmax(Double_t t) { fToFmax = t; }
    void SetTof2InvV_p0(Double_t p) { fTof2InvV_p0 = p; }
    void SetTof2InvV_p1(Double_t p) { fTof2InvV_p1 = p; }
    void SetFlightLength(Double_t L) { fL2 = L; }
    void SetPos_p0(Double_t p) { fPos_p0 = p; }
    void SetPos_p1(Double_t p) { fPos_p1 = p; }
    void SetDispersionS2(Double_t DS2) { fDispersionS2 = DS2; }
    void SetBrho0_S2toCC(Double_t Brho0) { fBrho0_S2toCC = Brho0; }
    void SetBetaCorrectionForZ(Double_t p0, Double_t p1, Double_t p2, Double_t Zprimary, Double_t Zoffset)
    {
        fP0 = p0;
        fP1 = p1;
        fP2 = p2;
        fZprimary = Zprimary;
        fZoffset = Zoffset;
    }

    void SetPosS2_limits(Double_t min, Double_t max)
    {
        fPos2min = min;
        fPos2max = max;
    }

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    TClonesArray* fHitItemsMus;
    TClonesArray* fFrsDataCA; /**< Array with FRS-output data. >*/

    Double_t fPos2min, fPos2max;

    TClonesArray* fTcalSci2; /**< Array with Tcal items. */

    enum DetectorInstances
    {
        DET_LOS,
        DET_SAMPLER,
        DET_MAX
    };

    const char* fDetectorNames[DET_MAX + 1] = { "Los", "Sampler", NULL };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    Int_t fSamp;
    Double_t fClockFreq; /**< Clock cycle in [ns]. */
    Int_t nLosEvents = 0;
    // Adjust the offset when a DAQ subsystem restarts
    Double_t fToFoffset;
    Bool_t fOnline; // Don't store data for online

    //   TClonesArray *fbmonMappedItems;
    Int_t fNofLosDetectors; /**< Number of LOS detectors. */

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
    Double_t fP0, fP1, fP2, fZprimary, fZoffset;

    unsigned long long time_V_mem = 0, time_start = 0, time = 0, time_mem = 0;
    std::vector<unsigned long long> time_prev;
    unsigned long long time_prev_read = 0, time_to_read = 0;
    unsigned long ic_mem = 0, see_mem = 0, tofdor_mem = 0;
    unsigned long ic_start = 0, see_start = 0, tofdor_start = 0;
    unsigned long long time_spill_start = 0, time_spill_end = 0;
    long samplerPrev = 0, samplerCurr = 0, samplerSpill = 0;
    unsigned long long time_previous_event = 0;

    Double_t time_clear = -1.;
    Double_t tdiff = 0.;
    Double_t fNorm = 1.;
    Double_t fEpileup;

    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

    // --- 1/ select the good hit from the ToF from S2 to Cave C:
    Double_t fToFmin;
    Double_t fToFmax;

    // --- 2/ calibration in beta and ToF
    Double_t fTof2InvV_p0;
    Double_t fTof2InvV_p1;
    Double_t fL2;

    // --- 3/ calibration in position
    Double_t fPos_p0;
    Double_t fPos_p1;

    // --- 4/ AoQ calculation
    Double_t fDispersionS2;
    Double_t fBrho0_S2toCC;

    Double_t ZMUSIC_cut;
    Double_t ZMUSIC_wcut;

    Double_t AoQ_cut;
    Double_t AoQ_wcut;

    /** Private method FrsData **/
    //** Adds a FrsData to the analysis
    R3BFrsData* AddData(Int_t StaId,
                        Int_t StoId,
                        Double_t z,
                        Double_t aq,
                        Double_t betaval,
                        Double_t brhoval,
                        Double_t xs2,
                        Double_t xc);

  public:
    ClassDef(R3BAnalysisIncomingFrs, 1)
};

#endif
