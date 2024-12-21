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
// ----------------------------------------------------------
// ----- Create histograms for parameters for Tofi      -----
// -----     Created July 2019 by L.Bott                -----
// ----------------------------------------------------------

#ifndef R3BTOFICAL2HISTO_H
#define R3BTOFICAL2HISTO_H

#define N_TOFI_HIT_PLANE_MAX 1
#define N_TOFI_HIT_PADDLE_MAX 24

#include <map>

#include "FairTask.h"
#include "R3BTCalEngine.h"
#include "TObject.h"

class R3BTofiHitPar;
class TClonesArray;
class R3BEventHeader;
class TH1F;
class TH2F;
#ifdef __CINT__
#pragma link C++ class R3BTofiHitModulePar + ;
#endif

class R3BTofiCal2Histo : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTofiCal2Histo();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTofiCal2Histo(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofiCal2Histo();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    virtual void SetParContainers();

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
     * Method for setting the update rate for control histograms
     * @param rate an update rate value (events).
     */
    inline void SetUpdateRate(Int_t rate) { fUpdateRate = rate; }

    /**
     * Method for setting minimum required statistics per module.
     * Only detector modules with number of entries in TDC
     * distribution greater than minimum statistics will be
     * calibrated.
     * @param minStats a value of minimum statistics required.
     */
    inline void SetMinStats(Int_t minStats) { fMinStats = minStats; }

    /**
     * Method for selecting events with certain trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Method for setting number of LOS detectors and channels.
     * @param nDets number of detectors.
     * @param nCh number of channels per detector (4+master trigger?)
     */
    inline void SetNofModules(Int_t nDets, Int_t nCh)
    {
        fNofPlanes = nDets;
        fPaddlesPerPlane = nCh;
        fNofModules = nDets * nCh * 1; // 1 edge per channel
    }
    /**
     * Method for setting the y-position of a horizonzal sweep run for calibration of effective velocity of light
     */
    inline void SetTofiPos(Double_t Y) { fTofiY = Y; }
    /**
     * Method for setting the nuclear charge of main beam
     */
    inline void SetTofiQ(Double_t Q) { fTofiQ = Q; }
    /**
     * Method for setting charge correction
     */
    inline void SetTofiZ(Bool_t Z) { fTofiZ = Z; }
    /**
     * Method for using smiley or double exponential charge correction
     */
    inline void SetTofiSmiley(Bool_t Smiley) { fTofiSmiley = Smiley; }
    /**
     *
     */
    inline void ReadParaFile(TString file) { fParaFile = file; }
    /**
     * old Method for walk calculation.
     */
    /// virtual Double_t walk(Double_t Q);
    /**
     * new Method for walk calculation.
     */
    virtual Double_t walk(Double_t Q, Double_t par1, Double_t par2, Double_t par3, Double_t par4, Double_t par5);
    /**
     * Method for creating histograms.
     */
    virtual void CreateHistograms(Int_t iPlane, Int_t iBar);

    /**
     * Method for calculation of saturation.
     */
    virtual Double_t saturation(Double_t x);

  private:
    Int_t fUpdateRate; /**< An update rate. */
    Int_t fMinStats;   /**< Minimum statistics required per module. */
    Int_t fTrigger;    /**< Trigger value. */
    Int_t fTpat;

    UInt_t fNofPlanes;       /**< Number of planes. */
    UInt_t fPaddlesPerPlane; /**< Number of bars per plane. */
    UInt_t fNofModules;      /**< Total number of modules (=edges) to calibrate */

    UInt_t fNEvents;                /**< Event counter. */
    R3BTofiHitPar* fCal_Par;        /**< Parameter container. */
    TClonesArray* fCalData;         /**< Array with mapped data - input data. */
    TClonesArray* fCalTriggerItems; /**< Array with trigger Cal items - input data. */
    R3BEventHeader* header;         /**< Event header - input data. */
    Double_t fClockFreq;            /**< Clock cycle in [ns]. */
    Double_t fTofiY;
    Double_t fTofiQ;
    Bool_t fTofiSmiley;
    Bool_t fTofiZ;
    TString fParaFile;
    UInt_t maxevent;

    // arrays of control histograms
    TH2F* fh_Tofi_TotPm[N_TOFI_HIT_PLANE_MAX];
    TH2F* fhTdiff[N_TOFI_HIT_PLANE_MAX];
    TH2F* fhTsync[N_TOFI_HIT_PLANE_MAX];
    TH2F* fhLogTot1vsLogTot2[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];
    TH2F* fhSqrtQvsPosToT[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];
    TH2F* fhQvsPos[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];
    TH2F* fhTot1vsPos[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];
    TH2F* fhTot2vsPos[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];

  public:
    ClassDef(R3BTofiCal2Histo, 1)
};

#endif
