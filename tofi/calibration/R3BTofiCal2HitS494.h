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
// -----                 R3BTofiCal2HitS494                   -----
// -----         Created May 30th 2016 by M.Heil          -----
// -----           Modified Mar 2020 by L.Bott            -----
// ------------------------------------------------------------

#ifndef R3BTOFICAL2HITS494
#define R3BTOFICAL2HITS494

#define N_TOFI_HIT_PLANE_MAX 1
#define N_TOFI_HIT_PADDLE_MAX 24

#include <map>

#include "FairTask.h"
#include "THnSparse.h"

class TClonesArray;
class R3BTofiHitModulePar;
class R3BTofiHitPar;
class R3BEventHeader;
class TH1F;
class TH2F;

/**
 * An analysis task to apply HIT calibration for ToFi.
 * it also applies walk corrections and Z-calibration
 */
class R3BTofiCal2HitS494 : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTofiCal2HitS494();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTofiCal2HitS494(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofiCal2HitS494();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    virtual InitStatus ReInit();

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

    virtual void SetParContainers();

    virtual void CreateHistograms(Int_t iPlane, Int_t iBar);

    /**
     * Method for setting the nuclear charge of main beam
     */
    inline void SetTofiQ(Double_t Q) { fTofiQ = Q; }

    /**
     * Method for setting histograms
     */
    inline void SetTofiHisto(Bool_t Histo) { fTofiHisto = Histo; }
    /**
     * Method for setting y calculation via ToT instead of tdiff
     */
    inline void SetTofiTotPos(Bool_t ToTy) { fTofiTotPos = ToTy; }
    /**
     * Method for setting gap between left and right half
     */
    inline void SetTofiGap(Float_t Gap) { fTofiGap = Gap; }

    inline void SetShowProgress(Float_t progress) { fShowProgress = progress; }

    inline void SetCalOpt(Bool_t cal) { cal_hand = cal; }
    /**
     * Old Method for walk calculation.
     */
    /// virtual Double_t walk(Double_t Q);
    /**
     * new Method for walk calculation.
     */
    virtual Double_t walk(Double_t Q, Double_t par1, Double_t par2, Double_t par3, Double_t par4, Double_t par5);
    /**
     * Method for beta correction.
     */
    virtual Double_t betaCorr(Double_t delta);
    /**
     * Method for calculation of saturation.
     */
    virtual Double_t saturation(Double_t x);
    /**

     * Method for selecting events with certain trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }
    /**
     * Methods for setting number of planes and paddles
     */
    inline void SetNofModules(Int_t planes, Int_t ppp)
    {
        fNofPlanes = planes;
        fPaddlesPerPlane = ppp;
    }

  private:
    TClonesArray* fCalItems;        /**< Array with Cal items - input data. */
    TClonesArray* fCalTriggerItems; /**< Array with trigger Cal items - input data. */
    TClonesArray* fHitItems;        /**< Array with Hit items - output data. */
    UInt_t fNofHitItems;            /**< Number of hit items for cur event. */
    R3BTofiHitPar* fHitPar;         /**< Hit parameter container. */
    UInt_t fNofHitPars;             /**< Number of modules in parameter file. */
    R3BEventHeader* header;         /**< Event header - input data. */
    Double_t fClockFreq;            /**< Clock cycle in [ns]. */
    Int_t fTrigger;                 /**< Trigger value. */
    Int_t fTpat1, fTpat2;
    Double_t fTofiQ;
    Bool_t cal_hand;
    Float_t fTofiGap;
    Bool_t fTofiHisto;
    Bool_t fTofiTotPos;
    Bool_t fShowProgress;
    UInt_t fnEvents;
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */
    UInt_t maxevent;
    UInt_t wrongtrigger;
    UInt_t wrongtpat;
    UInt_t headertpat;
    UInt_t events_in_cal_level;
    UInt_t inbarcoincidence;
    UInt_t eventstore;
    UInt_t singlehit;
    UInt_t multihit;
    UInt_t bars_with_multihit;
    UInt_t events_wo_tofi_hits;

    // arrays of control histograms
    TH2F* fhQvsPos[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];
    // TH2F* fhQvsTHit[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];
    // TH2F* fhTvsTHit[N_TOFI_HIT_PLANE_MAX][N_TOFI_HIT_PADDLE_MAX];
    TH2F* fhQ[N_TOFI_HIT_PLANE_MAX];
    TH2F* fhxy[N_TOFI_HIT_PLANE_MAX];
    TH2F* fhQvsEvent[N_TOFI_HIT_PLANE_MAX];
    TH2F* fhTdiff[N_TOFI_HIT_PLANE_MAX];
    TH2F* fhTsync[N_TOFI_HIT_PLANE_MAX];
    TH2F* fh_Tofi_TotPm[N_TOFI_HIT_PLANE_MAX];

  public:
    ClassDef(R3BTofiCal2HitS494, 1)
};

#endif
