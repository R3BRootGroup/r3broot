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
// -----                 R3BTofdCal2HitS454                   -----
// -----         Created May 30th 2016 by M.Heil          -----
// -----           Modified Mar 2020 by L.Bott            -----
// ------------------------------------------------------------

#ifndef R3BTOFDCAL2HITS454
#define R3BTOFDCAL2HITS454

#define N_TOFD_HIT_PLANE_MAX 4
#define N_TOFD_HIT_PADDLE_MAX 44

#include <map>

#include "FairTask.h"
#include "THnSparse.h"

class TClonesArray;
class R3BTofdHitModulePar;
class R3BTofdHitPar;
class R3BEventHeader;
class TH1F;
class TH2F;

/**
 * An analysis task to apply HIT calibration for ToFd.
 * it also applies walk corrections and Z-calibration
 */
class R3BTofdCal2HitS454 : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTofdCal2HitS454();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTofdCal2HitS454(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofdCal2HitS454();

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
    inline void SetTofdQ(Double_t Q) { fTofdQ = Q; }

    /**
     * Method for setting histograms
     */
    inline void SetTofdHisto(Bool_t Histo) { fTofdHisto = Histo; }
    /**
     * Method for setting y calculation via ToT instead of tdiff
     */
    inline void SetTofdTotPos(Bool_t ToTy) { fTofdTotPos = ToTy; }

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
     * Method for insert.
     */
    virtual Double_t* insertX(Int_t n, Double_t arr[], Double_t x, Int_t pos);

    /**
     * Method for selecting events with certain trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }
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
    R3BTofdHitPar* fHitPar;         /**< Hit parameter container. */
    UInt_t fNofHitPars;             /**< Number of modules in parameter file. */
    R3BEventHeader* header;         /**< Event header - input data. */
    Double_t fClockFreq;            /**< Clock cycle in [ns]. */
    Int_t fTrigger;                 /**< Trigger value. */
    Int_t fTpat;
    Double_t fTofdQ;
    Bool_t fTofdHisto;
    Bool_t fTofdTotPos;
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
    UInt_t events_wo_tofd_hits;

    // arrays of control histograms
    TH2F* fhQvsPos[N_TOFD_HIT_PLANE_MAX][N_TOFD_HIT_PADDLE_MAX];
    // TH2F* fhQvsTHit[N_TOFD_HIT_PLANE_MAX][N_TOFD_HIT_PADDLE_MAX];
    // TH2F* fhTvsTHit[N_TOFD_HIT_PLANE_MAX][N_TOFD_HIT_PADDLE_MAX];
    TH2F* fhQ[N_TOFD_HIT_PLANE_MAX];
    TH2F* fhxy[N_TOFD_HIT_PLANE_MAX];
    TH2F* fhQvsEvent[N_TOFD_HIT_PLANE_MAX];
    TH2F* fhTdiff[N_TOFD_HIT_PLANE_MAX];
    TH2F* fhTsync[N_TOFD_HIT_PLANE_MAX];

  public:
    ClassDef(R3BTofdCal2HitS454, 1)
};

#endif
