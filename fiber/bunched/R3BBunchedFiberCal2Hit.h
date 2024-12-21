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

#ifndef R3BBUNCHEDFIBERCAL2HIT_H
#define R3BBUNCHEDFIBERCAL2HIT_H 1

#include "FairTask.h"
#include <TClonesArray.h>

#include <R3BTCalEngine.h>

#include <list>

class TH1F;
class TH2F;

class R3BFiberMappingPar;
class R3BCoarseTimeStitch;
class R3BBunchedFiberCalData;
class R3BBunchedFiberHitPar;
class R3BBunchedFiberHitModulePar;
class R3BEventHeader;

class R3BBunchedFiberCal2Hit : public FairTask
{
  public:
    // This defines the fiber direction
    enum Direction
    {
        HORIZONTAL,
        VERTICAL
    };
    // This tells us where the first fiber is located:
    // VERTICAL direction ---------------------------------
    // STANDARD: fiber 1 wixhausen -> fiber 512 messel
    // INVERTED: fiber 1 messel -> fiber 512 wixhausen
    // HORIZONTAL direction -------------------------------
    // STANDARD: fiber 1 down -> fiber 512 up
    // INVERTED: fiber 1 up -> fiber 512 down
    enum Orientation
    {
        STANDARD,
        INVERTED
    };
    struct ToT
    {
        ToT(R3BBunchedFiberCalData const*, R3BBunchedFiberCalData const*, Double_t, Double_t, Double_t);
        R3BBunchedFiberCalData const* lead;
        R3BBunchedFiberCalData const* trail;
        Double_t lead_ns, tail_ns, tot_ns;
    };
    struct Channel
    {
        std::list<R3BBunchedFiberCalData const*> lead_list;
        std::list<ToT> tot_list;
    };

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BBunchedFiberCal2Hit(const char*,
                           Int_t,
                           enum R3BTCalEngine::CTDCVariant,
                           Direction,
                           UInt_t,
                           UInt_t,
                           UInt_t,
                           Bool_t,  // is calib
                           Bool_t,  // is gain
                           Bool_t); // is tsync

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BBunchedFiberCal2Hit();

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

    virtual void SetParContainers();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t*);

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

    void SetOrientation(Orientation opt) { fOrientation = opt; }

    void SetExpId(Int_t opt)
    {
        fExpId = opt;
    } // Mutator to set fExpId manually. It should be globally defined by EventHeader.

  private:
    void Standard();
    void S515();

    R3BEventHeader* header; /**< Event header. */

    TString fName;
    Int_t fnEvents;
    Int_t maxevent;
    Int_t fnEventsfill = 0;
    Int_t fNumFibers;
    Int_t fExpId;
    int fDetId;

    Int_t multi = 0;
    Double_t energy[2048];
    Int_t counts[2048];
    Double_t tsync_temp[2048];
    Double_t gain_temp[2048];
    Bool_t tofdin;

    double fClockFreq;
    Direction fDirection;
    Orientation fOrientation;
    R3BCoarseTimeStitch* fTimeStitch;
    UInt_t fSubNum;
    UInt_t fChPerSub[2];
    Bool_t fIsCalibrator;
    Bool_t fIsGain;
    Bool_t fIsTsync;
    TClonesArray* fCalItems;
    TClonesArray* fMAPMTCalTriggerItems;
    TClonesArray* fSPMTCalTriggerItems;
    TClonesArray* fHitItems;
    TClonesArray* fTofdHitItems;
    unsigned const* fMAPMTTriggerMap;
    unsigned const* fSPMTTriggerMap;
    R3BFiberMappingPar* fMapPar;
    R3BBunchedFiberHitPar* fCalPar; /**< Parameter container. */
    R3BBunchedFiberHitPar* fHitPar; /**< Hit parameter container. */
    Int_t fNofHitPars;              /**< Number of modules in parameter file. */
    Int_t fNofHitItems;
    // [0=MAPMT,1=SPMT][Channel].
    std::vector<Channel> fChannelArray[2];

    // histograms for gain matching
    TH2F* fh_ToT_MA_Fib;
    TH2F* fh_ToT_Single_Fib;
    TH2F* fh_ToT_s_Fib[4];
    TH2F* fh_time_s_Fib;
    TH2F* fh_ToT_ToT;
    TH2F* fh_dt_Fib;
    TH2F* fh_Fib_ToF;
    TH2F* fh_Test;
    TH1F* fh_multi;
    TH2F* fh_time_Fib;

  public:
    ClassDef(R3BBunchedFiberCal2Hit, 3)
};

#endif
