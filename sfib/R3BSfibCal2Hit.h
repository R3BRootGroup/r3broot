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

#ifndef R3BSFIBCAL2HIT
#define R3BSFIBCAL2HIT

#include "FairTask.h"
#include <TClonesArray.h>

#include <R3BTCalEngine.h>

#include <list>

class TH1F;
class TH2F;

class R3BSfibCalData;
class R3BSfibHitPar;

/**
 * Transforms bunched fiber Cal level data to Hit level.
 *
 * A note about the UCESB unpacker:
 * All channels in a full plane of fibers should be mapped in sequence. For
 * example, Fib4 is built out of two "sub-detectors" (left and right side),
 * where each subdet has 256 CTDC/Top 2 TAMEX2/Bot channels, which should be
 * mapped as:
 *  FIB4_Top1 .. FIB4_Top512
 *  FIB4_Bot1 .. FIB4_Bot4
 * This class will then figure out how to map the actual fibers (except for
 * non-ambiguous sorting errors).
 */
class R3BSfibCal2Hit : public FairTask
{
  public:
    enum Direction
    {
        HORIZONTAL,
        VERTICAL
    };
    struct ToT
    {
        ToT(R3BSfibCalData const*, R3BSfibCalData const*, Double_t, Double_t, Double_t);
        R3BSfibCalData const* lead;
        R3BSfibCalData const* trail;
        Double_t lead_ns, tail_ns, tot_ns;
    };
    struct Channel
    {
        std::list<R3BSfibCalData const*> lead_list;
        std::list<ToT> tot_list;
    };

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BSfibCal2Hit(Int_t, enum R3BTCalEngine::CTDCVariant = R3BTCalEngine::CTDCVariant::CTDC_8_12_FWD_250);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSfibCal2Hit();

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

  private:
    Int_t fnEvents;

    double fClockFreq;
    TClonesArray* fCalItems;
    TClonesArray* fTopCalTriggerItems;
    TClonesArray* fBotCalTriggerItems;
    TClonesArray* fHitItems;
    unsigned const* fTopTriggerMap;
    unsigned const* fBotTriggerMap;
    Int_t fNofHitItems;
    // [0=Bot,1=Top][Channel].
    std::vector<Channel> fChannelArray[2];

  public:
    ClassDef(R3BSfibCal2Hit, 1)
};

#endif
