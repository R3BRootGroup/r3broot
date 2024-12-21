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
// -----           R3BTwimvsFootOnlineSpectra             -----
// -----    Created 24/02/22 by J.L. Rodriguez-Sanchez    -----
// -----         Fill Twim-Foot online histograms         -----
// ------------------------------------------------------------

#ifndef R3BTwimvsFootOnlineSpectra_H
#define R3BTwimvsFootOnlineSpectra_H 1

#include "FairTask.h"
#include "R3BLogger.h"
#include "TCanvas.h"
#include "TH2F.h"
#include <vector>

class TClonesArray;
class R3BEventHeader;

class R3BTwimvsFootOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTwimvsFootOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTwimvsFootOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTwimvsFootOnlineSpectra();

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
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

    void SetFootDetIDs(UInt_t id1, UInt_t id2)
    {
        if (id1 < 1)
        {
            R3BLOG(fatal, "id1 < 1, it should be at least 1");
        }
        if (id2 < 1)
        {
            R3BLOG(fatal, "id2 < 1, it should be at least 1");
        }
        fFootDetId1 = id1;
        fFootDetId2 = id2;
    }

  private:
    Int_t fNbFootDet;
    UInt_t fFootDetId1;
    UInt_t fFootDetId2;
    Bool_t fMusli;
    TClonesArray* fHitItemsTwim; /* Array with hit-Twim items. */
    TClonesArray* fHitItemsFoot; /* Array with hit-Foot items. */

    R3BEventHeader* fEventHeader; /* Event header.      */
    Int_t fNEvents;               /* Event counter.     */

    // Histograms for Hit data
    std::vector<TH2F*> fh2_hit_e;
    std::vector<TH2F*> fh2_hit_x;

  public:
    ClassDef(R3BTwimvsFootOnlineSpectra, 1)
};

#endif /* R3BTwimvsFootOnlineSpectra_H */
