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

// ---------------------------------------------------------------
// -----  R3BFiberMAPMTCorrelationOnlineSpectra source file  -----
// -----    Created 27/04/22 by J.L. Rodriguez-Sanchez       -----
// ---------------------------------------------------------------

#pragma once

#include <FairTask.h>

class TClonesArray;
class TH2F;
class R3BEventHeader;
class R3BFiberMappingPar;

class R3BFiberMAPMTCorrelationOnlineSpectra : public FairTask
{
  public:
    R3BFiberMAPMTCorrelationOnlineSpectra();

    R3BFiberMAPMTCorrelationOnlineSpectra(const TString, const TString, Int_t iVerbose = 1);

    virtual ~R3BFiberMAPMTCorrelationOnlineSpectra();

    virtual InitStatus Init();

    virtual InitStatus ReInit();

    virtual void Exec(Option_t* option);

    virtual void FinishEvent();

    virtual void FinishTask();

    virtual void SetParContainers();

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }

    void Reset_Histo();

  private:
    TString fName1, fName2;
    UInt_t fNbfibers1, fNbfibers2;
    R3BFiberMappingPar* fMapPar1;
    R3BFiberMappingPar* fMapPar2;
    TClonesArray* fHitItems1;
    TClonesArray* fHitItems2;

    R3BEventHeader* header;
    Int_t fTrigger;
    Int_t fTpat1, fTpat2;
    unsigned long fNEvents;

    TH2F* fh_Fib_posX;
    TH2F* fh_Fib_posY;
    TH2F* fh_Fib_tot;

  public:
    ClassDef(R3BFiberMAPMTCorrelationOnlineSpectra, 1)
};
