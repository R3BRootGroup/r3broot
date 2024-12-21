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
// -----      R3BFibervsTofDOnlineSpectra source file        -----
// -----    Created 29/04/22 by J.L. Rodriguez-Sanchez       -----
// ---------------------------------------------------------------

#ifndef R3BFibervsTofDOnlineSpectra_H
#define R3BFibervsTofDOnlineSpectra_H 1

#include "FairTask.h"
#include <vector>

class TClonesArray;
class TH2F;
class R3BEventHeader;
class R3BFiberMappingPar;

class R3BFibervsTofDOnlineSpectra : public FairTask
{
  public:
    R3BFibervsTofDOnlineSpectra();

    R3BFibervsTofDOnlineSpectra(const TString, Int_t iVerbose = 1);

    virtual ~R3BFibervsTofDOnlineSpectra();

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
    TString fName;
    UInt_t fNbfibers;
    UInt_t fNbTofdPlanes;
    UInt_t fNbTofdPaddlesPerPlane;
    R3BFiberMappingPar* fMapPar;
    TClonesArray* fHitItems;
    TClonesArray* fHitTofdItems;

    R3BEventHeader* header;
    Int_t fTrigger;
    Int_t fTpat1, fTpat2;
    unsigned long fNEvents;

    std::vector<TH2F*> fh2_cor_posX;
    std::vector<TH2F*> fh2_cor_posY;

  public:
    ClassDef(R3BFibervsTofDOnlineSpectra, 1)
};

#endif /* R3BFibervsTofDOnlineSpectra_H */
