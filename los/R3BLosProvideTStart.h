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

#ifndef R3BROOT_R3BLOSPROVIDETSTART_H
#define R3BROOT_R3BLOSPROVIDETSTART_H 1

#include "FairTask.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosTCalData.h"
#include "TCAConnector.h"

class R3BEventHeader;
class R3BCoarseTimeStitch;

class R3BLosProvideTStart : public FairTask
{
  public:
    R3BLosProvideTStart();

    InitStatus Init() override;
    void Exec(Option_t*) override;
    // Setter to fix trigger peak window
    void SetWindow(Double_t left, Double_t right)
    {
        edgeL = left;
        edgeR = right;
        fUseTrigHit = kTRUE;
    }

  private:
    TCAOptionalInputConnector<R3BLosCalData> fLosCalData;
    TCAOptionalInputConnector<R3BLosCalData> fLosTriggerCalData;
    TCAOptionalInputConnector<R3BLosHitData> fLosHitData;
    TCAOptionalInputConnector<R3BLosTCalData> fLosTriggerData;
    R3BEventHeader* fEventHeader;
    R3BCoarseTimeStitch* fTimeStitch;

    bool IsBeam() const;
    Double_t GetTStart() const;
    Double_t GetTStartTrigHit() const;
    Double_t GetTStart_without_trigger() const;

    Bool_t fUseTrigHit;

    Double_t edgeL, edgeR;

    ClassDefOverride(R3BLosProvideTStart, 0)
};

#endif // R3BROOT_R3BLOSPROVIDETSTART_H
