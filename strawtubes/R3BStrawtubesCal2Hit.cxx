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

#include "R3BStrawtubesCal2Hit.h"
#include "FairLogger.h"
#include "R3BStrawtubesCalData.h"
#include "R3BStrawtubesHitData.h"
#include "TClonesArray.h"
#include <FairRootManager.h>

R3BStrawtubesCal2Hit::R3BStrawtubesCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems()
    , fHitItems(new TClonesArray("R3BStrawtubesHitData"))
    , fNofHitItems(0)
{
}

R3BStrawtubesCal2Hit::~R3BStrawtubesCal2Hit()
{
    if (fHitItems)
    {
        delete fHitItems;
    }
}

void R3BStrawtubesCal2Hit::Exec(Option_t* option)
{
    Int_t nDets = fCalItems->GetEntriesFast();
    for (Int_t i = 0; i < nDets; i++)
    {
        auto calItem = dynamic_cast<R3BStrawtubesCalData*>(fCalItems->At(i));
        if (!calItem)
        {
            continue;
        }
        new ((*fHitItems)[fNofHitItems]) R3BStrawtubesHitData(*calItem);
        ++fNofHitItems;
    }
}

void R3BStrawtubesCal2Hit::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
    }
    fNofHitItems = 0;
}

void R3BStrawtubesCal2Hit::FinishTask() {}

InitStatus R3BStrawtubesCal2Hit::Init()
{
    FairRootManager* mgr = FairRootManager::Instance();
    if (!mgr)
    {
        LOG(fatal) << "FairRootManager not found";
    }
    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("StrawtubesCal"));
    if (!fCalItems)
    {
        LOG(fatal) << "Branch StrawtubesCal not found";
    }

    mgr->Register("StrawtubesHit", "Land", fHitItems, kTRUE);

    return kSUCCESS;
}

InitStatus R3BStrawtubesCal2Hit::ReInit() { return kSUCCESS; }

ClassImp(R3BStrawtubesCal2Hit)
