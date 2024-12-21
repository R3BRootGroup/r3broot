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

// --------------------------------------------------------------
// -----             R3BWhiterabbitPropagator               -----
// -----   Created May 28th 2021 by J.L. Rodriguez-Sanchez  -----
// --------------------------------------------------------------

#include "FairLogger.h"
#include "FairRootManager.h"

#include "R3BWhiterabbitPropagator.h"

#include "TClonesArray.h"

R3BWhiterabbitPropagator::R3BWhiterabbitPropagator()
    : R3BWhiterabbitPropagator("R3BWhiterabbitPropagator", 1, "WRMaster")
{
}

R3BWhiterabbitPropagator::R3BWhiterabbitPropagator(const TString& name, Int_t iVerbose, const TString& namewr)
    : FairTask(name, iVerbose)
    , fInputItem(NULL)
    , fNameWR(namewr)
{
}

R3BWhiterabbitPropagator::~R3BWhiterabbitPropagator()
{
    LOG(debug) << "R3BWhiterabbitPropagator::Delete instance";
    if (fInputItem)
    {
        delete fInputItem;
        fInputItem = NULL;
    }
}

InitStatus R3BWhiterabbitPropagator::Init()
{
    LOG(info) << "R3BWhiterabbitPropagator::Init()";
    FairRootManager* frm = FairRootManager::Instance();
    if (NULL == frm)
    {
        LOG(fatal) << "R3BWhiterabbitPropagator::Init() FairRootManager not found";
        return kFATAL;
    }

    // get access to input data
    fInputItem = dynamic_cast<TClonesArray*>(frm->GetObject(fNameWR + "Data"));
    if (!fInputItem)
    {
        LOG(fatal) << "R3BWhiterabbitPropagator::Init() " + fNameWR + "Data not found";
        return kFATAL;
    }

    frm->Register(fNameWR + "Data", fNameWR, fInputItem, kTRUE);

    return kSUCCESS;
}

void R3BWhiterabbitPropagator::Exec(Option_t* option) { return; }

ClassImp(R3BWhiterabbitPropagator);
