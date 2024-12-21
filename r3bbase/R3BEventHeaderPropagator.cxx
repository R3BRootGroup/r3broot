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

#include <FairRootManager.h>

#include "R3BEventHeaderPropagator.h"
#include "R3BFileSource.h"
#include "R3BLogger.h"

R3BEventHeaderPropagator::R3BEventHeaderPropagator()
    : R3BEventHeaderPropagator("R3BEventHeaderPropagator", 1, "EventHeader.")
{
}

R3BEventHeaderPropagator::R3BEventHeaderPropagator(const TString& name, Int_t iVerbose, std::string_view nameheader)
    : FairTask(name, iVerbose)
    , fNameHeader(nameheader)
{
}

InitStatus R3BEventHeaderPropagator::Init()
{
    R3BLOG(info, "");
    FairRootManager* frm = FairRootManager::Instance();
    fHeader = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader."));
    R3BLOG_IF(fatal, !fHeader, "EventHeader. not found.");
    R3BLOG_IF(info, fHeader, "EventHeader. found.");

    frm->Register(fNameHeader.data(), "EventHeader", fHeader, kTRUE);

    fSource = frm->GetSource();
    R3BLOG_IF(fatal, !fSource, "R3BFileSource not found.");

    return kSUCCESS;
}

void R3BEventHeaderPropagator::Exec(Option_t* /*option*/)
{
    if (fSource != nullptr)
    {
        fHeader->SetRunId(fSource->GetRunId());
    }
}

ClassImp(R3BEventHeaderPropagator);
