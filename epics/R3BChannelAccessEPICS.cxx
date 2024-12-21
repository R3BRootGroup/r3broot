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

/* R3BChannelAccessEPICS.cxx
 * R3BROOT
 *
 * Author: Hans Törnqvist <h.toernqvist@gsi.de>
 * Date: 06/12/17
 * */

#include "R3BChannelAccessEPICS.h"
#include <FairLogger.h>

R3BChannelAccessEPICS::R3BChannelAccessEPICS(TString const& name)
    : fName(name)
    , fValue(0)
    , fId()
    , fType(DBR_DOUBLE)
{
    ca_search(name.Data(), &fId);
}

R3BChannelAccessEPICS::~R3BChannelAccessEPICS() {}

double R3BChannelAccessEPICS::Get() { return fValue; }

double* R3BChannelAccessEPICS::GetPointer() { return &fValue; }

chid R3BChannelAccessEPICS::GetID() const { return fId; }

chtype R3BChannelAccessEPICS::GetType() const { return fType; }

void R3BChannelAccessEPICS::Set(double value) { fValue = value; }

ClassImp(R3BChannelAccessEPICS)

    R3BChannelAccessGroupEPICS::R3BChannelAccessGroupEPICS()
    : fChannelArray("R3BChannelAccessEPICS")
    , fNeedSearch(false)
{
}

R3BChannelAccessGroupEPICS::~R3BChannelAccessGroupEPICS() {}

R3BChannelAccess* R3BChannelAccessGroupEPICS::CreateChannel(TString const& name)
{
    fNeedSearch = true;
    return new (fChannelArray[fChannelArray.GetEntriesFast()]) R3BChannelAccessEPICS(name);
}

bool R3BChannelAccessGroupEPICS::Commit()
{
    bool ret = true;
    ret &= Search();
    for (Int_t i = 0; i < fChannelArray.GetEntriesFast(); ++i)
    {
        auto ca = (R3BChannelAccessEPICS*)fChannelArray[i];
        ca_put(ca->GetType(), ca->GetID(), ca->GetPointer());
    }
    ret &= PendIO("get");
    return ret;
}

bool R3BChannelAccessGroupEPICS::Fetch()
{
    bool ret = true;
    ret &= Search();
    for (Int_t i = 0; i < fChannelArray.GetEntriesFast(); ++i)
    {
        auto ca = (R3BChannelAccessEPICS*)fChannelArray[i];
        ca_get(ca->GetType(), ca->GetID(), ca->GetPointer());
    }
    ret &= PendIO("get");
    return ret;
}

bool R3BChannelAccessGroupEPICS::PendIO(char const* op)
{
    auto status = ca_pend_io(5.0);
    switch (status)
    {
        case ECA_NORMAL:
            return true;
        case ECA_TIMEOUT:
            LOG(error) << "R3BChannelAccessGroupEPICS::Get : CA " << op << " failed.";
            return false;
        default:
            LOG(error) << "R3BChannelAccessGroupEPICS::Get : Unexpected CA " << op << " error.";
            return false;
    }
}

bool R3BChannelAccessGroupEPICS::Search()
{
    bool ret = true;
    if (fNeedSearch)
    {
        ret = PendIO("search");
        fNeedSearch = false;
    }
    return ret;
}

ClassImp(R3BChannelAccessGroupEPICS)

    R3BChannelAccessMasterEPICS::R3BChannelAccessMasterEPICS()
    : fGroupArray("R3BChannelAccessGroupEPICS")
{
    ca_task_initialize();
}

R3BChannelAccessMasterEPICS::~R3BChannelAccessMasterEPICS() { ca_task_exit(); }

R3BChannelAccessGroup* R3BChannelAccessMasterEPICS::CreateGroup()
{
    return new (fGroupArray[fGroupArray.GetEntriesFast()]) R3BChannelAccessGroupEPICS;
}

ClassImp(R3BChannelAccessMasterEPICS)
