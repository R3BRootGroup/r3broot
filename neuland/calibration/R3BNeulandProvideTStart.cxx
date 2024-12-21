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

#include "R3BNeulandProvideTStart.h"
#include "FairRootManager.h"
#include "R3BEventHeader.h"

R3BNeulandProvideTStart::R3BNeulandProvideTStart()
    : FairTask("R3BNeulandProvideTStart", 0)
    , fNeulandCalData("NeulandCalData")
    , fEventHeader(nullptr)
{
}

InitStatus R3BNeulandProvideTStart::Init()
{
    fNeulandCalData.Init();

    auto ioman = FairRootManager::Instance();
    if (ioman == nullptr)
    {
        throw std::runtime_error("R3BNeulandProvideTStart: No FairRootManager");
    }

    fEventHeader = dynamic_cast<R3BEventHeader*>(ioman->GetObject("EventHeader."));
    if (fEventHeader == nullptr)
    {
        throw std::runtime_error("R3BNeulandProvideTStart: No R3BEventHeader");
    }

    return kSUCCESS;
}

void R3BNeulandProvideTStart::Exec(Option_t*) { fEventHeader->SetTStart(GetTStart()); }

Double_t R3BNeulandProvideTStart::GetTStart() const
{
    const auto calData = fNeulandCalData.Retrieve();

    double tref[2] = { NAN, NAN };
    double eref[2] = { NAN, NAN };
    int fref = 0;

    for (const auto& data : calData)
    {
        const auto side = data->GetSide() - 1; // [1,2] -> [0,1]
        const auto bar = data->GetBarId();

        // proton beam 2020, neuland bar 701 used as a start detector
        if (bar == 701)
        {
            tref[side] = data->GetTime();
            eref[side] = data->GetQdc();
            fref++;
        }
    }

    double pstart = NAN;
    if (fref == 2)
    {
        if (tref[0] - tref[1] < -0.5 * 5. * 2048)
            tref[1] = tref[1] - 5. * 2048;
        if (tref[0] - tref[1] > 0.5 * 5. * 2048)
            tref[0] = tref[0] - 5. * 2048;

        pstart = 0.5 * (tref[0] + tref[1]);
    }

    return pstart;
}

bool R3BNeulandProvideTStart::IsBeam() const { return !std::isnan(GetTStart()); }

ClassImp(R3BNeulandProvideTStart)
