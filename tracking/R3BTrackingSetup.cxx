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

#include "R3BTrackingSetup.h"

#include "FairLogger.h"

#include <math.h>

using namespace std;

R3BTrackingSetup::R3BTrackingSetup() {}

R3BTrackingSetup::~R3BTrackingSetup()
{
    for (auto const& x : fDetectors)
    {
        delete x;
    }
    fDetectors.clear();
    fMapIndex.clear();
}

void R3BTrackingSetup::AddDetector(const string& name,
                                   EDetectorType type,
                                   const string& geoParName,
                                   const string& dataName)
{
    int index = fDetectors.size();
    if (dataName.length() > 0)
    {
        fDetectors.push_back(new R3BTrackingDetector(name.c_str(), type, geoParName.c_str(), dataName.c_str()));
    }
    else
    {
        fDetectors.push_back(new R3BTrackingDetector(name.c_str(), type, geoParName.c_str(), dataName.c_str()));
    }
    fMapIndex[name] = index;
}

R3BTrackingDetector* R3BTrackingSetup::GetByName(const string& name)
{
    if (fMapIndex.find(name) == fMapIndex.end())
    {
        LOG(error) << "Detector " << name << " was not found in setup.";
        return nullptr;
    }

    return fDetectors.at(fMapIndex[name]);
}

R3BTrackingDetector* R3BTrackingSetup::GetFirstByType(const EDetectorType& type)
{
    for (auto const& x : fDetectors)
    {
        if (type == x->section)
        {
            return x;
        }
    }

    return nullptr;
}

Double_t R3BTrackingSetup::GetAfterGladResolution()
{
    Double_t res = 0.;
    for (auto const& x : fDetectors)
    {
        if (kAfterGlad == x->section)
        {
            res += pow(x->res_x, 2);
        }
    }
    res = sqrt(res);
    return res;
}

void R3BTrackingSetup::Init()
{
    for (auto const& x : fDetectors)
    {
        x->Init();
    }
}

void R3BTrackingSetup::SetParContainers()
{
    for (auto const& x : fDetectors)
    {
        x->SetParContainers();
    }
}

void R3BTrackingSetup::CopyHits()
{
    for (auto const& x : fDetectors)
    {
        x->CopyHits();
    }
}

ClassImp(R3BTrackingSetup)
