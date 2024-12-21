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

#include "R3BDetector.h"
#include "R3BLogger.h"
#include <TGeoManager.h>

R3BDetector::R3BDetector()
    : FairDetector()
    , fCutE(1e-3)
{
}

R3BDetector::R3BDetector(const TString& name, const Int_t detId)
    : FairDetector(name, kTRUE, detId)
    , fCutE(1e-3)
{
    SetVerboseLevel(1);
}

R3BDetector::R3BDetector(const TString& name,
                         const Int_t detId,
                         const TString& geoFile,
                         const TGeoTranslation& trans,
                         const TGeoRotation& rot)
    : R3BDetector(name, detId, geoFile, { trans, rot })
{
}

R3BDetector::R3BDetector(const TString& name, const Int_t detId, const TString& geoFile, const TGeoCombiTrans& combi)
    : FairDetector(name, kTRUE, detId)
    , fCutE(1e-3)
    , fCombiTrans(combi)
{
    SetVerboseLevel(1);
    SetGeometryFileName(geoFile);
}

void R3BDetector::ConstructGeometry()
{
    if (!GetGeometryFileName().EndsWith(".root"))
    {
        R3BLOG(fatal, GetName() << " (which is a " << ClassName() << ") geometry file is not specified");
    }
    ConstructRootGeometry();
}

void R3BDetector::ConstructRootGeometry(TGeoMatrix*)
{
    R3BLOG(info,
           "Constructing " << GetName() << " (which is a " << ClassName() << ") geometry from ROOT file "
                           << GetGeometryFileName() << " ...");
    FairModule::ConstructRootGeometry();

    if (!fCombiTrans.IsIdentity())
    {
        auto n = gGeoManager->GetTopNode()->GetDaughter(gGeoManager->GetTopNode()->GetNdaughters() - 1);
        (dynamic_cast<TGeoNodeMatrix*>(n))->SetMatrix(fCombiTrans.MakeClone());
    }
}

ClassImp(R3BDetector)
