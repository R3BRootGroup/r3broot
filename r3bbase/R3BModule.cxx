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

#include "R3BModule.h"
#include "FairLogger.h"
#include "TGeoManager.h"

R3BModule::R3BModule()
    : FairModule()
{
}

R3BModule::R3BModule(const TString& name, const TString& title, const Bool_t active)
    : FairModule(name, title, active)
{
}

R3BModule::R3BModule(const TString& name,
                     const TString& title,
                     const Bool_t active,
                     const TString& geoFile,
                     const TGeoTranslation& trans,
                     const TGeoRotation& rot)
    : R3BModule(name, title, active, geoFile, { trans, rot })
{
}

R3BModule::R3BModule(const TString& name,
                     const TString& title,
                     const Bool_t active,
                     const TString& geoFile,
                     const TGeoCombiTrans& combi)
    : FairModule(name, title, active)
    , fCombiTrans(combi)
{
    SetGeometryFileName(geoFile);
}

void R3BModule::ConstructGeometry()
{
    if (!GetGeometryFileName().EndsWith(".root"))
    {
        LOG(fatal) << GetName() << " (which is a " << ClassName() << ") geometry file is not specified!";
    }
    ConstructRootGeometry();
}

void R3BModule::ConstructRootGeometry(TGeoMatrix*)
{
    LOG(info) << "R3BModule: Constructing " << GetName() << " (which is a " << ClassName()
              << ") geometry from ROOT file " << GetGeometryFileName() << " ...";
    FairModule::ConstructRootGeometry();

    if (!fCombiTrans.IsIdentity())
    {
        auto n = gGeoManager->GetTopNode()->GetDaughter(gGeoManager->GetTopNode()->GetNdaughters() - 1);
        (dynamic_cast<TGeoNodeMatrix*>(n))->SetMatrix(fCombiTrans.MakeClone());
    }
}

ClassImp(R3BModule)
