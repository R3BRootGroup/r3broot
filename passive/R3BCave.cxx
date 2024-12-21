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

#include "R3BCave.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairGeoVolume.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"

#include "TObjArray.h"
#include <TList.h>
#include <stddef.h>

#include "R3BGeoCave.h"
#include "R3BGeoPassivePar.h"

R3BCave::R3BCave()
    : FairModule()
{
}

R3BCave::R3BCave(const char* name, const char* Title)
    : FairModule(name, Title)
{
    world[0] = 0;
    world[1] = 0;
    world[2] = 0;
}

R3BCave::R3BCave(const R3BCave& right)
    : FairModule(right)
{
    world[0] = right.world[0];
    world[1] = right.world[1];
    world[2] = right.world[2];
}

R3BCave::~R3BCave() {}

void R3BCave::ConstructGeometry()
{
    FairGeoLoader* loader = FairGeoLoader::Instance();
    FairGeoInterface* GeoInterface = loader->getGeoInterface();
    R3BGeoCave* MGeo = new R3BGeoCave();
    MGeo->setGeomFile(GetGeometryFileName());
    GeoInterface->addGeoModule(MGeo);
    Bool_t rc = GeoInterface->readSet(MGeo);
    if (rc)
    {
        MGeo->create(loader->getGeoBuilder());
    }

    TList* volList = MGeo->getListOfVolumes();
    // store geo parameter

    FairRun* fRun = FairRun::Instance();
    FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();

    R3BGeoPassivePar* par = dynamic_cast<R3BGeoPassivePar*>(rtdb->getContainer("R3BGeoPassivePar"));
    TObjArray* fSensNodes = par->GetGeoSensitiveNodes();
    TObjArray* fPassNodes = par->GetGeoPassiveNodes();

    TListIter iter(volList);
    FairGeoNode* node = NULL;
    FairGeoVolume* aVol = NULL;

    while ((node = dynamic_cast<FairGeoNode*>(iter.Next())))
    {
        aVol = dynamic_cast<FairGeoVolume*>(node);
        if (node->isSensitive())
        {
            fSensNodes->AddLast(aVol);
        }
        else
        {
            fPassNodes->AddLast(aVol);
        }
    }
    par->setChanged();
    par->setInputVersion(fRun->GetRunId(), 1);
}

FairModule* R3BCave::CloneModule() const { return new R3BCave(*this); }

ClassImp(R3BCave);
