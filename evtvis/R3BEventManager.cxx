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

#include "R3BEventManager.h"
#include "FairLogger.h"
#include "R3BIonName.h"
#include "TDatabasePDG.h"
#include "TEveGeoNode.h"
#include "TEveManager.h"
#include "TGeoManager.h"

#include <iostream>
using namespace std;

R3BEventManager* R3BEventManager::fgRinstance = 0;

R3BEventManager* R3BEventManager::Instance() { return fgRinstance; }

R3BEventManager::R3BEventManager()
    : FairEventManager()
{
    LOG(info) << "calling ctor Event Manager";
    fgRinstance = this;
}

void R3BEventManager::AddParticlesToPdgDataBase(Int_t pdgCode)
{

    TDatabasePDG* pdgDB = TDatabasePDG::Instance();

    if (!pdgDB->GetParticle(pdgCode))
    {
        int temp = pdgCode / 1E8;
        temp = pdgCode - temp * 1E8;
        int element = temp / 1E4;
        int mass = (temp - element * 1E4) / 10;

        bool particleRecognised = true;
        char name[20];

        particleRecognised = ((R3BIonName*)fIonName)->GetIonName(element, name);

        if (particleRecognised)
        {
            char title[30];
            sprintf(title, "%s%i", name, mass);
            pdgDB->AddParticle(name, title, mass, kTRUE, 0, 0, "Ion", pdgCode);
        }
    }
}

ClassImp(R3BEventManager)
