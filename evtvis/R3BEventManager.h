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

#pragma once

#include <FairEventManager.h>

class R3BIonName;

class R3BEventManager : public FairEventManager
{
  public:
    static R3BEventManager* Instance();
    R3BEventManager();
    virtual ~R3BEventManager() { ; }
    virtual void AddParticlesToPdgDataBase(Int_t pdgCode);
    virtual void SetScaleByEnergy(Bool_t scale) { fScaleByEnergy = scale; }
    virtual Bool_t IsScaleByEnergy() { return fScaleByEnergy; }

  protected:
    Bool_t fScaleByEnergy; //!

    R3BIonName* fIonName;

  private:
    static R3BEventManager* fgRinstance; //!
    ClassDef(R3BEventManager, 1);
};
