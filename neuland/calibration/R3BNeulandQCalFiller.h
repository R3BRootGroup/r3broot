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

#ifndef R3BNEULANDQCALFILLER_H
#define R3BNEULANDQCALFILLER_H

#include "FairTask.h"

class TClonesArray;
class R3BEventHeader;

class R3BNeulandQCalFiller : public FairTask
{
  public:
    R3BNeulandQCalFiller();
    ~R3BNeulandQCalFiller();

    virtual InitStatus Init();
    virtual void Exec(Option_t* option);

  private:
    TClonesArray* data;
    R3BEventHeader* header;
    Int_t i = 0;

  public:
    ClassDef(R3BNeulandQCalFiller, 1);
};

#endif
