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

#ifndef R3BPDCHITDATA_H
#define R3BPDCHITDATA_H

#include "R3BHit.h"

// for the data analysis of the Pdc detectors.
// Introduced by M.Heil, June 2020

class R3BPdcHitData : public R3BHit
{
  public:
    // Default Constructor
    R3BPdcHitData();

    // Standard Constructor
    R3BPdcHitData(Double_t t, Double_t x, Double_t y, Double_t ELoss = 0, Double_t ID = 0, Int_t wire = 0);

    // Destructor
    virtual ~R3BPdcHitData() {}
    Int_t GetWireId() const;

  private:
    Int_t fWireId;

    ClassDef(R3BPdcHitData, 2)
};

#endif
