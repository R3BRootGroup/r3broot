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

// ---------------------------------------------------------------------------
// -----                                                                 -----
// -----                      R3BAmsMappedData                           -----
// -----                  Created 12/05/2018 by J.L. Rodriguez           -----
// -----                                                                 -----
// ---------------------------------------------------------------------------

#ifndef R3BAmsMappedData_H
#define R3BAmsMappedData_H
#include "TObject.h"

class R3BAmsMappedData : public TObject
{

  public:
    // Default Constructor
    R3BAmsMappedData();

    /** Standard Constructor
     *@param detId       Detector unique identifier
     *@param stripId     Strip unique identifier
     *@param energy      Total energy deposited in the strip
     **/
    R3BAmsMappedData(Int_t detId, Int_t stripId, Int_t energy);

    // Destructor
    virtual ~R3BAmsMappedData() {}

    // Getters
    inline const Int_t& GetDetectorId() const { return fDetId; }
    inline const Int_t& GetStripId() const { return fStripId; }
    inline const Int_t& GetEnergy() const { return fEnergy; }

  protected:
    Int_t fDetId, fStripId; // detector and strip unique identifiers
    Int_t fEnergy;          // total energy in the strip

  public:
    ClassDef(R3BAmsMappedData, 1)
};

#endif
