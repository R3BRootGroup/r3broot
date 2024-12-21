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

#ifndef R3BPSPXHITDATA_H
#define R3BPSPXHITDATA_H

#include "TObject.h"

/**
 * Class containing PSPX detector data on Hit level.
 * @author Ina Syndikus
 * @since April 2016
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxHitData : public TObject
{
  public:
    /** Default Constructor **/
    R3BPspxHitData();

    /** Standard Constructor **/
    R3BPspxHitData(Float_t energy, Float_t pos);

    /** Destructor **/
    virtual ~R3BPspxHitData() {}

    // Getters
    inline const Float_t& GetEnergy() const { return fEnergy; }
    inline const Float_t& GetPos() const { return fPos; }

  private:
    Float_t fEnergy; // Energy/Collected charge in MeV
    Float_t fPos;    // Position in x cm in the range [-l/2:l/2]

  public:
    ClassDef(R3BPspxHitData, 7)
};

#endif
