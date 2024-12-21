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

#ifndef R3BPSPXCALDATA_H
#define R3BPSPXCALDATA_H

#include "TObject.h"

/**
 * Class containing PSPX detector data on Cal level.
 * @author Ralf Plag, Ina Syndikus
 * @since January 2016
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxCalData : public TObject
{
  public:
    /** Default Constructor **/
    R3BPspxCalData();

    /** Standard Constructor **/
    R3BPspxCalData(Int_t strip, Float_t energy, Float_t pos);

    /** Destructor **/
    virtual ~R3BPspxCalData() {}

    // Getters
    inline const Int_t& GetStrip() const { return fStrip; }
    inline const Float_t& GetEnergy() const { return fEnergy; }
    inline const Float_t& GetPos() const { return fPos; }

  private:
    Int_t fStrip;    // Strip number, counting from 1
    Float_t fEnergy; // Total energy/collected charge from both sides of hit strip. This value is gainmatched.
    Float_t fPos;    // Position calculated from energy difference from both sides. Internal coordinates from -1 to +1

  public:
    ClassDef(R3BPspxCalData, 7)
};

#endif
