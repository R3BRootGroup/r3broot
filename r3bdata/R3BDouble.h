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

// Copy of RooDouble to store a double value in a TNamed to store it in a tree.

#ifndef R3B_Double
#define R3B_Double

#include "Rtypes.h"

class R3BDouble
{
  public:
    R3BDouble()
        : fValue(0)
    {
    }
    R3BDouble(Double_t value);
    R3BDouble(const R3BDouble& other)
        : fValue(other.fValue)
    {
    }
    virtual ~R3BDouble() {}

    inline operator Double_t() const { return fValue; }
    R3BDouble& operator=(Double_t value)
    {
        fValue = value;
        return *this;
    }

  protected:
    Double_t fValue; // Value payload
};

#endif
