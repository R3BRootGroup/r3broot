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

// -------------------------------------------------------------
// -----                   R3BTrloiiData                   -----
// -----       Created 10/11/2021 by J.L. Rodriguez        -----
// -------------------------------------------------------------

#ifndef R3BTrloiiData_H
#define R3BTrloiiData_H

#include "TObject.h"

class R3BTrloiiData : public TObject
{
  public:
    // Default Constructor
    R3BTrloiiData();

    /** Standard Constructor
     *@param type       Type of data
     *@param ch         Scaler channel
     *@param counts     Number of counts
     **/
    R3BTrloiiData(UInt_t type, Int_t ch, uint32_t counts);

    // Destructor
    virtual ~R3BTrloiiData() {}

    // Getters
    inline const UInt_t& GetType() const { return fType; }
    inline const Int_t& GetCh() const { return fCh; }
    inline const uint32_t& GetCounts() const { return fCounts; }

  protected:
    UInt_t fType;     // Type of data
    Int_t fCh;        // Scaler channel
    uint32_t fCounts; // Number of counts

  public:
    ClassDef(R3BTrloiiData, 1)
};

#endif /* R3BTrloiiData */
