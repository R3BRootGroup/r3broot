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

// ------------------------------------------------------------------------
// -----                      R3BWRData                               -----
// -----          Created 28/02/2019 by J.L. Rodriguez                -----
// ------------------------------------------------------------------------

#ifndef R3BWRData_H
#define R3BWRData_H 1

#include "TObject.h"
#include <Rtypes.h>

class R3BWRData : public TObject
{
  public:
    // Default Constructor
    R3BWRData();

    /** Standard Constructor
     *@param timestamp
     *@param id
     **/
    R3BWRData(uint64_t timestamp, UInt_t id = 1);

    // Destructor
    virtual ~R3BWRData() {}

    // Getters
    inline const uint64_t& GetTimeStamp() const { return fTimeStamp; }
    inline const UInt_t& GetId() const { return fId; }

  protected:
    UInt_t fId;
    uint64_t fTimeStamp; // timestamp for master

  public:
    ClassDef(R3BWRData, 1)
};

#endif /* R3BWRData_H */
