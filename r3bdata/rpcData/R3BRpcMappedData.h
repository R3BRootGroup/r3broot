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

#ifndef R3BRpcMappedData_H
#define R3BRpcMappedData_H

#include "TObject.h"
#include <stdint.h>

class R3BRpcMappedData : public TObject
{

  public:
    // Default Constructor
    R3BRpcMappedData();

    /** Standard Constructor
     *@param detId       Detector Id
     *@param ChannelId   Channel Id
     *@param fineTime    Fine Time
     *@param coarseTime  Coarse Time
     *@param side        Side
     *@param edge        Leading or Trailing
     **/
    R3BRpcMappedData(UShort_t DetId,
                     UShort_t ChannelId,
                     uint64_t FineTime,
                     uint64_t CoarseTime,
                     UShort_t Edge,
                     UShort_t Side);

    // Destructor
    virtual ~R3BRpcMappedData() {}

    // Getters
    inline UShort_t GetDetId() const { return fDetId; }
    inline UShort_t GetChannelId() const { return fChannelId; }
    inline uint64_t GetCoarseTime() const { return fCoarseTime; }
    inline uint64_t GetFineTime() const { return fFineTime; }
    inline UShort_t GetEdge() const { return fEdge; }
    inline UShort_t GetSide() const { return fSide; }

  protected:
    UShort_t fDetId;      // Detector unique identifier
    UShort_t fChannelId;  // Channel unique identifier
    uint64_t fFineTime;   // Fine time
    uint64_t fCoarseTime; // Coarse time
    UShort_t fEdge;       // Leading or Trailing
    UShort_t fSide;       // Left Or Right

  public:
    ClassDef(R3BRpcMappedData, 1)
};

#endif
