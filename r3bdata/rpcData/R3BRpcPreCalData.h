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

#ifndef R3BRPCPRECALDATA_H
#define R3BRPCPRECALDATA_H

#include "FairMultiLinkedData.h"
#include "TObject.h"

class R3BRpcPreCalData : public FairMultiLinkedData
{
  public:
    /** Default constructor **/
    R3BRpcPreCalData();

    /** Standard Constructor
     *@param DetId       Detector unique identifier
     *@param channelId   Channel unique identifier
     *@param Time        Time
     *@param Tot         Coarse Time
     *@param Side        left or right
     **/
    R3BRpcPreCalData(UShort_t DetId, UShort_t channelId, double Time, double Tot, UShort_t Side);

    /** Copy constructor **/
    R3BRpcPreCalData(const R3BRpcPreCalData&);

    R3BRpcPreCalData& operator=(const R3BRpcPreCalData&) { return *this; }

    /** Destructor **/
    virtual ~R3BRpcPreCalData() {}

    /** Accessors **/
    inline UShort_t GetDetId() const { return fDetId; }
    inline UShort_t GetChannelId() const { return fChannelId; }
    inline double GetTime() const { return fTime; }
    inline double GetTot() const { return fTot; }
    inline UShort_t GetSide() const { return fSide; }

  protected:
    UShort_t fDetId;     // Channel unique identifier
    UShort_t fChannelId; // Channel unique identifier
    double fTime;        // Fine time
    double fTot;         // Coarse time
    Short_t fSide;       // Top Or Bottom

  public:
    ClassDef(R3BRpcPreCalData, 1)
};

#endif
