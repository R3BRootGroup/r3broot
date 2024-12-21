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

#ifndef R3BRPCCALDATA_H
#define R3BRPCCALDATA_H

#include "FairMultiLinkedData.h"
#include "TObject.h"

class R3BRpcCalData : public FairMultiLinkedData
{
  public:
    /** Default constructor **/
    R3BRpcCalData();

    /** Standard Constructor
     *@param DetId       Detector identifier
     *@param channelId   Channel unique identifier
     *@param Time        Time (R/L,B/T)
     *@param Tot         Time over threshold (R/L,B/T)
     **/
    R3BRpcCalData(UShort_t DetId, UShort_t channelId, double TimeR_B, double TimeL_T, double TotR_B, double TotL_T);

    /** Copy constructor **/
    R3BRpcCalData(const R3BRpcCalData&);

    R3BRpcCalData& operator=(const R3BRpcCalData&) { return *this; }

    /** Destructor **/
    virtual ~R3BRpcCalData() {}

    /** Accessors **/
    inline UShort_t GetDetId() const { return fDetId; }
    inline UShort_t GetChannelId() const { return fChannelId; }
    inline double GetTimeL_T() const { return fTimeL_T; }
    inline double GetTimeR_B() const { return fTimeR_B; }
    inline double GetTotL_T() const { return fTotL_T; }
    inline double GetTotR_B() const { return fTotR_B; }

  protected:
    UShort_t fDetId;     // Detector identifier
    UShort_t fChannelId; // Channel unique identifier
    double fTimeR_B;     // Time from Right strip or Bottom Pmt channel
    double fTimeL_T;     // Time from Left or Top Pmt channel
    double fTotR_B;      // Tot from the Right or Bottom Pmt channel
    double fTotL_T;      // Tot from the Left or Top Pmt channel

  public:
    ClassDef(R3BRpcCalData, 1)
};

#endif
