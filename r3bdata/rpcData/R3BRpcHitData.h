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

#ifndef R3BRPCSTRIPHITDATA_H
#define R3BRPCSTRIPHITDATA_H

#include "TObject.h"

#include "FairMultiLinkedData.h"

class R3BRpcHitData : public FairMultiLinkedData
{

  public:
    /** Default constructor **/
    R3BRpcHitData();

    /** Constructor with arguments
     *@param fDetId       Detector unique identifier
     *@param fChannelId   Channel unique identifier
     *@param fPos         Position
     *@param fTime        Time
     *@param fCharge      Charge
     *@param Tof          Time of flight with respct to LOS detector
     **/
    R3BRpcHitData(UInt_t detId, UInt_t channel, double time, double pos, double charge, double Tof);

    /** Copy constructor **/
    R3BRpcHitData(const R3BRpcHitData&);
    R3BRpcHitData& operator=(const R3BRpcHitData&) { return *this; }

    /** Destructor **/
    virtual ~R3BRpcHitData() {}

    /** Accessors **/
    UInt_t GetDetId() const { return fDetId; }
    UInt_t GetChannelId() const { return fChannelId; }
    double GetTime() const { return fTime; }
    double GetPos() const { return fPos; }
    double GetCharge() const { return fCharge; }
    double GetTof() const { return fTof; }

    /** Modifiers **/
    void SetDetId(UInt_t detId) { fDetId = detId; }
    void SetChannelId(UInt_t ch) { fChannelId = ch; }
    void SetPos(double pos) { fPos = pos; }
    void SetTime(double tim) { fTime = tim; }
    void SetCharge(double Q) { fCharge = Q; }

  protected:
    // Basic Hit information
    UInt_t fDetId;
    UInt_t fChannelId;
    double fPos;
    double fTime;
    double fCharge;
    double fTof; // Time-of-flight with respect to LOS detector

    ClassDef(R3BRpcHitData, 1)
};

#endif
