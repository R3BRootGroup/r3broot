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

#ifndef R3BRPCPOINT_H
#define R3BRPCPOINT_H

#include "TObject.h"
#include "TVector3.h"

#include "FairMCPoint.h"

class R3BRpcPoint : public FairMCPoint
{

  public:
    /** Default constructor **/
    R3BRpcPoint();

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param ident    channel ID
     *@param posIn    Ccoordinates at entrance to active volume [cm]
     *@param momIn    Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     *@param EventId  Event Identifier
     **/
    R3BRpcPoint(Int_t trackID,
                Int_t detID,
                Int_t ident,
                TVector3 posIn,
                TVector3 momIn,
                Double_t tof,
                Double_t length,
                Double_t eLoss,
                UInt_t EventId);

    /** Copy constructor **/
    R3BRpcPoint(const R3BRpcPoint&);

    R3BRpcPoint& operator=(const R3BRpcPoint&) { return *this; }

    /** Destructor **/
    virtual ~R3BRpcPoint();

    /** Accessors **/
    Int_t GetChannelId() const { return fChannelId; }
    Double_t GetXIn() const { return fX; }
    Double_t GetYIn() const { return fY; }
    Double_t GetZIn() const { return fZ; }
    void PositionIn(TVector3& pos) { pos.SetXYZ(fX, fY, fZ); }
    void MomentumIn(TVector3& mom) { mom.SetXYZ(fPx, fPy, fPz); }

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

  protected:
    Int_t fChannelId; ///< channel index
    // MODIFY ME!!!!!!!!!!!!!!!!!

    ClassDef(R3BRpcPoint, 1)
};

#endif
