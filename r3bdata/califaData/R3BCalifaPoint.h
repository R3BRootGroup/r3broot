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

#pragma once

#include <FairMCPoint.h>
#include <TObject.h>
#include <TVector3.h>

class R3BCalifaPoint : public FairMCPoint
{
  public:
    // Default constructor
    R3BCalifaPoint() = default;

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param trackPID Particle PID
     *@param cryId    Crystal ID
     *@param posIn    Ccoordinates at entrance to active volume [cm]
     *@param momIn    Momentum of track at entrance [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [MeV]
     **/
    R3BCalifaPoint(int trackID,
                   int detID,
                   int trackPID,
                   int cryId,
                   TVector3 posIn,
                   TVector3 momIn,
                   double tof,
                   double length,
                   double eLoss);

    // Copy constructor
    R3BCalifaPoint(const R3BCalifaPoint&);
    R3BCalifaPoint& operator=(const R3BCalifaPoint&) { return *this; }

    // Destructor
    ~R3BCalifaPoint() override = default;

    /** Accessors **/
    int GetCrystalId() const { return fCrystalId; }
    int GetTrackPid() const { return fTrackPID; }
    double GetXIn() const { return fX; }
    double GetYIn() const { return fY; }
    double GetZIn() const { return fZ; }
    void PositionIn(TVector3& pos) { pos.SetXYZ(fX, fY, fZ); }

  protected:
    int fCrystalId = 0; // Crystal index
    int fTrackPID = 0;  // PID

    ClassDefOverride(R3BCalifaPoint, 2)
};
