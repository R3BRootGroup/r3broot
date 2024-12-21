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

// -------------------------------------------------------------------------
// -----                      R3BDchFullPoint header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BDchFullPoint.h
 **/

#ifndef R3BDCHFULLPOINT_H
#define R3BDCHFULLPOINT_H

#include "TObject.h"
#include "TVector3.h"

#include "FairMCPoint.h"

class R3BDchFullPoint : public FairMCPoint
{
  public:
    /** Default constructor **/
    R3BDchFullPoint();

    /** Constructor with arguments
     *@param trackID  Index of MCTrack
     *@param detID    Detector ID
     *@param posIn    Ccoordinates at entrance to active volume [cm]
     *@param posOut   Coordinates at exit of active volume [cm]
     *@param momIn    Momentum of track at entrance [GeV]
     *@param momOut   Momentum of track at exit [GeV]
     *@param tof      Time since event start [ns]
     *@param length   Track length since creation [cm]
     *@param eLoss    Energy deposit [GeV]
     **/
    R3BDchFullPoint(Int_t trackId,
                    Int_t module,
                    Int_t layer,
                    Int_t cell,
                    TVector3 pos,
                    TVector3 lpos,
                    TVector3 mom,
                    TVector3 lmom,
                    Double_t tof,
                    Double_t length,
                    Double_t eLoss);

    /** Destructor **/
    virtual ~R3BDchFullPoint();

    /** Accessors **/
    Double_t GetLocalX() const { return fLocalX; }
    Double_t GetLocalY() const { return fLocalY; }
    Double_t GetLocalZ() const { return fLocalZ; }

    Double_t GetLocalPx() const { return fLocalPx; }
    Double_t GetLocalPy() const { return fLocalPy; }
    Double_t GetLocalPz() const { return fLocalPz; }

  protected:
    Double32_t fLocalX, fLocalY, fLocalZ;
    Double32_t fLocalPx, fLocalPy, fLocalPz;
    Int_t fModule;
    Int_t fLayer;
    Int_t fCell;

    ClassDef(R3BDchFullPoint, 1)
};

#endif
