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

/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Hector Alvarez-Pol     hapol@fpddux.usc.es

//*-- Date: 08/2009
//*-- Last Update: 05/08/09 by <D.Bertini@gsi.de>
// --------------------------------------------------------------
// Description:
//   The information from the beam parameters used for the reaction
//   vertex. Information to be accessed in the ROOT file per event
//
// --------------------------------------------------------------
// Comments:
//     - 05/08/09 changes for Root compat.
//
//     - 07/05/07 Data structure for beam parameters
//
/////////////////////////////////////////////////////////////////

#ifndef R3BBeamInfo_H
#define R3BBeamInfo_H 1

#include "TROOT.h" //for including Rtypes.h

#include "TObject.h"

class R3BBeamInfo : public TObject
{

  private:
    Double_t theta; // theta emission angle in vertex
    Double_t phi;   // phi emission angle in vertex

    Double_t xPos; // vertex position
    Double_t yPos;
    Double_t zPos;

  public:
    R3BBeamInfo();
    ~R3BBeamInfo();

    inline Double_t GetTheta() const { return theta; }
    inline Double_t GetPhi() const { return phi; }

    inline Double_t GetXPos() const { return xPos; }
    inline Double_t GetYPos() const { return yPos; }
    inline Double_t GetZPos() const { return zPos; }

    inline void SetTheta(Double_t t) { theta = t; }
    inline void SetPhi(Double_t p) { phi = p; }

    inline void SetXPos(UInt_t x) { xPos = x; }
    inline void SetYPos(UInt_t y) { yPos = y; }
    inline void SetZPos(UInt_t z) { zPos = z; }

    void SetVertexPosition(Double_t, Double_t, Double_t);
    void SetAngles(Double_t, Double_t);

    void print(void);

    ClassDef(R3BBeamInfo, 1) // ROOT CINT
};
#endif
