/******************************************************************************
 **   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 **   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
 **                                                                            *
 **             This software is distributed under the terms of the            *
 **                 GNU General Public Licence (GPL) version 3,                *
 **                    copied verbatim in the file "LICENSE".                  *
 **                                                                            *
 ** In applying this license GSI does not waive the privileges and immunities  *
 ** granted to it by virtue of its status as an Intergovernmental Organization *
 ** or submit itself to any jurisdiction.                                      *
 *******************************************************************************/

// -------------------------------------------------------------------------
// // -----                  R3BMwpcHitData header file                   -----
// // -----            Created 09/10/19 by JL Rodriguez-Sanchez           -----
// // -------------------------------------------------------------------------
#ifndef R3BMwpcHitData_H
#define R3BMwpcHitData_H 1
#include "TObject.h"

class R3BMwpcHitData : public TObject
{

  public:
    /** Default constructor **/
    R3BMwpcHitData();

    /** Constructor with arguments
     *@param fX  Position X in [mm]
     *@param fY  Position Y in [mm]
     **/
    R3BMwpcHitData(Double_t X, Double_t Y, Int_t Plane = 0);

    /** Destructor **/
    virtual ~R3BMwpcHitData() {}

    /** Accessors **/
    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetY() const { return fY; }
    inline const Int_t& GetPlane() const { return fPlane; }

    /** Modifiers **/
    void SetX(Double_t x) { fX = x; };
    void SetY(Double_t y) { fY = y; };
    void SetPlane(Int_t plane) { fPlane = plane; };

  protected:
    Double_t fX, fY;
    Int_t fPlane;

  public:
    ClassDef(R3BMwpcHitData, 1)
};

#endif
