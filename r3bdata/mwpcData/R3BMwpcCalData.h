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
// -----                  R3BMwpcCalData header file                   -----
// -----            Created 09/10/19 by JL Rodriguez-Sanchez           -----
// -------------------------------------------------------------------------

#ifndef R3BMwpcCalData_H
#define R3BMwpcCalData_H 1

#include "TObject.h"

// This data structure must be filled after pedestal subtraction

class R3BMwpcCalData : public TObject
{
  public:
    /** Default constructor **/
    R3BMwpcCalData();

    /** Constructor with arguments
     *@param plane   1 for Xup, 2 for Xdown, and 3 for Y
     *@param pad     Pad number
     *@param charge  Charge of each pad [channels]
     **/
    R3BMwpcCalData(UInt_t plane, UInt_t pad, Float_t charge);

    /** Destructor **/
    virtual ~R3BMwpcCalData() {}

    /** Accessors **/
    inline const UInt_t& GetPlane() const { return fPlane; }
    inline const UInt_t& GetPad() const { return fPad; }
    inline const Float_t& GetQ() const { return fCharge; }

    /** Modifiers **/
    void SetPlane(UInt_t plane) { fPlane = plane; };
    void SetPad(UInt_t pad) { fPad = pad; };
    void SetQ(Float_t charge) { fCharge = charge; };

  protected:
    UInt_t fPlane;
    UInt_t fPad;
    Float_t fCharge;

  public:
    ClassDef(R3BMwpcCalData, 1)
};

#endif
