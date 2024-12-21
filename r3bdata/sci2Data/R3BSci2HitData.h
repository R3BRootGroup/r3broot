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

#ifndef R3BSCI2HITITEM_H
#define R3BSCI2HITITEM_H

#include "TObject.h"

// for the data analysis of the SCI2 detectors.

class R3BSci2HitData : public TObject
{
  public:
    // Default Constructor
    R3BSci2HitData();

    // Standard Constructor
    R3BSci2HitData(Int_t detid, Double_t xpos, Double_t time, Double_t time_w_tref);

    // Destructor
    virtual ~R3BSci2HitData() {}

    /** Accessors **/
    inline const Int_t& GetSciId() const { return fSciId; }
    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetTime() const { return fTime; }
    inline const Double_t& GetTimeWithTref() const { return fTimeWithTref; }

  protected:
    Int_t fSciId;
    Double_t fX, fTime, fTimeWithTref;

  public:
    ClassDef(R3BSci2HitData, 1)
};

#endif
