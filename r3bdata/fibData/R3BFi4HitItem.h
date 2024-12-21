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

#ifndef R3BFi4HITITEM_H
#define R3BFi4HITITEM_H

#include "R3BHit.h"

// for the data analysis of the s438b PSPX detectors.
// Introduced by Ralf, Jan 2016
// adapted by Joachim March 2016

class R3BFi4HitItem : public R3BHit
{
  public:
    // Default Constructor
    R3BFi4HitItem();

    // Standard Constructor
    R3BFi4HitItem(UInt_t detector,
                  UInt_t fiber,
                  Float_t xpos,
                  Float_t ypos,
                  Float_t zpos, // pos
                  Float_t energy_top,
                  Float_t energy_bottom, // energy
                  Float_t time_top,
                  Float_t time_bottom); // time

    // Destructor
    virtual ~R3BFi4HitItem() {}

    // Getters
    inline const UChar_t& GetDetector() const { return fDetector; }
    inline const UInt_t& GetFiber() const { return fFiber_n; }
    inline const Float_t& GetEnergy_top() const { return fEnergy_top; }
    inline const Float_t& GetEnergy_bottom() const { return fEnergy_bottom; }
    inline const Float_t& GetTime_top() const { return fTime_top; }
    inline const Float_t& GetTime_bottom() const { return fTime_bottom; }

  private:
    UChar_t fDetector;
    UInt_t fFiber_n;
    Float_t fEnergy_top; // 0xfff is energy data, Anything in 0xfffff000 indicates an error or overflow
    Float_t fEnergy_bottom;
    Float_t fTime_top;
    Float_t fTime_bottom;

  public:
    ClassDef(R3BFi4HitItem, 3)
};

#endif
