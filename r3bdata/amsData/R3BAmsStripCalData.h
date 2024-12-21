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
// -----              R3BAmsStripCalData header file                   -----
// -----             Created 29/05/18  by J.L. Rodriguez-Sanchez       -----
// -------------------------------------------------------------------------

/**  R3BAmsStripCalData.h
 **  A strip CalData is the representation of the
 **  information obtained from each strip in the silicon detector.
 **/

#ifndef R3BAmsStripCalData_H
#define R3BAmsStripCalData_H

#include "FairMultiLinkedData.h"
#include "TObject.h"

class R3BAmsStripCalData : public FairMultiLinkedData
{
  public:
    /** Default constructor **/
    R3BAmsStripCalData();

    /** Constructor with arguments
     *@param fDetId    Detector unique identifier
     *@param fSideId   Side unique identifier
     *@param fStripId  Strip unique identifier
     *@param fEnergy   Total energy deposited on the strip ([GeV] in sim)
     **/
    R3BAmsStripCalData(Int_t detid, Int_t sideid, Int_t stripid, Double32_t energy);

    /** Copy constructor **/
    R3BAmsStripCalData(const R3BAmsStripCalData&);

    R3BAmsStripCalData& operator=(const R3BAmsStripCalData&) { return *this; }

    /** Destructor **/
    virtual ~R3BAmsStripCalData() {}

    /** Accessors **/
    inline const Int_t& GetDetId() const { return fDetId; }
    inline const Int_t& GetSideId() const { return fSideId; }
    inline const Int_t& GetStripId() const { return fStripId; }
    inline const Double32_t& GetEnergy() const { return fEnergy; }

    /** Modifiers **/
    void SetDetId(Int_t detid) { fDetId = detid; }
    void SetSideId(Int_t sideid) { fSideId = sideid; }
    void SetStripId(Int_t stripid) { fStripId = stripid; }
    void SetEnergy(Double32_t energy) { fEnergy = energy; }

    void AddMoreEnergy(Double32_t moreEnergy) { fEnergy += moreEnergy; }

  protected:
    Int_t fDetId;       // detector unique identifier
    Int_t fSideId;      // side unique identifier
    Int_t fStripId;     // strip unique identifier
    Double32_t fEnergy; // total energy in the strip

    ClassDef(R3BAmsStripCalData, 1)
};

#endif
