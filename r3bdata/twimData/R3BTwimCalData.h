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

// -------------------------------------------------------------------
// -----           R3BTwimMappedData header file                 -----
// -----       Created 29/09/19 by J.L. Rodriguez-Sanchez        -----
// -------------------------------------------------------------------

#ifndef R3BTwimCalData_H
#define R3BTwimCalData_H 1

#include "TObject.h"

class R3BTwimCalData : public TObject
{
  public:
    /** Default constructor **/
    R3BTwimCalData();

    /** Constructor with arguments
     *@param secID    Section ID
     *@param anodeID  Anode ID
     *@param dtime    Drift time calibrated in [mm] for the position
     *@param eLoss    Energy deposit [channels], for simulation in [MeV]
     **/
    R3BTwimCalData(UInt_t secID, UInt_t anodeID, Double_t dtime, Double_t energy);

    /** Destructor **/
    virtual ~R3BTwimCalData() {}

    /** Accessors **/
    inline const UInt_t& GetSecID() const { return fSecID; }
    inline const UInt_t& GetAnodeID() const { return fAnodeID; }
    inline const Double_t& GetDTime() const { return fDT; }
    inline const Double_t& GetEnergy() const { return fEnergy; }

    /** Modifiers **/
    void SetSecID(UInt_t id) { fSecID = id; };
    void SetAnodeID(UInt_t id) { fAnodeID = id; };
    void SetDTime(Double_t dt) { fDT = dt; };
    void SetEnergy(Double_t energy) { fEnergy = energy; };

  protected:
    UInt_t fSecID, fAnodeID;
    Double_t fDT, fEnergy;

  public:
    ClassDef(R3BTwimCalData, 1)
};

#endif
