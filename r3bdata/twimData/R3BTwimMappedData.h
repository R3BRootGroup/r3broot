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
// -----       Created 13/10/19 by J.L. Rodriguez-Sanchez        -----
// -------------------------------------------------------------------

#ifndef R3BTwimMappedData_H
#define R3BTwimMappedData_H 1

#include "TObject.h"

class R3BTwimMappedData : public TObject
{
  public:
    /** Default constructor **/
    R3BTwimMappedData();

    /** Constructor with arguments
     *@param secID    Section ID attention 0-based !
     *@param anodeID  Anode ID attention 0-based !
     *@param dtime    Time [channels]
     *@param eLoss    Energy deposit [channels]
     **/
    R3BTwimMappedData(UInt_t secID, UInt_t anodeID, Int_t time, Int_t energy, Bool_t pu, Bool_t ov);

    /** Destructor **/
    virtual ~R3BTwimMappedData() {}

    /** Accessors **/
    inline const UInt_t& GetSecID() const { return fSecID; }
    inline const UInt_t& GetAnodeID() const { return fAnodeID; }
    inline const Int_t& GetTime() const { return fTime; }
    inline const Int_t& GetEnergy() const { return fEnergy; }
    inline const Bool_t& GetPileupStatus() const { return fPileup; }
    inline const Bool_t& GetOverflowStatus() const { return fOverflow; }

    /** Modifiers **/
    void SetSecID(UInt_t id) { fSecID = id; };
    void SetAnodeID(UInt_t id) { fAnodeID = id; };
    void SetDTime(Int_t time) { fTime = time; };
    void SetEnergy(Int_t energy) { fEnergy = energy; };
    void SetPileup(Bool_t pu) { fPileup = pu; }
    void SetOverflow(Bool_t ov) { fOverflow = ov; }

  protected:
    UInt_t fSecID;
    UInt_t fAnodeID;
    Int_t fTime;
    Int_t fEnergy;
    Bool_t fPileup;
    Bool_t fOverflow;

  public:
    ClassDef(R3BTwimMappedData, 1)
};

#endif
