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

// ----------------------------------------------------------------
// -----              R3BFibMappedData                -----
// -----         Created Jan 2018 by M.Heil        -----
// ----------------------------------------------------------------

#ifndef R3BFIBMAPPEDITEM_H
#define R3BFIBMAPPEDITEM_H

#include "TObject.h"

class R3BFibMappedData : public TObject
{
  public:
    /* Default Constructor */
    R3BFibMappedData();

    /* Standard Constructor */
    R3BFibMappedData(Int_t planeId, Int_t fiberId);

    // Destructor
    virtual ~R3BFibMappedData() {}

    // Getters
    inline const Int_t& GetPlaneId() const { return fPlane; }
    inline const Int_t& GetFiberId() const { return fFiber; }

    inline const Int_t& GetCoarseTime1LE() const { return fCoarseTime1LE; }
    inline const Int_t& GetFineTime1LE() const { return fFineTime1LE; }
    inline const Int_t& GetCoarseTime1TE() const { return fCoarseTime1TE; }
    inline const Int_t& GetFineTime1TE() const { return fFineTime1TE; }

    inline const Int_t& GetCoarseTime2LE() const { return fCoarseTime2LE; }
    inline const Int_t& GetFineTime2LE() const { return fFineTime2LE; }
    inline const Int_t& GetCoarseTime2TE() const { return fCoarseTime2TE; }
    inline const Int_t& GetFineTime2TE() const { return fFineTime2TE; }

    inline const Int_t& GetFineTime(int t, int e) const
    {
        return t ? (e ? fFineTime2TE : fFineTime2LE) : (e ? fFineTime1TE : fFineTime1LE);
    }

    inline const Int_t& GetCoarseTime(int t, int e) const
    {
        return t ? (e ? fCoarseTime2TE : fCoarseTime2LE) : (e ? fCoarseTime1TE : fCoarseTime1LE);
    }

    // To be compatible with Tacquila data:
    inline Bool_t Is17() const { return false; }

  protected:
    Int_t fPlane; //... number of plane 1..n
    Int_t fFiber; //... number of Fiber   1..n

  public:
    // PM1:
    Int_t fCoarseTime1LE; //... coarse time of leading edge
    Int_t fFineTime1LE;   //... fine time of leading edge
    Int_t fCoarseTime1TE; //... coarse time of trailing edge
    Int_t fFineTime1TE;   //... fine time of trailing edge
    // PM2:
    Int_t fCoarseTime2LE; //... coarse time of leading edge
    Int_t fFineTime2LE;   //... fine time of leading edge
    Int_t fCoarseTime2TE; //... coarse time of trailing edge
    Int_t fFineTime2TE;   //... fine time of trailing edge

  public:
    ClassDef(R3BFibMappedData, 1)
};

#endif
