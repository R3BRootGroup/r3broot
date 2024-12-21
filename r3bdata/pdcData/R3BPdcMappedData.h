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

#ifndef R3BPDCMAPPEDITEM_H
#define R3BPDCMAPPEDITEM_H

#include "TObject.h"

class R3BPdcMappedData : public TObject
{
  public:
    R3BPdcMappedData();
    R3BPdcMappedData(UInt_t, UInt_t, UInt_t, UInt_t, UInt_t);

    UInt_t GetPlaneId() const;
    UInt_t GetWireId() const;
    UInt_t GetEdgeId() const;
    UInt_t GetTimeCoarse() const;
    UInt_t GetTimeFine() const;

  private:
    UInt_t fPlane; // 1..n
    UInt_t fWire;  // 1..n
    UInt_t fEdge;  // 1 = leading, 2 = trailing
    UInt_t fTimeCoarse;
    UInt_t fTimeFine;

  public:
    ClassDef(R3BPdcMappedData, 1)
};

#endif
