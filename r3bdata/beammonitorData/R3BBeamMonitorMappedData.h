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

#ifndef R3BBEAMMONITORMAPPEDITEM_H
#define R3BBEAMMONITORMAPPEDITEM_H

#include "TObject.h"

class R3BBeamMonitorMappedData : public TObject
{
  public:
    // Default Constructor
    R3BBeamMonitorMappedData();

    // Standard Constructor
    R3BBeamMonitorMappedData(UInt_t, UInt_t, UInt_t);

    UInt_t GetIC() const;
    UInt_t GetSEETRAM() const;
    UInt_t GetTOFDOR() const;

  public:
    UInt_t fIC;
    UInt_t fSEETRAM;
    UInt_t fTOFDOR;

  public:
    ClassDef(R3BBeamMonitorMappedData, 1)
};

#endif
