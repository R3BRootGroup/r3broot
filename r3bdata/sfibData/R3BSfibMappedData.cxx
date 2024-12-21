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
// -----              R3BSfibMappedData                -----
// -----         Created Nov 2020 by A. Falduto        -----
// ----------------------------------------------------------------

#include "R3BSfibMappedData.h"

R3BSfibMappedData::R3BSfibMappedData()
    : fIsBottom(-1)
    , fChannel(-1)
    , fIsLeading()
    , fCoarse(-1)
    , fFine(-1)
{
}

R3BSfibMappedData::R3BSfibMappedData(Int_t a_is_bottom,
                                     Int_t a_channel,
                                     Bool_t a_is_leading,
                                     Int_t a_coarse,
                                     Int_t a_fine)
    : fIsBottom(a_is_bottom)
    , fChannel(a_channel)
    , fIsLeading(a_is_leading)
    , fCoarse(a_coarse)
    , fFine(a_fine)
{
}

R3BSfibMappedData::~R3BSfibMappedData() {}

Int_t R3BSfibMappedData::GetChannel() const { return fChannel; }

Int_t R3BSfibMappedData::GetCoarse() const { return fCoarse; }

Int_t R3BSfibMappedData::GetFine() const { return fFine; }

Bool_t R3BSfibMappedData::IsTop() const { return !fIsBottom; }

Bool_t R3BSfibMappedData::IsBottom() const { return fIsBottom; }

Bool_t R3BSfibMappedData::IsLeading() const { return fIsLeading; }

Bool_t R3BSfibMappedData::IsTrailing() const { return !fIsLeading; }

ClassImp(R3BSfibMappedData)
