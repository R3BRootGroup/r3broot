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

// ----------------------------------------------------------------------
// -----              R3BFiberMappedData header file                -----
// -----         Created 13/04/22 by J.L. Rodriguez-Sanchez         -----
// ----------------------------------------------------------------------

#include "R3BFiberMappedData.h"

R3BFiberMappedData::R3BFiberMappedData()
    : fSide(0)
    , fChannel(0)
    , fIsLeading()
    , fCoarse(-1)
    , fFine(-1)
{
}

R3BFiberMappedData::R3BFiberMappedData(UInt_t side, UInt_t channel, Bool_t is_leading, Int_t coarse, Int_t fine)
    : fSide(side)
    , fChannel(channel)
    , fIsLeading(is_leading)
    , fCoarse(coarse)
    , fFine(fine)
{
}

ClassImp(R3BFiberMappedData);
