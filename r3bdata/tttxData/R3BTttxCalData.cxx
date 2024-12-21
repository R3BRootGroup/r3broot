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
// -----                  R3BTttxCalData source file                   -----
// -------------------------------------------------------------------------

#include "R3BTttxCalData.h"
#include <fmt/core.h>

R3BTttxCalData::R3BTttxCalData(uint8_t DetID, uint8_t StripID, double time, double energy)
    : fDetID(DetID)
    , fStripID(StripID)
    , fTime(time)
    , fEnergy(energy)
{
}

std::string R3BTttxCalData::toString() const
{
    return fmt::format("DetID: {}, StripID: {}, Time: {:.2f} ns, Energy: {:.2f} keV",
                       GetDetID(),
                       GetStripID(),
                       GetTime(),
                       GetEnergy());
}

void R3BTttxCalData::Print(const Option_t*) const { std::cout << *this << std::endl; }

// Sort in energy with descending order
Int_t R3BTttxCalData::Compare(const TObject* obj) const
{
    const R3BTttxCalData* other = dynamic_cast<const R3BTttxCalData*>(obj);
    if (!other)
        return 0;
    if (fEnergy > other->GetEnergy())
        return -1;
    if (fEnergy < other->GetEnergy())
        return 1;
    return 0;
}

std::ostream& operator<<(std::ostream& os, const R3BTttxCalData& data)
{
    os << data.toString();
    return os;
}

ClassImp(R3BTttxCalData)
