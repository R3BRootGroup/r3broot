/******************************************************************************
 *   Copyright (C) 2025 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2025 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BActafMappedData.h"
#include <fmt/core.h>

R3BActafMappedData::R3BActafMappedData(UInt_t chn, std::vector<UInt_t> trace)
    : fCh(chn)
    , fTrace(trace)
{
}

std::string R3BActafMappedData::toString() const
{
    std::string trace_str = "[";
    for (size_t i = 0; i < fTrace.size(); ++i)
    {
        trace_str += std::to_string(fTrace[i]);
        if (i < fTrace.size())
            trace_str += ", ";
    }
    trace_str += "]";

    return fmt::format("Channel: {}, Trace: {}", GetChannel(), trace_str);
}

void R3BActafMappedData::Print(const Option_t*) const { std::cout << *this << std::endl; }

std::ostream& operator<<(std::ostream& os, const R3BActafMappedData& data)
{
    os << data.toString();
    return os;
}

ClassImp(R3BActafMappedData)
