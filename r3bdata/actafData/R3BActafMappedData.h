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

#pragma once

#include <TObject.h>
#include <cstdint>
#include <iostream>
#include <string>

class R3BActafMappedData : public TObject
{
  public:
    // Default Constructor
    R3BActafMappedData() = default;

    /** Standard Constructor
     *@param chn       Channel number
     *@param trace     Trace signal for each channel
     **/
    explicit R3BActafMappedData(UInt_t chn, std::vector<UInt_t> trace);

    // Destructor
    virtual ~R3BActafMappedData() = default;

    // Accessors with [[nodiscard]]
    [[nodiscard]] inline const UInt_t& GetChannel() const { return fCh; }
    [[nodiscard]] inline const std::vector<UInt_t>& GetTrace() const { return fTrace; }

    // Support for printing
    [[nodiscard]] std::string toString() const;
    void Print(const Option_t*) const override;

  protected:
    UInt_t fCh;                 // Channel number
    std::vector<UInt_t> fTrace; // Trace signal for each channel

  public:
    ClassDefOverride(R3BActafMappedData, 1)
};

// Operator overloading for printing R3BActafMappedData
std::ostream& operator<<(std::ostream& os, const R3BActafMappedData& data);
