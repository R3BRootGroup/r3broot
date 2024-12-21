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

#pragma once

#include <FairEventHeader.h>
#include <stdexcept>

class R3BEventHeader : public FairEventHeader
{
  public:
    R3BEventHeader() = default;

    ~R3BEventHeader() override = default;

    inline void SetExpId(const int expid) { fExpId = expid; }
    inline void SetEventno(const uint64_t eventno) { fEventno = eventno; }
    inline void SetTrigger(const int trigger) { fTrigger = trigger; }
    inline void SetTimeStamp(const uint64_t timeStamp) { fTimeStamp = timeStamp; }
    inline void SetTpat(const int tpat) { fTpat = tpat; }
    inline void SetTStart(const double tStart) { fTStart = tStart; }
    inline void SetTStartSimple(const double tStart) { fTStartSimple = tStart; }
    inline void SetTprev(const double tPrev) { fTprev = tPrev; }
    inline void SetTnext(const double tNext) { fTnext = tNext; }

    [[nodiscard]] inline int GetExpId() const { return fExpId; }
    [[nodiscard]] inline uint64_t GetEventno() const { return fEventno; }
    [[nodiscard]] inline int GetTrigger() const { return fTrigger; }
    [[nodiscard]] inline uint64_t GetTimeStamp() const { return fTimeStamp; }
    [[nodiscard]] inline int GetTpat() const { return fTpat; }

    static constexpr uint32_t MakeTpatBit(uint8_t trigNo)
    {
        return (1 <= trigNo && trigNo <= 16) ? (1 << (trigNo - 1)) : throw std::runtime_error("Bad trigNo.");
    }

    bool HasTpatTrig(int trigNo) const { return fTpat & MakeTpatBit(trigNo); }

    [[nodiscard]] inline double GetTStart() const { return fTStart; }

    [[nodiscard]] inline double GetTStartSimple() const { return fTStartSimple; }

    [[nodiscard]] inline double GetTprev() const { return fTprev; }

    [[nodiscard]] inline double GetTnext() const { return fTnext; }

    void Register(bool Persistance = true) override{};

  private:
    int fExpId = 0;
    uint64_t fEventno = 0;
    int fTrigger = 0;
    uint64_t fTimeStamp = 0;
    int fTpat = 0;
    double fTStart = 0;
    double fTStartSimple = 0;
    double fTprev = 0;
    double fTnext = 0;

  public:
    ClassDefOverride(R3BEventHeader, 10)
};
