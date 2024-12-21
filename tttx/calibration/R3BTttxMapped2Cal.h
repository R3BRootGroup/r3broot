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
// -----             R3BTttxMapped2Cal source file                     -----
// -----               Created 18/01/24  by XXX                        -----
// -------------------------------------------------------------------------

#pragma once

#include <FairTask.h>

#include "R3BTttxCalData.h"
#include "R3BTttxMappedData.h"
// #include "R3BLogger.h"

#include <Rtypes.h>
#include <TArrayF.h>
#include <TRandom.h>

class TClonesArray;
class R3BTttxStripCalPar;

class R3BTttxMapped2Cal : public FairTask
{
  public:
    // Default constructor
    R3BTttxMapped2Cal();

    // Standard constructor
    R3BTttxMapped2Cal(const TString& name, int iVerbose = 1);

    // Destructor
    ~R3BTttxMapped2Cal();

    // Method Exec
    void Exec(Option_t* /*option*/) override;

    // Method SetParContainers
    void SetParContainers() override;

    // Method Init
    InitStatus Init() override;

    // Method ReInit
    InitStatus ReInit() override;

    inline void SetOnline(Bool_t option) { fOnline = option; }

    [[nodiscard]] inline const uint8_t& GetTrefCh() { return fch_tref; }
    [[nodiscard]] inline const uint8_t& GetTrigCh() { return fch_trig; }
    [[nodiscard]] inline const double& GetTimeResolution() { return fTimeResolution; }

    inline void SetTimeWindow(double min, double max)
    {
        fTimeMin = min;
        fTimeMax = max;
    }
    inline void SetTrefCh(uint8_t val) { fch_tref = val; }
    inline void SetTrigCh(uint8_t val) { fch_trig = val; }
    inline void SetTimeResolution(double val) { fTimeResolution = val; }

  private:
    void Reset();
    void SetParameter();

    int NumDets = 0;
    int NumStrips = 0;
    int NumParams = 0;
    uint8_t fch_tref = 33; //(1-base)
    uint8_t fch_trig = 34; // (1-base)
    double fTimeMin = -1e5;
    double fTimeMax = 1e5;               // Maximum time of +/-100us
    double fTimeResolution = 25. / 128.; // For setting 3. in pico-second
    TArrayF* CalParams;

    bool fOnline = false; // Don't store data for online

    R3BTttxStripCalPar* fCal_Par = nullptr;
    TClonesArray* fTttxMapData = nullptr;
    TClonesArray* fTttxCalData = nullptr;

  public:
    class CalStrip
    {
      public:
        CalStrip() = default;
        void Init()
        {
            fmult = 0;
            fEraw.clear();
            fTime.clear();
        };
        void SetRawValues(R3BTttxMappedData* hit)
        {
            fEraw.push_back(static_cast<int>(hit->GetEnergy()));
            fTime.push_back(static_cast<int>(hit->GetTime()));
            fmult++;
        };
        [[nodiscard]] inline const auto& GetMult() { return fmult; };
        auto GetE(int mult) { return fEraw.at(mult); };
        auto GetT(int mult) { return fTime.at(mult); };

      private:
        int fmult{ 0 };
        std::vector<int> fEraw;
        std::vector<int> fTime;
    };

    class CalDetector
    {
      public:
        CalDetector() { ReInit(); };
        void ReInit(int NumS = 32, int NumTref = 1, int NumTrig = 1)
        {
            fStrip.clear();
            fTref.clear();
            fTrig.clear();
            fStrip.resize(NumS);
            fTref.resize(NumTref);
            fTrig.resize(NumTrig);
            Init();
        };
        void Init()
        {
            for (auto& val : fStrip)
            {
                val.Init();
            }
            for (auto& val : fTref)
            {
                val.Init();
            }
            for (auto& val : fTrig)
            {
                val.Init();
            }
        };
        auto& GetStrip(int val) { return fStrip.at(val); };
        auto& GetTref(int val = 0) { return fTref.at(val); };
        auto& GetTrig(int val = 0) { return fTrig.at(val); };

      private:
        std::vector<CalStrip> fStrip;
        std::vector<CalStrip> fTref;
        std::vector<CalStrip> fTrig;
    };

  private:
    std::vector<CalDetector> fTttx;
    void CalculateStrip(int idet, int istrip, CalStrip& fStrip, CalStrip& fTrig);

    // Private method AddCalData
    R3BTttxCalData* AddCalData(uint8_t detid, uint8_t stripid, double time, double energy);

  public:
    ClassDefOverride(R3BTttxMapped2Cal, 1)
};
