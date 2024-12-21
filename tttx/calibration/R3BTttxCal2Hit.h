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
// -----             R3BTttxCal2Hit source file                     -----
// -----               Created 18/01/24  by XXX                        -----
// -------------------------------------------------------------------------

#pragma once

#include <FairTask.h>

#include <Rtypes.h>
#include <TArrayF.h>
#include <TRandom.h>

#include "R3BLogger.h"
#include "R3BTttxCalData.h"
#include "R3BTttxHitData.h"
constexpr double WStrip = 95.97 / 32.; // mm
constexpr int NStrips = 32;

class TClonesArray;
class R3BTttxHitPar;

class R3BTttxCal2Hit : public FairTask
{
  public:
    // Default constructor
    R3BTttxCal2Hit();

    // Standard constructor
    R3BTttxCal2Hit(const TString& name, int iVerbose = 1);

    // Destructor
    ~R3BTttxCal2Hit();

    // Method Exec
    void Exec(Option_t* /*option*/) override;

    // Method SetParContainers
    void SetParContainers() override;

    // Method Init
    InitStatus Init() override;

    // Method ReInit
    InitStatus ReInit() override;

    inline void SetOnline(Bool_t option) { fOnline = option; }

  public:
    class HitDetector
    {
      public:
        HitDetector() = default;
        HitDetector(R3BTttxCalData* hit);
        ~HitDetector();
        void AddHit(R3BTttxCalData* hit);

        [[nodiscard]] inline const auto& GetX() { return fX; };
        [[nodiscard]] inline const auto& GetE() { return fE; };
        [[nodiscard]] inline const auto& GetT() { return fT; };
        [[nodiscard]] inline const auto& GetDetID() { return fDetID; };
        [[nodiscard]] inline const auto& IsTrefIn() { return fTrefIncluded; };

      private:
        double StripID2X(int ID)
        { // 1-oriented
            return (static_cast<double>(ID) - static_cast<double>(NStrips) / 2. - 0.5) * WStrip;
        }

        double fX{ NAN };
        double fE{ NAN };
        double fE2{ NAN }; // Squared energy used in AddHit()
        double fT{ NAN };
        int fDetID{ -1 };
        bool fTrefIncluded{ false };
    };

    class HitsDetector
    {
      public:
        HitsDetector() = default;

        void Init();
        void NewCluster(R3BTttxCalData* hit);
        void AddHitToLast(R3BTttxCalData* hit);
        auto GetMult() { return fHits.size(); };
        auto& GetHit(int mult) { return fHits.at(mult); }
        /*
        [[nodiscard]] inline const auto& GetMult() { return fHits.size(); };
        auto GetX(int mult) { return fHits.at(mult).fX; };
        auto GetE(int mult) { return fHits.at(mult).fE; };
        auto GetT(int mult) { return fHits.at(mult).fT; };
        */
      private:
        std::vector<HitDetector*> fHits;
        int fDetID;
    };

    [[nodiscard]] inline const bool& GetRequireTref() { return fRequireTref; }
    [[nodiscard]] inline const bool& GetStoreOnlyTwoDetHitEvent() { return fStoreOnlyTwoDetHitEvent; }
    [[nodiscard]] inline const double& GetFixedBeta() { return ffixedbeta; }
    [[nodiscard]] inline const double& GetTimeWindow() { return ftimewindow; }
    [[nodiscard]] inline const uint& GetNumNeighbours() { return fneighbours; }
    void SetRequireTref(bool val) { fRequireTref = val; }
    void SetStoreOnlyTwoDetHitEvent(bool val) { fStoreOnlyTwoDetHitEvent = val; }
    void SetFixedBeta(double val) { ffixedbeta = val; }
    void SetTimeWindow(double val) { ftimewindow = val; }
    void SetNumNeighbours(uint val) { fneighbours = val; }

    // For the later update
    // Method to access from external process for the event-by-event reconstruction.
    // void RecalculateZ(double fbeta = -1.);

    // Return calculated atomic number Z with the given velocity
    double CalculateZ(HitDetector* hit, double fbeta = -1.);

  private:
    bool IsInWindow(R3BTttxCalData* hit1, R3BTttxCalData* hit2);

    void Reset();
    void InitVHitsDetector();
    void ResetVHitsDetector();
    void SetParameter();

    int NumDets = 2;
    int NumPars = 2;

    bool fOnline = false;                  // Don't store data for online
    bool fRequireTref = false;             // Requires Tref in the cluster
    bool fStoreOnlyTwoDetHitEvent = false; // Store only events with two detectors have hits
    double ffixedbeta = 0.7;               // Velocity for the fixed analysis
    double ftimewindow = 1.e4;             // in nano seconds. Default = 10 us.
    uint fneighbours = 2;                  // Num of neighbouring strips, if 0, no clustering

    R3BTttxHitPar* fHit_Par = nullptr;
    TArrayF* HitParams = nullptr;
    TClonesArray* fTttxCalData = nullptr;
    TClonesArray* fTttxHitData = nullptr;

    // For clustering
    std::vector<HitsDetector*> vHitsDet;
    std::vector<bool> vUsed;

  private:
    // Private method AddHitData
    R3BTttxHitData* AddHitData(int8_t idet, double xpos, double energy, double ang, double charge, double time);
    R3BTttxHitData* AddHitData(HitDetector* hit1);
    R3BTttxHitData* AddHitData(HitDetector* hit1, HitDetector* hit2);

  public:
    ClassDefOverride(R3BTttxCal2Hit, 1)
};
