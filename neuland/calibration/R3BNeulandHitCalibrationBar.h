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

#ifndef R3BNEULANDHITCALIBRATIONBAR_H
#define R3BNEULANDHITCALIBRATIONBAR_H

#include <array>
#include <vector>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"

#include "R3BNeulandHitModulePar.h"

class TGraph;
class R3BNeulandHitPar;
class TDirectory;

namespace Neuland
{
    namespace Calibration
    {
        class TSyncer;

        class HitCalibrationBar
        {
            using DPair = std::array<Double_t, 2>;
            using IPair = std::array<Int_t, 2>;

          public:
            enum class CalibrationStatus
            {
                noData = 0,
                completeFail = 1,
                success = 2,
                energyGainFail = 3,
                timeSyncFail = 4,
                timeSyncError = 5,
                strangeParameters = 6,
                timeJump = 7
            };

            HitCalibrationBar(const Int_t id = 0);
            ~HitCalibrationBar();

            void Update(const R3BNeulandHitModulePar* par);

            void Set(const Int_t side, const Double_t time, const Int_t qdc);
            Bool_t Add(const Double_t timeOffset,
                       const Double_t entryPosition,
                       const Double_t exitPosition,
                       const Double_t energy,
                       const UInt_t eventNumber);
            void Reset();
            void Calibrate(TDirectory* histogramDir = nullptr);
            R3BNeulandHitModulePar GetParameters();
            CalibrationStatus GetCalibrationStatus() const;
            Bool_t IsValid() const;
            inline Double_t GetTime() const { return 0.5 * (CurrentHit.Time[0] + CurrentHit.Time[1]); }
            inline Double_t GetPosition() const { return CurrentHit.EntryPosition; }

            void SetGlobalTSync(const Double_t value, const Double_t error);

            inline static const char* const GetCalibrationStatusDescription(const CalibrationStatus status)
            {
                return CalibrationStatusDescription[static_cast<int>(status)];
            }
            inline static const char* const GetCalibrationStatusAbbreviation(const CalibrationStatus status)
            {
                return CalibrationStatusAbbreviation[static_cast<int>(status)];
            }

          private:
            static const char* CalibrationStatusDescription[];
            static const char* CalibrationStatusAbbreviation[];

            Int_t ID;
            Int_t Validity;
            UInt_t FailCounter;
            UInt_t LastEventNumber;

            struct Hit
            {
                DPair Time;
                IPair QDC;
                Double_t EntryPosition;
                Double_t ExitPosition;
                Double_t Energy;
            } CurrentHit;
            std::vector<Hit> LastHits; //!

            Double_t TimeDifference;
            Double_t EffectiveSpeed;
            DPair Gain;
            Double_t InvLightAttenuationLength;
            IPair Pedestal;
            DPair Saturation;
            DPair Threshold;
            DPair TimeSync;

            struct
            {
                TGraphErrors TimeDifference;
                TGraphErrors EffectiveSpeed;
                std::array<TGraphErrors, 2> Gain;
                TGraphErrors LightAttenuationLength;
                std::array<TH1F, 2> TotalHits;
                std::array<TH1F, 2> MissedHits;
                std::array<TGraph, 2> MissRatio;
                std::array<TH1F, 2> Energy;
                std::array<TH2F, 2> Saturation;
                std::array<TH1F, 2> Pedestal;
            } Log; //!

            // Some members so we do not have to create and delete them all the time
            TGraphErrors FitGraph;

            void positionCalibration(int firstHit, int nHits);
            void energyCalibration(int firstHit, int nHits);
            void pedestalCalibration();
            void thresholdCalibration();

            Int_t cleanupFit(TGraph& graph, TF1& fit, Double_t maxDifference) const;

            /**
             * @brief removes points from a graph in an efficient way
             *
             * @param points pointer to the array containing the point indexes
             * @param nPoints number of points which will be removed
             * @param graph the graph
             */
            void removePoints(Int_t* points, Int_t nPoints, TGraph& graph) const;

            /**
             * @brief created and writes histograms to the directory
             *
             * @param histogramDir pointer to the directory. If nullptr nothing is created or stored.
             */
            void createHistograms(TDirectory* histogramDir);

            Double_t getMean(const TGraphErrors& graph, Double_t expectedValue = 0.);
        };
    } // namespace Calibration
} // namespace Neuland
#endif
