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

#include "R3BNeulandHitCalibrationBar.h"
#include "R3BNeulandHitModulePar.h"
#include "R3BNeulandHitPar.h"

#include "FairLogger.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TPad.h"

#include <limits>
#include <numeric>
#include <string>

using DPair = std::array<Double_t, 2>;

// FIXME namespace Neuland::Calibration { with c++17
namespace Neuland
{
    namespace Calibration
    {
        bool GetJumps(double* data,
                      unsigned int size,
                      double_t threshold,
                      unsigned int averageSize,
                      std::vector<unsigned int>* output)
        {
            if (size < 2 * averageSize)
                return false;

            const auto scale = 1. / averageSize;

            auto lastAverage = std::accumulate(data, data + averageSize, 0.) * scale;
            const auto last = data + size;
            auto currentPos = data + averageSize;

            while (currentPos + averageSize < last)
            {
                const auto currentAverage = std::accumulate(currentPos, currentPos + averageSize, 0.) * scale;
                if (fabs(currentAverage - lastAverage) > threshold)
                {
                    if (output == nullptr)
                        return true;
                    output->push_back(currentPos - data);
                }
                lastAverage = currentAverage;
                currentPos += averageSize;
            }

            return false;
        }

        template <Int_t iterations = 8>
        const double FastExp(const double val)
        {
            auto exp = 1. + val / (1 << iterations);
            for (auto i = 0; i < iterations; ++i)
            {
                exp *= exp;
            }
            return exp;
        }

        const double GetLightAtPMT(const double energyLoss,
                                   const double position,
                                   const double invLightAttenuationLength,
                                   const int side)
        {
            return energyLoss *
                   FastExp(-invLightAttenuationLength * (TotalBarLength * 0.5 + (1 - 2 * side) * position));
        }

        constexpr auto RightSide = 0;
        constexpr auto LeftSide = 1;

        constexpr auto RightSideHitBit = RightSide;
        constexpr auto LeftSideHitBit = LeftSide;
        constexpr auto PosCalibrationBit = 2;
        constexpr auto TSyncCalibrationBit = 3;
        constexpr auto EnergyCalibrationBit = 4;
        constexpr auto ThresholdCalibrationBit = 5;
        constexpr auto PedestalCalibrationBit = 6;
        constexpr auto TimeJumpedBit = 20;
        constexpr auto LogCompleteBit = 30;
        constexpr auto StrangeBit = 31;

        constexpr auto PositionCalibrationSize = 1024;
        constexpr auto EnergyCalibrationSize = 4 * 1024;
        constexpr auto ThresholdCalibrationSize = 1024;
        constexpr auto CalibrationLogSize = 16 * 1024;
        constexpr auto LogInitialSize = 128;

        constexpr auto MinEnergyDeposit = 1.;     // MeV
        constexpr auto MaxTSyncError = 0.05;      // ns
        constexpr auto TimeJumpThreshold = 0.075; // ns
        constexpr auto MaxFastTDiffError = 0.05;

        constexpr auto MaxNumberOfFails = 10U;

        const void SetStatus(Int_t& var, Int_t statusBit) { var |= (1 << statusBit); }
        const void ClearStatus(Int_t& var, Int_t statusBit) { var &= ~(1 << statusBit); }
        const bool IsStatus(Int_t var, Int_t statusBit) { return (var & (1 << statusBit)); }

        const char* HitCalibrationBar::CalibrationStatusDescription[] = { "No Data recorded",
                                                                          "Calibration Failed",
                                                                          "Calibration was Successful",
                                                                          "Energy Calibration Failed",
                                                                          "Time Synchronisation Failed",
                                                                          "Time Synchronisation Error large",
                                                                          "Parameters look strange",
                                                                          "Time Jump" };
        const char* HitCalibrationBar::CalibrationStatusAbbreviation[] = { "XX", "FA", "  ", "EC",
                                                                           "TS", "SE", "??", "TJ" };

        HitCalibrationBar::HitCalibrationBar(const Int_t id)
            : ID(id)
            , Validity(0)
            , FailCounter(0)
            , LastEventNumber(0)
            , TimeDifference(NaN)
            , EffectiveSpeed(NaN)
            , Gain({ NaN, NaN })
            , InvLightAttenuationLength(NaN)
            , Pedestal({ 0, 0 })
            , Saturation({ 0., 0. })
            , TimeSync({ NaN, NaN })
            , Threshold({ NaN, NaN })
        {

            Reset();
            LastHits.reserve(CalibrationLogSize);

            Log.TimeDifference.Expand(LogInitialSize);
            Log.TimeDifference.SetTitle("Time Difference Parameter;Event Number;Parameter Value / ns");

            Log.EffectiveSpeed.Expand(LogInitialSize);
            Log.EffectiveSpeed.SetTitle("Effective Speed Parameter;Event Number;Parameter Value / cm/ns");

            Log.LightAttenuationLength.Expand(LogInitialSize);
            Log.LightAttenuationLength.SetTitle("Attenuation Length Parameter;Event Number;Parameter Value / cm");

            for (auto side = 0; side < 2; ++side)
            {
                Log.Gain[side].Expand(LogInitialSize);
                Log.Gain[side].SetTitle("Gain Parameter;Event Number;Parameter Value / ns/MeVle");

                Log.TotalHits[side] = TH1F("", "Total Hits; Light at PMT / MeVle; #Hits", 64, 0., 32.);
                Log.TotalHits[side].SetDirectory(nullptr);
                Log.TotalHits[side].SetStats(false);

                Log.MissedHits[side] = TH1F("", "Missed Hits; Light at PMT / MeVle; #Hits", 64, 0., 32.);
                Log.MissedHits[side].SetDirectory(nullptr);
                Log.MissedHits[side].SetStats(false);
                Log.MissedHits[side].SetLineColor(kGreen);

                Log.Energy[side] = TH1F("", "Energy Calibration; Gain / ns/MeVle; Counts", 100, 0., 30.);
                Log.Energy[side].SetDirectory(nullptr);
                Log.Energy[side].SetStats(false);

                Log.Saturation[side] =
                    TH2F("", "Saturation; Light at PMT / MeVle; QDC / ns", 128, 0., 128., 256, 0., 512.);
                Log.Saturation[side].SetDirectory(nullptr);
                Log.Saturation[side].SetStats(false);

                Log.Pedestal[side] = TH1F("", "Pedestal; QDC / ns", 512, 0., 512.);
                Log.Pedestal[side].SetDirectory(nullptr);
                Log.Pedestal[side].SetStats(false);
            }

            FitGraph.Expand(PositionCalibrationSize);
        }

        HitCalibrationBar::~HitCalibrationBar() {}

        void HitCalibrationBar::Update(const R3BNeulandHitModulePar* par)
        {
            if (!std::isnan(par->GetTimeOffset(1)))
            {
                TimeDifference = par->GetTDiff();
                EffectiveSpeed = par->GetEffectiveSpeed();
                SetStatus(Validity, PosCalibrationBit);
            }

            if (!std::isnan(par->GetEnergyGain(1)))
            {
                Gain[0] = par->GetEnergyGain(1);
                Gain[1] = par->GetEnergyGain(2);

                InvLightAttenuationLength = 1. / par->GetLightAttenuationLength();

                SetStatus(Validity, EnergyCalibrationBit);
            }

            if (!std::isnan(par->GetPMTThreshold(1)))
            {
                Threshold[0] = par->GetPMTThreshold(1);
                Threshold[1] = par->GetPMTThreshold(2);
                SetStatus(Validity, ThresholdCalibrationBit);
            }

            if (!std::isnan(par->GetPMTSaturation(1)))
            {
                Saturation[0] = par->GetPMTSaturation(1);
                Saturation[1] = par->GetPMTSaturation(2);
            }

            if (!std::isnan(par->GetPedestal(1)))
            {
                Pedestal[0] = par->GetPedestal(1);
                Pedestal[1] = par->GetPedestal(2);
                SetStatus(Validity, PedestalCalibrationBit);
            }

            LOG(debug) << "Loaded Parameters for Bar " << ID;
            return;
        }

        void HitCalibrationBar::Set(const Int_t side, const Double_t time, const Int_t qdc)
        {
            SetStatus(Validity, side);

            CurrentHit.Time[side] = time;
            CurrentHit.QDC[side] = std::max(qdc - Pedestal[side], 1);
            Log.Pedestal[side].Fill(qdc);

            if (IsStatus(Validity, LeftSideHitBit) && IsStatus(Validity, RightSideHitBit))
            {
                if (fabs(CurrentHit.Time[1] - CurrentHit.Time[0]) > 0.5 * MaxCalTime)
                {
                    if (CurrentHit.Time[1] < CurrentHit.Time[0])
                        CurrentHit.Time[1] += MaxCalTime;
                    else
                        CurrentHit.Time[0] += MaxCalTime;
                }

                if (IsStatus(Validity, PosCalibrationBit))
                {
                    CurrentHit.EntryPosition =
                        (CurrentHit.Time[1] - CurrentHit.Time[0] - TimeDifference) * EffectiveSpeed;

                    if (fabs(CurrentHit.EntryPosition) > TotalBarLength)
                    {
                        ++FailCounter;
                        if (FailCounter > 10)
                        {
                            LOG(debug) << "Bar " << ID << ":"
                                       << " Positions were too often too strange. Might be a Time-Jump in one PMT. "
                                          "Removing Position-Calibration.";
                            FailCounter = 0;
                            ClearStatus(Validity, PosCalibrationBit);
                        }
                    }
                    else
                    {
                        FailCounter = 0;
                    }
                }
            }
        }
        void HitCalibrationBar::SetGlobalTSync(const Double_t value, const Double_t error)
        {
            TimeSync = { value, error };
            if (std::isnan(value))
                ClearStatus(Validity, TSyncCalibrationBit);
            else
                SetStatus(Validity, TSyncCalibrationBit);
        }

        HitCalibrationBar::CalibrationStatus HitCalibrationBar::GetCalibrationStatus() const
        {
            if (LastHits.size() == 0 && !IsStatus(Validity, LogCompleteBit))
                return CalibrationStatus::noData;
            if (!IsStatus(Validity, PosCalibrationBit))
                return CalibrationStatus::completeFail;
            if (!IsStatus(Validity, EnergyCalibrationBit))
                return CalibrationStatus::energyGainFail;
            if (!IsStatus(Validity, TSyncCalibrationBit))
                return CalibrationStatus::timeSyncFail;
            if (IsStatus(Validity, TimeJumpedBit))
                return CalibrationStatus::timeJump;
            if (TimeSync[1] > MaxTSyncError)
                return CalibrationStatus::timeSyncError;
            if (IsStatus(Validity, StrangeBit))
                return CalibrationStatus::strangeParameters;
            return CalibrationStatus::success;
        }

        Bool_t HitCalibrationBar::Add(const Double_t timeOffset,
                                      const Double_t entryPosition,
                                      const Double_t exitPosition,
                                      const Double_t energy,
                                      const UInt_t eventNumber)
        {
            // If we have already a gain calibration, we can log the hits.
            // This way we might get some information about strange behaviour.
            if (IsStatus(Validity, EnergyCalibrationBit))
            {
                const auto centerPosition = (entryPosition + exitPosition) * 0.5;
                const DPair lightAtPMT = { GetLightAtPMT(energy, centerPosition, InvLightAttenuationLength, RightSide),
                                           GetLightAtPMT(energy, centerPosition, InvLightAttenuationLength, LeftSide) };

                if (!IsStatus(Validity, RightSideHitBit))
                {
                    Log.MissedHits[RightSide].Fill(lightAtPMT[RightSide]);
                }
                if (!IsStatus(Validity, LeftSideHitBit))
                {
                    Log.MissedHits[LeftSide].Fill(lightAtPMT[LeftSide]);
                }

                for (auto side = 0; side < 2; ++side)
                    Log.TotalHits[side].Fill(lightAtPMT[side]);

                if (IsStatus(Validity, LeftSideHitBit) && IsStatus(Validity, RightSideHitBit))
                {
                    for (auto side = 0; side < 2; ++side)
                    {
                        Log.Saturation[side].Fill(lightAtPMT[side], CurrentHit.QDC[side]);
                        const auto gain = CurrentHit.QDC[side] /
                                          (lightAtPMT[side] * (1. - CurrentHit.QDC[side] * SaturationCoefficient));
                        Log.Energy[side].Fill(gain);
                    }
                }
            }

            // We require both sides and reasonable energy deposit.
            if (!IsStatus(Validity, LeftSideHitBit) || !IsStatus(Validity, RightSideHitBit) ||
                energy < MinEnergyDeposit)
                return false;

            // Dont take the bar, if the cosmic was somewhere close to the lightguide:
            // Strange things happen there ;).
            const auto maxPos = BarLength * 0.5 - 1.;
            if (fabs(entryPosition) > maxPos || fabs(exitPosition) > maxPos)
                return false;

            LastEventNumber = eventNumber;

            for (auto side = 0; side < 2; ++side)
                CurrentHit.Time[side] -= timeOffset * 0.5;

            CurrentHit.EntryPosition = entryPosition;
            CurrentHit.ExitPosition = exitPosition;
            CurrentHit.Energy = energy;

            LastHits.push_back(CurrentHit);

            if (LastHits.size() % PositionCalibrationSize == 0)
                positionCalibration(LastHits.size() - PositionCalibrationSize, PositionCalibrationSize);

            if (LastHits.size() % EnergyCalibrationSize == 0)
                energyCalibration(LastHits.size() - EnergyCalibrationSize, EnergyCalibrationSize);

            // If we reached our limit, clear the vector.
            // The content itself is still there. We will set a (hacky) bit,
            // so we know when the vector actually contains more data than its size.
            if (LastHits.size() == CalibrationLogSize)
            {
                SetStatus(Validity, LogCompleteBit);
                LastHits.clear();
            }

            // We can say only that this was a somewhat good event, if we have a position calibration
            // and the position of the track matches with the one from the calibration
            if (IsStatus(Validity, PosCalibrationBit) &&
                fabs(0.5 * (entryPosition + exitPosition) - GetPosition()) > 25.)
                return true;

            return true;
        }

        void HitCalibrationBar::Reset()
        {
            ClearStatus(Validity, LeftSideHitBit);
            ClearStatus(Validity, RightSideHitBit);
            CurrentHit.EntryPosition = NaN;
        }

        void HitCalibrationBar::Calibrate(TDirectory* histogramDir)
        {
            const auto nHits = (IsStatus(Validity, LogCompleteBit) ? CalibrationLogSize : LastHits.size());

            // Check if we have any Hits at all.
            // If not, this bar is either not connected or has some problems.
            // In this case do not bother creating plots and just return.
            if (nHits == 0)
                return;

            // Position Calibration
            if (Log.TimeDifference.GetN() > 0)
            {
                TimeDifference = getMean(Log.TimeDifference);

                EffectiveSpeed = getMean(Log.EffectiveSpeed);

                // check for timejumps
                if (GetJumps(Log.TimeDifference.GetY(), Log.TimeDifference.GetN(), TimeJumpThreshold, 3, nullptr))
                    SetStatus(Validity, TimeJumpedBit);

                SetStatus(Validity, PosCalibrationBit);
            }
            else if (nHits > 0.5 * PositionCalibrationSize)
            {
                // There were not enough hits for the usual calibration.
                // If we have at least half the chunk size, we can try to calibrate anyway
                positionCalibration(0, nHits);
            }
            else
            {
                // We cannot have a reliable calibration.
                // Give up at this point.
                ClearStatus(Validity, PosCalibrationBit);
            }

            // Gain Calibration
            if (Log.Gain[0].GetN() > 0 && Log.Gain[1].GetN() > 0)
            {
                for (auto side = 0; side < 2; ++side)
                {
                    Gain[side] = getMean(Log.Gain[side], Gain[side]);
                    Saturation[side] = SaturationCoefficient * Gain[side];
                }

                InvLightAttenuationLength = 1. / getMean(Log.LightAttenuationLength, 1. / InvLightAttenuationLength);

                SetStatus(Validity, EnergyCalibrationBit);
            }
            else if (nHits > 0.5 * EnergyCalibrationSize)
            {
                // There were not enough hits for the usual calibration.
                // If we have at least half the chunk size, we can try to calibrate anyway
                energyCalibration(0, nHits);
            }
            else
            {
                // We cannot have a reliable calibration.
                // Give up at this point.
                ClearStatus(Validity, EnergyCalibrationBit);
            }

            // Threshold Calibration
            thresholdCalibration();

            createHistograms(histogramDir);
        }

        R3BNeulandHitModulePar HitCalibrationBar::GetParameters()
        {
            R3BNeulandHitModulePar parameter;
            parameter.SetModuleId(ID);

            parameter.SetTDiff(TimeDifference);
            parameter.SetTSync(TimeSync[0]);
            parameter.SetEffectiveSpeed(EffectiveSpeed);
            parameter.SetLightAttenuationLength(1. / InvLightAttenuationLength);

            for (auto side = 0; side < 2; ++side)
            {
                parameter.SetEnergyGain(Gain[side], side + 1);
                parameter.SetPMTThreshold(Threshold[side], side + 1);
                parameter.SetPMTSaturation(Saturation[side], side + 1);
                parameter.SetPedestal(Pedestal[side], side + 1);
            }

            LOG(debug) << ID << " Bar Parameters: ";
            LOG(debug) << "   Time Offset       : " << parameter.GetTimeOffset(1) << " " << parameter.GetTimeOffset(2);
            LOG(debug) << "   Effective Speed   : " << parameter.GetEffectiveSpeed();
            LOG(debug) << "   Energy Gain       : " << parameter.GetEnergyGain(1) << " " << parameter.GetEnergyGain(2);
            LOG(debug) << "   Attenuationlength : " << parameter.GetLightAttenuationLength();
            LOG(debug) << "   Threshold         : " << parameter.GetPMTThreshold(1) << " "
                       << parameter.GetEnergyGain(2);
            LOG(debug) << "   Pedestal          : " << parameter.GetPedestal(1) << " " << parameter.GetPedestal(2);
            LOG(debug) << "   Saturation        : " << parameter.GetPMTSaturation(1) << " "
                       << parameter.GetPMTSaturation(2);

            return parameter;
        }

        Bool_t HitCalibrationBar::IsValid() const
        {
            return IsStatus(Validity, LeftSideHitBit) && IsStatus(Validity, RightSideHitBit);
        }

        void HitCalibrationBar::positionCalibration(int firstHit, int nHits)
        {
            FitGraph.Set(0);
            const auto loopSize = (LastHits.size() >= (firstHit + nHits)) ? nHits : (LastHits.size() - firstHit);
            for (auto index = 0; index < loopSize; ++index)
            {
                const auto& hit = LastHits[firstHit + index];

                const auto meanPosition = 0.5 * (hit.EntryPosition + hit.ExitPosition);
                const auto tdiff = hit.Time[1] - hit.Time[0];
                FitGraph.SetPoint(index, meanPosition, tdiff);
            }

            TF1 linearFit("", "pol1", -TotalBarLength * 0.5, TotalBarLength * 0.5);
            FitGraph.Fit(&linearFit, "NQ");
            if (linearFit.GetParError(0) > MaxFastTDiffError)
                cleanupFit(FitGraph, linearFit, 2.5);

            // Use this calibration
            TimeDifference = linearFit.GetParameter(0);
            EffectiveSpeed = 1. / linearFit.GetParameter(1);

            // Write parameters to the Log.
            const auto nPoints = Log.TimeDifference.GetN();
            Log.TimeDifference.SetPoint(nPoints, LastEventNumber, TimeDifference);
            Log.TimeDifference.SetPointError(nPoints, 0, linearFit.GetParError(0));
            Log.EffectiveSpeed.SetPoint(nPoints, LastEventNumber, EffectiveSpeed);
            Log.EffectiveSpeed.SetPointError(nPoints, 0, linearFit.GetParError(1) * Sqr(EffectiveSpeed));

            SetStatus(Validity, PosCalibrationBit);
        }

        void HitCalibrationBar::energyCalibration(int firstHit, int nHits)
        {
            const auto loopSize = (LastHits.size() >= (firstHit + nHits)) ? nHits : (LastHits.size() - firstHit);
            if (!IsStatus(Validity, PedestalCalibrationBit))
            {
                // This is the first calibration
                // if we do not have pedestals yet, get them and subtract them from the hits

                pedestalCalibration();
                for (auto h_index = 0; h_index < loopSize; ++h_index)
                {
                    auto& hit = LastHits[firstHit + h_index];
                    for (auto side = 0; side < 2; ++side)
                    {
                        hit.QDC[side] = std::max(hit.QDC[side] - Pedestal[side], 1);
                    }
                }
            }

            FitGraph.Set(0);
            for (auto h_index = 0; h_index < loopSize; ++h_index)
            {
                const auto& hit = LastHits[firstHit + h_index];
                const auto centerPosition = (hit.EntryPosition + hit.ExitPosition) * 0.5;
                FitGraph.SetPoint(h_index, centerPosition, log(hit.QDC[1] * 1. / hit.QDC[0]));
            }

            TF1 linearFit("", "pol1", -TotalBarLength * 0.5, TotalBarLength * 0.5);
            // this has usually always outliers, so do not bother doing a normal fit first.
            cleanupFit(FitGraph, linearFit, 1.);

            InvLightAttenuationLength = 0.5 * linearFit.GetParameter(1);

            const auto logLightAttLenPoints = Log.LightAttenuationLength.GetN();
            Log.LightAttenuationLength.SetPoint(logLightAttLenPoints, LastEventNumber, 1. / InvLightAttenuationLength);
            Log.LightAttenuationLength.SetPointError(
                logLightAttLenPoints, 0, 0.5 * linearFit.GetParError(1) / Sqr(InvLightAttenuationLength));

            for (auto h_index = 0; h_index < loopSize; ++h_index)
            {
                const auto& hit = LastHits[firstHit + h_index];
                const auto centerPosition = (hit.EntryPosition + hit.ExitPosition) * 0.5;
                for (auto side = 0; side < 2; ++side)
                {
                    const auto lightAtPMT = GetLightAtPMT(hit.Energy, centerPosition, InvLightAttenuationLength, side);
                    const auto gain = hit.QDC[side] / (lightAtPMT * (1. - hit.QDC[side] * SaturationCoefficient));
                    Log.Energy[side].Fill(gain);
                }
            }

            std::array<TH1F, 2> hGainCal = { TH1F("", "", 200, 0., 40.), TH1F("", "", 200, 0., 40.) };

            for (auto h_index = 0; h_index < loopSize; ++h_index)
            {
                const auto& hit = LastHits[firstHit + h_index];
                const auto centerPosition = (hit.EntryPosition + hit.ExitPosition) * 0.5;
                for (auto side = 0; side < 2; ++side)
                {
                    const auto lightAtPMT = GetLightAtPMT(hit.Energy, centerPosition, InvLightAttenuationLength, side);
                    const auto distance = 0.5 * TotalBarLength + (1 - 2 * side) * centerPosition;
                    const auto gain = hit.QDC[side] / (lightAtPMT * (1. - hit.QDC[side] * SaturationCoefficient));
                    hGainCal[side].Fill(gain);
                }
            }

            TF1 gausFit("", "gaus", 0., 30.);
            for (auto side = 0; side < 2; ++side)
            {
                const auto maxPos = hGainCal[side].GetBinCenter(hGainCal[side].GetMaximumBin());
                gausFit.SetParameter(0, hGainCal[side].GetMaximum());
                gausFit.SetParameter(1, maxPos);
                hGainCal[side].Fit(&gausFit, "NQ", "", maxPos - 1.5, maxPos + 1.5);

                Gain[side] = gausFit.GetParameter(1);

                const auto nPoints = Log.Gain[side].GetN();
                Log.Gain[side].SetPoint(nPoints, LastEventNumber, Gain[side]);
                Log.Gain[side].SetPointError(nPoints, 0, gausFit.GetParError(1));
            }

            SetStatus(Validity, EnergyCalibrationBit);
        }

        void HitCalibrationBar::pedestalCalibration()
        {
            for (auto side = 0; side < 2; ++side)
                for (auto bin = 0; bin < Log.Pedestal[side].GetNbinsX(); ++bin)
                {
                    if (Log.Pedestal[side].GetBinContent(bin + 1) > 3.)
                    {
                        Pedestal[side] = bin;
                        break;
                    }
                }

            SetStatus(Validity, PedestalCalibrationBit);
        }

        void HitCalibrationBar::thresholdCalibration()
        {
            if (Log.TotalHits[0].GetEntries() < ThresholdCalibrationSize ||
                Log.TotalHits[1].GetEntries() < ThresholdCalibrationSize)
                return;

            // We have 50% of the maximum at x=[1]
            // For the bar at the edges it is hard to check if the cosmic was stopped.
            // Therefore we usually have [2] != 0.
            TF1 missFit("missFit", "(1. - [2]) * 0.5 * (1. - TMath::Erf([0] * (x - [1]))) + [2]", 0., 20.);

            for (auto side = 0; side < 2; ++side)
            {
                const auto nPoints = Log.TotalHits[side].GetNbinsX();
                Log.MissRatio[side].Expand(nPoints);
                Log.MissRatio[side].Set(0);
                for (auto bin = 1; bin <= nPoints; ++bin)
                {
                    if (Log.MissedHits[side].GetBinContent(bin) > 0)
                    {
                        const auto nextPoint = Log.MissRatio[side].GetN();
                        const auto binCenter = Log.MissedHits[side].GetBinCenter(bin);
                        const auto ratio =
                            Log.MissedHits[side].GetBinContent(bin) / Log.TotalHits[side].GetBinContent(bin);
                        Log.MissRatio[side].SetPoint(nextPoint, binCenter, ratio);
                    }
                }

                if (Log.MissRatio[side].GetN() > 5)
                {
                    missFit.SetParameter(0, 0.25);
                    missFit.SetParameter(1, AvgThreshold);
                    missFit.SetParLimits(1, 0., 10.);
                    missFit.SetParameter(2, 0.);
                    missFit.SetParLimits(2, 0., 0.5);
                    Log.MissRatio[side].Fit(&missFit, "QNE");
                    Threshold[side] = missFit.GetParameter(1);
                }
            }
            SetStatus(Validity, ThresholdCalibrationBit);
        }

        Int_t HitCalibrationBar::cleanupFit(TGraph& graph, TF1& fit, Double_t maxDifference) const
        {
            std::array<Int_t, 256> remove;
            auto totalRemoved = 0;

            graph.Fit(&fit, "QN ROB=0.90");

            auto removePosition = 0;
            auto x = graph.GetX();
            auto y = graph.GetY();
            auto nPoints = graph.GetN();
            for (auto p = 0; p < nPoints; ++p)
            {
                if (fabs(y[p] - fit.Eval(x[p])) < maxDifference)
                    continue;

                remove[removePosition] = p;
                if (++removePosition == remove.size())
                {
                    removePoints(&remove[0], removePosition, graph);
                    totalRemoved += removePosition;
                    nPoints -= removePosition;
                    p -= removePosition;
                    removePosition = 0;

                    x = graph.GetX();
                    y = graph.GetY();
                }
            }

            removePoints(&remove[0], removePosition, graph);
            totalRemoved += removePosition;

            graph.Fit(&fit, "QN");

            return totalRemoved;
        }

        void HitCalibrationBar::removePoints(Int_t* points, Int_t nPoints, TGraph& graph) const
        {
            if (nPoints == 0)
                return;

            auto offset = 1;
            const auto totalPoints = graph.GetN() - nPoints;
            const auto x = graph.GetX();
            const auto y = graph.GetY();

            for (auto p = points[0]; p < totalPoints; ++p)
            {

                if (offset < nPoints && p + offset == points[offset])
                {
                    ++offset;
                    --p;
                    continue;
                }

                x[p] = x[p + offset];
                y[p] = y[p + offset];
            }

            graph.Set(totalPoints);
        }

        void HitCalibrationBar::createHistograms(TDirectory* histogramDir)
        {
            if (!histogramDir)
                return;

            // First creat all the histograms we do not have yet.
            const auto nHits = (IsStatus(Validity, LogCompleteBit) ? CalibrationLogSize : LastHits.size());

            TGraph posCal;
            posCal.SetTitle("Position Calibration;Position / cm;Time Difference / ns");
            posCal.Expand(nHits);

            const auto loopSize = (LastHits.size() > nHits) ? nHits : LastHits.size();
            for (auto h_index = 0; h_index < loopSize; ++h_index)
            {
                const auto& hit = LastHits[h_index];
                const auto tdiff = hit.Time[1] - hit.Time[0];
                posCal.SetPoint(posCal.GetN(), hit.EntryPosition, tdiff);
            }

            std::array<TF1, 2> fSatCal = { TF1("", "[0]*x/(1+[1]*x)", 0., TotalBarLength),
                                           TF1("", "[0]*x/(1+[1]*x)", 0., TotalBarLength) };

            for (auto side = 0; side < 2; ++side)
            {
                fSatCal[side].SetParameter(0, Gain[side]);
                fSatCal[side].FixParameter(1, Saturation[side]);
            }

            // Go into the directory, create a Canvas, fill it with the plots and write it into the directory.
            histogramDir->cd();
            const auto canvname = TString::Format("Bar %d", ID);
            TCanvas Canvas(canvname, canvname);
            Canvas.SetWindowSize(1920, 1200);
            Canvas.Divide(1, 3);
            auto row = 0;
            TVirtualPad* currentPad;

            currentPad = Canvas.cd(++row);
            currentPad->Divide(3, 1);

            currentPad->cd(1);
            if (posCal.GetN() > 0)
                posCal.Draw("AP");

            currentPad->cd(2);
            if (Log.TimeDifference.GetN() > 0)
                Log.TimeDifference.Draw("A*");

            currentPad->cd(3);
            if (Log.EffectiveSpeed.GetN() > 0)
                Log.EffectiveSpeed.Draw("A*");

            currentPad = Canvas.cd(++row);
            currentPad->Divide(4, 1);

            currentPad->cd(1);
            Log.Energy[0].Draw("");
            Log.Energy[1].SetLineColor(kRed);
            Log.Energy[1].Draw("SAME");

            currentPad->cd(2);
            if (Log.LightAttenuationLength.GetN() > 0)
                Log.LightAttenuationLength.Draw("A*");

            currentPad->cd(3);
            if (Log.Gain[0].GetN() > 0)
                Log.Gain[0].Draw("A*");

            currentPad->cd(4);
            if (Log.Gain[1].GetN() > 0)
                Log.Gain[1].Draw("A*");

            currentPad = Canvas.cd(++row);
            currentPad->Divide(4, 1);

            TPad* mainPads[2];
            TPad* overlayPads[2];
            TGaxis* overlayAxis[2];
            for (auto side = 0; side < 2; ++side)
            {
                // Threshold
                auto pad = currentPad->cd(1 + side);
                mainPads[side] = new TPad("", "", 0, 0, 1, 1);
                mainPads[side]->Draw();
                mainPads[side]->cd();
                Log.TotalHits[side].Draw();
                Log.MissedHits[side].Draw("SAME");
                mainPads[side]->Update();
                pad->cd();

                overlayPads[side] = new TPad("", "", 0, 0, 1, 1);
                overlayPads[side]->SetFillStyle(4000);
                double ymin = 0.;
                double ymax = 1.;
                double dy = (ymax - ymin) / 0.8; // 10 per cent margins top and bottom
                double xmin = mainPads[side]->GetUxmin();
                double xmax = mainPads[side]->GetUxmax();
                double dx = (xmax - xmin) / 0.8; // 10 per cent margins left and right
                overlayPads[side]->Range(xmin - 0.1 * dx, ymin - 0.1 * dy, xmax + 0.1 * dx, ymax + 0.1 * dy);
                overlayPads[side]->Draw();
                overlayPads[side]->cd();

                if (Log.MissRatio[side].GetN() > 0)
                {
                    Log.MissRatio[side].SetLineColor(kRed);
                    Log.MissRatio[side].Draw("IL");
                }
                overlayPads[side]->Update();
                overlayAxis[side] = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 50510, "+L");
                overlayAxis[side]->SetLabelColor(kRed);
                overlayAxis[side]->Draw();

                overlayPads[side]->Update();

                pad->Modified();

                // Saturation

                currentPad->cd(3 + side);
                Log.Saturation[side].Draw("COLZ");
                fSatCal[side].Draw("SAME");
            }

            Canvas.Write();

            for (auto side = 0; side < 2; ++side)
            {
                delete mainPads[side];
                delete overlayPads[side];
                delete overlayAxis[side];
            }
        }

        Double_t HitCalibrationBar::getMean(const TGraphErrors& graph, Double_t expectedValue)
        {
            const auto nPoints = graph.GetN();
            if (nPoints == 0)
                return NaN;

            if (nPoints == 1)
                return graph.GetY()[0];

            TF1 constantFit("", "pol0", 0, nPoints);

            FitGraph.Set(0);
            for (auto point = 0; point < nPoints; ++point)
            {
                FitGraph.SetPoint(point, point, graph.GetY()[point]);
                FitGraph.SetPointError(point, 0, graph.GetEY()[point]);
            }
            constantFit.SetParameter(0, expectedValue);
            FitGraph.Fit(&constantFit, "QN");
            return constantFit.GetParameter(0);
        }
    } // namespace Calibration
} // namespace Neuland
