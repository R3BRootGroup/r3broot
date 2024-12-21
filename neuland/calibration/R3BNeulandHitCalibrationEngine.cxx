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

#include "R3BNeulandHitCalibrationEngine.h"
#include "R3BNeulandHitModulePar.h"
#include "R3BNeulandHitPar.h"

#include "FairLogger.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include <cmath>

using DPair = std::array<Double_t, 2>;
using CalibrationStatus = Neuland::Calibration::HitCalibrationBar::CalibrationStatus;

namespace Neuland
{

    // NaN2Value maps NaNs and Infs to a finite value.
    // Usefull if you want to fill a possible nan variable into a histogram
    inline double NaN2Value(const double val, const double valIfNaN = 0.)
    {
        return (std::isfinite(val) ? val : valIfNaN);
    }

    namespace Calibration
    {
        HitCalibrationEngine::HitCalibrationEngine() {}

        void HitCalibrationEngine::Init(const R3BNeulandHitPar* hitpar)
        {
            const auto nPlanes = (!hitpar ? MaxNumberOfPlanes : hitpar->GetNumberOfPlanes());
            const auto totalBars = nPlanes * BarsPerPlane;
            fHitMask.resize(nPlanes, 0L);

            fBars.reserve(totalBars);
            for (auto plane = 0; plane < nPlanes; ++plane)
                for (auto bar = 0; bar < BarsPerPlane; ++bar)
                    fBars.emplace_back(plane * BarsPerPlane + bar + 1);

            const auto nParameter = hitpar->GetNumModulePar();

            for (auto i = 0; i < nParameter; i++)
            {
                const auto modulePar = hitpar->GetModuleParAt(i);
                const auto barID = modulePar->GetModuleId() - 1;
                fBars[barID].Update(modulePar);
            }

            fBarDistribution = TH1F(
                "HitCalibrationEngine::BarDistribution", "Bar Distribution; Bar ID", totalBars, 0.5, totalBars + 0.5);
            fBarDistribution.SetStats(kFALSE);
            fBarDistribution.SetDirectory(nullptr);

            fStoppedDistribution = TH1F("HitCalibrationEngine::StoppedDistribution",
                                        "Stopped Bar Distribution; Bar ID",
                                        totalBars,
                                        0.5,
                                        totalBars + 0.5);
            fStoppedDistribution.SetStats(kFALSE);
            fStoppedDistribution.SetDirectory(nullptr);

            fInteractionsDistribution = TH1F("HitCalibrationEngine::InteractionsDistribution",
                                             "Interactions Distribution; Number of Interactions",
                                             250,
                                             0.5,
                                             250 + 0.5);
            fInteractionsDistribution.SetStats(kFALSE);
            fInteractionsDistribution.SetDirectory(nullptr);

            fStoppedInteractionsDistribution = TH1F("HitCalibrationEngine::StoppedInteractionsDistribution",
                                                    "Stopped Interactions Distribution; Number of Interactions",
                                                    250,
                                                    0.5,
                                                    250 + 0.5);
            fStoppedInteractionsDistribution.SetStats(kFALSE);
            fStoppedInteractionsDistribution.SetDirectory(nullptr);
            fStoppedInteractionsDistribution.SetLineColor(kRed);

            fTrackLengthDistribution = TH1F("HitCalibrationEngine::TrackLengthDistribution",
                                            "Track Length Distribution; Track Length / cm",
                                            5000,
                                            0.,
                                            100.);
            fTrackLengthDistribution.SetStats(kFALSE);
            fTrackLengthDistribution.SetDirectory(nullptr);

            fTotalTrackLengthDistribution = TH1F("HitCalibrationEngine::TotalTrackLengthDistribution",
                                                 "Total Track Length Distribution; Track Length / cm",
                                                 1000,
                                                 0.,
                                                 500.);
            fTotalTrackLengthDistribution.SetStats(kFALSE);
            fTotalTrackLengthDistribution.SetDirectory(nullptr);

            fTotalStoppedTrackLengthDistribution = TH1F("HitCalibrationEngine::fTotalStoppedTrackLengthDistribution",
                                                        "Total Stopped Track Length Distribution; Track Length / cm",
                                                        1000,
                                                        0.,
                                                        500.);
            fTotalStoppedTrackLengthDistribution.SetStats(kFALSE);
            fTotalStoppedTrackLengthDistribution.SetDirectory(nullptr);
            fTotalStoppedTrackLengthDistribution.SetLineColor(kRed);

            fCorrelationMatrix = TH2F("HitCalibrationEngine::CorrelationMatrix",
                                      "Correlation Matrix;Bar ID;Bar ID",
                                      totalBars,
                                      0.5,
                                      0.5 + totalBars,
                                      totalBars,
                                      0.5,
                                      0.5 + totalBars);
            fCorrelationMatrix.SetStats(kFALSE);
            fCorrelationMatrix.SetDirectory(nullptr);

            fTrackEntryPointDistribution = TH3F("HitCalibrationEngine::TrackEntryPointDistribution",
                                                "Entry Point Distribution; Z / cm; X / cm; Y / cm",
                                                150,
                                                0,
                                                300,
                                                130,
                                                -130,
                                                130,
                                                130,
                                                -130,
                                                130);
            fTrackEntryPointDistribution.SetStats(kFALSE);
            fTrackEntryPointDistribution.SetDirectory(nullptr);

            fTrackDirectionDistribution = TH3F("HitCalibrationEngine::TrackDirectionDistribution",
                                               "Track Direction Distribution; Z / cm/ns; X / cm/ns; Y / cm/ns",
                                               100,
                                               -50,
                                               50,
                                               100,
                                               -50,
                                               50,
                                               100,
                                               -50,
                                               50);
            fTrackDirectionDistribution.SetStats(kFALSE);
            fTrackDirectionDistribution.SetDirectory(nullptr);
        }

        void HitCalibrationEngine::Set(const Int_t id, const Int_t side, const Double_t time, const Int_t qdc)
        {
            const auto plane = id / BarsPerPlane;
            const auto bar = id % BarsPerPlane;

            fBars[id].Set(side, time, qdc);

            fHitMask[plane] |= 1UL << bar;
        }

        void HitCalibrationEngine::Add(const R3BNeulandCosmicTrack& track, const UInt_t eventNumber)
        {
            fTotalTrackLengthDistribution.Fill(track.TotalTrackLength);
            if (track.Stopped)
                fTotalStoppedTrackLengthDistribution.Fill(track.TotalTrackLength);

            fTrackEntryPointDistribution.Fill(track.EntryPoint.Z(), track.EntryPoint.X(), track.EntryPoint.Y());
            fTrackDirectionDistribution.Fill(track.Direction.Z(), track.Direction.X(), track.Direction.Y());

            const auto nPlanes = fHitMask.size();

            fInteractionsDistribution.Fill(track.Interactions.size());
            if (track.Stopped)
            {
                fStoppedDistribution.Fill(track.Interactions.back().BarID + 1);
                fStoppedInteractionsDistribution.Fill(track.Interactions.size());
            }

            // first add correlations
            for (auto plane = 0; plane < nPlanes; ++plane)
            {
                if (!fHitMask[plane])
                    continue;

                for (auto bar = 0; bar < BarsPerPlane; ++bar)
                {
                    if (!(fHitMask[plane] & (1UL << bar)))
                        continue;

                    for (auto otherPlane = 0; otherPlane < nPlanes; ++otherPlane)
                    {
                        if (!fHitMask[otherPlane])
                            continue;

                        for (auto otherBar = 0; otherBar < BarsPerPlane; ++otherBar)
                        {
                            if ((fHitMask[otherPlane] & (1UL << otherBar)) && (plane != otherPlane || bar != otherBar))
                                fCorrelationMatrix.Fill(BarsPerPlane * plane + bar + 1,
                                                        BarsPerPlane * otherPlane + otherBar + 1);
                        }
                    }
                }
            }

            // if the track stopped in NeuLAND, do not take the last interaction
            const auto nInteractions = track.Interactions.size() - (track.Stopped ? 1 : 0);
            for (auto i = 0; i < nInteractions; ++i)
            {
                const auto& interaction = track.Interactions[i];
                const auto barID = interaction.BarID;
                const auto bar = barID % BarsPerPlane;
                const auto plane = GetPlaneNumber(barID);

                fBarDistribution.Fill(barID + 1);
                fTrackLengthDistribution.Fill(interaction.TrackLength);

                if (fBars[barID].Add(interaction.EntryTime,
                                     interaction.EntryPosition,
                                     interaction.ExitPosition,
                                     interaction.Energy,
                                     eventNumber))
                {
                    fTSyncer.AddBarData(barID, fBars[barID].GetTime());
                }

                fBars[barID].Reset();
                fHitMask[plane] &= ~(1UL << bar);
            }

            fTSyncer.DoEvent();

            for (auto plane = 0; plane < nPlanes; ++plane)
            {
                if (fHitMask[plane] == 0UL)
                    continue;

                for (auto bar = 0; bar < BarsPerPlane; ++bar)
                {
                    if (fHitMask[plane] & (1UL << bar))
                        fBars[plane * BarsPerPlane + bar].Reset();
                }
                fHitMask[plane] = 0UL;
            }
        }

        void HitCalibrationEngine::Reset()
        {
            for (auto p = 0; p < fHitMask.size(); ++p)
            {
                if (!fHitMask[p])
                    continue;

                for (auto b = 0; b < BarsPerPlane; ++b)
                {
                    if (fHitMask[p] & (1UL << b))
                        fBars[p * BarsPerPlane + b].Reset();
                }
                fHitMask[p] = 0UL;
            }
        }

        std::vector<R3BNeulandHitModulePar> HitCalibrationEngine::Calibrate(TDirectory* histoDir)
        {
            const auto nPlanes = fHitMask.size();
            const auto nBars = nPlanes * BarsPerPlane;
            std::vector<R3BNeulandHitModulePar> allParameters;
            if (histoDir)
            {
                histoDir->cd();
                TCanvas canvasTracks("Tracks", "Tracks");
                canvasTracks.SetWindowSize(1920, 1200);
                canvasTracks.Divide(1, 2);
                auto pad = canvasTracks.cd(1);
                pad->Divide(5, 1);

                pad->cd(1);
                fBarDistribution.Draw();
                pad->cd(2);
                fStoppedDistribution.Draw();
                pad->cd(3);
                fInteractionsDistribution.Draw();
                fStoppedInteractionsDistribution.Draw("SAME");
                pad->cd(4);
                fTotalTrackLengthDistribution.Draw();
                fTotalStoppedTrackLengthDistribution.Draw("SAME");
                pad->cd(5);
                fTrackLengthDistribution.Draw();

                pad = canvasTracks.cd(2);
                pad->Divide(3, 1);
                pad->cd(1);
                fCorrelationMatrix.Draw("COLZ");
                pad->cd(2);
                fTrackEntryPointDistribution.Draw();
                pad->cd(3);
                fTrackDirectionDistribution.Draw();
                canvasTracks.Write("Tracking");
            }

            TDirectory* planeDir = nullptr;
            for (Int_t plane = 0; plane < nPlanes; ++plane)
            {
                if (histoDir)
                    planeDir = histoDir->mkdir(TString::Format("Plane_%d", plane + 1));

                for (Int_t bar = 0; bar < BarsPerPlane; ++bar)
                {
                    std::cout << "Calibrating Bar " << BarsPerPlane * plane + bar << "/" << BarsPerPlane * nPlanes
                              << "\r" << std::flush;
                    fBars[BarsPerPlane * plane + bar].Calibrate(planeDir);
                }
            }
            if (histoDir)
                histoDir->cd();

            std::cout << "Syncing NeuLAND Bars...                                 \r" << std::flush;
            const auto tsync = fTSyncer.GetTSync(nPlanes);

            for (Int_t plane = 0; plane < nPlanes; ++plane)
            {
                for (Int_t bar = 0; bar < BarsPerPlane; ++bar)
                {
                    const auto barID = plane * BarsPerPlane + bar;
                    fBars[barID].SetGlobalTSync(tsync[barID].Value, tsync[barID].Error);

                    const auto calibStatus = fBars[barID].GetCalibrationStatus();

                    if (calibStatus != CalibrationStatus::completeFail && calibStatus != CalibrationStatus::noData)
                        allParameters.push_back(fBars[barID].GetParameters());
                }
            }

            TH1F hTDiff("R3BNeulandCal2HitPar::TDiff", "TDiff;BarID;TDiff / ns", nBars, 0.5, 0.5 + nBars);
            hTDiff.SetStats(false);
            TH1F hTSync("R3BNeulandCal2HitPar::TSync", "TSync;BarID;TSync / ns", nBars, 0.5, 0.5 + nBars);
            hTSync.SetStats(false);
            TH1F hVEff("R3BNeulandCal2HitPar::VEff", "VEff;BarID;VEff / cm/ns", nBars, 0.5, 0.5 + nBars);
            hVEff.SetStats(false);
            TH1F hAtt("R3BNeulandCal2HitPar::Att",
                      "Attenuation Length;BarID;Attenuation Length / cm",
                      nBars,
                      0.5,
                      0.5 + nBars);
            hAtt.SetStats(false);
            TH1F hPosErr(
                "R3BNeulandCal2HitPar::PosError", "Position Error;BarID;\\Delta Pos / cm", nBars, 0.5, 0.5 + nBars);
            hPosErr.SetStats(false);
            TH1F hEPosErr("R3BNeulandCal2HitPar::EnergyPosError",
                          "Energy Position Error;BarID;\\Delta Pos / cm",
                          nBars,
                          0.5,
                          0.5 + nBars);
            hEPosErr.SetStats(false);
            std::array<TH1F, 2> hGain = {
                TH1F("R3BNeulandCal2HitPar::Gain1", "Gain;BarID;Gain / ns/MeV", nBars, 0.5, 0.5 + nBars),

                TH1F("R3BNeulandCal2HitPar::Gain2", "Gain;BarID;Gain / ns/MeV", nBars, 0.5, 0.5 + nBars)
            };
            hGain[0].SetLineColor(kRed);
            hGain[1].SetLineColor(kGreen);
            hGain[0].SetStats(false);
            hGain[1].SetStats(false);

            std::array<TH1F, 2> hSat = {
                TH1F("R3BNeulandCal2HitPar::Sat1", "Saturation;BarID;Saturation / 1/MeV", nBars, 0.5, 0.5 + nBars),
                TH1F("R3BNeulandCal2HitPar::Sat2", "Saturation;BarID;Saturation / 1/MeV", nBars, 0.5, 0.5 + nBars)
            };
            hSat[0].SetLineColor(kRed);
            hSat[1].SetLineColor(kGreen);
            hSat[0].SetStats(false);
            hSat[1].SetStats(false);

            std::array<TH1F, 2> hThr = {
                TH1F("R3BNeulandCal2HitPar::Thr1", "Threshold;BarID;Threshold / MeV", nBars, 0.5, 0.5 + nBars),
                TH1F("R3BNeulandCal2HitPar::Thr2", "Threshold;BarID;Threshold / MeV", nBars, 0.5, 0.5 + nBars)
            };
            hThr[0].SetLineColor(kRed);
            hThr[1].SetLineColor(kGreen);
            hThr[0].SetStats(false);
            hThr[1].SetStats(false);

            for (auto& parameter : allParameters)
            {
                const auto barID = parameter.GetModuleId();
                const auto attLength = NaN2Value(parameter.GetLightAttenuationLength());

                hTDiff.Fill(barID, parameter.GetTDiff());
                hTSync.Fill(barID, NaN2Value(parameter.GetTSync()));
                hVEff.Fill(barID, NaN2Value(parameter.GetEffectiveSpeed()));
                hAtt.Fill(barID, attLength);

                for (Int_t side = 0; side < 2; ++side)
                {
                    hGain[side].Fill(barID, NaN2Value(parameter.GetEnergyGain(side + 1)));
                    hSat[side].Fill(barID, NaN2Value(parameter.GetPMTSaturation(side + 1)));
                    hThr[side].Fill(barID, NaN2Value(parameter.GetPMTThreshold(side + 1)));
                }
            }

            draw();

            if (histoDir)
            {
                TCanvas canvasOverview("Overview", "Neuland Hit Calibration");
                canvasOverview.Divide(3, 2);

                canvasOverview.cd(1);
                hTDiff.Draw("HIST");
                canvasOverview.cd(2);
                hVEff.Draw("HIST");
                canvasOverview.cd(3);
                hTSync.Draw("HIST");
                canvasOverview.cd(4);
                hGain[0].Draw("HIST");
                hGain[1].Draw("HIST SAME");
                canvasOverview.cd(5);
                hAtt.Draw("HIST");
                canvasOverview.cd(6);
                hThr[0].Draw("HIST");
                hThr[1].Draw("HIST SAME");
                canvasOverview.Write("Overview");
            }
            return allParameters;
        }

        void HitCalibrationEngine::draw() const
        {
            const auto nPlanes = fHitMask.size();
            int stats = 0;

            std::vector<CalibrationStatus> statuses;
            statuses.reserve(16);

            std::cout << std::endl << std::endl << "   |";
            for (auto plane = 0; plane < nPlanes; ++plane)
            {
                std::cout << TString::Format("%2d|", plane + 1);
            }
            std::cout << std::endl << "   ";
            for (auto plane = 0; plane < nPlanes; ++plane)
            {
                std::cout << " ꓕꓕ";
            }
            std::cout << std::endl;
            for (auto bar = BarsPerPlane - 1; bar >= 0; --bar)
            {
                std::cout << TString::Format("%2d |", bar + 1);
                for (auto plane = 0; plane < nPlanes; ++plane)
                {
                    const auto cstatus = fBars[plane * BarsPerPlane + bar].GetCalibrationStatus();

                    if (cstatus != CalibrationStatus::success &&
                        std::find(statuses.begin(), statuses.end(), cstatus) == statuses.end())
                    {
                        statuses.push_back(cstatus);
                    }
                    std::cout << HitCalibrationBar::GetCalibrationStatusAbbreviation(cstatus);
                    std::cout << "|";
                }
                std::cout << std::endl;
            }
            std::cout << "   ";
            for (auto plane = 0; plane < nPlanes; ++plane)
            {
                std::cout << " TT";
            }
            std::cout << std::endl << "   |";
            for (auto plane = 0; plane < nPlanes; ++plane)
            {
                std::cout << TString::Format("%2d|", plane + 1);
            }
            std::cout << std::endl << std::endl;
            if (statuses.size() > 0)
            {
                for (const auto status : statuses)
                {
                    std::cout << "   '" << HitCalibrationBar::GetCalibrationStatusAbbreviation(status)
                              << "' = " << HitCalibrationBar::GetCalibrationStatusDescription(status) << std::endl;
                }
                std::cout << std::endl;
            }
        }

    } // namespace Calibration
} // namespace Neuland