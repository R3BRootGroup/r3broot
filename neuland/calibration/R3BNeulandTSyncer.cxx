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

#include "R3BNeulandTSyncer.h"
#include "LSQR.h"

#include "FairLogger.h"

#include "TF1.h"

#include <numeric>

constexpr auto NextBarLogSize = 128;
constexpr auto NextPlaneLogSize = 64;
constexpr auto NBins = 256;

using namespace Neuland;

namespace Neuland
{
    namespace Calibration
    {
        TSyncer::TSyncer()
            : SamplingHistogram("", "", MaxCalTime / 2, -MaxCalTime, MaxCalTime)
        {
            SamplingHistogram.SetDirectory(nullptr);
            for (auto& m : HitMask)
                m = 0UL;

            for (auto& bar : Data)
            {
                bar.NextBarLog.reserve(NextBarLogSize);
                for (auto b = 0; b < BarsPerPlane; ++b)
                    bar.NextPlaneLog[b].reserve(NextPlaneLogSize);
            }
        }

        void TSyncer::AddBarData(const Int_t barID, const Double_t time)
        {
            const auto bar = barID % BarsPerPlane;
            const auto plane = GetPlaneNumber(barID);
            HitMask[plane] |= 1UL << bar;
            EventData[barID] = time;
        }

        void TSyncer::ClearBarData(const Int_t barID)
        {
            Data[barID].NextBar.Clear();
            Data[barID].NextBarLog.clear();
            for (auto b = 0; b < BarsPerPlane; ++b)
            {
                Data[barID].NextPlane[b].Clear();
                Data[barID].NextPlaneLog[b].clear();
            }

            const auto bar = barID % BarsPerPlane;
            const auto plane = GetPlaneNumber(barID);

            if (bar != 0)
            {
                // This is not the first bar in the plane
                // Clear the Data from the bar before
                Data[barID - 1].NextBar.Clear();
                Data[barID - 1].NextBarLog.clear();
            }

            if (plane != 0)
            {
                // This is not the first plane
                // Clear the Data of this bar from the previous plane
                for (auto b = (plane - 1) * BarsPerPlane; b < plane * BarsPerPlane; ++b)
                {
                    Data[b].NextPlane[bar].Clear();
                    Data[b].NextPlaneLog[bar].clear();
                }
            }
        }

        void TSyncer::DoEvent()
        {
            for (auto plane = 0; plane < MaxNumberOfPlanes; ++plane)
            {
                if (HitMask[plane] == 0UL)
                    continue;

                for (auto bar = 0; bar < BarsPerPlane; ++bar)
                {
                    if ((HitMask[plane] & (1UL << bar)) == 0UL)
                        continue;

                    const auto barID = plane * BarsPerPlane + bar;
                    auto& barData = Data[barID];

                    if ((bar < BarsPerPlane - 1) && (HitMask[plane] & (1UL << (bar + 1))))
                    {
                        auto tdiff = EventData[barID + 1] - EventData[barID];
                        if (fabs(tdiff) > 0.5 * MaxCalTime)
                        {
                            if (tdiff < 0)
                                tdiff += MaxCalTime;
                            else
                                tdiff -= MaxCalTime;
                        }

                        if (barData.NextBarLog.size() < NextBarLogSize)
                        {
                            barData.NextBarLog.push_back(tdiff);
                            if (barData.NextBarLog.size() == NextBarLogSize)
                            {
                                SamplingHistogram.Reset();
                                for (auto t : barData.NextBarLog)
                                    SamplingHistogram.Fill(t);
                                const auto max = SamplingHistogram.GetBinCenter(SamplingHistogram.GetMaximumBin());
                                barData.NextBar = TH1F("", "", NBins, max - 10, max + 10);
                                barData.NextBar.SetDirectory(nullptr);
                                for (auto t : barData.NextBarLog)
                                    barData.NextBar.Fill(t);
                            }
                        }
                        else
                        {
                            barData.NextBar.Fill(tdiff);
                        }
                    }

                    const auto nextPlane = plane + 1;
                    if (nextPlane >= MaxNumberOfPlanes || !HitMask[nextPlane])
                        continue;

                    for (auto barInNextPlane = 0; barInNextPlane < BarsPerPlane; ++barInNextPlane)
                    {
                        if ((HitMask[nextPlane] & (1UL << barInNextPlane)) == 0UL)
                            continue;

                        const auto nextPlaneBarID = nextPlane * BarsPerPlane + barInNextPlane;
                        auto tdiff = EventData[nextPlaneBarID] - EventData[barID];
                        if (fabs(tdiff) > 0.5 * MaxCalTime)
                        {
                            if (tdiff < 0)
                                tdiff += MaxCalTime;
                            else
                                tdiff -= MaxCalTime;
                        }

                        if (barData.NextPlaneLog[barInNextPlane].size() < NextPlaneLogSize)
                        {
                            barData.NextPlaneLog[barInNextPlane].push_back(tdiff);
                            if (barData.NextPlaneLog[barInNextPlane].size() == NextPlaneLogSize)
                            {
                                SamplingHistogram.Reset();
                                for (auto t : barData.NextPlaneLog[barInNextPlane])
                                    SamplingHistogram.Fill(t);
                                const auto max = SamplingHistogram.GetBinCenter(SamplingHistogram.GetMaximumBin());
                                barData.NextPlane[barInNextPlane] = TH1F("", "", NBins, max - 25, max + 25);
                                barData.NextPlane[barInNextPlane].SetDirectory(nullptr);
                                for (auto t : barData.NextPlaneLog[barInNextPlane])
                                    barData.NextPlane[barInNextPlane].Fill(t);
                            }
                        }
                        else
                        {
                            barData.NextPlane[barInNextPlane].Fill(tdiff);
                        }
                    }
                }
                HitMask[plane] = 0UL;
            }
        }

        std::vector<TSyncer::ValueErrorPair> TSyncer::GetTSync(UInt_t nPlanes)
        {
            calcTSyncs();
            const auto nBars = nPlanes * BarsPerPlane;
            std::vector<ValueErrorPair> solution(nBars, { NaN, NaN });
            std::vector<Double_t> scaleFactors(nBars, 0.);

            // First get the number of Equations and the scaling Factor.
            // A bit slower but less Memory hungry
            UInt_t numberOfEquations = 0;
            for (auto b = 0; b < nBars; ++b)
            {
                if (!std::isnan(Data[b].TSyncNextBar.Value))
                {
                    scaleFactors[b] += 1. / Sqr(Data[b].TSyncNextBar.Error);
                    scaleFactors[b + 1] += 1. / Sqr(Data[b].TSyncNextBar.Error);
                    ++numberOfEquations;
                }

                const auto plane = GetPlaneNumber(b);
                if (plane == nPlanes - 1)
                    continue;

                for (auto ob = 0; ob < BarsPerPlane; ++ob)
                {
                    if (!std::isnan(Data[b].TSyncNextPlane[ob].Value))
                    {
                        scaleFactors[b] += 1. / Sqr(Data[b].TSyncNextPlane[ob].Error);
                        scaleFactors[plane * BarsPerPlane + ob] += 1. / Sqr(Data[b].TSyncNextPlane[ob].Error);
                        ++numberOfEquations;
                    }
                }
            }

            // we have less Equations than bars,
            // seems like we do not have enough statistics in most bars
            if (numberOfEquations < nBars)
            {
                LOG(info) << "Can not synchronize NeuLAND. Not enough equations (" << numberOfEquations << ").";
                return solution;
            }

            for (auto b = 0; b < nBars; ++b)
                scaleFactors[b] = 1. / sqrt(scaleFactors[b]);

            // Now setup for LSQR
            // --- Legacy Code Warning ---

            LSQR_INPUTS* input;
            LSQR_OUTPUTS* output;
            LSQR_WORK* work;

            struct Element
            {
                std::array<UInt_t, 2> Index;
                std::array<Double_t, 2> Value;
            };

            std::vector<Element> lhs(numberOfEquations);
            const auto numberOfVariables = BarsPerPlane * nPlanes;
            auto function =
                [&lhs, numberOfVariables](long mode, LSQR_DOUBLE_VECTOR* xVec, LSQR_DOUBLE_VECTOR* yVec, void* data)
            {
                double* x = xVec->elements;
                double* y = yVec->elements;

                if (mode == 0)
                {
                    for (uint r = 0; r < lhs.size(); ++r)
                    {
                        *y += x[lhs[r].Index[0]] * lhs[r].Value[0] + x[lhs[r].Index[1]] * lhs[r].Value[1];

                        ++y;
                    }
                }
                else
                {
                    for (uint r = 0; r < lhs.size(); ++r)
                    {
                        x[lhs[r].Index[0]] += (*y) * lhs[r].Value[0];
                        x[lhs[r].Index[1]] += (*y) * lhs[r].Value[1];

                        ++y;
                    }
                }

                x = xVec->elements;
                y = yVec->elements + lhs.size();

                if (mode == 0)
                {
                    for (auto i = 0; i < numberOfVariables; i++)
                        *y += x[i];
                }
                else
                {
                    for (auto i = 0; i < numberOfVariables; i++)
                        x[i] += *y;
                }
            };

            alloc_lsqr_mem(&input, &output, &work, numberOfEquations + 1, nBars);

            input->num_rows = numberOfEquations;
            input->num_cols = nBars;
            input->damp_val = 0.;         // we want damping (i.e. in this case average of all solution vars = 0)
            input->rel_mat_err = 1.0e-10; // TODO: this should be set to something reasonable
            input->rel_rhs_err = 1.0e-10; // TODO: this should be set to something reasonable
            input->cond_lim = 0.;         // 10.0 * act_mat_cond_num;
            input->max_iter = input->num_rows + input->num_cols + 50;
            input->lsqr_fp_out = nullptr;

            auto nEQ = 0;
            for (UInt_t id = 0; id < Data.size(); ++id)
            {
                const auto plane = GetPlaneNumber(id);
                if (!std::isnan(Data[id].TSyncNextBar.Value))
                {
                    const auto weight = 1. / Data[id].TSyncNextBar.Error;
                    lhs[nEQ] = { { id, id + 1U }, { -weight * scaleFactors[id], weight * scaleFactors[id + 1U] } };
                    input->rhs_vec->elements[nEQ] = Data[id].TSyncNextBar.Value * weight;
                    ++nEQ;
                }
                for (UInt_t barInNextPlane = 0; barInNextPlane < BarsPerPlane; ++barInNextPlane)
                {
                    if (!std::isnan(Data[id].TSyncNextPlane[barInNextPlane].Value))
                    {
                        const auto barInNextPlaneID = BarsPerPlane * (plane + 1) + barInNextPlane;
                        const auto weight = 1. / Data[id].TSyncNextPlane[barInNextPlane].Error;
                        lhs[nEQ] = { { id, barInNextPlaneID },
                                     { -weight * scaleFactors[id], weight * scaleFactors[barInNextPlaneID] } };
                        input->rhs_vec->elements[nEQ] = Data[id].TSyncNextPlane[barInNextPlane].Value * weight;
                        ++nEQ;
                    }
                }
            }

            for (auto bar = 0; bar < nBars; ++bar)
                input->sol_vec->elements[bar] = 0.;

            input->rhs_vec->elements[numberOfEquations] = 0.; // we will use this one the put the mean value to 0

            LOG(debug) << "Syncing Neuland with " << numberOfEquations << " equations...";

            lsqr(input, output, work, function, &lhs);

            for (auto id = 0; id < nBars; ++id)
            {
                solution[id] = { output->sol_vec->elements[id] * scaleFactors[id],
                                 output->std_err_vec->elements[id] * scaleFactors[id] };
            }
            free_lsqr_mem(input, output, work);

            return solution;
        }

        void TSyncer::calcTSyncs()
        {
            for (auto& bar : Data)
            {
                bar.TSyncNextBar = { NaN, NaN };
                if (bar.NextBar.Integral() > bar.NextBar.GetEntries() * 0.5 &&
                    bar.NextBar.Integral() > 0.5 * NextBarLogSize)
                {
                    bar.TSyncNextBar = { bar.NextBar.GetMean(), bar.NextBar.GetStdDev() };
                }

                for (auto nb = 0; nb < BarsPerPlane; ++nb)
                {
                    bar.TSyncNextPlane[nb] = { NaN, NaN };
                    if (bar.NextPlane[nb].Integral() > bar.NextPlane[nb].GetEntries() * 0.5 &&
                        bar.NextPlane[nb].Integral() > 0.5 * NextPlaneLogSize)
                    {
                        bar.TSyncNextPlane[nb] = { bar.NextPlane[nb].GetMean(), bar.NextPlane[nb].GetStdDev() };
                    }
                }
            }
        }
    } // namespace Calibration
} // namespace Neuland
