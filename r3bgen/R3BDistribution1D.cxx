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

#include "R3BDistribution1D.h"

#include <stdexcept>
#include <string>

const Int_t Dim = 1;
using Arr = std::array<Double_t, Dim>;
const Int_t nSigma = 5;

TGraph createLookupGraph(std::function<Double_t(Double_t)> distribution,
                         const Double_t lower_bound,
                         const Double_t upper_bound,
                         const Int_t samples = 1000)
{
    const Double_t step = (upper_bound - lower_bound) / samples;

    TGraph steppingGraph;
    steppingGraph.Expand(samples + 2);
    steppingGraph.SetPoint(steppingGraph.GetN(), lower_bound, 0);
    for (int i = 0; i <= samples; ++i)
    {
        const auto x = lower_bound + i * step;
        const auto y = distribution(x);
        if (y < 0.)
            throw std::underflow_error("Found negative value inside data!");

        steppingGraph.SetPoint(steppingGraph.GetN(), x, y);
    }
    steppingGraph.SetPoint(steppingGraph.GetN(), upper_bound, 0);

    const auto invInt = 1. / steppingGraph.Integral();

    TGraph integralGraph;
    integralGraph.Expand(samples);
    for (int i = 1; i < steppingGraph.GetN() - 2; ++i)
    {
        Double_t x, y;
        steppingGraph.GetPoint(i + 1, x, y);
        steppingGraph.SetPoint(i + 1, x - step, 0.); // For Integral Method

        integralGraph.SetPoint(integralGraph.GetN(), x - step, steppingGraph.Integral(0, i + 1) * invInt);
        steppingGraph.SetPoint(i + 1, x, y); // Reset
    }
    integralGraph.SetPoint(integralGraph.GetN(), upper_bound, 1);

    // inverse Graph
    auto x = integralGraph.GetX();
    auto y = integralGraph.GetY();
    for (int i = 0; i < integralGraph.GetN(); ++i)
        integralGraph.SetPoint(i, y[i], x[i]);

    return integralGraph;
}

TGraph createLookupGraph(const TH1& distribution, const Double_t lower_bound, const Double_t upper_bound)
{
    TGraph integralGraph;
    integralGraph.Expand(distribution.GetNbinsX() + 1);
    Double_t integral = 0;
    integralGraph.SetPoint(integralGraph.GetN(), lower_bound, 0);
    Int_t startbin = distribution.GetXaxis()->FindBin(lower_bound),
          endbin = distribution.GetXaxis()->FindBin(upper_bound);

    if (distribution.GetBinContent(startbin) < 0 || distribution.GetBinContent(endbin) < 0)
        throw std::underflow_error("Found negative value inside data!");

    // handle first bin --> might be "splitted"
    integral += (distribution.GetXaxis()->GetBinUpEdge(startbin) - lower_bound) /
                distribution.GetXaxis()->GetBinWidth(startbin) * distribution.GetBinContent(startbin);
    integralGraph.SetPoint(integralGraph.GetN(), distribution.GetXaxis()->GetBinUpEdge(startbin), integral);

    for (int i = startbin + 1; i < endbin; ++i)
    {
        if (distribution.GetBinContent(i) < 0)
            throw std::underflow_error("Found negative value inside data!");

        integral += distribution.GetBinContent(i);
        integralGraph.SetPoint(integralGraph.GetN(), distribution.GetXaxis()->GetBinUpEdge(i), integral);
    }

    // handle last bin --> might be "splitted"
    integral += (upper_bound - distribution.GetXaxis()->GetBinLowEdge(endbin)) /
                distribution.GetXaxis()->GetBinWidth(endbin) * distribution.GetBinContent(endbin);
    integralGraph.SetPoint(integralGraph.GetN(), upper_bound, integral);

    TGraph targetGraph;
    Double_t x, y, invIntegral = 1 / integral;
    for (int i = 0; i < integralGraph.GetN(); ++i)
    {
        integralGraph.GetPoint(i, x, y);
        targetGraph.SetPoint(targetGraph.GetN(), y * invIntegral, x);
    }

    return targetGraph;
}

TGraph createLogLogLookupGraph(const TGraph& distribution,
                               const Double_t lower_bound,
                               const Double_t upper_bound,
                               const Int_t samples = 1000)
{
    const auto nPoints = distribution.GetN();
    const auto logLowerBound = log(lower_bound);
    const auto logUpperBound = log(upper_bound);
    const auto xValues = distribution.GetX();
    const auto yValues = distribution.GetY();
    const auto xLogStep = log(upper_bound / lower_bound) / samples;

    for (auto p = 0; p < nPoints - 1; ++p)
        if (*(xValues + p) > *(xValues + p + 1))
            throw std::range_error(std::string(__func__) + " : X Values are not sorted! Run sort method on TGraph.");

    TGraph logGraph(nPoints);
    logGraph.SetBit(TGraph::kIsSortedX);

    for (auto p = 0; p < nPoints; ++p)
        logGraph.SetPoint(p, log(*(xValues + p)), log(*(yValues + p)));

    TGraph integralGraph;
    integralGraph.SetBit(TGraph::kIsSortedX);

    Double_t integral = 0;

    integralGraph.SetPoint(integralGraph.GetN(), lower_bound, 0);

    for (auto p = 0; p < samples; ++p)
    {
        const auto x0 = logLowerBound + xLogStep * p;
        const auto x1 = logLowerBound + xLogStep * (p + 1);
        const auto y0 = logGraph.Eval(x0, 0, "S");
        const auto y1 = logGraph.Eval(x1, 0, "S");

        const auto slope = (y1 - y0) / (x1 - x0);
        const auto offset = exp(y0 - x0 * slope);
        if (slope == -1.) // Just to be sure
            integral += offset * (x1 - x0);
        else
            integral += offset * (pow(exp(x1), slope + 1.) - pow(exp(x0), slope + 1.)) / (slope + 1);

        integralGraph.SetPoint(integralGraph.GetN(), exp(x1), integral);
    }

    const auto intXValues = integralGraph.GetX();
    const auto intYValues = integralGraph.GetY();
    const auto gPoints = integralGraph.GetN();
    const auto invIntegral = 1 / integral;

    TGraph targetGraph;
    for (auto p = 0; p < gPoints; ++p)
        targetGraph.SetPoint(p, (*(intYValues + p)) * invIntegral, *(intXValues + p));

    targetGraph.SetBit(TGraph::kIsSortedX);

    return targetGraph;
}

R3BDistribution<Dim> R3BDistribution1D::Delta(const Double_t value) { return R3BDistribution<Dim>({ value }); }

R3BDistribution<Dim> R3BDistribution1D::Flat(const Double_t lower_value, const Double_t upper_value)
{
    return R3BDistribution<Dim>({ lower_value }, { upper_value });
}

R3BDistribution<Dim> R3BDistribution1D::Gaussian(const Double_t mean, const Double_t sigma)
{
    const auto invSigma = 1 / sigma;
    auto g = createLookupGraph([mean, invSigma](const Double_t value)
                               { return TMath::Exp(-0.5 * (value - mean) * (value - mean) * invSigma); },
                               mean - nSigma * sigma,
                               mean + nSigma * sigma);
    return R3BDistribution<Dim>([g](Arr values) -> Arr { return { g.Eval(values[0]) }; });
}

R3BDistribution<Dim> R3BDistribution1D::Function(const std::function<Double_t(const Double_t)> func,
                                                 const Double_t lower_bound,
                                                 const Double_t upper_bound)
{
    auto g = createLookupGraph([func](Double_t val) { return func(val); }, lower_bound, upper_bound);
    return R3BDistribution<1>([g](Arr values) -> Arr { return { g.Eval(values[0]) }; });
}

R3BDistribution<Dim> R3BDistribution1D::Data(const TF1& data)
{
    Double_t lower_bound, upper_bound;
    data.GetRange(lower_bound, upper_bound);
    return R3BDistribution1D::Data(data, lower_bound, upper_bound);
}

R3BDistribution<Dim> R3BDistribution1D::Data(const TF1& data, const Double_t lower_bound, const Double_t upper_bound)
{
    if (lower_bound < data.GetXmin() || upper_bound > data.GetXmax())
        throw std::range_error(std::string(__func__) + " : bounds outside of data-range");

    auto g = createLookupGraph([&data](Double_t val) { return data.Eval(val); }, lower_bound, upper_bound);
    return R3BDistribution<1>([g](Arr values) -> Arr { return { g.Eval(values[0]) }; });
}

R3BDistribution<Dim> R3BDistribution1D::Data(const TH1& data)
{
    return R3BDistribution1D::Data(
        data, data.GetXaxis()->GetBinLowEdge(1), data.GetXaxis()->GetBinUpEdge(data.GetNbinsX()));
}

R3BDistribution<Dim> R3BDistribution1D::Data(const TH1& data, const Double_t lower_bound, const Double_t upper_bound)
{
    if (lower_bound < data.GetXaxis()->GetXmin() || upper_bound > data.GetXaxis()->GetXmax())
        throw std::range_error(std::string(__func__) + " : bounds outsie of data-range");

    auto g = createLookupGraph(data, lower_bound, upper_bound);
    return R3BDistribution<1>([g](Arr values) -> Arr { return { g.Eval(values[0]) }; });
}

R3BDistribution<Dim> R3BDistribution1D::Data(const TGraph& data)
{
    Double_t lower_bound, upper_bound, ymin, ymax;
    data.ComputeRange(lower_bound, ymin, upper_bound, ymax);
    return R3BDistribution1D::Data(data, lower_bound, upper_bound);
}

R3BDistribution<Dim> R3BDistribution1D::Data(const TGraph& data, const Double_t lower_bound, const Double_t upper_bound)
{
    Double_t xmin, xmax, ymin, ymax;
    data.ComputeRange(xmin, ymin, xmax, ymax);
    if (lower_bound < xmin || upper_bound > xmax)
        throw std::range_error(std::string(__func__) + " : bounds outsie of data-range");

    auto g = createLookupGraph([&data](Double_t val) { return data.Eval(val); }, lower_bound, upper_bound);
    return R3BDistribution<1>([g](Arr values) -> Arr { return { g.Eval(values[0]) }; });
}

R3BDistribution<Dim> R3BDistribution1D::DataLogLog(const TGraph& data)
{
    Double_t lower_bound, upper_bound, ymin, ymax;
    data.ComputeRange(lower_bound, ymin, upper_bound, ymax);
    return R3BDistribution1D::DataLogLog(data, lower_bound, upper_bound);
}

R3BDistribution<Dim> R3BDistribution1D::DataLogLog(const TGraph& data,
                                                   const Double_t lower_bound,
                                                   const Double_t upper_bound)
{
    Double_t xmin, xmax, ymin, ymax;
    data.ComputeRange(xmin, ymin, xmax, ymax);
    if (lower_bound < xmin || upper_bound > xmax)
        throw std::range_error(std::string(__func__) + " : bounds outsie of data-range");

    auto g = createLogLogLookupGraph(data, lower_bound, upper_bound);
    return R3BDistribution<1>([g](Arr values) -> Arr { return { g.Eval(values[0]) }; });
}
