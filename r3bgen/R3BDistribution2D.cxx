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

#include "R3BDistribution2D.h"

#include "TGraph.h"

#include <stdexcept>
#include <string>

const Int_t Dim = 2;
using Arr = std::array<Double_t, Dim>;
const Int_t nSigma = 5;

TGraph createLookupGraph(std::function<Double_t(Double_t)> distribution,
                         const Double_t lowerBound,
                         const Double_t upperBound,
                         const Int_t samples = 1000);

R3BDistribution<Dim> R3BDistribution2D::Delta(const Arr values) { return R3BDistribution<Dim>(values); }

R3BDistribution<Dim> R3BDistribution2D::Flat(const Arr lower_Values, const Arr upper_Values)
{
    return R3BDistribution<Dim>(lower_Values, upper_Values);
}

R3BDistribution<Dim> R3BDistribution2D::Gaussian(const Double_t mean, const Double_t sigma)
{
    const auto invSigma2 = 1 / (sigma * sigma);
    auto g = createLookupGraph([mean, invSigma2](const Double_t value)
                               { return TMath::Exp(-0.5 * (value - mean) * (value - mean) * invSigma2); },
                               mean - nSigma * sigma,
                               mean + nSigma * sigma);
    return R3BDistribution<Dim>(
        [g](Arr values) -> Arr
        {
            auto r = g.Eval(values[0]);
            auto phi = 2 * TMath::Pi() * values[1];
            return { r * TMath::Cos(phi), r * TMath::Sin(phi) };
        });
}

R3BDistribution<Dim> R3BDistribution2D::Gaussian(const Arr means, const Arr sigmas)
{
    std::array<TGraph, Dim> graphs;
    for (int i = 0; i < Dim; ++i)
    {
        const auto invSigma2 = 1 / (sigmas[i] * sigmas[i]), mean = means[i];
        graphs[i] = createLookupGraph([mean, invSigma2](const Double_t value)
                                      { return TMath::Exp(-0.5 * (value - mean) * (value - mean) * invSigma2); },
                                      means[0] - nSigma * sigmas[0],
                                      means[0] + nSigma * sigmas[0]);
    }

    return R3BDistribution<Dim>(
        [graphs](Arr values) -> Arr
        {
            Arr ret;
            for (int i = 0; i < Dim; ++i)
                ret[i] = graphs[i].Eval(values[i]);
            return ret;
        });
}

R3BDistribution<Dim> R3BDistribution2D::Square(const Arr center, const Double_t edgeLength)
{
    return R3BDistribution<Dim>(
        [center, edgeLength](const Arr values) -> Arr {
            return { center[0] + edgeLength * (values[0] - 0.5), center[1] + edgeLength * (values[1] - 0.5) };
        });
}

R3BDistribution<Dim> R3BDistribution2D::Rectangle(const Arr center, const Arr edgeLengths)
{
    return R3BDistribution<Dim>(
        [center, edgeLengths](const Arr values) -> Arr {
            return { center[0] + edgeLengths[0] * (values[0] - 0.5), center[1] + edgeLengths[1] * (values[1] - 0.5) };
        });
}

R3BDistribution<Dim> R3BDistribution2D::Circle(const Arr center, const Double_t radius)
{
    return R3BDistribution<Dim>(
        [center, radius](const Arr values) -> Arr
        {
            const auto r = radius * sqrt(values[0]);
            const auto phi = 2 * TMath::Pi() * values[1];

            const auto x = center[0] + r * cos(phi);
            const auto y = center[1] + r * sin(phi);
            return { x, y };
        });
}
/*
R3BDistribution<Dim> R3BDistribution2D::Data(const TH2& data){
    //TODO
}*/
/*
R3BDistribution<Dim> R3BDistribution2D::Data(const TGraph2D& data){
    //TODO
}*/
/*
R3BDistribution<Dim> R3BDistribution2D::Data(const TF2& data){
    //TODO
}*/