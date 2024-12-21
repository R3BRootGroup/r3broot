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

#include "R3BDistribution3D.h"
#include "R3BDistribution1D.h"

#include "TGraph.h"

#include <stdexcept>
#include <string>

TGraph createLookupGraph(std::function<Double_t(Double_t)> distribution,
                         const Double_t lowerBound,
                         const Double_t upperBound,
                         const Int_t samples = 1000);

const Int_t Dim = 3;
using Arr = std::array<Double_t, Dim>;
const Int_t nSigma = 5;

R3BDistribution<Dim> R3BDistribution3D::Delta(const Double_t value1, const Double_t value2, const Double_t value3)
{
    return R3BDistribution<Dim>({ value1, value2, value3 });
}

R3BDistribution<Dim> R3BDistribution3D::Constant(const Arr lower_values, const Arr upper_values)
{
    return R3BDistribution<Dim>(lower_values, upper_values);
}

R3BDistribution<Dim> R3BDistribution3D::Gaussian(const Arr means, const Arr sigmas)
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

R3BDistribution<Dim> R3BDistribution3D::Cube(const Arr center, const Double_t edgeLength)
{
    return R3BDistribution<Dim>(
        [center, edgeLength](Arr values) -> Arr
        {
            return { center[0] + edgeLength * (values[0] - 0.5),
                     center[1] + edgeLength * (values[1] - 0.5),
                     center[2] + edgeLength * (values[2] - 0.5) };
        });
}

R3BDistribution<Dim> R3BDistribution3D::Sphere(const Arr center, const Double_t radius)
{
    return R3BDistribution<Dim>(
        [center, radius](Arr values) -> Arr
        {
            const auto r = radius * cbrt(values[0]);
            const auto phi = 2 * TMath::Pi() * values[1];
            const auto theta = TMath::Pi() * values[2];

            const auto x = center[0] + r * sin(theta) * cos(phi);
            const auto y = center[1] + r * sin(theta) * sin(phi);
            const auto z = center[2] + r * cos(theta);

            return { x, y, z };
        });
}

R3BDistribution<Dim> R3BDistribution3D::SphereSurface(const Arr center, const Double_t radius)
{
    return R3BDistribution<Dim>(
        [center, radius](Arr values) -> Arr
        {
            const auto phi = 2 * TMath::Pi() * values[1];
            const auto theta = TMath::Pi() * values[2];

            const auto x = center[0] + radius * sin(theta) * cos(phi);
            const auto y = center[1] + radius * sin(theta) * sin(phi);
            const auto z = center[2] + radius * cos(theta);

            return { x, y, z };
        });
}

R3BDistribution<Dim> R3BDistribution3D::Prism(R3BDistribution<2> xydist, R3BDistribution<1> zdist)
{
    return R3BDistribution<Dim>(
        [xydist, zdist](const Arr values) mutable -> Arr
        {
            Arr ret;
            const auto xy = xydist.GetRandomValues({ values[0], values[1] });
            const auto z = zdist.GetRandomValues({ values[2] })[0];
            return { xy[0], xy[1], z };
        });
}