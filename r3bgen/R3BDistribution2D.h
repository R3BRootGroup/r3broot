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

#ifndef R3BDISTRIBUTION2D_H
#define R3BDISTRIBUTION2D_H

#include "TF2.h"
#include "TGraph2D.h"
#include "TH2.h"

#include "R3BDistribution.h"

class R3BDistribution2D : public R3BDistribution<2>
{
    using Arr = std::array<Double_t, 2>;

  public:
    static R3BDistribution<2> Delta(const Double_t value1, const Double_t value2) { return Delta({ value1, value2 }); }
    static R3BDistribution<2> Delta(const Arr values);
    static R3BDistribution<2> Flat(const Arr lower_values, const Arr upper_values);
    static R3BDistribution<2> Gaussian(const Double_t mean, const Double_t sigma);
    static R3BDistribution<2> Gaussian(const Arr means, const Arr sigmas);
    static R3BDistribution<2> Square(const Arr center, const Double_t edgeLength);
    static R3BDistribution<2> Rectangle(const Arr center, const Arr edgeLengths);
    static R3BDistribution<2> Circle(const Arr center, const Double_t radius);
    // static R3BDistribution<2> Data(const TH2& data); todo
    // static R3BDistribution<2> Data(const TGraph2D& data); todo
    // static R3BDistribution<2> Data(const TF2& data); todo
};
#endif
