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

#ifndef R3BDISTRIBUTION1D_H
#define R3BDISTRIBUTION1D_H

#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"

#include "R3BDistribution.h"

class R3BDistribution1D : public R3BDistribution<1>
{
  public:
    static R3BDistribution<1> Delta(const Double_t value);
    static R3BDistribution<1> Flat(const Double_t lower_value, const Double_t upper_value);
    static R3BDistribution<1> Gaussian(const Double_t mean, const Double_t sigma);
    static R3BDistribution<1> Function(const std::function<Double_t(const Double_t)> func,
                                       const Double_t lower_bound,
                                       const Double_t upper_bound);
    static R3BDistribution<1> Data(const TF1& data);
    static R3BDistribution<1> Data(const TF1& data, const Double_t lower_bound, const Double_t upper_bound);
    static R3BDistribution<1> Data(const TH1& data);
    static R3BDistribution<1> Data(const TH1& data, const Double_t lower_bound, const Double_t upper_bound);
    static R3BDistribution<1> Data(const TGraph& data);
    static R3BDistribution<1> Data(const TGraph& data, const Double_t lower_bound, const Double_t upper_bound);
    static R3BDistribution<1> DataLogLog(const TGraph& data);
    static R3BDistribution<1> DataLogLog(const TGraph& data, const Double_t lower_bound, const Double_t upper_bound);
};
#endif
