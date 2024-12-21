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

#ifndef R3BDISTRIBUTION_H
#define R3BDISTRIBUTION_H

#include <array>
#include <cmath>
#include <functional>
#include <memory>

#include "R3BDouble.h"

template <Int_t dimension>
class R3BDistribution
{
    template <Int_t dim>
    friend class R3BDistribution;

    using Array = std::array<Double_t, dimension>;

  public:
    R3BDistribution()
        : fLookupFunction([](const Array vals) { return vals; })
    {
        CreateValues();
    } // Required by Root Dict

    R3BDistribution(std::array<R3BDistribution<1>, dimension> dists)
    {
        fLookupFunction = [dists](const Array rndvalues) mutable
        {
            Array retArray;
            for (int i = 0; i < dimension; ++i)
                retArray[i] = dists[i].GetRandomValues({ rndvalues[i] })[0];
            return retArray;
        };

        for (int i = 0; i < dimension; ++i)
            fValues[i] = dists[i].fValues[0];
    }
    R3BDistribution(const std::function<Array(const Array)> lookupFunction)
        : fLookupFunction(lookupFunction)
    {
        CreateValues();
    }
    R3BDistribution(const Array& values)
        : fLookupFunction([values](const Array) { return values; })
    {
        CreateValues();
    }
    R3BDistribution(const Array& lvalues, const Array& uvalues)
        : fLookupFunction(
              [lvalues, uvalues](const Array rndvalues)
              {
                  Array v;
                  for (int i = 0; i < dimension; ++i)
                      v[i] = rndvalues[i] * (uvalues[i] - lvalues[i]) + lvalues[i];
                  return v;
              })
    {
        CreateValues();
    }

    inline operator Array() const
    {
        Array v;
        for (int i = 0; i < dimension; ++i)
            v[i] = fValues[i];
        return v;
    }

    Array GetRandomValues(const Array rnd)
    {
        auto v = fLookupFunction(rnd);
        Array retArr;
        for (int i = 0; i < dimension; ++i)
        {
            *fValues[i] = v[i];
            retArr[i] = v[i];
        }
        return retArr;
    }

    std::array<R3BDouble*, dimension> GetValueAddresses() const
    {
        std::array<R3BDouble*, dimension> retArr;
        for (int i = 0; i < dimension; ++i)
            retArr[i] = fValues[i].get();
        return retArr;
    }

  private:
    std::function<Array(const Array)> fLookupFunction;
    std::array<std::shared_ptr<R3BDouble>, dimension> fValues;

    void CreateValues()
    {
        for (int i = 0; i < dimension; ++i)
            fValues[i] = std::shared_ptr<R3BDouble>(new R3BDouble());
    }
};

#endif
