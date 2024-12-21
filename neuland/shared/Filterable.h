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

#ifndef FILTERABLE
#define FILTERABLE

#include <functional>
#include <vector>

template <typename T>
class Filterable
{
  public:
    using Filter = std::function<bool(const T)>;

  private:
    std::vector<Filter> filters;

  public:
    Filterable() = default;
    explicit Filterable(const std::vector<Filter>& vf)
        : filters(vf)
    {
    }
    inline void Add(const Filter& f) { filters.push_back(f); }
    inline bool IsValid(const T t) const
    {
        for (const auto& filter : filters)
        {
            if (filter(t) == false)
            {
                return false;
            }
        }
        return true;
    }
};

#endif // FILTERABLE
