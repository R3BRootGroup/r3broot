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

#ifndef NEULANDCLUSTERINGENGINEH
#define NEULANDCLUSTERINGENGINEH

#include <algorithm>
#include <functional>
#include <vector>

namespace Neuland
{

    template <typename T>
    class ClusteringEngine
    {
        using Tit = typename std::vector<T>::iterator;
        using BinaryPredicate = std::function<bool(const T&, const T&)>;

      private:
        BinaryPredicate f;

        /* Partitions an already divided range by clustering both the items in the first part AND the items added to the
         * first part during this process with the items in the second part.
         * Takes the begin of a range of interest that acts as reference/seed, the end of that range which is the
         * divider to the range to be analyzed, the end of the container, and the clustering function.
         * Returns the new divider between cluster and non-cluster */
        Tit moving_partition(const Tit begin, Tit moving_divider, const Tit end) const
        {
            // The end of the iteration will change if matching objects are found, such that these will be iterated over
            // as well
            for (Tit a = begin; a != moving_divider; a++)
            {
                // Partition returns iterator to the first element of the non-clustered rest of the container
                moving_divider = std::partition(moving_divider, end, [&](const T& b) { return f(*a, b); });
            }
            return moving_divider;
        }

      public:
        /* Default Constructor. Note: If the clustering condition is not set, a "bad_function_call" will be thrown upon
         * calling clusterize. This seems better than providing a default function which might produce hard-to-track
         * unwanted results. */
        ClusteringEngine(){};
        ClusteringEngine(const BinaryPredicate& _f)
            : f(_f)
        {
        }

        void SetClusteringCondition(const BinaryPredicate& _f) { f = _f; }

        bool SatisfiesClusteringCondition(const T& a, const T& b) const { return f(a, b); }

        std::vector<std::vector<T>> Clusterize(std::vector<T>& from) const
        {
            std::vector<std::vector<T>> out;

            /* Three iterators (read: markers for positions in the input vector) are required:
             * Begin is the start of the clustered part of the vector
             * Divider separates the clustered part from the unclustered part, and
             * end is the end of the vector*/
            Tit begin;
            Tit divider = from.begin();
            const Tit end = from.end();

            /* While the end of the input vector has not been reached, start a new cluster containing one element from
             * the unclustered part and then move all objects that match the clustering condition in that range.
             * This part of the input vector can then be move-constructed into a new vector = cluster */
            while (divider != end)
            {
                begin = divider;
                divider++;
                divider = moving_partition(begin, divider, end);
                out.push_back(
                    std::move(std::vector<T>(std::make_move_iterator(begin), std::make_move_iterator(divider))));
            }
            return out;
        }
    };

}; // namespace Neuland

#endif // NEULANDCLUSTERINGENGINEH
