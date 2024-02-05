#pragma once
/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// TODO: use C++20 std::source_location
#include <boost/assert/source_location.hpp>

#include <FairRootManager.h>
#include <R3BException.h>
#include <fmt/format.h>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

// TODO: Use C++20 Concept to put more constains on the template parameters
namespace R3B
{
    template <typename InputType>
    class InputConnector
    {
      public:
        using RawDataType = std::remove_const_t<std::remove_cv_t<InputType>>;
        explicit InputConnector(std::string_view branchName)
            : branch_name_{ branchName }
        {
        }

        // rule of 5
        ~InputConnector() = default;
        InputConnector(const InputConnector&) = default;
        InputConnector(InputConnector&&) = delete;
        InputConnector& operator=(const InputConnector&) = default;
        InputConnector& operator=(InputConnector&&) = delete;

        void init(const boost::source_location& loc = BOOST_CURRENT_LOCATION)
        {
            if (auto* ioman = FairRootManager::Instance(); ioman != nullptr)
            {
                data_ = ioman->InitObjectAs<const RawDataType*>(branch_name_.c_str());

                if (data_ == nullptr)
                {
                    throw R3B::runtime_error(
                        fmt::format("Initialisation of the input data with the branch name \"{}\" failed!",
                                    branch_name_),
                        loc);
                }
            }
            else
            {
                throw R3B::runtime_error(fmt::format("FairRootManager is nullptr during the initialisation of the "
                                                     "input data with  the branch name \"{}\"",
                                                     branch_name_),
                                         loc);
            }
        }

        [[nodiscard]] inline auto get(const boost::source_location& loc = BOOST_CURRENT_LOCATION) const
            -> const RawDataType&
        {
            check_init(loc);
            return *data_;
        }

        auto size(const boost::source_location& loc = BOOST_CURRENT_LOCATION) const
        {
            check_init(loc);
            return data_->size();
        }

        // implement range-based for loop:
        // TODO: is there a simpler way?
        auto begin(const boost::source_location& loc = BOOST_CURRENT_LOCATION)
        {
            check_init(loc);
            return data_->cbegin();
        }
        auto begin(const boost::source_location& loc = BOOST_CURRENT_LOCATION) const
        {
            check_init(loc);
            return data_->cbegin();
        }
        auto end(const boost::source_location& loc = BOOST_CURRENT_LOCATION)
        {
            check_init(loc);
            return data_->cend();
        }
        auto end(const boost::source_location& loc = BOOST_CURRENT_LOCATION) const
        {
            check_init(loc);
            return data_->cend();
        }

      private:
        std::string branch_name_;
        const RawDataType* data_ = nullptr;
        void check_init(const boost::source_location& loc) const
        {
            if (data_ == nullptr)
            {
                throw R3B::runtime_error(
                    fmt::format("Input data with the branch name \"{}\" cannot be queried without an initialisation!",
                                branch_name_),
                    loc);
            }
        }
    };

    template <typename OutputType>
    class OutputConnector
    {
      public:
        using RawDataType = std::remove_const_t<std::remove_cv_t<OutputType>>;
        explicit OutputConnector(std::string_view branchName)
            : branch_name_{ branchName }
        {
        }

        // rule of 5
        ~OutputConnector() = default;
        OutputConnector(const OutputConnector&) = delete;
        OutputConnector(OutputConnector&&) = delete;
        OutputConnector& operator=(const OutputConnector& other) = delete;
        OutputConnector& operator=(OutputConnector&&) = delete;

        void init(bool persistance = true, const boost::source_location& loc = BOOST_CURRENT_LOCATION)
        {
            if (auto* ioman = FairRootManager::Instance(); ioman != nullptr)
            {
                ioman->RegisterAny(branch_name_.c_str(), data_ptr_, persistance);
            }
            else
            {
                throw R3B::runtime_error(fmt::format("FairRootManager is nullptr during the initialisation of the "
                                                     "output data with the branch name \"{}\"",
                                                     branch_name_),
                                         loc);
            }
        }

        [[nodiscard]] inline auto get() -> RawDataType& { return data_; }
        [[nodiscard]] inline auto get_constref() const -> const RawDataType& { return data_; }

        inline void clear() { data_.clear(); }

        template <typename ResetOp>
        inline void clear(ResetOp&& opn)
        {
            opn(data_);
        }

      private:
        std::string branch_name_;
        RawDataType data_;
        RawDataType* data_ptr_ = &data_;
    };

    template <typename ElementType>
    using InputVectorConnector = InputConnector<std::vector<ElementType>>;

    template <typename ElementType>
    using OutputVectorConnector = OutputConnector<std::vector<ElementType>>;

    template <typename KeyType, typename ValueType>
    using InputMapConnector = InputConnector<std::map<KeyType, ValueType>>;

    template <typename KeyType, typename ValueType>
    using OutputMapConnector = OutputConnector<std::map<KeyType, ValueType>>;

    template <typename KeyType, typename ValueType>
    using InputHashConnector = InputConnector<std::unordered_map<KeyType, ValueType>>;

    template <typename KeyType, typename ValueType>
    using OutputHashConnector = OutputConnector<std::unordered_map<KeyType, ValueType>>;
} // namespace R3B
