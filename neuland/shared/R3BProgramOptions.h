#pragma once

#include <boost/program_options.hpp>
#include <fmt/format.h>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <utility>

namespace R3B
{
    namespace po = boost::program_options;

    // TODO: C++20 use concepts
    template <typename>
    struct is_bool : std::false_type
    {
    };

    template <>
    struct is_bool<bool> : std::true_type
    {
    };

    template <typename T, typename = void>
    struct IsString : std::false_type
    {
    };

    template <typename T>
    struct IsString<T, std::void_t<decltype(std::string_view{ std::declval<T>() })>> : std::true_type
    {
    };

    template <typename T>
    constexpr bool IsString_v = IsString<T>::value;

    class OptionConcept;
    template <typename Type>
    class Option;

    template <typename Type>
    class OptionHandle;

    class ProgramOptions
    {
      public:
        ProgramOptions() = default;
        explicit ProgramOptions(const std::string& desc)
            : desc_{ desc }
        {
        }

        template <typename OptionType>
        auto create_option(const std::string& optionName,
                           const std::string& option_desc,
                           std::optional<OptionType> defaultValue = {}) -> decltype(auto)
        {
            if (auto search = registries_.find(optionName); search != registries_.end())
            {
                std::cerr << "ERROR: option has been already defined!" << std::endl;
                exit(1);
            }
            auto option = std::make_unique<Option<OptionType>>(optionName, defaultValue.value_or(OptionType{}), this);
            option->add(option_desc, not defaultValue.has_value());
            registries_.emplace(optionName, option.get());
            return OptionHandle{ std::move(option) };
        }

        auto verify(int argc, const char** argv) -> bool;
        void delete_option(const std::string& optionName) { registries_.erase(optionName); }
        auto get_posDescRef() -> auto& { return pos_desc_; }
        auto get_desc_ref() -> auto& { return desc_; }

      private:
        std::unordered_map<std::string, OptionConcept*> registries_;
        po::positional_options_description pos_desc_;
        po::options_description desc_;
        po::variables_map varMap_;
    };

    class OptionConcept
    {
      public:
        OptionConcept(const OptionConcept&) = delete;
        OptionConcept(OptionConcept&&) = delete;
        auto operator=(const OptionConcept&) -> OptionConcept& = delete;
        auto operator=(OptionConcept&&) -> OptionConcept& = delete;
        virtual ~OptionConcept() = default;
        OptionConcept() = default;
        virtual void retrieve(const po::variables_map& varMap) = 0;
    };

    template <typename Type>
    class Option : public OptionConcept
    {
      public:
        using type = Type;
        Option(const Option&) = delete;
        Option(Option&&) = delete;
        auto operator=(const Option&) -> Option& = delete;
        auto operator=(Option&&) -> Option& = delete;
        Option(std::string name, Type defaultValue, ProgramOptions* program)
            : name_{ std::move(name) }
            , value_{ std::move(defaultValue) }
            , program_{ program }
        {
            if (auto end = name_.find(','))
            {
                key_ = name_.substr(0, end);
            }
            else
            {
                key_ = name_;
            }
        }
        ~Option() override { program_->delete_option(name_); }

        void add(const std::string& desc, bool is_requried = false)
        {
            auto& po_desc = program_->get_desc_ref();
            is_required_ = is_requried;

            auto desc_full = fmt::format("{} [ required ]", desc);
            if (not is_required_)
            {
                desc_full = fmt::format("{} [ = {} ]", desc, value_);
                if constexpr (IsString_v<Type>)
                {
                    desc_full = fmt::format("{} [ = {:?} ]", desc, value_);
                }
            }

            if constexpr (is_bool<Type>::value)
            {
                po_desc.add_options()(name_.c_str(), desc_full.c_str());
            }
            else
            {
                po_desc.add_options()(name_.c_str(), po::value<Type>(), desc_full.c_str());
            }
        }

        void as_positional(int option)
        {
            program_->get_posDescRef().add(key_.c_str(), option);
            is_positional_ = true;
        }

        void retrieve(const po::variables_map& varMap) override
        {
            if (varMap.count(key_) != 0U)
            {

                if constexpr (is_bool<Type>::value)
                {
                    value_ = true;
                }
                else
                {
                    value_ = varMap[key_].template as<Type>();
                }
            }
            else if (is_required_)
            {
                const auto error_msg = fmt::format(R"(Program option "--{}" is required! )", name_);
                throw std::runtime_error(error_msg);
            }
        }

        void set_required(bool p_rq = true) { is_required_ = p_rq; }

        [[nodiscard]] auto value() const { return value_; }
        [[nodiscard]] auto is_required() const -> bool { return is_required_; }
        [[nodiscard]] auto is_positional() const -> bool { return is_positional_; }

      private:
        bool is_positional_ = false;
        bool is_required_ = false;
        std::string name_;
        std::string key_;
        std::string desc_;
        Type value_{};
        ProgramOptions* program_;
    };

    template <typename Type>
    class OptionHandle
    {
      public:
        explicit OptionHandle(std::unique_ptr<Option<Type>> option)
            : option_{ std::move(option) }
        {
        }

        auto value() const { return option_->value(); }
        auto operator()() const { return option_->value(); }
        auto make_positional(int pos) -> OptionHandle<Type>&&
        {
            option_->as_positional(pos);
            return std::move(*this);
        }
        auto operator->() -> Option<Type>* { return option_.get(); }

      private:
        std::unique_ptr<Option<Type>> option_ = nullptr;
    };

}; // namespace R3B
