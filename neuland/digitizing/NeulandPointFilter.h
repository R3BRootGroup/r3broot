#pragma once
#include "R3BNeulandPoint.h"
#include <bitset>

namespace R3B::Neuland
{
    enum class BitSetParticle : uint32_t
    {
        none = 0x0000,
        proton = 0x0001,
        neutron = 0x0002,
        electron = 0x0004,
        positron = 0x0008,
        alpha = 0x0010,
        gamma = 0x0020,
        meson = 0x0040,
        other = 0x40000000
    };

    const std::unordered_map<int, BitSetParticle> PidToBitSetParticleHash = {
        { 2212, BitSetParticle::proton },  { 2112, BitSetParticle::neutron },     { 11, BitSetParticle::electron },
        { -11, BitSetParticle::positron }, { 1000040020, BitSetParticle::alpha }, { 22, BitSetParticle::gamma },
        { 0, BitSetParticle::none }
    };

    constexpr auto ParticleBitsetSize = 32U;

    using ParticleUType = std::underlying_type_t<BitSetParticle>;

    constexpr auto ParticleToBitSet(BitSetParticle particle)
    {
        return std::bitset<ParticleBitsetSize>{ static_cast<ParticleUType>(particle) };
    }

    inline auto BitSetToParticle(std::bitset<ParticleBitsetSize> bits) -> BitSetParticle
    {
        return static_cast<BitSetParticle>(static_cast<uint32_t>(bits.to_ulong()));
    }

    inline auto CheckCriteria(BitSetParticle particle, BitSetParticle criteria) -> bool
    {
        return (ParticleToBitSet(particle) & ParticleToBitSet(criteria)) == ParticleToBitSet(particle);
    }

    inline auto operator|(BitSetParticle left, BitSetParticle right) -> BitSetParticle
    {
        auto left_bitset = ParticleToBitSet(left);
        auto right_bitset = ParticleToBitSet(right);
        return BitSetToParticle(left_bitset | right_bitset);
    }

    inline auto operator&(BitSetParticle left, BitSetParticle right) -> BitSetParticle
    {
        auto left_bitset = ParticleToBitSet(left);
        auto right_bitset = ParticleToBitSet(right);
        return BitSetToParticle(left_bitset & right_bitset);
    }

    inline auto operator~(BitSetParticle particle) -> BitSetParticle
    {
        return BitSetToParticle(~ParticleToBitSet(particle));
    }

    inline auto PidToBitSetParticle(int pid) -> BitSetParticle
    {
        if (pid > 99 and pid < 1000) // NOLINT mesons have three digit pdgs
        {
            return BitSetParticle::meson;
        }
        auto pid_to_bitset_hash_iterator = PidToBitSetParticleHash.find(pid);

        if (pid_to_bitset_hash_iterator == PidToBitSetParticleHash.end())
        {
            return BitSetParticle::other;
        }

        return pid_to_bitset_hash_iterator->second;
    }

} // namespace R3B::Neuland
class NeulandPointFilter
{
  public:
    NeulandPointFilter() = default;
    void SetFilter(R3B::Neuland::BitSetParticle filtered_particles);
    void SetFilter(R3B::Neuland::BitSetParticle filtered_particles, double minimum_allowed_energy);
    [[nodiscard]] auto GetFilter() const -> R3B::Neuland::BitSetParticle { return filtered_particles_; }
    [[nodiscard]] auto GetMinimumAllowedEnergy() const -> double { return minimum_allowed_energy_; }
    auto CheckFiltered(const R3BNeulandPoint& neuland_point) -> bool;

  private:
    R3B::Neuland::BitSetParticle filtered_particles_ = R3B::Neuland::BitSetParticle::none;
    double minimum_allowed_energy_ = 0; // engergy in GeV
};
