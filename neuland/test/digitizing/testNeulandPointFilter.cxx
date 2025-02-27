/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "NeulandPointFilter.h"
#include "R3BNeulandPoint.h"
#include <bitset>
#include <gtest/gtest.h>

namespace
{

    namespace neuland = ::R3B::Neuland;

    using BitSetParticle = R3B::Neuland::BitSetParticle;
    constexpr auto ParticleBitsetSize = R3B::Neuland::ParticleBitsetSize;
    // Test for ParticleToBitSet function

    TEST(TestNeulandPointFilter, ParticleToBitSetTest)
    {
        auto bitset = neuland::ParticleToBitSet(BitSetParticle::proton);
        EXPECT_EQ(bitset.to_ulong(), 0x0001);

        bitset = ParticleToBitSet(BitSetParticle::neutron);
        EXPECT_EQ(bitset.to_ulong(), 0x0002);

        bitset = ParticleToBitSet(BitSetParticle::gamma);
        EXPECT_EQ(bitset.to_ulong(), 0x0020);
    }

    // Test for BitSetToParticle function
    TEST(TestNeulandPointFilter, BitSetToParticleTest)
    {
        constexpr auto proton_bitsetid = 0x0001;
        auto particle = neuland::BitSetToParticle(std::bitset<ParticleBitsetSize>{ proton_bitsetid });
        EXPECT_EQ(particle, BitSetParticle::proton);

        constexpr auto neutron_bitsetid = 0x0002;
        particle = neuland::BitSetToParticle(std::bitset<ParticleBitsetSize>{ neutron_bitsetid });
        EXPECT_EQ(particle, BitSetParticle::neutron);

        constexpr auto gamma_bitsetid = 0x0020;
        particle = neuland::BitSetToParticle(std::bitset<ParticleBitsetSize>{ gamma_bitsetid });
        EXPECT_EQ(particle, BitSetParticle::gamma);
    }

    // Test for CheckCriteria function
    TEST(TestNeulandPointFilter, CheckCriteriaTest)
    {
        EXPECT_TRUE(CheckCriteria(BitSetParticle::proton, BitSetParticle::proton));
        EXPECT_FALSE(CheckCriteria(BitSetParticle::proton, BitSetParticle::neutron));
        EXPECT_TRUE(CheckCriteria(BitSetParticle::proton, BitSetParticle::neutron | BitSetParticle::proton));
    }

    // Test for PidToBitSetParticle function
    TEST(TestNeulandPointFilter, PidToBitSetParticleTest)
    {
        EXPECT_EQ(neuland::PidToBitSetParticle(2212), BitSetParticle::proton);
        EXPECT_EQ(neuland::PidToBitSetParticle(2112), BitSetParticle::neutron);
        EXPECT_EQ(neuland::PidToBitSetParticle(22), BitSetParticle::gamma);
    }

    // Test for NeulandPointFilter class
    TEST(TestNeulandPointFilter, SetFilterTest)
    {
        NeulandPointFilter filter;
        filter.SetFilter(BitSetParticle::proton);
        EXPECT_EQ(filter.GetFilter(), BitSetParticle::proton);

        filter.SetFilter(BitSetParticle::neutron, 0.5);
        EXPECT_EQ(filter.GetFilter(), BitSetParticle::neutron);
        EXPECT_EQ(filter.GetMinimumAllowedEnergy(), 0.5);
    }

    // Test for NeulandPointFilter::ShouldNeulandPointBeFiltered function
    TEST(NeulandPointFilterTest, ShouldNeulandPointBeFilteredTest)
    {
        constexpr auto proton_pid = 2212;
        constexpr auto neutron_pid = 2112;
        constexpr auto eLoss_proton = 0.7;
        constexpr auto eLoss_neutron = 0.3;
        const TVector3 sample_vector;

        NeulandPointFilter filter;

        const R3BNeulandPoint protonPoint(
            0, 0, sample_vector, sample_vector, 0., 0., eLoss_proton, 0, 0., proton_pid, 0);
        const R3BNeulandPoint neutronPoint(
            0, 0, sample_vector, sample_vector, 0., 0., eLoss_neutron, 0, 0., neutron_pid, 0);

        // Test filtering criteria for protons
        filter.SetFilter(BitSetParticle::proton);
        EXPECT_TRUE(filter.CheckFiltered(protonPoint));
        EXPECT_FALSE(filter.CheckFiltered(neutronPoint));

        // Test minimum energy filter
        filter.SetFilter(BitSetParticle::neutron, 0.5);
        EXPECT_TRUE(filter.CheckFiltered(neutronPoint));
        EXPECT_FALSE(filter.CheckFiltered(protonPoint));
    }
} // namespace
