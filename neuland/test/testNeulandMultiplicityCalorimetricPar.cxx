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

#include "R3BNeulandMultiplicityCalorimetricPar.h"
#include "TCutG.h"
#include "gtest/gtest.h"
#include <map>

namespace
{
    TEST(testNMultiplicityCalorimetricPar, HandlesConversionBetweenMaps)
    {
        std::map<UInt_t, TCutG*> m;
        m[1] = new TCutG("cut1", 4);
        m[2] = new TCutG("cut2", 4);

        R3BNeulandMultiplicityCalorimetricPar par;
        par.SetNeutronCuts(m);

        EXPECT_STREQ(par.GetNeutronCuts().at(2)->GetName(), "cut2");
    }

    TEST(testMultiplicityCalorimetricPar, GetNeutronMultiplicity)
    {
        std::map<UInt_t, TCutG*> m;
        m[0] = new TCutG("cut0", 4);
        m[0]->SetPoint(0, 0, 0);
        m[0]->SetPoint(1, 0, 10);
        m[0]->SetPoint(2, 10, 0);
        m[0]->SetPoint(3, 0, 0);
        m[1] = new TCutG("cut1", 4);
        m[1]->SetPoint(0, 0, 10);
        m[1]->SetPoint(1, 0, 20);
        m[1]->SetPoint(2, 20, 0);
        m[1]->SetPoint(3, 10, 0);
        m[2] = new TCutG("cut2", 4);
        m[2]->SetPoint(0, 0, 20);
        m[2]->SetPoint(1, 0, 30);
        m[2]->SetPoint(2, 30, 0);
        m[2]->SetPoint(3, 20, 0);

        R3BNeulandMultiplicityCalorimetricPar par;
        par.SetNeutronCuts(m);

        EXPECT_EQ(par.GetNeutronMultiplicity(4, 4), 0u);
        EXPECT_EQ(par.GetNeutronMultiplicity(9, 9), 1u);
        EXPECT_EQ(par.GetNeutronMultiplicity(14, 14), 2u);
        EXPECT_EQ(par.GetNeutronMultiplicity(19, 19), 3u);
    }
} // namespace
