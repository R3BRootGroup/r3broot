
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

#include "R3BCalifaCrystalCalPar.h"
#include "R3BCalifaMappingPar.h"
#include "gtest/gtest.h"
#include <map>

namespace
{
    TEST(testCalifaMappingPar, GetNumCrystals)
    {
        R3BCalifaMappingPar par;
        par.SetNumCrystals(2000);

        EXPECT_EQ(par.GetNumCrystals(), 2000);
    }

    TEST(testCalifaMappingPar, GetInUse)
    {
        R3BCalifaMappingPar par;
        par.SetInUse(1, 1000);

        EXPECT_EQ(par.GetInUse(1000), 1);
    }

    TEST(testCalifaCrystalCalPar, GetNumCrystals)
    {
        R3BCalifaCrystalCalPar par;
        par.SetNumCrystals(100);

        EXPECT_EQ(par.GetNumCrystals(), 100);
    }

    TEST(testCalifaCrystalCalPar, GetNumParametersFit)
    {
        R3BCalifaCrystalCalPar par;
        par.SetNumParametersFit(10);

        EXPECT_EQ(par.GetNumParametersFit(), 10);
    }

} // namespace
