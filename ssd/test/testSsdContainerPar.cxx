
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

#include "R3BAmsMappingPar.h"
#include "R3BAmsStripCalPar.h"
#include "R3BFootCalPar.h"
#include "R3BFootMappingPar.h"
#include "gtest/gtest.h"
#include <map>

namespace
{
    TEST(testAmsMappingPar, GetNumDets)
    {
        R3BAmsMappingPar par;
        par.SetNumDets(6);

        EXPECT_EQ(par.GetNumDets(), 6);
    }

    TEST(testAmsStripCalPar, GetNumStrips)
    {
        R3BAmsStripCalPar par;
        par.SetNumStrips(300);

        EXPECT_EQ(par.GetNumStrips(), 300);
    }

    TEST(testFootMappingPar, GetNumDets)
    {
        R3BFootMappingPar par;
        par.SetNumDets(8);

        EXPECT_EQ(par.GetNumDets(), 8);
    }

    TEST(testFootStripCalPar, GetNumStrips)
    {
        R3BFootCalPar par;
        par.SetNumStrips(500);

        EXPECT_EQ(par.GetNumStrips(), 500);
    }
} // namespace
