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

#include "R3BMusicCalPar.h"
#include "R3BMusicHitPar.h"
#include "gtest/gtest.h"
#include <map>

namespace
{
    TEST(testMusicCalPar, GetInUse)
    {
        R3BMusicCalPar par;
        par.SetInUse(1, 4);

        EXPECT_EQ(par.GetInUse(4), 1);
    }

    TEST(testMusicCalPar, GetNumAnodes)
    {
        R3BMusicCalPar par;
        par.SetNumAnodes(9);

        EXPECT_EQ(par.GetNumAnodes(), 9);
    }

    TEST(testMusicHitPar, GetInUse)
    {
        R3BMusicHitPar par;
        par.SetInUse(1, 4);

        EXPECT_EQ(par.GetInUse(4), 1);
    }

    TEST(testMusicHitPar, GetNumAnodes)
    {
        R3BMusicHitPar par;
        par.SetNumAnodes(8);

        EXPECT_EQ(par.GetNumAnodes(), 8);
    }
} // namespace
