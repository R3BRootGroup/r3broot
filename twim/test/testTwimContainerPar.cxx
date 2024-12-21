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

#include "R3BTwimCalPar.h"
#include "R3BTwimHitPar.h"
#include "gtest/gtest.h"
#include <map>

namespace
{
    TEST(testTwimCalPar, GetNumSec)
    {
        R3BTwimCalPar par;
        par.SetNumSec(4);

        EXPECT_EQ(par.GetNumSec(), 4);
    }

    TEST(testTwimCalPar, GetNumAnodes)
    {
        R3BTwimCalPar par;
        par.SetNumAnodes(16);

        EXPECT_EQ(par.GetNumAnodes(), 16);
    }

    TEST(testTwimHitPar, GetNumSec)
    {
        R3BTwimHitPar par;
        par.SetNumSec(4);

        EXPECT_EQ(par.GetNumSec(), 4);
    }
} // namespace
