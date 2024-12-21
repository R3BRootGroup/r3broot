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

#include "R3BRpcHitPar.h"
#include "R3BRpcTotCalPar.h"
#include "gtest/gtest.h"
#include <map>

namespace
{
    TEST(testRpcTotCalPar, GetNumChannels)
    {
        R3BRpcTotCalPar par;
        par.SetNumChannels(44);

        EXPECT_EQ(par.GetNumChannels(), 44);
    }

    TEST(testRpcHitPar, GetNumChannels)
    {
        R3BRpcHitPar par;
        par.SetNumChannels(88);

        EXPECT_EQ(par.GetNumChannels(), 88);
    }
} // namespace
