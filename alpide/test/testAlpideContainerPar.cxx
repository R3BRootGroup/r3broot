
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

#include "R3BAlpideMappingPar.h"
#include "gtest/gtest.h"
#include <map>

namespace
{
    TEST(testAlpideMappingPar, GetNbSensors)
    {
        auto par = new R3BAlpideMappingPar();
        par->SetNbSensors(2);

        EXPECT_EQ(par->GetNbSensors(), 2);
    }
} // namespace
