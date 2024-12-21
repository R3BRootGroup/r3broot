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

// ----------------------------------------------------------------------
// -----                          R3BLogger                         -----
// -----             Created 15/12/21 by J.L. Rodriguez-Sanchez     -----
// ----------------------------------------------------------------------

#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <FairLogger.h>

// NOLINTBEGIN
class R3BLogger;

class R3BLogger : public FairLogger
{
  public:
#define R3BLOG(severity, x)                                                                            \
    if (true)                                                                                          \
    {                                                                                                  \
        std::string fN(__FILE__);                                                                      \
        std::stringstream ss;                                                                          \
        ss << fN.substr(fN.find_last_of("/") + 1) << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
        LOG(severity) << ss.str() << x;                                                                \
    }                                                                                                  \
    else                                                                                               \
        (void)0

#define R3BLOG_IF(severity, condition, x)                                                                    \
    if (true)                                                                                                \
    {                                                                                                        \
        std::string fNif(__FILE__);                                                                          \
        std::stringstream ssif;                                                                              \
        ssif << fNif.substr(fNif.find_last_of("/") + 1) << ":" << __LINE__ << ":" << __FUNCTION__ << "(): "; \
        LOG_IF(severity, condition) << ssif.str() << x;                                                      \
    }                                                                                                        \
    else                                                                                                     \
        (void)0

  private:
    R3BLogger();
    ~R3BLogger();

  public:
    ClassDefOverride(R3BLogger, 0)
};
// NOLINTEND
