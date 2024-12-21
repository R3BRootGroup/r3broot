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

#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>

#include "R3BMacroCompiler.h"

using namespace std;

#define XSTR(_d) #_d
#define STR(_d) XSTR(_d)

void R3BMacroCompiler::SetIncludeDirectories()
{
    TString dirstr(STR(INCLUDE_DIRECTORIES));
    TObjArray* dirs = dirstr.Tokenize(";");

    for (int i = 0; i < dirs->GetEntries(); i++)
    {
        TObjString* dir = dynamic_cast<TObjString*>(dirs->At(i));
        if (dir)
        {
            gSystem->AddIncludePath("-I" + dir->GetString());
            cout << "  " << dir->GetString() << endl;
        }
    }

    delete dirs;
}
