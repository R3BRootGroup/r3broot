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

// --------------------------------------------------------------
// -----                   R3BTcutPar                       -----
// -----   Created Jan 23th 2022 by J.L. Rodriguez-Sanchez  -----
// --------------------------------------------------------------

#include "R3BTcutPar.h"
#include "R3BLogger.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayF.h"
#include "TMath.h"

// -----   Default constructor   -----------------------------------------------
R3BTcutPar::R3BTcutPar(const TString& cutname)
    : TNamed(cutname, cutname)
    , fMaxPoints(40)
    , fCut(NULL)
{
}

// ----  Method putParams ------------------------------------------------------
void R3BTcutPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called for " << fCut->GetN() << " points");

    TArrayF* p = new TArrayF(2);
    for (int n = 0; n < fCut->GetN(); n++)
    {
        Double_t a[2];
        fCut->GetPoint(n, a[0], a[1]);

        TString name = GetNameObj() + "Point" + n;
        p->AddAt(a[0], 0);
        p->AddAt(a[1], 1);
        list->add(name, *p);
    }
    delete p;
}

// ----  Method getParams ------------------------------------------------------
TCutG* R3BTcutPar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");

    TArrayF* p = new TArrayF(2);
    for (int n = 0; n < fMaxPoints; n++)
    {
        TString name = GetNameObj() + "Point" + n;
        if (list->fill(name, p))
        {
            if (fCut == NULL)
                fCut = new TCutG(GetNameObj(), 1);
            fCut->SetPoint(n, p->GetAt(0), p->GetAt(1));
        }
        else
        {
            break;
        }
    }
    delete p;
    return fCut;
}

// ----  Method print ----------------------------------------------------------
void R3BTcutPar::print()
{
    if (fCut)
    {
        R3BLOG(info, "for " << fCut->GetName());
        for (int n = 0; n < fCut->GetN(); n++)
        {
            Double_t a[2];
            fCut->GetPoint(n, a[0], a[1]);
            LOG(info) << "CutPoint" << n << ": " << a[0] << " ; " << a[1];
        }
    }
    else
    {
        R3BLOG(warn, "Cut parameters not found");
    }
}

ClassImp(R3BTcutPar);
