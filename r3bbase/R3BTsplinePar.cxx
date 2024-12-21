
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
// -----                   R3BTsplinePar                    -----
// -----   Created Jan 23th 2022 by J.L. Rodriguez-Sanchez  -----
// --------------------------------------------------------------

#include "R3BTsplinePar.h"
#include "R3BLogger.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayF.h"
#include "TF1.h"
#include "TMath.h"

// -----   Default constructor   -----------------------------------------------
R3BTsplinePar::R3BTsplinePar(const TString& name)
    : TNamed(name, name)
    , fSpline(NULL)
    , fMaxPoints(4000) // 100
    , fxmin(-100.50)
    , fxmax(1030)
{
    // fSpline = new TSpline3();
}

// ----  Method putParams ------------------------------------------------------
void R3BTsplinePar::putParams(FairParamList* list)
{
    R3BLOG(info, "called for " << fSpline->GetNp() << " points");

    TArrayF* p = new TArrayF(5);
    // TArrayF* q = new TArrayF(3);
    for (int n = 0; n < fSpline->GetNp(); n++)
    {
        Double_t a[2];
        Double_t b[6];
        fSpline->GetKnot(n, a[0], a[1]);
        fSpline->GetCoeff(n, b[0], b[1], b[2], b[3], b[4]);
        TString name = GetNameObj() + "Point" + n;
        p->AddAt(a[0], 0);
        p->AddAt(a[1], 1);
        p->AddAt(b[2], 2);
        p->AddAt(b[3], 3);
        p->AddAt(b[4], 4);
        list->add(name, *p);
        /*TString namecoeff = GetNameObj() + "PointCoeff" + n;
        q->AddAt(b[0], 0);
        q->AddAt(b[1], 1);
        q->AddAt(b[2], 2);
        list->add(namecoeff, *q);*/
    }
    delete p;
}

// ----  Method getParams ------------------------------------------------------
TSpline3* R3BTsplinePar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");
    // if (fSpline == NULL)
    fSpline = new TSpline3("spline_func", fxmin, fxmax, (TF1*)0, fMaxPoints);
    TArrayF* p = new TArrayF(5);
    // TArrayF* q = new TArrayF(3);
    for (int n = 0; n < fMaxPoints; n++)
    {
        TString name = GetNameObj() + "Point" + n;
        // TString namecoeff = GetNameObj() + "PointCoeff" + n;
        if (list->fill(name, p))
        {
            fSpline->SetPoint(n, p->GetAt(0), p->GetAt(1));
            fSpline->SetPointCoeff(n, p->GetAt(2), p->GetAt(3), p->GetAt(4));
        }
        else
        {
            break;
        }
    }
    delete p;
    // std::cout << "Np spline = " << fSpline->GetNp() << std::endl;
    return fSpline;
}

// ----  Method print ----------------------------------------------------------
void R3BTsplinePar::print()
{
    if (fSpline)
    {
        R3BLOG(info, "for " << fSpline->GetName());

        for (int n = 0; n < fSpline->GetNp(); n++)
        {
            Double_t a[2];
            Double_t b[6];
            fSpline->GetKnot(n, a[0], a[1]);
            fSpline->GetCoeff(n, b[0], b[1], b[2], b[3], b[4]);
            LOG(info) << "CutPoint" << n << ": " << a[0] << " ; " << a[1];
            LOG(info) << "CutPoint" << n << ": " << b[0] << " ; " << b[1] << " ; " << b[2];
        }
    }
    else
    {
        R3BLOG(warn, "Spline parameters not found");
    }
}

ClassImp(R3BTsplinePar);
