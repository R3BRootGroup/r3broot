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

// -------------------------------------------------------------
// -----       R3BTwimVertexReconstruction source file     -----
// -----      Created on 18/07/22 by A. Grana Gonzalez     -----
// -------------------------------------------------------------

// ROOT headers
#include "TClonesArray.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVectorD.h"

// FAIR headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

// Twim headers
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTwimHitData.h"
#include "R3BTwimVertexReconstruction.h"

// R3BTwimVertexReconstruction: Default Constructor --------------------------
R3BTwimVertexReconstruction::R3BTwimVertexReconstruction()
    : R3BTwimVertexReconstruction("R3BTwimVertexReconstruction", 1)
{
}

// R3BTwimVertexReconstructionPar: Standard Constructor --------------------------
R3BTwimVertexReconstruction::R3BTwimVertexReconstruction(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTwimHitDataCA(NULL)
    , f1(NULL)
    , f2(NULL)
    , fTpat(-1)
{
}

// Virtual R3BTwimVertexReconstruction: Destructor
R3BTwimVertexReconstruction::~R3BTwimVertexReconstruction() { R3BLOG(debug1, ""); }

// -----   Public method Init   --------------------------------------------
InitStatus R3BTwimVertexReconstruction::Init()
{
    R3BLOG(info, "");
    FairRootManager* rootManager = FairRootManager::Instance();
    R3BLOG_IF(fatal, !rootManager, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
    R3BLOG_IF(warn, !header, "EventHeadder. not found");

    // INPUT DATA
    // get access to twim hit data
    fTwimHitDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("TwimHitData"));
    R3BLOG_IF(fatal, !fTwimHitDataCA, "TwimHitData not found");

    return kSUCCESS;
}

// -----   Public method Execution   --------------------------------------------
void R3BTwimVertexReconstruction::Exec(Option_t* option)
{
    if ((fTpat >= 0) && (header) && ((header->GetTpat() & fTpat) != fTpat))
    {
        return;
    }

    Int_t nHitTwim = fTwimHitDataCA->GetEntriesFast();
    // The two fission fragments are needed to reconstruct the vertex
    if (nHitTwim != 2)
    {
        return;
    }

    TVector3 trajparams[2];
    trajparams[0].SetXYZ(0., 0., 0.);
    trajparams[1].SetXYZ(0., 0., 0.);

    auto HitDat = new R3BTwimHitData*[nHitTwim];

    for (Int_t i = 0; i < nHitTwim; i++)
    {
        HitDat[i] = dynamic_cast<R3BTwimHitData*>(fTwimHitDataCA->At(i));
        trajparams[i].SetXYZ(HitDat[i]->GetSecID(), HitDat[i]->GetOffset(), HitDat[i]->GetTheta());
    }

    f1 = new TF1("f1", "pol1", -2800, 1000);
    f2 = new TF1("f2", "pol1", -2800, 1000);
    f1->SetParameter(0, trajparams[0].Y());
    f1->SetParameter(1, trajparams[0].Z());
    f2->SetParameter(0, trajparams[1].Y());
    f2->SetParameter(1, trajparams[1].Z());
    double min = 2000.; // in mm
    double xmin = NAN, zmin = NAN;
    for (int i = -2800; i < 1000; i++) // i is in mm
    {
        auto mintemp = TMath::Abs(f1->Eval(i) - f2->Eval(i));
        if (min > mintemp)
        {
            min = mintemp;
            zmin = i;
            xmin = (f1->Eval(i) + f2->Eval(i)) / 2.0; // Mean value for X
        }
    }

    if (!TMath::IsNaN(xmin) && !TMath::IsNaN(zmin))
    {
        HitDat[0]->SetVertexZ(zmin);
        HitDat[1]->SetVertexZ(zmin);
        HitDat[0]->SetVertexX(xmin);
        HitDat[1]->SetVertexX(xmin);
    }

    if (HitDat)
    {
        delete[] HitDat;
    }
    return;
}

ClassImp(R3BTwimVertexReconstruction);
