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

#include "R3BIncomingIDPar.h"
#include "R3BLogger.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayF.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>

// ---- Standard Constructor ---------------------------------------------------
R3BIncomingIDPar::R3BIncomingIDPar(const TString& name, const TString& title, const TString& context)
    : FairParGenericSet(name, title, context)
    , fNumDet(1)
{
    fToFoffset = new TArrayF(fNumDet);
    fPosS2Left = new TArrayF(fNumDet);
    fPosS2Right = new TArrayF(fNumDet);
    fBrho0_S2toCC = new TArrayF(fNumDet);
    fDispersionS2 = new TArrayF(fNumDet);
    fTof2InvV_p0 = new TArrayF(fNumDet);
    fTof2InvV_p1 = new TArrayF(fNumDet);
    // TCut to skip charge states at S2
    fTcutparS2 = new R3BTcutPar("CutS2");
    // TCut to select incoming projectile at the entrance of Cave-C
    fTcutparCave = new R3BTcutPar("CutCave");
}

// ----  Destructor ------------------------------------------------------------
R3BIncomingIDPar::~R3BIncomingIDPar()
{
    R3BLOG(debug, "");
    clear();
    if (fToFoffset)
        delete fToFoffset;
    if (fPosS2Left)
        delete fPosS2Left;
    if (fPosS2Right)
        delete fPosS2Right;
    if (fBrho0_S2toCC)
        delete fBrho0_S2toCC;
    if (fDispersionS2)
        delete fDispersionS2;
    if (fTof2InvV_p0)
        delete fTof2InvV_p1;
    if (fDispersionS2)
        delete fTof2InvV_p0;
}

// ----  Method clear ----------------------------------------------------------
void R3BIncomingIDPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BIncomingIDPar::putParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(fatal, "Could not find FairParamList");
        return;
    }
    fToFoffset->Set(fNumDet);
    list->add("TofOffset", *fToFoffset);
    fPosS2Left->Set(fNumDet);
    list->add("PosS2Left", *fPosS2Left);
    fPosS2Right->Set(fNumDet);
    list->add("PosS2Right", *fPosS2Right);
    fBrho0_S2toCC->Set(fNumDet);
    list->add("Brho0_S2toCC", *fBrho0_S2toCC);
    fDispersionS2->Set(fNumDet);
    list->add("DispersionS2", *fDispersionS2);
    fTof2InvV_p0->Set(fNumDet);
    list->add("Tof2InvV_p0", *fTof2InvV_p0);
    fTof2InvV_p1->Set(fNumDet);
    list->add("DTof2InvV_p1", *fTof2InvV_p1);
    list->add("x0_Aq", fx0_Aq);
    list->add("y0_Aq", fy0_Aq);
    list->add("ang_Aq", fang_Aq);
    list->add("Beta_min", fBeta_min);
    list->add("Beta_max", fBeta_max);

    fTcutparS2->putParams(list);
    fTcutparCave->putParams(list);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BIncomingIDPar::getParams(FairParamList* list)
{
    R3BLOG(info, "called");
    if (!list)
    {
        R3BLOG(fatal, "Could not initialize FairParamList");
        return kFALSE;
    }

    if (!(list->fill("TofOffset", fToFoffset)))
    {
        LOG(error) << "Could not initialize incomingIDTofOffsetPar";
        return kFALSE;
    }

    if (!(list->fill("PosS2Left", fPosS2Left)))
    {
        LOG(error) << "Could not initialize incomingIDPosS2LeftPar";
        return kFALSE;
    }

    if (!(list->fill("PosS2Right", fPosS2Right)))
    {
        LOG(error) << "Could not initialize incomingIDPosS2RightPar";
        return kFALSE;
    }

    if (!(list->fill("Brho0_S2toCC", fBrho0_S2toCC)))
    {
        LOG(error) << "Could not initialize incomingIDBrho0_S2toCCPar";
        return kFALSE;
    }

    if (!(list->fill("DispersionS2", fDispersionS2)))
    {
        LOG(error) << "---Could not initialize incomingIDDispersionS2Par";
        return kFALSE;
    }

    if (!(list->fill("Tof2InvV_p0", fTof2InvV_p0)))
    {
        LOG(error) << "Could not initialize incomingIDTof2InvV_p0Par";
        return kFALSE;
    }

    if (!(list->fill("Tof2InvV_p1", fTof2InvV_p1)))
    {
        LOG(error) << "Could not initialize incomingIDTof2InvV_p1Par";
        return kFALSE;
    }

    if (!(list->fill("x0_Aq", &fx0_Aq)))
    {
        LOG(error) << "Could not initialize incomingIDx0_AqPar";
        return kFALSE;
    }

    if (!(list->fill("y0_Aq", &fy0_Aq)))
    {
        LOG(error) << "Could not initialize incomingIDy0_AqPar";
        return kFALSE;
    }

    if (!(list->fill("ang_Aq", &fang_Aq)))
    {
        LOG(error) << "Could not initialize incomingIDang_AqPar";
        return kFALSE;
    }

    if (!(list->fill("Beta_min", &fBeta_min)))
    {
        LOG(error) << "Could not initialize incomingID Beta_min Par";
        return kFALSE;
    }

    if (!(list->fill("Beta_max", &fBeta_max)))
    {
        LOG(error) << "Could not initialize incomingID Beta_max Par";
        return kFALSE;
    }

    // TCuts
    fTcutparS2->getParams(list);
    fTcutparCave->getParams(list);

    return kTRUE;
}

// ----  Method print ----------------------------------------------------------
void R3BIncomingIDPar::print() { printParams(); }

// ----  Method printParams ----------------------------------------------------
void R3BIncomingIDPar::printParams()
{
    R3BLOG(info, "Incoming ID parameters");

    LOG(info) << "Brho: " << fBrho0_S2toCC->GetAt(0) << ", DispersionS2toCC: " << fDispersionS2->GetAt(0);

    for (Int_t d = 0; d < fNumDet; d++)
        LOG(info) << "Sci: " << d + 1 << ", TofOffset: " << fToFoffset->GetAt(d) << " PosS2Left" << fPosS2Left->GetAt(d)
                  << ", PosS2Right: " << fPosS2Right->GetAt(d);

    fTcutparS2->print();
    fTcutparCave->print();
}

ClassImp(R3BIncomingIDPar);
