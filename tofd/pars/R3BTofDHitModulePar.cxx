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

#include "R3BTofDHitModulePar.h"

#include "FairLogger.h"
#include "FairParamList.h"

// #include "TF1.h"
// #include "TH1F.h"

R3BTofDHitModulePar::R3BTofDHitModulePar(const char* name, const char* title, const char* context, Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fPlane(0)
    , fPaddle(0)
    , fOffset1(0)
    , fOffset2(0)
    , fToTOffset1(0)
    , fToTOffset2(0)
    , fVeff(1.)
    , fLambda(1.)
    , fSync(0)
    , fTofSyncOffset(0)
    , fPar1a(0)
    , fPar1b(0)
    , fPar1c(0)
    , fPar1d(0)
    , fPar2a(0)
    , fPar2b(0)
    , fPar2c(0)
    , fPar2d(0)
    , fPola(1.)
    , fPolb(0)
    , fPolc(0)
    , fPold(0)
    , fPar1za(0)
    , fPar1zb(0)
    , fPar1zc(0)
    , fPar1walk(0)
    , fPar2walk(0)
    , fPar3walk(0)
    , fPar4walk(0)
    , fPar5walk(0)
{
}

R3BTofDHitModulePar::~R3BTofDHitModulePar() { clear(); }

void R3BTofDHitModulePar::putParams(FairParamList* list)
{
    LOG(info) << "R3BTofDHitModulePar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("plane", fPlane);
    list->add("paddle", fPaddle);
}

Bool_t R3BTofDHitModulePar::getParams(FairParamList* list)
{
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("plane", &fPlane))
    {
        return kFALSE;
    }
    if (!list->fill("paddle", &fPaddle))
    {
        return kFALSE;
    }
    return kTRUE;
}

void R3BTofDHitModulePar::clear() { fPlane = fPaddle = 0; }

void R3BTofDHitModulePar::printParams()
{
    LOG(info) << "   TOFd HIT Calibration Parameters: ";
    LOG(info) << "   fPlane: " << fPlane;
    LOG(info) << "   fPaddle: " << fPaddle;
    LOG(info) << "   fOffset1: " << fOffset1;
    LOG(info) << "   fOffset2: " << fOffset2;
    LOG(info) << "   fToTOffset1: " << fToTOffset1;
    LOG(info) << "   fToTOffset2: " << fToTOffset2;
    LOG(info) << "   fSync: " << fSync;
    LOG(info) << "   fTofSyncOffset: " << fTofSyncOffset;
    LOG(info) << "   fVeff: " << fVeff;
    LOG(info) << "   fLambda: " << fLambda;
    LOG(info) << "   fDoubleExp1: " << fPar1a << "  " << fPar1b << "  " << fPar1c << "  " << fPar1d;
    LOG(info) << "   fDoubleExp2: " << fPar2a << "  " << fPar2b << "  " << fPar2c << "  " << fPar2d;
    LOG(info) << "   fPol3: " << fPola << "  " << fPolb << "  " << fPolc << "  " << fPold;
    LOG(info) << "   fParaZ: " << fPar1za << "  " << fPar1zb << "  " << fPar1zc;
    LOG(info) << "   fParaWalk: " << fPar1walk << "  " << fPar2walk << "  " << fPar3walk << "  " << fPar4walk << "  "
              << fPar5walk;
}

void R3BTofDHitModulePar::DrawParams()
{
    /*    Int_t type = 2; // VFTX
        if (fNofChannels > 0)
        {
            if (fSlope[0] > 0)
            {
                type = 1; // Tacquila
            }
        }
        else
        {
            return;
        }
        TH1F* h1 = new TH1F(Form("h1_%d_%d_%d", fPlane, fPaddle, fSide), "", 4096, -0.5, 4096.5);
        h1->Draw();
        for (Int_t i = 0; i < fNofChannels; i++)
        {
            if (1 == type)
            {
                TF1* f1 = new TF1(Form("f1_%d", i), "[0] + [1]*(x - [2])", fBinLow[i], fBinUp[i]);
                f1->SetParameter(0, fOffset[i]);
                f1->SetParameter(1, fSlope[i]);
                f1->SetParameter(2, fBinLow[i]);
                f1->Draw("same");
            }
            else if (2 == type)
            {
                h1->SetBinContent(fBinLow[i], fOffset[i]);
            }
        }

        if (1 == type)
        {
            h1->GetYaxis()->SetRangeUser(0., 1.2 * fOffset[fNofChannels - 1]);
        }
        else if (2 == type)
        {
            h1->Draw();
        }
        gPad->Update();*/
}

ClassImp(R3BTofDHitModulePar);
