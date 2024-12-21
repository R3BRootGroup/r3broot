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

#include "R3BFiberMAPMTHitModulePar.h"

#include "FairLogger.h"
#include "FairParamList.h"

// #include "TF1.h"
// #include "TH1F.h"
// #include "TPad.h"

R3BFiberMAPMTHitModulePar::R3BFiberMAPMTHitModulePar(const char* name,
                                                     const char* title,
                                                     const char* context,
                                                     Bool_t own)
    : FairParGenericSet(name, title, context, own)
    , fFiber(0)
    , fOffsetUp(0)
    , fOffsetDown(0)
    , fVeff(0)
    , fSync(0)
    , fGainUp(0)
    , fGainDown(0)
{
}

R3BFiberMAPMTHitModulePar::~R3BFiberMAPMTHitModulePar()
{
    // Reset all parameters
    clear();
}

void R3BFiberMAPMTHitModulePar::putParams(FairParamList* list)
{
    LOG(info) << "R3BFiberMAPMTHitModulePar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("fiber", fFiber);
    /*
        list->add("OffsetUp", fOffsetUp);
        list->add("OffsetDown", fOffsetDown);
        list->add("Veff", fOffset2);
    */
}

Bool_t R3BFiberMAPMTHitModulePar::getParams(FairParamList* list)
{
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("fiber", &fFiber))
    {
        return kFALSE;
    }
    return kTRUE;
}

void R3BFiberMAPMTHitModulePar::clear() { fFiber = 0; }

void R3BFiberMAPMTHitModulePar::printParams()
{
    LOG(info) << "   R3BFiberMAPMTHitModulePar: FIBER HIT Calibration Parameters: ";
    LOG(info) << "   fFiber: " << fFiber;
    LOG(info) << "   fOffsetUp: " << fOffsetUp;
    LOG(info) << "   fOffsetDown: " << fOffsetDown;
    LOG(info) << "   fVeff: " << fVeff;
    LOG(info) << "   fGainUp: " << fGainUp;
    LOG(info) << "   fGainDown: " << fGainDown;
    LOG(info) << "   fSync: " << fSync;
}

void R3BFiberMAPMTHitModulePar::DrawParams()
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

ClassImp(R3BFiberMAPMTHitModulePar);
