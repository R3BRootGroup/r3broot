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

// ------------------------------------------------------------
// -----                  R3BPtofCal2Hit                  -----
// -----          Created Mar 1st 2017 by V.Wagner        -----
// ------------------------------------------------------------

#include "R3BPtofCal2Hit.h"
#include "R3BPaddleCalData.h"
#include "R3BPtofHit.h"
#include "R3BPtofHitPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TArrayD.h"
#include "TClonesArray.h"
#include "TMath.h"

R3BPtofCal2Hit::R3BPtofCal2Hit()
    : FairTask("TofdTcal", 1)
    , fCalItems(NULL)
    , fHitItems(new TClonesArray("R3BPtofHit"))
    , fHitPar(NULL)
    , fNEvents(0)
{
    fTOffset1.Set(2 * PtofPaddlesPerPlane);
    fTOffset2.Set(2 * PtofPaddlesPerPlane);
    fZScale.Set(2 * PtofPaddlesPerPlane);
    fVEff.Set(2 * PtofPaddlesPerPlane);
}

R3BPtofCal2Hit::R3BPtofCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(NULL)
    , fHitItems(new TClonesArray("R3BPtofHit"))
    , fHitPar(NULL)
    , fNEvents(0)
{
    fTOffset1.Set(2 * PtofPaddlesPerPlane);
    fTOffset2.Set(2 * PtofPaddlesPerPlane);
    fZScale.Set(2 * PtofPaddlesPerPlane);
    fVEff.Set(2 * PtofPaddlesPerPlane);
}

R3BPtofCal2Hit::~R3BPtofCal2Hit()
{
    delete fHitItems;
    fHitItems = NULL;
}

InitStatus R3BPtofCal2Hit::Init()
{
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(fatal) << "FairRootManager not found";
    }

    fCalItems = dynamic_cast<TClonesArray*>(fMan->GetObject("PtofCal"));
    if (!fCalItems)
    {
        LOG(fatal) << "PtofCalData not found";
    }

    fMan->Register("PtofHit", "Ptof", fHitItems, kTRUE);

    SetParameter();
    return kSUCCESS;
}

void R3BPtofCal2Hit::SetParContainers()
{
    fHitPar = dynamic_cast<R3BPtofHitPar*>(FairRuntimeDb::instance()->getContainer("PtofHitPar"));
    ;
}

void R3BPtofCal2Hit::SetParameter()
{

    for (Int_t paddle = 0; paddle < 2 * PtofPaddlesPerPlane; paddle++)
    {
        fTOffset1[paddle] = fHitPar->GetTOffsetAt(paddle + 1, 1);
        fTOffset2[paddle] = fHitPar->GetTOffsetAt(paddle + 1, 2);
        fZScale[paddle] = fHitPar->GetZScaletAt(paddle + 1);
        fVEff[paddle] = fHitPar->GetVEfftAt(paddle + 1);
        LOG(error) << "Offset1: " << fTOffset1[paddle] << "  offset2: " << fTOffset2[paddle]
                   << "  ZScale: " << fZScale[paddle] << "  Veff: " << fVEff[paddle];
    }
}

InitStatus R3BPtofCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BPtofCal2Hit::Exec(Option_t* option)
{
    const Double_t paddle_width = 5.0;    // cm
    const Double_t paddle_height = 105.0; // cm
    const Double_t paddle_depth = 0.5;    // cm
    const Double_t plane_offset = 0.0;    // cm

    Double_t charge[2][PtofPaddlesPerPlane];
    Double_t time[2][PtofPaddlesPerPlane];
    Double_t xx[2][PtofPaddlesPerPlane];
    Double_t yy[2][PtofPaddlesPerPlane];
    Double_t zz[2][PtofPaddlesPerPlane];

    for (Int_t i = 0; i < PtofPaddlesPerPlane; i++)
    {
        charge[0][i] = 0.;
        charge[1][i] = 0.;
        time[0][i] = 0.;
        time[1][i] = 0.;
        xx[0][i] = 0.;
        xx[1][i] = 0.;
        yy[0][i] = 0.;
        yy[1][i] = 0.;
        zz[0][i] = 0.;
        zz[1][i] = 0.;
    }

    Int_t nItems = fCalItems->GetEntriesFast();
    Int_t nHitData = 0;

    for (Int_t i = 0; i < nItems; i++)
    {
        R3BPaddleCalData* caldata = dynamic_cast<R3BPaddleCalData*>(fCalItems->At(i));
        Int_t iPlane = caldata->GetPlane() - 1;
        Int_t iBar = caldata->GetBar() - 1;
        assert(iBar >= 0 && iBar < PtofPaddlesPerPlane);
        Int_t id = iPlane * PtofPaddlesPerPlane + iBar;

        Double_t t1 = caldata->fTime1L_ns + fTOffset1[id];
        Double_t t2 = caldata->fTime2L_ns + fTOffset2[id];
        time[iPlane][iBar] = (t1 + t2) * 0.5;

        charge[iPlane][iBar] =
            sqrt((caldata->fTime1T_ns - caldata->fTime1L_ns) * (caldata->fTime2T_ns - caldata->fTime2L_ns)) *
            fZScale[id];

        xx[iPlane][iBar] = (id % PtofPaddlesPerPlane + (id < PtofPaddlesPerPlane ? 0 : -0.5)) * paddle_width;
        yy[iPlane][iBar] = (t2 - t1) * fVEff[id];
        zz[iPlane][iBar] = paddle_depth * 0.5 + (id < PtofPaddlesPerPlane ? 0 : paddle_depth + plane_offset);

        //		For KVI experiment, bar 8 was a horizontal bar
        if (iBar + 1 == 8)
        {
            yy[iPlane][iBar] = 0;
            xx[iPlane][iBar] = (t2 - t1) * 6.5;
            zz[iPlane][iBar] = 10.;
        }

        if (PtofPaddlesPerPlane >= 40 && iBar + 1 == 40)
        {
            LOG(error) << "cal2hit: Plane: " << iPlane + 1 << "  Bar: " << iBar + 1 << "  t1= " << t1 << "  t2= " << t2;
            LOG(error) << "cal2hit: xx: " << xx[iPlane][iBar] << "  charge: " << charge[iPlane][iBar];
        }

        //		new ((*fHitItems)[nHitData++]) R3BPtofHit(id+1, charge, time, xx, yy, zz);
    }
    for (Int_t iBar = 0; iBar < PtofPaddlesPerPlane; iBar++)
    {
        if (PtofPaddlesPerPlane >= 40 && iBar + 1 == 40)
        {
            LOG(error) << "Planne 1  Bar: " << iBar + 1 << "  charge: " << charge[0][iBar] << "  x: " << xx[0][iBar];
            LOG(error) << "Planne 2  Bar: " << iBar + 1 << "  charge: " << charge[1][iBar] << "  x: " << xx[1][iBar];
        }
        /*
                if(time[0][iBar]>0. && time[1][iBar]>0.){
                    Double_t charge_avg = (charge[0][iBar] + charge[1][iBar])/2.;
                    Double_t time_avg = (time[0][iBar] + time[1][iBar])/2.;
                    Double_t xx_avg = (xx[0][iBar] + xx[1][iBar])/2.;
                    Double_t yy_avg = (yy[0][iBar] + yy[1][iBar])/2.;
                    Double_t zz_avg = (zz[0][iBar] + zz[1][iBar])/2.;

        //			LOG(error)<<"Bar: "<<iBar+1<<"  charge: "<<charge_avg<<"  time: "<<time_avg<<"  x: "<<xx_avg;

                    new ((*fHitItems)[nHitData++]) R3BPtofHit(iBar+1, charge_avg, time_avg, xx_avg, yy_avg, zz_avg);
                }
        */
        if (charge[0][iBar] > 0. || charge[1][iBar] > 1.)
        {
            if (abs(charge[0][iBar] - charge[1][iBar]) / (charge[0][iBar] + charge[1][iBar]) < 0.1)
            {
                Double_t charge_avg = (charge[0][iBar] + charge[1][iBar]) / 2.;
                Double_t time_avg = (time[0][iBar] + time[1][iBar]) / 2.;
                Double_t xx_avg = (xx[0][iBar] + xx[1][iBar]) / 2.;
                Double_t yy_avg = (yy[0][iBar] + yy[1][iBar]) / 2.;
                Double_t zz_avg = (zz[0][iBar] + zz[1][iBar]) / 2.;

                new ((*fHitItems)[nHitData++]) R3BPtofHit(iBar + 1, charge_avg, time_avg, xx_avg, yy_avg, zz_avg);
            }
            else if (charge[0][iBar] > charge[1][iBar])
            {
                Double_t charge_avg = charge[0][iBar];
                Double_t time_avg = time[0][iBar];
                Double_t xx_avg = xx[0][iBar];
                Double_t yy_avg = yy[0][iBar];
                Double_t zz_avg = zz[0][iBar];
            }
            else if (charge[1][iBar] > charge[0][iBar])
            {
                Double_t charge_avg = charge[1][iBar];
                Double_t time_avg = time[1][iBar];
                Double_t xx_avg = xx[1][iBar];
                Double_t yy_avg = yy[1][iBar];
                Double_t zz_avg = zz[1][iBar];

                //		    	if (iBar==4 && fNEvents>500000)
                //				LOG(error)<<"Event: "<<fNEvents<<"  Bar: "<<iBar+1<<"  time1: "<<time[0][iBar]<<"  time:
                //"<<time[1][iBar];

                new ((*fHitItems)[nHitData++]) R3BPtofHit(iBar + 1, charge_avg, time_avg, xx_avg, yy_avg, zz_avg);
            }
        }
    }
    fNEvents += 1;
}

void R3BPtofCal2Hit::FinishEvent() { fHitItems->Clear(); }

void R3BPtofCal2Hit::FinishTask() {}

ClassImp(R3BPtofCal2Hit)
