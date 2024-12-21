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
// -----                  R3BLosTCal2Hit                   -----
// -----      Created January 17th 2023 by N.Mozumdar      -----
// -------------------------------------------------------------

#include "R3BLosTCal2Hit.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BLosHitData.h"
#include "R3BLosHitPar.h"
#include "R3BLosMapped2TCal.h"
#include "R3BLosMappedData.h"
#include "R3BLosTCalData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <cstdlib>
#include <iterator>

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

R3BLosTCal2Hit::R3BLosTCal2Hit()
    : R3BLosTCal2Hit("LosTCal2Hit", 1)
{
}

R3BLosTCal2Hit::R3BLosTCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTCalItems(NULL)
    , fTCalTriggerItems(NULL)
    , fHitItems(new TClonesArray("R3BLosHitData"))
    , fTimeStitch(NULL)
    , fLEMatchParams(NULL)
    , fTEMatchParams(NULL)
    , fNofHitItems(0)
    , fNofDetectors(1)
    , fTrigger(-1)
    , fTpat(-1)
    , fNumParamsTamexLE(2)
    , fNumParamsTamexTE(2)
    , fWindowV(10.)
    , flosVeffX(1.)
    , flosVeffY(1.)
    , flosOffsetX(0.)
    , flosOffsetY(0.)
    , flosVeffXT(1.)
    , flosVeffYT(1.)
    , flosOffsetXT(0.)
    , flosOffsetYT(0.)
    , flosVeffXQ(1.)
    , flosVeffYQ(1.)
    , flosOffsetXQ(0.)
    , flosOffsetYQ(0.)
    , fOnline(kFALSE)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{
}

R3BLosTCal2Hit::~R3BLosTCal2Hit()
{
    LOG(debug) << "R3BLosTCal2Hit::Destructor";
    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

void R3BLosTCal2Hit::SetParContainers()
{
    LOG(info) << "R3BLosTcal2Hit::SetParContainers()";
    // Parameter Container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }
    fLosHit_Par = dynamic_cast<R3BLosHitPar*>(rtdb->getContainer("LosHitPar"));
    if (!fLosHit_Par)
    {
        LOG(error) << "R3BLosTcal2Hit:: Couldn't get handle on R3BLosHitPar container";
    }
    else
    {
        LOG(info) << "R3BLosTcal2Hit:: R3BLosHitPar container open";
    }
}

void R3BLosTCal2Hit::SetParameter()
{
    //--- Parameter Container ---
    fp0 = fLosHit_Par->Getp0();
    fp1 = fLosHit_Par->Getp1();
    fNumParamsTamexLE = fLosHit_Par->GetNumParamsTamexLE();
    fNumParamsTamexTE = fLosHit_Par->GetNumParamsTamexTE();

    flosOffsetX = fLosHit_Par->Getxoffset_MCFD();
    flosOffsetY = fLosHit_Par->Getyoffset_MCFD();
    flosVeffX = fLosHit_Par->Getxveff_MCFD();
    flosVeffY = fLosHit_Par->Getyveff_MCFD();

    Int_t array_LE = 8 * fNumParamsTamexLE;
    fLEMatchParams = new TArrayF();
    fLEMatchParams->Set(array_LE);
    fLEMatchParams = fLosHit_Par->GetLEMatchParams();

    Int_t array_TE = 8 * fNumParamsTamexTE;
    fTEMatchParams = new TArrayF();
    fTEMatchParams->Set(array_TE);
    fTEMatchParams = fLosHit_Par->GetTEMatchParams();
}

InitStatus R3BLosTCal2Hit::Init()
{
    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(error) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    fTCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("LosTCal"));

    fTCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("LosTriggerTCal"));
    if (NULL == fTCalItems)
    {
        LOG(fatal) << "Branch LosTCal not found";
        return kFATAL;
    }

    // request storage of Hit data in output tree
    mgr->Register("LosHit", "LosHitData", fHitItems, !fOnline);

    fHitItems->Clear();

    fTimeStitch = new R3BCoarseTimeStitch();
    // file with walk-correction parameters
    ifstream infile(fwalk_param_file.c_str());
    if (infile.is_open())
    {
        for (Int_t ivec = 0; ivec < 16; ivec++)
        {
            /*
             ************* Parameters 0-7 MCFD  **********************
             ************* Parameters 8-15 TAMEX *********************
             */

            infile >> walk_par[ivec][0] >> walk_par[ivec][1] >> walk_par[ivec][2] >> walk_par[ivec][3] >>
                walk_par[ivec][4] >> walk_par[ivec][5] >> walk_par[ivec][6] >> walk_par[ivec][7] >> walk_par[ivec][8] >>
                walk_par[ivec][9] >> walk_par[ivec][10];
        }
    }
    else
    {
        LOG(warn) << "UNABLE TO OPEN FILE WITH WALK PARAMETERS! Parameters set to zero!";
        for (Int_t ivec = 0; ivec < 16; ivec++)
        {
            walk_par[ivec][0] = 10.;
            walk_par[ivec][1] = 1000.;
            for (Int_t iv = 2; iv < 11; iv++)
            {
                walk_par[ivec][iv] = 0.;
                ;
            }
        }
    }
    // file with tot-correction parameters
    ifstream infile1(ftot_param_file.c_str());

    if (infile1.is_open())
    {
        for (Int_t ivec = 0; ivec < 8; ivec++)
        {
            infile1 >> tot_par[ivec][0] >> tot_par[ivec][1] >> tot_par[ivec][2] >> tot_par[ivec][3];
        }
    }
    else
    {
        LOG(warn) << "UNABLE TO OPEN FILE WITH ToT PARAMETERS! Parameters set to zero!";
        for (Int_t ivec = 0; ivec < 8; ivec++)
        {
            tot_par[ivec][0] = 0.;
            tot_par[ivec][1] = 0.;
            tot_par[ivec][2] = 0.;
            tot_par[ivec][3] = 1.; // Normalization factor
        }
    }

    SetParameter();
    return kSUCCESS;
}

InitStatus R3BLosTCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BLosTCal2Hit::Exec(Option_t* option)
{
    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = header->GetTpat();
        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if ((header) && (tpatvalue == 0))
            return;
    }

    Int_t trigHits = fTCalTriggerItems->GetEntries();
    Double_t trigTime[trigHits];
    if (trigHits == 0)
        return;
    for (Int_t ihit = 0; ihit < trigHits; ihit++)
    {
        R3BLosTCalData* hit = dynamic_cast<R3BLosTCalData*>(fTCalTriggerItems->At(ihit));
        Int_t typ = hit->GetType();
        trigTime[typ] = hit->GetRawTimeNs();
    }

    Int_t nHits = fTCalItems->GetEntries();
    Int_t hits[fNofDetectors], v_npmt[fNofDetectors][nHits], tle_npmt[fNofDetectors][nHits],
        tte_npmt[fNofDetectors][nHits];
    Double_t vTime[fNofDetectors * 8][nHits], TLeTime[fNofDetectors * 8][nHits], TTeTime[fNofDetectors * 8][nHits],
        avg_v[fNofDetectors][nHits], avg_t[fNofDetectors][nHits];
    Bool_t tHit[fNofDetectors * 8 * 3][nHits];
    for (Int_t d = 0; d < fNofDetectors; d++)
    {
        hits[d] = 0;
        for (Int_t i = 0; i < nHits; i++)
        {
            v_npmt[d][i] = 0;
            tle_npmt[d][i] = 0;
            tte_npmt[d][i] = 0;
            avg_v[d][i] = 0.;
            avg_t[d][i] = 0.;
            for (Int_t ch = 0; ch < 8; ch++)
            {
                vTime[d * 8 + ch][i] = 0.;
                TLeTime[d * 8 + ch][i] = 0.;
                TTeTime[d * 8 + ch][i] = 0.;
                for (Int_t typ = 0; typ < 3; typ++)
                    tHit[d * 8 * 3 + ch * 3 + typ][i] = false;
            }
        }
    }

    if (nHits == 0)
        return;

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BLosTCalData* hit = dynamic_cast<R3BLosTCalData*>(fTCalItems->At(ihit));
        Int_t ch = hit->GetChannel();
        Int_t typ = hit->GetType();
        Int_t det = hit->GetDetector();
        Double_t rawTime = hit->GetRawTimeNs();
        Bool_t inHit = false;

        Int_t hitNo = 0, detNo = 0;
        for (Int_t d = 0; d < fNofDetectors; d++)
        {
            if ((det - 1) != d)
                continue;
            for (Int_t i = 0; i < hits[d]; i++)
            {
                if (typ == 0)
                {
                    Double_t t_ref = 0.;
                    t_ref = fTimeStitch->GetTime(rawTime - trigTime[typ], "vftx", "vftx");
                    if (fabs((avg_v[d][i] / v_npmt[d][i]) - t_ref) < fWindowV)
                    {
                        v_npmt[d][i]++;
                        if (tHit[(det - 1) * 8 * 3 + (ch - 1) * 3 + typ][i] || v_npmt[d][i] > 8)
                        {
                            LOG(warn) << "Possible pileup in LOS VFTX. Skipping event.";
                            return;
                        }
                        vTime[d * 8 + ch - 1][i] = rawTime;
                        avg_v[d][i] += fTimeStitch->GetTime(rawTime - trigTime[typ], "vftx", "vftx");
                        tHit[d * 8 * 3 + (ch - 1) * 3 + typ][i] = true;
                        inHit = true;
                    }
                }
                else if (typ == 1)
                {
                    if (v_npmt[d][i] != 8)
                        continue;
                    Double_t t_ref = 0., tdiff = 0.;
                    t_ref = fTimeStitch->GetTime(rawTime - trigTime[typ], "tamex", "tamex");
                    tdiff = fTimeStitch->GetTime(t_ref - avg_v[d][i] / 8., "tamex", "vftx");
                    if (tdiff > fLEMatchParams->GetAt(fNumParamsTamexLE * (ch - 1) + 0) &&
                        tdiff < fLEMatchParams->GetAt(fNumParamsTamexLE * (ch - 1) + 1) && tle_npmt[d][i] < 8)
                    {
                        tle_npmt[d][i]++;
                        if (tHit[(det - 1) * 8 * 3 + (ch - 1) * 3 + typ][i] || tle_npmt[d][i] > 8)
                        {
                            LOG(warn) << "Possible pileup in LOS TAMEX LE. Ignoring TAMEX hit.";
                            tle_npmt[d][i] = 32768;
                            break;
                        }
                        if (inHit)
                        {
                            LOG(warn) << "Possible pileup in LOS TAMEX LE. Ignoring TAMEX hit.";
                            tle_npmt[d][i] = 32768;
                            tle_npmt[detNo][hitNo] = 32768;
                            break;
                        }
                        TLeTime[d * 8 + ch - 1][i] = rawTime;
                        avg_t[d][i] += fTimeStitch->GetTime(rawTime - trigTime[typ], "tamex", "tamex");
                        tHit[d * 8 * 3 + (ch - 1) * 3 + typ][i] = true;
                        hitNo = i;
                        detNo = d;
                        inHit = true;
                    }
                }
                else if (typ == 2)
                {
                    if (tle_npmt[d][i] != 8)
                        continue;
                    Double_t t_ref = 0., tdiff = 0.;
                    t_ref = fTimeStitch->GetTime(rawTime - trigTime[1], "tamex", "tamex");
                    tdiff = fTimeStitch->GetTime(t_ref - avg_t[d][i] / 8., "tamex", "tamex");
                    if (tdiff > fTEMatchParams->GetAt(fNumParamsTamexTE * (ch - 1) + 0) &&
                        tdiff < fTEMatchParams->GetAt(fNumParamsTamexTE * (ch - 1) + 1) && tte_npmt[d][i] < 8)
                    {
                        tte_npmt[d][i]++;
                        if (tHit[(det - 1) * 8 * 3 + (ch - 1) * 3 + typ][i] || tte_npmt[d][i] > 8)
                        {
                            //   LOG(warn) << "Possible pileup in LOS TAMEX TE. Ignoring TAMEX hit."; commented to avoid
                            //   too many messages
                            tte_npmt[d][i] = 32768;
                            break;
                        }
                        if (inHit)
                        {
                            LOG(warn) << "Possible pileup in LOS TAMEX TE. Ignoring TAMEX hit.";
                            std::cout << "\nor here\t";
                            tte_npmt[d][i] = 32768;
                            tte_npmt[detNo][hitNo] = 32768;
                            break;
                        }
                        TTeTime[d * 8 + ch - 1][i] = rawTime;
                        tHit[d * 8 * 3 + (ch - 1) * 3 + typ][i] = true;
                        hitNo = i;
                        detNo = d;
                        inHit = true;
                    }
                }
            }
        }
        if (!inHit && typ == 0 && fNofDetectors > 0)
        {
            vTime[(det - 1) * 8 + ch - 1][hits[det - 1]] = rawTime;
            tHit[(det - 1) * 8 * 3 + (ch - 1) * 3 + typ][hits[det - 1]] = true;
            v_npmt[(det - 1)][hits[det - 1]]++;
            avg_v[(det - 1)][hits[det - 1]] += fTimeStitch->GetTime(rawTime - trigTime[typ], "vftx", "vftx");
            hits[det - 1]++;
        }
    }
    for (Int_t d = 0; d < fNofDetectors; d++)
    {
        for (Int_t i = 0; i < hits[d]; i++)
        {

            if (v_npmt[d][i] != 8)
                continue;
            Double_t t_hit = 0., Z = 0., x_cm = 0., y_cm = 0.;
            for (Int_t pm = 0; pm < 8; pm++)
            {
                t_hit += vTime[d * 8 + pm][i];
            }
            t_hit = t_hit / 8.;
            x_cm =
                ((vTime[d * 8 + 5][i] + vTime[d * 8 + 6][i]) / 2. - (vTime[d * 8 + 1][i] + vTime[d * 8 + 2][i]) / 2.) *
                (-1.);
            y_cm =
                ((vTime[d * 8 + 7][i] + vTime[d * 8 + 0][i]) / 2. - (vTime[d * 8 + 3][i] + vTime[d * 8 + 4][i]) / 2.) *
                (-1.);
            x_cm = (x_cm - flosOffsetX) * flosVeffX;
            y_cm = (y_cm - flosOffsetY) * flosVeffY;

            if (tle_npmt[d][i] == 8 && tte_npmt[d][i] == 8)
            {
                Double_t tot = 0.;
                for (Int_t pm = 0; pm < 8; pm++)
                {
                    Double_t tle_ref = 0., tte_ref = 0.;
                    tle_ref = fTimeStitch->GetTime(TLeTime[d * 8 + pm][i] - trigTime[1], "tamex", "tamex");
                    tte_ref = fTimeStitch->GetTime(TTeTime[d * 8 + pm][i] - trigTime[1], "tamex", "tamex");
                    tot += fTimeStitch->GetTime(tte_ref - tle_ref, "tamex", "tamex");
                }
                tot = tot / 8.;
                Z = tot * fp1 + fp0;
            }
            else
                Z = 0.;

            new ((*fHitItems)[fHitItems->GetEntriesFast()]) R3BLosHitData(d + 1, t_hit, x_cm, y_cm, Z);
        }
    }
}

void R3BLosTCal2Hit::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
}

Double_t R3BLosTCal2Hit::walk(Int_t inum, Double_t tot)
{

    Double_t y = 0. / 0., ysc = 0. / 0., term[8] = { 0. };
    Double_t x;

    x = tot;
    term[0] = x;
    for (Int_t i = 0; i < 7; i++)
    {
        term[i + 1] = term[i] * x;
    }

    ysc = walk_par[inum][2] + walk_par[inum][3] * term[0] + walk_par[inum][4] * term[1] + walk_par[inum][5] * term[2] +
          walk_par[inum][6] * term[3] + walk_par[inum][7] * term[4] + walk_par[inum][8] * term[5] +
          walk_par[inum][9] * term[6] + walk_par[inum][10] * term[7];

    if (tot < walk_par[inum][0] || tot > walk_par[inum][1])
        ysc = 0.0 / 0.0;

    return ysc;
}

Double_t R3BLosTCal2Hit::satu(Int_t inum, Double_t tot, Double_t dt)
{

    Double_t ysc = 0. / 0.;

    // if(tot_par[inum][0] > 0.)
    // ysc  = (tot_par[inum][0]*tot+tot_par[inum][1])/(tot_par[inum][2]-tot)*tot_par[inum][3] ;

    ysc = tot_par[inum][0] + tot_par[inum][1] * (1. - 1. / (exp((dt - tot_par[inum][2]) / tot_par[inum][3]) + 1.));
    ysc = tot / ysc;
    ysc = ysc * (tot_par[inum][0] + tot_par[inum][1]);

    return ysc;
}

ClassImp(R3BLosTCal2Hit);
