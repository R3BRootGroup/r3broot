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

/*
 * This task should fill histograms with MWPC online data correlations
 */

#include "R3BMusliVsMwpcOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BMusliCalData.h"
#include "R3BMusliHitData.h"
#include "R3BMusliMappedData.h"
#include "R3BMwpcHitData.h"
#include "R3BTGeoPar.h"
#include "THttpServer.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "R3BLogger.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"

#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;

R3BMusliVsMwpcOnlineSpectra::R3BMusliVsMwpcOnlineSpectra()
    : R3BMusliVsMwpcOnlineSpectra("MusicVsMwpcOnlineSpectra", 1, "Mwpc1", "Mwpc0")
{
}

R3BMusliVsMwpcOnlineSpectra::R3BMusliVsMwpcOnlineSpectra(const TString& name,
                                                         Int_t iVerbose,
                                                         const TString& namedet1,
                                                         const TString& namedet2)
    : FairTask(name, iVerbose)
    , fMapMusli(NULL)
    , fCalMusli(NULL)
    , fHitMusli(NULL)
    , fHitMwpcDet1(NULL)
    , fHitMwpcDet2(NULL)
    , fMw1Geo_Par(NULL)
    , fMw2Geo_Par(NULL)
    , fNameDet1(namedet1)
    , fNameDet2(namedet2)
    , fNEvents(0)
{
}

R3BMusliVsMwpcOnlineSpectra::~R3BMusliVsMwpcOnlineSpectra()
{
    LOG(info) << "R3BMusicVsMwpcOnlineSpectra::Delete instance";
    if (fMapMusli)
        delete fMapMusli;
    if (fCalMusli)
        delete fCalMusli;
    if (fHitMusli)
        delete fHitMusli;
    if (fHitMwpcDet1)
        delete fHitMwpcDet1;
    if (fHitMwpcDet2)
        delete fHitMwpcDet2;
}

void R3BMusliVsMwpcOnlineSpectra::SetParContainers()
{
    // Parameter Container

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "R3BMusliVsMwpcOnlineSpectra::SetParContainers(), FairRuntimeDb not found");

    fMw1Geo_Par = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer(fNameDet1 + "GeoPar"));
    R3BLOG_IF(error,
              !fMw1Geo_Par,
              "R3BMusliVsMwpcOnlineSpectra::SetParContainers() Couldn´t access to " + fNameDet1 + "GeoPar container.");

    fMw2Geo_Par = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer(fNameDet2 + "GeoPar"));
    R3BLOG_IF(error,
              !fMw2Geo_Par,
              "R3BMusliVsMwpcOnlineSpectra::SetParContainers() Couldn´t access to " + fNameDet2 + "GeoPar container.");
}

InitStatus R3BMusliVsMwpcOnlineSpectra::Init()
{

    LOG(info) << "R3B" + fNameDet1 + "vs" + fNameDet2 + "CorrelationOnlineSpectra::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BMusliVsMwpcOnlineSpectra::Init FairRootManager not found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to mapped data of the musli detectors
    fMapMusli = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliMappedData"));
    if (!fMapMusli)
    {
        return kFATAL;
    }

    // get access to cal data of the musli detectors
    fCalMusli = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliCalData"));
    if (!fCalMusli)
    {
        LOG(info) << "R3BMusliOnlineSpectra::Init() No MusliCalData found";
    }

    // get access to hit data of the musli detectors
    fHitMusli = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliHitData"));
    if (!fHitMusli)
    {
        LOG(info) << "R3BMusliOnlineSpectra::Init() No MusliHitData found";
    }

    // get access to hit data of mwpcs
    fHitMwpcDet1 = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet1 + "HitData"));
    if (!fHitMwpcDet1)
        LOG(warn) << "R3BMusliVsMwpcOnlineSpectra: " + fNameDet1 + "HitData not found";

    fHitMwpcDet2 = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet2 + "HitData"));
    if (!fHitMwpcDet2)
        LOG(warn) << "R3BMusliVsMwpcOnlineSpectra: " + fNameDet2 + "HitData not found";

    // Create histograms for detectors
    char Name1[255];
    char Name2[255];
    char* cstrNameDet1 = new char[fNameDet1.Length() + 1];
    char* cstrNameDet2 = new char[fNameDet2.Length() + 1];
    strcpy(cstrNameDet1, fNameDet1);
    strcpy(cstrNameDet2, fNameDet2);

    // === ======================== === //
    // === MAPPED-MUSLI vs HIT-MWPC === //
    // === ======================== === //

    cMusliMap_DTvsX1 = new TCanvas("mapMusVsMw_DTvsX1", "", 10, 10, 800, 700);
    cMusliMap_DTvsX1->Divide(5, 3);
    cMusliMap_DTvsX2 = new TCanvas("mapMusVsMw_DTvsX2", "", 10, 10, 800, 700);
    cMusliMap_DTvsX2->Divide(5, 3);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh2_MusMapDT_vs_X1_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "DTmap A%02d and A%02d vs X1", 2 * j + 1, 2 * j + 2);
        fh2_mapDTvsX1[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_mapDTvsX1[j]->GetXaxis()->SetTitle(Name1);
        fh2_mapDTvsX1[j]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_mapDTvsX1[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_mapDTvsX1[j]->GetXaxis()->CenterTitle(true);
        fh2_mapDTvsX1[j]->GetYaxis()->CenterTitle(true);
        fh2_mapDTvsX1[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX1[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX1[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX1[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX1[j]->SetFillColor(31);
        cMusliMap_DTvsX1->cd(j + 1);
        fh2_mapDTvsX1[j]->Draw("COL");

        sprintf(Name1, "fh2_MusMapDT_vs_X2_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "DTmap A%02d and A%02d vs X2", 2 * j + 1, 2 * j + 2);
        fh2_mapDTvsX2[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_mapDTvsX2[j]->GetXaxis()->SetTitle(Name2);
        fh2_mapDTvsX2[j]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_mapDTvsX2[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_mapDTvsX2[j]->GetXaxis()->CenterTitle(true);
        fh2_mapDTvsX2[j]->GetYaxis()->CenterTitle(true);
        fh2_mapDTvsX2[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX2[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX2[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX2[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX2[j]->SetFillColor(31);
        cMusliMap_DTvsX2->cd(j + 1);
        fh2_mapDTvsX2[j]->Draw("COL");
    }
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh2_MusMapDT_vs_X1_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "DTmap A%02d to A%02d vs X1", j * 4 + 1, j * 4 + 5);
        fh2_mapDTvsX1[j + 8] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_mapDTvsX1[j + 8]->GetXaxis()->SetTitle(Name1);
        fh2_mapDTvsX1[j + 8]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_mapDTvsX1[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh2_mapDTvsX1[j + 8]->GetXaxis()->CenterTitle(true);
        fh2_mapDTvsX1[j + 8]->GetYaxis()->CenterTitle(true);
        fh2_mapDTvsX1[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX1[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX1[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX1[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX1[j + 8]->SetFillColor(31);
        cMusliMap_DTvsX1->cd(j + 9);
        fh2_mapDTvsX1[j + 8]->Draw("COL");

        sprintf(Name1, "fh2_MusMapDT_vs_X2_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "DTmap A%02d to A%02d vs X2", j * 4 + 1, j * 4 + 5);
        fh2_mapDTvsX2[j + 8] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_mapDTvsX2[j + 8]->GetXaxis()->SetTitle(Name2);
        fh2_mapDTvsX2[j + 8]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_mapDTvsX2[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh2_mapDTvsX2[j + 8]->GetXaxis()->CenterTitle(true);
        fh2_mapDTvsX2[j + 8]->GetYaxis()->CenterTitle(true);
        fh2_mapDTvsX2[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX2[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX2[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX2[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX2[j + 8]->SetFillColor(31);
        cMusliMap_DTvsX2->cd(j + 9);
        fh2_mapDTvsX2[j + 8]->Draw("COL");
    }
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh2_MusMapDT_vs_X1_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "DTmap A%02d to A%02d vs X1", j * 8 + 1, j * 8 + 8);
        fh2_mapDTvsX1[j + 12] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_mapDTvsX1[j + 12]->GetXaxis()->SetTitle(Name1);
        fh2_mapDTvsX1[j + 12]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_mapDTvsX1[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh2_mapDTvsX1[j + 12]->GetXaxis()->CenterTitle(true);
        fh2_mapDTvsX1[j + 12]->GetYaxis()->CenterTitle(true);
        fh2_mapDTvsX1[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX1[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX1[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX1[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX1[j + 12]->SetFillColor(31);
        cMusliMap_DTvsX1->cd(j + 13);
        fh2_mapDTvsX1[j + 12]->Draw("COL");

        sprintf(Name1, "fh2_MusMapDT_vs_X2_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "DTmap A%02d to A%02d vs X2", j * 8 + 1, j * 8 + 8);
        fh2_mapDTvsX2[j + 12] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_mapDTvsX2[j + 12]->GetXaxis()->SetTitle(Name2);
        fh2_mapDTvsX2[j + 12]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_mapDTvsX2[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh2_mapDTvsX2[j + 12]->GetXaxis()->CenterTitle(true);
        fh2_mapDTvsX2[j + 12]->GetYaxis()->CenterTitle(true);
        fh2_mapDTvsX2[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX2[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX2[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh2_mapDTvsX2[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh2_mapDTvsX2[j + 12]->SetFillColor(31);
        cMusliMap_DTvsX2->cd(j + 13);
        fh2_mapDTvsX2[j + 12]->Draw("COL");
    }
    sprintf(Name1, "fh2_MusMapDT_vs_X1_a1_to_a16");
    sprintf(Name2, "DTmap A01 to A16 vs X1");
    fh2_mapDTvsX1[14] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
    sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
    fh2_mapDTvsX1[14]->GetXaxis()->SetTitle(Name1);
    fh2_mapDTvsX1[14]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
    fh2_mapDTvsX1[14]->GetYaxis()->SetTitleOffset(1.1);
    fh2_mapDTvsX1[14]->GetXaxis()->CenterTitle(true);
    fh2_mapDTvsX1[14]->GetYaxis()->CenterTitle(true);
    fh2_mapDTvsX1[14]->GetXaxis()->SetLabelSize(0.045);
    fh2_mapDTvsX1[14]->GetXaxis()->SetTitleSize(0.045);
    fh2_mapDTvsX1[14]->GetYaxis()->SetLabelSize(0.045);
    fh2_mapDTvsX1[14]->GetYaxis()->SetTitleSize(0.045);
    fh2_mapDTvsX1[14]->SetFillColor(31);
    cMusliMap_DTvsX1->cd(15);
    fh2_mapDTvsX1[14]->Draw("COL");

    sprintf(Name1, "fh2_MusMapDT_vs_X2_a1_to_a16");
    sprintf(Name2, "DTmap A01 to A16 vs X2");
    fh2_mapDTvsX2[14] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, 10000, 30000);
    sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
    fh2_mapDTvsX2[14]->GetXaxis()->SetTitle(Name2);
    fh2_mapDTvsX2[14]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
    fh2_mapDTvsX2[14]->GetYaxis()->SetTitleOffset(1.1);
    fh2_mapDTvsX2[14]->GetXaxis()->CenterTitle(true);
    fh2_mapDTvsX2[14]->GetYaxis()->CenterTitle(true);
    fh2_mapDTvsX2[14]->GetXaxis()->SetLabelSize(0.045);
    fh2_mapDTvsX2[14]->GetXaxis()->SetTitleSize(0.045);
    fh2_mapDTvsX2[14]->GetYaxis()->SetLabelSize(0.045);
    fh2_mapDTvsX2[14]->GetYaxis()->SetTitleSize(0.045);
    fh2_mapDTvsX2[14]->SetFillColor(31);
    cMusliMap_DTvsX2->cd(15);
    fh2_mapDTvsX2[14]->Draw("COL");

    // === ===================== === //
    // === CAL-MUSLI vs HIT-MWPC === //
    // === ===================== === //

    cMusliCal_DTvsX1 = new TCanvas("calMusVsMw_DTvsX1", "", 10, 10, 800, 700);
    cMusliCal_DTvsX1->Divide(5, 3);
    cMusliCal_DTvsX2 = new TCanvas("calMusVsMw_DTvsX2", "", 10, 10, 800, 700);
    cMusliCal_DTvsX2->Divide(5, 3);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh2_MusCalDT_vs_X1_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "DT A%02d and A%02d vs X1", 2 * j + 1, 2 * j + 2);
        fh2_calDTvsX1[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_calDTvsX1[j]->GetXaxis()->SetTitle(Name1);
        fh2_calDTvsX1[j]->GetYaxis()->SetTitle("Calibrated Drift Time [mm]");
        fh2_calDTvsX1[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_calDTvsX1[j]->GetXaxis()->CenterTitle(true);
        fh2_calDTvsX1[j]->GetYaxis()->CenterTitle(true);
        fh2_calDTvsX1[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_calDTvsX1[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_calDTvsX1[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_calDTvsX1[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_calDTvsX1[j]->SetFillColor(31);
        cMusliCal_DTvsX1->cd(j + 1);
        fh2_calDTvsX1[j]->Draw("COL");

        sprintf(Name1, "fh2_MusCalDT_vs_X2_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "DT A%02d and A%02d vs X2", 2 * j + 1, 2 * j + 2);
        fh2_calDTvsX2[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_calDTvsX2[j]->GetXaxis()->SetTitle(Name2);
        fh2_calDTvsX2[j]->GetYaxis()->SetTitle("Drift Time [100ps]");
        fh2_calDTvsX2[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_calDTvsX2[j]->GetXaxis()->CenterTitle(true);
        fh2_calDTvsX2[j]->GetYaxis()->CenterTitle(true);
        fh2_calDTvsX2[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_calDTvsX2[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_calDTvsX2[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_calDTvsX2[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_calDTvsX2[j]->SetFillColor(31);
        cMusliCal_DTvsX2->cd(j + 1);
        fh2_calDTvsX2[j]->Draw("COL");
    }
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh2_MusCalDT_vs_X1_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "DTcal A%02d to A%02d vs X1", j * 4 + 1, j * 4 + 5);
        fh2_calDTvsX1[j + 8] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_calDTvsX1[j + 8]->GetXaxis()->SetTitle(Name1);
        fh2_calDTvsX1[j + 8]->GetYaxis()->SetTitle("Calibrated Drift Time [mm]");
        fh2_calDTvsX1[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh2_calDTvsX1[j + 8]->GetXaxis()->CenterTitle(true);
        fh2_calDTvsX1[j + 8]->GetYaxis()->CenterTitle(true);
        fh2_calDTvsX1[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh2_calDTvsX1[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh2_calDTvsX1[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh2_calDTvsX1[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh2_calDTvsX1[j + 8]->SetFillColor(31);
        cMusliCal_DTvsX1->cd(j + 9);
        fh2_calDTvsX1[j + 8]->Draw("COL");

        sprintf(Name1, "fh2_MusCalDT_vs_X2_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "DTcal A%02d to A%02d vs X2", j * 4 + 1, j * 4 + 5);
        fh2_calDTvsX2[j + 8] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_calDTvsX2[j + 8]->GetXaxis()->SetTitle(Name2);
        fh2_calDTvsX2[j + 8]->GetYaxis()->SetTitle("Calibrated Drift Time [mm]");
        fh2_calDTvsX2[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh2_calDTvsX2[j + 8]->GetXaxis()->CenterTitle(true);
        fh2_calDTvsX2[j + 8]->GetYaxis()->CenterTitle(true);
        fh2_calDTvsX2[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh2_calDTvsX2[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh2_calDTvsX2[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh2_calDTvsX2[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh2_calDTvsX2[j + 8]->SetFillColor(31);
        cMusliCal_DTvsX2->cd(j + 9);
        fh2_calDTvsX2[j + 8]->Draw("COL");
    }
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh2_MusCalDT_vs_X1_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "DTcal A%02d to A%02d vs X1", j * 8 + 1, j * 8 + 8);
        fh2_calDTvsX1[j + 12] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_calDTvsX1[j + 12]->GetXaxis()->SetTitle(Name1);
        fh2_calDTvsX1[j + 12]->GetYaxis()->SetTitle("Calibrated Drift Time [mm]");
        fh2_calDTvsX1[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh2_calDTvsX1[j + 12]->GetXaxis()->CenterTitle(true);
        fh2_calDTvsX1[j + 12]->GetYaxis()->CenterTitle(true);
        fh2_calDTvsX1[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh2_calDTvsX1[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh2_calDTvsX1[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh2_calDTvsX1[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh2_calDTvsX1[j + 12]->SetFillColor(31);
        cMusliCal_DTvsX1->cd(j + 13);
        fh2_calDTvsX1[j + 12]->Draw("COL");

        sprintf(Name1, "fh2_MusCalDT_vs_X2_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "DTcal A%02d to A%02d vs X2", j * 8 + 1, j * 8 + 8);
        fh2_calDTvsX2[j + 12] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_calDTvsX2[j + 12]->GetXaxis()->SetTitle(Name2);
        fh2_calDTvsX2[j + 12]->GetYaxis()->SetTitle("Calibrated Drift Time [mm]");
        fh2_calDTvsX2[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh2_calDTvsX2[j + 12]->GetXaxis()->CenterTitle(true);
        fh2_calDTvsX2[j + 12]->GetYaxis()->CenterTitle(true);
        fh2_calDTvsX2[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh2_calDTvsX2[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh2_calDTvsX2[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh2_calDTvsX2[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh2_calDTvsX2[j + 12]->SetFillColor(31);
        cMusliCal_DTvsX2->cd(j + 13);
        fh2_calDTvsX2[j + 12]->Draw("COL");
    }
    sprintf(Name1, "fh2_MusCalDT_vs_X1_a1_to_a16");
    sprintf(Name2, "DTcal A01 to A16 vs X1");
    fh2_calDTvsX1[14] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
    sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
    fh2_calDTvsX1[14]->GetXaxis()->SetTitle(Name1);
    fh2_calDTvsX1[14]->GetYaxis()->SetTitle("Calibrated Drift Time [mm]");
    fh2_calDTvsX1[14]->GetYaxis()->SetTitleOffset(1.1);
    fh2_calDTvsX1[14]->GetXaxis()->CenterTitle(true);
    fh2_calDTvsX1[14]->GetYaxis()->CenterTitle(true);
    fh2_calDTvsX1[14]->GetXaxis()->SetLabelSize(0.045);
    fh2_calDTvsX1[14]->GetXaxis()->SetTitleSize(0.045);
    fh2_calDTvsX1[14]->GetYaxis()->SetLabelSize(0.045);
    fh2_calDTvsX1[14]->GetYaxis()->SetTitleSize(0.045);
    fh2_calDTvsX1[14]->SetFillColor(31);
    cMusliCal_DTvsX1->cd(15);
    fh2_calDTvsX1[14]->Draw("COL");

    sprintf(Name1, "fh2_MusCalDT_vs_X2_a1_to_a16");
    sprintf(Name2, "DTcal A01 to A16 vs X2");
    fh2_calDTvsX2[14] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
    sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
    fh2_calDTvsX2[14]->GetXaxis()->SetTitle(Name2);
    fh2_calDTvsX2[14]->GetYaxis()->SetTitle("Calibrated Drift Time [mm]");
    fh2_calDTvsX2[14]->GetYaxis()->SetTitleOffset(1.1);
    fh2_calDTvsX2[14]->GetXaxis()->CenterTitle(true);
    fh2_calDTvsX2[14]->GetYaxis()->CenterTitle(true);
    fh2_calDTvsX2[14]->GetXaxis()->SetLabelSize(0.045);
    fh2_calDTvsX2[14]->GetXaxis()->SetTitleSize(0.045);
    fh2_calDTvsX2[14]->GetYaxis()->SetLabelSize(0.045);
    fh2_calDTvsX2[14]->GetYaxis()->SetTitleSize(0.045);
    fh2_calDTvsX2[14]->SetFillColor(31);
    cMusliCal_DTvsX2->cd(15);
    fh2_calDTvsX2[14]->Draw("COL");

    // === ===================== === //
    // === HIT-MUSLI vs HIT-MWPC === //
    // === ===================== === //

    cMusliHit_EvsX1 = new TCanvas("hitMusVsMw_EvsX1", "", 10, 10, 800, 700);
    cMusliHit_EvsX1->Divide(2, 2);

    cMusliHit_ZvsX1 = new TCanvas("hitMusVsMw_ZvsX1", "", 10, 10, 800, 700);
    cMusliHit_ZvsX1->Divide(2, 2);

    cMusliHit_XvsX1 = new TCanvas("hitMusVsMw_XvsX1", "", 10, 10, 800, 700);
    cMusliHit_XvsX1->Divide(2, 2);

    cMusliHit_TvsX1 = new TCanvas("hitMusVsMw_THETAvsX1", "", 10, 10, 800, 700);

    cMusliHit_EvsX2 = new TCanvas("hitMusVsMw_EvsX2", "", 10, 10, 800, 700);
    cMusliHit_EvsX2->Divide(2, 2);

    cMusliHit_ZvsX2 = new TCanvas("hitMusVsMw_ZvsX2", "", 10, 10, 800, 700);
    cMusliHit_ZvsX2->Divide(2, 2);

    cMusliHit_XvsX2 = new TCanvas("hitMusVsMw_XvsX2", "", 10, 10, 800, 700);
    cMusliHit_XvsX2->Divide(2, 2);

    cMusliHit_TvsX2 = new TCanvas("hitMusVsMw_THETAvsX2", "", 10, 10, 800, 700);

    Int_t nb_anodes[4] = { 2, 4, 8, 16 };

    for (Int_t j = 0; j < 4; j++)
    {

        sprintf(Name1, "fh2_MusHitE%i_vs_X1", nb_anodes[j]);
        sprintf(Name2, " Average E for data with nb anodes %02d  vs X1", nb_anodes[j]);
        fh2_hitEvsX1[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1600, 0, 64000);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_hitEvsX1[j]->GetXaxis()->SetTitle(Name1);
        fh2_hitEvsX1[j]->GetYaxis()->SetTitle("Average Energy [channels]");
        fh2_hitEvsX1[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_hitEvsX1[j]->GetXaxis()->CenterTitle(true);
        fh2_hitEvsX1[j]->GetYaxis()->CenterTitle(true);
        fh2_hitEvsX1[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_hitEvsX1[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_hitEvsX1[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_hitEvsX1[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_hitEvsX1[j]->SetFillColor(31);
        cMusliHit_EvsX1->cd(j + 1);
        fh2_hitEvsX1[j]->Draw("COL");

        sprintf(Name1, "fh2_MusHitE%i_vs_X2", nb_anodes[j]);
        sprintf(Name2, " Average E for data with nb anodes %02d  vs X2", nb_anodes[j]);
        fh2_hitEvsX2[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1600, 0, 64000);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_hitEvsX2[j]->GetXaxis()->SetTitle("Position of the second Mwpc [mm]");
        fh2_hitEvsX2[j]->GetYaxis()->SetTitle("Average Energy [channels]");
        fh2_hitEvsX2[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_hitEvsX2[j]->GetXaxis()->CenterTitle(true);
        fh2_hitEvsX2[j]->GetYaxis()->CenterTitle(true);
        fh2_hitEvsX2[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_hitEvsX2[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_hitEvsX2[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_hitEvsX2[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_hitEvsX2[j]->SetFillColor(31);
        cMusliHit_EvsX2->cd(j + 1);
        fh2_hitEvsX2[j]->Draw("COL");

        sprintf(Name1, "fh2_MusHitZ%i_vs_X1", nb_anodes[j]);
        sprintf(Name2, " Z charge for data with nb anodes %02d  vs X1", nb_anodes[j]);
        fh2_hitZvsX1[j] = new TH2D(Name1, Name2, 1000, -100, 100, 800, 2, 10);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_hitZvsX1[j]->GetXaxis()->SetTitle(Name1);
        fh2_hitZvsX1[j]->GetYaxis()->SetTitle("Atomic charge (Z)");
        fh2_hitZvsX1[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_hitZvsX1[j]->GetXaxis()->CenterTitle(true);
        fh2_hitZvsX1[j]->GetYaxis()->CenterTitle(true);
        fh2_hitZvsX1[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_hitZvsX1[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_hitZvsX1[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_hitZvsX1[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_hitZvsX1[j]->SetFillColor(31);
        cMusliHit_ZvsX1->cd(j + 1);
        fh2_hitZvsX1[j]->Draw("COL");

        sprintf(Name1, "fh2_MusHitZ%i_vs_X2", nb_anodes[j]);
        sprintf(Name2, " Z charge for data with nb anodes %02d  vs X2", nb_anodes[j]);
        fh2_hitZvsX2[j] = new TH2D(Name1, Name2, 1000, -100, 100, 800, 2, 10);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_hitZvsX2[j]->GetXaxis()->SetTitle("Position of the second Mwpc [mm]");
        fh2_hitZvsX2[j]->GetYaxis()->SetTitle("Atomic charge (Z)");
        fh2_hitZvsX2[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_hitZvsX2[j]->GetXaxis()->CenterTitle(true);
        fh2_hitZvsX2[j]->GetYaxis()->CenterTitle(true);
        fh2_hitZvsX2[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_hitZvsX2[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_hitZvsX2[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_hitZvsX2[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_hitZvsX2[j]->SetFillColor(31);
        cMusliHit_ZvsX2->cd(j + 1);
        fh2_hitZvsX2[j]->Draw("COL");

        sprintf(Name1, "fh2_MusHitX%i_vs_X1", nb_anodes[j]);
        sprintf(Name2, " X middle of MusLI with nb anodes %02d  vs X1", nb_anodes[j]);
        fh2_hitXvsX1[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
        fh2_hitXvsX1[j]->GetXaxis()->SetTitle(Name1);
        fh2_hitXvsX1[j]->GetYaxis()->SetTitle("Rough position in the middle of MusLI [mm]");
        fh2_hitXvsX1[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_hitXvsX1[j]->GetXaxis()->CenterTitle(true);
        fh2_hitXvsX1[j]->GetYaxis()->CenterTitle(true);
        fh2_hitXvsX1[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_hitXvsX1[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_hitXvsX1[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_hitXvsX1[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_hitXvsX1[j]->SetFillColor(31);
        cMusliHit_XvsX1->cd(j + 1);
        fh2_hitXvsX1[j]->Draw("COL");

        sprintf(Name1, "fh2_MusHitX%i_vs_X2", nb_anodes[j]);
        sprintf(Name2, " X middle of MusLI for data with nb anodes %02d  vs X2", nb_anodes[j]);
        fh2_hitXvsX2[j] = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -100, 100);
        sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
        fh2_hitXvsX2[j]->GetXaxis()->SetTitle("Position of the second Mwpc [mm]");
        fh2_hitXvsX2[j]->GetYaxis()->SetTitle("Rough position in the middle of MusLI [mm]");
        fh2_hitXvsX2[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_hitXvsX2[j]->GetXaxis()->CenterTitle(true);
        fh2_hitXvsX2[j]->GetYaxis()->CenterTitle(true);
        fh2_hitXvsX2[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_hitXvsX2[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_hitXvsX2[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_hitXvsX2[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_hitXvsX2[j]->SetFillColor(31);
        cMusliHit_XvsX2->cd(j + 1);
        fh2_hitXvsX2[j]->Draw("COL");
    }

    sprintf(Name1, "fh2_MusHitTheta_vs_X1");
    sprintf(Name2, " Theta from MusLI  vs X1");
    fh2_hitTvsX1 = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -0.1, 0.1);
    sprintf(Name1, "X position in the first Mwpc [mm] (%s)", cstrNameDet1);
    fh2_hitTvsX1->GetXaxis()->SetTitle(Name1);
    fh2_hitTvsX1->GetYaxis()->SetTitle("Theta from MusLI [mrad]");
    fh2_hitTvsX1->GetYaxis()->SetTitleOffset(1.1);
    fh2_hitTvsX1->GetXaxis()->CenterTitle(true);
    fh2_hitTvsX1->GetYaxis()->CenterTitle(true);
    fh2_hitTvsX1->GetXaxis()->SetLabelSize(0.045);
    fh2_hitTvsX1->GetXaxis()->SetTitleSize(0.045);
    fh2_hitTvsX1->GetYaxis()->SetLabelSize(0.045);
    fh2_hitTvsX1->GetYaxis()->SetTitleSize(0.045);
    fh2_hitTvsX1->SetFillColor(31);
    cMusliHit_TvsX1->cd();
    fh2_hitTvsX1->Draw("COL");

    sprintf(Name1, "fh2_MusHitTheta_vs_X2");
    sprintf(Name2, " Theta from MusLI  vs X2");
    fh2_hitTvsX2 = new TH2D(Name1, Name2, 1000, -100, 100, 1000, -0.1, 0.1);
    sprintf(Name2, "X position in the second Mwpc [mm] (%s)", cstrNameDet2);
    fh2_hitTvsX2->GetXaxis()->SetTitle(Name2);
    fh2_hitTvsX2->GetYaxis()->SetTitle("Theta from MusLI [mrad]");
    fh2_hitTvsX2->GetYaxis()->SetTitleOffset(1.1);
    fh2_hitTvsX2->GetXaxis()->CenterTitle(true);
    fh2_hitTvsX2->GetYaxis()->CenterTitle(true);
    fh2_hitTvsX2->GetXaxis()->SetLabelSize(0.045);
    fh2_hitTvsX2->GetXaxis()->SetTitleSize(0.045);
    fh2_hitTvsX2->GetYaxis()->SetLabelSize(0.045);
    fh2_hitTvsX2->GetYaxis()->SetTitleSize(0.045);
    fh2_hitTvsX2->SetFillColor(31);
    cMusliHit_TvsX2->cd();
    fh2_hitTvsX2->Draw("COL");

    // === ======================= === //
    // === MAIN FOLDER-MusliVsMwpc === //
    // === ======================= === //

    TFolder* mainfolMusliVsMw = new TFolder("MusliVsMw", "Music Light Ion vs Mwpc info");

    mainfolMusliVsMw->Add(cMusliMap_DTvsX1);
    mainfolMusliVsMw->Add(cMusliMap_DTvsX2);

    mainfolMusliVsMw->Add(cMusliCal_DTvsX1);
    mainfolMusliVsMw->Add(cMusliCal_DTvsX2);

    mainfolMusliVsMw->Add(cMusliHit_EvsX1);
    mainfolMusliVsMw->Add(cMusliHit_EvsX2);

    mainfolMusliVsMw->Add(cMusliHit_ZvsX1);
    mainfolMusliVsMw->Add(cMusliHit_ZvsX2);

    mainfolMusliVsMw->Add(cMusliHit_XvsX1);
    mainfolMusliVsMw->Add(cMusliHit_XvsX2);

    mainfolMusliVsMw->Add(cMusliHit_TvsX1);
    mainfolMusliVsMw->Add(cMusliHit_TvsX2);

    run->AddObject(mainfolMusliVsMw);

    // === ==================================== === //
    // === Register command to reset histograms === //
    // === ==================================== === //

    run->GetHttpServer()->RegisterCommand("Reset_MUSLIvsMWPC_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BMusliVsMwpcOnlineSpectra::Reset_Histo()
{

    LOG(info) << "R3BMusliVsMwpcOnlineSpectra::Reset_Histo";
    for (Int_t j = 0; j < 15; j++)
    {
        fh2_mapDTvsX1[j]->Reset();
        fh2_mapDTvsX2[j]->Reset();
        fh2_calDTvsX1[j]->Reset();
        fh2_calDTvsX2[j]->Reset();
    }
    for (Int_t j = 0; j < 4; j++)
    {
        fh2_hitEvsX1[j]->Reset();
        fh2_hitEvsX2[j]->Reset();
        fh2_hitZvsX1[j]->Reset();
        fh2_hitZvsX2[j]->Reset();
        fh2_hitXvsX1[j]->Reset();
        fh2_hitXvsX2[j]->Reset();
    }
    fh2_hitTvsX1->Reset();
    fh2_hitTvsX2->Reset();
}

void R3BMusliVsMwpcOnlineSpectra::Exec(Option_t* option)
{
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BMusliVsMwpcOnlineSpectra::Exec FairRootManager not found";

    Int_t nHits;
    UInt_t rank;

    // --- --------------- --- //
    // --- HIT MWPC 1 DATA --- //
    // --- --------------- --- //
    Double_t fX1 = -100.;
    if (fHitMwpcDet1 && fHitMwpcDet1->GetEntriesFast() == 1)
    {
        R3BMwpcHitData* hitDataMw1 = dynamic_cast<R3BMwpcHitData*>(fHitMwpcDet1->At(0));
        fX1 = hitDataMw1->GetX() + fMw1Geo_Par->GetPosX();
    }

    // --- --------------- --- //
    // --- HIT MWPC 2 DATA --- //
    // --- --------------- --- //
    Double_t fX2 = -100.;
    if (fHitMwpcDet2 && fHitMwpcDet2->GetEntriesFast() == 1)
    {
        R3BMwpcHitData* hitDataMw2 = dynamic_cast<R3BMwpcHitData*>(fHitMwpcDet2->At(0));
        fX2 = hitDataMw2->GetX() + fMw2Geo_Par->GetPosX();
    }

    // --- ----------------- --- //
    // --- MAPPED MUSLI DATA --- //
    // --- ----------------- --- //
    if (fMapMusli && fMapMusli->GetEntriesFast() > 0)
    {
        // Local variables for mapped data
        UInt_t mult[18], t[18];
        UInt_t tref = 0;
        for (UShort_t i = 0; i < 15; i++)
        {
            mult[i] = 0;
            t[i] = 0;
        }
        for (UShort_t i = 15; i < 18; i++)
        {
            mult[i] = 0;
            t[i] = 0;
        }
        // read mapped data
        nHits = fMapMusli->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusliMappedData* hit = dynamic_cast<R3BMusliMappedData*>(fMapMusli->At(ihit));
            if (!hit)
                continue;
            rank = hit->GetSignal() - 1;
            mult[rank]++;
            t[rank] = hit->GetTime();
        }
        // calculate drift time and fill mapped data histo
        if (mult[16] == 1)
            tref = t[16];
        for (UShort_t i = 0; i < 15; i++)
        {
            if (mult[i] == 1 && tref > 0)
            {
                fh2_mapDTvsX1[i]->Fill(fX1, t[i] - tref);
                fh2_mapDTvsX2[i]->Fill(fX2, t[i] - tref);
            }
        }

    } // end of mapped data

    // --- -------------- --- //
    // --- CAL MUSLI DATA --- //
    // --- -------------- --- //
    if (fCalMusli && fCalMusli->GetEntriesFast() > 0)
    {
        // Local variables for cal data
        Double_t multcal[15], dtcal[18];
        for (UShort_t i = 0; i < 15; i++)
        {
            multcal[i] = 0;
            dtcal[i] = 0;
        }
        // read cal data
        nHits = fCalMusli->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusliCalData* hit = dynamic_cast<R3BMusliCalData*>(fCalMusli->At(ihit));
            if (!hit)
                continue;
            rank = hit->GetSignal() - 1;
            multcal[rank]++;
            dtcal[rank] = hit->GetDT();
        }
        // implement the spectra only if mult[signal] == 1
        for (Int_t i = 0; i < 15; i++)
        {
            if (multcal[i] == 1)
            {
                fh2_calDTvsX1[i]->Fill(fX1, dtcal[i]);
                fh2_calDTvsX2[i]->Fill(fX2, dtcal[i]);
            }
        }
    } // end of cal data

    // --- -------------- --- //
    // --- HIT MUSLI DATA --- //
    // --- -------------- --- //
    if (fHitMusli && fHitMusli->GetEntriesFast() > 0)
    {
        // Local variables for hit data
        Double_t multhit[4], ehit[4], zhit[4], xhit[4], theta[4];
        for (UShort_t i = 0; i < 4; i++)
        {
            multhit[i] = 0;
            ehit[i] = 0;
            zhit[i] = 0;
            xhit[i] = -1000.;
            theta[i] = -1000.;
        }
        // read  data
        nHits = fHitMusli->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusliHitData* hit = dynamic_cast<R3BMusliHitData*>(fHitMusli->At(ihit));
            if (!hit)
                continue;
            rank = hit->GetType() - 1;
            multhit[rank]++;
            ehit[rank] = hit->GetEave();
            zhit[rank] = hit->GetZcharge();
            xhit[rank] = hit->GetX();
            theta[rank] = hit->GetTheta();
        }
        // implement the spectra only if mult[signal] == 1
        for (Int_t i = 0; i < 4; i++)
        {
            if (multhit[i] == 1)
            {
                fh2_hitEvsX1[i]->Fill(fX1, ehit[i]);
                fh2_hitEvsX2[i]->Fill(fX2, ehit[i]);
                fh2_hitZvsX1[i]->Fill(fX1, zhit[i]);
                fh2_hitZvsX2[i]->Fill(fX2, zhit[i]);
                fh2_hitXvsX1[i]->Fill(fX1, xhit[i]);
                fh2_hitXvsX2[i]->Fill(fX2, xhit[i]);
                if (i == 0)
                {
                    fh2_hitTvsX1->Fill(fX1, theta[i]);
                    fh2_hitTvsX2->Fill(fX2, theta[i]);
                }
            }
        }

    } // end of hit data

    fNEvents += 1;
}

void R3BMusliVsMwpcOnlineSpectra::FinishEvent()
{
    if (fMapMusli)
    {
        fMapMusli->Clear();
    }
    if (fCalMusli)
    {
        fCalMusli->Clear();
    }
    if (fHitMusli)
    {
        fHitMusli->Clear();
    }
    if (fHitMwpcDet1)
    {
        fHitMwpcDet1->Clear();
    }
    if (fHitMwpcDet2)
    {
        fHitMwpcDet2->Clear();
    }
}

void R3BMusliVsMwpcOnlineSpectra::FinishTask()
{

    if (fHitMwpcDet1 && fHitMwpcDet2)
    {

        if (fMapMusli)
        {
            cMusliMap_DTvsX1->Clear();
            cMusliMap_DTvsX2->Clear();
        }

        if (fCalMusli)
        {
            cMusliCal_DTvsX1->Clear();
            cMusliCal_DTvsX2->Clear();
        }

        if (fHitMusli)
        {
            cMusliHit_EvsX1->Clear();
            cMusliHit_ZvsX1->Clear();
            cMusliHit_XvsX1->Clear();
            cMusliHit_TvsX1->Clear();

            cMusliHit_EvsX2->Clear();
            cMusliHit_ZvsX2->Clear();
            cMusliHit_XvsX2->Clear();
            cMusliHit_TvsX2->Clear();
        }
    }
}

ClassImp(R3BMusliVsMwpcOnlineSpectra);
