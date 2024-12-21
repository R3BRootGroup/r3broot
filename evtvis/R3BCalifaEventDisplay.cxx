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

// --------------------------------------------------------------------
// -----            R3BCalifaEventDisplay source file             -----
// --------------------------------------------------------------------

#include "R3BCalifaEventDisplay.h"

#include "FairEventManager.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairTask.h"
#include "TClonesArray.h"
#include "TObjArray.h"

#include "R3BCalifaCrystalCalData.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TMath.h"
#include "TVector3.h"

#include "TEveBrowser.h"
#include "TEveCalo.h"
#include "TEveCaloData.h"
#include "TEveCaloLegoOverlay.h"
#include "TEveLegoEventHandler.h"
#include "TEveManager.h"
#include "TEveProjectionAxes.h"
#include "TEveProjectionManager.h"
#include "TEveProjections.h"
#include "TEveScene.h"
#include "TEveTrans.h"
#include "TEveViewer.h"
#include "TEveWindow.h"
#include "TGLSAViewer.h"
#include "TGLViewer.h"
#include "TGTab.h"

#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include <iostream>

#include <assert.h>

using std::cout;
using std::endl;

// -----   Default constructor   -------------------------------------------
R3BCalifaEventDisplay::R3BCalifaEventDisplay()
    : FairTask("R3BCalifaEventDisplay", 0)
    , fEventManager(NULL)
    , fGeometryVersion(10)
    , fDataHist(NULL)
    , fCalo3d(NULL)
    , fCalo2d(NULL)
    , fLego(NULL)
    , slotLeftTop(NULL)
    , slotLeftBottom(NULL)
    , slotRightTop(NULL)
    , slotRightBottom(NULL)
    , viewerLeftTop(NULL)
    , sceneLeftTop(NULL)
    , viewerRightTop(NULL)
    , sceneRightTop(NULL)
    , viewerLeftBottom(NULL)
    , sceneLeftBottom(NULL)
    , viewerRightBottom(NULL)
    , sceneRightBottom(NULL)
    , legoSlot(NULL)
    , legoScene(NULL)
    , legoViewer(NULL)
    , fProjManager1(NULL)
    , fProjManager2(NULL)
{
}
// -----------------------------------------------------------------------

// -----   Standard constructor   -------------------------------------------
R3BCalifaEventDisplay::R3BCalifaEventDisplay(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fEventManager(NULL)
    , fGeometryVersion(10)
    , fDataHist(NULL)
    , fCalo3d(NULL)
    , fCalo2d(NULL)
    , fLego(NULL)
    , slotLeftTop(NULL)
    , slotLeftBottom(NULL)
    , slotRightTop(NULL)
    , slotRightBottom(NULL)
    , viewerLeftTop(NULL)
    , sceneLeftTop(NULL)
    , viewerRightTop(NULL)
    , sceneRightTop(NULL)
    , viewerLeftBottom(NULL)
    , sceneLeftBottom(NULL)
    , viewerRightBottom(NULL)
    , sceneRightBottom(NULL)
    , legoSlot(NULL)
    , legoScene(NULL)
    , legoViewer(NULL)
    , fProjManager1(NULL)
    , fProjManager2(NULL)
{
}
// -----------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
R3BCalifaEventDisplay::~R3BCalifaEventDisplay() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------
InitStatus R3BCalifaEventDisplay::Init()
{

    FairRootManager* ioManager = FairRootManager::Instance();
    if (!ioManager)
        LOG(fatal) << "Init: No FairRootManager";
    if (!fCrystalHitCA)
        fCrystalHitCA = dynamic_cast<TClonesArray*>(ioManager->GetObject("CalifaCrystalCalData"));

    fEventManager = FairEventManager::Instance();

    CreateHistograms();
    MakeSlots();
    MakeViewerScene();

    return kSUCCESS;
}

// -----   Public method ReInit   --------------------------------------------
InitStatus R3BCalifaEventDisplay::ReInit() { return kSUCCESS; }

// -----  Public method SelectGeometryVersion  ----------------------------------
void R3BCalifaEventDisplay::SelectGeometryVersion(Int_t version) { fGeometryVersion = version; }

// -----   Public method Exec   --------------------------------------------
void R3BCalifaEventDisplay::Exec(Option_t* opt)
{

    if (IsActive())
    {

        Reset();

        R3BCalifaCrystalCalData* crystalHit = nullptr;

        Int_t crystalHits; // Nb of CrystalHits in current event
        crystalHits = fCrystalHitCA->GetEntriesFast();

        Int_t binx = -1, biny = -1;
        Double_t theta = 0., phi = 0., rho = 0., eta = 0.;

        // Loop in Crystal Hits
        for (Int_t i = 0; i < crystalHits; i++)
        {
            crystalHit = dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalHitCA->At(i));
            GetAngles(fGeometryVersion, crystalHit->GetCrystalId(), &theta, &phi, &rho);
        }

        eta = -TMath::Log(TMath::Tan(theta * 0.5f));

        // Matching hit with histogram bin
        for (Int_t b = 1; b <= hcalo->GetNbinsX(); b++)
        {
            if (eta < hcalo->GetXaxis()->GetBinLowEdge(b))
            {
                binx = b - 1;
                break;
            }
        }
        for (Int_t b = 1; b <= hcalo->GetNbinsY(); b++)
        {
            if (phi < hcalo->GetYaxis()->GetBinLowEdge(b))
            {
                biny = b - 1;
                break;
            }
        }

        assert(binx != -1 && biny != -1);

        // crystalHit seems to point to whatever came last when looping over the TCA.
        // Not sure why we want that, but at least make sure it exists.
        assert(crystalHit);

        // Filling histograms

        if (hcalo->GetBinContent(binx, biny) == 0)
        {
            hcalo->SetBinContent(binx, biny, crystalHit->GetEnergy() * 1000);
        }
        else
        {
            hcalo->SetBinContent(binx, biny, hcalo->GetBinContent(binx, biny) + crystalHit->GetEnergy() * 1000);
        }

        // ------- Data Histogram --------------------------------
        fDataHist = new TEveCaloDataHist();
        fDataHist->AddHistogram(hcalo);
        fDataHist->SetSliceColor(0, 2);

        // ------- Calo3D -------------------------------------
        fCalo3d = new TEveCalo3D(fDataHist);
        fCalo3d->SetBarrelRadius(47.00);
        fCalo3d->SetEndCapPos(50.10);
        sceneRightTop->AddElement(fCalo3d);
        Float_t maxH = 60;
        fCalo3d->SetMaxTowerH(maxH);
        fCalo3d->SetMainColor(5);

        fDataHist->IncDenyDestroy();
        fCalo3d->IncDenyDestroy();

        // ------- Projections -------------------------------------
        // --- note.- just two ways of doing the same... -----------
        fProjManager1 = new TEveProjectionManager(TEveProjection::kPT_RPhi);
        sceneLeftBottom->AddElement(fProjManager1);
        TEveProjectionAxes* axes = new TEveProjectionAxes(fProjManager1);
        fProjManager1->ImportElements(fCalo3d);
        sceneLeftBottom->AddElement(axes);

        fProjManager2 = new TEveProjectionManager(TEveProjection::kPT_RhoZ);
        fCalo2d = dynamic_cast<TEveCalo2D*>(fProjManager2->ImportElements(fCalo3d));
        sceneRightBottom->AddElement(fCalo2d);
        TEveProjectionAxes* axes2 = new TEveProjectionAxes(fProjManager2);
        sceneRightBottom->AddElement(axes2);
        fCalo2d->SetMaxTowerH(maxH);

        // ------ Lego ------------------------------------------
        gStyle->SetPalette(1, 0);
        fLego = new TEveCaloLego((TEveCaloData*)fDataHist);
        legoScene->AddElement(fLego);
        fLego->InitMainTrans();
        fLego->RefMainTrans().SetScale(TMath::TwoPi(), TMath::TwoPi(), TMath::Pi());
        // set event handler to move from perspective to orthographic view.
        legoViewer->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
        legoViewer->GetGLViewer()->SetEventHandler(new TEveLegoEventHandler(
            (TGWindow*)legoViewer->GetGLViewer()->GetGLWidget(), legoViewer->GetGLViewer(), fLego));
        fLego->SetGridColor(5);

        gEve->FullRedraw3D(kFALSE);
    }
}

// ---- Public method Reset   --------------------------------------------------
void R3BCalifaEventDisplay::Reset()
{

    hcalo->Reset();

    if (fCalo3d)
        delete fCalo3d;
    if (fLego)
        delete fLego;
}

// ---- Public method Finish   --------------------------------------------------
void R3BCalifaEventDisplay::Finish() {}

// -----   Private method CreateHistograms  -------------------------------------
void R3BCalifaEventDisplay::CreateHistograms()
{

    // Variable x-bin size: defined by crystals polar angle edges
    // Virtual backwards end-cap assumed to be equals to forward end-cap
    Double_t xedges[64] = { 8.71,   11.76,  14.74,  17.21,  19.68,  21.83,  23.97,  26.02,  28.06,  30.00,  31.94,
                            33.85,  35.75,  37.53,  39.30,  41.23,  43.16,  45.09,  47.02,  49.08,  51.14,  53.34,
                            55.53,  57.87,  60.21,  62.70,  65.18,  67.81,  70.43,  73.18,  75.93,  78.78,  81.63,
                            84.59,  87.55,  90.97,  94.39,  97.78,  101.17, 105.67, 110.71, 114.34, 118.51, 122.24,
                            125.97, 129.27, 132.56, 135.35, 138.36, 141.37, 143.84, 146.31, 148.46, 150.60, 152.65,
                            154.69, 156.63, 158.57, 160.48, 162.38, 164.16, 165.93, 167.86, 169.79 };

    Double_t xedgesrad[64], xedgeseta[64];

    for (Int_t i = 0; i < 64; i++)
        xedgesrad[i] = xedges[i] * 3.14 / 180.;

    // conversion of Theta into Eta (default CaloLego x-axis angle)
    for (Int_t i = 0; i < 64; i++)
        xedgeseta[i] = TMath::Log(TMath::Tan(xedgesrad[i] * 0.5f));

    hcalo = new TH2F("hcalo", "hcalo", 63, xedgeseta, 64, -3.14, 3.14);
    hcalo->GetXaxis()->SetTitle("Eta [rad]");
    hcalo->GetYaxis()->SetTitle("Phi [rad]");
}

// -----   Private method MakeSlots  -------------------------------------
void R3BCalifaEventDisplay::MakeSlots()
{

    // frames
    TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    TEveWindowPack* packH = slot->MakePack();
    packH->SetElementName("CrystalHit Projections");
    packH->SetHorizontal();
    packH->SetShowTitleBar(kFALSE);

    slot = packH->NewSlot();
    TEveWindowPack* pack0 = slot->MakePack();
    pack0->SetShowTitleBar(kFALSE);
    slotLeftTop = pack0->NewSlot();
    slotLeftBottom = pack0->NewSlot();

    slot = packH->NewSlot();
    TEveWindowPack* pack1 = slot->MakePack();
    pack1->SetShowTitleBar(kFALSE);
    slotRightTop = pack1->NewSlot();
    slotRightBottom = pack1->NewSlot();

    legoSlot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
}

// -----   Private method MakeViewerScene  ---------------------------------
void R3BCalifaEventDisplay::MakeViewerScene()
{

    slotLeftTop->MakeCurrent();
    viewerLeftTop = gEve->SpawnNewViewer("Main", "Main");
    viewerLeftTop->AddScene(gEve->GetGlobalScene());
    viewerLeftTop->AddScene(gEve->GetEventScene());

    slotRightTop->MakeCurrent();
    viewerRightTop = gEve->SpawnNewViewer("Calo3D", "Calo3D");
    sceneRightTop = gEve->SpawnNewScene("Calo 3D Scene");
    viewerRightTop->AddScene(sceneRightTop);

    slotLeftBottom->MakeCurrent();
    viewerLeftBottom = gEve->SpawnNewViewer("RPhi Projection", "RPhi");
    sceneLeftBottom = gEve->SpawnNewScene("RPhi Scene");
    viewerLeftBottom->AddScene(sceneLeftBottom);
    viewerLeftBottom->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);

    slotRightBottom->MakeCurrent();
    viewerRightBottom = gEve->SpawnNewViewer("RhoZ Projection", "RhoZ");
    sceneRightBottom = gEve->SpawnNewScene("RhoZ Scene");
    viewerRightBottom->AddScene(sceneRightBottom);
    viewerRightBottom->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);

    legoSlot->MakeCurrent();
    legoViewer = gEve->SpawnNewViewer("CrystalHit LegoView", "Lego");
    legoScene = gEve->SpawnNewScene("Lego Scene");
    legoViewer->AddScene(legoScene);
    TEveCaloLegoOverlay* overlay = new TEveCaloLegoOverlay();
    legoViewer->GetGLViewer()->AddOverlayElement(overlay);
    overlay->SetCaloLego(fLego);
}

// -----   Private method GetAngles  ---------------------------------
void R3BCalifaEventDisplay::GetAngles(Int_t geoVersion, Int_t iD, Double_t* polar, Double_t* azimuthal, Double_t* rho)
{

    // Scheme for 8.11 BARREL geometry

    Double_t local[3] = { 0, 0, 0 };
    Double_t master[3] = { 0, 0, 0 };
    Int_t crystalType = 0;
    Int_t crystalCopy = 0;
    Int_t alveolusCopy = 0;
    Int_t crystalInAlveolus = 0;
    Int_t alveoliType30[30] = { 1, 1, 2, 2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,
                                8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15 };
    Int_t alveoliType15[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    Char_t nameVolume[200];

    if (iD > 0)
    {

        if (iD < 3000)
        { // BARREL

            if (iD < 33)
                crystalType = 1; // Alv type 1
            else
                crystalType = (Int_t)((iD - 33) / 128) + 2; // Alv type (2 to 16)
            if (iD < 33)
                crystalCopy = iD; // for Alv type 1
            else
                crystalCopy = ((iD - 33) % 128) + 1; // CrystalCopy (1 to 128)
            if (iD < 33)
                alveolusCopy = iD; // Alv type 1
            else
                alveolusCopy = (Int_t)(((iD - 33) % 128) / 4) + 1; // Alveolus copy (1 to 32)
            if (iD < 33)
                crystalInAlveolus = 1; // Alv type 1
            else
                crystalInAlveolus = (iD - 33) % 4 + 1; // Crystal number in alveolus (1 to 4)

            Int_t alveoliType[16] = { 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6 };

            sprintf(
                nameVolume,
                "/cave_1/CalifaWorld_0/Alveolus_%i_%i/AlveolusInner_%i_1/CrystalWithWrapping_%i_%i_%i/Crystal_%i_%i_1",
                crystalType,
                alveolusCopy - 1,
                crystalType,
                alveoliType[crystalType - 1],
                crystalInAlveolus,
                crystalInAlveolus - 1,
                alveoliType[crystalType - 1],
                crystalInAlveolus);

            // the origin for each crystal is the alveoli corner
            if (crystalType == 1)
            {
                local[0] = 27.108 / 8;
                local[1] = -28.0483 / 8;
                local[2] = 0;
            }
            else if (crystalType == 2 || crystalType == 3 || crystalType == 4 || crystalType == 5)
            {
                if (crystalInAlveolus == 1)
                {
                    local[0] = 37.4639 / 8;
                    local[1] = -8.57573 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 2)
                {
                    local[0] = 37.4639 / 8;
                    local[1] = -31.1043 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 3)
                {
                    local[0] = 9.52012 / 8;
                    local[1] = -8.57573 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 4)
                {
                    local[0] = 9.52012 / 8;
                    local[1] = -31.1043 / 8;
                    local[2] = 0;
                }
            }
            else if (crystalType == 6 || crystalType == 7)
            {
                if (crystalInAlveolus == 1)
                {
                    local[0] = 38.3282 / 8;
                    local[1] = -5.49819 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 2)
                {
                    local[0] = 38.3282 / 8;
                    local[1] = -23.0538 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 3)
                {
                    local[0] = 8.66384 / 8;
                    local[1] = -5.49819 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 4)
                {
                    local[0] = 8.66384 / 8;
                    local[1] = -23.0538 / 8;
                    local[2] = 0;
                }
            }
            else if (crystalType == 8 || crystalType == 9 || crystalType == 10)
            {
                if (crystalInAlveolus == 1)
                {
                    local[0] = 38.3683 / 8;
                    local[1] = -4.71618 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 2)
                {
                    local[0] = 38.3683 / 8;
                    local[1] = -19.8438 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 3)
                {
                    local[0] = 8.43569 / 8;
                    local[1] = -4.71618 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 4)
                {
                    local[0] = 8.43569 / 8;
                    local[1] = -19.8438 / 8;
                    local[2] = 0;
                }
            }
            else if (crystalType == 11 || crystalType == 12 || crystalType == 13)
            {
                if (crystalInAlveolus == 1)
                {
                    local[0] = 38.3495 / 8;
                    local[1] = -4.70373 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 2)
                {
                    local[0] = 38.3495 / 8;
                    local[1] = -19.8403 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 3)
                {
                    local[0] = 8.66654 / 8;
                    local[1] = -4.70373 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 4)
                {
                    local[0] = 8.66654 / 8;
                    local[1] = -19.8403 / 8;
                    local[2] = 0;
                }
            }
            else if (crystalType == 14 || crystalType == 15 || crystalType == 16)
            {
                if (crystalInAlveolus == 1)
                {
                    local[0] = 37.9075 / 8;
                    local[1] = -4.66458 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 2)
                {
                    local[0] = 37.9075 / 8;
                    local[1] = -19.8474 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 3)
                {
                    local[0] = 9.07247 / 8;
                    local[1] = -19.8474 / 8;
                    local[2] = 0;
                }
                else if (crystalInAlveolus == 4)
                {
                    local[0] = 9.07247 / 8;
                    local[1] = -4.66458 / 8;
                    local[2] = 0;
                }
            }

            gGeoManager->cd(nameVolume);
            TGeoNode* currentNode = gGeoManager->GetCurrentNode();
            currentNode->LocalToMaster(local, master);

            sprintf(nameVolume,
                    "/cave_1/CalifaWorld_0/Alveolus_%i_%i/AlveolusInner_%i_1/CrystalWithWrapping_%i_%i_%i",
                    crystalType,
                    alveolusCopy - 1,
                    crystalType,
                    alveoliType[crystalType - 1],
                    crystalInAlveolus,
                    crystalInAlveolus - 1);
            gGeoManager->cd(nameVolume);
            currentNode = gGeoManager->GetCurrentNode();
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            currentNode->LocalToMaster(local, master);

            sprintf(nameVolume,
                    "/cave_1/CalifaWorld_0/Alveolus_%i_%i/AlveolusInner_%i_1",
                    crystalType,
                    alveolusCopy - 1,
                    crystalType);
            gGeoManager->cd(nameVolume);
            currentNode = gGeoManager->GetCurrentNode();
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            currentNode->LocalToMaster(local, master);

            sprintf(nameVolume, "/cave_1/CalifaWorld_0/Alveolus_%i_%i", crystalType, alveolusCopy - 1);
            gGeoManager->cd(nameVolume);
            currentNode = gGeoManager->GetCurrentNode();
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            currentNode->LocalToMaster(local, master);

            sprintf(nameVolume, "/cave_1/CalifaWorld_0");
            gGeoManager->cd(nameVolume);
            currentNode = gGeoManager->GetCurrentNode();
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            currentNode->LocalToMaster(local, master);
        }
        else
        { // ENDCAP

            if (geoVersion == 11)
            {
                crystalType = ((iD - 3000) % 30) + 1;
                crystalCopy = (iD - 3000) / 30 + 1;
                // alveoliType30[30]={1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,
                //  10,10,11,11,12,12,13,13,14,14,15,15};
                sprintf(nameVolume,
                        "/cave_1/CalifaWorld_0/Alveolus_EC_%i_%i/CrystalWithWrapping_%i_1/Crystal_%i_1",
                        alveoliType30[crystalType - 1],
                        crystalCopy - 1,
                        crystalType,
                        crystalType);
            }
            else
            {
                crystalType = ((iD - 3000) % 15) + 1;
                crystalCopy = (iD - 3000) / 15 + 1;
                // alveoliType15[15]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
                sprintf(nameVolume,
                        "/cave_1/CalifaWorld_0/Alveolus_EC_%i_%i/CrystalWithWrapping_%i_1/Crystal_%i_1",
                        alveoliType15[crystalType - 1],
                        crystalCopy - 1,
                        crystalType,
                        crystalType);
            }

            gGeoManager->cd(nameVolume);
            TGeoNode* currentNode = gGeoManager->GetCurrentNode();
            currentNode->LocalToMaster(local, master);

            if (geoVersion == 11)
            {
                sprintf(nameVolume,
                        "/cave_1/CalifaWorld_0/Alveolus_EC_%i_%i/CrystalWithWrapping_%i_1",
                        alveoliType30[crystalType - 1],
                        crystalCopy - 1,
                        crystalType);
            }
            else
            {
                sprintf(nameVolume,
                        "/cave_1/CalifaWorld_0/Alveolus_EC_%i_%i/CrystalWithWrapping_%i_1",
                        alveoliType15[crystalType - 1],
                        crystalCopy - 1,
                        crystalType);
            }

            gGeoManager->cd(nameVolume);
            currentNode = gGeoManager->GetCurrentNode();
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            currentNode->LocalToMaster(local, master);

            if (geoVersion == 11)
            {
                sprintf(nameVolume,
                        "/cave_1/CalifaWorld_0/Alveolus_EC_%i_%i",
                        alveoliType30[crystalType - 1],
                        crystalCopy - 1);
            }
            else
            {
                sprintf(nameVolume,
                        "/cave_1/CalifaWorld_0/Alveolus_EC_%i_%i",
                        alveoliType15[crystalType - 1],
                        crystalCopy - 1);
            }

            gGeoManager->cd(nameVolume);
            currentNode = gGeoManager->GetCurrentNode();
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            currentNode->LocalToMaster(local, master);

            sprintf(nameVolume, "/cave_1/CalifaWorld_0");
            gGeoManager->cd(nameVolume);
            currentNode = gGeoManager->GetCurrentNode();
            local[0] = master[0];
            local[1] = master[1];
            local[2] = master[2];
            currentNode->LocalToMaster(local, master);
        }
    }

    TVector3 masterV(master[0], master[1], master[2]);
    *polar = masterV.Theta();
    *azimuthal = masterV.Phi();
    *rho = masterV.Mag();
}

ClassImp(R3BCalifaEventDisplay)
