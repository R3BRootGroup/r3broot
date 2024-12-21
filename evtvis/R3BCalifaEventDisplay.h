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

// -------------------------------------------------------------------------
// -----            R3BCalifaEventDisplay header file                  -----
// -----            Created 01/10/2012  by P.Cabanelas                 -----
// -------------------------------------------------------------------------

#ifndef R3BCALIFAEVENTDISPLAY_H
#define R3BCALIFAEVENTDISPLAY_H

#include "FairTask.h"

#include "TEveCalo.h"
#include "TEveCaloData.h"
#include "TEveProjectionManager.h"
#include "TEveProjections.h"
#include "TEveScene.h"
#include "TEveViewer.h"
#include "TEveWindow.h"

class TClonesArray;
class FairEventManager;
class FairRootManager;

class TH2F;

class R3BCalifaEventDisplay : public FairTask
{
  public:
    /** Default constructor **/
    R3BCalifaEventDisplay();

    /** Standard Constructor **/
    R3BCalifaEventDisplay(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    ~R3BCalifaEventDisplay();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Finish **/
    virtual void Finish();

    /** Public method SelectGeometryVersion
     **
     ** Defines the geometry
     *@param version  Integer parameter used to select the geometry:
     ** (see documentation /r3broot/cal/perlScripts/README))
     **/
    void SelectGeometryVersion(Int_t version);

    /** Accessors **/

  protected:
    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    TClonesArray* fCrystalHitCA;

    FairEventManager* fEventManager;
    FairRootManager* fManager;

    // Selecting the geometry of the CALIFA calorimeter
    Int_t fGeometryVersion;

    TEveCaloDataHist* fDataHist;
    TEveCalo3D* fCalo3d;
    TEveCalo2D* fCalo2d;
    TEveCaloLego* fLego;

    TEveWindowSlot* slotLeftTop;
    TEveWindowSlot* slotLeftBottom;
    TEveWindowSlot* slotRightTop;
    TEveWindowSlot* slotRightBottom;
    TEveViewer* viewerLeftTop;
    TEveScene* sceneLeftTop;
    TEveViewer* viewerRightTop;
    TEveScene* sceneRightTop;
    TEveViewer* viewerLeftBottom;
    TEveScene* sceneLeftBottom;
    TEveViewer* viewerRightBottom;
    TEveScene* sceneRightBottom;
    TEveWindowSlot* legoSlot;
    TEveScene* legoScene;
    TEveViewer* legoViewer;

    TEveProjectionManager* fProjManager1;
    TEveProjectionManager* fProjManager2;

  private:
    /** Data members **/

    TH2F* hcalo;

    /** Private Methods **/

    void CreateHistograms();
    void MakeSlots();
    void MakeViewerScene();

    void GetAngles(Int_t geoVersion, Int_t iD, Double_t* polar, Double_t* azimuthal, Double_t* rho);

    ClassDef(R3BCalifaEventDisplay, 2);
};

#endif
