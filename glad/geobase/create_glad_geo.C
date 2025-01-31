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

// Last update: Solved incorrect geometry subtractions 31/01/23 (Jose Luis)

#include "TGeoManager.h"
#include "TMath.h"
#include <iomanip>
#include <iostream>

void ConstructSubPartEcrans(TGeoVolume* pWorld);
void ConstructDemiEcransTh(TGeoVolume* pWorld);
void ConstructEnceinteI(TGeoVolume* pWorld);
void ConstructEnceinteE(TGeoVolume* pWorld);
void ConstructFonfE(TGeoVolume* pWorld);
void ConstructFondS(TGeoVolume* pWorld);
void ConstructGToles(TGeoVolume* pWorld);
void ConstructFlange(TGeoVolume* pWorld);
void ConstructFlangeWindows(TGeoVolume* pWorld);

void create_glad_geo(bool create_flange_windows = false, const char* geoTag = "v2023.1", bool create_flange = true)
{
    TGeoManager* gGeoMan = NULL;

    // Global positioning definition
    Double_t glad_angle = 14.08;       // deg
    Double_t distanceToTarget = 308.8; // cm
    Double_t zeroLineOffset = -42.0;   // cm

    // -------   Load media from media file   -----------------------------------
    FairGeoLoader* geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();
    TString geoPath = gSystem->Getenv("VMCWORKDIR");
    TString medFile = geoPath + "/geometry/media_r3b.geo";
    geoFace->setMediaFile(medFile);
    geoFace->readMedia();
    gGeoMan = gGeoManager;
    // --------------------------------------------------------------------------

    // -------   Geometry file name (output)   ----------------------------------
    TString geoFileName = geoPath + "/geometry/glad_";
    if (create_flange_windows)
    {
        geoFileName = geoFileName + "s455_" + geoTag + ".geo.root";
    }
    else
    {
        geoFileName = geoFileName + geoTag + ".geo.root";
    }
    // --------------------------------------------------------------------------

    // -----------------   Get and create the required media    -----------------
    FairGeoMedia* geoMedia = geoFace->getMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

    FairGeoMedium* mCopper = geoMedia->getMedium("copper");
    if (!mCopper)
        Fatal("Main", "FairMedium copper not found");
    geoBuild->createMedium(mCopper);
    TGeoMedium* pMedCu = gGeoMan->GetMedium("copper");
    if (!pMedCu)
        Fatal("Main", "Medium copper not found");

    FairGeoMedium* mIron = geoMedia->getMedium("iron");
    if (!mIron)
        Fatal("Main", "FairMedium iron not found");
    geoBuild->createMedium(mIron);
    TGeoMedium* pMedFe = gGeoMan->GetMedium("iron");
    if (!pMedFe)
        Fatal("Main", "Medium iron not found");

    FairGeoMedium* mAlu = geoMedia->getMedium("aluminium");
    if (!mAlu)
        Fatal("Main", "FairMedium aluminium not found");
    geoBuild->createMedium(mAlu);
    TGeoMedium* pMedAl = gGeoMan->GetMedium("aluminium");
    if (!pMedAl)
        Fatal("Main", "Medium aluminium not found");

    FairGeoMedium* mAir = geoMedia->getMedium("Air");
    if (!mAir)
        Fatal("Main", "FairMedium Air not found");
    geoBuild->createMedium(mAir);
    TGeoMedium* pMed2 = gGeoMan->GetMedium("Air");
    if (!pMed2)
        Fatal("Main", "Medium Air not found");

    FairGeoMedium* mMylar = geoMedia->getMedium("mylar");
    if (!mMylar)
        Fatal("Main", "FairMedium mylar not found");
    geoBuild->createMedium(mMylar);
    TGeoMedium* pMed4 = gGeoMan->GetMedium("mylar");
    if (!pMed4)
        Fatal("Main", "Medium mylar not found");
    // --------------------------------------------------------------------------

    // --------------   Create geometry and top volume  -------------------------
    gGeoMan = (TGeoManager*)gROOT->FindObject("FAIRGeom");
    gGeoMan->SetName("GLADgeom");
    TGeoVolume* top = new TGeoVolumeAssembly("TOP");
    gGeoMan->SetTopVolume(top);
    // --------------------------------------------------------------------------

    // Main Glad Enveloppe
    TGeoVolumeAssembly* pWorld = new TGeoVolumeAssembly("Glad_box");
    TGeoRotation* rot_glob = new TGeoRotation();
    rot_glob->RotateY(90.0);
    rot_glob->RotateZ(-90.0);
    rot_glob->RotateY(glad_angle);
    rot_glob->RotateZ(180.0);
    TGeoCombiTrans* t0 = new TGeoCombiTrans(zeroLineOffset, 0., distanceToTarget, rot_glob);
    top->AddNode(pWorld, 0, t0);

    ConstructSubPartEcrans(pWorld);
    ConstructDemiEcransTh(pWorld);
    ConstructEnceinteI(pWorld);
    ConstructEnceinteE(pWorld);
    ConstructFonfE(pWorld);
    ConstructFondS(pWorld);
    ConstructGToles(pWorld);

    if (create_flange)
    {
        ConstructFlange(pWorld);
    }

    if (create_flange_windows)
    {
        ConstructFlangeWindows(pWorld);
    }

    // ---------------   Finish   -----------------------------------------------
    gGeoMan->CloseGeometry();
    gGeoMan->CheckOverlaps(0.001);
    gGeoMan->PrintOverlaps();
    gGeoMan->Test();

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    // top->Draw();
    top->Write();
    geoFile->Close();

    std::cout << "\033[34m Creating geometry:\033[0m "
              << "\033[33m" << geoFileName << " \033[0m" << std::endl;
    std::cout << "Macro finished successfully." << std::endl;
}

void ConstructFlange(TGeoVolume* pWorld)
{
    Float_t flange_z_offset = -198.;
    Float_t disc_radius = 56.;
    Float_t disc_half_thickn = 1.5;
    Float_t hole_radius = 13.5338;
    Float_t hole_radius_2 = 12.7769;
    Float_t hole_x_offset = -10.408;
    Float_t tube_thickn = 0.3193;
    Float_t tube_half_length = 5.0114 + disc_half_thickn;
    Float_t segm_xy_thickn = 2.4604;
    Float_t segm_z_thickn = 2.3696;

    TGeoVolumeAssembly* flange_world = new TGeoVolumeAssembly("Flange");

    // Convertion to GLAD coordinate system
    TGeoRotation* rot1 = new TGeoRotation();
    rot1->RotateZ(-180.0);
    rot1->RotateZ(90.0);
    rot1->RotateY(-90.0);

    // Zero rotation
    TGeoRotation* zero_rot = new TGeoRotation();

    // Disc and hole
    TGeoShape* f_disc = new TGeoTube(0., disc_radius, disc_half_thickn);
    TGeoShape* f_hole = new TGeoTube(0., hole_radius, 10. * disc_half_thickn);

    // Subtract hole from disc and create new shape
    TGeoRotation* rot_tube = new TGeoRotation();
    rot_tube->RotateY(14.);
    TGeoCombiTrans* m_hole = new TGeoCombiTrans(hole_x_offset, 0., 0., rot_tube);
    TGeoSubtraction* f_subtr_tube = new TGeoSubtraction(f_disc, f_hole, 0, m_hole);
    TGeoShape* f_tube = new TGeoCompositeShape("FlangeDisc", f_subtr_tube);

    TGeoShape* f_tube_2 = new TGeoTube(hole_radius, hole_radius + tube_thickn, tube_half_length);
    TGeoShape* f_segm = new TGeoTube(hole_radius_2, hole_radius_2 + segm_xy_thickn, segm_z_thickn / 2.);
    TGeoCombiTrans* m_segm = new TGeoCombiTrans(0., 0., -tube_half_length - segm_z_thickn / 2., zero_rot);
    TGeoUnion* u_secondary = new TGeoUnion(f_tube_2, f_segm, 0, m_segm);
    TGeoShape* f_secondary = new TGeoCompositeShape("SecondPart", u_secondary);

    TGeoCombiTrans* m_tube_2 = new TGeoCombiTrans(hole_x_offset, 0., 0., rot_tube);
    TGeoUnion* flange_u = new TGeoUnion(f_tube, f_secondary, 0, m_tube_2);
    TGeoShape* flange_u_shape = new TGeoCompositeShape("FlangeDiscTube", flange_u);

    TGeoShape* f_cut = new TGeoTube(0., disc_radius + 0.1, tube_half_length);
    TGeoCombiTrans* m_cut = new TGeoCombiTrans(0., 0., disc_half_thickn + tube_half_length, zero_rot);
    TGeoSubtraction* flange_u_c = new TGeoSubtraction(flange_u_shape, f_cut, 0, m_cut);
    TGeoShape* flange_u_c_shape = new TGeoCompositeShape("FlangeDiscTubeCut", flange_u_c);

    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");
    TGeoVolume* vol_tube = new TGeoVolume("Flange_disc", flange_u_c_shape, pMedFe);
    vol_tube->SetLineColor(921);

    TGeoCombiTrans* m_disc = new TGeoCombiTrans(0., 0., 0., zero_rot);
    flange_world->AddNode(vol_tube, 0, m_disc);

    TGeoCombiTrans* m1 = new TGeoCombiTrans(-1. * flange_z_offset, -4., 0., rot1);
    pWorld->AddNode(flange_world, 0, m1);
}

void ConstructFlangeWindows(TGeoVolume* pWorld)
{
    Float_t z_offset = -198. - 10.408 + 1.8;
    Float_t disc_radius = 13.5338;
    Float_t disc_half_thickn = 0.005;
    TGeoMedium* pMed4 = gGeoManager->GetMedium("mylar");

    TGeoVolumeAssembly* wd_world = new TGeoVolumeAssembly("EntranceWindow");

    // Convertion to GLAD coordinate system
    TGeoRotation* rot1 = new TGeoRotation();
    rot1->RotateZ(-180.0);
    rot1->RotateZ(90.0);
    rot1->RotateY(-90.0);

    // Zero rotation
    TGeoRotation* rote = new TGeoRotation();
    rote->RotateY(14.);

    TGeoShape* f_disc = new TGeoTube(0., disc_radius, disc_half_thickn);

    // Volume:
    TGeoVolume* wd_log = new TGeoVolume("Entrance_WD", f_disc, pMed4);
    wd_log->SetVisLeaves(kTRUE);
    wd_log->SetLineColor(11);

    TGeoCombiTrans* m_disc = new TGeoCombiTrans(0., 0., 0., rote);
    wd_world->AddNode(wd_log, 0, m_disc);

    TGeoCombiTrans* m1 = new TGeoCombiTrans(-1. * z_offset + disc_half_thickn, 9., 0., rot1);
    pWorld->AddNode(wd_world, 0, m1);

    z_offset = -210.;
    float half_thickn = 0.005;

    TGeoVolumeAssembly* wd_world2 = new TGeoVolumeAssembly("ExitWindow");

    // Convertion to GLAD coordinate system
    TGeoRotation* rot2 = new TGeoRotation();
    rot2->RotateZ(-180.0);
    rot2->RotateZ(90.0);
    rot2->RotateY(-90.0);

    TGeoShape* ExitWDBox = new TGeoBBox("Exit_WD", 229., 59., half_thickn);
    // Volume:
    TGeoVolume* ewd_log = new TGeoVolume("Exit_WD_log", ExitWDBox, pMed4);
    ewd_log->SetVisLeaves(kTRUE);

    // Zero rotation
    TGeoRotation* zero_rot = new TGeoRotation();
    m_disc = new TGeoCombiTrans(0., 0., 0., zero_rot);
    wd_world2->AddNode(ewd_log, 0, m_disc);

    m1 = new TGeoCombiTrans(z_offset, 0., 0., rot2);
    pWorld->AddNode(wd_world2, 0, m1);
}

void ConstructSubPartEcrans(TGeoVolume* pWorld)
{
    // AllEcran part description
    Double_t dx, dy, dz;
    Double_t rmin, rmax, phi1, phi2;
    Double_t thx, thy, thz, phx, phy, phz;

    // -- Mixture: Air
    TGeoMedium* pMed2 = gGeoManager->GetMedium("Air");
    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");
    TGeoMedium* pMedC = gGeoManager->GetMedium("copper");

    dx = 114.000000;
    dy = 40.000000;
    dz = 0.150000;
    TGeoShape* pGEcrans_42_1 = new TGeoBBox("GEcrans_42", dx, dy, dz);
    // Shape: GEcrans_50 type: TGeoBBox
    dx = 114.000000;
    dy = 4.210618;
    dz = 0.150000;
    TGeoShape* pGEcrans_50_2 = new TGeoBBox("GEcrans_50", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -44.210600;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix3 = new TGeoCombiTrans("");
    pMatrix3->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode1 = new TGeoUnion(pGEcrans_42_1, pGEcrans_50_2, 0, pMatrix3);
    // Shape: GEcrans_288 type: TGeoCompositeShape
    TGeoShape* pGEcrans_288_89 = new TGeoCompositeShape("GEcrans_288", pBoolNode1);
    // Shape: GEcrans_46 type: TGeoTubeSeg
    rmin = 150.000000;
    rmax = 150.300000;
    dz = 114.000000;
    phi1 = 342.374600;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_46_15 = new TGeoTubeSeg("GEcrans_46", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -48.421200;
    dz = 150.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix6 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix5 = new TGeoCombiTrans("", dx, dy, dz, pMatrix6);
    TGeoUnion* pBoolNode2 = new TGeoUnion(pGEcrans_288_89, pGEcrans_46_15, 0, pMatrix5);
    // Shape: GEcrans_287 type: TGeoCompositeShape
    TGeoShape* pGEcrans_287_90 = new TGeoCompositeShape("GEcrans_287", pBoolNode2);
    // Shape: GEcrans_48 type: TGeoTubeSeg
    rmin = 449.000000;
    rmax = 449.300000;
    dz = 114.000000;
    phi1 = 333.575610;
    phi2 = 342.397460;
    TGeoShape* pGEcrans_48_16 = new TGeoTubeSeg("GEcrans_48", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 42.000000;
    dz = 435.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix9 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix8 = new TGeoCombiTrans("", dx, dy, dz, pMatrix9);
    TGeoUnion* pBoolNode3 = new TGeoUnion(pGEcrans_287_90, pGEcrans_48_16, 0, pMatrix8);
    // Shape: GEcrans_286 type: TGeoCompositeShape
    TGeoShape* pGEcrans_286_91 = new TGeoCompositeShape("GEcrans_286", pBoolNode3);
    // Shape: GEcrans_52 type: TGeoTubeSeg
    rmin = 114.000000;
    rmax = 114.300000;
    dz = 114.000000;
    phi1 = 308.624830;
    phi2 = 333.575610;
    TGeoShape* pGEcrans_52_17 = new TGeoTubeSeg("GEcrans_52", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix12 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix11 = new TGeoCombiTrans("", dx, dy, dz, pMatrix12);
    TGeoUnion* pBoolNode4 = new TGeoUnion(pGEcrans_286_91, pGEcrans_52_17, 0, pMatrix11);
    // Shape: GEcrans_285 type: TGeoCompositeShape
    TGeoShape* pGEcrans_285_92 = new TGeoCompositeShape("GEcrans_285", pBoolNode4);
    // Shape: GEcrans_44 type: TGeoTubeSeg
    rmin = 50.000000;
    rmax = 50.300000;
    dz = 114.000000;
    phi1 = 270.000000;
    phi2 = 308.624830;
    TGeoShape* pGEcrans_44_18 = new TGeoTubeSeg("GEcrans_44", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -157.080500;
    dz = 95.200000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix15 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix14 = new TGeoCombiTrans("", dx, dy, dz, pMatrix15);
    TGeoUnion* pBoolNode5 = new TGeoUnion(pGEcrans_285_92, pGEcrans_44_18, 0, pMatrix14);
    // Shape: GEcrans_290 type: TGeoCompositeShape
    TGeoShape* pGEcrans_290_93 = new TGeoCompositeShape("GEcrans_290", pBoolNode5);
    // Shape: GEcrans_53 type: TGeoBBox
    dx = 114.000000;
    dy = 0.150000;
    dz = 19.474985;
    TGeoShape* pGEcrans_53_3 = new TGeoBBox("GEcrans_53", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -207.230500;
    dz = 114.675000;
    TGeoCombiTrans* pMatrix17 = new TGeoCombiTrans("");
    pMatrix17->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode6 = new TGeoUnion(pGEcrans_290_93, pGEcrans_53_3, 0, pMatrix17);
    // Shape: GEcrans_289 type: TGeoCompositeShape
    TGeoShape* pGEcrans_289_94 = new TGeoCompositeShape("GEcrans_289", pBoolNode6);
    // Shape: GEcrans_303 type: TGeoBBox
    dx = 114.000000;
    dy = 40.000000;
    dz = 0.150000;
    TGeoShape* pGEcrans_303_4 = new TGeoBBox("GEcrans_303", dx, dy, dz);
    // Shape: GEcrans_302 type: TGeoBBox
    dx = 114.000000;
    dy = 4.210618;
    dz = 0.150000;
    TGeoShape* pGEcrans_302_5 = new TGeoBBox("GEcrans_302", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -44.210600;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix19 = new TGeoCombiTrans("");
    pMatrix19->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode7 = new TGeoUnion(pGEcrans_303_4, pGEcrans_302_5, 0, pMatrix19);
    // Shape: GEcrans_301 type: TGeoCompositeShape
    TGeoShape* pGEcrans_301_95 = new TGeoCompositeShape("GEcrans_301", pBoolNode7);
    // Shape: GEcrans_300 type: TGeoTubeSeg
    rmin = 150.000000;
    rmax = 150.300000;
    dz = 114.000000;
    phi1 = 342.374600;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_300_19 = new TGeoTubeSeg("GEcrans_300", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -48.421200;
    dz = 150.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix22 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix21 = new TGeoCombiTrans("", dx, dy, dz, pMatrix22);
    TGeoUnion* pBoolNode8 = new TGeoUnion(pGEcrans_301_95, pGEcrans_300_19, 0, pMatrix21);
    // Shape: GEcrans_299 type: TGeoCompositeShape
    TGeoShape* pGEcrans_299_96 = new TGeoCompositeShape("GEcrans_299", pBoolNode8);
    // Shape: GEcrans_298 type: TGeoTubeSeg
    rmin = 449.000000;
    rmax = 449.300000;
    dz = 114.000000;
    phi1 = 333.575610;
    phi2 = 342.397460;
    TGeoShape* pGEcrans_298_20 = new TGeoTubeSeg("GEcrans_298", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 42.000000;
    dz = 435.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix25 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix24 = new TGeoCombiTrans("", dx, dy, dz, pMatrix25);
    TGeoUnion* pBoolNode9 = new TGeoUnion(pGEcrans_299_96, pGEcrans_298_20, 0, pMatrix24);
    // Shape: GEcrans_297 type: TGeoCompositeShape
    TGeoShape* pGEcrans_297_97 = new TGeoCompositeShape("GEcrans_297", pBoolNode9);
    // Shape: GEcrans_296 type: TGeoTubeSeg
    rmin = 114.000000;
    rmax = 114.300000;
    dz = 114.000000;
    phi1 = 308.624830;
    phi2 = 333.575610;
    TGeoShape* pGEcrans_296_21 = new TGeoTubeSeg("GEcrans_296", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix28 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix27 = new TGeoCombiTrans("", dx, dy, dz, pMatrix28);
    TGeoUnion* pBoolNode10 = new TGeoUnion(pGEcrans_297_97, pGEcrans_296_21, 0, pMatrix27);
    // Shape: GEcrans_295 type: TGeoCompositeShape
    TGeoShape* pGEcrans_295_98 = new TGeoCompositeShape("GEcrans_295", pBoolNode10);
    // Shape: GEcrans_294 type: TGeoTubeSeg
    rmin = 50.000000;
    rmax = 50.300000;
    dz = 114.000000;
    phi1 = 270.000000;
    phi2 = 308.624830;
    TGeoShape* pGEcrans_294_22 = new TGeoTubeSeg("GEcrans_294", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -157.080500;
    dz = 95.200000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix31 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix30 = new TGeoCombiTrans("", dx, dy, dz, pMatrix31);
    TGeoUnion* pBoolNode11 = new TGeoUnion(pGEcrans_295_98, pGEcrans_294_22, 0, pMatrix30);
    // Shape: GEcrans_293 type: TGeoCompositeShape
    TGeoShape* pGEcrans_293_99 = new TGeoCompositeShape("GEcrans_293", pBoolNode11);
    // Shape: GEcrans_292 type: TGeoBBox
    dx = 114.000000;
    dy = 0.150000;
    dz = 19.474985;
    TGeoShape* pGEcrans_292_6 = new TGeoBBox("GEcrans_292", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -207.230500;
    dz = 114.675000;
    TGeoCombiTrans* pMatrix33 = new TGeoCombiTrans("");
    pMatrix33->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode12 = new TGeoUnion(pGEcrans_293_99, pGEcrans_292_6, 0, pMatrix33);
    // Shape: GEcrans_291 type: TGeoCompositeShape
    TGeoShape* pGEcrans_291_100 = new TGeoCompositeShape("GEcrans_291", pBoolNode12);
    // Combi transformation:
    dx = 0.000000;
    dy = 84.000000;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix150 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix149 = new TGeoCombiTrans("", dx, dy, dz, pMatrix150);
    TGeoUnion* pBoolNode13 = new TGeoUnion(pGEcrans_289_94, pGEcrans_291_100, 0, pMatrix149);
    // Shape: GEcrans_304 type: TGeoCompositeShape
    TGeoShape* pGEcrans_304_144 = new TGeoCompositeShape("GEcrans_304", pBoolNode13);
    // Shape: GEcrans_331 type: TGeoBBox
    dx = 114.000000;
    dy = 40.000000;
    dz = 0.150000;
    TGeoShape* pGEcrans_331_7 = new TGeoBBox("GEcrans_331", dx, dy, dz);
    // Shape: GEcrans_330 type: TGeoBBox
    dx = 114.000000;
    dy = 4.210618;
    dz = 0.150000;
    TGeoShape* pGEcrans_330_8 = new TGeoBBox("GEcrans_330", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -44.210600;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix35 = new TGeoCombiTrans("");
    pMatrix35->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode14 = new TGeoUnion(pGEcrans_331_7, pGEcrans_330_8, 0, pMatrix35);
    // Shape: GEcrans_329 type: TGeoCompositeShape
    TGeoShape* pGEcrans_329_101 = new TGeoCompositeShape("GEcrans_329", pBoolNode14);
    // Shape: GEcrans_328 type: TGeoTubeSeg
    rmin = 150.000000;
    rmax = 150.300000;
    dz = 114.000000;
    phi1 = 342.374600;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_328_23 = new TGeoTubeSeg("GEcrans_328", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -48.421200;
    dz = 150.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix38 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix37 = new TGeoCombiTrans("", dx, dy, dz, pMatrix38);
    TGeoUnion* pBoolNode15 = new TGeoUnion(pGEcrans_329_101, pGEcrans_328_23, 0, pMatrix37);
    // Shape: GEcrans_327 type: TGeoCompositeShape
    TGeoShape* pGEcrans_327_102 = new TGeoCompositeShape("GEcrans_327", pBoolNode15);
    // Shape: GEcrans_326 type: TGeoTubeSeg
    rmin = 449.000000;
    rmax = 449.300000;
    dz = 114.000000;
    phi1 = 333.575610;
    phi2 = 342.397460;
    TGeoShape* pGEcrans_326_24 = new TGeoTubeSeg("GEcrans_326", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 42.000000;
    dz = 435.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix41 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix40 = new TGeoCombiTrans("", dx, dy, dz, pMatrix41);
    TGeoUnion* pBoolNode16 = new TGeoUnion(pGEcrans_327_102, pGEcrans_326_24, 0, pMatrix40);
    // Shape: GEcrans_325 type: TGeoCompositeShape
    TGeoShape* pGEcrans_325_103 = new TGeoCompositeShape("GEcrans_325", pBoolNode16);
    // Shape: GEcrans_324 type: TGeoTubeSeg
    rmin = 114.000000;
    rmax = 114.300000;
    dz = 114.000000;
    phi1 = 308.624830;
    phi2 = 333.575610;
    TGeoShape* pGEcrans_324_25 = new TGeoTubeSeg("GEcrans_324", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix44 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix43 = new TGeoCombiTrans("", dx, dy, dz, pMatrix44);
    TGeoUnion* pBoolNode17 = new TGeoUnion(pGEcrans_325_103, pGEcrans_324_25, 0, pMatrix43);
    // Shape: GEcrans_323 type: TGeoCompositeShape
    TGeoShape* pGEcrans_323_104 = new TGeoCompositeShape("GEcrans_323", pBoolNode17);
    // Shape: GEcrans_322 type: TGeoTubeSeg
    rmin = 50.000000;
    rmax = 50.300000;
    dz = 114.000000;
    phi1 = 270.000000;
    phi2 = 308.624830;
    TGeoShape* pGEcrans_322_26 = new TGeoTubeSeg("GEcrans_322", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -157.080500;
    dz = 95.200000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix47 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix46 = new TGeoCombiTrans("", dx, dy, dz, pMatrix47);
    TGeoUnion* pBoolNode18 = new TGeoUnion(pGEcrans_323_104, pGEcrans_322_26, 0, pMatrix46);
    // Shape: GEcrans_321 type: TGeoCompositeShape
    TGeoShape* pGEcrans_321_105 = new TGeoCompositeShape("GEcrans_321", pBoolNode18);
    // Shape: GEcrans_320 type: TGeoBBox
    dx = 114.000000;
    dy = 0.150000;
    dz = 19.474985;
    TGeoShape* pGEcrans_320_9 = new TGeoBBox("GEcrans_320", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -207.230500;
    dz = 114.675000;
    TGeoCombiTrans* pMatrix49 = new TGeoCombiTrans("");
    pMatrix49->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode19 = new TGeoUnion(pGEcrans_321_105, pGEcrans_320_9, 0, pMatrix49);
    // Shape: GEcrans_319 type: TGeoCompositeShape
    TGeoShape* pGEcrans_319_106 = new TGeoCompositeShape("GEcrans_319", pBoolNode19);
    // Shape: GEcrans_318 type: TGeoBBox
    dx = 114.000000;
    dy = 40.000000;
    dz = 0.150000;
    TGeoShape* pGEcrans_318_10 = new TGeoBBox("GEcrans_318", dx, dy, dz);
    // Shape: GEcrans_317 type: TGeoBBox
    dx = 114.000000;
    dy = 4.210618;
    dz = 0.150000;
    TGeoShape* pGEcrans_317_11 = new TGeoBBox("GEcrans_317", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -44.210600;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix51 = new TGeoCombiTrans("");
    pMatrix51->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode20 = new TGeoUnion(pGEcrans_318_10, pGEcrans_317_11, 0, pMatrix51);
    // Shape: GEcrans_316 type: TGeoCompositeShape
    TGeoShape* pGEcrans_316_107 = new TGeoCompositeShape("GEcrans_316", pBoolNode20);
    // Shape: GEcrans_315 type: TGeoTubeSeg
    rmin = 150.000000;
    rmax = 150.300000;
    dz = 114.000000;
    phi1 = 342.374600;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_315_27 = new TGeoTubeSeg("GEcrans_315", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -48.421200;
    dz = 150.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix54 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix53 = new TGeoCombiTrans("", dx, dy, dz, pMatrix54);
    TGeoUnion* pBoolNode21 = new TGeoUnion(pGEcrans_316_107, pGEcrans_315_27, 0, pMatrix53);
    // Shape: GEcrans_314 type: TGeoCompositeShape
    TGeoShape* pGEcrans_314_108 = new TGeoCompositeShape("GEcrans_314", pBoolNode21);
    // Shape: GEcrans_313 type: TGeoTubeSeg
    rmin = 449.000000;
    rmax = 449.300000;
    dz = 114.000000;
    phi1 = 333.575610;
    phi2 = 342.397460;
    TGeoShape* pGEcrans_313_28 = new TGeoTubeSeg("GEcrans_313", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 42.000000;
    dz = 435.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix57 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix56 = new TGeoCombiTrans("", dx, dy, dz, pMatrix57);
    TGeoUnion* pBoolNode22 = new TGeoUnion(pGEcrans_314_108, pGEcrans_313_28, 0, pMatrix56);
    // Shape: GEcrans_312 type: TGeoCompositeShape
    TGeoShape* pGEcrans_312_109 = new TGeoCompositeShape("GEcrans_312", pBoolNode22);
    // Shape: GEcrans_311 type: TGeoTubeSeg
    rmin = 114.000000;
    rmax = 114.300000;
    dz = 114.000000;
    phi1 = 308.624830;
    phi2 = 333.575610;
    TGeoShape* pGEcrans_311_29 = new TGeoTubeSeg("GEcrans_311", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix60 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix59 = new TGeoCombiTrans("", dx, dy, dz, pMatrix60);
    TGeoUnion* pBoolNode23 = new TGeoUnion(pGEcrans_312_109, pGEcrans_311_29, 0, pMatrix59);
    // Shape: GEcrans_310 type: TGeoCompositeShape
    TGeoShape* pGEcrans_310_110 = new TGeoCompositeShape("GEcrans_310", pBoolNode23);
    // Shape: GEcrans_309 type: TGeoTubeSeg
    rmin = 50.000000;
    rmax = 50.300000;
    dz = 114.000000;
    phi1 = 270.000000;
    phi2 = 308.624830;
    TGeoShape* pGEcrans_309_30 = new TGeoTubeSeg("GEcrans_309", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -157.080500;
    dz = 95.200000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix63 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix62 = new TGeoCombiTrans("", dx, dy, dz, pMatrix63);
    TGeoUnion* pBoolNode24 = new TGeoUnion(pGEcrans_310_110, pGEcrans_309_30, 0, pMatrix62);
    // Shape: GEcrans_308 type: TGeoCompositeShape
    TGeoShape* pGEcrans_308_111 = new TGeoCompositeShape("GEcrans_308", pBoolNode24);
    // Shape: GEcrans_307 type: TGeoBBox
    dx = 114.000000;
    dy = 0.150000;
    dz = 19.474985;
    TGeoShape* pGEcrans_307_12 = new TGeoBBox("GEcrans_307", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -207.230500;
    dz = 114.675000;
    TGeoCombiTrans* pMatrix65 = new TGeoCombiTrans("");
    pMatrix65->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode25 = new TGeoUnion(pGEcrans_308_111, pGEcrans_307_12, 0, pMatrix65);
    // Shape: GEcrans_306 type: TGeoCompositeShape
    TGeoShape* pGEcrans_306_112 = new TGeoCompositeShape("GEcrans_306", pBoolNode25);
    // Combi transformation:
    dx = 0.000000;
    dy = 84.000000;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix153 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix152 = new TGeoCombiTrans("", dx, dy, dz, pMatrix153);
    TGeoUnion* pBoolNode26 = new TGeoUnion(pGEcrans_319_106, pGEcrans_306_112, 0, pMatrix152);
    // Shape: GEcrans_305 type: TGeoCompositeShape
    TGeoShape* pGEcrans_305_145 = new TGeoCompositeShape("GEcrans_305", pBoolNode26);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 270.300000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix156 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix155 = new TGeoCombiTrans("", dx, dy, dz, pMatrix156);
    TGeoUnion* pBoolNode27 = new TGeoUnion(pGEcrans_304_144, pGEcrans_305_145, 0, pMatrix155);
    // Shape: GEcrans_332 type: TGeoCompositeShape
    TGeoShape* pGEcrans_332_146 = new TGeoCompositeShape("GEcrans_332", pBoolNode27);
    // Shape: GEcrans_264 type: TGeoBBox
    dx = 114.300000;
    dy = 114.300000;
    dz = 74.300000;
    TGeoShape* pGEcrans_264_85 = new TGeoBBox("GEcrans_264", dx, dy, dz);
    // Shape: GEcrans_263 type: TGeoBBox
    dx = 114.000000;
    dy = 114.000000;
    dz = 74.000000;
    TGeoShape* pGEcrans_263_86 = new TGeoBBox("GEcrans_263", dx, dy, dz);
    TGeoSubtraction* pSubNode = new TGeoSubtraction(pGEcrans_264_85, pGEcrans_263_86, 0, 0);
    // Shape: GEcrans_334 type: TGeoCompositeShape
    TGeoShape* pGEcrans_334_134 = new TGeoCompositeShape("GEcrans_334", pSubNode);
    // Shape: GEcrans_261 type: TGeoTubeSeg
    rmin = 40.000000;
    rmax = 115.000000;
    dz = 37.500000;
    phi1 = 206.424390;
    phi2 = 333.575610;
    TGeoShape* pGEcrans_261_31 = new TGeoTubeSeg("GEcrans_261", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -37.500000;
    TGeoCombiTrans* pMatrix133 = new TGeoCombiTrans("");
    pMatrix133->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode28 = new TGeoIntersection(pGEcrans_334_134, pGEcrans_261_31, 0, pMatrix133);
    // Shape: GEcrans_335 type: TGeoCompositeShape
    TGeoShape* pGEcrans_335_136 = new TGeoCompositeShape("GEcrans_335", pBoolNode28);
    // Combi transformation:
    dx = 115.000000;
    dy = -107.080700;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix159 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix158 = new TGeoCombiTrans("", dx, dy, dz, pMatrix159);
    TGeoUnion* pBoolNode29 = new TGeoUnion(pGEcrans_332_146, pGEcrans_335_136, 0, pMatrix158);
    // Shape: GEcrans_424 type: TGeoCompositeShape
    TGeoShape* pGEcrans_424_147 = new TGeoCompositeShape("GEcrans_424", pBoolNode29);
    // Shape: GEcrans_340 type: TGeoBBox
    dx = 114.300000;
    dy = 114.300000;
    dz = 74.300000;
    TGeoShape* pGEcrans_340_87 = new TGeoBBox("GEcrans_340", dx, dy, dz);
    // Shape: GEcrans_339 type: TGeoBBox
    dx = 114.000000;
    dy = 114.000000;
    dz = 74.000000;
    TGeoShape* pGEcrans_339_88 = new TGeoBBox("GEcrans_339", dx, dy, dz);
    TGeoSubtraction* pBoolNode30 = new TGeoSubtraction(pGEcrans_340_87, pGEcrans_339_88, 0, 0);
    // Shape: GEcrans_338 type: TGeoCompositeShape
    TGeoShape* pGEcrans_338_135 = new TGeoCompositeShape("GEcrans_338", pBoolNode30);
    // Shape: GEcrans_337 type: TGeoTubeSeg
    rmin = 40.000000;
    rmax = 115.000000;
    dz = 37.500000;
    phi1 = 26.424000;
    phi2 = 153.576000;
    TGeoShape* pGEcrans_337_32 = new TGeoTubeSeg("GEcrans_337", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -37.500000;
    TGeoCombiTrans* pMatrix135 = new TGeoCombiTrans("");
    pMatrix135->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode31 = new TGeoIntersection(pGEcrans_338_135, pGEcrans_337_32, 0, pMatrix135);
    // Shape: GEcrans_336 type: TGeoCompositeShape
    TGeoShape* pGEcrans_336_137 = new TGeoCompositeShape("GEcrans_336", pBoolNode31);
    // Combi transformation:
    dx = 115.000000;
    dy = 191.080700;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix162 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix161 = new TGeoCombiTrans("", dx, dy, dz, pMatrix162);
    TGeoUnion* pBoolNode32 = new TGeoUnion(pGEcrans_424_147, pGEcrans_336_137, 0, pMatrix161);
    // Shape: GEcrans_425 type: TGeoCompositeShape
    TGeoShape* pGEcrans_425_148 = new TGeoCompositeShape("GEcrans_425", pBoolNode32);
    // Shape: GEcrans_259 type: TGeoTorus
    Double_t r = 375.000000;
    rmin = 74.000000;
    rmax = 74.300000;
    phi1 = 333.575610;
    Double_t dphi = 52.848780;
    TGeoShape* pGEcrans_259_79 = new TGeoTorus("GEcrans_259", r, rmin, rmax, phi1, dphi);
    // Shape: GEcrans_258 type: TGeoTubeSeg
    rmin = 375.000000;
    rmax = 450.000000;
    dz = 37.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_258_33 = new TGeoTubeSeg("GEcrans_258", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 37.500000;
    TGeoCombiTrans* pMatrix137 = new TGeoCombiTrans("");
    pMatrix137->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode33 = new TGeoIntersection(pGEcrans_259_79, pGEcrans_258_33, 0, pMatrix137);
    // Shape: GEcrans_341 type: TGeoCompositeShape
    TGeoShape* pGEcrans_341_138 = new TGeoCompositeShape("GEcrans_341", pBoolNode33);
    // Combi transformation:
    dx = 115.000000;
    dy = 42.000000;
    dz = 435.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix165 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix164 = new TGeoCombiTrans("", dx, dy, dz, pMatrix165);
    TGeoUnion* pBoolNode34 = new TGeoUnion(pGEcrans_425_148, pGEcrans_341_138, 0, pMatrix164);
    // Shape: GEcrans_426 type: TGeoCompositeShape
    TGeoShape* pGEcrans_426_149 = new TGeoCompositeShape("GEcrans_426", pBoolNode34);
    // Shape: GEcrans_344 type: TGeoTorus
    r = 375.000000;
    rmin = 74.000000;
    rmax = 74.300000;
    phi1 = 153.576000;
    dphi = 52.848000;
    TGeoShape* pGEcrans_344_80 = new TGeoTorus("GEcrans_344", r, rmin, rmax, phi1, dphi);
    // Shape: GEcrans_343 type: TGeoTubeSeg
    rmin = 375.000000;
    rmax = 450.000000;
    dz = 37.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_343_34 = new TGeoTubeSeg("GEcrans_343", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 37.500000;
    TGeoCombiTrans* pMatrix139 = new TGeoCombiTrans("");
    pMatrix139->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode35 = new TGeoIntersection(pGEcrans_344_80, pGEcrans_343_34, 0, pMatrix139);
    // Shape: GEcrans_342 type: TGeoCompositeShape
    TGeoShape* pGEcrans_342_139 = new TGeoCompositeShape("GEcrans_342", pBoolNode35);
    // Combi transformation:
    dx = 115.000000;
    dy = 42.000000;
    dz = -164.850000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix168 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix167 = new TGeoCombiTrans("", dx, dy, dz, pMatrix168);
    TGeoUnion* pBoolNode36 = new TGeoUnion(pGEcrans_426_149, pGEcrans_342_139, 0, pMatrix167);
    // Shape: GEcrans_427 type: TGeoCompositeShape
    TGeoShape* pGEcrans_427_150 = new TGeoCompositeShape("GEcrans_427", pBoolNode36);
    // Shape: GEcrans_276 type: TGeoTorus
    r = 80.000000;
    rmin = 34.000000;
    rmax = 34.300000;
    phi1 = 206.424410;
    dphi = 127.151180;
    TGeoShape* pGEcrans_276_81 = new TGeoTorus("GEcrans_276", r, rmin, rmax, phi1, dphi);
    // Shape: GEcrans_275 type: TGeoTubeSeg
    rmin = 80.000000;
    rmax = 115.000000;
    dz = 35.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_275_35 = new TGeoTubeSeg("GEcrans_275", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -35.000000;
    TGeoCombiTrans* pMatrix141 = new TGeoCombiTrans("");
    pMatrix141->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode37 = new TGeoIntersection(pGEcrans_276_81, pGEcrans_275_35, 0, pMatrix141);
    // Shape: GEcrans_274 type: TGeoCompositeShape
    TGeoShape* pGEcrans_274_140 = new TGeoCompositeShape("GEcrans_274", pBoolNode37);
    // Combi transformation:
    dx = -115.000000;
    dy = -107.080700;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix171 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix170 = new TGeoCombiTrans("", dx, dy, dz, pMatrix171);
    TGeoUnion* pBoolNode38 = new TGeoUnion(pGEcrans_427_150, pGEcrans_274_140, 0, pMatrix170);
    // Shape: GEcrans_428 type: TGeoCompositeShape
    TGeoShape* pGEcrans_428_151 = new TGeoCompositeShape("GEcrans_428", pBoolNode38);
    // Shape: GEcrans_269 type: TGeoTorus
    r = 415.000000;
    rmin = 34.000000;
    rmax = 34.300000;
    phi1 = 333.575610;
    dphi = 52.848780;
    TGeoShape* pGEcrans_269_82 = new TGeoTorus("GEcrans_269", r, rmin, rmax, phi1, dphi);
    // Shape: GEcrans_270 type: TGeoTubeSeg
    rmin = 415.000000;
    rmax = 450.000000;
    dz = 35.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_270_36 = new TGeoTubeSeg("GEcrans_270", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 35.000000;
    TGeoCombiTrans* pMatrix143 = new TGeoCombiTrans("");
    pMatrix143->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode39 = new TGeoIntersection(pGEcrans_269_82, pGEcrans_270_36, 0, pMatrix143);
    // Shape: GEcrans_268 type: TGeoCompositeShape
    TGeoShape* pGEcrans_268_141 = new TGeoCompositeShape("GEcrans_268", pBoolNode39);
    // Combi transformation:
    dx = -115.000000;
    dy = 42.000000;
    dz = -164.850000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix174 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix173 = new TGeoCombiTrans("", dx, dy, dz, pMatrix174);
    TGeoUnion* pBoolNode40 = new TGeoUnion(pGEcrans_428_151, pGEcrans_268_141, 0, pMatrix173);
    // Shape: GEcrans_429 type: TGeoCompositeShape
    TGeoShape* pGEcrans_429_152 = new TGeoCompositeShape("GEcrans_429", pBoolNode40);
    // Shape: GEcrans_347 type: TGeoTorus
    r = 80.000000;
    rmin = 34.000000;
    rmax = 34.300000;
    phi1 = 26.424000;
    dphi = 127.152000;
    TGeoShape* pGEcrans_347_83 = new TGeoTorus("GEcrans_347", r, rmin, rmax, phi1, dphi);
    // Shape: GEcrans_346 type: TGeoTubeSeg
    rmin = 80.000000;
    rmax = 115.000000;
    dz = 35.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_346_37 = new TGeoTubeSeg("GEcrans_346", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -35.000000;
    TGeoCombiTrans* pMatrix145 = new TGeoCombiTrans("");
    pMatrix145->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode41 = new TGeoIntersection(pGEcrans_347_83, pGEcrans_346_37, 0, pMatrix145);
    // Shape: GEcrans_345 type: TGeoCompositeShape
    TGeoShape* pGEcrans_345_142 = new TGeoCompositeShape("GEcrans_345", pBoolNode41);
    // Combi transformation:
    dx = -115.000000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix177 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix176 = new TGeoCombiTrans("", dx, dy, dz, pMatrix177);
    TGeoUnion* pBoolNode42 = new TGeoUnion(pGEcrans_429_152, pGEcrans_345_142, 0, pMatrix176);
    // Shape: GEcrans_430 type: TGeoCompositeShape
    TGeoShape* pGEcrans_430_153 = new TGeoCompositeShape("GEcrans_430", pBoolNode42);
    // Shape: GEcrans_350 type: TGeoTorus
    r = 415.000000;
    rmin = 34.000000;
    rmax = 34.300000;
    phi1 = 153.576000;
    dphi = 52.848000;
    TGeoShape* pGEcrans_350_84 = new TGeoTorus("GEcrans_350", r, rmin, rmax, phi1, dphi);
    // Shape: GEcrans_349 type: TGeoTubeSeg
    rmin = 415.000000;
    rmax = 450.000000;
    dz = 35.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEcrans_349_38 = new TGeoTubeSeg("GEcrans_349", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 35.000000;
    TGeoCombiTrans* pMatrix147 = new TGeoCombiTrans("");
    pMatrix147->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode43 = new TGeoIntersection(pGEcrans_350_84, pGEcrans_349_38, 0, pMatrix147);
    // Shape: GEcrans_348 type: TGeoCompositeShape
    TGeoShape* pGEcrans_348_143 = new TGeoCompositeShape("GEcrans_348", pBoolNode43);
    // Combi transformation:
    dx = -115.000000;
    dy = 42.000000;
    dz = 435.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix180 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix179 = new TGeoCombiTrans("", dx, dy, dz, pMatrix180);
    TGeoUnion* pBoolNode44 = new TGeoUnion(pGEcrans_430_153, pGEcrans_348_143, 0, pMatrix179);
    // Shape: GEcrans_431 type: TGeoCompositeShape
    TGeoShape* pGEcrans_431_154 = new TGeoCompositeShape("GEcrans_431", pBoolNode44);
    // Shape: GEcrans_96 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_96_39 = new TGeoTubeSeg("GEcrans_96", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_95 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_95_40 = new TGeoTubeSeg("GEcrans_95", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix68 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix67 = new TGeoCombiTrans("", dx, dy, dz, pMatrix68);
    TGeoUnion* pBoolNode45 = new TGeoUnion(pGEcrans_96_39, pGEcrans_95_40, 0, pMatrix67);
    // Shape: GEcrans_94 type: TGeoCompositeShape
    TGeoShape* pGEcrans_94_113 = new TGeoCompositeShape("GEcrans_94", pBoolNode45);
    // Combi transformation:
    dx = 114.500000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix183 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix182 = new TGeoCombiTrans("", dx, dy, dz, pMatrix183);
    TGeoUnion* pBoolNode46 = new TGeoUnion(pGEcrans_431_154, pGEcrans_94_113, 0, pMatrix182);
    // Shape: GEcrans_432 type: TGeoCompositeShape
    TGeoShape* pGEcrans_432_155 = new TGeoCompositeShape("GEcrans_432", pBoolNode46);
    // Shape: GEcrans_354 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_354_41 = new TGeoTubeSeg("GEcrans_354", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_353 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_353_42 = new TGeoTubeSeg("GEcrans_353", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix71 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix70 = new TGeoCombiTrans("", dx, dy, dz, pMatrix71);
    TGeoUnion* pBoolNode47 = new TGeoUnion(pGEcrans_354_41, pGEcrans_353_42, 0, pMatrix70);
    // Shape: GEcrans_352 type: TGeoCompositeShape
    TGeoShape* pGEcrans_352_114 = new TGeoCompositeShape("GEcrans_352", pBoolNode47);
    // Combi transformation:
    dx = 53.500000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix186 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix185 = new TGeoCombiTrans("", dx, dy, dz, pMatrix186);
    TGeoUnion* pBoolNode48 = new TGeoUnion(pGEcrans_432_155, pGEcrans_352_114, 0, pMatrix185);
    // Shape: GEcrans_433 type: TGeoCompositeShape
    TGeoShape* pGEcrans_433_156 = new TGeoCompositeShape("GEcrans_433", pBoolNode48);
    // Shape: GEcrans_357 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_357_43 = new TGeoTubeSeg("GEcrans_357", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_356 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_356_44 = new TGeoTubeSeg("GEcrans_356", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix74 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix73 = new TGeoCombiTrans("", dx, dy, dz, pMatrix74);
    TGeoUnion* pBoolNode49 = new TGeoUnion(pGEcrans_357_43, pGEcrans_356_44, 0, pMatrix73);
    // Shape: GEcrans_355 type: TGeoCompositeShape
    TGeoShape* pGEcrans_355_115 = new TGeoCompositeShape("GEcrans_355", pBoolNode49);
    // Combi transformation:
    dx = -8.000000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix189 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix188 = new TGeoCombiTrans("", dx, dy, dz, pMatrix189);
    TGeoUnion* pBoolNode50 = new TGeoUnion(pGEcrans_433_156, pGEcrans_355_115, 0, pMatrix188);
    // Shape: GEcrans_434 type: TGeoCompositeShape
    TGeoShape* pGEcrans_434_157 = new TGeoCompositeShape("GEcrans_434", pBoolNode50);
    // Shape: GEcrans_360 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_360_45 = new TGeoTubeSeg("GEcrans_360", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_359 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_359_46 = new TGeoTubeSeg("GEcrans_359", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix77 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix76 = new TGeoCombiTrans("", dx, dy, dz, pMatrix77);
    TGeoUnion* pBoolNode51 = new TGeoUnion(pGEcrans_360_45, pGEcrans_359_46, 0, pMatrix76);
    // Shape: GEcrans_358 type: TGeoCompositeShape
    TGeoShape* pGEcrans_358_116 = new TGeoCompositeShape("GEcrans_358", pBoolNode51);
    // Combi transformation:
    dx = -65.000000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix192 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix191 = new TGeoCombiTrans("", dx, dy, dz, pMatrix192);
    TGeoUnion* pBoolNode52 = new TGeoUnion(pGEcrans_434_157, pGEcrans_358_116, 0, pMatrix191);
    // Shape: GEcrans_435 type: TGeoCompositeShape
    TGeoShape* pGEcrans_435_158 = new TGeoCompositeShape("GEcrans_435", pBoolNode52);
    // Shape: GEcrans_363 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_363_47 = new TGeoTubeSeg("GEcrans_363", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_362 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_362_48 = new TGeoTubeSeg("GEcrans_362", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix80 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix79 = new TGeoCombiTrans("", dx, dy, dz, pMatrix80);
    TGeoUnion* pBoolNode53 = new TGeoUnion(pGEcrans_363_47, pGEcrans_362_48, 0, pMatrix79);
    // Shape: GEcrans_361 type: TGeoCompositeShape
    TGeoShape* pGEcrans_361_117 = new TGeoCompositeShape("GEcrans_361", pBoolNode53);
    // Combi transformation:
    dx = -114.500000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix195 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix194 = new TGeoCombiTrans("", dx, dy, dz, pMatrix195);
    TGeoUnion* pBoolNode54 = new TGeoUnion(pGEcrans_435_158, pGEcrans_361_117, 0, pMatrix194);
    // Shape: GEcrans_436 type: TGeoCompositeShape
    TGeoShape* pGEcrans_436_159 = new TGeoCompositeShape("GEcrans_436", pBoolNode54);
    // Shape: GEcrans_370 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_370_49 = new TGeoTubeSeg("GEcrans_370", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_369 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_369_50 = new TGeoTubeSeg("GEcrans_369", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix83 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix82 = new TGeoCombiTrans("", dx, dy, dz, pMatrix83);
    TGeoUnion* pBoolNode55 = new TGeoUnion(pGEcrans_370_49, pGEcrans_369_50, 0, pMatrix82);
    // Shape: GEcrans_364 type: TGeoCompositeShape
    TGeoShape* pGEcrans_364_118 = new TGeoCompositeShape("GEcrans_364", pBoolNode55);
    // Combi transformation:
    dx = 114.500000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix198 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix197 = new TGeoCombiTrans("", dx, dy, dz, pMatrix198);
    TGeoUnion* pBoolNode56 = new TGeoUnion(pGEcrans_436_159, pGEcrans_364_118, 0, pMatrix197);
    // Shape: GEcrans_437 type: TGeoCompositeShape
    TGeoShape* pGEcrans_437_160 = new TGeoCompositeShape("GEcrans_437", pBoolNode56);
    // Shape: GEcrans_372 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_372_51 = new TGeoTubeSeg("GEcrans_372", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_371 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_371_52 = new TGeoTubeSeg("GEcrans_371", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix86 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix85 = new TGeoCombiTrans("", dx, dy, dz, pMatrix86);
    TGeoUnion* pBoolNode57 = new TGeoUnion(pGEcrans_372_51, pGEcrans_371_52, 0, pMatrix85);
    // Shape: GEcrans_365 type: TGeoCompositeShape
    TGeoShape* pGEcrans_365_119 = new TGeoCompositeShape("GEcrans_365", pBoolNode57);
    // Combi transformation:
    dx = 53.500000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix201 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix200 = new TGeoCombiTrans("", dx, dy, dz, pMatrix201);
    TGeoUnion* pBoolNode58 = new TGeoUnion(pGEcrans_437_160, pGEcrans_365_119, 0, pMatrix200);
    // Shape: GEcrans_438 type: TGeoCompositeShape
    TGeoShape* pGEcrans_438_161 = new TGeoCompositeShape("GEcrans_438", pBoolNode58);
    // Shape: GEcrans_374 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_374_53 = new TGeoTubeSeg("GEcrans_374", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_373 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_373_54 = new TGeoTubeSeg("GEcrans_373", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix89 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix88 = new TGeoCombiTrans("", dx, dy, dz, pMatrix89);
    TGeoUnion* pBoolNode59 = new TGeoUnion(pGEcrans_374_53, pGEcrans_373_54, 0, pMatrix88);
    // Shape: GEcrans_366 type: TGeoCompositeShape
    TGeoShape* pGEcrans_366_120 = new TGeoCompositeShape("GEcrans_366", pBoolNode59);
    // Combi transformation:
    dx = -8.000000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix204 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix203 = new TGeoCombiTrans("", dx, dy, dz, pMatrix204);
    TGeoUnion* pBoolNode60 = new TGeoUnion(pGEcrans_438_161, pGEcrans_366_120, 0, pMatrix203);
    // Shape: GEcrans_439 type: TGeoCompositeShape
    TGeoShape* pGEcrans_439_162 = new TGeoCompositeShape("GEcrans_439", pBoolNode60);
    // Shape: GEcrans_376 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_376_55 = new TGeoTubeSeg("GEcrans_376", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_375 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_375_56 = new TGeoTubeSeg("GEcrans_375", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix92 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix91 = new TGeoCombiTrans("", dx, dy, dz, pMatrix92);
    TGeoUnion* pBoolNode61 = new TGeoUnion(pGEcrans_376_55, pGEcrans_375_56, 0, pMatrix91);
    // Shape: GEcrans_367 type: TGeoCompositeShape
    TGeoShape* pGEcrans_367_121 = new TGeoCompositeShape("GEcrans_367", pBoolNode61);
    // Combi transformation:
    dx = -65.000000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix207 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix206 = new TGeoCombiTrans("", dx, dy, dz, pMatrix207);
    TGeoUnion* pBoolNode62 = new TGeoUnion(pGEcrans_439_162, pGEcrans_367_121, 0, pMatrix206);
    // Shape: GEcrans_440 type: TGeoCompositeShape
    TGeoShape* pGEcrans_440_163 = new TGeoCompositeShape("GEcrans_440", pBoolNode62);
    // Shape: GEcrans_378 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_378_57 = new TGeoTubeSeg("GEcrans_378", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_377 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_377_58 = new TGeoTubeSeg("GEcrans_377", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix95 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix94 = new TGeoCombiTrans("", dx, dy, dz, pMatrix95);
    TGeoUnion* pBoolNode63 = new TGeoUnion(pGEcrans_378_57, pGEcrans_377_58, 0, pMatrix94);
    // Shape: GEcrans_368 type: TGeoCompositeShape
    TGeoShape* pGEcrans_368_122 = new TGeoCompositeShape("GEcrans_368", pBoolNode63);
    // Combi transformation:
    dx = -114.500000;
    dy = 191.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix210 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix209 = new TGeoCombiTrans("", dx, dy, dz, pMatrix210);
    TGeoUnion* pBoolNode64 = new TGeoUnion(pGEcrans_440_163, pGEcrans_368_122, 0, pMatrix209);
    // Shape: GEcrans_441 type: TGeoCompositeShape
    TGeoShape* pGEcrans_441_164 = new TGeoCompositeShape("GEcrans_441", pBoolNode64);
    // Shape: GEcrans_390 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_390_59 = new TGeoTubeSeg("GEcrans_390", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_389 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_389_60 = new TGeoTubeSeg("GEcrans_389", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix98 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix97 = new TGeoCombiTrans("", dx, dy, dz, pMatrix98);
    TGeoUnion* pBoolNode65 = new TGeoUnion(pGEcrans_390_59, pGEcrans_389_60, 0, pMatrix97);
    // Shape: GEcrans_379 type: TGeoCompositeShape
    TGeoShape* pGEcrans_379_123 = new TGeoCompositeShape("GEcrans_379", pBoolNode65);
    // Combi transformation:
    dx = 114.500000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix213 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix212 = new TGeoCombiTrans("", dx, dy, dz, pMatrix213);
    TGeoUnion* pBoolNode66 = new TGeoUnion(pGEcrans_441_164, pGEcrans_379_123, 0, pMatrix212);
    // Shape: GEcrans_442 type: TGeoCompositeShape
    TGeoShape* pGEcrans_442_165 = new TGeoCompositeShape("GEcrans_442", pBoolNode66);
    // Shape: GEcrans_392 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_392_61 = new TGeoTubeSeg("GEcrans_392", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_391 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_391_62 = new TGeoTubeSeg("GEcrans_391", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix101 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix100 = new TGeoCombiTrans("", dx, dy, dz, pMatrix101);
    TGeoUnion* pBoolNode67 = new TGeoUnion(pGEcrans_392_61, pGEcrans_391_62, 0, pMatrix100);
    // Shape: GEcrans_380 type: TGeoCompositeShape
    TGeoShape* pGEcrans_380_124 = new TGeoCompositeShape("GEcrans_380", pBoolNode67);
    // Combi transformation:
    dx = 53.500000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix216 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix215 = new TGeoCombiTrans("", dx, dy, dz, pMatrix216);
    TGeoUnion* pBoolNode68 = new TGeoUnion(pGEcrans_442_165, pGEcrans_380_124, 0, pMatrix215);
    // Shape: GEcrans_443 type: TGeoCompositeShape
    TGeoShape* pGEcrans_443_166 = new TGeoCompositeShape("GEcrans_443", pBoolNode68);
    // Shape: GEcrans_394 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_394_63 = new TGeoTubeSeg("GEcrans_394", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_393 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_393_64 = new TGeoTubeSeg("GEcrans_393", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix104 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix103 = new TGeoCombiTrans("", dx, dy, dz, pMatrix104);
    TGeoUnion* pBoolNode69 = new TGeoUnion(pGEcrans_394_63, pGEcrans_393_64, 0, pMatrix103);
    // Shape: GEcrans_381 type: TGeoCompositeShape
    TGeoShape* pGEcrans_381_125 = new TGeoCompositeShape("GEcrans_381", pBoolNode69);
    // Combi transformation:
    dx = -8.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix219 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix218 = new TGeoCombiTrans("", dx, dy, dz, pMatrix219);
    TGeoUnion* pBoolNode100 = new TGeoUnion(pGEcrans_443_166, pGEcrans_381_125, 0, pMatrix218);
    // Shape: GEcrans_444 type: TGeoCompositeShape
    TGeoShape* pGEcrans_444_167 = new TGeoCompositeShape("GEcrans_444", pBoolNode100);
    // Shape: GEcrans_396 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_396_65 = new TGeoTubeSeg("GEcrans_396", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_395 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_395_66 = new TGeoTubeSeg("GEcrans_395", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix107 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix106 = new TGeoCombiTrans("", dx, dy, dz, pMatrix107);
    TGeoUnion* pBoolNode70 = new TGeoUnion(pGEcrans_396_65, pGEcrans_395_66, 0, pMatrix106);
    // Shape: GEcrans_382 type: TGeoCompositeShape
    TGeoShape* pGEcrans_382_126 = new TGeoCompositeShape("GEcrans_382", pBoolNode70);
    // Combi transformation:
    dx = -65.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix222 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix221 = new TGeoCombiTrans("", dx, dy, dz, pMatrix222);
    TGeoUnion* pBoolNode71 = new TGeoUnion(pGEcrans_444_167, pGEcrans_382_126, 0, pMatrix221);
    // Shape: GEcrans_445 type: TGeoCompositeShape
    TGeoShape* pGEcrans_445_168 = new TGeoCompositeShape("GEcrans_445", pBoolNode71);
    // Shape: GEcrans_398 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_398_67 = new TGeoTubeSeg("GEcrans_398", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_397 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_397_68 = new TGeoTubeSeg("GEcrans_397", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix110 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix109 = new TGeoCombiTrans("", dx, dy, dz, pMatrix110);
    TGeoUnion* pBoolNode72 = new TGeoUnion(pGEcrans_398_67, pGEcrans_397_68, 0, pMatrix109);
    // Shape: GEcrans_383 type: TGeoCompositeShape
    TGeoShape* pGEcrans_383_127 = new TGeoCompositeShape("GEcrans_383", pBoolNode72);
    // Combi transformation:
    dx = -114.500000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix225 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix224 = new TGeoCombiTrans("", dx, dy, dz, pMatrix225);
    TGeoUnion* pBoolNode73 = new TGeoUnion(pGEcrans_445_168, pGEcrans_383_127, 0, pMatrix224);
    // Shape: GEcrans_446 type: TGeoCompositeShape
    TGeoShape* pGEcrans_446_169 = new TGeoCompositeShape("GEcrans_446", pBoolNode73);
    // Shape: GEcrans_415 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_415_69 = new TGeoTubeSeg("GEcrans_415", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_414 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_414_70 = new TGeoTubeSeg("GEcrans_414", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix113 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix112 = new TGeoCombiTrans("", dx, dy, dz, pMatrix113);
    TGeoUnion* pBoolNode74 = new TGeoUnion(pGEcrans_415_69, pGEcrans_414_70, 0, pMatrix112);
    // Shape: GEcrans_409 type: TGeoCompositeShape
    TGeoShape* pGEcrans_409_128 = new TGeoCompositeShape("GEcrans_409", pBoolNode74);
    // Combi transformation:
    dx = -114.500000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 360.000000;
    TGeoRotation* pMatrix228 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix227 = new TGeoCombiTrans("", dx, dy, dz, pMatrix228);
    TGeoUnion* pBoolNode75 = new TGeoUnion(pGEcrans_446_169, pGEcrans_409_128, 0, pMatrix227);
    // Shape: GEcrans_447 type: TGeoCompositeShape
    TGeoShape* pGEcrans_447_170 = new TGeoCompositeShape("GEcrans_447", pBoolNode75);
    // Shape: GEcrans_417 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_417_71 = new TGeoTubeSeg("GEcrans_417", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_416 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_416_72 = new TGeoTubeSeg("GEcrans_416", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix116 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix115 = new TGeoCombiTrans("", dx, dy, dz, pMatrix116);
    TGeoUnion* pBoolNode76 = new TGeoUnion(pGEcrans_417_71, pGEcrans_416_72, 0, pMatrix115);
    // Shape: GEcrans_410 type: TGeoCompositeShape
    TGeoShape* pGEcrans_410_129 = new TGeoCompositeShape("GEcrans_410", pBoolNode76);
    // Combi transformation:
    dx = -65.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 360.000000;
    TGeoRotation* pMatrix231 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix230 = new TGeoCombiTrans("", dx, dy, dz, pMatrix231);
    TGeoUnion* pBoolNode77 = new TGeoUnion(pGEcrans_447_170, pGEcrans_410_129, 0, pMatrix230);
    // Shape: GEcrans_448 type: TGeoCompositeShape
    TGeoShape* pGEcrans_448_171 = new TGeoCompositeShape("GEcrans_448", pBoolNode77);
    // Shape: GEcrans_419 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_419_73 = new TGeoTubeSeg("GEcrans_419", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_418 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_418_74 = new TGeoTubeSeg("GEcrans_418", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix119 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix118 = new TGeoCombiTrans("", dx, dy, dz, pMatrix119);
    TGeoUnion* pBoolNode78 = new TGeoUnion(pGEcrans_419_73, pGEcrans_418_74, 0, pMatrix118);
    // Shape: GEcrans_411 type: TGeoCompositeShape
    TGeoShape* pGEcrans_411_130 = new TGeoCompositeShape("GEcrans_411", pBoolNode78);
    // Combi transformation:
    dx = -8.000000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 360.000000;
    TGeoRotation* pMatrix234 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix233 = new TGeoCombiTrans("", dx, dy, dz, pMatrix234);
    TGeoUnion* pBoolNode79 = new TGeoUnion(pGEcrans_448_171, pGEcrans_411_130, 0, pMatrix233);
    // Shape: GEcrans_449 type: TGeoCompositeShape
    TGeoShape* pGEcrans_449_172 = new TGeoCompositeShape("GEcrans_449", pBoolNode79);
    // Shape: GEcrans_421 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_421_75 = new TGeoTubeSeg("GEcrans_421", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_420 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_420_76 = new TGeoTubeSeg("GEcrans_420", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix122 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix121 = new TGeoCombiTrans("", dx, dy, dz, pMatrix122);
    TGeoUnion* pBoolNode80 = new TGeoUnion(pGEcrans_421_75, pGEcrans_420_76, 0, pMatrix121);
    // Shape: GEcrans_412 type: TGeoCompositeShape
    TGeoShape* pGEcrans_412_131 = new TGeoCompositeShape("GEcrans_412", pBoolNode80);
    // Combi transformation:
    dx = 53.500000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 360.000000;
    TGeoRotation* pMatrix237 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix236 = new TGeoCombiTrans("", dx, dy, dz, pMatrix237);
    TGeoUnion* pBoolNode81 = new TGeoUnion(pGEcrans_449_172, pGEcrans_412_131, 0, pMatrix236);
    // Shape: GEcrans_450 type: TGeoCompositeShape
    TGeoShape* pGEcrans_450_173 = new TGeoCompositeShape("GEcrans_450", pBoolNode81);
    // Shape: GEcrans_423 type: TGeoTubeSeg
    rmin = 99.000000;
    rmax = 115.000000;
    dz = 0.500000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEcrans_423_77 = new TGeoTubeSeg("GEcrans_423", rmin, rmax, dz, phi1, phi2);
    // Shape: GEcrans_422 type: TGeoTubeSeg
    rmin = 434.000000;
    rmax = 450.000000;
    dz = 0.500000;
    phi1 = 153.575610;
    phi2 = 180.000000;
    TGeoShape* pGEcrans_422_78 = new TGeoTubeSeg("GEcrans_422", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.006700;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix125 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix124 = new TGeoCombiTrans("", dx, dy, dz, pMatrix125);
    TGeoUnion* pBoolNode811 = new TGeoUnion(pGEcrans_423_77, pGEcrans_422_78, 0, pMatrix124);
    // Shape: GEcrans_413 type: TGeoCompositeShape
    TGeoShape* pGEcrans_413_132 = new TGeoCompositeShape("GEcrans_413", pBoolNode811);
    // Combi transformation:
    dx = 114.500000;
    dy = -107.080500;
    dz = 135.150000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 90.000000;
    phz = 360.000000;
    TGeoRotation* pMatrix240 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix239 = new TGeoCombiTrans("", dx, dy, dz, pMatrix240);
    TGeoUnion* pBoolNode82 = new TGeoUnion(pGEcrans_450_173, pGEcrans_413_132, 0, pMatrix239);
    // Shape: GEcrans_451 type: TGeoCompositeShape
    TGeoShape* pGEcrans_451_174 = new TGeoCompositeShape("GEcrans_451", pBoolNode82);
    // Shape: GEcrans_139 type: TGeoBBox
    dx = 220.000000;
    dy = 300.000000;
    dz = 1.000000;
    TGeoShape* pGEcrans_139_13 = new TGeoBBox("GEcrans_139", dx, dy, dz);
    // Shape: GEcrans_138 type: TGeoBBox
    dx = 220.000000;
    dy = 2.000000;
    dz = 200.000000;
    TGeoShape* pGEcrans_138_14 = new TGeoBBox("GEcrans_138", dx, dy, dz);
    TGeoUnion* pBoolNode83 = new TGeoUnion(pGEcrans_139_13, pGEcrans_138_14, 0, 0);
    // Shape: GEcrans_140 type: TGeoCompositeShape
    TGeoShape* pGEcrans_140_133 = new TGeoCompositeShape("GEcrans_140", pBoolNode83);
    // Combi transformation:
    dx = 0.000000;
    dy = 42.000000;
    dz = 135.150000;
    TGeoCombiTrans* pMatrix242 = new TGeoCombiTrans("");
    pMatrix242->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode84 = new TGeoSubtraction(pGEcrans_451_174, pGEcrans_140_133, 0, pMatrix242);
    // Shape: GEcrans_PartBody type: TGeoCompositeShape
    TGeoShape* pGEcrans_PartBody_175 = new TGeoCompositeShape("GEcrans_PartBody", pBoolNode84);
    // Volume: GEcrans
    TGeoVolume* pGEcrans = new TGeoVolume("GEcrans", pGEcrans_PartBody_175, pMedC);
    pGEcrans->SetLineColor(920);
    pGEcrans->SetVisLeaves(kTRUE);

    // HMatrix:
    Double_t tr[3];
    Double_t rot[9];

    tr[0] = 0.000000;
    tr[1] = -42.000000;
    tr[2] = -130.000000;

    rot[0] = 1.000000;
    rot[1] = 0.000000;
    rot[2] = 0.000000;
    rot[3] = 0.000000;
    rot[4] = 1.000000;
    rot[5] = 0.000000;
    rot[6] = 0.000000;
    rot[7] = 0.000000;
    rot[8] = 1.000000;

    TGeoHMatrix* pMatrix243 = new TGeoHMatrix("");
    pMatrix243->SetTranslation(tr);
    pMatrix243->SetRotation(rot);
    pMatrix243->SetBit(TGeoMatrix::kGeoTranslation);
    pMatrix243->SetBit(TGeoMatrix::kGeoRotation);

    // Add the top Ecran volume
    pWorld->AddNode(pGEcrans, 1, pMatrix243);
}

void ConstructDemiEcransTh(TGeoVolume* pWorld)
{
    // Demi Ecrans part
    Double_t dx, dy, dz;
    Double_t dx1, dy1, dz1;
    Double_t dx2, dy2, dz2;
    Double_t rmin, rmax, phi1, phi2;
    Double_t thx, thy, thz, phx, phy, phz;

    // -- Mixture: Air
    TGeoMedium* pMed2 = gGeoManager->GetMedium("Air");
    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");

    Double_t tr[3], rot[9];

    // Shape: G2202001_Demi_Ecran_thermique_interne_5 type: TGeoBBox
    dx = 4.119954;
    dy = 24.750000;
    dz = 18.500000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_5_2 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_5", dx, dy, dz);
    // Shape: G2202001_Demi_Ecran_thermique_interne_4 type: TGeoTrd2
    dx1 = 18.505235;
    dx2 = 22.207822;
    dy1 = 24.750000;
    dy2 = 24.750000;
    dz = 21.190413;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_4_10 =
        new TGeoTrd2("G2202001_Demi_Ecran_thermique_interne_4", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -25.310300;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix4 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix3 = new TGeoCombiTrans("", dx, dy, dz, pMatrix4);
    TGeoUnion* pBoolNode1 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_5_2, pG2202001_Demi_Ecran_thermique_interne_4_10, 0, pMatrix3);
    // Shape: G2202001_Demi_Ecran_thermique_interne_6 type: TGeoCompositeShape
    // D. Kresan: skip entrance box - use Trd2 shape directly
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_6_16 = pG2202001_Demi_Ecran_thermique_interne_4_10;
    //    new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_6", pBoolNode1);
    // Shape: G2202001_Demi_Ecran_thermique_interne_3 type: TGeoTrd2
    dx1 = 22.207822;
    dx2 = 40.111503;
    dy1 = 24.750000;
    dy2 = 93.621885;
    dz = 102.320005;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_3_11 =
        new TGeoTrd2("G2202001_Demi_Ecran_thermique_interne_3", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -148.820700;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix7 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix6 = new TGeoCombiTrans("", dx, dy, dz, pMatrix7);
    // D. Kresan: skip entrance box - use Trd2 shape directly - use pMatrix3 instead of 0
    TGeoUnion* pBoolNode2 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_6_16, pG2202001_Demi_Ecran_thermique_interne_3_11, pMatrix3, pMatrix6);
    // Shape: G2202001_Demi_Ecran_thermique_interne_7 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_7_17 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_7", pBoolNode2);
    // Shape: G2202001_Demi_Ecran_thermique_interne_2 type: TGeoTrd2
    dx1 = 40.111503;
    dx2 = 47.787693;
    dy1 = 93.621885;
    dy2 = 167.243865;
    dz = 43.869630;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_2_12 =
        new TGeoTrd2("G2202001_Demi_Ecran_thermique_interne_2", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -295.010400;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix10 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix9 = new TGeoCombiTrans("", dx, dy, dz, pMatrix10);
    TGeoUnion* pBoolNode3 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_7_17, pG2202001_Demi_Ecran_thermique_interne_2_12, 0, pMatrix9);
    // Shape: G2202001_Demi_Ecran_thermique_interne_8 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_8_18 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_8", pBoolNode3);
    // Shape: G2202001_Demi_Ecran_thermique_interne_15 type: TGeoBBox
    dx = 4.130858;
    dy = 24.250000;
    dz = 18.005711;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_15_3 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_15", dx, dy, dz);
    // Shape: G2202001_Demi_Ecran_thermique_interne_14 type: TGeoTrd2
    dx1 = 18.005711;
    dx2 = 21.713076;
    dy1 = 24.250000;
    dy2 = 24.250000;
    dz = 21.220449;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_14_13 =
        new TGeoTrd2("G2202001_Demi_Ecran_thermique_interne_14", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -25.351300;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix13 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix12 = new TGeoCombiTrans("", dx, dy, dz, pMatrix13);
    TGeoUnion* pBoolNode4 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_15_3, pG2202001_Demi_Ecran_thermique_interne_14_13, 0, pMatrix12);
    // Shape: G2202001_Demi_Ecran_thermique_interne_13 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_13_19 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_13", pBoolNode4);
    // Shape: G2202001_Demi_Ecran_thermique_interne_12 type: TGeoTrd2
    dx1 = 21.713076;
    dx2 = 39.631380;
    dy1 = 24.250000;
    dy2 = 93.178137;
    dz = 102.403576;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_12_14 =
        new TGeoTrd2("G2202001_Demi_Ecran_thermique_interne_12", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -148.975300;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix16 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix15 = new TGeoCombiTrans("", dx, dy, dz, pMatrix16);
    TGeoUnion* pBoolNode5 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_13_19, pG2202001_Demi_Ecran_thermique_interne_12_14, 0, pMatrix15);
    // Shape: G2202001_Demi_Ecran_thermique_interne_11 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_11_20 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_11", pBoolNode5);
    // Shape: G2202001_Demi_Ecran_thermique_interne_10 type: TGeoTrd2
    dx1 = 39.631380;
    dx2 = 47.285784;
    dy1 = 93.178137;
    dy2 = 166.591161;
    dz = 43.80000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_10_15 =
        new TGeoTrd2("G2202001_Demi_Ecran_thermique_interne_10", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -295.134000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix19 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix18 = new TGeoCombiTrans("", dx, dy, dz, pMatrix19);
    TGeoUnion* pBoolNode6 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_11_20, pG2202001_Demi_Ecran_thermique_interne_10_15, 0, pMatrix18);
    // Shape: G2202001_Demi_Ecran_thermique_interne_9 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_9_21 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_9", pBoolNode6);
    TGeoSubtraction* pBoolNode7 = new TGeoSubtraction(
        pG2202001_Demi_Ecran_thermique_interne_8_18, pG2202001_Demi_Ecran_thermique_interne_9_21, 0, 0);
    // Shape: G2202001_Demi_Ecran_thermique_interne_1 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_1_22 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_1", pBoolNode7);
    // Shape: G2202001_Demi_Ecran_thermique_interne_26 type: TGeoBBox
    dx = 0.250000;
    dy = 6.000000;
    dz = 43.500000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_26_4 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_26", dx, dy, dz);
    // Combi transformation:
    dx = -339.130000;
    dy = 172.500000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix23 = new TGeoCombiTrans("");
    pMatrix23->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode8 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_1_22, pG2202001_Demi_Ecran_thermique_interne_26_4, 0, pMatrix23);
    // Shape: G2202001_Demi_Ecran_thermique_interne_28 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_28_23 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_28", pBoolNode8);
    // Shape: G2202001_Demi_Ecran_thermique_interne_25 type: TGeoBBox
    dx = 0.250000;
    dy = 82.250000;
    dz = 6.000000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_25_5 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_25", dx, dy, dz);
    // Combi transformation:
    dx = -339.130000;
    dy = 83.250000;
    dz = -52.500000;
    TGeoCombiTrans* pMatrix25 = new TGeoCombiTrans("");
    pMatrix25->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode9 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_28_23, pG2202001_Demi_Ecran_thermique_interne_25_5, 0, pMatrix25);
    // Shape: G2202001_Demi_Ecran_thermique_interne_29 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_29_24 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_29", pBoolNode9);
    // Shape: G2202001_Demi_Ecran_thermique_interne_23 type: TGeoBBox
    dx = 0.250000;
    dy = 82.250000;
    dz = 6.000000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_23_6 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_23", dx, dy, dz);
    // Combi transformation:
    dx = -339.130000;
    dy = 83.250000;
    dz = 52.500000;
    TGeoCombiTrans* pMatrix27 = new TGeoCombiTrans("");
    pMatrix27->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode10 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_29_24, pG2202001_Demi_Ecran_thermique_interne_23_6, 0, pMatrix27);
    // Shape: G2202001_Demi_Ecran_thermique_interne_30 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_30_25 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_30", pBoolNode10);
    // Shape: G2202001_Demi_Ecran_thermique_interne_19 type: TGeoBBox
    dx = 0.250000;
    dy = 6.000000;
    dz = 43.500000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_19_7 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_19", dx, dy, dz);
    // Combi transformation:
    dx = -339.130000;
    dy = -172.500000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix29 = new TGeoCombiTrans("");
    pMatrix29->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode11 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_30_25, pG2202001_Demi_Ecran_thermique_interne_19_7, 0, pMatrix29);
    // Shape: G2202001_Demi_Ecran_thermique_interne_31 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_31_26 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_31", pBoolNode11);
    // Shape: G2202001_Demi_Ecran_thermique_interne_18 type: TGeoBBox
    dx = 0.250000;
    dy = 82.250000;
    dz = 6.000000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_18_8 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_18", dx, dy, dz);
    // Combi transformation:
    dx = -339.130000;
    dy = -83.250000;
    dz = -52.500000;
    TGeoCombiTrans* pMatrix31 = new TGeoCombiTrans("");
    pMatrix31->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode12 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_31_26, pG2202001_Demi_Ecran_thermique_interne_18_8, 0, pMatrix31);
    // Shape: G2202001_Demi_Ecran_thermique_interne_32 type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_32_27 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_32", pBoolNode12);
    // Shape: G2202001_Demi_Ecran_thermique_interne_17 type: TGeoBBox
    dx = 0.250000;
    dy = 82.250000;
    dz = 6.000000;
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_17_9 =
        new TGeoBBox("G2202001_Demi_Ecran_thermique_interne_17", dx, dy, dz);
    // Combi transformation:
    dx = -339.130000;
    dy = -83.250000;
    dz = 52.500000;
    TGeoCombiTrans* pMatrix33 = new TGeoCombiTrans("");
    pMatrix33->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode13 = new TGeoUnion(
        pG2202001_Demi_Ecran_thermique_interne_32_27, pG2202001_Demi_Ecran_thermique_interne_17_9, 0, pMatrix33);
    // Shape: G2202001_Demi_Ecran_thermique_interne_PartBody type: TGeoCompositeShape
    TGeoShape* pG2202001_Demi_Ecran_thermique_interne_PartBody_28 =
        new TGeoCompositeShape("G2202001_Demi_Ecran_thermique_interne_PartBody", pBoolNode13);
    // Volume: G2202001_Demi_Ecran_thermique_interne
    TGeoVolume* p_Demi_Ecran_thermique_interne = new TGeoVolume(
        "G2202001_Demi_Ecran_thermique_interne", pG2202001_Demi_Ecran_thermique_interne_PartBody_28, pMedFe);
    p_Demi_Ecran_thermique_interne->SetLineColor(920);
    p_Demi_Ecran_thermique_interne->SetVisLeaves(kTRUE);

    // Add to top World volume
    tr[0] = 200.500000;
    tr[1] = 0.000000;
    tr[2] = 0.000000;
    rot[0] = 1.000000;
    rot[1] = 0.000000;
    rot[2] = 0.000000;
    rot[3] = 0.000000;
    rot[4] = 1.000000;
    rot[5] = 0.000000;
    rot[6] = 0.000000;
    rot[7] = 0.000000;
    rot[8] = 1.000000;
    TGeoHMatrix* pMatrix34 = new TGeoHMatrix("");
    pMatrix34->SetTranslation(tr);
    pMatrix34->SetRotation(rot);
    pMatrix34->SetBit(TGeoMatrix::kGeoTranslation);
    pMatrix34->SetBit(TGeoMatrix::kGeoRotation);

    pWorld->AddNode(p_Demi_Ecran_thermique_interne, 1, pMatrix34);
}

void ConstructEnceinteI(TGeoVolume* pWorld)
{
    // Enceinte part
    Double_t dx, dy, dz;
    Double_t dx1, dy1, dz1;
    Double_t dx2, dy2, dz2;
    Double_t rmin, rmax, phi1, phi2;
    Double_t thx, thy, thz, phx, phy, phz;

    // -- Mixture: Air
    TGeoMedium* pMed2 = gGeoManager->GetMedium("Air");
    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");
    TGeoMedium* pMedAl = gGeoManager->GetMedium("aluminium");

    // local transformation
    Double_t tr[3], rot[9];

    // Shape: G2402001_Enceinte_interne_2 type: TGeoBBox
    dx = 7.424479;
    dy = 22.250000;
    dz = 16.000000;
    TGeoShape* pG2402001_Enceinte_interne_2_2 = new TGeoBBox("G2402001_Enceinte_interne_2", dx, dy, dz);
    // Shape: G2402001_Enceinte_interne_3 type: TGeoTrd2
    dx1 = 16.000000;
    dx2 = 19.994047;
    dy1 = 22.250000;
    dy2 = 22.250000;
    dz = 22.826083;
    TGeoShape* pG2402001_Enceinte_interne_3_9 = new TGeoTrd2("G2402001_Enceinte_interne_3", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -30.250500;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix4 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix3 = new TGeoCombiTrans("", dx, dy, dz, pMatrix4);
    TGeoUnion* pBoolNode1 = new TGeoUnion(pG2402001_Enceinte_interne_2_2, pG2402001_Enceinte_interne_3_9, 0, pMatrix3);
    // Shape: G2402001_Enceinte_interne_26 type: TGeoCompositeShape
    // D. Kresan: skip entrance box - use Trd2 shape directly
    TGeoShape* pG2402001_Enceinte_interne_26_15 =
        pG2402001_Enceinte_interne_3_9; // new TGeoCompositeShape("G2402001_Enceinte_interne_26", pBoolNode1);
    // Shape: G2402001_Enceinte_interne_4 type: TGeoTrd2
    dx1 = 19.994047;
    dx2 = 37.710887;
    dy1 = 23.250000;
    dy2 = 91.403144;
    dz = 101.252203;
    TGeoShape* pG2402001_Enceinte_interne_4_10 = new TGeoTrd2("G2402001_Enceinte_interne_4", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -154.328800;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix7 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix6 = new TGeoCombiTrans("", dx, dy, dz, pMatrix7);
    // D. Kresan: skip entrance box - use Trd2 shape directly - use pMatrix3 instead of 0
    TGeoUnion* pBoolNode2 =
        new TGeoUnion(pG2402001_Enceinte_interne_26_15, pG2402001_Enceinte_interne_4_10, pMatrix3, pMatrix6);
    // Shape: G2402001_Enceinte_interne_27 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_27_16 = new TGeoCompositeShape("G2402001_Enceinte_interne_27", pBoolNode2);
    // Shape: G2402001_Enceinte_interne_5 type: TGeoTrd2
    dx1 = 37.710887;
    dx2 = 46.328008;
    dy1 = 91.403144;
    dy2 = 174.049542;
    dz = 49.247071;
    TGeoShape* pG2402001_Enceinte_interne_5_11 = new TGeoTrd2("G2402001_Enceinte_interne_5", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -304.828100;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix10 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix9 = new TGeoCombiTrans("", dx, dy, dz, pMatrix10);
    TGeoUnion* pBoolNode3 =
        new TGeoUnion(pG2402001_Enceinte_interne_27_16, pG2402001_Enceinte_interne_5_11, 0, pMatrix9);
    // Shape: G2402001_Enceinte_interne_28 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_28_17 = new TGeoCompositeShape("G2402001_Enceinte_interne_28", pBoolNode3);
    // Shape: G2402001_Enceinte_interne_6 type: TGeoBBox
    dx = 16.250000;
    dx = 16.3;
    dy = 207.116855;
    dz = 47.852199;
    // dz= 47.89;
    TGeoShape* pG2402001_Enceinte_interne_6_3 = new TGeoBBox("G2402001_Enceinte_interne_6", dx, dy, dz);
    // Combi transformation:
    dx = -370.325200;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix12 = new TGeoCombiTrans("");
    pMatrix12->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode4 =
        new TGeoUnion(pG2402001_Enceinte_interne_28_17, pG2402001_Enceinte_interne_6_3, 0, pMatrix12);
    // Shape: G2402001_Enceinte_interne_29 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_29_18 = new TGeoCompositeShape("G2402001_Enceinte_interne_29", pBoolNode4);

    // Shape: G2402001_Enceinte_interne_7 type: TGeoBBox
    dx = 4.000000;
    dy = 234.116855;
    dz = 74.852199;
    dz = 74.89;
    TGeoShape* pG2402001_Enceinte_interne_7_4 = new TGeoBBox("G2402001_Enceinte_interne_7", dx, dy, dz);
    // Combi transformation:
    dx = -390.575200;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix14 = new TGeoCombiTrans("");
    pMatrix14->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode5 =
        new TGeoUnion(pG2402001_Enceinte_interne_29_18, pG2402001_Enceinte_interne_7_4, 0, pMatrix14);
    // Shape: G2402001_Enceinte_interne_30 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_30_19 = new TGeoCompositeShape("G2402001_Enceinte_interne_30", pBoolNode5);
    // Shape: G2402001_Enceinte_interne_8 type: TGeoBBox
    dx = 2.500000;
    dy = 222.116855;
    dz = 62.852199;
    dz = 62.89;
    TGeoShape* pG2402001_Enceinte_interne_8_5 = new TGeoBBox("G2402001_Enceinte_interne_8", dx, dy, dz);
    // Combi transformation:
    dx = -368.575500;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix16 = new TGeoCombiTrans("");
    pMatrix16->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode6 =
        new TGeoUnion(pG2402001_Enceinte_interne_30_19, pG2402001_Enceinte_interne_8_5, 0, pMatrix16);
    // Shape: G2402001_Enceinte_interne_31 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_31_20 = new TGeoCompositeShape("G2402001_Enceinte_interne_31", pBoolNode6);
    // Shape: G2402001_Enceinte_interne_15 type: TGeoBBox
    dx = 7.450000;
    dy = 20.250000;
    dz = 11.000000;
    TGeoShape* pG2402001_Enceinte_interne_15_6 = new TGeoBBox("G2402001_Enceinte_interne_15", dx, dy, dz);
    // Shape: G2402001_Enceinte_interne_16 type: TGeoTrd2
    dx1 = 11.000000;
    dx2 = 15.017930;
    dy1 = 20.250000;
    dy2 = 20.250000;
    dz = 23.246371;
    TGeoShape* pG2402001_Enceinte_interne_16_12 = new TGeoTrd2("G2402001_Enceinte_interne_16", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -30.496400;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix19 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix18 = new TGeoCombiTrans("", dx, dy, dz, pMatrix19);
    TGeoUnion* pBoolNode7 =
        new TGeoUnion(pG2402001_Enceinte_interne_15_6, pG2402001_Enceinte_interne_16_12, 0, pMatrix18);
    // Shape: G2402001_Enceinte_interne_24 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_24_21 = new TGeoCompositeShape("G2402001_Enceinte_interne_24", pBoolNode7);
    // Shape: G2402001_Enceinte_interne_17 type: TGeoTrd2
    dx1 = 15.017930;
    dx2 = 32.822508;
    dy1 = 20.250000;
    dy2 = 88.740656;
    dz = 101.753630;
    TGeoShape* pG2402001_Enceinte_interne_17_13 = new TGeoTrd2("G2402001_Enceinte_interne_17", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -155.496400;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix22 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix21 = new TGeoCombiTrans("", dx, dy, dz, pMatrix22);
    TGeoUnion* pBoolNode8 =
        new TGeoUnion(pG2402001_Enceinte_interne_24_21, pG2402001_Enceinte_interne_17_13, 0, pMatrix21);
    // Shape: G2402001_Enceinte_interne_21 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_21_22 = new TGeoCompositeShape("G2402001_Enceinte_interne_21", pBoolNode8);
    // Shape: G2402001_Enceinte_interne_18 type: TGeoTrd2
    dx1 = 32.822508;
    dx2 = 41.308909;
    dy1 = 88.740656;
    dy2 = 170.133320;
    dz = 48.500000;
    TGeoShape* pG2402001_Enceinte_interne_18_14 = new TGeoTrd2("G2402001_Enceinte_interne_18", dx1, dx2, dy1, dy2, dz);
    // Combi transformation:
    dx = -305.750000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix25 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix24 = new TGeoCombiTrans("", dx, dy, dz, pMatrix25);
    TGeoUnion* pBoolNode9 =
        new TGeoUnion(pG2402001_Enceinte_interne_21_22, pG2402001_Enceinte_interne_18_14, 0, pMatrix24);
    // Shape: G2402001_Enceinte_interne_22 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_22_23 = new TGeoCompositeShape("G2402001_Enceinte_interne_22", pBoolNode9);
    // Shape: G2402001_Enceinte_interne_19 type: TGeoBBox
    dx = 1.500000;
    dy = 170.133320;
    dz = 41.308909;
    TGeoShape* pG2402001_Enceinte_interne_19_7 = new TGeoBBox("G2402001_Enceinte_interne_19", dx, dy, dz);
    // Combi transformation:
    dx = -355.750000;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix27 = new TGeoCombiTrans("");
    pMatrix27->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode10 =
        new TGeoUnion(pG2402001_Enceinte_interne_22_23, pG2402001_Enceinte_interne_19_7, 0, pMatrix27);
    // Shape: G2402001_Enceinte_interne_23 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_23_24 = new TGeoCompositeShape("G2402001_Enceinte_interne_23", pBoolNode10);
    // Shape: G2402001_Enceinte_interne_20 type: TGeoBBox
    dx = 22.1;
    dy = 204.116855;
    dz = 44.852199;
    TGeoShape* pG2402001_Enceinte_interne_20_8 = new TGeoBBox("G2402001_Enceinte_interne_20", dx, dy, dz);
    // Combi transformation:
    dx = -376.250000;
    //   dx = -380.0;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix29 = new TGeoCombiTrans("");
    pMatrix29->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode11 =
        new TGeoUnion(pG2402001_Enceinte_interne_23_24, pG2402001_Enceinte_interne_20_8, 0, pMatrix29);
    // Shape: G2402001_Enceinte_interne_25 type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_25_25 = new TGeoCompositeShape("G2402001_Enceinte_interne_25", pBoolNode11);
    TGeoSubtraction* pBoolNode12 =
        new TGeoSubtraction(pG2402001_Enceinte_interne_31_20, pG2402001_Enceinte_interne_25_25, 0, 0);
    // Shape: G2402001_Enceinte_interne_PartBody type: TGeoCompositeShape
    TGeoShape* pG2402001_Enceinte_interne_PartBody_26 =
        new TGeoCompositeShape("G2402001_Enceinte_interne_PartBody", pBoolNode12);
    // Volume: G2402001_Enceinte_interne
    TGeoVolume* pG2402001_Enceinte_interne =
        new TGeoVolume("G2402001_Enceinte_interne", pG2402001_Enceinte_interne_PartBody_26, pMedAl);
    pG2402001_Enceinte_interne->SetLineColor(920);
    pG2402001_Enceinte_interne->SetVisLeaves(kTRUE);

    // Add the Nurbs like node to the world

    tr[0] = 185.000000 + 18.5;
    tr[1] = 0.000000;
    tr[2] = 0.000000;
    rot[0] = 1.000000;
    rot[1] = 0.000000;
    rot[2] = 0.000000;
    rot[3] = 0.000000;
    rot[4] = 1.000000;
    rot[5] = 0.000000;
    rot[6] = 0.000000;
    rot[7] = 0.000000;
    rot[8] = 1.000000;
    TGeoHMatrix* pMatrix32 = new TGeoHMatrix("");
    pMatrix32->SetTranslation(tr);
    pMatrix32->SetRotation(rot);
    pMatrix32->SetBit(TGeoMatrix::kGeoTranslation);
    pMatrix32->SetBit(TGeoMatrix::kGeoRotation);

    pWorld->AddNode(pG2402001_Enceinte_interne, 1, pMatrix32);
}

void ConstructEnceinteE(TGeoVolume* pWorld)
{
    // Enceinte Extern part
    Double_t dx, dy, dz;
    Double_t dx1, dy1, dz1;
    Double_t dx2, dy2, dz2;
    Double_t rmin, rmax, phi1, phi2;
    Double_t thx, thy, thz, phx, phy, phz;

    // -- Mixture: Air
    TGeoMedium* pMed2 = gGeoManager->GetMedium("Air");
    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");
    TGeoMedium* pMedAl = gGeoManager->GetMedium("aluminium");

    Double_t tr[3], rot[9];

    // Geometry Description
    rmin = 120.000000;
    rmax = 121.500000;
    dz = 115.000000;
    phi1 = 26.424390;
    phi2 = 90.000000;
    TGeoShape* pGEnceinte_externe_131_2 = new TGeoTubeSeg("GEnceinte_externe_131", rmin, rmax, dz, phi1, phi2);
    // Shape: GEnceinte_externe_132 type: TGeoTubeSeg
    rmin = 455.000000;
    rmax = 456.500000;
    dz = 115.000000;
    phi1 = 333.575610;
    phi2 = 386.424390;
    TGeoShape* pGEnceinte_externe_132_3 = new TGeoTubeSeg("GEnceinte_externe_132", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.000000;
    dy = -149.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix3 = new TGeoCombiTrans("");
    pMatrix3->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode1 = new TGeoUnion(pGEnceinte_externe_131_2, pGEnceinte_externe_132_3, 0, pMatrix3);
    // Shape: GEnceinte_externe_133 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_133_49 = new TGeoCompositeShape("GEnceinte_externe_133", pBoolNode1);
    // Shape: GEnceinte_externe_130 type: TGeoTubeSeg
    rmin = 120.000000;
    rmax = 121.500000;
    dz = 115.000000;
    phi1 = 270.000000;
    phi2 = 333.575610;
    TGeoShape* pGEnceinte_externe_130_4 = new TGeoTubeSeg("GEnceinte_externe_130", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -298.161000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix5 = new TGeoCombiTrans("");
    pMatrix5->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode2 = new TGeoUnion(pGEnceinte_externe_133_49, pGEnceinte_externe_130_4, 0, pMatrix5);
    // Shape: GEnceinte_externe_134 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_134_50 = new TGeoCompositeShape("GEnceinte_externe_134", pBoolNode2);
    // Shape: GEnceinte_externe_128 type: TGeoTorus
    Double_t r = 375.000000;
    rmin = 80.000000;
    rmax = 81.500000;
    phi1 = 333.575610;
    Double_t dphi = 52.848780;
    TGeoShape* pGEnceinte_externe_128_39 = new TGeoTorus("GEnceinte_externe_128", r, rmin, rmax, phi1, dphi);
    // Shape: GEnceinte_externe_129 type: TGeoTubeSeg
    rmin = 375.000000;
    rmax = 456.500000;
    dz = 40.750000;
    phi1 = 330.000000;
    phi2 = 390.000000;
    TGeoShape* pGEnceinte_externe_129_5 = new TGeoTubeSeg("GEnceinte_externe_129", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 40.750000;
    TGeoCombiTrans* pMatrix19 = new TGeoCombiTrans("");
    pMatrix19->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode3 =
        new TGeoIntersection(pGEnceinte_externe_128_39, pGEnceinte_externe_129_5, 0, pMatrix19);
    // Shape: GEnceinte_externe_135 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_135_56 = new TGeoCompositeShape("GEnceinte_externe_135", pBoolNode3);
    // Combi transformation:
    dx = -300.000000;
    dy = -149.080500;
    dz = 115.000000;
    TGeoCombiTrans* pMatrix35 = new TGeoCombiTrans("");
    pMatrix35->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode4 = new TGeoUnion(pGEnceinte_externe_134_50, pGEnceinte_externe_135_56, 0, pMatrix35);
    // Shape: GEnceinte_externe_136 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_136_62 = new TGeoCompositeShape("GEnceinte_externe_136", pBoolNode4);
    // Shape: GEnceinte_externe_158 type: TGeoBBox
    dx = 86.200000;
    dy = 121.500000;
    dz = 108.805970;
    TGeoShape* pGEnceinte_externe_158_41 = new TGeoBBox("GEnceinte_externe_158", dx, dy, dz);
    // Shape: GEnceinte_externe_157 type: TGeoTubeSeg
    rmin = 40.000000;
    rmax = 122.500000;
    dz = 45.000000;
    phi1 = 26.424390;
    phi2 = 90.000000;
    TGeoShape* pGEnceinte_externe_157_6 = new TGeoTubeSeg("GEnceinte_externe_157", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 45.000000;
    dy = -0.003800;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix22 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix21 = new TGeoCombiTrans("", dx, dy, dz, pMatrix22);
    TGeoIntersection* pBoolNode5 =
        new TGeoIntersection(pGEnceinte_externe_158_41, pGEnceinte_externe_157_6, 0, pMatrix21);
    // Shape: GEnceinte_externe_156 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_156_57 = new TGeoCompositeShape("GEnceinte_externe_156", pBoolNode5);
    // Shape: GEnceinte_externe_155 type: TGeoBBox
    dx = 84.900000;
    dy = 120.000000;
    dz = 122.000000;
    TGeoShape* pGEnceinte_externe_155_42 = new TGeoBBox("GEnceinte_externe_155", dx, dy, dz);
    TGeoSubtraction* pBoolNode6 = new TGeoSubtraction(pGEnceinte_externe_156_57, pGEnceinte_externe_155_42, 0, 0);
    // Shape: GEnceinte_externe_154 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_154_66 = new TGeoCompositeShape("GEnceinte_externe_154", pBoolNode6);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 115.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix75 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix74 = new TGeoCombiTrans("", dx, dy, dz, pMatrix75);
    TGeoUnion* pBoolNode7 = new TGeoUnion(pGEnceinte_externe_136_62, pGEnceinte_externe_154_66, 0, pMatrix74);
    // Shape: GEnceinte_externe_159 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_159_76 = new TGeoCompositeShape("GEnceinte_externe_159", pBoolNode7);
    // Shape: GEnceinte_externe_153 type: TGeoBBox
    dx = 86.200000;
    dy = 121.500000;
    dz = 108.805970;
    TGeoShape* pGEnceinte_externe_153_43 = new TGeoBBox("GEnceinte_externe_153", dx, dy, dz);
    // Shape: GEnceinte_externe_152 type: TGeoTubeSeg
    rmin = 40.000000;
    rmax = 122.500000;
    dz = 45.000000;
    phi1 = 270.000000;
    phi2 = 333.575610;
    TGeoShape* pGEnceinte_externe_152_7 = new TGeoTubeSeg("GEnceinte_externe_152", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 45.000000;
    dy = 0.003800;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix25 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix24 = new TGeoCombiTrans("", dx, dy, dz, pMatrix25);
    TGeoIntersection* pBoolNode8 =
        new TGeoIntersection(pGEnceinte_externe_153_43, pGEnceinte_externe_152_7, 0, pMatrix24);
    // Shape: GEnceinte_externe_151 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_151_58 = new TGeoCompositeShape("GEnceinte_externe_151", pBoolNode8);
    // Shape: GEnceinte_externe_150 type: TGeoBBox
    dx = 84.900000;
    dy = 120.000000;
    dz = 122.000000;
    TGeoShape* pGEnceinte_externe_150_44 = new TGeoBBox("GEnceinte_externe_150", dx, dy, dz);
    TGeoSubtraction* pBoolNode9 = new TGeoSubtraction(pGEnceinte_externe_151_58, pGEnceinte_externe_150_44, 0, 0);
    // Shape: GEnceinte_externe_149 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_149_67 = new TGeoCompositeShape("GEnceinte_externe_149", pBoolNode9);
    // Combi transformation:
    dx = 0.000000;
    dy = -298.161000;
    dz = 115.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix78 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix77 = new TGeoCombiTrans("", dx, dy, dz, pMatrix78);
    TGeoUnion* pBoolNode10 = new TGeoUnion(pGEnceinte_externe_159_76, pGEnceinte_externe_149_67, 0, pMatrix77);
    // Shape: GEnceinte_externe_160 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_160_77 = new TGeoCompositeShape("GEnceinte_externe_160", pBoolNode10);
    // Shape: GEnceinte_externe_163 type: TGeoBBox
    dx = 19.500000;
    dy = 30.000000;
    dz = 25.500000;
    TGeoShape* pGEnceinte_externe_163_1 = new TGeoBBox("GEnceinte_externe_163", dx, dy, dz);
    // Shape: GEnceinte_externe_162 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 19.500000;
    dz = 30.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_162_8 = new TGeoTubeSeg("GEnceinte_externe_162", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -25.500000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 0.000000;
    phy = 0.000000;
    thz = 90.000000;
    phz = 270.000000;
    TGeoRotation* pMatrix8 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix7 = new TGeoCombiTrans("", dx, dy, dz, pMatrix8);
    TGeoUnion* pBoolNode11 = new TGeoUnion(pGEnceinte_externe_163_1, pGEnceinte_externe_162_8, 0, pMatrix7);
    // Shape: GEnceinte_externe_161 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_161_51 = new TGeoCompositeShape("GEnceinte_externe_161", pBoolNode11);
    // Combi transformation:
    dx = 24.500000;
    dy = 90.000000;
    dz = 158.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix81 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix80 = new TGeoCombiTrans("", dx, dy, dz, pMatrix81);
    TGeoSubtraction* pBoolNode12 =
        new TGeoSubtraction(pGEnceinte_externe_160_77, pGEnceinte_externe_161_51, 0, pMatrix80);
    // Shape: GEnceinte_externe_167 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_167_78 = new TGeoCompositeShape("GEnceinte_externe_167", pBoolNode12);
    // Shape: GEnceinte_externe_3 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 5.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_3_9 = new TGeoTubeSeg("GEnceinte_externe_3", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 133.500000;
    dy = -8.080500;
    dz = 78.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix84 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix83 = new TGeoCombiTrans("", dx, dy, dz, pMatrix84);
    TGeoSubtraction* pBoolNode13 = new TGeoSubtraction(pGEnceinte_externe_167_78, pGEnceinte_externe_3_9, 0, pMatrix83);
    // Shape: GEnceinte_externe_166 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_166_79 = new TGeoCompositeShape("GEnceinte_externe_166", pBoolNode13);
    // Shape: GEnceinte_externe_5 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 15.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_5_10 = new TGeoTubeSeg("GEnceinte_externe_5", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 100.000000;
    dy = 70.919500;
    dz = -37.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix87 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix86 = new TGeoCombiTrans("", dx, dy, dz, pMatrix87);
    TGeoSubtraction* pBoolNode14 =
        new TGeoSubtraction(pGEnceinte_externe_166_79, pGEnceinte_externe_5_10, 0, pMatrix86);
    // Shape: GEnceinte_externe_165 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_165_80 = new TGeoCompositeShape("GEnceinte_externe_165", pBoolNode14);
    // Shape: GEnceinte_externe_7 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 5.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_7_11 = new TGeoTubeSeg("GEnceinte_externe_7", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 133.000000;
    dy = -290.080500;
    dz = 78.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix90 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix89 = new TGeoCombiTrans("", dx, dy, dz, pMatrix90);
    TGeoSubtraction* pBoolNode15 =
        new TGeoSubtraction(pGEnceinte_externe_165_80, pGEnceinte_externe_7_11, 0, pMatrix89);
    // Shape: GEnceinte_externe_168 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_168_81 = new TGeoCompositeShape("GEnceinte_externe_168", pBoolNode15);
    // Shape: GEnceinte_externe_9 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 15.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_9_12 = new TGeoTubeSeg("GEnceinte_externe_9", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 100.000000;
    dy = -369.080500;
    dz = -37.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix93 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix92 = new TGeoCombiTrans("", dx, dy, dz, pMatrix93);
    TGeoSubtraction* pBoolNode16 =
        new TGeoSubtraction(pGEnceinte_externe_168_81, pGEnceinte_externe_9_12, 0, pMatrix92);
    // Shape: GEnceinte_externe_169 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_169_82 = new TGeoCompositeShape("GEnceinte_externe_169", pBoolNode16);
    // Shape: GEnceinte_externe_11 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 28.500000;
    dz = 50.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_11_13 = new TGeoTubeSeg("GEnceinte_externe_11", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 120.000000;
    dy = -149.080500;
    dz = 158.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix96 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix95 = new TGeoCombiTrans("", dx, dy, dz, pMatrix96);
    TGeoSubtraction* pBoolNode17 =
        new TGeoSubtraction(pGEnceinte_externe_169_82, pGEnceinte_externe_11_13, 0, pMatrix95);
    // Shape: GEnceinte_externe_170 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_170_83 = new TGeoCompositeShape("GEnceinte_externe_170", pBoolNode17);
    // Shape: GEnceinte_externe_13 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 28.500000;
    dz = 14.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_13_14 = new TGeoTubeSeg("GEnceinte_externe_13", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 125.000000;
    dy = -310.080500;
    dz = -37.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix99 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix98 = new TGeoCombiTrans("", dx, dy, dz, pMatrix99);
    TGeoSubtraction* pBoolNode18 =
        new TGeoSubtraction(pGEnceinte_externe_170_83, pGEnceinte_externe_13_14, 0, pMatrix98);
    // Shape: GEnceinte_externe_171 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_171_84 = new TGeoCompositeShape("GEnceinte_externe_171", pBoolNode18);
    // Shape: GEnceinte_externe_15 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 28.500000;
    dz = 30.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_15_15 = new TGeoTubeSeg("GEnceinte_externe_15", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 130.000000;
    dy = 11.919500;
    dz = -37.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix102 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix101 = new TGeoCombiTrans("", dx, dy, dz, pMatrix102);
    TGeoSubtraction* pBoolNode19 =
        new TGeoSubtraction(pGEnceinte_externe_171_84, pGEnceinte_externe_15_15, 0, pMatrix101);
    // Shape: GEnceinte_externe_173 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_173_85 = new TGeoCompositeShape("GEnceinte_externe_173", pBoolNode19);
    // Shape: GEnceinte_externe_17 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 13.050000;
    dz = 150.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_17_16 = new TGeoTubeSeg("GEnceinte_externe_17", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -37.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 150.000000;
    thz = 90.000000;
    phz = 240.000000;
    TGeoRotation* pMatrix105 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix104 = new TGeoCombiTrans("", dx, dy, dz, pMatrix105);
    TGeoSubtraction* pBoolNode20 =
        new TGeoSubtraction(pGEnceinte_externe_173_85, pGEnceinte_externe_17_16, 0, pMatrix104);
    // Shape: GEnceinte_externe_172 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_172_86 = new TGeoCompositeShape("GEnceinte_externe_172", pBoolNode20);
    // Shape: GEnceinte_externe_97 type: TGeoTubeSeg
    rmin = 365.000000;
    rmax = 375.000000;
    dz = 0.750000;
    phi1 = 0.000000;
    phi2 = 26.424390;
    TGeoShape* pGEnceinte_externe_97_17 = new TGeoTubeSeg("GEnceinte_externe_97", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.000000;
    dy = -149.080500;
    dz = 195.750000;
    TGeoCombiTrans* pMatrix107 = new TGeoCombiTrans("");
    pMatrix107->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode21 = new TGeoUnion(pGEnceinte_externe_172_86, pGEnceinte_externe_97_17, 0, pMatrix107);
    // Shape: GEnceinte_externe_174 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_174_87 = new TGeoCompositeShape("GEnceinte_externe_174", pBoolNode21);
    // Shape: GEnceinte_externe_96 type: TGeoTubeSeg
    rmin = 365.000000;
    rmax = 375.000000;
    dz = 0.750000;
    phi1 = 333.575610;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_96_18 = new TGeoTubeSeg("GEnceinte_externe_96", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.000000;
    dy = -149.080500;
    dz = 195.750000;
    TGeoCombiTrans* pMatrix109 = new TGeoCombiTrans("");
    pMatrix109->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode22 = new TGeoUnion(pGEnceinte_externe_174_87, pGEnceinte_externe_96_18, 0, pMatrix109);
    // Shape: GEnceinte_externe_175 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_175_88 = new TGeoCompositeShape("GEnceinte_externe_175", pBoolNode22);
    // Shape: GEnceinte_externe_94 type: TGeoTubeSeg
    rmin = 30.000000;
    rmax = 40.000000;
    dz = 0.750000;
    phi1 = 26.424390;
    phi2 = 90.000000;
    TGeoShape* pGEnceinte_externe_94_19 = new TGeoTubeSeg("GEnceinte_externe_94", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 195.750000;
    TGeoCombiTrans* pMatrix111 = new TGeoCombiTrans("");
    pMatrix111->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode23 = new TGeoUnion(pGEnceinte_externe_175_88, pGEnceinte_externe_94_19, 0, pMatrix111);
    // Shape: GEnceinte_externe_176 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_176_89 = new TGeoCompositeShape("GEnceinte_externe_176", pBoolNode23);
    // Shape: GEnceinte_externe_92 type: TGeoTubeSeg
    rmin = 30.000000;
    rmax = 40.000000;
    dz = 0.750000;
    phi1 = 270.000000;
    phi2 = 333.575610;
    TGeoShape* pGEnceinte_externe_92_20 = new TGeoTubeSeg("GEnceinte_externe_92", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -298.161000;
    dz = 195.750000;
    TGeoCombiTrans* pMatrix113 = new TGeoCombiTrans("");
    pMatrix113->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode24 = new TGeoUnion(pGEnceinte_externe_176_89, pGEnceinte_externe_92_20, 0, pMatrix113);
    // Shape: GEnceinte_externe_177 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_177_90 = new TGeoCompositeShape("GEnceinte_externe_177", pBoolNode24);
    // Shape: GEnceinte_externe_87 type: TGeoTubeSeg
    rmin = 30.000000;
    rmax = 40.000000;
    dz = 0.750000;
    phi1 = 26.424390;
    phi2 = 90.000000;
    TGeoShape* pGEnceinte_externe_87_21 = new TGeoTubeSeg("GEnceinte_externe_87", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -298.161000;
    dz = 195.750000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix116 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix115 = new TGeoCombiTrans("", dx, dy, dz, pMatrix116);
    TGeoUnion* pBoolNode25 = new TGeoUnion(pGEnceinte_externe_177_90, pGEnceinte_externe_87_21, 0, pMatrix115);
    // Shape: GEnceinte_externe_178 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_178_91 = new TGeoCompositeShape("GEnceinte_externe_178", pBoolNode25);
    // Shape: GEnceinte_externe_85 type: TGeoTubeSeg
    rmin = 30.000000;
    rmax = 40.000000;
    dz = 0.750000;
    phi1 = 270.000000;
    phi2 = 333.575610;
    TGeoShape* pGEnceinte_externe_85_22 = new TGeoTubeSeg("GEnceinte_externe_85", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 195.750000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix119 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix118 = new TGeoCombiTrans("", dx, dy, dz, pMatrix119);
    TGeoUnion* pBoolNode26 = new TGeoUnion(pGEnceinte_externe_178_91, pGEnceinte_externe_85_22, 0, pMatrix118);
    // Shape: GEnceinte_externe_179 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_179_92 = new TGeoCompositeShape("GEnceinte_externe_179", pBoolNode26);
    // Shape: GEnceinte_externe_90 type: TGeoTubeSeg
    rmin = 365.000000;
    rmax = 375.000000;
    dz = 0.750000;
    phi1 = 0.000000;
    phi2 = 26.424390;
    TGeoShape* pGEnceinte_externe_90_23 = new TGeoTubeSeg("GEnceinte_externe_90", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 300.000000;
    dy = -149.080500;
    dz = 195.750000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix122 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix121 = new TGeoCombiTrans("", dx, dy, dz, pMatrix122);
    TGeoUnion* pBoolNode27 = new TGeoUnion(pGEnceinte_externe_179_92, pGEnceinte_externe_90_23, 0, pMatrix121);
    // Shape: GEnceinte_externe_180 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_180_93 = new TGeoCompositeShape("GEnceinte_externe_180", pBoolNode27);
    // Shape: GEnceinte_externe_89 type: TGeoTubeSeg
    rmin = 365.000000;
    rmax = 375.000000;
    dz = 0.750000;
    phi1 = 333.575610;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_89_24 = new TGeoTubeSeg("GEnceinte_externe_89", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 300.000000;
    dy = -149.080500;
    dz = 195.750000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix125 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix124 = new TGeoCombiTrans("", dx, dy, dz, pMatrix125);
    TGeoUnion* pBoolNode28 = new TGeoUnion(pGEnceinte_externe_180_93, pGEnceinte_externe_89_24, 0, pMatrix124);
    // Shape: GEnceinte_externe_181 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_181_94 = new TGeoCompositeShape("GEnceinte_externe_181", pBoolNode28);
    // Shape: GEnceinte_externe_80 type: TGeoTubeSeg
    rmin = 120.000000;
    rmax = 121.500000;
    dz = 115.000000;
    phi1 = 90.000000;
    phi2 = 153.576000;
    TGeoShape* pGEnceinte_externe_80_25 = new TGeoTubeSeg("GEnceinte_externe_80", rmin, rmax, dz, phi1, phi2);
    // Shape: GEnceinte_externe_79 type: TGeoTubeSeg
    rmin = 455.000000;
    rmax = 456.500000;
    dz = 115.000000;
    phi1 = 333.575610;
    phi2 = 386.424390;
    TGeoShape* pGEnceinte_externe_79_26 = new TGeoTubeSeg("GEnceinte_externe_79", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 300.000000;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix11 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix10 = new TGeoCombiTrans("", dx, dy, dz, pMatrix11);
    TGeoUnion* pBoolNode29 = new TGeoUnion(pGEnceinte_externe_80_25, pGEnceinte_externe_79_26, 0, pMatrix10);
    // Shape: GEnceinte_externe_184 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_184_52 = new TGeoCompositeShape("GEnceinte_externe_184", pBoolNode29);
    // Shape: GEnceinte_externe_77 type: TGeoTubeSeg
    rmin = 120.000000;
    rmax = 121.500000;
    dz = 115.000000;
    phi1 = 206.424000;
    phi2 = 270.000000;
    TGeoShape* pGEnceinte_externe_77_27 = new TGeoTubeSeg("GEnceinte_externe_77", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -298.100000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix13 = new TGeoCombiTrans("");
    pMatrix13->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode30 = new TGeoUnion(pGEnceinte_externe_184_52, pGEnceinte_externe_77_27, 0, pMatrix13);
    // Shape: GEnceinte_externe_185 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_185_53 = new TGeoCompositeShape("GEnceinte_externe_185", pBoolNode30);
    // Shape: GEnceinte_externe_82 type: TGeoTubeSeg
    rmin = 375.000000;
    rmax = 456.500000;
    dz = 40.750000;
    phi1 = 330.000000;
    phi2 = 390.000000;
    TGeoShape* pGEnceinte_externe_82_28 = new TGeoTubeSeg("GEnceinte_externe_82", rmin, rmax, dz, phi1, phi2);
    // Shape: GEnceinte_externe_83 type: TGeoTorus
    r = 375.000000;
    rmin = 80.000000;
    rmax = 81.500000;
    phi1 = 333.575610;
    dphi = 52.848780;
    TGeoShape* pGEnceinte_externe_83_40 = new TGeoTorus("GEnceinte_externe_83", r, rmin, rmax, phi1, dphi);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -40.750000;
    TGeoCombiTrans* pMatrix27 = new TGeoCombiTrans("");
    pMatrix27->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode31 =
        new TGeoIntersection(pGEnceinte_externe_82_28, pGEnceinte_externe_83_40, 0, pMatrix27);
    // Shape: GEnceinte_externe_183 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_183_59 = new TGeoCompositeShape("GEnceinte_externe_183", pBoolNode31);
    // Combi transformation:
    dx = 300.000000;
    dy = -149.080500;
    dz = 155.750000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 0.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix38 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix37 = new TGeoCombiTrans("", dx, dy, dz, pMatrix38);
    TGeoUnion* pBoolNode32 = new TGeoUnion(pGEnceinte_externe_185_53, pGEnceinte_externe_183_59, 0, pMatrix37);
    // Shape: GEnceinte_externe_186 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_186_63 = new TGeoCompositeShape("GEnceinte_externe_186", pBoolNode32);
    // Shape: GEnceinte_externe_74 type: TGeoBBox
    dx = 86.200000;
    dy = 121.500000;
    dz = 108.805970;
    TGeoShape* pGEnceinte_externe_74_45 = new TGeoBBox("GEnceinte_externe_74", dx, dy, dz);
    // Shape: GEnceinte_externe_71 type: TGeoBBox
    dx = 84.900000;
    dy = 120.000000;
    dz = 122.000000;
    TGeoShape* pGEnceinte_externe_71_46 = new TGeoBBox("GEnceinte_externe_71", dx, dy, dz);
    TGeoSubtraction* pBoolNode33 = new TGeoSubtraction(pGEnceinte_externe_74_45, pGEnceinte_externe_71_46, 0, 0);
    // Shape: GEnceinte_externe_187 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_187_54 = new TGeoCompositeShape("GEnceinte_externe_187", pBoolNode33);
    // Shape: GEnceinte_externe_73 type: TGeoTubeSeg
    rmin = 40.000000;
    rmax = 122.500000;
    dz = 45.000000;
    phi1 = 270.000000;
    phi2 = 333.575610;
    TGeoShape* pGEnceinte_externe_73_29 = new TGeoTubeSeg("GEnceinte_externe_73", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 45.000000;
    dy = 0.003800;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix30 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix29 = new TGeoCombiTrans("", dx, dy, dz, pMatrix30);
    TGeoIntersection* pBoolNode34 =
        new TGeoIntersection(pGEnceinte_externe_187_54, pGEnceinte_externe_73_29, 0, pMatrix29);
    // Shape: GEnceinte_externe_188 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_188_60 = new TGeoCompositeShape("GEnceinte_externe_188", pBoolNode34);
    // Combi transformation:
    dx = 0.000000;
    dy = -298.164800;
    dz = 115.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix41 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix40 = new TGeoCombiTrans("", dx, dy, dz, pMatrix41);
    TGeoUnion* pBoolNode35 = new TGeoUnion(pGEnceinte_externe_186_63, pGEnceinte_externe_188_60, 0, pMatrix40);
    // Shape: GEnceinte_externe_194 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_194_64 = new TGeoCompositeShape("GEnceinte_externe_194", pBoolNode35);
    // Shape: GEnceinte_externe_193 type: TGeoBBox
    dx = 86.200000;
    dy = 121.500000;
    dz = 108.805970;
    TGeoShape* pGEnceinte_externe_193_47 = new TGeoBBox("GEnceinte_externe_193", dx, dy, dz);
    // Shape: GEnceinte_externe_192 type: TGeoBBox
    dx = 84.900000;
    dy = 120.000000;
    dz = 122.000000;
    TGeoShape* pGEnceinte_externe_192_48 = new TGeoBBox("GEnceinte_externe_192", dx, dy, dz);
    TGeoSubtraction* pBoolNode36 = new TGeoSubtraction(pGEnceinte_externe_193_47, pGEnceinte_externe_192_48, 0, 0);
    // Shape: GEnceinte_externe_191 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_191_55 = new TGeoCompositeShape("GEnceinte_externe_191", pBoolNode36);
    // Shape: GEnceinte_externe_190 type: TGeoTubeSeg
    rmin = 40.000000;
    rmax = 122.500000;
    dz = 45.000000;
    phi1 = 26.424000;
    phi2 = 90.000000;
    TGeoShape* pGEnceinte_externe_190_30 = new TGeoTubeSeg("GEnceinte_externe_190", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 45.000000;
    dy = 0.003800;
    dz = 0.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix33 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix32 = new TGeoCombiTrans("", dx, dy, dz, pMatrix33);
    TGeoIntersection* pBoolNode37 =
        new TGeoIntersection(pGEnceinte_externe_191_55, pGEnceinte_externe_190_30, 0, pMatrix32);
    // Shape: GEnceinte_externe_189 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_189_61 = new TGeoCompositeShape("GEnceinte_externe_189", pBoolNode37);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 115.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix44 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix43 = new TGeoCombiTrans("", dx, dy, dz, pMatrix44);
    TGeoUnion* pBoolNode38 = new TGeoUnion(pGEnceinte_externe_194_64, pGEnceinte_externe_189_61, 0, pMatrix43);
    // Shape: GEnceinte_externe_195 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_195_65 = new TGeoCompositeShape("GEnceinte_externe_195", pBoolNode38);
    // Shape: GEnceinte_externe_48 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 300.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_48_31 = new TGeoTubeSeg("GEnceinte_externe_48", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -33.000000;
    dy = -150.000000;
    dz = -37.000000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 0.000000;
    phy = 0.000000;
    thz = 90.000000;
    phz = 270.000000;
    TGeoRotation* pMatrix51 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix50 = new TGeoCombiTrans("", dx, dy, dz, pMatrix51);
    TGeoSubtraction* pBoolNode39 =
        new TGeoSubtraction(pGEnceinte_externe_195_65, pGEnceinte_externe_48_31, 0, pMatrix50);
    // Shape: GEnceinte_externe_200 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_200_68 = new TGeoCompositeShape("GEnceinte_externe_200", pBoolNode39);
    // Shape: GEnceinte_externe_50 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 300.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_50_32 = new TGeoTubeSeg("GEnceinte_externe_50", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -33.000000;
    dy = -150.000000;
    dz = 38.000000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 0.000000;
    phy = 0.000000;
    thz = 90.000000;
    phz = 270.000000;
    TGeoRotation* pMatrix54 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix53 = new TGeoCombiTrans("", dx, dy, dz, pMatrix54);
    TGeoSubtraction* pBoolNode40 =
        new TGeoSubtraction(pGEnceinte_externe_200_68, pGEnceinte_externe_50_32, 0, pMatrix53);
    // Shape: GEnceinte_externe_199 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_199_69 = new TGeoCompositeShape("GEnceinte_externe_199", pBoolNode40);
    // Shape: GEnceinte_externe_52 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 10.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_52_33 = new TGeoTubeSeg("GEnceinte_externe_52", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -130.000000;
    dy = 11.919500;
    dz = -37.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix57 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix56 = new TGeoCombiTrans("", dx, dy, dz, pMatrix57);
    TGeoSubtraction* pBoolNode41 =
        new TGeoSubtraction(pGEnceinte_externe_199_69, pGEnceinte_externe_52_33, 0, pMatrix56);
    // Shape: GEnceinte_externe_198 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_198_70 = new TGeoCompositeShape("GEnceinte_externe_198", pBoolNode41);
    // Shape: GEnceinte_externe_54 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 10.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_54_34 = new TGeoTubeSeg("GEnceinte_externe_54", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -130.000000;
    dy = 11.919500;
    dz = 38.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix60 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix59 = new TGeoCombiTrans("", dx, dy, dz, pMatrix60);
    TGeoSubtraction* pBoolNode42 =
        new TGeoSubtraction(pGEnceinte_externe_198_70, pGEnceinte_externe_54_34, 0, pMatrix59);
    // Shape: GEnceinte_externe_197 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_197_71 = new TGeoCompositeShape("GEnceinte_externe_197", pBoolNode42);
    // Shape: GEnceinte_externe_56 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 10.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_56_35 = new TGeoTubeSeg("GEnceinte_externe_56", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -130.000000;
    dy = -310.080500;
    dz = -37.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix63 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix62 = new TGeoCombiTrans("", dx, dy, dz, pMatrix63);
    TGeoSubtraction* pBoolNode43 =
        new TGeoSubtraction(pGEnceinte_externe_197_71, pGEnceinte_externe_56_35, 0, pMatrix62);
    // Shape: GEnceinte_externe_196 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_196_72 = new TGeoCompositeShape("GEnceinte_externe_196", pBoolNode43);
    // Shape: GEnceinte_externe_58 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 10.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_58_36 = new TGeoTubeSeg("GEnceinte_externe_58", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -130.000000;
    dy = -310.080500;
    dz = 38.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix66 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix65 = new TGeoCombiTrans("", dx, dy, dz, pMatrix66);
    TGeoSubtraction* pBoolNode44 =
        new TGeoSubtraction(pGEnceinte_externe_196_72, pGEnceinte_externe_58_36, 0, pMatrix65);
    // Shape: GEnceinte_externe_203 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_203_73 = new TGeoCompositeShape("GEnceinte_externe_203", pBoolNode44);
    // Shape: GEnceinte_externe_60 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 19.980000;
    dz = 5.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_60_37 = new TGeoTubeSeg("GEnceinte_externe_60", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -155.000000;
    dy = -114.080500;
    dz = 55.500000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix69 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix68 = new TGeoCombiTrans("", dx, dy, dz, pMatrix69);
    TGeoSubtraction* pBoolNode45 =
        new TGeoSubtraction(pGEnceinte_externe_203_73, pGEnceinte_externe_60_37, 0, pMatrix68);
    // Shape: GEnceinte_externe_202 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_202_74 = new TGeoCompositeShape("GEnceinte_externe_202", pBoolNode45);
    // Shape: GEnceinte_externe_62 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 25.000000;
    dz = 5.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGEnceinte_externe_62_38 = new TGeoTubeSeg("GEnceinte_externe_62", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -155.000000;
    dy = -184.080500;
    dz = 55.500000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix72 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix71 = new TGeoCombiTrans("", dx, dy, dz, pMatrix72);
    TGeoSubtraction* pBoolNode46 =
        new TGeoSubtraction(pGEnceinte_externe_202_74, pGEnceinte_externe_62_38, 0, pMatrix71);
    // Shape: GEnceinte_externe_201 type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_201_75 = new TGeoCompositeShape("GEnceinte_externe_201", pBoolNode46);
    TGeoUnion* pBoolNode48 = new TGeoUnion(pGEnceinte_externe_181_94, pGEnceinte_externe_201_75, 0, 0);
    // Shape: GEnceinte_externe_PartBody type: TGeoCompositeShape
    TGeoShape* pGEnceinte_externe_PartBody_95 = new TGeoCompositeShape("GEnceinte_externe_PartBody", pBoolNode48);
    // Volume: GEnceinte_externe
    TGeoVolume* pGEnceinte_externe = new TGeoVolume("GEnceinte_externe", pGEnceinte_externe_PartBody_95, pMedAl);
    pGEnceinte_externe->SetLineColor(920);
    pGEnceinte_externe->SetVisLeaves(kTRUE);

    tr[0] = 0.000000;
    tr[1] = 150.000000;
    tr[2] = -9.000000;
    rot[0] = 1.000000;
    rot[1] = 0.000000;
    rot[2] = 0.000000;
    rot[3] = 0.000000;
    rot[4] = 1.000000;
    rot[5] = 0.000000;
    rot[6] = 0.000000;
    rot[7] = 0.000000;
    rot[8] = 1.000000;
    TGeoHMatrix* pMatrix128 = new TGeoHMatrix("");
    pMatrix128->SetTranslation(tr);
    pMatrix128->SetRotation(rot);
    pMatrix128->SetBit(TGeoMatrix::kGeoTranslation);
    pMatrix128->SetBit(TGeoMatrix::kGeoRotation);
    pMatrix128->RotateY(+90.0);
    pWorld->AddNode(pGEnceinte_externe, 1, pMatrix128);
}

void ConstructFonfE(TGeoVolume* pWorld)
{
    // Enceinte Extern part
    Double_t dx, dy, dz;
    Double_t dx1, dy1, dz1;
    Double_t dx2, dy2, dz2;
    Double_t rmin, rmax, phi1, phi2;
    Double_t thx, thy, thz, phx, phy, phz;

    // -- Mixture: Air
    TGeoMedium* pMed2 = gGeoManager->GetMedium("Air");
    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");
    TGeoMedium* pMedCu = gGeoManager->GetMedium("copper");
    Double_t tr[3], rot[9];

    // Shape: G2403002_Fonf_cote_sortie_42 type: TGeoBBox
    dx = 2.500000;
    dy = 240.500000;
    dz = 73.500000;
    TGeoShape* pG2403002_Fonf_cote_sortie_42_2 = new TGeoBBox("G2403002_Fonf_cote_sortie_42", dx, dy, dz);
    // Shape: G2403002_Fonf_cote_sortie_41 type: TGeoBBox
    dx = 6.250000;
    dy = 223.500000;
    dz = 56.500000;
    TGeoShape* pG2403002_Fonf_cote_sortie_41_3 = new TGeoBBox("G2403002_Fonf_cote_sortie_41", dx, dy, dz);
    // Combi transformation:
    dx = 8.750000;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix3 = new TGeoCombiTrans("");
    pMatrix3->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode1 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_42_2, pG2403002_Fonf_cote_sortie_41_3, 0, pMatrix3);
    // Shape: G2403002_Fonf_cote_sortie_40 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_40_40 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_40", pBoolNode1);
    // Shape: G2403002_Fonf_cote_sortie_39 type: TGeoBBox
    dx = 10.000000;
    dy = 220.500000;
    dz = 53.500000;
    TGeoShape* pG2403002_Fonf_cote_sortie_39_4 = new TGeoBBox("G2403002_Fonf_cote_sortie_39", dx, dy, dz);
    // Combi transformation:
    dx = 6.500000;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix21 = new TGeoCombiTrans("");
    pMatrix21->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode2 =
        new TGeoSubtraction(pG2403002_Fonf_cote_sortie_40_40, pG2403002_Fonf_cote_sortie_39_4, 0, pMatrix21);
    // Shape: G2403002_Fonf_cote_sortie_2 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_2_49 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_2", pBoolNode2);
    // Shape: G2403002_Fonf_cote_sortie_36 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 80.000000;
    dz = 0.750000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_36_6 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_36", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 15.750000;
    dy = -149.080515;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix32 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix31 = new TGeoCombiTrans("", dx, dy, dz, pMatrix32);
    TGeoUnion* pBoolNode3 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_2_49, pG2403002_Fonf_cote_sortie_36_6, 0, pMatrix31);
    // Shape: G2403002_Fonf_cote_sortie_82 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_82_54 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_82", pBoolNode3);
    // Shape: G2403002_Fonf_cote_sortie_35 type: TGeoTubeSeg
    rmin = 350.000000;
    rmax = 415.000000;
    dz = 0.750000;
    phi1 = 333.575610;
    phi2 = 386.424390;
    TGeoShape* pG2403002_Fonf_cote_sortie_35_7 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_35", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 15.750000;
    dy = 0.000000;
    dz = -300.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix35 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix34 = new TGeoCombiTrans("", dx, dy, dz, pMatrix35);
    TGeoUnion* pBoolNode4 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_82_54, pG2403002_Fonf_cote_sortie_35_7, 0, pMatrix34);
    // Shape: G2403002_Fonf_cote_sortie_81 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_81_55 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_81", pBoolNode4);
    // Shape: G2403002_Fonf_cote_sortie_34 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 80.000000;
    dz = 0.750000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_34_8 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_34", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 15.750000;
    dy = 149.080515;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix38 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix37 = new TGeoCombiTrans("", dx, dy, dz, pMatrix38);
    TGeoUnion* pBoolNode5 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_81_55, pG2403002_Fonf_cote_sortie_34_8, 0, pMatrix37);
    // Shape: G2403002_Fonf_cote_sortie_80 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_80_56 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_80", pBoolNode5);
    // Shape: G2403002_Fonf_cote_sortie_32 type: TGeoTubeSeg
    rmin = 350.000000;
    rmax = 415.000000;
    dz = 0.750000;
    phi1 = 333.575610;
    phi2 = 386.424390;
    TGeoShape* pG2403002_Fonf_cote_sortie_32_9 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_32", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 15.750000;
    dy = 0.000000;
    dz = 300.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix41 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix40 = new TGeoCombiTrans("", dx, dy, dz, pMatrix41);
    TGeoUnion* pBoolNode6 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_80_56, pG2403002_Fonf_cote_sortie_32_9, 0, pMatrix40);
    // Shape: G2403002_Fonf_cote_sortie_79 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_79_57 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_79", pBoolNode6);
    // Shape: G2403002_Fonf_cote_sortie_29 type: TGeoTorus
    Double_t r = 80.000000;
    rmin = 40.000000;
    rmax = 41.500000;
    phi1 = 26.424000;
    Double_t dphi = 127.152000;
    TGeoShape* pG2403002_Fonf_cote_sortie_29_36 =
        new TGeoTorus("G2403002_Fonf_cote_sortie_29", r, rmin, rmax, phi1, dphi);
    // Shape: G2403002_Fonf_cote_sortie_28 type: TGeoTubeSeg
    rmin = 80.000000;
    rmax = 122.000000;
    dz = 42.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_28_10 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_28", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -42.000000;
    TGeoCombiTrans* pMatrix23 = new TGeoCombiTrans("");
    pMatrix23->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode7 =
        new TGeoIntersection(pG2403002_Fonf_cote_sortie_29_36, pG2403002_Fonf_cote_sortie_28_10, 0, pMatrix23);
    // Shape: G2403002_Fonf_cote_sortie_78 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_78_50 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_78", pBoolNode7);
    // Combi transformation:
    dx = 56.500000;
    dy = 149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix44 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix43 = new TGeoCombiTrans("", dx, dy, dz, pMatrix44);
    TGeoUnion* pBoolNode8 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_79_57, pG2403002_Fonf_cote_sortie_78_50, 0, pMatrix43);
    // Shape: G2403002_Fonf_cote_sortie_92 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_92_58 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_92", pBoolNode8);
    // Shape: G2403002_Fonf_cote_sortie_22 type: TGeoTorus
    r = 415.000000;
    rmin = 40.000000;
    rmax = 41.500000;
    phi1 = 333.575610;
    dphi = 52.848780;
    TGeoShape* pG2403002_Fonf_cote_sortie_22_37 =
        new TGeoTorus("G2403002_Fonf_cote_sortie_22", r, rmin, rmax, phi1, dphi);
    // Shape: G2403002_Fonf_cote_sortie_21 type: TGeoTubeSeg
    rmin = 415.000000;
    rmax = 457.000000;
    dz = 42.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_21_11 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_21", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -42.000000;
    TGeoCombiTrans* pMatrix25 = new TGeoCombiTrans("");
    pMatrix25->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode9 =
        new TGeoIntersection(pG2403002_Fonf_cote_sortie_22_37, pG2403002_Fonf_cote_sortie_21_11, 0, pMatrix25);
    // Shape: G2403002_Fonf_cote_sortie_77 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_77_51 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_77", pBoolNode9);
    // Combi transformation:
    dx = 56.500000;
    dy = 0.000000;
    dz = 300.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix47 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix46 = new TGeoCombiTrans("", dx, dy, dz, pMatrix47);
    TGeoUnion* pBoolNode10 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_92_58, pG2403002_Fonf_cote_sortie_77_51, 0, pMatrix46);
    // Shape: G2403002_Fonf_cote_sortie_91 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_91_59 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_91", pBoolNode10);
    // Shape: G2403002_Fonf_cote_sortie_88 type: TGeoTorus
    r = 415.000000;
    rmin = 40.000000;
    rmax = 41.500000;
    phi1 = 153.576000;
    dphi = 52.848000;
    TGeoShape* pG2403002_Fonf_cote_sortie_88_38 =
        new TGeoTorus("G2403002_Fonf_cote_sortie_88", r, rmin, rmax, phi1, dphi);
    // Shape: G2403002_Fonf_cote_sortie_87 type: TGeoTubeSeg
    rmin = 415.000000;
    rmax = 457.000000;
    dz = 42.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_87_12 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_87", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -42.000000;
    TGeoCombiTrans* pMatrix27 = new TGeoCombiTrans("");
    pMatrix27->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode11 =
        new TGeoIntersection(pG2403002_Fonf_cote_sortie_88_38, pG2403002_Fonf_cote_sortie_87_12, 0, pMatrix27);
    // Shape: G2403002_Fonf_cote_sortie_86 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_86_52 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_86", pBoolNode11);
    // Combi transformation:
    dx = 56.500000;
    dy = 0.000000;
    dz = -300.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix50 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix49 = new TGeoCombiTrans("", dx, dy, dz, pMatrix50);
    TGeoUnion* pBoolNode12 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_91_59, pG2403002_Fonf_cote_sortie_86_52, 0, pMatrix49);
    // Shape: G2403002_Fonf_cote_sortie_90 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_90_60 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_90", pBoolNode12);
    // Shape: G2403002_Fonf_cote_sortie_85 type: TGeoTorus
    r = 80.000000;
    rmin = 40.000000;
    rmax = 41.500000;
    phi1 = 206.424000;
    dphi = 127.152000;
    TGeoShape* pG2403002_Fonf_cote_sortie_85_39 =
        new TGeoTorus("G2403002_Fonf_cote_sortie_85", r, rmin, rmax, phi1, dphi);
    // Shape: G2403002_Fonf_cote_sortie_84 type: TGeoTubeSeg
    rmin = 80.000000;
    rmax = 122.000000;
    dz = 42.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_84_13 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_84", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -42.000000;
    TGeoCombiTrans* pMatrix29 = new TGeoCombiTrans("");
    pMatrix29->SetTranslation(dx, dy, dz);
    TGeoIntersection* pBoolNode13 =
        new TGeoIntersection(pG2403002_Fonf_cote_sortie_85_39, pG2403002_Fonf_cote_sortie_84_13, 0, pMatrix29);
    // Shape: G2403002_Fonf_cote_sortie_83 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_83_53 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_83", pBoolNode13);
    // Combi transformation:
    dx = 56.500000;
    dy = -149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix53 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix52 = new TGeoCombiTrans("", dx, dy, dz, pMatrix53);
    TGeoUnion* pBoolNode14 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_90_60, pG2403002_Fonf_cote_sortie_83_53, 0, pMatrix52);
    // Shape: G2403002_Fonf_cote_sortie_89 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_89_61 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_89", pBoolNode14);
    // Shape: G2403002_Fonf_cote_sortie_13 type: TGeoBBox
    dx = 20.000000;
    dy = 220.500000;
    dz = 53.500000;
    TGeoShape* pG2403002_Fonf_cote_sortie_13_5 = new TGeoBBox("G2403002_Fonf_cote_sortie_13", dx, dy, dz);
    TGeoSubtraction* pBoolNode15 =
        new TGeoSubtraction(pG2403002_Fonf_cote_sortie_89_61, pG2403002_Fonf_cote_sortie_13_5, 0, 0);
    // Shape: G2403002_Fonf_cote_sortie_93 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_93_62 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_93", pBoolNode15);
    // Shape: G2403002_Fonf_cote_sortie_11 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 20.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_11_14 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_11", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -55.000000;
    dz = 100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix58 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix57 = new TGeoCombiTrans("", dx, dy, dz, pMatrix58);
    TGeoSubtraction* pBoolNode16 =
        new TGeoSubtraction(pG2403002_Fonf_cote_sortie_93_62, pG2403002_Fonf_cote_sortie_11_14, 0, pMatrix57);
    // Shape: G2403002_Fonf_cote_sortie_98 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_98_63 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_98", pBoolNode16);
    // Shape: G2403002_Fonf_cote_sortie_9 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 20.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_9_15 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_9", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 55.000000;
    dz = 100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix61 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix60 = new TGeoCombiTrans("", dx, dy, dz, pMatrix61);
    TGeoSubtraction* pBoolNode17 =
        new TGeoSubtraction(pG2403002_Fonf_cote_sortie_98_63, pG2403002_Fonf_cote_sortie_9_15, 0, pMatrix60);
    // Shape: G2403002_Fonf_cote_sortie_97 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_97_64 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_97", pBoolNode17);
    // Shape: G2403002_Fonf_cote_sortie_7 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 20.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_7_16 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_7", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = -55.000000;
    dz = -100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix64 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix63 = new TGeoCombiTrans("", dx, dy, dz, pMatrix64);
    TGeoSubtraction* pBoolNode18 =
        new TGeoSubtraction(pG2403002_Fonf_cote_sortie_97_64, pG2403002_Fonf_cote_sortie_7_16, 0, pMatrix63);
    // Shape: G2403002_Fonf_cote_sortie_96 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_96_65 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_96", pBoolNode18);
    // Shape: G2403002_Fonf_cote_sortie_5 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.600000;
    dz = 20.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_5_17 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_5", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 55.000000;
    dz = -100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix67 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix66 = new TGeoCombiTrans("", dx, dy, dz, pMatrix67);
    TGeoSubtraction* pBoolNode19 =
        new TGeoSubtraction(pG2403002_Fonf_cote_sortie_96_65, pG2403002_Fonf_cote_sortie_5_17, 0, pMatrix66);
    // Shape: G2403002_Fonf_cote_sortie_95 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_95_66 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_95", pBoolNode19);
    // Shape: G2403002_Fonf_cote_sortie_3 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 25.000000;
    dz = 25.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_3_18 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_3", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -109.400000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix70 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix69 = new TGeoCombiTrans("", dx, dy, dz, pMatrix70);
    TGeoSubtraction* pBoolNode20 =
        new TGeoSubtraction(pG2403002_Fonf_cote_sortie_95_66, pG2403002_Fonf_cote_sortie_3_18, 0, pMatrix69);
    // Shape: G2403002_Fonf_cote_sortie_94 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_94_67 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_94", pBoolNode20);
    // Shape: G2403002_Fonf_cote_sortie_113 type: TGeoTubeSeg
    rmin = 11.400000;
    rmax = 15.500000;
    dz = 3.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_113_19 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_113", rmin, rmax, dz, phi1, phi2);
    // Shape: G2403002_Fonf_cote_sortie_112 type: TGeoTubeSeg
    rmin = 15.500000;
    rmax = 18.100000;
    dz = 11.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_112_20 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_112", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.250000;
    TGeoCombiTrans* pMatrix5 = new TGeoCombiTrans("");
    pMatrix5->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode21 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_113_19, pG2403002_Fonf_cote_sortie_112_20, 0, pMatrix5);
    // Shape: G2403002_Fonf_cote_sortie_111 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_111_41 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_111", pBoolNode21);
    // Shape: G2403002_Fonf_cote_sortie_110 type: TGeoTubeSeg
    rmin = 18.100000;
    rmax = 19.000000;
    dz = 12.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_110_21 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_110", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.750000;
    TGeoCombiTrans* pMatrix7 = new TGeoCombiTrans("");
    pMatrix7->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode22 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_111_41, pG2403002_Fonf_cote_sortie_110_21, 0, pMatrix7);
    // Shape: G2403002_Fonf_cote_sortie_200 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_200_42 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_200", pBoolNode22);
    // Combi transformation:
    dx = 11.750000;
    dy = 55.000000;
    dz = 100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix73 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix72 = new TGeoCombiTrans("", dx, dy, dz, pMatrix73);
    TGeoUnion* pBoolNode23 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_94_67, pG2403002_Fonf_cote_sortie_200_42, 0, pMatrix72);
    // Shape: G2403002_Fonf_cote_sortie_117 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_117_68 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_117", pBoolNode23);
    // Shape: G2403002_Fonf_cote_sortie_109 type: TGeoTubeSeg
    rmin = 11.400000;
    rmax = 15.500000;
    dz = 3.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_109_22 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_109", rmin, rmax, dz, phi1, phi2);
    // Shape: G2403002_Fonf_cote_sortie_108 type: TGeoTubeSeg
    rmin = 15.500000;
    rmax = 18.100000;
    dz = 11.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_108_23 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_108", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.250000;
    TGeoCombiTrans* pMatrix9 = new TGeoCombiTrans("");
    pMatrix9->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode24 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_109_22, pG2403002_Fonf_cote_sortie_108_23, 0, pMatrix9);
    // Shape: G2403002_Fonf_cote_sortie_107 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_107_43 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_107", pBoolNode24);
    // Shape: G2403002_Fonf_cote_sortie_106 type: TGeoTubeSeg
    rmin = 18.100000;
    rmax = 19.000000;
    dz = 12.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_106_24 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_106", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.750000;
    TGeoCombiTrans* pMatrix11 = new TGeoCombiTrans("");
    pMatrix11->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode25 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_107_43, pG2403002_Fonf_cote_sortie_106_24, 0, pMatrix11);
    // Shape: G2403002_Fonf_cote_sortie_201 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_201_44 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_201", pBoolNode25);
    // Combi transformation:
    dx = 11.750000;
    dy = -55.000000;
    dz = 100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix76 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix75 = new TGeoCombiTrans("", dx, dy, dz, pMatrix76);
    TGeoUnion* pBoolNode26 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_117_68, pG2403002_Fonf_cote_sortie_201_44, 0, pMatrix75);
    // Shape: G2403002_Fonf_cote_sortie_116 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_116_69 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_116", pBoolNode26);
    // Shape: G2403002_Fonf_cote_sortie_103 type: TGeoTubeSeg
    rmin = 11.400000;
    rmax = 15.500000;
    dz = 3.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_103_25 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_103", rmin, rmax, dz, phi1, phi2);
    // Shape: G2403002_Fonf_cote_sortie_102 type: TGeoTubeSeg
    rmin = 15.500000;
    rmax = 18.100000;
    dz = 11.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_102_26 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_102", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.250000;
    TGeoCombiTrans* pMatrix13 = new TGeoCombiTrans("");
    pMatrix13->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode27 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_103_25, pG2403002_Fonf_cote_sortie_102_26, 0, pMatrix13);
    // Shape: G2403002_Fonf_cote_sortie_101 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_101_45 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_101", pBoolNode27);
    // Shape: G2403002_Fonf_cote_sortie_100 type: TGeoTubeSeg
    rmin = 18.100000;
    rmax = 19.000000;
    dz = 12.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_100_27 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_100", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.750000;
    TGeoCombiTrans* pMatrix15 = new TGeoCombiTrans("");
    pMatrix15->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode28 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_101_45, pG2403002_Fonf_cote_sortie_100_27, 0, pMatrix15);
    // Shape: G2403002_Fonf_cote_sortie_202 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_202_46 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_202", pBoolNode28);
    // Combi transformation:
    dx = 11.750000;
    dy = -55.000000;
    dz = -100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix79 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix78 = new TGeoCombiTrans("", dx, dy, dz, pMatrix79);
    TGeoUnion* pBoolNode29 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_116_69, pG2403002_Fonf_cote_sortie_202_46, 0, pMatrix78);
    // Shape: G2403002_Fonf_cote_sortie_115 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_115_70 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_115", pBoolNode29);
    // Shape: G2403002_Fonf_cote_sortie_64 type: TGeoTubeSeg
    rmin = 11.400000;
    rmax = 15.500000;
    dz = 3.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_64_28 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_64", rmin, rmax, dz, phi1, phi2);
    // Shape: G2403002_Fonf_cote_sortie_63 type: TGeoTubeSeg
    rmin = 15.500000;
    rmax = 18.100000;
    dz = 11.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_63_29 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_63", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.250000;
    TGeoCombiTrans* pMatrix17 = new TGeoCombiTrans("");
    pMatrix17->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode30 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_64_28, pG2403002_Fonf_cote_sortie_63_29, 0, pMatrix17);
    // Shape: G2403002_Fonf_cote_sortie_62 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_62_47 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_62", pBoolNode30);
    // Shape: G2403002_Fonf_cote_sortie_61 type: TGeoTubeSeg
    rmin = 18.100000;
    rmax = 19.000000;
    dz = 12.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_61_30 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_61", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -8.750000;
    TGeoCombiTrans* pMatrix19 = new TGeoCombiTrans("");
    pMatrix19->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode31 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_62_47, pG2403002_Fonf_cote_sortie_61_30, 0, pMatrix19);
    // Shape: G2403002_Fonf_cote_sortie_203 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_203_48 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_203", pBoolNode31);
    // Combi transformation:
    dx = 11.750000;
    dy = 55.000000;
    dz = -100.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix82 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix81 = new TGeoCombiTrans("", dx, dy, dz, pMatrix82);
    TGeoUnion* pBoolNode32 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_115_70, pG2403002_Fonf_cote_sortie_203_48, 0, pMatrix81);
    // Shape: G2403002_Fonf_cote_sortie_114 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_114_71 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_114", pBoolNode32);
    // Shape: G2403002_Fonf_cote_sortie_45 type: TGeoTubeSeg
    rmin = 25.000000;
    rmax = 25.400000;
    dz = 1.200000;
    phi1 = 323.000000;
    phi2 = 397.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_45_31 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_45", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 18.800000;
    dy = 0.000000;
    dz = -109.450000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix85 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix84 = new TGeoCombiTrans("", dx, dy, dz, pMatrix85);
    TGeoUnion* pBoolNode33 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_114_71, pG2403002_Fonf_cote_sortie_45_31, 0, pMatrix84);
    // Shape: G2403002_Fonf_cote_sortie_119 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_119_72 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_119", pBoolNode33);
    // Shape: G2403002_Fonf_cote_sortie_49 type: TGeoTubeSeg
    rmin = 25.000000;
    rmax = 25.400000;
    dz = 7.200000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_49_32 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_49", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 9.200000;
    dy = 0.000000;
    dz = -109.450000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix88 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix87 = new TGeoCombiTrans("", dx, dy, dz, pMatrix88);
    TGeoUnion* pBoolNode34 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_119_72, pG2403002_Fonf_cote_sortie_49_32, 0, pMatrix87);
    // Shape: G2403002_Fonf_cote_sortie_118 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_118_73 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_118", pBoolNode34);
    // Shape: G2403002_Fonf_cote_sortie_47 type: TGeoTubeSeg
    rmin = 25.000000;
    rmax = 25.400000;
    dz = 0.600000;
    phi1 = 307.000000;
    phi2 = 413.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_47_33 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_47", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 17.000000;
    dy = 0.000000;
    dz = -109.450000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix91 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix90 = new TGeoCombiTrans("", dx, dy, dz, pMatrix91);
    TGeoUnion* pBoolNode35 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_118_73, pG2403002_Fonf_cote_sortie_47_33, 0, pMatrix90);
    // Shape: G2403002_Fonf_cote_sortie_121 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_121_74 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_121", pBoolNode35);
    // Shape: G2403002_Fonf_cote_sortie_48 type: TGeoTubeSeg
    rmin = 25.000000;
    rmax = 30.500000;
    dz = 1.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_48_34 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_48", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 1.000000;
    dy = 0.000000;
    dz = -109.450000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix94 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix93 = new TGeoCombiTrans("", dx, dy, dz, pMatrix94);
    TGeoUnion* pBoolNode36 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_121_74, pG2403002_Fonf_cote_sortie_48_34, 0, pMatrix93);
    // Shape: G2403002_Fonf_cote_sortie_120 type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_120_75 = new TGeoCompositeShape("G2403002_Fonf_cote_sortie_120", pBoolNode36);
    // Shape: G2403002_Fonf_cote_sortie_1 type: TGeoTubeSeg
    rmin = 21.000000;
    rmax = 25.000000;
    dz = 0.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403002_Fonf_cote_sortie_1_35 =
        new TGeoTubeSeg("G2403002_Fonf_cote_sortie_1", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 10.500000;
    dy = 0.000000;
    dz = -109.450000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix97 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix96 = new TGeoCombiTrans("", dx, dy, dz, pMatrix97);
    TGeoUnion* pBoolNode37 =
        new TGeoUnion(pG2403002_Fonf_cote_sortie_120_75, pG2403002_Fonf_cote_sortie_1_35, 0, pMatrix96);
    // Shape: G2403002_Fonf_cote_sortie_PartBody type: TGeoCompositeShape
    TGeoShape* pG2403002_Fonf_cote_sortie_PartBody_76 =
        new TGeoCompositeShape("G2403002_Fonf_cote_sortie_PartBody", pBoolNode37);
    // Volume: G2403002_Fonf_cote_sortie

    TGeoVolume* pG2403002_Fonf_cote_sortie =
        new TGeoVolume("G2403002_Fonf_cote_sortie", pG2403002_Fonf_cote_sortie_PartBody_76, pMedCu);
    pG2403002_Fonf_cote_sortie->SetLineColor(920);
    pG2403002_Fonf_cote_sortie->SetVisLeaves(kTRUE);

    tr[0] = -221.500000 + 41;
    tr[1] = 0.000000;
    tr[2] = 0.000000;
    rot[0] = 1.000000;
    rot[1] = 0.000000;
    rot[2] = 0.000000;
    rot[3] = 0.000000;
    rot[4] = 1.000000;
    rot[5] = 0.000000;
    rot[6] = 0.000000;
    rot[7] = 0.000000;
    rot[8] = 1.000000;
    TGeoHMatrix* pMatrix98 = new TGeoHMatrix("");

    // pMatrix98->SetRotation(rot);
    // pMatrix98->RotateZ(180.0);
    pMatrix98->SetTranslation(tr);
    pMatrix98->SetBit(TGeoMatrix::kGeoTranslation);
    pMatrix98->SetBit(TGeoMatrix::kGeoRotation);

    pWorld->AddNode(pG2403002_Fonf_cote_sortie, 1, pMatrix98);
}

void ConstructFondS(TGeoVolume* pWorld)
{
    // Enceinte Extern part
    Double_t dx, dy, dz;
    Double_t dx1, dy1, dz1;
    Double_t dx2, dy2, dz2;
    Double_t rmin, rmax, phi1, phi2;
    Double_t thx, thy, thz, phx, phy, phz;

    // -- Mixture: Air
    TGeoMedium* pMed2 = gGeoManager->GetMedium("Air");
    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");

    Double_t tr[3], rot[9];

    // Shape: G2403001_Fond_cote_entree_92 type: TGeoBBox
    dx = 3.500000;
    dy = 60.750000;
    dz = 21.500000;
    TGeoShape* pG2403001_Fond_cote_entree_92_1 = new TGeoBBox("G2403001_Fond_cote_entree_92", dx, dy, dz);
    // Shape: G2403001_Fond_cote_entree_96 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 365.000000;
    dz = 1.500000;
    phi1 = 333.575610;
    phi2 = 386.424390;
    TGeoShape* pG2403001_Fond_cote_entree_96_22 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_96", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -5.000000;
    dy = 0.000000;
    dz = -300.000000;
    // Rotation:
    thx = 0.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 180.000000;
    TGeoRotation* pMatrix4 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix3 = new TGeoCombiTrans("", dx, dy, dz, pMatrix4);
    TGeoUnion* pBoolNode1 =
        new TGeoUnion(pG2403001_Fond_cote_entree_92_1, pG2403001_Fond_cote_entree_96_22, 0, pMatrix3);
    // Shape: G2403001_Fond_cote_entree_100 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_100_36 = new TGeoCompositeShape("G2403001_Fond_cote_entree_100", pBoolNode1);
    // Shape: G2403001_Fond_cote_entree_95 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 30.000000;
    dz = 1.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_95_23 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_95", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -5.000000;
    dy = 149.080515;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix7 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix6 = new TGeoCombiTrans("", dx, dy, dz, pMatrix7);
    TGeoUnion* pBoolNode2 =
        new TGeoUnion(pG2403001_Fond_cote_entree_100_36, pG2403001_Fond_cote_entree_95_23, 0, pMatrix6);
    // Shape: G2403001_Fond_cote_entree_99 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_99_37 = new TGeoCompositeShape("G2403001_Fond_cote_entree_99", pBoolNode2);
    // Shape: G2403001_Fond_cote_entree_94 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 30.000000;
    dz = 1.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_94_24 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_94", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -5.000000;
    dy = -149.080515;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix10 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix9 = new TGeoCombiTrans("", dx, dy, dz, pMatrix10);
    TGeoUnion* pBoolNode3 =
        new TGeoUnion(pG2403001_Fond_cote_entree_99_37, pG2403001_Fond_cote_entree_94_24, 0, pMatrix9);
    // Shape: G2403001_Fond_cote_entree_98 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_98_38 = new TGeoCompositeShape("G2403001_Fond_cote_entree_98", pBoolNode3);
    // Shape: G2403001_Fond_cote_entree_93 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 365.000000;
    dz = 1.500000;
    phi1 = 333.575610;
    phi2 = 386.424390;
    TGeoShape* pG2403001_Fond_cote_entree_93_25 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_93", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -5.000000;
    dy = 0.000000;
    dz = 300.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix13 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix12 = new TGeoCombiTrans("", dx, dy, dz, pMatrix13);
    TGeoUnion* pBoolNode4 =
        new TGeoUnion(pG2403001_Fond_cote_entree_98_38, pG2403001_Fond_cote_entree_93_25, 0, pMatrix12);
    // Shape: G2403001_Fond_cote_entree_97 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_97_39 = new TGeoCompositeShape("G2403001_Fond_cote_entree_97", pBoolNode4);
    // Shape: G2403001_Fond_cote_entree_101 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 59.000000;
    dz = 3.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_101_26 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_101", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix16 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix15 = new TGeoCombiTrans("", dx, dy, dz, pMatrix16);
    TGeoUnion* pBoolNode5 =
        new TGeoUnion(pG2403001_Fond_cote_entree_97_39, pG2403001_Fond_cote_entree_101_26, 0, pMatrix15);
    // Shape: G2403001_Fond_cote_entree_102 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_102_40 = new TGeoCompositeShape("G2403001_Fond_cote_entree_102", pBoolNode5);
    // Shape: G2403001_Fond_cote_entree_103 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 47.500000;
    dz = 10.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_103_27 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_103", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 3.000000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix25 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix24 = new TGeoCombiTrans("", dx, dy, dz, pMatrix25);
    TGeoSubtraction* pBoolNode6 =
        new TGeoSubtraction(pG2403001_Fond_cote_entree_102_40, pG2403001_Fond_cote_entree_103_27, 0, pMatrix24);
    // Shape: G2403001_Fond_cote_entree_104 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_104_43 = new TGeoCompositeShape("G2403001_Fond_cote_entree_104", pBoolNode6);
    // Shape: G2403001_Fond_cote_entree_106 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 56.000000;
    dz = 5.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_106_28 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_106", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 5.500000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix28 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix27 = new TGeoCombiTrans("", dx, dy, dz, pMatrix28);
    TGeoSubtraction* pBoolNode7 =
        new TGeoSubtraction(pG2403001_Fond_cote_entree_104_43, pG2403001_Fond_cote_entree_106_28, 0, pMatrix27);
    // Shape: G2403001_Fond_cote_entree_107 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_107_44 = new TGeoCompositeShape("G2403001_Fond_cote_entree_107", pBoolNode7);
    // Shape: G2403001_Fond_cote_entree_108 type: TGeoTubeSeg
    rmin = 32.500000;
    rmax = 47.500000;
    dz = 1.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_108_29 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_108", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -5.000000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix37 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix36 = new TGeoCombiTrans("", dx, dy, dz, pMatrix37);
    TGeoUnion* pBoolNode8 =
        new TGeoUnion(pG2403001_Fond_cote_entree_107_44, pG2403001_Fond_cote_entree_108_29, 0, pMatrix36);
    // Shape: G2403001_Fond_cote_entree_109 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_109_47 = new TGeoCompositeShape("G2403001_Fond_cote_entree_109", pBoolNode8);
    // Shape: G2403001_Fond_cote_entree_110 type: TGeoBBox
    dx = 3.500000;
    dy = 19.000000;
    dz = 19.000000;
    TGeoShape* pG2403001_Fond_cote_entree_110_2 = new TGeoBBox("G2403001_Fond_cote_entree_110", dx, dy, dz);
    // Shape: G2403001_Fond_cote_entree_111 type: TGeoTubeSeg
    rmin = 11.400000;
    rmax = 18.000000;
    dz = 8.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_111_30 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_111", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -11.750000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix19 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix18 = new TGeoCombiTrans("", dx, dy, dz, pMatrix19);
    TGeoUnion* pBoolNode9 =
        new TGeoUnion(pG2403001_Fond_cote_entree_110_2, pG2403001_Fond_cote_entree_111_30, 0, pMatrix18);
    // Shape: G2403001_Fond_cote_entree_113 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_113_41 = new TGeoCompositeShape("G2403001_Fond_cote_entree_113", pBoolNode9);
    // Shape: G2403001_Fond_cote_entree_112 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 18.100000;
    dz = 2.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_112_31 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_112", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 5.000000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix31 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix30 = new TGeoCombiTrans("", dx, dy, dz, pMatrix31);
    TGeoSubtraction* pBoolNode10 =
        new TGeoSubtraction(pG2403001_Fond_cote_entree_113_41, pG2403001_Fond_cote_entree_112_31, 0, pMatrix30);
    // Shape: G2403001_Fond_cote_entree_114 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_114_45 = new TGeoCompositeShape("G2403001_Fond_cote_entree_114", pBoolNode10);
    // Combi transformation:
    dx = 0.000000;
    dy = -133.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix39 = new TGeoCombiTrans("");
    pMatrix39->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode11 =
        new TGeoUnion(pG2403001_Fond_cote_entree_109_47, pG2403001_Fond_cote_entree_114_45, 0, pMatrix39);
    // Shape: G2403001_Fond_cote_entree_115 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_115_48 = new TGeoCompositeShape("G2403001_Fond_cote_entree_115", pBoolNode11);
    // Shape: G2403001_Fond_cote_entree_116 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 15.500000;
    dz = 15.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_116_32 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_116", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.400000;
    dy = -133.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix42 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix41 = new TGeoCombiTrans("", dx, dy, dz, pMatrix42);
    TGeoSubtraction* pBoolNode12 =
        new TGeoSubtraction(pG2403001_Fond_cote_entree_115_48, pG2403001_Fond_cote_entree_116_32, 0, pMatrix41);
    // Shape: G2403001_Fond_cote_entree_117 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_117_49 = new TGeoCompositeShape("G2403001_Fond_cote_entree_117", pBoolNode12);
    // Shape: G2403001_Fond_cote_entree_122 type: TGeoBBox
    dx = 3.500000;
    dy = 19.000000;
    dz = 19.000000;
    TGeoShape* pG2403001_Fond_cote_entree_122_3 = new TGeoBBox("G2403001_Fond_cote_entree_122", dx, dy, dz);
    // Shape: G2403001_Fond_cote_entree_121 type: TGeoTubeSeg
    rmin = 11.400000;
    rmax = 18.000000;
    dz = 8.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_121_33 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_121", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -11.750000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix22 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix21 = new TGeoCombiTrans("", dx, dy, dz, pMatrix22);
    TGeoUnion* pBoolNode13 =
        new TGeoUnion(pG2403001_Fond_cote_entree_122_3, pG2403001_Fond_cote_entree_121_33, 0, pMatrix21);
    // Shape: G2403001_Fond_cote_entree_120 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_120_42 = new TGeoCompositeShape("G2403001_Fond_cote_entree_120", pBoolNode13);
    // Shape: G2403001_Fond_cote_entree_119 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 18.100000;
    dz = 2.500000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_119_34 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_119", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 5.000000;
    dy = 0.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix34 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix33 = new TGeoCombiTrans("", dx, dy, dz, pMatrix34);
    TGeoSubtraction* pBoolNode14 =
        new TGeoSubtraction(pG2403001_Fond_cote_entree_120_42, pG2403001_Fond_cote_entree_119_34, 0, pMatrix33);
    // Shape: G2403001_Fond_cote_entree_118 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_118_46 = new TGeoCompositeShape("G2403001_Fond_cote_entree_118", pBoolNode14);
    // Combi transformation:
    dx = 0.000000;
    dy = 133.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix44 = new TGeoCombiTrans("");
    pMatrix44->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode =
        new TGeoUnion(pG2403001_Fond_cote_entree_117_49, pG2403001_Fond_cote_entree_118_46, 0, pMatrix44);
    // Shape: G2403001_Fond_cote_entree_124 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_124_50 = new TGeoCompositeShape("G2403001_Fond_cote_entree_124", pBoolNode);
    // Shape: G2403001_Fond_cote_entree_123 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 15.500000;
    dz = 15.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pG2403001_Fond_cote_entree_123_35 =
        new TGeoTubeSeg("G2403001_Fond_cote_entree_123", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.400000;
    dy = 133.000000;
    dz = 0.000000;
    // Rotation:
    thx = 180.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 90.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix47 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix46 = new TGeoCombiTrans("", dx, dy, dz, pMatrix47);
    TGeoSubtraction* pBoolNode15 =
        new TGeoSubtraction(pG2403001_Fond_cote_entree_124_50, pG2403001_Fond_cote_entree_123_35, 0, pMatrix46);
    // Shape: G2403001_Fond_cote_entree_125 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_125_51 = new TGeoCompositeShape("G2403001_Fond_cote_entree_125", pBoolNode15);
    // Shape: G2403001_Fond_cote_entree_41 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 11.750000;
    TGeoShape* pG2403001_Fond_cote_entree_41_4 = new TGeoBBox("G2403001_Fond_cote_entree_41", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 120.000000;
    dz = 30.750000;
    TGeoCombiTrans* pMatrix49 = new TGeoCombiTrans("");
    pMatrix49->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode16 =
        new TGeoUnion(pG2403001_Fond_cote_entree_125_51, pG2403001_Fond_cote_entree_41_4, 0, pMatrix49);
    // Shape: G2403001_Fond_cote_entree_144 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_144_52 = new TGeoCompositeShape("G2403001_Fond_cote_entree_144", pBoolNode16);
    // Shape: G2403001_Fond_cote_entree_43 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 11.750000;
    TGeoShape* pG2403001_Fond_cote_entree_43_5 = new TGeoBBox("G2403001_Fond_cote_entree_43", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 120.000000;
    dz = -30.750000;
    TGeoCombiTrans* pMatrix51 = new TGeoCombiTrans("");
    pMatrix51->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode17 =
        new TGeoUnion(pG2403001_Fond_cote_entree_144_52, pG2403001_Fond_cote_entree_43_5, 0, pMatrix51);
    // Shape: G2403001_Fond_cote_entree_145 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_145_53 = new TGeoCompositeShape("G2403001_Fond_cote_entree_145", pBoolNode17);
    // Shape: G2403001_Fond_cote_entree_45 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 7.750000;
    TGeoShape* pG2403001_Fond_cote_entree_45_6 = new TGeoBBox("G2403001_Fond_cote_entree_45", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 140.000000;
    dz = 26.750000;
    TGeoCombiTrans* pMatrix53 = new TGeoCombiTrans("");
    pMatrix53->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode18 =
        new TGeoUnion(pG2403001_Fond_cote_entree_145_53, pG2403001_Fond_cote_entree_45_6, 0, pMatrix53);
    // Shape: G2403001_Fond_cote_entree_146 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_146_54 = new TGeoCompositeShape("G2403001_Fond_cote_entree_146", pBoolNode18);
    // Shape: G2403001_Fond_cote_entree_46 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 7.750000;
    TGeoShape* pG2403001_Fond_cote_entree_46_7 = new TGeoBBox("G2403001_Fond_cote_entree_46", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 140.000000;
    dz = -26.750000;
    TGeoCombiTrans* pMatrix55 = new TGeoCombiTrans("");
    pMatrix55->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode19 =
        new TGeoUnion(pG2403001_Fond_cote_entree_146_54, pG2403001_Fond_cote_entree_46_7, 0, pMatrix55);
    // Shape: G2403001_Fond_cote_entree_147 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_147_55 = new TGeoCompositeShape("G2403001_Fond_cote_entree_147", pBoolNode19);
    // Shape: G2403001_Fond_cote_entree_47 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 11.750000;
    TGeoShape* pG2403001_Fond_cote_entree_47_8 = new TGeoBBox("G2403001_Fond_cote_entree_47", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -120.000000;
    dz = -30.750000;
    TGeoCombiTrans* pMatrix57 = new TGeoCombiTrans("");
    pMatrix57->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode20 =
        new TGeoUnion(pG2403001_Fond_cote_entree_147_55, pG2403001_Fond_cote_entree_47_8, 0, pMatrix57);
    // Shape: G2403001_Fond_cote_entree_148 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_148_56 = new TGeoCompositeShape("G2403001_Fond_cote_entree_148", pBoolNode20);
    // Shape: G2403001_Fond_cote_entree_48 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 11.750000;
    TGeoShape* pG2403001_Fond_cote_entree_48_9 = new TGeoBBox("G2403001_Fond_cote_entree_48", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -120.000000;
    dz = 30.750000;
    TGeoCombiTrans* pMatrix59 = new TGeoCombiTrans("");
    pMatrix59->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode21 =
        new TGeoUnion(pG2403001_Fond_cote_entree_148_56, pG2403001_Fond_cote_entree_48_9, 0, pMatrix59);
    // Shape: G2403001_Fond_cote_entree_149 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_149_57 = new TGeoCompositeShape("G2403001_Fond_cote_entree_149", pBoolNode21);
    // Shape: G2403001_Fond_cote_entree_49 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 7.750000;
    TGeoShape* pG2403001_Fond_cote_entree_49_10 = new TGeoBBox("G2403001_Fond_cote_entree_49", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -140.000000;
    dz = -26.750000;
    TGeoCombiTrans* pMatrix61 = new TGeoCombiTrans("");
    pMatrix61->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode22 =
        new TGeoUnion(pG2403001_Fond_cote_entree_149_57, pG2403001_Fond_cote_entree_49_10, 0, pMatrix61);
    // Shape: G2403001_Fond_cote_entree_150 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_150_58 = new TGeoCompositeShape("G2403001_Fond_cote_entree_150", pBoolNode22);
    // Shape: G2403001_Fond_cote_entree_50 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 7.750000;
    TGeoShape* pG2403001_Fond_cote_entree_50_11 = new TGeoBBox("G2403001_Fond_cote_entree_50", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -140.000000;
    dz = 26.750000;
    TGeoCombiTrans* pMatrix63 = new TGeoCombiTrans("");
    pMatrix63->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode23 =
        new TGeoUnion(pG2403001_Fond_cote_entree_150_58, pG2403001_Fond_cote_entree_50_11, 0, pMatrix63);
    // Shape: G2403001_Fond_cote_entree_151 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_151_59 = new TGeoCompositeShape("G2403001_Fond_cote_entree_151", pBoolNode23);
    // Shape: G2403001_Fond_cote_entree_51 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 25.500000;
    TGeoShape* pG2403001_Fond_cote_entree_51_12 = new TGeoBBox("G2403001_Fond_cote_entree_51", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 160.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix65 = new TGeoCombiTrans("");
    pMatrix65->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode24 =
        new TGeoUnion(pG2403001_Fond_cote_entree_151_59, pG2403001_Fond_cote_entree_51_12, 0, pMatrix65);
    // Shape: G2403001_Fond_cote_entree_152 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_152_60 = new TGeoCompositeShape("G2403001_Fond_cote_entree_152", pBoolNode24);
    // Shape: G2403001_Fond_cote_entree_52 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 18.250000;
    TGeoShape* pG2403001_Fond_cote_entree_52_13 = new TGeoBBox("G2403001_Fond_cote_entree_52", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -60.000000;
    dz = 39.750000;
    TGeoCombiTrans* pMatrix67 = new TGeoCombiTrans("");
    pMatrix67->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode25 =
        new TGeoUnion(pG2403001_Fond_cote_entree_152_60, pG2403001_Fond_cote_entree_52_13, 0, pMatrix67);
    // Shape: G2403001_Fond_cote_entree_153 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_153_61 = new TGeoCompositeShape("G2403001_Fond_cote_entree_153", pBoolNode25);
    // Shape: G2403001_Fond_cote_entree_53 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 18.250000;
    TGeoShape* pG2403001_Fond_cote_entree_53_14 = new TGeoBBox("G2403001_Fond_cote_entree_53", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -60.000000;
    dz = -39.750000;
    TGeoCombiTrans* pMatrix69 = new TGeoCombiTrans("");
    pMatrix69->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode26 =
        new TGeoUnion(pG2403001_Fond_cote_entree_153_61, pG2403001_Fond_cote_entree_53_14, 0, pMatrix69);
    // Shape: G2403001_Fond_cote_entree_154 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_154_62 = new TGeoCompositeShape("G2403001_Fond_cote_entree_154", pBoolNode26);
    // Shape: G2403001_Fond_cote_entree_54 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 18.250000;
    TGeoShape* pG2403001_Fond_cote_entree_54_15 = new TGeoBBox("G2403001_Fond_cote_entree_54", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 60.000000;
    dz = 39.750000;
    TGeoCombiTrans* pMatrix71 = new TGeoCombiTrans("");
    pMatrix71->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode27 =
        new TGeoUnion(pG2403001_Fond_cote_entree_154_62, pG2403001_Fond_cote_entree_54_15, 0, pMatrix71);
    // Shape: G2403001_Fond_cote_entree_155 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_155_63 = new TGeoCompositeShape("G2403001_Fond_cote_entree_155", pBoolNode27);
    // Shape: G2403001_Fond_cote_entree_55 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 54.000000;
    TGeoShape* pG2403001_Fond_cote_entree_55_16 = new TGeoBBox("G2403001_Fond_cote_entree_55", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -80.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix73 = new TGeoCombiTrans("");
    pMatrix73->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode28 =
        new TGeoUnion(pG2403001_Fond_cote_entree_155_63, pG2403001_Fond_cote_entree_55_16, 0, pMatrix73);
    // Shape: G2403001_Fond_cote_entree_156 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_156_64 = new TGeoCompositeShape("G2403001_Fond_cote_entree_156", pBoolNode28);
    // Shape: G2403001_Fond_cote_entree_56 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 18.250000;
    TGeoShape* pG2403001_Fond_cote_entree_56_17 = new TGeoBBox("G2403001_Fond_cote_entree_56", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 60.000000;
    dz = -39.750000;
    TGeoCombiTrans* pMatrix75 = new TGeoCombiTrans("");
    pMatrix75->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode29 =
        new TGeoUnion(pG2403001_Fond_cote_entree_156_64, pG2403001_Fond_cote_entree_56_17, 0, pMatrix75);
    // Shape: G2403001_Fond_cote_entree_157 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_157_65 = new TGeoCompositeShape("G2403001_Fond_cote_entree_157", pBoolNode29);
    // Shape: G2403001_Fond_cote_entree_57 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 25.500000;
    TGeoShape* pG2403001_Fond_cote_entree_57_18 = new TGeoBBox("G2403001_Fond_cote_entree_57", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -160.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix77 = new TGeoCombiTrans("");
    pMatrix77->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode30 =
        new TGeoUnion(pG2403001_Fond_cote_entree_157_65, pG2403001_Fond_cote_entree_57_18, 0, pMatrix77);
    // Shape: G2403001_Fond_cote_entree_158 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_158_66 = new TGeoCompositeShape("G2403001_Fond_cote_entree_158", pBoolNode30);
    // Shape: G2403001_Fond_cote_entree_58 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 49.000000;
    TGeoShape* pG2403001_Fond_cote_entree_58_19 = new TGeoBBox("G2403001_Fond_cote_entree_58", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 100.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix79 = new TGeoCombiTrans("");
    pMatrix79->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode31 =
        new TGeoUnion(pG2403001_Fond_cote_entree_158_66, pG2403001_Fond_cote_entree_58_19, 0, pMatrix79);
    // Shape: G2403001_Fond_cote_entree_159 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_159_67 = new TGeoCompositeShape("G2403001_Fond_cote_entree_159", pBoolNode31);
    // Shape: G2403001_Fond_cote_entree_59 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 49.000000;
    TGeoShape* pG2403001_Fond_cote_entree_59_20 = new TGeoBBox("G2403001_Fond_cote_entree_59", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = -100.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix81 = new TGeoCombiTrans("");
    pMatrix81->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode32 =
        new TGeoUnion(pG2403001_Fond_cote_entree_159_67, pG2403001_Fond_cote_entree_59_20, 0, pMatrix81);
    // Shape: G2403001_Fond_cote_entree_160 type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_160_68 = new TGeoCompositeShape("G2403001_Fond_cote_entree_160", pBoolNode32);
    // Shape: G2403001_Fond_cote_entree_60 type: TGeoBBox
    dx = 3.500000;
    dy = 0.750000;
    dz = 54.000000;
    TGeoShape* pG2403001_Fond_cote_entree_60_21 = new TGeoBBox("G2403001_Fond_cote_entree_60", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 80.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix83 = new TGeoCombiTrans("");
    pMatrix83->SetTranslation(dx, dy, dz);
    pMatrix83->RotateY(180.0);
    TGeoUnion* pBoolNode33 =
        new TGeoUnion(pG2403001_Fond_cote_entree_160_68, pG2403001_Fond_cote_entree_60_21, 0, pMatrix83);
    // Shape: G2403001_Fond_cote_entree_PartBody type: TGeoCompositeShape
    TGeoShape* pG2403001_Fond_cote_entree_PartBody_69 =
        new TGeoCompositeShape("G2403001_Fond_cote_entree_PartBody", pBoolNode33);
    // Volume: G2403001_Fond_cote_entree

    tr[0] = +184.000000;
    tr[1] = 0.000000;
    tr[2] = 0.000000;
    rot[0] = 1.000000;
    rot[1] = 0.000000;
    rot[2] = 0.000000;
    rot[3] = 0.000000;
    rot[4] = 1.000000;
    rot[5] = 0.000000;
    rot[6] = 0.000000;
    rot[7] = 0.000000;
    rot[8] = 1.000000;
    TGeoHMatrix* pMatrix84 = new TGeoHMatrix("");

    pMatrix84->RotateY(180.0);
    pMatrix84->SetTranslation(tr);
    pMatrix84->SetBit(TGeoMatrix::kGeoTranslation);
    pMatrix84->SetBit(TGeoMatrix::kGeoRotation);

    TGeoVolume* pG2403001_Fond_cote_entree =
        new TGeoVolume("G2403001_Fond_cote_entree", pG2403001_Fond_cote_entree_PartBody_69, pMedFe);
    pG2403001_Fond_cote_entree->SetLineColor(920);
    pG2403001_Fond_cote_entree->SetVisLeaves(kTRUE);
    pWorld->AddNode(pG2403001_Fond_cote_entree, 1, pMatrix84);
}

void ConstructGToles(TGeoVolume* pWorld)
{
    // Enceinte Extern part
    Double_t dx, dy, dz;
    Double_t dx1, dy1, dz1;
    Double_t dx2, dy2, dz2;
    Double_t rmin, rmax, phi1, phi2;
    Double_t thx, thy, thz, phx, phy, phz;

    // -- Mixture: Air
    TGeoMedium* pMed2 = gGeoManager->GetMedium("Air");
    TGeoMedium* pMedFe = gGeoManager->GetMedium("iron");

    Double_t tr[3], rot[9];

    // SHAPES, VOLUMES AND GEOMETRICAL HIERARCHY
    // Shape: GToles_28 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 38.000000;
    dz = 0.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_28_28 = new TGeoTubeSeg("GToles_28", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_30 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 373.000000;
    dz = 0.250000;
    phi1 = 180.000000;
    phi2 = 206.424390;
    TGeoShape* pGToles_30_29 = new TGeoTubeSeg("GToles_30", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 300.000000;
    dy = 149.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix3 = new TGeoCombiTrans("");
    pMatrix3->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode1 = new TGeoUnion(pGToles_28_28, pGToles_30_29, 0, pMatrix3);
    // Shape: GToles_89 type: TGeoCompositeShape
    TGeoShape* pGToles_89_46 = new TGeoCompositeShape("GToles_89", pBoolNode1);
    // Shape: GToles_31 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 373.000000;
    dz = 0.250000;
    phi1 = 180.000000;
    phi2 = 206.424390;
    TGeoShape* pGToles_31_30 = new TGeoTubeSeg("GToles_31", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.000000;
    dy = 149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix6 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix5 = new TGeoCombiTrans("", dx, dy, dz, pMatrix6);
    TGeoUnion* pBoolNode2 = new TGeoUnion(pGToles_89_46, pGToles_31_30, 0, pMatrix5);
    // Shape: GToles_27 type: TGeoCompositeShape
    TGeoShape* pGToles_27_47 = new TGeoCompositeShape("GToles_27", pBoolNode2);
    // Shape: GToles_26 type: TGeoBBox
    dx = 2.500000;
    dy = 10.000000;
    dz = 1.000000;
    TGeoShape* pGToles_26_2 = new TGeoBBox("GToles_26", dx, dy, dz);
    // Combi transformation:
    dx = -30.000000;
    dy = 24.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix39 = new TGeoCombiTrans("");
    pMatrix39->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode3 = new TGeoSubtraction(pGToles_27_47, pGToles_26_2, 0, pMatrix39);
    // Shape: GToles_120 type: TGeoCompositeShape
    TGeoShape* pGToles_120_62 = new TGeoCompositeShape("GToles_120", pBoolNode3);
    // Shape: GToles_24 type: TGeoBBox
    dx = 2.500000;
    dy = 10.000000;
    dz = 1.000000;
    TGeoShape* pGToles_24_3 = new TGeoBBox("GToles_24", dx, dy, dz);
    // Combi transformation:
    dx = 30.000000;
    dy = 24.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix41 = new TGeoCombiTrans("");
    pMatrix41->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode4 = new TGeoSubtraction(pGToles_120_62, pGToles_24_3, 0, pMatrix41);
    // Shape: GToles_123 type: TGeoCompositeShape
    TGeoShape* pGToles_123_63 = new TGeoCompositeShape("GToles_123", pBoolNode4);
    // Shape: GToles_22 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 21.000000;
    dz = 1.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_22_31 = new TGeoTubeSeg("GToles_22", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 16.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix43 = new TGeoCombiTrans("");
    pMatrix43->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode5 = new TGeoSubtraction(pGToles_123_63, pGToles_22_31, 0, pMatrix43);
    // Shape: GToles_124 type: TGeoCompositeShape
    TGeoShape* pGToles_124_64 = new TGeoCompositeShape("GToles_124", pBoolNode5);
    // Shape: GToles_20 type: TGeoBBox
    dx = 74.000000;
    dy = 2.000000;
    dz = 2.000000;
    TGeoShape* pGToles_20_4 = new TGeoBBox("GToles_20", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 148.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix45 = new TGeoCombiTrans("");
    pMatrix45->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode6 = new TGeoSubtraction(pGToles_124_64, pGToles_20_4, 0, pMatrix45);
    // Shape: GToles_122 type: TGeoCompositeShape
    TGeoShape* pGToles_122_65 = new TGeoCompositeShape("GToles_122", pBoolNode6);
    // Shape: GToles_18 type: TGeoBBox
    dx = 19.500000;
    dy = 11.875000;
    dz = 1.000000;
    TGeoShape* pGToles_18_5 = new TGeoBBox("GToles_18", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 135.205500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix47 = new TGeoCombiTrans("");
    pMatrix47->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode7 = new TGeoSubtraction(pGToles_122_65, pGToles_18_5, 0, pMatrix47);
    // Shape: GToles_121 type: TGeoCompositeShape
    TGeoShape* pGToles_121_66 = new TGeoCompositeShape("GToles_121", pBoolNode7);
    // Shape: GToles_139 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 38.000000;
    dz = 0.150000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_139_32 = new TGeoTubeSeg("GToles_139", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_138 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 373.000000;
    dz = 0.150000;
    phi1 = 180.000000;
    phi2 = 206.424000;
    TGeoShape* pGToles_138_33 = new TGeoTubeSeg("GToles_138", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 300.000000;
    dy = 149.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix8 = new TGeoCombiTrans("");
    pMatrix8->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode8 = new TGeoUnion(pGToles_139_32, pGToles_138_33, 0, pMatrix8);
    // Shape: GToles_137 type: TGeoCompositeShape
    TGeoShape* pGToles_137_48 = new TGeoCompositeShape("GToles_137", pBoolNode8);
    // Shape: GToles_136 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 373.000000;
    dz = 0.150000;
    phi1 = 180.000000;
    phi2 = 206.424000;
    TGeoShape* pGToles_136_34 = new TGeoTubeSeg("GToles_136", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.000000;
    dy = 149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix11 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix10 = new TGeoCombiTrans("", dx, dy, dz, pMatrix11);
    TGeoUnion* pBoolNode9 = new TGeoUnion(pGToles_137_48, pGToles_136_34, 0, pMatrix10);
    // Shape: GToles_135 type: TGeoCompositeShape
    TGeoShape* pGToles_135_49 = new TGeoCompositeShape("GToles_135", pBoolNode9);
    // Shape: GToles_134 type: TGeoBBox
    dx = 2.500000;
    dy = 10.000000;
    dz = 1.000000;
    TGeoShape* pGToles_134_6 = new TGeoBBox("GToles_134", dx, dy, dz);
    // Combi transformation:
    dx = -30.000000;
    dy = 24.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix49 = new TGeoCombiTrans("");
    pMatrix49->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode10 = new TGeoSubtraction(pGToles_135_49, pGToles_134_6, 0, pMatrix49);
    // Shape: GToles_133 type: TGeoCompositeShape
    TGeoShape* pGToles_133_67 = new TGeoCompositeShape("GToles_133", pBoolNode10);
    // Shape: GToles_132 type: TGeoBBox
    dx = 2.500000;
    dy = 10.000000;
    dz = 1.000000;
    TGeoShape* pGToles_132_7 = new TGeoBBox("GToles_132", dx, dy, dz);
    // Combi transformation:
    dx = 30.000000;
    dy = 24.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix51 = new TGeoCombiTrans("");
    pMatrix51->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode11 = new TGeoSubtraction(pGToles_133_67, pGToles_132_7, 0, pMatrix51);
    // Shape: GToles_131 type: TGeoCompositeShape
    TGeoShape* pGToles_131_68 = new TGeoCompositeShape("GToles_131", pBoolNode11);
    // Shape: GToles_130 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 21.000000;
    dz = 1.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_130_35 = new TGeoTubeSeg("GToles_130", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 0.000000;
    dy = 16.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix53 = new TGeoCombiTrans("");
    pMatrix53->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode = new TGeoSubtraction(pGToles_131_68, pGToles_130_35, 0, pMatrix53);
    // Shape: GToles_129 type: TGeoCompositeShape
    TGeoShape* pGToles_129_69 = new TGeoCompositeShape("GToles_129", pBoolNode);
    // Shape: GToles_128 type: TGeoBBox
    dx = 74.000000;
    dy = 2.000000;
    dz = 2.000000;
    TGeoShape* pGToles_128_8 = new TGeoBBox("GToles_128", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 148.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix55 = new TGeoCombiTrans("");
    pMatrix55->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode12 = new TGeoSubtraction(pGToles_129_69, pGToles_128_8, 0, pMatrix55);
    // Shape: GToles_127 type: TGeoCompositeShape
    TGeoShape* pGToles_127_70 = new TGeoCompositeShape("GToles_127", pBoolNode12);
    // Shape: GToles_126 type: TGeoBBox
    dx = 19.500000;
    dy = 11.875000;
    dz = 1.000000;
    TGeoShape* pGToles_126_9 = new TGeoBBox("GToles_126", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 135.205500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix57 = new TGeoCombiTrans("");
    pMatrix57->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode13 = new TGeoSubtraction(pGToles_127_70, pGToles_126_9, 0, pMatrix57);
    // Shape: GToles_125 type: TGeoCompositeShape
    TGeoShape* pGToles_125_71 = new TGeoCompositeShape("GToles_125", pBoolNode13);
    // Combi transformation:
    dx = 0.000000;
    dy = 298.161000;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix90 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix89 = new TGeoCombiTrans("", dx, dy, dz, pMatrix90);
    TGeoUnion* pBoolNode14 = new TGeoUnion(pGToles_121_66, pGToles_125_71, 0, pMatrix89);
    // Shape: GToles_140 type: TGeoCompositeShape
    TGeoShape* pGToles_140_86 = new TGeoCompositeShape("GToles_140", pBoolNode14);
    // Shape: GToles_57 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 78.000000;
    dz = 0.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_57_36 = new TGeoTubeSeg("GToles_57", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_55 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 413.000000;
    dz = 0.250000;
    phi1 = 333.575610;
    phi2 = 360.000000;
    TGeoShape* pGToles_55_37 = new TGeoTubeSeg("GToles_55", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 300.000000;
    dy = 149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix14 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix13 = new TGeoCombiTrans("", dx, dy, dz, pMatrix14);
    TGeoUnion* pBoolNode15 = new TGeoUnion(pGToles_57_36, pGToles_55_37, 0, pMatrix13);
    // Shape: GToles_141 type: TGeoCompositeShape
    TGeoShape* pGToles_141_50 = new TGeoCompositeShape("GToles_141", pBoolNode15);
    // Shape: GToles_58 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 413.000000;
    dz = 0.250000;
    phi1 = 333.575610;
    phi2 = 360.000000;
    TGeoShape* pGToles_58_38 = new TGeoTubeSeg("GToles_58", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.000000;
    dy = 149.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix16 = new TGeoCombiTrans("");
    pMatrix16->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode16 = new TGeoUnion(pGToles_141_50, pGToles_58_38, 0, pMatrix16);
    // Shape: GToles_142 type: TGeoCompositeShape
    TGeoShape* pGToles_142_51 = new TGeoCompositeShape("GToles_142", pBoolNode16);
    // Shape: GToles_149 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.000000;
    dz = 1.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_149_39 = new TGeoTubeSeg("GToles_149", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_148 type: TGeoBBox
    dx = 8.000000;
    dy = 11.000000;
    dz = 1.000000;
    TGeoShape* pGToles_148_10 = new TGeoBBox("GToles_148", dx, dy, dz);
    // Combi transformation:
    dx = -8.000000;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix18 = new TGeoCombiTrans("");
    pMatrix18->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode17 = new TGeoUnion(pGToles_149_39, pGToles_148_10, 0, pMatrix18);
    // Shape: GToles_147 type: TGeoCompositeShape
    TGeoShape* pGToles_147_52 = new TGeoCompositeShape("GToles_147", pBoolNode17);
    // Shape: GToles_146 type: TGeoBBox
    dx = 2.500000;
    dy = 5.000000;
    dz = 1.000000;
    TGeoShape* pGToles_146_11 = new TGeoBBox("GToles_146", dx, dy, dz);
    // Combi transformation:
    dx = 12.996800;
    dy = 15.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix20 = new TGeoCombiTrans("");
    pMatrix20->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode18 = new TGeoUnion(pGToles_147_52, pGToles_146_11, 0, pMatrix20);
    // Shape: GToles_145 type: TGeoCompositeShape
    TGeoShape* pGToles_145_53 = new TGeoCompositeShape("GToles_145", pBoolNode18);
    // Combi transformation:
    dx = 100.000000;
    dy = 94.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix60 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix59 = new TGeoCombiTrans("", dx, dy, dz, pMatrix60);
    TGeoSubtraction* pBoolNode19 = new TGeoSubtraction(pGToles_142_51, pGToles_145_53, 0, pMatrix59);
    // Shape: GToles_150 type: TGeoCompositeShape
    TGeoShape* pGToles_150_72 = new TGeoCompositeShape("GToles_150", pBoolNode19);
    // Shape: GToles_39 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.000000;
    dz = 1.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_39_40 = new TGeoTubeSeg("GToles_39", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_38 type: TGeoBBox
    dx = 8.000000;
    dy = 11.000000;
    dz = 1.000000;
    TGeoShape* pGToles_38_12 = new TGeoBBox("GToles_38", dx, dy, dz);
    // Combi transformation:
    dx = -8.000000;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix22 = new TGeoCombiTrans("");
    pMatrix22->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode20 = new TGeoUnion(pGToles_39_40, pGToles_38_12, 0, pMatrix22);
    // Shape: GToles_143 type: TGeoCompositeShape
    TGeoShape* pGToles_143_54 = new TGeoCompositeShape("GToles_143", pBoolNode20);
    // Shape: GToles_51 type: TGeoBBox
    dx = 2.500000;
    dy = 5.000000;
    dz = 1.000000;
    TGeoShape* pGToles_51_13 = new TGeoBBox("GToles_51", dx, dy, dz);
    // Combi transformation:
    dx = 12.996800;
    dy = 15.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix24 = new TGeoCombiTrans("");
    pMatrix24->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode21 = new TGeoUnion(pGToles_143_54, pGToles_51_13, 0, pMatrix24);
    // Shape: GToles_144 type: TGeoCompositeShape
    TGeoShape* pGToles_144_55 = new TGeoCompositeShape("GToles_144", pBoolNode21);
    // Combi transformation:
    dx = -100.000000;
    dy = 94.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix62 = new TGeoCombiTrans("");
    pMatrix62->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode22 = new TGeoSubtraction(pGToles_150_72, pGToles_144_55, 0, pMatrix62);
    // Shape: GToles_151 type: TGeoCompositeShape
    TGeoShape* pGToles_151_73 = new TGeoCompositeShape("GToles_151", pBoolNode22);
    // Shape: GToles_49 type: TGeoBBox
    dx = 115.000000;
    dy = 2.000000;
    dz = 2.000000;
    TGeoShape* pGToles_49_14 = new TGeoBBox("GToles_49", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 148.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix64 = new TGeoCombiTrans("");
    pMatrix64->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode23 = new TGeoSubtraction(pGToles_151_73, pGToles_49_14, 0, pMatrix64);
    // Shape: GToles_152 type: TGeoCompositeShape
    TGeoShape* pGToles_152_74 = new TGeoCompositeShape("GToles_152", pBoolNode23);
    // Shape: GToles_41 type: TGeoBBox
    dx = 10.000000;
    dy = 10.000000;
    dz = 1.000000;
    TGeoShape* pGToles_41_15 = new TGeoBBox("GToles_41", dx, dy, dz);
    // Combi transformation:
    dx = 104.000000;
    dy = 137.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix66 = new TGeoCombiTrans("");
    pMatrix66->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode24 = new TGeoSubtraction(pGToles_152_74, pGToles_41_15, 0, pMatrix66);
    // Shape: GToles_153 type: TGeoCompositeShape
    TGeoShape* pGToles_153_75 = new TGeoCompositeShape("GToles_153", pBoolNode24);
    // Shape: GToles_43 type: TGeoBBox
    dx = 5.000000;
    dy = 12.500000;
    dz = 1.000000;
    TGeoShape* pGToles_43_16 = new TGeoBBox("GToles_43", dx, dy, dz);
    // Combi transformation:
    dx = 89.000000;
    dy = 134.580500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix68 = new TGeoCombiTrans("");
    pMatrix68->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode25 = new TGeoSubtraction(pGToles_153_75, pGToles_43_16, 0, pMatrix68);
    // Shape: GToles_154 type: TGeoCompositeShape
    TGeoShape* pGToles_154_76 = new TGeoCompositeShape("GToles_154", pBoolNode25);
    // Shape: GToles_45 type: TGeoBBox
    dx = 50.000000;
    dy = 84.500000;
    dz = 2.000000;
    TGeoShape* pGToles_45_17 = new TGeoBBox("GToles_45", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 63.580500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix70 = new TGeoCombiTrans("");
    pMatrix70->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode66 = new TGeoSubtraction(pGToles_154_76, pGToles_45_17, 0, pMatrix70);
    // Shape: GToles_155 type: TGeoCompositeShape
    TGeoShape* pGToles_155_77 = new TGeoCompositeShape("GToles_155", pBoolNode66);
    // Shape: GToles_47 type: TGeoBBox
    dx = 7.500000;
    dy = 6.500000;
    dz = 1.000000;
    TGeoShape* pGToles_47_18 = new TGeoBBox("GToles_47", dx, dy, dz);
    // Combi transformation:
    dx = -93.500000;
    dy = 140.580500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix72 = new TGeoCombiTrans("");
    pMatrix72->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode26 = new TGeoSubtraction(pGToles_155_77, pGToles_47_18, 0, pMatrix72);
    // Shape: GToles_156 type: TGeoCompositeShape
    TGeoShape* pGToles_156_78 = new TGeoCompositeShape("GToles_156", pBoolNode26);
    // Combi transformation:
    dx = 0.000000;
    dy = 0.000000;
    dz = -339.100000;
    TGeoCombiTrans* pMatrix92 = new TGeoCombiTrans("");
    pMatrix92->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode27 = new TGeoUnion(pGToles_140_86, pGToles_156_78, 0, pMatrix92);
    // Shape: GToles_184 type: TGeoCompositeShape
    TGeoShape* pGToles_184_87 = new TGeoCompositeShape("GToles_184", pBoolNode27);
    // Shape: GToles_183 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 78.000000;
    dz = 0.250000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_183_41 = new TGeoTubeSeg("GToles_183", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_182 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 413.000000;
    dz = 0.250000;
    phi1 = 333.575610;
    phi2 = 360.000000;
    TGeoShape* pGToles_182_42 = new TGeoTubeSeg("GToles_182", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = 300.000000;
    dy = 149.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix27 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix26 = new TGeoCombiTrans("", dx, dy, dz, pMatrix27);
    TGeoUnion* pBoolNode28 = new TGeoUnion(pGToles_183_41, pGToles_182_42, 0, pMatrix26);
    // Shape: GToles_181 type: TGeoCompositeShape
    TGeoShape* pGToles_181_56 = new TGeoCompositeShape("GToles_181", pBoolNode28);
    // Shape: GToles_180 type: TGeoTubeSeg
    rmin = 300.000000;
    rmax = 413.000000;
    dz = 0.250000;
    phi1 = 333.575610;
    phi2 = 360.000000;
    TGeoShape* pGToles_180_43 = new TGeoTubeSeg("GToles_180", rmin, rmax, dz, phi1, phi2);
    // Combi transformation:
    dx = -300.000000;
    dy = 149.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix29 = new TGeoCombiTrans("");
    pMatrix29->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode29 = new TGeoUnion(pGToles_181_56, pGToles_180_43, 0, pMatrix29);
    // Shape: GToles_179 type: TGeoCompositeShape
    TGeoShape* pGToles_179_57 = new TGeoCompositeShape("GToles_179", pBoolNode29);
    // Shape: GToles_178 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.000000;
    dz = 1.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_178_44 = new TGeoTubeSeg("GToles_178", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_177 type: TGeoBBox
    dx = 8.000000;
    dy = 11.000000;
    dz = 1.000000;
    TGeoShape* pGToles_177_19 = new TGeoBBox("GToles_177", dx, dy, dz);
    // Combi transformation:
    dx = -8.000000;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix31 = new TGeoCombiTrans("");
    pMatrix31->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode30 = new TGeoUnion(pGToles_178_44, pGToles_177_19, 0, pMatrix31);
    // Shape: GToles_176 type: TGeoCompositeShape
    TGeoShape* pGToles_176_58 = new TGeoCompositeShape("GToles_176", pBoolNode30);
    // Shape: GToles_175 type: TGeoBBox
    dx = 2.500000;
    dy = 5.000000;
    dz = 1.000000;
    TGeoShape* pGToles_175_20 = new TGeoBBox("GToles_175", dx, dy, dz);
    // Combi transformation:
    dx = 12.996800;
    dy = 15.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix33 = new TGeoCombiTrans("");
    pMatrix33->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode31 = new TGeoUnion(pGToles_176_58, pGToles_175_20, 0, pMatrix33);
    // Shape: GToles_174 type: TGeoCompositeShape
    TGeoShape* pGToles_174_59 = new TGeoCompositeShape("GToles_174", pBoolNode31);
    // Combi transformation:
    dx = 100.000000;
    dy = 94.080500;
    dz = 0.000000;
    // Rotation:
    thx = 90.000000;
    phx = 180.000000;
    thy = 90.000000;
    phy = 90.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix75 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix74 = new TGeoCombiTrans("", dx, dy, dz, pMatrix75);
    TGeoSubtraction* pBoolNode32 = new TGeoSubtraction(pGToles_179_57, pGToles_174_59, 0, pMatrix74);
    // Shape: GToles_173 type: TGeoCompositeShape
    TGeoShape* pGToles_173_79 = new TGeoCompositeShape("GToles_173", pBoolNode32);
    // Shape: GToles_172 type: TGeoTubeSeg
    rmin = 0.000000;
    rmax = 11.000000;
    dz = 1.000000;
    phi1 = 0.000000;
    phi2 = 360.000000;
    TGeoShape* pGToles_172_45 = new TGeoTubeSeg("GToles_172", rmin, rmax, dz, phi1, phi2);
    // Shape: GToles_171 type: TGeoBBox
    dx = 8.000000;
    dy = 11.000000;
    dz = 1.000000;
    TGeoShape* pGToles_171_21 = new TGeoBBox("GToles_171", dx, dy, dz);
    // Combi transformation:
    dx = -8.000000;
    dy = 0.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix35 = new TGeoCombiTrans("");
    pMatrix35->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode33 = new TGeoUnion(pGToles_172_45, pGToles_171_21, 0, pMatrix35);
    // Shape: GToles_170 type: TGeoCompositeShape
    TGeoShape* pGToles_170_60 = new TGeoCompositeShape("GToles_170", pBoolNode33);
    // Shape: GToles_169 type: TGeoBBox
    dx = 2.500000;
    dy = 5.000000;
    dz = 1.000000;
    TGeoShape* pGToles_169_22 = new TGeoBBox("GToles_169", dx, dy, dz);
    // Combi transformation:
    dx = 12.996800;
    dy = 15.000000;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix37 = new TGeoCombiTrans("");
    pMatrix37->SetTranslation(dx, dy, dz);
    TGeoUnion* pBoolNode34 = new TGeoUnion(pGToles_170_60, pGToles_169_22, 0, pMatrix37);
    // Shape: GToles_168 type: TGeoCompositeShape
    TGeoShape* pGToles_168_61 = new TGeoCompositeShape("GToles_168", pBoolNode34);
    // Combi transformation:
    dx = -100.000000;
    dy = 94.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix77 = new TGeoCombiTrans("");
    pMatrix77->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode35 = new TGeoSubtraction(pGToles_173_79, pGToles_168_61, 0, pMatrix77);
    // Shape: GToles_167 type: TGeoCompositeShape
    TGeoShape* pGToles_167_80 = new TGeoCompositeShape("GToles_167", pBoolNode35);
    // Shape: GToles_166 type: TGeoBBox
    dx = 114.000000;
    dy = 2.000000;
    dz = 2.000000;
    TGeoShape* pGToles_166_23 = new TGeoBBox("GToles_166", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 148.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix79 = new TGeoCombiTrans("");
    pMatrix79->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode36 = new TGeoSubtraction(pGToles_167_80, pGToles_166_23, 0, pMatrix79);
    // Shape: GToles_165 type: TGeoCompositeShape
    TGeoShape* pGToles_165_81 = new TGeoCompositeShape("GToles_165", pBoolNode36);
    // Shape: GToles_164 type: TGeoBBox
    dx = 10.000000;
    dy = 10.000000;
    dz = 1.000000;
    TGeoShape* pGToles_164_24 = new TGeoBBox("GToles_164", dx, dy, dz);
    // Combi transformation:
    dx = 104.000000;
    dy = 137.080500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix81 = new TGeoCombiTrans("");
    pMatrix81->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode37 = new TGeoSubtraction(pGToles_165_81, pGToles_164_24, 0, pMatrix81);
    // Shape: GToles_163 type: TGeoCompositeShape
    TGeoShape* pGToles_163_82 = new TGeoCompositeShape("GToles_163", pBoolNode37);
    // Shape: GToles_162 type: TGeoBBox
    dx = 5.000000;
    dy = 12.500000;
    dz = 1.000000;
    TGeoShape* pGToles_162_25 = new TGeoBBox("GToles_162", dx, dy, dz);
    // Combi transformation:
    dx = 89.000000;
    dy = 134.580500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix83 = new TGeoCombiTrans("");
    pMatrix83->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode38 = new TGeoSubtraction(pGToles_163_82, pGToles_162_25, 0, pMatrix83);
    // Shape: GToles_161 type: TGeoCompositeShape
    TGeoShape* pGToles_161_83 = new TGeoCompositeShape("GToles_161", pBoolNode38);
    // Shape: GToles_160 type: TGeoBBox
    dx = 50.000000;
    dy = 83.500000;
    dz = 2.000000;
    TGeoShape* pGToles_160_26 = new TGeoBBox("GToles_160", dx, dy, dz);
    // Combi transformation:
    dx = 0.000000;
    dy = 63.580500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix85 = new TGeoCombiTrans("");
    pMatrix85->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode39 = new TGeoSubtraction(pGToles_161_83, pGToles_160_26, 0, pMatrix85);
    // Shape: GToles_159 type: TGeoCompositeShape
    TGeoShape* pGToles_159_84 = new TGeoCompositeShape("GToles_159", pBoolNode39);
    // Shape: GToles_158 type: TGeoBBox
    dx = 7.500000;
    dy = 6.500000;
    dz = 1.000000;
    TGeoShape* pGToles_158_27 = new TGeoBBox("GToles_158", dx, dy, dz);
    // Combi transformation:
    dx = -93.500000;
    dy = 140.580500;
    dz = 0.000000;
    TGeoCombiTrans* pMatrix87 = new TGeoCombiTrans("");
    pMatrix87->SetTranslation(dx, dy, dz);
    TGeoSubtraction* pBoolNode40 = new TGeoSubtraction(pGToles_159_84, pGToles_158_27, 0, pMatrix87);
    // Shape: GToles_157 type: TGeoCompositeShape
    TGeoShape* pGToles_157_85 = new TGeoCompositeShape("GToles_157", pBoolNode40);
    // Combi transformation:
    dx = 0.000000;
    dy = 298.161000;
    dz = -339.100000;
    // Rotation:
    thx = 90.000000;
    phx = 0.000000;
    thy = 90.000000;
    phy = 270.000000;
    thz = 180.000000;
    phz = 0.000000;
    TGeoRotation* pMatrix95 = new TGeoRotation("", thx, phx, thy, phy, thz, phz);
    TGeoCombiTrans* pMatrix94 = new TGeoCombiTrans("", dx, dy, dz, pMatrix95);
    TGeoUnion* pBoolNode41 = new TGeoUnion(pGToles_184_87, pGToles_157_85, 0, pMatrix94);
    // Shape: GToles_PartBody type: TGeoCompositeShape
    TGeoShape* pGToles_PartBody_88 = new TGeoCompositeShape("GToles_PartBody", pBoolNode41);

    // HMatrix:
    tr[0] = 0.000000;
    tr[1] = -150.000000;
    tr[2] = 185.000000;
    rot[0] = 1.000000;
    rot[1] = 0.000000;
    rot[2] = 0.000000;
    rot[3] = 0.000000;
    rot[4] = 1.000000;
    rot[5] = 0.000000;
    rot[6] = 0.000000;
    rot[7] = 0.000000;
    rot[8] = 1.000000;
    TGeoHMatrix* pMatrix96 = new TGeoHMatrix("");

    pMatrix96->SetTranslation(tr);
    pMatrix96->SetRotation(rot);
    pMatrix96->RotateY(90.0);
    pMatrix96->SetBit(TGeoMatrix::kGeoTranslation);
    pMatrix96->SetBit(TGeoMatrix::kGeoRotation);
    // Shape: World_PartBody type: TGeoBBox

    // Volume: GToles
    TGeoVolume* pGToles = new TGeoVolume("GToles", pGToles_PartBody_88, pMedFe);
    pGToles->SetLineColor(920);
    pGToles->SetVisLeaves(kTRUE);
    pWorld->AddNode(pGToles, 1, pMatrix96);
}
