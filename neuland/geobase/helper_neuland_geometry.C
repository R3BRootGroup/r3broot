/******************************************************************************
 *   Copyright (C) 2023 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2023 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

const Int_t gPaddlesPerPlane = 50;

const Double_t gPaddleDistance = 2.5;
const Double_t gPaddleBaseLength = 125.0;
const Double_t gPaddleConeLength = 5.;

const Double_t gBC408ConeRadius = 1.2;
const Double_t gBC408Thickness = 2.4;
const Double_t gAlThickness = 0.02;
const Double_t gTapeThickness = 0.05;

TGeoMedium* BuildMaterial(const TString material, FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild)
{
    FairGeoMedium* fairMedium = geoMedia->getMedium(material);
    if (!fairMedium)
    {
        std::cout << "FairGeoMedium " << material << " not found" << std::endl;
    }

    geoBuild->createMedium(fairMedium);
    TGeoMedium* med = gGeoManager->GetMedium(material);
    if (!med)
    {
        std::cout << "TGeoMedium " << material << " not found" << std::endl;
    }
    return med;
}

/** Building the shape of a paddle
    The paddle end cones are created by intersection ("*") of a small box
   (width*width*coneLength) with a cone The intersected structures are moved,
   rotated and then to a base box (length*width*width) */
TGeoShape* BuildPaddleShape(const TString& name,
                            const Double_t length,
                            const Double_t width,
                            const Double_t coneRadius,
                            const Double_t coneLength)
{
    new TGeoBBox(name + "Box", length, width, width);
    new TGeoCone(name + "Cone", coneLength + 0.001, 0., coneRadius, 0., width * TMath::Sqrt(2.));
    new TGeoBBox(name + "Conebox", width, width, coneLength);
    TGeoShape* shape = new TGeoCompositeShape(
        name,
        name + "Box + ((" + name + "Conebox*" + name + "Cone):trc1) + ((" + name + "Conebox*" + name + "Cone):trc2)");
    return shape;
}

/** A paddle consists of the scintillator volume, wrapped with Al and Tape.
    The shapes for the wrappings are created by subtracting the solid inner part
   from a larger solid volume */
TGeoVolume* BuildPaddleVolume()
{
    // Load Interfaces to build materials
    FairGeoLoader* geoLoad = FairGeoLoader::Instance();
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();
    geoFace->setMediaFile(TString(gSystem->Getenv("VMCWORKDIR")) + "/geometry/media_r3b.geo");
    geoFace->readMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();
    FairGeoMedia* geoMedia = geoFace->getMedia();
    const TGeoMedium* medBC408 = BuildMaterial("BC408", geoMedia, geoBuild);
    const TGeoMedium* medCH2 = BuildMaterial("polyethylene", geoMedia, geoBuild);
    const TGeoMedium* medAl = BuildMaterial("aluminium", geoMedia, geoBuild);

    // Prepare Transformations for cones
    TGeoRotation* r1 = new TGeoRotation();
    r1->RotateY(90);
    TGeoCombiTrans* trc1 = new TGeoCombiTrans(TGeoTranslation(-(gPaddleBaseLength + gPaddleConeLength), 0., 0.), *r1);
    trc1->SetName("trc1");
    trc1->RegisterYourself();

    TGeoRotation* r2 = new TGeoRotation();
    r2->RotateY(-90);
    TGeoCombiTrans* trc2 = new TGeoCombiTrans(TGeoTranslation(+(gPaddleBaseLength + gPaddleConeLength), 0., 0.), *r2);
    trc2->SetName("trc2");
    trc2->RegisterYourself();

    // Build shapes
    const TGeoShape* shapeBC408 =
        BuildPaddleShape("shapeBC408", gPaddleBaseLength, gBC408Thickness, gBC408ConeRadius, gPaddleConeLength);

    BuildPaddleShape("shapeAlWrappingSolid",
                     gPaddleBaseLength,
                     gBC408Thickness + gAlThickness,
                     gBC408ConeRadius + gAlThickness,
                     gPaddleConeLength);
    const TGeoShape* shapeAlWrapping = new TGeoCompositeShape("shapeAlWrapping", "shapeAlWrappingSolid - shapeBC408");

    BuildPaddleShape("shapeTapeWrappingSolid",
                     gPaddleBaseLength,
                     gBC408Thickness + gAlThickness + gTapeThickness,
                     gBC408ConeRadius + gAlThickness + gTapeThickness,
                     gPaddleConeLength);
    const TGeoShape* shapeTapeWrapping =
        new TGeoCompositeShape("shapeTapeWrapping", "shapeTapeWrappingSolid - shapeAlWrappingSolid");

    // Build Volumes. Note that the volume names are important and reappear in
    // R3BNeuland()
    TGeoVolume* volBC408 = new TGeoVolume("volBC408", shapeBC408, medBC408);
    volBC408->SetLineColor(33); // greyish plue
    volBC408->SetTransparency(30);
    TGeoVolume* volAlWrapping = new TGeoVolume("volAlWrapping", shapeAlWrapping, medAl);
    volAlWrapping->SetLineColor(17); // grey/silver
    TGeoVolume* volTapeWrapping = new TGeoVolume("volTapeWrapping", shapeTapeWrapping, medCH2);
    volTapeWrapping->SetLineColor(1); // black

    // Make the elementary assembly
    TGeoVolume* volPaddle = new TGeoVolumeAssembly("volPaddle");
    volPaddle->AddNode(volBC408, 1);
    volPaddle->AddNode(volAlWrapping, 1);
    volPaddle->AddNode(volTapeWrapping, 1);

    return volPaddle;
}
