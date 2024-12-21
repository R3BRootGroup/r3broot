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

// Includes from C
#include <fstream>
#include <iomanip>
#include <iostream>

// Includes from ROOT
#include "Math/RotationZYX.h"
#include "TArrayD.h"
#include "TArrayF.h"
#include "TFile.h"
#include "TMath.h"

#include "R3BGladFieldMap.h"
#include "R3BLogger.h"

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::ofstream;
using std::right;
using std::setw;
using std::showpoint;
using TMath::Nint;

// -------------   Default constructor  ----------------------------------
//
R3BGladFieldMap::R3BGladFieldMap()
{
    fPosX = fPosY = fPosZ = 0.;
    fXmin = fYmin = fZmin = 0.;
    fXmax = fYmax = fZmax = 0.;
    fXstep = fYstep = fZstep = 0.;
    fNx = fNy = fNz = 0;
    fScale = 1.;
    fBx = fBy = fBz = NULL;
    fPosX = fPosY = fPosZ = 0.;
    fName = "";
    fFileName = "";
    fType = 2;
}

// -------------   Standard constructor   ---------------------------------
//
R3BGladFieldMap::R3BGladFieldMap(const TString& mapName, const TString& fileType)
    : FairField(mapName.Data())
{
    // Default field positions (in cm) in lab:
    // between target position (0,0,0) and GLAD rotation point (field origin)
    // Override these values by calling SetPosition(x,y,z) before Init()
    fPosX = 0.;
    fPosY = 1.75;
    fPosZ = 163.4;

    // Default translation vector of the local filed coordiantes
    gTrans = new TVector3(-fPosX, -fPosY, -fPosZ);

    // Default Euler rotations of the local field (in degrees)
    // Override these values by calling SetXAngle(), SetYAngle(), SetZAngle() before Init()
    fXAngle = 0.;
    fYAngle = -14.;
    fZAngle = 0.;

    fXmin = 0;
    fYmin = 0;
    fZmin = 0.;
    fXmax = 0;
    fYmax = 0;
    fZmax = 0.;
    fXstep = 0;
    fYstep = 0;
    fZstep = 0.;
    fNx = 0;
    fNy = 0;
    fNz = 0;
    fScale = 1.;
    fBx = nullptr;
    fBy = nullptr;
    fBz = nullptr;
    fName = mapName;
    TString dir = getenv("VMCWORKDIR");
    fFileName = dir + "/field/magField/R3B/" + mapName;
    if (fileType[0] == 'R')
    {
        fFileName += ".root";
    }
    else
    {
        fFileName += ".dat";
    }

    fType = 2;
}

void R3BGladFieldMap::SetFieldfromCurrent(double current)
{
    R3BLOG_IF(fatal, TMath::Abs(current) > 3583.81, "GLAD current cannot be larger than 3583.81 A.");
    fScale = current / 3583.81;
    R3BLOG(info, "GLAD current set to " << current << " A, which corresponds to a scaling factor of " << fScale);
    return;
}

// ------------   Constructor from R3BGladFieldPar   -----------------------
//
R3BGladFieldMap::R3BGladFieldMap(R3BFieldPar* fieldPar)
{
    fType = 2;
    fPosX = 0;
    fPosY = 0;
    fPosZ = 0.;
    fXAngle = 0;
    fYAngle = 0;
    fZAngle = 0.;
    fXmin = 0;
    fYmin = 0;
    fZmin = 0.;
    fXmax = 0;
    fYmax = 0;
    fZmax = 0.;
    fXstep = 0;
    fYstep = 0;
    fZstep = 0.;
    fNx = 0;
    fNy = 0;
    fNz = 0;
    fScale = 1.;
    fBx = fBy = fBz = NULL;
    if (!fieldPar)
    {
        R3BLOG(warn, "empty parameter container!");
        fName = "";
        fFileName = "";
        fType = 1;
    }
    else
    {
        fieldPar->MapName(fName);
        fFileName = fieldPar->GetFileName();
        fPosX = fieldPar->GetPositionX();
        fPosY = fieldPar->GetPositionY();
        fPosZ = fieldPar->GetPositionZ();
        fXAngle = fieldPar->GetXAngle();
        fYAngle = fieldPar->GetYAngle();
        fZAngle = fieldPar->GetZAngle();
        fScale = fieldPar->GetScale();
        fType = fieldPar->GetType();
    }
}

// ------------   Destructor   --------------------------------------------
//
R3BGladFieldMap::~R3BGladFieldMap()
{
    if (fBx)
        delete fBx;
    if (fBy)
        delete fBy;
    if (fBz)
        delete fBz;
}

// -----------   Intialisation   ------------------------------------------
//
void R3BGladFieldMap::Init()
{
    if (!gTrans)
        gTrans = new TVector3(-fPosX, -fPosY, -fPosZ);
    if (fFileName.EndsWith(".dat"))
        ReadAsciiFile(fFileName);
    else if (fFileName.EndsWith(".root"))
        ReadRootFile(fFileName);
    else
    {
        R3BLOG(fatal, "No proper file name defined! (" << fFileName.Data() << ")");
    }
    Print();
}
// -----------   Get x component of the field   ---------------------------
//
Double_t R3BGladFieldMap::GetBx(Double_t x, Double_t y, Double_t z)
{
    TVector3 B = GetBtrans(x, y, z);
    Double_t val = B.X() * 10.0;
    return (val); // should be in kGaus units
}
// -----------   Get y component of the field   ---------------------------
//
Double_t R3BGladFieldMap::GetBy(Double_t x, Double_t y, Double_t z)
{
    TVector3 B = GetBtrans(x, y, z);
    Double_t val = B.Y() * 10.0;
    return (val); // should be in kGaus units
}

// -----------   Get z component of the field   ---------------------------
//
Double_t R3BGladFieldMap::GetBz(Double_t x, Double_t y, Double_t z)
{
    TVector3 B = GetBtrans(x, y, z);
    Double_t val = B.Z() * 10.0;
    return (val); // should be in kGaus units
}

/* -----------   Get total B vector  of the field   -----------------------
 *
 * General function to get properly transformed B-field vector
 * The function translates and rotates the original field map using current values
 * of translation vector gTrans and Euler angles: fXAngle, fYAngle and fZAngle.
 * Input lab coordinates where the field should be calculated.
 * Returns B field vector with rotated components (in Tesla units)
 */
TVector3 R3BGladFieldMap::GetBtrans(Double_t x, Double_t y, Double_t z)
{
    TVector3 localPoint(x, y, z);
    localPoint += (*gTrans);

    TRotation r;
    // First Euler rotation around Y axis (default: -14 deg)
    r.RotateY(fYAngle * TMath::DegToRad());

    // get local X axis after first rotation
    TVector3 v3_localX(1, 0, 0);
    v3_localX.SetTheta(r.ThetaX());
    v3_localX.SetPhi(r.PhiX());

    // Second Euler rotation around local X axis
    r.Rotate(fXAngle * TMath::DegToRad(), v3_localX);

    // get local Z axis after second rotation
    TVector3 v3_localZ(0, 0, 1);
    v3_localZ.SetTheta(r.ThetaZ());
    v3_localZ.SetPhi(r.PhiZ());

    // final rotation around local Z axis
    r.Rotate(fZAngle * TMath::DegToRad(), v3_localZ);

    auto rr = r.Inverse();
    localPoint.Transform(rr);

    Int_t ix = 0;
    Int_t iy = 0;
    Int_t iz = 0;
    Double_t dx = 0.;
    Double_t dy = 0.;
    Double_t dz = 0.;

    TVector3 B(0., 0., 0.);

    if (!IsInside(localPoint.X(), localPoint.Y(), localPoint.Z(), ix, iy, iz, dx, dy, dz))
    {
        return B;
    }

    Double_t Bx = 0;
    Double_t By = 0;
    Double_t Bz = 0;

    // Get Bx field values at grid cell corners
    fHa[0][0][0] = fBx->At(ix * fNy * fNz + iy * fNz + iz);
    fHa[1][0][0] = fBx->At((ix + 1) * fNy * fNz + iy * fNz + iz);
    fHa[0][1][0] = fBx->At(ix * fNy * fNz + (iy + 1) * fNz + iz);
    fHa[1][1][0] = fBx->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + iz);
    fHa[0][0][1] = fBx->At(ix * fNy * fNz + iy * fNz + (iz + 1));
    fHa[1][0][1] = fBx->At((ix + 1) * fNy * fNz + iy * fNz + (iz + 1));
    fHa[0][1][1] = fBx->At(ix * fNy * fNz + (iy + 1) * fNz + (iz + 1));
    fHa[1][1][1] = fBx->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + (iz + 1));
    Bx = Interpolate(dx, dy, dz);

    // Get By field values at grid cell corners
    fHa[0][0][0] = fBy->At(ix * fNy * fNz + iy * fNz + iz);
    fHa[1][0][0] = fBy->At((ix + 1) * fNy * fNz + iy * fNz + iz);
    fHa[0][1][0] = fBy->At(ix * fNy * fNz + (iy + 1) * fNz + iz);
    fHa[1][1][0] = fBy->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + iz);
    fHa[0][0][1] = fBy->At(ix * fNy * fNz + iy * fNz + (iz + 1));
    fHa[1][0][1] = fBy->At((ix + 1) * fNy * fNz + iy * fNz + (iz + 1));
    fHa[0][1][1] = fBy->At(ix * fNy * fNz + (iy + 1) * fNz + (iz + 1));
    fHa[1][1][1] = fBy->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + (iz + 1));
    By = Interpolate(dx, dy, dz);

    // Get Bz field values at grid cell corners
    fHa[0][0][0] = fBz->At(ix * fNy * fNz + iy * fNz + iz);
    fHa[1][0][0] = fBz->At((ix + 1) * fNy * fNz + iy * fNz + iz);
    fHa[0][1][0] = fBz->At(ix * fNy * fNz + (iy + 1) * fNz + iz);
    fHa[1][1][0] = fBz->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + iz);
    fHa[0][0][1] = fBz->At(ix * fNy * fNz + iy * fNz + (iz + 1));
    fHa[1][0][1] = fBz->At((ix + 1) * fNy * fNz + iy * fNz + (iz + 1));
    fHa[0][1][1] = fBz->At(ix * fNy * fNz + (iy + 1) * fNz + (iz + 1));
    fHa[1][1][1] = fBz->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + (iz + 1));
    Bz = Interpolate(dx, dy, dz);

    // Set total B vector and transform accordingly
    B.SetXYZ(Bx * fScale, By * fScale, Bz * fScale);
    B.Transform(r);

    return B;
}

// -----------   Check whether a point is inside the map   ----------------
//
Bool_t R3BGladFieldMap::IsInside(Double_t x,
                                 Double_t y,
                                 Double_t z,
                                 Int_t& ix,
                                 Int_t& iy,
                                 Int_t& iz,
                                 Double_t& dx,
                                 Double_t& dy,
                                 Double_t& dz)
{

    // --- Transform into local coordinate system
    Double_t xl = x;
    Double_t yl = y;
    Double_t zl = z;

    // ---  Check for being outside the map range
    if (!(xl >= fXmin && xl < fXmax && yl >= fYmin && yl < fYmax && zl >= fZmin && zl < fZmax))
    {
        ix = iy = iz = 0;
        dx = dy = dz = 0.;
        return kFALSE;
    }

    // --- Determine grid cell
    ix = Int_t((xl - fXmin) / fXstep);
    iy = Int_t((yl - fYmin) / fYstep);
    iz = Int_t((zl - fZmin) / fZstep);

    // Relative distance from grid point (in units of cell size)
    dx = (xl - fXmin) / fXstep - Double_t(ix);
    dy = (yl - fYmin) / fYstep - Double_t(iy);
    dz = (zl - fZmin) / fZstep - Double_t(iz);

    // cout << "-I- isInside true " << endl;
    return kTRUE;
}

// ----------   Write the map to an ASCII file   --------------------------
//
void R3BGladFieldMap::WriteAsciiFile(const TString& fileName)
{
    // Open file
    R3BLOG(info, "Writing field map to ASCII file " << fileName.Data());
    ofstream mapFile(fileName);
    if (!mapFile.is_open())
    {
        R3BLOG(error, "Could not open file!");
        return;
    }

    // Write field map grid parameters
    mapFile.precision(6);
    mapFile << showpoint;
    if (fType == 1)
        mapFile << "nosym" << endl;
    if (fType == 2)
        mapFile << "sym2" << endl;
    if (fType == 3)
        mapFile << "sym3" << endl;
    mapFile << fXmin << " " << fXmax << " " << fNx << endl;
    mapFile << fYmin << " " << fYmax << " " << fNy << endl;
    mapFile << fZmin << " " << fZmax << " " << fNz << endl;

    // Write field values
    Double_t factor = 10. * fScale; // Takes out scaling and converts kG->T
    cout << right;
    Int_t nTot = fNx * fNy * fNz;
    R3BLOG(info, fNx * fNy * fNz << " entries to write... " << setw(3) << 0 << " % ");
    Int_t index = 0;
    div_t modul;
    Int_t iDiv = TMath::Nint(nTot / 100.);
    for (Int_t ix = 0; ix < fNx; ix++)
    {
        for (Int_t iy = 0; iy < fNy; iy++)
        {
            for (Int_t iz = 0; iz < fNz; iz++)
            {
                index = ix * fNy * fNz + iy * fNz + iz;
                modul = div(index, iDiv);
                if (modul.rem == 0)
                {
                    Double_t perc = TMath::Nint(100. * index / nTot);
                    cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
                }
                mapFile << fBx->At(index) / factor << " " << fBy->At(index) / factor << " " << fBz->At(index) / factor
                        << endl;
            } // z-Loop
        }     // y-Loop
    }         // x-Loop
    R3BLOG(info, "   " << index + 1 << " written");
    mapFile.close();
}

// -----  Set the position of the field centre in global coordinates  -----
//
void R3BGladFieldMap::SetPosition(Double_t x, Double_t y, Double_t z)
{
    fPosX = x;
    fPosY = y;
    fPosZ = z;
    gTrans->SetXYZ(-fPosX, -fPosY, -fPosZ);
}

// ---------   Screen output   --------------------------------------------
//
void R3BGladFieldMap::Print(Option_t*) const
{
    TString type = "Map";
    if (fType == 2)
    {
        type = "Map sym2";
    }
    else if (fType == 3)
    {
        type = "Map sym3";
    }

    std::stringstream sprint;
    sprint << endl;
    sprint << "======================================================" << endl;
    sprint.precision(4);
    sprint << showpoint;
    sprint << "----  " << fTitle << " : " << fName << endl;
    sprint << "----" << endl;
    sprint << "----  Field type     : " << type << endl;
    sprint << "----" << endl;
    sprint << "----  Field map grid : " << endl;
    sprint << "----  x = " << setw(4) << fXmin << " to " << setw(4) << fXmax << " cm, " << fNx
           << " grid points, dx = " << fXstep << " cm" << endl;
    sprint << "----  y = " << setw(4) << fYmin << " to " << setw(4) << fYmax << " cm, " << fNy
           << " grid points, dy = " << fYstep << " cm" << endl;
    sprint << "----  z = " << setw(4) << fZmin << " to " << setw(4) << fZmax << " cm, " << fNz
           << " grid points, dz = " << fZstep << " cm" << endl;
    sprint << endl;
    sprint << "----  Field centre position: ( " << setw(6) << fPosX << ", " << setw(6) << fPosY << ", " << setw(6)
           << fPosZ << ") cm" << endl;
    sprint << "----  Field rotation X: " << setw(6) << fXAngle << " deg" << endl;
    sprint << "----  Field rotation Y: " << setw(6) << fYAngle << " deg" << endl;
    sprint << "----  Field rotation Z: " << setw(6) << fZAngle << " deg" << endl;
    sprint << "----  Field scaling factor: " << fScale << endl;
    sprint << "======================================================" << endl;

    R3BLOG(info, sprint.str());
}

// ---------    Reset parameters and data (private)  ----------------------
//
void R3BGladFieldMap::Reset()
{
    fPosX = 0.;
    fPosY = 0.;
    fPosZ = 0.;
    fXmin = 0.;
    fYmin = 0.;
    fZmin = 0.;
    fXmax = 0.;
    fYmax = 0.;
    fZmax = 0.;
    fXstep = 0.;
    fYstep = 0.;
    fZstep = 0.;
    fNx = 0;
    fNy = 0;
    fNz = 0;
    fScale = 1.;
    if (fBx)
    {
        delete fBx;
        fBx = NULL;
    }
    if (fBy)
    {
        delete fBy;
        fBy = NULL;
    }
    if (fBz)
    {
        delete fBz;
        fBz = NULL;
    }
}

// -----   Read field map from ASCII file (private)   ---------------------
//
void R3BGladFieldMap::ReadAsciiFile(const TString& fileName)
{
    Double_t bx = 0., by = 0., bz = 0.;
    Double_t ax = 0., ay = 0., az = 0.;

    // Open file
    R3BLOG(info, "Reading field map from ASCII file " << fileName.Data());
    ifstream mapFile(fileName);
    if (!mapFile.is_open())
    {
        R3BLOG(fatal, "Could not open file! ");
    }
    // Read map type
    TString type;
    mapFile >> type;
    Int_t iType = 0;
    if (type == "nosym")
    {
        iType = 1;
    }
    else if (type == "sym2")
    {
        iType = 2;
    }
    else if (type == "sym3")
    {
        iType = 3;
    }

    // Read grid parameters
    mapFile >> fXmin >> fXmax >> fNx;
    mapFile >> fYmin >> fYmax >> fNy;
    mapFile >> fZmin >> fZmax >> fNz;
    fXstep = (fXmax - fXmin) / Double_t(fNx);
    fYstep = (fYmax - fYmin) / Double_t(fNy);
    fZstep = (fZmax - fZmin) / Double_t(fNz);

    // Create field arrays
    fNx += 1;
    fNy += 1;
    fNz += 1;
    fBx = new TArrayD(fNx * fNy * fNz);
    fBy = new TArrayD(fNx * fNy * fNz);
    fBz = new TArrayD(fNx * fNy * fNz);

    // Read the field values
    Double_t factor = fScale * 10.; // Factor 10 for T -> kG
    cout << right;
    Int_t nTot = fNx * fNy * fNz;
    R3BLOG(info, nTot << " entries to read... " << setw(3) << 0 << " % ");
    Int_t index = 0;
    div_t modul;
    Int_t iDiv = TMath::Nint(nTot / 100.);
    for (Int_t ix = 0; ix < fNx; ix++)
    {
        for (Int_t iy = 0; iy < fNy; iy++)
        {
            for (Int_t iz = 0; iz < fNz; iz++)
            {
                if (!mapFile.good())
                {
                    R3BLOG(error, "I/O Error at " << ix << " " << iy << " " << iz);
                }
                index = ix * fNy * fNz + iy * fNz + iz;
                modul = div(index, iDiv);
                if (modul.rem == 0)
                {
                    Double_t perc = TMath::Nint(100. * index / nTot);
                    cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
                }
                mapFile >> ax >> ay >> az >> bx >> by >> bz;

                // Implement mirroring ------------
                Double_t x = ax;
                Double_t y = ay;
                Double_t z = az;

                Int_t index1 = Int_t((x - fXmin) / fXstep) * fNy * fNz + Int_t((y - fYmin) / fYstep) * fNz +
                               Int_t((z - fZmin) / fZstep);

                TVector3 B(bx, by, bz);

                fBx->AddAt(B.X(), index1);
                fBy->AddAt(B.Y(), index1);
                fBz->AddAt(B.Z(), index1);
                // ----------------------------------

                //  cout << "-I- " << bx << " : " << by << " : "  << bz  << " : " << endl;
                if (mapFile.eof())
                {
                    R3BLOG(error,
                           " EOF"
                               << " reached at " << ix << " " << iy << " " << iz);
                    mapFile.close();
                    break;
                }
            } // z-Loop
        }     // y-Loop
    }         // x-Loop

    R3BLOG(info, "   " << index + 1 << " read");
    mapFile.close();
}

// -----   Read field map from ROOT file (private)   ---------------------
//
void R3BGladFieldMap::ReadRootFile(const TString& fileName)
{
    // Opening root file
    R3BLOG(info, "Reading field map from ROOT file " << fileName.Data());
    fFile = new TFile(fileName, "READ");
    if (!(fFile->IsOpen()))
    {
        R3BLOG(fatal, "Cannot read from file!");
    }

    TTree* fTreeMap = NULL;
    fTreeMap = dynamic_cast<TTree*>(fFile->Get("tree"));

    if (!fTreeMap)
    {
        R3BLOG(error, "No TTree named 'tree' found in the file " << fileName.Data());
        R3BLOG(fatal, "No field map data");
    }

    Double_t tBx, tBy, tBz; // branches
    fTreeMap->SetBranchAddress("Bx", &tBx);
    fTreeMap->SetBranchAddress("By", &tBy);
    fTreeMap->SetBranchAddress("Bz", &tBz);

    // Read map info from the first three entries in the tree
    fTreeMap->GetEntry(0);
    fXmin = tBx;
    fXmax = tBy;
    fNx = tBz;

    fTreeMap->GetEntry(1);
    fYmin = tBx;
    fYmax = tBy;
    fNy = tBz;

    fTreeMap->GetEntry(2);
    fZmin = tBx;
    fZmax = tBy;
    fNz = tBz;

    fXstep = (fXmax - fXmin) / Double_t(fNx - 1);
    fYstep = (fYmax - fYmin) / Double_t(fNy - 1);
    fZstep = (fZmax - fZmin) / Double_t(fNz - 1);

    // Create field arrays
    fBx = new TArrayD(fNx * fNy * fNz);
    fBy = new TArrayD(fNx * fNy * fNz);
    fBz = new TArrayD(fNx * fNy * fNz);

    Long64_t Nentries = fTreeMap->GetEntries();

    R3BLOG(info, "Reading GLAD field data from root tree");

    TVector3 fBvec;

    for (Long64_t ev = 3; ev < Nentries; ev++) // first 3 entries are map info
    {
        if (ev % 100000 == 0)
            cout << "\rProcessed " << ev << " entries..." << flush;

        fTreeMap->GetEntry(ev);

        fBvec.SetXYZ(tBx, tBy, tBz);

        fBx->AddAt(fBvec.X(), ev - 3);
        fBy->AddAt(fBvec.Y(), ev - 3);
        fBz->AddAt(fBvec.Z(), ev - 3);
    }
    R3BLOG(info, "Finished reading root tree");

    return;
}

// ------------   Interpolation in a grid cell (private)  -----------------
Double_t R3BGladFieldMap::Interpolate(Double_t dx, Double_t dy, Double_t dz)
{

    // Interpolate in x coordinate
    fHb[0][0] = fHa[0][0][0] + (fHa[1][0][0] - fHa[0][0][0]) * dx;
    fHb[1][0] = fHa[0][1][0] + (fHa[1][1][0] - fHa[0][1][0]) * dx;
    fHb[0][1] = fHa[0][0][1] + (fHa[1][0][1] - fHa[0][0][1]) * dx;
    fHb[1][1] = fHa[0][1][1] + (fHa[1][1][1] - fHa[0][1][1]) * dx;

    // Interpolate in y coordinate
    fHc[0] = fHb[0][0] + (fHb[1][0] - fHb[0][0]) * dy;
    fHc[1] = fHb[0][1] + (fHb[1][1] - fHb[0][1]) * dy;

    // Interpolate in z coordinate
    return fHc[0] + (fHc[1] - fHc[0]) * dz;
}
// ------------------------------------------------------------------------

ClassImp(R3BGladFieldMap)
