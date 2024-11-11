/*   Copyright (C) 2019 Members of R3B Collaboration                          *
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
#include "TArrayD.h"
#include "TArrayF.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"

#include "R3BGladFieldMap.h"

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
    fTrackerCorr = 1.;
}
// ------------------------------------------------------------------------

// -------------   Standard constructor   ---------------------------------
R3BGladFieldMap::R3BGladFieldMap(const char* mapName, const char* fileType)
    : FairField(mapName)
{
    fPosX = fPosY = fPosZ = 0.;
    fXmin = fYmin = fZmin = 0.;
    fXmax = fYmax = fZmax = 0.;
    fXstep = fYstep = fZstep = 0.;
    fNx = fNy = fNz = 0;
    fScale = 1.;
    fBx = fBy = fBz = NULL;
    fName = mapName;
    TString dir = getenv("VMCWORKDIR");
    fFileName = dir + "/field/magField/R3B/" + mapName;
    if (fileType[0] == 'R')
        fFileName += ".root";
    else
        fFileName += ".dat";

    fType = 2;
    fTrackerCorr = 1.;
}

// ------------   Constructor from R3BGladFieldPar   --------------------------

R3BGladFieldMap::R3BGladFieldMap(R3BFieldPar* fieldPar)
{
    fType = 2;
    fPosX = fPosY = fPosZ = 0.;
    fXmin = fYmin = fZmin = 0.;
    fXmax = fYmax = fZmax = 0.;
    fXstep = fYstep = fZstep = 0.;
    fNx = fNy = fNz = 0;
    fScale = 1.;
    fBx = fBy = fBz = NULL;
    if (!fieldPar)
    {
        cerr << "-W- R3BGladFieldConst::R3BGladFieldMap: empty parameter container!" << endl;
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
        fScale = fieldPar->GetScale();
        // TString dir = getenv("VMCWORKDIR");
        // fFileName = dir + "/field/magField/R3B/" + fName;
        // fFileName += ".dat";
        fType = fieldPar->GetType();
    }
    fTrackerCorr = 1.;
}

// ------------   Destructor   --------------------------------------------
R3BGladFieldMap::~R3BGladFieldMap()
{
    if (fBx)
        delete fBx;
    if (fBy)
        delete fBy;
    if (fBz)
        delete fBz;
}
// ------------------------------------------------------------------------

// -----------   Intialisation   ------------------------------------------
void R3BGladFieldMap::Init()
{

    Bool_t newField = true;

    // Distance in z: turning point to flanch:   539.5 mm
    //				       GLAD flanch target:  1210 mm
    //										    1749.5 mm
    //

    if (newField)
    {
        // 2022
        /*
        fPosX = 0.86128;
        fPosY = -1.543865 + 1.5;
        fPosZ = 174.5681; // from the fit of the magnetic field measurements
        fXAngle = -0.1446;
        fYAngle = -14.19588;
        fZAngle = -0.7113;
        */
        // 2022+2023
        /*fPosX = 0.2443951;
        fPosY = -0.97518 + 1.5;
        fPosZ = 174.5638;
        fXAngle = 0.1388603;
        fYAngle = -14.20966;
        fZAngle = -0.1381274;*/

        fPosX = 0.4942526;
        fPosY = -1.022616 + 1.5;
        fPosZ = 174.5851;
        fXAngle = 0.04190000; // fixed to the measured value by Michael
        fYAngle = -14.18602;
        fZAngle = -0.1464000; // fixed to the measured value by Michael
    }

    if (!newField)
    {
        fPosX = 0.0;
        fPosY = 1.5;
        fPosZ = 174.95; // from the drawings
        fXAngle = 0.0;
        fYAngle = -14.;
        fZAngle = 0.0;
    }

    gTrans = new TVector3(-fPosX, -fPosY, -fPosZ);
    if (fFileName.EndsWith(".dat"))
        ReadAsciiFile(fFileName);
    else if (fFileName.EndsWith(".root"))
        ReadRootFile(fFileName);
    else
    {
        cerr << "-E- R3BGladFieldMap::Init: No proper file name defined! (" << fFileName << ")" << endl;
        LOG(fatal) << "Init: No proper file name";
    }
}
// ------------------------------------------------------------------------

// -----------   Get x component of the field   ---------------------------
Double_t R3BGladFieldMap::GetBx(Double_t x, Double_t y, Double_t z)
{

    // cout << "-I- get Bx called " << endl;
    // transform to local coordinates
    // x, y, z are global coordinates
    TVector3 localPoint(x, y, z);

    // gTrans = new TVector3(-fPosX, -fPosY, -fPosZ);
    // localPoint = localPoint + (*gTrans);

    // fPosX, fPosY, fPosZ origin of GLAD system
    TVector3 gTrans1;
    gTrans1.SetX(-fPosX);
    gTrans1.SetY(-fPosY);
    gTrans1.SetZ(-fPosZ);
    localPoint = localPoint + gTrans1;

    localPoint.RotateX(-fXAngle * TMath::DegToRad());
    localPoint.RotateY(-fYAngle * TMath::DegToRad());
    localPoint.RotateZ(-fZAngle * TMath::DegToRad());

    Int_t ix = 0;
    Int_t iy = 0;
    Int_t iz = 0;
    Double_t dx = 0.;
    Double_t dy = 0.;
    Double_t dz = 0.;

    if (IsInside(localPoint.X(), localPoint.Y(), localPoint.Z(), ix, iy, iz, dx, dy, dz))
    {
        // Get Bx field values at grid cell corners
        fHa[0][0][0] = fBx->At(ix * fNy * fNz + iy * fNz + iz);
        fHa[1][0][0] = fBx->At((ix + 1) * fNy * fNz + iy * fNz + iz);
        fHa[0][1][0] = fBx->At(ix * fNy * fNz + (iy + 1) * fNz + iz);
        fHa[1][1][0] = fBx->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + iz);
        fHa[0][0][1] = fBx->At(ix * fNy * fNz + iy * fNz + (iz + 1));
        fHa[1][0][1] = fBx->At((ix + 1) * fNy * fNz + iy * fNz + (iz + 1));
        fHa[0][1][1] = fBx->At(ix * fNy * fNz + (iy + 1) * fNz + (iz + 1));
        fHa[1][1][1] = fBx->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + (iz + 1));

        // Return interpolated field value
        Double_t val = Interpolate(dx, dy, dz);
        // cout << " (X) interpolated " << val << endl;
        //        val = val + 0.5;
        // val = val*0.8;

        // cout<<"xx: "<<localPoint.X()<<", yy: "<<localPoint.Y()<<", zz: "<<localPoint.Z()<<", Bx: "<<fTrackerCorr *
        // val<<endl;

        return (fTrackerCorr * val);
    }
    return 0;
}
// ------------------------------------------------------------------------

// -----------   Get y component of the field   ---------------------------
Double_t R3BGladFieldMap::GetBy(Double_t x, Double_t y, Double_t z)
{

    // cout << "-I- get By called " << endl;
    // transform to local coordinates
    // local to global
    TVector3 localPoint(x, y, z);

    // gTrans = new TVector3(-fPosX, -fPosY, -fPosZ);
    // localPoint = localPoint + (*gTrans);

    TVector3 gTrans1;
    gTrans1.SetX(-fPosX);
    gTrans1.SetY(-fPosY);
    gTrans1.SetZ(-fPosZ);
    localPoint = localPoint + gTrans1;
    localPoint.RotateX(-fXAngle * TMath::DegToRad());
    localPoint.RotateY(-fYAngle * TMath::DegToRad());
    localPoint.RotateZ(-fZAngle * TMath::DegToRad());

    Int_t ix = 0;
    Int_t iy = 0;
    Int_t iz = 0;
    Double_t dx = 0.;
    Double_t dy = 0.;
    Double_t dz = 0.;

    if (IsInside(localPoint.X(), localPoint.Y(), localPoint.Z(), ix, iy, iz, dx, dy, dz))
    {
        // Get By field values at grid cell corners
        fHa[0][0][0] = fBy->At(ix * fNy * fNz + iy * fNz + iz);
        fHa[1][0][0] = fBy->At((ix + 1) * fNy * fNz + iy * fNz + iz);
        fHa[0][1][0] = fBy->At(ix * fNy * fNz + (iy + 1) * fNz + iz);
        fHa[1][1][0] = fBy->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + iz);
        fHa[0][0][1] = fBy->At(ix * fNy * fNz + iy * fNz + (iz + 1));
        fHa[1][0][1] = fBy->At((ix + 1) * fNy * fNz + iy * fNz + (iz + 1));
        fHa[0][1][1] = fBy->At(ix * fNy * fNz + (iy + 1) * fNz + (iz + 1));
        fHa[1][1][1] = fBy->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + (iz + 1));

        // Return interpolated field value
        Double_t val = Interpolate(dx, dy, dz);
        // cout << " (Y) interpolated " << val << endl;
        // if(localPoint.Z() > 300.) val = val * 1.2;
        // if(localPoint.Z() > 200. && localPoint.Z() < 300.) val = val *1.02;

        /*

                Double_t fcor = 1.;
                if(localPoint.Z()>36.67 && localPoint.Z()<294.){
        //			if(localPoint.Y() == 0. ){  // correction for y=0 plane
                        Double_t
        par[8]={1.0157655,-5.491E-04,8.7452E-06,-8.849E-08,5.6868E-10,-2.001E-12,2.8387E-15,-3.150E-19}; Double_t zz =
        localPoint.Z();

                        fcor = par[0] + par[1]* TMath::Power(zz,1)
                                      + par[2]* TMath::Power(zz,2)
                                      + par[3]* TMath::Power(zz,3)
                                      + par[4]* TMath::Power(zz,4)
                                      + par[5]* TMath::Power(zz,5)
                                      + par[6]* TMath::Power(zz,6)
                                      + par[7]* TMath::Power(zz,7);
        //			}
                }
            //	fcor = 1.;
                val = val*fcor;
        */
        Double_t xx = localPoint.X();
        Double_t yy = localPoint.Y();
        Double_t zz = localPoint.Z();
        /*
          Double_t p0,p1,p2,fitparfit[7];
          p0=-0.00830474-0.0043102*yy-1.149e-4*yy*yy;
          p1=0.00690782+5.522e-4*yy-3.835e-4*yy*yy;
          p2=-0.00142735+2.2243e-5*yy;
          fitparfit[0]=p0+p1*xx+p2*xx*xx;

          p0=-0.01915779+0.00120021;
          p1=-0.00180119+3.8234e-5*yy+7.8526e-6*yy*yy;
          fitparfit[1]=p0+p1*xx;

          p0=7.8322e-4-8.265e-5*yy;
          p1=3.9034e-5-1.169e-6*yy+4.2791e-7*yy*yy;
          fitparfit[2]=p0+p1*xx;

          p0=-2.22e-6+2.3711e-6*yy;
          p1=-6.132e-7+3.6396e-8*yy-9.395e-9*yy*yy;
          fitparfit[3]=p0+p1*xx;

          p0=1.6986e-7;
          p1=-1.236e-8-6.819e-10*yy;
          fitparfit[4]=p0+p1*xx;

          p0=-1.597e-8-2.079e-9*yy-1.827e-10*yy*yy;
          p1=4.9526e-10-4.448e-11*yy;
          fitparfit[5]=p0+p1*xx;

          p0=2.809e-10+3.9839e-11*yy+3.3516e-12*yy*yy;
          p1=6.4899e-13+1.4236e-12*yy;
          fitparfit[6]=p0+p1*xx;

          Double_t fcorr=0.;
          if(zz < 25.)
          {
              for(Int_t ipar=0; ipar<7; ipar++)
              {
                  fcorr=fcorr+fitparfit[ipar]*TMath::Power(zz,ipar);
              }
          }
          //val = val/(1.-fcorr/100.);
          */
        /*
         Double_t fcorr = 0.0;
         if(zz < 5.75 && zz > -44.25)
          {
              zz = zz + 121. +  53.25;
              fcorr=fcorr + 1417.05
                          - 33.2292 * TMath::Power(zz,1)
                          + 0.292215 * TMath::Power(zz,2)
                          - 0.00114258 * TMath::Power(zz,3)
                          + 1.67625e-06* TMath::Power(zz,4);
          }
          val = val/(1.-fcorr/100.);
       */

        gRandom->SetSeed(0);
        /*
          //if(zz > -84.25 && zz < -43.25) //fi23-10cm
         // if(zz < -84.25 && zz > -174.25)  // target and fi23
          //if(zz < -53.25 && zz > -174.25)  // target and entrance
          if(zz < -43.25 && zz > -174.25)  // target and 10cm
          //if(zz > 400.)
          {
                val = val * 1.;
           }
         */
        //   cout<<"xx: "<<xx<<", yy: "<<yy<<", zz: "<<zz<<", By: "<<fTrackerCorr * val<<endl;

        return (fTrackerCorr * val);
    }
    return 0.;
}
// ------------------------------------------------------------------------

// -----------   Get z component of the field   ---------------------------
Double_t R3BGladFieldMap::GetBz(Double_t x, Double_t y, Double_t z)
{

    // cout << "-I- get Bz called " << endl;
    // transform to local coordinates
    // local to global
    TVector3 localPoint(x, y, z);

    // gTrans = new TVector3(-fPosX, -fPosY, -fPosZ);
    // localPoint = localPoint + (*gTrans);

    TVector3 gTrans1;
    gTrans1.SetX(-fPosX);
    gTrans1.SetY(-fPosY);
    gTrans1.SetZ(-fPosZ);
    localPoint = localPoint + gTrans1;

    localPoint.RotateX(-fXAngle * TMath::DegToRad());
    localPoint.RotateY(-fYAngle * TMath::DegToRad());
    localPoint.RotateZ(-fZAngle * TMath::DegToRad());

    Int_t ix = 0;
    Int_t iy = 0;
    Int_t iz = 0;
    Double_t dx = 0.;
    Double_t dy = 0.;
    Double_t dz = 0.;

    if (IsInside(localPoint.X(), localPoint.Y(), localPoint.Z(), ix, iy, iz, dx, dy, dz))
    {
        // Get Bz field values at grid cell corners
        fHa[0][0][0] = fBz->At(ix * fNy * fNz + iy * fNz + iz);
        fHa[1][0][0] = fBz->At((ix + 1) * fNy * fNz + iy * fNz + iz);
        fHa[0][1][0] = fBz->At(ix * fNy * fNz + (iy + 1) * fNz + iz);
        fHa[1][1][0] = fBz->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + iz);
        fHa[0][0][1] = fBz->At(ix * fNy * fNz + iy * fNz + (iz + 1));
        fHa[1][0][1] = fBz->At((ix + 1) * fNy * fNz + iy * fNz + (iz + 1));
        fHa[0][1][1] = fBz->At(ix * fNy * fNz + (iy + 1) * fNz + (iz + 1));
        fHa[1][1][1] = fBz->At((ix + 1) * fNy * fNz + (iy + 1) * fNz + (iz + 1));

        // Return interpolated field value
        Double_t val = Interpolate(dx, dy, dz);
        // cout << " (Z) interpolated " << val << endl;
        //        val = val + 0.5;
        // val = val*1.2;

        //   cout<<"xx: "<<localPoint.X()<<", yy: "<<localPoint.Y()<<", zz: "<<localPoint.Z()<<", Bz: "<<fTrackerCorr *
        //   val<<endl;

        return (fTrackerCorr * val);
    }

    return 0.;
}
// ------------------------------------------------------------------------

// -----------   Check whether a point is inside the map   ----------------
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
// ------------------------------------------------------------------------

// ----------   Write the map to an ASCII file   --------------------------
void R3BGladFieldMap::WriteAsciiFile(const char* fileName)
{

    // Open file
    cout << "-I- R3BGladFieldMap: Writing field map to ASCII file " << fileName << endl;
    ofstream mapFile(fileName);
    if (!mapFile.is_open())
    {
        cerr << "-E- R3BGladFieldMap:WriteAsciiFile: Could not open file! " << endl;
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
    cout << "-I- R3BGladFieldMap: " << fNx * fNy * fNz << " entries to write... " << setw(3) << 0 << " % ";
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
    cout << "   " << index + 1 << " written" << endl;
    mapFile.close();
}
// ------------------------------------------------------------------------

// -------   Write field map to a ROOT file   -----------------------------
/*
   void R3BGladFieldMap::WriteRootFile(const char* fileName,
   const char* mapName) {

   R3BGladFieldMapData* data = new R3BGladFieldMapData(mapName, *this);
   TFile* oldFile = gFile;
   TFile* file = new TFile(fileName, "RECREATE");
   data->Write();
   file->Close();
   if(oldFile) oldFile->cd();

   }
   */

// ------------------------------------------------------------------------

// -----  Set the position of the field centre in global coordinates  -----
void R3BGladFieldMap::SetPosition(Double_t x, Double_t y, Double_t z)
{
    fPosX = x;
    fPosY = y;
    fPosZ = z;
}
// ------------------------------------------------------------------------

// ---------   Screen output   --------------------------------------------
void R3BGladFieldMap::Print(Option_t* option) const
{
    TString type = "Map";
    if (fType == 2)
        type = "Map sym2";
    if (fType == 3)
        type = "Map sym3";
    cout << "======================================================" << endl;
    cout.precision(4);
    cout << showpoint;
    cout << "----  " << fTitle << " : " << fName << endl;
    cout << "----" << endl;
    cout << "----  Field type     : " << type << endl;
    cout << "----" << endl;
    cout << "----  Field map grid : " << endl;
    cout << "----  x = " << setw(4) << fXmin << " to " << setw(4) << fXmax << " cm, " << fNx
         << " grid points, dx = " << fXstep << " cm" << endl;
    cout << "----  y = " << setw(4) << fYmin << " to " << setw(4) << fYmax << " cm, " << fNy
         << " grid points, dy = " << fYstep << " cm" << endl;
    cout << "----  z = " << setw(4) << fZmin << " to " << setw(4) << fZmax << " cm, " << fNz
         << " grid points, dz = " << fZstep << " cm" << endl;
    cout << endl;
    cout << "----  Field centre position: ( " << setw(6) << fPosX << ", " << setw(6) << fPosY << ", " << setw(6)
         << fPosZ << ") cm" << endl;
    cout << "----  Field rotation X: " << setw(6) << fXAngle << " deg" << endl;
    cout << "----  Field rotation Y: " << setw(6) << fYAngle << " deg" << endl;
    cout << "----  Field rotation Z: " << setw(6) << fZAngle << " deg" << endl;
    cout << "----  Field scaling factor: " << fScale << endl;
    //  Double_t bx = GetBx(0.,0.,0.);
    //  Double_t by = GetBy(0.,0.,0.);
    //  Double_t bz = GetBz(0.,0.,0.);
    //  cout << "----" << endl;
    //  cout << "----  Field at origin is ( " << setw(6) << bx << ", " << setw(6)
    //       << by << ", " << setw(6) << bz << ") kG" << endl;
    cout << "======================================================" << endl;
}
// ------------------------------------------------------------------------

// ---------    Reset parameters and data (private)  ----------------------
void R3BGladFieldMap::Reset()
{
    fPosX = fPosY = fPosZ = 0.;
    fXmin = fYmin = fZmin = 0.;
    fXmax = fYmax = fZmax = 0.;
    fXstep = fYstep = fZstep = 0.;
    fNx = fNy = fNz = 0;
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
// ------------------------------------------------------------------------

// -----   Read field map from ASCII file (private)   ---------------------
void R3BGladFieldMap::ReadAsciiFile(const char* fileName)
{

    Double_t bx = 0., by = 0., bz = 0.;
    Double_t ax = 0., ay = 0., az = 0.;
    // Open file
    cout << "-I- R3BGladFieldMap: Reading field map from ASCII file " << fileName << endl;
    ifstream mapFile(fileName);
    if (!mapFile.is_open())
    {
        cerr << "-E- R3BGladFieldMap:ReadAsciiFile: Could not open file! " << endl;
        LOG(fatal) << "ReadAsciiFile: Could not open file";
    }

    // Read map type
    TString type;
    mapFile >> type;
    Int_t iType = 0;
    if (type == "nosym")
        iType = 1;
    if (type == "sym2")
        iType = 2;
    if (type == "sym3")
        iType = 3;
    // if (fType != iType)
    //{
    //    cout << "-E- R3BGladFieldMap::ReadAsciiFile: Incompatible map types!" << endl;
    //    cout << "    Field map is of type " << fType << " but map on file is of type " << iType << endl;
    //    LOG(fatal) << "ReadAsciiFile: Incompatible map types";
    //}

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
    cout << "-I- R3BGladFieldMap: " << nTot << " entries to read... " << setw(3) << 0 << " % ";
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
                    cerr << "-E- R3BGladFieldMap::ReadAsciiFile: "
                         << "I/O Error at " << ix << " " << iy << " " << iz << endl;
                index = ix * fNy * fNz + iy * fNz + iz;
                modul = div(index, iDiv);
                if (modul.rem == 0)
                {
                    Double_t perc = TMath::Nint(100. * index / nTot);
                    cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
                }
                mapFile >> ax >> ay >> az >> bx >> by >> bz;

                // Implement mirroring ----------------------------------------------------------------------
                Double_t x = ax;
                Double_t y = ay;
                Double_t z = az;

                Int_t index1 = Int_t((x - fXmin) / fXstep) * fNy * fNz + Int_t((y - fYmin) / fYstep) * fNz +
                               Int_t((z - fZmin) / fZstep);

                TVector3 B(bx, by, bz);
                B.RotateX(fXAngle * TMath::DegToRad());
                B.RotateY(fYAngle * TMath::DegToRad());
                B.RotateZ(fZAngle * TMath::DegToRad());

                fBx->AddAt(factor * B.X(), index1);
                fBy->AddAt(factor * B.Y(), index1);
                fBz->AddAt(factor * B.Z(), index1);
                // ------------------------------------------------------------------------------------------

                //  cout << "-I- " << bx << " : " << by << " : "  << bz  << " : " << endl;
                if (mapFile.eof())
                {
                    cerr << endl
                         << "-E- R3BGladFieldMap::ReadAsciiFile: EOF"
                         << " reached at " << ix << " " << iy << " " << iz << endl;
                    mapFile.close();
                    break;
                }
            } // z-Loop
        }     // y-Loop0)
    }         // x-Loop

    cout << "   " << index + 1 << " read" << endl;

    mapFile.close();
    //  exit(0);
}
// ------------------------------------------------------------------------

void R3BGladFieldMap::ReadRootFile(const char* fileName)
{
    // Opening root file
    cout << "-I- R3BGladFieldMap: Reading field map from ROOT file " << fileName;
    fFile = new TFile(fileName, "READ");
    if (!(fFile->IsOpen()))
    {
        cerr << "-E- R3BGladFieldMap::ReadRootFile: Cannot read from file! " << endl;
        LOG(fatal) << "ReadRootFile: Cannot read from file";
    }

    TTree* fTreeMap = NULL;
    fTreeMap = (TTree*)fFile->Get("tree");
    // fTreeMap = (TTree*)fFile->Get("tree")->Clone("fTreeMap");

    if (!fTreeMap)
    {
        cerr << "-E- R3BGladFieldMap::ReadRootFile: no TTree named 'tree' found in the file" << fileName;
        LOG(fatal) << "No field map data";
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

    Double_t factor = fScale * 10.; // Factor 10 for T -> kG

    // Create field arrays
    fBx = new TArrayD(fNx * fNy * fNz);
    fBy = new TArrayD(fNx * fNy * fNz);
    fBz = new TArrayD(fNx * fNy * fNz);

    Long64_t Nentries = fTreeMap->GetEntries();

    cout << "\n-I- Reading GLAD field data from root tree" << endl;

    TVector3 fBvec;

    for (Long64_t ev = 3; ev < Nentries; ev++) // first 3 entries aremap info
    {
        if (ev % 100000 == 0)
            cout << "\rProcessed " << ev << " entries..." << flush;

        fTreeMap->GetEntry(ev);

        fBvec.SetXYZ(tBx * factor, tBy * factor, tBz * factor);
        fBvec.RotateX(fXAngle * TMath::DegToRad());
        fBvec.RotateY(fYAngle * TMath::DegToRad());
        fBvec.RotateZ(fZAngle * TMath::DegToRad());

        fBx->AddAt(fBvec.X(), ev - 3);
        fBy->AddAt(fBvec.Y(), ev - 3);
        fBz->AddAt(fBvec.Z(), ev - 3);
    }
    cout << "\n-I- Finished reading root tree" << endl;

    return;
}

/*
// Get the field parameters
SetField(data);

*/

// ------------------------------------------------------------------------

// ------------   Set field parameters and data (private)  ----------------
/*
   void R3BGladFieldMap::SetField(const R3BGladFieldMapData* data) {

// Check compatibility
if ( data->GetType() != fType ) {
if (!((data->GetType()==3)&&(fType==5)))                   // E.Litvinenko
{
cout << "-E- R3BGladFieldMap::SetField: Incompatible map types!"
<< endl;
cout << "    Field map is of type " << fType
<< " but map on file is of type " << data->GetType() << endl;
LOG(fatal) << "SetField: Incompatible map types";
}
else
cout << "   R3BGladFieldMap::SetField: Warning:  You are using PosDepScaled map (original map type = 3)" << endl;
}


*/

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
