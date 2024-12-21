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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Includes from ROOT
#include "TArrayF.h"
#include "TArrayI.h"
#include "TFile.h"
#include "TMath.h"
#include <assert.h>

#include "FairLogger.h"

#include "R3BAladinFieldMap.h"

// Local Macros
#define SQR(x) ((x) * (x))

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::right;
using std::setw;
using std::showpoint;
using TMath::Nint;

// Static local init

/*
 static const TVector3 af_box[2][2] = { { {  123.300, 0.00,  -10.    },
 { -123.224, 0.00,  -10.    }  },
 { {  123.300, 0.00,   10.    },
 { -123.197, 0.00,   10.    }  } };
 static const TVector3 af_mag[2][2] = { { {  123.319, -.261,  -9.408 },
 { -123.022, -.261,   0.092 }  },
 { {  123.259, -.210,  10.411 },
 { -123.022, -.210,   0.092 }  } };

 */

map_fields_ALADiN R3BAladinFieldMap::gMapIFieldOrig;
map_fields_ALADiN R3BAladinFieldMap::gMapIField;

coords_ALADiN R3BAladinFieldMap::gCoords[2];
Bool_t R3BAladinFieldMap::gInitialized = kFALSE;

R3BAladinFieldMap::R3BAladinFieldMap()
    : fBx(NULL)
    , fBy(NULL)
    , fBz(NULL)
{
    fType = 1;
}

R3BAladinFieldMap::R3BAladinFieldMap(const char* mapName, const char* fileType)
    : FairField(mapName)
    , fBx(NULL)
    , fBy(NULL)
    , fBz(NULL)
{
    fType = 1;
}

R3BAladinFieldMap::R3BAladinFieldMap(R3BFieldPar* fieldPar)
    : fBx(NULL)
    , fBy(NULL)
    , fBz(NULL)
{
    fType = 1;
    fCurrent = fieldPar->GetCurrent();
    fScale = fieldPar->GetScale();
}

// ------------   Destructor   --------------------------------------------
R3BAladinFieldMap::~R3BAladinFieldMap()
{
    if (fBx)
        delete fBx;
    if (fBy)
        delete fBy;
    if (fBz)
        delete fBz;
}
// -----------   Intialisation   ------------------------------------------
void R3BAladinFieldMap::Init()
{

    if (gInitialized)
        return;

    gFringeField = kTRUE;

    af_box[0][0].SetXYZ(123.300, 0.00, -10.0);
    af_box[0][1].SetXYZ(-123.224, 0.00, -10.0);
    af_box[1][0].SetXYZ(123.300, 0.00, 10.0);
    af_box[1][1].SetXYZ(-123.197, 0.00, 10.0);

    af_mag[0][0].SetXYZ(123.319, -.261, -9.408);
    af_mag[0][1].SetXYZ(-123.022, -.261, 0.092);
    af_mag[1][0].SetXYZ(123.259, -.210, 10.411);
    af_mag[1][1].SetXYZ(-123.022, -.210, 0.092);

    memcpy((void*)gCoords[0].fBox_pt, (void*)af_box[0], sizeof(af_box[0]));
    memcpy((void*)gCoords[1].fBox_pt, (void*)af_box[1], sizeof(af_box[0]));
    memcpy((void*)gCoords[0].fMag_pt, (void*)af_mag[0], sizeof(af_mag[0]));
    memcpy((void*)gCoords[1].fMag_pt, (void*)af_mag[1], sizeof(af_mag[0]));

    char str[500];
    for (Int_t i = 0; i < 2; i++)
    {
        // The y's are aligned
        assert(gCoords[i].fBox_pt[0].Y() == gCoords[i].fBox_pt[1].Y());
        assert(gCoords[i].fMag_pt[0].Y() == gCoords[i].fMag_pt[1].Y());
        // Calculate the angle between the two coordinate systems
        TVector3 v1, v2;

        v1 = gCoords[i].fBox_pt[1] - gCoords[i].fBox_pt[0];
        v2 = gCoords[i].fMag_pt[1] - gCoords[i].fMag_pt[0];

        Double_t cos_angle = v1.Dot(v2) / TMath::Sqrt(v1.Dot(v1) * v2.Dot(v2));

        gCoords[i].fCosa = cos_angle;
        gCoords[i].fSina = TMath::Sqrt(1 - SQR(cos_angle)) * (i ? -1 : 1);

        sprintf(str, "R3BAladinFieldMap ---> %d: cos %7.4f sin %7.4f", i, gCoords[i].fCosa, gCoords[i].fSina);

        LOG(info) << str;
    }

    // Read the field maps
    // Measured Ic(A) points
    // int measured_I[] = { 500, 1100, 1300, 1500, 1700, 1900, 2100, 2300, 2500 };

    TArrayI measured_I;
    measured_I.Set(9);
    measured_I.AddAt(500, 0);
    measured_I.AddAt(1100, 1);
    measured_I.AddAt(1300, 2);
    measured_I.AddAt(1500, 3);
    measured_I.AddAt(1700, 4);
    measured_I.AddAt(1900, 5);
    measured_I.AddAt(2100, 6);
    measured_I.AddAt(2300, 7);
    measured_I.AddAt(2500, 8);

    for (unsigned int i = 0; i < measured_I.GetSize(); i++)
    {
        FILE* fin;
        Char_t filename[256];
        sprintf(filename, "ala_%04d.dat", measured_I[i]);
        TString mapName(filename);
        TString dir = getenv("VMCWORKDIR");
        TString fMapFileName = dir + "/field/magField/Aladin/newmap/" + mapName;
        LOG(info) << "R3BAladinFieldMap opening Field Map file : " << fMapFileName;

        fin = fopen(fMapFileName, "r");

        if (!fin)
        {
            LOG(error) << "Failure opening field map : " << fMapFileName;
        }

        // free(fMapFileName.Data());

        fields_ALADiN* field = new fields_ALADiN;

        for (int j = 0; j < 3; j++)
        {
            field->f[0][j]._np[0] = 90; //  1 .. 90
            field->f[0][j]._np[1] = 19; // -2 .. 16
            field->f[0][j]._np[2] = 21; //  0 .. 20
            field->f[0][j].init();
            field->f[1][j]._np[0] = 85; //  1 .. 85
            field->f[1][j]._np[1] = 19; // -2 .. 16
            field->f[1][j]._np[2] = 17; //  1 .. 17
            field->f[1][j].init();
        }

        static const int off[2][3] = { { 1, -2, 0 }, { 1, -2, 1 } };

        for (int line = 0; !feof(fin); line++)
        {
            int I, rl, ixyz[3];
            float bxyz[3], bdummy;

            int n = fscanf(fin,
                           " %d %d %d %d %d %f %f %f %f\n",
                           &I,
                           &rl,
                           &ixyz[2],
                           &ixyz[1],
                           &ixyz[0],
                           &bxyz[0],
                           &bxyz[1],
                           &bxyz[2],
                           &bdummy);

            if (n != 9)
                LOG(error) << "Failure parsing field from map: " << filename << " @ line: % " << line;

            if (I != measured_I[i])
                LOG(error) << "Wrong current " << I << " when parsing field from map " << filename
                           << " @ line: " << line;

            if (rl != 0 && rl != 1)
                LOG(error) << "Wrong box " << rl << " when parsing field from map " << filename << " @line: " << line;

            for (int j = 0; j < 3; j++)
            {
                ixyz[j] -= off[rl][j];

                if (ixyz[j] < 0 || ixyz[j] >= field->f[rl][0]._np[j])
                {
                    sprintf(str,
                            "Wrong coordinate(%d) (%d -> %d) >= %d when parsing field from map %s, line? %d.",
                            j,
                            ixyz[j] + off[rl][j],
                            ixyz[j],
                            field->f[rl][0]._np[j],
                            filename,
                            line);
                    LOG(error) << str;
                }
            }

            for (int j = 0; j < 3; j++)
            {
                field->f[rl][j].set_data_pt(ixyz[0], ixyz[1], ixyz[2], bxyz[j]);
            }
        }

        fclose(fin);

        for (int j = 0; j < 3; j++)
        {
            for (int rl = 0; rl < 2; rl++)
            {
                while (field->f[rl][j].expand())
                {
                    // expand as long as there are nodes to be expanded
                }
            }
        }

        gMapIFieldOrig.insert(map_fields_ALADiN::value_type(measured_I[i], field));

        LOG(info) << "R3BAladinFieldMap: Reading field map: " << filename;
    }

    // Has to ben changed somehow using
    // parameters here
    // <DB> check me !!
    Double_t DistanceToTarget = 350.0; // cm
    // Double_t Correction = -117.5; // cm
    // Double_t Correction = -350.5; // cm
    Double_t Correction = -95.0; // cm
    // Double_t Glad_angle = +7.3; // degree
    Double_t Glad_angle = -7.0; // degree
    Double_t DistanceFromtargetToAladinCenter = DistanceToTarget + Correction;
    // Transformations inverse
    gRot = new TRotation();
    gRot->RotateY(-1. * Glad_angle);
    gTrans = new TVector3(0.0, 0.0, -1. * DistanceFromtargetToAladinCenter);

    LOG(info) << "R3BAladinFieldMap::Init() called";
    InitField();

    // Init has been called
    gInitialized = kTRUE;
}

void R3BAladinFieldMap::CalcFieldDiv(R3BFieldInterp f[3], double d[3])
{
    // Assume that the measurement values are at the same locations (not
    // completely true)

    int i[3];
    char str[100];

    for (i[0] = 0; i[0] < f[0]._np[0] - 1; i[0]++)
    {
        // printf ("-- %d --\n",i[0]);
        for (i[1] = 0; i[1] < f[0]._np[1] - 1; i[1]++)
        {
            for (i[2] = 0; i[2] < f[0]._np[2] - 1; i[2]++)
            {
                // divergence is dFx/dx + dFy/dy + dFz/dz

                double divF = 0;
                double divF_var = 0;

                for (int c = 0; c < 3; c++)
                {
                    int ai[3];

                    int abc[3][2] = {
                        { 1, 2 },
                        { 2, 0 },
                        { 0, 1 },
                    };

                    int a = abc[c][0];
                    int b = abc[c][1];

                    double sum_dFa_da = 0, sum_dFa_da_sqr = 0;

                    for (ai[a] = 0; ai[a] < 2; ai[a]++)
                        for (ai[b] = 0; ai[b] < 2; ai[b]++)
                        {
                            ai[c] = 0;
                            double v1 = f[c].get_data_pt(i[0] + ai[0], i[1] + ai[1], i[2] + ai[2]);
                            ai[c] = 1;
                            double v2 = f[c].get_data_pt(i[0] + ai[0], i[1] + ai[1], i[2] + ai[2]);

                            double dFa = v2 - v1;
                            double da = d[c];

                            double dFa_da = dFa / da;

                            sum_dFa_da += dFa_da;
                            sum_dFa_da_sqr += dFa_da * dFa_da;
                        }

                    double dFa_da = sum_dFa_da / 4;
                    double dFa_da_var = (sum_dFa_da_sqr - sum_dFa_da * sum_dFa_da / 4) / 3;
                    ;

                    divF += dFa_da;
                    divF_var += dFa_da_var;
                }
                if (divF == 0)
                {
                    LOG(debug) << "  0   ";
                }
                else
                {
                    sprintf(str, "%5.2f ", divF / sqrt(divF_var));
                    LOG(debug) << str;
                }
            }
            LOG(debug);
        }
    }
}

void R3BAladinFieldMap::InitField()
{
    fCurField = NULL;

    map_fields_ALADiN::iterator iter;

    fFieldSign = fCurrent >= 0 ? 1 : -1;
    double current = TMath::Abs(fCurrent);

    // First look in the original maps

    iter = gMapIFieldOrig.find(current);

    if (iter != gMapIFieldOrig.end())
    {
        fCurField = iter->second;
        return;
    }

    // Then try any interpolation we may already have done

    iter = gMapIField.find(current);

    if (iter != gMapIField.end())
    {
        fCurField = iter->second;
        return;
    }

    // If no interpolation found, then we need to create a new one

    map_fields_ALADiN::iterator iter1;
    map_fields_ALADiN::iterator iter2;

    iter2 = gMapIFieldOrig.lower_bound(current);

    if (iter2 == gMapIFieldOrig.end())
    {
        // Too high current requested...  Extrapolate
        --iter2;
    }

    iter1 = iter2;

    if (iter2 == gMapIFieldOrig.begin())
    {
        LOG(error) << "R3BAladinFieldMap: Cannot interpolate ALADiN field for current " << current << " A, too low.";
    }
    --iter1;

    // printf ("# Interpolate ALADiN %7.1f A (%7.1f A, %7.1f A)\n",
    //        current,iter1->first,iter2->first);

    fCurField = new fields_ALADiN;

    if (!fCurField)
        LOG(debug) << "R3BAladinFieldMap: ALADiN field interpolation, memory allocation failure.";

    double w2 = (current - iter1->first) / (iter2->first - iter1->first);
    double w1 = 1 - w2;

    for (int rl = 0; rl < 2; rl++)
        for (int j = 0; j < 3; j++)
        {
            fCurField->f[rl][j].interpolate(iter1->second->f[rl][j], w1, iter2->second->f[rl][j], w2);
        }

    gMapIField.insert(map_fields_ALADiN::value_type(current, fCurField));

    // cout << "-I- R3BAladinFieldMap Init Field called ********* " << endl;
}

void R3BAladinFieldMap::GetFieldValue(const Double_t point[3], Double_t* bField)
{

    //<D.B> this is a modified version to work
    // in GEANT simulation environment

    TVector3 B(0.0, 0.0, 0.0);
    TVector3 p(point[0], point[1], point[2]);

    // first we need to decompose the position into the magnet coordinates.
    // printf ("Vectorized  p: %10.5f %10.5f %10.5f : ",p.X(),p.Y(),p.Z());

    Double_t DistanceToTarget = 350.0; // cm
    //  Double_t Correction = -114.5; // cm
    Double_t Correction = -95.0; // cm
    //  Double_t Magnet_angle = +7.3; // degree
    Double_t Magnet_angle = -7.0; // degree

    Double_t xx = 0.;
    Double_t yy = 0.;
    Double_t zz = 0.;

    Double_t angle = -1. * Magnet_angle * TMath::Pi() / 180.;

    // Translation
    Double_t zt = point[2] - (DistanceToTarget + Correction);

    // Rotation
    xx = zt * TMath::Sin(angle) + point[0] * TMath::Cos(angle);
    yy = point[1];
    zz = zt * TMath::Cos(angle) - point[0] * TMath::Sin(angle);

    // Assign to c
    TVector3 c(xx, yy, zz);

    // printf ("c: %10.5f %10.5f %10.5f\n",c.X(),c.Y(),c.Z());
    // c is in the magnet coordinate system.

    double z_abs = TMath::Abs(c.Z());

    if ((gFringeField == kFALSE) && (z_abs > (160.0 / 1.)))
    {
        bField[0] = 0.0;
        bField[1] = 0.0;
        bField[2] = 0.0;
        return; // the fringe field is cut back/forward
    }

    // for each field map, there is one mar transformation to do
    // to get into each field map's coordinate system.

    // We want the field to be continous also outside the field maps.
    // So the maps will give their boundary values when we happen to
    // reach outside.  In Y and X directions (i.e. up and to the sides)
    // the means to give non-zero values.  In the Z direction, this
    // requires the field maps to vanish.

    // One issue left is when switching maps from left to right.  In
    // region where both maps are valid (note: this differs on the x y
    // and z compononets) we can take an average.  To not get a sharp
    // discontinuity at the edge where the maps end, we'll do the
    // averaging with a decreasing function, such that it does not
    // contribute to the average at all 1 unit away from the end of the
    // field map.

    double Bi[3];
    double Bp[3];
    double wsum = 0.0;

    for (int i = 0; i < 3; i++)
        Bi[i] = 0.0;

    for (int rl = 0; rl < 2; rl++)
    {
        // coords_ALADiN  gCoords;

        // swap axis
        double magx = c.Z();
        double magy = c.Y();
        double magz = -c.X();

        // move to origin at pt in magnet coord system
        double boxxp = magx - gCoords[rl].fMag_pt[0].X();
        double boxy = magy - gCoords[rl].fMag_pt[0].Y();
        double boxzp = magz - gCoords[rl].fMag_pt[0].Z();

        // rotate into nbox coordinate system
        double boxx = gCoords[rl].fCosa * boxxp - gCoords[rl].fSina * boxzp;
        double boxz = gCoords[rl].fSina * boxxp + gCoords[rl].fCosa * boxzp;

        // mova away from rotation origin to real box origin
        boxx += gCoords[rl].fBox_pt[0].X();
        boxy += gCoords[rl].fBox_pt[0].Y();
        boxz += gCoords[rl].fBox_pt[0].Z();

        // printf ("rl %d: %10.5f %10.5f %10.5f : ",rl,
        //        boxx,boxy,boxz);

        // So box coordinates are now fine, need to calculate our
        // interpolation coordinates

        static const double box_grid_off[2][5] = {
            { 140. - 1 * 3, 139. - 1 * 3, 138. - 1 * 3, 19. + 2 * 3, 66. },
            { 125. - 1 * 3, 124. - 1 * 3, 123. - 1 * 3, 19. + 2 * 3, 6. - 1 * 4 }
        };

        double dc[3];
        int ic[3];

        dc[1] = (boxy + box_grid_off[rl][3]) * (1 / 3.);
        dc[2] = (boxz + box_grid_off[rl][4]) * (1 / 4.);

        ic[1] = (int)floor(dc[1]);
        dc[1] -= ic[1];

        ic[2] = (int)floor(dc[2]);
        dc[2] -= ic[2];
        /*
         printf ("(1:%2d/%7.5f 2:%2d/%7.5f ",
         ic[1],dc[1],ic[2],dc[2]);
         */
        double w = 1;

        if (rl == 0)
        {
            if (ic[2] > 21)
            {
                continue; // this is outside range
            }
            else if (ic[2] > 20)
            {
                w = 1 - dc[2];
            }
        }
        else
        {
            if (ic[2] < -1)
            {
                continue; // this is outside range
            }
            else if (ic[2] < 0)
            {
                w = dc[2];
            }
        }

        wsum += w;

        double Bbi[3];

        for (int i = 0; i < 3; i++)
        {
            dc[0] = (boxx + box_grid_off[rl][i]) * (1 / 3.);

            ic[0] = (int)floor(dc[0]);
            dc[0] -= ic[0];

            // printf (" -I- Interpolation parameters ---->> 0:%2d/%7.5f",
            //         ic[0],dc[0]);

            Bbi[i] = fCurField->f[rl][i].interp(ic, dc);
        }

        // The field must be rotated into the ALADiN FRAME (from the box
        // frame)

        Bi[0] += w * (gCoords[rl].fCosa * Bbi[0] + gCoords[rl].fSina * Bbi[2]);
        Bi[1] += w * Bbi[1];
        Bi[2] += w * (gCoords[rl].fSina * Bbi[0] - gCoords[rl].fCosa * Bbi[2]); // * -1
        /*
         printf (")\n");

         printf ("%8.5f %8.5f %8.5f : ",Bbi[0],Bbi[1],Bbi[2]);
         printf ("%8.5f %8.5f %8.5f\n",Bi[0],Bi[1],Bi[2]);
         */
        Bp[0] = -gCoords[rl].fCosa * Bbi[2] + gCoords[rl].fSina * Bbi[0];
        Bp[1] = Bbi[1];
        Bp[2] = -gCoords[rl].fSina * Bbi[2] + gCoords[rl].fCosa * Bbi[0];

        continue;
        // printf (")\n");
    }

    double wsuminv = 0.0;
    if (wsum != 0.)
        wsuminv = fFieldSign / wsum;

    if ((z_abs > 120.0) && (z_abs <= 160.))
    {
        // Between 150 and 160 cm, we force the field down to 0
        double factor = 0.025 * (160. - z_abs);
        wsuminv *= factor;
    }
    if (z_abs > 160.)
        wsuminv = 0.0;

    for (int i = 0; i < 3; i++)
        Bi[i] *= wsuminv;

    // And then transformed into world coordinates
    // rotate_field(B,Bi[0],Bi[1],Bi[2]);
    TVector3 BLab(Bi[0], Bi[1], Bi[2]);
    //   BLab.Transform(*gRot);
    // copy value @ end /old/standard version
    // bField[0] = -1.*BLab.X()*10.; // [kGauss]
    // bField[1] = -1.*BLab.Y()*10.; // [kGauss]
    // bField[2] = -1.*BLab.Z()*10.; // [kGauss]
    // swap field like in tracker, according to denis and justyna
    // git/justyna version:
    bField[2] = -1. * BLab.X() * 10.; // [kGauss]
    bField[1] = -1. * BLab.Y() * 10.; // [kGauss]
    bField[0] = +1. * BLab.Z() * 10.; // [kGauss]

    //  cout << "-I- ALADIN --> X: " << p.X() << " Y: " << p.Y() << "Z:" << p.Z() << endl;
    //  cout << "-I- ALADIN --> Bx: " << bField[0] << " By: " << bField[1] << "Bz:" << bField[2] << endl;

    return;
}

//-------------------------- Standard R3BROOT API ----------------------//
Double_t R3BAladinFieldMap::GetBx(Double_t x, Double_t y, Double_t z) { return 0.; }
Double_t R3BAladinFieldMap::GetBy(Double_t x, Double_t y, Double_t z) { return 0.; }

Double_t R3BAladinFieldMap::GetBz(Double_t x, Double_t y, Double_t z) { return 0.; }
Bool_t R3BAladinFieldMap::IsInside(Double_t x,
                                   Double_t y,
                                   Double_t z,
                                   Int_t& ix,
                                   Int_t& iy,
                                   Int_t& iz,
                                   Double_t& dx,
                                   Double_t& dy,
                                   Double_t& dz)
{
    return kFALSE;
}
void R3BAladinFieldMap::WriteAsciiFile(const char* fileName) {}

void R3BAladinFieldMap::SetPosition(Double_t x, Double_t y, Double_t z)
{
    fPosX = x;
    fPosY = y;
    fPosZ = z;
}

void R3BAladinFieldMap::Print(Option_t* option) const {}

void R3BAladinFieldMap::Reset() {}

void R3BAladinFieldMap::ReadAsciiFile(const char* fileName) {}

Double_t R3BAladinFieldMap::Interpolate(Double_t dx, Double_t dy, Double_t dz) { return 0.; }

ClassImp(R3BAladinFieldMap);
