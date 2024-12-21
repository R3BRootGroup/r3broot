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
// -----                   R3BTraHitFinder source file                 -----
// -----                  Created 01/06/11  by N.Ashwood               -----
// -------------------------------------------------------------------------
#include "R3BTraHitFinder.h"
#include "TGeoMatrix.h"
#include "TMath.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"
#include "TRandom.h"

#include "TGeoManager.h"

#include "R3BTraPoint.h"
#include "R3BTrackerHit.h"

R3BTraHitFinder::R3BTraHitFinder()
    : FairTask("R3B Tracker Hit Finder ")
{
    fThreshold = 0.;         // no threshold
    fTrackerResolution = 0.; // perfect resolution
}

R3BTraHitFinder::~R3BTraHitFinder() {}

// -----   Public method Init   --------------------------------------------
InitStatus R3BTraHitFinder::Init()
{
    FairRootManager* ioManager = FairRootManager::Instance();
    if (!ioManager)
        LOG(fatal) << "Init: No FairRootManager";
    fTrackerHitCA = dynamic_cast<TClonesArray*>(ioManager->GetObject("TraPoint"));

    // Register output array TraHit
    fTraHitCA = new TClonesArray("R3BTrackerHit", 1000);
    ioManager->Register("TrackerHit", "Tracker Hit", fTraHitCA, kTRUE);

    return kSUCCESS;
}

// -----   Public method ReInit   --------------------------------------------
InitStatus R3BTraHitFinder::ReInit() { return kSUCCESS; }

// -----   Public method Exec   --------------------------------------------
void R3BTraHitFinder::Exec(Option_t* opt)
{

    Reset();

    // Si Geometrical parameter:
    // Inner layer
    Double_t Length1 = 19.03;     // cm
    Double_t WidthMax1 = 7.945;   // cm
    Double_t WidthMin1 = 2.25;    // cm
    Double_t StripPitch1 = 0.005; // = 50 um
    // Double_t AngRangeMin1=7.;
    // Double_t InclAng1=14.9;
    Double_t AngTrap1 = atan((WidthMax1 / 2 - WidthMin1 / 2) / Length1);
    // Double_t WidthHalf1=WidthMax1 - (Length1/cos(AngTrap1))*sin(AngTrap1);
    Double_t StepZ1 = StripPitch1 / sin(AngTrap1); // step along the z axis of the detector (in xz plan)
    Double_t StepX1 = StripPitch1 / cos(AngTrap1); // step along the x axis of the detector (in xz plan)
    Int_t NbStrip1 = int(WidthMax1 / StepX1);      //
    // LOG(info) << "NbStrip1= " << NbStrip1;
    // Double_t ShiftalongZ;
    Double_t Proj;      // cm
    Double_t ProjStrip; // cm
    Int_t strip = 0;
    Int_t StripA_Id = 0;
    Int_t StripB_Id = 0;
    Double_t SlopA, SlopB, OffsetA = nan(""), OffsetB = nan("");

    // Middle layer
    Double_t Length2 = 30.06;                      // cm
    Double_t WidthMax2 = 10.4;                     // cm
    Double_t WidthMin2 = 1.3;                      // cm
    Double_t StripPitch2 = 0.005;                  // = 50 um
    Double_t AngTrap2 = 0.14853;                   // in rad = 17/2 degrees
    Double_t StepZ2 = StripPitch2 / sin(AngTrap2); // step along the z axis of the detector (in xz plan)
    Double_t StepX2 = StripPitch2 / cos(AngTrap2); // step along the x axis of the detector (in xz plan)
    Int_t NbStrip2 = int(WidthMax2 / StepX2);      //
    // LOG(info) << "NbStrip2= " << NbStrip2;

    // Outer layer
    Double_t Length3 = 30.06;                      // cm
    Double_t WidthMax3 = 10.4;                     // cm
    Double_t WidthMin3 = 1.3;                      // cm
    Double_t StripPitch3 = 0.005;                  // = 50 um
    Double_t AngTrap3 = 0.14853;                   // in rad = 17/2 degrees
    Double_t StepZ3 = StripPitch3 / sin(AngTrap3); // step along the z axis of the detector (in xz plan)
    Double_t StepX3 = StripPitch3 / cos(AngTrap3); // step along the x axis of the detector (in xz plan)
    Int_t NbStrip3 = int(WidthMax3 / StepX3);      //
    // LOG(info) << "NbStrip3= " << NbStrip3;

    Double_t M_Inner[6][4][4] = { {
                                      { 1, 0, 0, 0 },                                // Matrice 1 row 0
                                      { 0, 0.966376079, -0.257132793, 6.991301064 }, // Matrice 1 row 1
                                      { 0, 0.257132793, 0.966376079, -7.093624461 }, // Matrice 1 row 2
                                      { 0, 0, 0, 1 }                                 // Matrice 1 row 3
                                  },

                                  {
                                      { 0.5, -0.866025404, 0, 0 },                             // Matrice 2 row 0
                                      { 0.836906234, 0.483188040, -0.257132793, 6.991301064 }, // Matrice 2 row 1
                                      { 0.222683531, 0.128566397, 0.966376079, -7.093624461 }, // Matrice 2 row 2
                                      { 0, 0, 0, 1 }                                           // Matrice 2 row 3
                                  },

                                  {
                                      { -0.5, -0.866025404, 0, 0 },                             // Matrice 3 row 0
                                      { 0.836906234, -0.483188040, -0.257132793, 6.991301064 }, // Matrice 3 row 1
                                      { 0.222683531, -0.128566397, 0.966376079, -7.093624461 }, // Matrice 3 row 2
                                      { 0, 0, 0, 1 }                                            // Matrice 3 row 3
                                  },

                                  {
                                      { -1., 0., 0, 0 },                               // Matrice 4 row 0
                                      { 0., -0.966376079, -0.257132793, 6.991301064 }, // Matrice 4 row 1
                                      { 0., -0.257132793, 0.966376079, -7.093624461 }, // Matrice 4 row 2
                                      { 0., 0, 0, 1 }                                  // Matrice 4 row 3
                                  },

                                  {
                                      { -0.5, 0.866025404, 0, 0 },                               // Matrice 5 row 0
                                      { -0.836906234, -0.483188040, -0.257132793, 6.991301064 }, // Matrice 5 row 1
                                      { -0.222683531, -0.128566397, 0.966376079, -7.093624461 }, // Matrice 5 row 2
                                      { 0, 0, 0, 1 }                                             // Matrice 5 row 3
                                  },

                                  {
                                      { 0.5, 0.866025404, 0, 0 },                               // Matrice 6 row 0
                                      { -0.836906234, 0.483188040, -0.257132793, 6.991301064 }, // Matrice 6 row 1
                                      { -0.222683531, 0.128566397, 0.966376079, -7.093624461 }, // Matrice 6 row 2
                                      { 0, 0, 0, 1 }                                            // Matrice 6 row 3
                                  } };

    /*
     LOG(info) << M_Inner[0][0][0] << " " << M_Inner[0][0][1] << " " << M_Inner[0][0][2] << " " << M_Inner[0][0][3] << "
     "; LOG(info) << M_Inner[0][1][0] << " " << M_Inner[0][1][1] << " " << M_Inner[0][1][2] << " " << M_Inner[0][1][3]
     << " "; LOG(info) << M_Inner[0][2][0] << " " << M_Inner[0][2][1] << " " << M_Inner[0][2][2] << " " <<
     M_Inner[0][2][3] << " "; LOG(info) << M_Inner[0][3][0] << " " << M_Inner[0][3][1] << " " << M_Inner[0][3][2] << " "
     << M_Inner[0][3][3] << " ";
     */

    // Middle layer:
    Double_t M_Mid[12][4][4] = { {
                                     { 1, 0, 0, 0 },                                // Matrice 1 row 0
                                     { 0, 0.831954122, -0.554844427, 19.71955601 }, // Matrice 1 row 1
                                     { 0, 0.554844427, 0.831954122, -7.631968001 }, // Matrice 1 row 2
                                     { 0, 0, 0, 1 }                                 // Matrice 1 row 3
                                 },

                                 {
                                     { 0.866025404, -0.5, 0, 0 },                             // Matrice 2 row 0
                                     { 0.415977061, 0.720493405, -0.554844427, 19.71955601 }, // Matrice 2 row 1
                                     { 0.277422214, 0.480509369, 0.831954122, -7.631968001 }, // Matrice 2 row 2
                                     { 0, 0, 0, 1 }                                           // Matrice 2 row 3
                                 },

                                 {
                                     { 0.5, -0.866025404, 0, 0 },                             // Matrice 3 row 0
                                     { 0.720493405, 0.415977061, -0.554844427, 19.71955601 }, // Matrice 3 row 1
                                     { 0.480509369, 0.277422214, 0.831954122, -7.631968001 }, // Matrice 3 row 2
                                     { 0, 0, 0, 1 }                                           // Matrice 3 row 3
                                 },

                                 {
                                     { 0, -1, 0, 0 },                               // Matrice 4 row 0
                                     { 0.831954122, 0, -0.554844427, 19.71955601 }, // Matrice 4 row 1
                                     { 0.554844427, 0, 0.831954122, -7.631968001 }, // Matrice 4 row 2
                                     { 0, 0, 0, 1 }                                 // Matrice 4 row 3
                                 },

                                 {
                                     { -0.5, -0.866025404, 0, 0 },                             // Matrice 5 row 0
                                     { 0.720493405, -0.415977061, -0.554844427, 19.71955601 }, // Matrice 5 row 1
                                     { 0.480509369, -0.277422214, 0.831954122, -7.631968001 }, // Matrice 5 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 5 row 3
                                 },

                                 {
                                     { -0.866025404, -0.5, 0, 0 },                             // Matrice 6 row 0
                                     { 0.415977061, -0.720493405, -0.554844427, 19.71955601 }, // Matrice 6 row 1
                                     { 0.277422214, -0.480509369, 0.831954122, -7.631968001 }, // Matrice 6 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 6 row 3
                                 },

                                 {
                                     { -1, 0, 0, 0 },                                // Matrice 7 row 0
                                     { 0, -0.831954122, -0.554844427, 19.71955601 }, // Matrice 7 row 1
                                     { 0, -0.554844427, 0.831954122, -7.631968001 }, // Matrice 7 row 2
                                     { 0, 0, 0, 1 }                                  // Matrice 7 row 3
                                 },

                                 {
                                     { -0.866025404, 0.5, 0, 0 },                               // Matrice 8 row 0
                                     { -0.415977061, -0.720493405, -0.554844427, 19.71955601 }, // Matrice 8 row 1
                                     { -0.277422214, -0.480509369, 0.831954122, -7.631968001 }, // Matrice 8 row 2
                                     { 0, 0, 0, 1 }                                             // Matrice 8 row 3
                                 },

                                 {
                                     { -0.5, 0.866025404, 0, 0 },                               // Matrice 9 row 0
                                     { -0.720493405, -0.415977061, -0.554844427, 19.71955601 }, // Matrice 9 row 1
                                     { -0.480509369, -0.277422214, 0.831954122, -7.631968001 }, // Matrice 9 row 2
                                     { 0, 0, 0, 1 }                                             // Matrice 9 row 3
                                 },

                                 {
                                     { 0, 1, 0, 0 },                                 // Matrice 10 row 0
                                     { -0.831954122, 0, -0.554844427, 19.71955601 }, // Matrice 10 row 1
                                     { -0.554844427, 0, 0.831954122, -7.631968001 }, // Matrice 10 row 2
                                     { 0, 0, 0, 1 }                                  // Matrice 10 row 3
                                 },

                                 {
                                     { 0.5, 0.866025404, 0, 0 },                               // Matrice 11 row 0
                                     { -0.720493405, 0.415977061, -0.554844427, 19.71955601 }, // Matrice 11 row 1
                                     { -0.480509369, 0.277422214, 0.831954122, -7.631968001 }, // Matrice 11 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 11 row 3
                                 },

                                 {
                                     { 0.866025404, 0.5, 0, 0 },                               // Matrice 12 row 0
                                     { -0.415977061, 0.720493405, -0.554844427, 19.71955601 }, // Matrice 12 row 1
                                     { -0.277422214, 0.480509369, 0.831954122, -7.631968001 }, // Matrice 12 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 12 row 3
                                 }

    };

    // Outer layer:
    Double_t M_Out[12][4][4] = { {
                                     { 1, 0, 0, 0 },                                // Matrice 1 row 0
                                     { 0, 0.831954122, -0.554844427, 19.90455601 }, // Matrice 1 row 1
                                     { 0, 0.554844427, 0.831954122, -9.138672091 }, // Matrice 1 row 2
                                     { 0, 0, 0, 1 }                                 // Matrice 1 row 3
                                 },

                                 {
                                     { 0.866025404, -0.5, 0, 0 },                             // Matrice 2 row 0
                                     { 0.415977061, 0.720493405, -0.554844427, 19.90455601 }, // Matrice 2 row 1
                                     { 0.277422214, 0.480509369, 0.831954122, -9.138672091 }, // Matrice 2 row 2
                                     { 0, 0, 0, 1 }                                           // Matrice 2 row 3
                                 },

                                 {
                                     { 0.5, -0.866025404, 0, 0 },                             // Matrice 3 row 0
                                     { 0.720493405, 0.415977061, -0.554844427, 19.90455601 }, // Matrice 3 row 1
                                     { 0.480509369, 0.277422214, 0.831954122, -9.138672091 }, // Matrice 3 row 2
                                     { 0, 0, 0, 1 }                                           // Matrice 3 row 3
                                 },

                                 {
                                     { 0, -1, 0, 0 },                               // Matrice 4 row 0
                                     { 0.831954122, 0, -0.554844427, 19.90455601 }, // Matrice 4 row 1
                                     { 0.554844427, 0, 0.831954122, -9.138672091 }, // Matrice 4 row 2
                                     { 0, 0, 0, 1 }                                 // Matrice 4 row 3
                                 },

                                 {
                                     { -0.5, -0.866025404, 0, 0 },                             // Matrice 5 row 0
                                     { 0.720493405, -0.415977061, -0.554844427, 19.90455601 }, // Matrice 5 row 1
                                     { 0.480509369, -0.277422214, 0.831954122, -9.138672091 }, // Matrice 5 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 5 row 3
                                 },

                                 {
                                     { -0.866025404, -0.5, 0, 0 },                             // Matrice 6 row 0
                                     { 0.415977061, -0.720493405, -0.554844427, 19.90455601 }, // Matrice 6 row 1
                                     { 0.277422214, -0.480509369, 0.831954122, -9.138672091 }, // Matrice 6 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 6 row 3
                                 },

                                 {
                                     { -1, 0, 0, 0 },                                // Matrice 7 row 0
                                     { 0, -0.831954122, -0.554844427, 19.90455601 }, // Matrice 7 row 1
                                     { 0, -0.554844427, 0.831954122, -9.138672091 }, // Matrice 7 row 2
                                     { 0, 0, 0, 1 }                                  // Matrice 7 row 3
                                 },

                                 {
                                     { 0.866025404, -0.5, 0, 0 },                               // Matrice 8 row 0
                                     { -0.415977061, -0.720493405, -0.554844427, 19.90455601 }, // Matrice 8 row 1
                                     { -0.277422214, -0.480509369, 0.831954122, -9.138672091 }, // Matrice 8 row 2
                                     { 0, 0, 0, 1 }                                             // Matrice 8 row 3
                                 },

                                 {
                                     { -0.5, 0.866025404, 0, 0 },                               // Matrice 9 row 0
                                     { -0.720493405, -0.415977061, -0.554844427, 19.90455601 }, // Matrice 9 row 1
                                     { -0.480509369, -0.277422214, 0.831954122, -9.138672091 }, // Matrice 9 row 2
                                     { 0, 0, 0, 1 }                                             // Matrice 9 row 3
                                 },

                                 {
                                     { 0, 1, 0, 0 },                                 // Matrice 10 row 0
                                     { -0.831954122, 0, -0.554844427, 19.90455601 }, // Matrice 10 row 1
                                     { -0.554844427, 0, 0.831954122, -9.138672091 }, // Matrice 10 row 2
                                     { 0, 0, 0, 1 }                                  // Matrice 10 row 3
                                 },

                                 {
                                     { 0.5, 0.866025404, 0, 0 },                               // Matrice 11 row 0
                                     { -0.720493405, 0.415977061, -0.554844427, 19.90455601 }, // Matrice 11 row 1
                                     { -0.480509369, 0.277422214, 0.831954122, -9.138672091 }, // Matrice 11 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 11 row 3
                                 },

                                 {
                                     { 0.866025404, 0.5, 0, 0 },                               // Matrice 12 row 0
                                     { -0.415977061, 0.720493405, -0.554844427, 19.90455601 }, // Matrice 12 row 1
                                     { -0.277422214, 0.480509369, 0.831954122, -9.138672091 }, // Matrice 12 row 2
                                     { 0, 0, 0, 1 }                                            // Matrice 12 row 3
                                 }

    };

    // Transformation inverse matrices: (ie transformation from det coord. system to lab)

    // Inner layer:
    Double_t M_INV_Inner[6][4][4] = { {
                                          { 1, 0, 0, 0 },                                // Matrice 1 row 0
                                          { 0, 0.966376079, 0.257132793, -4.932222641 }, // Matrice 1 row 1
                                          { 0, -0.257132793, 0.966376079, 8.652801765 }, // Matrice 1 row 2
                                          { 0, 0, 0, 1 },                                // Matrice 1 row 3
                                      },

                                      {
                                          { 0.5, 0.836906234, 0.222683531, -4.271430104 },        // Matrice 2 row 0
                                          { -0.866025404, 0.48318804, 0.128566397, -2.46611132 }, // Matrice 2 row 1
                                          { 0, -0.257132793, 0.966376079, 8.652801765 },          // Matrice 2 row 2
                                          { 0, 0, 0, 1 },                                         // Matrice 2 row 3

                                      },

                                      {
                                          { -0.5, 0.836906234, 0.222683531, -4.271430104 },          // Matrice 3 row 0
                                          { -0.866025404, -0.483188040, -0.128566397, 2.466111320 }, // Matrice 3 row 1
                                          { 0, -0.257132793, 0.966376079, 8.652801765 },             // Matrice 3 row 2
                                          { 0, 0, 0, 1 },                                            // Matrice 3 row 3
                                      },

                                      {
                                          { -1., 0., 0., 0. },                             // Matrice 4 row 0
                                          { 0., -0.966376079, -0.257132793, 4.932222641 }, // Matrice 4 row 1
                                          { 0., -0.257132793, 0.966376079, 8.652801765 },  // Matrice 4 row 2
                                          { 0., 0., 0., 1. },                              // Matrice 4 row 3
                                      },

                                      {
                                          { -0.5, -0.836906234, -0.222683531, 4.271430104 },        // Matrice 5 row 0
                                          { 0.866025404, -0.483188040, -0.128566397, 2.466111320 }, // Matrice 5 row 1
                                          { 0., -0.257132793, 0.966376079, 8.652801765 },           // Matrice 5 row 2
                                          { 0., 0, 0, 1 },                                          // Matrice 5 row 3
                                      },

                                      {
                                          { 0.5, -0.836906234, -0.222683531, 4.271430104 },        // Matrice 6 row 0
                                          { 0.866025404, 0.483188040, 0.128566397, -2.466111320 }, // Matrice 6 row 1
                                          { 0., -0.257132793, 0.966376079, 8.652801765 },          // Matrice 6 row 2
                                          { 0., 0, 0, 1 },                                         // Matrice 6 row 3
                                      } };

    // Middle layer:
    Double_t M_INV_Mid[12][4][4] = { {
                                         { 1, 0, 0, 0 },                                // Matrice 1 row 0
                                         { 0, 0.831954122, 0.554844427, -12.17121099 }, // Matrice 1 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },   // Matrice 1 row 2
                                         { 0, 0, 0, 1 }                                 // Matrice 1 row 3
                                     },

                                     {
                                         { 0.866025404, 0.415977061, 0.277422214, -6.085605496 }, // Matrice 2 row 0
                                         { -0.5, 0.720493405, 0.480509369, -10.54057791 },        // Matrice 2 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },             // Matrice 2 row 2
                                         { 0, 0, 0, 1 }                                           // Matrice 2 row 3
                                     },

                                     {
                                         { 0.5, 0.720493405, 0.480509369, -10.54057791 },          // Matrice 3 row 0
                                         { -0.866025404, 0.415977061, 0.277422214, -6.085605496 }, // Matrice 3 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },              // Matrice 3 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 3 row 3
                                     },

                                     {
                                         { 0., 0.831954122, 0.554844427, -12.17121099 }, // Matrice 4 row 0
                                         { -1., 0., 0., 0. },                            // Matrice 4 row 1
                                         { 0., -0.554844427, 0.831954122, 17.290733 },   // Matrice 4 row 2
                                         { 0., 0, 0, 1 }                                 // Matrice 4 row 3
                                     },

                                     {
                                         { -0.5, 0.720493405, 0.480509369, -10.54057791 },          // Matrice 5 row 0
                                         { -0.866025404, -0.415977061, -0.277422214, 6.085605496 }, // Matrice 5 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },               // Matrice 5 row 2
                                         { 0, 0, 0, 1 }                                             // Matrice 5 row 3
                                     },

                                     {
                                         { -0.866025404, 0.415977061, 0.277422214, -6.085605496 }, // Matrice 6 row 0
                                         { -0.5, -0.720493405, -0.480509369, 10.54057791 },        // Matrice 6 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },              // Matrice 6 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 6 row 3
                                     },

                                     {
                                         { -1, 0, 0, 0 },                                // Matrice 7 row 0
                                         { 0, -0.831954122, -0.554844427, 12.17121099 }, // Matrice 7 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },    // Matrice 7 row 2
                                         { 0, 0, 0, 1 }                                  // Matrice 7 row 3
                                     },

                                     {
                                         { -0.866025404, -0.415977061, -0.277422214, 6.085605496 }, // Matrice 8 row 0
                                         { 0.5, -0.720493405, -0.480509369, 10.54057791 },          // Matrice 8 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },               // Matrice 8 row 2
                                         { 0, 0, 0, 1 }                                             // Matrice 8 row 3
                                     },

                                     {
                                         { -0.5, -0.720493405, -0.480509369, 10.54057791 },        // Matrice 9 row 0
                                         { 0.866025404, -0.415977061, -0.277422214, 6.085605496 }, // Matrice 9 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },              // Matrice 9 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 9 row 3
                                     },

                                     {
                                         { 0., -0.831954122, -0.554844427, 12.17121099 }, // Matrice 10 row 0
                                         { 1., 0., 0., 0. },                              // Matrice 10 row 1
                                         { 0., -0.554844427, 0.831954122, 17.290733 },    // Matrice 10 row 2
                                         { 0., 0, 0, 1 }                                  // Matrice 10 row 3
                                     },

                                     {
                                         { 0.5, -0.720493405, -0.480509369, 10.54057791 },        // Matrice 11 row 0
                                         { 0.866025404, 0.415977061, 0.277422214, -6.085605496 }, // Matrice 11 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },             // Matrice 11 row 2
                                         { 0, 0, 0, 1 }                                           // Matrice 11 row 3
                                     },

                                     {
                                         { 0.866025404, -0.415977061, -0.277422214, 6.085605496 }, // Matrice 12 row 0
                                         { 0.5, 0.720493405, 0.480509369, -10.54057791 },          // Matrice 12 row 1
                                         { 0, -0.554844427, 0.831954122, 17.290733 },              // Matrice 12 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 12 row 3
                                     } };

    // Outer layer:
    Double_t M_INV_Out[12][4][4] = { {
                                         { 1, 0, 0, 0 },                                // Matrice 1 row 0
                                         { 0, 0.831954122, 0.554844427, -11.48913614 }, // Matrice 1 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },  // Matrice 1 row 2
                                         { 0, 0, 0, 1 }                                 // Matrice 1 row 3
                                     },

                                     {
                                         { 0.866025404, 0.415977061, 0.277422214, -5.744568068 }, // Matrice 2 row 0
                                         { -0.5, 0.720493405, 0.480509369, -9.949883762 },        // Matrice 2 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },            // Matrice 2 row 2
                                         { 0, 0, 0, 1 }                                           // Matrice 2 row 3
                                     },

                                     {
                                         { 0.5, 0.720493405, 0.480509369, -9.949883762 },          // Matrice 3 row 0
                                         { -0.866025404, 0.415977061, 0.277422214, -5.744568068 }, // Matrice 3 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },             // Matrice 3 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 3 row 3
                                     },

                                     {
                                         { 0., 0.831954122, 0.554844427, -11.48913614 }, // Matrice 4 row 0
                                         { -1., 0., 0., 0. },                            // Matrice 4 row 1
                                         { 0., -0.554844427, 0.831954122, 18.6468879 },  // Matrice 4 row 2
                                         { 0., 0, 0, 1 }                                 // Matrice 4 row 3
                                     },

                                     {
                                         { -0.5, 0.720493405, 0.480509369, -9.949883762 },          // Matrice 5 row 0
                                         { -0.866025404, -0.415977061, -0.277422214, 5.744568068 }, // Matrice 5 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },              // Matrice 5 row 2
                                         { 0, 0, 0, 1 }                                             // Matrice 5 row 3
                                     },

                                     {
                                         { -0.866025404, 0.415977061, 0.277422214, -5.744568068 }, // Matrice 6 row 0
                                         { -0.5, -0.720493405, -0.480509369, 9.949883762 },        // Matrice 6 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },             // Matrice 6 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 6 row 3
                                     },

                                     {
                                         { -1, 0, 0, 0 },                                // Matrice 7 row 0
                                         { 0, -0.831954122, -0.554844427, 11.48913614 }, // Matrice 7 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },   // Matrice 7 row 2
                                         { 0, 0, 0, 1 }                                  // Matrice 7 row 3
                                     },

                                     {
                                         { 0.866025404, -0.415977061, -0.277422214, 5.744568068 }, // Matrice 8 row 0
                                         { -0.5, -0.720493405, -0.480509369, 9.949883762 },        // Matrice 8 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },             // Matrice 8 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 8 row 3
                                     },

                                     {
                                         { -0.5, -0.720493405, -0.480509369, 9.949883762 },        // Matrice 9 row 0
                                         { 0.866025404, -0.415977061, -0.277422214, 5.744568068 }, // Matrice 9 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },             // Matrice 9 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 9 row 3
                                     },

                                     {
                                         { 0., -0.831954122, -0.554844427, 11.48913614 }, // Matrice 10 row 0
                                         { 1., 0., 0., 0. },                              // Matrice 10 row 1
                                         { 0., -0.554844427, 0.831954122, 18.6468879 },   // Matrice 10 row 2
                                         { 0., 0, 0, 1 }                                  // Matrice 10 row 3
                                     },

                                     {
                                         { 0.5, -0.720493405, -0.480509369, 9.949883762 },        // Matrice 11 row 0
                                         { 0.866025404, 0.415977061, 0.277422214, -5.744568068 }, // Matrice 11 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },            // Matrice 11 row 2
                                         { 0, 0, 0, 1 }                                           // Matrice 11 row 3
                                     },

                                     {
                                         { 0.866025404, -0.415977061, -0.277422214, 5.744568068 }, // Matrice 12 row 0
                                         { 0.5, 0.720493405, 0.480509369, -9.949883762 },          // Matrice 12 row 1
                                         { 0, -0.554844427, 0.831954122, 18.6468879 },             // Matrice 12 row 2
                                         { 0, 0, 0, 1 }                                            // Matrice 12 row 3
                                     } };

    Double_t Energy = 0.;
    Double_t X_track = 0;
    Double_t Y_track = 0;
    Double_t Z_track = 0;
    Double_t X_Hit = 0;
    Double_t Y_Hit = 0;
    Double_t Z_Hit = 0;
    Double_t Px = 0.;
    Double_t Py = 0.;
    Double_t Pz = 0.;
    Double_t X_track_det = 0;
    Double_t Y_track_det = 0;
    Double_t Z_track_det = 0;
    Double_t X_intersect = 0.; // Position X in Detector frame
    Double_t Y_intersect = 0.; // Position Y in Detector frame
    Double_t Z_intersect = 0.; // Position Z in Detector frame
    Double_t Theta = 0.;       // Theta from (0,0,0)
    Double_t Phi = 0.;         // Phi (0,0,0)

    Int_t Detector;

    R3BTraPoint** traHit;

    Int_t traHitsPerEvent = 0;
    traHitsPerEvent = fTrackerHitCA->GetEntries();

    // LOG(info) << " NEW EVENT ";

    // LOG(info) << traHitsPerEvent;

    // LOG(info) << "";

    if (traHitsPerEvent > 0)
    {
        traHit = new R3BTraPoint*[traHitsPerEvent];
        for (Int_t i = 0; i < traHitsPerEvent; i++)
        {
            traHit[i] = new R3BTraPoint;
            traHit[i] = dynamic_cast<R3BTraPoint*>(fTrackerHitCA->At(i));
            Energy = ExpResSmearing(traHit[i]->GetEnergyLoss());
            // Energy = traHit[i]->GetEnergyLoss();
            Detector = traHit[i]->GetDetCopyID();

            X_track = traHit[i]->GetXIn();
            Y_track = traHit[i]->GetYIn();
            Z_track = traHit[i]->GetZIn();

            Px = traHit[i]->GetPxOut();
            Py = traHit[i]->GetPyOut();
            Pz = traHit[i]->GetPzOut();

            if (Detector <= 6)
            { // = inner layer

                // LOG(info) << Detector << " " << j << " " << k;

                // LOG(info) <<  X_track << " " <<   Y_track << " " <<  Z_track << " ";

                // apply transformation M_INV_inner[x][4][4]
                X_track_det = X_track * M_Inner[Detector - 1][0][0] + Y_track * M_Inner[Detector - 1][0][1] +
                              Z_track * M_Inner[Detector - 1][0][2] + M_Inner[Detector - 1][0][3];

                // LOG(info) << M_Inner[Detector-1][0][0] << " " << M_Inner[Detector-1][0][1] << " " <<
                // M_Inner[Detector-1][0][2] << " " << M_Inner[Detector-1][0][3] ;

                Y_track_det = X_track * M_Inner[Detector - 1][1][0] + Y_track * M_Inner[Detector - 1][1][1] +
                              Z_track * M_Inner[Detector - 1][1][2] + M_Inner[Detector - 1][1][3];

                // LOG(info) << M_Inner[Detector-1][1][0] << " " << M_Inner[Detector-1][1][1] << " " <<
                // M_Inner[Detector-1][1][2] << " " << M_Inner[Detector-1][1][3] ;

                Z_track_det = X_track * M_Inner[Detector - 1][2][0] + Y_track * M_Inner[Detector - 1][2][1] +
                              Z_track * M_Inner[Detector - 1][2][2] + M_Inner[Detector - 1][2][3];

                // LOG(info) << M_Inner[Detector-1][2][0] << " " << M_Inner[Detector-1][2][1] << " " <<
                // M_Inner[Detector-1][2][2] << " " << M_Inner[Detector-1][2][3] ;

                // LOG(info) <<  X_track_det << " " <<   Y_track_det << " " <<  Z_track_det << " ";

                // find 1st strip hit:
                // step 1 : projection parallel to the 1st longitudinal side of the detector (ie: offset at z=0 of this
                // straight line in plane xz)
                // SlopA= (2*Length1)/(WidthMax1-WidthMin1);
                SlopA = (2 * Length1) / (WidthMin1 - WidthMax1);
                Proj = Z_track_det - SlopA * X_track_det;

                // step 2 : compare with projection of middle line of each strips.
                strip = 0;
                ProjStrip = 0.;
                // ProjStrip= (-Length1/2-(StepZ1/2)*(2*strip+1)) - SlopA*(-WidthMax1/2);  // 1st strip proj at z=0
                ProjStrip =
                    (-Length1 / 2 - 1.1) -
                    SlopA *
                        (WidthMax1 / 2 -
                         (StepX1 / 2) *
                             (2 * strip +
                              1)); // 1st strip proj at x=0. -1.1 is a shift along z axis (TODO: find analytic formula).

                // LOG(info) << "SlopA= " << SlopA << " Length1= " << Length1 << " WidthMax1= " << WidthMax1 << "
                // WidthMin1= " << WidthMin1; LOG(info) << "Proj= " << Proj  << " ProjStrip= " << ProjStrip; LOG(info)
                // << "StepZ1/2 = " << StepZ1/2  << " ProjStrip - Proj= " << ProjStrip-Proj;

                while (strip < NbStrip1)
                {
                    if ((ProjStrip - Proj) < StepZ1 / 2)
                    {
                        // first strip hit is:
                        StripA_Id = strip;
                        // LOG(info) << "StripA_Id= " << StripA_Id;
                        OffsetA = ProjStrip;
                        strip = NbStrip1; // to end the loop
                    }
                    else
                    {
                        strip++;
                        // ProjStrip= (-Length1/2-(StepZ1/2)*(2*strip+1)) - SlopA*(-WidthMax1/2);
                        ProjStrip = (-Length1 / 2 - 1.1) - SlopA * (WidthMax1 / 2 - (StepX1 / 2) * (2 * strip + 1));
                    }
                }

                // find 2nd strip hit:
                // step 1 : projection parallel to the 2nd longitudinal side of the detector (ie: offset at z=0 of this
                // straight line in plane xz)
                SlopB = -SlopA;
                Proj = Z_track_det - SlopB * X_track_det;
                // step 2 : compare with projection of middle line of each strips.
                strip = 0;
                ProjStrip = 0.;
                // ProjStrip= (-Length1/2-(StepZ1/2)*(2*strip+1)) - SlopB*(WidthMax1/2);  // 1st strip proj at z=0
                ProjStrip = (-Length1 / 2 - 1.1) - SlopB * (-WidthMax1 / 2 + (StepX1 / 2) * (2 * strip + 1));

                // LOG(info) << "SlopB= " << SlopB << " Length1= " << Length1 << " WidthMax1= " << WidthMax1 << "
                // WidthMin1= " << WidthMin1; LOG(info) << "Proj= " << Proj  << " ProjStrip= " << ProjStrip; LOG(info)
                // << "StepZ1/2 = " << StepZ1/2  << " ProjStrip - Proj= " << ProjStrip-Proj;

                while (strip < NbStrip1)
                {
                    if ((ProjStrip - Proj) < StepZ1 / 2)
                    {
                        // 2nd strip hit:
                        StripB_Id = strip;
                        // LOG(info) << "StripB_Id= " << StripB_Id;
                        OffsetB = ProjStrip;
                        strip = NbStrip1; // to end the loop
                    }
                    else
                    {
                        strip++;
                        // ProjStrip= (-Length1/2-(StepZ1/2)*(2*strip+1)) - SlopB*(WidthMax1/2);
                        ProjStrip = (-Length1 / 2 - 1.1) - SlopB * (-WidthMax1 / 2 + (StepX1 / 2) * (2 * strip + 1));
                    }
                }

                // find intersection of the 2 hit strips:
                X_intersect = (OffsetB - OffsetA) / (SlopA - SlopB);
                // X_intersect= X_track_det[j];
                Y_intersect = Y_track_det;
                Z_intersect = SlopA * X_intersect + OffsetA;
                // Z_intersect= Z_track_det[j];

                // then transform back in Lab frame:

                X_Hit = X_intersect * M_INV_Inner[Detector - 1][0][0] + Y_intersect * M_INV_Inner[Detector - 1][0][1] +
                        Z_intersect * M_INV_Inner[Detector - 1][0][2] + M_INV_Inner[Detector - 1][0][3];
                Y_Hit = X_intersect * M_INV_Inner[Detector - 1][1][0] + Y_intersect * M_INV_Inner[Detector - 1][1][1] +
                        Z_intersect * M_INV_Inner[Detector - 1][1][2] + M_INV_Inner[Detector - 1][1][3];
                Z_Hit = X_intersect * M_INV_Inner[Detector - 1][2][0] + Y_intersect * M_INV_Inner[Detector - 1][2][1] +
                        Z_intersect * M_INV_Inner[Detector - 1][2][2] + M_INV_Inner[Detector - 1][2][3];

                // LOG(info) <<  X_Hit << " " <<   Y_Hit << " " <<  Z_Hit << " ";

            } // end of Det=7 ( inner layer)

            if (Detector >= 7 && Detector <= 18)
            { // = Middle layer

                // LOG(info) << Detector;

                // apply transformation M_INV_Mid[x][4][4]
                X_track_det = X_track * M_Mid[Detector - 7][0][0] + Y_track * M_Mid[Detector - 7][0][1] +
                              Z_track * M_Mid[Detector - 7][0][2] + M_Mid[Detector - 7][0][3];
                Y_track_det = X_track * M_Mid[Detector - 7][1][0] + Y_track * M_Mid[Detector - 7][1][1] +
                              Z_track * M_Mid[Detector - 7][1][2] + M_Mid[Detector - 7][1][3];
                Z_track_det = X_track * M_Mid[Detector - 7][2][0] + Y_track * M_Mid[Detector - 7][2][1] +
                              Z_track * M_Mid[Detector - 7][2][2] + M_Mid[Detector - 7][2][3];

                // LOG(info) <<  X_track << " " <<   Y_track << " " <<  Z_track << " ";
                // LOG(info) <<  X_track_det << " " <<   Y_track_det << " " <<  Z_track_det << " ";

                // find 1st strip hit:
                // step 1 : projection parallel to the 1st longitudinal side of the detector (ie: offset at z=0 of this
                // straight line in plane xz)
                // SlopA= (2*Length2)/(WidthMax2-WidthMin2);
                SlopA = (2 * Length2) / (WidthMin2 - WidthMax2);
                Proj = Z_track_det - SlopA * X_track_det;

                // step 2 : compare with projection of middle line of each strips.
                strip = 0;
                ProjStrip = 0.;
                // ProjStrip= (-Length2/2-(StepZ2/2)*(2*strip+1)) - SlopA*(-WidthMax2/2);  // 1st strip proj at z=0
                ProjStrip = (-Length2 / 2 - 5.5) -
                            SlopA * (WidthMax2 / 2 -
                                     (StepX2 / 2) * (2 * strip + 1)); // 1st strip proj at x=0, -5.5 ia a shift along Z
                                                                      // axis (TODO: find the analytic formula)

                // LOG(info) << "SlopA= " << SlopA << " Length2= " << Length2 << " WidthMax2= " << WidthMax2 << "
                // WidthMin2= " << WidthMin2; LOG(info) << "Proj= " << Proj  << " ProjStrip= " << ProjStrip; LOG(info)
                // << "StepZ2/2 = " << StepZ2/2  << " ProjStrip - Proj= " << ProjStrip-Proj;

                while (strip < NbStrip2)
                {
                    if ((ProjStrip - Proj) < StepZ2 / 2)
                    {
                        // first strip hit is:
                        StripA_Id = strip;
                        // LOG(info) << "StripA_Id= " << StripA_Id;
                        OffsetA = ProjStrip;
                        strip = NbStrip2; // to end the loop
                    }
                    else
                    {
                        strip++;
                        // ProjStrip= (-Length2/2-(StepZ2/2)*(2*strip+1)) - SlopA*(-WidthMax2/2);
                        ProjStrip = (-Length2 / 2 - 5.5) - SlopA * (WidthMax2 / 2 - (StepX2 / 2) * (2 * strip + 1));
                    }
                }

                // find 2nd strip hit:
                // step 1 : projection parallel to the 2nd longitudinal side of the detector (ie: offset at z=0 of this
                // straight line in plane xz)
                SlopB = -SlopA;
                Proj = Z_track_det - SlopB * X_track_det;
                // step 2 : compare with projection of middle line of each strips.
                strip = 0;
                ProjStrip = 0.;
                // ProjStrip= (-Length2/2-(StepZ2/2)*(2*strip+1)) - SlopB*(WidthMax2/2);  // 1st strip proj at z=0
                ProjStrip = (-Length2 / 2 - 5.5) - SlopB * (-WidthMax2 / 2 + (StepX2 / 2) * (2 * strip + 1));

                while (strip < NbStrip2)
                {
                    if ((ProjStrip - Proj) < StepZ2 / 2)
                    {
                        // 2nd strip hit:
                        StripB_Id = strip;
                        // LOG(info) << "StripB_Id= " << StripB_Id;
                        OffsetB = ProjStrip;
                        strip = NbStrip2; // to end the loop
                    }
                    else
                    {
                        strip++;
                        // ProjStrip= (-Length2/2-(StepZ2/2)*(2*strip+1)) - SlopB*(WidthMax2/2);
                        ProjStrip = (-Length2 / 2 - 5.5) - SlopB * (-WidthMax2 / 2 + (StepX2 / 2) * (2 * strip + 1));
                    }
                }

                // find intersection of the 2 hit strips:
                X_intersect = (OffsetB - OffsetA) / (SlopA - SlopB);
                // X_intersect= X_track_det[j];
                Y_intersect = Y_track_det;
                Z_intersect = SlopA * X_intersect + OffsetA;
                // Z_intersect= Z_track_det[j];

                // then transform back in Lab frame:

                X_Hit = X_intersect * M_INV_Mid[Detector - 7][0][0] + Y_intersect * M_INV_Mid[Detector - 7][0][1] +
                        Z_intersect * M_INV_Mid[Detector - 7][0][2] + M_INV_Mid[Detector - 7][0][3];
                Y_Hit = X_intersect * M_INV_Mid[Detector - 7][1][0] + Y_intersect * M_INV_Mid[Detector - 7][1][1] +
                        Z_intersect * M_INV_Mid[Detector - 7][1][2] + M_INV_Mid[Detector - 7][1][3];
                Z_Hit = X_intersect * M_INV_Mid[Detector - 7][2][0] + Y_intersect * M_INV_Mid[Detector - 7][2][1] +
                        Z_intersect * M_INV_Mid[Detector - 7][2][2] + M_INV_Mid[Detector - 7][2][3];

                // LOG(info) <<  X_Hit << " " <<   Y_Hit << " " <<  Z_Hit << " ";
            }

            if (Detector >= 19)
            { // = outer layer

                // LOG(info) << Detector;

                // apply transformation M_INV_Out[x][4][4]
                X_track_det = X_track * M_Out[Detector - 19][0][0] + Y_track * M_Out[Detector - 19][0][1] +
                              Z_track * M_Out[Detector - 19][0][2] + M_Out[Detector - 19][0][3];
                Y_track_det = X_track * M_Out[Detector - 18][1][0] + Y_track * M_Out[Detector - 19][1][1] +
                              Z_track * M_Out[Detector - 19][1][2] + M_Out[Detector - 19][1][3];
                Z_track_det = X_track * M_Out[Detector - 19][2][0] + Y_track * M_Out[Detector - 19][2][1] +
                              Z_track * M_Out[Detector - 19][2][2] + M_Out[Detector - 19][2][3];

                // LOG(info) <<  X_track << " " <<   Y_track << " " <<  Z_track << " ";
                // LOG(info) <<  X_track_det << " " <<   Y_track_det << " " <<  Z_track_det << " ";

                // find 1st strip hit:
                // step 1 : projection parallel to the 1st longitudinal side of the detector (ie: offset at z=0 of this
                // straight line in plane xz)
                SlopA = (2 * Length3) / (WidthMin3 - WidthMax3);
                // SlopA= (2*Length3)/(WidthMax3-WidthMin3);
                Proj = Z_track_det - SlopA * X_track_det;

                // step 2 : compare with projection of middle line of each strips.
                strip = 0;
                ProjStrip = 0.;
                // ProjStrip= (-Length3/2-(StepZ3/2)*(2*strip+1)) - SlopA*(-WidthMax3/2);  // 1st strip proj at z=0
                ProjStrip = (-Length3 / 2 - 5) -
                            SlopA * (WidthMax3 / 2 - (StepX3 / 2) * (2 * strip + 1)); // 1st strip proj at x=0

                // LOG(info) << "SlopA= " << SlopA << " Length3= " << Length3 << " WidthMax3= " << WidthMax3 << "
                // WidthMin3= " << WidthMin3; LOG(info) << "Proj= " << Proj  << " ProjStrip= " << ProjStrip; LOG(info)
                // << "StepZ3/2 = " << StepZ3/2  << " ProjStrip - Proj= " << ProjStrip-Proj;

                while (strip < NbStrip3)
                {
                    if ((ProjStrip - Proj) < StepZ3 / 2)
                    {
                        // first strip hit is:
                        StripA_Id = strip;
                        // LOG(info) << "StripA_Id= " << StripA_Id;
                        OffsetA = ProjStrip;
                        strip = NbStrip3; // to end the loop
                    }
                    else
                    {
                        strip++;
                        // ProjStrip= (-Length3/2-(StepZ3/2)*(2*strip+1)) - SlopA*(-WidthMax3/2);
                        ProjStrip = (-Length3 / 2 - 5) - SlopA * (WidthMax3 / 2 - (StepX3 / 2) * (2 * strip + 1));
                    }
                }

                // find 2nd strip hit:
                // step 1 : projection parallel to the 2nd longitudinal side of the detector (ie: offset at z=0 of this
                // straight line in plane xz)
                SlopB = -SlopA;
                Proj = Z_track_det - SlopB * X_track_det;
                // step 2 : compare with projection of middle line of each strips.
                strip = 0;
                ProjStrip = 0.;
                // ProjStrip= (-Length3/2-(StepZ3/2)*(2*strip+1)) - SlopB*(WidthMax3/2);  // 1st strip proj at z=0
                ProjStrip = (-Length3 / 2 - 5) -
                            SlopB * (-WidthMax3 / 2 + (StepX3 / 2) * (2 * strip + 1)); // 1st strip proj at x=0

                while (strip < NbStrip3)
                {
                    if ((ProjStrip - Proj) < StepZ3 / 2)
                    {
                        // 2nd strip hit:
                        StripB_Id = strip;
                        // LOG(info) << "StripB_Id= " << StripB_Id;
                        OffsetB = ProjStrip;
                        strip = NbStrip3; // to end the loop
                    }
                    else
                    {
                        strip++;
                        // ProjStrip= (-Length3/2-(StepZ3/2)*(2*strip+1)) - SlopB*(WidthMax3/2);
                        ProjStrip = (-Length3 / 2 - 5) -
                                    SlopB * (-WidthMax3 / 2 + (StepX3 / 2) * (2 * strip + 1)); // 1st strip proj at x=0
                    }
                }

                // find intersection of the 2 hit strips:
                X_intersect = (OffsetB - OffsetA) / (SlopA - SlopB);
                // X_intersect= X_track_det[j];
                Y_intersect = Y_track_det;
                Z_intersect = SlopA * X_intersect + OffsetA;
                // Z_intersect= Z_track_det[j];

                // then transform back in Lab frame:

                X_Hit = X_intersect * M_INV_Out[Detector - 19][0][0] + Y_intersect * M_INV_Out[Detector - 19][0][1] +
                        Z_intersect * M_INV_Out[Detector - 19][0][2] + M_INV_Out[Detector - 19][0][3];
                Y_Hit = X_intersect * M_INV_Out[Detector - 19][1][0] + Y_intersect * M_INV_Out[Detector - 19][1][1] +
                        Z_intersect * M_INV_Out[Detector - 19][1][2] + M_INV_Out[Detector - 19][1][3];
                Z_Hit = X_intersect * M_INV_Out[Detector - 19][2][0] + Y_intersect * M_INV_Out[Detector - 19][2][1] +
                        Z_intersect * M_INV_Out[Detector - 19][2][2] + M_INV_Out[Detector - 19][2][3];

                // LOG(info) <<  X_Hit << " " <<   Y_Hit << " " <<  Z_Hit << " ";
            }

            Theta = GetThetaScatZero(X_Hit, Y_Hit, Z_Hit);
            Phi = GetPhiScatZero(X_Hit, Y_Hit, Z_Hit);

            /*if ((X_track/X_Hit)<0 || (Y_track/Y_Hit)<0) {

                LOG(info) << "Detector" << " " << "Smeared Energy" << " " << "Raw Energy";
                LOG(info) << Detector << " " << Energy << " " << traHit[i]->GetEnergyLoss();
                LOG(info) << "X_in" << " " << "Y_in" << " " << "Z_in";
                LOG(info) << X_track << " " << Y_track << " " << Z_track;
                LOG(info) << "X_out" << " " << "Y_out" << " " << "Z_out";
                LOG(info) << X_Hit << " " << Y_Hit << " " << Z_Hit;
                LOG(info) << "Theta " << "Phi";
                LOG(info) << Theta << " " << Phi;
                LOG(info) << " ";

            }*/

            if (Energy >= fThreshold)
                AddHit(Energy, Detector, X_Hit, Y_Hit, Z_Hit, Px, Py, Pz, Theta, Phi);
            // AddHit(Energy, Detector);
        }
    }

    // LOG(info) << " ";
}

// ---- Public method Reset   --------------------------------------------------
void R3BTraHitFinder::Reset()
{
    // Clear the CA structure
    // LOG(info) << " -I- R3BTraHitFinder:Reset() called ";

    if (fTraHitCA)
        fTraHitCA->Clear();
}

// ---- Public method Finish   --------------------------------------------------
void R3BTraHitFinder::Finish()
{
    // here event. write histos
    LOG(info) << "";
    LOG(info) << " -I- Digit Finish() called ";
    // Write control histograms
}

// -----  Public method SetExperimentalResolution  ----------------------------------
void R3BTraHitFinder::SetExperimentalResolution(Double_t trackerRes)
{
    fTrackerResolution = trackerRes;
    LOG(info) << "-I- R3BTraHitFinder::SetExperimentalResolution to " << fTrackerResolution << " GeV";
}

// -----  Public method SetDetectionThreshold  ----------------------------------
void R3BTraHitFinder::SetDetectionThreshold(Double_t thresholdEne)
{
    fThreshold = thresholdEne;
    LOG(info) << "-I- R3BTraHitFinder::SetDetectionThreshold to " << fThreshold << " GeV";
}

// -----   Private method ExpResSmearing  --------------------------------------------
Double_t R3BTraHitFinder::ExpResSmearing(Double_t inputEnergy)
{

    if (fTrackerResolution == 0)
        return inputEnergy;
    else
    {
        Double_t energy = gRandom->Gaus(inputEnergy, fTrackerResolution);
        // LOG(info) << "energy " << energy << " for and Energy of "<< inputEnergy ;
        return energy;
    }
}

// -----   Private method GetThetaScatZero  --------------------------------------------
Double_t R3BTraHitFinder::GetThetaScatZero(Double_t X, Double_t Y, Double_t Z)
{

    Double_t Theta;

    Theta = atan(sqrt(pow(X, 2) + pow(Y, 2)) / Z);

    return Theta;
}

// -----   Private method GetPhiScatZero  --------------------------------------------
Double_t R3BTraHitFinder::GetPhiScatZero(Double_t X, Double_t Y, Double_t Z) { return atan2(Y, X); }

// -----   Private method AddHit  --------------------------------------------
R3BTrackerHit* R3BTraHitFinder::AddHit(Double_t ene,
                                       Int_t det,
                                       Double_t x,
                                       Double_t y,
                                       Double_t z,
                                       Double_t px,
                                       Double_t py,
                                       Double_t pz,
                                       Double_t th,
                                       Double_t phi)
{
    // R3BTrackerHit* R3BTraHitFinder::AddHit(Double_t ene,Int_t det){
    TClonesArray& clref = *fTraHitCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTrackerHit(ene, det, x, y, z, px, py, pz, th, phi);
    // return new(clref[size]) R3BTrackerHit(ene, det);
}

ClassImp(R3BTraHitFinder)
