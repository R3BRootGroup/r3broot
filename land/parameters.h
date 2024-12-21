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

static const Int_t NofUnitLayers = 25;               // number of detector units
static const Int_t NofDetLayers = 2 * NofUnitLayers; // number of detector layers
static const Int_t NofSegs = 4;                      // number of detector segments
static const Int_t NofStrips = 19;                   // number of readout strips in one unit
static const Int_t NofCells = 6;                     // number of gas cells in one layer
static const Double_t StripWidth = 2.5;              // strip width in cm
static const Double_t AnGapWidth = 0.15;             // anode gap width in cm
static const Double_t SegGapWidth = 5.0;             // segment gap width in cm
static const Double_t DetH = (NofStrips * StripWidth + (NofStrips - 1) * AnGapWidth); // height of detector (x)
static const Double_t LayerH = NofSegs * (DetH + SegGapWidth);                        // height of layer (x)
static const Double_t DetW = 200.0;                                                   // width of detector (y) in cm
static const Double_t St1T = 0.2;                                                     // cm
static const Double_t Gas1T = 0.32;                                                   // cm
static const Double_t Myl1T = 0.09;                                                   // cm
static const Double_t Myl2T = 0.04;                                                   // cm
static const Double_t GlassT = 0.10;                                                  // cm
static const Double_t Gas2T = 0.03;                                                   // cm
static const Double_t St2T = 0.4;                                                     // cm
static const Double_t DetT = 2 * (St1T + Myl1T + (NofCells / 2) * (GlassT + Gas2T) + GlassT + Myl2T) + Gas1T + St2T;
static const Double_t UnitT = 2 * DetT;
