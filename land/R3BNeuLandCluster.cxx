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

// -----------------------------------------------------------------------------
// -----                       R3BNeuLandCluster                           -----
// -----                 Created 12-07-2011 by D.Kresan                    -----
// -----------------------------------------------------------------------------

#include "R3BNeuLandCluster.h"

ClassImp(R3BNeuLandCluster)

    // -----------------------------------------------------------------------------
    R3BNeuLandCluster::R3BNeuLandCluster()
{
    fStartX = 0.;
    fStartY = 0.;
    fStartZ = 0.;
    fStartT = 0.;
    fStopX = 0.;
    fStopY = 0.;
    fStopZ = 0.;
    fStopT = 0.;
    fE = 0.;
    fSize = 0.;
    fMarked = kFALSE;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BNeuLandCluster::R3BNeuLandCluster(Double_t startX,
                                     Double_t startY,
                                     Double_t startZ,
                                     Double_t startT,
                                     Double_t stopX,
                                     Double_t stopY,
                                     Double_t stopZ,
                                     Double_t stopT,
                                     Double_t e,
                                     Double_t size)
{
    fStartX = startX;
    fStartY = startY;
    fStartZ = startZ;
    fStartT = startT;
    fStopX = stopX;
    fStopY = stopY;
    fStopZ = stopZ;
    fStopT = stopT;
    fE = e;
    fSize = size;
    fMarked = kFALSE;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BNeuLandCluster::~R3BNeuLandCluster()
{
    fStartX = 0.;
    fStartY = 0.;
    fStartZ = 0.;
    fStartT = 0.;
    fStopX = 0.;
    fStopY = 0.;
    fStopZ = 0.;
    fStopT = 0.;
    fE = 0.;
    fSize = 0.;
    fMarked = kFALSE;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::SetStartPosition(Double_t startX, Double_t startY, Double_t startZ)
{
    fStartX = startX;
    fStartY = startY;
    fStartZ = startZ;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::SetStopPosition(Double_t stopX, Double_t stopY, Double_t stopZ)
{
    fStopX = stopX;
    fStopY = stopY;
    fStopZ = stopZ;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::SetStartT(Double_t startT) { fStartT = startT; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::SetStopT(Double_t stopT) { fStopT = stopT; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::SetE(Double_t e) { fE = e; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::SetSize(Double_t size) { fSize = size; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::Mark() { fMarked = kTRUE; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::StartPosition(TVector3& pos) const { pos.SetXYZ(fStartX, fStartY, fStartZ); }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeuLandCluster::StopPosition(TVector3& pos) const { pos.SetXYZ(fStopX, fStopY, fStopZ); }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Double_t R3BNeuLandCluster::GetStartT() const { return fStartT; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Double_t R3BNeuLandCluster::GetStopT() const { return fStopT; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Double_t R3BNeuLandCluster::GetE() const { return fE; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Double_t R3BNeuLandCluster::GetSize() const { return fSize; }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Bool_t R3BNeuLandCluster::IsMarked() const { return fMarked; }
// -----------------------------------------------------------------------------
