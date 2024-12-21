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

#include "R3BNeulandGeoPar.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TVector3.h"
#include <algorithm>
#include <iostream>

#include "FairParamList.h"

R3BNeulandGeoPar::R3BNeulandGeoPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNeulandGeoNode(nullptr)
{
}

R3BNeulandGeoPar::~R3BNeulandGeoPar()
{
    // Note: Deleting stuff here or in clear() causes segfaults?
}

void R3BNeulandGeoPar::clear() {}

void R3BNeulandGeoPar::putParams(FairParamList* l)
{
    if (!l)
    {
        return;
    }
    l->addObject("NeulandGeoNode", fNeulandGeoNode);
}

Bool_t R3BNeulandGeoPar::getParams(FairParamList* l)
{
    if (!l)
    {
        return kFALSE;
    }
    if (!l->fillObject("NeulandGeoNode", fNeulandGeoNode))
    {
        return kFALSE;
    }
    BuildPaddleLookup();
    return kTRUE;
}

void R3BNeulandGeoPar::printParams()
{
    std::cout << "R3BNeulandGeoPar: Neuland Paddle Positions ..." << std::endl;

    fNeulandGeoNode->GetMatrix()->Print();

    for (Int_t i = 0; i < fNeulandGeoNode->GetNdaughters(); i++)
    {
        TGeoNode* node = fNeulandGeoNode->GetDaughter(i);
        std::cout << node->GetNumber() << std::endl;
        node->GetMatrix()->Print();
    }
}

void R3BNeulandGeoPar::SetNeulandGeoNode(const TGeoNode* const p)
{
    fNeulandGeoNode = dynamic_cast<TGeoNode*>(p->Clone());
    BuildPaddleLookup();
}

Double_t R3BNeulandGeoPar::GetPaddleHalfLength() const
{
    // All paddles have to have the same length
    return (dynamic_cast<TGeoBBox*>(fNeulandGeoNode->GetDaughter(0)->GetVolume()->GetShape()))->GetDX();
}

// Convert positions of e.g. points to the local coordinate of the respective paddle [(-135,135),(-2.5,2.5),(-2.5,2.5)]
TVector3 R3BNeulandGeoPar::ConvertToLocalCoordinates(const TVector3& position, const Int_t paddleID) const
{
    Double_t pos_in[3] = { position.X(), position.Y(), position.Z() };
    Double_t pos_tmp[3];
    Double_t pos_out[3];

    // First, convert to Neuland-local coordinates (consisting of all paddles)
    fNeulandGeoNode->GetMatrix()->MasterToLocal(pos_in, pos_tmp);
    // Second, convert to the repective paddle
    fPaddleGeoNodes.at(paddleID)->MasterToLocal(pos_tmp, pos_out);

    return TVector3(pos_out[0], pos_out[1], pos_out[2]);
}

TVector3 R3BNeulandGeoPar::ConvertToGlobalCoordinates(const TVector3& position, const Int_t paddleID) const
{
    Double_t pos_in[3] = { position.X(), position.Y(), position.Z() };
    Double_t pos_tmp[3];
    Double_t pos_out[3];

    // Note reverse order of Global->Local
    fPaddleGeoNodes.at(paddleID)->LocalToMaster(pos_in, pos_tmp);

    fNeulandGeoNode->GetMatrix()->LocalToMaster(pos_tmp, pos_out);

    return TVector3(pos_out[0], pos_out[1], pos_out[2]);
}

TVector3 R3BNeulandGeoPar::ConvertGlobalToPixel(const TVector3& position) const
{
    const Int_t nPixels = 50;
    const Double_t sizePixel = 5;
    const Int_t nPaddles = fNeulandGeoNode->GetNdaughters();
    const Int_t nPlanes = nPaddles / nPixels;

    Double_t pos_in[3] = { position.X(), position.Y(), position.Z() };
    Double_t pos_tmp[3];

    // First, convert to Neuland-local coordinates (consisting of all paddles)
    fNeulandGeoNode->GetMatrix()->MasterToLocal(pos_in, pos_tmp);

    // Note: PaddleHalfLength is 135 (light guides)
    // Map x and y values with [-125.:125.] float to [0:nPixels-1] int
    const Int_t x = std::min(std::max<Int_t>(0, pos_tmp[0] / sizePixel + nPixels / 2), nPixels - 1);
    const Int_t y = std::min(std::max<Int_t>(0, pos_tmp[1] / sizePixel + nPixels / 2), nPixels - 1);

    // Map z to [0:nPlanes-1]
    const Int_t z = std::min(std::max<Int_t>(0, pos_tmp[2] / sizePixel + nPlanes / 2), nPlanes - 1);

    return TVector3(x, y, z);
}

void R3BNeulandGeoPar::BuildPaddleLookup()
{
    for (Int_t i = 0; i < fNeulandGeoNode->GetNdaughters(); i++)
    {
        TGeoNode* node = fNeulandGeoNode->GetDaughter(i);
        fPaddleGeoNodes[node->GetNumber()] = node;
    }
}
ClassImp(R3BNeulandGeoPar);
