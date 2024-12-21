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

#ifndef R3BNEULANDGEOPAR_H
#define R3BNEULANDGEOPAR_H

#include "FairParGenericSet.h"
#include "TGeoNode.h"
#include <map>
class FairParamList;
class TVector3;

/**
 * NeuLAND geometry parameter storage
 * @author Jan Mayer
 *
 * Stores the full Neuland geo node used in the simulation for later reference, especially for coordinate
 * transformation from and to local and global coordinates.
 */

class R3BNeulandGeoPar : public FairParGenericSet
{
  public:
    // needs to be public?
    TGeoNode* fNeulandGeoNode;

    R3BNeulandGeoPar(const char* name = "R3BNeulandGeoPar",
                     const char* title = "Neuland Geometry Parameters",
                     const char* context = "TestDefaultContext");
    ~R3BNeulandGeoPar() override;

    void clear() override;
    void putParams(FairParamList*) override;
    Bool_t getParams(FairParamList*) override;
    void printParams() override;

    TGeoNode* GetNeulandGeoNode() { return fNeulandGeoNode; } // FIXME: const?
    void SetNeulandGeoNode(const TGeoNode* const p);

    Double_t GetPaddleHalfLength() const;
    TVector3 ConvertToLocalCoordinates(const TVector3& position, const Int_t paddleID) const;
    TVector3 ConvertToGlobalCoordinates(const TVector3& position, const Int_t paddleID) const;
    TVector3 ConvertGlobalToPixel(const TVector3& position) const;

  private:
    std::map<Int_t, TGeoNode*> fPaddleGeoNodes;
    void BuildPaddleLookup();

    R3BNeulandGeoPar(const R3BNeulandGeoPar&);
    R3BNeulandGeoPar& operator=(const R3BNeulandGeoPar&);

    ClassDefOverride(R3BNeulandGeoPar, 1)
};

#endif // R3BNEULANDGEOPAR_H
