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

#ifndef R3BCALIFAHITPAR_H
#define R3BCALIFAHITPAR_H

#include <TObjString.h>
#include <TVector3.h>

#include "FairParGenericSet.h"
#include "FairParamList.h"

class R3BCalifaClusterPar : public FairParGenericSet
{
  public:
    R3BCalifaClusterPar(const char* name = "R3BCalifaClusterPar",
                        const char* title = "The Cluster Finder Parameters",
                        const char* context = "TestDefaultContext");
    ~R3BCalifaClusterPar(void){};
    void clear(void){};
    void putParams(FairParamList* list);
    Bool_t getParams(FairParamList* list);

    void Print(Option_t* option = "") const;
    /** Accessor functions **/
    const Double_t GetThreshold() { return fThreshold; }
    const Double_t GetDRThreshold() { return fDRThreshold; }
    const Double_t GetDeltaPolar() { return fDeltaPolar; }
    const Double_t GetDeltaAzimuthal() { return fDeltaAzimuthal; }
    const Double_t GetDeltaAngleClust() { return fDeltaAngleClust; }
    const Int_t GetClusteringAlgorithmSelector() { return fClusterAlgorithmSelector; }
    const Double_t GetParCluster1() { return fParCluster1; }

    void SetThreshold(Double_t value) { fThreshold = value; }
    void SetDRThreshold(Double_t value) { fDRThreshold = value; }
    void SetDeltaPolar(Double_t value) { fDeltaPolar = value; }
    void SetDeltaAzimuthal(Double_t value) { fDeltaAzimuthal = value; }
    void SetDeltaAngleClust(Double_t value) { fDeltaAngleClust = value; }
    void SetClusteringAlgorithmSelector(Int_t value) { fClusterAlgorithmSelector = value; }
    void SetParCluster1(Double_t value) { fParCluster1 = value; }

  private:
    // Minimum energy requested in a crystal to be considered in a calorimeter Cluster
    Double_t fThreshold;
    // Threshold for selecting gamma or proton branch in double reading channels
    Double_t fDRThreshold;
    // Angular window (polar)
    Double_t fDeltaPolar;
    // Angular window (azimuthal)
    Double_t fDeltaAzimuthal;
    // Angular opening used for the clustering condition
    Double_t fDeltaAngleClust;
    // Clustering algorithm selector
    Int_t fClusterAlgorithmSelector;
    // Clustering parameter 1
    Double_t fParCluster1;

    ClassDef(R3BCalifaClusterPar, 2); //
};

#endif /* !R3BCALIFACLUSTERPAR_H*/
