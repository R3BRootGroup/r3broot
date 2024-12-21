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

#ifndef R3BMDFWRAPPER_H
#define R3BMDFWRAPPER_H

#include "FairLogger.h"
#include "R3BLogger.h"
#include "Rtypes.h"
#include "TMath.h"
#include "TObject.h"

using namespace std;

class R3BMDFWrapper : public TObject
{
  public:
    R3BMDFWrapper();                     // default constructor
    R3BMDFWrapper(const char* mdf_file); // standard constructor
    virtual ~R3BMDFWrapper();
    //~R3BMDFWrapper() override=default;

    void InitMDF(const char* mdf_file);
    void InitPCA(const char* pca_file);
    void PrintMDF();
    void PrintPCA();
    void X2P(Double_t* x, Double_t* p);
    void P2X(Double_t* p, Double_t* x, Int_t nTest);
    int GetNVariables() { return mdf_NVariables; }
    Double_t MDF(Double_t* x);

  private:
    // MDF data
    int mdf_NVariables;
    int mdf_NCoefficients;
    int mdf_NMaxFunctions;
    double mdf_DMean;
    int mdf_PolyType;
    double* mdf_XMean;
    double* mdf_XMin;
    double* mdf_XMax;
    double* mdf_Coefficient;
    double* mdf_CoefficientRMS;
    int* mdf_Power;

    // PCA data
    int pca_NVariables;
    double* pca_EigenVectors;
    double* pca_EigenValues;
    double* pca_MeanValues;
    double* pca_SigmaValues;

    bool is_PCA;

    ClassDef(R3BMDFWrapper, 1);
};

#endif
