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
#include "R3BMDFWrapper.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

R3BMDFWrapper::R3BMDFWrapper()
    : mdf_NVariables(0)
    , pca_NVariables(0)
    , is_PCA(false)
{
}

R3BMDFWrapper::R3BMDFWrapper(const char* mdf_file)
{
    InitMDF(mdf_file);
    is_PCA = false;
}

R3BMDFWrapper::~R3BMDFWrapper() {}

// This method will read a special txt file with MDF parameters, don't change
void R3BMDFWrapper::InitMDF(const char* mdf_file)
{
    int i;
    std::ifstream file(mdf_file);
    if (file.is_open())
    {
        file >> mdf_NVariables >> mdf_NCoefficients >> mdf_NMaxFunctions >> mdf_DMean >> mdf_PolyType;

        mdf_XMean = new double[mdf_NVariables];
        for (i = 0; i < mdf_NVariables; i++)
            file >> mdf_XMean[i];

        mdf_XMin = new double[mdf_NVariables];
        for (i = 0; i < mdf_NVariables; i++)
            file >> mdf_XMin[i];

        mdf_XMax = new double[mdf_NVariables];
        for (i = 0; i < mdf_NVariables; i++)
            file >> mdf_XMax[i];

        mdf_Coefficient = new double[mdf_NCoefficients];
        for (i = 0; i < mdf_NCoefficients; i++)
            file >> mdf_Coefficient[i];

        mdf_CoefficientRMS = new double[mdf_NCoefficients];
        for (i = 0; i < mdf_NCoefficients; i++)
            file >> mdf_CoefficientRMS[i];

        mdf_Power = new int[mdf_NCoefficients * mdf_NVariables];
        for (i = 0; i < mdf_NCoefficients * mdf_NVariables; i++)
            file >> mdf_Power[i];
        file.close();
    }

    else
        R3BLOG(fatal, "\nUnable to open MDF parameter file!!");
    return;
}

void R3BMDFWrapper::PrintMDF()
{
    R3BLOG(info, "\n\n======== Initialized MDF from the R3BMDFWrapper:");
    cout << std::setprecision(9) << "\nmdf_NVariables = " << mdf_NVariables
         << "\nmdf_NCoefficients = " << mdf_NCoefficients << "\nmdf_NMaxFunctions = " << mdf_NMaxFunctions
         << "\nmdf_DMean = " << mdf_DMean << "\nmdf_PolyType = " << mdf_PolyType;

    cout << "\nmdf_Xmean = ";
    for (int i = 0; i < mdf_NVariables; i++)
        cout << mdf_XMean[i] << ", ";

    cout << "\nmdf_Xmin = ";
    for (int i = 0; i < mdf_NVariables; i++)
        cout << mdf_XMin[i] << ", ";

    cout << "\nmdf_Xmax = ";
    for (int i = 0; i < mdf_NVariables; i++)
        cout << mdf_XMax[i] << ", ";

    cout << "\nmdf_Coefficient = \n";
    for (int i = 0; i < mdf_NCoefficients; i++)
        cout << mdf_Coefficient[i] << ", ";

    cout << "\nmdf_CoefficientRMS = \n";
    for (int i = 0; i < mdf_NCoefficients; i++)
        cout << mdf_CoefficientRMS[i] << ", ";

    cout << "\nmdf_Power = \n";
    for (int i = 0; i < mdf_NCoefficients * mdf_NVariables; i++)
        cout << mdf_Power[i] << ", ";
}

// This method will read a special txt file with PCA parameter, don't change
void R3BMDFWrapper::InitPCA(const char* pca_file)
{
    int i;
    std::ifstream file(pca_file);
    if (file.is_open())
    {
        file >> pca_NVariables;

        pca_EigenVectors = new double[pca_NVariables * pca_NVariables];
        for (i = 0; i < pca_NVariables * pca_NVariables; i++)
            file >> pca_EigenVectors[i];

        pca_EigenValues = new double[pca_NVariables];
        for (i = 0; i < pca_NVariables; i++)
            file >> pca_EigenValues[i];

        pca_MeanValues = new double[pca_NVariables];
        for (i = 0; i < pca_NVariables; i++)
            file >> pca_MeanValues[i];

        pca_SigmaValues = new double[pca_NVariables];
        for (i = 0; i < pca_NVariables; i++)
            file >> pca_SigmaValues[i];

        file.close();
    }
    else
        R3BLOG(fatal, "\nUnable to open PCA parameter file!!");
    return;
}

void R3BMDFWrapper::X2P(Double_t* x, Double_t* p)
{
    for (Int_t i = 0; i < pca_NVariables; i++)
    {
        p[i] = 0;
        for (Int_t j = 0; j < pca_NVariables; j++)
            p[i] += (x[j] - pca_MeanValues[j]) * pca_EigenVectors[j * pca_NVariables + i] / pca_SigmaValues[j];
    }

    return;
}

void R3BMDFWrapper::P2X(Double_t* p, Double_t* x, Int_t nTest)
{
    for (Int_t i = 0; i < pca_NVariables; i++)
    {
        x[i] = pca_MeanValues[i];
        for (Int_t j = 0; j < nTest; j++)
            x[i] += p[j] * pca_SigmaValues[i] * pca_EigenVectors[i * pca_NVariables + j];
    }

    return;
}

Double_t R3BMDFWrapper::MDF(Double_t* x)
{
    Double_t returnValue = mdf_DMean;
    int i = 0, j = 0, k = 0;
    for (i = 0; i < mdf_NCoefficients; i++)
    {
        // Evaluate the ith term in the expansion
        double term = mdf_Coefficient[i];
        for (j = 0; j < mdf_NVariables; j++)
        {
            // Evaluate the polynomial in the jth variable.
            int power = mdf_Power[mdf_NVariables * i + j];
            double p1 = 1, p2 = 0, p3 = 0, r = 0;
            double v = 1 + 2. / (mdf_XMax[j] - mdf_XMin[j]) * (x[j] - mdf_XMax[j]);
            // what is the power to use!
            switch (power)
            {
                case 1:
                    r = 1;
                    break;
                case 2:
                    r = v;
                    break;
                default:
                    p2 = v;

                    for (k = 3; k <= power; k++)
                    {
                        p3 = p2 * v;

                        if (mdf_PolyType == 2) // kLegendre
                            p3 = ((2 * i - 3) * p2 * v - (i - 2) * p1) / (i - 1);

                        if (mdf_PolyType == 1) // kChebyshev
                            p3 = 2 * v * p2 - p1;

                        p1 = p2;
                        p2 = p3;
                    }
                    r = p3;
            }
            // multiply this term by the poly in the jth var
            term *= r;
        }
        // Add this term to the final result
        returnValue += term;
    }
    return returnValue;
}

ClassImp(R3BMDFWrapper)
