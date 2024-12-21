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

#include "R3BNeulandNeutronReconstructionStatistics.h"
#include "FairLogger.h"
#include <FairRootManager.h>
#include <TFile.h>
#include <iostream>
#include <numeric>

template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type almost_equal(T x, T y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           // unless the result is subnormal
           || std::abs(x - y) < std::numeric_limits<T>::min();
}

R3BNeulandNeutronReconstructionStatistics::R3BNeulandNeutronReconstructionStatistics(const TString primary,
                                                                                     const TString secondary,
                                                                                     const TString predicted,
                                                                                     std::ostream& out)
    : fPrimaryClusters(primary)
    , fSecondaryClusters(secondary)
    , fPredictedNeutrons(predicted)
    , fPredictedName(predicted)
    , fhTP(nullptr)
    , fhFP(nullptr)
    , fhFN(nullptr)
    , fhTN(nullptr)
    , fTP(0)
    , fFP(0)
    , fFN(0)
    , fTN(0)
    , fOut(out)
    , fMult(100)
{
}

InitStatus R3BNeulandNeutronReconstructionStatistics::Init()
{
    try
    {
        fPrimaryClusters.Init();
        fSecondaryClusters.Init();
        fPredictedNeutrons.Init();
    }
    catch (const std::exception& e)
    {
        LOG(fatal) << "R3BNeulandNeutronReconstruction" << e.what();
    }

    TH1::AddDirectory(kFALSE);
    fhTP = new TH1D("fhTP", "True Positive", 10, 0., 10.);
    fhFN = new TH1D("fhFN", "False Negative", 10, 0., 10.);
    fhFP = new TH1D("fhFP", "False Positive", 10, 0., 10.);
    fhTN = new TH1D("fhTN", "True Negative", 30, 0., 30.);
    fhF1 = new TH1D("fhF1", "F1 Value", 11, 0., 1.1);

    return kSUCCESS;
}

void R3BNeulandNeutronReconstructionStatistics::Exec(Option_t*)
{
    const auto actualPositive = fPrimaryClusters.Retrieve();
    const auto actualNegative = fSecondaryClusters.Retrieve();
    const auto predictedPositive = fPredictedNeutrons.Retrieve();

    auto comp = [](const R3BNeulandNeutron* n, const R3BNeulandCluster* c)
    {
        return almost_equal(c->GetT(), n->GetT(), 2) && almost_equal(c->GetPosition().X(), n->GetPosition().X(), 2) &&
               almost_equal(c->GetPosition().Y(), n->GetPosition().Y(), 2) &&
               almost_equal(c->GetPosition().Z(), n->GetPosition().Z(), 2);
    };

    fMult[predictedPositive.size()]++;

    double truePositives = 0.;
    for (const auto& c : actualPositive)
    {
        for (const auto& n : predictedPositive)
        {
            if (comp(n, c))
            {
                truePositives++;
                break;
            }
        }
    }
    const double falseNegatives = (double)actualPositive.size() - truePositives;
    fTP += truePositives;
    fhTP->Fill(truePositives);
    fFN += falseNegatives;
    fhFN->Fill(falseNegatives);

    double falsePositives = 0.;
    for (const auto& c : actualNegative)
    {
        for (const auto& n : predictedPositive)
        {
            if (comp(n, c))
            {
                falsePositives++;
                break;
            }
        }
    }
    const double trueNegatives = (double)actualNegative.size() - falsePositives;
    fFP += falsePositives;
    fhFP->Fill(falsePositives);
    fTN += trueNegatives;
    fhTN->Fill(trueNegatives);

    const double precision = truePositives / (truePositives + falsePositives);
    const double recall = truePositives / (truePositives + falseNegatives);

    auto Fbeta = [&](double beta)
    {
        if (precision == 0. && recall == 0.)
            return 0.;
        return (1. + beta * beta) * (precision * recall) / ((beta * beta * precision) + recall);
    };
    fhF1->Fill(Fbeta(1.));
}

void R3BNeulandNeutronReconstructionStatistics::Finish()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();

    gDirectory->mkdir(fPredictedName + "Statistics");
    gDirectory->cd(fPredictedName + "Statistics");

    fhTP->Write();
    fhFN->Write();
    fhFP->Write();
    fhTN->Write();
    fhF1->Write();

    gDirectory = tmp;

    const double precision = (double)fTP / ((double)fTP + (double)fFP);
    const double recall = (double)fTP / ((double)fTP + (double)fFN);
    auto Fbeta = [&](double beta)
    { return (1. + beta * beta) * (precision * recall) / ((beta * beta * precision) + recall); };
    const double accuracy = (double)(fTP + fTN) / (double)(fTP + fTN + fFP + fFN);

    fOut << "PREC DATA"
         << "\n";
    fOut << fPredictedName << "\tTruePositive \t" << fTP << "\t" << (double)fTP / ((double)fTP + (double)fFN) << "\n"
         << fPredictedName << "\tFalsePositive\t" << fFP << "\t" << (double)fFP / ((double)fTN + (double)fFP) << "\n"
         << fPredictedName << "\tFalseNegative\t" << fFN << "\t" << (double)fFN / ((double)fTP + (double)fFN) << "\n"
         << fPredictedName << "\tTrueNegative \t" << fTN << "\t" << (double)fTN / ((double)fTN + (double)fFP) << "\n"
         << fPredictedName << "\tAccuracy\t" << accuracy << "\n"
         << fPredictedName << "\tPrecision\t" << precision << "\n"
         << fPredictedName << "\tSensitivity\t" << recall << "\n"
         << fPredictedName << "\tF1\t" << Fbeta(1.) << "\n"
         << "\n";

    fOut << "MULT DATA"
         << "\n";
    const int sum = std::accumulate(fMult.begin(), fMult.end(), 0);
    for (unsigned int i = 0; i < fMult.size(); i++)
    {
        fOut << i << "\t" << fMult[i] << "\t" << (double)fMult[i] / (double)sum << "\n";
    }
    fOut << "END" << std::endl;
}
