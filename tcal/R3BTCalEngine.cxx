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

#include <TH1F.h>
#include <TMath.h>
#include <string>

#include <FairLogger.h>

#include "R3BLogger.h"
#include "R3BTCalEngine.h"

R3BTCalEngine::R3BTCalEngine(R3BTCalPar* param, Int_t minStats)
    : fMinStats(minStats)
    , fCal_Par(param)
    , fClockFreq(0.)
{
    for (Int_t i = 0; i < N_PLANE_MAX; i++)
    {
        for (Int_t j = 0; j < N_PADDLE_MAX; j++)
        {
            for (Int_t k = 0; k < N_SIDE_MAX; k++)
            {
                fhData[i][j][k] = NULL;
                fhTime[i][j][k] = NULL;
            }
        }
    }
}

void R3BTCalEngine::Fill(Int_t plane, Int_t paddle, Int_t side, Int_t tdc)
{
    if (plane < 1 || plane > N_PLANE_MAX || paddle < 1 || paddle > N_PADDLE_MAX || side < 1 || side > N_SIDE_MAX)
    {
        R3BLOG(error, "index out of max range " << plane << " / " << paddle << " / " << side);
        R3BLOG(error, "ranges: " << N_PLANE_MAX << " / " << N_PADDLE_MAX << " / " << N_SIDE_MAX);
        return;
    }
    if (NULL == fhData[plane - 1][paddle - 1][side - 1])
    {
        char strName[255];
        sprintf(strName, "%s_tcaldata_%d_%d_%d", fCal_Par->GetName(), plane, paddle, side);
        fhData[plane - 1][paddle - 1][side - 1] = new TH1F(strName, "", 4097, -0.5, 4096.5);
        sprintf(strName, "%s_time_%d_%d_%d", fCal_Par->GetName(), plane, paddle, side);
        fhTime[plane - 1][paddle - 1][side - 1] = new TH1F(strName, "", 4097, -0.5, 4096.5);
    }
    fhData[plane - 1][paddle - 1][side - 1]->Fill(tdc);

    //  LOG(info) << "R3BTCalEngine:: " << plane << " " << paddle<< " "<<side<<" "<< fhData[plane - 1][paddle - 1][side
    //  - 1]->GetEntries();
}

void R3BTCalEngine::CalculateParamClockTDC(enum CTDCVariant a_variant)
{
    switch (a_variant)
    {
        case CTDC_8_12_FWD_250:
            fClockFreq = 1000. / 250;
            break;
        case CTDC_16_BWD_150:
            fClockFreq = 1000. / CTDC_16_CLOCK_MHZ;
            break;
        default:
            assert(0 && "Invalid CTDC variant!");
    }

    for (Int_t i = 0; i < N_PLANE_MAX; i++)
    {
        for (Int_t j = 0; j < N_PADDLE_MAX; j++)
        {
            for (Int_t k = 0; k < N_SIDE_MAX; k++)
            {
                if (NULL == fhData[i][j][k])
                {
                    continue;
                }
                if (fhData[i][j][k]->GetEntries() < fMinStats)
                {
                    continue;
                }

                // Define range of channels
                Int_t ic, iMin, iMax;
                FindRange(fhData[i][j][k], ic, iMin, iMax);
                if (iMin < 0 || iMax > 4097)
                {
                    return;
                }
                R3BLOG(info, "Range of channels: " << iMin << " - " << iMax);

                Int_t nparam = 0;
                auto pTCal = new R3BTCalModulePar;
                pTCal->SetPlane(i + 1);
                pTCal->SetPaddle(j + 1);
                pTCal->SetSide(k + 1);

                Double_t total = fhData[i][j][k]->Integral(iMin, iMax);
                for (Int_t ii = iMin; ii < iMax; ii++)
                {
                    auto bin_mid = fhData[i][j][k]->Integral(iMin, ii) + fhData[i][j][k]->GetBinContent(1 + ii) * 0.5;
                    auto time_ns = bin_mid / total * fClockFreq;

                    fhTime[i][j][k]->SetBinContent(1 + ii, time_ns);

                    pTCal->SetBinLowAt(ii, nparam);
                    pTCal->SetOffsetAt(time_ns, nparam);
                    pTCal->IncrementNofChannels();
                    nparam++;
                }
                fCal_Par->AddModulePar(pTCal);

                R3BLOG(info, "Number of parameters: " << nparam);

                fhData[i][j][k]->Write();
                fhTime[i][j][k]->Write();

                R3BLOG(info, "Module: " << (i + 1) << " / " << (j + 1) << " / " << (k + 1) << " is calibrated.");
            }
        }
    }

    fCal_Par->setChanged();
}

void R3BTCalEngine::CalculateParamTacquila()
{
    fClockFreq = 1. / TACQUILA_CLOCK_MHZ * 1000.;

    for (Int_t i = 0; i < N_PLANE_MAX; i++)
    {
        for (Int_t j = 0; j < N_PADDLE_MAX; j++)
        {
            for (Int_t k = 0; k < N_SIDE_MAX; k++)
            {
                if (NULL == fhData[i][j][k])
                {
                    continue;
                }
                if (fhData[i][j][k]->GetEntries() < fMinStats)
                {
                    continue;
                }

                // Define range of channels
                Int_t ic, iMin, iMax;
                FindRange(fhData[i][j][k], ic, iMin, iMax);
                if (iMin < 0 || iMax > 4097)
                {
                    return;
                }
                R3BLOG(info, "Range of channels: " << iMin << " - " << iMax);

                Double_t total = fhData[i][j][k]->Integral(iMin, iMax);
                for (Int_t ii = iMin; ii <= iMax; ii++)
                {
                    fhTime[i][j][k]->SetBinContent(ii, fhData[i][j][k]->Integral(iMin, ii) / total * fClockFreq);
                }

                Int_t nparam = 0;

                R3BTCalModulePar* pTCal = NULL;
                pTCal = new R3BTCalModulePar();
                pTCal->SetPlane(i + 1);
                pTCal->SetPaddle(j + 1);
                pTCal->SetSide(k + 1);

                Int_t il = ic - 10 + 1;
                Int_t ih = ic;
                while (il > iMin)
                {
                    Double_t slope = 0, offset = 0;
                    LinearDown(fhData[i][j][k], iMin, iMax, il, ih, slope, offset);

                    pTCal->SetBinLowAt(il, nparam);
                    pTCal->SetBinUpAt(ih, nparam);
                    pTCal->SetSlopeAt(slope, nparam);
                    pTCal->SetOffsetAt(offset, nparam);
                    pTCal->IncrementNofChannels();

                    nparam += 1;

                    ih = il;
                    il = ih - 10 + 1;
                }

                if (ih > iMin)
                {
                    Double_t slope = 0, offset = 0;
                    Double_t tot = fhData[i][j][k]->Integral(iMin, iMax);
                    Double_t t1 = 0.;
                    Double_t t2 = fhData[i][j][k]->Integral(iMin, ih) / tot * fClockFreq;
                    slope = (t2 - t1) / (Double_t)(ih - iMin);
                    offset = t1;

                    pTCal->SetBinLowAt(iMin, nparam);
                    pTCal->SetBinUpAt(ih, nparam);
                    pTCal->SetSlopeAt(slope, nparam);
                    pTCal->SetOffsetAt(offset, nparam);
                    pTCal->IncrementNofChannels();

                    nparam += 1;
                }

                il = ic;
                ih = ic + 10 - 1;
                while (ih <= iMax)
                {
                    Double_t slope = 0, offset = 0;
                    LinearUp(fhData[i][j][k], iMin, iMax, il, ih, slope, offset);

                    pTCal->SetBinLowAt(il, nparam);
                    pTCal->SetBinUpAt(ih, nparam);
                    pTCal->SetSlopeAt(slope, nparam);
                    pTCal->SetOffsetAt(offset, nparam);
                    pTCal->IncrementNofChannels();

                    nparam += 1;

                    il = ih;
                    if ((iMax - ih) < 100)
                    {
                        ih = il + 5 - 1;
                    }
                    else
                    {
                        ih = il + 10 - 1;
                    }
                }

                if (il < iMax)
                {
                    Double_t slope = 0, offset = 0;
                    Double_t tot = fhData[i][j][k]->Integral(iMin, iMax);
                    Double_t t1 = fhData[i][j][k]->Integral(iMin, il) / tot * fClockFreq;
                    Double_t t2 = fClockFreq;
                    slope = (t2 - t1) / (Double_t)(iMax - il);
                    offset = t1;

                    pTCal->SetBinLowAt(il, nparam);
                    pTCal->SetBinUpAt(iMax, nparam);
                    pTCal->SetSlopeAt(slope, nparam);
                    pTCal->SetOffsetAt(offset, nparam);
                    pTCal->IncrementNofChannels();

                    nparam += 1;
                }

                fCal_Par->AddModulePar(pTCal);

                R3BLOG(info, "Number of parameters: " << nparam);

                fhData[i][j][k]->Write();
                fhTime[i][j][k]->Write();

                R3BLOG(info, "Module: " << (i + 1) << " / " << (j + 1) << " / " << (k + 1) << " is calibrated.");
            }
        }
    }
    fCal_Par->setChanged();
}

void R3BTCalEngine::CalculateParamVFTX()
{
    fClockFreq = 1. / VFTX_CLOCK_MHZ * 1000.;
    R3BLOG(info, "Num Det: " << N_PLANE_MAX << ", Num Paddles: " << N_PADDLE_MAX);

    for (Int_t i = 0; i < N_PLANE_MAX; i++)
    {
        //  LOG(info) << "R3BTCalEngine::CalculateParamVFTX() : Plane: " << i;
        for (Int_t j = 0; j < N_PADDLE_MAX; j++)
        {
            //	 LOG(info) << "R3BTCalEngine::CalculateParamVFTX() : Paddle: " << j;

            for (Int_t k = 0; k < N_SIDE_MAX; k++)
            {

                // if(i == 0)   LOG(info) << "R3BTCalEngine::CalculateParamVFTX() : Detector: "<<i<<", Channel: "<<j<<",
                // Type: " << k;

                if (NULL == fhData[i][j][k])
                {
                    //			LOG(info) << "R3BTCalEngine::CalculateParamVFTX() : NULL: " << fhData[i][j][k];
                    continue;
                }
                if (fhData[i][j][k]->GetEntries() < fMinStats)
                {
                    R3BLOG(info, " fMinStatus: " << fhData[i][j][k]->GetEntries() << ", " << fMinStats);
                    continue;
                }

                // Define range of channels
                Int_t ic, iMin, iMax;
                FindRange(fhData[i][j][k], ic, iMin, iMax);
                if (iMin < 0 || iMax > 4097)
                {
                    R3BLOG(warn, "Out of range, Min: " << iMin << " , Max: " << iMax);
                    return;
                }
                R3BLOG(info, "Range of channels: " << iMin << " - " << iMax);

                Double_t total = fhData[i][j][k]->Integral(iMin, iMax);
                for (Int_t ii = iMin; ii <= iMax; ii++)
                {
                    fhTime[i][j][k]->SetBinContent(ii, fhData[i][j][k]->Integral(iMin, ii) / total * fClockFreq);
                }

                Int_t nparam = 0;

                R3BTCalModulePar* pTCal = NULL;
                pTCal = new R3BTCalModulePar();
                pTCal->SetPlane(i + 1);
                pTCal->SetPaddle(j + 1);
                pTCal->SetSide(k + 1);

                for (Int_t ibin = iMin; ibin <= iMax; ibin++)
                {
                    Double_t time = fhData[i][j][k]->Integral(iMin, ibin) / total;
                    if (time > 1.)
                    {
                        LOG(fatal) << "Integration error.";
                    }
                    time *= fClockFreq;

                    pTCal->SetBinLowAt(ibin, nparam);
                    pTCal->SetOffsetAt(time, nparam);
                    pTCal->IncrementNofChannels();
                    nparam += 1;
                }

                fCal_Par->AddModulePar(pTCal);

                R3BLOG(info, "Number of parameters: " << nparam);

                fhData[i][j][k]->Write();
                fhTime[i][j][k]->Write();

                R3BLOG(info, "Module: " << (i + 1) << " / " << (j + 1) << " / " << (k + 1) << " is calibrated.");
            }
        }
    }
    fCal_Par->setChanged();
}

// iMin == left side of fine times.
// iMax == right side of fine times.
// I.e. iMin <= fine-time <= iMax-1.
void R3BTCalEngine::FindRange(TH1F* h1, Int_t& ic, Int_t& iMin, Int_t& iMax)
{
    Double_t mean = h1->GetMean();
    ic = (Int_t)(mean + 0.5);

    for (Int_t i = 1; i <= 4097; i++)
    {
        if (h1->GetBinContent(i) > 0)
        {
            iMin = i - 1;
            break;
        }
    }

    for (Int_t i = 4097; i >= 1; i--)
    {
        if (h1->GetBinContent(i) > 0)
        {
            iMax = i;
            break;
        }
    }
}

void R3BTCalEngine::LinearUp(TH1F* h1, Int_t iMin, Int_t iMax, Int_t& il, Int_t& ih, Double_t& slope, Double_t& offset)
{
    Double_t tot = h1->Integral(iMin, iMax);
    Double_t t1 = h1->Integral(iMin, il) / tot; // * fClockFreq;
    Double_t t2 = h1->Integral(iMin, ih) / tot; // * fClockFreq;
    if (t1 > 1. || t2 > 1.)
    {
        R3BLOG(fatal, "Integration error");
    }
    t1 *= fClockFreq;
    t2 *= fClockFreq;
    slope = (t2 - t1) / (Double_t)(ih - il);
    offset = t1;

    Double_t prec = 3. / TMath::Sqrt(h1->GetEntries());

    Double_t slope1;

    Int_t iter = 1;
    while (ih <= iMax)
    {
        Int_t ih_next = ih + iter;
        if (ih_next > iMax)
        {
            break;
        }
        Double_t t21 = h1->Integral(iMin, ih_next) / tot * fClockFreq;
        slope1 = (t21 - t1) / (Double_t)(ih_next - il);

        Double_t dev = TMath::Abs(slope1 - slope) / TMath::Abs(slope);

        ih = ih_next;
        iter += 1;

        slope = slope1;
        if (dev > prec)
        {
            break;
        }
    }
}

void R3BTCalEngine::LinearDown(TH1F* h1,
                               Int_t iMin,
                               Int_t iMax,
                               Int_t& il,
                               Int_t& ih,
                               Double_t& slope,
                               Double_t& offset)
{
    Double_t tot = h1->Integral(iMin, iMax);
    Double_t t1 = h1->Integral(iMin, il) / tot * fClockFreq;
    Double_t t2 = h1->Integral(iMin, ih) / tot * fClockFreq;
    slope = (t2 - t1) / (Double_t)(ih - il);
    offset = t1;

    Double_t prec = 3. / TMath::Sqrt(h1->GetEntries());

    Double_t slope1;
    Double_t offset1;

    Int_t iter = 1;
    while (il >= iMin)
    {
        Int_t il_next = il - iter, ih_next = ih;
        if (il_next < iMin)
        {
            break;
        }
        Double_t t11 = h1->Integral(iMin, il_next) / tot * fClockFreq;
        Double_t t21 = h1->Integral(iMin, ih_next) / tot * fClockFreq;
        slope1 = (t21 - t11) / (Double_t)(ih_next - il_next);
        offset1 = t11;

        Double_t dev = TMath::Abs(slope1 - slope) / TMath::Abs(slope);

        il = il_next;
        ih = ih_next;
        iter += 1;

        slope = slope1;
        offset = offset1;
        if (dev > prec)
        {
            break;
        }
    }
}

ClassImp(R3BTCalEngine)
