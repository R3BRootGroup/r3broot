/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ------------------------------------------------------------
// -----                  R3BLosCal2Hit                   -----
// -----          Created Mar 10th 2016 by R.Plag         -----
// ------------------------------------------------------------

#include "R3BLosCal2Hit.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BLosHitPar.h"
#include "R3BLosMapped2Cal.h"
#include "R3BLosMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <FairRootManager.h>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <cstdlib>
#include <iterator>

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

R3BLosCal2Hit::R3BLosCal2Hit()
    : R3BLosCal2Hit("LosCal2Hit", 1)
{
}

R3BLosCal2Hit::R3BLosCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(NULL)
    , fHitItems(new TClonesArray("R3BLosHitData"))
    , fNofHitItems(0)
    , fTrigger(-1)
    , fTpat(-1)
    , flosVeffX(1.)
    , flosVeffY(1.)
    , flosOffsetX(0.)
    , flosOffsetY(0.)
    , flosVeffXT(1.)
    , flosVeffYT(1.)
    , flosOffsetXT(0.)
    , flosOffsetYT(0.)
    , flosVeffXQ(1.)
    , flosVeffYQ(1.)
    , flosOffsetXQ(0.)
    , flosOffsetYQ(0.)
    , fOptHisto(kFALSE)
    , fOnline(kFALSE)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{
    fhTres_M = NULL;
    fhTres_T = NULL;
    fhTres_T_corr = NULL;
    fhTres_M_corr = NULL;
    fhQ_L = NULL;
    fhQ_B = NULL;
    fhQ_R = NULL;
    fhQ_T = NULL;
    fhQ_L_corr = NULL;
    fhQ_B_corr = NULL;
    fhQ_R_corr = NULL;
    fhQ_T_corr = NULL;
    fhQ_LT = NULL;
    fhQ_LB = NULL;
    fhQ_RB = NULL;
    fhQ_RT = NULL;
    fhQ_LT_corr = NULL;
    fhQ_LB_corr = NULL;
    fhQ_RB_corr = NULL;
    fhQ_RT_corr = NULL;
    fhQ = NULL;
    fhQ_vs_X = NULL;
    fhQ_vs_Y = NULL;
    fhQ_vs_X_corr = NULL;
    fhQ_vs_Y_corr = NULL;
    fhTM_vs_Q = NULL;
    fhTT_vs_Q = NULL;
    fhTM_vs_Q_corr = NULL;
    fhTT_vs_Q_corr = NULL;
    fhXY = NULL;
    fhXYmean = NULL;
    fhXY_ToT = NULL;
    fhXYproj = NULL;
    fhXYT = NULL;
    fhQ1_vs_Q5 = NULL;
    fhQ1_vs_Q5_corr = NULL;
    fhQ2_vs_Q6 = NULL;
    fhQ2_vs_Q6_corr = NULL;
    fhQ3_vs_Q7 = NULL;
    fhQ3_vs_Q7_corr = NULL;
    fhQ4_vs_Q8 = NULL;
    fhQ4_vs_Q8_corr = NULL;
    fhQtest = NULL;
    fhTresX_M = NULL;
    fhTresY_M = NULL;
    fhTresX_M_corr = NULL;
    fhTresY_M_corr = NULL;
    fhTresX_T = NULL;
    fhTresY_T = NULL;
    fhTresX_T_corr = NULL;
    fhTresY_T_corr = NULL;
    fhTresMvsIcount = NULL;
    fhTresTvsIcount = NULL;
    fhTreswcMvsIcount = NULL;
    fhTreswcTvsIcount = NULL;
    fh_los_dt_hits_ToT_corr = NULL;
    fh_los_ihit_ToTcorr = NULL;
    for (Int_t j = 0; j < 8; j++)
    {
        fhQvsdt[j] = NULL;
        fhQcorrvsIcount[j] = NULL;
        fhQvsIcount[j] = NULL;
    }
}

R3BLosCal2Hit::~R3BLosCal2Hit()
{
    LOG(debug) << "R3BLosCal2Hit::Destructor";
    if (fhTres_M)
        delete (fhTres_M);
    if (fhTres_T)
        delete (fhTres_T);
    if (fhTres_T_corr)
        delete (fhTres_T_corr);
    if (fhTres_M_corr)
        delete (fhTres_M_corr);
    if (fhQ_L)
        delete (fhQ_L);
    if (fhQ_B)
        delete (fhQ_B);
    if (fhQ_R)
        delete (fhQ_R);
    if (fhQ_T)
        delete (fhQ_T);
    if (fhQ_L_corr)
        delete (fhQ_L_corr);
    if (fhQ_B_corr)
        delete (fhQ_B_corr);
    if (fhQ_R_corr)
        delete (fhQ_R_corr);
    if (fhQ_T_corr)
        delete (fhQ_T_corr);
    if (fhQ_LT)
        delete (fhQ_LT);
    if (fhQ_LB)
        delete (fhQ_LB);
    if (fhQ_RB)
        delete (fhQ_RB);
    if (fhQ_RT)
        delete (fhQ_RT);
    if (fhQ_LT_corr)
        delete (fhQ_LT_corr);
    if (fhQ_LB_corr)
        delete (fhQ_LB_corr);
    if (fhQ_RB_corr)
        delete (fhQ_RB_corr);
    if (fhQ_RT_corr)
        delete (fhQ_RT_corr);
    if (fhQ)
        delete (fhQ);
    if (fhQtest)
        delete (fhQtest);
    if (fhQ_vs_X)
        delete (fhQ_vs_X);
    if (fhQ_vs_Y)
        delete (fhQ_vs_Y);
    if (fhQ_vs_X_corr)
        delete (fhQ_vs_X_corr);
    if (fhQ_vs_Y_corr)
        delete (fhQ_vs_Y_corr);
    if (fhTM_vs_Q)
        delete (fhTM_vs_Q);
    if (fhTT_vs_Q)
        delete (fhTT_vs_Q);
    if (fhTM_vs_Q_corr)
        delete (fhTM_vs_Q_corr);
    if (fhTT_vs_Q_corr)
        delete (fhTT_vs_Q_corr);
    if (fhXY)
        delete (fhXY);
    if (fhXYmean)
        delete (fhXYmean);
    if (fhXY_ToT)
        delete (fhXY_ToT);
    if (fhXYproj)
        delete (fhXYproj);
    if (fhXYT)
        delete (fhXYT);
    if (fhQ1_vs_Q5)
        delete (fhQ1_vs_Q5);
    if (fhQ1_vs_Q5_corr)
        delete (fhQ1_vs_Q5_corr);
    if (fhQ2_vs_Q6)
        delete (fhQ2_vs_Q6);
    if (fhQ2_vs_Q6_corr)
        delete (fhQ2_vs_Q6_corr);
    if (fhQ3_vs_Q7)
        delete (fhQ3_vs_Q7);
    if (fhQ3_vs_Q7_corr)
        delete (fhQ3_vs_Q7_corr);
    if (fhQ4_vs_Q8)
        delete (fhQ4_vs_Q8);
    if (fhQ4_vs_Q8_corr)
        delete (fhQ4_vs_Q8_corr);
    if (fhTresX_M)
        delete (fhTresX_M);
    if (fhTresY_M)
        delete (fhTresY_M);
    if (fhTresX_M_corr)
        delete (fhTresX_M_corr);
    if (fhTresY_M_corr)
        delete (fhTresY_M_corr);
    if (fhTresX_T)
        delete (fhTresX_T);
    if (fhTresY_T)
        delete (fhTresY_T);
    if (fhTresX_T_corr)
        delete (fhTresX_T_corr);
    if (fhTresY_T_corr)
        delete (fhTresY_T_corr);
    if (fh_los_ihit_ToTcorr)
        delete (fh_los_ihit_ToTcorr);
    if (fhTresMvsIcount)
        delete (fhTresMvsIcount);
    if (fhTresTvsIcount)
        delete (fhTresTvsIcount);
    if (fhTreswcMvsIcount)
        delete (fhTreswcMvsIcount);
    if (fhTreswcTvsIcount)
        delete (fhTreswcTvsIcount);
    if (fh_los_dt_hits_ToT_corr)
        delete (fh_los_dt_hits_ToT_corr);
    for (Int_t j = 0; j < 8; j++)
    {
        if (fhQvsdt[j])
            delete (fhQvsdt[j]);
        if (fhQcorrvsIcount[j])
            delete (fhQcorrvsIcount[j]);
        if (fhQvsIcount[j])
            delete (fhQvsIcount[j]);
    }

    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

void R3BLosCal2Hit::SetParContainers()
{
    LOG(info) << "R3BLosTcal2Hit::SetParContainers()";
    // Parameter Container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }
    fLosHit_Par = dynamic_cast<R3BLosHitPar*>(rtdb->getContainer("LosHitPar"));
    if (!fLosHit_Par)
    {
        LOG(error) << "R3BLosTcal2Hit:: Couldn't get handle on R3BLosHitPar container";
    }
    else
    {
        LOG(info) << "R3BLosTcal2Hit:: R3BLosHitPar container open";
    }
}

void R3BLosCal2Hit::SetParameter()
{
    //--- Parameter Container ---
    fp0 = fLosHit_Par->Getp0();
    fp1 = fLosHit_Par->Getp1();

    flosOffsetX = fLosHit_Par->Getxoffset_MCFD();
    flosOffsetY = fLosHit_Par->Getyoffset_MCFD();
    flosVeffX = fLosHit_Par->Getxveff_MCFD();
    flosVeffY = fLosHit_Par->Getyveff_MCFD();
}

InitStatus R3BLosCal2Hit::Init()
{
    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(error) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("LosCal"));
    if (NULL == fCalItems)
    {
        LOG(fatal) << "Branch LosCal not found";
        return kFATAL;
    }

    // request storage of Hit data in output tree
    mgr->Register("LosHit", "LosHitData", fHitItems, !fOnline);

    fHitItems->Clear();

    Icount = 0;

    // file with walk-correction parameters
    ifstream infile(fwalk_param_file.c_str());
    if (infile.is_open())
    {
        for (Int_t ivec = 0; ivec < 16; ivec++)
        {
            /*
             ************* Parameters 0-7 MCFD  **********************
             ************* Parameters 8-15 TAMEX *********************
             */

            infile >> walk_par[ivec][0] >> walk_par[ivec][1] >> walk_par[ivec][2] >> walk_par[ivec][3] >>
                walk_par[ivec][4] >> walk_par[ivec][5] >> walk_par[ivec][6] >> walk_par[ivec][7] >> walk_par[ivec][8] >>
                walk_par[ivec][9] >> walk_par[ivec][10];

            /*
                 cout<<setprecision(10)<< ivec<<", "<< walk_par[ivec][0] <<", "<< walk_par[ivec][1] <<", "<<
               walk_par[ivec][2] <<", "<< walk_par[ivec][3] <<", "<< walk_par[ivec][4] <<", "<< walk_par[ivec][5] <<",
               "<< walk_par[ivec][6] <<", "<< walk_par[ivec][7] <<", "<< walk_par[ivec][8] <<", "<< walk_par[ivec][9]
               <<", "<< walk_par[ivec][10] <<endl;
             */
        }
    }
    else // cerr << "Unable to open file \""<<fwalk_param_file<<"\" with walk parameters!"<<endl;
    {
        cout << "*****************************************************************" << endl;
        cout << "UNABLE TO OPEN FILE WITH WALK PARAMETERS! Parameters set to zero!" << endl;
        cout << "*****************************************************************" << endl;
        for (Int_t ivec = 0; ivec < 16; ivec++)
        {
            walk_par[ivec][0] = 10.;
            walk_par[ivec][1] = 1000.;
            for (Int_t iv = 2; iv < 11; iv++)
            {
                walk_par[ivec][iv] = 0.;
                ;
            }
        }
    }
    // file with tot-correction parameters
    ifstream infile1(ftot_param_file.c_str());

    if (infile1.is_open())
    {
        for (Int_t ivec = 0; ivec < 8; ivec++)
        {
            infile1 >> tot_par[ivec][0] >> tot_par[ivec][1] >> tot_par[ivec][2] >> tot_par[ivec][3];
        }
    }
    else // cerr << "Unable to open file \""<<ftot_param_file<<"\" with tot parameters!"<<endl;
    {
        cout << "****************************************************************" << endl;
        cout << "UNABLE TO OPEN FILE WITH ToT PARAMETERS! Parameters set to zero!" << endl;
        cout << "****************************************************************" << endl;
        for (Int_t ivec = 0; ivec < 8; ivec++)
        {
            tot_par[ivec][0] = 0.;
            tot_par[ivec][1] = 0.;
            tot_par[ivec][2] = 0.;
            tot_par[ivec][3] = 1.; // Normalization factor
        }
    }

    // if fOptHisto = true histograms are created
    if (fOptHisto)
        CreateHisto();

    // cout << "R3BLosCal2Hit::Init END" << endl;
    SetParameter();
    return kSUCCESS;
}

InitStatus R3BLosCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

/* Calculate a single hit time for each LOS detector
 *
 * Remember: The times of individual channels depend on the position of
 * the particle on the scintillator. To obtain a precise time of the
 * particle, we need to average either over all four signals (right, top,
 * left, bottom) or over two opposite signals (left+right or top+bottom).
 */
void R3BLosCal2Hit::Exec(Option_t* option)
{
    // cout << "R3BLosCal2Hit::Exec BEGIN: " << Icount << endl;

    // ofstream myFile("data_s473_run197.dat",ios_base::out|ios_base::app);

    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = header->GetTpat();
        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if ((header) && (tpatvalue == 0))
            return;
    }

    Int_t nHits = fCalItems->GetEntries();
    if (nHits == 0)
        return;

    // missing times are NAN, hence other times will also
    // be NAN if one time is missing.
    //  cout<<nHits<<"; "<<narr<<endl;

    Double_t time_V[nHits][8]; // [multihit][pm]
    Double_t time_V_temp[nHits][8];
    Double_t time_L[nHits][8];
    Double_t time_T[nHits][8];
    Double_t time_V_corr[nHits][8]; // [multihit][pm]
    Double_t time_L_corr[nHits][8];
    Double_t timeLosM[nHits];
    Double_t LosTresM[nHits];
    Double_t timeLosT[nHits];
    Double_t LosTresT[nHits];
    Double_t timeLos[nHits];
    Double_t timeLosM_corr[nHits];
    Double_t LosTresM_corr[nHits];
    Double_t timeLosT_corr[nHits];
    Double_t LosTresT_corr[nHits];
    Double_t timeLos_corr[nHits];
    Double_t totsum[nHits];
    Double_t tot[nHits][8];
    Double_t totsum_corr[nHits];
    Double_t tot_corr[nHits][8];
    Double_t xT_cm[nHits];
    Double_t yT_cm[nHits];
    Double_t xV_cm[nHits];
    Double_t yV_cm[nHits];
    Double_t xToT_cm[nHits];
    Double_t yToT_cm[nHits];
    Double_t x_cm[nHits];
    Double_t y_cm[nHits];
    Double_t Z[nHits];
    Double_t t_hit[nHits]; // NAN

    memset(time_V, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_V_temp, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_L, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_T, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_V_corr, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(time_L_corr, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(timeLosM, 0., nHits * sizeof(Double_t));
    memset(LosTresM, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLosT, 0., nHits * sizeof(Double_t));
    memset(LosTresT, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLos, 0., nHits * sizeof(Double_t));
    memset(timeLosM_corr, 0., nHits * sizeof(Double_t));
    memset(LosTresM_corr, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLosT_corr, 0., nHits * sizeof(Double_t));
    memset(LosTresT_corr, 0. / 0., nHits * sizeof(Double_t));
    memset(timeLos_corr, 0., nHits * sizeof(Double_t));
    memset(totsum, 0., nHits * sizeof(Double_t));
    memset(tot, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(totsum_corr, 0., nHits * sizeof(Double_t));
    memset(tot_corr, 0. / 0., nHits * 8 * sizeof(Double_t));
    memset(xT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(yT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(xV_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(yV_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(xToT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(yToT_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(x_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(y_cm, 0. / 0., nHits * sizeof(Double_t));
    memset(Z, 0., nHits * sizeof(Double_t));
    memset(t_hit, 0. / 0., nHits * sizeof(Double_t));

    Int_t nDet = 0;

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {

        auto calItem = dynamic_cast<R3BLosCalData*>(fCalItems->At(ihit));

        nDet = calItem->GetDetector();

        // lt=0, l=1,lb=2,b=3,rb=4,r=5,rt=6,t=7
        for (Int_t iCha = 0; iCha < 8; iCha++)
        {
            time_V[ihit][iCha] = 0. / 0.;
            if (!(IS_NAN(calItem->GetTimeV_ns(iCha))))
            { // VFTX
                time_V[ihit][iCha] = calItem->GetTimeV_ns(iCha);
            }
            time_L[ihit][iCha] = 0. / 0.;
            if (!(IS_NAN(calItem->GetTimeL_ns(iCha))))
            { // TAMEX leading
                time_L[ihit][iCha] = calItem->GetTimeL_ns(iCha);
            }
            time_T[ihit][iCha] = 0. / 0.;
            if (!(IS_NAN(calItem->GetTimeT_ns(iCha))))
            { // TAMEX trailing
                time_T[ihit][iCha] = calItem->GetTimeT_ns(iCha);
            }
        }
    }

    // Sorting VFTX data:

    std::qsort(time_V,
               nHits,
               sizeof(*time_V),
               [](const void* arg1, const void* arg2) -> int
               {
                   double const* lhs = static_cast<double const*>(arg1);
                   double const* rhs = static_cast<double const*>(arg2);

                   return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
               });
    // End sorting
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        Bool_t iLOSTypeMCFD = false;
        Bool_t iLOSTypeTAMEX = false;
        Bool_t iLOSType = false;
        UInt_t Igood_event = 0;

        if (time_V[ihit][0] > 0. && !(IS_NAN(time_V[ihit][0])) && time_V[ihit][1] > 0. && !(IS_NAN(time_V[ihit][1])) &&
            time_V[ihit][2] > 0. && !(IS_NAN(time_V[ihit][2])) && time_V[ihit][3] > 0. && !(IS_NAN(time_V[ihit][3])) &&
            time_V[ihit][4] > 0. && !(IS_NAN(time_V[ihit][4])) && time_V[ihit][5] > 0. && !(IS_NAN(time_V[ihit][5])) &&
            time_V[ihit][6] > 0. && !(IS_NAN(time_V[ihit][6])) && time_V[ihit][7] > 0. && !(IS_NAN(time_V[ihit][7])))
        {
            iLOSTypeMCFD = true; // all 8 MCFD times
        }

        if (time_L[ihit][0] > 0. && !(IS_NAN(time_L[ihit][0])) && time_L[ihit][1] > 0. && !(IS_NAN(time_L[ihit][1])) &&
            time_L[ihit][2] > 0. && !(IS_NAN(time_L[ihit][2])) && time_L[ihit][3] > 0. && !(IS_NAN(time_L[ihit][3])) &&
            time_L[ihit][4] > 0. && !(IS_NAN(time_L[ihit][4])) && time_L[ihit][5] > 0. && !(IS_NAN(time_L[ihit][5])) &&
            time_L[ihit][6] > 0. && !(IS_NAN(time_L[ihit][6])) && time_L[ihit][7] > 0. && !(IS_NAN(time_L[ihit][7])) &&

            time_T[ihit][0] > 0. && !(IS_NAN(time_T[ihit][0])) && time_T[ihit][1] > 0. && !(IS_NAN(time_T[ihit][1])) &&
            time_T[ihit][2] > 0. && !(IS_NAN(time_T[ihit][2])) && time_T[ihit][3] > 0. && !(IS_NAN(time_T[ihit][3])) &&
            time_T[ihit][4] > 0. && !(IS_NAN(time_T[ihit][4])) && time_T[ihit][5] > 0. && !(IS_NAN(time_T[ihit][5])) &&
            time_T[ihit][6] > 0. && !(IS_NAN(time_T[ihit][6])) && time_T[ihit][7] > 0. && !(IS_NAN(time_T[ihit][7])))
        {
            iLOSTypeTAMEX = true; // all 8 leading and trailing times
        }

        // We will consider only events in which booth MCFD and TAMEX see same number of channels:
        if (iLOSTypeTAMEX && iLOSTypeMCFD)
            iLOSType = true;
        //   if(iLOSTypeMCFD ) LOSType = true;

        if (iLOSType)
        {
            // calculate time over threshold and check if clock counter went out of range

            Igood_event = 1;

            int nPMT = 0;
            int nPMV = 0;
            for (int ipm = 0; ipm < 8; ipm++)
            {
                tot[ihit][ipm] = 0. / 0.;
                if (time_T[ihit][ipm] > 0. && time_L[ihit][ipm] > 0. && !(IS_NAN(time_T[ihit][ipm])) &&
                    !(IS_NAN(time_L[ihit][ipm])))
                {
                    while (time_T[ihit][ipm] - time_L[ihit][ipm] <= 0.)
                    {
                        time_T[ihit][ipm] = time_T[ihit][ipm] + 2048. * fClockFreq;
                    }

                    nPMT = nPMT + 1;
                    tot[ihit][ipm] = time_T[ihit][ipm] - time_L[ihit][ipm];
                }

                if (tot[ihit][ipm] != 0. && !(IS_NAN(tot[ihit][ipm])))
                    totsum[ihit] += tot[ihit][ipm];

                // Corrected for saturation and absorption
                if (ihit > 0)
                {
                    Double_t dthit = time_V[ihit][ipm] - time_V[ihit - 1][ipm];
                    tot_corr[ihit][ipm] = satu(ipm, tot[ihit][ipm], dthit);
                }
                else
                    tot_corr[ihit][ipm] = tot[ihit][ipm];

                totsum_corr[ihit] += tot_corr[ihit][ipm];

                if (time_L[ihit][ipm] > 0. && !(IS_NAN(time_L[ihit][ipm])))
                    timeLosT[ihit] += time_L[ihit][ipm];

                if (time_V[ihit][ipm] > 0. && !(IS_NAN(time_V[ihit][ipm])))
                {
                    timeLosM[ihit] += time_V[ihit][ipm];
                    nPMV = nPMV + 1;
                }
            }

            timeLosM[ihit] = timeLosM[ihit] / static_cast<double>(nPMV);
            timeLosT[ihit] = timeLosT[ihit] / static_cast<double>(nPMT);

            totsum[ihit] = totsum[ihit] / static_cast<double>(nPMT);
            totsum_corr[ihit] = totsum_corr[ihit] / static_cast<double>(nPMT);

            // Time resolution TAMEX
            LosTresT[ihit] = ((time_L[ihit][0] + time_L[ihit][2] + time_L[ihit][4] + time_L[ihit][6]) -
                              (time_L[ihit][1] + time_L[ihit][3] + time_L[ihit][5] + time_L[ihit][7])) /
                             4.;

            // Time resolution MCFD
            LosTresM[ihit] = ((time_V[ihit][0] + time_V[ihit][2] + time_V[ihit][4] + time_V[ihit][6]) -
                              (time_V[ihit][1] + time_V[ihit][3] + time_V[ihit][5] + time_V[ihit][7])) /
                             4.;

            // Position TAMEX:
            xT_cm[ihit] = ((time_L[ihit][5] + time_L[ihit][6]) / 2. - (time_L[ihit][1] + time_L[ihit][2]) / 2.) * (-1.);
            yT_cm[ihit] = ((time_L[ihit][7] + time_L[ihit][0]) / 2. - (time_L[ihit][3] + time_L[ihit][4]) / 2.) * (-1.);
            xT_cm[ihit] = (xT_cm[ihit] - flosOffsetXT) * flosVeffXT;
            yT_cm[ihit] = (yT_cm[ihit] - flosOffsetYT) * flosVeffYT;

            // Position MCFD:
            xV_cm[ihit] = ((time_V[ihit][5] + time_V[ihit][6]) / 2. - (time_V[ihit][1] + time_V[ihit][2]) / 2.) * (-1.);
            yV_cm[ihit] = ((time_V[ihit][7] + time_V[ihit][0]) / 2. - (time_V[ihit][3] + time_V[ihit][4]) / 2.) * (-1.);
            xV_cm[ihit] = (xV_cm[ihit] - flosOffsetX) * flosVeffX;
            yV_cm[ihit] = (yV_cm[ihit] - flosOffsetY) * flosVeffY;

            /*
                    for(int ipm=0; ipm<8; ipm++)
                    {
                              if(Icount < 500000)
                     {
                        myFile<<setprecision(10)<<ihit<<" "<<ipm<<" "<<time_V[ihit][ipm]<<" "<<time_L[ihit][ipm]<<"
               "<<tot[ihit][ipm]<<endl;
                      }
                    }
            */

            // Walk correction for MCFD and TAMEX
            for (int ipm = 0; ipm < 8; ipm++)
            {
                time_V_corr[ihit][ipm] = time_V[ihit][ipm] - walk(ipm, tot[ihit][ipm]);
                time_L_corr[ihit][ipm] = time_L[ihit][ipm] - walk(ipm + 8, tot[ihit][ipm]);
                timeLosM_corr[ihit] += time_V_corr[ihit][ipm];
                timeLosT_corr[ihit] += time_L_corr[ihit][ipm];
            }

            // Walk-corrected time-properties

            // TAMEX:
            LosTresT_corr[ihit] =
                ((time_L_corr[ihit][0] + time_L_corr[ihit][2] + time_L_corr[ihit][4] + time_L_corr[ihit][6]) -
                 (time_L_corr[ihit][1] + time_L_corr[ihit][3] + time_L_corr[ihit][5] + time_L_corr[ihit][7])) /
                4.;
            timeLosT_corr[ihit] = timeLosT_corr[ihit] / static_cast<double>(nPMT);

            // MCFD
            LosTresM_corr[ihit] = ((time_V_corr[ihit][0] + time_V[ihit][2] + time_V[ihit][4] + time_V[ihit][6]) -
                                   (time_V_corr[ihit][1] + time_V[ihit][3] + time_V[ihit][5] + time_V[ihit][7])) /
                                  4.;
            timeLosM_corr[ihit] = timeLosM_corr[ihit] / static_cast<double>(nPMV);

            // Position from ToT:
            xToT_cm[ihit] = (((tot[ihit][5] + tot[ihit][6]) / 2. - (tot[ihit][1] + tot[ihit][2]) / 2.) /
                             ((tot[ihit][1] + tot[ihit][2] + tot[ihit][5] + tot[ihit][6]) / 4.));

            yToT_cm[ihit] = (((tot[ihit][0] + tot[ihit][7]) / 2. - (tot[ihit][3] + tot[ihit][4]) / 2.) /
                             ((tot[ihit][7] + tot[ihit][0] + tot[ihit][3] + tot[ihit][4]) / 4.));

            xToT_cm[ihit] = (xToT_cm[ihit] - flosOffsetXQ) * flosVeffXQ;
            yToT_cm[ihit] = (yToT_cm[ihit] - flosOffsetYQ) * flosVeffYQ;

            x_cm[ihit] = xV_cm[ihit];
            y_cm[ihit] = yV_cm[ihit];
            Z[ihit] = totsum_corr[ihit] * fp1 + fp0;
            // Z[ihit] = totsum_corr[ihit];
            t_hit[ihit] = timeLosM_corr[ihit];

            if (fOptHisto && nPMV == 8 && nPMT == 8 && Igood_event)
            {
                // MCFD times:
                fhTres_M->Fill(LosTresM[ihit]);
                // fhTresMvsIcount->Fill(Icount,LosTresM[ihit]);
                fhTresMvsIcount->Fill(Icount, LosTresM[ihit]);
                // MCFD walk corrected times:
                fhTres_M_corr->Fill(LosTresM_corr[ihit]);
                if (ihit > 0)
                    fhTreswcMvsIcount->Fill(timeLosM[ihit] - timeLosM[ihit - 1], LosTresM[ihit]);
                // TAMEX times:
                fhTres_T->Fill(LosTresT[ihit]);
                // fhTresTvsIcount->Fill(Icount,LosTresT[ihit]);
                fhTresTvsIcount->Fill(Icount, LosTresT[ihit]);
                // TAMEX walk corrected times:
                fhTres_T_corr->Fill(LosTresT_corr[ihit]);
                if (ihit > 0)
                    fhTreswcTvsIcount->Fill(timeLosM[ihit] - timeLosM[ihit - 1], LosTresT_corr[ihit]);
                // x,y from MCFD:
                fhXY->Fill(xV_cm[ihit], yV_cm[ihit]);
                // x,y from TAMEX:
                fhXYT->Fill(xT_cm[ihit], yT_cm[ihit]);
                // x,y from ToT:
                fhXY_ToT->Fill(xToT_cm[ihit], yToT_cm[ihit]);
                // Average of three:
                fhXYproj->Fill(((xToT_cm[ihit] + xV_cm[ihit] + xT_cm[ihit]) / 3.),
                               (yToT_cm[ihit] + yV_cm[ihit] + yT_cm[ihit]) / 3.);
                // Mean value between TAMEX and MCFD:
                fhXYmean->Fill(((xV_cm[ihit] + xT_cm[ihit]) / 2.), (yV_cm[ihit] + yT_cm[ihit]) / 2.);

                // ToT
                fhQ->Fill(totsum[ihit]);
                fhQtest->Fill(totsum_corr[ihit]);
                fh_los_ihit_ToTcorr->Fill(ihit, totsum_corr[ihit]);

                for (int ipm = 0; ipm < 8; ipm++)
                {
                    //  fhQcorrvsIcount[ipm]->Fill(Icount,tot_corr[ihit][ipm]);
                    fhQvsIcount[ipm]->Fill(Icount, tot[ihit][ipm]);
                    if (ihit > 0)
                    {
                        fh_los_dt_hits_ToT_corr->Fill(time_V[ihit][ipm] - time_V[ihit - 1][ipm], tot_corr[ihit][ipm]);
                        fhQvsdt[ipm]->Fill(time_V[ihit][ipm] - time_V[ihit - 1][ipm], tot[ihit][ipm]);
                        fhQcorrvsIcount[ipm]->Fill(time_V[ihit][ipm] - time_V[ihit - 1][ipm], tot_corr[ihit][ipm]);
                    }
                }

                fhQ1_vs_Q5->Fill(tot[ihit][0], tot[ihit][4]);
                fhQ3_vs_Q7->Fill(tot[ihit][2], tot[ihit][6]);
                fhQ2_vs_Q6->Fill(tot[ihit][1], tot[ihit][5]);
                fhQ4_vs_Q8->Fill(tot[ihit][3], tot[ihit][7]);

                fhQ1_vs_Q5_corr->Fill(tot_corr[ihit][0], tot_corr[ihit][4]);
                fhQ3_vs_Q7_corr->Fill(tot_corr[ihit][2], tot_corr[ihit][6]);
                fhQ2_vs_Q6_corr->Fill(tot_corr[ihit][1], tot_corr[ihit][5]);
                fhQ4_vs_Q8_corr->Fill(tot_corr[ihit][3], tot_corr[ihit][7]);

                fhQ_T->Fill(y_cm[ihit], tot[ihit][7]);
                fhQ_L->Fill(x_cm[ihit], tot[ihit][1]);
                fhQ_B->Fill(y_cm[ihit], tot[ihit][3]);
                fhQ_R->Fill(x_cm[ihit], tot[ihit][5]);

                fhQ_RT->Fill(x_cm[ihit], tot[ihit][6]);
                fhQ_LT->Fill(y_cm[ihit], tot[ihit][0]);
                fhQ_LB->Fill(x_cm[ihit], tot[ihit][2]);
                fhQ_RB->Fill(y_cm[ihit], tot[ihit][4]);

                fhQ_T_corr->Fill(y_cm[ihit], tot_corr[ihit][7]);
                fhQ_L_corr->Fill(x_cm[ihit], tot_corr[ihit][1]);
                fhQ_B_corr->Fill(y_cm[ihit], tot_corr[ihit][3]);
                fhQ_R_corr->Fill(x_cm[ihit], tot_corr[ihit][4]);

                fhQ_RT_corr->Fill(x_cm[ihit], tot_corr[ihit][6]);
                fhQ_LT_corr->Fill(y_cm[ihit], tot_corr[ihit][0]);
                fhQ_LB_corr->Fill(y_cm[ihit], tot_corr[ihit][2]);
                fhQ_RB_corr->Fill(y_cm[ihit], tot_corr[ihit][4]);

                // Correlations:
                fhTresX_M->Fill(x_cm[ihit], LosTresM[ihit]);
                fhTresY_M->Fill(y_cm[ihit], LosTresM[ihit]);
                fhTresX_M_corr->Fill(x_cm[ihit], LosTresM_corr[ihit]);
                fhTresY_M_corr->Fill(y_cm[ihit], LosTresM_corr[ihit]);

                fhTresX_T->Fill(x_cm[ihit], LosTresT[ihit]);
                fhTresY_T->Fill(y_cm[ihit], LosTresT[ihit]);
                fhTresX_T_corr->Fill(x_cm[ihit], LosTresT_corr[ihit]);
                fhTresY_T_corr->Fill(y_cm[ihit], LosTresT_corr[ihit]);

                fhTM_vs_Q->Fill(totsum[ihit], LosTresM[ihit]);
                fhTM_vs_Q_corr->Fill(totsum[ihit], LosTresM_corr[ihit]);

                fhTT_vs_Q->Fill(totsum[ihit], LosTresT[ihit]);
                fhTT_vs_Q_corr->Fill(totsum[ihit], LosTresT_corr[ihit]);

                fhQ_vs_X->Fill(x_cm[ihit], totsum[ihit]);
                fhQ_vs_Y->Fill(y_cm[ihit], totsum[ihit]);

                fhQ_vs_X_corr->Fill(x_cm[ihit], totsum_corr[ihit]);
                fhQ_vs_Y_corr->Fill(y_cm[ihit], totsum_corr[ihit]);
            }

        } // end of LOSTYPE = true
        new ((*fHitItems)[fHitItems->GetEntriesFast()])
            R3BLosHitData(nDet, t_hit[ihit], x_cm[ihit], y_cm[ihit], Z[ihit]);

        Icount++;
    }
    //  myFile.close();
    // cout << "R3BLosCal2Hit::Exec END: " << Icount << endl;
}

void R3BLosCal2Hit::CreateHisto()
{
    Double_t fhdTxmin = 0.;
    Double_t fhdTxmax = 100.;
    Int_t fhdTxbin = 20;
    // min,max,Nbins for time spectra
    Double_t fhTmin = -5.; //-5.; //-10
    Double_t fhTmax = 5.;  // 5.;  // 10
    Int_t fhTbin = 10000;  // 10000; // 20000
    // min,max,Nbins for ToT spectra
    Double_t fhQmin = 0.;
    Double_t fhQmax = 200; // 300.; //150
    Int_t fhQbin = 2000;   // 0; //3000; //1500
    // min,max,Nbins for X and Y spectra
    Double_t fhXmin = -5.; //-5
    Double_t fhXmax = 5.;  // 5
    Int_t fhXbin = 1000;   // 1000
                           // min,max,Nbins for radius spectra
    Double_t fhRmin = -5.; // 0.;
    Double_t fhRmax = 5.;  // 10.;
    Int_t fhRbin = 1000;

    if (NULL == fhTres_M)
    {
        char strName[255];
        sprintf(strName, "LOS_dt_MCFD");
        fhTres_M = new TH1F(strName, "", fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTres_M_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_dt_MCFD_wc");
        fhTres_M_corr = new TH1F(strName, "", fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTres_T)
    {
        char strName[255];
        sprintf(strName, "LOS_dt_TAMEX");
        fhTres_T = new TH1F(strName, "", fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTres_T_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_dt_TAMEX_wc");
        fhTres_T_corr = new TH1F(strName, "", fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhQ_L)
    {
        char strName[255];
        sprintf(strName, "LOS_QL_vs_X");
        fhQ_L = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_T)
    {
        char strName[255];
        sprintf(strName, "LOS_QT_vs_Y");
        fhQ_T = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_B)
    {
        char strName[255];
        sprintf(strName, "LOS_QB_vs_Y");
        fhQ_B = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_R)
    {
        char strName[255];
        sprintf(strName, "LOS_QR_vs_X");
        fhQ_R = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_L_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QL_vs_X_corr");
        fhQ_L_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_T_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QT_vs_Y_corr");
        fhQ_T_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_B_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QB_vs_Y_corr");
        fhQ_B_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_R_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QR_vs_X_corr");
        fhQ_R_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_LB)
    {
        char strName[255];
        sprintf(strName, "LOS_QLB_vs_X");
        fhQ_LB = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_LT)
    {
        char strName[255];
        sprintf(strName, "LOS_QLT_vs_Y");
        fhQ_LT = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_RB)
    {
        char strName[255];
        sprintf(strName, "LOS_QRB_vs_Y");
        fhQ_RB = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_RT)
    {
        char strName[255];
        sprintf(strName, "LOS_QRT_vs_X");
        fhQ_RT = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_LB_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QLB_vs_X_corr");
        fhQ_LB_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_LT_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QLT_vs_Y_corr");
        fhQ_LT_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_RB_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QRB_vs_Y_corr");
        fhQ_RB_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_RT_corr)
    {
        char strName[255];
        sprintf(strName, "LOS_QRT_vs_X_corr");
        fhQ_RT_corr = new TH2F(strName, "", fhRbin, fhRmin, fhRmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ)
    {
        char strName[255];
        sprintf(strName, "LOS_Q");
        fhQ = new TH1F(strName, "", 10000, fhQmin, fhQmax);
    }

    if (NULL == fhQtest)
    {
        char strName[255];
        sprintf(strName, "LOS_Q_test");
        fhQtest = new TH1F(strName, "", 10000, fhQmin, fhQmax);
    }

    if (NULL == fhTM_vs_Q)
    {
        char strName[255];
        sprintf(strName, "TMCFD_vs_Q");
        fhTM_vs_Q = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTT_vs_Q)
    {
        char strName[255];
        sprintf(strName, "TTAMEX_vs_Q");
        fhTT_vs_Q = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTM_vs_Q_corr)
    {
        char strName[255];
        sprintf(strName, "TMCFD_vs_Q_wc");
        fhTM_vs_Q_corr = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTT_vs_Q_corr)
    {
        char strName[255];
        sprintf(strName, "TTAMEX_vs_Q_wc");
        fhTT_vs_Q_corr = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresX_M)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_X_MCFD");
        fhTresX_M = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresY_M)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_Y_MCFD");
        fhTresY_M = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresX_M_corr)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_X_MCFD_wc");
        fhTresX_M_corr = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresY_M_corr)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_Y_MCFD_wc");
        fhTresY_M_corr = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresX_T)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_X_TAMEX");
        fhTresX_T = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresY_T)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_Y_TAMEX");
        fhTresY_T = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresX_T_corr)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_X_TAMEX_wc");
        fhTresX_T_corr = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhTresY_T_corr)
    {
        char strName[255];
        sprintf(strName, "Tres_vs_Y_TAMEX_wc");
        fhTresY_T_corr = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhTbin, fhTmin, fhTmax);
    }

    if (NULL == fhXY)
    {
        char strName[255];
        sprintf(strName, "X_vs_Y");
        fhXY = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhXbin, fhXmin, fhXmax);
    }

    if (NULL == fhXYT)
    {
        char strName[255];
        sprintf(strName, "X_vs_Y_TAMEX");
        fhXYT = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhXbin, fhXmin, fhXmax);
    }

    if (NULL == fhXYmean)
    {
        char strName[255];
        sprintf(strName, "X_vs_Y_mean");
        fhXYmean = new TH2F(strName, "", 4000, fhXmin, fhXmax, 4000, fhXmin, fhXmax);
    }

    if (NULL == fhXY_ToT)
    {
        char strName[255];
        sprintf(strName, "X_vs_Y_ToT");
        fhXY_ToT = new TH2F(strName, "", 4000, fhXmin, fhXmax, 4000, fhXmin, fhXmax);
        fhXY_ToT->GetXaxis()->SetTitle("X position / cm");
        fhXY_ToT->GetYaxis()->SetTitle("Y position / cm");
    }

    if (NULL == fhXYproj)
    {
        char strName[255];
        sprintf(strName, "X2_vs_Y2_proj");
        fhXYproj = new TH2F(strName, "", 4000, fhXmin, fhXmax, 4000, fhXmin, fhXmax);
    }

    if (NULL == fhQ_vs_X)
    {
        char strName[255];
        sprintf(strName, "Q_vs_X");
        fhQ_vs_X = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_vs_Y)
    {
        char strName[255];
        sprintf(strName, "Q_vs_Y");
        fhQ_vs_Y = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_vs_X_corr)
    {
        char strName[255];
        sprintf(strName, "Qcorr_vs_X");
        fhQ_vs_X_corr = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ_vs_Y_corr)
    {
        char strName[255];
        sprintf(strName, "Qcorr_vs_Y");
        fhQ_vs_Y_corr = new TH2F(strName, "", fhXbin, fhXmin, fhXmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ1_vs_Q5)
    {
        char strName[255];
        sprintf(strName, "QPM1_vs_QPM5");
        fhQ1_vs_Q5 = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ1_vs_Q5_corr)
    {
        char strName[255];
        sprintf(strName, "QPM1_vs_QPM5_corr");
        fhQ1_vs_Q5_corr = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ2_vs_Q6)
    {
        char strName[255];
        sprintf(strName, "QPM2_vs_QPM6");
        fhQ2_vs_Q6 = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ2_vs_Q6_corr)
    {
        char strName[255];
        sprintf(strName, "QPM2_vs_QPM6_corr");
        fhQ2_vs_Q6_corr = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ3_vs_Q7)
    {
        char strName[255];
        sprintf(strName, "QPM3_vs_QPM7");
        fhQ3_vs_Q7 = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ3_vs_Q7_corr)
    {
        char strName[255];
        sprintf(strName, "QPM3_vs_QPM7_corr");
        fhQ3_vs_Q7_corr = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ4_vs_Q8)
    {
        char strName[255];
        sprintf(strName, "QPM4_vs_QPM8");
        fhQ4_vs_Q8 = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    if (NULL == fhQ4_vs_Q8_corr)
    {
        char strName[255];
        sprintf(strName, "QPM4_vs_QPM8_corr");
        fhQ4_vs_Q8_corr = new TH2F(strName, "", fhQbin, fhQmin, fhQmax, fhQbin, fhQmin, fhQmax);
    }

    for (Int_t j = 0; j < 8; j++)
    {
        if (NULL == fhQvsdt[j])
        {
            char strName[255];
            sprintf(strName, "Q_vs_dt_ch%d", j + 1);
            fhQvsdt[j] = new TH2F(strName, "", 6000, 0, 3000, fhQbin, fhQmin, fhQmax);
            fhQvsdt[j]->GetXaxis()->SetTitle("time between two hits / ns");
            char strYname[255];
            sprintf(strYname, "ToT / ns of ch. %d", j + 1);
            fhQvsdt[j]->GetYaxis()->SetTitle(strYname);
        }
        if (NULL == fhQcorrvsIcount[j])
        {
            char strName[255];
            sprintf(strName, "Qcorr_vs_dt_ch%d", j + 1);
            fhQcorrvsIcount[j] = new TH2F(strName, "", 6000, 0, 3000, fhQbin, fhQmin, fhQmax);
            fhQcorrvsIcount[j]->GetXaxis()->SetTitle("Event number");
            char strYname[255];
            sprintf(strYname, "ToTcorr / ns of ch. %d", j + 1);
            fhQcorrvsIcount[j]->GetYaxis()->SetTitle(strYname);
        }

        if (NULL == fhQvsIcount[j])
        {
            char strName[255];
            sprintf(strName, "Q_vs_Icount_ch%d", j + 1);
            fhQvsIcount[j] = new TH2F(strName, "", 10000, 0, 10000000, fhQbin, fhQmin, fhQmax);
            fhQvsIcount[j]->GetXaxis()->SetTitle("Event number");
            char strYname[255];
            sprintf(strYname, "ToT / ns of ch. %d", j + 1);
            fhQvsIcount[j]->GetYaxis()->SetTitle(strYname);
        }
    }

    if (NULL == fh_los_ihit_ToTcorr)
    {
        char strName[255];
        sprintf(strName, "los_totcorr_ihit");
        fh_los_ihit_ToTcorr = new TH2F(strName, "", 10, 0, 10, 3100, -10., 300.);
        fh_los_ihit_ToTcorr->GetXaxis()->SetTitle("iHit");
        fh_los_ihit_ToTcorr->GetYaxis()->SetTitle("ToT / ns");
    }

    if (NULL == fhTresMvsIcount)
    {
        char strName[255];
        sprintf(strName, "TresM_vs_Icount");
        fhTresMvsIcount = new TH2F(strName, "", 10000, 0, 10000000, fhTbin, fhTmin, fhTmax);
        fhTresMvsIcount->GetXaxis()->SetTitle("Event number");
        fhTresMvsIcount->GetYaxis()->SetTitle("MCFD time precision / ns");
    }

    if (NULL == fhTreswcMvsIcount)
    {
        char strName[255];
        sprintf(strName, "TreswcM_vs_dt");
        fhTreswcMvsIcount = new TH2F(strName, "", 6000, 0, 3000, fhTbin, fhTmin, fhTmax);
        fhTreswcMvsIcount->GetXaxis()->SetTitle("time between two hits / ns");
        fhTreswcMvsIcount->GetYaxis()->SetTitle("MCFD time precision / ns");
    }

    if (NULL == fhTresTvsIcount)
    {
        char strName[255];
        sprintf(strName, "TresT_vs_Icount");
        fhTresTvsIcount = new TH2F(strName, "", 10000, 0, 10000000, fhTbin, fhTmin, fhTmax);
        fhTresTvsIcount->GetXaxis()->SetTitle("Event number");
        fhTresTvsIcount->GetYaxis()->SetTitle("MCFD time precision / ns");
    }

    if (NULL == fhTreswcTvsIcount)
    {
        char strName[255];
        sprintf(strName, "TreswcT_vs_dt");
        fhTreswcTvsIcount = new TH2F(strName, "", 6000, 0, 3000, fhTbin, fhTmin, fhTmax);
        fhTreswcTvsIcount->GetXaxis()->SetTitle("time between two hits / ns");
        fhTreswcTvsIcount->GetYaxis()->SetTitle("MCFD time precision / ns");
    }

    if (NULL == fh_los_dt_hits_ToT_corr)
    {
        char strName[255];
        sprintf(strName, "los_dt_hits_ToT_corr");
        fh_los_dt_hits_ToT_corr = new TH2F(strName, "", 6000, 0, 3000, fhQbin, fhQmin, fhQmax);
    }
}

void R3BLosCal2Hit::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
}

void R3BLosCal2Hit::FinishTask()
{
    if (fOptHisto)
    {
        if (fhTres_M)
            fhTres_M->Write();
        if (fhTres_T)
            fhTres_T->Write();
        if (fhQ)
            fhQ->Write();
        if (fhQtest)
            fhQtest->Write();
        if (fhQ_vs_X)
            fhQ_vs_X->Write();
        if (fhQ_vs_Y)
            fhQ_vs_Y->Write();
        if (fhTM_vs_Q)
            fhTM_vs_Q->Write();
        if (fhTT_vs_Q)
            fhTT_vs_Q->Write();
        if (fhXY)
            fhXY->Write();
        if (fhXYT)
            fhXYT->Write();
        if (fhXYmean)
            fhXYmean->Write();
        if (fhXY_ToT)
            fhXY_ToT->Write();
        if (fhXYproj)
            fhXYproj->Write();
        if (fhQ1_vs_Q5)
            fhQ1_vs_Q5->Write();
        if (fhQ2_vs_Q6)
            fhQ2_vs_Q6->Write();
        if (fhQ3_vs_Q7)
            fhQ3_vs_Q7->Write();
        if (fhQ4_vs_Q8)
            fhQ4_vs_Q8->Write();
        if (fhTresX_M)
            fhTresX_M->Write();
        if (fhTresY_M)
            fhTresY_M->Write();
        if (fhTresX_T)
            fhTresX_T->Write();
        if (fhTresY_T)
            fhTresY_T->Write();
        if (fhTres_T_corr)
            fhTres_T_corr->Write();
        if (fhTres_M_corr)
            fhTres_M_corr->Write();
        if (fhQ_L)
            fhQ_L->Write();
        if (fhQ_B)
            fhQ_B->Write();
        if (fhQ_R)
            fhQ_R->Write();
        if (fhQ_T)
            fhQ_T->Write();
        if (fhQ_L_corr)
            fhQ_L_corr->Write();
        if (fhQ_B_corr)
            fhQ_B_corr->Write();
        if (fhQ_R_corr)
            fhQ_R_corr->Write();
        if (fhQ_T_corr)
            fhQ_T_corr->Write();
        if (fhQ_LT)
            fhQ_LT->Write();
        if (fhQ_LB)
            fhQ_LB->Write();
        if (fhQ_RB)
            fhQ_RB->Write();
        if (fhQ_RT)
            fhQ_RT->Write();
        if (fhQ_LT_corr)
            fhQ_LT_corr->Write();
        if (fhQ_LB_corr)
            fhQ_LB_corr->Write();
        if (fhQ_RB_corr)
            fhQ_RB_corr->Write();
        if (fhQ_RT_corr)
            fhQ_RT_corr->Write();
        if (fhQ_vs_X_corr)
            fhQ_vs_X_corr->Write();
        if (fhQ_vs_Y_corr)
            fhQ_vs_Y_corr->Write();
        if (fhTM_vs_Q_corr)
            fhTM_vs_Q_corr->Write();
        if (fhTT_vs_Q_corr)
            fhTT_vs_Q_corr->Write();
        if (fhQ1_vs_Q5_corr)
            fhQ1_vs_Q5_corr->Write();
        if (fhQ2_vs_Q6_corr)
            fhQ2_vs_Q6_corr->Write();
        if (fhQ3_vs_Q7_corr)
            fhQ3_vs_Q7_corr->Write();
        if (fhQ4_vs_Q8_corr)
            fhQ4_vs_Q8_corr->Write();

        if (fhTresX_M_corr)
            fhTresX_M_corr->Write();
        if (fhTresY_M_corr)
            fhTresY_M_corr->Write();
        if (fhTresX_T_corr)
            fhTresX_T_corr->Write();
        if (fhTresY_T_corr)
            fhTresY_T_corr->Write();

        if (fhTresMvsIcount)
            fhTresMvsIcount->Write();
        if (fhTresTvsIcount)
            fhTresTvsIcount->Write();
        if (fhTreswcMvsIcount)
            fhTreswcMvsIcount->Write();
        if (fhTreswcTvsIcount)
            fhTreswcTvsIcount->Write();

        if (fh_los_ihit_ToTcorr)
            fh_los_ihit_ToTcorr->Write();
        if (fh_los_dt_hits_ToT_corr)
            fh_los_dt_hits_ToT_corr->Write();

        for (Int_t j = 0; j < 8; j++)
        {
            if (fhQvsdt[j])
                fhQvsdt[j]->Write();
            if (fhQcorrvsIcount[j])
                fhQcorrvsIcount[j]->Write();
            if (fhQvsIcount[j])
                fhQvsIcount[j]->Write();
        }
    }
}

Double_t R3BLosCal2Hit::walk(Int_t inum, Double_t tot)
{

    Double_t y = 0. / 0., ysc = 0. / 0., term[8] = { 0. };
    Double_t x;

    x = tot;
    term[0] = x;
    for (Int_t i = 0; i < 7; i++)
    {
        term[i + 1] = term[i] * x;
    }

    ysc = walk_par[inum][2] + walk_par[inum][3] * term[0] + walk_par[inum][4] * term[1] + walk_par[inum][5] * term[2] +
          walk_par[inum][6] * term[3] + walk_par[inum][7] * term[4] + walk_par[inum][8] * term[5] +
          walk_par[inum][9] * term[6] + walk_par[inum][10] * term[7];

    if (tot < walk_par[inum][0] || tot > walk_par[inum][1])
        ysc = 0.0 / 0.0;

    /*
       cout<< inum<<", "<< walk_par[inum][0] <<", "<< walk_par[inum][1] <<", "<< walk_par[inum][2] <<", "<<
                           walk_par[inum][3] <<", "<< walk_par[inum][4] <<", "<< walk_par[inum][5] <<", "<<
                           walk_par[inum][6] <<", "<< walk_par[inum][7] <<", "<< walk_par[inum][8] <<", "<<
                           walk_par[inum][9] <<", "<< walk_par[inum][10] <<endl;
    */
    return ysc;
}

Double_t R3BLosCal2Hit::satu(Int_t inum, Double_t tot, Double_t dt)
{

    Double_t ysc = 0. / 0.;

    // if(tot_par[inum][0] > 0.)
    // ysc  = (tot_par[inum][0]*tot+tot_par[inum][1])/(tot_par[inum][2]-tot)*tot_par[inum][3] ;

    ysc = tot_par[inum][0] + tot_par[inum][1] * (1. - 1. / (exp((dt - tot_par[inum][2]) / tot_par[inum][3]) + 1.));
    ysc = tot / ysc;
    ysc = ysc * (tot_par[inum][0] + tot_par[inum][1]);

    /*
     cout<< inum<<", "<< tot_par[inum][0] <<", "<< tot_par[inum][1] <<", "<< tot_par[inum][2]<<", "<< tot_par[inum][3]
     <<endl; cout<<tot<<", "<<ysc<<endl;
  */
    return ysc;
}

ClassImp(R3BLosCal2Hit);
