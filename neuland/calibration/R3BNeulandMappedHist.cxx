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

#include "R3BNeulandMappedHist.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "R3BEventHeader.h"
#include "R3BLosMappedData.h"
#include "R3BNeulandTacquilaMappedData.h"
#include "R3BPaddleTamexMappedData.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"

R3BNeulandMappedHist::R3BNeulandMappedHist()
    : FairTask("NeulandMappedHist", 1)
    , fnEvents(0)
    , fNItemsTotal(0)
    , fHeader(NULL)
    , fLandMappedData(NULL)
    , fNeulandTamexHitMapped(NULL)
    , fLosMappedData(NULL)
{
}

R3BNeulandMappedHist::R3BNeulandMappedHist(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fnEvents(0)
    , fNItemsTotal(0)
    , fHeader(NULL)
    , fLandMappedData(NULL)
    , fNeulandTamexHitMapped(NULL)
    , fLosMappedData(NULL)
{
}

R3BNeulandMappedHist::~R3BNeulandMappedHist() {}

InitStatus R3BNeulandMappedHist::Init()
{
    FairRootManager* fMan = FairRootManager::Instance();
    fHeader = dynamic_cast<R3BEventHeader*>(fMan->GetObject("EventHeader."));
    fLandMappedData = dynamic_cast<TClonesArray*>(fMan->GetObject("NeulandTacquilaMappedData"));
    fLosMappedData = dynamic_cast<TClonesArray*>(fMan->GetObject("LosMapped"));
    fNeulandTamexHitMapped = dynamic_cast<TClonesArray*>(fMan->GetObject("NeulandMappedData"));
    CreateHistos();

    return kSUCCESS;
}

void R3BNeulandMappedHist::Exec(Option_t* option)
{
    if (fHeader)
    {
        fh_trigger->Fill(fHeader->GetTrigger());
    }

    if (fLandMappedData)
    {
        Int_t nLandMapped = fLandMappedData->GetEntries();
        fNItemsTotal += nLandMapped;
        R3BNeulandTacquilaMappedData* hitmapped;
        for (Int_t i = 0; i < nLandMapped; i++)
        {
            hitmapped = dynamic_cast<R3BNeulandTacquilaMappedData*>(fLandMappedData->At(i));
            fh_land_mapped_barid->Fill((hitmapped->GetPlane() - 1) * 50 + hitmapped->GetPaddle());
            fh_land_mapped_side->Fill(hitmapped->GetSide());
            fh_land_mapped_clock->Fill(hitmapped->GetClock());
            fh_land_mapped_tac->Fill(hitmapped->GetTacData());
            fh_land_mapped_qdc->Fill(hitmapped->GetQdcData());
        }
    }

    if (fNeulandTamexHitMapped)
    {
        Int_t nNeulandTamexHitsMapped = fNeulandTamexHitMapped->GetEntries();
        //        fNItemsTotal += nLandRawHitsMapped;
        R3BPaddleTamexMappedData* hitmapped;
        for (Int_t i = 0; i < nNeulandTamexHitsMapped; i++)
        {
            hitmapped = dynamic_cast<R3BPaddleTamexMappedData*>(fNeulandTamexHitMapped->At(i));
            fh_neuland_mapped_is17->Fill(hitmapped->Is17());
            if (!hitmapped->Is17())
            {
                fh_neuland_mapped_planeid->Fill(hitmapped->GetPlaneId());
                fh_neuland_mapped_barid->Fill(hitmapped->GetBarId());
                // fh_neuland_mapped_side->Fill(hitmapped->GetSide());

                fh_neuland_mapped_cle->Fill(hitmapped->GetCoarseTime(0, 0)); // PM1 LE
                fh_neuland_mapped_cte->Fill(hitmapped->GetCoarseTime(0, 1)); // PM1 TE
                fh_neuland_mapped_fle->Fill(hitmapped->GetFineTime(0, 0));   // PM1 LE
                fh_neuland_mapped_fte->Fill(hitmapped->GetFineTime(0, 1));   // PM1 TE

                fh_neuland_mapped_cle->Fill(hitmapped->GetCoarseTime(1, 0)); // PM2 LE
                fh_neuland_mapped_cte->Fill(hitmapped->GetCoarseTime(1, 1)); // PM2 TE
                fh_neuland_mapped_fle->Fill(hitmapped->GetFineTime(1, 0));   // PM2 LE
                fh_neuland_mapped_fte->Fill(hitmapped->GetFineTime(1, 1));   // PM2 TE
                //        cout<<"in ana: "<<hitmapped->GetPlaneId()<<"  "<<hitmapped->GetBarId()<<endl;
            }
        }
    }

    if (fLosMappedData)
    {
        Int_t nLosMapped = fLosMappedData->GetEntries();
        R3BLosMappedData* loshit;
        for (Int_t i = 0; i < nLosMapped; i++)
        {
            loshit = dynamic_cast<R3BLosMappedData*>(fLosMappedData->At(i));
            fh_los_det->Fill(loshit->GetDetector());
            fh_los_ch->Fill(loshit->GetChannel());
            fh_los_tcoarse->Fill(loshit->GetTimeCoarse());
            fh_los_tfine->Fill(loshit->GetTimeFine());
        }
    }

    fnEvents += 1;
    // if (0 == (fnEvents % 1000))
    //{
    //    LOG(info) << "R3BNeulandMappedHist : " << fnEvents << " events collected";
    //}
}

void R3BNeulandMappedHist::FinishTask() { WriteHistos(); }

void R3BNeulandMappedHist::CreateHistos()
{
    fh_trigger = new TH1F("h_trigger", "Trigger", 10, -0.5, 9.5);

    fh_land_mapped_barid = new TH1F("h_land_mapped_barid", "Bar ID", 500, -0.5, 499.5);
    fh_land_mapped_side = new TH1F("h_land_mapped_side", "Side", 5, -0.5, 4.5);
    fh_land_mapped_clock = new TH1F("h_land_mapped_clock", "Clock count", 70, -0.5, 69.5);
    fh_land_mapped_tac = new TH1F("h_land_mapped_tac", "TAC data", 500, 0., 5000.);
    fh_land_mapped_qdc = new TH1F("h_land_mapped_qdc", "QDC data", 100, 0., 1000.);

    fh_neuland_mapped_is17 = new TH1F("h_neuland_mapped_is17", "Is 17", 4, -0.5, 3.5);
    fh_neuland_mapped_planeid = new TH1F("h_neuland_mapped_planeid", "Plane ID", 30, -0.5, 29.5);
    fh_neuland_mapped_barid = new TH1F("h_neuland_mapped_barid", "Bar ID", 50, -0.5, 49.5);
    fh_neuland_mapped_side = new TH1F("h_neuland_mapped_side", "Side", 5, -0.5, 4.5);
    fh_neuland_mapped_cle = new TH1F("h_neuland_mapped_cle", "Coarse time LE", 8200, 0., 8200.);
    fh_neuland_mapped_cte = new TH1F("h_neuland_mapped_cte", "Coarse time TE", 8200, 0., 8200.);
    fh_neuland_mapped_fle = new TH1F("h_neuland_mapped_fle", "Fine time LE", 600, 0., 600.);
    fh_neuland_mapped_fte = new TH1F("h_neuland_mapped_fte", "Fine time TE", 600, 0., 600.);

    fh_los_det = new TH1F("h_los_det", "Detector", 20, -0.5, 19.5);
    fh_los_ch = new TH1F("h_los_ch", "Channel", 20, -0.5, 19.5);
    fh_los_tcoarse = new TH1F("h_los_tcoarse", "Time coarse", 1000, 0., 10000.);
    fh_los_tfine = new TH1F("h_los_tfine", "Time fine", 200, 0., 2000.);

    FairRunOnline* run = FairRunOnline::Instance();

    run->AddObject(fh_trigger);

    run->AddObject(fh_land_mapped_barid);
    run->AddObject(fh_land_mapped_side);
    run->AddObject(fh_land_mapped_clock);
    run->AddObject(fh_land_mapped_tac);
    run->AddObject(fh_land_mapped_qdc);

    run->AddObject(fh_neuland_mapped_is17);
    run->AddObject(fh_neuland_mapped_planeid);
    run->AddObject(fh_neuland_mapped_barid);
    run->AddObject(fh_neuland_mapped_side);
    run->AddObject(fh_neuland_mapped_cle);
    run->AddObject(fh_neuland_mapped_cte);
    run->AddObject(fh_neuland_mapped_fle);
    run->AddObject(fh_neuland_mapped_fte);

    run->AddObject(fh_los_det);
    run->AddObject(fh_los_ch);
    run->AddObject(fh_los_tcoarse);
    run->AddObject(fh_los_tfine);

    run->RegisterHttpCommand("/Reset_h_land_mapped_barid", "/h_land_mapped_barid/->Reset()");
}

void R3BNeulandMappedHist::WriteHistos()
{
    fh_neuland_mapped_is17->Write();
    fh_neuland_mapped_planeid->Write();
    fh_neuland_mapped_barid->Write();
    fh_neuland_mapped_side->Write();
    fh_neuland_mapped_cle->Write();
    fh_neuland_mapped_cte->Write();
    fh_neuland_mapped_fle->Write();
    fh_neuland_mapped_fte->Write();
}

ClassImp(R3BNeulandMappedHist)
