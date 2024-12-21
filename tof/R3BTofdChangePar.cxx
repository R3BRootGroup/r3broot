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
// --------------------------------------------------
// -----            R3BTofdChangePar            -----
// -----       Created Sep 2019 by L.Bott       -----
// --------------------------------------------------

/* Some notes:
 *
 *
 */

#include "R3BTofdChangePar.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"
#include "R3BTofdCalData.h"
#include "R3BTofdHitPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <stdlib.h>

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
} // namespace

R3BTofdChangePar::R3BTofdChangePar()
    : FairTask("R3BTofdChangePar", 1)
    , fCalItemsLos(NULL)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fParaFile("")
{
}

R3BTofdChangePar::R3BTofdChangePar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItemsLos(NULL)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fParaFile("")
{
}

R3BTofdChangePar::~R3BTofdChangePar()
{
    if (fCal_Par)
    {
        delete fCal_Par;
    }
}

InitStatus R3BTofdChangePar::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(rm->GetObject("R3BEventHeader"));
    // may be = NULL!

    fCalData = dynamic_cast<TClonesArray*>(rm->GetObject("TofdCal"));
    if (!fCalData)
    {
        return kFATAL;
    }

    if (!fNofModules)
    {
        LOG(error) << "R3BTofdChangePar::Init() Number of modules not set. ";
        return kFATAL;
    }

    return kSUCCESS;
}

void R3BTofdChangePar::SetParContainers()
{
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCal_Par = dynamic_cast<R3BTofdHitPar*>(FairRuntimeDb::instance()->getContainer("TofdHitPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BTofdChangePar::Init() Couldn't get handle on TofdHitPar. ";
    }
    //	    fCal_Par->setChanged();
}

void R3BTofdChangePar::Exec(Option_t* option) {}

void R3BTofdChangePar::FinishEvent() {}

void R3BTofdChangePar::FinishTask()
{
    // set new parameters
    Double_t pars[5];
    Double_t parameter;
    Int_t plane, pplane, bar, pbar, pm, ppm, p, pp, i;
    if (fParaFile != "")
    {
        ifstream a_file(fParaFile);
        if (!a_file)
        {
            cout << "Could not open file!!!" << endl;
            exit(1);
        }
        i = 0;
        while (a_file.good())
        {
            for (Int_t a = 0; a < 5; a++)
            {
                pars[a] = 0.;
            }
            pplane = plane;
            pbar = bar;
            ppm = pm;
            pp = p;
            a_file >> plane >> bar >> pm >> p >> parameter;
            if (plane != pplane || bar != pbar || pm != ppm || p != pp)
                i = 0;
            pars[i] = parameter;
            if (a_file.eof())
                break;
            cout << "Reading Parameter, plane " << plane << " bar " << bar << " PM " << pm << " #" << p << " : "
                 << "i " << i << " = " << pars[i] << endl;
            i++;

            if (p == 0)
                changeAll0(plane, bar, pm, pars);
            if (p == 1)
                changeOffset(plane, bar, pm, pars);
            if (p == 2)
                changeToTOffset(plane, bar, pm, pars);
            if (p == 3)
                changeSync(plane, bar, pm, pars);
            if (p == 4)
                changeVeff(plane, bar, pm, pars);
            if (p == 5)
                changeLambda(plane, bar, pm, pars);
            if (p == 6)
                changeDoubleExp(plane, bar, pm, pars);
            if (p == 7)
                changeZ(plane, bar, pm, pars);
            if (p == 8)
                changeWalk(plane, bar, pm, pars);
        }
        a_file.close();
    }
    fCal_Par->setChanged();
}

void R3BTofdChangePar::changeAll0(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    par->SetOffset1(0.);
    par->SetOffset2(0.);
    par->SetToTOffset1(0.);
    par->SetToTOffset2(0.);
    par->SetSync(0.);
    par->SetVeff(0.);
    par->SetLambda(0.);
    par->SetPar1a(0.);
    par->SetPar1b(0.);
    par->SetPar1c(0.);
    par->SetPar1d(0.);
    par->SetPar1za(0.);
    par->SetPar1zb(0.);
    par->SetPar1zc(0.);
    par->SetPar2a(0.);
    par->SetPar2b(0.);
    par->SetPar2c(0.);
    par->SetPar2d(0.);
    par->SetPar1za(0.);
    par->SetPar1zb(0.);
    par->SetPar1zc(0.);
    par->SetPar1Walk(0.);
    par->SetPar2Walk(0.);
    par->SetPar3Walk(0.);
    par->SetPar4Walk(0.);
    par->SetPar5Walk(0.);
}

void R3BTofdChangePar::changeOffset(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[4];
    Double_t para2[3];
    Double_t para3[5];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    if (pm == 1)
    {
        if (pars[0] == 0.)
            pars[0] = par->GetOffset1();
        par->SetOffset1(pars[0]);
    }
    if (pm == 2)
    {
        if (pars[0] == 0.)
            pars[0] = par->GetOffset2();
        par->SetOffset2(pars[0]);
    }
    Double_t ToTOffset1 = par->GetToTOffset1();
    Double_t ToTOffset2 = par->GetToTOffset2();
    Double_t sync = par->GetSync();
    Double_t veff = par->GetVeff();
    Double_t lambda = par->GetLambda();
    if (pm == 1)
    {
        para1[0] = par->GetPar1a();
        para1[1] = par->GetPar1b();
        para1[2] = par->GetPar1c();
        para1[3] = par->GetPar1d();
    }
    else
    {
        para1[0] = par->GetPar2a();
        para1[1] = par->GetPar2b();
        para1[2] = par->GetPar2c();
        para1[3] = par->GetPar2d();
    }
    para2[0] = par->GetPar1za();
    para2[1] = par->GetPar1zb();
    para2[2] = par->GetPar1zc();
    para3[0] = par->GetPar1Walk();
    para3[1] = par->GetPar2Walk();
    para3[2] = par->GetPar3Walk();
    para3[3] = par->GetPar4Walk();
    para3[4] = par->GetPar5Walk();
}

void R3BTofdChangePar::changeToTOffset(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[4];
    Double_t para2[3];
    Double_t para3[5];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    Double_t offset1 = par->GetOffset1();
    Double_t offset2 = par->GetOffset2();
    if (pm == 1)
    {
        if (pars[0] == 0.)
            pars[0] = par->GetToTOffset1();
        par->SetToTOffset1(pars[0]);
    }
    if (pm == 2)
    {
        if (pars[0] == 0.)
            pars[0] = par->GetToTOffset2();
        par->SetToTOffset2(pars[0]);
    }
    Double_t sync = par->GetSync();
    Double_t veff = par->GetVeff();
    Double_t lambda = par->GetLambda();
    if (pm == 1)
    {
        para1[0] = par->GetPar1a();
        para1[1] = par->GetPar1b();
        para1[2] = par->GetPar1c();
        para1[3] = par->GetPar1d();
    }
    else
    {
        para1[0] = par->GetPar2a();
        para1[1] = par->GetPar2b();
        para1[2] = par->GetPar2c();
        para1[3] = par->GetPar2d();
    }
    para2[0] = par->GetPar1za();
    para2[1] = par->GetPar1zb();
    para2[2] = par->GetPar1zc();
    para3[0] = par->GetPar1Walk();
    para3[1] = par->GetPar2Walk();
    para3[2] = par->GetPar3Walk();
    para3[3] = par->GetPar4Walk();
    para3[4] = par->GetPar5Walk();
}

void R3BTofdChangePar::changeSync(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[4];
    Double_t para2[3];
    Double_t para3[5];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    Double_t offset1 = par->GetOffset1();
    Double_t offset2 = par->GetOffset2();
    Double_t ToTOffset1 = par->GetToTOffset1();
    Double_t ToTOffset2 = par->GetToTOffset2();
    if (pars[0] == 0.)
        pars[0] = par->GetSync();
    par->SetSync(pars[0]);
    Double_t veff = par->GetVeff();
    Double_t lambda = par->GetLambda();
    if (pm == 1)
    {
        para1[0] = par->GetPar1a();
        para1[1] = par->GetPar1b();
        para1[2] = par->GetPar1c();
        para1[3] = par->GetPar1d();
    }
    else
    {
        para1[0] = par->GetPar2a();
        para1[1] = par->GetPar2b();
        para1[2] = par->GetPar2c();
        para1[3] = par->GetPar2d();
    }
    para2[0] = par->GetPar1za();
    para2[1] = par->GetPar1zb();
    para2[2] = par->GetPar1zc();
    para3[0] = par->GetPar1Walk();
    para3[1] = par->GetPar2Walk();
    para3[2] = par->GetPar3Walk();
    para3[3] = par->GetPar4Walk();
    para3[4] = par->GetPar5Walk();
}

void R3BTofdChangePar::changeVeff(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[4];
    Double_t para2[3];
    Double_t para3[5];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    Double_t offset1 = par->GetOffset1();
    Double_t offset2 = par->GetOffset2();
    Double_t ToTOffset1 = par->GetToTOffset1();
    Double_t ToTOffset2 = par->GetToTOffset2();
    Double_t sync = par->GetSync();
    if (pars[0] == 0.)
        pars[0] = par->GetVeff();
    par->SetVeff(pars[0]);
    Double_t lambda = par->GetLambda();
    if (pm == 1)
    {
        para1[0] = par->GetPar1a();
        para1[1] = par->GetPar1b();
        para1[2] = par->GetPar1c();
        para1[3] = par->GetPar1d();
    }
    else
    {
        para1[0] = par->GetPar2a();
        para1[1] = par->GetPar2b();
        para1[2] = par->GetPar2c();
        para1[3] = par->GetPar2d();
    }
    para2[0] = par->GetPar1za();
    para2[1] = par->GetPar1zb();
    para2[2] = par->GetPar1zc();
    para3[0] = par->GetPar1Walk();
    para3[1] = par->GetPar2Walk();
    para3[2] = par->GetPar3Walk();
    para3[3] = par->GetPar4Walk();
    para3[4] = par->GetPar5Walk();
}

void R3BTofdChangePar::changeLambda(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[4];
    Double_t para2[3];
    Double_t para3[5];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    Double_t offset1 = par->GetOffset1();
    Double_t offset2 = par->GetOffset2();
    Double_t ToTOffset1 = par->GetToTOffset1();
    Double_t ToTOffset2 = par->GetToTOffset2();
    Double_t sync = par->GetSync();
    Double_t veff = par->GetVeff();
    if (pars[0] == 0.)
        pars[0] = par->GetLambda();
    par->SetLambda(pars[0]);
    if (pm == 1)
    {
        para1[0] = par->GetPar1a();
        para1[1] = par->GetPar1b();
        para1[2] = par->GetPar1c();
        para1[3] = par->GetPar1d();
    }
    else
    {
        para1[0] = par->GetPar2a();
        para1[1] = par->GetPar2b();
        para1[2] = par->GetPar2c();
        para1[3] = par->GetPar2d();
    }
    para2[0] = par->GetPar1za();
    para2[1] = par->GetPar1zb();
    para2[2] = par->GetPar1zc();
    para3[0] = par->GetPar1Walk();
    para3[1] = par->GetPar2Walk();
    para3[2] = par->GetPar3Walk();
    para3[3] = par->GetPar4Walk();
    para3[4] = par->GetPar5Walk();
}

void R3BTofdChangePar::changeDoubleExp(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[3];
    Double_t para2[5];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    Double_t offset1 = par->GetOffset1();
    Double_t offset2 = par->GetOffset2();
    Double_t ToTOffset1 = par->GetToTOffset1();
    Double_t ToTOffset2 = par->GetToTOffset2();
    Double_t sync = par->GetSync();
    Double_t veff = par->GetVeff();
    Double_t lambda = par->GetLambda();
    if (pars[0] == 0.)
        pars[0] = par->GetPar1a();
    if (pars[1] == 0.)
        pars[1] = par->GetPar1b();
    if (pars[2] == 0.)
        pars[2] = par->GetPar1c();
    if (pars[3] == 0.)
        pars[3] = par->GetPar1d();
    if (pm == 1)
    {
        par->SetPar1a(pars[0]);
        par->SetPar1b(pars[1]);
        par->SetPar1c(pars[2]);
        par->SetPar1d(pars[3]);
    }
    else
    {
        par->SetPar2a(pars[0]);
        par->SetPar2b(pars[1]);
        par->SetPar2c(pars[2]);
        par->SetPar2d(pars[3]);
    }
    para1[0] = par->GetPar1za();
    para1[1] = par->GetPar1zb();
    para1[2] = par->GetPar1zc();
    para2[0] = par->GetPar1Walk();
    para2[1] = par->GetPar2Walk();
    para2[2] = par->GetPar3Walk();
    para2[3] = par->GetPar4Walk();
    para2[4] = par->GetPar5Walk();
}

void R3BTofdChangePar::changeZ(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[4];
    Double_t para2[5];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    Double_t offset1 = par->GetOffset1();
    Double_t offset2 = par->GetOffset2();
    Double_t ToTOffset1 = par->GetToTOffset1();
    Double_t ToTOffset2 = par->GetToTOffset2();
    Double_t sync = par->GetSync();
    Double_t veff = par->GetVeff();
    Double_t lambda = par->GetLambda();
    if (pars[0] == 0.)
        pars[0] = par->GetPar1za();
    if (pars[1] == 0.)
        pars[1] = par->GetPar1zb();
    if (pars[2] == 0.)
        pars[2] = par->GetPar1zc();
    para1[0] = par->GetPar1a();
    para1[1] = par->GetPar1b();
    para1[2] = par->GetPar1c();
    para1[3] = par->GetPar1d();
    par->SetPar1za(pars[0]);
    par->SetPar1zb(pars[1]);
    par->SetPar1zc(pars[2]);
    para2[0] = par->GetPar1Walk();
    para2[1] = par->GetPar2Walk();
    para2[2] = par->GetPar3Walk();
    para2[3] = par->GetPar4Walk();
    para2[4] = par->GetPar5Walk();
}

void R3BTofdChangePar::changeWalk(Int_t plane, Int_t bar, Int_t pm, Double_t* pars)
{
    Double_t para1[4];
    Double_t para2[3];
    R3BTofdHitModulePar* par = fCal_Par->GetModuleParAt(plane, bar);
    Double_t offset1 = par->GetOffset1();
    Double_t offset2 = par->GetOffset2();
    Double_t ToTOffset1 = par->GetToTOffset1();
    Double_t ToTOffset2 = par->GetToTOffset2();
    Double_t sync = par->GetSync();
    Double_t veff = par->GetVeff();
    Double_t lambda = par->GetLambda();
    if (pm == 1)
    {
        para1[0] = par->GetPar1a();
        para1[1] = par->GetPar1b();
        para1[2] = par->GetPar1c();
        para1[3] = par->GetPar1d();
    }
    else
    {
        para1[0] = par->GetPar2a();
        para1[1] = par->GetPar2b();
        para1[2] = par->GetPar2c();
        para1[3] = par->GetPar2d();
    }
    para2[0] = par->GetPar1za();
    para2[1] = par->GetPar1zb();
    para2[2] = par->GetPar1zc();
    if (pars[0] == 0.)
        pars[0] = par->GetPar1Walk();
    if (pars[1] == 0.)
        pars[1] = par->GetPar2Walk();
    if (pars[2] == 0.)
        pars[2] = par->GetPar3Walk();
    if (pars[3] == 0.)
        pars[3] = par->GetPar4Walk();
    if (pars[4] == 0.)
        pars[4] = par->GetPar5Walk();
    par->SetPar1Walk(pars[0]);
    par->SetPar2Walk(pars[1]);
    par->SetPar3Walk(pars[2]);
    par->SetPar4Walk(pars[3]);
    par->SetPar5Walk(pars[4]);
}

ClassImp(R3BTofdChangePar)
