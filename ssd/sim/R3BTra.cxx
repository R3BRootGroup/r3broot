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

#include "R3BTra.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "R3BLogger.h"
#include "R3BMCStack.h"
#include "R3BTraPoint.h"
#include "TClonesArray.h"
#include "TGeoMCGeometry.h"
#include "TGeoManager.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TVirtualMC.h"

R3BTra::R3BTra()
    : R3BTra("")
{
}

R3BTra::R3BTra(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BTra(geoFile, { trans, rot })
{
}

R3BTra::R3BTra(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BTra", kTRA, geoFile, combi)
    , fTraCollection(new TClonesArray("R3BTraPoint"))
    , fPosIndex(0)
    , kGeoSaved(kFALSE)
    , flGeoPar(new TList())
    , fNbDet(20)
{
    flGeoPar->SetName(GetName());
    ResetParameters();
}

R3BTra::~R3BTra()
{
    if (flGeoPar)
    {
        delete flGeoPar;
    }
    if (fTraCollection)
    {
        fTraCollection->Delete();
        delete fTraCollection;
    }
}

void R3BTra::Initialize()
{
    FairDetector::Initialize();

    LOG(info) << "R3BTra: initialisation";
    LOG(debug) << "R3BTra: Sens. Vol. (McId) " << gMC->VolId("TraLog");
    LOG(debug) << "R3BTra: Sens. Vol. (McId) " << gMC->VolId("Strip");
    LOG(debug) << "R3BTra: Sens. Vol. (McId) " << gMC->VolId("Foot_");

    Char_t buffer[126];
    for (Int_t i = 0; i < fNbDet; i++)
    {
        sprintf(buffer, "Det%i", i + 1);
        LOG(debug) << "Det " << i << " connected to (McId) ---> " << buffer << "  " << gMC->VolId(buffer);
        fdetid[i] = gMC->VolId(buffer);
    }
}

void R3BTra::SetSpecialPhysicsCuts()
{
    LOG(info) << "R3BTra: Adding customized Physics cut ... ";

    if (gGeoManager)
    {
        TGeoMedium* pSi = gGeoManager->GetMedium("silicon");
        if (pSi)
        {
            // Setting processes for Si only
            gMC->Gstpar(pSi->GetId(), "LOSS", 3);
            gMC->Gstpar(pSi->GetId(), "STRA", 1.0);
            gMC->Gstpar(pSi->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pSi->GetId(), "COMP", 1.0);
            gMC->Gstpar(pSi->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pSi->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pSi->GetId(), "BREM", 1.0);
            gMC->Gstpar(pSi->GetId(), "HADR", 1.0);
            gMC->Gstpar(pSi->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pSi->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pSi->GetId(), "MULS", 1.0);
            gMC->Gstpar(pSi->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Si Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: silicon Medium Id " << pSi->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Si
            gMC->Gstpar(pSi->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pSi->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pSi->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pSi->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }
        // <DB> trick to remove too much internal
        // tracking in the Aladin magnet yoke
        TGeoMedium* pFe = gGeoManager->GetMedium("iron");

        if (pFe)
        {
            Double_t cutM = fCutE;
            gMC->Gstpar(pFe->GetId(), "CUTELE", cutM);
            gMC->Gstpar(pFe->GetId(), "DRAY", 0.0);
        }

        TGeoMedium* pLiH = gGeoManager->GetMedium("H2");
        if (pLiH)
        {
            // Setting processes for LiH only
            gMC->Gstpar(pLiH->GetId(), "LOSS", 3);
            gMC->Gstpar(pLiH->GetId(), "STRA", 1.0);
            gMC->Gstpar(pLiH->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pLiH->GetId(), "COMP", 1.0);
            gMC->Gstpar(pLiH->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pLiH->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pLiH->GetId(), "BREM", 1.0);
            gMC->Gstpar(pLiH->GetId(), "HADR", 1.0);
            gMC->Gstpar(pLiH->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pLiH->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pLiH->GetId(), "MULS", 1.0);
            gMC->Gstpar(pLiH->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for LiH Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: LiH Medium Id " << pLiH->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Si
            gMC->Gstpar(pLiH->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pLiH->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pLiH->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pVac = gGeoManager->GetMedium("vacuum");
        if (pVac)
        {
            // Setting processes for Vac only
            gMC->Gstpar(pVac->GetId(), "LOSS", 3);
            gMC->Gstpar(pVac->GetId(), "STRA", 1.0);
            gMC->Gstpar(pVac->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pVac->GetId(), "COMP", 1.0);
            gMC->Gstpar(pVac->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pVac->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pVac->GetId(), "BREM", 1.0);
            gMC->Gstpar(pVac->GetId(), "HADR", 1.0);
            gMC->Gstpar(pVac->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pVac->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pVac->GetId(), "MULS", 1.0);
            gMC->Gstpar(pVac->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Vac Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: Vac Medium Id " << pVac->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Vac
            gMC->Gstpar(pVac->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pVac->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pVac->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pVac->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pVac->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pVac->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pVac->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pVac->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pVac->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pVac->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pGold = gGeoManager->GetMedium("gold");
        if (pGold)
        {
            // Setting processes for Vac only
            gMC->Gstpar(pGold->GetId(), "LOSS", 3);
            gMC->Gstpar(pGold->GetId(), "STRA", 1.0);
            gMC->Gstpar(pGold->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pGold->GetId(), "COMP", 1.0);
            gMC->Gstpar(pGold->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pGold->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pGold->GetId(), "BREM", 1.0);
            gMC->Gstpar(pGold->GetId(), "HADR", 1.0);
            gMC->Gstpar(pGold->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pGold->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pGold->GetId(), "MULS", 1.0);
            gMC->Gstpar(pGold->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Vac Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: Gold Medium Id " << pGold->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Gold
            gMC->Gstpar(pGold->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pGold->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pGold->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pGold->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pGold->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pGold->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pGold->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pGold->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pGold->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pGold->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pM = gGeoManager->GetMedium("mylar");
        if (pM)
        {
            // Setting processes for Mylar only
            gMC->Gstpar(pM->GetId(), "LOSS", 3);
            gMC->Gstpar(pM->GetId(), "STRA", 1.0);
            gMC->Gstpar(pM->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pM->GetId(), "COMP", 1.0);
            gMC->Gstpar(pM->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pM->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pM->GetId(), "BREM", 1.0);
            gMC->Gstpar(pM->GetId(), "HADR", 1.0);
            gMC->Gstpar(pM->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pM->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pM->GetId(), "MULS", 1.0);
            gMC->Gstpar(pM->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Mylar Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: Mylar Medium Id " << pM->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Mylar
            gMC->Gstpar(pM->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pM->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pM->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pM->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pM->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pM->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pM->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pM->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pM->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pM->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pAl = gGeoManager->GetMedium("aluminium");
        if (pAl)
        {
            // Setting processes for Mylar only
            gMC->Gstpar(pAl->GetId(), "LOSS", 3);
            gMC->Gstpar(pAl->GetId(), "STRA", 1.0);
            gMC->Gstpar(pAl->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pAl->GetId(), "COMP", 1.0);
            gMC->Gstpar(pAl->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pAl->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pAl->GetId(), "BREM", 1.0);
            gMC->Gstpar(pAl->GetId(), "HADR", 1.0);
            gMC->Gstpar(pAl->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pAl->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pAl->GetId(), "MULS", 1.0);
            gMC->Gstpar(pAl->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for aluminium Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: Aluminium Medium Id " << pAl->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Al
            gMC->Gstpar(pAl->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pAl->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pAl->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pAl->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pAl->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pAl->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pAl->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pAl->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pAl->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pAl->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pC = gGeoManager->GetMedium("carbon");
        if (pC)
        {
            // Setting processes for Carbon only
            gMC->Gstpar(pC->GetId(), "LOSS", 3);
            gMC->Gstpar(pC->GetId(), "STRA", 1.0);
            gMC->Gstpar(pC->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pC->GetId(), "COMP", 1.0);
            gMC->Gstpar(pC->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pC->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pC->GetId(), "BREM", 1.0);
            gMC->Gstpar(pC->GetId(), "HADR", 1.0);
            gMC->Gstpar(pC->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pC->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pC->GetId(), "MULS", 1.0);
            gMC->Gstpar(pC->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for carbon Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: Carbon Medium Id " << pC->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // C
            gMC->Gstpar(pC->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pC->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pC->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pC->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pC->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pC->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pC->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pC->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pC->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pC->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pHe = gGeoManager->GetMedium("helium");
        if (pHe)
        {
            // Setting processes for Helium only
            gMC->Gstpar(pHe->GetId(), "LOSS", 3);
            gMC->Gstpar(pHe->GetId(), "STRA", 1.0);
            gMC->Gstpar(pHe->GetId(), "PAIR", 1.0);
            gMC->Gstpar(pHe->GetId(), "COMP", 1.0);
            gMC->Gstpar(pHe->GetId(), "PHOT", 1.0);
            gMC->Gstpar(pHe->GetId(), "ANNI", 1.0);
            gMC->Gstpar(pHe->GetId(), "BREM", 1.0);
            gMC->Gstpar(pHe->GetId(), "HADR", 1.0);
            gMC->Gstpar(pHe->GetId(), "DRAY", 1.0);
            gMC->Gstpar(pHe->GetId(), "DCAY", 1.0);
            gMC->Gstpar(pHe->GetId(), "MULS", 1.0);
            gMC->Gstpar(pHe->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for carbon Only
            Double_t cutE = fCutE;

            LOG(info) << "R3BTra: Helium Medium Id " << pHe->GetId() << " Energy Cut-Off : " << cutE << " GeV";

            // Helium
            gMC->Gstpar(pHe->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            gMC->Gstpar(pHe->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            gMC->Gstpar(pHe->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            gMC->Gstpar(pHe->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }
    } //! gGeoManager
}

Int_t R3BTra::GetDetId(Int_t volID)
{
    Int_t detid = -1;

    for (Int_t i = 0; i < fNbDet; i++)
    {
        if (volID == fdetid[i] + 1 || volID == fdetid[i] + 2)
        {
            detid = i;
            return (detid);
        }
    }
    return detid;
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BTra::ProcessHits(FairVolume* vol)
{
    if (gMC->IsTrackEntering())
    {
        fELoss = 0.;
        fTime = gMC->TrackTime() * 1.0e09;
        fLength = gMC->TrackLength();
        gMC->TrackPosition(fPosIn);
        gMC->TrackMomentum(fMomIn);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += gMC->Edep();

    // Set additional parameters at exit of active volume. Create R3BTraPoint.
    if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
    {
        fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
        fVolumeID = vol->getMCid();
        fVolumeID = GetDetId(fVolumeID);
        fDetCopyID = vol->getCopyNo();
        gMC->TrackPosition(fPosOut);
        gMC->TrackMomentum(fMomOut);
        if (fELoss == 0.)
            return kFALSE;

        if (gMC->IsTrackExiting())
        {
            const Double_t* oldpos;
            const Double_t* olddirection;
            Double_t newpos[3];
            Double_t newdirection[3];
            Double_t safety;

            gGeoManager->FindNode(fPosOut.X(), fPosOut.Y(), fPosOut.Z());
            oldpos = gGeoManager->GetCurrentPoint();
            olddirection = gGeoManager->GetCurrentDirection();

            for (Int_t i = 0; i < 3; i++)
            {
                newdirection[i] = -1 * olddirection[i];
            }

            gGeoManager->SetCurrentDirection(newdirection);
            //   TGeoNode *bla = gGeoManager->FindNextBoundary(2);
            safety = gGeoManager->GetSafeDistance();

            gGeoManager->SetCurrentDirection(-newdirection[0], -newdirection[1], -newdirection[2]);

            for (Int_t i = 0; i < 3; i++)
            {
                newpos[i] = oldpos[i] - (3 * safety * olddirection[i]);
            }

            fPosOut.SetX(newpos[0]);
            fPosOut.SetY(newpos[1]);
            fPosOut.SetZ(newpos[2]);
        }

        AddHit(fTrackID,
               fVolumeID,
               fDetCopyID,
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss,
               gMC->TrackPid());

        // Increment number of TraPoints for this track
        R3BStack* stack = dynamic_cast<R3BStack*>(gMC->GetStack());
        stack->AddPoint(kTRA);

        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BTra::BeginEvent() {}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BTra::EndOfEvent()
{
    if (fVerboseLevel)
        Print();
    fTraCollection->Clear();

    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BTra::Register() { FairRootManager::Instance()->Register("TraPoint", GetName(), fTraCollection, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BTra::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
        return fTraCollection;
    else
        return NULL;
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BTra::Print(Option_t* option) const
{
    Int_t nHits = fTraCollection->GetEntriesFast();
    LOG(info) << "R3BTra: " << nHits << " points registered in this event";
}
// ----------------------------------------------------------------------------

// -----   Public method Reset   ----------------------------------------------
void R3BTra::Reset()
{
    fTraCollection->Clear();
    ResetParameters();
}
// ----------------------------------------------------------------------------

// -----   Public method CopyClones   -----------------------------------------
void R3BTra::CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset)
{
    Int_t nEntries = cl1->GetEntriesFast();
    LOG(info) << "R3BTra: " << nEntries << " entries to add";
    TClonesArray& clref = *cl2;
    R3BTraPoint* oldpoint = NULL;
    for (Int_t i = 0; i < nEntries; i++)
    {
        oldpoint = dynamic_cast<R3BTraPoint*>(cl1->At(i));
        Int_t index = oldpoint->GetTrackID() + offset;
        oldpoint->SetTrackID(index);
        new (clref[fPosIndex]) R3BTraPoint(*oldpoint);
        fPosIndex++;
    }
    LOG(info) << "R3BTra: " << cl2->GetEntriesFast() << " merged entries";
}

// -----   Private method AddHit   --------------------------------------------
R3BTraPoint* R3BTra::AddHit(Int_t trackID,
                            Int_t detID,
                            Int_t detCopyID, // added by Marc
                            TVector3 posIn,
                            TVector3 posOut,
                            TVector3 momIn,
                            TVector3 momOut,
                            Double_t time,
                            Double_t length,
                            Double_t eLoss,
                            Int_t pdgcode)
{
    TClonesArray& clref = *fTraCollection;
    Int_t size = clref.GetEntriesFast();
    if (fVerboseLevel > 1)
        LOG(info) << "R3BTra: Adding Point at (" << posIn.X() << ", " << posIn.Y() << ", " << posIn.Z()
                  << ") cm,  detector " << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";
    return new (clref[size]) R3BTraPoint(trackID,
                                         detID,
                                         detCopyID,
                                         posIn,
                                         posOut, // detCopyID added by Marc
                                         momIn,
                                         momOut,
                                         time,
                                         length,
                                         eLoss,
                                         pdgcode);
}

Bool_t R3BTra::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("TraLog") || TString(name).Contains("Strip") || TString(name).Contains("Alpide") ||
        TString(name).Contains("Foot_"))
    {
        // LOG(info) << "Found TRA geometry from ROOT file: " << name;
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BTra)
