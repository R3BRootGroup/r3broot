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

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairVolume.h"

#include "R3BAlpide.h"
#include "R3BAlpideGeometry.h"
#include "R3BAlpideMappingPar.h"
#include "R3BAlpidePoint.h"
#include "R3BLogger.h"
#include "R3BMCStack.h"

#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TVirtualMC.h"

R3BAlpide::R3BAlpide()
    : R3BAlpide("")
{
}

R3BAlpide::R3BAlpide(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot)
    : R3BAlpide(geoFile, { trans, rot })
{
}

R3BAlpide::R3BAlpide(const TString& geoFile, const TGeoCombiTrans& combi)
    : R3BDetector("R3BAlpide", kTRA, geoFile, combi)
    , fAlpidePoint(new TClonesArray("R3BAlpidePoint"))
    , fMap_Par(nullptr)
    , fAlpideGeo(nullptr)
    , fGeoversion(2024)
{
    ResetParameters();
}

R3BAlpide::~R3BAlpide()
{
    R3BLOG(debug1, "");
    if (fAlpidePoint)
    {
        fAlpidePoint->Delete();
        delete fAlpidePoint;
    }
}

void R3BAlpide::SetParameter()
{
    // Parameter Container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    fMap_Par = dynamic_cast<R3BAlpideMappingPar*>(rtdb->getContainer("alpideMappingPar"));
    R3BLOG_IF(warning, !fMap_Par, "Container alpideMappingPar not found");
    //--- Parameter Container ---
    if (fMap_Par)
    {
        fGeoversion = fMap_Par->GetGeoVersion();
    }
    R3BLOG(info, "Geometry version: " << fGeoversion);
    return;
}

void R3BAlpide::Initialize()
{
    FairDetector::Initialize();

    R3BLOG(info, " ");
    R3BLOG(debug, "Sens. Vol. (McId) " << TVirtualMC::GetMC()->VolId("Alpide"));

    SetParameter();
    fAlpideGeo = R3BAlpideGeometry::Instance();
    R3BLOG_IF(error, !fAlpideGeo->Init(fGeoversion), "Alpide geometry " << fGeoversion << " not found");
}

void R3BAlpide::SetSpecialPhysicsCuts()
{
    R3BLOG(info, "Adding customized Physics cut.");

    if (gGeoManager)
    {
        TGeoMedium* pSi = static_cast<TGeoMedium*>(gGeoManager->GetMedium("silicon"));
        if (pSi)
        {
            // Setting processes for Si only
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Si Only
            Double_t cutE = fCutE;

            R3BLOG(info, "silicon Medium Id " << pSi->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // Si
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pSi->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }
        // <DB> trick to remove too much internal
        // tracking in the Aladin magnet yoke
        TGeoMedium* pFe = static_cast<TGeoMedium*>(gGeoManager->GetMedium("iron"));

        if (pFe)
        {
            Double_t cutM = fCutE;
            TVirtualMC::GetMC()->Gstpar(pFe->GetId(), "CUTELE", cutM);
            TVirtualMC::GetMC()->Gstpar(pFe->GetId(), "DRAY", 0.0);
        }

        TGeoMedium* pLiH = static_cast<TGeoMedium*>(gGeoManager->GetMedium("H2"));
        if (pLiH)
        {
            // Setting processes for LiH only
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for LiH Only
            Double_t cutE = fCutE;

            R3BLOG(info, "LiH Medium Id " << pLiH->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // Si
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pLiH->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pVac = static_cast<TGeoMedium*>(gGeoManager->GetMedium("vacuum"));
        if (pVac)
        {
            // Setting processes for Vac only
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Vac Only
            Double_t cutE = fCutE;

            R3BLOG(info, "Vacuum Medium Id " << pVac->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // Vac
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pVac->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pGold = static_cast<TGeoMedium*>(gGeoManager->GetMedium("gold"));
        if (pGold)
        {
            // Setting processes for Vac only
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Vac Only
            Double_t cutE = fCutE;

            R3BLOG(info, "Gold Medium Id " << pGold->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // Gold
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pGold->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pM = static_cast<TGeoMedium*>(gGeoManager->GetMedium("mylar"));
        if (pM)
        {
            // Setting processes for Mylar only
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for Mylar Only
            Double_t cutE = fCutE;

            R3BLOG(info, "Mylar Medium Id " << pM->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // Mylar
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pM->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pAl = static_cast<TGeoMedium*>(gGeoManager->GetMedium("aluminium"));
        if (pAl)
        {
            // Setting processes for Mylar only
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for aluminium Only
            Double_t cutE = fCutE;

            R3BLOG(info, "Aluminium Medium Id " << pAl->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // Al
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pAl->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pC = static_cast<TGeoMedium*>(gGeoManager->GetMedium("carbon"));
        if (pC)
        {
            // Setting processes for Carbon only
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for carbon Only
            Double_t cutE = fCutE;

            R3BLOG(info, "Carbon Medium Id " << pC->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // C
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pC->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }

        TGeoMedium* pHe = static_cast<TGeoMedium*>(gGeoManager->GetMedium("helium"));
        if (pHe)
        {
            // Setting processes for Helium only
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "LOSS", 3);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "STRA", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "PAIR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "COMP", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "PHOT", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "ANNI", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "BREM", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "HADR", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "DRAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "DCAY", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "MULS", 1.0);
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "RAYL", 1.0);

            // Setting Energy-CutOff for carbon Only
            Double_t cutE = fCutE;

            R3BLOG(info, "Helium Medium Id " << pHe->GetId() << " Energy Cut-Off : " << cutE << " GeV");

            // Helium
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "CUTGAM", cutE); /** gammas (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "CUTELE", cutE); /** electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "CUTNEU", cutE); /** neutral hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "CUTHAD", cutE); /** charged hadrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "CUTMUO", cutE); /** muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "BCUTE", cutE);  /** electron bremsstrahlung (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "BCUTM", cutE);  /** muon and hadron bremsstrahlung(GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "DCUTE", cutE);  /** delta-rays by electrons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "DCUTM", cutE);  /** delta-rays by muons (GeV)*/
            TVirtualMC::GetMC()->Gstpar(pHe->GetId(), "PPCUTM", -1.);  /** direct pair production by muons (GeV)*/
        }
    } //! gGeoManager
}

// -----   Public method ProcessHits  --------------------------------------
Bool_t R3BAlpide::ProcessHits(FairVolume*)
{
    if (TVirtualMC::GetMC()->IsTrackEntering())
    {
        fELoss = 0.;
        fTime = TVirtualMC::GetMC()->TrackTime() * 1.0e09;
        fLength = TVirtualMC::GetMC()->TrackLength();
        TVirtualMC::GetMC()->TrackPosition(fPosIn);
        TVirtualMC::GetMC()->TrackMomentum(fMomIn);
    }

    // Sum energy loss for all steps in the active volume
    fELoss += TVirtualMC::GetMC()->Edep();

    // Set additional parameters at exit of active volume. Create R3BAlpidePoint.
    if (TVirtualMC::GetMC()->IsTrackExiting() || TVirtualMC::GetMC()->IsTrackStop() ||
        TVirtualMC::GetMC()->IsTrackDisappeared())
    {
        fTrackID = TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();
        R3BLOG(debug, TVirtualMC::GetMC()->CurrentVolPath());
        fBarrelID = fAlpideGeo->GetBarrelId(TVirtualMC::GetMC()->CurrentVolPath());
        fSensorID = fAlpideGeo->GetSensorId(TVirtualMC::GetMC()->CurrentVolPath());

        TVirtualMC::GetMC()->TrackPosition(fPosOut);
        TVirtualMC::GetMC()->TrackMomentum(fMomOut);
        if (fELoss == 0.)
        {
            return kFALSE;
        }

        AddHit(fTrackID,
               fBarrelID,
               fSensorID,
               TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
               TVector3(fPosOut.X(), fPosOut.Y(), fPosOut.Z()),
               TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()),
               TVector3(fMomOut.Px(), fMomOut.Py(), fMomOut.Pz()),
               fTime,
               fLength,
               fELoss,
               TVirtualMC::GetMC()->TrackPid());

        // Increment number of AlpidePoints for this track
        R3BStack* stack = static_cast<R3BStack*>(TVirtualMC::GetMC()->GetStack());
        stack->AddPoint(kTRA);
        ResetParameters();
    }

    return kTRUE;
}

// -----   Public method EndOfEvent   -----------------------------------------
void R3BAlpide::EndOfEvent()
{
    if (fVerboseLevel)
    {
        Print();
    }
    Reset();
}
// ----------------------------------------------------------------------------

// -----   Public method Register   -------------------------------------------
void R3BAlpide::Register() { FairRootManager::Instance()->Register("AlpidePoint", GetName(), fAlpidePoint, kTRUE); }
// ----------------------------------------------------------------------------

// -----   Public method GetCollection   --------------------------------------
TClonesArray* R3BAlpide::GetCollection(Int_t iColl) const
{
    if (iColl == 0)
    {
        return fAlpidePoint;
    }
    else
    {
        return nullptr;
    }
}
// ----------------------------------------------------------------------------

// -----   Public method Print   ----------------------------------------------
void R3BAlpide::Print(Option_t*) const
{
    Int_t nHits = fAlpidePoint->GetEntriesFast();
    R3BLOG(info, nHits << " points registered in this event");
}

// -----   Public method Reset   ----------------------------------------------
void R3BAlpide::Reset()
{
    fAlpidePoint->Clear();
    ResetParameters();
}

// -----   Private method AddHit   --------------------------------------------
R3BAlpidePoint* R3BAlpide::AddHit(Int_t trackID,
                                  Int_t detID,
                                  Int_t detCopyID,
                                  TVector3 posIn,
                                  TVector3 posOut,
                                  TVector3 momIn,
                                  TVector3 momOut,
                                  Double_t time,
                                  Double_t length,
                                  Double_t eLoss,
                                  Int_t pdgcode)
{
    TClonesArray& clref = *fAlpidePoint;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size])
        R3BAlpidePoint(trackID, detID, detCopyID, posIn, posOut, momIn, momOut, time, length, eLoss, pdgcode);
}

// -----   Public method CheckIfSensitive   -----------------------------------
Bool_t R3BAlpide::CheckIfSensitive(std::string name)
{
    if (TString(name).Contains("Alpide"))
    {
        // LOG(info) << "Found geometry from ROOT file: " << name;
        return kTRUE;
    }
    return kFALSE;
}

ClassImp(R3BAlpide)
