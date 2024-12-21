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

#include <R3BCalifaGeometry.h>
#include <R3BLogger.h>

#include <FairLogger.h>

#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoNavigator.h>
#include <TGeoVolume.h>
#include <TMath.h>
#include <TSystem.h>
#include <vector>

#include <boost/regex.hpp>

R3BCalifaGeometry* R3BCalifaGeometry::Instance()
{
    // Returns singleton instance
    static thread_local R3BCalifaGeometry instance;
    return &instance;
}

R3BCalifaGeometry::R3BCalifaGeometry()
    : TObject()
{
}

bool R3BCalifaGeometry::Init(int version)
{
    if (!IsInitialize)
    {
        IsInitialize = kTRUE;
    }
    else
    {
        return kTRUE;
    }

    R3BLOG(info, "");
    TString geoPath = gSystem->Getenv("VMCWORKDIR");
    geoPath += "/geometry/";

    switch (version)
    {
        case 2019:
            // BARREL: demonstrator from 2019
            geoPath += "califa_s444.geo.root";
            fNumCrystals = 10000;
            fGeometryVersion = version;
            break;

        case 2020:
            // Half BARREL+ 6 IPHOS sectors
            geoPath += "califa_v2019.11.geo.root";
            fNumCrystals = 4864;
            fGeometryVersion = version;
            break;

        case 2021:
            // s455, S515, S509, S522 experiments: Half Barrel + Full IPHOS
            geoPath += "califa_v2021.3.geo.root";
            fNumCrystals = 4864;
            fGeometryVersion = version;
            break;

        case 2024:
            // S118, S091 experiments: Half Barrel (extended) + Full IPHOS + Full CEPA
            geoPath += "califa_v2024.1.geo.root";
            fNumCrystals = 5088;
            fGeometryVersion = version;
            break;

        default:
            // Full Barrel + Full IPHOS
            geoPath += "califa_full.geo.root";
            R3BLOG(warn,
                   "Unsupported geometry version: " << version << ", so standard full configuration will be used.");
            // return kFALSE;
    }

    if (gGeoManager && strcmp(gGeoManager->GetTopVolume()->GetName(), "cave") == 0)
    {
        // Already set up (MC mode)
        R3BLOG(info, "Using existing geometry");
        // LOG(info) << "R3BCalifaGeometry::Opened geometry file " << geoPath;
        fIsSimulation = kTRUE;
        return kTRUE;
    }

    // Stand alone mode
    R3BLOG(info, "Open geometry file " << geoPath << " for analysis.");
    f = new TFile(geoPath, "READ");
    TGeoVolume* v = dynamic_cast<TGeoVolume*>(f->Get("TOP"));
    if (!v)
    {
        R3BLOG(error, "Could not open geometry file, No TOP volume");
        return kFALSE;
    }

    v->SetName("cave");
    if (!gGeoManager)
    {
        gGeoManager = new TGeoManager();
    }
    gGeoManager->SetTopVolume(v);
    fIsSimulation = kFALSE;
    return kTRUE;
}

R3BCalifaGeometry::~R3BCalifaGeometry()
{
    R3BLOG(debug1, "");
    if (gGeoManager)
    {
        delete gGeoManager;
    }
    if (f)
    {
        f->Close();
    }
}

const TVector3& R3BCalifaGeometry::GetAngles(int iD)
{
    static std::map<int, TVector3> cache;
    Double_t local[3] = { 0, 0, 0 };
    Double_t master[3];
    const static TVector3 invalid(NAN, NAN, NAN);
    std::string nameVolume;
    if (cache.count(iD))
    {
        return cache[iD];
    }

    // SOLUTION FOR DOUBLE READING CHANNELS
    if (iD > fNumCrystals / 2 && iD <= fNumCrystals)
        iD = iD - fNumCrystals / 2; // for double reading crystals (crystals from 1 to max.[2544])

    if (iD >= 1 && iD <= fNumCrystals / 2)
    {
        nameVolume = GetCrystalVolumePath(iD);

        gGeoManager->CdTop();

        if (gGeoManager->CheckPath(nameVolume.c_str()))
            gGeoManager->cd(nameVolume.c_str());
        else
        {
            R3BLOG(error, "Invalid crystal path: " << nameVolume);
            return invalid;
        }
        gGeoManager->LocalToMaster(local, master);
    }
    else
    {
        R3BLOG(error, "Invalid crystalId: " << iD);
        return invalid;
    }

    cache[iD] = master - fRefPoint;

    return cache[iD];
}

void R3BCalifaGeometry::GetAngles(int iD, double* polar, double* azimuthal, double* rho)
{
    auto& masterV = this->GetAngles(iD);
    *polar = masterV.Theta();
    *azimuthal = masterV.Phi();
    *rho = masterV.Mag();
    if (std::isnan(*polar) || std::isnan(*azimuthal) || std::isnan(*rho))
    {
        R3BLOG(error, " returns NaN");
    }
}

std::string R3BCalifaGeometry::GetCrystalVolumePath(int iD)
{
    int alveolusCopy = -1;
    int alvType = -1;
    int cryType = -1;

    // SOLUTION FOR DOUBLE READING CHANNELS
    if (iD > fNumCrystals / 2 && iD <= fNumCrystals)
    {
        iD = iD - fNumCrystals / 2; // for double reading crystals (crystals from 1 to max.[2544])
    }

    std::ostringstream oss;
    static constexpr int kMaxNameAlv = 23;
    std::string name_Alv[kMaxNameAlv] = { "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
                                          "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23" };

    if (iD >= 1 && iD <= fNumCrystals / 2)
    {
        if (iD < 33)
        {                          // First ring in BARREL (single crystal per alveolus)
            alvType = 1;           // Alveolus type 1
            alveolusCopy = iD - 1; // copy from 0 to 31
            cryType = 1;           // Only one crystal type (1)
        }
        else if (iD < 2337)
        {                                                                // All 4-crystals alveoli in BARREL and iPhos
            alvType = (Int_t)((iD - 33) / 128) + 2;                      // Alveolus type (2, ..., 19)
            alveolusCopy = (Int_t)((iD - 33 - (alvType - 2) * 128) / 4); // copy from 0 to 31
            cryType = iD - 33 - (alvType - 2) * 128 - alveolusCopy * 4 + 1; // Four crystal types (1,2,3,4)
        }
        else if (iD < 2433)
        {                                                                     // 3-crystals alveoli in iPhos
            alvType = (Int_t)((iD - 2337) / 24) + 20;                         // Alveolus type (20, 21, 22, 23)
            alveolusCopy = (Int_t)((iD - 2337 - (alvType - 20) * 24) / 3);    // copy from 0 to 7
            cryType = iD - 2337 - (alvType - 20) * 24 - alveolusCopy * 3 + 1; // Three crystal types (1,2,3)
        }
        else if (iD < 2465) // CEPA_CsI, four crystals alveoli (iD from 2433 to 2464)
        {
            alvType = 1;
            alveolusCopy = (Int_t)((iD - 2433) / 4);    // copy from 0 to 7
            cryType = iD - 2433 - alveolusCopy * 4 + 1; // Four crystal types (1,2,3,4)
        }
        else if (iD < 2513) // CEPA_CsI, six crystals alveoli (iD from 2465 to 2512)
        {
            alvType = 2;
            alveolusCopy = (Int_t)((iD - 2465) / 6);    // copy from 0 to 7
            cryType = iD - 2465 - alveolusCopy * 6 + 1; // Four crystal types (1,2,3,4,5,6)
        }
        else if (iD < 2545) // CEPA_CsI, four crystals alveoli (iD from 2513 to 2544)
        {
            alvType = 3;
            alveolusCopy = (Int_t)((iD - 2513) / 4);    // copy from 0 to 7
            cryType = iD - 2513 - alveolusCopy * 4 + 1; // Four crystal types (1,2,3,4)
        }

        if (iD < 2433)
        {
            oss << "/cave_1/CalifaWorld_0/Alveolus_" << name_Alv[alvType - 1] << "_" << alveolusCopy << "/InnerAlv_"
                << name_Alv[alvType - 1] << "_0/WrapCry_" << name_Alv[alvType - 1] << "_" << cryType << "_0/Crystal_"
                << name_Alv[alvType - 1] << "_" << cryType << "_0";
        }
        else
        {
            oss << "/cave_1/CalifaWorld_0/Alveolus_CCSI_" << name_Alv[alvType - 1] << "_" << alveolusCopy
                << "/InnerAlv_CCSI_" << name_Alv[alvType - 1] << "_0/WrapCry_CCSI_" << name_Alv[alvType - 1] << "_"
                << cryType << "_0/Crystal_CCSI_" << name_Alv[alvType - 1] << "_" << cryType << "_0";
        }
    }
    else
    {
        R3BLOG(error, "Invalid crystalId: " << iD);
    }

    return oss.str();
}

double R3BCalifaGeometry::GetDistanceThroughCrystals(TVector3& startVertex,
                                                     TVector3& direction,
                                                     TVector3* hitPos,
                                                     int* numCrystals,
                                                     int* crystalIds)
{
    int maxNumCrystals = 0;

    if (numCrystals != NULL && crystalIds != NULL)
    {
        maxNumCrystals = *numCrystals;
        *numCrystals = 0;
    }

    TGeoNode* n;

    gGeoManager->InitTrack(startVertex.X(),
                           startVertex.Y(),
                           startVertex.Z(),
                           direction.X() / direction.Mag(),
                           direction.Y() / direction.Mag(),
                           direction.Z() / direction.Mag());

    Double_t distance = 0.;
    const Double_t* pos;
    bool inCrystal = false, wasInCrystal = false;
    TString nodeName;

    // find the distance to the next boundary and then extrapolate the current point/direction
    // with this distance making sure that the boundary was crossed. Finally the goal would be to find the next state
    // after crossing the boundary
    while ((n = gGeoManager->FindNextBoundaryAndStep()))
    {
        nodeName = n->GetName();

        if (inCrystal)
            distance += gGeoManager->GetStep();

        inCrystal = nodeName.BeginsWith("WrapCry_"); // can't enter to Crystal volume

        if (inCrystal && maxNumCrystals != 0)
        {
            Int_t cid = GetCrystalId(gGeoManager->GetPath());
            if (cid != -1 && (*numCrystals == 0 || cid != crystalIds[(*numCrystals) - 1]))
            {
                crystalIds[(*numCrystals)++] = cid;
                maxNumCrystals--;
            }
        }

        if (hitPos != NULL && inCrystal && !wasInCrystal)
        {
            pos = gGeoManager->GetCurrentPoint();
            hitPos->SetXYZ(pos[0], pos[1], pos[2]);
            wasInCrystal = true;
        }
    }

    return distance;
}

// ROOT (or pyroot) is unable to automatically load dependencies of so's for whatever reason.
// so we try to force loading the lib during var initialisation.
// yes, this is ugly.
// static bool _dummy_load_libboost_regex = []() {
//    char buf[100];
//    snprintf(buf, 100, "%s/lib/libboost_regex.so", getenv("SIMPATH"));
//    return gSystem->Load(buf);
//}();

int R3BCalifaGeometry::GetCrystalId(const std::string volumePath)
{ /*
   * Yes, regex matching in the stepping action is probably slow.
   * But if we cared about speed, we would not use string processing here
   * in the first place.
   * At least, this way we have well defined behaviour if something fails,
   * instead of possible out of bounds access consequences (nasal demons etc).
   * Also, it hopefully works with the 2019/s444 geo file too, now.
   */
    static auto restr = "(Alveolus|Alveolus_CCSI)_([0-9]+)_([0-9]+).*(Crystal|Crystal_CCSI)_[^_]+_([0-9]+)_";
    static auto re = boost::regex(restr, boost::regex::extended);
    boost::cmatch m;
    if (!boost::regex_search(volumePath.c_str(), m, re))
    {
        R3BLOG(error,
               "\"" << volumePath
                    << "\"\n"
                       "does not match RE \""
                    << restr << "\".\n");
        return 0;
    }

    int crystalId = 0;
    bool isCCSI = m[1].str() == "Alveolus_CCSI"; // Adding CEPA CSI
    int alvType = std::stoi(m[2].str());         // converting to int the alveolus type
    int alveolusCopy = std::stoi(m[3].str());    // converting to int the alveolus copy
    int cryType = std::stoi(m[5].str());         // converting to int the crystal type
    bool invalid = kFALSE;

    // cryType runs from 1 to 4 while alvType runs from 1 to 23, otherwise invaled
    if (cryType < 1 || cryType > 4 || alvType < 1 || alvType > 23)
        invalid = kTRUE;
    // EXCEPTION: cryType runs up to 6 in CEPA CSI alveoli 2
    if (isCCSI && alvType == 2 && (cryType == 5 || cryType == 6))
        invalid = kFALSE;
    if (invalid)
    {
        R3BLOG(error, "Wrong crystal numbers (1)");
        LOG(info) << "---- cryType: " << cryType << "   alvType: " << alvType;
        LOG(info) << "path=" << volumePath;
        return 0;
    }

    if (isCCSI)
    {
        if (alvType == 1)
            crystalId = 2433 + alveolusCopy * 4 + (cryType - 1); // four crystal per alveolus (2433 to 2464)
        if (alvType == 2)
            crystalId = 2465 + alveolusCopy * 6 + (cryType - 1); // six crystal per alveolus (2465 to 2512)
        if (alvType == 3)
            crystalId = 2513 + alveolusCopy * 4 + (cryType - 1); // four crystal per alveolus (2513 to 2544)
    }
    else
    {
        if (alvType == 1)
            crystalId = 1 + alveolusCopy; // first alveoli ring, one crystal per alveolus
        else if (alvType < 20)
            crystalId = 33 + (alvType - 2) * 128 + alveolusCopy * 4 + (cryType - 1); // four crystal per alveolus
        else
            crystalId = 2337 + (alvType - 20) * 24 + alveolusCopy * 3 + (cryType - 1); // three crystal per alveolus
    }

    if (crystalId < 1 || crystalId > fNumCrystals / 2)
    { // crystalId runs from 1 to fNumCrystals/2
        R3BLOG(error, "Wrong crystal numbers (2)");
        LOG(info) << "---- crystalId: " << crystalId;
        return 0;
    }

    return crystalId;
}

ClassImp(R3BCalifaGeometry)
