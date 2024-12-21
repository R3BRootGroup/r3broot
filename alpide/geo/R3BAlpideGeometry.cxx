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

#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoNavigator.h>
#include <TGeoVolume.h>
#include <TMath.h>
#include <TRotation.h>
#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>
#include <iostream>

#include "R3BAlpideGeometry.h"
#include "R3BLogger.h"

#include <boost/regex.hpp>

R3BAlpideGeometry* R3BAlpideGeometry::Instance()
{
    // Returns singleton instance
    static thread_local R3BAlpideGeometry instance;
    return &instance;
}

R3BAlpideGeometry::R3BAlpideGeometry()
    : TObject()
    , IsInitialize(kFALSE)
    , fGeometryVersion(2024)
    , fNbSensor(108)
{
}

bool R3BAlpideGeometry::Init(Int_t version)
{
    if (!IsInitialize)
        IsInitialize = kTRUE;
    else
        return kTRUE;

    R3BLOG(info, "");
    TString geoPath = gSystem->Getenv("VMCWORKDIR");
    geoPath += "/geometry/";

    switch (version)
    {
        case 2028:
            // Two barrels + one bending barrel
            geoPath += "target_area_alpide_barrels_bending_v28.geo.root";
            fNbSensor = 363;
            fGeometryVersion = version;
            break;

        case 2026:
            // Two barrels
            geoPath += "target_area_alpide_barrel_v26.geo.root";
            fNbSensor = 363;
            fGeometryVersion = version;
            break;

        case 2024:
            // Two arms
            geoPath += "target_area_alpide_twoarms_v2024.2.geo.root";
            fNbSensor = 108;
            fGeometryVersion = version;
            break;

        case 202402:
            // Two arms
            geoPath += "target_area_alpide_singlearms_v202402.geo.root";
            fNbSensor = 24;
            fGeometryVersion = version;
            break;

        case 202210:
            // 6 ALPIDEs in a telescope configuration
            geoPath += "tracking_alpide_cern_202210.geo.root";
            fNbSensor = 31;
            fGeometryVersion = version;
            break;

        default:
            R3BLOG(error, "Unsupported geometry version: " << version);
            return kFALSE;
    }

    if (gGeoManager && strcmp(gGeoManager->GetTopVolume()->GetName(), "cave") == 0)
    {
        // Already set up (MC mode)
        R3BLOG(info, "Using existing geometry");
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
        gGeoManager = new TGeoManager();
    gGeoManager->SetTopVolume(v);
    fIsSimulation = kFALSE;
    return kTRUE;
}

R3BAlpideGeometry::~R3BAlpideGeometry()
{
    R3BLOG(debug, "");
    if (gGeoManager)
        delete gGeoManager;
    if (f)
    {
        f->Close();
    }
}

const TVector3& R3BAlpideGeometry::GetAngles(Int_t iD)
{
    static std::map<int, TVector3> cache;
    Double_t local[3] = { 0, 0, 0 };
    Double_t master[3];

    const static TVector3 invalid(NAN, NAN, NAN);
    const char* nameVolume;
    if (cache.count(iD))
        return cache[iD];

    if (iD >= 1 && iD <= fNbSensor)
    {
        nameVolume = GetSensorVolumePath(iD);
        gGeoManager->CdTop();

        if (gGeoManager->CheckPath(nameVolume))
            gGeoManager->cd(nameVolume);
        else
        {
            R3BLOG(error, "Invalid sensor path: " << nameVolume);
            return invalid;
        }
        gGeoManager->LocalToMaster(local, master);
    }
    else
    {
        R3BLOG(error, "Invalid sensorId: " << iD);
        return invalid;
    }
    return cache[iD] = master;
}

// Rotation matrix
const TRotation R3BAlpideGeometry::GetRotation(Int_t iD)
{
    static std::map<int, TRotation> cache;
    Double_t local[3] = { 0, 0, 0 };
    Double_t master[3];

    TRotation tt;

    static TRotation invalid;
    invalid.SetToIdentity();

    R3BLOG(debug, "SensorId: " << iD);

    const char* nameVolume;
    if (cache.count(iD))
        return cache[iD];

    if (iD >= 1 && iD <= fNbSensor)
    {
        nameVolume = GetSensorVolumePath(iD);
        gGeoManager->CdTop();

        if (gGeoManager->CheckPath(nameVolume))
            gGeoManager->cd(nameVolume);
        else
        {
            R3BLOG(error, "Invalid sensor path: " << nameVolume);
            return invalid;
        }
        gGeoManager->LocalToMaster(local, master);

        Double_t* rMatrix = gGeoManager->GetCurrentMatrix()->GetRotationMatrix();
        //  std::cout <<"rMatrix "<< rMatrix[0] <<" "<< rMatrix[1] <<" "<< rMatrix[2] << std::endl;
        //  std::cout <<"rMatrix "<< rMatrix[3] <<" "<< rMatrix[4] <<" "<< rMatrix[5] << std::endl;
        //  std::cout <<"rMatrix "<< rMatrix[6] <<" "<< rMatrix[7] <<" "<< rMatrix[8] << std::endl;

        TVector3 rx;
        rx.SetX(rMatrix[0]);
        rx.SetY(rMatrix[3]);
        rx.SetZ(rMatrix[6]);

        TVector3 ry;
        ry.SetX(rMatrix[1]);
        ry.SetY(rMatrix[4]);
        ry.SetZ(rMatrix[7]);

        TVector3 rz;
        rz.SetX(rMatrix[2]);
        rz.SetY(rMatrix[5]);
        rz.SetZ(rMatrix[8]);

        tt.RotateAxes(rx, ry, rz);

        R3BLOG(debug, "Rx " << tt.XX() << ", " << tt.XY() << ", " << tt.XZ());
        R3BLOG(debug, "Ry " << tt.YX() << ", " << tt.YY() << ", " << tt.YZ());
        R3BLOG(debug, "Rz " << tt.ZX() << ", " << tt.ZY() << ", " << tt.ZZ());
    }
    else
    {
        R3BLOG(error, "Invalid sensorId: " << iD);
        return invalid;
    }
    return cache[iD] = tt;
}

// Translation vector
const TVector3& R3BAlpideGeometry::GetTranslation(Int_t iD)
{
    static std::map<int, TVector3> cache;
    Double_t local[3] = { 0, 0, 0 };
    Double_t master[3];
    Double_t* trans;

    R3BLOG(debug, "SensorId: " << iD);

    const static TVector3 invalid(NAN, NAN, NAN);
    const char* nameVolume;
    if (cache.count(iD))
        return cache[iD];

    if (iD >= 1 && iD <= fNbSensor)
    {
        nameVolume = GetSensorVolumePath(iD);
        gGeoManager->CdTop();

        if (gGeoManager->CheckPath(nameVolume))
            gGeoManager->cd(nameVolume);
        else
        {
            R3BLOG(error, "Invalid sensor path: " << nameVolume);
            return invalid;
        }
        gGeoManager->LocalToMaster(local, master);
        trans = gGeoManager->GetCurrentMatrix()->GetTranslation();
        // std::cout << master[0] <<" "<< master[1] <<" "<< master[2] << std::endl;
        // std::cout << trans[0] <<" "<< trans[1] <<" "<< trans[2] << std::endl;
    }
    else
    {
        R3BLOG(error, "Invalid sensorId: " << iD);
        return invalid;
    }

    return cache[iD] = trans;
}

void R3BAlpideGeometry::GetAngles(Int_t iD, Double_t* polar, Double_t* azimuthal, Double_t* rho)
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

const char* R3BAlpideGeometry::GetSensorVolumePath(Int_t iD)
{
    static char nameVolume[300];
    Int_t sid = 0;
    Int_t bartype = 0;
    Int_t layertype = 0;

    R3BLOG(debug, "SensorId: " << iD);

    if (iD >= 1 && iD <= fNbSensor)
    {
        if (fGeometryVersion == 2026 || fGeometryVersion == 2028)
        {
            if (iD <= 153)
            {
                bartype = 1;
                layertype = (iD - 1) / 9;
                sid = iD - layertype * 9;
                layertype++;
            }
            else
            {
                bartype = 2;
                sid = iD - 154;
                layertype = sid / 10;
                sid = iD - 153 - layertype * 10;
                layertype++;
            }
        }
        else
        {

            bartype = 1;
            layertype = (iD - 1) / 6;
            sid = iD - layertype * 6;
            layertype++;
        }

        sprintf(nameVolume, "/cave_1/TargetChamberWorld_0/Multilayer_%i_%i/Alpide_%i", bartype, layertype, sid);
    }
    else
    {
        R3BLOG(error, "Invalid sensorId: " << iD);
    }

    return nameVolume;
}

int R3BAlpideGeometry::GetBarrelId(const char* volumePath)
{
    Int_t barID = 0;
    static auto restr = "Multilayer_([0-9]+)_([0-9]+)/Alpide_([0-9]+)";
    static auto re = boost::regex(restr, boost::regex::extended);
    boost::cmatch m;
    if (!boost::regex_search(volumePath, m, re))
    {
        R3BLOG(error,
               "\"" << volumePath
                    << "\"\n"
                       "does not match RE \""
                    << restr << "\".\n");
        return barID;
    }

    barID = std::stoi(m[1].str()); // converting to int the barrel type

    R3BLOG(debug, "Barrel ID: " << barID);

    return barID;
}

int R3BAlpideGeometry::GetSensorId(const char* volumePath)
{
    Int_t sensorId = 0;
    Int_t barID = 0;
    Int_t layerID = 0;
    Int_t alpideID = 0;
    static auto restr = "Multilayer_([0-9]+)_([0-9]+)/Alpide_([0-9]+)";
    static auto re = boost::regex(restr, boost::regex::extended);
    boost::cmatch m;
    if (!boost::regex_search(volumePath, m, re))
    {
        R3BLOG(error,
               "\"" << volumePath
                    << "\"\n"
                       "does not match RE \""
                    << restr << "\".\n");
        return sensorId;
    }

    barID = std::stoi(m[1].str());    // converting to int the barrel type
    layerID = std::stoi(m[2].str());  // converting to int the layer type
    alpideID = std::stoi(m[3].str()); // converting to int the alpide type

    if (fGeometryVersion == 2026 || fGeometryVersion == 2028)
    {
        if (barID == 2)
        {
            sensorId = (layerID - 1) * 10 + alpideID;
            sensorId += 153;
        }
        else
        {
            sensorId = (layerID - 1) * 9 + alpideID;
        }
    }
    else
    {
        sensorId = (layerID - 1) * 6 + alpideID;
    }

    R3BLOG(debug,
           "Barrel ID: " << barID << ", multiflex ID: " << layerID << ", sensorID in multiflex: " << alpideID
                         << ", sensorID: " << sensorId);

    return sensorId;
}

ClassImp(R3BAlpideGeometry)
