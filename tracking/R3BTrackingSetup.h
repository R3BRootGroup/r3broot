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

#ifndef R3B_TRACKING_SETUP
#define R3B_TRACKING_SETUP

#include "R3BTrackingDetector.h"
#include <map>
#include <string>
#include <vector>

class R3BTrackingSetup
{
  public:
    R3BTrackingSetup();
    virtual ~R3BTrackingSetup();

    void AddDetector(const std::string& name,
                     EDetectorType type,
                     const std::string& geoParName,
                     const std::string& dataName = "");

    R3BTrackingDetector* GetByName(const std::string& name);

    R3BTrackingDetector* GetFirstByType(const EDetectorType& type);

    void Init();

    void SetParContainers();

    void CopyHits();

    std::vector<R3BTrackingDetector*>& GetArray() { return fDetectors; }

    R3BHit* GetHit(const std::string& detName, const Int_t& hitId) { return GetByName(detName)->hits[hitId]; }

    Double_t GetAfterGladResolution();

  private:
    std::vector<R3BTrackingDetector*> fDetectors;
    std::map<std::string, int> fMapIndex;

    ClassDef(R3BTrackingSetup, 1)
};

#endif
