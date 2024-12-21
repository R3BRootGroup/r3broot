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

#pragma once

#ifndef R3BCALIFAGEOMETRY_H
#define R3BCALIFAGEOMETRY_H 1

#include <Rtypes.h>
#include <TFile.h>
#include <TObject.h>
#include <TVector3.h>
#include <iostream>
#include <sstream>

class TGeoNavigator;

/**
 * Geometrical queries to CALIFA
 */
class R3BCalifaGeometry : public TObject
{
  public:
    /** Default constructor */
    R3BCalifaGeometry();

    /** Destructor **/
    ~R3BCalifaGeometry();

    // Initialization with a Root version
    bool Init(int fGeo);

    // Setter and getter for a CALIFA reference point, by default (0,0,0)
    void SetReferencePoint(const TVector3 refpoint) { fRefPoint = refpoint; }
    const TVector3& SetReferencePoint() { return fRefPoint; }

    /**
     * Gets position in polar coordinates of crystal with given ID.
     * On error, the x,y and z component of the TVector3 are set to NAN.
     * @param iD crystal ID (depending on geometry version)
     */
    const TVector3& GetAngles(int iD);

    /**
     * Legacy: Gets position in polar coordinates of crystal with given ID.
     *
     * Before this comment, the last three parameters were untouched on error.
     * Now, they should be set to NAN;
     *
     * @param iD crystal ID (depending on geometry version)
     * @param polar [out] Will be filled with polar angle (radians) of crystal center
     * @param azimuthal [out] Will be filled with azimuthal angle (radians) of crystal center
     * @param rho [out] Will be filled with distance (cm) of crystal center to target position (0,0,0)
     */
    void GetAngles(int iD, double* polar, double* azimuthal, double* rho);

    /**
     * Gets volume path of crystal with given ID.
     *
     * @param iD crystal ID (depending on geometry version)
     * @return Volume path
     */
    std::string GetCrystalVolumePath(int iD);

    /**
     * Gets crystal ID for given volume path.
     *
     * @param volumePath Volume path
     * @return Crystal ID
     */
    int GetCrystalId(const std::string volumePath);

    /**
     * Calculate the distance of a given straight track through the active detector volume (crystal(s)). Usefull for
     * iPhos.
     *
     * @param startVertex Start position of track (e.g. (0,0,0) for center of target)
     * @param direction Direction of the track (does not need to be normalized)
     * @param hitPos [out] (optional) If given, this TVector3 will be filled with the position of the first contact
     * between the track and the first crystal along the path.
     * @param numCrystals [in,out]
     * \parblock
     *  (optional) If given together with crystalIDs,
     *  this function will fill numCrystals with the number of crystals along the straight
     *  and the array crystalIds with the IDs of the crystals along the straight.
     *  The array crystalIds MUST be allocated by the caller.
     *  The maximum size of the array has to be given with the content of numCrystals.
     *  The function will stop filling the array after reaching the specified, maximum capacity.
     * \endparblock
     * @param crystalIds [out] (optional, array) See numCrystals.
     * @return Distance (cm) through crystal(s)
     */
    double GetDistanceThroughCrystals(TVector3& startVertex,
                                      TVector3& direction,
                                      TVector3* hitPos = NULL,
                                      int* numCrystals = NULL,
                                      int* crystalIds = NULL);

    /**
     * @return if we are running the simulation or data analysis
     */
    const bool IsSimulation() { return fIsSimulation; }

    const int GetNbCrystals() { return fNumCrystals; }

    const int GetGeoVersion() { return fGeometryVersion; }

    /**
     * Returns singleton instance of R3BCalifaGeometry for given geometry version.
     * Warning: Switching the geometry version at run time will lead to undefined beheaviour!
     *
     * @param version Geometry version to use. If in doubt, use 2020.
     * @return Instance of R3BCalifaGeometry
     */
    static R3BCalifaGeometry* Instance();

  private:
    int fGeometryVersion = 2024;
    int fNumCrystals = 5088;
    bool fIsSimulation = false;
    bool IsInitialize = false;
    TVector3 fRefPoint;
    TFile* f;

  public:
    ClassDefOverride(R3BCalifaGeometry, 9);
};

#endif /* R3BCALIFAGEOMETRY_H */
