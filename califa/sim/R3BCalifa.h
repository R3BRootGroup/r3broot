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

#ifndef R3BCALIFA_H
#define R3BCALIFA_H 1

#include <R3BDetector.h>
#include <Rtypes.h>
#include <TLorentzVector.h>
#include <map>
#include <string>

class TClonesArray;
class R3BCalifaPoint;
class FairVolume;
class TGeoRotation;
class R3BCalifaGeometry;

class R3BCalifa : public R3BDetector
{
  public:
    /** Default constructor */
    R3BCalifa();

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param trans   position
     *@param rot     rotation
     */
    explicit R3BCalifa(const std::string& geoFile,
                       const TGeoTranslation& trans,
                       const TGeoRotation& rot = TGeoRotation());

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param combi   position + rotation
     */
    explicit R3BCalifa(const std::string& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    /** Destructor **/
    ~R3BCalifa();

    /** Virtual method ProcessHits
     **
     ** Defines the action to be taken when a step is inside the active
     ** volume. Creates a R3BCalifaPoint and adds it to the collection.
     *@param vol  Pointer to the active volume
     **/
    Bool_t ProcessHits(FairVolume* vol = 0) override;

    /** Virtual method EndOfEvent
     **
     ** If verbosity level is set, print hit collection at the
     ** end of the event and resets it afterwards.
     **/
    void EndOfEvent() override;

    /** Virtual method Register
     **
     ** Registers the hit collection in the ROOT manager.
     **/
    void Register() override;

    /** Accessor to the hit collection **/
    TClonesArray* GetCollection(Int_t iColl) const override;

    /** Virtual method Print
     **
     ** Screen output of hit collection.
     **/
    void Print(Option_t* option = "") const override;

    /** Virtual method Reset
     **
     ** Clears the hit collection
     **/
    void Reset() override;

    /** Virtual method CopyClones
     **
     ** Copies the hit collection with a given track index offset
     *@param cl1     Origin
     *@param cl2     Target
     *@param offset  Index offset
     **/
    void CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset) override;

    Bool_t CheckIfSensitive(std::string name) override;

    /** Select the version of the CALIFA geometry
     **
     *@param version Geometry version (final BARREL+iPhos: 2020)
     **/
    void SelectGeometryVersion(Int_t version);

    void Initialize() override;

  private:
    /** Track information to be stored until the track leaves the
    active volume. **/
    int fTrackID = 0;  //!  track index
    int fVolumeID = 0; //!  volume id
    int fTrackPID = 0;
    TLorentzVector fPosIn;                 //!  position
    TLorentzVector fMomIn;                 //!  momentum
    double fTime = 0.;                     //!  time
    double fLength = 0., fLengthzero = 0.; //!  length
    double fELoss = 0.;                    //!  energy loss
    int fPosIndex = 0;                     //!
    int fNSteps = 0;                       //!  Number of steps in the active volume

    TClonesArray* fCalifaCollection; //!  The point collection

    // Selecting the geometry of the CALIFA calorimeter (final BARREL+iPhos: 2021)
    int fGeometryVersion = 2021;

    /** Private method AddPoint
     **
     ** Adds a CalifaPoint to the HitCollection
     **/
    R3BCalifaPoint* AddPoint(int trackID,
                             int detID,
                             int trackPID,
                             int ident,
                             TVector3 posIn,
                             TVector3 momIn,
                             double time,
                             double length,
                             double eLoss);

    /** Private method ResetParameters
     **
     ** Resets the private members for the track parameters
     **/
    void ResetParameters();

  public:
    ClassDefOverride(R3BCalifa, 9);
};

inline void R3BCalifa::ResetParameters()
{
    fTrackID = fVolumeID = fTrackPID = 0;
    fPosIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fMomIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fTime = fLength = fLengthzero = fELoss = 0.0;
    fPosIndex = 0;
    fNSteps = 0;
};

#endif /* R3BCALIFA_H */
