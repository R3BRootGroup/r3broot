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

#ifndef R3BMusic_H
#define R3BMusic_H 1

#include "R3BDetector.h"
#include "TLorentzVector.h"
#include <map>

class TClonesArray;
class R3BMusicPoint;
class FairVolume;
class TGeoRotation;

class R3BMusic : public R3BDetector
{
  public:
    /** Default constructor **/
    R3BMusic();

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param trans   position
     *@param rot     rotation
     */
    R3BMusic(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot = TGeoRotation());

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param combi   position + rotation
     */
    R3BMusic(const TString& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    /** Destructor **/
    ~R3BMusic() override;

    /** Method ProcessHits **/
    Bool_t ProcessHits(FairVolume* vol = 0) override;

    /** Virtual method EndOfEvent **/
    void EndOfEvent() override;

    /** Method Register **/
    void Register() override;

    /** Accessor to the hit collection **/
    TClonesArray* GetCollection(Int_t iColl) const override;

    /** Method Print **/
    void Print(Option_t* option = "") const override;

    /** Method Reset **/
    void Reset() override;

    Bool_t CheckIfSensitive(std::string name) override;

    void Initialize() override;

  private:
    /** Track information to be stored until the track leaves the
        active volume. **/
    Int_t fTrackID;                 //!  track index
    Int_t fTrackPID;                //!  particle identification
    Int_t fVolumeID;                //!  volume id
    Int_t fDetCopyID;               //!  Det volume id
    Int_t fParentTrackID;           //!  parent track index
    TLorentzVector fPosIn, fPosOut; //!  position
    TLorentzVector fMomIn, fMomOut; //!  momentum
    Double32_t fTime;               //!  time
    Double32_t fLength;             //!  length
    Double32_t fELoss;              //!  energy loss
    Int_t fPosIndex;

    TClonesArray* fMusicCollection; //!  The point collection

    /** Private method AddPoint  **/
    R3BMusicPoint* AddPoint(Int_t trackID,
                            Int_t detID,
                            Int_t volid,
                            Double_t Z,
                            Double_t A,
                            TVector3 posIn,
                            TVector3 pos_out,
                            TVector3 momIn,
                            TVector3 momOut,
                            Double_t time,
                            Double_t length,
                            Double_t eLoss);

    /** Private method ResetParameters **/
    void ResetParameters();

  public:
    ClassDefOverride(R3BMusic, 1);
};

#endif
