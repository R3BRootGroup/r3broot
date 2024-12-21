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

#ifndef R3BRPC_H
#define R3BRPC_H 1

#include "R3BDetector.h"
#include "Rtypes.h"
#include "TLorentzVector.h"

class TClonesArray;
class R3BRpcPoint;
class FairVolume;

class R3BRpc : public R3BDetector
{
  public:
    /** Default constructor */
    R3BRpc();

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param trans   position
     *@param rot     rotation
     */
    R3BRpc(const TString& geoFile, const TGeoTranslation& trans, const TGeoRotation& rot = TGeoRotation());

    /** Standard constructor.
     *@param geoFile name of the ROOT geometry file
     *@param combi   position + rotation
     */
    R3BRpc(const TString& geoFile, const TGeoCombiTrans& combi = TGeoCombiTrans());

    /** Destructor **/
    virtual ~R3BRpc();

    /** Virtual method ProcessHits
     **
     ** Defines the action to be taken when a step is inside the active
     ** volume. Creates a R3BRpcPoint and adds it to the collection.
     *@param vol  Pointer to the active volume
     **/
    virtual Bool_t ProcessHits(FairVolume* vol = 0);

    /** Virtual method EndOfEvent
     **
     ** If verbosity level is set, print hit collection at the
     ** end of the event and resets it afterwards.
     **/
    virtual void EndOfEvent();

    /** Virtual method Register
     **
     ** Registers the hit collection in the ROOT manager.
     **/
    virtual void Register();

    /** Accessor to the hit collection **/
    virtual TClonesArray* GetCollection(Int_t iColl) const;

    /** Virtual method Print
     **
     ** Screen output of hit collection.
     **/
    virtual void Print(Option_t* option = "") const;

    /** Virtual method Reset
     **
     ** Clears the hit collection
     **/
    virtual void Reset();

    /** Virtual method CopyClones
     **
     ** Copies the hit collection with a given track index offset
     *@param cl1     Origin
     *@param cl2     Target
     *@param offset  Index offset
     **/
    virtual void CopyClones(TClonesArray* cl1, TClonesArray* cl2, Int_t offset);

    virtual Bool_t CheckIfSensitive(std::string name);

    virtual void Initialize();
    virtual void SetSpecialPhysicsCuts() {}

  private:
    /** Track information to be stored until the track leaves the
    active volume. **/
    Int_t fTrackID;        //!  track index
    Int_t fTrackPID;       //!  particle identification
    Int_t fVolumeID;       //!  volume id
    Int_t fStripID;        //!  strip id
    TLorentzVector fPosIn; //!  position
    TLorentzVector fMomIn; //!  momentum
    Double32_t fTime;      //!  time
    Double32_t fLength;    //!  length
    Double32_t fELoss;     //!  energy loss
    Int_t fPosIndex;       //!

    TClonesArray* fRpcCollection; //!  The point collection

    /** Private method AddPoint
     **
     ** Adds a RpcPoint to the HitCollection
     **/
    R3BRpcPoint* AddPoint(Int_t trackID,
                          Int_t detID,
                          Int_t ident,
                          TVector3 posIn,
                          TVector3 momIn,
                          Double_t time,
                          Double_t length,
                          Double_t eLoss,
                          UInt_t EventId);

    /** Private method ResetParameters
     **
     ** Resets the private members for the track parameters
     **/
    void ResetParameters();

  public:
    ClassDef(R3BRpc, 1);
};

inline void R3BRpc::ResetParameters()
{
    fTrackID = fVolumeID = fTrackPID = fStripID = 0;
    fPosIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fMomIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fTime = fLength = fELoss = 0;
    fPosIndex = 0;
};

#endif /* R3BRPC_H */
