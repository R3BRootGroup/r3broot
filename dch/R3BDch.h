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

#ifndef R3BDCH_H
#define R3BDCH_H

#include "R3BDetector.h"
#include "TLorentzVector.h"
#include <string>

class TClonesArray;
class R3BDchFullPoint;
class R3BDchPoint;
class FairVolume;

class R3BDch : public R3BDetector
{

  public:
    /** Default constructor **/
    R3BDch(const TString& geoFile = "");

    /** Standard constructor.
     *@param name    detector name
     *@param geoFile name of the geometry version
     *@param active  sensitivity flag
     *@param x,y,z   position of station 1 and 2
     *@param rot_y,rot_z  rotation of station 1 and 2
     **/
    R3BDch(const TString& geoFile,
           const Double_t x1,
           const Double_t y1,
           const Double_t z1,
           const Double_t rot_y1,
           const Double_t rot_z1,
           const Double_t x2,
           const Double_t y2,
           const Double_t z2,
           const Double_t rot_y2,
           const Double_t rot_z2);

    /** Destructor **/
    ~R3BDch();

    /** Virtual method ProcessHits
     **
     ** Defines the action to be taken when a step is inside the
     ** active volume. Creates a R3BDchPoint and adds it to the
     ** collection.
     *@param vol  Pointer to the active volume
     **/
    virtual Bool_t ProcessHits(FairVolume* vol = 0);

    /** Virtual method BeginEvent
     **
     ** If verbosity level is set, print hit collection at the
     ** end of the event and resets it afterwards.
     **/

    virtual void BeginEvent();

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

    /** Virtaul method Construct geometry
     **
     ** Constructs the STS geometry
     **/
    virtual void ConstructGeometry();

    virtual Bool_t CheckIfSensitive(std::string name);

    virtual void Initialize();
    virtual void SetSpecialPhysicsCuts();
    void SetHeliumBag(Bool_t bag) { kHelium = bag; }
    void SetDynamicStepSize(Bool_t step) { fDynamicStepSize = step; }

    // renormalization of MC step size
    Double_t BetheBloch(Double_t bg);
    void SetStepToNextCollision();
    void SetVerbosity(Bool_t verb) { fVerbose = verb; }

    void PrintInfo();
    void RecordFullMcHit();
    void RecordPartialMcHit();
    void FindNodePath(TObjArray* arr);

  private:
    /** Track information to be stored until the track leaves the
    active volume. **/
    Int_t fTrackID;                 //!  track index
    Int_t fVolumeID;                //!  volume id
    TLorentzVector fPosIn, fPosOut; //!  position
    TLorentzVector fMomIn, fMomOut; //!  momentum
    Double32_t fTime_in;            //!  time when entering active volume
    Double32_t fTime_out;           //!  time when exiting active volume
    Double32_t fTime;               //!  time
    Double32_t fLength_in;          //!  length when entering active volume
    Double32_t fLength_out;         //!  length when exiting active volume
    Double32_t fLength;             //!  length
    Double32_t fELoss;              //!  energy loss
    Int_t fPosIndex;                //!
    TClonesArray* fDchCollection;   //!  The hit collection
    Bool_t kGeoSaved;               //!
    Bool_t kHelium;                 //!  Helium bag part
    TList* flGeoPar;                //!
    Bool_t fVerbose;                //!
    Bool_t fDynamicStepSize;        //!
    // Sensitive Ref Node
    TGeoMatrix* refMatrix; //! Trans reference Node

    TVector3 fPos1;      //!
    TVector3 fPos2;      //!
    TGeoRotation* fRot1; //!
    TGeoRotation* fRot2; //!

    /** Private method AddHit
     **
     ** Adds a DchPoint to the HitCollection
     **/
    R3BDchFullPoint* AddFullHit(Int_t trackId,
                                Int_t mod,
                                Int_t layer,
                                Int_t cell,
                                TVector3 pos,
                                TVector3 lpos,
                                TVector3 mom,
                                TVector3 lmom,
                                Double_t time,
                                Double_t length,
                                Double_t eLoss);

    R3BDchPoint* AddHit(Int_t trackId,
                        Int_t mod,
                        Int_t layer,
                        Int_t cell,
                        TVector3 posIn,
                        TVector3 pos_out,
                        TVector3 momIn,
                        TVector3 momOut,
                        TVector3 lpos1,
                        TVector3 lmom1,
                        TVector3 lpos2,
                        TVector3 lmom2,
                        Double_t time,
                        Double_t length,
                        Double_t eLoss);

    /** Private method ResetParameters
     **
     ** Resets the private members for the track parameters
     **/
    void ResetParameters();

    ClassDef(R3BDch, 3);
};

inline void R3BDch::ResetParameters()
{
    fTrackID = fVolumeID = 0;
    fPosIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fPosOut.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fMomIn.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fMomOut.SetXYZM(0.0, 0.0, 0.0, 0.0);
    fTime = fLength = fELoss = 0;
    fPosIndex = 0;
};

#endif
