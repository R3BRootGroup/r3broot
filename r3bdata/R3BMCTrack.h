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

#ifndef R3BMCTRACK_H
#define R3BMCTRACK_H

#include "Math/Vector4D.h"
#include "R3BDetectorList.h"
#include "TLorentzVector.h"
#include "TObject.h"
#include "TParticle.h"
#include "TVector3.h"
#include <array>

class R3BMCTrack : public TObject
{
  public:
    /**  Default constructor  **/
    R3BMCTrack();

    /**  Standard constructor  **/
    R3BMCTrack(int pdgCode,
               int motherID,
               ROOT::Math::XYZTVector xyzt,
               ROOT::Math::PxPyPzMVector pm,
               std::array<int, kLAST + 1> nPoints);

    /**  Constructor from TParticle  **/
    R3BMCTrack(TParticle* particle, std::array<int, kLAST + 1> nPoints, int fMC);

    /**  Destructor  **/
    ~R3BMCTrack() override = default;

    /**  Output to screen  **/
    void Print(Option_t* option = "") const override;

    /**  Accessors  **/
    int GetPdgCode() const { return fPdgCode; }
    int GetMotherId() const { return fMotherId; }

    double GetStartX() const { return fStartVertex.X(); }
    double GetStartY() const { return fStartVertex.Y(); }
    double GetStartZ() const { return fStartVertex.Z(); }
    double GetStartT() const { return fStartVertex.T(); }
    ROOT::Math::XYZTVector GetStartVertex() const { return fStartVertex; }

    double GetPx() const { return fMomentumMass.Px(); }
    double GetPy() const { return fMomentumMass.Py(); }
    double GetPz() const { return fMomentumMass.Pz(); }
    double GetMass() const { return fMomentumMass.M(); }
    ROOT::Math::PxPyPzMVector GetMomentumMass() const { return fMomentumMass; };
    ROOT::Math::PxPyPzEVector GetFourMomentum() const { return ROOT::Math::PxPyPzEVector(fMomentumMass); };
    double GetEnergy() const { return fMomentumMass.E(); }
    double GetPt() const { return fMomentumMass.Pt(); }
    double GetP() const { return fMomentumMass.P(); }
    double GetRapidity() const { return fMomentumMass.Rapidity(); };

    /** Accessors to the number of MCPoints in the detectors **/
    int GetNPoints(DetectorId detId) const;

    /**  Modifiers  **/
    void SetMotherId(int id) { fMotherId = id; } // Used by R3BMCStack UpdateTrackIndex. Needed? -> Could make const!

  private:
    /**  PDG particle code  **/
    const int fPdgCode;

    /**  Index of mother track. -1 for primary particles.  **/
    int fMotherId;

    /** Coordinates of start vertex [cm, ns]  **/
    const ROOT::Math::XYZTVector fStartVertex;

    /** Momentum and mass at start vertex [GeV]  **/
    const ROOT::Math::PxPyPzMVector fMomentumMass;

    /**  Array representing the number of MCPoints for this track in each subdetector. **/
    const std::array<int, kLAST + 1> fNPoints;

    ClassDefOverride(R3BMCTrack, 3);
};

#endif
