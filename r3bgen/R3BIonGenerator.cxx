/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// -------------------------------------------------------------------------
// -----            Based on FairIonGenerator source file              -----
// -----            Created 18/03/11  by M. Labiche                    -----
// -------------------------------------------------------------------------
#include "R3BIonGenerator.h"

#include "FairLogger.h"
#include "FairPrimaryGenerator.h"

#include "FairIon.h"
#include "FairRunSim.h"

#include "TDatabasePDG.h"
#include "TObjArray.h"
#include "TParticlePDG.h"

#include "TMath.h"
#include "TRandom.h"

#include "FairIon.h"
#include "FairRunSim.h"
#include "TParticle.h"
#include <iostream>
using std::cout;
using std::endl;

// -----   Initialsisation of static variables   --------------------------
Int_t R3BIonGenerator::fgNIon = 0;
// ------------------------------------------------------------------------

// -----   Default constructor   ------------------------------------------
R3BIonGenerator::R3BIonGenerator()
    : fMult(0)
    , fPx(0.)
    , fPy(0.)
    , fPz(0.)
    , fR(0.)
    , fz(0.)
    , fOffsetX(0.)
    , fOffsetY(0.)
    , fSigmaP(0.)
    , fAngleX(0.)
    , fAngleY(0.)
    , fVx(0.)
    , fVy(0.)
    , fVz(0.)
    , fIon(NULL)
    , fQ(0)
    , fBeamSpotIsSet(kFALSE)
{
    //  cout << "-W- R3BIonGenerator: "
    //      << " Please do not use the default constructor! " << endl;
}
// ------------------------------------------------------------------------

R3BIonGenerator::R3BIonGenerator(const Char_t* ionName, Int_t mult, Double_t px, Double_t py, Double_t pz)
    : fMult(0)
    , fPx(0.)
    , fPy(0.)
    , fPz(0.)
    , fR(0.)
    , fz(0.)
    , fOffsetX(0.)
    , fOffsetY(0.)
    , fSigmaP(0.)
    , fAngleX(0.)
    , fAngleY(0.)
    , fVx(0.)
    , fVy(0.)
    , fVz(0.)
    , fIon(NULL)
    , fQ(0)
    , fBeamSpotIsSet(kFALSE)
{

    FairRunSim* fRun = FairRunSim::Instance();
    TObjArray* UserIons = fRun->GetUserDefIons();
    TObjArray* UserParticles = fRun->GetUserDefParticles();
    FairParticle* part = 0;
    fIon = (FairIon*)UserIons->FindObject(ionName);
    if (fIon)
    {
        fgNIon++;
        fMult = mult;
        fPx = Double_t(fIon->GetA()) * px;
        fPy = Double_t(fIon->GetA()) * py;
        fPz = Double_t(fIon->GetA()) * pz;
        cout << "Testing ion mass: " << fIon->GetA() << endl;
        // fVx   = vx;
        // fVy   = vy;
        // fVz   = vz;
        //}
    }
    else
    {
        part = (FairParticle*)UserParticles->FindObject(ionName);
        if (part)
        {
            fgNIon++;
            TParticle* particle = part->GetParticle();
            fMult = mult;
            fPx = Double_t(particle->GetMass() / 0.92827231) * px;
            fPy = Double_t(particle->GetMass() / 0.92827231) * py;
            fPz = Double_t(particle->GetMass() / 0.92827231) * pz;
            // fVx   = vx;
            // fVy   = vy;
            // fVz   = vz;
        }
    }
    if (fIon == 0 && part == 0)
    {
        cout << "-E- R3BIonGenerator: Ion or Particle is not defined !" << endl;
        LOG(fatal) << "R3BIonGenerator: No FairRun instantised!";
    }
}
// ------------------------------------------------------------------------

// -----   Default constructor   ------------------------------------------
R3BIonGenerator::R3BIonGenerator(Int_t z, Int_t a, Int_t q, Int_t mult, Double_t px, Double_t py, Double_t pz)
    : fMult(0)
    , fPx(0.)
    , fPy(0.)
    , fPz(0.)
    , fR(0.)
    , fz(0.)
    , fOffsetX(0.)
    , fOffsetY(0.)
    , fSigmaP(0.)
    , fAngleX(0.)
    , fAngleY(0.)
    , fVx(0.)
    , fVy(0.)
    , fVz(0.)
    , fIon(NULL)
    , fQ(0)
    , fBeamSpotIsSet(kFALSE)
{
    fgNIon++;
    fMult = mult;
    fPx = Double_t(a) * px;
    fPy = Double_t(a) * py;
    fPz = Double_t(a) * pz;
    // fVx   = vx;
    // fVy   = vy;
    // fVz   = vz;
    char buffer[20];
    sprintf(buffer, "FairIon%d", fgNIon);
    fIon = new FairIon(buffer, z, a, q);

    cout << "New Ion with Z = " << z << " A = " << a << " q = " << q << endl;

    FairRunSim* run = FairRunSim::Instance();
    if (!run)
    {
        cout << "-E- FairIonGenerator: No FairRun instantised!" << endl;
        LOG(fatal) << "FairIonGenerator: No FairRun instantised!";
    }
    run->AddNewIon(fIon);
}
//_________________________________________________________________________

R3BIonGenerator::R3BIonGenerator(const R3BIonGenerator& right)
    : fMult(right.fMult)
    , fPx(right.fPx)
    , fPy(right.fPy)
    , fPz(right.fPz)
    , fR(right.fR)
    , fz(right.fz)
    , fOffsetX(right.fOffsetX)
    , fOffsetY(right.fOffsetY)
    , fSigmaP(right.fSigmaP)
    , fAngleX(right.fAngleX)
    , fAngleY(right.fAngleY)
    , fVx(right.fVx)
    , fVy(right.fVy)
    , fVz(right.fVz)
    , fIon(right.fIon)
    , fQ(right.fQ)
    , fBeamSpotIsSet(right.fBeamSpotIsSet)
{
}

// -----   Destructor   ---------------------------------------------------
R3BIonGenerator::~R3BIonGenerator()
{
    // if (fIon) delete fIon;
}
//_________________________________________________________________________

// -----   Public method SetExcitationEnergy   ----------------------------
void R3BIonGenerator::SetExcitationEnergy(Double_t eExc) { fIon->SetExcEnergy(eExc); }
//_________________________________________________________________________

// -----   Public method SetMass   ----------------------------------------
void R3BIonGenerator::SetMass(Double_t mass) { fIon->SetMass(mass); }
//_________________________________________________________________________

// -----   Public method ReadEvent   --------------------------------------
Bool_t R3BIonGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{

    Double_t Phi, SpotR;

    // if ( ! fIon ) {
    //   cout << "-W- FairIonGenerator: No ion defined! " << endl;
    //   return kFALSE;
    // }

    TParticlePDG* thisPart = TDatabasePDG::Instance()->GetParticle(fIon->GetName());
    if (!thisPart)
    {
        cout << "-W- FairIonGenerator: Ion " << fIon->GetName() << " not found in database!" << endl;
        return kFALSE;
    }

    int pdgType = thisPart->PdgCode();

    // cout << "fR=" << fR << " fz=" << fz <<endl;

    Phi = gRandom->Uniform(0, 360) * TMath::DegToRad();
    SpotR = sqrt(gRandom->Uniform(0, fR * fR));

    if (fBeamSpotIsSet)
    {
        fVx = fOffsetX + SpotR * cos(Phi); // gRandom->Uniform(-fx,fx);
        fVy = fOffsetY + SpotR * sin(Phi); // gRandom->Uniform(-fy,fy);
        fVz = fz;

		//cout << "x, y :" << fVx << "  "  << fVy << endl;

		// Special for experiment s494. we make a focus on the hole of fib23

        Double_t p = sqrt(fPx * fPx + fPy * fPy + fPz * fPz);
        fPx = -atan(fVx / 245.05) * fPz;
        fPy = -atan(fVy / 245.05) * fPz;
        // cout << "Px, Py :" << fPx << "  "  << fPy << endl;
        fPz = sqrt(p * p - fPx * fPx - fPy * fPy);
    }
    else
    {
        fVx = 0.0;
        fVy = 0.0;
        fVz = 0.0;
    }

    if (fSigmaP > 0.)
    {
        Double_t p0 = gRandom->Uniform(fPz - fSigmaP, fPz + fSigmaP);
        fPz = p0;
    }

    Double_t thetaY, thetaX;

    if (fAngleX > 0. && fAngleY == 0.)
    {
        Double_t p = sqrt(fPx * fPx + fPy * fPy + fPz * fPz);

        thetaX = gRandom->Uniform(-fAngleX, fAngleX); // max angle in rad
        /*if(thetaX < 0)
        {
            thetaX=-fAngleX;
        }
        else
        {
            thetaX=fAngleX;
        }*/

        // Double_t thetaX = fAngle;
        fPx = sin(thetaX) * p;

        fPz = sqrt(p * p - fPx * fPx - fPy * fPy);

        // cout<< " theta "<< thetaX << "  " << thetaY << "  "<<endl;
        // cout << "p alt: "<< p << " p neu: " << sqrt(fPx*fPx+fPy*fPy+fPz*fPz) << endl;
    }
    if (fAngleX == 0. && fAngleY > 0.)
    {
        Double_t p = sqrt(fPx * fPx + fPy * fPy + fPz * fPz);

        thetaY = gRandom->Uniform(-fAngleY, fAngleY); // max angle in rad
        // Double_t thetaY = 0.;
        // fPy = tan(thetaY) * fPz;
        /*if(thetaY < 0)
        {
            thetaY=-fAngleY;
        }
        else
        {
            thetaY=fAngleY;
        }*/
        fPy = sin(thetaY) * p;

        // Double_t thetaY = 0.04;
        // fPy = tan(thetaY) * fPz;

        fPz = sqrt(p * p - fPx * fPx - fPy * fPy);

        // cout<< " theta "<< thetaX << "  " << thetaY << "  "<<endl;
        // cout << "p alt: "<< p << " p neu: " << sqrt(fPx*fPx+fPy*fPy+fPz*fPz) << endl;
    }
    if (fAngleX > 0. && fAngleY > 0.)
    {
        Double_t p = sqrt(fPx * fPx + fPy * fPy + fPz * fPz);

        thetaY = gRandom->Uniform(-fAngleY, fAngleY); // max angle in rad
        // Double_t thetaY = 0.;
        /*if(thetaY < 0)
        {
            thetaY=-fAngleY;
        }
        else
        {
            thetaY=fAngleY;
        }*/
        fPy = sin(thetaY) * p;

        thetaX = gRandom->Uniform(-fAngleX, fAngleX); // max angle in rad
        // Double_t thetaX = fAngle;

        fPx = cos(thetaY) * sin(thetaX) * p;

        fPz = sqrt(p * p - fPx * fPx - fPy * fPy);

        // cout<< " theta "<< thetaX << "  " << thetaY << "  "<<endl;
        // cout << "p alt: "<< p << " p neu: " << sqrt(fPx*fPx+fPy*fPy+fPz*fPz) << endl;
    }

//      cout << "-I- FairIonGenerator: Generating " << fMult << " ions of type "
//           << fIon->GetName() << " (PDG code " << pdgType << ")" << endl;
//      cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz
//           << ") Gev from vertex (" << fVx << ", " << fVy
//           << ", " << fVz << ") cm" << endl;

    for (Int_t i = 0; i < fMult; i++)
        primGen->AddTrack(pdgType, fPx, fPy, fPz, fVx, fVy, fVz);

    return kTRUE;

    return kTRUE;
}

//_____________________________________________________________________________

ClassImp(R3BIonGenerator)
