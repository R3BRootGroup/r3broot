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

// -----------------------------------------------------------------------------
// -----                      R3BNeutronTracker2D                          -----
// -----                Created 29-02-2012 by D.Kresan                     -----
// -----            Based on algorithm developed by M.Heil                 -----
// -----------------------------------------------------------------------------

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>

#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BLandDigi.h"
#include "R3BLandDigiPar.h"
#include "R3BLandFirstHits.h"
#include "R3BLandPoint.h"
#include "R3BMCTrack.h"
#include "R3BNeutHit.h"
#include "R3BNeutronTracker2D.h"

#include <assert.h>

using std::cout;
using std::endl;
using std::ifstream;

ClassImp(R3BNeutronTracker2D)

    Double_t c = 29.9792458;
Double_t gBeamBeta;
bool AuxSortClustersBeta(R3BNeuLandCluster*, R3BNeuLandCluster*);

// -----------------------------------------------------------------------------
R3BNeutronTracker2D::R3BNeutronTracker2D()
    : FairTask("R3B NeuLAND Neutron Tracker")
    , f2DCutEnabled(kTRUE)
    , fNNeutrons(0)
{
    dio = 10.6; // 3 times half the diogonal of a paddle
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BNeutronTracker2D::~R3BNeutronTracker2D() {}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::SetParContainers()
{
    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fLandDigiPar = dynamic_cast<R3BLandDigiPar*>(rtdb->getContainer("R3BLandDigiPar"));

    if (fLandDigiPar)
    {
        cout << "-I- R3BLandDigitizer::SetParContainers() " << endl;
        cout << "-I- Container R3BLandDigiPar  loaded " << endl;
    }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
InitStatus R3BNeutronTracker2D::Init()
{
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fLandPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("LandPoint"));
    fLandMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    fLandDigi = dynamic_cast<TClonesArray*>(ioman->GetObject("LandDigi"));
    fLandFirstHits = dynamic_cast<TClonesArray*>(ioman->GetObject("LandFirstHits"));
    fArrayCluster = dynamic_cast<TClonesArray*>(ioman->GetObject("NeuLandCluster"));

    // New structure created by the Neutron Tracker
    fNeutHits = new TClonesArray("R3BNeutHit");
    ioman->Register("LandNeutHits", "Neutron Hits", fNeutHits, kTRUE);

    npaddles = fLandDigiPar->GetMaxPaddle() + 1;
    nplanes = fLandDigiPar->GetMaxPlane();
    cout << "# paddles: " << npaddles - 1 << "  # planes: " << nplanes << endl;
    amu = 931.494028;
    mNeutron = 1.0086649156 * amu;
    cMedia = 14.; // speed of light in material in [cm/ns]
    eventNo = 0;
    printing = 0;

    PRIM_part = new R3BPrimPart*[10];
    PRIM_prot = new R3BPrimPart*[10];
    PRIM_frag = new R3BPrimPart*[10];
    PRIM_gamma = new R3BPrimPart*[10];

    CreateHistograms();

    return kSUCCESS;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::ReadCalibrFile(char* name)
{
    // Open calibration file
    ifstream cFile;
    cFile.open(name);
    // Check if opened succesfully
    if (!cFile.is_open())
    {
        LOG(fatal) << "ReadCalibrFile: Calibration file does not exist! Aborting...";
    }

    // Read parameters
    cFile >> fKappa;
    Int_t i;
    for (i = 0; i < 5; i++)
    {
        cFile >> fCuts[i];
        if (cFile.eof())
        {
            break;
        }
    }
    for (Int_t j = i + 1; j < 5; j++)
    {
        fCuts[j] = 0.;
    }

    // Close file
    cFile.close();

    // Control output
    cout << "-I- R3BNeutronTracker2D : calibration : " << endl;
    cout << "                                kappa : " << fKappa << endl;
    for (i = 0; i < 5; i++)
    {
        cout << "                                      : " << fCuts[i] << endl;
    }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::Exec(Option_t* opt)
{
    if (0 == (eventNo % 100))
    {
        cout << "-I- R3BNeutronTracker2D : Event #: " << eventNo << "." << endl;
    }
    eventNo += 1;

    Reset();

    R3BMCTrack* aTrack1;
    Int_t prim;
    Int_t particleID;
    for (Int_t it = 0; it < fLandMCTrack->GetEntriesFast(); it++)
    {
        aTrack1 = dynamic_cast<R3BMCTrack*>(fLandMCTrack->At(it));

        prim = aTrack1->GetMotherId();
        if (-1 != prim)
        {
            break;
        }

        particleID = aTrack1->GetPdgCode();

        if (2112 == particleID)
        {
            // neutron
            PRIM_part[fNPrimNeut] = new R3BPrimPart(aTrack1->GetPdgCode(),
                                                    aTrack1->GetPx() * 1000.,
                                                    aTrack1->GetPy() * 1000.,
                                                    aTrack1->GetPz() * 1000.,
                                                    aTrack1->GetStartX(),
                                                    aTrack1->GetStartY(),
                                                    aTrack1->GetStartZ(),
                                                    aTrack1->GetStartT(),
                                                    1.0086649,
                                                    mNeutron);
            fNPrimNeut += 1;
        }
        else if (2212 == particleID)
        {
            // proton
            PRIM_prot[fNPrimProt] = new R3BPrimPart(aTrack1->GetPdgCode(),
                                                    aTrack1->GetPx() * 1000.,
                                                    aTrack1->GetPy() * 1000.,
                                                    aTrack1->GetPz() * 1000.,
                                                    aTrack1->GetStartX(),
                                                    aTrack1->GetStartY(),
                                                    aTrack1->GetStartZ(),
                                                    aTrack1->GetStartT(),
                                                    1.00782503207,
                                                    1.00782503207 * amu);
            fNPrimProt += 1;
        }
        else if (22 == particleID)
        {
            // gamma
            PRIM_gamma[fNPrimGamma] = new R3BPrimPart(aTrack1->GetPdgCode(),
                                                      aTrack1->GetPx() * 1000.,
                                                      aTrack1->GetPy() * 1000.,
                                                      aTrack1->GetPz() * 1000.,
                                                      aTrack1->GetStartX(),
                                                      aTrack1->GetStartY(),
                                                      aTrack1->GetStartZ(),
                                                      aTrack1->GetStartT(),
                                                      0.,
                                                      0.);
            fNPrimGamma += 1;
        }
        else
        {
            // fragment
            Double_t A = aTrack1->GetMass();
            PRIM_frag[fNofFrag] = new R3BPrimPart(aTrack1->GetPdgCode(),
                                                  aTrack1->GetPx() * 1000.,
                                                  aTrack1->GetPy() * 1000.,
                                                  aTrack1->GetPz() * 1000.,
                                                  aTrack1->GetStartX(),
                                                  aTrack1->GetStartY(),
                                                  aTrack1->GetStartZ(),
                                                  aTrack1->GetStartT(),
                                                  A,
                                                  A * amu);
            fNofFrag += 1;
        }
    }

    fNofClusters = fArrayCluster->GetEntries();

    // Determine amount of neutrons
    // Loop over digis
    sumTotalEnergy = 0;
    R3BLandDigi* digi;
    for (Int_t id = 0; id < fLandDigi->GetEntries(); id++)
    {
        digi = dynamic_cast<R3BLandDigi*>(fLandDigi->At(id));
        fEtot += digi->GetQdc();
        sumTotalEnergy += digi->GetQdc();
    }
    // Double_t neutmult=0.0098379928527*fEtot/beamEnergy*600.;
    Double_t neutmult = 0.00960093705145622332 * fEtot / beamEnergy * 600.;
    nNeut = (Int_t)(neutmult + 0.5);
    hNeutmult->Fill(nNeut);

    h_nofclusters->Fill(fNofClusters);
    h_etot->Fill(fEtot);
    h_ncl_etot->Fill(fEtot, fNofClusters);
    h_ndigi_etot->Fill(fEtot, nentries);

    if (f2DCutEnabled)
    {
        if (fNofClusters >= ((0. - fKappa * fCuts[4]) / (fCuts[4] - 0.) * (fEtot - 0.) + fKappa * fCuts[4]))
        {
            nNeut = 5;
        }
        else if (fNofClusters >= ((0. - fKappa * fCuts[3]) / (fCuts[3] - 0.) * (fEtot - 0.) + fKappa * fCuts[3]))
        {
            nNeut = 4;
        }
        else if (fNofClusters >= ((0. - fKappa * fCuts[2]) / (fCuts[2] - 0.) * (fEtot - 0.) + fKappa * fCuts[2]))
        {
            nNeut = 3;
        }
        else if (fNofClusters >= ((0. - fKappa * fCuts[1]) / (fCuts[1] - 0.) * (fEtot - 0.) + fKappa * fCuts[1]))
        {
            nNeut = 2;
        }
        else if (fNofClusters >= ((0. - fKappa * fCuts[0]) / (fCuts[0] - 0.) * (fEtot - 0.) + fKappa * fCuts[0]))
        {
            nNeut = 1;
        }
        else
        {
            nNeut = 0;
        }
    }
    else
    {
        nNeut = fNNeutrons;
    }

    h_ntracks->Fill(nNeut);

    if (nNeut > 0)
    {
        // Find neutrons
        nAboveThresh = 0;
        AdvancedMethod();
        // if(nAboveThresh > 3) cout<<"Neutrons above threshold: "<<nAboveThresh<<"  event: "<<eventNo<< "  Neutrons: "
        // <<nNeut<<endl;
        hClusters1->Fill(fNofClusters1);
        h_ncl1_etot->Fill(fEtot, fNofClusters1);
    }

    // Calculate invariant mass of neutrons
    CalculateMassInv();
    //   CalculateExce();

    hMinv0->Fill(fMinvTrue);
    hExce0->Fill(fExceTrue);
    if (nNeut != 0 && fNofTracks == fNPrimNeut)
    {
        hMinv->Fill(fMinv);
        hExce->Fill(fExce);
    }
    /*
      if(fMinv>0.45 && fMinv<0.55 && nNeut != 0 && fNofTracks == fNPrimNeut){
          cout<<" Good event: "<<eventNo<< "  Neutrons: " <<nNeut<<endl;
          h_nofclusters->Fill(fNofClusters);
          h_etot->Fill(fEtot);
          h_ncl_etot->Fill(fEtot, fNofClusters);
          h_ndigi_etot->Fill(fEtot, nentries);
      }
    */
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Int_t R3BNeutronTracker2D::AdvancedMethod()
{
    // Loop over clusters
    R3BNeuLandCluster* cluster;
    for (Int_t ic = 0; ic < fNofClusters; ic++)
    {
        // Get pointer to the cluster
        cluster = dynamic_cast<R3BNeuLandCluster*>(fArrayCluster->At(ic));
        // Process starting from this cluster
        NextIteration(ic, cluster);
    } // clusters

    // Sort clusters
    SortClustersBeta();

    Int_t nOutput = 0;
    Bool_t mapUsed[3000];
    for (Int_t i = 0; i < 3000; i++)
    {
        mapUsed[i] = kFALSE;
    }
    for (Int_t i = 0; i < nNeut; i++)
    {
        for (Int_t ic = 0; ic < fNofClusters1; ic++)
        {
            cluster = fVectorCluster.at(ic);

            // Skip found clusters
            if (mapUsed[ic])
            {
                continue;
            }

            // Check energy of the cluster
            if (cluster->GetE() < 2.5 && ic > 0)
            {
                continue;
            }

            // Check if in the beta-window
            TVector3 pos;
            cluster->StartPosition(pos);
            Double_t beta = pos.Z() / cluster->GetStartT() / c;
            if (TMath::Abs(beta - beamBeta) > (0.05 * 600. / beamEnergy) && ic > 0)
            {
                continue;
            }

            /*
                  // Check if cluster has more energy than half the neutron energy
                  if(cluster->GetE() > 150. ) {
                      nAboveThresh++;
                   //cout<<"Cluster energy: "<<	cluster->GetE() <<"  event: "<<eventNo<<endl;
                  }
            */
            /*
                  // Check cluster size
                  if(cluster->GetSize() < 8. ) {
                    //  nAboveThresh++;
                   //cout<<"Cluster energy: "<<	cluster->GetE() <<"  event: "<<eventNo<<endl;
                 continue;
                  }
            */
            /*
                  if(cluster->GetSize()*cluster->GetE()>2000. ) {
            //		  cout<<"Cluster energy: "<<	cluster->GetE()*cluster->GetSize() <<"  event: "<<eventNo<<endl;
                      nAboveThresh++;
                  }
            */

            // Create neutron track
            fTracks[fNofTracks][0] = ic;
            mapUsed[ic] = kTRUE;
            fNofHits[fNofTracks] += 1;
            fMapTracks[fNofTracks] = kTRUE;
            new ((*fNeutHits)[fNofTracks]) R3BNeutHit(pos.X(), pos.Y(), pos.Z(), cluster->GetStartT());
            nOutput += 1;
            fNofTracks += 1;
            break;
        } // clusters
    }     // nNeut

    return nOutput;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::SortClustersBeta()
{
    if (0 == fNofClusters)
    {
        return;
    }

    R3BNeuLandCluster* cluster;
    for (Int_t ic = 1; ic < fNofClusters; ic++)
    {
        cluster = dynamic_cast<R3BNeuLandCluster*>(fArrayCluster->At(ic));
        if (cluster->IsMarked())
        {
            continue;
        }
        fVectorCluster.push_back(cluster);
    } // Clusters

    fNofClusters1 = fVectorCluster.size();
    if (fNofClusters1 > 1)
    {
        std::sort(fVectorCluster.begin(), fVectorCluster.end(), AuxSortClustersBeta);
    }

    cluster = dynamic_cast<R3BNeuLandCluster*>(fArrayCluster->At(0));
    fVectorCluster.insert(fVectorCluster.begin(), cluster);
    fNofClusters1 = fVectorCluster.size();
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::NextIteration(Int_t curIndex, R3BNeuLandCluster* curClus)
{
    if (curClus->GetSize() < 2)
    {
        return;
    }

    // Loop over clusters
    R3BNeuLandCluster* cluster;
    for (Int_t ic = (curIndex + 1); ic < fNofClusters; ic++)
    {
        cluster = dynamic_cast<R3BNeuLandCluster*>(fArrayCluster->At(ic));
        // Check if elastic scattering
        if (IsElastic(curClus, cluster))
        {
            // Elastic, so mark the cluster
            cluster->Mark();
            // Stop the loop
            break;
        }
    } // clusters
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
Bool_t R3BNeutronTracker2D::IsElastic(R3BNeuLandCluster* c1, R3BNeuLandCluster* c2)
{
    Double_t x0 = 0.;
    Double_t y0 = 0.;
    Double_t z0 = 0.;
    Double_t t0 = 0.;

    //  Double_t v1xmin,v1ymin,v1zmin,v1xmax,v1ymax,v1zmax;
    Double_t v4xmin, v4ymin, v4zmin; //,v4xmax,v4ymax,v4zmax;
    Double_t v3xmin, v3ymin, v3zmin; //,v3xmax,v3ymax,v3zmax;
                                     //  Double_t v6xmin,v6ymin,v6zmin,v6xmax,v6ymax,v6zmax;

    //    search for scattering
    //    elastic scattering of particle 1 (neutron) on a particle 2 (proton) at rest.
    //    Outgoing particles are 3 (scattered neutron) and 4 (recoil proton).
    //    angle of particles after scattering is theta3 and theta4
    //    K is kinetic energy, p is momentum

    TVector3 pos1;
    TVector3 pos2;

    c1->StartPosition(pos1);

    //    incoming particle
    // vector from previous interaction to present interaction
    Double_t v1x = pos1.X() - x0;
    Double_t v1y = pos1.Y() - y0;
    Double_t v1z = pos1.Z() - z0;
    // time difference and distance between previous interaction to present interaction
    Double_t dt = c1->GetStartT() - t0;
    Double_t dr = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);

    Double_t beta1 = dr / dt / c;
    Double_t beta1max = (dr + dio) / dt / c;
    Double_t beta1min = (dr - dio) / dt / c;
    if (beta1max > 0.99)
        beta1max = 0.99;
    if (beta1min < 0.)
        beta1min = 0.;

    Double_t gamma1 = 1. / sqrt(1. - beta1 * beta1);
    Double_t gamma1min = 1. / sqrt(1. - beta1min * beta1min);
    Double_t gamma1max = 1. / sqrt(1. - beta1max * beta1max);
    Double_t p1 = beta1 * gamma1 * 1. * amu;
    Double_t p1min = beta1min * gamma1min * 1. * amu;
    Double_t p1max = beta1max * gamma1max * 1. * amu;
    Double_t En1 = sqrt(p1 * p1 + amu * amu);
    Double_t En1min = sqrt(p1min * p1min + amu * amu);
    Double_t En1max = sqrt(p1max * p1max + amu * amu);
    Double_t K1 = En1 - amu;
    Double_t K1min = En1min - amu;
    Double_t K1max = En1max - amu;

    // particle 4 is proton
    c1->StopPosition(pos2);
    Double_t v4x = (pos2 - pos1).X();
    Double_t v4y = (pos2 - pos1).Y();
    Double_t v4z = (pos2 - pos1).Z();
    // vector perpendicular to scattering plane
    Double_t v5x = v1y * v4z - v1z * v4y;
    Double_t v5y = v1z * v4x - v1x * v4z;
    Double_t v5z = v1x * v4y - v1y * v4x;

    Double_t tempAngle;
    dt = c1->GetStopT() - c1->GetStartT();
    dr = sqrt(v4x * v4x + v4y * v4y + v4z * v4z);

    Double_t theta4Measured = acos((v1x * v4x + v1y * v4y + v1z * v4z) / sqrt(v1x * v1x + v1y * v1y + v1z * v1z) /
                                   sqrt(v4x * v4x + v4y * v4y + v4z * v4z));
    Double_t theta4Measuredmin = theta4Measured;
    Double_t theta4Measuredmax = theta4Measured;
    if ((v4x * v4x + v4y * v4y + v4z * v4z) == 0.)
    {
        theta4Measured = 1.55;
        theta4Measuredmin = 0.;
        theta4Measuredmax = 1.55;
    }
    for (Int_t k = 1; k < 3; k++)
    {
        if (k == 1)
        {
            v4xmin = v4x - dio;
        }
        else
        {
            v4xmin = v4x + dio;
        }
        for (Int_t l = 1; l < 3; l++)
        {
            if (l == 1)
            {
                v4ymin = v4y - dio;
            }
            else
            {
                v4ymin = v4y + dio;
            }
            for (Int_t m = 1; m < 3; m++)
            {
                if (m == 1)
                {
                    v4zmin = v4z - dio;
                }
                else
                {
                    v4zmin = v4z + dio;
                }
                tempAngle =
                    acos((v1x * v4xmin + v1y * v4ymin + v1z * v4zmin) / sqrt(v1x * v1x + v1y * v1y + v1z * v1z) /
                         sqrt(v4xmin * v4xmin + v4ymin * v4ymin + v4zmin * v4zmin));
                if (tempAngle > theta4Measuredmax)
                    theta4Measuredmax = tempAngle;
                if (tempAngle < theta4Measuredmin)
                    theta4Measuredmin = tempAngle;
            }
        }
    }
    if (theta4Measuredmax > 1.55)
        theta4Measuredmax = 1.55;

    // calculate velocity of neutron after scattering
    TVector3 pos3;

    // Vector of scattered neutron
    c2->StartPosition(pos3);
    Double_t v3x = (pos3 - pos1).X();
    Double_t v3y = (pos3 - pos1).Y();
    Double_t v3z = (pos3 - pos1).Z();

    Double_t v6x = v1y * v3z - v1z * v3y;
    Double_t v6y = v1z * v3x - v1x * v3z;
    Double_t v6z = v1x * v3y - v1y * v3x;

    Double_t theta56 = acos((v5x * v6x + v5y * v6y + v5z * v6z) / sqrt(v5x * v5x + v5y * v5y + v5z * v5z) /
                            sqrt(v6x * v6x + v6y * v6y + v6z * v6z));
    if ((v5x * v5x + v5y * v5y + v5z * v5z) == 0. || (v6x * v6x + v6y * v6y + v6z * v6z) == 0.)
    {
        theta56 = 3.14;
    }

    Double_t theta3 = acos((v1x * v3x + v1y * v3y + v1z * v3z) / sqrt(v1x * v1x + v1y * v1y + v1z * v1z) /
                           sqrt(v3x * v3x + v3y * v3y + v3z * v3z));
    Double_t theta3min = theta3;
    Double_t theta3max = theta3;

    for (Int_t k = 1; k < 3; k++)
    {
        if (k == 1)
        {
            v3xmin = v3x - dio;
        }
        else
        {
            v3xmin = v3x + dio;
        }
        for (Int_t l = 1; l < 3; l++)
        {
            if (l == 1)
            {
                v3ymin = v3y - dio;
            }
            else
            {
                v3ymin = v3y + dio;
            }
            for (Int_t m = 1; m < 3; m++)
            {
                if (m == 1)
                {
                    v3zmin = v3z - dio;
                }
                else
                {
                    v3zmin = v3z + dio;
                }
                tempAngle =
                    acos((v1x * v3xmin + v1y * v3ymin + v1z * v3zmin) / sqrt(v1x * v1x + v1y * v1y + v1z * v1z) /
                         sqrt(v3xmin * v3xmin + v3ymin * v3ymin + v3zmin * v3zmin));
                if (tempAngle > theta3max)
                    theta3max = tempAngle;
                if (tempAngle < theta3min)
                    theta3min = tempAngle;
            }
        }
    }

    if (theta3max > 1.55)
        theta3max = 1.55;

    dt = c2->GetStartT() - c1->GetStartT();
    dr = sqrt(v3x * v3x + v3y * v3y + v3z * v3z);

    Double_t beta3 = dr / dt / c;
    Double_t beta3max = (dr + dio) / dt / c;
    Double_t beta3min = (dr - dio) / dt / c;

    if (beta3max > 0.99)
        beta3max = 0.99;
    if (beta3min < 0.)
        beta3min = 0.;

    Double_t gamma3 = 1. / sqrt(1. - beta3 * beta3);
    Double_t gamma3min = 1. / sqrt(1. - beta3min * beta3min);
    Double_t gamma3max = 1. / sqrt(1. - beta3max * beta3max);
    Double_t p3 = beta3 * gamma3 * 1. * amu;
    Double_t p3min = beta3min * gamma3min * 1. * amu;
    Double_t p3max = beta3max * gamma3max * 1. * amu;
    if (p3 > p1)
        p3 = p1;
    //  Double_t En3=sqrt(p3*p3+amu*amu);
    //  Double_t En3min=sqrt(p3min*p3min+amu*amu);
    //  Double_t En3max=sqrt(p3max*p3max+amu*amu);
    //  Double_t K3=En3-amu;
    //  Double_t K3min=En3min-amu;
    //  Double_t K3max=En3max-amu;

    Double_t Ma, Mb, Mc, Md, Ka, Thc, Ei, Pa, AA, BB, a, b, cc, Pc1, Pc2;
    Double_t Pd1, Thd1, Ec1, Ed1, Kc1, Kd1, Qsqr1;
    //  Double_t Pd2,Thd2,Ec2,Ed2,Kc2,Kd2,Qsqr2;
    Ma = 1.0087 * amu;
    Mb = 1.0073 * amu;
    Mc = Ma;
    Md = Mb;
    // calculate momentum of scattered proton and neutron
    Ka = K1;
    Thc = theta3;
    Ei = Ma + Mb + Ka;
    Pa = sqrt(Ka * Ka + 2. * Ka * Ma);
    AA = Ei * Ei - Md * Md + Mc * Mc - Pa * Pa;
    BB = AA * AA - 4. * Ei * Ei * Mc * Mc;
    a = 4. * Pa * Pa * cos(Thc) * cos(Thc) - 4. * Ei * Ei;
    b = 4. * AA * Pa * cos(Thc);
    cc = BB;
    Pc1 = -b / (2. * a) + sqrt((b * b) / (4. * a * a) - (cc / a));
    if (Pc1 < 0.)
        Pc1 = 0.;
    Pc2 = -b / (2. * a) - sqrt((b * b) / (4. * a * a) - (cc / a));
    Pd1 = sqrt(Pc1 * Pc1 + Pa * Pa - 2. * Pc1 * Pa * cos(Thc));
    Thd1 = acos((Pc1 * Pc1 - Pd1 * Pd1 - Pa * Pa) / (-2. * Pd1 * Pa));
    Ec1 = sqrt(Pc1 * Pc1 + Ma * Ma);
    Ed1 = sqrt(Pd1 * Pd1 + Mb * Mb);
    Kc1 = Ec1 - Ma;
    Kd1 = Ed1 - Mb;
    Qsqr1 = (-(Ka - Kc1) * (Ka - Kc1) + (Pa * Pa + Pc1 * Pc1 - 2. * Pa * Pc1 * cos(Thc))) / 197.327 / 197.327;
    Double_t p3b = Pc1;
    //  Double_t K3b=Kc1;
    //  Double_t E3b=Ec1;
    Double_t theta4 = Thd1;
    Double_t p4b = Pd1;
    //  Double_t K4b=Kd1;
    //  Double_t E4b=Ed1;

    Double_t p3bmin = p3b;
    Double_t p3bmax = p3b;
    Double_t p4bmin = p4b;
    Double_t p4bmax = p4b;
    Double_t theta4min = theta4;
    Double_t theta4max = theta4;

    // calculate minimum and maximum(within errors) momentum of scattered proton and neutron
    for (Int_t m = 1; m < 5; m++)
    {
        if (m == 1)
        {
            Ka = K1min;
            Thc = theta3min;
        }
        if (m == 2)
        {
            Ka = K1min;
            Thc = theta3max;
        }
        if (m == 3)
        {
            Ka = K1max;
            Thc = theta3min;
        }
        if (m == 4)
        {
            Ka = K1max;
            Thc = theta3max;
        }

        Ei = Ma + Mb + Ka;
        Pa = sqrt(Ka * Ka + 2. * Ka * Ma);
        AA = Ei * Ei - Md * Md + Mc * Mc - Pa * Pa;
        BB = AA * AA - 4. * Ei * Ei * Mc * Mc;
        a = 4. * Pa * Pa * cos(Thc) * cos(Thc) - 4. * Ei * Ei;
        b = 4. * AA * Pa * cos(Thc);
        cc = BB;
        Pc1 = -b / (2. * a) + sqrt((b * b) / (4. * a * a) - (cc / a));
        if (((b * b) / (4. * a * a) - (cc / a)) < 0.)
            Pc1 = 0.;
        if (Pc1 < 0.)
            Pc1 = 0.;
        Pc2 = -b / (2. * a) - sqrt((b * b) / (4. * a * a) - (cc / a));
        Pd1 = sqrt(Pc1 * Pc1 + Pa * Pa - 2. * Pc1 * Pa * cos(Thc));
        Thd1 = acos((Pc1 * Pc1 - Pd1 * Pd1 - Pa * Pa) / (-2. * Pd1 * Pa));

        if (Pc1 < p3bmin)
            p3bmin = Pc1;
        if (Pc1 > p3bmax)
            p3bmax = Pc1;
        if (Pd1 < p4bmin)
            p4bmin = Pd1;
        if (Pd1 > p4bmax)
            p4bmax = Pd1;
        if (Thd1 < theta4min)
            theta4min = Thd1;
        if (Thd1 > theta4max)
            theta4max = Thd1;
    }

    Ec1 = sqrt(p3bmin * p3bmin + Ma * Ma);
    Kc1 = Ec1 - Ma;
    Ed1 = sqrt(p4bmin * p4bmin + Mb * Mb);
    Kd1 = Ed1 - Mb;
    Qsqr1 = (-(Ka - Kc1) * (Ka - Kc1) + (Pa * Pa + Pc1 * Pc1 - 2. * Pa * Pc1 * cos(Thc))) / 197.327 / 197.327;
    //  Double_t K3bmin=Kc1;
    //  Double_t E3bmin=Ec1;
    Double_t K4bmin = Kd1;
    //  Double_t E4bmin=Ed1;
    Ec1 = sqrt(p3bmax * p3bmax + Ma * Ma);
    Kc1 = Ec1 - Ma;
    Ed1 = sqrt(p4bmax * p4bmax + Mb * Mb);
    Kd1 = Ed1 - Mb;
    Qsqr1 = (-(Ka - Kc1) * (Ka - Kc1) + (Pa * Pa + Pc1 * Pc1 - 2. * Pa * Pc1 * cos(Thc))) / 197.327 / 197.327;

    //  Double_t K3bmax=Kc1;
    //  Double_t E3bmax=Ec1;
    Double_t K4bmax = Kd1;
    //  Double_t E4bmax=Ed1;

    // decide if Cluster comes from scattered neutron or another neutron is needed!
    Double_t protonEnergy = 8.76839 + 4.13858 * c1->GetE() - 0.00337368 * c1->GetE() * c1->GetE();
    if (p3bmax > p3min && p3bmin < p3max && theta4max > theta4Measuredmin && theta4min < theta4Measuredmax &&
        K4bmax > protonEnergy && K4bmin < protonEnergy && theta56 * 180. / 3.14 > 120.)
    {
        return kTRUE;
    }
    else
    {
        return kFALSE;
    }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::CalculateMassInv()
{
    R3BLandFirstHits* fh = dynamic_cast<R3BLandFirstHits*>(fLandFirstHits->At(0));
    TVector3 fhv[6];
    fhv[0].SetXYZ(fh->GetX0(), fh->GetY0(), fh->GetZ0());
    fhv[1].SetXYZ(fh->GetX1(), fh->GetY1(), fh->GetZ1());
    fhv[2].SetXYZ(fh->GetX2(), fh->GetY2(), fh->GetZ2());
    fhv[3].SetXYZ(fh->GetX3(), fh->GetY3(), fh->GetZ3());
    // fhv[4].SetXYZ(fh->GetX4(), fh->GetY4(), fh->GetZ4());
    // fhv[5].SetXYZ(fh->GetX5(), fh->GetY5(), fh->GetZ5());

    Double_t sum_energy_tr = 0.;
    Double_t sum_momentumX_tr = 0.;
    Double_t sum_momentumY_tr = 0.;
    Double_t sum_momentumZ_tr = 0.;
    Double_t sum_masses_tr = 0.;
    Double_t sum_energy = 0.;
    Double_t sum_momentumX = 0.;
    Double_t sum_momentumY = 0.;
    Double_t sum_momentumZ = 0.;
    Double_t sum_masses = 0.;

    for (Int_t i = 0; i < fNPrimNeut; i++)
    {
        TVector3 mom;
        PRIM_part[i]->Momentum(mom);
        Double_t momentumZ = mom.Z();
        Double_t momentumX = mom.X();
        Double_t momentumY = mom.Y();
        Double_t energy = sqrt(mom.Mag2() + mNeutron * mNeutron);

        sum_energy_tr += energy;
        sum_momentumX_tr += momentumX;
        sum_momentumY_tr += momentumY;
        sum_momentumZ_tr += momentumZ;
        sum_masses_tr += mNeutron;
    }

    // for (Int_t i = 0; i < 1/*fNPrimProt*/; i++) { // !!!!!!!!!!!!! FIXME
    //   TVector3 mom;
    //   PRIM_prot[i]->Momentum(mom);
    //   Double_t momentumZ = mom.Z();
    //   Double_t momentumX = mom.X();
    //   Double_t momentumY = mom.Y();
    //   Double_t energy = sqrt(mom.Mag2()+PRIM_prot[i]->GetM()*PRIM_prot[i]->GetM());

    //   sum_energy_tr += energy;
    //   sum_momentumX_tr += momentumX;
    //   sum_momentumY_tr += momentumY;
    //   sum_momentumZ_tr += momentumZ;
    //   sum_masses_tr += PRIM_prot[i]->GetM();

    //   sum_energy += energy;
    //   sum_momentumX += momentumX;
    //   sum_momentumY += momentumY;
    //   sum_momentumZ += momentumZ;
    //   sum_masses += PRIM_prot[i]->GetM();
    // }

    R3BNeutHit* hit;
    for (Int_t i = 0; i < fNeutHits->GetEntriesFast(); i++)
    {
        // add up momentum of neutrons
        // neutron: beta is calculated with position and time of first hit

        hit = dynamic_cast<R3BNeutHit*>(fNeutHits->At(i));

        Double_t momentumT = hit->GetP();
        Double_t momentumZ =
            momentumT *
            cos(TMath::ATan2(TMath::Sqrt(hit->GetX() * hit->GetX() + hit->GetY() * hit->GetY()), hit->GetZ()));
        Double_t momentumX = momentumZ * hit->GetX() / hit->GetZ();
        Double_t momentumY = momentumZ * hit->GetY() / hit->GetZ();
        Double_t energy = TMath::Sqrt(momentumT * momentumT + mNeutron * mNeutron);

        sum_energy += energy;
        sum_momentumX += momentumX;
        sum_momentumY += momentumY;
        sum_momentumZ += momentumZ;
        sum_masses += mNeutron;

        R3BLandPoint* point;
        TVector3 posPoint;
        Double_t d;
        Int_t index = -1;
        Double_t dmin = 1e6;
        for (Int_t ip = 0; ip < fLandPoints->GetEntriesFast(); ip++)
        {
            point = dynamic_cast<R3BLandPoint*>(fLandPoints->At(ip));
            point->PositionIn(posPoint);
            d = TMath::Sqrt(TMath::Power(hit->GetX() - posPoint.X(), 2) + TMath::Power(hit->GetY() - posPoint.Y(), 2) +
                            TMath::Power(hit->GetZ() - posPoint.Z(), 2));
            if (d < dmin)
            {
                dmin = d;
                index = ip;
            }
        }
        assert(index != -1);
        point = dynamic_cast<R3BLandPoint*>(fLandPoints->At(index));
        Int_t trackID = point->GetTrackID();
        R3BMCTrack* track = dynamic_cast<R3BMCTrack*>(fLandMCTrack->At(trackID));
        Int_t mothID = track->GetMotherId();
        Int_t ind;
        while (-1 != mothID)
        {
            track = dynamic_cast<R3BMCTrack*>(fLandMCTrack->At(mothID));
            ind = mothID;
            mothID = track->GetMotherId();
        }
        hDeltaPx2->Fill(momentumX - track->GetPx() * 1000.);
        hDeltaPy2->Fill(momentumY - track->GetPy() * 1000.);
        hDeltaPz2->Fill(momentumZ - track->GetPz() * 1000.);

        Double_t deltaP = momentumT - track->GetP() * 1000.;
        if (0 == i)
        {
            hDeltaP1->Fill(deltaP);
        }
        else if (1 == i)
        {
            hDeltaP2->Fill(deltaP);
        }
        else if (2 == i)
        {
            hDeltaP3->Fill(deltaP);
        }
        else if (3 == i)
        {
            hDeltaP4->Fill(deltaP);
        }

        // TVector3 mom;
        // PRIM_part[0]->Momentum(mom);
        // hDeltaPx2->Fill(momentumX-mom.X());
        // hDeltaPy2->Fill(momentumY-mom.Y());
        // hDeltaPz2->Fill(momentumZ-mom.Z());

        dmin = 1e6;
        for (Int_t ifh = 0; ifh < 4; ifh++)
        {
            d = TMath::Sqrt(TMath::Power(hit->GetX() - fhv[ifh].X(), 2) + TMath::Power(hit->GetY() - fhv[ifh].Y(), 2) +
                            TMath::Power(hit->GetZ() - fhv[ifh].Z(), 2));
            if (d < dmin)
            {
                dmin = d;
                index = ifh;
            }
        } // first hits

        // hDeltaX->Fill((pos-fhv[1]).X());
        // hDeltaY->Fill((pos-fhv[1]).Y());
        // hDeltaZ->Fill((pos-fhv[1]).Z());

        hDeltaX->Fill(hit->GetX() - fhv[index].X());
        hDeltaY->Fill(hit->GetY() - fhv[index].Y());
        hDeltaZ->Fill(hit->GetZ() - fhv[index].Z());
    } // tracks

    // fragment: This information have to come later from the Tracker
    if (fNofFrag > 0)
    {
        TVector3 momF;
        PRIM_frag[0]->Momentum(momF);
        Double_t beta_frag = momF.Mag() / sqrt(momF.Mag2() + PRIM_frag[0]->GetM() * PRIM_frag[0]->GetM());
        Double_t gamma_frag = 1. / sqrt(1. - beta_frag * beta_frag);
        Double_t energy_frag = gamma_frag * PRIM_frag[0]->GetM();
        sum_energy += energy_frag;
        sum_momentumX += momF.X();
        sum_momentumY += momF.Y();
        sum_momentumZ += momF.Z();
        sum_masses += PRIM_frag[0]->GetM();
        sum_energy_tr += energy_frag;
        sum_momentumX_tr += momF.X();
        sum_momentumY_tr += momF.Y();
        sum_momentumZ_tr += momF.Z();
        sum_masses_tr += PRIM_frag[0]->GetM();
    }

    // Gamma
    TVector3 momGamma;
    for (Int_t ig = 0; ig < fNPrimGamma; ig++)
    {
        PRIM_gamma[ig]->Momentum(momGamma);
        sum_energy += momGamma.Mag();
        sum_momentumX += momGamma.X();
        sum_momentumY += momGamma.Y();
        sum_momentumZ += momGamma.Z();
        sum_energy_tr += momGamma.Mag();
        sum_momentumX_tr += momGamma.X();
        sum_momentumY_tr += momGamma.Y();
        sum_momentumZ_tr += momGamma.Z();
    }

    // calculate invariant mass
    // calculate relative energy
    Double_t totMom2_tr =
        sum_momentumX_tr * sum_momentumX_tr + sum_momentumY_tr * sum_momentumY_tr + sum_momentumZ_tr * sum_momentumZ_tr;
    Double_t totMom2 = sum_momentumX * sum_momentumX + sum_momentumY * sum_momentumY + sum_momentumZ * sum_momentumZ;
    fMinvTrue = (sqrt(sum_energy_tr * sum_energy_tr - totMom2_tr) - sum_masses_tr);
    fMinv = (sqrt(sum_energy * sum_energy - totMom2) - sum_masses);
    // cout << fMinvTrue << "   " << fMinv << endl;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::CalculateExce()
{
    TVector3 mom1;
    TVector3 mom2;
    Double_t en_g = 0;
    for (Int_t i = 0; i < fNPrimGamma; i++)
    {
        PRIM_gamma[i]->Momentum(mom1);
        en_g += mom1.Mag();
    }

    Int_t N = fNofFrag + fNPrimProt + fNPrimNeut - 1; // !!!!!!!!!!!!! FIXME
    R3BPrimPart* p[30];

    Int_t ii = 0;
    for (Int_t i = 0; i < fNofFrag; i++)
    {
        p[ii] = PRIM_frag[i];
        ii += 1;
    }
    for (Int_t i = 0; i < 1 /*fNPrimProt*/; i++)
    { // !!!!!!!!!!!!! FIXME
        p[ii] = PRIM_prot[i];
        ii += 1;
    }
    for (Int_t i = 0; i < fNPrimNeut; i++)
    {
        p[ii] = PRIM_part[i];
        ii += 1;
    }

    // MC truth - ideal neutron reconstruction -------------------------
    Double_t sum_mc = 0.;
    for (Int_t i = 0; i < N; i++)
    {
        p[i]->Momentum(mom1);
        for (Int_t j = 0; j < N; j++)
        {
            p[j]->Momentum(mom2);
            sum_mc += p[i]->GetGamma() * p[j]->GetGamma() * p[i]->GetM() * p[j]->GetM() *
                      (1. - p[i]->GetBeta() * p[j]->GetBeta() * cos(mom1.Angle(mom2)));
        }
    }
    // Double_t m_proj = 135.939300*amu;
    Double_t m_proj = 2.0141017778 * amu;
    fExceTrue = sqrt(sum_mc) + en_g - m_proj;
    // -----------------------------------------------------------------

    // RECO spectrum ---------------------------------------------------
    if (!fNofTracks)
    {
        fExce = -1.;
        return;
    }

    N = fNofFrag + fNPrimProt + fNofTracks - 1; // !!!!!!!!!!!!! FIXME
    ii = 0;
    for (Int_t i = 0; i < fNofFrag; i++)
    {
        p[ii] = PRIM_frag[i];
        ii += 1;
    }
    for (Int_t i = 0; i < 1 /*fNPrimProt*/; i++)
    { // !!!!!!!!!!!!! FIXME
        p[ii] = PRIM_prot[i];
        ii += 1;
    }
    Int_t k = ii;
    for (Int_t i = 0; i < fNofTracks; i++)
    {
        R3BNeuLandCluster* c1 = fVectorCluster.at(fTracks[i][0]);
        TVector3 pos;
        c1->StartPosition(pos);
        Double_t s = pos.Mag();
        Double_t rr = pos.Perp();

        Double_t beta = s / c1->GetStartT() / c;
        Double_t gamma = 1. / sqrt(1. - beta * beta);
        Double_t momT = beta * gamma * mNeutron;
        Double_t momZ = momT * cos(TMath::ATan2(rr, pos.Z()));
        Double_t momX = momZ * pos.X() / pos.Z();
        Double_t momY = momZ * pos.Y() / pos.Z();
        p[k] = new R3BPrimPart(2112, momX, momY, momZ, 0., 0., 0., 0., 1.0086649, mNeutron);
        k += 1;
    }
    Double_t sum = 0.;
    for (Int_t i = 0; i < N; i++)
    {
        p[i]->Momentum(mom1);
        for (Int_t j = 0; j < N; j++)
        {
            p[j]->Momentum(mom2);
            sum += p[i]->GetGamma() * p[j]->GetGamma() * p[i]->GetM() * p[j]->GetM() *
                   (1. - p[i]->GetBeta() * p[j]->GetBeta() * cos(mom1.Angle(mom2)));
        }
    }
    fExce = sqrt(sum) + en_g - m_proj;
    for (Int_t i = ii; i < N; i++)
    {
        delete p[i];
    }
    // -----------------------------------------------------------------
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::Reset()
{
    for (Int_t i = 0; i < fNPrimNeut; i++)
    {
        delete PRIM_part[i];
    }
    for (Int_t i = 0; i < fNPrimProt; i++)
    {
        delete PRIM_prot[i];
    }
    for (Int_t i = 0; i < fNPrimGamma; i++)
    {
        delete PRIM_gamma[i];
    }
    for (Int_t i = 0; i < fNofFrag; i++)
    {
        delete PRIM_frag[i];
    }

    fMinvTrue = 0.;
    fMinv = 0.;
    fExceTrue = 0.;
    fExce = 0.;
    fEtot = 0.;

    fNPrimNeut = 0;
    fNPrimProt = 0;
    fNPrimGamma = 0;
    fNofFrag = 0;
    fNofTracks = 0;
    fNofClusters = 0;
    fNofClusters1 = 0;

    if (fVectorCluster.size() > 0)
    {
        fVectorCluster.clear();
    }

    // Reset the neutron tracks
    for (Int_t i = 0; i < 1000; i++)
    {
        fMapTracks[i] = kTRUE;
        for (Int_t k = 0; k < 300; k++)
        {
            fTracks[i][k] = -1;
        }
        fNofHits[i] = 0;
    }

    fNeutHits->Clear();
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::Finish()
{
    hClusters1->Write();

    hNeutmult->Write();

    hMinv->Write();
    hMinv0->Write();
    hMinv1->Write();
    hMinv2->Write();

    hExce->Write();
    hExce0->Write();

    hErel1->Write();
    hErel2->Write();
    hErel3->Write();
    hErel4->Write();

    hDeltaPx1->Write();
    hDeltaPy1->Write();
    hDeltaPz1->Write();
    hDeltaPx2->Write();
    hDeltaPy2->Write();
    hDeltaPz2->Write();

    hDeltaX->Write();
    hDeltaY->Write();
    hDeltaZ->Write();
    hDeltaT->Write();

    hDeltaP1->Write();
    hDeltaP2->Write();
    hDeltaP3->Write();
    hDeltaP4->Write();
    hDeltaP5->Write();
    hDeltaP6->Write();

    hDelta->Write();

    hFirstHitZ->Write();

    h_theta_cand->Write();
    h_beta_min->Write();
    h_beta_max->Write();
    h_ntracks->Write();
    h_nhits->Write();
    h_theta->Write();

    h_nofclusters->Write();
    h_etot->Write();
    h_ncl_etot->Write();
    h_ncl_etot_1->Write();
    h_ndigi_etot->Write();
    h_ncl1_etot->Write();
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::UseBeam(Double_t beam_energy, Double_t beam_beta)
{
    beamEnergy = beam_energy;
    beamBeta = beam_beta;
    gBeamBeta = beam_beta;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void R3BNeutronTracker2D::CreateHistograms()
{
    hNeutmult = new TH1F("Neutmult", "Neutron multiplicity from energy considerations", 10, -0.5, 9.5);
    hNeutmult->GetXaxis()->SetTitle("Number of Neutrons");
    hNeutmult->GetYaxis()->SetTitle("Counts");

    hErel1 = new TH1F("Erel1", "Erel for 1 neutron", 400, 0., 40.);
    hErel1->GetXaxis()->SetTitle("Erel (MeV)");
    hErel1->GetYaxis()->SetTitle("Counts");

    hErel2 = new TH1F("Erel2", "Erel for 2 neutron", 400, 0., 40.);
    hErel2->GetXaxis()->SetTitle("Erel (MeV)");
    hErel2->GetYaxis()->SetTitle("Counts");

    hErel3 = new TH1F("Erel3", "Erel for 3 neutron", 400, 0., 40.);
    hErel3->GetXaxis()->SetTitle("Erel (MeV)");
    hErel3->GetYaxis()->SetTitle("Counts");

    hErel4 = new TH1F("Erel4", "Erel for 4 neutron", 400, 0., 40.);
    hErel4->GetXaxis()->SetTitle("Erel (MeV)");
    hErel4->GetYaxis()->SetTitle("Counts");

    hMinv = new TH1F("Minv", "Minv for reconstructed hits", 4000, -0.01, 39.99);
    // hMinv = new TH1F("Minv","Minv for reconstructed hits",1000,0.,1000.);
    hMinv->GetXaxis()->SetTitle("Erel (MeV)");
    hMinv->GetYaxis()->SetTitle("Counts");

    hExce = new TH1F("Exce", "Reconstructed Exce", 100, 0., 600.);

    hMinv1 = new TH1F("Minv1", "Minv for first hits", 4000, 0., 40.);
    hMinv1->GetXaxis()->SetTitle("Erel (MeV)");
    hMinv1->GetYaxis()->SetTitle("Counts");

    hMinv2 = new TH1F("Minv2", "Minv for selected hits", 4000, 0., 40.);
    hMinv2->GetXaxis()->SetTitle("Erel (MeV)");
    hMinv2->GetYaxis()->SetTitle("Counts");

    hMinv0 = new TH1F("Minv0", "Minv for ideal hits", 2000, -0.01, 39.99);
    // hMinv0 = new TH1F("Minv0","Minv for ideal hits",1000,0.,1000.);
    hMinv0->GetXaxis()->SetTitle("Erel (MeV)");
    hMinv0->GetYaxis()->SetTitle("Counts");

    hExce0 = new TH1F("Exce0", "Input Exce", 100, 0., 600.);

    hDeltaX = new TH1F("DeltaX", "error in x determination", 300, -150., 150.);
    hDeltaX->GetXaxis()->SetTitle("x position (cm)");
    hDeltaX->GetYaxis()->SetTitle("Counts");

    hDeltaY = new TH1F("DeltaY", "error in y determination", 300, -150., 150.);
    hDeltaY->GetXaxis()->SetTitle("y position (cm)");
    hDeltaY->GetYaxis()->SetTitle("Counts");

    hDeltaZ = new TH1F("DeltaZ", "error in z determination", 500, -250., 250.);
    hDeltaZ->GetXaxis()->SetTitle("z position (cm)");
    hDeltaZ->GetYaxis()->SetTitle("Counts");

    hDeltaT = new TH1F("DeltaT", "error in time determination", 2000, -10., 10.);
    hDeltaT->GetXaxis()->SetTitle("time (ns)");
    hDeltaT->GetYaxis()->SetTitle("Counts");

    hDeltaP1 = new TH1F("DeltaP1", "difference in reconstucted momenta for 1st neutron)", 1000, -150., 150.);
    hDeltaP1->GetXaxis()->SetTitle("Delta P (MeV/c)");
    hDeltaP1->GetYaxis()->SetTitle("Counts");

    hDeltaP2 = new TH1F("DeltaP2", "difference in reconstucted momenta for 2nd neutron)", 1000, -150., 150.);
    hDeltaP2->GetXaxis()->SetTitle("Delta P (MeV/c)");
    hDeltaP2->GetYaxis()->SetTitle("Counts");

    hDeltaP3 = new TH1F("DeltaP3", "difference in reconstucted momenta for 3rd neutron)", 1000, -150., 150.);
    hDeltaP3->GetXaxis()->SetTitle("Delta P (MeV/c)");
    hDeltaP3->GetYaxis()->SetTitle("Counts");

    hDeltaP4 = new TH1F("DeltaP4", "difference in reconstucted momenta for 4th neutron)", 1000, -150., 150.);
    hDeltaP4->GetXaxis()->SetTitle("Delta P (MeV/c)");
    hDeltaP4->GetYaxis()->SetTitle("Counts");

    hDeltaP5 = new TH1F("DeltaP5", "difference in reconstucted momenta for 5th neutron)", 1000, -150., 150.);
    hDeltaP5->GetXaxis()->SetTitle("Delta P (MeV/c)");
    hDeltaP5->GetYaxis()->SetTitle("Counts");

    hDeltaP6 = new TH1F("DeltaP6", "difference in reconstucted momenta for 6th neutron)", 1000, -150., 150.);
    hDeltaP6->GetXaxis()->SetTitle("Delta P (MeV/c)");
    hDeltaP6->GetYaxis()->SetTitle("Counts");

    hDeltaPx1 = new TH1F("DeltaPx1", "difference in reconstucted momenta px (ideal case)", 1000, -50., 50.);
    hDeltaPx1->GetXaxis()->SetTitle("Delta Px (MeV/c)");
    hDeltaPx1->GetYaxis()->SetTitle("Counts");

    hDeltaPy1 = new TH1F("DeltaPy1", "difference in reconstucted momenta py (ideal case)", 1000, -50., 50.);
    hDeltaPy1->GetXaxis()->SetTitle("Delta Py (MeV/c)");
    hDeltaPy1->GetYaxis()->SetTitle("Counts");

    hDeltaPz1 = new TH1F("DeltaPz1", "difference in reconstucted momenta pz (ideal case)", 1000, -150., 150.);
    hDeltaPz1->GetXaxis()->SetTitle("Delta Pz (MeV/c)");
    hDeltaPz1->GetYaxis()->SetTitle("Counts");

    hDeltaPx2 = new TH1F("DeltaPx2", "difference in reconstucted momenta px (exp case)", 1000, -50., 50.);
    hDeltaPx2->GetXaxis()->SetTitle("Delta Px (MeV/c)");
    hDeltaPx2->GetYaxis()->SetTitle("Counts");

    hDeltaPy2 = new TH1F("DeltaPy2", "difference in reconstucted momenta py (exp case)", 1000, -50., 50.);
    hDeltaPy2->GetXaxis()->SetTitle("Delta Py (MeV/c)");
    hDeltaPy2->GetYaxis()->SetTitle("Counts");

    hDeltaPz2 = new TH1F("DeltaPz2", "difference in reconstucted momenta pz (exp case)", 1000, -150., 150.);
    hDeltaPz2->GetXaxis()->SetTitle("Delta Pz (MeV/c)");
    hDeltaPz2->GetYaxis()->SetTitle("Counts");

    hClusters1 = new TH1F("Cluster1", "Number of clusters after eliminating elastic scattering", 100, -0.5, 99.5);
    hClusters1->GetXaxis()->SetTitle("number of clusters");
    hClusters1->GetYaxis()->SetTitle("Counts");

    hDelta = new TH1F("Delta", "distance between two primary interactions", 300, -150., 150.);
    hDelta->GetXaxis()->SetTitle("distance (cm)");
    hDelta->GetYaxis()->SetTitle("Counts");

    hFirstHitZ = new TH1F("FirstHitZ", "z positions of first hits", 200, 1000., 2000.);
    hFirstHitZ->GetXaxis()->SetTitle("z position (cm)");
    hFirstHitZ->GetYaxis()->SetTitle("Counts");

    h_theta_cand = new TH1F("h_theta_cand", "Polar angle of a track candidate", 180, 0., 180.);
    h_theta_cand->GetXaxis()->SetTitle("#theta (deg)");
    h_theta_cand->GetYaxis()->SetTitle("counts");

    h_beta_min = new TH1F("h_beta_min", "velocity of the potential neutron", 400, -2., 2.);
    h_beta_min->GetXaxis()->SetTitle("#beta_{min}");
    h_beta_min->GetYaxis()->SetTitle("counts");

    h_beta_max = new TH1F("h_beta_max", "velocity of the potential neutron", 1200, -2., 10.);
    h_beta_max->GetXaxis()->SetTitle("#beta_{max}");
    h_beta_max->GetYaxis()->SetTitle("counts");

    h_ntracks = new TH1F("h_ntracks", "Number of reconstructed neutron tracks", 50, -0.5, 49.5);
    h_ntracks->GetXaxis()->SetTitle("N_{tracks}");
    h_ntracks->GetYaxis()->SetTitle("events");

    h_nhits = new TH1F("h_nhits", "Number of hits on a neutron track", 50, -0.5, 49.5);
    h_nhits->GetXaxis()->SetTitle("N_{hits} / track");
    h_nhits->GetYaxis()->SetTitle("all tracks");

    h_theta = new TH1F("h_theta", "Polar angle", 180, 0., 180.);
    h_theta->GetXaxis()->SetTitle("#theta (deg)");
    h_theta->GetYaxis()->SetTitle("counts");

    h_nofclusters = new TH1F("h_nofclusters", "Number of clusters", 100, -0.5, 99.5);

    h_etot = new TH1F("h_etot", "Total light", 100, 0., 1000.);

    if (beamEnergy < 300)
    {
        h_ncl_etot = new TH2F("h_ncl_etot", "Number of clusters vs. total light", 1500, 0., 1500., 150, -0.5, 149.5);
    }
    else
    {
        h_ncl_etot = new TH2F("h_ncl_etot", "Number of clusters vs. total light", 150, 0., 1500., 150, -0.5, 149.5);
    }

    h_ncl_etot_1 =
        new TH2F("h_ncl_etot_1", "Number of clusters vs. total light after cut", 100, 0., 1000., 100, -0.5, 99.5);

    h_ndigi_etot = new TH2F("h_ndigi_etot", "", 100, 0., 1000., 500, -0.5, 499.5);

    h_ncl1_etot = new TH2F("h_ncl1_etot", "", 100, 0., 1000., 100, -0.5, 99.5);
}
// -----------------------------------------------------------------------------

bool AuxSortClustersBeta(R3BNeuLandCluster* c1, R3BNeuLandCluster* c2)
{
    TVector3 pos1;
    TVector3 pos2;
    c1->StartPosition(pos1);
    c2->StartPosition(pos2);
    Double_t beta1 = pos1.Z() / c1->GetStartT() / c;
    Double_t beta2 = pos2.Z() / c2->GetStartT() / c;
    Double_t val1 = TMath::Abs(beta1 - gBeamBeta) / c1->GetE();
    Double_t val2 = TMath::Abs(beta2 - gBeamBeta) / c2->GetE();
    if (val1 < val2)
    {
        return true;
    }
    else
    {
        return false;
    }
}
