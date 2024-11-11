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

#include "R3BTPropagator.h"
#include "FairLogger.h"
#include "FairRKPropagator.h"
#include "FairRunAna.h"
#include "R3BGladFieldMap.h"
#include "R3BTGeoPar.h"
#include "R3BTrackingDetector.h"
#include "R3BTrackingParticle.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMath.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

R3BTPropagator::R3BTPropagator(R3BGladFieldMap* field, Bool_t vis)
    : fFairProp(new FairRKPropagator(field))
    , fField(field)
    , fmTofGeo(NULL)
    , fVis(vis)
{
    // Define magnetic field boundaries ------------------------------------
    TVector3 pos(field->GetPositionX(), field->GetPositionY(), field->GetPositionZ());
    Double_t angle = field->GetYAngle() * TMath::DegToRad();
    fPlane1[0] = TVector3(0., 0., field->GetZmin());
    fPlane1[1] = TVector3(field->GetXmax(), field->GetYmax(), field->GetZmin());
    fPlane1[2] = TVector3(-field->GetXmax(), field->GetYmax(), field->GetZmin());
    fPlane2[0] = TVector3(0., 0., field->GetZmax());
    fPlane2[1] = TVector3(field->GetXmax(), field->GetYmax(), field->GetZmax());
    fPlane2[2] = TVector3(-field->GetXmax(), field->GetYmax(), field->GetZmax());

    for (Int_t i = 0; i < 3; i++)
    {
        fPlane1[i].RotateY(angle);
        fPlane1[i] += pos;
        fPlane2[i].RotateY(angle);
        fPlane2[i] += pos;
    }
    fNorm1 = ((fPlane1[1] - fPlane1[0]).Cross(fPlane1[2] - fPlane1[0])).Unit();
    fNorm2 = ((fPlane2[1] - fPlane2[0]).Cross(fPlane2[2] - fPlane2[0])).Unit();
    //----------------------------------------------------------------------

    //  fVis = kTRUE;
    if (fVis)
    {
        fc4 = new TCanvas("c4", "", 50, 50, 800, 800);
        TH2F* h3 = new TH2F("h3", "", 180, -450., 450., 180, -100., 800.);
        h3->SetStats(kFALSE);
        h3->Draw();

        TLine* l1 = new TLine(-fPlane1[1].X(), fPlane1[1].Z(), -fPlane1[2].X(), fPlane1[2].Z());
        l1->SetLineColor(2);
        l1->Draw();
        TLine* l2 = new TLine(-fPlane2[1].X(), fPlane2[1].Z(), -fPlane2[2].X(), fPlane2[2].Z());
        l2->SetLineColor(2);
        l2->Draw();
        TLine* l3 = new TLine(-fPlane1[1].X(), fPlane1[1].Z(), -fPlane2[1].X(), fPlane2[1].Z());
        l3->SetLineColor(2);
        l3->Draw();
        TLine* l4 = new TLine(-fPlane1[2].X(), fPlane1[2].Z(), -fPlane2[2].X(), fPlane2[2].Z());
        l4->SetLineColor(2);
        l4->Draw();
    }
}

R3BTPropagator::~R3BTPropagator() {}

Bool_t R3BTPropagator::PropagateToDetector(R3BTrackingParticle* particle, R3BTrackingDetector* detector)
{
    // std::cout<<"In PropagateToDetector: "<<detector->GetDetectorName()<<std::endl;
    // detector->pos0.Print();
    // detector->pos1.Print();
    // detector->pos2.Print();

    if (fVis)
    {
        detector->Draw();
    }

    return PropagateToPlane(particle, detector->pos0, detector->pos1, detector->pos2);
}

Bool_t R3BTPropagator::PropagateToDetectorBackward(R3BTrackingParticle* particle, R3BTrackingDetector* detector)
{

    if (fVis)
        detector->Draw();
    return PropagateToPlaneBackward(particle, detector->pos0, detector->pos1, detector->pos2);
}

Bool_t R3BTPropagator::PropagateToPlane(R3BTrackingParticle* particle,
                                        const TVector3& v1,
                                        const TVector3& v2,
                                        const TVector3& v3)
{
    TVector3 intersect;
    Bool_t crossed;
    Double_t step = 0.;
    Bool_t result = false;

    LOG(debug2) << "Entered PropagateToPlane!" << endl;
    LOG(debug2) << "Magnetic field: " << endl;
    //  fPlane1[0].Print();
    //  fPlane2[1].Print();
    //  fPlane2[0].Print();

    TVector3 norm = ((v2 - v1).Cross(v3 - v1)).Unit(); // normalized vector of vector product of (v2-1) and (v3-v1)
    // Check if particle is already on plane
    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
    TVector3 dist = (particle->GetPosition() - intersect);
    Double_t diff = dist.Mag();

    // if (TMath::Abs((particle->GetPosition() - v1).Dot(norm)) < 1e-6)
    if (TMath::Abs(diff) < 1.e-6)
    {
        LOG(debug2) << "Already at plane!" << endl;
        return kTRUE;
    }

    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), fPlane1[0], fNorm1, intersect);

    if (crossed)
    {
        LOG(debug2) << "crossed at glad entrance: " << crossed << " at:" << endl;
        //  intersect.Print();
        //   particle->GetPosition().Print();
        //   particle->GetMomentum().Print();

        LOG(debug2) << "Starting upstream of magnetic field boundaries...";
        if ((v1 - particle->GetPosition()).Mag() < (fPlane1[0] - particle->GetPosition()).Mag())
        {
            LOG(debug2) << "Propagating to end-plane and stop.";
            crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
            step = (intersect - particle->GetPosition()).Mag();
            //   cout<<"Step 1: "<<step<<", intersect & particle: "<<endl;
            // intersect.Print();
            // particle->GetPosition().Print();
            // particle->GetMomentum().Print();

            particle->SetPosition(intersect);
            particle->AddStep(step);
            return kTRUE;
        }

        LOG(debug2) << "Propagating to entrance of magnetic field.";

        // fVis = kTRUE;
        if (fVis)
        {
            //	cout<<"Line 1 to be plotted: xpart= "<<-particle->GetX()<<"; zpart= "<<particle->GetZ()<<", xinter=
            //"<<-intersect.X()<<
            //	", zinter= "<<intersect.Z()<<endl;
            TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
            l1->SetLineColor(1);
            l1->Draw();
        }

        step = (intersect - particle->GetPosition()).Mag();
        // cout<<"Step 2: "<<step<<", intersect: "<<endl;
        // intersect.Print();
        // particle->GetPosition().Print();
        // particle->GetMomentum().Print();

        particle->SetPosition(intersect);
        particle->AddStep(step);
        LOG(debug2) << intersect.X() << " " << intersect.Y() << " " << intersect.Z();
        // particle->GetPosition().Print();
        // particle->GetMomentum().Print();
    }
    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), fPlane2[0], fNorm2, intersect);

    if (crossed)
    {
        LOG(debug2) << "will cross at the end of Bfeld: " << crossed << " at:";
        //	intersect.Print();
        //	particle->GetPosition().Print();
        //    particle->GetMomentum().Print();

        LOG(debug2) << "####### Propagating inside of field using RK4..." << endl;
        TVector3 tpos;
        // cout << "Test1: dist to det: " << (v1 - particle->GetPosition()).Mag() << ", dist to exit: " << (fPlane2[0] -
        // particle->GetPosition()).Mag() << ", xpart: "<<particle->GetX()<<", ypart: "<<particle->GetY()<<", zpart: "<<
        // particle->GetZ()<<endl;
        //  v1.Print();
        //  v2.Print();
        //   v3.Print();
        if ((v1 - particle->GetPosition()).Mag() < (fPlane2[0] - particle->GetPosition()).Mag())
        {
            LOG(debug2) << "Propagating to end-plane using RK4 and stop. Particle before: " << endl;
            tpos = particle->GetPosition();
            // particle->GetPosition().Print();
            // particle->GetMomentum().Print();
            // LOG(info) << particle->GetCharge();
            result = PropagateToPlaneRK(particle, v1, v2, v3);
            //  cout<<"particle after propagaging: "<<endl;
            // particle->GetPosition().Print();
            // particle->GetMomentum().Print();

            //  fVis = kTRUE;

            if (fVis)
            {
                TLine* l1 = new TLine(-tpos.X(), tpos.Z(), -particle->GetX(), particle->GetZ());
                l1->SetLineColor(1);
                l1->Draw();
            }
            // cout<<"Result: "<<result<<endl;
            return result;
        }

        LOG(debug2) << "Propagating to exit from magnetic field.";
        tpos = particle->GetPosition();
        result = PropagateToPlaneRK(particle, fPlane2[0], fPlane2[1], fPlane2[2]);

        //  fVis = true; // kTRUE;
        if (fVis)
        {
            TLine* l1 = new TLine(-tpos.X(), tpos.Z(), -particle->GetX(), particle->GetZ());
            l1->Draw();
        }
        if (!result)
        {
            return result;
        }
    }
    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
    // cout << "crossed 3rd: " << crossed << " at:";
    // intersect.Print();
    if (crossed)
    {

        LOG(debug2) << "Propagating to end plane. Finish." << endl;

        //  fVis = kTRUE;
        if (fVis)
        {
            TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
            l1->Draw();
        }

        step = (intersect - particle->GetPosition()).Mag(); // * 0.1;

        // cout<<"Step 3: "<<step<<", interesct: "<<endl;
        // intersect.Print();
        particle->SetPosition(intersect);
        particle->AddStep(step);
        LOG(debug2);
        return kTRUE;
    }
    // LOG(error) << "!!! Failed !!!";
    return kFALSE;
}

/*
Bool_t R3BTPropagator::PropagateToPlaneBackward(R3BTrackingParticle* particle,
                                                const TVector3& v1,
                                                const TVector3& v2,
                                                const TVector3& v3)
{

    TVector3 intersect;
    Bool_t crossed;
    Double_t step = 0.;
    Bool_t result = false;

    TVector3 norm = ((v2 - v1).Cross(v3 - v1)).Unit(); // normalized vector of vector product of (v2-1) and (v3-v1)
    // Check if particle is already on plane
    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
    TVector3 dist = (particle->GetPosition() - intersect);
    Double_t diff = dist.Mag();

    // if (TMath::Abs((particle->GetPosition() - v1).Dot(norm)) < 1e-6)
    if (TMath::Abs(diff) < 1.e-6)
    {
        LOG(debug2) <<"Already at plane!"<<endl;
        return kTRUE;
    }

    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), fPlane2[0], -fNorm2, intersect);
    if (crossed)
    {
        LOG(debug2) << "Starting downstream of magnetic field boundaries..."<<endl;

        if ((v1 - particle->GetPosition()).Mag() < (fPlane2[0] - particle->GetPosition()).Mag())
        {
            LOG(debug2) << "Propagating to end-plane and stop."<<endl;
            crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
          //  fVis = kTRUE;
            if (fVis)
            {
                TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
                l1->Draw();
            }

            step = (intersect - particle->GetPosition()).Mag();
            particle->SetPosition(intersect);
            particle->AddStep(step);

            return kTRUE;
        }
        LOG(debug2) << "Propagating to exit of magnetic field."<<endl;
     //   fVis = kTRUE;
        if (fVis)
        {
            TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
            l1->Draw();
        }
        step = (intersect - particle->GetPosition()).Mag();
        particle->SetPosition(intersect);
        particle->AddStep(step);

        LOG(debug2) << intersect.X() << " " << intersect.Y() << " " << intersect.Z();
    }

    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), fPlane1[0], -fNorm1, intersect);
    if (crossed)
    {
       LOG(debug2) << "Propagating inside of field using RK4..."<<endl;

        TVector3 tpos;

        if ((v1 - particle->GetPosition()).Mag() < (fPlane1[0] - particle->GetPosition()).Mag())
        {
            LOG(debug2) << "Propagating to end-plane using RK4 and stop. v1/particle"<<endl;
            tpos = particle->GetPosition();
            result = PropagateToPlaneRK(particle, v1, v3, v2);

         //   fVis = kTRUE;
            if (fVis)
            {
                TLine* l1 = new TLine(-tpos.X(), tpos.Z(), -particle->GetX(), particle->GetZ());
                l1->Draw();
            }
            LOG(debug2) << particle->GetX() << ", " << particle->GetY() << ", " << particle->GetZ();
            return result;
        }
        LOG(debug2) << "Propagating to entrance of magnetic field. Finish."<<endl;
        tpos = particle->GetPosition();
        result = PropagateToPlaneRK(particle, fPlane1[0], fPlane1[1], fPlane1[2]);

      //  fVis = kTRUE;
        if (fVis)
        {
            TLine* l1 = new TLine(-tpos.X(), tpos.Z(), -particle->GetX(), particle->GetZ());
            l1->Draw();
        }

        LOG(debug2) << particle->GetX() << ", " << particle->GetY() << ", " << particle->GetZ();
        if (!result)
        {
            return result;
        }
    }

    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
 / *   if (crossed)
    {
       // LOG(debug2)
        cout<< "Propagating to end plane. Finish.";
     //   fVis = kTRUE;
        if (fVis)
        {
            TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
            l1->Draw();
        }
        step = (intersect - particle->GetPosition()).Mag();
        particle->SetPosition(intersect);
        particle->AddStep(step);
        LOG(debug2);
        return kTRUE;
    } * /

    // LOG(error) << "!!! Failed !!!";
    return kFALSE;
}
*/
Bool_t R3BTPropagator::PropagateToPlaneBackward(R3BTrackingParticle* particle,
                                                const TVector3& v1,
                                                const TVector3& v2,
                                                const TVector3& v3)
{
    TVector3 norm = ((v3 - v1).Cross(v2 - v1)).Unit();

    // Check if particle is already on plane
    if (TMath::Abs((particle->GetPosition() - v1).Dot(norm)) < 1e-6)
    {
        return kTRUE;
    }

    TVector3 intersect;
    Bool_t crossed;
    Double_t step = 0.;
    Bool_t result;
    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), fPlane2[0], -fNorm2, intersect);
    if (crossed)
    {
        LOG(debug2) << "Starting downstream of magnetic field boundaries...";
        if ((v1 - particle->GetPosition()).Mag() < (fPlane2[0] - particle->GetPosition()).Mag())
        {
            LOG(debug2) << "Propagating to end-plane and stop.";
            crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
            // fVis = kFALSE;
            if (fVis)
            {
                TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
                l1->Draw();
            }
            step = (intersect - particle->GetPosition()).Mag();
            particle->SetPosition(intersect);
            particle->AddStep(step);
            return kTRUE;
        }
        LOG(debug2) << "Propagating to exit of magnetic field.";
        fVis = kFALSE;
        if (fVis)
        {
            TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
            l1->Draw();
        }
        step = (intersect - particle->GetPosition()).Mag();
        particle->SetPosition(intersect);
        particle->AddStep(step);
        LOG(debug2) << intersect.X() << " " << intersect.Y() << " " << intersect.Z();
    }

    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), fPlane1[0], -fNorm1, intersect);
    if (crossed)
    {
        LOG(debug2) << "Propagating inside of field using RK4...";
        TVector3 tpos;
        if ((v1 - particle->GetPosition()).Mag() < (fPlane1[0] - particle->GetPosition()).Mag())
        {
            LOG(debug2) << "Propagating to end-plane using RK4 and stop.";
            tpos = particle->GetPosition();
            result = PropagateToPlaneRK(particle, v1, v3, v2);
            fVis = kFALSE;
            if (fVis)
            {
                TLine* l1 = new TLine(-tpos.X(), tpos.Z(), -particle->GetX(), particle->GetZ());
                l1->Draw();
            }
            LOG(debug2) << particle->GetX() << ", " << particle->GetY() << ", " << particle->GetZ();
            return result;
        }
        LOG(debug2) << "Propagating to entrance of magnetic field.";
        tpos = particle->GetPosition();
        result = PropagateToPlaneRK(particle, fPlane1[0], fPlane1[2], fPlane1[1]);
        fVis = kFALSE;
        if (fVis)
        {
            TLine* l1 = new TLine(-tpos.X(), tpos.Z(), -particle->GetX(), particle->GetZ());
            l1->Draw();
        }
        LOG(debug2) << particle->GetX() << ", " << particle->GetY() << ", " << particle->GetZ();
        if (!result)
        {
            return result;
        }
    }

    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
    if (crossed)
    {
        LOG(debug2) << "Propagating to end plane. Finish.";
        fVis = kFALSE;
        if (fVis)
        {
            TLine* l1 = new TLine(-particle->GetX(), particle->GetZ(), -intersect.X(), intersect.Z());
            l1->Draw();
        }
        step = (intersect - particle->GetPosition()).Mag();
        particle->SetPosition(intersect);
        particle->AddStep(step);
        LOG(debug2);
        return kTRUE;
    }

    // LOG(error) << "!!! Failed !!!";
    return kFALSE;
}

Bool_t R3BTPropagator::PropagateToPlaneRK(R3BTrackingParticle* particle,
                                          const TVector3& v1,
                                          const TVector3& v2,
                                          const TVector3& v3)
{
    Int_t nStep = 0;

    Double_t vecRKIn[7];
    Double_t vecTemp[7];
    Double_t vecOut[7];
    TVector3 dist;
    Double_t diff;
    TVector3 intersect;
    Bool_t crossed;
    Double_t distance;

    particle->GetPosition(vecRKIn);

    // cout<<"In R3BTPropagator::PropagateToPlaneRK: "<<endl;
    // particle->GetPosition().Print();

    particle->GetCosines(&vecRKIn[3]);

    TVector3 norm = ((v2 - v1).Cross(v3 - v1)).Unit();

    // dist = particle->GetPosition() - v1;
    // diff = dist.Dot(norm);
    crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
    dist = (particle->GetPosition() - intersect);
    diff = dist.Mag();
    // cout << "diff: " << diff << endl;

    // Double_t step = TMath::Abs(diff) / 400.;// / 100.;
    Double_t step = 0.1;
    Double_t length = 0.;
    Double_t res = 10000.;
    Double_t res_old = 10000.;

    // cout << "Step: " << step << endl;

    while (kTRUE)
    {

        // cout<<"TEST0: "<<particle->GetPosition().Z()<<", "<<step<<"; "<<diff<<endl;

        if (step > diff)
            step = diff;
        length = fFairProp->OneStepRungeKutta(particle->GetCharge(), step, vecRKIn, vecOut);
        // length = length / 100.;
        // cout << "Length: " << length << endl;

        for (Int_t ii = 0; ii < 7; ii++)
        {
            vecRKIn[ii] = vecOut[ii];
        }

        particle->SetPosition(vecOut);
        particle->SetCosines(&vecOut[3]);
        particle->AddStep(length);

        //   particle->GetPosition().Print();

        nStep += 1;

        // dist = particle->GetPosition() - v1;
        // distance = (TVector3(dist.X() * norm.X(), dist.Y() * norm.Y(), dist.Z() * norm.Z())).Mag();

        crossed = LineIntersectPlane(particle->GetPosition(), particle->GetMomentum(), v1, norm, intersect);
        // cout<<"In propagatRK intersect/v1/particle: "<<endl;
        // intersect.Print();
        // v1.Print();
        // particle->GetPosition().Print();

        dist = (particle->GetPosition() - intersect);
        distance = dist.Mag();
        // cout<<"Distance: "<<distance<<endl;
        // dist.Print();

        res = TMath::Abs(distance); // / diff);

        if (res < 0.1 || res > res_old)
        {

            // cout << "break############" << endl;
            // cout << "Res: " << res << "  " << res_old << "  " << nStep << "  " << distance << endl;
            // cout<<"Intersection point: "<<intersect.X()<<", "<<intersect.Y()<<", "<<intersect.Z()<<endl;

            break;
        }
        else
        {
            res_old = res;

            // step = 1.;

            TVector3 pos = particle->GetPosition();
            Double_t field = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetBy(pos.X(), pos.Y(), pos.Z());

            // particle->GetPosition().Print();
            // cout<<"By: "<<field<<endl;

            step = 20. / TMath::Abs(field); // 10.
            // if(step < 0.1) step = 0.1;

            //  cout << "x: " << pos.X() << " y: " << pos.Y() << " z: " << pos.Z() <<
            // " Field:" << field << " step from Bfield: " << 10./ TMath::Abs(field) << endl;

            //    cout <<  pos.X()<<" " <<  pos.Y()<<" " << pos.Z()<<" " << field <<  endl;

            if (step > TMath::Abs(distance))
            {
                step = TMath::Abs(distance) / 2.;
            }
            if (step > 50.) // 20
            {
                step = 50.; // 20
            }
            /*
            if (step < 0.2) //0.1
            {
                step = 0.2; //0.1
            }
            * */
            if (TMath::Abs(distance) < 0.25) // 0.1
            {
                step = TMath::Abs(distance); // 0.1
            }
            //   step = 0.1;

            //  cout<<"Step final: "<<step<<endl;

            //  cout<<"TEST1: "<<particle->GetPosition().Z()<<", "<<step<<"; "<<diff<<endl;
        }

        if (nStep > 100000)
        {
            return kFALSE;
        }
    }

    return kTRUE;
}

Bool_t R3BTPropagator::LineIntersectPlane(const TVector3& pos,
                                          const TVector3& mom,
                                          const TVector3& v1,
                                          const TVector3& normal,
                                          TVector3& intersect)
{
    TVector3 delta = v1 - pos;
    Double_t pn = mom.Dot(normal);
    if (pn == 0)
    {
        return kFALSE;
    }
    Double_t t = delta.Dot(normal) / pn;
    if (t < 0)
    {
        return kFALSE;
    }
    else
    {
        intersect = pos + mom * t;
    }
    return kTRUE;
}

ClassImp(R3BTPropagator)
