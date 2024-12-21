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

#include "IsElastic.h"
#include "TVector3.h"

// Code extracted from R3BNeutronTracker2D. It is absolutely atrocious. Don't look at it. Nothing too see here. Go away.

bool Neuland::IsElastic(const R3BNeulandCluster* cl1, const R3BNeulandCluster* cl2)
{
    const R3BNeulandCluster* c1;
    const R3BNeulandCluster* c2;

    // if (cl1->GetT() < cl2->GetT())
    // {
    c1 = cl1;
    c2 = cl2;
    /*}
    else
    {
        c1 = cl2;
        c2 = cl1;
    }*/

    Double_t c = 29.9792458;
    Double_t dio = 10.6; // 3 times half the diogonal of a paddle
    Double_t amu = 931.494028;
    // Double_t mNeutron = 1.0086649156 * amu;

    Double_t x0 = 0.;
    Double_t y0 = 0.;
    Double_t z0 = 0.;
    Double_t t0 = 0.;

    //  Double_t v1xmin,v1ymin,v1zmin,v1xmax,v1ymax,v1zmax;
    Double_t v4xmin, v4ymin, v4zmin; //,v4xmax,v4ymax,v4zmax;
    Double_t v3xmin, v3ymin, v3zmin; //,v3xmax,v3ymax,v3zmax;
    //  Double_t v6xmin,v6ymin,v6zmin,v6xmax,v6ymax,v6zmax;

    //    search for scattering
    //    elastic scattering of particle 1 (neutron) on a particle 2 (proton) at
    //    rest. Outgoing particles are 3 (scattered neutron) and 4 (recoil
    //    proton). angle of particles after scattering is theta3 and theta4 K is
    //    kinetic energy, p is momentum

    // TVector3 pos1;
    // TVector3 pos2;

    const TVector3 pos1 = c1->GetFirstHit().GetPosition();

    //    incoming particle
    // vector from previous interaction to present interaction
    Double_t v1x = pos1.X() - x0;
    Double_t v1y = pos1.Y() - y0;
    Double_t v1z = pos1.Z() - z0;
    // time difference and distance between previous interaction to present
    // interaction
    Double_t dt = c1->GetFirstHit().GetT() - t0;
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
    const TVector3 pos2 = c1->GetLastHit().GetPosition();
    Double_t v4x = (pos2 - pos1).X();
    Double_t v4y = (pos2 - pos1).Y();
    Double_t v4z = (pos2 - pos1).Z();
    // vector perpendicular to scattering plane
    Double_t v5x = v1y * v4z - v1z * v4y;
    Double_t v5y = v1z * v4x - v1x * v4z;
    Double_t v5z = v1x * v4y - v1y * v4x;

    Double_t tempAngle;
    dt = c1->GetLastHit().GetT() - c1->GetFirstHit().GetT();
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
    const TVector3 pos3 = c2->GetFirstHit().GetPosition();
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

    dt = c2->GetLastHit().GetT() - c2->GetFirstHit().GetT();
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

    Double_t Ma, Mb, Mc, Md, Ka, Thc, Ei, Pa, AA, BB, a, b, cc, Pc1 /*, Pc2*/;
    Double_t Pd1, Thd1, /*Ec1,*/ Ed1, /*Kc1,*/ Kd1 /*, Qsqr1*/;
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
    // Pc2 = -b / (2. * a) - sqrt((b * b) / (4. * a * a) - (cc / a));
    Pd1 = sqrt(Pc1 * Pc1 + Pa * Pa - 2. * Pc1 * Pa * cos(Thc));
    Thd1 = acos((Pc1 * Pc1 - Pd1 * Pd1 - Pa * Pa) / (-2. * Pd1 * Pa));
    // Ec1 = sqrt(Pc1 * Pc1 + Ma * Ma);
    Ed1 = sqrt(Pd1 * Pd1 + Mb * Mb);
    // Kc1 = Ec1 - Ma;
    Kd1 = Ed1 - Mb;
    // Qsqr1 = (-(Ka - Kc1) * (Ka - Kc1) + (Pa * Pa + Pc1 * Pc1 - 2. * Pa * Pc1 * cos(Thc))) / 197.327 / 197.327;
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

    // calculate minimum and maximum(within errors) momentum of scattered proton
    // and neutron
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
        // Pc2 = -b / (2. * a) - sqrt((b * b) / (4. * a * a) - (cc / a));
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

    // Ec1 = sqrt(p3bmin * p3bmin + Ma * Ma);
    // Kc1 = Ec1 - Ma;
    Ed1 = sqrt(p4bmin * p4bmin + Mb * Mb);
    Kd1 = Ed1 - Mb;
    // Qsqr1 = (-(Ka - Kc1) * (Ka - Kc1) + (Pa * Pa + Pc1 * Pc1 - 2. * Pa * Pc1 * cos(Thc))) / 197.327 / 197.327;
    //  Double_t K3bmin=Kc1;
    //  Double_t E3bmin=Ec1;
    Double_t K4bmin = Kd1;
    //  Double_t E4bmin=Ed1;
    // Ec1 = sqrt(p3bmax * p3bmax + Ma * Ma);
    // Kc1 = Ec1 - Ma;
    Ed1 = sqrt(p4bmax * p4bmax + Mb * Mb);
    Kd1 = Ed1 - Mb;
    // Qsqr1 = (-(Ka - Kc1) * (Ka - Kc1) + (Pa * Pa + Pc1 * Pc1 - 2. * Pa * Pc1 * cos(Thc))) / 197.327 / 197.327;

    //  Double_t K3bmax=Kc1;
    //  Double_t E3bmax=Ec1;
    Double_t K4bmax = Kd1;
    //  Double_t E4bmax=Ed1;

    // decide if Cluster comes from scattered neutron or another neutron is
    // needed!
    Double_t protonEnergy = 8.76839 + 4.13858 * c1->GetE() - 0.00337368 * c1->GetE() * c1->GetE();

    return p3bmax > p3min && p3bmin < p3max && theta4max > theta4Measuredmin && theta4min < theta4Measuredmax &&
           K4bmax > protonEnergy && K4bmin < protonEnergy && theta56 * 180. / 3.14 > 120.;
}
