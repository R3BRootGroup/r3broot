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

#include "R3Bp2pevtGenerator.h"

#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <FairLogger.h>
#include <FairPrimaryGenerator.h>

using namespace std;

R3Bp2pevtGenerator::R3Bp2pevtGenerator()
{
    this->SetValues(0, 0, 0, 0, false, false);
    T_LIMIT = 10000000.;
    fRandom.SetSeed(0);
}

Bool_t R3Bp2pevtGenerator::Init()
{
    fRandom.SetSeed(0);
    return kTRUE;
}

void R3Bp2pevtGenerator::SetBeamEnergy(double E)
{
    ENERGY = E;
    return;
}

void R3Bp2pevtGenerator::SetHeavyNucleus(int A, double M_def, double M_res)
{
    A_BEAM = A;
    double UNIT = 931.494061;
    MA = A * UNIT + M_def;
    MB = (A - 1) * UNIT + M_res;
    return;
}
void R3Bp2pevtGenerator::SetLightNucleus(double ma, double mi)
{
    Ma = ma;
    Mi = mi;
    return;
}

void R3Bp2pevtGenerator::SetExcitation(double exe)
{
    MB += exe;
    return;
}

void R3Bp2pevtGenerator::SetInverse(bool inverse)
{
    INVERSE = inverse;
    return;
}

void R3Bp2pevtGenerator::SetIsotropic(bool isotropic)
{
    ISOTROPIC = isotropic;
    return;
}

void R3Bp2pevtGenerator::SetMomDistrib(double mom)
{
    MOM_SIGMA = mom;
    return;
}

void R3Bp2pevtGenerator::Print()
{
    cout << "***** Print generator values *****" << endl;
    cout << "Energy: \t" << ENERGY << endl;
    cout << "A Beam: \t" << A_BEAM << endl;
    cout << "Mass A: \t" << MA << endl;
    cout << "Mass B: \t" << MB << endl;
    cout << "Mass a: \t" << Ma << endl;
    cout << "Mass i: \t" << Mi << endl;
    cout << "T_Limit: \t" << T_LIMIT << endl;
    cout << "Sigma Fermi: \t" << MOM_SIGMA << endl;
    cout << "WQ Isotropic: \t" << ISOTROPIC << endl;
    cout << "Inverse kinematics: \t" << INVERSE << endl;
}

Bool_t R3Bp2pevtGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
    double E = ENERGY;
    // GenerateEvent(ENERGY,primGen);

    if (E <= 0)
    {
        LOG(error) << "R3Bp2pevtGenerator: E < 0!";
        return kFALSE;
    }
    fRandom.SetSeed(0);
    //  	gRandom = new TRandom3(); //default
    //  	gRandom->SetSeed(0);//using computer time
    //  	TRandom2 r1;
    //  	r1.SetSeed(0);
    R3Bp2pevtGenerator::Print();
    bool evt = false;
    if (INVERSE)
    {
        //-------------- The reaction 11C(p,pn)10C --------------------
        const double Tkin = E * A_BEAM; // Total kinetic energy (MeV) of the projectile nucleus
        //--------------------- Beam parameters -----------------------------
        const double PA = sqrt(Tkin * (Tkin + 2 * MA));         // Total 3-momentum of the beam
        const double EA = sqrt(MA * MA + PA * PA);              // Total energy of the beam
        const double bA = -PA / EA;                             // Beta of the beam
        const double gA = 1 / sqrt(1 - bA * bA);                // Gamma of the beam
        const double S_first = (EA + Mi) * (EA + Mi) - PA * PA; // Invariant mass (Mandelstam S-variable)
        const double sigma = MOM_SIGMA; // Internal momentum spread of the cluster "a" inside "A"

        LOG(debug) << "\n****** Beam parameters ********";
        LOG(debug) << "\nMA:\t" << MA << " MeV";
        LOG(debug) << "\nTotal momentum:\t" << PA << " MeV";
        LOG(debug) << "\nTotal energy:\t" << EA << " MeV";
        LOG(debug) << "\nBeta (beam):\t" << (-bA) << "\nGamma (beam):\t" << gA << "\n\n";

        while (!evt) // eventloop
        {

            //------------ Internal momentum of a cluster -------------------
            TVector3 Pa;
            /*Pa.SetX(input_xy->GetRandom()); ### it is for momentum distribution from carlos calculations
             *input_filename Pa.SetY(input_xy->GetRandom()); Pa.SetZ(input_z->GetRandom()); ### */
            Pa.SetX(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
            Pa.SetY(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
            Pa.SetZ(fRandom.Gaus(0, sigma)); // ### fixed experimental width of mom distribution in info.hh
            //------------ Internal momentum of the residual-----------------
            TVector3 PB;
            PB.SetX(-Pa.X());
            PB.SetY(-Pa.Y());
            PB.SetZ(-Pa.Z());

            // From the energy conservation in the virtual dissociation A->B+a
            double rrtt = MA * MA + MB * MB - 2 * MA * sqrt(MB * MB + Pa.Mag2());

            if (rrtt <= 0)
            {
                // cout<<"\nerror off-shell mass!!";//non-zero and real off-shell mass
                // cout<<"\nP:"<< PBx << "\t" << PBy << "\t" << PBz_rf << "\n";//non-zero and real off-shell mass
                continue;
            }

            // Off-shell mass of the bound cluster
            double Ma_off = sqrt(rrtt);

            // Total energies of "a" and "B" in the restframe of "A"
            double Ea = sqrt(Ma_off * Ma_off + Pa.Mag2());
            double EB = sqrt(MB * MB + PB.Mag2());

            //------- Lorentz transformations into laboratory system ---------
            std::pair<double, double> lora = R3Bp2pevtGenerator::Lorentz(gA, bA, Ea, Pa.Z());
            double EaL = lora.first; // cluster energy in lab
            Pa.SetZ(lora.second);    // cluster Pz in lab

            std::pair<double, double> lorB = R3Bp2pevtGenerator::Lorentz(gA, bA, EB, PB.Z());
            double EBL = lorB.first; // energy of the residual B in lab
            PB.SetZ(lorB.second);    // Pz of the residual B in lab

            //---------- Generating CM scattering process ----------
            double S = Ma_off * Ma_off + Mi * Mi + 2 * Mi * EaL; // Mandelstam invariant

            // Now generate CM scattering kinematics
            cm_values CM = R3Bp2pevtGenerator::CENMASS(S, Ma_off, Mi, Ma, ISOTROPIC);
            if (!CM.good)
                continue; // non-physical output

            TVector3 P1cm(0., 0., 1.), P2cm(0., 0., 1.);
            double phi_rand = fRandom.Uniform(-PI, PI);

            P2cm.SetMag(CM.p_clust);
            P2cm.SetTheta(CM.theta_clust);
            P2cm.SetPhi(phi_rand);

            P1cm.SetX(-P2cm.X());
            P1cm.SetY(-P2cm.Y());
            P1cm.SetZ(-P2cm.Z());

            //------- Calculate realtive to the direction of the quasi-particle (cluster) --------
            double beta_cm = -Pa.Mag() / (EaL + Mi);
            double gamma_cm = 1 / sqrt(1 - beta_cm * beta_cm);

            std::pair<double, double> lora1 = R3Bp2pevtGenerator::Lorentz(gamma_cm, beta_cm, CM.e_scat, P1cm.Z());
            std::pair<double, double> lora2 = R3Bp2pevtGenerator::Lorentz(gamma_cm, beta_cm, CM.e_clust, P2cm.Z());

            P1cm.SetZ(lora1.second);
            P2cm.SetZ(lora2.second);
            //-------- Rotating back to the beam direction -----------
            TVector3 P1L = R3Bp2pevtGenerator::DREHUNG(P1cm, Pa);
            TVector3 P2L = R3Bp2pevtGenerator::DREHUNG(P2cm, Pa);

            evt = true;
            LOG(debug) << "R3Bp2pevtGenerator: Sending p2pevt: P1 : " << P1L.Px() << " , " << P1L.Py() << " , "
                       << P1L.Pz() << "\n P2 : " << P2L.Px() << " , " << P2L.Py() << " , " << P2L.Pz() << " ";
            primGen->AddTrack(2212, P1L.Px() / 1000., P1L.Py() / 1000., P1L.Pz() / 1000., 0, 0, 0);
            primGen->AddTrack(2212, P2L.Px() / 1000., P2L.Py() / 1000., P2L.Pz() / 1000., 0, 0, 0);

            //--------- Filling in the ROOTTree variables------------
        }
    }
    else if (!INVERSE)
    {
        //--------------------- Beam parameters -----------------------------
        const Double_t Tkin = E;                          // Total kinetic energy (MeV) of the projectile
        const Double_t Pi = sqrt(Tkin * (Tkin + 2 * Mi)); // Total 3-momentum of the beam
        const Double_t Ei = sqrt(Mi * Mi + Pi * Pi);      // Total energy of the beam
        const Double_t bi = -Pi / Ei;                     // Beta of the beam
        const Double_t gi = 1 / sqrt(1 - bi * bi);        // Gamma of the beam

        std::cout << "\n****** Beam parameters ********";
        std::cout << "\nTotal momentum:\t" << Pi << " MeV/c";
        std::cout << "\nTotal energy:\t" << Ei << " MeV";
        std::cout << "\nBeta (beam):\t" << (-bi) << "\nGamma (beam):\t" << gi << "\n\n";

        TLorentzVector LVa; // Lorentz vector of the target proton
        TLorentzVector LVi; // Lorentz vector of the beam
        LVi.SetPxPyPzE(0.0, 0.0, Pi, Ei);
        TVector3 Pa;
        TVector3 P1cm(1e-12, 0.0, 0.0);
        TVector3 P2cm(1e-12, 0.0, 0.0);

        while (!evt) // eventloop
        {
            //------------ Internal momentum of a cluster -------------------
            Pa.SetX(fRandom.Gaus(0., MOM_SIGMA));
            Pa.SetY(fRandom.Gaus(0., MOM_SIGMA));
            Pa.SetZ(fRandom.Gaus(0., MOM_SIGMA));

            //------------ Internal momentum of the residual-----------------
            TVector3 PB;
            PB.SetXYZ((-Pa.X()), (-Pa.Y()), (-Pa.Z()));

            double rrtt = MA * MA + MB * MB - 2 * MA * sqrt(MB * MB + Pa.Mag2());
            if (rrtt <= 0)
            {
                // cout<<"\nerror off-shell mass!!\n";//non-zero and real off-shell mass
                continue;
            }
            double Ma_off = sqrt(rrtt);

            // Total energies of "a" and "B" in the restframe of "A"
            double Ea = sqrt(Ma_off * Ma_off + Pa.Mag2());
            double EB = sqrt(MB * MB + PB.Mag2());
            LVa.SetPxPyPzE(Pa.X(), Pa.Y(), Pa.Z(), Ea);

            TLorentzVector LVstart = LVa + LVi;
            double S = LVstart.Mag2(); // Mandelstam invariant

            cm_values CM = R3Bp2pevtGenerator::CENMASS(S, Ma_off, Mi, Ma, ISOTROPIC);
            if (!CM.good)
            {
                // cout << "non physical" << endl;
                continue; // non-physical output
            }
            double phi_rand = fRandom.Uniform(-PI, PI);

            P2cm.SetMag(CM.p_clust);
            P2cm.SetTheta(CM.theta_clust);
            P2cm.SetPhi(phi_rand);

            P1cm.SetX(-P2cm.X());
            P1cm.SetY(-P2cm.Y());
            P1cm.SetZ(-P2cm.Z());

            //------- Calculate realtive to the direction of the CM motion ----------
            double beta_cm = 0.00000 - LVstart.Beta();
            double gamma_cm = 1 / sqrt(1 - beta_cm * beta_cm);

            std::pair<double, double> lora1 = R3Bp2pevtGenerator::Lorentz(gamma_cm, beta_cm, CM.e_scat, P1cm.Z());
            std::pair<double, double> lora2 = R3Bp2pevtGenerator::Lorentz(gamma_cm, beta_cm, CM.e_clust, P2cm.Z());

            P1cm.SetZ(lora1.second);
            P2cm.SetZ(lora2.second);

            //-------- Rotating back to the beam direction -----------
            // TVector3 P1L = DREHUNG(P1cm,(LVstart.Vect()));
            // TVector3 P2L = DREHUNG(P2cm,(LVstart.Vect()));

            TVector3 direction = LVstart.Vect();
            direction = direction.Unit();
            P1cm.RotateUz(direction);
            P2cm.RotateUz(direction);

            TVector3 P1L = P1cm;
            TVector3 P2L = P2cm;

            // TVector3 P1L = DREHUNG(P1cm,(LVstart.Vect()));
            // TVector3 P2L = DREHUNG(P2cm,(LVstart.Vect()));
            evt = true;

            LOG(debug) << "R3Bp2pevtGenerator: Sending p2pevt: P1 : " << P1L.Px() << " , " << P1L.Py() << " , "
                       << P1L.Pz() << "\n P2 : " << P2L.Px() << " , " << P2L.Py() << " , " << P2L.Pz() << " ";
            primGen->AddTrack(2212, P1cm.Px() / 1000., P1cm.Py() / 1000., P1cm.Pz() / 1000., 0, 0, 0);
            primGen->AddTrack(2212, P2cm.Px() / 1000., P2cm.Py() / 1000., P2cm.Pz() / 1000., 0, 0, 0);
        }
    }
    return kTRUE;
}

void R3Bp2pevtGenerator::SetValues(double E, int A, double MOM, double exe, bool invert, bool iso)
{
    ENERGY = E;
    MOM_SIGMA = MOM;
    A_BEAM = A;
    double UNIT = 931.494061;
    MA = A * UNIT;
    MB = (A - 1) * UNIT + exe;
    Ma = 938.727;
    Mi = 938.727;
    T_LIMIT = 10000000.;
    INVERSE = invert;
    ISOTROPIC = iso;
    return;
}

TVector3 R3Bp2pevtGenerator::DREHUNG(TVector3 v1, TVector3 v2)
{
    double CT = v2.Z() / v2.Mag(); // cos(theta) of v2 wrt. Z-axis
    double ST = sqrt(1 - CT * CT); // sin(theta)
    double CF = v2.X() / v2.Mag() / ST;
    double SF = v2.Y() / v2.Mag() / ST;

    TVector3 v3;
    double _v3x = v1.X() * CT * CF - v1.Y() * SF + v1.Z() * ST * CF;
    double _v3y = v1.X() * CT * SF + v1.Y() * CF + v1.Z() * ST * SF;
    double _v3z = -v1.X() * ST + v1.Z() * CT;
    v3.SetXYZ(_v3x, _v3y, _v3z);
    return v3;
}

// Kinematical function
double R3Bp2pevtGenerator::CINEMA(double x, double y, double z)
{
    double lambda = x * x + y * y + z * z - 2 * x * y - 2 * x * z - 2 * y * z;
    return lambda;
}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
cm_values R3Bp2pevtGenerator::CENMASS(double s, double m2off, double m1, double m2, bool isotropic)
{
    cm_values output;
    output.good = false;
    double X = s;
    double Y = m2off * m2off;
    double Z = m1 * m1;
    double sqrs = sqrt(s);

    // Kinematics before the scattering process
    // (with one off-shell mass)
    double p2_off = sqrt(R3Bp2pevtGenerator::CINEMA(X, Y, Z)) / 2 / sqrs;
    double p1_off = p2_off;
    // CM energies
    double e1_off = (s + Z - Y) / 2 / sqrs;
    double e2_off = (s + Y - Z) / 2 / sqrs;

    // Now take the real masses (after scattering)
    Y = m2 * m2;
    Z = m1 * m1;
    // And check whether the kinematical function is ok
    // for this specific kinematical case
    double error_CI = R3Bp2pevtGenerator::CINEMA(X, Y, Z);
    if (error_CI <= 0.)
    {
        // cout << "\nerror!!! Kinematical function is negative!";
        return output;
    }

    // Kinematics after the scattering process
    // (with all real masses)
    double p2 = sqrt(R3Bp2pevtGenerator::CINEMA(X, Y, Z)) / 2 / sqrs;
    double p1 = p2;
    double e1 = (s + Z - Y) / 2 / sqrs;
    double e2 = (s + Y - Z) / 2 / sqrs;

    // Let's consider momentum transfer <t> from the
    // target particle 1 to the cluster 2
    double tmax = 2 * (m1 * m1 - e1_off * e1 - p1_off * p1); // COSINE=(-1)
    double tmin = 2 * (m1 * m1 - e1_off * e1 + p1_off * p1); // COSINE=(1)
    // cout << "\n\n Tmax = " << tmax;
    // cout << "\n Tmin = " << tmin;
    // cout << "\n Mandels = " << X;

    double t;
    // Generate random momentum transfer for this kinematical case
    if (!isotropic)
    {
        t = R3Bp2pevtGenerator::get_T(s, tmax);
    } // Using parameterized cross sections
    else
    {
        t = fRandom.Uniform(0., T_LIMIT) * (-1);
    } // Isotropic scattering

    // double COSINE = (t - m2off*m2off - m2*m2 + 2*e2_off*e2)/(2*p2_off*p2);
    double COSINE = (t - 2 * m1 * m1 + 2 * e1_off * e1) / (2 * p1_off * p1);
    if (fabs(COSINE) >= 1)
    { // momentum transfer out of range
        // cout << "\nerror! Scattering cosine is larger than 1";
        return output;
    }

    // CM scattering angles
    double theta1 = acos(COSINE);
    double theta2 = PI - theta1;

    output.e_clust = e2;
    output.p_clust = p2;
    output.theta_clust = theta2;

    output.e_scat = e1;
    output.p_scat = p1;
    output.theta_scat = theta1;

    output.T = t;
    output.good = true;

    return output;
}

// Calculate 3-momentum in CM system of two particles M1 and M2
// when M1 has kinetic energy TLAB and M2 is at rest
double R3Bp2pevtGenerator::momentum_CM(double TLAB, double M1, double M2)
{
    // Particle M2 is assumed to be in rest
    double PLAB = sqrt(TLAB * (TLAB + 2 * M1));      //  Total 3-momentum of an incident particle in Lab
    double ELAB = sqrt(PLAB * PLAB + M1 * M1);       //  Total energy of an incident particle in lab
    double SLAB = M1 * M1 + M2 * M2 + 2 * M2 * ELAB; // Mandelstam invariant S in lab
    double PCM = PLAB * M2 / sqrt(SLAB);             // Momentum of both particles in CM frame
    return PCM;
}

std::pair<double, double> R3Bp2pevtGenerator::Lorentz(double g, double b, double e, double p)
{
    double eL = g * e - g * b * p;
    double pL = g * p - g * b * e;
    return std::make_pair(eL, pL);
}

// Returns a random value of mandelstam T (in (MeV/c)� units)
// distributed according to the parameterized proton-proton
// invarant cross section. Pass as a parameter "sm" the
// Mandelstam variable S (in MeV�)
// and the maximum possible momentum transfer
double R3Bp2pevtGenerator::get_T(double sm, double max)
{
    // TRandom1 rand;
    // rand.SetSeed(0);

    double Tmax = max * 0.000001; // convert to GeV� units
    // double Tmin = min*0.000001;

    // double Tmax = -2*pCM*pCM*(1 - cos(PI))*0.000001; //in (GeV/c)�
    Double_t rr = fRandom.Uniform(-1., 1.); // to randomize wrt 90 degrees
    double mandels = sm * 0.000001;         // in GeV�
    // cout << "\nMandelstam S = " << mandels << "\t Tmax/2 = " << Tmax/2 << "\t Random: " << rr;

    // Probability function from the parameterization
    TF1* foo = new TF1("foo", "[0]*exp(x*[1])*(1+0.02*exp((-6)*x))", Tmax / 2, 0);
    double c = 0.;
    if (mandels <= 4.79)
        c = -3283.75 + 3064.11 * mandels - 1068.44 * mandels * mandels + 164.844 * pow(mandels, 3) -
            9.48152 * pow(mandels, 4);
    else if (mandels > 4.79)
        c = -776.822 + 586.016 * mandels - 175.347 * mandels * mandels + 26.1823 * pow(mandels, 3) -
            1.94889 * pow(mandels, 4) + 0.0578352 * pow(mandels, 5);

    foo->FixParameter(0, 25.); // normalization constant (could be anything)
    foo->FixParameter(1, c);

    double Trand = foo->GetRandom(Tmax / 2, 0.); // from 90 to 0 degrees
    if (rr > 0)
        Trand = Tmax - Trand; // symmetrization relative to 90 degrees
    // cout << "\n Tmax/2 = " << Tmax/2 *  1000000;
    // cout << "\n Random T = " << Trand*1000000;
    delete foo;
    return (Trand * 1000000); // returning value in MeV�
}

ClassImp(R3Bp2pevtGenerator)
