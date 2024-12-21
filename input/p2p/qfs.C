/* Original code was writen in FORTRAN by Leonid Chulkov l.chulkov@gsi.de
 *
 * The program is simulating inverse kinematics of a quasi-elastic scattering
 * of a proton "i" (being at rest) on a cluster/nucleon "a" inside an incident
 * nucleus A. The kinematical variables of the residual fragment "B" and the
 * two scattered particles "1" (scattered) and "2" (knocked-out) are stored
 * in the output ROOT tree.
 *
 *  A --> i  ==> B + (i -> a) =>  B + 1 + 2
 *
 * Valerii Panin October 2012 (valerii.panin@riken.jp)
 *
 *  Modified to reproduce original Chulkov output
 *  for use in R3Bp2pGenerator (B. Pietras)
 */

#include "./info/headers.hh"
#include "./info/info.hh"
#include "libs.hh"

using namespace std;

void run(TString fname)
{
    //------- Output tree ----------
    TFile file("quasi.root", "RECREATE");
    TTree* tree = new TTree("tree", "Tree with simulated QFS kinematics");
    //        #include "tree.hh"
    //--------------------- Beam parameters -----------------------------
    const double Tkin = ENERGY * A; // Total kinetic energy (MeV) of the projectile
    //
    const double PA = sqrt(Tkin * (Tkin + 2 * MA));         // Total 3-momentum of the beam
    const double EA = sqrt(MA * MA + PA * PA);              // Total energy of the beam
    const double bA = -PA / EA;                             // Beta of the beam
    const double gA = 1 / sqrt(1 - bA * bA);                // Gamma of the beam
    const double S_first = (EA + Mi) * (EA + Mi) - PA * PA; // Invariant mass (A+i) (Mandelstam S)

    double target_depth;

    cout << "\n****** Beam parameters ********" << endl;
    cout << "Total momentum:\t\t" << PA << " MeV" << endl;
    cout << "Total energy:\t\t" << EA << " MeV" << endl;
    cout << "Excitation energy:\t" << Exe << " MeV" << endl;
    cout << "Beta (beam):\t\t" << (-bA) << "\nGamma (beam):\t" << gA << endl;

    cout << "Beta (CM):\t\t" << PA / (EA + Mi) << endl;
    cout << "Processing " << MAX_STORY << " events........\n" << endl;

    // Random number generators of ROOT
    auto gRandom = new TRandom3(); // default
    gRandom->SetSeed(0);           // using computer time
    TRandom3 r1;
    r1.SetSeed(0);

    // Output text file for R3BROOT event generator
    ofstream outfile, outfile_ascii;
    outfile.open(fname);
    outfile_ascii.open("quasi_ascii.txt");
    char tooutfile[300];

    int events = 0;            // generated event counter
    while (events < MAX_STORY) // eventloop
    {
        //------------ Internal momentum of a cluster -------------------
        TVector3 Pa;
        Pa.SetX(r1.Gaus(0., MOM_SIGMA));
        Pa.SetY(r1.Gaus(0., MOM_SIGMA));
        Pa.SetZ(r1.Gaus(0., MOM_SIGMA));

        //------------ Internal momentum of the residual-----------------
        TVector3 PB;
        PB.SetXYZ((-Pa.X()), (-Pa.Y()), (-Pa.Z()));
        // Tree variables for the fragmnent
        auto PBx = PB.X();
        auto PBy = PB.Y();
        auto PBz_rf = PB.Z();

        // Off-shell mass of the bound nucleon from the energy conservation
        // in the virtual dissociation A->B+a
        double rrtt = MA * MA + MB * MB - 2 * MA * sqrt(MB * MB + Pa.Mag2());
        if (rrtt <= 0)
        {
            cout << "\nerror off-shell mass!!\n"; // non-zero and real off-shell mass
            continue;
        }
        double Ma_off = sqrt(rrtt);

        // Total energies of "a" and "B" in the restframe of "A"
        double Ea = sqrt(Ma_off * Ma_off + Pa.Mag2());
        double EB = sqrt(MB * MB + PB.Mag2());

        //------- Lorentz transformations into laboratory system ---------
        std::pair<double, double> lora = Lorentz(gA, bA, Ea, Pa.Z());
        double EaL = lora.first; // cluster energy in lab
        Pa.SetZ(lora.second);    // cluster Pz in lab

        std::pair<double, double> lorB = Lorentz(gA, bA, EB, PB.Z());
        double EBL = lorB.first; // energy of the residual B in lab
        PB.SetZ(lorB.second);    // Pz of the residual B in lab

        //---------- Generating CM scattering process ----------
        double S = Ma_off * Ma_off + Mi * Mi + 2 * Mi * EaL; // Mandelstam invariant
        auto Mandelstam_S = S;                               // filling tree variable
        // Now generate CM scattering kinematics
        cm_values CM = CENMASS(S, Ma_off, Mi, Ma, ISOTROPIC);
        if (!CM.good)
            continue; // non-physical output

        TVector3 P1cm(0., 0., 1.), P2cm(0., 0., 1.);
        double phi_rand = r1.Uniform(-PI, PI);

        P2cm.SetMag(CM.p_clust);
        P2cm.SetTheta(CM.theta_clust);
        P2cm.SetPhi(phi_rand);

        P1cm.SetX(-P2cm.X());
        P1cm.SetY(-P2cm.Y());
        P1cm.SetZ(-P2cm.Z());

        //------- Calculate realtive to the direction of the quasi-particle (cluster) --------
        double beta_cm = -Pa.Mag() / (EaL + Mi);
        double gamma_cm = 1 / sqrt(1 - beta_cm * beta_cm);

        std::pair<double, double> lora1 = Lorentz(gamma_cm, beta_cm, CM.e_scat, P1cm.Z());
        std::pair<double, double> lora2 = Lorentz(gamma_cm, beta_cm, CM.e_clust, P2cm.Z());

        P1cm.SetZ(lora1.second);
        P2cm.SetZ(lora2.second);

        //-------- Rotating back to the beam direction -----------
        TVector3 P1L = DREHUNG(P1cm, Pa);
        TVector3 P2L = DREHUNG(P2cm, Pa);

        //--------- Filling in the ROOTTree variables------------
        auto theta_1 = P1L.Theta();
        auto theta_2 = P2L.Theta();
        auto theta_B = PB.Theta();

        auto phi_1 = P1L.Phi();
        auto phi_2 = P2L.Phi();
        auto phi_B = PB.Phi();

        auto P1x = P1L.X();
        auto P1y = P1L.Y();
        auto P1z = P1L.Z();

        auto P2x = P2L.X();
        auto P2y = P2L.Y();
        auto P2z = P2L.Z();

        auto PBz_lab = PB.Z();

        auto E1 = sqrt(Mi * Mi + P1L.Mag2()) - Mi;
        auto E2 = sqrt(Ma * Ma + P2L.Mag2()) - Ma;
        EB = sqrt(MB * MB + PB.Mag2()) - MB;

        auto th1_cm = CM.theta_scat;
        auto th2_cm = CM.theta_clust;
        auto P1_cm = CM.p_scat;
        auto P2_cm = CM.p_clust;
        auto Moff = Ma_off;
        auto Mandelstam_T = CM.T;
        auto Opang = acos(sin(P1L.Theta()) * sin(P2L.Theta()) * cos(P1L.Phi() - P2L.Phi()) +
                          cos(P1L.Theta()) * cos(P2L.Theta()));

        double df = fabs(phi_1 - phi_2);
        double Dif_phi;

        if (df > 0 && df <= PI)
            Dif_phi = df;
        else
            Dif_phi = 2 * PI - df;

        // if(events%10000==0) cout<< events <<" of "<<MAX_STORY<<" ("<<(float)events/MAX_STORY*100<<"%)"<<endl;

        // the new R3BAsciiGenerator wants this as an input:
        // eventId >> nTracks >> pBeam;
        // iPid  >> iZ >> iA >> px >> py >> pz >> vx >> vy >> vz;
        sprintf(tooutfile,
                "%d\t2\t%0.5E\n2212\t1\t1\t%0.5E\t%0.5E\t%0.5E\t0\t0\t0\n2212\t1\t1\t%0.5E\t%0.5E\t%0.5E\t0\t0\t0\n",
                events,
                PA / 1000,
                P2x / 1000,
                P2y / 1000,
                P2z / 1000,
                P1x / 1000,
                P1y / 1000,
                P1z / 1000);

        outfile_ascii << tooutfile; // create an input for the R3BAsciiGenerator

        Double_t P3F = sqrt(PA * PA + P1L.Mag2() - 2 * PA * P1L.Z());
        Double_t A = MA * MA - 2 * EA * sqrt(Mi * Mi + P1L.Mag2()) + 2 * PA * P1L.Z();
        Double_t B = A + 2 * Mi * Mi;
        Double_t QUATER = B + sqrt(B * B + 4 * Mi * Mi * P3F * P3F - A * A);
        Double_t ARAMISS = sqrt(QUATER) - MA;

        sprintf(tooutfile,
                "\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%"
                "0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\t%0.5E\n",
                P2x,
                P2y,
                P2z,
                P1x,
                P1y,
                P1z,
                PBx,
                PBy,
                PBz_lab,
                MA - Ma - MB,
                E1,
                theta_1,
                phi_1,
                E2,
                theta_2,
                phi_2,
                EB,
                theta_B,
                phi_B,
                ARAMISS,
                Mandelstam_T);

        outfile << tooutfile; // Emulate Chulkov's old format for the R3Bp2pgen input

        tree->Fill();
        events++;

        printf(
            "Processed: \033[32m %d \033[0m of \033[34m %d \033[0m (\033[33m %.2f \033[0m of 100) \033[31m \033[0m \r",
            events,
            MAX_STORY,
            100. * events / (double)(MAX_STORY));
        fflush(stdout);
    }

    outfile.close();
    outfile_ascii.close();
    // tree->Print();
    tree->AutoSave();
    file.Close();

    cout << "\nDone!\n";
    return;
}

//---- Two consecutive rotations
// first around Z on <phi>, then around new X' on <theta> (1=Pcm, 2=Pa in lab)
TVector3 DREHUNG(TVector3 v1, TVector3 v2)
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
double CINEMA(double x, double y, double z)
{
    double lambda = x * x + y * y + z * z - 2 * x * y - 2 * x * z - 2 * y * z;
    return lambda;
}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
cm_values CENMASS(double s, double m2off, double m1, double m2, bool isotropic)
{
    cm_values output;
    output.good = false;
    double X = s;
    double Y = m2off * m2off;
    double Z = m1 * m1;
    double sqrs = sqrt(s);

    // Kinematics before the scattering process
    // (with one off-shell mass)
    double p2_off = sqrt(CINEMA(X, Y, Z)) / 2 / sqrs;
    double p1_off = p2_off;
    // CM energies
    double e1_off = (s + Z - Y) / 2 / sqrs;
    double e2_off = (s + Y - Z) / 2 / sqrs;

    // Now take the real masses (after scattering)
    Y = m2 * m2;
    Z = m1 * m1;
    // And check whether the kinematical function is ok
    // for this specific kinematical case
    double error_CI = CINEMA(X, Y, Z);
    if (error_CI <= 0.)
    {
        // cout << "\nerror!!! Kinematical function is negative!";
        return output;
    }

    // Kinematics after the scattering process
    // (with all real masses)
    double p2 = sqrt(CINEMA(X, Y, Z)) / 2 / sqrs;
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
        t = get_T(s, tmax);
    } // Using parameterized cross sections
    else
    {
        t = gRandom->Uniform(0., 3000000.) * (-1);
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
double momentum_CM(double TLAB, double M1, double M2)
{
    // Particle M2 is assumed to be in rest
    double PLAB = sqrt(TLAB * (TLAB + 2 * M1));      //  Total 3-momentum of an incident particle in Lab
    double ELAB = sqrt(PLAB * PLAB + M1 * M1);       //  Total energy of an incident particle in lab
    double SLAB = M1 * M1 + M2 * M2 + 2 * M2 * ELAB; // Mandelstam invariant S in lab
    double PCM = PLAB * M2 / sqrt(SLAB);             // Momentum of both particles in CM frame
    return PCM;
}

std::pair<double, double> Lorentz(double g, double b, double e, double p)
{
    double eL = g * e - g * b * p;
    double pL = g * p - g * b * e;
    return std::make_pair(eL, pL);
}

// Returns a random value of mandelstam T (in (MeV/c) units)
// distributed according to the parameterized proton-proton
// invarant cross section. Pass as a parameter "sm" the
// Mandelstam variable S (in MeV)
// and the maximum possible momentum transfer
double get_T(double sm, double max)
{
    // TRandom1 rand;
    // rand.SetSeed(0);

    double Tmax = max * 0.000001; // convert to GeV units
    // double Tmin = min*0.000001;

    // double Tmax = -2*pCM*pCM*(1 - cos(PI))*0.000001; //in (GeV/c)
    Double_t rr = gRandom->Uniform(-1., 1.); // to randomize wrt 90 degrees
    double mandels = sm * 0.000001;          // in GeV
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
    return (Trand * 1000000); // returning value in MeV
}

int main()
{
    run("quasi.out");
    return 1;
}
