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

/******************************************************************************
 *****                   Task for fitting a fragment                      *****
 *****     Created Nov 2016 by D. Kresan                                  *****
 *****     Comments to the algorithm implementation from R. Plag          *****
 *****     Implemented ideal fitter, July 2017, D. Kresan                 *****
 ******************************************************************************/

#include "R3BOptimizeMagfieldS494.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "G4NistManager.hh"
#include "R3BGladFieldMap.h"
#include "R3BTGeoPar.h"

#include "TArc.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>

using namespace std;

namespace GlobalChi2
{
    R3BOptimizeMagfieldS494* gThisTask;
} // namespace GlobalChi2

using namespace GlobalChi2;

R3BOptimizeMagfieldS494::R3BOptimizeMagfieldS494(const char* name, Bool_t vis, Int_t verbose)
    : FairTask(name, verbose)
    , fFieldPar(NULL)
{
}

R3BOptimizeMagfieldS494::~R3BOptimizeMagfieldS494() {}

InitStatus R3BOptimizeMagfieldS494::Init()
{
    FairRootManager* man = FairRootManager::Instance();
    gThisTask = this;

    return kSUCCESS;
}

void R3BOptimizeMagfieldS494::Exec(const Option_t*) { cout << "In Exec" << endl; }

double Chi2AllEvents(const double* xx)
{

    Double_t chi2_all = 0;
    Int_t nof = 0;
    cout.precision(7);
    //   cout << "new correction: " << xx[0] << "  " << xx[1] << "  " << xx[2] << "  " << xx[3] << "  " << xx[4] << "  "
    //      << xx[5] << "  " << xx[6] << endl;

    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(xx[0]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(xx[1]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(xx[2]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(xx[3]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(xx[4]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(xx[5]);

    Double_t scale = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetScale();
    //  cout << "Field scale: " << scale << endl;
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(xx[6]);

    for (Int_t iev = 0; iev < 1; iev++)
    {

        Double_t chi2 = gThisTask->Chi2();
        // cout << "Chi2 value: " << chi2 << endl;
        if (chi2 > 0)
        {
            chi2_all += chi2;
            nof++;
        }
    }
    chi2_all = chi2_all / nof;
    // cout << "Chi2 all: " << chi2_all << ", nof: "<<nof<<endl;
    return chi2_all;
}

double R3BOptimizeMagfieldS494::Chi2()
{

    Double_t mChi2 = 0.;
    Double_t Chi2 = 0.;
    Double_t Chi2x = 0.;
    Double_t Chi2y = 0.;
    Double_t Chi2z = 0.;
    Double_t weightX = 1., weightY = 1., weightZ = 1.;
    Double_t sum_weightX = 0., sum_weightY = 0., sum_weightZ = 0.;
    // calculate chi2
    for (Int_t i = 0; i < numPoints; i++)
    // for (Int_t i = 0; i < 1596; i++)
    // for (Int_t i = 1596; i < numPoints; i++)
    {

        Double_t bx = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetBx(d[i][0], d[i][1], d[i][2]);
        Double_t by = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetBy(d[i][0], d[i][1], d[i][2]);
        Double_t bz = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetBz(d[i][0], d[i][1], d[i][2]);

        Chi2x = Chi2x + (bx - d[i][3]) * (bx - d[i][3]);
        weightX = (d[i][4] * 0.001);
        weightY = (d[i][4] * 0.001);
        weightZ = (d[i][5] * 0.001);
        sum_weightX = sum_weightX + weightX * weightX;
        sum_weightY = sum_weightY + weightY * weightY;
        sum_weightZ = sum_weightZ + weightZ * weightZ;

        /*      mChi2 = mChi2 +
                      ((bx - d[i][3]) * (bx - d[i][3]) * weightX * weightX +
                       (by - d[i][4]) * (by - d[i][4]) * weightY * weightY +
                       (bz - d[i][5]) * (bz - d[i][5]) * weightZ * weightZ) /
                          (1.00 * 1.00);*/

        mChi2 =
            mChi2 +
            (sqrt(bx * bx + by * by + bz * bz) - sqrt(d[i][3] * d[i][3] + d[i][4] * d[i][4] + d[i][5] * d[i][5])) *
                (sqrt(bx * bx + by * by + bz * bz) - sqrt(d[i][3] * d[i][3] + d[i][4] * d[i][4] + d[i][5] * d[i][5]));

        Chi2y = Chi2y + (by - d[i][4]) * (by - d[i][4]); //* (weightY * weightY);

        // cout<<"TESTING: "<<i<<"; "<<by - d[i][4]<<", "<<weight<<", "<<(by - d[i][4]) * (by - d[i][4]) / (weight *
        // weight)<<", "<<Chi2y<<endl;

        Chi2z = Chi2z + (bz - d[i][5]) * (bz - d[i][5]);

        if (i < 0)
        {
            // cout << "x: " << i << "  " << bx << "  " << d[i][3] << "  " << Chi2x << endl;
            // cout << "y: " << i << "  " << by << "  " << d[i][4] << "  " << Chi2y << endl;
            // cout << "z: " << i << "  " << bz << "  " << d[i][5] << "  " << Chi2z << endl;
            cout << "Position: " << d[i][0] << "  " << d[i][1] << "  " << d[i][2] << endl;
            cout << "Field theory: " << bx << "  " << by << "  " << bz << endl;
            cout << "Field exp   : " << d[i][3] << "  " << d[i][4] << "  " << d[i][5] << endl;
            // cout << "chi2: " << mChi2 << endl;
        }
    }
    // Chi2y = Chi2y / ((double)numPoints-7.);
    // mChi2 = mChi2 / (sum_weightX+sum_weightY+sum_weightZ) / ((double)numPoints-7.);
    // cout << "Chi2x: " << Chi2x << endl;
    // cout << "Chi2y: " << Chi2y << endl;
    // cout << "Chi2z: " << Chi2z << endl;
    // cout << "Chi2: " << Chi2 << endl;
    // cout << "mChi2 " << mChi2 << endl;

    // return Chi2y;

    // return Chi2z;

    return mChi2;
}

void R3BOptimizeMagfieldS494::Finish()
{
    // cout << "data" << endl;
    Double_t x, y, z, bx, by, bz, current;

    // Set start values of field
    /*
    // in Glad system for starting parameters
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(-6.77);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1.);
    */

    // in Glad system for optimized parameters
    /*
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0.4352514);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(-0.9668069);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(-6.717881);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(0.1382386);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-0.4046498);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(0.002756169);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(0.9987661);
    */
    /*
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(0);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(0.0419);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(0);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(-0.1464);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1);
    */

    /*
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(0.);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(174.95);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(0.0419);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-14.2);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(-0.1464);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1.);
    */

    /*
        // exp_field_map_2022_2023_2000A_0deg_labsyst.dat
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0.3845166);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(-0.9629298);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(174.6417);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(0.09071250);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-14.11291);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(0.1001246);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(0.9999184);
    */
    /*
        // exp_field_map_2023_2000A_0deg_labsyst.dat
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0.2584010);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(-0.9657198);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(174.6048);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(0.1701205);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-14.35249);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(-0.003491582);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(0.9986935);
    */
    /*
        //exp_field_map_longarm.dat
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0.8409866);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(-1.445724);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(174.5730);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(-0.08488216);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-14.18750);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(-0.5596786);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(1.000468);
    */
    /*
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(0.3342123);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(-0.9594673);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(-6.718458);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(0.1505736);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(-0.3753170);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(0.03437903);
        ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(0.9987568);
    */

    string line;
    // ifstream myfile("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_2022_2023_2000A_0deg_labsyst.dat");
    // ifstream myfile("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_2023_2000A_0deg_labsyst.dat");
    // ifstream myfile("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_longarm.dat");
    // ifstream myfile ("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_shortarm0deg.dat");
    // ifstream myfile
    // ("/u/kelic/glad/field_measurements_2023/GLAD_field_measurements_202307/20220613A01/compare/exp_field_map_2022_2023_2000A_0deg_labsyst_v2.dat");
    // ifstream myfile
    // ("/u/kelic/glad/field_measurements_2023/GLAD_field_measurements_202307/20220613A01/compare/exp_field_map_2022_2000A_0deg_labsyst.dat");

    // with measured GLAD angles (Michael): 0.0419,-14.,-0.1464, and TP (Daniel): 53.25cm
    ifstream myfile("/u/kelic/glad/field_measurements_2023/GLAD_field_measurements_202307/20220613A01/compare/"
                    "exp_field_map_2022_2023_2000A_0deg_labsyst_default.dat");
    // with 1stoptim:   {0.05588722, -1.092624, 174.3953, 0.005748210, -14.07997, -0.7210894, 1.001873}; //
    // chi2= 1.248268 1st optim
    //  ifstream myfile
    //  ("/u/kelic/glad/field_measurements_2023/GLAD_field_measurements_202307/20220613A01/compare/exp_field_map_2022_2023_2000A_0deg_labsyst_1stoptim.dat");

    // ofstream myfile1("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_2022_2023_2000A_0deg_labsyst_out.dat");
    // ofstream myfile1("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_2023_2000A_0deg_labsyst_out.dat");
    // ofstream myfile1("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_longarm_out.dat");
    // ofstream myfile1("/u/mheil/R3BRoot/macros/r3b/tracking/s494/exp_field_map_longarm_out_test.dat");
    ofstream myfile1("/u/kelic/glad/field_measurements_2023/GLAD_field_measurements_202307/20220613A01/compare/"
                     "exp_field_map_2022_2023_2000A_0deg_labsyst_optim_tmp.dat");

    getline(myfile, line);
    myfile1 << line << endl;
    // cout << "Text: " << line << endl;
    getline(myfile, line);
    myfile1 << line << endl;
    // cout << "Text: " << line << endl;

    for (Int_t i = 0; i < numPoints; i++)
    {
        myfile >> current >> d0[i][0] >> d0[i][1] >> d0[i][2] >> d0[i][3] >> d0[i][4] >> d0[i][5];
        d0[i][2] = d0[i][2] + 121.0;
        d0[i][3] = d0[i][3] * 10.;
        d0[i][4] = d0[i][4] * 10.;
        d0[i][5] = d0[i][5] * 10.;

        /*
               if(i > 1595)
               {
               / *	d[i][3] = d[i][3] * 1.0015;
                   d[i][4] = d[i][4] * 1.0015;
                   d[i][5] = d[i][5] * 1.0015;* /
               // 2022 and 2023 data (lab system) fitted separately, all params except scale fixed to deafult values;
           ratio between two scale values
               // used to normalize 2023 data to 2022 data
                   / *d[i][3] = d[i][3] * 1.0063853;
                   d[i][4] = d[i][4] * 1.0063853;
                   d[i][5] = d[i][5] * 1.0063853;* /

                   d[i][3] = d[i][3] * 1.0021;
                   d[i][4] = d[i][4] * 1.0021;
                   d[i][5] = d[i][5] * 1.0021;
               }
        */
    }

    // Place here Minuit part for
    // minimization of Chi2AllEvents()

    // ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
    // ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
    // ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
    // Minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    // Migrad, Simplex, Scan,
    // set tolerance , etc...
    mini->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    mini->SetMaxIterations(100000);     // for GSL
    mini->SetTolerance(1);
    mini->SetPrintLevel(1);
    //	mini->SetPrecision(1.e-15);
    // set minimization quality level for low (0), medium (1) and high (2) quality
    mini->SetStrategy(0);

    // create funciton wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor fm1(&Chi2AllEvents, 7);
    mini->SetFunction(fm1);

    //						 x    y   z       thetax thetay thetaz
    // Start values with Michael's measurements for anglesand Daniel's for turning point:
    Double_t variable_default[7] = { 0., 0., 174.25, 0.0419, -14., -0.1464, 1. };

    //  Double_t variable_default[7] = {0. , 0., 174.95, 0,-14.,0, 1.};
    // 1st optim
    // Double_t variable_default[7] = {0.4672199, -1.070673, 174.5847, -0.1062992, -14.18890, -0.6624220, 1.000436};

    Double_t variable[7];
    Double_t variable_mem[7];
    Double_t variable_mem_global[7];
    Double_t chi2mem = 1e6;
    Double_t chi2mem_global = 1e6;

    // Double_t step[7] = { 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.01 };
    Double_t step[7] = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001 };
    Double_t fdata_scale = 1, fdata_scale_mem = 1, fdata_scale_mem_global = 1;

    Int_t nroundsmax = 500;
    Int_t datascalemax = 100;
    for (Int_t ij = 0; ij < datascalemax; ij++)
    {
        fdata_scale = 1.000 + ((float)ij - 50.) / 10000.;
        for (Int_t ik = 0; ik < numPoints; ik++)
        {
            d[ik][0] = d0[ik][0];
            d[ik][1] = d0[ik][1];
            d[ik][2] = d0[ik][2];
            if (ik > 1595)
            {
                d[ik][3] = d0[ik][3] * fdata_scale;
                d[ik][4] = d0[ik][4] * fdata_scale;
                d[ik][5] = d0[ik][5] * fdata_scale;
            }
            else
            {
                d[ik][3] = d0[ik][3];
                d[ik][4] = d0[ik][4];
                d[ik][5] = d0[ik][5];
            }
        }
        for (Int_t i = 0; i < nroundsmax; i++)
        {

            gRandom->SetSeed(0);
            variable[0] = gRandom->Gaus(variable_default[0], 2.); // x
            variable[1] = gRandom->Gaus(variable_default[1], 2.); // y
            variable[2] = gRandom->Gaus(variable_default[2], 3.); // z
            // variable[3] = gRandom->Gaus(variable_default[3], 0.1);	// aX
            variable[4] = gRandom->Gaus(variable_default[4], 1.); // aY
            // variable[5] = gRandom->Gaus(variable_default[5], 0.1);	// aZ
            variable[6] = gRandom->Gaus(variable_default[6], 0.1); // scale

            // variable[0] = variable_default[0];	// x
            // variable[1] = variable_default[1];	// y
            // variable[2] = variable_default[2];	// z
            variable[3] = variable_default[3]; // aX
            // variable[4] = variable_default[4];	// aY
            variable[5] = variable_default[5]; // aZ
            // variable[6] = variable_default[6];	// scale

            /*
            for(Int_t ip = 0; ip < 7; ip++)
            {
                variable[ip] = variable_default[ip];
            }*/

            cout << "** Step if fitpar: " << i << ", step in data scale: " << ij << endl;
            cout << "Start values in this step: " << endl;
            cout.precision(7);
            cout << "{x, y, z, aX, aY, aZ, scale}:  {" << variable[0] << ", " << variable[1] << ", " << variable[2]
                 << ", " << variable[3] << ", " << variable[4] << ", " << variable[5] << ", " << variable[6] << endl;

            Double_t delta_aX = 0., delta_aZ = 0.;
            delta_aX = variable[3] * 0.1;
            delta_aZ = variable[4] * 0.1;

            // Set the free variables to be minimized!
            mini->SetLimitedVariable(0, "xpos", variable[0], step[0], variable[0] - 200., variable[0] + 200.);
            mini->SetLimitedVariable(1, "ypos", variable[1], step[1], variable[1] - 200., variable[1] + 200.);
            mini->SetLimitedVariable(2, "zpos", variable[2], step[2], variable[2] - 200., variable[2] + 200.);
            mini->SetLimitedVariable(3, "angleX", variable[3], step[3], variable[3] - delta_aX, variable[3] + delta_aX);
            mini->SetLimitedVariable(4, "angleY", variable[4], step[4], variable[4] - 20., variable[4] + 20.);
            mini->SetLimitedVariable(5, "angleZ", variable[5], step[5], variable[5] - delta_aZ, variable[5] + delta_aZ);
            mini->SetLimitedVariable(6, "scale", variable[6], step[6], variable[6] - 0.5, variable[6] + 0.5);

            // mini->FixVariable(0);
            // mini->FixVariable(1);
            // mini->FixVariable(2);
            mini->FixVariable(3);
            // mini->FixVariable(4);
            mini->FixVariable(5);
            // mini->FixVariable(6);

            Int_t status = 0;

            // do the minimization
            mini->Minimize();

            status = mini->Status();
            cout << "fdata_scale: " << fdata_scale << endl;
            cout << "Status:  " << status << endl;

            if (mini->MinValue() < chi2mem)
            {
                for (Int_t ip = 0; ip < 7; ip++)
                {
                    variable_mem[ip] = mini->X()[ip];
                }
                chi2mem = mini->MinValue();
                fdata_scale_mem = fdata_scale;
            }

            mini->Clear();
        } // nroundsmax
        if (chi2mem < chi2mem_global)
        {
            for (Int_t ip = 0; ip < 7; ip++)
            {
                variable_mem_global[ip] = variable_mem[ip];
            }
            chi2mem_global = chi2mem;
            fdata_scale_mem_global = fdata_scale_mem;
        }

    } // end over field_scale

    cout << "**** Final results: " << endl;
    cout.precision(7);
    cout << "{x, y, z, aX, aY, aZ, scale}:  " << endl;
    cout << "{" << variable_mem_global[0] << ", " << variable_mem_global[1] << ", " << variable_mem_global[2] << ", "
         << variable_mem_global[3] << ", " << variable_mem_global[4] << ", " << variable_mem_global[5] << ", "
         << variable_mem_global[6] << "}; // chi2= " << chi2mem_global << ", fdata_scale = " << fdata_scale_mem_global
         << endl;
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionX(variable_mem_global[0]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionY(variable_mem_global[1]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionZ(variable_mem_global[2]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleX(variable_mem_global[3]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleY(variable_mem_global[4]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionAngleZ(variable_mem_global[5]);
    ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->SetTrackerCorrectionScale(variable_mem_global[6]);

    Double_t datascale = 1;
    for (Int_t i = 0; i < numPoints; i++)
    {

        Double_t Bx = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetBx(d[i][0], d[i][1], d[i][2]);
        Double_t By = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetBy(d[i][0], d[i][1], d[i][2]);
        Double_t Bz = ((R3BGladFieldMap*)FairRunAna::Instance()->GetField())->GetBz(d[i][0], d[i][1], d[i][2]);
        if (i > 1595)
        {
            datascale = fdata_scale_mem_global;
        }

        myfile1 << current << "  " << d[i][0] << "  " << d[i][1] << "  " << d[i][2] << "  " << d0[i][3] * datascale
                << "  " << d0[i][4] * datascale << "  " << d0[i][5] * datascale << "  " << Bx << "  " << By << "  "
                << Bz << endl;

        //   cout << d[i][0] << "  " << d[i][1] << "  " << d[i][2] << "  " << d[i][3] << "  " << d[i][4] << "  " <<
        //   d[i][5]
        //      << endl;
    }

    myfile.close();
    myfile1.close();
}

Bool_t R3BOptimizeMagfieldS494::InitPropagator()
{
    FairField* fairField = FairRunAna::Instance()->GetField();
    R3BGladFieldMap* gladField = NULL;
    fFieldPar->printParams();
    if (2 == fFieldPar->GetType())
    {
        gladField = (R3BGladFieldMap*)fairField;
    }
    else
    {
        LOG(error) << "Unsupported type of field.";
        return kFALSE;
    }
    return kTRUE;
}

ClassImp(R3BOptimizeMagfieldS494)
