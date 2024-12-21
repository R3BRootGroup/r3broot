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

/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Oleg Kiselev		O.Kiselev@gsi.de
//             Denis Bertini            <D.Bertini@gsi.de>
//*-- Date: 08/2009
//*-- Last Update: 08/08/09 Denis Bertini
// --------------------------------------------------------------
// Description:
//   Event generator kinematics reader
//
// --------------------------------------------------------------
// Comments:
//          - 05/08/09  <D.Bertini@gsi.de>
//
// --------------------------------------------------------------
//
/////////////////////////////////////////////////////////////////

#include "R3BReadKinematics.h"
#include "TString.h"
#include "TSystem.h"

ClassImp(R3BReadKinematics)

    using namespace std;

R3BReadKinematics::R3BReadKinematics()
{

    TString WorkDir = gSystem->Getenv("VMCWORKDIR");

    //
    // Constructor, all job here
    //

    Double_t a, b, c, d, e, f, j, h, k, l, m, n;
    Double_t Integral = 0.;
    Int_t i;
    Int_t nbins = 0;
    const Int_t MAX = 140;
    Char_t buffer[MAX];

    // zero initial values
    for (i = 0; i < 1000; i++)
    {
        LabAngle[i] = 0.;
        LabEnergy[i] = 0.;
    }

    for (i = 0; i < 5000; i++)
    {
        T_xsec[i] = 0.;
        CrossSection_t[i] = 0.;
        Icross_section[i] = 0.;
    }

    // Read kinematical file
    TString kine1 = WorkDir + "/input/kine/KINR.OUT.12c_p_p_0MeV_400MeV_long.dat";
    ifstream infile(kine1.Data());

    //      if (!infile);
    for (i = 0; i < 17; i++)
    {
        infile.getline(buffer, MAX); // readout header
        //	  cout << "i=" << i << " buffer" << infile.rdbuf() <<  G4endl;
    }
    //	cout << "read file" << i << G4endl;
    i = -1;
    while (!infile.eof())
    {
        i++;

        // readout of data - elastic case
        infile >> a >> b >> c >> d >> e >> f >> j >> LabAngle[i] >> LabEnergy[i] >> h >> k >> l >> m >> n;
        //	cout << "read file a=" << a << G4endl;
        // readout of data - inelastic case
        //	infile >> a >> b >> c >> d >> e >> f >> j >>
        //	LabAngle[i] >> LabEnergy[i] >> h >> k >> l >> m;

        // readout of data - transfer case
        //	infile >> LabAngle[i] >> LabEnergy[i] >> a >> b;
    }

    // Read cross section file

    TString kine2 = WorkDir + "/input/kine/c12_400_long.dat";
    ifstream infile_xsec(kine2.Data());

    i = -1;
    while (!infile_xsec.eof())
    {
        i++;
        infile_xsec >> T_xsec[i] >> CrossSection_t[i]; // readout of data

        //  cout << i << " " << T_xsec[i] << " " << CrossSection_t[i] <<  endl;
        nbins = i;
    }
    // Int_Tegral of cross section
    for (i = 0; i < nbins; i++)
    {
        Integral = Integral + CrossSection_t[i];
    }
    // cout << " Integral " << Integral << endl;

    // Integrated cross section - probability
    Icross_section[0] = CrossSection_t[0] / Integral;
    //  cout << "Cross_section[0] " << CrossSection_t[0]  << endl;
    //  cout << "Icross_section[0] " <<  Icross_section[0] << endl;

    for (i = 1; i < nbins; i++)
    {
        Icross_section[i] = Icross_section[i - 1] + CrossSection_t[i] / Integral;
        //       cout << "i" << i << "Icross_section[i]" <<  Icross_section[i]   << endl;
    }
}

R3BReadKinematics::~R3BReadKinematics()
{
    //
    // Destructor
    //
}
