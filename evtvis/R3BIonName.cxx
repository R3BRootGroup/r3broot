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

// --------------------------------------------------------------------
// -----                    R3BIonName source file                -----
// --------------------------------------------------------------------

#include "R3BIonName.h"

#include <iostream>
using namespace std;

// -----   Default constructor   -------------------------------------------
R3BIonName::R3BIonName()
    : TNamed()
{
}
// -----------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
R3BIonName::~R3BIonName(){};
// -------------------------------------------------------------------------

Bool_t R3BIonName::GetIonName(Int_t element, Char_t* name)
{

    Bool_t flag = true;

    switch (element)
    {
        case 2:
            sprintf(name, "He");
            break;
        case 3:
            sprintf(name, "Li");
            break;
        case 4:
            sprintf(name, "Be");
            break;
        case 5:
            sprintf(name, "B");
            break;
        case 6:
            sprintf(name, "C");
            break;
        case 7:
            sprintf(name, "N");
            break;
        case 8:
            sprintf(name, "O");
            break;
        case 9:
            sprintf(name, "F");
            break;
        case 10:
            sprintf(name, "Ne");
            break;
        case 11:
            sprintf(name, "Na");
            break;
        case 12:
            sprintf(name, "Mg");
            break;
        case 13:
            sprintf(name, "Al");
            break;
        case 14:
            sprintf(name, "Si");
            break;
        case 15:
            sprintf(name, "P");
            break;
        case 16:
            sprintf(name, "S");
            break;
        case 17:
            sprintf(name, "Cl");
            break;
        case 18:
            sprintf(name, "Ar");
            break;
        case 19:
            sprintf(name, "K");
            break;
        case 20:
            sprintf(name, "Ca");
            break;
        case 21:
            sprintf(name, "Sc");
            break;
        case 22:
            sprintf(name, "Ti");
            break;
        case 23:
            sprintf(name, "V");
            break;
        case 24:
            sprintf(name, "Cr");
            break;
        case 25:
            sprintf(name, "Mn");
            break;
        case 26:
            sprintf(name, "Fe");
            break;
        case 27:
            sprintf(name, "Co");
            break;
        case 28:
            sprintf(name, "Ni");
            break;
        case 29:
            sprintf(name, "Cu");
            break;
        case 30:
            sprintf(name, "Zn");
            break;
        case 31:
            sprintf(name, "Ga");
            break;
        case 32:
            sprintf(name, "Ge");
            break;
        case 33:
            sprintf(name, "As");
            break;
        case 34:
            sprintf(name, "Se");
            break;
        case 35:
            sprintf(name, "Br");
            break;
        case 36:
            sprintf(name, "Kr");
            break;
        case 37:
            sprintf(name, "Rb");
            break;
        case 38:
            sprintf(name, "Sr");
            break;
        case 39:
            sprintf(name, "Y");
            break;
        case 40:
            sprintf(name, "Zr");
            break;
        case 41:
            sprintf(name, "Nb");
            break;
        case 42:
            sprintf(name, "Mo");
            break;
        case 43:
            sprintf(name, "Tc");
            break;
        case 44:
            sprintf(name, "Ru");
            break;
        case 45:
            sprintf(name, "Rh");
            break;
        case 46:
            sprintf(name, "Pd");
            break;
        case 47:
            sprintf(name, "Ag");
            break;
        case 48:
            sprintf(name, "Cd");
            break;
        case 49:
            sprintf(name, "In");
            break;
        case 50:
            sprintf(name, "Sn");
            break;
        case 51:
            sprintf(name, "Sb");
            break;
        case 52:
            sprintf(name, "Te");
            break;
        case 53:
            sprintf(name, "I");
            break;
        case 54:
            sprintf(name, "Xe");
            break;
        case 55:
            sprintf(name, "Cs");
            break;
        case 56:
            sprintf(name, "Ba");
            break;
        case 57:
            sprintf(name, "La");
            break;
        case 58:
            sprintf(name, "Ce");
            break;
        case 59:
            sprintf(name, "Pr");
            break;
        case 60:
            sprintf(name, "Nd");
            break;
        case 61:
            sprintf(name, "Pm");
            break;
        case 62:
            sprintf(name, "Sm");
            break;
        case 63:
            sprintf(name, "Eu");
            break;
        case 64:
            sprintf(name, "Gd");
            break;
        case 65:
            sprintf(name, "Tb");
            break;
        case 66:
            sprintf(name, "Dy");
            break;
        case 67:
            sprintf(name, "Ho");
            break;
        case 68:
            sprintf(name, "Er");
            break;
        case 69:
            sprintf(name, "Tm");
            break;
        case 70:
            sprintf(name, "Yb");
            break;
        case 71:
            sprintf(name, "Lu");
            break;
        case 72:
            sprintf(name, "Hf");
            break;
        case 73:
            sprintf(name, "Ta");
            break;
        case 74:
            sprintf(name, "W");
            break;
        case 75:
            sprintf(name, "Re");
            break;
        case 76:
            sprintf(name, "Os");
            break;
        case 77:
            sprintf(name, "Ir");
            break;
        case 78:
            sprintf(name, "Pt");
            break;
        case 79:
            sprintf(name, "Au");
            break;
        case 80:
            sprintf(name, "Hg");
            break;
        case 81:
            sprintf(name, "Tl");
            break;
        case 82:
            sprintf(name, "Pb");
            break;
        case 83:
            sprintf(name, "Bi");
            break;
        case 84:
            sprintf(name, "Po");
            break;
        case 85:
            sprintf(name, "At");
            break;
        case 86:
            sprintf(name, "Rn");
            break;
        case 87:
            sprintf(name, "Fr");
            break;
        case 88:
            sprintf(name, "Ra");
            break;
        case 89:
            sprintf(name, "Ac");
            break;
        case 90:
            sprintf(name, "Th");
            break;
        case 91:
            sprintf(name, "Pa");
            break;
        case 92:
            sprintf(name, "U");
            break;
        case 93:
            sprintf(name, "Np");
            break;
        case 94:
            sprintf(name, "Pu");
            break;
        case 95:
            sprintf(name, "Am");
            break;
        case 96:
            sprintf(name, "Cm");
            break;
        case 97:
            sprintf(name, "Bk");
            break;
        case 98:
            sprintf(name, "Cf");
            break;
        case 99:
            sprintf(name, "Es");
            break;
        case 100:
            sprintf(name, "Fm");
            break;
        case 101:
            sprintf(name, "Md");
            break;
        case 102:
            sprintf(name, "No");
            break;
        case 103:
            sprintf(name, "Lr");
            break;
        case 104:
            sprintf(name, "Rf");
            break;
        case 105:
            sprintf(name, "Db");
            break;
        case 106:
            sprintf(name, "Sg");
            break;
        case 107:
            sprintf(name, "Bh");
            break;
        case 108:
            sprintf(name, "Hs");
            break;
        case 109:
            sprintf(name, "Mt");
            break;
        case 110:
            sprintf(name, "Ds");
            break;
        case 111:
            sprintf(name, "Uuu");
            break;
        case 112:
            sprintf(name, "Uub");
            break;
        case 113:
            sprintf(name, "Uut");
            break;
        case 114:
            sprintf(name, "Uuq");
            break;
        case 115:
            sprintf(name, "Uup");
            break;
        case 116:
            sprintf(name, "Uuh");
            break;
        default:
            flag = false;
            break;
    }

    return flag;
}

ClassImp(R3BIonName)
