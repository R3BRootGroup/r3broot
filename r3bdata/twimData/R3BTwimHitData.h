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

// -------------------------------------------------------------------
// -----           R3BTwimMappedData header file                 -----
// -----       Created 13/10/19 by J.L. Rodriguez-Sanchez        -----
// -------------------------------------------------------------------

#pragma once

#include <TObject.h>

class R3BTwimHitData : public TObject
{

  public:
    /** Default constructor **/
    R3BTwimHitData();

    /** Constructor with arguments
     *@param secID    Section ID
     *@param theta    Angle theta in [mrad]
     *@param z        Atomic number Z in charge units
     *@param x        Position of a reference anode
     *@param offset   Offset of the fission fragment trajectory fit[mm]
     *@param vertex   vertex of the fission fragments trajectory (fission position)
     **/
    R3BTwimHitData(UInt_t secID,
                   Double_t theta,
                   Double_t z,
                   Double_t x = 0.,
                   Double_t Offset = 0.,
                   Double_t ene = 0.,
                   Double_t VertexX = 0.,
                   Double_t VertexZ = 0.);

    /** Destructor **/
    virtual ~R3BTwimHitData() {}

    /** Accessors **/
    inline const UInt_t& GetSecID() const { return fSecID; }
    inline const Double_t& GetTheta() const { return fTheta; }
    inline const Double_t& GetZcharge() const { return fZ; }
    inline const Double_t& GetEave() const { return fE; }
    inline const Double_t& GetX() const { return fX; }
    inline const Double_t& GetOffset() const { return fOffset; }
    inline const Double_t& GetVertexX() const { return fVertexX; }
    inline const Double_t& GetVertexZ() const { return fVertexZ; }

    /** Modifiers **/
    void SetSecID(UInt_t id) { fSecID = id; };
    void SetTheta(Double_t theta) { fTheta = theta; };
    void SetZcharge(Double_t z) { fZ = z; };
    void SetEave(Double_t ene) { fE = ene; };
    void SetX(Double_t x) { fX = x; };
    void SetOffset(Double_t offset) { fOffset = offset; };
    void SetVertexX(Double_t vertexX) { fVertexX = vertexX; };
    void SetVertexZ(Double_t vertexZ) { fVertexZ = vertexZ; };

  protected:
    UInt_t fSecID;
    Double_t fTheta, fZ, fX, fOffset, fE, fVertexX, fVertexZ;

  public:
    ClassDef(R3BTwimHitData, 1)
};
