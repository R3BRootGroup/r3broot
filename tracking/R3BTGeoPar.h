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

#ifndef R3BTGEOPAR_H
#define R3BTGEOPAR_H

#include "FairParGenericSet.h"

using namespace std;

class FairParamList;

/**
 * @author D. Kresan
 * @since November 11, 2016
 */
class R3BTGeoPar : public FairParGenericSet
{

  public:
    /**
     * Standard constructor. Creates instance of this class.
     * @param name a name of container.
     * @param title a title of container.
     * @param context context/purpose for parameters and conditions.
     * @param own class ownership, if flag is kTRUE FairDB has the par.
     */
    R3BTGeoPar(const char* name = "TGeoPar",
               const char* title = "Geometry parameters",
               const char* context = "TestDefaultContext",
               Bool_t own = kTRUE);
    /**
     * Destructor. Cleares the memory used by the object.
     */
    virtual ~R3BTGeoPar(void);

    /**
     * Method to reset the values of parameters.
     */
    virtual void clear(void);

    /**
     * Method to store parameters using FairRuntimeDB.
     * @param list a list of parameters.
     */
    virtual void putParams(FairParamList* list);

    /**
     * Method to retrieve parameters using FairRuntimeDB.
     * @param list a list of parameters.
     * @return kTRUE if successful, else kFALSE.
     */
    Bool_t getParams(FairParamList* list);

    /**
     * Method to print value of parameters to the standard output.
     * Calls printParams() for every module container.
     * Accessible from Context menu in TBrowser.
     */
    void printParams(); // *MENU*

    inline const Double_t& GetPosX() const { return fPosX; }
    inline const Double_t& GetPosY() const { return fPosY; }
    inline const Double_t& GetPosZ() const { return fPosZ; }

    inline const Double_t& GetRotX() const { return fRotX; }
    inline const Double_t& GetRotY() const { return fRotY; }
    inline const Double_t& GetRotZ() const { return fRotZ; }

    inline const Double_t& GetDimX() const { return fDimX; }
    inline const Double_t& GetDimY() const { return fDimY; }
    inline const Double_t& GetDimZ() const { return fDimZ; }

    inline const Double_t& GetSigmaX() const { return fSigmaX; }
    inline const Double_t& GetSigmaY() const { return fSigmaY; }

    inline void SetPosXYZ(const Double_t& posX, const Double_t& posY, const Double_t& posZ)
    {
        fPosX = posX;
        fPosY = posY;
        fPosZ = posZ;
    }
    inline void SetRotXYZ(const Double_t& rotX, const Double_t& rotY, const Double_t& rotZ)
    {
        fRotX = rotX;
        fRotY = rotY;
        fRotZ = rotZ;
    }
    inline void SetDimXYZ(const Double_t& dimX, const Double_t& dimY, const Double_t& dimZ)
    {
        fDimX = dimX;
        fDimY = dimY;
        fDimZ = dimZ;
    }

    inline void SetSigmaXY(const Double_t& sigmaX, const Double_t& sigmaY)
    {
        fSigmaX = sigmaX;
        fSigmaY = sigmaY;
    }

    void SetMaterial(const Double_t& Z, const Double_t& A, const Double_t& density, const Double_t& I);

    inline const Double_t& GetZ() const { return fZ; }
    inline const Double_t& GetA() const { return fA; }
    inline const Double_t& GetDensity() const { return fDensity; }
    inline const Double_t& GetI() const { return fI; }

  private:
    const R3BTGeoPar& operator=(const R3BTGeoPar&); /**< an assignment operator */
    R3BTGeoPar(const R3BTGeoPar&);                  /**< a copy constructor */

    Double_t fPosX;
    Double_t fPosY;
    Double_t fPosZ;

    Double_t fRotX;
    Double_t fRotY;
    Double_t fRotZ;

    Double_t fDimX;
    Double_t fDimY;
    Double_t fDimZ;

    Double_t fSigmaX;
    Double_t fSigmaY;

    Double_t fZ;
    Double_t fA;
    Double_t fDensity;
    Double_t fI;

    ClassDef(R3BTGeoPar, 2);
};

#endif /* !R3BTGEOPAR_H*/
