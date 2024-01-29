#include "R3BFrsSciTofCalData.h"

R3BFrsSciTofCalData::R3BFrsSciTofCalData()
    : fRank(1)
    , fDetIdSta(1)
    , fDetIdSto(2)
    , fCalPosStaMm(-1000)
    , fCalPosStoMm(-1000)
    , fRawTofNs(0)
    , fCalTofNs(0)
    , fBeta(0)
{
}

R3BFrsSciTofCalData::R3BFrsSciTofCalData(UShort_t rank,
                                         UShort_t detIdSta,
                                         UShort_t detIdSto,
                                         Float_t calPosSta,
                                         Float_t calPosSto,
                                         Double_t rawTof,
                                         Double_t calTof,
                                         Double_t beta)
    : fRank(rank)
    , fDetIdSta(detIdSta)
    , fDetIdSto(detIdSto)
    , fCalPosStaMm(calPosSta)
    , fCalPosStoMm(calPosSto)
    , fRawTofNs(rawTof)
    , fCalTofNs(calTof)
    , fBeta(beta)
{
}
ClassImp(R3BFrsSciTofCalData)
