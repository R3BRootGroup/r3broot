#include "R3BFrsSciPosCalData.h"

R3BFrsSciPosCalData::R3BFrsSciPosCalData()
    : fDetector(0)
    , fRawTimeNs(0)
    , fRawTimeNs_wTref(0)
    , fRawPosNs(0)
    , fCalPosMm(0)
{
}

R3BFrsSciPosCalData::R3BFrsSciPosCalData(UShort_t det,
                                         Double_t rawT,
                                         Double_t rawT_wTref,
                                         Float_t rawPos,
                                         Float_t calPos)
    : fDetector(det)
    , fRawTimeNs(rawT)
    , fRawTimeNs_wTref(rawT_wTref)
    , fRawPosNs(rawPos)
    , fCalPosMm(calPos)
{
}
ClassImp(R3BFrsSciPosCalData)
