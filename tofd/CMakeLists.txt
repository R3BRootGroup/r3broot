##############################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2019-2025 Members of R3B Collaboration                     #
#                                                                            #
#             This software is distributed under the terms of the            #
#                 GNU General Public Licence (GPL) version 3,                #
#                    copied verbatim in the file "LICENSE".                  #
#                                                                            #
# In applying this license GSI does not waive the privileges and immunities  #
# granted to it by virtue of its status as an Intergovernmental Organization #
# or submit itself to any jurisdiction.                                      #
##############################################################################

set(SRCS
    sim/R3BTofD.cxx
    digi/R3BTofDDigitizer.cxx
    digi/R3BTofDDigitizerCal.cxx
    pars/R3BTofDContFact.cxx
    pars/R3BTofDHitPar.cxx
    pars/R3BTofDHitModulePar.cxx
    pars/R3BTofDMappingPar.cxx
    pars/R3BTofDMapped2CalPar.cxx
    pars/R3BTofDCal2HitPar.cxx
    calibration/R3BTofDMapped2Cal.cxx
    calibration/R3BTofDCal2Hit.cxx
    online/R3BTofDOnlineSpectra.cxx
    #R3BTofDCal2Histo.cxx
    #R3BTofDCal2HistoPar.cxx
    #R3BTofDChangePar.cxx
    #R3BTofDCal2Hit.cxx
    #R3BGeoTofD.cxx
    #R3BGeoTofDPar.cxx
    #R3BTofDContFact.cxx
    #R3BTofDDigiPar.cxx
)

set(HEADERS
    sim/R3BTofD.h
    digi/R3BTofDDigitizer.h
    digi/R3BTofDDigitizerCal.h
    pars/R3BTofDContFact.h
    pars/R3BTofDHitPar.h
    pars/R3BTofDHitModulePar.h
    pars/R3BTofDMappingPar.h
    pars/R3BTofDMapped2CalPar.h
    pars/R3BTofDCal2HitPar.h
    calibration/R3BTofDMapped2Cal.h
    calibration/R3BTofDCal2Hit.h
    online/R3BTofDOnlineSpectra.h
    #R3BTofDCal2Histo.h
    #R3BTofDCal2HistoPar.h
    #R3BTofDChangePar.h
    #R3BTofDCal2Hit.h
    #R3BGeoTofD.h
    #R3BGeoTofDPar.h
    #R3BTofDContFact.h
    #R3BTofDDigiPar.h
)

add_library_with_dictionary(
    LIBNAME
    R3BTofD
    LINKDEF
    TofDLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    calibration
    digi
    online
    pars
    sim
    DEPENDENCIES
    R3BTracking
    R3BTCal)
