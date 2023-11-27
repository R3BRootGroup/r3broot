##############################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2019-2024 Members of R3B Collaboration                     #
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
    sim/R3BTra.cxx
    sim/R3BTraHitFinder.cxx
    pars/R3BGeoTra.cxx
    pars/R3BGeoTraPar.cxx
    pars/R3BAmsMappingPar.cxx
    pars/R3BAmsStripCalPar.cxx
    pars/R3BSsdContFact.cxx
    digi/R3BTraDigiPar.cxx
    digi/R3BTraDigitizer.cxx
    digi/R3BTargetDigiPar.cxx
    digi/R3BTargetDigitizer.cxx
    digi/R3BAmsDigitizer.cxx
    digi/R3BTra2pDigiPar.cxx
    digi/R3BTra2pDigitizer.cxx
    digi/R3BTraFraDigiPar.cxx
    digi/R3BTraFraDigitizer.cxx
    digi/R3BTarget2pDigiPar.cxx
    digi/R3BTarget2pDigitizer.cxx
    calibration/R3BAmsMapped2StripCalPar.cxx
    calibration/R3BAmsMapped2StripCal.cxx
    calibration/R3BAmsStripCal2Hit.cxx
    online/R3BAmsOnlineSpectra.cxx
    online/R3BAmsCalifaCorrelatedOnlineSpectra.cxx
    online/R3BFootOnlineSpectra.cxx
    pars/R3BFootMappingPar.cxx
    pars/R3BFootCalPar.cxx
    calibration/R3BFootMapped2StripCal.cxx
    calibration/R3BFootStripCal2Hit.cxx)

set(HEADERS
    sim/R3BTra.h
    sim/R3BTraHitFinder.h
    pars/R3BGeoTra.h
    pars/R3BGeoTraPar.h
    pars/R3BAmsMappingPar.h
    pars/R3BAmsStripCalPar.h
    pars/R3BSsdContFact.h
    digi/R3BTraDigiPar.h
    digi/R3BTraDigitizer.h
    digi/R3BTargetDigiPar.h
    digi/R3BTargetDigitizer.h
    digi/R3BAmsDigitizer.h
    digi/R3BTra2pDigiPar.h
    digi/R3BTra2pDigitizer.h
    digi/R3BTraFraDigiPar.h
    digi/R3BTraFraDigitizer.h
    digi/R3BTarget2pDigiPar.h
    digi/R3BTarget2pDigitizer.h
    calibration/R3BAmsMapped2StripCalPar.h
    calibration/R3BAmsMapped2StripCal.h
    calibration/R3BAmsStripCal2Hit.h
    online/R3BAmsOnlineSpectra.h
    online/R3BAmsCalifaCorrelatedOnlineSpectra.h
    online/R3BFootOnlineSpectra.h
    pars/R3BFootMappingPar.h
    pars/R3BFootCalPar.h
    calibration/R3BFootMapped2StripCal.h
    calibration/R3BFootStripCal2Hit.h)

add_library_with_dictionary(
    LIBNAME
    R3BSsd
    LINKDEF
    SsdLinkDef.h
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
    R3BTracking R3BTCal)

add_subdirectory(test)
