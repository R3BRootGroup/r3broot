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

# Create a library called "R3BBunchedFiber" which includes the source files
# given in the array. The extension is already found. Any number of
# sources could be listed here.

set(SRCS
    sim/R3BFiber.cxx
    digi/R3BFiberDigitizer.cxx
    bunched/R3BBunchedFiberMapped2CalPar.cxx
    bunched/R3BBunchedFiberMapped2Cal.cxx
    bunched/R3BBunchedFiberCal2Hit.cxx
    bunched/R3BBunchedFiberCal2Hit_s494.cxx
    bunched/R3BBunchedFiberSPMTTrigMapped2CalPar.cxx
    bunched/R3BBunchedFiberSPMTTrigMapped2Cal.cxx
    mapmt/R3BFiberMAPMTMapped2Cal.cxx
    mapmt/R3BFiberMAPMTCal2Hit.cxx
    mapmt/R3BFiberMAPMTMapped2CalPar.cxx
    pars/R3BFiberContFact.cxx
    pars/R3BBunchedFiberHitModulePar.cxx
    pars/R3BBunchedFiberHitPar.cxx
    pars/R3BFiberMAPMTHitPar.cxx
    pars/R3BFiberMAPMTHitModulePar.cxx
    pars/R3BFiberMappingPar.cxx
    online/R3BFiberMAPMTOnlineSpectra.cxx
    online/R3BFiberMAPMTCorrelationOnlineSpectra.cxx)

set(HEADERS
    sim/R3BFiber.h
    digi/R3BFiberDigitizer.h
    bunched/R3BBunchedFiberMapped2CalPar.h
    bunched/R3BBunchedFiberMapped2Cal.h
    bunched/R3BBunchedFiberCal2Hit.h
    bunched/R3BBunchedFiberCal2Hit_s494.h
    bunched/R3BBunchedFiberSPMTTrigMapped2CalPar.h
    bunched/R3BBunchedFiberSPMTTrigMapped2Cal.h
    mapmt/R3BFiberMAPMTMapped2Cal.h
    mapmt/R3BFiberMAPMTCal2Hit.h
    mapmt/R3BFiberMAPMTMapped2CalPar.h
    pars/R3BFiberContFact.h
    pars/R3BBunchedFiberHitModulePar.h
    pars/R3BBunchedFiberHitPar.h
    pars/R3BFiberMAPMTHitPar.h
    pars/R3BFiberMAPMTHitModulePar.h
    pars/R3BFiberMappingPar.h
    online/R3BFiberMAPMTOnlineSpectra.h
    online/R3BFiberMAPMTCorrelationOnlineSpectra.h)

add_library_with_dictionary(
    LIBNAME
    R3BBunchedFiber
    LINKDEF
    FiberLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    bunched
    digi
    mapmt
    online
    pars
    sim
    DEPENDENCIES
    R3BTCal
    R3BTracking)
