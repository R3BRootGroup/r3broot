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

# Create a library called "libR3BTofi" which includes the source files given
# in the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    #Put here your sourcefiles
    sim/R3BTofi.cxx
    sim/R3BGeoTofi.cxx
    pars/R3BGeoTofiPar.cxx
    pars/R3BTofiDigiPar.cxx
    digi/R3BTofiDigitizerCal.cxx
    pars/R3BTofiHisto2HitPar.cxx
    pars/R3BTofiHitPar.cxx
    pars/R3BTofiHitModulePar.cxx
    pars/R3BTofiChangePar.cxx
    pars/R3BTofiContFact.cxx
    calibration/R3BTofiCal2Histo.cxx
    calibration/R3BTofiMapped2Cal.cxx
    calibration/R3BTofiMapped2CalPar.cxx
    calibration/R3BTofiCal2HitS494.cxx)

set(HEADERS
    #Put here your sourcefiles
    sim/R3BTofi.h
    sim/R3BGeoTofi.h
    pars/R3BGeoTofiPar.h
    pars/R3BTofiDigiPar.h
    digi/R3BTofiDigitizerCal.h
    pars/R3BTofiHisto2HitPar.h
    pars/R3BTofiHitPar.h
    pars/R3BTofiHitModulePar.h
    pars/R3BTofiChangePar.h
    pars/R3BTofiContFact.h
    calibration/R3BTofiCal2Histo.h
    calibration/R3BTofiMapped2Cal.h
    calibration/R3BTofiMapped2CalPar.h
    calibration/R3BTofiCal2HitS494.h)

add_library_with_dictionary(
    LIBNAME
    R3BTofi
    LINKDEF
    TofiLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    calibration
    digi
    pars
    sim
    DEPENDENCIES
    R3BTCal
    R3BTracking)
