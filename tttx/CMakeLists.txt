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
    #Put here your sourcefiles
    #pars/R3BTttxMappingPar.cxx
    pars/R3BTttxStripCalPar.cxx
    pars/R3BTttxHitPar.cxx
    pars/R3BTttxContFact.cxx
    #calibration/R3BTttxMapped2CalPar.cxx
    calibration/R3BTttxMapped2Cal.cxx
    calibration/R3BTttxCal2Hit.cxx
    online/R3BTttxOnlineSpectra.cxx)

set(HEADERS
    #Put here your sourcefiles
    #pars/R3BTttxMappingPar.h
    pars/R3BTttxStripCalPar.h
    pars/R3BTttxHitPar.h
    pars/R3BTttxContFact.h
    #calibration/R3BTttxMapped2CalPar.h
    calibration/R3BTttxMapped2Cal.h
    calibration/R3BTttxCal2Hit.h
    online/R3BTttxOnlineSpectra.h)

add_library_with_dictionary(
    LIBNAME
    R3BTttx
    LINKDEF
    TttxLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    calibration
    online
    pars
    DEPENDENCIES
    R3BTracking)
