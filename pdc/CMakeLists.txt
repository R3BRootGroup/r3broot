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
    R3BGeoPdcPar.cxx
    R3BPdcContFact.cxx
    R3BPdcMapped2CalPar.cxx
    R3BPdcMapped2Cal.cxx
    R3BPdcCal2Hit.cxx)

set(HEADERS
    R3BGeoPdcPar.h
    R3BPdcContFact.h
    R3BPdcMapped2CalPar.h
    R3BPdcMapped2Cal.h
    R3BPdcCal2Hit.h)

add_library_with_dictionary(
    LIBNAME
    R3BPdc
    LINKDEF
    PdcLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BTracking
    R3BTCal
    R3BSource)
