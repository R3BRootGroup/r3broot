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
    R3BPsp.cxx
    R3BGeoPsp.cxx
    # R3BGeoPspPar.cxx
    # R3BPspContFact.cxx
    R3BPspDigiPar.cxx
    R3BPspDigitizer.cxx
    # R3BPspxMapped2Cal.cxx
    R3BPspxMapped2Precal.cxx
    R3BPspxPrecal2Cal.cxx
    R3BPspxCal2Hit.cxx
    R3BPspxReader.cxx
    R3BPspxOnlineSpectra.cxx
    par/R3BPspxMappedPar.cxx
    par/R3BPspxPrecalPar.cxx
    par/R3BPspxCalPar.cxx
    par/R3BPspxHitPar.cxx
    par/R3BPspxContFact.cxx)

set(HEADERS
    R3BPsp.h
    R3BGeoPsp.h
    # R3BGeoPspPar.h
    # R3BPspContFact.h
    R3BPspDigiPar.h
    R3BPspDigitizer.h
    # R3BPspxMapped2Cal.h
    R3BPspxMapped2Precal.h
    R3BPspxPrecal2Cal.h
    R3BPspxCal2Hit.h
    R3BPspxReader.h
    R3BPspxOnlineSpectra.h
    par/R3BPspxMappedPar.h
    par/R3BPspxPrecalPar.h
    par/R3BPspxCalPar.h
    par/R3BPspxHitPar.h
    par/R3BPspxContFact.h)

add_library_with_dictionary(
    LIBNAME
    R3BPsp
    LINKDEF
    PspLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    par
    DEPENDENCIES
    R3BTracking R3BSource)
