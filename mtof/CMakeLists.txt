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
    R3BmTof.cxx
    R3BGeomTof.cxx
    R3BGeomTofPar.cxx
    R3BmTofContFact.cxx
    R3BmTofDigiPar.cxx
    R3BmTofDigitizer.cxx)

set(HEADERS
    R3BmTof.h
    R3BGeomTof.h
    R3BGeomTofPar.h
    R3BmTofContFact.h
    R3BmTofDigiPar.h
    R3BmTofDigitizer.h)

add_library_with_dictionary(
    LIBNAME
    R3BmTof
    LINKDEF
    mTofLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BTracking)
