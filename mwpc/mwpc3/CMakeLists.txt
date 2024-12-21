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
    R3BMwpc3.cxx
    R3BMwpc3ContFact.cxx
    R3BMwpc3CalPar.cxx
    R3BMwpc3Mapped2Cal.cxx
    R3BMwpc3Mapped2CalPar.cxx
    R3BMwpc3Cal2Hit.cxx)

set(HEADERS
    R3BMwpc3.h
    R3BMwpc3ContFact.h
    R3BMwpc3CalPar.h
    R3BMwpc3Mapped2Cal.h
    R3BMwpc3Mapped2CalPar.h
    R3BMwpc3Cal2Hit.h)

add_library_with_dictionary(
    LIBNAME
    R3BMwpc3
    LINKDEF
    Mwpc3LinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BSofData
    R3BTracking)
