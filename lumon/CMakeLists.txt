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

set(SRCS ELILuMon.cxx ELIGeoLuMon.cxx ELIGeoLuMonPar.cxx ELILuMonContFact.cxx)

set(HEADERS ELILuMon.h ELIGeoLuMon.h ELIGeoLuMonPar.h ELILuMonContFact.h)

add_library_with_dictionary(
    LIBNAME
    ELILuMon
    LINKDEF
    LuMonLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData)
