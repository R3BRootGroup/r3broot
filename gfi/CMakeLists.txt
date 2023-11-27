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

set(SRCS R3BGfi.cxx R3BGeoGfi.cxx R3BGeoGfiPar.cxx R3BGfiContFact.cxx R3BGfiDigiPar.cxx
         R3BGfiDigitizer.cxx)

set(HEADERS R3BGfi.h R3BGeoGfi.h R3BGeoGfiPar.h R3BGfiContFact.h R3BGfiDigiPar.h R3BGfiDigitizer.h)

add_library_with_dictionary(
    LIBNAME
    R3BGfi
    LINKDEF
    GfiLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData)
