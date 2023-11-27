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

set(SRCS R3BXBall.cxx R3BGeoXBall.cxx R3BGeoXBallPar.cxx R3BXBallContFact.cxx)
set(HEADERS R3BXBall.h R3BGeoXBall.h R3BGeoXBallPar.h R3BXBallContFact.h)

add_library_with_dictionary(
    LIBNAME
    R3BXBall
    LINKDEF
    ${LINKDEF}
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData)
