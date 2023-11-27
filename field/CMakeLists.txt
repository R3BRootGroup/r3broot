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
    R3BFieldConst.cxx
    R3BFieldContFact.cxx
    R3BFieldPar.cxx
    R3BFieldCreator.cxx
    R3BGladFieldMap.cxx
    R3BFieldInterp.cxx
    R3BAladinFieldMap.cxx)

set(HEADERS
    R3BFieldConst.h
    R3BFieldContFact.h
    R3BFieldPar.h
    R3BFieldCreator.h
    R3BGladFieldMap.h
    R3BFieldInterp.h
    R3BAladinFieldMap.h)

add_library_with_dictionary(
    LIBNAME
    Field
    LINKDEF
    FieldLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BBase)
