##############################################################################
#   Copyright (C) 2022 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2022-2023 Members of R3B Collaboration                     #
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
    sim/R3BRpc.cxx
    pars/R3BRpcContFact.cxx
    #R3BRpcDigitizer.cxx
    calibration/R3BRpcMapped2PreCalPar.cxx
    calibration/R3BRpcMapped2PreCal.cxx
    calibration/R3BRpcPreCal2CalPar.cxx
    calibration/R3BRpcPreCal2Cal.cxx
    pars/R3BRpcTotCalPar.cxx
    calibration/R3BRpcCal2Hit.cxx
    calibration/R3BRpcCal2HitPar.cxx
    pars/R3BRpcHitPar.cxx
    pars/R3BRpcPars4Sim.cxx
    online/R3BRpcOnlineSpectra.cxx)

set(HEADERS
    sim/R3BRpc.h
    pars/R3BRpcContFact.h
    #R3BRpcDigitizer.h
    calibration/R3BRpcMapped2PreCalPar.h
    calibration/R3BRpcMapped2PreCal.h
    calibration/R3BRpcPreCal2CalPar.h
    calibration/R3BRpcPreCal2Cal.h
    pars/R3BRpcTotCalPar.h
    calibration/R3BRpcCal2Hit.h
    calibration/R3BRpcCal2HitPar.h
    pars/R3BRpcHitPar.h
    pars/R3BRpcPars4Sim.h
    online/R3BRpcOnlineSpectra.h)

add_library_with_dictionary(
    LIBNAME
    R3BRpc
    LINKDEF
    RpcLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    calibration
    digi
    online
    pars
    sim
    DEPENDENCIES
    R3BTracking
    R3BTCal)

if(BUILD_GEOMETRY)
    add_subdirectory(geobase)
endif()
add_subdirectory(test)
