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


if(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_full.geo.root")
message(STATUS "Starting to generate the standard root geometry for CALIFA (FULL)")
execute_process(COMMAND root -l -q -x ${R3BROOT_SOURCE_DIR}/califa/geobase/create_califa_geo.C)
message(STATUS "${BYellow}Full CALIFA standard root geometry has been created${CR} ")
else(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_full.geo.root")
message(STATUS "${BYellow}Full CALIFA standard root geometry already exists${CR} ")
endif(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_full.geo.root")
