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

generate_root_test_script(${R3BROOT_SOURCE_DIR}/califa/geobase/create_califa_geo.C)
add_test(NAME CalifaRootGeometry COMMAND ${R3BROOT_BINARY_DIR}/califa/geobase/create_califa_geo.sh)
set_tests_properties(CalifaRootGeometry PROPERTIES TIMEOUT "2000")
set_tests_properties(CalifaRootGeometry PROPERTIES PASS_REGULAR_EXPRESSION
                                                  "Macro finished successfully.")
                                                  
if(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2019.11.geo.root")
message(STATUS "Starting to generate the standard root geometry v2019.11 for CALIFA")
execute_process(COMMAND root -l -q -x "${R3BROOT_SOURCE_DIR}/califa/geobase/create_califa_geo.C(1)")
message(STATUS "${BYellow}CALIFA standard root geometry v2019.11 has been created${CR} ")
else(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2019.11.geo.root")
message(STATUS "${BYellow}CALIFA standard root geometry v2019.11 already exists${CR} ")
endif(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2019.11.geo.root")

if(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2021.3.geo.root")
message(STATUS "Starting to generate the standard root geometry v2021.3 for CALIFA")
execute_process(COMMAND root -l -q -x "${R3BROOT_SOURCE_DIR}/califa/geobase/create_califa_geo.C(2)")
message(STATUS "${BYellow}CALIFA standard root geometry v2021.3 has been created${CR} ")
else(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2021.3.geo.root")
message(STATUS "${BYellow}CALIFA standard root geometry v2021.3 already exists${CR} ")
endif(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2021.3.geo.root")

if(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2024.1.geo.root")
message(STATUS "Starting to generate the standard root geometry v2024.1 for CALIFA")
execute_process(COMMAND root -l -q -x "${R3BROOT_SOURCE_DIR}/califa/geobase/create_califa_geo.C(3)")
message(STATUS "${BYellow}CALIFA standard root geometry v2024.1 has been created${CR} ")
else(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2024.1.geo.root")
message(STATUS "${BYellow}CALIFA standard root geometry v2024.1 already exists${CR} ")
endif(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_v2024.1.geo.root")

if(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_full.geo.root")
message(STATUS "Starting to generate the standard root geometry for CALIFA (FULL)")
execute_process(COMMAND root -l -q -x ${R3BROOT_SOURCE_DIR}/califa/geobase/create_califa_geo.C)
message(STATUS "${BYellow}Full CALIFA standard root geometry has been created${CR} ")
else(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_full.geo.root")
message(STATUS "${BYellow}Full CALIFA standard root geometry already exists${CR} ")
endif(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/califa_full.geo.root")
