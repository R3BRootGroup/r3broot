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

if(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/glad_v2023.1.geo.root")
message(STATUS "Starting to generate the standard root geometry for GLAD")
execute_process(COMMAND root -l -q -x ${R3BROOT_SOURCE_DIR}/glad/geobase/create_glad_geo.C)
message(STATUS "${BYellow}GLAD standard root geometry has been created${CR} ")
else(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/glad_v2023.1.geo.root")
message(STATUS "${BYellow}GLAD standard root geometry already created${CR} ")
endif(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/glad_v2023.1.geo.root")

if(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/glad_s455_v2023.1.geo.root")
message(STATUS "Starting to generate the standard root geometry of GLAD for S455")
execute_process(COMMAND root -l -q -x "${R3BROOT_SOURCE_DIR}/glad/geobase/create_glad_geo.C(true)")
message(STATUS "${BYellow}GLAD standard root geometry for S455 has been created${CR} ")
else(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/glad_s455_v2023.1.geo.root")
message(STATUS "${BYellow}GLAD standard root geometry for S455 already created${CR} ")
endif(NOT EXISTS "${R3BROOT_SOURCE_DIR}/geometry/glad_s455_v2023.1.geo.root")

generate_root_test_script(${R3BROOT_SOURCE_DIR}/glad/geobase/create_glad_geo.C)
add_test(NAME GladRootGeometry COMMAND ${R3BROOT_BINARY_DIR}/glad/geobase/create_glad_geo.sh)
set_tests_properties(GladRootGeometry PROPERTIES TIMEOUT "2000")
set_tests_properties(GladRootGeometry PROPERTIES PASS_REGULAR_EXPRESSION
                                                  "Macro finished successfully.")
