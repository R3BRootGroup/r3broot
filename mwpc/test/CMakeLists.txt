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


generate_root_test_script(${R3BROOT_SOURCE_DIR}/mwpc/test/testMwpcSimulation.C)
add_test(NAME MwpcSimulation COMMAND ${R3BROOT_BINARY_DIR}/mwpc/test/testMwpcSimulation.sh)
set_tests_properties(MwpcSimulation PROPERTIES TIMEOUT "2000")
set_tests_properties(MwpcSimulation PROPERTIES PASS_REGULAR_EXPRESSION
                                                  "Macro finished successfully.")
