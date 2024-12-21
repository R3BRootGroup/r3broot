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

# Defines the following variables:
#
#   ClangFormat_FOUND - Found clang-format and git-clang-format
#   CLANG_FORMAT_BIN - clang-format executable 
#   GIT_CLANG_FORMAT_BIN - git-clang-format executable 

find_program(CLANG_FORMAT_BIN
  NAMES clang-format-11
        clang-format-10
        clang-format-9
        clang-format-8
)

find_program(GIT_CLANG_FORMAT_BIN
  NAMES git-clang-format-11
        git-clang-format-10
        git-clang-format-9
        git-clang-format-8
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(R3BClangFormat
  REQUIRED_VARS CLANG_FORMAT_BIN GIT_CLANG_FORMAT_BIN
)
