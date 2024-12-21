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

# this file is automatically symlinked to GNUmakefile when
# the configuration script ${CMAKE_SOURCE_DIR}/lustretar.sh exists.

$(info You are using the c3w makefile wrapper)

SRC_UTIL_MAKE=$(shell readlink -f GNUmakefile)
SRC_UTIL=$(shell dirname $(SRC_UTIL_MAKE) )

c3w: all
	@echo "Trying to create c3w tarball. "
	$(SRC_UTIL)/build_c3w.sh

include Makefile # include the cmake generated makefile
