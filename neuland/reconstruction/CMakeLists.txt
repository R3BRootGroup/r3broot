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

# link_libraries(/home/jmayer/.pyenv/versions/3.7.6/lib/libpython3.7m.so.1.0)
set(SRCS
    multiplicity/R3BNeulandMultiplicityBayes.cxx
    multiplicity/R3BNeulandMultiplicityBayesPar.cxx
    multiplicity/R3BNeulandMultiplicityBayesTrain.cxx
    multiplicity/R3BNeulandMultiplicityCalorimetric.cxx
    multiplicity/R3BNeulandMultiplicityCalorimetricPar.cxx
    multiplicity/R3BNeulandMultiplicityCalorimetricTrain.cxx
    multiplicity/R3BNeulandMultiplicityCheat.cxx
    multiplicity/R3BNeulandMultiplicityFixed.cxx
    # multiplicity/R3BNeulandMultiplicityScikit.cxx
    neutrons/R3BNeulandNeutronsCheat.cxx
    neutrons/R3BNeulandNeutronsRValue.cxx
    # neutrons/R3BNeulandNeutronsScikit.cxx
    # neutrons/R3BNeulandNeutronsKeras.cxx
    R3BNeulandReconstructionContFact.cxx
    R3BNeulandNeutronReconstructionMon.cxx
    R3BNeulandNeutronReconstructionStatistics.cxx)
set(HEADERS
    multiplicity/R3BNeulandMultiplicityBayes.h
    multiplicity/R3BNeulandMultiplicityBayesPar.h
    multiplicity/R3BNeulandMultiplicityBayesTrain.h
    multiplicity/R3BNeulandMultiplicityCalorimetric.h
    multiplicity/R3BNeulandMultiplicityCalorimetricPar.h
    multiplicity/R3BNeulandMultiplicityCalorimetricTrain.h
    multiplicity/R3BNeulandMultiplicityCheat.h
    multiplicity/R3BNeulandMultiplicityFixed.h
    # multiplicity/R3BNeulandMultiplicityScikit.h
    neutrons/R3BNeulandNeutronsCheat.h
    neutrons/R3BNeulandNeutronsRValue.h
    # neutrons/R3BNeulandNeutronsScikit.h
    # neutrons/R3BNeulandNeutronsKeras.h
    R3BNeulandReconstructionContFact.h
    R3BNeulandNeutronReconstructionMon.h
    R3BNeulandNeutronReconstructionStatistics.h)

add_library(R3BNeulandReconstruction SHARED ${SRCS})
target_include_directories(R3BNeulandReconstruction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} multiplicity neutrons)
target_link_libraries(R3BNeulandReconstruction PUBLIC R3BNeulandShared)
root_generate_dictionary(R3BNeulandReconstruction_dict ${HEADERS} MODULE R3BNeulandReconstruction LINKDEF
                         NeulandReconstructionLinkDef.h)
