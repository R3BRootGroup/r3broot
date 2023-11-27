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

#########################################################
# temporary workaround for testR3BPhaseSpaceGenerator.cxx
# while FairVersion.h is not working
execute_process(COMMAND ${FAIRROOTPATH}/bin/fairroot-config --major_version
                OUTPUT_VARIABLE fairroot_major_version)
execute_process(COMMAND ${FAIRROOTPATH}/bin/fairroot-config --minor_version
                OUTPUT_VARIABLE fairroot_minor_version)

if((${fairroot_major_version} EQUAL 17 AND ${fairroot_minor_version} EQUAL 10)
   OR ${fairroot_major_version} GREATER 17)
    message("using new FairPrimaryGenerator interface")
    add_definitions(-DFairPrimaryGeneratorAddTrackNewInterface)
endif((${fairroot_major_version} EQUAL 17 AND ${fairroot_minor_version} EQUAL 10)
      OR ${fairroot_major_version} GREATER 17)
#########################################################

set(SRCS
    R3BSpecificGenerator.cxx
    R3BReadKinematics.cxx
    R3BCDGenerator.cxx
    R3BBeamInfo.cxx
    R3BBackTracking.cxx
    R3BBackTrackingStorageState.cxx
    R3BAsciiGenerator.cxx
    R3BCosmicGenerator.cxx
    R3BCryAsciiGenerator.cxx
    R3Bp2pGenerator.cxx
    R3BLandGenerator.cxx
    R3BCALIFATestGenerator.cxx
    R3BParticleGenerator.cxx
    R3BIonGenerator.cxx
    R3BGammaGenerator.cxx
    R3BPhaseSpaceGenerator.cxx
    R3BDistribution1D.cxx
    R3BDistribution2D.cxx
    R3BDistribution3D.cxx
    R3Bp2pevtGenerator.cxx
    R3BParticleSelector.cxx
    R3BBeamProperties.cxx
    R3BINCLRootGenerator.cxx)

set(HEADERS
    R3BSpecificGenerator.h
    R3BReadKinematics.h
    R3BCDGenerator.h
    R3BBeamInfo.h
    R3BBackTracking.h
    R3BBackTrackingStorageState.h
    R3BAsciiGenerator.h
    R3BCosmicGenerator.h
    R3BCryAsciiGenerator.h
    R3Bp2pGenerator.h
    R3BLandGenerator.h
    R3BCALIFATestGenerator.h
    R3BParticleGenerator.h
    R3BIonGenerator.h
    R3BGammaGenerator.h
    R3BPhaseSpaceGenerator.h
    R3BDistribution1D.h
    R3BDistribution2D.h
    R3BDistribution3D.h
    R3Bp2pevtGenerator.h
    R3BParticleSelector.h
    R3BBeamProperties.h
    R3BINCLRootGenerator.h)

add_library_with_dictionary(
    LIBNAME
    R3BGen
    LINKDEF
    R3BGenLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData
    PRIVATE_DEPENDENCIES
    Geant4::G4materials
)

add_subdirectory(test)
