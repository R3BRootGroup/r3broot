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

# Create a library called "libR3BCal" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    R3BTof.cxx
    R3BGeoTof.cxx
    R3BGeoTofPar.cxx
    R3BTofContFact.cxx
    R3BTofDigiPar.cxx
    R3BTofDigitizer.cxx
    R3BTof2pDigiPar.cxx
    R3BTof2pDigitizer.cxx
    R3BTofdMapped2Cal.cxx
    R3BTofdMapped2CalPar.cxx
    R3BTofdCal2HitPar.cxx
    R3BTofdCal2Histo.cxx
    R3BTofdCal2HistoPar.cxx
    R3BTofdChangePar.cxx
    R3BTofdCal2Hit.cxx
    R3BTofdCal2HitS454.cxx
    R3BTofdCal2HitS494.cxx
    R3BTofdHitModulePar.cxx
    R3BTofdHitPar.cxx
    R3BTofd.cxx
    R3BTofdDigitizer.cxx
    R3BTofdDigitizerCal.cxx
    R3BGeoTofd.cxx
    R3BGeoTofdPar.cxx
    R3BTofdContFact.cxx
    R3BTofdDigiPar.cxx
    R3BPtofMapped2CalPar.cxx
    R3BPtofMapped2Cal.cxx
    calibration/R3BPtofCal2HitPar.cxx
    calibration/R3BPtofCal2Hit.cxx
    calibration/R3BPtofHitPar.cxx
    calibration/R3BPtofParFact.cxx)

set(HEADERS
    R3BTof.h
    R3BGeoTof.h
    R3BGeoTofPar.h
    R3BTofContFact.h
    R3BTofDigiPar.h
    R3BTofDigitizer.h
    R3BTof2pDigiPar.h
    R3BTof2pDigitizer.h
    R3BTofdMapped2Cal.h
    R3BTofdMapped2CalPar.h
    R3BTofdCal2HitPar.h
    R3BTofdCal2Histo.h
    R3BTofdCal2HistoPar.h
    R3BTofdChangePar.h
    R3BTofdCal2Hit.h
    R3BTofdCal2HitS454.h
    R3BTofdCal2HitS494.h
    R3BTofdHitModulePar.h
    R3BTofdHitPar.h
    R3BTofd.h
    R3BTofdDigitizer.h
    R3BTofdDigitizerCal.h
    R3BGeoTofd.h
    R3BGeoTofdPar.h
    R3BTofdContFact.h
    R3BTofdDigiPar.h
    R3BPtofMapped2CalPar.h
    R3BPtofMapped2Cal.h
    calibration/R3BPtofCal2HitPar.h
    calibration/R3BPtofCal2Hit.h
    calibration/R3BPtofHitPar.h
    calibration/R3BPtofParFact.h)

add_library_with_dictionary(
    LIBNAME
    R3BTof
    LINKDEF
    TofLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    calibration
    DEPENDENCIES
    R3BTracking
    R3BTCal)
