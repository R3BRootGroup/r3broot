add_executable(neulandSim neulandSim.cxx)
target_link_libraries(neulandSim PRIVATE R3BNeulandSimulation ${Geant4_LIBRARIES}
                                        ${Geant4VMC_LIBRARIES} ExTestDetector)

add_executable(neulandAna neulandAna.cxx)
target_link_libraries(neulandAna PRIVATE R3BNeulandDigitizing)
