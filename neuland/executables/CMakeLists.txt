add_executable(neulandSim neulandSim.cxx)
target_link_libraries(
    neulandSim
    PRIVATE R3BNeulandSimulation
            R3BGen
            ${Geant4_LIBRARIES}
            ${Geant4VMC_LIBRARIES})

add_executable(neulandAna neulandAna.cxx)
target_link_libraries(neulandAna PRIVATE R3BNeulandDigitizing)

if(CONAN_ENABLED)
    add_executable(neuland_lmd_to_map neuland_lmd_to_map.cxx)
    target_link_libraries(neuland_lmd_to_map PRIVATE R3BNeulandCalibration R3BSource -pthread)

    add_executable(neuland_map_to_cal neuland_map_to_cal.cxx)
    target_link_libraries(neuland_map_to_cal PRIVATE R3BNeulandCalibration R3BSource)

    add_executable(neuland_cal_to_hit neuland_cal_to_hit.cxx)
    target_link_libraries(neuland_cal_to_hit PRIVATE R3BNeulandCalibration R3BSource)

    add_executable(neuland_online_monitor neuland_online_monitor.cxx)
    target_link_libraries(
        neuland_online_monitor
        PRIVATE R3BNeulandOnline
                R3BLos
                R3BSource
                -pthread)

    add_executable(neuland_par_merger neuland_par_merger.cxx)
    target_link_libraries(neuland_par_merger PRIVATE R3BNeulandShared)

    add_executable(neuland_create_base_par neuland_create_base_par.cxx)
    target_link_libraries(neuland_create_base_par PRIVATE R3BNeulandCalibration R3BNeulandShared)
endif()
