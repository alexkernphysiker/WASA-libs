        if(NOT $ENV{PLUTOSYS} STREQUAL "")
        link_directories($ENV{PLUTOSYS})
                include_directories(
                          $ENV{PLUTOSYS}/src
                          $ENV{PLUTOSYS}/plugins/brems
                          $ENV{PLUTOSYS}/plugins/dalitz_mod
                          $ENV{PLUTOSYS}/plugins/elementary
                          $ENV{PLUTOSYS}/plugins/eta_decays
                          $ENV{PLUTOSYS}/plugins/fairroot
                          $ENV{PLUTOSYS}/plugins/hades
                          $ENV{PLUTOSYS}/plugins/nucleus_fermi
                          $ENV{PLUTOSYS}/plugins/pdg_unigen
                          $ENV{PLUTOSYS}/plugins/pion_beam
                          $ENV{PLUTOSYS}/plugins/scatter_mod
                          $ENV{PLUTOSYS}/plugins/strangeness
                          $ENV{PLUTOSYS}/plugins/w_decays
                )
                add_executable(run_plutto pluto.cpp)
                target_link_libraries(run_plutto Pluto ${ROOT_LIBRARIES} ${libs} experiment kinematics)
        endif()
