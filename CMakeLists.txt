cmake_minimum_required(VERSION 2.6)
project(WASAatCOSY_etamesic_3He_search_libs)
add_subdirectory(FitGen)
include_directories(FitGen/include/)
include_directories(.)
set(CMAKE_CXX_FLAGS "-pthread -std=gnu++11")
if(NOT $ENV{ROOTSYS} STREQUAL "")
	include_directories($ENV{ROOTSYS}/include/)
	link_directories($ENV{ROOTSYS}/lib/)
	set(ROOT_LIBRARIES Core Gui Gpad Graf Hist RIO Physics Tree)
	
	file(GLOB EXPERIMENT Experiment/*.cpp)
	add_library(experiment SHARED ${EXPERIMENT})
	target_link_libraries(experiment ${ROOT_LIBRARIES})
	
	file(GLOB KINEMATICS Kinematics/*.cpp)
	add_library(kinematics SHARED ${KINEMATICS})

	file(GLOB RECONSTRUCTION ReconstructionFit/*.cpp)
	add_library(reconstructionfit SHARED ${RECONSTRUCTION})
	target_link_libraries(reconstructionfit FitGen)

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
endif()
