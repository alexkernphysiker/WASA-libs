cmake_minimum_required(VERSION 3.0)
project(WASAatCOSY_etamesic_3He_search_libs)
add_subdirectory(FitGen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/FitGen/include/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS "-pthread -std=gnu++11")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
foreach(mode QUIET REQUIRED)
  find_package(ROOT 5 ${mode} COMPONENTS Core Gui Gpad Graf Hist RIO Physics Tree)
  if(ROOT_USE_FILE)
    include(${ROOT_USE_FILE})
  endif()
  if(ROOT_FOUND)
    break()
  endif()
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/fallback)
endforeach()
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_DEFINITIONS})
message(STATUS "ROOT version: ${ROOT_VERSION}")
if(ROOT_FOUND AND NOT ROOT_VERSION VERSION_LESS "6.0")
  message(FATAL_ERROR "ROOT 6.0 is not compatible")
endif()

get_directory_property(has_parent PARENT_DIRECTORY)
if(has_parent)
  set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} PARENT_SCOPE)
  set(ROOT_DEFINITIONS ${ROOT_DEFINITIONS} PARENT_SCOPE)
  if(NOT MSVC)
    set(ROOT_DEFINITIONS -std=c++11 ${ROOT_DEFINITIONS} PARENT_SCOPE)
  endif()
  set(ROOT_rootcint_CMD   ${ROOT_rootcint_CMD}   PARENT_SCOPE)
  set(rootcint_cmd        ${ROOT_rootcint_CMD}   PARENT_SCOPE)
  set(ROOTCINT_EXECUTABLE ${ROOTCINT_EXECUTABLE} PARENT_SCOPE)
  set(ROOT_LIBRARY_DIR    ${ROOT_LIBRARY_DIR}    PARENT_SCOPE)
endif()

link_directories(${ROOT_LIBRARY_DIR})

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
