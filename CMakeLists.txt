cmake_minimum_required(VERSION 3.0)
project(WASAatCOSY_etamesic_3He_search_libs)
set(CMAKE_CXX_FLAGS "-pthread -std=gnu++11")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/star-cmake")
find_package(ROOT)
if(NOT ROOT_FOUND)
	message(FATAL_ERROR "Fatal error: ROOT package not found")
endif()

add_subdirectory(FitGen)
include_directories(${FITGEN_INC})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIR})

get_directory_property(has_parent PARENT_DIRECTORY)
if(has_parent)
  set(ROOT_LIBRARY_DIR ${ROOT_LIBRARY_DIR} PARENT_SCOPE)
  set(WASA_LIBS_INC ${CMAKE_CURRENT_SOURCE_DIR} ${FITGEN_INC} ${ROOT_INCLUDE_DIR} PARENT_SCOPE)
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
