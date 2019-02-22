# Set default project to unknown
if(NOT PROJECT)
  message(STATUS "Setting CCPP project to 'unknown' as none was specified.")
  set(PROJECT "Unknown")
endif (NOT PROJECT)

#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.11)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
    project(micmchem VERSION 1.0.0)
else(POLICY CMP0048)
    project(micmchem)
    set(PROJECT_VERSION 1.0.0)
    set(PROJECT_VERSION_MAJOR 1)
    set(PROJECT_VERSION_MINOR 0)
    set(PROJECT_VERSION_PATCH 0)
endif(POLICY CMP0048)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

#------------------------------------------------------------------------------
set(PACKAGE "micm-chemistry")
set(AUTHORS  "Cheryl Craig" "Francis Vitt" "Andrew Conley" "Stacy Walters")

#------------------------------------------------------------------------------
# Enable Fortran
enable_language(Fortran)

# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../ccpp-framework/cmake")

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
else(OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# By default we want a static library
option(BUILD_SHARED_LIBS "Build a shared library" OFF)


#------------------------------------------------------------------------------
# Set the sources: physics schemes
message (STATUS "MODEL is " ${PROJECT})
INCLUDE_DIRECTORIES(generated/${PROJECT})

set(SCHEMES
      ./src/json_loader.F90  
      ./src/const_props_mod.F90
      ./src/kinetics/kinetics.F90
      ./src/kinetics/k_rateConst.F90
      ./src/kinetics/kinetics_module.F90
      ./src/json_lib/json_file_module.F90
      ./src/json_lib/json_kinds.F90
      ./src/json_lib/json_module.F90
      ./src/json_lib/json_parameters.F90
      ./src/json_lib/json_string_utilities.F90
      ./src/json_lib/json_value_module.F90
  )

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${GFS_PHYS_SRC}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

#apply general fortran tags to all fortran source files
if(${CMAKE_VERSION} LESS 3.3)
  string (REPLACE ";" " " f_flags_str "${f_flags}")
  SET_PROPERTY(SOURCE ${SCHEMES} ${CCPP_CAPS}
    APPEND_STRING PROPERTY COMPILE_FLAGS " ${f_flags_str}")
else(${CMAKE_VERSION} LESS 3.3)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
endif (${CMAKE_VERSION} LESS 3.3)

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/MICM_chemistry)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------

add_library(micmchem ${SCHEMES})
target_link_libraries(micmchem LINK_PUBLIC ${LIBS} )
set_target_properties(base PROPERTIES VERSION ${PROJECT_VERSION}
                           SOVERSION ${PROJECT_VERSION_MAJOR}
                           COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                           LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
set_target_properties(micmchem PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR}
                                     COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                                     LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
