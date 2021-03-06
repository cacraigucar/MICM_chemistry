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

#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# The Fortran compiler/linker flag inserted by cmake to create shared libraries
# with the Intel compiler is deprecated (-i_dynamic), correct here.
# CMAKE_Fortran_COMPILER_ID = {"Intel", "PGI", "GNU", "Clang", "MSVC", ...}
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
endif()

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# By default we want a shared library
option(BUILD_SHARED_LIBS "Build a shared library" ON)

#------------------------------------------------------------------------------
# Add the CCPP include/module directory
set(CCPP_INCLUDE_DIRS "" CACHE FILEPATH "Path to ccpp includes")
set_property(DIRECTORY PROPERTY INCLUDE_DIRECTORIES ${CCPP_INCLUDE_DIRS})

#------------------------------------------------------------------------------
# Add the CCPP library
set(CCPP_LIB_DIRS "" CACHE FILEPATH "Path to ccpp library")
link_directories(${CCPP_LIB_DIRS})
list(APPEND LIBS "ccpp")

# Add the required preprocessor flags so that cmake can sort out the dependencies
ADD_DEFINITIONS(-DNEMS_GSM)

#------------------------------------------------------------------------------
# Set the sources: physics schemes
message (STATUS "MODEL is " ${PROJECT})
INCLUDE_DIRECTORIES(generated/${PROJECT})

set(SCHEMES
      ./src/machine.F
      ./src/const_props_mod.F90
      ./src/json_loader.F90
      ./src/ODE_solver.F90
      ./src/mozart_solver.F90
      ./src/rosenbrock_solver.F90
      ./src/kinetics/kinetics.F90
      ./src/kinetics/k_rateConst.F90
      ./src/kinetics/prepare_chemistry.F90
      ./src/kinetics/kinetics_module.F90
      ./src/chemistry_driver_ros.F90
      ./src/chemistry_driver_moz.F90
      ./src/tuv/la_srb_mod.F90
      ./src/tuv/module_params.F90
      ./src/tuv/module_prates_tuv.F90
      ./src/tuv/module_rxn.F90
      ./src/tuv/module_xsections.F90
      ./src/tuv/numer.F90
      ./src/tuv/params_mod.F90
      ./src/tuv/rad_trans.F90
      ./src/tuv/tuv_error_mod.F90
      ./src/tuv/tuv_photolysis.F90
      ./src/tuv/tuv_subs.F90
      ./src/photolysis_interstitial.F90
      ./src/mass_quantities_util.F90
      ./src/json_lib/json_file_module.F90
      ./src/json_lib/json_kinds.F90
      ./src/json_lib/json_module.F90
      ./src/json_lib/json_parameters.F90
      ./src/json_lib/json_string_utilities.F90
      ./src/json_lib/json_value_module.F90
  )
  # Set the sources: physics scheme caps
set(CAPS
      ./caps/chemistry_driver_ros_cap.F90
      ./caps/chemistry_driver_moz_cap.F90
      ./caps/k_rateConst_cap.F90
      ./caps/photolysis_interstitial_cap.F90
      ./caps/tuv_photolysis_cap.F90
      ./caps/mass_quantities_util_cap.F90
  )
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  # Basic settings
  set (CMAKE_C_FLAGS       "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS     "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS "-O0 -g -fPIC" CACHE STRING "" FORCE)
  # Compiler-dependent settings
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fcheck=all -fbacktrace -ffpe-trap=zero -finit-real=nan -ggdb -ffree-line-length-none -cpp -fdefault-real-8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fpe0 -warn -traceback -debug extended -fpp -r8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mnoipa -traceback -Mpreprocess -r8")
  else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
  endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
  set (CMAKE_C_FLAGS       "-O2 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS     "-O2 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS "-O2 -g -fPIC" CACHE STRING "" FORCE)
  # Compiler-dependent settings
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -cpp -fdefault-real-8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ftz -debug minimal -fpp -r8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mipa=noconst -Mpreprocess -r8")
  else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
  endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")

#apply general fortran tags to all fortran source files
if(${CMAKE_VERSION} LESS 3.3)
  string (REPLACE ";" " " f_flags_str "${f_flags}")
  SET_PROPERTY(SOURCE ${SCHEMES} ${CAPS} APPEND_STRING PROPERTY COMPILE_FLAGS " ${f_flags_str}")
else(${CMAKE_VERSION} LESS 3.3)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
endif (${CMAKE_VERSION} LESS 3.3)


#this should call the custom FindNetCDF.cmake file in /src/cmake/modules and populate the NETCDF_INCLUDES and NETCDF_LIBRARIES variables
SET(NETCDF_F90 "YES")
FIND_PACKAGE(NetCDF)

INCLUDE_DIRECTORIES(${NETCDF_INCLUDES})

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp-framework/src)

#------------------------------------------------------------------------------

add_library(micmchem ${SCHEMES} ${CAPS})
target_link_libraries(micmchem LINK_PUBLIC ${LIBS} )
set_target_properties(micmchem PROPERTIES VERSION ${PROJECT_VERSION}
                                     SOVERSION ${PROJECT_VERSION_MAJOR}
                                     COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                                     LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
# DH* hack for PGI compiler: rename objects in scheme cap object files for ISO_C compliancy
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
set(CAPOBJS)
foreach(cap ${CAPS})
    string(REPLACE "_cap.F90" "_cap.F90.o" capobj "./${CMAKE_FILES_DIRECTORY}/micmchem.dir/${cap}")
    list(APPEND CAPOBJS ${capobj})
endforeach(cap)

add_custom_command(TARGET micmchem
                   PRE_LINK
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/pgifix.py --cmake ${CAPOBJS}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   COMMENT "Running pgifix_wrapper.py over all scheme caps")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
# *DH end hack for PGI compiler

