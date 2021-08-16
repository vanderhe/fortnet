cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(FortnetUtils)

fnet_ensure_out_of_source_build()
fnet_load_build_settings()

set(FORTNET_VERSION "0.2")

project(fortnet VERSION ${FORTNET_VERSION} LANGUAGES Fortran)

fnet_setup_build_type()
fnet_load_toolchain_settings()
fnet_setup_global_compiler_flags()
fnet_get_release_name(RELEASE)

find_package(HDF5 REQUIRED COMPONENTS Fortran_HL)

find_package(MPI QUIET)
if(WITH_MPI AND NOT MPI_FORTRAN_FOUND)
  message(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER} is not MPI capable but is specified for "
    "a WITH_MPI=TRUE build")
elseif(NOT WITH_MPI AND MPI_FORTRAN_FOUND)
  message(WARNING "MPI enabled compiler ${CMAKE_Fortran_COMPILER} found for a non-MPI build. Your "
    "build will NOT be MPI-parallelised. Set WITH_MPI=TRUE in order to obtain an MPI-parallelised "
    "build.")
endif()

#
# Git repository check
#
find_package(Git)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  set(GIT_WORKING_COPY True)
endif()

#
# Check external libraries
#
set(PKG_CONFIG_REQUIRES)
set(PKG_CONFIG_LIBS)
set(PKG_CONFIG_LIBS_PRIVATE)


#
# Preprocessor details
#
set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")
fnet_add_fypp_defines(FYPP_FLAGS)

set(FYPP_CONFIG_FLAGS "${FYPP_FLAGS}")
# Make sure, the line-marker option is not set
list(REMOVE_ITEM FYPP_CONFIG_FLAGS "-n")
set(FYPP_BUILD_FLAGS "${FYPP_FLAGS}" "$<IF:$<CONFIG:Debug>,-DDEBUG=1,-DDEBUG=0>")

set(PYTHON_INTERPRETER "python3" CACHE STRING
  "Python interpreter to use for installing and test python components")


#
# Add optional external components
#

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Note: GIT_TAG hashes below must be updated with the utils/test/check_submodule_commits script!

if(WITH_MPI)
  set(MPIFX_GIT_REPOSITORY "https://github.com/dftbplus/mpifx.git")
  set(MPIFX_GIT_TAG "0cb07ee08cbb20f3f7bb2527152a4ec317c579ad")  # do not change manually!
  fnet_config_hybrid_dependency(MpiFx MpiFx::MpiFx "${HYBRID_CONFIG_METHODS}" "QUIET" 
    external/mpifx "${exclude}" "${MPIFX_GIT_REPOSITORY}" "${MPIFX_GIT_TAG}")
endif()

add_subdirectory(external/xmlf90 EXCLUDE_FROM_ALL)


#
# Add internal components
#
add_subdirectory(prog)


#
# Testing
#
string(CONFIGURE "${TEST_RUNNER_TEMPLATE}" TEST_RUNNER)
enable_testing()
add_subdirectory(test)
