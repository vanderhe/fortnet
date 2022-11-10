cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(FortnetUtils)

fnet_ensure_out_of_source_build()
fnet_load_build_settings()

set(FORTNET_VERSION "0.7")

project(fortnet VERSION ${FORTNET_VERSION} LANGUAGES Fortran C)

fnet_setup_build_type()
fnet_load_toolchain_settings()
fnet_setup_global_compiler_flags()
fnet_ensure_config_consistency()
fnet_get_release_name(RELEASE)

find_package(HDF5 REQUIRED COMPONENTS Fortran_HL)

if(WITH_MPI)
  find_package(MPI REQUIRED)
  if(NOT MPI_FORTRAN_FOUND)
    message(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER} is not MPI capable but is specified "
      "for a WITH_MPI=TRUE build")
  endif()
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

find_package(CustomBlas REQUIRED)
list(APPEND PKG_CONFIG_LIBS_PRIVATE ${BLAS_LIBRARY})


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
if(BUILD_SHARED_LIBS)
  set(exclude)
  option(INSTALL_INCLUDE_FILES "Whether include files should be installed" FALSE)
else()
  set(exclude EXCLUDE_FROM_ALL)
  option(INSTALL_INCLUDE_FILES "Whether include files should be installed" FALSE)
endif()


# Follow GNU conventions for installing directories
include(GNUInstallDirs)

add_subdirectory(external/xmlf90 EXCLUDE_FROM_ALL)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets EXCLUDE_FROM_ALL)
endif()

# If INCLUDE_INDIRECT_DEPS is non-empty, indirect dependencies must also be explicitely treated
string(REGEX MATCH "(^|;)[Ss]ubmodule(^|;)" INCLUDE_INDIRECT_DEPS "${HYBRID_CONFIG_METHODS}")

# Note: GIT_TAG hashes below must be updated with the utils/test/check_submodule_commits script!

if(WITH_MPI)
  set(MPIFX_GIT_REPOSITORY "https://github.com/dftbplus/mpifx.git")
  set(MPIFX_GIT_TAG "da51073aa87831a91ae756a9773e39ea000d1c3a")  # do not change manually!
  fnet_config_hybrid_dependency(MpiFx MpiFx::MpiFx "${HYBRID_CONFIG_METHODS}" "QUIET"
    external/mpifx "${exclude}" "${MPIFX_GIT_REPOSITORY}" "${MPIFX_GIT_TAG}")
endif()


#
# Add internal components
#
add_subdirectory(prog)


#
# Testing
#
string(CONFIGURE "${TEST_RUNNER_TEMPLATE}" TEST_RUNNER)
enable_testing()
if(NOT BUILD_EXPORTED_TARGETS_ONLY)
  add_subdirectory(test)
endif()
