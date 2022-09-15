#
# Global architecture independent build settings
#

#set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (Release|RelWithDebInfo|Debug|MinSizeRel)")
# CMAKE_BUILD_TYPE is commented out in order to allow for multi-configuration builds. It will
# automatically default to RelWithDebInfo if used in a single configuration build. Uncomment or
# override it only if you want a non-default single configuration build.

option(WITH_MPI "Whether Fortnet should support MPI-parallelism" TRUE)

option(WITH_SOCKETS "Whether socket communication should be allowed for" TRUE)

option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" FALSE)
# Turn this on, if the Fortnet library (and other compiled libraries) should be shared libraries and
# dynamically linked to their applications. This results in smaller applications, but the libraries
# must be present at run-time (and the correct LD_LIBRARY_PATH environment variable must be set, so
# that they can be found by the operating system).


#
# Test environment settings
#
set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of MPI processes used for testing")

# Command line used to launch the test code.
# The escaped variables (\${VARIABLE}) will be substituted by the corresponding CMake variables.
if(WITH_MPI)
  set(TEST_RUNNER_TEMPLATE "mpiexec -n \${TEST_MPI_PROCS}"
    CACHE STRING "How to run the tests")
else()
  set(TEST_RUNNER_TEMPLATE " "
    CACHE STRING "How to run the tests")
endif()


#
# Installation options
#
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_install" CACHE STRING
  "Directory to install the compiled code into")

set(INSTALL_INCLUDEDIR "fortnet" CACHE PATH
  "Name of the project specific sub-folder within the install folder for include files")

set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within the install folder for include files)")

set(PKGCONFIG_LANGUAGE "Fortran" CACHE STRING
  "Compiler and Linker language to assume when creating the pkg-config export file (C or Fortran)")
# The pkg-config export file (lib/pkgconfig/dftbplus.pc) contains the compiler and linker options
# needed to link the Fortnet library to an application. (It can be queried with the pkg-config tool.)
# Depending on the language setting ("C" or "Fortran") you would get the flags for the case of using
# that compiler for the linking.


#
# Advanced options (e.g. for developers and packagers)
#

#set(TOOLCHAIN "gnu" CACHE STRING "Prefix of the toolchain file to be read from the sys/ folder")
# Uncomment and set it if you want to override the automatic, compiler based toolchain file
# selection.


set(HYBRID_CONFIG_METHODS "Submodule;Find;Fetch" CACHE STRING
  "Configuration methods to try in order to satisfy hybrid dependencies")
#
# This list can be used to control how hybrid dependencies (external dependencies which can
# optionally be built during the build process) are configured. The listed methods are applied in
# the following order:
#
# Submodule: Use the source in external/*/origin and build the dependency as part of the build
#     process. If the source is not present, try to retrieve it via the 'git submodule' command
#     (provided the source tree is a git repository and git is available)
#
# Find: Find the dependency as an already installed package in the system.
#
# Fetch: Fetch the source into the build folder and build the dependency as part of the build
#     process (works also in cases where the source tree is not a Git repository)
