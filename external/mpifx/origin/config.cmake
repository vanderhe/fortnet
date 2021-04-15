#
# Build options
#

# CMAKE_BUILD_TYPE is commented out in order to allow for multi-configuration builds. It will
# automatically default to RelWithDebInfo if used in a single configuration build. Uncomment or
# override it only if you want a non-default single configuration build.
#
#set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (Release|RelWithDebInfo|Debug|MinSizeRel)")

option(BUILD_SHARED_LIBS "Whether the library should be a shared one" FALSE)

#
# Installation options
#

option(INSTALL_INCLUDE_FILES "Whether include / module files should be installed" TRUE)

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_install" CACHE STRING
  "Directory to install the compiled code into")

set(INSTALL_INCLUDEDIR "mpifx" CACHE PATH
  "Installation directory for header and include files (within standard include folder)")

set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within standard include folder)")
