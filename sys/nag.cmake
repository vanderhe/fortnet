#
# Toolchain file for
#
# NAG Fortran compiler, GNU C compiler
#
# Notes:
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually


#
# Fortran compiler settings
#
set(Fortran_FLAGS "-ieee=full -maxcontin=512 ${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O2"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "-g -f2008 -nan -C=all"
  CACHE STRING "Fortran compiler flags for Debug build")


#
# External libraries
#

# LAPACK and BLAS
# (if the BLAS library contains the LAPACK functions, set LAPACK_LIBRARY to "NONE")
#set(BLAS_LIBRARY "openblas" CACHE STRING "BLAS libraries to link")
#set(BLAS_LIBRARY_DIR "" CACHE STRING "Directories where BLAS libraries can be found")
#set(LAPACK_LIBRARY "NONE" CACHE STRING "LAPACK libraries to link")
#set(LAPACK_LIBRARY_DIR "" CACHE STRING "Directories where LAPACK libraries can be found")
