#
# Toolchain file for
#
# Intel compiler, MKL library
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


#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O1 -fp-model precise -ip -heap-arrays 10"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "-g -warn all -stand f08 -check -diag-error-limit 1 -traceback"
  CACHE STRING "Fortran compiler flags for Debug build")


#
# External libraries
#

# LAPACK and BLAS
# (if the BLAS library contains the LAPACK functions, set LAPACK_LIBRARY to "NONE")

# set(BLAS_LIBRARY "mkl_intel_lp64;mkl_sequential;mkl_core" CACHE STRING "BLAS libraries to link")
# set(BLAS_LIBRARY_DIR "$ENV{MKLROOT}/lib/intel64" CACHE STRING
#     "Directories where BLAS libraries can be found")

# set(LAPACK_LIBRARY_DIR "$ENV{MKLROOT}/lib/intel64" CACHE STRING
#    "Directories where LAPACK libraries can be found")
set(LAPACK_LIBRARY "NONE")
