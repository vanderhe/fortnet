file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../" projectdir)

#
# General options for all targets
#
set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${RELEASE}'")

#
# Compile and install library
#
set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)

add_subdirectory(lib_analysis)
add_subdirectory(lib_common)
add_subdirectory(lib_descriptors)
add_subdirectory(lib_dftbp)
add_subdirectory(lib_extlibs)
add_subdirectory(lib_fortnet)
add_subdirectory(lib_hdf5fx)
add_subdirectory(lib_io)
add_subdirectory(lib_nn)
add_subdirectory(lib_optimizers)
add_subdirectory(lib_types)
add_subdirectory(lib_utils)

fnet_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

#
# Library components
#
add_library(fortnet ${ALL-SOURCES-F90} ${all-sources-f90-preproc})
target_sources(fortnet PRIVATE $<TARGET_OBJECTS:xmlf90_objlib>)
set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(fortnet PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
target_include_directories(fortnet PUBLIC
  $<BUILD_INTERFACE:${includedir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}>)

#
# External libraries
#
target_include_directories(fortnet PUBLIC
  $<BUILD_INTERFACE:$<TARGET_PROPERTY:xmlf90_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)

if(WITH_MPI)
  target_link_libraries(fortnet PUBLIC MPI::MPI_Fortran)
  target_link_libraries(fortnet PUBLIC MpiFx::MpiFx)
endif()

# cmake 3.19.x and 3.20.x feature broken targets :(
# target_link_libraries(fortnet PUBLIC hdf5::hdf5_hl_fortran hdf5::hdf5_fortran)

target_link_libraries(fortnet PUBLIC ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_Fortran_LIBRARIES})
target_include_directories(fortnet PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
target_link_libraries(fortnet PRIVATE ${OTHER_LIBRARIES})
target_include_directories(fortnet PUBLIC ${OTHER_INCLUDE_DIRS})

#
# Installation
#
list(APPEND PKG_CONFIG_LIBS fortnet)
set(PKG_CONFIG_LIBS "${PKG_CONFIG_LIBS}" PARENT_SCOPE)

#
# Compile and install program
#
set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)

add_subdirectory(prg_fnet)

fnet_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}"
  all-sources-f90-preproc)

add_executable(fnet ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

target_link_libraries(fnet fortnet)

install(TARGETS fnet DESTINATION ${CMAKE_INSTALL_BINDIR})
